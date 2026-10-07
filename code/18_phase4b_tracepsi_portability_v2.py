#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE-Ψ Phase 4B v2 — weighted LOTO portability model
=====================================================

This is the optimized Phase4B model layer. It consumes the frozen Phase4A
benchmark and DOES NOT rebuild loci.

Why v2
------
Phase4A revealed two modeling details that should be fixed before formal
development evaluation:

1) source_breadth = source_BID + source_BACS + source_ELAP exactly.
   The old P1/P2 therefore contained exact linear redundancy.

2) The same locus contributes 2-3 target-task rows. Without weighting,
   broadly-supported loci receive more training weight simply because they
   generate more tasks.

v2 fixes both:
- ELAP is used as the provenance reference category:
      P1 = breadth + source_BID + source_BACS
- training weights make each locus contribute equal TOTAL weight.
- hyperparameter tuning maximizes macro target-specific AUROC, so a target
  task with more rows does not dominate tuning.
- inner CV remains group-disjoint by locus_id.

Scientific estimand
-------------------
Rank, among source-positive loci, which reported evidence is more likely to be
re-observed by another profiling technology.

This is NOT P(true Ψ).

Frozen model hierarchy
----------------------
P0_breadth:
    source_breadth

P1_provenance:
    source_breadth + source_BID + source_BACS
    (ELAP is the reference source; no exact collinearity)

A_context_only (secondary ablation only):
    P1 + additive flanking-base context

P2_TRACEpsi:
    P1
    + frozen cross-cell BID-vs-ELAP fingerprint
    + fingerprint magnitude
    + source/fingerprint alignment
    + fingerprint-seen flag
    + additive flanking-base context

At each flanking position, U is the reference base; U dummy columns are
omitted to avoid dummy-variable redundancy.

Validation
----------
A) HeLa held-target-role LOTO:
   BID / BACS / ELAP target roles are held out one at a time.
   This is NOT claimed as novel-locus external validation because loci can
   participate in other target roles. It tests target-role transport.

B) HEK293T cross-cell secondary benchmark:
   models fit only in HeLa, score HEK tasks.

C) HeLa DRS development benchmark:
   DRS is never used for fitting or hyperparameter tuning.

D) Primary utility challenge:
   HeLa DRS, source_breadth=1.

A549 DRS remains completely untouched for Phase4C lockbox validation.

Practical utility
-----------------
Besides AUROC/AP, report:
- Lift@10% / Lift@20%
- Precision among fixed experimental budgets K=50, 100, 200
- global and source-pattern-specific portability percentiles

Model selection for the future A549 lockbox is DEVELOPMENT selection, not a
claim of statistical truth. The later A549 lockbox is the decisive validation.
"""

from __future__ import annotations

import json
import math
import traceback
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import joblib
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.model_selection import GroupKFold
except Exception as e:
    raise RuntimeError("Phase4B v2 requires scikit-learn and joblib.") from e


ROOT = Path(r"D:\RNA\Trace")
TRACEPSI = ROOT / "07_tracepsi"
META = ROOT / "00_meta"
LOGS = ROOT / "logs"

for p in [TRACEPSI, META, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

LOTO = TRACEPSI / "phase4a_HeLa_LOTO.tsv.gz"
DRS_DEV = TRACEPSI / "phase4a_HeLa_DRS_development.tsv.gz"
HEK_SECONDARY = TRACEPSI / "phase4a_HEK_crosscell_tasks.tsv.gz"
PHASE4A_SUMMARY = META / "phase4a_summary.json"

SEED = 20261004
C_GRID = np.array([1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0], float)
N_BOOT = 3000

# Drop U at each position as reference category.
CONTEXT_FEATURES = [
    f"motif_{off:+d}_{base}"
    for off in [-2, -1, 1, 2]
    for base in ["A", "C", "G"]
]

P0 = ["source_breadth"]

# ELAP is reference source. Together with breadth this spans provenance
# without exact linear redundancy.
P1 = [
    "source_breadth",
    "source_BID",
    "source_BACS",
]

FINGERPRINT_FEATURES = [
    "fp_crosscell_ELAPvsBID",
    "fp_abs",
    "fp_source_alignment",
    "fp_seen_in_HEK",
]

P2 = P1 + FINGERPRINT_FEATURES + CONTEXT_FEATURES
A_CONTEXT = P1 + CONTEXT_FEATURES

PRIMARY_MODELS = {
    "P0_breadth": P0,
    "P1_provenance": P1,
    "P2_TRACEpsi": P2,
}
SECONDARY_ABLATIONS = {
    "A_context_only": A_CONTEXT,
}
ALL_MODELS = {**PRIMARY_MODELS, **SECONDARY_ABLATIONS}


# =============================================================================
# Helpers
# =============================================================================

def safe_auc(y, score) -> float:
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    if len(np.unique(y)) != 2:
        return np.nan
    return float(roc_auc_score(y, score))


def safe_ap(y, score) -> float:
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    if y.sum() == 0:
        return np.nan
    return float(average_precision_score(y, score))


def equal_locus_weights(df: pd.DataFrame) -> np.ndarray:
    """
    Every unique locus contributes equal TOTAL weight across its task rows.
    """
    counts = df.groupby("locus_id")["locus_id"].transform("size").astype(float)
    w = 1.0 / counts.to_numpy()
    return w / np.mean(w)


def weighted_standardizer(X: np.ndarray, w: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mu = np.average(X, axis=0, weights=w)
    var = np.average((X - mu) ** 2, axis=0, weights=w)
    sd = np.sqrt(np.maximum(var, 0.0))
    sd = np.where(sd > 1e-12, sd, 1.0)
    return mu, sd


def fit_weighted_logit(
    train: pd.DataFrame,
    features: Sequence[str],
    C: float,
) -> Dict:
    X = train[list(features)].to_numpy(float)
    y = train["target_reported"].astype(int).to_numpy()
    w = equal_locus_weights(train)

    mu, sd = weighted_standardizer(X, w)
    Xs = (X - mu) / sd

    model = LogisticRegression(
        C=float(C),
        penalty="l2",
        solver="lbfgs",
        max_iter=5000,
        random_state=SEED,
    )
    model.fit(Xs, y, sample_weight=w)

    return {
        "features": list(features),
        "C": float(C),
        "mean": mu,
        "scale": sd,
        "model": model,
    }


def decision_function(bundle: Dict, df: pd.DataFrame) -> np.ndarray:
    X = df[bundle["features"]].to_numpy(float)
    Xs = (X - bundle["mean"]) / bundle["scale"]
    return bundle["model"].decision_function(Xs)


def grouped_splits(df: pd.DataFrame, n_splits: int = 5):
    n_groups = df["locus_id"].nunique()
    n_splits = min(n_splits, n_groups)
    if n_splits < 3:
        raise RuntimeError("Insufficient locus groups for tuning.")
    splitter = GroupKFold(n_splits=n_splits)
    return list(
        splitter.split(
            df,
            df["target_reported"].astype(int),
            groups=df["locus_id"].astype(str),
        )
    )


def macro_target_auc(df: pd.DataFrame, score: np.ndarray) -> float:
    vals = []
    score = np.asarray(score, float)

    for target, idx_labels in df.groupby("target_assay").groups.items():
        # Convert original index labels to positional indices.
        pos = df.index.get_indexer(idx_labels)
        y = df.iloc[pos]["target_reported"].astype(int).to_numpy()
        if len(np.unique(y)) == 2:
            vals.append(safe_auc(y, score[pos]))

    return float(np.mean(vals)) if vals else np.nan


def tune_C(
    train: pd.DataFrame,
    features: Sequence[str],
) -> Tuple[float, pd.DataFrame]:
    train = train.reset_index(drop=True)
    splits = grouped_splits(train, 5)
    rows = []

    for C in C_GRID:
        fold_scores = []

        for tr_idx, va_idx in splits:
            tr = train.iloc[tr_idx].copy()
            va = train.iloc[va_idx].copy()

            if tr["target_reported"].nunique() < 2:
                continue

            bundle = fit_weighted_logit(tr, features, float(C))
            score = decision_function(bundle, va)
            auc = macro_target_auc(va.reset_index(drop=True), score)

            if np.isfinite(auc):
                fold_scores.append(auc)

        rows.append({
            "C": float(C),
            "mean_macro_target_AUROC": (
                float(np.mean(fold_scores)) if fold_scores else np.nan
            ),
            "folds_used": int(len(fold_scores)),
        })

    tab = pd.DataFrame(rows)
    valid = tab.dropna(subset=["mean_macro_target_AUROC"]).copy()
    if len(valid) == 0:
        raise RuntimeError("No valid grouped tuning folds.")

    best = valid["mean_macro_target_AUROC"].max()
    tied = valid[
        valid["mean_macro_target_AUROC"] >= best - 1e-10
    ].sort_values("C")

    # Tie -> stronger regularization.
    best_C = float(tied.iloc[0]["C"])
    tab["selected"] = tab["C"].eq(best_C)
    return best_C, tab


def fit_tuned(
    train: pd.DataFrame,
    features: Sequence[str],
) -> Tuple[Dict, pd.DataFrame]:
    C, tuning = tune_C(train, features)
    bundle = fit_weighted_logit(train, features, C)
    return bundle, tuning


# =============================================================================
# Metrics
# =============================================================================

def bootstrap_auc_ci(y, score, seed: int) -> Tuple[float, float, int]:
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    rng = np.random.default_rng(seed)
    vals = []

    for _ in range(N_BOOT):
        idx = rng.integers(0, len(y), size=len(y))
        yy = y[idx]
        if len(np.unique(yy)) != 2:
            continue
        vals.append(safe_auc(yy, score[idx]))

    if not vals:
        return np.nan, np.nan, 0
    return (
        float(np.quantile(vals, .025)),
        float(np.quantile(vals, .975)),
        len(vals),
    )


def top_fraction(y, score, fraction: float) -> Dict[str, float]:
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    k = max(1, int(math.ceil(len(y) * fraction)))
    order = np.argsort(-score, kind="mergesort")
    precision = float(y[order[:k]].mean())
    base = float(y.mean())
    return {
        "k": k,
        "precision": precision,
        "lift": precision / base if base > 0 else np.nan,
    }


def top_k(y, score, k: int) -> Dict[str, float]:
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    k = min(int(k), len(y))
    order = np.argsort(-score, kind="mergesort")
    precision = float(y[order[:k]].mean())
    base = float(y.mean())
    return {
        "k": k,
        "positives_in_topK": int(y[order[:k]].sum()),
        "precision": precision,
        "lift": precision / base if base > 0 else np.nan,
    }


def bootstrap_lift_ci(y, score, fraction: float, seed: int) -> Tuple[float, float]:
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    rng = np.random.default_rng(seed)
    vals = []

    for _ in range(N_BOOT):
        idx = rng.integers(0, len(y), size=len(y))
        yy, ss = y[idx], score[idx]
        if yy.sum() == 0:
            continue
        v = top_fraction(yy, ss, fraction)["lift"]
        if np.isfinite(v):
            vals.append(v)

    if not vals:
        return np.nan, np.nan
    return (
        float(np.quantile(vals, .025)),
        float(np.quantile(vals, .975)),
    )


def evaluate(
    df: pd.DataFrame,
    score: np.ndarray,
    label: str,
    seed: int,
) -> Dict:
    y = df["target_reported"].astype(int).to_numpy()
    auc = safe_auc(y, score)
    ap = safe_ap(y, score)
    lo, hi, nb = bootstrap_auc_ci(y, score, seed)

    f10 = top_fraction(y, score, .10)
    f20 = top_fraction(y, score, .20)
    l10lo, l10hi = bootstrap_lift_ci(y, score, .10, seed + 1)
    l20lo, l20hi = bootstrap_lift_ci(y, score, .20, seed + 2)

    rec = {
        "analysis": label,
        "n": int(len(y)),
        "positives": int(y.sum()),
        "prevalence": float(y.mean()),
        "AUROC": auc,
        "AUROC_CI_low": lo,
        "AUROC_CI_high": hi,
        "average_precision": ap,
        "Lift10": f10["lift"],
        "Precision10": f10["precision"],
        "K10": f10["k"],
        "Lift10_CI_low": l10lo,
        "Lift10_CI_high": l10hi,
        "Lift20": f20["lift"],
        "Precision20": f20["precision"],
        "K20": f20["k"],
        "Lift20_CI_low": l20lo,
        "Lift20_CI_high": l20hi,
        "bootstrap_replicates": nb,
    }

    for k in [50, 100, 200]:
        if len(y) >= k:
            x = top_k(y, score, k)
            rec[f"top{k}_positives"] = x["positives_in_topK"]
            rec[f"top{k}_precision"] = x["precision"]
            rec[f"top{k}_lift"] = x["lift"]

    return rec


def paired_delta_auc(
    y,
    score_a,
    score_b,
    seed: int,
) -> Dict:
    y = np.asarray(y, int)
    a = np.asarray(score_a, float)
    b = np.asarray(score_b, float)
    observed = safe_auc(y, a) - safe_auc(y, b)

    rng = np.random.default_rng(seed)
    vals = []

    for _ in range(N_BOOT):
        idx = rng.integers(0, len(y), size=len(y))
        yy = y[idx]
        if len(np.unique(yy)) != 2:
            continue
        vals.append(
            safe_auc(yy, a[idx]) - safe_auc(yy, b[idx])
        )

    vals = np.asarray(vals, float)
    return {
        "delta_AUROC": float(observed),
        "CI_low": float(np.quantile(vals, .025)),
        "CI_high": float(np.quantile(vals, .975)),
        "bootstrap_P_delta_gt_0": float(np.mean(vals > 0)),
        "bootstrap_replicates": int(len(vals)),
    }


# =============================================================================
# Diagnostics / scoring surfaces
# =============================================================================

def provenance_rates(df: pd.DataFrame) -> pd.DataFrame:
    z = df[df["source_breadth"].eq(1)].copy()
    return (
        z.groupby(["target_assay", "source_pattern"], as_index=False)
        .agg(
            n=("locus_id", "size"),
            target_reported_n=("target_reported", "sum"),
            reobservation_rate=("target_reported", "mean"),
        )
        .sort_values(["target_assay", "reobservation_rate"], ascending=[True, False])
    )


def alignment_diagnostic(df: pd.DataFrame) -> pd.DataFrame:
    """
    Predefined Aim1->portability bridge:
    higher positive source alignment was hypothesized to reflect greater
    source-assay dependence; use -alignment as portability score.
    """
    rows = []
    z0 = df[df["source_breadth"].eq(1)].copy()

    for i, (target, z) in enumerate(z0.groupby("target_assay")):
        score = -z["fp_source_alignment"].to_numpy(float)
        rec = evaluate(
            z,
            score,
            label=f"breadth1_minus_source_alignment_{target}",
            seed=SEED + 9000 + i,
        )
        rec["target_assay"] = target
        rows.append(rec)

    return pd.DataFrame(rows)


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 108)
    print("TRACE-Ψ PHASE 4B v2 — WEIGHTED LOTO PORTABILITY MODEL")
    print("=" * 108)

    for p in [LOTO, DRS_DEV, HEK_SECONDARY, PHASE4A_SUMMARY]:
        if not p.exists():
            raise RuntimeError(f"Missing Phase4A input: {p}")

    loto = pd.read_csv(LOTO, sep="\t", compression="gzip", low_memory=False)
    drs = pd.read_csv(DRS_DEV, sep="\t", compression="gzip", low_memory=False)
    hek = pd.read_csv(HEK_SECONDARY, sep="\t", compression="gzip", low_memory=False)

    with open(PHASE4A_SUMMARY, "r", encoding="utf-8") as f:
        s4a = json.load(f)

    if s4a.get("status") != "PHASE4A_PORTABILITY_BENCHMARK_FROZEN":
        raise RuntimeError("Phase4A is not frozen/ready.")

    for name, df in [("LOTO", loto), ("DRS", drs), ("HEK", hek)]:
        required = {
            "locus_id", "target_assay", "target_reported",
            "source_breadth", "source_pattern",
            *P2,
        }
        missing = sorted(required - set(df.columns))
        if missing:
            raise RuntimeError(f"{name}: missing required columns {missing}")
        if df[list(P2)].isna().any().any():
            raise RuntimeError(f"{name}: model features contain missing values.")

    # Exact redundancy hard gate for chosen model matrix.
    if not np.array_equal(
        loto["source_breadth"].astype(int).to_numpy(),
        (
            loto["source_BID"].astype(int)
            + loto["source_BACS"].astype(int)
            + loto["source_ELAP"].astype(int)
        ).to_numpy(),
    ):
        raise RuntimeError("source_breadth no longer equals source support sum.")

    # -------------------------------------------------------------------------
    # Transparent 4A-derived diagnostics
    # -------------------------------------------------------------------------
    prov = provenance_rates(loto)
    prov.to_csv(
        TRACEPSI / "phase4b_v2_breadth1_provenance_rates.tsv",
        sep="\t", index=False
    )

    align = alignment_diagnostic(loto)
    align.to_csv(
        TRACEPSI / "phase4b_v2_breadth1_alignment_diagnostic.tsv",
        sep="\t", index=False
    )

    # -------------------------------------------------------------------------
    # A. Held-target-role HeLa LOTO
    # -------------------------------------------------------------------------
    internal_metrics = []
    internal_predictions = []
    internal_tuning = []

    targets = sorted(loto["target_assay"].unique())
    if targets != ["BACS", "BID", "ELAP"]:
        raise RuntimeError(f"Unexpected target assays: {targets}")

    for ti, held in enumerate(targets):
        train = loto[loto["target_assay"].ne(held)].copy().reset_index(drop=True)
        test = loto[loto["target_assay"].eq(held)].copy().reset_index(drop=True)

        for mi, (model_name, features) in enumerate(ALL_MODELS.items()):
            bundle, tuning = fit_tuned(train, features)

            tuning["held_target"] = held
            tuning["model"] = model_name
            internal_tuning.append(tuning)

            score = decision_function(bundle, test)
            rec = evaluate(
                test,
                score,
                label=f"HeLa_held_target_{held}_{model_name}",
                seed=SEED + 100*ti + mi,
            )
            rec.update({
                "held_target": held,
                "model": model_name,
                "selected_C": bundle["C"],
                "train_rows": int(len(train)),
                "train_loci": int(train["locus_id"].nunique()),
                "test_rows": int(len(test)),
            })
            internal_metrics.append(rec)

            q = test[
                [
                    "locus_id", "target_assay", "target_reported",
                    "source_breadth", "source_pattern", "motif",
                ]
            ].copy()
            q["model"] = model_name
            q["score"] = score
            internal_predictions.append(q)

    internal_metrics = pd.DataFrame(internal_metrics)
    internal_metrics.to_csv(
        TRACEPSI / "phase4b_v2_internal_metrics.tsv",
        sep="\t", index=False
    )
    pd.concat(internal_tuning, ignore_index=True).to_csv(
        TRACEPSI / "phase4b_v2_internal_tuning.tsv",
        sep="\t", index=False
    )
    pd.concat(internal_predictions, ignore_index=True).to_csv(
        TRACEPSI / "phase4b_v2_internal_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    internal_macro = (
        internal_metrics[
            internal_metrics["model"].isin(PRIMARY_MODELS)
        ]
        .groupby("model", as_index=False)
        .agg(
            macro_AUROC=("AUROC", "mean"),
            macro_AP=("average_precision", "mean"),
            macro_Lift10=("Lift10", "mean"),
            min_target_AUROC=("AUROC", "min"),
        )
    )

    # -------------------------------------------------------------------------
    # B. Fit all HeLa LOTO development models (still without DRS labels)
    # -------------------------------------------------------------------------
    fitted = {}
    fitall_tuning = []
    coef_rows = []

    for model_name, features in ALL_MODELS.items():
        bundle, tuning = fit_tuned(loto.reset_index(drop=True), features)
        fitted[model_name] = bundle

        tuning["model"] = model_name
        fitall_tuning.append(tuning)

        coef = bundle["model"].coef_[0]
        for feature, b, mu, sd in zip(
            features, coef, bundle["mean"], bundle["scale"]
        ):
            coef_rows.append({
                "model": model_name,
                "feature": feature,
                "standardized_coefficient": float(b),
                "training_weighted_mean": float(mu),
                "training_weighted_scale": float(sd),
                "selected_C": float(bundle["C"]),
            })

    pd.concat(fitall_tuning, ignore_index=True).to_csv(
        TRACEPSI / "phase4b_v2_fitall_tuning.tsv",
        sep="\t", index=False
    )
    pd.DataFrame(coef_rows).to_csv(
        TRACEPSI / "phase4b_v2_coefficients.tsv",
        sep="\t", index=False
    )

    # -------------------------------------------------------------------------
    # C. HEK cross-cell secondary
    # -------------------------------------------------------------------------
    hek_rows = []
    hek_pred = hek[
        [
            "locus_id", "target_assay", "target_reported",
            "source_breadth", "source_pattern", "motif",
        ]
    ].copy()

    for mi, (model_name, bundle) in enumerate(fitted.items()):
        score = decision_function(bundle, hek)
        hek_pred[f"score_{model_name}"] = score

        # Combined
        if hek["target_reported"].nunique() == 2:
            rec = evaluate(
                hek, score,
                label=f"HEK_all_{model_name}",
                seed=SEED + 2000 + mi,
            )
            rec.update({"model": model_name, "target": "combined"})
            hek_rows.append(rec)

        # Per target
        for tj, target in enumerate(sorted(hek["target_assay"].unique())):
            mask = hek["target_assay"].eq(target).to_numpy()
            z = hek.loc[mask].copy()
            if z["target_reported"].nunique() != 2:
                continue
            rec = evaluate(
                z, score[mask],
                label=f"HEK_{target}_{model_name}",
                seed=SEED + 2100 + 10*mi + tj,
            )
            rec.update({"model": model_name, "target": target})
            hek_rows.append(rec)

    pd.DataFrame(hek_rows).to_csv(
        TRACEPSI / "phase4b_v2_HEK_crosscell_metrics.tsv",
        sep="\t", index=False
    )
    hek_pred.to_csv(
        TRACEPSI / "phase4b_v2_HEK_crosscell_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # -------------------------------------------------------------------------
    # D. HeLa DRS development benchmark
    # -------------------------------------------------------------------------
    drs_rows = []
    drs_pred = drs[
        [
            "locus_id", "target_reported", "source_breadth",
            "source_pattern", "motif",
            "source_BID", "source_BACS", "source_ELAP",
        ]
    ].copy()

    scores = {}
    for mi, (model_name, bundle) in enumerate(fitted.items()):
        score = decision_function(bundle, drs)
        scores[model_name] = score

        drs_pred[f"score_{model_name}"] = score
        drs_pred[f"percentile_{model_name}"] = (
            pd.Series(score)
            .rank(method="average", pct=True)
            .to_numpy()
        )

        # Global
        rec = evaluate(
            drs, score,
            label=f"HeLa_DRS_all_{model_name}",
            seed=SEED + 3000 + mi,
        )
        rec.update({"model": model_name, "subset": "all"})
        drs_rows.append(rec)

        # Core breadth=1 challenge
        mask = drs["source_breadth"].eq(1).to_numpy()
        z = drs.loc[mask].copy()
        rec = evaluate(
            z, score[mask],
            label=f"HeLa_DRS_breadth1_{model_name}",
            seed=SEED + 3100 + mi,
        )
        rec.update({"model": model_name, "subset": "breadth1"})
        drs_rows.append(rec)

    drs_metrics = pd.DataFrame(drs_rows)
    drs_metrics.to_csv(
        TRACEPSI / "phase4b_v2_HeLa_DRS_metrics.tsv",
        sep="\t", index=False
    )
    drs_pred.to_csv(
        TRACEPSI / "phase4b_v2_HeLa_DRS_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # -------------------------------------------------------------------------
    # E. Paired breadth=1 comparisons
    # -------------------------------------------------------------------------
    b1 = drs["source_breadth"].eq(1).to_numpy()
    y = drs.loc[b1, "target_reported"].astype(int).to_numpy()

    comparisons = []
    for ai, (a, b) in enumerate([
        ("P1_provenance", "P0_breadth"),
        ("A_context_only", "P1_provenance"),
        ("P2_TRACEpsi", "P1_provenance"),
        ("P2_TRACEpsi", "A_context_only"),
    ]):
        rec = paired_delta_auc(
            y,
            scores[a][b1],
            scores[b][b1],
            seed=SEED + 4000 + ai,
        )
        rec.update({
            "subset": "HeLa_DRS_breadth1",
            "model_A": a,
            "model_B": b,
        })
        comparisons.append(rec)

    comparisons = pd.DataFrame(comparisons)
    comparisons.to_csv(
        TRACEPSI / "phase4b_v2_HeLa_DRS_breadth1_deltaAUC.tsv",
        sep="\t", index=False
    )

    # -------------------------------------------------------------------------
    # F. Development model selection for the future A549 LOCKBOX
    # -------------------------------------------------------------------------
    b1m = (
        drs_metrics[drs_metrics["subset"].eq("breadth1")]
        .set_index("model")
    )
    p1 = b1m.loc["P1_provenance"]
    p2 = b1m.loc["P2_TRACEpsi"]

    d21 = comparisons[
        (comparisons["model_A"] == "P2_TRACEpsi")
        & (comparisons["model_B"] == "P1_provenance")
    ].iloc[0]

    # Development thresholds are deliberately practical rather than framed as
    # definitive hypothesis tests; A549 is the untouched confirmatory lockbox.
    p2_minimum_utility = bool(
        p2["AUROC"] >= 0.58
        and p2["AUROC_CI_low"] > 0.50
        and p2["Lift10"] >= 1.15
    )
    p2_clear_gain = bool(
        d21["delta_AUROC"] >= 0.02
        and d21["bootstrap_P_delta_gt_0"] >= 0.80
    )

    p1_minimum_utility = bool(
        p1["AUROC"] >= 0.58
        and p1["AUROC_CI_low"] > 0.50
        and p1["Lift10"] >= 1.15
    )

    if p2_minimum_utility and p2_clear_gain:
        status = "TRACEPSI_P2_READY_FOR_A549_LOCKBOX"
        selected = "P2_TRACEpsi"
    elif p2_minimum_utility:
        status = "TRACEPSI_P2_USEFUL_BUT_GAIN_OVER_PROVENANCE_MODEST"
        selected = "P2_TRACEpsi_qualified"
    elif p1_minimum_utility:
        status = "TRACEPSI_PROVENANCE_ONLY_SIGNAL"
        selected = "P1_provenance"
    else:
        status = "TRACEPSI_NO_USEFUL_BREADTH1_RANKING_YET"
        selected = "none"

    # -------------------------------------------------------------------------
    # G. Prototype tool scorecard using selected P2 scores (always output P2;
    # selection status controls claims, not file availability)
    # -------------------------------------------------------------------------
    scorecard = drs_pred.copy()
    scorecard["TRACEpsi_score"] = scorecard["score_P2_TRACEpsi"]
    scorecard["TRACEpsi_percentile_global"] = (
        scorecard["TRACEpsi_score"]
        .rank(method="average", pct=True)
    )
    scorecard["TRACEpsi_percentile_within_source_pattern"] = (
        scorecard.groupby("source_pattern")["TRACEpsi_score"]
        .rank(method="average", pct=True)
    )
    scorecard.to_csv(
        TRACEPSI / "phase4b_v2_TRACEpsi_HeLa_scorecard.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # Freeze fitted models BEFORE opening A549.
    model_artifact = {
        "phase": "4B_v2_development",
        "selected_complexity_for_lockbox": selected,
        "status": status,
        "primary_models": PRIMARY_MODELS,
        "secondary_ablation": SECONDARY_ABLATIONS,
        "fitted_models": fitted,
        "score_semantics": (
            "relative cross-platform evidence-portability score; not P(true Ψ)"
        ),
        "ELAP_reference_for_provenance": True,
        "equal_total_training_weight_per_locus": True,
        "A549_DRS_seen": False,
    }
    joblib.dump(
        model_artifact,
        TRACEPSI / "phase4b_v2_TRACEpsi_development_model.joblib"
    )

    model_contract = {
        "phase": "4B_v2",
        "purpose": "TRACE-Ψ development model selection before A549 lockbox",
        "estimand": (
            "Ranking of source-positive loci by re-observation propensity on "
            "another profiling technology; not biological truth probability."
        ),
        "models": {
            **PRIMARY_MODELS,
            **SECONDARY_ABLATIONS,
        },
        "redundancy_fix": (
            "ELAP is provenance reference; U is flanking-base reference."
        ),
        "training_weight": (
            "Each unique locus contributes equal total weight across its target-task rows."
        ),
        "hyperparameter_selection": (
            "5-fold GroupKFold by locus_id; select C by mean target-specific AUROC "
            "across validation targets."
        ),
        "HeLa_internal_interpretation": (
            "held-target-role transport, not novel-locus external validation"
        ),
        "DRS_policy": (
            "HeLa DRS used only after models are fitted/tuned on HeLa LOTO; "
            "DRS labels do not affect model coefficients or C."
        ),
        "development_selection_rule": {
            "P2_minimum_utility": (
                "breadth1 DRS AUROC>=0.58, bootstrap lower CI>0.50, Lift10>=1.15"
            ),
            "P2_clear_gain_over_P1": (
                "delta AUROC>=0.02 and bootstrap P(delta>0)>=0.80"
            ),
            "note": (
                "These are practical development criteria; independent A549 lockbox "
                "is the final validation."
            ),
        },
        "A549_DRS_lockbox": {
            "status": "UNTOUCHED",
            "read": False,
            "used_for_selection": False,
        },
    }
    (META / "phase4b_v2_analysis_contract.json").write_text(
        json.dumps(model_contract, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    summary = {
        "status": status,
        "selected_complexity_for_A549_lockbox": selected,
        "HeLa_internal_macro": internal_macro.to_dict(orient="records"),
        "HeLa_internal_per_target": internal_metrics.to_dict(orient="records"),
        "breadth1_provenance_rates": prov.to_dict(orient="records"),
        "breadth1_alignment_diagnostic": align.to_dict(orient="records"),
        "HeLa_DRS_metrics": drs_metrics.to_dict(orient="records"),
        "HeLa_DRS_breadth1_delta_AUC": comparisons.to_dict(orient="records"),
        "P2_minimum_utility": p2_minimum_utility,
        "P2_clear_gain_over_P1": p2_clear_gain,
        "P1_minimum_utility": p1_minimum_utility,
        "prototype_scorecard": str(
            TRACEPSI / "phase4b_v2_TRACEpsi_HeLa_scorecard.tsv.gz"
        ),
        "model_artifact": str(
            TRACEPSI / "phase4b_v2_TRACEpsi_development_model.joblib"
        ),
        "next_gate": (
            "If P2 is selected, freeze it and run one A549 BID->DRS lockbox "
            "without changing features/hyperparameters. If P1 only is selected, "
            "do not pretend it can rank within a BID-only A549 list; reconsider "
            "whether source-strength should become a separately justified v2 feature."
        ),
    }
    (META / "phase4b_v2_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )

    print()
    print("STATUS:", status)
    print("Selected complexity:", selected)
    print()
    print("HELA INTERNAL MACRO")
    print(internal_macro.to_string(index=False))
    print()
    print("BREADTH=1 PROVENANCE RATES")
    print(prov.to_string(index=False))
    print()
    print("AIM1->PORTABILITY ALIGNMENT DIAGNOSTIC")
    print(
        align[
            [
                "target_assay", "n", "positives",
                "AUROC", "AUROC_CI_low", "AUROC_CI_high",
                "Lift10", "Lift10_CI_low", "Lift10_CI_high",
            ]
        ].to_string(index=False)
    )
    print()
    print("HELA DRS BREADTH=1")
    print(
        b1m[
            [
                "n", "positives", "AUROC", "AUROC_CI_low", "AUROC_CI_high",
                "average_precision", "Lift10", "Precision10",
                "top50_positives", "top100_positives", "top200_positives",
            ]
        ].to_string()
    )
    print()
    print("BREADTH=1 PAIRED DELTA AUC")
    print(comparisons.to_string(index=False))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase4b_v2_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase4b_v2_analysis_contract.json")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4b_v2_internal_metrics.tsv")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4b_v2_HeLa_DRS_metrics.tsv")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4b_v2_HeLa_DRS_breadth1_deltaAUC.tsv")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4b_v2_TRACEpsi_HeLa_scorecard.tsv.gz")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase4b_v2_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase4b_v2_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
