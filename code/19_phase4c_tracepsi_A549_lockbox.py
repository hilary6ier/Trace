#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE-Ψ Phase 4C — A549 BID -> DRS external lockbox
===================================================

This phase opens the previously untouched A549 cell+platform lockbox exactly
once. It DOES NOT fit, tune, recalibrate, or alter TRACE-Ψ.

Frozen source list
------------------
A549 BID-seq locally-single-U exact calls.

Frozen target
-------------
A549 native direct-RNA (DRS/Mod-p ID) positive calls from the already-normalized
Supplementary Table S6.

Frozen model
------------
Phase4B-v2 P2_TRACEpsi, trained/tuned only on HeLa LOTO development tasks.
The frozen A_context_only model and the predefined -source-alignment score are
evaluated only as secondary ablations; they cannot replace P2 after seeing A549.

Primary scientific question
---------------------------
Within a BID-only source list (all loci have breadth=1 and identical source
provenance), can the frozen technology-aware TRACE-Ψ score rank loci by
re-observation on an unseen native-RNA platform in another cell context?

Outcome semantics
-----------------
DRS_reported=1 means exact coordinate + same 5-mer was reported by the A549
DRS/Mod-p ID published call table.
DRS_reported=0 means not reported by that published target table, NOT
biological absence / false Ψ.

Matching
--------
Same conservative rule as frozen Aim2B:
- source = locally-single-U exact BID locus
- target = A549 DRS psi-positive call
- match requires exact chromosome + genomic position + identical 5-mer
- no +/-k window
- source/target motif conflicts are audited and excluded rather than labeled 0

No raw sequencing. No threshold optimization.

Primary metrics
---------------
- AUROC with bootstrap CI
- average precision
- tie-aware Lift@10% and Lift@20%
- expected precision/lift at budgets K=50/100/200
- score-quintile re-observation gradient

Lockbox interpretation (frozen before reading A549)
----------------------------------------------------
STRONG external validation:
    positives >=15
    P2 AUROC >=0.60
    bootstrap 95% CI lower >0.50
    tie-aware Lift@20% >=1.20

SUPPORTIVE external validation:
    positives >=10
    P2 AUROC >=0.55
    bootstrap support fraction P*(AUROC>0.5) >=0.90
    and (Lift@10% >=1.15 OR Lift@20% >=1.15)

LOW POWER:
    positives <10

A549 is not used to choose a new model. If P2 fails, report domain-shift
failure/qualification; do not modify P2 and re-open the same lockbox.
"""

from __future__ import annotations

import importlib.util
import json
import math
import re
import traceback
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import joblib
    from sklearn.metrics import average_precision_score, roc_auc_score
except Exception as e:
    raise RuntimeError("Phase4C requires scikit-learn and joblib.") from e


# =============================================================================
# Paths
# =============================================================================

ROOT = Path(r"D:\RNA\Trace")
CODE = ROOT / "code"
NORM = ROOT / "03_harmonized" / "normalized_sources"
TRACEPSI = ROOT / "07_tracepsi"
META = ROOT / "00_meta"
LOGS = ROOT / "logs"

for p in [TRACEPSI, META, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

PHASE4A_CODE = CODE / "16_phase4a_tracepsi_benchmark.py"
BID_A549 = NORM / "BID_A549.tsv"
DRS = NORM / "DRS_S6_psi_calls.tsv"

FP_HEK = TRACEPSI / "phase4a_fingerprint_HEK_for_HeLa.tsv"
MODEL_ARTIFACT = TRACEPSI / "phase4b_v2_TRACEpsi_development_model.joblib"
PHASE4B_SUMMARY = META / "phase4b_v2_summary.json"
PHASE4B_CONTRACT = META / "phase4b_v2_analysis_contract.json"

SEED = 20261004
N_BOOT = 5000


# =============================================================================
# Frozen Phase4A helper import
# =============================================================================

def load_phase4a_module():
    if not PHASE4A_CODE.exists():
        raise RuntimeError(
            f"Missing frozen Phase4A code: {PHASE4A_CODE}"
        )
    spec = importlib.util.spec_from_file_location(
        "tracepsi_phase4a_frozen", PHASE4A_CODE
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not import frozen Phase4A helper module.")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# =============================================================================
# Generic helpers
# =============================================================================

def ss(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    return str(x).strip()


def truthy(v: Any) -> bool:
    z = ss(v).lower()
    if z in {"", "0", "0.0", "false", "no", "na", "n/a", "nan", "none"}:
        return False
    try:
        return float(z) > 0
    except Exception:
        return True


def norm_chr(x: Any) -> str:
    z = ss(x)
    if not z:
        return ""
    if z.lower().startswith("chr"):
        tail = z[3:]
        return "chr" + ("M" if tail.upper() == "MT" else tail)
    if re.fullmatch(r"(?:\d+|X|Y|M|MT)", z, flags=re.I):
        return "chr" + ("M" if z.upper() == "MT" else z.upper())
    return z


def motif5(x: Any) -> str:
    z = ss(x).upper().replace("Ψ", "U").replace("T", "U")
    z = re.sub(r"[^ACGU]", "", z)
    if len(z) == 5 and z[2] == "U":
        return z
    return ""


def locally_single_u(m: str) -> bool:
    return bool(
        len(m) == 5
        and m[2] == "U"
        and m[1] != "U"
        and m[3] != "U"
    )


# =============================================================================
# DRS A549 parser
# =============================================================================

def find_a549_rows(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict]:
    # Preferred: normalized source_sheet explicitly tells which original sheet.
    if "source_sheet" in df.columns:
        mask = (
            df["source_sheet"].astype(str)
            .str.contains(r"(?i)\bA549\b", regex=True, na=False)
        )
        if mask.sum() > 0:
            return df.loc[mask].copy(), {
                "mode": "source_sheet",
                "column": "source_sheet",
                "rows": int(mask.sum()),
            }

    # Conservative fallback: column values explicitly containing A549.
    candidates = []
    for c in df.columns:
        if c == "source_file":
            continue
        vals = df[c].astype(str)
        mask = vals.str.contains(r"(?i)\bA549\b", regex=True, na=False)
        n = int(mask.sum())
        if n > 0:
            candidates.append((n, c, mask))

    if not candidates:
        raise RuntimeError(
            "Could not identify A549 rows in normalized DRS table."
        )

    candidates.sort(reverse=True, key=lambda x: x[0])
    n, c, mask = candidates[0]
    return df.loc[mask].copy(), {
        "mode": "value_search",
        "column": c,
        "rows": int(n),
    }


def parse_a549_drs() -> Tuple[pd.DataFrame, Dict]:
    raw = pd.read_csv(
        DRS,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )
    a549, row_report = find_a549_rows(raw)

    required = {
        "raw__chr",
        "raw__position",
        "raw__kmer",
        "raw__psi",
    }
    missing = sorted(required - set(a549.columns))
    if missing:
        raise RuntimeError(
            f"A549 DRS standardized schema missing {missing}. "
            f"Available columns={list(a549.columns)}"
        )

    x = pd.DataFrame({
        "chrom": a549["raw__chr"].map(norm_chr),
        "pos1": pd.to_numeric(a549["raw__position"], errors="coerce"),
        "motif": a549["raw__kmer"].map(motif5),
        "psi_detected": a549["raw__psi"].map(truthy).astype(int),
    })

    x = x[
        x["chrom"].ne("")
        & x["pos1"].notna()
        & x["motif"].ne("")
    ].copy()
    x["pos1"] = x["pos1"].astype(int)

    x["single_u"] = x["motif"].map(locally_single_u)

    # Positive-call target set only.
    pos = x[
        x["single_u"] & x["psi_detected"].eq(1)
    ].copy()

    pos["match_key"] = (
        pos["chrom"]
        + ":"
        + pos["pos1"].astype(str)
        + ":"
        + pos["motif"]
    )

    pos = pos.drop_duplicates("match_key").copy()

    report = {
        "normalized_DRS_rows_total": int(len(raw)),
        "A549_row_detection": row_report,
        "A549_rows_schema_resolved": int(len(x)),
        "A549_singleU_rows": int(x["single_u"].sum()),
        "A549_singleU_psi_positive_unique_pairs": int(len(pos)),
        "schema": {
            "chr": "raw__chr",
            "position": "raw__position",
            "kmer": "raw__kmer",
            "psi": "raw__psi",
        },
    }
    return pos, report


# =============================================================================
# Model scoring
# =============================================================================

def decision_function(bundle: Dict, df: pd.DataFrame) -> np.ndarray:
    features = list(bundle["features"])
    missing = sorted(set(features) - set(df.columns))
    if missing:
        raise RuntimeError(f"Lockbox missing frozen model features: {missing}")

    X = df[features].to_numpy(float)
    Xs = (X - np.asarray(bundle["mean"])) / np.asarray(bundle["scale"])
    return bundle["model"].decision_function(Xs)


# =============================================================================
# Tie-aware utility metrics
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


def tie_inclusive_top_fraction(
    y,
    score,
    fraction: float,
) -> Dict[str, float]:
    """
    Select all loci tied at the boundary score. Avoid arbitrary row-order
    behavior when several loci share the same motif-derived score.
    """
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    n = len(y)
    k_nominal = max(1, int(math.ceil(fraction * n)))

    ordered = np.sort(np.unique(score))[::-1]
    selected = np.zeros(n, dtype=bool)
    cum = 0

    for threshold in ordered:
        idx = score == threshold
        selected |= idx
        cum += int(idx.sum())
        if cum >= k_nominal:
            break

    precision = float(y[selected].mean())
    prevalence = float(y.mean())

    return {
        "fraction_nominal": float(fraction),
        "k_nominal": int(k_nominal),
        "k_tie_inclusive": int(selected.sum()),
        "threshold": float(threshold),
        "precision": precision,
        "lift": (
            precision / prevalence if prevalence > 0 else np.nan
        ),
    }


def expected_budget_k(
    y,
    score,
    k: int,
) -> Dict[str, float]:
    """
    Expected positives for an exact experimental budget K if the cutoff score
    falls inside a tie. Loci tied at the boundary are treated exchangeably.
    """
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    k = min(int(k), len(y))

    order_scores = np.sort(np.unique(score))[::-1]
    chosen_n = 0
    expected_pos = 0.0
    boundary_score = np.nan

    for s in order_scores:
        idx = score == s
        n_tie = int(idx.sum())
        p_tie = int(y[idx].sum())

        if chosen_n + n_tie <= k:
            chosen_n += n_tie
            expected_pos += p_tie
            boundary_score = s
            if chosen_n == k:
                break
        else:
            need = k - chosen_n
            expected_pos += need * (p_tie / n_tie)
            chosen_n = k
            boundary_score = s
            break

    precision = expected_pos / k
    prevalence = float(y.mean())
    return {
        "k": int(k),
        "expected_positives": float(expected_pos),
        "expected_precision": float(precision),
        "expected_lift": (
            float(precision / prevalence) if prevalence > 0 else np.nan
        ),
        "boundary_score": float(boundary_score),
    }


def bootstrap_auc(y, score, seed: int) -> Dict[str, float]:
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

    vals = np.asarray(vals, float)
    return {
        "CI_low": float(np.quantile(vals, .025)),
        "CI_high": float(np.quantile(vals, .975)),
        "bootstrap_support_fraction_AUC_gt_0.5": float(np.mean(vals > .5)),
        "bootstrap_replicates": int(len(vals)),
    }


def bootstrap_delta_auc(y, a, b, seed: int) -> Dict[str, float]:
    y = np.asarray(y, int)
    a = np.asarray(a, float)
    b = np.asarray(b, float)

    obs = safe_auc(y, a) - safe_auc(y, b)
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
        "delta_AUROC": float(obs),
        "CI_low": float(np.quantile(vals, .025)),
        "CI_high": float(np.quantile(vals, .975)),
        "bootstrap_support_fraction_delta_gt_0": float(np.mean(vals > 0)),
        "bootstrap_replicates": int(len(vals)),
    }


def evaluate_score(
    df: pd.DataFrame,
    score: np.ndarray,
    name: str,
    seed: int,
) -> Dict:
    y = df["DRS_reported"].astype(int).to_numpy()
    auc = safe_auc(y, score)
    ap = safe_ap(y, score)
    boot = bootstrap_auc(y, score, seed)

    f10 = tie_inclusive_top_fraction(y, score, .10)
    f20 = tie_inclusive_top_fraction(y, score, .20)

    rec = {
        "score": name,
        "n": int(len(y)),
        "positives": int(y.sum()),
        "prevalence": float(y.mean()),
        "AUROC": auc,
        "average_precision": ap,
        **{f"AUROC_{k}": v for k, v in boot.items()},
        "Lift10_tieaware": f10["lift"],
        "Precision10_tieaware": f10["precision"],
        "K10_nominal": f10["k_nominal"],
        "K10_tie_inclusive": f10["k_tie_inclusive"],
        "Lift20_tieaware": f20["lift"],
        "Precision20_tieaware": f20["precision"],
        "K20_nominal": f20["k_nominal"],
        "K20_tie_inclusive": f20["k_tie_inclusive"],
    }

    for k in [50, 100, 200]:
        if len(y) >= k:
            b = expected_budget_k(y, score, k)
            rec[f"budget{k}_expected_positives"] = b["expected_positives"]
            rec[f"budget{k}_expected_precision"] = b["expected_precision"]
            rec[f"budget{k}_expected_lift"] = b["expected_lift"]

    return rec


def score_gradient(df: pd.DataFrame) -> pd.DataFrame:
    z = df.copy()

    # Average-rank percentile keeps equal scores tied.
    z["score_percentile"] = (
        z["TRACEpsi_score"]
        .rank(method="average", pct=True)
    )

    z["quintile"] = np.minimum(
        5,
        np.maximum(
            1,
            np.ceil(z["score_percentile"] * 5).astype(int),
        ),
    )

    rows = []
    for q, g in z.groupby("quintile", sort=True):
        rows.append({
            "quintile": int(q),
            "n": int(len(g)),
            "DRS_reported_n": int(g["DRS_reported"].sum()),
            "DRS_reobservation_rate": float(g["DRS_reported"].mean()),
            "score_min": float(g["TRACEpsi_score"].min()),
            "score_max": float(g["TRACEpsi_score"].max()),
        })
    return pd.DataFrame(rows)


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 108)
    print("TRACE-Ψ PHASE 4C — A549 BID -> DRS EXTERNAL LOCKBOX")
    print("=" * 108)

    for p in [
        PHASE4A_CODE,
        BID_A549,
        DRS,
        FP_HEK,
        MODEL_ARTIFACT,
        PHASE4B_SUMMARY,
        PHASE4B_CONTRACT,
    ]:
        if not p.exists():
            raise RuntimeError(f"Missing required frozen input: {p}")

    with open(PHASE4B_SUMMARY, "r", encoding="utf-8") as f:
        s4b = json.load(f)

    if s4b.get("status") != "TRACEPSI_P2_READY_FOR_A549_LOCKBOX":
        raise RuntimeError(
            "Phase4B did not freeze P2 as ready for A549 lockbox."
        )
    if s4b.get("selected_complexity_for_A549_lockbox") != "P2_TRACEpsi":
        raise RuntimeError("Frozen selected model is not P2_TRACEpsi.")

    phase4a = load_phase4a_module()

    # -------------------------------------------------------------------------
    # 1. Frozen A549 BID source universe
    # -------------------------------------------------------------------------
    bid = phase4a.load_bid_elap_single_u(
        BID_A549,
        assay="BID",
        cell="A549",
    ).copy()

    bid["source_BID"] = 1
    bid["source_BACS"] = 0
    bid["source_ELAP"] = 0
    bid["source_breadth"] = 1
    bid["source_pattern"] = "BID"

    fp = pd.read_csv(FP_HEK, sep="\t", low_memory=False)
    lock = phase4a.add_sequence_features(bid, fp)

    # Development model artifact MUST remain untouched.
    artifact = joblib.load(MODEL_ARTIFACT)

    if artifact.get("selected_complexity_for_lockbox") != "P2_TRACEpsi":
        raise RuntimeError(
            "Model artifact does not freeze P2_TRACEpsi for lockbox."
        )
    if artifact.get("A549_DRS_seen", None) is not False:
        raise RuntimeError(
            "Model artifact does not assert A549_DRS_seen=False."
        )

    models = artifact["fitted_models"]
    if "P2_TRACEpsi" not in models:
        raise RuntimeError("Frozen P2 model missing from artifact.")

    p2_score = decision_function(models["P2_TRACEpsi"], lock)
    lock["TRACEpsi_score"] = p2_score
    lock["TRACEpsi_percentile"] = (
        pd.Series(p2_score, index=lock.index)
        .rank(method="average", pct=True)
    )

    # Secondary frozen ablation if present.
    if "A_context_only" in models:
        lock["context_only_score"] = decision_function(
            models["A_context_only"], lock
        )
    else:
        lock["context_only_score"] = np.nan

    # Predefined simple Aim1->portability bridge.
    lock["minus_source_alignment_score"] = (
        -lock["fp_source_alignment"].astype(float)
    )

    # -------------------------------------------------------------------------
    # 2. OPEN A549 DRS lockbox now
    # -------------------------------------------------------------------------
    drs_pos, drs_report = parse_a549_drs()

    lock["match_key"] = (
        lock["chrom"].map(norm_chr)
        + ":"
        + lock["pos1"].astype(int).astype(str)
        + ":"
        + lock["motif"].astype(str)
    )

    # Coordinate-level audit: if DRS has same coordinate but a different motif,
    # exclude that source locus from the lockbox rather than call it unsupported.
    drs_coord = set(
        zip(
            drs_pos["chrom"].astype(str),
            drs_pos["pos1"].astype(int),
        )
    )
    drs_keys = set(drs_pos["match_key"].astype(str))

    conflict = []
    support = []

    for r in lock.itertuples():
        coord = (str(r.chrom), int(r.pos1))
        key = str(r.match_key)

        if key in drs_keys:
            support.append(1)
            conflict.append(0)
        elif coord in drs_coord:
            support.append(np.nan)
            conflict.append(1)
        else:
            support.append(0)
            conflict.append(0)

    lock["DRS_reported"] = support
    lock["DRS_coordinate_motif_conflict"] = conflict

    conflicts = lock[lock["DRS_coordinate_motif_conflict"].eq(1)].copy()
    conflicts.to_csv(
        META / "phase4c_A549_coordinate_motif_conflicts.tsv",
        sep="\t", index=False
    )

    if len(conflicts) / max(len(lock), 1) > 0.01:
        raise RuntimeError(
            f"A549 coordinate/motif conflicts={len(conflicts)} "
            f"({len(conflicts)/len(lock):.2%}), unexpectedly high. "
            "Review genome build / motif orientation before interpreting lockbox."
        )

    evaluable = lock[lock["DRS_reported"].notna()].copy()
    evaluable["DRS_reported"] = evaluable["DRS_reported"].astype(int)

    # Basic power/label gate.
    n_pos = int(evaluable["DRS_reported"].sum())
    n_neg = int(len(evaluable) - n_pos)

    if n_pos == 0 or n_neg == 0:
        raise RuntimeError(
            f"A549 lockbox has unusable outcome classes: positives={n_pos}, negatives={n_neg}"
        )

    # -------------------------------------------------------------------------
    # 3. Frozen metrics
    # -------------------------------------------------------------------------
    metrics = []

    metrics.append(
        evaluate_score(
            evaluable,
            evaluable["TRACEpsi_score"].to_numpy(float),
            "P2_TRACEpsi_frozen",
            SEED + 1,
        )
    )

    if evaluable["context_only_score"].notna().all():
        metrics.append(
            evaluate_score(
                evaluable,
                evaluable["context_only_score"].to_numpy(float),
                "A_context_only_frozen",
                SEED + 2,
            )
        )

    metrics.append(
        evaluate_score(
            evaluable,
            evaluable["minus_source_alignment_score"].to_numpy(float),
            "minus_source_alignment_predefined",
            SEED + 3,
        )
    )

    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(
        TRACEPSI / "phase4c_A549_lockbox_metrics.tsv",
        sep="\t", index=False
    )

    # -------------------------------------------------------------------------
    # 4. Frozen P2 vs ablation comparisons
    # -------------------------------------------------------------------------
    comparisons = []
    y = evaluable["DRS_reported"].to_numpy(int)

    if evaluable["context_only_score"].notna().all():
        rec = bootstrap_delta_auc(
            y,
            evaluable["TRACEpsi_score"].to_numpy(float),
            evaluable["context_only_score"].to_numpy(float),
            SEED + 100,
        )
        rec.update({
            "model_A": "P2_TRACEpsi_frozen",
            "model_B": "A_context_only_frozen",
        })
        comparisons.append(rec)

    rec = bootstrap_delta_auc(
        y,
        evaluable["TRACEpsi_score"].to_numpy(float),
        evaluable["minus_source_alignment_score"].to_numpy(float),
        SEED + 101,
    )
    rec.update({
        "model_A": "P2_TRACEpsi_frozen",
        "model_B": "minus_source_alignment_predefined",
    })
    comparisons.append(rec)

    comparisons_df = pd.DataFrame(comparisons)
    comparisons_df.to_csv(
        TRACEPSI / "phase4c_A549_lockbox_deltaAUC.tsv",
        sep="\t", index=False
    )

    # -------------------------------------------------------------------------
    # 5. Score-gradient utility
    # -------------------------------------------------------------------------
    gradient = score_gradient(evaluable)
    gradient.to_csv(
        TRACEPSI / "phase4c_A549_TRACEpsi_quintile_gradient.tsv",
        sep="\t", index=False
    )

    evaluable.to_csv(
        TRACEPSI / "phase4c_A549_TRACEpsi_scorecard.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # -------------------------------------------------------------------------
    # 6. Lockbox interpretation
    # -------------------------------------------------------------------------
    p2 = metrics_df[
        metrics_df["score"].eq("P2_TRACEpsi_frozen")
    ].iloc[0]

    strong = bool(
        n_pos >= 15
        and p2["AUROC"] >= 0.60
        and p2["AUROC_CI_low"] > 0.50
        and p2["Lift20_tieaware"] >= 1.20
    )

    supportive = bool(
        n_pos >= 10
        and p2["AUROC"] >= 0.55
        and p2["AUROC_bootstrap_support_fraction_AUC_gt_0.5"] >= 0.90
        and (
            p2["Lift10_tieaware"] >= 1.15
            or p2["Lift20_tieaware"] >= 1.15
        )
    )

    if n_pos < 10:
        status = "TRACEPSI_A549_LOCKBOX_INCONCLUSIVE_LOW_POWER"
    elif strong:
        status = "TRACEPSI_EXTERNALLY_VALIDATED_A549"
    elif supportive:
        status = "TRACEPSI_EXTERNALLY_SUPPORTED_A549"
    else:
        status = "TRACEPSI_A549_DOMAIN_SHIFT_NOT_CONFIRMED"

    fingerprint_external_gain = None
    if len(comparisons_df) and "A_context_only_frozen" in comparisons_df["model_B"].values:
        q = comparisons_df[
            comparisons_df["model_B"].eq("A_context_only_frozen")
        ].iloc[0]
        fingerprint_external_gain = {
            "delta_AUROC_P2_minus_context": float(q["delta_AUROC"]),
            "CI_low": float(q["CI_low"]),
            "CI_high": float(q["CI_high"]),
            "bootstrap_support_fraction_delta_gt_0": float(
                q["bootstrap_support_fraction_delta_gt_0"]
            ),
            "interpretation": (
                "Secondary evidence for transfer of the learned technology-aware "
                "fingerprint beyond generic local sequence context."
            ),
        }

    contract = {
        "phase": "4C_A549_lockbox",
        "tool": "TRACE-Ψ",
        "frozen_model": "Phase4B-v2 P2_TRACEpsi",
        "source": "A549 BID locally-single-U exact calls",
        "target": "A549 DRS/Mod-p ID Supplementary Table S6 positive calls",
        "matching": (
            "exact chromosome + position + identical 5-mer; coordinate/motif "
            "conflicts excluded; no +/-k matching"
        ),
        "no_retraining": True,
        "no_recalibration": True,
        "no_feature_change": True,
        "primary_metric": "AUROC",
        "utility_metrics": [
            "average_precision",
            "tie-aware Lift@10%",
            "tie-aware Lift@20%",
            "expected precision/lift at K=50/100/200",
            "score quintile re-observation gradient",
        ],
        "strong_gate": (
            ">=15 positives; AUROC>=0.60; bootstrap CI lower>0.50; Lift20>=1.20"
        ),
        "supportive_gate": (
            ">=10 positives; AUROC>=0.55; bootstrap support(AUC>0.5)>=0.90; "
            "Lift10>=1.15 or Lift20>=1.15"
        ),
        "outcome_boundary": (
            "DRS non-report is not biological absence and TRACE-Ψ score is not P(true Ψ)."
        ),
        "after_lockbox": (
            "Do not retune on A549. If strong/supportive, freeze method and package "
            "the user-facing TRACE-Ψ scorer. If not confirmed, qualify domain transfer "
            "rather than re-optimizing on this lockbox."
        ),
    }
    (META / "phase4c_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    summary = {
        "status": status,
        "A549_BID_source_rows_total": int(
            pd.read_csv(BID_A549, sep="\t", low_memory=False).shape[0]
        ),
        "A549_BID_locally_singleU_exact_n": int(len(lock)),
        "A549_DRS_QC": drs_report,
        "coordinate_motif_conflicts_excluded": int(len(conflicts)),
        "evaluable_source_loci": int(len(evaluable)),
        "DRS_reobserved_n": n_pos,
        "DRS_reobservation_rate": float(evaluable["DRS_reported"].mean()),
        "lockbox_metrics": metrics_df.to_dict(orient="records"),
        "lockbox_delta_AUC": comparisons_df.to_dict(orient="records"),
        "TRACEpsi_quintile_gradient": gradient.to_dict(orient="records"),
        "fingerprint_external_gain": fingerprint_external_gain,
        "final_tool_interpretation": (
            "TRACE-Ψ ranks evidence portability, not biological truth. "
            "A549 provides an external cell+native-platform transfer test."
        ),
        "next_gate": (
            "If externally supported/validated, Phase4D should package the frozen "
            "scorer + scorecard + inference-audit interface. No further model tuning."
        ),
    }
    (META / "phase4c_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print("A549 BID locally-single-U source loci:", len(lock))
    print("Evaluable after motif-conflict audit:", len(evaluable))
    print("DRS re-observed:", n_pos, "/", len(evaluable))
    print()
    print("LOCKBOX METRICS")
    print(metrics_df.to_string(index=False))
    print()
    print("P2 vs FROZEN ABLATIONS")
    print(comparisons_df.to_string(index=False))
    print()
    print("TRACE-Ψ QUINTILE GRADIENT")
    print(gradient.to_string(index=False))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase4c_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase4c_analysis_contract.json")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4c_A549_lockbox_metrics.tsv")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4c_A549_lockbox_deltaAUC.tsv")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4c_A549_TRACEpsi_quintile_gradient.tsv")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4c_A549_TRACEpsi_scorecard.tsv.gz")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase4c_FATAL_ERROR.txt").write_text(
            err, encoding="utf-8"
        )
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase4c_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
