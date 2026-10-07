#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRACE Phase 3B-R — final robustness closure
===========================================

This is NOT a new scientific Aim.

It closes the only remaining methodological questions in the frozen Phase 3B:
1) the original Ridge alpha grid ended at 1e3 and every outer fold selected
   that upper boundary;
2) one deterministic 5-fold partition may be unstable when delta-R2 is small;
3) one common Ridge penalty may over-shrink the 1–2 low-dimensional Ψ
   variables together with ~325 sequence features.

The biological definitions are unchanged:
    M0 = sequence only
    M1 = M0 + naive Ψ burden
    M2 = M1 + extra cross-assay confirmations

No new outcome.
No new Ψ map.
No new model family.
No DRS inside M1/M2.
No threshold/model shopping after seeing results.

Primary robustness replication
------------------------------
Five independent deterministic gene-fold repeats.
Within every repeat:
- same fold assignment across all TE datasets and M0/M1/M2;
- each TE dataset is modeled independently;
- nested inner CV selects Ridge alpha using outer-training genes only;
- alpha grid spans 1e-4 ... 1e16 PLUS an explicit intercept-only (alpha=inf)
  candidate, so "stronger than maximum tested shrinkage" is identifiable.

Orthogonalized predictive sensitivity
--------------------------------------
For each outer test fold:
1) tune/final-fit sequence-only M0 on outer training genes;
2) create CROSS-FITTED M0 predictions for every outer-training gene using
   only the other training folds (including fold-excluded alpha tuning);
3) fit a LOW-DIMENSIONAL OLS model to those cross-fitted sequence residuals:
       R1: residual ~ naive Ψ
       R2: residual ~ naive Ψ + confirmation
4) add the residual prediction to the untouched outer-test M0 prediction.

This is a predictive partialling-out sensitivity, NOT a causal DML estimator.
Its sole purpose is to test whether the single shared Ridge penalty masked
low-dimensional Ψ information.

Final inference
---------------
Repeated CV partitions are NOT treated as independent biological samples.
For uncertainty:
- predictions are averaged across repeats for each gene;
- gene bootstrap gives dataset-level CIs;
- hierarchical TE-study + gene bootstrap gives macro CIs.

The original Phase3B remains frozen. This script can strengthen, qualify, or
reveal penalty-sensitivity of that result; it does not erase the original run.

Project root
------------
D:\\RNA\\Trace
"""

from __future__ import annotations

import importlib.util
import json
import math
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple, Union

import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr
except Exception as e:
    raise RuntimeError("Phase3B-R requires scipy.") from e


# =============================================================================
# Paths / frozen settings
# =============================================================================

ROOT = Path(r"D:\RNA\Trace")
CODE = ROOT / "code"
META = ROOT / "00_meta"
AIM3 = ROOT / "06_aim3"
LOGS = ROOT / "logs"

for p in [META, AIM3, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

ORIGINAL_SCRIPT = CODE / "14_phase3b_final_sequence_baseline.py"
FUNCTIONAL = AIM3 / "aim3a_functional_cohort.tsv"
ORIGINAL_FEATURE_CONTRACT = META / "phase3b_feature_contract.json"
ORIGINAL_SUMMARY = META / "phase3b_summary.json"
ORIGINAL_TUNING = AIM3 / "aim3b_nested_cv_tuning.tsv"

CACHE_FEATURES = AIM3 / "aim3bR_sequence_feature_matrix.tsv.gz"

N_REPEATS = 5
N_OUTER = 5
N_BOOT = 3000
BASE_SEED = 20261004

# A broad deterministic grid plus an explicit zero-coefficient/intercept-only model.
FINITE_ALPHAS = np.power(10.0, np.arange(-4, 17, dtype=float))  # 1e-4 ... 1e16
ALL_ALPHA_LABELS: List[Union[float, str]] = list(FINITE_ALPHAS) + ["INF"]

M0 = "M0_sequence"
M1 = "M1_sequence_plus_naivePsi"
M2 = "M2_plus_confirmation"

R0 = "R0_sequence"
R1 = "R1_sequence_plus_naivePsi_residual"
R2 = "R2_plus_confirmation_residual"


# =============================================================================
# Basic utilities
# =============================================================================

def r2(y: np.ndarray, pred: np.ndarray) -> float:
    y = np.asarray(y, float)
    pred = np.asarray(pred, float)
    den = float(np.sum((y - y.mean()) ** 2))
    if den <= 0:
        return np.nan
    return float(1.0 - np.sum((y - pred) ** 2) / den)


def rmse(y: np.ndarray, pred: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(y) - np.asarray(pred)) ** 2)))


def rho(y: np.ndarray, pred: np.ndarray) -> float:
    if len(y) < 3:
        return np.nan
    return float(spearmanr(y, pred).statistic)


def alpha_to_label(alpha: Union[float, str]) -> str:
    if alpha == "INF" or (isinstance(alpha, float) and np.isinf(alpha)):
        return "INF"
    return f"{float(alpha):.12g}"


def label_to_alpha(label: str) -> Union[float, str]:
    return "INF" if str(label).upper() == "INF" else float(label)


def balanced_gene_folds(genes: Sequence[str], repeat_id: int) -> Dict[str, int]:
    """
    Deterministic, balanced, repeat-specific global gene folds.
    Same map is used for every TE dataset and every model in the repeat.
    """
    g = np.array(sorted(set(map(str, genes))), dtype=object)
    rng = np.random.default_rng(BASE_SEED + 10007 * repeat_id)
    rng.shuffle(g)
    return {str(x): int(i % N_OUTER) for i, x in enumerate(g)}


# =============================================================================
# Import frozen Phase3B sequence-feature constructor
# =============================================================================

def load_phase3b_module():
    if not ORIGINAL_SCRIPT.exists():
        raise RuntimeError(
            f"Missing frozen Phase3B script: {ORIGINAL_SCRIPT}\n"
            "Keep 14_phase3b_final_sequence_baseline.py in D:\\RNA\\Trace\\code."
        )
    spec = importlib.util.spec_from_file_location(
        "trace_phase3b_frozen", ORIGINAL_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not import frozen Phase3B script.")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_or_load_feature_matrix(
    functional: pd.DataFrame,
    sequence_features: Sequence[str],
) -> Tuple[pd.DataFrame, Dict]:
    """
    Re-use the EXACT frozen Phase3B sequence feature constructor.
    Cache the resulting gene-level matrix for fast reruns.
    """
    required_gene_cols = [
        "gene_name", "transcript_id", "transcript_length",
        "naive_singleU_union_count", "breadth_sum_singleU",
    ]
    missing = [c for c in required_gene_cols if c not in functional.columns]
    if missing:
        raise RuntimeError(f"Functional cohort missing {missing}")

    gene_base = (
        functional[required_gene_cols]
        .drop_duplicates("gene_name")
        .copy()
    )

    gene_base["confirmation_excess"] = (
        pd.to_numeric(gene_base["breadth_sum_singleU"], errors="coerce")
        - pd.to_numeric(gene_base["naive_singleU_union_count"], errors="coerce")
    )
    if (gene_base["confirmation_excess"] < 0).any():
        raise RuntimeError("confirmation_excess < 0 for one or more genes.")

    gene_base["psi_naive_log"] = np.log1p(
        pd.to_numeric(gene_base["naive_singleU_union_count"], errors="coerce")
    )
    gene_base["psi_confirmation_log"] = np.log1p(
        pd.to_numeric(gene_base["confirmation_excess"], errors="coerce")
    )

    # A valid cache must contain the exact current transcript IDs and frozen features.
    if CACHE_FEATURES.exists():
        try:
            cached = pd.read_csv(CACHE_FEATURES, sep="\t", compression="gzip")
            needed = {
                "gene_name", "transcript_id",
                "psi_naive_log", "psi_confirmation_log",
                *sequence_features,
            }
            if needed.issubset(cached.columns):
                left = set(
                    zip(
                        gene_base["gene_name"].astype(str),
                        gene_base["transcript_id"].astype(str),
                    )
                )
                right = set(
                    zip(
                        cached["gene_name"].astype(str),
                        cached["transcript_id"].astype(str),
                    )
                )
                if left == right:
                    return cached, {
                        "mode": "validated_cache",
                        "cache": str(CACHE_FEATURES),
                        "genes": int(len(cached)),
                    }
        except Exception:
            pass

    frozen = load_phase3b_module()
    seq, seq_qc = frozen.extract_sequence_features(gene_base)

    gene = gene_base.merge(
        seq,
        on=["gene_name", "transcript_id"],
        how="inner",
        validate="one_to_one",
    )

    missing_features = [c for c in sequence_features if c not in gene.columns]
    if missing_features:
        raise RuntimeError(
            f"Frozen sequence feature constructor did not reproduce "
            f"{len(missing_features)} contract features: {missing_features[:20]}"
        )

    gene.to_csv(
        CACHE_FEATURES,
        sep="\t",
        index=False,
        compression="gzip",
    )
    return gene, {
        "mode": "recomputed_from_frozen_phase3b",
        "cache": str(CACHE_FEATURES),
        "genes": int(len(gene)),
        "sequence_qc": seq_qc,
    }


# =============================================================================
# Fast standardized Ridge using one SVD per train/validation split
# =============================================================================

def standardize_train_test(
    Xtr: np.ndarray,
    Xte: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0, ddof=0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-12), sd, 1.0)
    return (Xtr - mu) / sd, (Xte - mu) / sd


def ridge_grid_predictions(
    Xtr: np.ndarray,
    ytr: np.ndarray,
    Xva: np.ndarray,
) -> Dict[str, np.ndarray]:
    """
    Standardize from training data only, perform one compact SVD, and produce
    validation predictions for every finite alpha plus the intercept-only model.
    """
    Xs, Xv = standardize_train_test(
        np.asarray(Xtr, float),
        np.asarray(Xva, float),
    )
    y = np.asarray(ytr, float)
    ymean = float(y.mean())
    yc = y - ymean

    U, s, Vt = np.linalg.svd(Xs, full_matrices=False)
    uy = U.T @ yc

    out: Dict[str, np.ndarray] = {}
    for alpha in FINITE_ALPHAS:
        factors = (s * uy) / (s * s + float(alpha))
        beta = Vt.T @ factors
        out[alpha_to_label(alpha)] = ymean + Xv @ beta

    out["INF"] = np.full(len(Xv), ymean, dtype=float)
    return out


def ridge_single_prediction(
    Xtr: np.ndarray,
    ytr: np.ndarray,
    Xte: np.ndarray,
    alpha: Union[float, str],
) -> np.ndarray:
    Xs, Xt = standardize_train_test(
        np.asarray(Xtr, float),
        np.asarray(Xte, float),
    )
    y = np.asarray(ytr, float)
    ymean = float(y.mean())

    if alpha == "INF" or (isinstance(alpha, float) and np.isinf(alpha)):
        return np.full(len(Xt), ymean, dtype=float)

    yc = y - ymean
    U, s, Vt = np.linalg.svd(Xs, full_matrices=False)
    uy = U.T @ yc
    factors = (s * uy) / (s * s + float(alpha))
    beta = Vt.T @ factors
    return ymean + Xt @ beta


# =============================================================================
# Nested tuning with explicit alpha=inf
# =============================================================================

def tune_alpha(
    d: pd.DataFrame,
    features: Sequence[str],
    allowed_folds: Sequence[int],
    fold_col: str,
    cache: Dict,
    cache_key: Tuple,
) -> Tuple[Union[float, str], pd.DataFrame, bool]:
    """
    CV only inside allowed_folds.

    If the largest finite alpha is numerically equivalent to intercept-only,
    choose INF. If the largest finite alpha is genuinely better than INF,
    boundary_unresolved=True and the run will stop rather than silently accept
    another truncated search.
    """
    if cache_key in cache:
        return cache[cache_key]

    allowed = sorted(set(map(int, allowed_folds)))
    if len(allowed) < 3:
        raise RuntimeError(
            f"Need >=3 folds for tuning, received {allowed}."
        )

    sqerr: Dict[str, List[float]] = defaultdict(list)
    ns: Dict[str, List[int]] = defaultdict(list)

    for val_fold in allowed:
        tr = d[
            d[fold_col].isin(allowed)
            & d[fold_col].ne(val_fold)
        ]
        va = d[d[fold_col].eq(val_fold)]

        if len(tr) < 100 or len(va) < 20:
            raise RuntimeError(
                f"Insufficient inner split: train={len(tr)}, val={len(va)}, "
                f"allowed={allowed}, val={val_fold}"
            )

        pred_grid = ridge_grid_predictions(
            tr[list(features)].to_numpy(float),
            tr["TE_z"].to_numpy(float),
            va[list(features)].to_numpy(float),
        )
        yv = va["TE_z"].to_numpy(float)

        for label, pred in pred_grid.items():
            sqerr[label].append(float(np.sum((yv - pred) ** 2)))
            ns[label].append(len(yv))

    rows = []
    for label in [alpha_to_label(a) for a in FINITE_ALPHAS] + ["INF"]:
        total_n = int(np.sum(ns[label]))
        mse = float(np.sum(sqerr[label]) / total_n)
        rows.append({
            "alpha_label": label,
            "alpha_numeric": np.inf if label == "INF" else float(label),
            "inner_MSE": mse,
            "inner_n": total_n,
            "folds_used": "|".join(map(str, allowed)),
        })

    tab = pd.DataFrame(rows)
    min_mse = float(tab["inner_MSE"].min())
    tolerance = max(1e-12, abs(min_mse) * 1e-10)
    tied = tab[tab["inner_MSE"] <= min_mse + tolerance].copy()

    # Numeric ties are resolved toward stronger regularization.
    def strength(row):
        return np.inf if row["alpha_label"] == "INF" else float(row["alpha_numeric"])

    chosen_row = max(
        (r for _, r in tied.iterrows()),
        key=strength,
    )
    chosen_label = str(chosen_row["alpha_label"])
    chosen: Union[float, str] = label_to_alpha(chosen_label)

    max_label = alpha_to_label(FINITE_ALPHAS[-1])
    max_mse = float(
        tab.loc[tab["alpha_label"].eq(max_label), "inner_MSE"].iloc[0]
    )
    inf_mse = float(
        tab.loc[tab["alpha_label"].eq("INF"), "inner_MSE"].iloc[0]
    )

    # At alpha=1e16, standardized Ridge is already effectively intercept-only.
    # If max finite wins by more than numerical tolerance, do not extrapolate.
    relative_gap = (inf_mse - max_mse) / max(abs(inf_mse), 1e-12)
    boundary_unresolved = bool(
        chosen_label == max_label and relative_gap > 1e-8
    )

    tab["selected"] = tab["alpha_label"].eq(chosen_label)
    tab["boundary_unresolved"] = boundary_unresolved
    tab["maxfinite_vs_inf_relative_MSE_gain"] = relative_gap

    result = (chosen, tab, boundary_unresolved)
    cache[cache_key] = result
    return result


# =============================================================================
# Joint M0/M1/M2 repeated nested CV
# =============================================================================

def nested_joint_predictions(
    d: pd.DataFrame,
    features: Sequence[str],
    model_name: str,
    repeat_id: int,
    fold_col: str,
    tune_cache: Dict,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[int, Union[float, str]]]:
    preds = []
    tuning = []
    outer_alphas: Dict[int, Union[float, str]] = {}

    all_folds = sorted(d[fold_col].unique())
    if len(all_folds) != N_OUTER:
        raise RuntimeError(
            f"{d['TE_dataset'].iloc[0]} repeat {repeat_id}: "
            f"expected {N_OUTER} folds, got {all_folds}"
        )

    for outer in all_folds:
        train_folds = [f for f in all_folds if f != outer]

        key = (
            "joint",
            d["TE_dataset"].iloc[0],
            repeat_id,
            model_name,
            tuple(train_folds),
        )
        alpha, tab, unresolved = tune_alpha(
            d=d,
            features=features,
            allowed_folds=train_folds,
            fold_col=fold_col,
            cache=tune_cache,
            cache_key=key,
        )
        if unresolved:
            raise RuntimeError(
                f"Unresolved alpha upper boundary for {key}. "
                "Even alpha=1e16 was materially better than intercept-only."
            )

        q = tab.copy()
        q["TE_dataset"] = d["TE_dataset"].iloc[0]
        q["repeat"] = repeat_id
        q["outer_fold"] = outer
        q["model"] = model_name
        tuning.append(q)

        tr = d[d[fold_col].ne(outer)]
        te = d[d[fold_col].eq(outer)]
        pred = ridge_single_prediction(
            tr[list(features)].to_numpy(float),
            tr["TE_z"].to_numpy(float),
            te[list(features)].to_numpy(float),
            alpha,
        )

        p = te[
            ["gene_name", "TE_dataset", "TE_z", fold_col]
        ].copy()
        p["repeat"] = repeat_id
        p["model"] = model_name
        p["prediction"] = pred
        p["selected_alpha"] = alpha_to_label(alpha)
        preds.append(p)

        outer_alphas[int(outer)] = alpha

    return (
        pd.concat(preds, ignore_index=True),
        pd.concat(tuning, ignore_index=True),
        outer_alphas,
    )


# =============================================================================
# Cross-fitted residual-addition sensitivity
# =============================================================================

def ols_lowdim_predict(
    Xtr: np.ndarray,
    ytr: np.ndarray,
    Xte: np.ndarray,
) -> np.ndarray:
    """
    Low-dimensional OLS with training-only standardization and unpenalized intercept.
    """
    Xs, Xt = standardize_train_test(
        np.asarray(Xtr, float),
        np.asarray(Xte, float),
    )
    A = np.column_stack([np.ones(len(Xs)), Xs])
    At = np.column_stack([np.ones(len(Xt)), Xt])
    beta = np.linalg.pinv(A) @ np.asarray(ytr, float)
    return At @ beta


def partialling_predictions(
    d: pd.DataFrame,
    seq_features: Sequence[str],
    repeat_id: int,
    fold_col: str,
    joint_m0: pd.DataFrame,
    tune_cache: Dict,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Outer-test predictions remain untouched throughout.

    Outer-training sequence residuals are generated by an additional internal
    cross-fitting layer. Each residual observation is predicted by an M0 model
    fit without its own fold, and the nuisance alpha is tuned without that fold.
    """
    out = []
    tuning = []
    all_folds = sorted(d[fold_col].unique())

    m0lookup = joint_m0.set_index(["gene_name", "repeat"])[
        "prediction"
    ].to_dict()

    for outer in all_folds:
        outer_train_folds = [f for f in all_folds if f != outer]
        tr_outer = d[d[fold_col].isin(outer_train_folds)].copy()
        te_outer = d[d[fold_col].eq(outer)].copy()

        # Cross-fitted nuisance predictions for every OUTER-TRAIN gene.
        cf_parts = []

        for held in outer_train_folds:
            nuisance_fit_folds = [
                f for f in outer_train_folds if f != held
            ]

            key = (
                "partial_nuisance",
                d["TE_dataset"].iloc[0],
                repeat_id,
                M0,
                tuple(nuisance_fit_folds),
            )
            alpha, tab, unresolved = tune_alpha(
                d=d,
                features=seq_features,
                allowed_folds=nuisance_fit_folds,
                fold_col=fold_col,
                cache=tune_cache,
                cache_key=key,
            )
            if unresolved:
                raise RuntimeError(
                    f"Unresolved alpha boundary in residual nuisance {key}."
                )

            q = tab.copy()
            q["TE_dataset"] = d["TE_dataset"].iloc[0]
            q["repeat"] = repeat_id
            q["outer_fold"] = outer
            q["heldout_train_fold"] = held
            q["model"] = "M0_nuisance_for_partialling"
            tuning.append(q)

            nuisance_tr = d[d[fold_col].isin(nuisance_fit_folds)]
            nuisance_va = d[d[fold_col].eq(held)]

            pred = ridge_single_prediction(
                nuisance_tr[list(seq_features)].to_numpy(float),
                nuisance_tr["TE_z"].to_numpy(float),
                nuisance_va[list(seq_features)].to_numpy(float),
                alpha,
            )
            z = nuisance_va[
                [
                    "gene_name", "TE_z",
                    "psi_naive_log", "psi_confirmation_log",
                ]
            ].copy()
            z["m0_cf_prediction"] = pred
            z["sequence_residual"] = (
                z["TE_z"] - z["m0_cf_prediction"]
            )
            cf_parts.append(z)

        cf = pd.concat(cf_parts, ignore_index=True)
        if cf["gene_name"].nunique() != tr_outer["gene_name"].nunique():
            raise RuntimeError(
                "Cross-fitted residual universe does not match outer training genes."
            )

        # Outer-test M0 prediction comes from the already nested/tuned joint M0.
        m0_test = np.array(
            [
                m0lookup[(str(g), repeat_id)]
                for g in te_outer["gene_name"].astype(str)
            ],
            dtype=float,
        )

        residual1 = ols_lowdim_predict(
            cf[["psi_naive_log"]].to_numpy(float),
            cf["sequence_residual"].to_numpy(float),
            te_outer[["psi_naive_log"]].to_numpy(float),
        )
        residual2 = ols_lowdim_predict(
            cf[["psi_naive_log", "psi_confirmation_log"]].to_numpy(float),
            cf["sequence_residual"].to_numpy(float),
            te_outer[["psi_naive_log", "psi_confirmation_log"]].to_numpy(float),
        )

        base = te_outer[
            ["gene_name", "TE_dataset", "TE_z", fold_col]
        ].copy()
        base["repeat"] = repeat_id

        for model, pred in [
            (R0, m0_test),
            (R1, m0_test + residual1),
            (R2, m0_test + residual2),
        ]:
            q = base.copy()
            q["model"] = model
            q["prediction"] = pred
            out.append(q)

    return (
        pd.concat(out, ignore_index=True),
        pd.concat(tuning, ignore_index=True),
    )


# =============================================================================
# Metrics / repeated-CV aggregation
# =============================================================================

def performance_table(pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (ds, rep, model), z in pred.groupby(
        ["TE_dataset", "repeat", "model"],
        sort=True,
    ):
        rows.append({
            "TE_dataset": ds,
            "repeat": int(rep),
            "model": model,
            "n_genes": int(len(z)),
            "R2": r2(z["TE_z"], z["prediction"]),
            "RMSE": rmse(z["TE_z"], z["prediction"]),
            "Spearman": rho(z["TE_z"], z["prediction"]),
        })
    return pd.DataFrame(rows)


def repeat_delta_table(
    perf: pd.DataFrame,
    family: str,
) -> pd.DataFrame:
    if family == "joint":
        a, b, c = M0, M1, M2
    elif family == "partial":
        a, b, c = R0, R1, R2
    else:
        raise ValueError(family)

    p = perf.pivot_table(
        index=["TE_dataset", "repeat"],
        columns="model",
        values="R2",
    ).reset_index()

    for req in [a, b, c]:
        if req not in p.columns:
            raise RuntimeError(f"Missing {req} in {family} performance.")

    p = p.rename(
        columns={a: "R2_M0", b: "R2_M1", c: "R2_M2"}
    )
    p["delta_R2_M1_minus_M0"] = p["R2_M1"] - p["R2_M0"]
    p["delta_R2_M2_minus_M1"] = p["R2_M2"] - p["R2_M1"]
    p["family"] = family
    return p


def average_repeat_predictions(pred: pd.DataFrame) -> pd.DataFrame:
    """
    Repeat partitions are algorithmic resamples, not independent evidence.
    Average each gene's OOF prediction across repeats before inference.
    """
    ycheck = (
        pred.groupby(["TE_dataset", "gene_name"])["TE_z"]
        .nunique()
        .max()
    )
    if ycheck != 1:
        raise RuntimeError("TE_z changed across repeats for the same gene.")

    return (
        pred.groupby(
            ["TE_dataset", "gene_name", "model"],
            as_index=False,
        )
        .agg(
            TE_z=("TE_z", "first"),
            prediction=("prediction", "mean"),
            prediction_repeat_SD=("prediction", "std"),
        )
    )


def paired_table(
    averaged: pd.DataFrame,
    ds: str,
    models: Tuple[str, str, str],
) -> pd.DataFrame:
    a, b, c = models
    z = averaged[averaged["TE_dataset"].eq(ds)].copy()

    y = z[["gene_name", "TE_z"]].drop_duplicates("gene_name")
    p = z.pivot(
        index="gene_name",
        columns="model",
        values="prediction",
    ).reset_index()

    out = y.merge(p, on="gene_name", how="inner", validate="one_to_one")
    missing = [m for m in models if m not in out.columns]
    if missing:
        raise RuntimeError(f"{ds}: missing averaged predictions {missing}")

    return out


def deltas_from_paired(
    z: pd.DataFrame,
    models: Tuple[str, str, str],
) -> Dict[str, float]:
    a, b, c = models
    y = z["TE_z"].to_numpy(float)
    r0 = r2(y, z[a].to_numpy(float))
    r1 = r2(y, z[b].to_numpy(float))
    r2v = r2(y, z[c].to_numpy(float))
    return {
        "R2_M0": r0,
        "R2_M1": r1,
        "R2_M2": r2v,
        "delta_R2_M1_minus_M0": r1 - r0,
        "delta_R2_M2_minus_M1": r2v - r1,
    }


def bootstrap_dataset(
    z: pd.DataFrame,
    models: Tuple[str, str, str],
    seed: int,
) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    n = len(z)
    vals10 = []
    vals21 = []

    for _ in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        b = z.iloc[idx]
        if np.var(b["TE_z"].to_numpy(float)) <= 1e-12:
            continue
        d = deltas_from_paired(b, models)
        vals10.append(d["delta_R2_M1_minus_M0"])
        vals21.append(d["delta_R2_M2_minus_M1"])

    return {
        "delta10_CI_low": float(np.quantile(vals10, .025)),
        "delta10_CI_high": float(np.quantile(vals10, .975)),
        "delta21_CI_low": float(np.quantile(vals21, .025)),
        "delta21_CI_high": float(np.quantile(vals21, .975)),
        "bootstrap_replicates": int(len(vals10)),
    }


def averaged_delta_table(
    averaged: pd.DataFrame,
    family: str,
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    models = (M0, M1, M2) if family == "joint" else (R0, R1, R2)

    rows = []
    tables = {}
    for i, ds in enumerate(sorted(averaged["TE_dataset"].unique())):
        z = paired_table(averaged, ds, models)
        tables[ds] = z
        d = deltas_from_paired(z, models)
        b = bootstrap_dataset(
            z, models, BASE_SEED + 50000 + i + (0 if family == "joint" else 1000)
        )
        rows.append({
            "TE_dataset": ds,
            "family": family,
            "n_genes": len(z),
            **d,
            **b,
        })

    return pd.DataFrame(rows), tables


def hierarchical_macro_bootstrap(
    tables: Dict[str, pd.DataFrame],
    models: Tuple[str, str, str],
    seed: int,
) -> Dict[str, float]:
    datasets = sorted(tables)
    rng = np.random.default_rng(seed)

    vals10 = []
    vals21 = []

    for _ in range(N_BOOT):
        chosen = rng.choice(datasets, size=len(datasets), replace=True)
        d10 = []
        d21 = []

        for ds in chosen:
            z = tables[str(ds)]
            idx = rng.integers(0, len(z), size=len(z))
            b = z.iloc[idx]
            if np.var(b["TE_z"].to_numpy(float)) <= 1e-12:
                continue
            d = deltas_from_paired(b, models)
            d10.append(d["delta_R2_M1_minus_M0"])
            d21.append(d["delta_R2_M2_minus_M1"])

        if d10:
            vals10.append(float(np.mean(d10)))
            vals21.append(float(np.mean(d21)))

    return {
        "macro_delta10_CI_low": float(np.quantile(vals10, .025)),
        "macro_delta10_CI_high": float(np.quantile(vals10, .975)),
        "macro_delta21_CI_low": float(np.quantile(vals21, .025)),
        "macro_delta21_CI_high": float(np.quantile(vals21, .975)),
        "bootstrap_replicates": int(len(vals10)),
    }


def macro_summary(
    avg_deltas: pd.DataFrame,
    repeat_deltas: pd.DataFrame,
    tables: Dict[str, pd.DataFrame],
    family: str,
) -> Dict:
    models = (M0, M1, M2) if family == "joint" else (R0, R1, R2)

    rec = {
        "family": family,
        "macro_R2_M0": float(avg_deltas["R2_M0"].mean()),
        "macro_R2_M1": float(avg_deltas["R2_M1"].mean()),
        "macro_R2_M2": float(avg_deltas["R2_M2"].mean()),
        "macro_delta_R2_M1_minus_M0": float(
            avg_deltas["delta_R2_M1_minus_M0"].mean()
        ),
        "macro_delta_R2_M2_minus_M1": float(
            avg_deltas["delta_R2_M2_minus_M1"].mean()
        ),
        "datasets_delta10_positive": int(
            (avg_deltas["delta_R2_M1_minus_M0"] > 0).sum()
        ),
        "datasets_delta21_positive": int(
            (avg_deltas["delta_R2_M2_minus_M1"] > 0).sum()
        ),
        "repeat_dataset_cells_delta10_positive": int(
            (repeat_deltas["delta_R2_M1_minus_M0"] > 0).sum()
        ),
        "repeat_dataset_cells_delta21_positive": int(
            (repeat_deltas["delta_R2_M2_minus_M1"] > 0).sum()
        ),
        "repeat_dataset_cells_total": int(len(repeat_deltas)),
    }
    rec.update(
        hierarchical_macro_bootstrap(
            tables=tables,
            models=models,
            seed=BASE_SEED + (70000 if family == "joint" else 80000),
        )
    )

    rec["increment10_supported_by_frozen_rule"] = bool(
        rec["macro_delta_R2_M1_minus_M0"] > 0
        and rec["macro_delta10_CI_low"] > 0
        and rec["datasets_delta10_positive"] >= 3
    )
    rec["increment21_supported_by_frozen_rule"] = bool(
        rec["macro_delta_R2_M2_minus_M1"] > 0
        and rec["macro_delta21_CI_low"] > 0
        and rec["datasets_delta21_positive"] >= 3
    )
    return rec


# =============================================================================
# Boundary diagnostics
# =============================================================================

def tuning_diagnostics(
    tuning: pd.DataFrame,
    original_tuning: pd.DataFrame | None,
) -> Dict:
    selected = tuning[tuning["selected"].astype(bool)].copy()

    counts = Counter(selected["alpha_label"].astype(str))
    maxfinite = alpha_to_label(FINITE_ALPHAS[-1])
    unresolved = int(selected["boundary_unresolved"].astype(bool).sum())

    old_diag = None
    if original_tuning is not None and len(original_tuning):
        idx = original_tuning.groupby(
            ["TE_dataset", "model", "outer_fold"]
        )["inner_mean_MSE"].idxmin()
        oldsel = original_tuning.loc[idx]
        oldmax = float(original_tuning["alpha"].max())
        old_diag = {
            "old_tuning_contexts": int(len(oldsel)),
            "old_max_alpha": oldmax,
            "old_selected_max_count": int((oldsel["alpha"] == oldmax).sum()),
            "old_selected_max_fraction": float(
                (oldsel["alpha"] == oldmax).mean()
            ),
        }

    return {
        "new_selected_contexts": int(len(selected)),
        "new_selected_alpha_counts": dict(counts),
        "new_intercept_only_selected_count": int(
            (selected["alpha_label"].astype(str) == "INF").sum()
        ),
        "new_max_finite_alpha": float(FINITE_ALPHAS[-1]),
        "new_max_finite_selected_count": int(
            (selected["alpha_label"].astype(str) == maxfinite).sum()
        ),
        "new_boundary_unresolved_count": unresolved,
        "old_boundary_diagnostic": old_diag,
    }


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 108)
    print("TRACE PHASE 3B-R — FINAL ROBUSTNESS CLOSURE")
    print("=" * 108)

    for p in [
        FUNCTIONAL,
        ORIGINAL_FEATURE_CONTRACT,
        ORIGINAL_SUMMARY,
        ORIGINAL_SCRIPT,
    ]:
        if not p.exists():
            raise RuntimeError(f"Missing required frozen input: {p}")

    functional = pd.read_csv(FUNCTIONAL, sep="\t", low_memory=False)

    with open(ORIGINAL_FEATURE_CONTRACT, "r", encoding="utf-8") as f:
        feature_contract_json = json.load(f)

    seq_features = list(feature_contract_json["sequence_features_used"])
    if len(seq_features) < 250:
        raise RuntimeError(
            f"Frozen sequence feature contract unexpectedly small: {len(seq_features)}"
        )

    gene, feature_qc = build_or_load_feature_matrix(
        functional=functional,
        sequence_features=seq_features,
    )

    analysis = functional.merge(
        gene[
            [
                "gene_name", "transcript_id",
                "psi_naive_log", "psi_confirmation_log",
                *seq_features,
            ]
        ],
        on=["gene_name", "transcript_id"],
        how="inner",
        validate="many_to_one",
    )
    analysis = analysis.dropna(
        subset=["TE_z", "psi_naive_log", "psi_confirmation_log", *seq_features]
    ).copy()

    if analysis["gene_name"].nunique() < 700:
        raise RuntimeError(
            f"Only {analysis['gene_name'].nunique()} sequence-eligible genes."
        )
    if analysis["TE_dataset"].nunique() != 5:
        raise RuntimeError(
            f"Expected 5 TE datasets, found {analysis['TE_dataset'].nunique()}."
        )

    model_features = {
        M0: seq_features,
        M1: seq_features + ["psi_naive_log"],
        M2: seq_features + ["psi_naive_log", "psi_confirmation_log"],
    }

    tune_cache: Dict = {}
    joint_preds_all = []
    partial_preds_all = []
    joint_tuning_all = []
    partial_tuning_all = []
    fold_rows = []

    global_genes = sorted(analysis["gene_name"].astype(str).unique())

    for repeat_id in range(N_REPEATS):
        fold_map = balanced_gene_folds(global_genes, repeat_id)
        fold_col = f"fold_repeat_{repeat_id}"
        analysis[fold_col] = (
            analysis["gene_name"].astype(str).map(fold_map).astype(int)
        )

        for g, f0 in fold_map.items():
            fold_rows.append({
                "repeat": repeat_id,
                "gene_name": g,
                "outer_fold": f0,
            })

        print(f"\nRepeat {repeat_id + 1}/{N_REPEATS}")

        for ds, d0 in analysis.groupby("TE_dataset", sort=True):
            d = d0.copy()
            fold_counts = d[fold_col].value_counts().sort_index().to_dict()
            if len(fold_counts) != N_OUTER or min(fold_counts.values()) < 25:
                raise RuntimeError(
                    f"{ds}, repeat {repeat_id}: bad fold counts {fold_counts}"
                )

            print(f"  {ds}: n={len(d)}")

            ds_joint = []
            m0_joint = None

            for model_name in [M0, M1, M2]:
                pred, tuning, _ = nested_joint_predictions(
                    d=d,
                    features=model_features[model_name],
                    model_name=model_name,
                    repeat_id=repeat_id,
                    fold_col=fold_col,
                    tune_cache=tune_cache,
                )
                ds_joint.append(pred)
                joint_tuning_all.append(tuning)
                if model_name == M0:
                    m0_joint = pred.copy()

            joint_preds_all.extend(ds_joint)

            partial_pred, partial_tune = partialling_predictions(
                d=d,
                seq_features=seq_features,
                repeat_id=repeat_id,
                fold_col=fold_col,
                joint_m0=m0_joint,
                tune_cache=tune_cache,
            )
            partial_preds_all.append(partial_pred)
            partial_tuning_all.append(partial_tune)

    folds = pd.DataFrame(fold_rows)
    folds.to_csv(
        AIM3 / "aim3bR_repeated_global_gene_folds.tsv",
        sep="\t", index=False
    )

    joint_pred = pd.concat(joint_preds_all, ignore_index=True)
    partial_pred = pd.concat(partial_preds_all, ignore_index=True)
    joint_tuning = pd.concat(joint_tuning_all, ignore_index=True)
    partial_tuning = pd.concat(partial_tuning_all, ignore_index=True)

    joint_pred.to_csv(
        AIM3 / "aim3bR_joint_oof_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )
    partial_pred.to_csv(
        AIM3 / "aim3bR_partial_oof_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )
    joint_tuning.to_csv(
        AIM3 / "aim3bR_joint_tuning.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )
    partial_tuning.to_csv(
        AIM3 / "aim3bR_partial_tuning.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # Repeat-level algorithmic stability.
    joint_perf = performance_table(joint_pred)
    partial_perf = performance_table(partial_pred)
    joint_repeat_delta = repeat_delta_table(joint_perf, "joint")
    partial_repeat_delta = repeat_delta_table(partial_perf, "partial")

    joint_perf.to_csv(
        AIM3 / "aim3bR_joint_repeat_performance.tsv",
        sep="\t", index=False
    )
    partial_perf.to_csv(
        AIM3 / "aim3bR_partial_repeat_performance.tsv",
        sep="\t", index=False
    )
    joint_repeat_delta.to_csv(
        AIM3 / "aim3bR_joint_repeat_deltas.tsv",
        sep="\t", index=False
    )
    partial_repeat_delta.to_csv(
        AIM3 / "aim3bR_partial_repeat_deltas.tsv",
        sep="\t", index=False
    )

    # Repeat-averaged OOF predictions are the inferential surface.
    joint_avg = average_repeat_predictions(joint_pred)
    partial_avg = average_repeat_predictions(partial_pred)

    joint_avg.to_csv(
        AIM3 / "aim3bR_joint_repeat_averaged_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )
    partial_avg.to_csv(
        AIM3 / "aim3bR_partial_repeat_averaged_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    joint_avg_delta, joint_tables = averaged_delta_table(
        joint_avg, "joint"
    )
    partial_avg_delta, partial_tables = averaged_delta_table(
        partial_avg, "partial"
    )

    joint_avg_delta.to_csv(
        AIM3 / "aim3bR_joint_averaged_deltas.tsv",
        sep="\t", index=False
    )
    partial_avg_delta.to_csv(
        AIM3 / "aim3bR_partial_averaged_deltas.tsv",
        sep="\t", index=False
    )

    joint_macro = macro_summary(
        joint_avg_delta,
        joint_repeat_delta,
        joint_tables,
        "joint",
    )
    partial_macro = macro_summary(
        partial_avg_delta,
        partial_repeat_delta,
        partial_tables,
        "partial",
    )

    original_tuning = (
        pd.read_csv(ORIGINAL_TUNING, sep="\t")
        if ORIGINAL_TUNING.exists()
        else None
    )
    alpha_diag = tuning_diagnostics(
        tuning=pd.concat(
            [joint_tuning, partial_tuning],
            ignore_index=True,
        ),
        original_tuning=original_tuning,
    )

    if alpha_diag["new_boundary_unresolved_count"] > 0:
        raise RuntimeError(
            "At least one broadened-alpha tuning context remains unresolved. "
            "Do not interpret Phase3B-R."
        )

    # Final closure logic: original frozen rule is applied independently to
    # the joint Ridge replication and the low-dimensional partialling sensitivity.
    j10 = bool(joint_macro["increment10_supported_by_frozen_rule"])
    j21 = bool(joint_macro["increment21_supported_by_frozen_rule"])
    p10 = bool(partial_macro["increment10_supported_by_frozen_rule"])
    p21 = bool(partial_macro["increment21_supported_by_frozen_rule"])

    if j10 and p10:
        psi_status = "ROBUST_PSI_INCREMENT"
    elif (not j10) and (not p10):
        psi_status = "NO_ROBUST_PSI_INCREMENT_CONFIRMED"
    else:
        psi_status = "PSI_INCREMENT_PENALTY_METHOD_SENSITIVE"

    if j21 and p21:
        provenance_status = "ROBUST_PROVENANCE_INCREMENT"
    elif (not j21) and (not p21):
        provenance_status = "NO_ROBUST_PROVENANCE_INCREMENT_CONFIRMED"
    else:
        provenance_status = "PROVENANCE_INCREMENT_PENALTY_METHOD_SENSITIVE"

    if (
        psi_status == "NO_ROBUST_PSI_INCREMENT_CONFIRMED"
        and provenance_status == "NO_ROBUST_PROVENANCE_INCREMENT_CONFIRMED"
    ):
        status = "PHASE3B_ROBUST_NULL_CLOSED"
    elif (
        psi_status == "ROBUST_PSI_INCREMENT"
        or provenance_status == "ROBUST_PROVENANCE_INCREMENT"
    ):
        status = "PHASE3B_ROBUST_INCREMENT_DETECTED"
    else:
        status = "PHASE3B_RESULT_METHOD_SENSITIVE_REQUIRES_QUALIFIED_INTERPRETATION"

    contract = {
        "phase": "3B_R_final_robustness_closure",
        "scientific_story": "Calibration -> Measurement -> Evidence -> Inference",
        "frozen_primary_result": (
            "Original Phase3B remains the registered primary run; this closure "
            "addresses alpha-boundary, partition stability, and shared-penalty sensitivity."
        ),
        "unchanged_biology": {
            "M0": "same frozen MANE sequence features",
            "M1": "M0 + same naive BID/BACS/ELAP single-U union burden",
            "M2": "M1 + same extra cross-assay confirmation count",
            "DRS_in_M1_M2": False,
            "new_outcomes": False,
            "new_Psi_thresholds": False,
        },
        "alpha_grid": (
            "1e-4 through 1e16 in decade steps plus explicit intercept-only alpha=INF"
        ),
        "repeated_nested_CV": {
            "repeats": N_REPEATS,
            "outer_folds": N_OUTER,
            "same_gene_folds_across_models_and_TE_datasets_within_repeat": True,
            "repeat_partitions_not_treated_as_independent_samples": True,
        },
        "partialling_sensitivity": (
            "cross-fitted sequence-only residuals on outer-training genes; "
            "unpenalized low-dimensional OLS residual addition for naive Ψ and confirmation"
        ),
        "uncertainty": (
            "average each gene's OOF predictions across repeats, then dataset-level "
            "gene bootstrap and hierarchical TE-study+gene bootstrap"
        ),
        "frozen_positive_rule": (
            "macro delta R2 > 0, hierarchical 95% CI lower > 0, "
            "and positive delta in >=3/5 TE datasets"
        ),
        "interpretation_boundary": (
            "This is predictive incremental information conditional on the frozen "
            "sequence baseline, not causal evidence for Ψ effects on translation."
        ),
        "final_primary_analysis_after_this": False,
    }
    (META / "phase3bR_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    summary = {
        "status": status,
        "psi_increment_status": psi_status,
        "provenance_increment_status": provenance_status,
        "feature_matrix_QC": feature_qc,
        "functional_sequence_eligible_genes": int(
            analysis["gene_name"].nunique()
        ),
        "TE_datasets": sorted(
            analysis["TE_dataset"].astype(str).unique().tolist()
        ),
        "sequence_features_n": int(len(seq_features)),
        "alpha_boundary_diagnostics": alpha_diag,
        "joint_repeat_averaged_dataset_deltas": (
            joint_avg_delta.to_dict(orient="records")
        ),
        "partial_repeat_averaged_dataset_deltas": (
            partial_avg_delta.to_dict(orient="records")
        ),
        "joint_macro": joint_macro,
        "partial_macro": partial_macro,
        "final_scientific_gate": (
            "Freeze the main computational study after this closure. "
            "If both approaches remain null, conclude that generic Ψ burden and "
            "cross-assay confirmation do not provide robust sequence-independent "
            "TE predictive information in these five HeLa datasets. "
            "Do not interpret this as absence of context-specific Ψ function."
        ),
    }
    (META / "phase3bR_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print("Ψ increment:", psi_status)
    print("Provenance increment:", provenance_status)
    print()
    print("ALPHA DIAGNOSTICS")
    print(json.dumps(alpha_diag, ensure_ascii=False, indent=2))
    print()
    print("JOINT RIDGE — REPEAT-AVERAGED DELTAS")
    print(joint_avg_delta.to_string(index=False))
    print()
    print("PARTIALLED — REPEAT-AVERAGED DELTAS")
    print(partial_avg_delta.to_string(index=False))
    print()
    print("JOINT MACRO")
    print(json.dumps(joint_macro, indent=2))
    print()
    print("PARTIAL MACRO")
    print(json.dumps(partial_macro, indent=2))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase3bR_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase3bR_analysis_contract.json")
    print(r"  D:\RNA\Trace\06_aim3\aim3bR_joint_averaged_deltas.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3bR_partial_averaged_deltas.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3bR_joint_repeat_deltas.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3bR_partial_repeat_deltas.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3bR_joint_tuning.tsv.gz")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase3bR_FATAL_ERROR.txt").write_text(
            err, encoding="utf-8"
        )
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase3bR_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
