#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRACE Phase 3A v3 — heterogeneity closure before final M0/M1/M2
================================================================

Purpose
-------
Phase 3A v2 showed that changing the published Ψ map changes the apparent
Ψ-burden / translation-efficiency association on the same HeLa cohort.

However, the five HeLa TE datasets are independent public experimental
contexts, not technical replicates of one homogeneous outcome. Therefore a
single common burden slope must be stress-tested before it is interpreted.

This closure analysis does four things and NOTHING else:

1) Formal burden × TE-dataset interaction tests for each Ψ map.
2) Fixed-effect and random-effects meta-analysis of the five per-dataset slopes.
3) Leave-one-TE-dataset-out influence analysis.
4) The same per-dataset/meta-analysis audit for the four technology-aware
   secondary burden definitions from Phase 3A v2.

No thresholds are changed.
No maps are added or removed.
No functional outcome is added.
No sequence model is fit here.

Interpretation
--------------
If map effects differ across TE datasets, the paper must say that functional
inference is jointly conditioned by Ψ-map provenance AND TE-study context.
It must not present the stacked common slope as a universal HeLa effect.

Project root
------------
D:\\RNA\\Trace
"""

from __future__ import annotations

import json
import math
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

try:
    from scipy.optimize import minimize_scalar
    from scipy.stats import chi2, t as t_dist
except Exception as e:
    raise RuntimeError(
        "Phase3A v3 requires scipy (optimize + stats), which is standard in the "
        "current scientific Python environment."
    ) from e

ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
AIM3 = ROOT / "06_aim3"
LOGS = ROOT / "logs"

for p in [META, AIM3, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

FUNCTIONAL = AIM3 / "aim3a_functional_cohort.tsv"
PER_DATASET_V2 = AIM3 / "aim3a_TE_per_dataset_effects.tsv"
PAIRWISE_V2 = AIM3 / "aim3a_map_pairwise_bootstrap.tsv"

MAP_BURDENS = {
    "BID": "BID_count",
    "BACS": "BACS_count",
    "ELAP": "ELAP_count",
    "DRS": "DRS_count",
}

TECH_BURDENS = [
    "naive_singleU_union_count",
    "orthogonal_singleU_count",
    "triple_singleU_count",
    "breadth_sum_singleU",
]


# -----------------------------------------------------------------------------
# utilities
# -----------------------------------------------------------------------------

def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    mu = x.mean()
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index)
    return (x - mu) / sd


def normal_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2)) if np.isfinite(z) else np.nan


def bh_adjust(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    if len(vals) == 0:
        return out
    order = np.argsort(vals)
    ranked = vals[order]
    q = np.empty(len(vals))
    prev = 1.0
    n = len(vals)
    for i in range(n - 1, -1, -1):
        rank = i + 1
        v = min(prev, ranked[i] * n / rank)
        q[order[i]] = v
        prev = v
    out[np.where(ok)[0]] = np.minimum(q, 1.0)
    return out


def hc3_ols(y: np.ndarray, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    inv = np.linalg.pinv(X.T @ X)
    beta = inv @ (X.T @ y)
    resid = y - X @ beta
    h = np.sum((X @ inv) * X, axis=1)
    adj = resid / np.clip(1 - h, 1e-6, None)
    meat = X.T @ ((adj ** 2)[:, None] * X)
    cov = inv @ meat @ inv
    return beta, cov


def cluster_ols(
    y: np.ndarray,
    X: np.ndarray,
    clusters: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    inv = np.linalg.pinv(X.T @ X)
    beta = inv @ (X.T @ y)
    resid = y - X @ beta

    groups = defaultdict(list)
    for i, g in enumerate(clusters.astype(str)):
        groups[g].append(i)

    meat = np.zeros((X.shape[1], X.shape[1]))
    for ids in groups.values():
        sg = X[ids].T @ resid[ids]
        meat += np.outer(sg, sg)

    G = len(groups)
    N, K = X.shape
    correction = (G / max(G - 1, 1)) * ((N - 1) / max(N - K, 1))
    cov = inv @ meat @ inv * correction
    return beta, cov


# -----------------------------------------------------------------------------
# per-dataset effects
# -----------------------------------------------------------------------------

def per_dataset_effects(
    functional: pd.DataFrame,
    burden_col: str,
    feature_label: str,
) -> pd.DataFrame:
    rows = []

    for ds, z0 in functional.groupby("TE_dataset", sort=True):
        z = z0[
            ["gene_name", "TE_z", burden_col, "transcript_length"]
        ].dropna().copy()

        if len(z) < 100:
            continue

        z["burden_z"] = zscore(np.log1p(
            pd.to_numeric(z[burden_col], errors="coerce").clip(lower=0)
        ))
        z["length_z"] = zscore(np.log(
            pd.to_numeric(z["transcript_length"], errors="coerce").clip(lower=1)
        ))
        z = z.dropna()

        X = np.column_stack([
            np.ones(len(z)),
            z["burden_z"].to_numpy(float),
            z["length_z"].to_numpy(float),
        ])
        y = z["TE_z"].to_numpy(float)

        beta, cov = hc3_ols(y, X)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))

        b = float(beta[1])
        s0 = float(se[1])
        zz = b / s0 if s0 > 0 else np.nan

        rows.append({
            "feature": feature_label,
            "burden_column": burden_col,
            "TE_dataset": ds,
            "beta_standardized": b,
            "SE_HC3": s0,
            "CI_low": b - 1.96 * s0,
            "CI_high": b + 1.96 * s0,
            "p": normal_p(zz),
            "n_genes": int(len(z)),
        })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# random-effects meta-analysis
# -----------------------------------------------------------------------------

def reml_tau2(y: np.ndarray, v: np.ndarray) -> float:
    """
    Restricted maximum-likelihood tau^2 for a normal-normal meta-analysis.
    """
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    k = len(y)

    def nll(log_tau2):
        tau2 = math.exp(log_tau2)
        w = 1.0 / (v + tau2)
        mu = np.sum(w * y) / np.sum(w)
        resid = y - mu

        # REML negative log likelihood up to constants.
        return 0.5 * (
            np.sum(np.log(v + tau2))
            + np.log(np.sum(w))
            + np.sum(w * resid * resid)
        )

    # Permit values arbitrarily close to zero.
    result = minimize_scalar(
        nll,
        bounds=(math.log(1e-12), math.log(max(1.0, np.var(y) * 100 + 1e-6))),
        method="bounded",
        options={"xatol": 1e-12},
    )
    tau2 = float(math.exp(result.x))

    # If boundary likelihood is effectively equal, treat as zero.
    nll0 = nll(math.log(1e-12))
    if nll0 <= result.fun + 1e-8:
        tau2 = 0.0
    return tau2


def meta_one(z: pd.DataFrame) -> Dict:
    yi = z["beta_standardized"].to_numpy(float)
    sei = z["SE_HC3"].to_numpy(float)
    vi = sei ** 2
    k = len(yi)

    # Fixed effect.
    wf = 1.0 / vi
    mu_f = float(np.sum(wf * yi) / np.sum(wf))
    se_f = float(math.sqrt(1.0 / np.sum(wf)))
    zf = mu_f / se_f

    Q = float(np.sum(wf * (yi - mu_f) ** 2))
    q_df = k - 1
    q_p = float(chi2.sf(Q, q_df))
    I2 = float(max(0.0, (Q - q_df) / Q) * 100.0) if Q > 0 else 0.0

    # REML random effect + modified Hartung-Knapp CI.
    tau2 = reml_tau2(yi, vi)
    wr = 1.0 / (vi + tau2)
    mu_r = float(np.sum(wr * yi) / np.sum(wr))

    q_star = float(np.sum(wr * (yi - mu_r) ** 2) / max(k - 1, 1))
    hk_scale = max(1.0, q_star)  # modified KH: never narrower than unscaled RE
    se_r = float(math.sqrt(hk_scale / np.sum(wr)))

    crit = float(t_dist.ppf(0.975, df=k - 1))
    ci_r_lo = mu_r - crit * se_r
    ci_r_hi = mu_r + crit * se_r
    tstat = mu_r / se_r if se_r > 0 else np.nan
    p_r = float(2 * t_dist.sf(abs(tstat), df=k - 1)) if np.isfinite(tstat) else np.nan

    # Prediction interval: conservative t-based interval.
    pred_se = math.sqrt(max(0.0, tau2 + se_r ** 2))
    pred_lo = mu_r - crit * pred_se
    pred_hi = mu_r + crit * pred_se

    return {
        "k_TE_datasets": k,
        "fixed_beta": mu_f,
        "fixed_SE": se_f,
        "fixed_CI_low": mu_f - 1.96 * se_f,
        "fixed_CI_high": mu_f + 1.96 * se_f,
        "fixed_p": normal_p(zf),
        "Cochran_Q": Q,
        "Q_df": q_df,
        "Q_p": q_p,
        "I2_percent": I2,
        "tau2_REML": tau2,
        "random_beta_mKH": mu_r,
        "random_SE_mKH": se_r,
        "random_CI_low_mKH": ci_r_lo,
        "random_CI_high_mKH": ci_r_hi,
        "random_p_mKH": p_r,
        "prediction_interval_low": pred_lo,
        "prediction_interval_high": pred_hi,
    }


def meta_table(perds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature, z in perds.groupby("feature", sort=False):
        rec = {"feature": feature}
        rec.update(meta_one(z))
        rows.append(rec)

    out = pd.DataFrame(rows)
    out["q_BH_random_effect"] = bh_adjust(out["random_p_mKH"])
    out["heterogeneous_Q_FDR"] = bh_adjust(out["Q_p"])
    return out


# -----------------------------------------------------------------------------
# burden x TE-dataset interaction
# -----------------------------------------------------------------------------

def interaction_test(
    functional: pd.DataFrame,
    burden_col: str,
    feature_label: str,
) -> Dict:
    z = functional[
        ["gene_name", "TE_dataset", "TE_z", burden_col, "transcript_length"]
    ].dropna().copy()

    # Standardize gene-level predictors ONCE over the shared functional genes.
    gene = z[
        ["gene_name", burden_col, "transcript_length"]
    ].drop_duplicates("gene_name").copy()
    gene["burden_z"] = zscore(np.log1p(
        pd.to_numeric(gene[burden_col], errors="coerce").clip(lower=0)
    ))
    gene["length_z"] = zscore(np.log(
        pd.to_numeric(gene["transcript_length"], errors="coerce").clip(lower=1)
    ))

    z = (
        z.drop(columns=[burden_col, "transcript_length"])
         .merge(
             gene[["gene_name", "burden_z", "length_z"]],
             on="gene_name",
             how="inner",
             validate="many_to_one",
         )
         .dropna()
    )

    datasets = sorted(z["TE_dataset"].unique())
    if len(datasets) < 2:
        raise RuntimeError("Need >=2 TE datasets for interaction analysis.")

    base = datasets[0]

    cols = [
        np.ones(len(z)),
        z["burden_z"].to_numpy(float),
        z["length_z"].to_numpy(float),
    ]
    names = ["intercept", "burden_z", "length_z"]

    for ds in datasets[1:]:
        ind = (z["TE_dataset"] == ds).astype(float).to_numpy()
        cols.append(ind)
        names.append(f"dataset::{ds}")

    interaction_names = []
    for ds in datasets[1:]:
        ind = (z["TE_dataset"] == ds).astype(float).to_numpy()
        cols.append(ind * z["burden_z"].to_numpy(float))
        nm = f"burden_x::{ds}"
        names.append(nm)
        interaction_names.append(nm)

    X = np.column_stack(cols)
    y = z["TE_z"].to_numpy(float)
    clusters = z["gene_name"].astype(str).to_numpy()

    beta, cov = cluster_ols(y, X, clusters)

    ids = [names.index(nm) for nm in interaction_names]
    b = beta[ids]
    V = cov[np.ix_(ids, ids)]
    stat = float(b.T @ np.linalg.pinv(V) @ b)
    df = len(ids)
    p = float(chi2.sf(stat, df))

    # Recover dataset-specific slopes on the common gene-level scaling.
    slopes = {base: float(beta[names.index("burden_z")])}
    for ds, nm in zip(datasets[1:], interaction_names):
        slopes[ds] = float(
            beta[names.index("burden_z")] + beta[names.index(nm)]
        )

    return {
        "feature": feature_label,
        "burden_column": burden_col,
        "base_dataset": base,
        "interaction_Wald_chi2": stat,
        "interaction_df": df,
        "interaction_p": p,
        "n_stacked_rows": int(len(z)),
        "n_genes": int(z["gene_name"].nunique()),
        "dataset_slopes_common_scaling": json.dumps(slopes, ensure_ascii=False),
    }


# -----------------------------------------------------------------------------
# leave-one-study-out
# -----------------------------------------------------------------------------

def leave_one_out(perds: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for feature, z in perds.groupby("feature", sort=False):
        datasets = list(z["TE_dataset"])
        for omitted in datasets:
            zz = z[z["TE_dataset"] != omitted].copy()
            rec = {
                "feature": feature,
                "omitted_TE_dataset": omitted,
                "remaining_k": len(zz),
            }
            m = meta_one(zz)
            for key in [
                "fixed_beta", "fixed_CI_low", "fixed_CI_high",
                "I2_percent", "tau2_REML",
                "random_beta_mKH", "random_CI_low_mKH",
                "random_CI_high_mKH", "random_p_mKH"
            ]:
                rec[key] = m[key]
            rows.append(rec)

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------

def main():
    print("=" * 104)
    print("TRACE PHASE 3A v3 — TE-STUDY HETEROGENEITY CLOSURE")
    print("=" * 104)

    for p in [FUNCTIONAL, PER_DATASET_V2, PAIRWISE_V2]:
        if not p.exists():
            raise RuntimeError(f"Missing required Phase3A input: {p}")

    functional = pd.read_csv(FUNCTIONAL, sep="\t", low_memory=False)
    v2_per = pd.read_csv(PER_DATASET_V2, sep="\t", low_memory=False)
    pairwise = pd.read_csv(PAIRWISE_V2, sep="\t", low_memory=False)

    required = {
        "gene_name", "TE_dataset", "TE_z", "transcript_length",
        *MAP_BURDENS.values(), *TECH_BURDENS
    }
    missing = sorted(required - set(functional.columns))
    if missing:
        raise RuntimeError(f"Functional cohort missing columns: {missing}")

    if functional["TE_dataset"].nunique() != 5:
        raise RuntimeError(
            f"Expected 5 HeLa TE datasets, found {functional['TE_dataset'].nunique()}."
        )

    # Recalculate map per-study effects from the frozen cohort, rather than
    # trusting the previous TSV blindly.
    map_per = pd.concat([
        per_dataset_effects(functional, bcol, map_name)
        for map_name, bcol in MAP_BURDENS.items()
    ], ignore_index=True)

    # Cross-check the v2 values. Small numerical differences are allowed only
    # at floating-point precision.
    check = map_per.merge(
        v2_per[["TE_dataset", "map", "beta_standardized"]],
        left_on=["TE_dataset", "feature"],
        right_on=["TE_dataset", "map"],
        how="left",
        suffixes=("_v3", "_v2"),
    )
    max_diff = float(
        np.nanmax(np.abs(
            check["beta_standardized_v3"] - check["beta_standardized_v2"]
        ))
    )
    if max_diff > 1e-8:
        raise RuntimeError(
            f"Recomputed Phase3A map effects differ from v2 by up to {max_diff}. "
            "Do not reinterpret until cohort consistency is resolved."
        )

    tech_per = pd.concat([
        per_dataset_effects(functional, c, c)
        for c in TECH_BURDENS
    ], ignore_index=True)

    all_per = pd.concat([map_per, tech_per], ignore_index=True)
    all_per.to_csv(
        AIM3 / "aim3a_v3_per_dataset_all_features.tsv",
        sep="\t", index=False
    )

    meta = meta_table(all_per)
    meta.to_csv(
        AIM3 / "aim3a_v3_meta_analysis.tsv",
        sep="\t", index=False
    )

    interactions = pd.DataFrame([
        interaction_test(functional, bcol, label)
        for label, bcol in {
            **MAP_BURDENS,
            **{c: c for c in TECH_BURDENS},
        }.items()
    ])
    interactions["q_BH_interaction"] = bh_adjust(interactions["interaction_p"])
    interactions.to_csv(
        AIM3 / "aim3a_v3_dataset_interactions.tsv",
        sep="\t", index=False
    )

    loo = leave_one_out(all_per)
    loo.to_csv(
        AIM3 / "aim3a_v3_leave_one_dataset_out.tsv",
        sep="\t", index=False
    )

    # Concise map-only interpretation table.
    map_meta = meta[meta["feature"].isin(MAP_BURDENS)].copy()
    map_inter = interactions[
        interactions["feature"].isin(MAP_BURDENS)
    ][["feature", "interaction_p", "q_BH_interaction"]]

    map_closure = map_meta.merge(map_inter, on="feature", how="left")
    map_closure["stable_direction_all_5"] = map_closure["feature"].map(
        map_per.groupby("feature")["beta_standardized"].apply(
            lambda x: bool((x > 0).all() or (x < 0).all())
        )
    )
    map_closure.to_csv(
        AIM3 / "aim3a_v3_map_closure.tsv",
        sep="\t", index=False
    )

    # Existing pairwise map differences remain valid as same-cohort
    # marginal-effect comparisons, but must be interpreted alongside
    # TE-study heterogeneity.
    pairwise_sig = bool(
        (pd.to_numeric(pairwise["q_BH_pairwise"], errors="coerce") < .05).any()
    )

    hetero_maps = map_closure.loc[
        map_closure["q_BH_interaction"] < .05, "feature"
    ].tolist()

    high_I2_maps = map_closure.loc[
        map_closure["I2_percent"] >= 50, "feature"
    ].tolist()

    if pairwise_sig and (hetero_maps or high_I2_maps):
        status = "AIM3A_CLOSED_MAP_PROVENANCE_AND_TE_CONTEXT_BOTH_MATTER"
    elif pairwise_sig:
        status = "AIM3A_CLOSED_MAP_PROVENANCE_WITH_STABLE_TE_CONTEXT"
    else:
        status = "AIM3A_CLOSED_NO_ROBUST_MAP_DIFFERENCE"

    contract = {
        "phase": "3A_v3_closure",
        "reason": (
            "Phase3A v2 per-dataset results showed potential TE-study slope "
            "heterogeneity, so the common stacked slope is formally audited "
            "before the final sequence-baseline analysis."
        ),
        "frozen_data": (
            "Same Phase3A functional cohort, same four maps, same five HeLa TE "
            "datasets, same transcript-length covariate."
        ),
        "heterogeneity_tests": [
            "burden x TE-dataset joint Wald test with gene-cluster robust covariance",
            "Cochran Q and I2 across five per-dataset slopes",
            "REML random-effects meta-analysis with modified Hartung-Knapp CI",
            "leave-one-TE-dataset-out influence analysis",
        ],
        "interpretation_rule": (
            "When substantial heterogeneity is present, random-effects and "
            "per-dataset estimates supersede the common stacked slope as the "
            "generalizable functional-effect summary."
        ),
        "technology_aware_secondary": (
            "The same heterogeneity audit is applied to naive union, orthogonal, "
            "triple, and breadth-sum burdens; no positive secondary result is "
            "carried into Phase3B without this audit."
        ),
        "next_step": (
            "After closure, Phase3B must evaluate M0/M1/M2 separately within each "
            "TE dataset with common gene folds, then aggregate predictive increments "
            "across datasets. It must not fit one universal TE slope."
        ),
    }
    (META / "phase3a_v3_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    summary = {
        "status": status,
        "functional_rows": int(len(functional)),
        "functional_genes": int(functional["gene_name"].nunique()),
        "TE_datasets": sorted(functional["TE_dataset"].unique().tolist()),
        "v2_effect_recalculation_max_abs_difference": max_diff,
        "pairwise_map_difference_after_BH_exists": pairwise_sig,
        "maps_with_FDR_significant_dataset_interaction": hetero_maps,
        "maps_with_I2_at_least_50_percent": high_I2_maps,
        "map_closure": map_closure.to_dict(orient="records"),
        "technology_aware_meta": meta[
            meta["feature"].isin(TECH_BURDENS)
        ].to_dict(orient="records"),
        "technology_aware_interactions": interactions[
            interactions["feature"].isin(TECH_BURDENS)
        ].to_dict(orient="records"),
        "scientific_boundary": (
            "Phase3A supports map-dependent functional inference only to the "
            "extent supported by same-cohort pairwise map differences. A common "
            "pooled burden slope is not treated as universal when TE-study "
            "heterogeneity is present."
        ),
        "next_gate": (
            "If closure confirms heterogeneity, proceed to a per-TE-dataset "
            "nested-CV M0/M1/M2 design. Do not pool all TE datasets into one "
            "common-slope prediction model."
        ),
    }
    (META / "phase3a_v3_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print()
    print("MAP CLOSURE")
    print(map_closure.to_string(index=False))
    print()
    print("TECHNOLOGY-AWARE META")
    print(meta[meta["feature"].isin(TECH_BURDENS)].to_string(index=False))
    print()
    print("DATASET INTERACTIONS")
    print(interactions.to_string(index=False))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase3a_v3_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase3a_v3_analysis_contract.json")
    print(r"  D:\RNA\Trace\06_aim3\aim3a_v3_map_closure.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3a_v3_meta_analysis.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3a_v3_dataset_interactions.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3a_v3_leave_one_dataset_out.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase3a_v3_FATAL_ERROR.txt").write_text(
            err, encoding="utf-8"
        )
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase3a_v3_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
