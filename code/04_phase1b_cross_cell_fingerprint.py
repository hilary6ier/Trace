#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE Phase 1B — exact-locus-disjoint cross-cell 5-mer fingerprint
==================================================================

Primary scientific question
---------------------------
Does a sequence-context signature learned from BID-vs-ELAP reported-map
discordance in one cell line transfer to the other cell line?

Key design choices
------------------
- Single-U exact-locus universe inherited from Aim 1A.
- Remove every exact genomic locus observed in BOTH cell lines from BOTH
  training and test sets before primary transfer analysis.
- Train only on BID-exclusive-reported vs ELAP-exclusive-reported loci.
- Shared loci are never used for fitting; they are used only for a
  predeclared continuum check.
- Primary score is an interpretable Jeffreys-smoothed 5-mer log relative
  frequency ratio, not a black-box classifier.
- No feature/model tuning.
- Synthetic calibration decomposition is secondary/explanatory:
  it asks how much of the transferable motif score aligns with the
  BID-vs-ELAP synthetic calibration contrast, and whether a residual
  fingerprint remains transferable.

Project root:
    D:\RNA\Trace

Inputs:
    04_aim1\aim1a_site_table.tsv
    04_aim1\calibration_resolved.tsv

Outputs:
    00_meta\phase1b_analysis_contract.json
    00_meta\phase1b_summary.json
    04_aim1\aim1b_transfer_results.tsv
    04_aim1\aim1b_alpha_sensitivity.tsv
    04_aim1\aim1b_shared_continuum.tsv
    04_aim1\aim1b_motif_effects.tsv
    04_aim1\aim1b_calibration_decomposition.tsv
"""

from __future__ import annotations

import json
import math
import traceback
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
AIM1 = ROOT / "04_aim1"
LOGS = ROOT / "logs"

for p in [META, AIM1, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

SITE_PATH = AIM1 / "aim1a_site_table.tsv"
CAL_PATH = AIM1 / "calibration_resolved.tsv"

SEED = 20261001
RNG = np.random.default_rng(SEED)

PRIMARY_ALPHA = 0.5
SENSITIVITY_ALPHAS = [0.25, 0.5, 1.0, 2.0]
N_BOOT = 3000
N_PERM = 5000
MOTIF_UNIVERSE_SIZE = 256

BID_CLASS = "BID_exclusive_reported"
ELAP_CLASS = "ELAP_exclusive_reported"
SHARED_CLASS = "shared"


def auc_score(y, score) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    ok = np.isfinite(score)
    y = y[ok]
    score = score[ok]
    n1 = int((y == 1).sum())
    n0 = int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = pd.Series(score).rank(method="average").to_numpy()
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def learn_motif_score(train: pd.DataFrame, alpha: float) -> Tuple[Dict[str, float], float, pd.DataFrame]:
    """
    score(m) = log relative frequency of motif m in ELAP-exclusive reported
               minus log relative frequency of motif m in BID-exclusive reported.

    Jeffreys alpha=0.5 is primary. 256 is fixed by the NNUNN motif universe.
    """
    n_bid = int((train["y"] == 0).sum())
    n_elap = int((train["y"] == 1).sum())

    ct = (
        train.groupby(["motif", "y"])
        .size()
        .unstack(fill_value=0)
    )
    if 0 not in ct.columns:
        ct[0] = 0
    if 1 not in ct.columns:
        ct[1] = 0

    score = {}
    rows = []
    motifs = sorted(set(train["motif"]))
    for motif in motifs:
        nb = int(ct.loc[motif, 0]) if motif in ct.index else 0
        ne = int(ct.loc[motif, 1]) if motif in ct.index else 0
        value = (
            math.log((ne + alpha) / (n_elap + alpha * MOTIF_UNIVERSE_SIZE))
            - math.log((nb + alpha) / (n_bid + alpha * MOTIF_UNIVERSE_SIZE))
        )
        score[motif] = value
        rows.append({
            "motif": motif,
            "n_BID": nb,
            "n_ELAP": ne,
            "n_total": nb + ne,
            "human_motif_log_ratio": value,
        })

    unseen = (
        math.log(alpha / (n_elap + alpha * MOTIF_UNIVERSE_SIZE))
        - math.log(alpha / (n_bid + alpha * MOTIF_UNIVERSE_SIZE))
    )
    return score, unseen, pd.DataFrame(rows)


def motif_cluster_bootstrap_auc(test: pd.DataFrame, n_boot: int = N_BOOT) -> Tuple[float, float]:
    motifs = test["motif"].dropna().unique()
    groups = {m: test[test["motif"] == m] for m in motifs}
    vals = []
    for _ in range(n_boot):
        sampled = RNG.choice(motifs, size=len(motifs), replace=True)
        b = pd.concat([groups[m] for m in sampled], ignore_index=True)
        a = auc_score(b["y"], b["score"])
        if np.isfinite(a):
            vals.append(a)
    if not vals:
        return float("nan"), float("nan")
    return float(np.quantile(vals, .025)), float(np.quantile(vals, .975))


def motif_mapping_permutation_p(
    test: pd.DataFrame,
    learned_score: Dict[str, float],
    unseen_score: float,
    observed_auc: float,
    n_perm: int = N_PERM,
) -> float:
    """
    Permute learned score assignments across motif identities.
    This preserves the distribution of learned motif effects but destroys
    train->test motif alignment: a direct test of cross-cell reproducibility.
    """
    motifs = list(learned_score)
    values = np.array([learned_score[m] for m in motifs], dtype=float)
    row_motifs = test["motif"].to_numpy()
    y = test["y"].to_numpy()

    ge = 0
    for _ in range(n_perm):
        perm = RNG.permutation(values)
        mapping = dict(zip(motifs, perm))
        score = np.array([mapping.get(m, unseen_score) for m in row_motifs], dtype=float)
        a = auc_score(y, score)
        if a >= observed_auc:
            ge += 1
    return float((ge + 1) / (n_perm + 1))


def residual_mapping_permutation_p(
    test: pd.DataFrame,
    motif_table: pd.DataFrame,
    observed_auc: float,
    n_perm: int = N_PERM,
) -> float:
    mp = motif_table.set_index("motif")["residual_score"].to_dict()
    motifs = list(mp)
    values = np.array([mp[m] for m in motifs], float)
    row_motifs = test["motif"].to_numpy()
    y = test["y"].to_numpy()
    fallback = float(np.median(values))

    ge = 0
    for _ in range(n_perm):
        perm = RNG.permutation(values)
        pm = dict(zip(motifs, perm))
        score = np.array([pm.get(m, fallback) for m in row_motifs], float)
        a = auc_score(y, score)
        if a >= observed_auc:
            ge += 1
    return float((ge + 1) / (n_perm + 1))


def weighted_calibration_decomposition(
    motif_table: pd.DataFrame,
    calibration_contrast: Dict[str, float],
) -> Tuple[pd.DataFrame, Dict]:
    """
    Weighted least squares across motifs observed in training:
        human motif log-ratio ~ synthetic calibration contrast

    Weight = number of training loci carrying the motif.
    This makes the decomposition match the locus-weighted estimand while
    avoiding outcome-driven feature selection.

    Weighted R² is descriptive: fraction of training-cell motif-score variance
    aligned with the fixed synthetic calibration contrast.
    """
    all_motifs = sorted(calibration_contrast)
    base = pd.DataFrame({"motif": all_motifs})
    x = base.merge(motif_table, on="motif", how="left")
    x["n_total"] = x["n_total"].fillna(0)
    x["human_motif_log_ratio"] = x["human_motif_log_ratio"].fillna(np.nan)
    x["calibration_contrast"] = x["motif"].map(calibration_contrast)

    observed = x[(x["n_total"] > 0) & x["human_motif_log_ratio"].notna()].copy()
    X = np.column_stack([
        np.ones(len(observed)),
        observed["calibration_contrast"].to_numpy(float),
    ])
    y = observed["human_motif_log_ratio"].to_numpy(float)
    w = observed["n_total"].to_numpy(float)

    xtwx = X.T @ (w[:, None] * X)
    xtwy = X.T @ (w * y)
    beta = np.linalg.pinv(xtwx) @ xtwy

    pred_obs = X @ beta
    ybar = float(np.average(y, weights=w))
    sse = float(np.sum(w * (y - pred_obs) ** 2))
    sst = float(np.sum(w * (y - ybar) ** 2))
    wr2 = float(1 - sse / sst) if sst > 0 else float("nan")

    x["calibration_component"] = beta[0] + beta[1] * x["calibration_contrast"]

    # For motifs not observed in training, human score is the empirical prior
    # score. Caller fills this before invoking this function.
    x["residual_score"] = x["human_motif_log_ratio"] - x["calibration_component"]

    meta = {
        "intercept": float(beta[0]),
        "calibration_slope": float(beta[1]),
        "weighted_R2_training_motif_score": wr2,
        "n_observed_training_motifs": int(len(observed)),
    }
    return x, meta


def prepare(site: pd.DataFrame):
    required = {"cell_line", "locus", "class", "motif"}
    missing = required - set(site.columns)
    if missing:
        raise RuntimeError(f"aim1a_site_table missing columns: {sorted(missing)}")

    site = site.copy()
    site["motif"] = site["motif"].fillna("").astype(str)

    hek_loci = set(site.loc[site["cell_line"] == "HEK293T", "locus"])
    hela_loci = set(site.loc[site["cell_line"] == "HeLa", "locus"])
    cross_cell_overlap = hek_loci & hela_loci

    # Primary no-locus-leak universe.
    disjoint = site[~site["locus"].isin(cross_cell_overlap)].copy()

    return site, disjoint, cross_cell_overlap


def transfer_direction(
    disjoint: pd.DataFrame,
    train_cell: str,
    test_cell: str,
    alpha: float,
    calibration_contrast: Dict[str, float],
):
    train = disjoint[
        (disjoint["cell_line"] == train_cell)
        & disjoint["class"].isin([BID_CLASS, ELAP_CLASS])
        & disjoint["motif"].ne("")
    ].copy()
    test = disjoint[
        (disjoint["cell_line"] == test_cell)
        & disjoint["class"].isin([BID_CLASS, ELAP_CLASS])
        & disjoint["motif"].ne("")
    ].copy()

    train["y"] = (train["class"] == ELAP_CLASS).astype(int)
    test["y"] = (test["class"] == ELAP_CLASS).astype(int)

    learned, unseen, motif_table = learn_motif_score(train, alpha)

    # Add all 256 motifs with empirical-prior score for unseen motifs.
    full = pd.DataFrame({"motif": sorted(calibration_contrast)})
    full["human_motif_log_ratio"] = full["motif"].map(learned).fillna(unseen)

    counts = (
        train.groupby(["motif", "y"]).size().unstack(fill_value=0)
        if len(train) else pd.DataFrame()
    )
    if 0 not in counts.columns:
        counts[0] = 0
    if 1 not in counts.columns:
        counts[1] = 0
    count_rows = []
    for m in full["motif"]:
        nb = int(counts.loc[m, 0]) if m in counts.index else 0
        ne = int(counts.loc[m, 1]) if m in counts.index else 0
        count_rows.append((m, nb, ne, nb + ne))
    cdf = pd.DataFrame(count_rows, columns=["motif", "n_BID", "n_ELAP", "n_total"])
    full = full.merge(cdf, on="motif", how="left")

    test["score"] = test["motif"].map(learned).fillna(unseen)
    observed_auc = auc_score(test["y"], test["score"])
    ci_lo, ci_hi = motif_cluster_bootstrap_auc(test)
    perm_p = motif_mapping_permutation_p(test, learned, unseen, observed_auc)

    # Calibration decomposition, secondary/explanatory.
    decomposed, decmeta = weighted_calibration_decomposition(full, calibration_contrast)

    score_map = decomposed.set_index("motif")["human_motif_log_ratio"].to_dict()
    cal_map = decomposed.set_index("motif")["calibration_component"].to_dict()
    resid_map = decomposed.set_index("motif")["residual_score"].to_dict()

    test["calibration_component"] = test["motif"].map(cal_map)
    test["residual_score"] = test["motif"].map(resid_map)

    cal_auc = auc_score(test["y"], test["calibration_component"])
    resid_auc = auc_score(test["y"], test["residual_score"])

    test_resid = test.copy()
    test_resid["score"] = test_resid["residual_score"]
    resid_ci_lo, resid_ci_hi = motif_cluster_bootstrap_auc(test_resid)
    resid_p = residual_mapping_permutation_p(test, decomposed, resid_auc)

    result = {
        "train_cell": train_cell,
        "test_cell": test_cell,
        "alpha": alpha,
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "n_train_BID_exclusive": int((train["y"] == 0).sum()),
        "n_train_ELAP_exclusive": int((train["y"] == 1).sum()),
        "n_test_BID_exclusive": int((test["y"] == 0).sum()),
        "n_test_ELAP_exclusive": int((test["y"] == 1).sum()),
        "n_train_motifs_observed": int(train["motif"].nunique()),
        "n_test_motifs_observed": int(test["motif"].nunique()),
        "n_test_loci_with_unseen_train_motif": int((~test["motif"].isin(learned)).sum()),
        "AUC_cross_cell_5mer": float(observed_auc),
        "AUC_CI_low_motif_boot": float(ci_lo),
        "AUC_CI_high_motif_boot": float(ci_hi),
        "motif_alignment_permutation_p": float(perm_p),
        "calibration_component_AUC_test": float(cal_auc),
        "residual_fingerprint_AUC_test": float(resid_auc),
        "residual_AUC_CI_low_motif_boot": float(resid_ci_lo),
        "residual_AUC_CI_high_motif_boot": float(resid_ci_hi),
        "residual_motif_alignment_permutation_p": float(resid_p),
        **decmeta,
    }

    decomposed["train_cell"] = train_cell
    decomposed["test_cell"] = test_cell

    return result, decomposed, learned, unseen


def shared_continuum(
    disjoint: pd.DataFrame,
    train_cell: str,
    test_cell: str,
    learned: Dict[str, float],
    unseen: float,
):
    d = disjoint[
        (disjoint["cell_line"] == test_cell)
        & disjoint["class"].isin([BID_CLASS, SHARED_CLASS, ELAP_CLASS])
        & disjoint["motif"].ne("")
    ].copy()
    d["score"] = d["motif"].map(learned).fillna(unseen)

    order = [BID_CLASS, SHARED_CLASS, ELAP_CLASS]
    rows = []
    medians = {}
    for cls in order:
        x = d.loc[d["class"] == cls, "score"].astype(float)
        medians[cls] = float(x.median()) if len(x) else float("nan")
        rows.append({
            "train_cell": train_cell,
            "test_cell": test_cell,
            "class": cls,
            "n": int(len(x)),
            "median_score": medians[cls],
            "mean_score": float(x.mean()) if len(x) else float("nan"),
            "q25": float(x.quantile(.25)) if len(x) else float("nan"),
            "q75": float(x.quantile(.75)) if len(x) else float("nan"),
        })

    ordered = (
        np.isfinite(medians[BID_CLASS])
        and np.isfinite(medians[SHARED_CLASS])
        and np.isfinite(medians[ELAP_CLASS])
        and medians[BID_CLASS] <= medians[SHARED_CLASS] <= medians[ELAP_CLASS]
    )
    for r in rows:
        r["predeclared_median_order_holds"] = bool(ordered)
    return rows


def motif_cross_cell_correlation(
    effects_a: pd.DataFrame,
    effects_b: pd.DataFrame,
) -> Dict:
    a = effects_a[["motif", "human_motif_log_ratio", "n_total"]].rename(
        columns={"human_motif_log_ratio": "score_HEK", "n_total": "n_HEK"}
    )
    b = effects_b[["motif", "human_motif_log_ratio", "n_total"]].rename(
        columns={"human_motif_log_ratio": "score_HeLa", "n_total": "n_HeLa"}
    )
    m = a.merge(b, on="motif", how="inner")
    z = m[(m["n_HEK"] > 0) & (m["n_HeLa"] > 0)].copy()
    if len(z) < 3:
        return {"n_common_observed_motifs": len(z), "pearson": np.nan, "spearman": np.nan}
    return {
        "n_common_observed_motifs": int(len(z)),
        "pearson": float(z["score_HEK"].corr(z["score_HeLa"], method="pearson")),
        "spearman": float(z["score_HEK"].corr(z["score_HeLa"], method="spearman")),
    }


def main():
    print("=" * 96)
    print("TRACE PHASE 1B — EXACT-LOCUS-DISJOINT CROSS-CELL 5-MER FINGERPRINT")
    print("=" * 96)

    if not SITE_PATH.exists():
        raise RuntimeError(f"Missing: {SITE_PATH}")
    if not CAL_PATH.exists():
        raise RuntimeError(f"Missing: {CAL_PATH}")

    site = pd.read_csv(SITE_PATH, sep="\t")
    cal = pd.read_csv(CAL_PATH, sep="\t")

    # Fixed BID-vs-ELAP synthetic calibration contrast from Aim1A.
    c = cal[cal["assay"].isin(["BID", "ELAP"])].pivot(
        index="motif", columns="assay", values="response_percentile"
    )
    if set(c.columns) != {"BID", "ELAP"} or len(c) != 256:
        raise RuntimeError("Expected exactly 256 BID+ELAP calibration motifs.")
    c["contrast"] = c["ELAP"] - c["BID"]
    calibration_contrast = c["contrast"].to_dict()

    site_all, disjoint, overlap = prepare(site)

    contract = {
        "phase": "1B",
        "primary_question": (
            "Does a BID-vs-ELAP 5-mer sequence fingerprint learned in one cell line "
            "transfer to the other after exact genomic-locus leakage removal?"
        ),
        "primary_universe": "Aim1A locally single-U exact loci",
        "cross_cell_exact_locus_overlap_removed_from_both_cells": True,
        "n_exact_loci_removed": int(len(overlap)),
        "training_classes": [BID_CLASS, ELAP_CLASS],
        "shared_sites_used_for_training": False,
        "primary_score": (
            "Jeffreys-smoothed 5-mer log relative-frequency ratio; alpha=0.5; "
            "256-motif universe"
        ),
        "primary_metric": "held-out-cell AUROC",
        "confidence_interval": "motif-cluster bootstrap",
        "null_test": (
            "permute learned motif-score assignments across motif identities before "
            "application to held-out cell"
        ),
        "alpha_sensitivity": SENSITIVITY_ALPHAS,
        "secondary_calibration_decomposition": (
            "training-cell motif log-ratio ~ fixed synthetic calibration contrast, "
            "weighted by training motif counts; test calibration component and residual "
            "fingerprint on held-out cell"
        ),
        "interpretation_boundary": (
            "Cross-cell transfer supports a stable technology-associated sequence "
            "fingerprint, not proof that the fingerprint is entirely assay-caused. "
            "Stable biological sequence preferences may contribute and are addressed "
            "in subsequent biological/annotation controls."
        ),
    }
    (META / "phase1b_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    transfer_rows = []
    decomposition_tables = {}
    primary_models = {}

    for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
        result, effects, learned, unseen = transfer_direction(
            disjoint,
            train_cell,
            test_cell,
            PRIMARY_ALPHA,
            calibration_contrast,
        )
        transfer_rows.append(result)
        decomposition_tables[(train_cell, test_cell)] = effects
        primary_models[(train_cell, test_cell)] = (learned, unseen)

    transfer_df = pd.DataFrame(transfer_rows)
    transfer_df.to_csv(AIM1 / "aim1b_transfer_results.tsv", sep="\t", index=False)

    # Alpha sensitivity: report all values; never select the best alpha.
    sens_rows = []
    for alpha in SENSITIVITY_ALPHAS:
        for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
            result, _, _, _ = transfer_direction(
                disjoint,
                train_cell,
                test_cell,
                alpha,
                calibration_contrast,
            )
            sens_rows.append({
                "alpha": alpha,
                "train_cell": train_cell,
                "test_cell": test_cell,
                "AUC": result["AUC_cross_cell_5mer"],
                "CI_low": result["AUC_CI_low_motif_boot"],
                "CI_high": result["AUC_CI_high_motif_boot"],
            })
    pd.DataFrame(sens_rows).to_csv(
        AIM1 / "aim1b_alpha_sensitivity.tsv", sep="\t", index=False
    )

    # Shared-site continuum with fixed cross-cell score.
    continuum_rows = []
    for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
        learned, unseen = primary_models[(train_cell, test_cell)]
        continuum_rows.extend(
            shared_continuum(disjoint, train_cell, test_cell, learned, unseen)
        )
    continuum_df = pd.DataFrame(continuum_rows)
    continuum_df.to_csv(AIM1 / "aim1b_shared_continuum.tsv", sep="\t", index=False)

    # Motif effects / calibration decomposition.
    effects_all = []
    for (train_cell, test_cell), df in decomposition_tables.items():
        z = df.copy()
        z["train_cell"] = train_cell
        z["test_cell"] = test_cell
        effects_all.append(z)
    effects_df = pd.concat(effects_all, ignore_index=True)
    effects_df.to_csv(AIM1 / "aim1b_motif_effects.tsv", sep="\t", index=False)

    dec_cols = [
        "train_cell", "test_cell",
        "weighted_R2_training_motif_score",
        "calibration_slope",
        "calibration_component_AUC_test",
        "residual_fingerprint_AUC_test",
        "residual_AUC_CI_low_motif_boot",
        "residual_AUC_CI_high_motif_boot",
        "residual_motif_alignment_permutation_p",
    ]
    transfer_df[dec_cols].to_csv(
        AIM1 / "aim1b_calibration_decomposition.tsv", sep="\t", index=False
    )

    # Cross-cell motif-effect correlation using the two independently learned,
    # exact-locus-disjoint scores.
    e_hek = decomposition_tables[("HEK293T", "HeLa")]
    e_hela = decomposition_tables[("HeLa", "HEK293T")]
    correlation = motif_cross_cell_correlation(e_hek, e_hela)

    both_ci_above_half = bool(
        (transfer_df["AUC_CI_low_motif_boot"] > 0.5).all()
    )
    both_perm = bool(
        (transfer_df["motif_alignment_permutation_p"] <= 0.01).all()
    )
    alpha_stable = bool(
        pd.DataFrame(sens_rows)
        .groupby(["train_cell", "test_cell"])["AUC"]
        .min()
        .gt(0.5)
        .all()
    )
    continuum_holds_both = bool(
        continuum_df
        .groupby(["train_cell", "test_cell"])["predeclared_median_order_holds"]
        .first()
        .all()
    )

    status = (
        "AIM1B_STRONG_CROSS_CELL_FINGERPRINT"
        if both_ci_above_half and both_perm and alpha_stable
        else "AIM1B_RESULT_READY_FOR_INTERPRETATION"
    )

    summary = {
        "status": status,
        "n_exact_cross_cell_loci_removed_from_both_cells": int(len(overlap)),
        "primary_transfer": transfer_df.to_dict(orient="records"),
        "motif_effect_cross_cell_correlation": correlation,
        "shared_continuum_order_holds_both_directions": continuum_holds_both,
        "alpha_sensitivity_all_AUC_above_0.5": alpha_stable,
        "interpretation": (
            "If strong: the human BID-vs-ELAP reported maps contain a reproducible "
            "5-mer technology-associated fingerprint that transfers across biological "
            "contexts even after exact-locus leakage removal. The calibration decomposition "
            "then distinguishes the component aligned with measured synthetic assay chemistry "
            "from the remaining transferable sequence structure."
        ),
        "next_gate": (
            "Do not add larger models automatically. If the 5-mer fingerprint is strong, "
            "next test known biological/measurement explanations (writer motifs, region, "
            "confidence) and only then decide whether a common GRCh38 +/-10 nt sequence "
            "sensitivity analysis adds scientific value."
        ),
    }
    (META / "phase1b_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print("Exact cross-cell loci removed from BOTH cells:", len(overlap))
    print()
    print(transfer_df.to_string(index=False))
    print()
    print("Motif-effect cross-cell correlation:", correlation)
    print("Shared median order holds both directions:", continuum_holds_both)
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase1b_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase1b_analysis_contract.json")
    print(r"  D:\RNA\Trace\04_aim1\aim1b_transfer_results.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1b_alpha_sensitivity.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1b_shared_continuum.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1b_calibration_decomposition.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1b_motif_effects.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase1b_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase1b_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
