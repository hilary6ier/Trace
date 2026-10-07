#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE Phase 1D — Aim 1 closure
==============================

Purpose
-------
Close Aim 1 rigorously after the strong Phase 1B/1C results.

This script does NOT add a larger classifier. It:
1) fixes the invalid ELAP confidence sensitivity from Phase 1C;
2) uses ELAP's AUTHOR-DEFINED confidence categories rather than inventing
   cross-assay raw-signal thresholds;
3) formally re-tests the fixed 5-mer cross-cell fingerprint after excluding
   lower-confidence ELAP calls and, separately, using highest-confidence
   ELAP calls only;
4) recomputes within-CDS transfer with motif-cluster CI/permutation;
5) adds multiplicity-controlled positional-base effect interpretation;
6) emits a single Aim-1 closure decision.

Scientific boundaries
---------------------
- Published call absence is not a biological negative.
- No raw FASTQ/BAM.
- No model search.
- No cross-assay comparison of raw signal units.
- BID published mRNA sites already satisfy BID's published calling criteria
  and >10% estimated modification-fraction cutoff.
- ELAP confidence sensitivity uses the published confidence classes, whose
  definition combines stop ratio and replicate detection.

Project root:
    D:\RNA\Trace
"""

from __future__ import annotations

import json
import math
import re
import traceback
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
NORM = ROOT / "03_harmonized" / "normalized_sources"
AIM1 = ROOT / "04_aim1"
LOGS = ROOT / "logs"

for p in [META, AIM1, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

SITE = AIM1 / "aim1a_site_table.tsv"
ANN = AIM1 / "aim1c_locus_annotations.tsv"
POS = AIM1 / "aim1c_positional_base_effects.tsv"
PHASE1C_SUMMARY = META / "phase1c_summary.json"

ELAP_FILES = {
    "HEK293T": NORM / "ELAP_HEK293T.tsv",
    "HeLa": NORM / "ELAP_HeLa.tsv",
}

BID_CLASS = "BID_exclusive_reported"
ELAP_CLASS = "ELAP_exclusive_reported"
SHARED_CLASS = "shared"

ALPHA = 0.5
MOTIF_UNIVERSE = 256
SEED = 20261001
RNG = np.random.default_rng(SEED)
N_BOOT = 3000
N_PERM = 5000


# ---------------------------------------------------------------------
# Core cross-cell score
# ---------------------------------------------------------------------

def auc_score(y, score) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    ok = np.isfinite(score)
    y, score = y[ok], score[ok]
    n1 = int((y == 1).sum())
    n0 = int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = pd.Series(score).rank(method="average").to_numpy()
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def learn_motif_score(train: pd.DataFrame):
    n0 = int((train["y"] == 0).sum())
    n1 = int((train["y"] == 1).sum())
    ct = train.groupby(["motif", "y"]).size().unstack(fill_value=0)
    if 0 not in ct.columns:
        ct[0] = 0
    if 1 not in ct.columns:
        ct[1] = 0

    score = {}
    for m in ct.index:
        nb = int(ct.loc[m, 0])
        ne = int(ct.loc[m, 1])
        score[m] = (
            math.log((ne + ALPHA) / (n1 + ALPHA * MOTIF_UNIVERSE))
            - math.log((nb + ALPHA) / (n0 + ALPHA * MOTIF_UNIVERSE))
        )

    unseen = (
        math.log(ALPHA / (n1 + ALPHA * MOTIF_UNIVERSE))
        - math.log(ALPHA / (n0 + ALPHA * MOTIF_UNIVERSE))
    )
    return score, unseen


def motif_boot(test: pd.DataFrame):
    motifs = test["motif"].dropna().unique()
    if len(motifs) < 5:
        return np.nan, np.nan
    groups = {m: test[test["motif"] == m] for m in motifs}
    vals = []
    for _ in range(N_BOOT):
        picked = RNG.choice(motifs, size=len(motifs), replace=True)
        z = pd.concat([groups[m] for m in picked], ignore_index=True)
        a = auc_score(z["y"], z["score"])
        if np.isfinite(a):
            vals.append(a)
    if not vals:
        return np.nan, np.nan
    return float(np.quantile(vals, .025)), float(np.quantile(vals, .975))


def motif_perm(test: pd.DataFrame, learned: Dict[str, float], unseen: float, obs: float):
    motifs = list(learned)
    vals = np.asarray([learned[m] for m in motifs], float)
    row_motifs = test["motif"].to_numpy()
    y = test["y"].to_numpy()
    ge = 0
    for _ in range(N_PERM):
        pm = dict(zip(motifs, RNG.permutation(vals)))
        score = np.asarray([pm.get(m, unseen) for m in row_motifs], float)
        if auc_score(y, score) >= obs:
            ge += 1
    return float((ge + 1) / (N_PERM + 1))


def transfer(train: pd.DataFrame, test: pd.DataFrame) -> Dict:
    train = train.copy()
    test = test.copy()
    train["y"] = (train["class"] == ELAP_CLASS).astype(int)
    test["y"] = (test["class"] == ELAP_CLASS).astype(int)

    if train["y"].nunique() < 2 or test["y"].nunique() < 2:
        return {
            "n_train": len(train), "n_test": len(test),
            "AUC": np.nan, "CI_low": np.nan, "CI_high": np.nan,
            "permutation_p": np.nan,
        }

    learned, unseen = learn_motif_score(train)
    test["score"] = test["motif"].map(learned).fillna(unseen)
    obs = auc_score(test["y"], test["score"])
    lo, hi = motif_boot(test)
    p = motif_perm(test, learned, unseen, obs)

    return {
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "n_train_BID": int((train["y"] == 0).sum()),
        "n_train_ELAP": int((train["y"] == 1).sum()),
        "n_test_BID": int((test["y"] == 0).sum()),
        "n_test_ELAP": int((test["y"] == 1).sum()),
        "n_train_motifs": int(train["motif"].nunique()),
        "n_test_motifs": int(test["motif"].nunique()),
        "AUC": float(obs),
        "CI_low": float(lo),
        "CI_high": float(hi),
        "permutation_p": float(p),
    }


# ---------------------------------------------------------------------
# ELAP author-confidence repair
# ---------------------------------------------------------------------

def normalize_confidence(x) -> Tuple[str, float]:
    z = str(x).strip().lower()
    if z in {"", "nan", "none"}:
        return "", np.nan

    # Explicit labels: no inference from downstream results.
    if "highest" in z:
        return "highest", 3.0
    if "higher" in z:
        return "higher", 2.0
    if "lower" in z:
        return "lower", 1.0

    # Common compact spellings, still only if the semantic word is present.
    if re.search(r"\bhigh[-_ ]?est\b", z):
        return "highest", 3.0
    if re.search(r"\bhigh[-_ ]?er\b", z):
        return "higher", 2.0
    if re.search(r"\blow[-_ ]?er\b", z):
        return "lower", 1.0

    return f"unresolved:{z}", np.nan


def load_elap_confidence(cell: str):
    path = ELAP_FILES[cell]
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = {"base_locus_id", "confidence_raw"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"{path.name}: missing {sorted(missing)}")

    labels = df["confidence_raw"].map(normalize_confidence)
    df["confidence_class"] = [x[0] for x in labels]
    df["confidence_ordinal"] = [x[1] for x in labels]

    unresolved = sorted({
        x for x in df["confidence_class"].dropna().astype(str)
        if x.startswith("unresolved:")
    })
    if unresolved:
        raise RuntimeError(
            f"{cell}: unresolved ELAP confidence labels: {unresolved[:20]}"
        )

    # Optional numeric QC only; NOT used to choose thresholds.
    numeric_candidates = [
        c for c in df.columns
        if (
            "avg(stop_ratio_stopped_reads)" in c
            or "modification_estimation" in c
            or c == "signal_primary"
        )
    ]

    mapping = df.set_index("base_locus_id")[["confidence_class", "confidence_ordinal"]]
    qc = {
        "cell_line": cell,
        "n_source_rows": int(len(df)),
        "confidence_distribution": df["confidence_class"].value_counts(dropna=False).to_dict(),
        "n_confidence_resolved": int(df["confidence_ordinal"].notna().sum()),
        "numeric_candidate_columns_present": numeric_candidates,
    }
    return mapping, qc


def exact_locus_disjoint(site: pd.DataFrame):
    a = set(site.loc[site["cell_line"] == "HEK293T", "locus"])
    b = set(site.loc[site["cell_line"] == "HeLa", "locus"])
    overlap = a & b
    return site[~site["locus"].isin(overlap)].copy(), overlap


def attach_elap_confidence(site: pd.DataFrame):
    x = site.copy()
    qcs = []
    x["elap_confidence_class"] = ""
    x["elap_confidence_ordinal"] = np.nan

    for cell in ["HEK293T", "HeLa"]:
        mapping, qc = load_elap_confidence(cell)
        qcs.append(qc)

        mask = x["cell_line"].eq(cell)
        x.loc[mask, "elap_confidence_class"] = (
            x.loc[mask, "locus"].map(mapping["confidence_class"]).fillna("")
        )
        x.loc[mask, "elap_confidence_ordinal"] = (
            x.loc[mask, "locus"].map(mapping["confidence_ordinal"])
        )

        elap_reported = mask & x["elap_reported"].eq(1)
        frac = x.loc[elap_reported, "elap_confidence_ordinal"].notna().mean()
        qc["mapping_fraction_among_ELAP_reported_primary_loci"] = float(frac)
        qc["n_ELAP_reported_primary_loci"] = int(elap_reported.sum())

        if frac < 0.95:
            raise RuntimeError(
                f"{cell}: only {frac:.3f} of ELAP-reported Aim1 loci mapped "
                "to author confidence classes."
            )

    return x, qcs


def class_only(df):
    return df[
        df["class"].isin([BID_CLASS, ELAP_CLASS])
        & df["motif"].fillna("").astype(str).ne("")
    ].copy()


def confidence_sensitivity(disjoint: pd.DataFrame):
    """
    Primary concern: are abundant lower-confidence ELAP calls creating the
    fingerprint? BID published calls are kept unchanged because they already
    satisfy the published BID calling pipeline and >10% fraction criterion.

    Two predeclared ELAP restrictions:
      >= higher-confidence  (higher + highest)
      highest-confidence only
    """
    rows = []
    for threshold_name, min_ord in [
        ("ELAP_higher_or_highest_only", 2.0),
        ("ELAP_highest_only", 3.0),
    ]:
        z = class_only(disjoint).copy()
        keep = (
            z["class"].eq(BID_CLASS)
            |
            (
                z["class"].eq(ELAP_CLASS)
                & pd.to_numeric(z["elap_confidence_ordinal"], errors="coerce").ge(min_ord)
            )
        )
        z = z[keep].copy()

        for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
            tr = z[z["cell_line"] == train_cell]
            te = z[z["cell_line"] == test_cell]
            res = transfer(tr, te)
            res.update({
                "analysis": threshold_name,
                "train_cell": train_cell,
                "test_cell": test_cell,
                "ELAP_min_confidence_ordinal": min_ord,
            })
            rows.append(res)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# CDS-specific formal sensitivity
# ---------------------------------------------------------------------

def cds_sensitivity(disjoint: pd.DataFrame, ann: pd.DataFrame):
    x = disjoint.merge(
        ann[["cell_line", "locus", "mane_region"]].drop_duplicates(),
        on=["cell_line", "locus"],
        how="left",
    )
    x["mane_region"] = x["mane_region"].fillna("UNANNOTATED")
    rows = []

    for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
        # Train on all classed training-cell loci; evaluate only held-out CDS.
        tr = class_only(x[x["cell_line"] == train_cell])
        te = class_only(
            x[(x["cell_line"] == test_cell) & (x["mane_region"] == "CDS")]
        )
        res = transfer(tr, te)
        res.update({
            "analysis": "fixed_all_region_training_to_CDS_test",
            "train_cell": train_cell,
            "test_cell": test_cell,
        })
        rows.append(res)

        # Stricter: both training and test restricted to CDS.
        tr2 = class_only(
            x[(x["cell_line"] == train_cell) & (x["mane_region"] == "CDS")]
        )
        te2 = te
        res2 = transfer(tr2, te2)
        res2.update({
            "analysis": "CDS_to_CDS",
            "train_cell": train_cell,
            "test_cell": test_cell,
        })
        rows.append(res2)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Multiplicity-controlled positional effects
# ---------------------------------------------------------------------

def normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2))


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = min(prev, ranked[i] * n / rank)
        q[order[i]] = val
        prev = val
    return np.minimum(q, 1.0)


def positional_fdr(path: Path):
    df = pd.read_csv(path, sep="\t")
    rows = []
    for r in df.itertuples():
        a = r.ELAP_with_base + 0.5
        b = r.ELAP_without_base + 0.5
        c = r.BID_with_base + 0.5
        d = r.BID_without_base + 0.5
        se = math.sqrt(1/a + 1/b + 1/c + 1/d)
        z = r.logOR_ELAP_vs_BID / se
        rows.append({
            **r._asdict(),
            "SE_logOR": se,
            "z": z,
            "p": normal_two_sided_p(z),
        })
    out = pd.DataFrame(rows)

    # Adjust the 16 position/base tests separately in each cell.
    out["q_BH_within_cell"] = np.nan
    for cell, idx in out.groupby("cell_line").groups.items():
        out.loc[idx, "q_BH_within_cell"] = bh_adjust(out.loc[idx, "p"])

    a = out[out["cell_line"] == "HEK293T"][
        ["relative_position", "base", "logOR_ELAP_vs_BID", "q_BH_within_cell"]
    ].rename(columns={
        "logOR_ELAP_vs_BID": "logOR_HEK293T",
        "q_BH_within_cell": "q_HEK293T",
    })
    b = out[out["cell_line"] == "HeLa"][
        ["relative_position", "base", "logOR_ELAP_vs_BID", "q_BH_within_cell"]
    ].rename(columns={
        "logOR_ELAP_vs_BID": "logOR_HeLa",
        "q_BH_within_cell": "q_HeLa",
    })
    rep = a.merge(b, on=["relative_position", "base"], how="inner")
    rep["direction_consistent"] = (
        np.sign(rep["logOR_HEK293T"]) == np.sign(rep["logOR_HeLa"])
    )
    rep["FDR05_both_cells"] = (
        rep["q_HEK293T"].le(.05)
        & rep["q_HeLa"].le(.05)
        & rep["direction_consistent"]
    )
    rep["min_abs_logOR"] = rep[["logOR_HEK293T", "logOR_HeLa"]].abs().min(axis=1)
    rep = rep.sort_values(
        ["FDR05_both_cells", "min_abs_logOR"],
        ascending=[False, False],
    )
    return out, rep


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    print("=" * 98)
    print("TRACE PHASE 1D — AIM 1 CLOSURE")
    print("=" * 98)

    required = [SITE, ANN, POS, PHASE1C_SUMMARY, *ELAP_FILES.values()]
    for p in required:
        if not p.exists():
            raise RuntimeError(f"Missing required input: {p}")

    site = pd.read_csv(SITE, sep="\t", low_memory=False)
    ann = pd.read_csv(ANN, sep="\t", low_memory=False)

    site, conf_qc = attach_elap_confidence(site)
    disjoint, overlap = exact_locus_disjoint(site)

    conf = confidence_sensitivity(disjoint)
    conf.to_csv(
        AIM1 / "aim1d_corrected_confidence_sensitivity.tsv",
        sep="\t", index=False
    )

    cds = cds_sensitivity(disjoint, ann)
    cds.to_csv(
        AIM1 / "aim1d_CDS_sensitivity.tsv",
        sep="\t", index=False
    )

    pos_all, pos_rep = positional_fdr(POS)
    pos_all.to_csv(
        AIM1 / "aim1d_positional_effects_FDR.tsv",
        sep="\t", index=False
    )
    pos_rep.to_csv(
        AIM1 / "aim1d_replicated_positional_effects_FDR.tsv",
        sep="\t", index=False
    )

    phase1c = json.loads(PHASE1C_SUMMARY.read_text(encoding="utf-8-sig"))

    gene_ok = bool(
        phase1c["robustness_flags"]["gene_disjoint_both_CI_above_0.5"]
    )
    pus7_ok = bool(
        phase1c["robustness_flags"]["remove_PUS7_like_both_CI_above_0.5"]
    )

    higher = conf[conf["analysis"] == "ELAP_higher_or_highest_only"]
    higher_ok = bool(
        len(higher) == 2
        and higher["AUC"].notna().all()
        and higher["CI_low"].gt(.5).all()
    )

    highest = conf[conf["analysis"] == "ELAP_highest_only"]
    highest_direction_ok = bool(
        len(highest) == 2
        and highest["AUC"].notna().all()
        and highest["AUC"].gt(.5).all()
    )

    cds_to_cds = cds[cds["analysis"] == "CDS_to_CDS"]
    cds_ok = bool(
        len(cds_to_cds) == 2
        and cds_to_cds["AUC"].notna().all()
        and cds_to_cds["CI_low"].gt(.5).all()
    )

    # Highest-confidence-only can be underpowered; it is supportive rather than
    # required for closure. The author-defined higher+highest subset is the
    # primary confidence robustness gate.
    core_robust = gene_ok and pus7_ok and higher_ok and cds_ok

    status = (
        "AIM1_CLOSED_ROBUST_TECHNOLOGY_ASSOCIATED_FINGERPRINT"
        if core_robust
        else "AIM1_CLOSURE_REQUIRES_INTERPRETATION"
    )

    contract = {
        "phase": "1D",
        "role": "Aim 1 closure; no larger model added",
        "confidence_correction": {
            "Phase1C_problem": (
                "ELAP signal_primary was empty because Phase0C referenced the "
                "pre-cleaned column name avg(stop_ratio*stopped_reads); therefore "
                "the previous top50/top25 results were invalid and are superseded."
            ),
            "replacement_primary": (
                "Use ELAP author-defined higher/highest/lower confidence classes."
            ),
            "primary_confidence_sensitivity": (
                "Keep all published BID calls; remove ELAP lower-confidence calls "
                "in both training and test."
            ),
            "secondary_confidence_sensitivity": (
                "Keep all published BID calls; retain ELAP highest-confidence only."
            ),
        },
        "region_closure": (
            "Formal motif-cluster CI/permutation for held-out CDS; additionally "
            "CDS-to-CDS transfer."
        ),
        "positional_interpretation": (
            "16 base-by-position effects per cell; BH FDR within cell; replicated "
            "features require same direction and q<0.05 in both cells."
        ),
        "closure_gate": {
            "gene_disjoint_both_CI_above_0.5": "required",
            "remove_PUS7_like_both_CI_above_0.5": "required",
            "ELAP_higher_or_highest_both_CI_above_0.5": "required",
            "CDS_to_CDS_both_CI_above_0.5": "required",
            "ELAP_highest_only": "supportive; direction reported, not required",
        },
        "interpretation_boundary": (
            "A robust result supports a technology-associated sequence fingerprint "
            "that cannot be reduced to exact locus/gene reuse, broad transcript "
            "region, lower-confidence ELAP calls, or canonical PUS7-like motifs. "
            "It remains non-causal with respect to assay mechanism."
        ),
    }
    (META / "phase1d_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    summary = {
        "status": status,
        "n_exact_cross_cell_loci_removed": len(overlap),
        "ELAP_confidence_mapping_QC": conf_qc,
        "corrected_confidence_sensitivity": conf.to_dict(orient="records"),
        "CDS_sensitivity": cds.to_dict(orient="records"),
        "n_replicated_positional_features_FDR05_both_cells": int(
            pos_rep["FDR05_both_cells"].sum()
        ),
        "replicated_positional_features_FDR05": (
            pos_rep[pos_rep["FDR05_both_cells"]]
            .to_dict(orient="records")
        ),
        "closure_flags": {
            "gene_disjoint": gene_ok,
            "PUS7_removed": pus7_ok,
            "ELAP_higher_or_highest_confidence": higher_ok,
            "ELAP_highest_only_direction": highest_direction_ok,
            "CDS_to_CDS": cds_ok,
        },
        "scientific_interpretation": (
            "If closed: freeze Aim 1. The main result is a strong, reproducible "
            "BID-vs-ELAP technology-associated 5-mer fingerprint that transfers "
            "across cell types and persists under major biological/measurement "
            "stress tests; synthetic assay calibration explains only a minority "
            "of its variation."
        ),
        "next_step": (
            "If Aim1 is closed, do NOT add +/-10 nt elastic-net merely to improve "
            "prediction. Move to Aim2: chemistry-diverse support -> independent "
            "biological evidence using source-positive PUS perturbation and held-out DRS."
        ),
    }
    (META / "phase1d_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print()
    print("ELAP confidence mapping:")
    print(json.dumps(conf_qc, ensure_ascii=False, indent=2))
    print()
    print("CORRECTED CONFIDENCE SENSITIVITY")
    print(conf.to_string(index=False))
    print()
    print("CDS SENSITIVITY")
    print(cds.to_string(index=False))
    print()
    print("FDR-significant replicated positional features:")
    sig = pos_rep[pos_rep["FDR05_both_cells"]]
    print(sig.to_string(index=False) if len(sig) else "None")
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase1d_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase1d_analysis_contract.json")
    print(r"  D:\RNA\Trace\04_aim1\aim1d_corrected_confidence_sensitivity.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1d_CDS_sensitivity.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1d_replicated_positional_effects_FDR.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase1d_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase1d_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
