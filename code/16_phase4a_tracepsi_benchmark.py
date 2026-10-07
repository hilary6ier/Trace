#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE-Ψ Phase 4A — frozen LOTO evidence-portability benchmark
=============================================================

Purpose
-------
Construct the FINAL development benchmark for TRACE-Ψ without fitting the
portability model.

Scientific estimand
-------------------
For a locus already reported by >=1 source technology:

    Y_target = 1  if a held-out target technology also reports the locus
               0  if it is not present in that target technology's published
                  positive-call table.

Y_target=0 is NOT "unmodified" and NOT a biological negative. It means only
"not reported by that target assay/pipeline under the published experiment".

Primary locus universe
----------------------
Locally-single-U, exact genomic loci only.

This deliberately excludes consecutive-U ambiguity. The frozen HeLa universe
comes directly from Phase 2B:
    D:\RNA\Trace\05_aim2\aim2b_HeLa_source_union_DRS.tsv

Primary development technologies:
    BID, BACS, ELAP in HeLa

For each target technology t:
    candidate universe = loci reported by >=1 of the OTHER source technologies
    outcome            = whether t also reported that exact locus

Thus each target's evidence is unavailable to its own feature construction.

Feature hierarchy frozen here
-----------------------------
P0:
    source_breadth

P1:
    P0 +
    source_BID, source_BACS, source_ELAP

P2:
    P1 +
    four-position flanking-base one-hot features (-2,-1,+1,+2) +
    a frozen cross-cell BID-vs-ELAP 5-mer fingerprint learned ONLY in the
    opposite cell line (HEK293T -> HeLa development; HeLa -> HEK secondary) +
    fingerprint magnitude +
    source-alignment of the fingerprint

No source-strength/confidence feature is included in v1 because the published
assays expose non-equivalent signal definitions and prior project audits did
not establish a uniformly comparable confidence mapping for every technology.
This is a deliberate scientific restriction, not missing engineering.

Benchmarks created
------------------
1) HeLa BID/BACS/ELAP leave-one-technology-out task table.
2) HeLa native-DRS development benchmark (already frozen in Aim 2B).
3) HEK293T BID/ELAP cross-cell secondary transport task table.
4) Frozen cross-cell 5-mer fingerprint tables.

IMPORTANT LOCKBOX
-----------------
A549 DRS is NOT read, summarized, matched or scored in this phase.
It remains reserved for the later cell+platform lockbox validation.

Project root
------------
D:\RNA\Trace
"""

from __future__ import annotations

import json
import math
import re
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from sklearn.metrics import roc_auc_score
except Exception as e:
    raise RuntimeError("Phase4A requires scikit-learn.") from e


ROOT = Path(r"D:\RNA\Trace")
NORM = ROOT / "03_harmonized" / "normalized_sources"
AIM2 = ROOT / "05_aim2"
TRACEPSI = ROOT / "07_tracepsi"
META = ROOT / "00_meta"
LOGS = ROOT / "logs"

for p in [TRACEPSI, META, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

HELA_FROZEN = AIM2 / "aim2b_HeLa_source_union_DRS.tsv"

BID_HEK = NORM / "BID_HEK293T.tsv"
ELAP_HEK = NORM / "ELAP_HEK293T.tsv"
BID_HELA = NORM / "BID_HeLa.tsv"
ELAP_HELA = NORM / "ELAP_HeLa.tsv"

EXPECTED_HELA_UNION = 1331
EXPECTED_AIM1B = {
    "HEK293T": {"n": 623, "BID": 159, "ELAP": 464, "motifs": 124, "AUC_to_other": 0.7335020029517183},
    "HeLa":    {"n": 492, "BID": 186, "ELAP": 306, "motifs": 115, "AUC_to_other": 0.8203887443070917},
}
SMOOTHING_ALPHA = 0.5
BASES = "ACGU"
OFFSETS = [-2, -1, 1, 2]


# =============================================================================
# Generic helpers
# =============================================================================

def ss(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    return str(x).strip()


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
    return z if len(z) == 5 and z[2] == "U" else ""


def locally_single_u(m: str) -> bool:
    return bool(len(m) == 5 and m[2] == "U" and m[1] != "U" and m[3] != "U")


def exact_key(chrom: Any, strand: Any, pos1: Any) -> str:
    try:
        p = int(float(pos1))
    except Exception:
        return ""
    c = norm_chr(chrom)
    s = ss(strand)
    if not c or s not in {"+", "-"}:
        return ""
    return f"{c}:{s}:{p}"


def auc_safe(y, score) -> float:
    y = np.asarray(y, int)
    score = np.asarray(score, float)
    if len(np.unique(y)) != 2:
        return np.nan
    return float(roc_auc_score(y, score))


# =============================================================================
# Recover exact single-U BID / ELAP tables
# =============================================================================

def find_elap_seqcols(df: pd.DataFrame) -> List[str]:
    raw = [c for c in df.columns if c.startswith("raw__")]
    candidates = []
    for i in range(max(0, len(raw) - 4)):
        cols = raw[i:i+5]
        fracs = []
        for c in cols:
            z = df[c].fillna("").astype(str).str.upper().str.strip()
            fracs.append(z.isin(["A", "C", "G", "T", "U"]).mean())
        center = df[cols[2]].fillna("").astype(str).str.upper().str.strip()
        center_u = center.isin(["T", "U"]).mean()
        if min(fracs) > 0.90 and center_u > 0.90:
            candidates.append((float(np.mean(fracs) + center_u), cols))
    if not candidates:
        raise RuntimeError("Could not identify five ELAP sequence columns.")
    candidates.sort(reverse=True, key=lambda x: x[0])
    return list(candidates[0][1])


def load_bid_elap_single_u(path: Path, assay: str, cell: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)

    for c in ["chrom", "pos1", "strand"]:
        if c not in df.columns:
            raise RuntimeError(f"{assay}/{cell}: missing required column {c}")

    if assay == "BID":
        if "motif_5mer_reported" not in df.columns:
            raise RuntimeError(f"BID/{cell}: motif_5mer_reported missing.")
        motif = df["motif_5mer_reported"].map(motif5)
        motif_source = "motif_5mer_reported"
    elif assay == "ELAP":
        seqcols = find_elap_seqcols(df)
        motif = (
            df[seqcols]
            .fillna("")
            .astype(str)
            .agg("".join, axis=1)
            .map(motif5)
        )
        motif_source = "|".join(seqcols)
    else:
        raise ValueError(assay)

    x = pd.DataFrame({
        "chrom": df["chrom"].map(norm_chr),
        "pos1": pd.to_numeric(df["pos1"], errors="coerce"),
        "strand": df["strand"].astype(str).str.strip(),
        "motif": motif,
    })
    x = x[
        x["chrom"].ne("")
        & x["pos1"].notna()
        & x["strand"].isin(["+", "-"])
        & x["motif"].map(locally_single_u)
    ].copy()
    x["pos1"] = x["pos1"].astype(int)
    x["locus_id"] = [
        exact_key(c, s, p)
        for c, s, p in zip(x["chrom"], x["strand"], x["pos1"])
    ]
    x["assay"] = assay
    x["cell_line"] = cell

    dup = x.duplicated("locus_id", keep=False)
    if dup.any():
        # Exact duplicates must agree in motif; then collapse.
        bad = (
            x.loc[dup]
            .groupby("locus_id")["motif"]
            .nunique()
            .loc[lambda s: s > 1]
        )
        if len(bad):
            raise RuntimeError(
                f"{assay}/{cell}: duplicated loci disagree in motif: {bad.index[:10].tolist()}"
            )
        x = x.drop_duplicates("locus_id", keep="first").copy()

    x.attrs["motif_source"] = motif_source
    return x


def build_cell_bid_elap(cell: str) -> Tuple[pd.DataFrame, Dict]:
    if cell == "HEK293T":
        bid_path, elap_path = BID_HEK, ELAP_HEK
    elif cell == "HeLa":
        bid_path, elap_path = BID_HELA, ELAP_HELA
    else:
        raise ValueError(cell)

    bid = load_bid_elap_single_u(bid_path, "BID", cell)
    elap = load_bid_elap_single_u(elap_path, "ELAP", cell)

    b = bid[["locus_id", "chrom", "pos1", "strand", "motif"]].copy()
    b["BID_support"] = 1
    e = elap[["locus_id", "chrom", "pos1", "strand", "motif"]].copy()
    e["ELAP_support"] = 1

    u = b.merge(
        e[["locus_id", "motif", "ELAP_support"]],
        on="locus_id",
        how="outer",
        suffixes=("", "_ELAP"),
    )

    # Recover genomic fields for ELAP-only rows.
    e_index = e.set_index("locus_id")
    for c in ["chrom", "pos1", "strand"]:
        miss = u[c].isna()
        if miss.any():
            u.loc[miss, c] = u.loc[miss, "locus_id"].map(e_index[c])

    # Motif agreement gate on shared exact loci.
    shared = u["motif_ELAP"].notna() & u["motif"].notna()
    if shared.any():
        agree = (u.loc[shared, "motif"] == u.loc[shared, "motif_ELAP"]).mean()
        if agree < 0.999:
            bad = u.loc[shared & (u["motif"] != u["motif_ELAP"])]
            bad.to_csv(
                META / f"phase4a_{cell}_BID_ELAP_motif_disagreement.tsv",
                sep="\t",
                index=False,
            )
            raise RuntimeError(
                f"{cell}: BID/ELAP shared motif agreement={agree:.3%}; expected ~100%."
            )
    else:
        agree = np.nan

    u["motif"] = u["motif"].fillna(u["motif_ELAP"])
    u = u.drop(columns=["motif_ELAP"])
    u["BID_support"] = u["BID_support"].fillna(0).astype(int)
    u["ELAP_support"] = u["ELAP_support"].fillna(0).astype(int)
    u["cell_line"] = cell

    qc = {
        "cell_line": cell,
        "BID_singleU": int(len(bid)),
        "ELAP_singleU": int(len(elap)),
        "shared_exact": int(((u["BID_support"] == 1) & (u["ELAP_support"] == 1)).sum()),
        "union_exact": int(len(u)),
        "shared_motif_agreement": None if not np.isfinite(agree) else float(agree),
        "BID_motif_source": bid.attrs.get("motif_source", ""),
        "ELAP_motif_source": elap.attrs.get("motif_source", ""),
    }
    return u, qc


# =============================================================================
# Reconstruct the frozen Aim1B cross-cell 5-mer fingerprint
# =============================================================================

def aim1b_exclusive_after_crosscell_removal(
    hek_union: pd.DataFrame,
    hela_union: pd.DataFrame,
) -> Tuple[Dict[str, pd.DataFrame], Dict]:
    cross_keys = set(hek_union["locus_id"]) & set(hela_union["locus_id"])

    out = {}
    qc = {"cross_cell_exact_loci_removed": int(len(cross_keys))}

    for cell, u in [("HEK293T", hek_union), ("HeLa", hela_union)]:
        z = u[~u["locus_id"].isin(cross_keys)].copy()
        z = z[z["BID_support"] != z["ELAP_support"]].copy()
        z["label_ELAP"] = z["ELAP_support"].astype(int)
        z["class"] = np.where(z["label_ELAP"] == 1, "ELAP", "BID")
        out[cell] = z

        qc[cell] = {
            "n": int(len(z)),
            "BID": int((z["class"] == "BID").sum()),
            "ELAP": int((z["class"] == "ELAP").sum()),
            "motifs": int(z["motif"].nunique()),
        }

        exp = EXPECTED_AIM1B[cell]
        for key in ["n", "BID", "ELAP", "motifs"]:
            if qc[cell][key] != exp[key]:
                raise RuntimeError(
                    f"Aim1B reconstruction mismatch for {cell}/{key}: "
                    f"observed={qc[cell][key]}, expected={exp[key]}. "
                    "Do not build TRACE-Ψ fingerprint from a non-identical universe."
                )

    return out, qc


def fit_motif_logratio(train: pd.DataFrame, alpha: float = SMOOTHING_ALPHA) -> pd.DataFrame:
    """
    Positive score = ELAP-like; negative score = BID-like.

    Smoothed class-conditional log frequency ratio:
      log p(motif|ELAP) - log p(motif|BID)
    """
    motifs = sorted(set(train["motif"]))
    m = len(motifs)
    ne = int((train["class"] == "ELAP").sum())
    nb = int((train["class"] == "BID").sum())

    ce = Counter(train.loc[train["class"] == "ELAP", "motif"])
    cb = Counter(train.loc[train["class"] == "BID", "motif"])

    rows = []
    for motif in motifs:
        pe = (ce[motif] + alpha) / (ne + alpha * m)
        pb = (cb[motif] + alpha) / (nb + alpha * m)
        rows.append({
            "motif": motif,
            "n_ELAP": int(ce[motif]),
            "n_BID": int(cb[motif]),
            "fingerprint_logratio_ELAP_vs_BID": float(math.log(pe) - math.log(pb)),
        })
    return pd.DataFrame(rows)


def score_fingerprint(
    test: pd.DataFrame,
    fp: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray]:
    mp = fp.set_index("motif")["fingerprint_logratio_ELAP_vs_BID"].to_dict()
    score = test["motif"].map(mp)
    seen = score.notna().astype(int).to_numpy()
    return score.fillna(0.0).to_numpy(float), seen


def verify_crosscell_fingerprint(
    exclusives: Dict[str, pd.DataFrame],
) -> Tuple[Dict[str, pd.DataFrame], Dict]:
    fps = {
        "HEK293T": fit_motif_logratio(exclusives["HEK293T"]),
        "HeLa": fit_motif_logratio(exclusives["HeLa"]),
    }

    checks = []
    # Train HEK -> test HeLa; train HeLa -> test HEK.
    for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
        score, seen = score_fingerprint(exclusives[test_cell], fps[train_cell])
        auc = auc_safe(exclusives[test_cell]["label_ELAP"], score)

        expected = EXPECTED_AIM1B[train_cell]["AUC_to_other"]
        # Different but equivalent Laplace parameterizations can move AUC slightly
        # only through unseen-motif ties; require close reproduction, not bit identity.
        if not np.isfinite(auc) or abs(auc - expected) > 0.035:
            raise RuntimeError(
                f"Frozen Aim1B fingerprint reproduction failed {train_cell}->{test_cell}: "
                f"AUC={auc:.4f}, expected≈{expected:.4f}."
            )

        checks.append({
            "train_cell": train_cell,
            "test_cell": test_cell,
            "AUC": float(auc),
            "expected_Aim1B_AUC": float(expected),
            "abs_difference": float(abs(auc - expected)),
            "test_n": int(len(exclusives[test_cell])),
            "test_unseen_motif_n": int((seen == 0).sum()),
        })

    return fps, {"crosscell_checks": checks}


# =============================================================================
# Frozen TRACE-Ψ features
# =============================================================================

def add_sequence_features(
    df: pd.DataFrame,
    fp: pd.DataFrame,
) -> pd.DataFrame:
    z = df.copy()
    if not z["motif"].map(locally_single_u).all():
        raise RuntimeError("Phase4A primary table contains non-single-U motif.")

    fp_map = fp.set_index("motif")["fingerprint_logratio_ELAP_vs_BID"].to_dict()
    z["fp_crosscell_ELAPvsBID"] = z["motif"].map(fp_map)
    z["fp_seen_in_HEK"] = z["fp_crosscell_ELAPvsBID"].notna().astype(int)
    z["fp_crosscell_ELAPvsBID"] = z["fp_crosscell_ELAPvsBID"].fillna(0.0)
    z["fp_abs"] = z["fp_crosscell_ELAPvsBID"].abs()

    if not {"source_BID", "source_ELAP"}.issubset(z.columns):
        raise RuntimeError("source_BID/source_ELAP required before sequence features.")

    # Positive means motif is aligned with the positive source technology:
    # ELAP-like score for ELAP support; BID-like (negative fp) for BID support.
    z["fp_source_alignment"] = (
        z["fp_crosscell_ELAPvsBID"]
        * (z["source_ELAP"].astype(float) - z["source_BID"].astype(float))
    )

    # Four flanking positions. Center is always U and intentionally omitted.
    idx = {-2: 0, -1: 1, 1: 3, 2: 4}
    for off in OFFSETS:
        for base in BASES:
            z[f"motif_{off:+d}_{base}"] = (
                z["motif"].str[idx[off]].eq(base).astype(int)
            )

    return z


def make_loto_tasks(
    union: pd.DataFrame,
    cell_line: str,
    assays: Sequence[str],
    fp_for_cell: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for target in assays:
        target_col = f"{target}_support"
        if target_col not in union.columns:
            raise RuntimeError(f"{cell_line}: missing {target_col}")

        other = [a for a in assays if a != target]
        source_cols = [f"{a}_support" for a in other]

        cand = union[union[source_cols].sum(axis=1) >= 1].copy()
        cand["cell_line"] = cell_line
        cand["target_assay"] = target
        cand["target_reported"] = cand[target_col].astype(int)

        for a in ["BID", "BACS", "ELAP"]:
            if a == target:
                cand[f"source_{a}"] = 0
            elif f"{a}_support" in cand.columns:
                cand[f"source_{a}"] = cand[f"{a}_support"].astype(int)
            else:
                cand[f"source_{a}"] = 0

        cand["source_breadth"] = (
            cand["source_BID"] + cand["source_BACS"] + cand["source_ELAP"]
        )
        cand["source_pattern"] = cand.apply(
            lambda r: "+".join([
                a for a in ["BID", "BACS", "ELAP"]
                if int(r[f"source_{a}"]) == 1
            ]),
            axis=1,
        )

        keep = [
            "locus_id", "chrom", "pos1", "strand", "motif",
            "cell_line", "target_assay", "target_reported",
            "source_BID", "source_BACS", "source_ELAP",
            "source_breadth", "source_pattern",
        ]
        rows.append(cand[keep])

    out = pd.concat(rows, ignore_index=True)
    out = add_sequence_features(out, fp_for_cell)
    return out


def make_drs_benchmark(
    hela: pd.DataFrame,
    fp_hek: pd.DataFrame,
) -> pd.DataFrame:
    z = hela.copy()
    for a in ["BID", "BACS", "ELAP"]:
        z[f"source_{a}"] = z[f"{a}_support"].astype(int)
    z["source_breadth"] = (
        z["source_BID"] + z["source_BACS"] + z["source_ELAP"]
    )
    z["source_pattern"] = z.apply(
        lambda r: "+".join([
            a for a in ["BID", "BACS", "ELAP"]
            if int(r[f"source_{a}"]) == 1
        ]),
        axis=1,
    )
    z["cell_line"] = "HeLa"
    z["target_assay"] = "DRS"
    z["target_reported"] = z["DRS_reported_support"].astype(int)

    keep = [
        "locus_id", "chrom", "pos1", "strand", "motif",
        "cell_line", "target_assay", "target_reported",
        "source_BID", "source_BACS", "source_ELAP",
        "source_breadth", "source_pattern",
    ]
    out = add_sequence_features(z[keep], fp_hek)
    return out


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 104)
    print("TRACE-Ψ PHASE 4A — FROZEN EVIDENCE-PORTABILITY BENCHMARK")
    print("=" * 104)

    for p in [HELA_FROZEN, BID_HEK, ELAP_HEK, BID_HELA, ELAP_HELA]:
        if not p.exists():
            raise RuntimeError(f"Missing required input: {p}")

    # -------------------------------------------------------------------------
    # Frozen HeLa source universe from Aim2B
    # -------------------------------------------------------------------------
    hela = pd.read_csv(HELA_FROZEN, sep="\t", low_memory=False)

    required_hela = {
        "key", "chrom", "pos1", "strand", "motif",
        "BID_support", "BACS_support", "ELAP_support",
        "chemistry_breadth", "DRS_reported_support",
    }
    missing = sorted(required_hela - set(hela.columns))
    if missing:
        raise RuntimeError(f"Frozen HeLa Phase2B table missing: {missing}")

    if len(hela) != EXPECTED_HELA_UNION:
        raise RuntimeError(
            f"Frozen HeLa union rows={len(hela)}, expected {EXPECTED_HELA_UNION}."
        )

    hela["motif"] = hela["motif"].map(motif5)
    if not hela["motif"].map(locally_single_u).all():
        raise RuntimeError(
            "Frozen HeLa Phase2B table is no longer a pure locally-single-U universe."
        )

    hela["locus_id"] = [
        exact_key(c, s, p)
        for c, s, p in zip(hela["chrom"], hela["strand"], hela["pos1"])
    ]
    if hela["locus_id"].eq("").any() or hela["locus_id"].duplicated().any():
        raise RuntimeError("HeLa frozen union does not have unique exact locus IDs.")

    breadth_check = (
        hela[["BID_support", "BACS_support", "ELAP_support"]]
        .astype(int)
        .sum(axis=1)
    )
    if not np.array_equal(
        breadth_check.to_numpy(),
        pd.to_numeric(hela["chemistry_breadth"], errors="coerce").astype(int).to_numpy(),
    ):
        raise RuntimeError("Frozen chemistry_breadth disagrees with source support columns.")

    # -------------------------------------------------------------------------
    # Reconstruct exact BID/ELAP cross-cell universe and frozen fingerprint
    # -------------------------------------------------------------------------
    hek_union, hek_qc = build_cell_bid_elap("HEK293T")
    hela_be_union, hela_be_qc = build_cell_bid_elap("HeLa")

    exclusives, fp_universe_qc = aim1b_exclusive_after_crosscell_removal(
        hek_union, hela_be_union
    )
    fps, fp_qc = verify_crosscell_fingerprint(exclusives)

    fp_hek = fps["HEK293T"]  # primary HeLa feature: trained outside HeLa
    fp_hela = fps["HeLa"]    # secondary HEK feature: trained outside HEK

    fp_hek.to_csv(
        TRACEPSI / "phase4a_fingerprint_HEK_for_HeLa.tsv",
        sep="\t", index=False
    )
    fp_hela.to_csv(
        TRACEPSI / "phase4a_fingerprint_HeLa_for_HEK.tsv",
        sep="\t", index=False
    )

    # -------------------------------------------------------------------------
    # Primary HeLa LOTO table
    # -------------------------------------------------------------------------
    hela_loto = make_loto_tasks(
        union=hela,
        cell_line="HeLa",
        assays=["BID", "BACS", "ELAP"],
        fp_for_cell=fp_hek,
    )

    # Each (target,locus) row must be unique.
    if hela_loto.duplicated(["target_assay", "locus_id"]).any():
        raise RuntimeError("Duplicated (target,locus) rows in HeLa LOTO table.")

    # Every LOTO row must have source evidence.
    if (hela_loto["source_breadth"] < 1).any():
        raise RuntimeError("HeLa LOTO contains source_breadth < 1.")

    hela_loto.to_csv(
        TRACEPSI / "phase4a_HeLa_LOTO.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # -------------------------------------------------------------------------
    # HeLa DRS development benchmark; uses only the already-frozen Aim2B DRS label
    # and never opens the raw DRS table (protects A549 lockbox).
    # -------------------------------------------------------------------------
    drs_dev = make_drs_benchmark(hela, fp_hek)
    drs_dev.to_csv(
        TRACEPSI / "phase4a_HeLa_DRS_development.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # -------------------------------------------------------------------------
    # HEK cross-cell secondary transport table (BID <-> ELAP).
    # Fingerprint is trained only in HeLa.
    # -------------------------------------------------------------------------
    # Prepare support table in same schema.
    hek_for_tasks = hek_union.copy()
    hek_for_tasks["BACS_support"] = 0
    hek_tasks = make_loto_tasks(
        union=hek_for_tasks,
        cell_line="HEK293T",
        assays=["BID", "ELAP"],
        fp_for_cell=fp_hela,
    )
    hek_tasks.to_csv(
        TRACEPSI / "phase4a_HEK_crosscell_tasks.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # -------------------------------------------------------------------------
    # Frozen feature contract
    # -------------------------------------------------------------------------
    motif_features = [
        f"motif_{off:+d}_{base}"
        for off in OFFSETS for base in BASES
    ]
    feature_contract = {
        "P0_breadth": ["source_breadth"],
        "P1_provenance": [
            "source_breadth",
            "source_BID", "source_BACS", "source_ELAP",
        ],
        "P2_TRACEpsi": [
            "source_breadth",
            "source_BID", "source_BACS", "source_ELAP",
            "fp_crosscell_ELAPvsBID",
            "fp_abs",
            "fp_source_alignment",
            "fp_seen_in_HEK",
            *motif_features,
        ],
        "metadata_not_features": [
            "locus_id", "chrom", "pos1", "strand", "motif",
            "cell_line", "target_assay", "source_pattern",
        ],
        "outcome": "target_reported",
        "target_assay_is_feature": False,
        "source_confidence_in_primary_v1": False,
        "source_confidence_reason": (
            "Published assays expose heterogeneous, non-equivalent signal definitions; "
            "no uniformly comparable source-strength feature was frozen across BID/BACS/ELAP."
        ),
    }
    (META / "phase4a_feature_contract.json").write_text(
        json.dumps(feature_contract, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    # -------------------------------------------------------------------------
    # QC / contracts
    # -------------------------------------------------------------------------
    target_qc = []
    for target, z in hela_loto.groupby("target_assay"):
        target_qc.append({
            "target_assay": target,
            "n": int(len(z)),
            "target_reported_n": int(z["target_reported"].sum()),
            "target_reported_rate": float(z["target_reported"].mean()),
            "source_breadth_counts": {
                str(k): int(v)
                for k, v in z["source_breadth"].value_counts().sort_index().items()
            },
            "unique_loci": int(z["locus_id"].nunique()),
        })

    drs_qc = {
        "n": int(len(drs_dev)),
        "DRS_reported_n": int(drs_dev["target_reported"].sum()),
        "DRS_reported_rate": float(drs_dev["target_reported"].mean()),
        "breadth_counts": {
            str(k): int(v)
            for k, v in drs_dev["source_breadth"].value_counts().sort_index().items()
        },
        "breadth1_n": int((drs_dev["source_breadth"] == 1).sum()),
        "breadth1_DRS_reported_n": int(
            drs_dev.loc[drs_dev["source_breadth"] == 1, "target_reported"].sum()
        ),
    }

    if drs_qc["breadth1_n"] < 800:
        raise RuntimeError(
            f"DRS breadth=1 challenge unexpectedly small: {drs_qc['breadth1_n']}"
        )

    contract = {
        "phase": "4A",
        "framework": "TRACE-Ψ / leave-one-technology-out evidence portability",
        "estimand": (
            "Probability/ranking of re-observation by a held-out profiling technology "
            "conditional on prior source-positive evidence. It is NOT P(true Ψ)."
        ),
        "primary_universe": (
            "Frozen Phase2B HeLa locally-single-U exact genomic source union "
            "(BID/BACS/ELAP)."
        ),
        "development_targets": ["BID", "BACS", "ELAP"],
        "development_external_platform": (
            "HeLa DRS labels already frozen in Aim2B; used only as development benchmark."
        ),
        "secondary_cross_cell": "HEK293T BID<->ELAP task table",
        "feature_hierarchy": feature_contract,
        "fingerprint_policy": (
            "Primary HeLa fingerprint is trained only in HEK293T after reproducing the "
            "frozen Aim1B cross-cell-disjoint universe; HEK secondary uses HeLa-trained fingerprint."
        ),
        "coordinate_policy": (
            "Exact locally-single-U loci only; no arbitrary +/-k matching; consecutive-U "
            "ambiguity is excluded from the primary portability benchmark."
        ),
        "noncall_interpretation": (
            "target_reported=0 means absent from target published positive-call table, "
            "not biological absence."
        ),
        "A549_DRS_lockbox": {
            "status": "UNTOUCHED_IN_PHASE4A",
            "A549_DRS_read": False,
            "A549_DRS_matched": False,
            "A549_DRS_metrics_computed": False,
        },
    }
    (META / "phase4a_data_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    summary = {
        "status": "PHASE4A_PORTABILITY_BENCHMARK_FROZEN",
        "HeLa_frozen_union_n": int(len(hela)),
        "HeLa_BID_ELAP_reconstruction": hela_be_qc,
        "HEK_BID_ELAP_reconstruction": hek_qc,
        "Aim1B_reconstruction": fp_universe_qc,
        "fingerprint_QC": fp_qc,
        "HeLa_LOTO_rows": int(len(hela_loto)),
        "HeLa_LOTO_target_QC": target_qc,
        "HeLa_DRS_development_QC": drs_qc,
        "HEK_secondary_rows": int(len(hek_tasks)),
        "HEK_secondary_target_QC": [
            {
                "target_assay": t,
                "n": int(len(z)),
                "target_reported_n": int(z["target_reported"].sum()),
                "target_reported_rate": float(z["target_reported"].mean()),
            }
            for t, z in hek_tasks.groupby("target_assay")
        ],
        "feature_contract_file": str(META / "phase4a_feature_contract.json"),
        "next_gate": (
            "Run Phase4B exactly once: target-held-out P0/P1/P2 + HeLa DRS development "
            "benchmark + breadth=1 challenge. Do not inspect A549 DRS."
        ),
    }
    (META / "phase4a_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )

    print()
    print("STATUS: PHASE4A_PORTABILITY_BENCHMARK_FROZEN")
    print("HeLa frozen loci:", len(hela))
    print("HeLa LOTO rows:", len(hela_loto))
    print("HeLa DRS benchmark rows:", len(drs_dev))
    print("HeLa DRS breadth=1:", drs_qc["breadth1_n"])
    print("HEK secondary task rows:", len(hek_tasks))
    print()
    print("TARGET QC")
    print(pd.DataFrame(target_qc).to_string(index=False))
    print()
    print("FINGERPRINT QC")
    print(pd.DataFrame(fp_qc["crosscell_checks"]).to_string(index=False))
    print()
    print(r"Return/upload after Phase4A:")
    print(r"  D:\RNA\Trace\00_meta\phase4a_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase4a_data_contract.json")
    print(r"  D:\RNA\Trace\00_meta\phase4a_feature_contract.json")
    print(r"  D:\RNA\Trace\07_tracepsi\phase4a_HeLa_LOTO.tsv.gz")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase4a_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase4a_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
