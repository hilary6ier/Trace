#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE Phase 1C — explain and stress-test the cross-cell fingerprint
===================================================================

Scientific purpose
------------------
Aim 1B established a strong cross-cell BID-vs-ELAP 5-mer fingerprint.
Phase 1C asks whether that result can be reduced to obvious biological or
measurement-selection explanations, WITHOUT increasing model complexity.

Predeclared questions
1) Gene identity leakage:
   Does transfer persist when held-out loci come from genes never represented
   in the training cell?
2) Transcript-region composition:
   Can CDS/UTR composition alone explain the assay label, and does the fixed
   5-mer score still transfer within major transcript regions?
3) Source-signal strength:
   Does transfer persist among assay-internally stronger published calls?
4) Known writer sequence biology:
   Does transfer persist after removing canonical PUS7-like UNUAR contexts?

Important boundaries
--------------------
- Primary universe remains Aim1A/1B locally single-U exact loci.
- No raw FASTQ/BAM.
- No new sequence model or feature tuning.
- No missing published call is interpreted as a biological negative.
- MANE annotation is used only for covariates/sensitivity analysis.
- Source signal is ranked WITHIN assay/cell; raw BID and ELAP units are never
  compared to each other.

Project root:
    D:\RNA\Trace

Inputs
------
    04_aim1\aim1a_site_table.tsv
    02_reference\MANE.GRCh38.v1.5.refseq_genomic.gtf.gz
    03_harmonized\normalized_sources\BID_HEK293T.tsv
    03_harmonized\normalized_sources\BID_HeLa.tsv
    03_harmonized\normalized_sources\ELAP_HEK293T.tsv
    03_harmonized\normalized_sources\ELAP_HeLa.tsv

Outputs
-------
    00_meta\phase1c_analysis_contract.json
    00_meta\phase1c_summary.json
    04_aim1\aim1c_locus_annotations.tsv
    04_aim1\aim1c_gene_disjoint.tsv
    04_aim1\aim1c_region_analysis.tsv
    04_aim1\aim1c_confidence_sensitivity.tsv
    04_aim1\aim1c_writer_motif_sensitivity.tsv
    04_aim1\aim1c_positional_base_effects.tsv
"""

from __future__ import annotations

import gzip
import json
import math
import re
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import numpy as np
import pandas as pd

ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
REF = ROOT / "02_reference"
NORM = ROOT / "03_harmonized" / "normalized_sources"
AIM1 = ROOT / "04_aim1"
LOGS = ROOT / "logs"

for p in [META, AIM1, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

SITE = AIM1 / "aim1a_site_table.tsv"
GTF = REF / "MANE.GRCh38.v1.5.refseq_genomic.gtf.gz"

SOURCE_TABLES = {
    ("HEK293T", "BID"): NORM / "BID_HEK293T.tsv",
    ("HeLa", "BID"): NORM / "BID_HeLa.tsv",
    ("HEK293T", "ELAP"): NORM / "ELAP_HEK293T.tsv",
    ("HeLa", "ELAP"): NORM / "ELAP_HeLa.tsv",
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

# No tuning: these are declared before this script reads outcome results.
SIGNAL_QUANTILES = [0.50, 0.75]

# Human PUS7 canonical degenerate motif: UNΨAR / UNUAR.
PUS7_RE = re.compile(r"^U[ACGU]UA[AG]$")
# A stricter motif highlighted in 2026 work (USUAG; S=C/G).
PUS7_STRICT_RE = re.compile(r"^U[CG]UAG$")
# Canonical TRUB1 local 5-mer around target from RGUΨCN... is GUΨCN -> GUUCN.
TRUB1_5MER_RE = re.compile(r"^GUUC[ACGU]$")


# ---------------------------------------------------------------------
# Core utilities
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


def learn_motif_score(train: pd.DataFrame, alpha: float = ALPHA):
    n_bid = int((train["y"] == 0).sum())
    n_elap = int((train["y"] == 1).sum())
    ct = train.groupby(["motif", "y"]).size().unstack(fill_value=0)
    if 0 not in ct.columns:
        ct[0] = 0
    if 1 not in ct.columns:
        ct[1] = 0

    score = {}
    for motif in ct.index:
        nb = int(ct.loc[motif, 0])
        ne = int(ct.loc[motif, 1])
        score[motif] = (
            math.log((ne + alpha) / (n_elap + alpha * MOTIF_UNIVERSE))
            - math.log((nb + alpha) / (n_bid + alpha * MOTIF_UNIVERSE))
        )

    unseen = (
        math.log(alpha / (n_elap + alpha * MOTIF_UNIVERSE))
        - math.log(alpha / (n_bid + alpha * MOTIF_UNIVERSE))
    )
    return score, unseen


def motif_cluster_bootstrap(test: pd.DataFrame, n_boot: int = N_BOOT):
    motifs = test["motif"].dropna().unique()
    groups = {m: test[test["motif"] == m] for m in motifs}
    vals = []
    if len(motifs) < 5:
        return float("nan"), float("nan")
    for _ in range(n_boot):
        sampled = RNG.choice(motifs, size=len(motifs), replace=True)
        z = pd.concat([groups[m] for m in sampled], ignore_index=True)
        a = auc_score(z["y"], z["score"])
        if np.isfinite(a):
            vals.append(a)
    if not vals:
        return float("nan"), float("nan")
    return float(np.quantile(vals, .025)), float(np.quantile(vals, .975))


def motif_permutation_p(test: pd.DataFrame, learned: Dict[str, float], unseen: float,
                        obs: float, n_perm: int = N_PERM) -> float:
    motifs = list(learned)
    values = np.asarray([learned[m] for m in motifs], float)
    row_motifs = test["motif"].to_numpy()
    y = test["y"].to_numpy()
    ge = 0
    for _ in range(n_perm):
        pm = dict(zip(motifs, RNG.permutation(values)))
        score = np.asarray([pm.get(m, unseen) for m in row_motifs], float)
        if auc_score(y, score) >= obs:
            ge += 1
    return float((ge + 1) / (n_perm + 1))


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
    lo, hi = motif_cluster_bootstrap(test)
    p = motif_permutation_p(test, learned, unseen, obs)
    return {
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "n_train_BID": int((train["y"] == 0).sum()),
        "n_train_ELAP": int((train["y"] == 1).sum()),
        "n_test_BID": int((test["y"] == 0).sum()),
        "n_test_ELAP": int((test["y"] == 1).sum()),
        "n_train_motifs": int(train["motif"].nunique()),
        "n_test_motifs": int(test["motif"].nunique()),
        "n_test_unseen_motif_loci": int((~test["motif"].isin(learned)).sum()),
        "AUC": float(obs),
        "CI_low": float(lo),
        "CI_high": float(hi),
        "permutation_p": float(p),
    }


# ---------------------------------------------------------------------
# Locus / MANE annotation
# ---------------------------------------------------------------------

def parse_locus(x: str):
    # locus = chr:strand:pos ; chromosome itself does not contain ":"
    m = re.match(r"^([^:]+):([+-]):(\d+)$", str(x))
    if not m:
        return "", "", np.nan
    return m.group(1), m.group(2), int(m.group(3))


def parse_gtf_attrs(text: str) -> Dict[str, str]:
    out = {}
    for key, val in re.findall(r'(\S+)\s+"([^"]*)"', text):
        out[key] = val
    return out


def region_label(feature: str) -> str:
    f = feature.lower()
    if f == "cds" or f in {"start_codon", "stop_codon"}:
        return "CDS"
    if "five_prime" in f and "utr" in f:
        return "5UTR"
    if "three_prime" in f and "utr" in f:
        return "3UTR"
    if f == "utr":
        return "UTR"
    if f == "exon":
        return "EXON_OTHER"
    return ""


def annotate_with_mane(site: pd.DataFrame) -> pd.DataFrame:
    q = site[["locus"]].drop_duplicates().copy()
    parsed = q["locus"].map(parse_locus)
    q["chrom"] = [x[0] for x in parsed]
    q["strand"] = [x[1] for x in parsed]
    q["pos1"] = [x[2] for x in parsed]

    wanted_groups = set(zip(q["chrom"], q["strand"]))
    intervals = defaultdict(list)

    opener = gzip.open if str(GTF).endswith(".gz") else open
    with opener(GTF, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = parts
            if (chrom, strand) not in wanted_groups:
                continue
            reg = region_label(feature)
            if not reg:
                continue
            a = parse_gtf_attrs(attrs)
            gene_id = a.get("gene_id", "")
            gene_name = a.get("gene_name", a.get("gene", gene_id))
            transcript_id = a.get("transcript_id", "")
            intervals[(chrom, strand)].append(
                (int(start), int(end), gene_id, gene_name, transcript_id, reg)
            )

    ann_rows = []
    for key, subq in q.groupby(["chrom", "strand"], sort=False):
        ints = sorted(intervals.get(key, []), key=lambda z: z[0])
        positions = sorted(
            [(int(r.pos1), r.locus) for r in subq.itertuples() if np.isfinite(r.pos1)]
        )
        active = []
        j = 0
        for pos, locus in positions:
            while j < len(ints) and ints[j][0] <= pos:
                active.append(ints[j])
                j += 1
            active = [iv for iv in active if iv[1] >= pos]
            hits = [iv for iv in active if iv[0] <= pos <= iv[1]]

            genes = sorted({h[3] or h[2] for h in hits if h[3] or h[2]})
            gene_ids = sorted({h[2] for h in hits if h[2]})
            txs = sorted({h[4] for h in hits if h[4]})
            regs = {h[5] for h in hits if h[5]}

            # MANE should usually resolve one transcript/gene. If feature labels
            # overlap, use a transparent priority but retain all raw regions.
            if "CDS" in regs:
                primary = "CDS"
            elif "5UTR" in regs and "3UTR" not in regs:
                primary = "5UTR"
            elif "3UTR" in regs and "5UTR" not in regs:
                primary = "3UTR"
            elif "UTR" in regs:
                primary = "UTR"
            elif "EXON_OTHER" in regs:
                primary = "EXON_OTHER"
            elif regs:
                primary = "MULTI"
            else:
                primary = "UNANNOTATED"

            ann_rows.append({
                "locus": locus,
                "mane_gene": "|".join(genes),
                "mane_gene_id": "|".join(gene_ids),
                "mane_transcript_id": "|".join(txs),
                "mane_n_genes": len(genes),
                "mane_regions_all": "|".join(sorted(regs)),
                "mane_region": primary,
            })

    ann = pd.DataFrame(ann_rows)
    return q.merge(ann, on="locus", how="left")


# ---------------------------------------------------------------------
# Source signal ranks
# ---------------------------------------------------------------------

def load_signal_map(path: Path) -> Dict[str, float]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "base_locus_id" not in df.columns or "signal_primary" not in df.columns:
        raise RuntimeError(f"{path.name}: base_locus_id/signal_primary missing")
    df["signal_num"] = pd.to_numeric(df["signal_primary"], errors="coerce")
    if df["base_locus_id"].duplicated().any():
        # Exact duplicates are not expected in these normalized BID/ELAP maps.
        d = int(df["base_locus_id"].duplicated().sum())
        raise RuntimeError(f"{path.name}: duplicate base_locus_id rows={d}")
    return df.set_index("base_locus_id")["signal_num"].to_dict()


def attach_signal_ranks(site: pd.DataFrame) -> Tuple[pd.DataFrame, Dict]:
    x = site.copy()
    report = {}
    for cell in ["HEK293T", "HeLa"]:
        for assay in ["BID", "ELAP"]:
            mp = load_signal_map(SOURCE_TABLES[(cell, assay)])
            col = f"{assay.lower()}_signal"
            mask = x["cell_line"].eq(cell)
            x.loc[mask, col] = x.loc[mask, "locus"].map(mp)

            vals = pd.to_numeric(x.loc[mask, col], errors="coerce")
            ranks = vals.rank(method="average", pct=True)
            rank_col = f"{assay.lower()}_signal_rank"
            x.loc[mask, rank_col] = ranks

            report[f"{cell}_{assay}"] = {
                "n_site_rows": int(mask.sum()),
                "n_signal_mapped": int(vals.notna().sum()),
                "signal_mapping_fraction": float(vals.notna().mean()),
            }
    return x, report


# ---------------------------------------------------------------------
# Analyses
# ---------------------------------------------------------------------

def exact_locus_disjoint(site: pd.DataFrame) -> Tuple[pd.DataFrame, int]:
    a = set(site.loc[site["cell_line"] == "HEK293T", "locus"])
    b = set(site.loc[site["cell_line"] == "HeLa", "locus"])
    overlap = a & b
    return site[~site["locus"].isin(overlap)].copy(), len(overlap)


def class_only(df: pd.DataFrame) -> pd.DataFrame:
    return df[
        df["class"].isin([BID_CLASS, ELAP_CLASS])
        & df["motif"].fillna("").ne("")
    ].copy()


def run_gene_disjoint(disjoint: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
        train = class_only(disjoint[disjoint["cell_line"] == train_cell])
        test = class_only(disjoint[disjoint["cell_line"] == test_cell])

        # Strict sensitivity: require a unique MANE gene in both sets.
        train_u = train[train["mane_n_genes"] == 1].copy()
        test_u = test[test["mane_n_genes"] == 1].copy()

        train_genes = set(train_u["mane_gene"])
        strict_test = test_u[~test_u["mane_gene"].isin(train_genes)].copy()

        res = transfer(train_u, strict_test)
        res.update({
            "train_cell": train_cell,
            "test_cell": test_cell,
            "analysis": "heldout_test_genes_absent_from_training",
            "n_train_unique_gene_loci": len(train_u),
            "n_test_unique_gene_loci_before_filter": len(test_u),
            "n_test_gene_novel_loci": len(strict_test),
            "n_training_genes": len(train_genes),
            "n_test_gene_novel_genes": strict_test["mane_gene"].nunique(),
        })
        rows.append(res)
    return pd.DataFrame(rows)


def learn_region_score(train: pd.DataFrame):
    train = train.copy()
    train["y"] = (train["class"] == ELAP_CLASS).astype(int)
    n0 = (train["y"] == 0).sum()
    n1 = (train["y"] == 1).sum()
    cats = ["5UTR", "CDS", "3UTR", "UTR", "EXON_OTHER", "UNANNOTATED"]
    score = {}
    for r in cats:
        a = ((train["mane_region"] == r) & (train["y"] == 1)).sum()
        b = ((train["mane_region"] == r) & (train["y"] == 0)).sum()
        score[r] = (
            math.log((a + .5) / (n1 + .5 * len(cats)))
            - math.log((b + .5) / (n0 + .5 * len(cats)))
        )
    return score


def run_region_analysis(disjoint: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
        train = class_only(disjoint[disjoint["cell_line"] == train_cell])
        test = class_only(disjoint[disjoint["cell_line"] == test_cell])
        test["y"] = (test["class"] == ELAP_CLASS).astype(int)

        # Region-only cross-cell discrimination.
        rscore = learn_region_score(train)
        test["region_score"] = test["mane_region"].map(rscore).fillna(0)
        rows.append({
            "train_cell": train_cell,
            "test_cell": test_cell,
            "analysis": "region_only",
            "region": "ALL",
            "n_test": len(test),
            "n_BID": int((test["y"] == 0).sum()),
            "n_ELAP": int((test["y"] == 1).sum()),
            "AUC": auc_score(test["y"], test["region_score"]),
        })

        # Fixed 5-mer score learned on all training loci, then evaluated within
        # each held-out region. No refitting within region.
        train["y"] = (train["class"] == ELAP_CLASS).astype(int)
        learned, unseen = learn_motif_score(train)
        test["motif_score"] = test["motif"].map(learned).fillna(unseen)

        for region in ["CDS", "3UTR", "5UTR", "EXON_OTHER", "UNANNOTATED"]:
            z = test[test["mane_region"] == region]
            if len(z) < 20 or z["y"].nunique() < 2:
                continue
            rows.append({
                "train_cell": train_cell,
                "test_cell": test_cell,
                "analysis": "fixed_5mer_within_region",
                "region": region,
                "n_test": len(z),
                "n_BID": int((z["y"] == 0).sum()),
                "n_ELAP": int((z["y"] == 1).sum()),
                "AUC": auc_score(z["y"], z["motif_score"]),
            })
    return pd.DataFrame(rows)


def run_confidence_sensitivity(disjoint: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for q in SIGNAL_QUANTILES:
        for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
            train = class_only(disjoint[disjoint["cell_line"] == train_cell]).copy()
            test = class_only(disjoint[disjoint["cell_line"] == test_cell]).copy()

            def strong_mask(d):
                return (
                    ((d["class"] == BID_CLASS) & (pd.to_numeric(d["bid_signal_rank"], errors="coerce") >= q))
                    |
                    ((d["class"] == ELAP_CLASS) & (pd.to_numeric(d["elap_signal_rank"], errors="coerce") >= q))
                )

            tr = train[strong_mask(train)].copy()
            te = test[strong_mask(test)].copy()
            res = transfer(tr, te)
            res.update({
                "train_cell": train_cell,
                "test_cell": test_cell,
                "analysis": f"assay_internal_signal_rank_ge_{q:.2f}",
                "rank_threshold": q,
            })
            rows.append(res)
    return pd.DataFrame(rows)


def motif_flags(df: pd.DataFrame) -> pd.DataFrame:
    x = df.copy()
    x["pus7_like_UNUAR"] = x["motif"].astype(str).map(lambda m: bool(PUS7_RE.match(m)))
    x["pus7_strict_USUAG"] = x["motif"].astype(str).map(lambda m: bool(PUS7_STRICT_RE.match(m)))
    x["trub1_local_GUUCN"] = x["motif"].astype(str).map(lambda m: bool(TRUB1_5MER_RE.match(m)))
    x["plus1U"] = x["motif"].astype(str).str.len().eq(5) & x["motif"].astype(str).str[3].eq("U")
    x["minus1U"] = x["motif"].astype(str).str.len().eq(5) & x["motif"].astype(str).str[1].eq("U")
    return x


def run_writer_sensitivity(disjoint: pd.DataFrame) -> Tuple[pd.DataFrame, Dict]:
    x = motif_flags(disjoint)
    diagnostic = {
        "n_rows": len(x),
        "n_pus7_like_UNUAR": int(x["pus7_like_UNUAR"].sum()),
        "n_pus7_strict_USUAG": int(x["pus7_strict_USUAG"].sum()),
        "n_trub1_local_GUUCN": int(x["trub1_local_GUUCN"].sum()),
        "n_plus1U": int(x["plus1U"].sum()),
        "n_minus1U": int(x["minus1U"].sum()),
        "scientific_note": (
            "The primary single-U universe should contain zero +/-1U contexts, "
            "therefore it structurally excludes ELAP +1U preference, consecutive-U "
            "localization ambiguity, and the canonical local TRUB1 GUΨCN context."
        ),
    }

    rows = []
    for removal_name, removal_col in [
        ("remove_PUS7_like_UNUAR", "pus7_like_UNUAR"),
        ("remove_strict_PUS7_USUAG", "pus7_strict_USUAG"),
    ]:
        z = x[~x[removal_col]].copy()
        for train_cell, test_cell in [("HEK293T", "HeLa"), ("HeLa", "HEK293T")]:
            tr = class_only(z[z["cell_line"] == train_cell])
            te = class_only(z[z["cell_line"] == test_cell])
            res = transfer(tr, te)
            res.update({
                "train_cell": train_cell,
                "test_cell": test_cell,
                "analysis": removal_name,
                "n_removed_total": int(x[removal_col].sum()),
            })
            rows.append(res)
    return pd.DataFrame(rows), diagnostic


def positional_effects(disjoint: pd.DataFrame) -> pd.DataFrame:
    """
    Interpretable effect sizes, not a model:
    for each cell, position and base, log OR of ELAP-exclusive vs BID-exclusive.
    """
    rows = []
    positions = {-2: 0, -1: 1, +1: 3, +2: 4}
    for cell in ["HEK293T", "HeLa"]:
        z = class_only(disjoint[disjoint["cell_line"] == cell]).copy()
        for rel, idx in positions.items():
            for base in "ACGU":
                is_base = z["motif"].str[idx].eq(base)
                is_elap = z["class"].eq(ELAP_CLASS)
                a = int((is_base & is_elap).sum())
                b = int((~is_base & is_elap).sum())
                c = int((is_base & ~is_elap).sum())
                d = int((~is_base & ~is_elap).sum())
                # Haldane-Anscombe correction.
                aa, bb, cc, dd = a+.5, b+.5, c+.5, d+.5
                logor = math.log((aa*dd)/(bb*cc))
                se = math.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
                rows.append({
                    "cell_line": cell,
                    "relative_position": rel,
                    "base": base,
                    "ELAP_with_base": a,
                    "ELAP_without_base": b,
                    "BID_with_base": c,
                    "BID_without_base": d,
                    "logOR_ELAP_vs_BID": logor,
                    "OR": math.exp(logor),
                    "CI_low": math.exp(logor - 1.96*se),
                    "CI_high": math.exp(logor + 1.96*se),
                })

    out = pd.DataFrame(rows)
    wide = out.pivot_table(
        index=["relative_position", "base"],
        columns="cell_line",
        values="logOR_ELAP_vs_BID",
    ).reset_index()
    if {"HEK293T", "HeLa"}.issubset(wide.columns):
        wide["direction_replicates"] = np.sign(wide["HEK293T"]) == np.sign(wide["HeLa"])
        wide["min_abs_logOR"] = wide[["HEK293T", "HeLa"]].abs().min(axis=1)
        rep = wide[wide["direction_replicates"]].sort_values("min_abs_logOR", ascending=False)
        rep.to_csv(AIM1 / "aim1c_replicated_positional_features.tsv", sep="\t", index=False)
    return out


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------

def main():
    print("="*98)
    print("TRACE PHASE 1C — EXPLAIN / STRESS-TEST CROSS-CELL FINGERPRINT")
    print("="*98)

    for p in [SITE, GTF, *SOURCE_TABLES.values()]:
        if not p.exists():
            raise RuntimeError(f"Missing required input: {p}")

    site = pd.read_csv(SITE, sep="\t", low_memory=False)
    site["motif"] = site["motif"].fillna("").astype(str)

    ann = annotate_with_mane(site)
    site = site.merge(
        ann[["locus","mane_gene","mane_gene_id","mane_transcript_id",
             "mane_n_genes","mane_regions_all","mane_region"]],
        on="locus", how="left"
    )
    site["mane_n_genes"] = pd.to_numeric(site["mane_n_genes"], errors="coerce").fillna(0).astype(int)
    site["mane_region"] = site["mane_region"].fillna("UNANNOTATED")
    site["mane_gene"] = site["mane_gene"].fillna("")

    site, signal_report = attach_signal_ranks(site)
    disjoint, n_overlap = exact_locus_disjoint(site)

    # Annotation QC
    annotation_qc = {
        "n_site_rows": int(len(site)),
        "n_unique_loci": int(site["locus"].nunique()),
        "n_exact_cross_cell_loci_removed": int(n_overlap),
        "mane_gene_annotation_fraction": float(site["mane_n_genes"].gt(0).mean()),
        "mane_unique_gene_fraction": float(site["mane_n_genes"].eq(1).mean()),
        "mane_region_counts": site["mane_region"].value_counts().to_dict(),
        "source_signal_mapping": signal_report,
    }

    ann_out = site[[
        "cell_line","locus","class","motif",
        "mane_gene","mane_gene_id","mane_transcript_id","mane_n_genes",
        "mane_regions_all","mane_region",
        "bid_signal","bid_signal_rank","elap_signal","elap_signal_rank"
    ]].copy()
    ann_out.to_csv(AIM1 / "aim1c_locus_annotations.tsv", sep="\t", index=False)

    gene = run_gene_disjoint(disjoint)
    gene.to_csv(AIM1 / "aim1c_gene_disjoint.tsv", sep="\t", index=False)

    region = run_region_analysis(disjoint)
    region.to_csv(AIM1 / "aim1c_region_analysis.tsv", sep="\t", index=False)

    conf = run_confidence_sensitivity(disjoint)
    conf.to_csv(AIM1 / "aim1c_confidence_sensitivity.tsv", sep="\t", index=False)

    writer, writer_diag = run_writer_sensitivity(disjoint)
    writer.to_csv(AIM1 / "aim1c_writer_motif_sensitivity.tsv", sep="\t", index=False)

    pos = positional_effects(disjoint)
    pos.to_csv(AIM1 / "aim1c_positional_base_effects.tsv", sep="\t", index=False)

    # Gate interpretation is deliberately descriptive, not a hidden score.
    gene_informative = gene["n_test"].ge(50).all() and gene["AUC"].notna().all()
    gene_both_above = gene_informative and gene["CI_low"].gt(.5).all()

    high50 = conf[conf["rank_threshold"] == .50]
    conf50_both_above = len(high50) == 2 and high50["CI_low"].gt(.5).all()

    pus7b = writer[writer["analysis"] == "remove_PUS7_like_UNUAR"]
    pus7_both_above = len(pus7b) == 2 and pus7b["CI_low"].gt(.5).all()

    region_only = region[region["analysis"] == "region_only"]
    region_only_auc = {
        f"{r.train_cell}->{r.test_cell}": float(r.AUC)
        for r in region_only.itertuples()
    }

    within = region[region["analysis"] == "fixed_5mer_within_region"]
    major_within = within[within["n_test"] >= 50].copy()
    within_major_above_half = (
        len(major_within) > 0 and major_within["AUC"].gt(.5).all()
    )

    if gene_both_above and conf50_both_above and pus7_both_above:
        status = "AIM1C_FINGERPRINT_ROBUST_TO_MAJOR_EXPLANATIONS"
    else:
        status = "AIM1C_RESULTS_READY_FOR_INTERPRETATION"

    contract = {
        "phase": "1C",
        "purpose": (
            "Explain/stress-test the strong Aim1B cross-cell fingerprint without "
            "increasing model complexity."
        ),
        "frozen_primary_universe": "Aim1A locally single-U exact loci",
        "tests": {
            "gene_disjoint": (
                "Train on uniquely MANE-annotated training loci; evaluate only held-out "
                "loci whose MANE gene is absent from training."
            ),
            "region": (
                "Region-only cross-cell score plus fixed 5-mer score evaluated within "
                "held-out transcript regions."
            ),
            "confidence": (
                "Repeat transfer after retaining calls at or above 50% and 75% "
                "within-assay/cell source-signal ranks."
            ),
            "writer_sequence": (
                "Repeat after removing PUS7-like UNUAR and strict USUAG contexts."
            ),
            "positional_effects": (
                "Report -2/-1/+1/+2 base log-odds effects independently in HEK293T/HeLa."
            ),
        },
        "known_mechanism_exclusion_by_design": (
            "Single-U primary loci contain no +/-1 U, so the main Aim1B result cannot "
            "be directly driven by ELAP's +1U preference, consecutive-U localization "
            "ambiguity, or canonical local TRUB1 GUΨCN context."
        ),
        "no_model_selection": True,
    }
    (META / "phase1c_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    summary = {
        "status": status,
        "annotation_QC": annotation_qc,
        "known_motif_diagnostic": writer_diag,
        "gene_disjoint": gene.to_dict(orient="records"),
        "region_only_AUC": region_only_auc,
        "major_within_region_results": major_within.to_dict(orient="records"),
        "confidence_sensitivity": conf.to_dict(orient="records"),
        "writer_motif_sensitivity": writer.to_dict(orient="records"),
        "robustness_flags": {
            "gene_disjoint_both_CI_above_0.5": bool(gene_both_above),
            "top50_signal_both_CI_above_0.5": bool(conf50_both_above),
            "remove_PUS7_like_both_CI_above_0.5": bool(pus7_both_above),
            "all_major_within_region_AUC_above_0.5": bool(within_major_above_half),
        },
        "next_gate": (
            "If the fingerprint remains robust, freeze Aim1 as a completed result and "
            "move to Aim2 orthogonal biological evidence (PUS perturbation + held-out DRS). "
            "Only add GRCh38 +/-10-nt elastic-net if Phase1C exposes an unresolved sequence "
            "mechanism that 5-mer cannot localize; do not add it merely to increase AUROC."
        ),
    }
    (META / "phase1c_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print("Annotation QC:", json.dumps(annotation_qc, ensure_ascii=False))
    print()
    print("Known motif diagnostic:", writer_diag)
    print()
    print("GENE-DISJOINT")
    print(gene.to_string(index=False))
    print()
    print("REGION")
    print(region.to_string(index=False))
    print()
    print("SOURCE-SIGNAL SENSITIVITY")
    print(conf.to_string(index=False))
    print()
    print("WRITER-MOTIF SENSITIVITY")
    print(writer.to_string(index=False))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase1c_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase1c_analysis_contract.json")
    print(r"  D:\RNA\Trace\04_aim1\aim1c_gene_disjoint.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1c_region_analysis.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1c_confidence_sensitivity.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1c_writer_motif_sensitivity.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1c_positional_base_effects.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1c_replicated_positional_features.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase1c_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase1c_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
