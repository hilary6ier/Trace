#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE Phase 2A v5 — direct PRAISE perturbation validation
==========================================================

Why v4
------
The previous BID Fig.3k route is retired for the PRIMARY Aim-2 analysis:
Fig.3k is a gene-labelled visualization and is not a uniquely reversible
nucleotide-level outcome table.

PRAISE provides a cleaner site-level design:
- Supplementary Dataset 2: 2,209 PRAISE-reported HEK293T records
  (2,201 unique published chr_site loci).
- Supplementary Dataset 3: direct PUS-dependent site-level records with an
  `enzyme_dependency` field (1,021 records; 1,017 unique chr_site loci).

Scientific question
-------------------
Within a fixed PRAISE-positive HEK293T universe, does support from an
orthogonal chemistry (ELAP) carry independent biological evidence of
PUS dependency, beyond:
  1) PRAISE's own signal strength, and
  2) support from BID, which belongs to the SAME bisulfite/deletion family?

This directly tests the project's principle:
    orthogonal chemistry != repeated support by related chemistry.

Primary model
-------------
writer_assigned ~ ELAP_support + BID_support
                  + source_signal_rank + source_interval_flag

- ELAP_support = independent enzymatic-labeling / RT-stop chemistry.
- BID_support  = same bisulfite/deletion family as PRAISE.
- source_signal_rank = within-PRAISE rank of mean replicate deletion ratio.
- source_interval_flag controls PRAISE's published consecutive-U ambiguity.

SEs are gene-cluster robust.

Important interpretation
------------------------
writer_assigned=0 means only:
    "not listed in PRAISE Supplementary Dataset 3 as PUS1/PUS7/TRUB1/DKC1-
     dependent under the published perturbation experiment."
It is NOT a biological negative and does NOT imply the site is unmodified.

Matching
--------
PRAISE preserves its published `chr_site` semantics:
- point:    chrN_position
- interval: chrN_start-end

An external exact BID/ELAP genomic call supports PRAISE if:
- point source: exact same genomic position
- interval source: external call lies within the published PRAISE interval

No arbitrary +/-k matching is used.
PRAISE does not provide strand in Dataset 2; therefore:
- primary support is coordinate/interval based;
- a sensitivity analysis excludes source loci where an external assay has
  matching calls on >1 genomic strand within the source interval.

No FASTQ/BAM/raw sequencing.
No truth probability.
No outcome-derived threshold tuning.
"""

from __future__ import annotations

import bisect
import json
import math
import re
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
NORM = ROOT / "03_harmonized" / "normalized_sources"
AIM2 = ROOT / "05_aim2"
LOGS = ROOT / "logs"

for p in [META, AIM2, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

PRAISE_SITES = NORM / "PRAISE_HEK293T.tsv"
PRAISE_WRITERS = NORM / "PRAISE_writer_assignments.tsv"
BID = NORM / "BID_HEK293T.tsv"
ELAP = NORM / "ELAP_HEK293T.tsv"

EXPECTED = {
    "source_records": 2209,
    "source_unique_chr_site": 2201,
    "writer_records": 1021,
    "writer_unique_chr_site": 1017,
    "writer_counts": {
        "DKC1": 473,
        "TRUB1": 346,
        "PUS7": 165,
        "PUS1": 37,
    },
}

SEED = 20261003
RNG = np.random.default_rng(SEED)

PRIMARY_PREDICTORS = [
    "ELAP_support",
    "BID_support",
    "source_signal_rank",
    "source_interval_flag",
]


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------

def s(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    return str(x).strip()


def num(x: Any) -> float:
    try:
        z = float(str(x).replace("%", "").replace(",", ""))
        return z if np.isfinite(z) else np.nan
    except Exception:
        return np.nan


def norm_chr(x: Any) -> str:
    z = s(x)
    if not z:
        return ""
    if z.lower().startswith("chr"):
        return "chr" + z[3:]
    if re.fullmatch(r"(?:\d+|X|Y|M|MT)", z, flags=re.I):
        return "chr" + ("M" if z.upper() == "MT" else z)
    return z


def wilson(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if n <= 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z*z/n
    ctr = (p + z*z/(2*n)) / den
    rad = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return ctr-rad, ctr+rad


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    if len(vals) == 0:
        return q
    order = np.argsort(vals)
    ranked = vals[order]
    out = np.empty(len(vals))
    prev = 1.0
    n = len(vals)
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = min(prev, ranked[i] * n / rank)
        out[order[i]] = val
        prev = val
    q[np.where(ok)[0]] = np.minimum(out, 1.0)
    return q


# -----------------------------------------------------------------------------
# PRAISE source / outcome
# -----------------------------------------------------------------------------

# PRAISE chr_site is author-provided text. We parse from the RIGHT edge so
# chromosome/contig names may themselves contain underscores or version tokens.
_COORD_TAIL_PATTERNS = [
    # chr..._12345-12347 ; chr...:12345-12347
    re.compile(
        r"^(?P<chrom>chr.+?)[_:](?P<a>\d+)"
        r"(?:\s*[-–—~/:;,]\s*(?P<b>\d+))?$",
        flags=re.I
    ),
    # chr..._12345_12347 (interval encoded with a second underscore)
    re.compile(
        r"^(?P<chrom>chr.+?)[_:](?P<a>\d+)[_:](?P<b>\d+)$",
        flags=re.I
    ),
]


def parse_chr_site(x: Any):
    """
    Conservative PRAISE coordinate parser.

    Returns:
        chrom, start, end, interval_flag, parser_mode

    The parser never invents a coordinate. If no documented-looking tail can
    be recovered, the locus remains UNKNOWN and is excluded from coordinate-
    based evidence analyses rather than being treated as unsupported.
    """
    z0 = s(x)
    if not z0:
        return "", np.nan, np.nan, np.nan, "missing"

    z = (
        z0.replace("−", "-")
          .replace("–", "-")
          .replace("—", "-")
          .strip()
    )
    # Remove harmless outer quotes/brackets only.
    z = z.strip("\"'[]() ")

    for i, pat in enumerate(_COORD_TAIL_PATTERNS, start=1):
        m = pat.match(z)
        if not m:
            continue
        chrom = norm_chr(m.group("chrom"))
        a = int(m.group("a"))
        btxt = m.groupdict().get("b")
        b = int(btxt) if btxt else a
        lo, hi = min(a, b), max(a, b)
        return chrom, lo, hi, int(lo != hi), f"tail_pattern_{i}"

    # Last conservative fallback: split at the LAST '_' or ':' before an
    # all-numeric coordinate tail. This handles contig names containing many
    # underscores without consuming their version digits.
    m = re.match(
        r"^(?P<chrom>chr.+)[_:](?P<tail>\d+(?:[-–—~/:;,_]\d+)?)$",
        z,
        flags=re.I
    )
    if m:
        chrom = norm_chr(m.group("chrom"))
        nums = [int(v) for v in re.findall(r"\d+", m.group("tail"))]
        if len(nums) in {1, 2}:
            a = nums[0]
            b = nums[1] if len(nums) == 2 else a
            lo, hi = min(a, b), max(a, b)
            return chrom, lo, hi, int(lo != hi), "right_tail_fallback"

    return "", np.nan, np.nan, np.nan, "unparsed"


def require_cols(df: pd.DataFrame, cols: List[str], label: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"{label}: missing required columns {missing}")


def build_praise_source() -> Tuple[pd.DataFrame, Dict]:
    df = pd.read_csv(PRAISE_SITES, sep="\t", low_memory=False)

    require_cols(
        df,
        [
            "praise_chr_site", "gene",
            "raw__rep1_deletion_ratio",
            "raw__rep2_deletion_ratio",
        ],
        "PRAISE source"
    )

    if len(df) != EXPECTED["source_records"]:
        raise RuntimeError(
            f"PRAISE source rows={len(df)}, expected {EXPECTED['source_records']}."
        )

    df["praise_chr_site"] = df["praise_chr_site"].astype(str).str.strip()
    unique_n = df["praise_chr_site"].nunique()
    if unique_n != EXPECTED["source_unique_chr_site"]:
        raise RuntimeError(
            f"PRAISE source unique chr_site={unique_n}, "
            f"expected {EXPECTED['source_unique_chr_site']}."
        )

    parsed = df["praise_chr_site"].map(parse_chr_site)
    df["src_chrom"] = [x[0] for x in parsed]
    df["src_start"] = [x[1] for x in parsed]
    df["src_end"] = [x[2] for x in parsed]
    df["source_interval_flag"] = [x[3] for x in parsed]
    df["coordinate_parser_mode"] = [x[4] for x in parsed]
    df["coordinate_usable"] = (
        df["src_chrom"].ne("")
        & df["src_start"].notna()
        & df["src_end"].notna()
    )

    # Always write the parser audit; no hidden attrition.
    parser_audit = (
        df[[
            "praise_chr_site", "gene", "src_chrom", "src_start", "src_end",
            "source_interval_flag", "coordinate_parser_mode", "coordinate_usable"
        ]]
        .drop_duplicates("praise_chr_site")
        .copy()
    )
    parser_audit.to_csv(
        META / "phase2a_v5_praise_coordinate_audit.tsv",
        sep="\t", index=False
    )

    unique_parse_fraction = float(parser_audit["coordinate_usable"].mean())
    bad_unique = parser_audit[~parser_audit["coordinate_usable"]].copy()
    bad_unique.to_csv(
        META / "phase2a_v5_praise_unparsed_loci.tsv",
        sep="\t", index=False
    )

    # Reasonable manuscript-level gate:
    # tiny irreducible source-format attrition is allowed, but never silently.
    if unique_parse_fraction < 0.99:
        raise RuntimeError(
            f"Only {unique_parse_fraction:.3%} of unique PRAISE chr_site loci "
            "have safe genomic coordinates. This exceeds the predeclared 1% "
            "maximum coordinate attrition. "
            r"See 00_meta\phase2a_v5_praise_unparsed_loci.tsv"
        )

    r1 = pd.to_numeric(df["raw__rep1_deletion_ratio"], errors="coerce")
    r2 = pd.to_numeric(df["raw__rep2_deletion_ratio"], errors="coerce")
    df["source_signal"] = pd.concat([r1, r2], axis=1).mean(axis=1, skipna=True)

    spearman = float(r1.corr(r2, method="spearman"))
    if not np.isfinite(spearman) or spearman < 0.85:
        raise RuntimeError(
            f"PRAISE replicate deletion-ratio Spearman={spearman}; "
            "expected strong agreement from prior audit (~0.902)."
        )

    # Collapse transcript-level duplicate records to the published genomic locus.
    # Coordinate usability is explicit; unparsed loci remain in the source
    # universe for outcome bookkeeping but are not eligible for external matching.
    src = (
        df.groupby("praise_chr_site", as_index=False)
        .agg(
            src_chrom=("src_chrom", "first"),
            src_start=("src_start", "first"),
            src_end=("src_end", "first"),
            source_interval_flag=("source_interval_flag", "first"),
            coordinate_parser_mode=("coordinate_parser_mode", "first"),
            coordinate_usable=("coordinate_usable", "max"),
            source_signal=("source_signal", "mean"),
            gene=("gene", lambda z: "|".join(sorted({
                s(v) for v in z if s(v)
            }))),
            n_source_records=("praise_chr_site", "size"),
        )
    )
    src["source_signal_rank"] = src["source_signal"].rank(
        method="average", pct=True
    )

    qc = {
        "records": len(df),
        "unique_chr_site": len(src),
        "replicate_spearman": spearman,
        "coordinate_usable_unique_loci": int(src["coordinate_usable"].sum()),
        "coordinate_unusable_unique_loci": int((~src["coordinate_usable"]).sum()),
        "coordinate_parse_fraction_unique": unique_parse_fraction,
        "parser_mode_counts": (
            parser_audit["coordinate_parser_mode"].value_counts().to_dict()
        ),
        "interval_loci_among_usable": int(
            pd.to_numeric(
                src.loc[src["coordinate_usable"], "source_interval_flag"],
                errors="coerce"
            ).fillna(0).sum()
        ),
        "point_loci_among_usable": int(
            (
                pd.to_numeric(
                    src.loc[src["coordinate_usable"], "source_interval_flag"],
                    errors="coerce"
                ).fillna(-1) == 0
            ).sum()
        ),
        "duplicated_record_loci": int((src["n_source_records"] > 1).sum()),
    }
    return src, qc


def build_writer_outcome() -> Tuple[pd.DataFrame, Dict]:
    df = pd.read_csv(PRAISE_WRITERS, sep="\t", low_memory=False)
    require_cols(
        df,
        ["raw__chr_site", "raw__enzyme_dependency"],
        "PRAISE writer assignments"
    )

    if len(df) != EXPECTED["writer_records"]:
        raise RuntimeError(
            f"PRAISE writer rows={len(df)}, expected {EXPECTED['writer_records']}."
        )

    df["chr_site"] = df["raw__chr_site"].astype(str).str.strip()
    df["writer"] = (
        df["raw__enzyme_dependency"].astype(str).str.strip().str.upper()
    )

    counts = df["writer"].value_counts().to_dict()
    expected_counts = EXPECTED["writer_counts"]
    for w, n in expected_counts.items():
        if int(counts.get(w, 0)) != n:
            raise RuntimeError(
                f"PRAISE writer {w}: observed {counts.get(w,0)}, expected {n}."
            )

    unexpected = sorted(set(df["writer"]) - set(expected_counts))
    if unexpected:
        raise RuntimeError(f"Unexpected PRAISE writer labels: {unexpected}")

    unique_sites = df["chr_site"].nunique()
    if unique_sites != EXPECTED["writer_unique_chr_site"]:
        raise RuntimeError(
            f"PRAISE writer unique chr_site={unique_sites}, "
            f"expected {EXPECTED['writer_unique_chr_site']}."
        )

    out = (
        df.groupby("chr_site", as_index=False)
        .agg(
            writer_set=("writer", lambda z: ",".join(sorted(set(z)))),
            n_writer_assignments=("writer", "nunique"),
        )
    )
    out["writer_assigned"] = 1

    qc = {
        "records": len(df),
        "unique_chr_site": unique_sites,
        "writer_counts": {w: int(counts[w]) for w in expected_counts},
        "multi_writer_sites": int((out["n_writer_assignments"] > 1).sum()),
    }
    return out, qc


# -----------------------------------------------------------------------------
# External exact-call indexes
# -----------------------------------------------------------------------------

def build_external_index(path: Path, assay: str) -> Dict[str, Dict]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    require_cols(df, ["chrom", "strand", "pos1"], assay)

    z = pd.DataFrame({
        "chrom": df["chrom"].map(norm_chr),
        "strand": df["strand"].astype(str).str.strip(),
        "pos1": pd.to_numeric(df["pos1"], errors="coerce"),
    }).dropna(subset=["pos1"])
    z = z[z["chrom"].ne("") & z["strand"].isin(["+", "-"])].copy()
    z["pos1"] = z["pos1"].astype(int)

    by_chr = {}
    for chrom, g in z.groupby("chrom", sort=False):
        g = g.sort_values("pos1")
        positions = g["pos1"].astype(int).tolist()
        strands = g["strand"].astype(str).tolist()
        by_chr[chrom] = {"positions": positions, "strands": strands}

    return {
        "by_chr": by_chr,
        "n_calls": len(z),
        "n_unique_exact": len(set(zip(z["chrom"], z["strand"], z["pos1"]))),
    }


def interval_support(chrom: str, lo: int, hi: int, index: Dict) -> Tuple[int, int, str]:
    """
    O(log n + k) interval query using bisect.
    Returns support, strand ambiguity, and matching strand label.
    Source PRAISE has no strand, so strand ambiguity is explicitly recorded.
    """
    block = index["by_chr"].get(chrom)
    if not block:
        return 0, 0, ""
    positions = block["positions"]
    strands_all = block["strands"]
    left = bisect.bisect_left(positions, lo)
    right = bisect.bisect_right(positions, hi)
    if left >= right:
        return 0, 0, ""
    strands = sorted(set(strands_all[left:right]))
    return 1, int(len(strands) > 1), "|".join(strands)


def attach_external_support(src: pd.DataFrame, index: Dict, prefix: str) -> pd.DataFrame:
    """
    External support is undefined (NaN), not zero, when the PRAISE coordinate
    itself cannot be safely parsed.
    """
    x = src.copy()
    support = []
    ambig = []
    strands = []
    for r in x.itertuples():
        usable = bool(getattr(r, "coordinate_usable", False))
        if (
            not usable
            or pd.isna(r.src_start)
            or pd.isna(r.src_end)
            or not s(r.src_chrom)
        ):
            support.append(np.nan)
            ambig.append(np.nan)
            strands.append("")
            continue

        a, b, c = interval_support(
            r.src_chrom, int(r.src_start), int(r.src_end), index
        )
        support.append(float(a))
        ambig.append(float(b))
        strands.append(c)

    x[f"{prefix}_support"] = support
    x[f"{prefix}_strand_ambiguous"] = ambig
    x[f"{prefix}_matching_strands"] = strands
    return x


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------

def logistic_cluster(df: pd.DataFrame, outcome: str, predictors: List[str],
                     cluster: str = "gene"):
    d = df[[outcome, cluster] + predictors].copy()
    for c in [outcome] + predictors:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d[cluster] = df.loc[d.index, cluster].astype(str)
    d = d.dropna()

    y = d[outcome].to_numpy(float)
    X = np.column_stack([
        np.ones(len(d)),
        *[d[p].to_numpy(float) for p in predictors],
    ])
    names = ["intercept"] + predictors

    if len(np.unique(y)) < 2:
        raise RuntimeError(f"{outcome} has only one class in model universe.")

    beta = np.zeros(X.shape[1])
    for _ in range(200):
        eta = np.clip(X @ beta, -30, 30)
        pr = 1 / (1 + np.exp(-eta))
        W = np.clip(pr * (1-pr), 1e-8, None)
        H = X.T @ (W[:, None] * X)
        score = X.T @ (y-pr)
        step = np.linalg.pinv(H) @ score
        beta2 = beta + step
        if np.max(np.abs(step)) < 1e-9:
            beta = beta2
            break
        beta = beta2

    eta = np.clip(X @ beta, -30, 30)
    pr = 1 / (1 + np.exp(-eta))
    W = np.clip(pr * (1-pr), 1e-8, None)
    bread = np.linalg.pinv(X.T @ (W[:, None] * X))

    groups = defaultdict(list)
    for i, g in enumerate(d[cluster].astype(str)):
        groups[g].append(i)

    meat = np.zeros((X.shape[1], X.shape[1]))
    for ids in groups.values():
        sg = X[ids].T @ (y[ids] - pr[ids])
        meat += np.outer(sg, sg)

    G = len(groups)
    N, K = X.shape
    correction = (G / max(G-1, 1)) * ((N-1) / max(N-K, 1))
    cov = bread @ meat @ bread * correction
    se = np.sqrt(np.clip(np.diag(cov), 0, None))

    rows = []
    for i, (name, b, se0) in enumerate(zip(names, beta, se)):
        z = b / se0 if se0 > 0 else np.nan
        p = math.erfc(abs(z)/math.sqrt(2)) if np.isfinite(z) else np.nan
        rows.append({
            "term": name,
            "beta": float(b),
            "SE_gene_cluster": float(se0),
            "z": float(z) if np.isfinite(z) else np.nan,
            "p": float(p) if np.isfinite(p) else np.nan,
            "OR": math.exp(float(b)),
            "OR_CI_low": math.exp(float(b - 1.96*se0)),
            "OR_CI_high": math.exp(float(b + 1.96*se0)),
            "n": N,
            "n_gene_clusters": G,
        })

    return pd.DataFrame(rows), cov, names


def coefficient_contrast(table: pd.DataFrame, cov: np.ndarray, names: List[str],
                         term_a: str, term_b: str) -> Dict:
    ia, ib = names.index(term_a), names.index(term_b)
    beta_a = float(table.loc[table["term"] == term_a, "beta"].iloc[0])
    beta_b = float(table.loc[table["term"] == term_b, "beta"].iloc[0])
    diff = beta_a - beta_b
    var = cov[ia,ia] + cov[ib,ib] - 2*cov[ia,ib]
    se = math.sqrt(max(var, 0))
    z = diff/se if se > 0 else np.nan
    p = math.erfc(abs(z)/math.sqrt(2)) if np.isfinite(z) else np.nan
    return {
        "contrast": f"{term_a}_minus_{term_b}",
        "beta_difference": diff,
        "SE": se,
        "z": z,
        "p": p,
        "ratio_of_ORs": math.exp(diff),
        "CI_low": math.exp(diff - 1.96*se),
        "CI_high": math.exp(diff + 1.96*se),
    }


def rate_table(df: pd.DataFrame, group: str) -> pd.DataFrame:
    rows = []
    for label, z in df.groupby(group, dropna=False):
        n = len(z)
        k = int(z["writer_assigned"].sum())
        lo, hi = wilson(k, n)
        rows.append({
            group: label,
            "n": n,
            "writer_assigned_n": k,
            "writer_assignment_rate": k/n if n else np.nan,
            "wilson95_low": lo,
            "wilson95_high": hi,
        })
    return pd.DataFrame(rows)


def writer_specific_models(df: pd.DataFrame) -> pd.DataFrame:
    """
    Exploratory-but-predeclared writer-specific analysis.
    BH correction applies to ELAP_support tests across four writers.
    """
    rows = []
    for writer in ["DKC1", "TRUB1", "PUS7", "PUS1"]:
        z = df.copy()
        z["outcome"] = z["writer_set"].fillna("").str.split(",").map(
            lambda xs: int(writer in xs)
        )
        if z["outcome"].sum() < 20:
            continue
        tab, _, _ = logistic_cluster(
            z, "outcome",
            ["ELAP_support", "BID_support", "source_signal_rank", "source_interval_flag"]
        )
        effect = tab[tab["term"] == "ELAP_support"].iloc[0].to_dict()
        effect["writer"] = writer
        rows.append(effect)

    out = pd.DataFrame(rows)
    if len(out):
        out["q_BH_ELAP_across_writers"] = bh_adjust(out["p"])
    return out


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    print("=" * 100)
    print("TRACE PHASE 2A v5 — PRAISE DIRECT SITE-LEVEL PUS VALIDATION")
    print("=" * 100)

    for p in [PRAISE_SITES, PRAISE_WRITERS, BID, ELAP]:
        if not p.exists():
            raise RuntimeError(f"Missing required input: {p}")

    source, source_qc = build_praise_source()
    writers, writer_qc = build_writer_outcome()

    source = source.merge(
        writers[["chr_site", "writer_set", "n_writer_assignments", "writer_assigned"]],
        left_on="praise_chr_site", right_on="chr_site", how="left"
    )
    source["writer_assigned"] = source["writer_assigned"].fillna(0).astype(int)
    source["writer_set"] = source["writer_set"].fillna("")
    source["n_writer_assignments"] = (
        source["n_writer_assignments"].fillna(0).astype(int)
    )
    source = source.drop(columns=["chr_site"])

    overlap_n = int(source["writer_assigned"].sum())
    if overlap_n != EXPECTED["writer_unique_chr_site"]:
        missing = writers[~writers["chr_site"].isin(set(source["praise_chr_site"]))]
        missing.to_csv(
            META / "phase2a_v5_writer_sites_not_in_source.tsv",
            sep="\t", index=False
        )
        raise RuntimeError(
            f"Writer-dependent unique sites mapped to source={overlap_n}, "
            f"expected {EXPECTED['writer_unique_chr_site']}. "
            r"See 00_meta\phase2a_v5_writer_sites_not_in_source.tsv"
        )

    # Coordinate-unusable loci retain their biological-outcome annotation but
    # are excluded from external-support inference. They are NOT coded as zero support.
    full_source = source.copy()

    bid_idx = build_external_index(BID, "BID_HEK293T")
    elap_idx = build_external_index(ELAP, "ELAP_HEK293T")

    source = attach_external_support(source, bid_idx, "BID")
    source = attach_external_support(source, elap_idx, "ELAP")

    source["analysis_included"] = (
        source["coordinate_usable"].astype(bool)
        & source["BID_support"].notna()
        & source["ELAP_support"].notna()
        & source["source_signal_rank"].notna()
        & source["gene"].fillna("").astype(str).ne("")
    )

    analysis = source[source["analysis_included"]].copy()
    if len(analysis) < 0.98 * len(source):
        raise RuntimeError(
            f"Aim2A analyzable source universe fell to {len(analysis)}/{len(source)} "
            "(>2% total attrition). See coordinate audit before proceeding."
        )

    source["support_pattern"] = np.where(
        source["analysis_included"],
        np.select(
        [
            (source["BID_support"] == 0) & (source["ELAP_support"] == 0),
            (source["BID_support"] == 1) & (source["ELAP_support"] == 0),
            (source["BID_support"] == 0) & (source["ELAP_support"] == 1),
            (source["BID_support"] == 1) & (source["ELAP_support"] == 1),
        ],
        ["PRAISE_only", "PRAISE+BID_same_family",
         "PRAISE+ELAP_cross_chemistry", "PRAISE+BID+ELAP"],
            default="unexpected",
        ),
        "coordinate_unknown_excluded",
    )
    source["assay_breadth"] = np.where(
        source["analysis_included"],
        1 + source["BID_support"] + source["ELAP_support"],
        np.nan,
    )
    source["chemistry_breadth"] = np.where(
        source["analysis_included"],
        1 + source["ELAP_support"],
        np.nan,
    )

    # Main descriptive result uses only the explicitly analyzable universe.
    pattern_rates = rate_table(analysis.assign(
        support_pattern=np.select(
            [
                (analysis["BID_support"] == 0) & (analysis["ELAP_support"] == 0),
                (analysis["BID_support"] == 1) & (analysis["ELAP_support"] == 0),
                (analysis["BID_support"] == 0) & (analysis["ELAP_support"] == 1),
                (analysis["BID_support"] == 1) & (analysis["ELAP_support"] == 1),
            ],
            ["PRAISE_only", "PRAISE+BID_same_family",
             "PRAISE+ELAP_cross_chemistry", "PRAISE+BID+ELAP"],
            default="unexpected",
        )
    ), "support_pattern")
    pattern_rates.to_csv(
        AIM2 / "aim2a_v5_support_pattern_PUS_rates.tsv",
        sep="\t", index=False
    )

    analysis["chemistry_breadth"] = 1 + analysis["ELAP_support"]
    analysis["assay_breadth"] = 1 + analysis["BID_support"] + analysis["ELAP_support"]
    chemistry_rates = rate_table(analysis, "chemistry_breadth")
    chemistry_rates.to_csv(
        AIM2 / "aim2a_v5_chemistry_breadth_PUS_rates.tsv",
        sep="\t", index=False
    )

    # Primary adjusted model.
    primary, cov, names = logistic_cluster(
        analysis, "writer_assigned", PRIMARY_PREDICTORS
    )
    primary["analysis"] = "primary_full_interval_aware"
    primary.to_csv(
        AIM2 / "aim2a_v5_primary_model.tsv",
        sep="\t", index=False
    )

    contrast = coefficient_contrast(
        primary, cov, names, "ELAP_support", "BID_support"
    )
    pd.DataFrame([contrast]).to_csv(
        AIM2 / "aim2a_v5_ELAP_vs_BID_effect_contrast.tsv",
        sep="\t", index=False
    )

    # Sensitivity 1: point source sites only, eliminating PRAISE interval ambiguity.
    point = analysis[pd.to_numeric(analysis["source_interval_flag"], errors="coerce") == 0].copy()
    point_model, _, _ = logistic_cluster(
        point, "writer_assigned",
        ["ELAP_support", "BID_support", "source_signal_rank"]
    )
    point_model["analysis"] = "sensitivity_PRAISE_point_sites_only"

    # Sensitivity 2: external matches must not have both strands inside source interval.
    unambig = analysis[
        (source["BID_strand_ambiguous"] == 0)
        & (source["ELAP_strand_ambiguous"] == 0)
    ].copy()
    strand_model, _, _ = logistic_cluster(
        unambig, "writer_assigned", PRIMARY_PREDICTORS
    )
    strand_model["analysis"] = "sensitivity_no_external_strand_ambiguity"

    # Sensitivity 3: stronger half of PRAISE source signal.
    high = analysis[analysis["source_signal_rank"] >= 0.5].copy()
    high_model, _, _ = logistic_cluster(
        high, "writer_assigned", PRIMARY_PREDICTORS
    )
    high_model["analysis"] = "sensitivity_top_half_PRAISE_signal"

    sensitivities = pd.concat(
        [point_model, strand_model, high_model], ignore_index=True
    )
    sensitivities.to_csv(
        AIM2 / "aim2a_v5_sensitivity_models.tsv",
        sep="\t", index=False
    )

    # Writer-specific secondary biology.
    writer_specific = writer_specific_models(analysis)
    writer_specific.to_csv(
        AIM2 / "aim2a_v5_writer_specific_ELAP_effects.tsv",
        sep="\t", index=False
    )

    source.to_csv(
        AIM2 / "aim2a_v5_PRAISE_source_positive_evidence_table.tsv",
        sep="\t", index=False
    )

    elap_eff = primary[primary["term"] == "ELAP_support"].iloc[0].to_dict()
    bid_eff = primary[primary["term"] == "BID_support"].iloc[0].to_dict()

    # Closure flags are descriptive, not a hidden score.
    elap_positive = bool(
        elap_eff["beta"] > 0 and elap_eff["OR_CI_low"] > 1
    )
    bid_positive = bool(
        bid_eff["beta"] > 0 and bid_eff["OR_CI_low"] > 1
    )

    point_elap = point_model[point_model["term"] == "ELAP_support"].iloc[0]
    strand_elap = strand_model[strand_model["term"] == "ELAP_support"].iloc[0]
    high_elap = high_model[high_model["term"] == "ELAP_support"].iloc[0]
    sensitivity_direction_consistent = bool(
        point_elap["beta"] > 0
        and strand_elap["beta"] > 0
        and high_elap["beta"] > 0
    )

    status = (
        "AIM2A_ORTHOGONAL_SUPPORT_RESULT_READY"
        if sensitivity_direction_consistent
        else "AIM2A_RESULT_READY_FOR_INTERPRETATION"
    )

    contract = {
        "phase": "2A_v5",
        "scientific_story": "Calibration -> Measurement -> Evidence -> Inference",
        "source_positive_universe": (
            "2,201 unique PRAISE-reported HEK293T chr_site loci from "
            "Supplementary Dataset 2"
        ),
        "outcome": (
            "membership in PRAISE Supplementary Dataset 3 direct site-level "
            "PUS1/PUS7/TRUB1/DKC1-dependent lists; 1,017 unique chr_site loci"
        ),
        "outcome_boundary": (
            "writer_assigned=0 means not listed as dependent under these four "
            "published perturbations, not biological negative"
        ),
        "external_support": {
            "BID_HEK293T": "same bisulfite/deletion family as PRAISE",
            "ELAP_HEK293T": "orthogonal enzymatic-labeling / RT-stop chemistry",
        },
        "primary_estimand": (
            "association of ELAP_support with PUS dependency conditional on "
            "BID_support, PRAISE signal rank, and PRAISE interval status"
        ),
        "primary_model": (
            "writer_assigned ~ ELAP_support + BID_support + "
            "source_signal_rank + source_interval_flag; gene-cluster robust SE"
        ),
        "key_secondary_contrast": (
            "beta_ELAP_support - beta_BID_support; asks whether orthogonal "
            "chemistry association exceeds same-family replication"
        ),
        "matching": (
            "PRAISE published point -> exact coordinate; PRAISE published "
            "interval -> external exact call inside that interval; no +/-k. "
            "Source chr_site strings that cannot be safely parsed remain unknown "
            "and are excluded, never recoded as external non-support."
        ),
        "coordinate_attrition_policy": (
            "<=1% unique PRAISE loci may be transparently excluded for unresolved "
            "author coordinate text; >1% stops the analysis."
        ),
        "sensitivities": [
            "PRAISE point loci only",
            "exclude external strand-ambiguous interval matches",
            "top half of PRAISE source-signal rank",
        ],
        "DRS": (
            "not used in HEK293T primary validation because current held-out "
            "DRS dataset has no matched HEK293T cell line; deferred to Aim2B "
            "matched-cell platform-support analysis"
        ),
        "no_truth_score": True,
        "no_raw_sequencing": True,
    }
    (META / "phase2a_v5_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    summary = {
        "status": status,
        "source_QC": source_qc,
        "writer_QC": writer_qc,
        "external_QC": {
            "BID_HEK293T": {
                "n_calls": bid_idx["n_calls"],
                "n_unique_exact": bid_idx["n_unique_exact"],
            },
            "ELAP_HEK293T": {
                "n_calls": elap_idx["n_calls"],
                "n_unique_exact": elap_idx["n_unique_exact"],
            },
        },
        "n_source_unique_published": int(len(source)),
        "n_source_analyzable": int(len(analysis)),
        "n_coordinate_unknown_excluded": int((~source["analysis_included"]).sum()),
        "n_writer_assigned": int(source["writer_assigned"].sum()),
        "BID_support_n_analyzable": int(analysis["BID_support"].sum()),
        "ELAP_support_n_analyzable": int(analysis["ELAP_support"].sum()),
        "support_pattern_counts": source["support_pattern"].value_counts().to_dict(),
        "coordinate_attrition_outcome_check": {
            "writer_assignment_rate_included": float(analysis["writer_assigned"].mean()),
            "writer_assignment_rate_excluded": (
                float(source.loc[~source["analysis_included"], "writer_assigned"].mean())
                if (~source["analysis_included"]).any() else None
            ),
        },
        "support_pattern_rates": pattern_rates.to_dict(orient="records"),
        "primary_ELAP_effect": elap_eff,
        "primary_BID_effect": bid_eff,
        "ELAP_minus_BID_contrast": contrast,
        "sensitivity_direction_consistent_for_ELAP": sensitivity_direction_consistent,
        "primary_ELAP_CI_above_1": elap_positive,
        "primary_BID_CI_above_1": bid_positive,
        "writer_specific_ELAP_effects": writer_specific.to_dict(orient="records"),
        "next_gate": (
            "If orthogonal ELAP support shows a robust positive association with "
            "PUS dependency, freeze this as Aim2 biological validation and proceed "
            "to Aim2B held-out native-DRS support in matched HeLa/A549 evidence "
            "strata. If null, do not tune thresholds; interpret same-family vs "
            "cross-chemistry effects and proceed to the held-out platform test."
        ),
    }
    (META / "phase2a_v5_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print("PRAISE published source unique loci:", len(source))
    print("Analyzable coordinate-safe loci:", len(analysis))
    print("PUS-dependent unique loci:", int(source["writer_assigned"].sum()))
    print("BID same-family support:", int(analysis["BID_support"].sum()))
    print("ELAP cross-chemistry support:", int(analysis["ELAP_support"].sum()))
    print()
    print("SUPPORT PATTERN RATES")
    print(pattern_rates.to_string(index=False))
    print()
    print("PRIMARY MODEL")
    print(primary.to_string(index=False))
    print()
    print("ELAP vs BID coefficient contrast")
    print(pd.DataFrame([contrast]).to_string(index=False))
    print()
    print("SENSITIVITIES")
    print(sensitivities.to_string(index=False))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase2a_v5_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase2a_v5_analysis_contract.json")
    print(r"  D:\RNA\Trace\05_aim2\aim2a_v5_support_pattern_PUS_rates.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2a_v5_primary_model.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2a_v5_ELAP_vs_BID_effect_contrast.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2a_v5_sensitivity_models.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2a_v5_writer_specific_ELAP_effects.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase2a_v5_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase2a_v5_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
