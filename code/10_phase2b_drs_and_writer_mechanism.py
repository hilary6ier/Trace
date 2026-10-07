#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRACE Phase 2B — held-out native DRS + writer-conditioned measurement effects
=============================================================================

Scientific role
---------------
Aim 2A falsified the simple hypothesis that orthogonal ELAP support is
universally associated with stronger PUS-dependency evidence in PRAISE.
This script does NOT retune that result.

It now asks two narrower questions:

A. Predeclared held-out platform question
   Does broader pre-DRS support across BID/BACS/ELAP predict whether a HeLa
   locus is also reported by the independent native direct-RNA DRS platform?

B. Explanatory bridge back to Aim 1
   Is the surprising Aim2A heterogeneity writer-conditioned, and is it
   directionally coherent with ELAP's independent synthetic sequence calibration?

Interpretation boundary
-----------------------
- DRS support means "reported/detected by the published DRS pipeline", not truth.
- Absence from a source positive list is never called a biological negative.
- The PRAISE writer-composition analysis is explicitly explanatory/post-Aim2A,
  not a new confirmatory rescue analysis.
- No thresholds are tuned from downstream results.
- No raw sequencing is used.

Primary DRS design
------------------
Build a HeLa source-positive union from locally single-U exact loci:
    BID + BACS + ELAP
DRS is held out until this union and chemistry breadth are frozen.

Locally single-U means the reported five-mer has a central U and neither
adjacent base (-1,+1) is U. This removes U-run localization ambiguity and
allows exact coordinate comparison to DRS, which does not report strand in S6.

Primary outcome:
    DRS_reported_support = 1 if the exact genomic coordinate AND five-mer
    match a HeLa S6 row with author field psi=1; otherwise 0.

Primary inference:
    DRS report rate by chemistry breadth (1/2/3)
    + motif-stratified permutation test for monotonic association.

Qualified sensitivity:
    among source loci explicitly present in S6 with N_reads_Direct >= 10,
    model author psi (1/0) as a function of chemistry breadth and direct-read
    coverage, with motif-cluster robust SE.

Project root
------------
D:\\RNA\\Trace
"""

from __future__ import annotations

import json
import math
import re
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
NORM = ROOT / "03_harmonized" / "normalized_sources"
AIM1 = ROOT / "04_aim1"
AIM2 = ROOT / "05_aim2"
LOGS = ROOT / "logs"

for p in [META, AIM2, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

BID_HELA = NORM / "BID_HeLa.tsv"
BACS_HELA = NORM / "BACS_HeLa.tsv"
ELAP_HELA = NORM / "ELAP_HeLa.tsv"
DRS_S6 = NORM / "DRS_S6_psi_calls.tsv"

PRAISE_AIM2A = AIM2 / "aim2a_v5_PRAISE_source_positive_evidence_table.tsv"
CALIBRATION = AIM1 / "calibration_resolved.tsv"

SEED = 20261003
RNG = np.random.default_rng(SEED)
N_PERM = 20000

WRITERS = ["DKC1", "TRUB1", "PUS7", "PUS1"]


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

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
        return "chr" + z[3:]
    if re.fullmatch(r"(?:\d+|X|Y|M|MT)", z, flags=re.I):
        return "chr" + ("M" if z.upper() == "MT" else z)
    return z


def motif5(x: Any) -> str:
    z = ss(x).upper().replace("Ψ", "U").replace("T", "U")
    z = re.sub(r"[^ACGU]", "", z)
    if len(z) == 5 and z[2] == "U":
        return z
    return ""


def single_u(m: str) -> bool:
    return bool(len(m) == 5 and m[2] == "U" and m[1] != "U" and m[3] != "U")


def wilson(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if n <= 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z*z/n
    ctr = (p + z*z/(2*n)) / den
    rad = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return ctr-rad, ctr+rad


def bh_adjust(pvals: Sequence[float]) -> np.ndarray:
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
    for i in range(n-1, -1, -1):
        rank = i + 1
        val = min(prev, ranked[i] * n / rank)
        out[order[i]] = val
        prev = val
    q[np.where(ok)[0]] = np.minimum(out, 1.0)
    return q


def semantic_col(df: pd.DataFrame, candidates: Sequence[str], required=True) -> Optional[str]:
    """
    Resolve normalized provenance columns such as raw__N_reads_Direct without
    depending on exact punctuation/case.
    """
    canon = {}
    for c in df.columns:
        z = c.lower()
        if z.startswith("raw__"):
            z = z[5:]
        z = re.sub(r"[^a-z0-9]+", "", z)
        canon[c] = z

    for token in candidates:
        t = re.sub(r"[^a-z0-9]+", "", token.lower())
        exact = [c for c, z in canon.items() if z == t]
        if len(exact) == 1:
            return exact[0]

    for token in candidates:
        t = re.sub(r"[^a-z0-9]+", "", token.lower())
        partial = [c for c, z in canon.items() if t in z]
        if len(partial) == 1:
            return partial[0]

    if required:
        raise RuntimeError(
            f"Could not resolve one of {list(candidates)} from columns: {list(df.columns)}"
        )
    return None


# -----------------------------------------------------------------------------
# Source assay exact single-U calls
# -----------------------------------------------------------------------------

def find_elap_seqcols(df: pd.DataFrame) -> List[str]:
    raw = [c for c in df.columns if c.startswith("raw__")]
    candidates = []
    for i in range(len(raw)-4):
        cols = raw[i:i+5]
        fractions = []
        for c in cols:
            z = df[c].fillna("").astype(str).str.upper().str.strip()
            fractions.append(z.isin(["A","C","G","T","U"]).mean())
        center = df[cols[2]].fillna("").astype(str).str.upper().str.strip()
        center_u = center.isin(["T","U"]).mean()
        if min(fractions) > .90 and center_u > .90:
            candidates.append((float(np.mean(fractions)+center_u), cols))
    if not candidates:
        raise RuntimeError("Could not identify ELAP five-base sequence columns.")
    candidates.sort(reverse=True)
    return candidates[0][1]


def load_source_assay(path: Path, assay: str) -> Tuple[pd.DataFrame, Dict]:
    df = pd.read_csv(path, sep="\t", low_memory=False)

    for c in ["chrom", "pos1", "strand"]:
        if c not in df.columns:
            raise RuntimeError(f"{assay}: missing {c}")

    if assay in {"BID", "BACS"}:
        if "motif_5mer_reported" not in df.columns:
            raise RuntimeError(f"{assay}: motif_5mer_reported missing")
        motifs = df["motif_5mer_reported"].map(motif5)
        motif_source = "motif_5mer_reported"
    elif assay == "ELAP":
        cols = find_elap_seqcols(df)
        seq = df[cols].fillna("").astype(str).agg("".join, axis=1)
        motifs = seq.map(motif5)
        motif_source = "|".join(cols)
    else:
        raise ValueError(assay)

    x = pd.DataFrame({
        "chrom": df["chrom"].map(norm_chr),
        "pos1": pd.to_numeric(df["pos1"], errors="coerce"),
        "strand": df["strand"].astype(str).str.strip(),
        "motif": motifs,
    })
    x = x[
        x["chrom"].ne("")
        & x["pos1"].notna()
        & x["strand"].isin(["+","-"])
        & x["motif"].map(single_u)
    ].copy()
    x["pos1"] = x["pos1"].astype(int)
    x["key"] = x["chrom"] + ":" + x["pos1"].astype(str)
    x["assay"] = assay

    duplicate_exact = int(x.duplicated(["chrom","strand","pos1"]).sum())
    if duplicate_exact:
        raise RuntimeError(
            f"{assay}: {duplicate_exact} duplicated exact single-U calls."
        )

    qc = {
        "assay": assay,
        "source_rows": len(df),
        "single_U_exact_rows": len(x),
        "motif_source": motif_source,
    }
    return x, qc


def build_hela_union() -> Tuple[pd.DataFrame, Dict]:
    parts = []
    qcs = []
    for assay, path in [
        ("BID", BID_HELA),
        ("BACS", BACS_HELA),
        ("ELAP", ELAP_HELA),
    ]:
        x, qc = load_source_assay(path, assay)
        parts.append(x)
        qcs.append(qc)

    allx = pd.concat(parts, ignore_index=True)

    rows = []
    ambiguous = []
    for key, g in allx.groupby("key", sort=False):
        strands = sorted(set(g["strand"]))
        motifs = sorted(set(g["motif"]))
        assays = sorted(set(g["assay"]))

        if len(strands) != 1 or len(motifs) != 1:
            ambiguous.append({
                "key": key,
                "strands": "|".join(strands),
                "motifs": "|".join(motifs),
                "assays": "|".join(assays),
            })
            continue

        chrom, pos = key.rsplit(":", 1)
        row = {
            "key": key,
            "chrom": chrom,
            "pos1": int(pos),
            "strand": strands[0],
            "motif": motifs[0],
            "BID_support": int("BID" in assays),
            "BACS_support": int("BACS" in assays),
            "ELAP_support": int("ELAP" in assays),
        }
        row["chemistry_breadth"] = (
            row["BID_support"] + row["BACS_support"] + row["ELAP_support"]
        )
        present = [a for a in ["BID","BACS","ELAP"] if row[f"{a}_support"]]
        row["support_pattern"] = "+".join(present)
        rows.append(row)

    amb = pd.DataFrame(ambiguous)
    amb.to_csv(
        META / "phase2b_source_coordinate_ambiguities.tsv",
        sep="\t", index=False
    )

    union = pd.DataFrame(rows)

    if len(union) < 500:
        raise RuntimeError(f"HeLa single-U source union unexpectedly small: {len(union)}")

    qc = {
        "source_assays": qcs,
        "union_unambiguous_loci": int(len(union)),
        "ambiguous_coordinate_or_motif_loci_excluded": int(len(amb)),
        "breadth_counts": union["chemistry_breadth"].value_counts().sort_index().to_dict(),
        "pattern_counts": union["support_pattern"].value_counts().to_dict(),
    }
    return union, qc


# -----------------------------------------------------------------------------
# DRS S6
# -----------------------------------------------------------------------------

def load_drs_hela() -> Tuple[pd.DataFrame, Dict]:
    df = pd.read_csv(DRS_S6, sep="\t", low_memory=False, dtype=str)
    if "source_sheet" not in df.columns:
        raise RuntimeError("DRS S6 normalized table lacks source_sheet.")

    h = df[df["source_sheet"].astype(str).str.lower().eq("hela")].copy()
    if len(h) < 1000:
        raise RuntimeError(f"DRS HeLa S6 rows unexpectedly small: {len(h)}")

    chr_col = semantic_col(h, ["chr","chromosome"])
    pos_col = semantic_col(h, ["position","pos"])
    kmer_col = semantic_col(h, ["kmer","5mer"])
    psi_col = semantic_col(h, ["psi"])
    direct_reads_col = semantic_col(h, ["N_reads_Direct","reads_direct"], required=False)
    ivt_reads_col = semantic_col(h, ["N_reads_IVT","reads_ivt"], required=False)

    out = pd.DataFrame({
        "chrom": h[chr_col].map(norm_chr),
        "pos1": pd.to_numeric(h[pos_col], errors="coerce"),
        "motif": h[kmer_col].map(motif5),
        "psi": pd.to_numeric(h[psi_col], errors="coerce"),
    })
    out["N_reads_Direct"] = (
        pd.to_numeric(h[direct_reads_col], errors="coerce")
        if direct_reads_col else np.nan
    )
    out["N_reads_IVT"] = (
        pd.to_numeric(h[ivt_reads_col], errors="coerce")
        if ivt_reads_col else np.nan
    )

    out = out[
        out["chrom"].ne("")
        & out["pos1"].notna()
        & out["motif"].map(single_u)
        & out["psi"].isin([0,1])
    ].copy()
    out["pos1"] = out["pos1"].astype(int)
    out["psi"] = out["psi"].astype(int)
    out["key"] = out["chrom"] + ":" + out["pos1"].astype(str)

    # The match object is (genomic coordinate, five-mer), not coordinate alone.
    # If the same pair repeats, author detection=1 dominates 0 and coverage uses max.
    out = (
        out.groupby(["key","motif"], as_index=False)
        .agg(
            chrom=("chrom","first"),
            pos1=("pos1","first"),
            psi=("psi","max"),
            N_reads_Direct=("N_reads_Direct","max"),
            N_reads_IVT=("N_reads_IVT","max"),
            n_DRS_rows=("psi","size"),
        )
    )

    qc = {
        "HeLa_rows_raw": int(len(h)),
        "single_U_coordinate_motif_pairs": int(len(out)),
        "psi_detected_n": int(out["psi"].sum()),
        "direct_read_column": direct_reads_col,
        "ivt_read_column": ivt_reads_col,
        "schema": {
            "chr": chr_col, "position": pos_col, "kmer": kmer_col, "psi": psi_col
        },
    }
    return out, qc


def attach_drs(union: pd.DataFrame, drs: pd.DataFrame) -> pd.DataFrame:
    d = drs[[
        "key","motif","psi","N_reads_Direct","N_reads_IVT","n_DRS_rows"
    ]].copy()
    d["DRS_table_present"] = 1
    d = d.rename(columns={"psi":"DRS_psi"})

    x = union.merge(d, on=["key","motif"], how="left", validate="one_to_one")
    x["DRS_table_present"] = x["DRS_table_present"].fillna(0).astype(int)
    x["DRS_reported_support"] = x["DRS_psi"].fillna(0).astype(int)
    return x


# -----------------------------------------------------------------------------
# DRS inference
# -----------------------------------------------------------------------------

def rate_table(df: pd.DataFrame, group: str, outcome: str) -> pd.DataFrame:
    rows = []
    for label, z in df.groupby(group, dropna=False):
        n = len(z)
        k = int(z[outcome].sum())
        lo, hi = wilson(k,n)
        rows.append({
            group: label,
            "n": n,
            "DRS_supported_n": k,
            "DRS_support_rate": k/n if n else np.nan,
            "wilson95_low": lo,
            "wilson95_high": hi,
        })
    return pd.DataFrame(rows)


def logistic_cluster(df: pd.DataFrame, outcome: str, predictors: List[str],
                     cluster: str = "motif") -> pd.DataFrame:
    d = df[[outcome, cluster] + predictors].copy()
    for c in [outcome] + predictors:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d[cluster] = df.loc[d.index, cluster].astype(str)
    d = d.dropna()

    y = d[outcome].to_numpy(float)
    X = np.column_stack([np.ones(len(d))] + [
        d[p].to_numpy(float) for p in predictors
    ])
    names = ["intercept"] + predictors

    if len(np.unique(y)) < 2:
        raise RuntimeError(f"{outcome} has one class in model universe.")

    beta = np.zeros(X.shape[1])
    for _ in range(250):
        eta = np.clip(X @ beta, -30, 30)
        pr = 1/(1+np.exp(-eta))
        W = np.clip(pr*(1-pr), 1e-8, None)
        H = X.T @ (W[:,None]*X)
        score = X.T @ (y-pr)
        step = np.linalg.pinv(H) @ score
        beta2 = beta + step
        if np.max(np.abs(step)) < 1e-9:
            beta = beta2
            break
        beta = beta2

    eta = np.clip(X @ beta, -30, 30)
    pr = 1/(1+np.exp(-eta))
    W = np.clip(pr*(1-pr), 1e-8, None)
    bread = np.linalg.pinv(X.T @ (W[:,None]*X))

    groups = defaultdict(list)
    for i,g in enumerate(d[cluster].astype(str)):
        groups[g].append(i)
    meat = np.zeros((X.shape[1], X.shape[1]))
    for ids in groups.values():
        sg = X[ids].T @ (y[ids]-pr[ids])
        meat += np.outer(sg,sg)

    G = len(groups)
    N,K = X.shape
    correction = (G/max(G-1,1))*((N-1)/max(N-K,1))
    cov = bread @ meat @ bread * correction
    se = np.sqrt(np.clip(np.diag(cov),0,None))

    rows = []
    for name,b,s0 in zip(names,beta,se):
        z = b/s0 if s0 > 0 else np.nan
        p = math.erfc(abs(z)/math.sqrt(2)) if np.isfinite(z) else np.nan
        rows.append({
            "term": name,
            "beta": float(b),
            "SE_motif_cluster": float(s0),
            "z": float(z) if np.isfinite(z) else np.nan,
            "p": float(p) if np.isfinite(p) else np.nan,
            "OR": math.exp(float(b)),
            "OR_CI_low": math.exp(float(b-1.96*s0)),
            "OR_CI_high": math.exp(float(b+1.96*s0)),
            "n": N,
            "n_motif_clusters": G,
        })
    return pd.DataFrame(rows)


def motif_stratified_permutation(df: pd.DataFrame) -> Dict:
    """
    Statistic = difference in mean chemistry breadth between DRS-supported
    and DRS-not-reported source loci. Labels are permuted only within the
    same five-mer motif, preserving local sequence composition.
    """
    d = df[["motif","chemistry_breadth","DRS_reported_support"]].copy()
    y = d["DRS_reported_support"].to_numpy(int)
    b = d["chemistry_breadth"].to_numpy(float)

    def stat(y0):
        if (y0==1).sum() == 0 or (y0==0).sum() == 0:
            return np.nan
        return float(b[y0==1].mean() - b[y0==0].mean())

    obs = stat(y)
    motif_groups = {
        m: idx.to_numpy()
        for m, idx in d.groupby("motif").groups.items()
    }

    ge = 0
    le = 0
    for _ in range(N_PERM):
        yp = y.copy()
        for ids in motif_groups.values():
            yp[ids] = RNG.permutation(yp[ids])
        st = stat(yp)
        if st >= obs:
            ge += 1
        if st <= obs:
            le += 1

    p_upper = (ge+1)/(N_PERM+1)
    p_lower = (le+1)/(N_PERM+1)
    p_two = min(1.0, 2*min(p_upper,p_lower))
    return {
        "statistic": "mean_breadth_DRS1_minus_DRS0",
        "observed": obs,
        "permutations": N_PERM,
        "p_upper": p_upper,
        "p_lower": p_lower,
        "p_two_sided": p_two,
    }


# -----------------------------------------------------------------------------
# Explanatory Aim2A writer-conditioned follow-up
# -----------------------------------------------------------------------------

def writer_composition_followup() -> pd.DataFrame:
    df = pd.read_csv(PRAISE_AIM2A, sep="\t", low_memory=False)
    req = [
        "writer_assigned","writer_set","ELAP_support","BID_support",
        "source_signal_rank","source_interval_flag","gene"
    ]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise RuntimeError(f"Aim2A evidence table missing {missing}")

    z = df[df["writer_assigned"] == 1].copy()
    z["writer_n"] = z["writer_set"].fillna("").astype(str).map(
        lambda q: len([v for v in q.split(",") if v])
    )
    z = z[z["writer_n"] == 1].copy()
    z["writer"] = z["writer_set"].astype(str)

    rows = []
    # Full adjusted and point-only adjusted are both reported; no best one selected.
    for scope_name, q in [
        ("all_coordinate_resolved", z),
        ("PRAISE_point_sites_only", z[z["source_interval_flag"] == 0].copy()),
    ]:
        for w in WRITERS:
            if (q["writer"] == w).sum() < 10:
                continue
            q2 = q.copy()
            q2["writer_identity_outcome"] = (q2["writer"] == w).astype(int)
            predictors = ["ELAP_support","BID_support","source_signal_rank"]
            if scope_name == "all_coordinate_resolved":
                predictors.append("source_interval_flag")
            tab = logistic_cluster(
                q2, "writer_identity_outcome", predictors, cluster="gene"
            )
            eff = tab[tab["term"] == "ELAP_support"].iloc[0].to_dict()
            eff.update({
                "writer": w,
                "scope": scope_name,
                "writer_cases": int((q2["writer"] == w).sum()),
                "writer_assigned_universe": int(len(q2)),
            })
            rows.append(eff)

    out = pd.DataFrame(rows)
    out["q_BH_within_scope"] = np.nan
    for scope, idx in out.groupby("scope").groups.items():
        out.loc[idx, "q_BH_within_scope"] = bh_adjust(out.loc[idx, "p"])
    return out


def calibration_writer_context() -> pd.DataFrame:
    cal = pd.read_csv(CALIBRATION, sep="\t")
    req = {"motif","response_percentile","assay"}
    if not req.issubset(cal.columns):
        raise RuntimeError("calibration_resolved.tsv schema unexpected")

    wide = cal[cal["assay"].isin(["BID","ELAP"])].pivot(
        index="motif", columns="assay", values="response_percentile"
    ).dropna()
    if len(wide) != 256:
        raise RuntimeError(f"Expected 256 BID/ELAP contexts, observed {len(wide)}")
    wide["ELAP_minus_BID"] = wide["ELAP"] - wide["BID"]

    sets = {
        "PUS7_UNUAR": lambda m: bool(re.match(r"^U[ACGU]UA[AG]$", m)),
        "TRUB1_GUUCN": lambda m: bool(re.match(r"^GUUC[ACGU]$", m)),
    }

    rows = []
    all_elap = wide["ELAP"].to_numpy(float)
    all_diff = wide["ELAP_minus_BID"].to_numpy(float)

    for name, fn in sets.items():
        sub = wide[[fn(m) for m in wide.index]]
        k = len(sub)
        if k == 0:
            continue

        obs_elap = float(sub["ELAP"].mean())
        obs_diff = float(sub["ELAP_minus_BID"].mean())

        # Two-sided Monte Carlo random-set tests because this is explanatory
        # follow-up after observing Aim2A heterogeneity.
        nperm = 100000
        elap_vals = np.empty(nperm)
        diff_vals = np.empty(nperm)
        for i in range(nperm):
            ids = RNG.choice(len(wide), size=k, replace=False)
            elap_vals[i] = all_elap[ids].mean()
            diff_vals[i] = all_diff[ids].mean()

        def two_sided(vals, obs):
            pl = (np.sum(vals <= obs)+1)/(len(vals)+1)
            pu = (np.sum(vals >= obs)+1)/(len(vals)+1)
            return min(1.0, 2*min(pl,pu))

        rows.append({
            "motif_set": name,
            "n_contexts": k,
            "mean_ELAP_percentile": obs_elap,
            "median_ELAP_percentile": float(sub["ELAP"].median()),
            "mean_BID_percentile": float(sub["BID"].mean()),
            "mean_ELAP_minus_BID_percentile": obs_diff,
            "two_sided_random_set_p_ELAP": two_sided(elap_vals, obs_elap),
            "two_sided_random_set_p_ELAP_minus_BID": two_sided(diff_vals, obs_diff),
            "motifs": "|".join(sub.index.tolist()),
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    print("="*104)
    print("TRACE PHASE 2B — HELD-OUT DRS + WRITER-CONDITIONED MEASUREMENT EFFECTS")
    print("="*104)

    for p in [
        BID_HELA, BACS_HELA, ELAP_HELA, DRS_S6,
        PRAISE_AIM2A, CALIBRATION
    ]:
        if not p.exists():
            raise RuntimeError(f"Missing required input: {p}")

    union, source_qc = build_hela_union()
    drs, drs_qc = load_drs_hela()
    evidence = attach_drs(union, drs)

    evidence.to_csv(
        AIM2 / "aim2b_HeLa_source_union_DRS.tsv",
        sep="\t", index=False
    )

    breadth_rates = rate_table(
        evidence, "chemistry_breadth", "DRS_reported_support"
    )
    breadth_rates.to_csv(
        AIM2 / "aim2b_DRS_breadth_rates.tsv",
        sep="\t", index=False
    )

    pattern_rates = rate_table(
        evidence, "support_pattern", "DRS_reported_support"
    )
    pattern_rates.to_csv(
        AIM2 / "aim2b_DRS_pattern_rates.tsv",
        sep="\t", index=False
    )

    trend_model = logistic_cluster(
        evidence,
        "DRS_reported_support",
        ["chemistry_breadth"],
        cluster="motif"
    )
    trend_model["analysis"] = "all_source_positive_union_reported_support"
    trend_model.to_csv(
        AIM2 / "aim2b_DRS_breadth_model.tsv",
        sep="\t", index=False
    )

    pattern_model = logistic_cluster(
        evidence,
        "DRS_reported_support",
        ["BID_support","BACS_support","ELAP_support"],
        cluster="motif"
    )
    pattern_model["analysis"] = "assay_specific_source_support"
    pattern_model.to_csv(
        AIM2 / "aim2b_DRS_assay_specific_model.tsv",
        sep="\t", index=False
    )

    perm = motif_stratified_permutation(evidence)
    pd.DataFrame([perm]).to_csv(
        AIM2 / "aim2b_DRS_motif_stratified_permutation.tsv",
        sep="\t", index=False
    )

    # Qualified DRS-evaluable sensitivity.
    evaluable = evidence[
        (evidence["DRS_table_present"] == 1)
        & pd.to_numeric(evidence["N_reads_Direct"], errors="coerce").ge(10)
        & evidence["DRS_psi"].isin([0,1])
    ].copy()
    evaluable["log1p_DRS_direct_reads"] = np.log1p(
        pd.to_numeric(evaluable["N_reads_Direct"], errors="coerce")
    )

    if len(evaluable) >= 100 and evaluable["DRS_psi"].nunique() == 2:
        eval_model = logistic_cluster(
            evaluable,
            "DRS_psi",
            ["chemistry_breadth","log1p_DRS_direct_reads"],
            cluster="motif"
        )
        eval_model["analysis"] = "DRS_table_evaluable_direct_reads_ge_10"
    else:
        eval_model = pd.DataFrame([{
            "term": "NOT_RUN",
            "analysis": "DRS_table_evaluable_direct_reads_ge_10",
            "n": len(evaluable),
            "reason": "insufficient evaluable rows or outcome classes",
        }])
    eval_model.to_csv(
        AIM2 / "aim2b_DRS_evaluable_sensitivity.tsv",
        sep="\t", index=False
    )

    # Explanatory bridge to Aim1/Aim2A.
    writer_follow = writer_composition_followup()
    writer_follow.to_csv(
        AIM2 / "aim2b_writer_conditioned_ELAP_effects.tsv",
        sep="\t", index=False
    )

    cal_writer = calibration_writer_context()
    cal_writer.to_csv(
        AIM2 / "aim2b_writer_motif_calibration_context.tsv",
        sep="\t", index=False
    )

    breadth_eff = trend_model[
        trend_model["term"] == "chemistry_breadth"
    ].iloc[0].to_dict()

    rates_sorted = breadth_rates.sort_values("chemistry_breadth")
    monotonic = bool(
        rates_sorted["DRS_support_rate"].is_monotonic_increasing
        and len(rates_sorted) >= 2
    )
    positive_model = bool(
        breadth_eff["beta"] > 0 and breadth_eff["OR_CI_low"] > 1
    )
    motif_perm_positive = bool(
        perm["observed"] > 0 and perm["p_upper"] < 0.05
    )

    if monotonic and positive_model and motif_perm_positive:
        status = "AIM2B_DRS_BREADTH_SUPPORTS_CROSS_PLATFORM_REPRODUCIBILITY"
    else:
        status = "AIM2B_DRS_RESULT_READY_FOR_INTERPRETATION"

    contract = {
        "phase": "2B",
        "scientific_story": "Calibration -> Measurement -> Evidence -> Inference",
        "Aim2A_frozen_result": (
            "No robust positive generic association between orthogonal ELAP "
            "support and PRAISE PUS-dependency. No threshold/model retuning."
        ),
        "primary_Aim2B_question": (
            "Among HeLa loci already reported by >=1 of BID/BACS/ELAP, does "
            "pre-DRS chemistry breadth predict exact single-U support by held-out DRS?"
        ),
        "source_assays": ["BID_HeLa","BACS_HeLa","ELAP_HeLa"],
        "heldout_platform": "DRS HeLa Supplementary Table S6",
        "primary_locus_rule": (
            "locally single-U exact genomic coordinate; source strand/motif conflicts "
            "excluded; DRS match requires exact coordinate + same five-mer"
        ),
        "primary_outcome": (
            "DRS_reported_support; absence means not reported/detected in the "
            "published DRS table, not biological negative"
        ),
        "primary_statistics": [
            "DRS support rate by chemistry breadth 1/2/3",
            "motif-cluster robust logistic breadth effect",
            "motif-stratified permutation of DRS labels"
        ],
        "qualified_sensitivity": (
            "source loci explicitly present in S6 with N_reads_Direct>=10; "
            "author psi 1/0 modeled with chemistry breadth + direct-read coverage"
        ),
        "explanatory_writer_followup": (
            "post-Aim2A, explicitly exploratory: writer identity among already "
            "PUS-dependent PRAISE sites vs ELAP/BID support, plus independent "
            "synthetic calibration context for PUS7/TRUB1 motifs"
        ),
        "no_truth_score": True,
        "no_raw_sequencing": True,
    }
    (META / "phase2b_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    summary = {
        "status": status,
        "source_union_QC": source_qc,
        "DRS_QC": drs_qc,
        "n_source_union": int(len(evidence)),
        "DRS_reported_support_n": int(evidence["DRS_reported_support"].sum()),
        "DRS_table_present_n": int(evidence["DRS_table_present"].sum()),
        "DRS_evaluable_direct_reads_ge10_n": int(len(evaluable)),
        "breadth_rates": breadth_rates.to_dict(orient="records"),
        "pattern_rates": pattern_rates.to_dict(orient="records"),
        "breadth_effect": breadth_eff,
        "motif_stratified_permutation": perm,
        "monotonic_rate_pattern": monotonic,
        "writer_conditioned_followup": writer_follow.to_dict(orient="records"),
        "writer_motif_calibration_context": cal_writer.to_dict(orient="records"),
        "next_gate": (
            "If DRS breadth is positively supported, Aim2 concludes that chemistry "
            "breadth predicts cross-platform reproducibility but not generic PUS "
            "dependency. If DRS breadth is also null/non-monotonic, freeze Aim2 as "
            "evidence that simple consensus breadth is not a universal confidence "
            "measure and proceed directly to Aim3 map-choice functional inference."
        ),
    }
    (META / "phase2b_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print("HeLa source union:", len(evidence))
    print("Held-out DRS reported support:", int(evidence["DRS_reported_support"].sum()))
    print("DRS S6 exact+motif present:", int(evidence["DRS_table_present"].sum()))
    print()
    print("DRS BREADTH RATES")
    print(breadth_rates.to_string(index=False))
    print()
    print("DRS BREADTH MODEL")
    print(trend_model.to_string(index=False))
    print()
    print("MOTIF-STRATIFIED PERMUTATION")
    print(pd.DataFrame([perm]).to_string(index=False))
    print()
    print("WRITER-CONDITIONED FOLLOW-UP")
    print(writer_follow.to_string(index=False))
    print()
    print("SYNTHETIC CALIBRATION / WRITER MOTIF CONTEXT")
    print(cal_writer.to_string(index=False))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase2b_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase2b_analysis_contract.json")
    print(r"  D:\RNA\Trace\05_aim2\aim2b_DRS_breadth_rates.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2b_DRS_pattern_rates.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2b_DRS_breadth_model.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2b_DRS_evaluable_sensitivity.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2b_writer_conditioned_ELAP_effects.tsv")
    print(r"  D:\RNA\Trace\05_aim2\aim2b_writer_motif_calibration_context.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase2b_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase2b_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
