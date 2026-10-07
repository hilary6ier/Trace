#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE Phase 1A — calibration resolution + first scientific test
===============================================================

Purpose
-------
1. Resolve synthetic calibration using published assay semantics.
2. Amend the PRAISE interval parser discrepancy without allowing it to block Aim 1.
3. Reconstruct 5-mer sequence context from normalized BID / ELAP source tables.
4. Run the first preregistered scientific experiment:

    synthetic assay-context sensitivity
        -> BID-vs-ELAP human-map discordance

Primary human analysis is restricted to exact, locally single-U loci.
Shared sites are NOT used as training labels and are evaluated separately
as a continuum check.

No FASTQ/BAM/raw sequencing.
No "true Psi" inference.
No missing-call-as-biological-negative interpretation.
No model/feature hunting.

Project root:
    D:\RNA\Trace

Requirements:
    pandas
    numpy
    openpyxl

Outputs:
    00_meta/phase1a_*.{json,tsv,txt}
    04_aim1/calibration_resolved.tsv
    04_aim1/aim1a_site_table.tsv
    04_aim1/aim1a_cell_results.tsv
    04_aim1/aim1a_pooled_logit.tsv
    04_aim1/aim1a_continuum.tsv
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
import openpyxl

ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
SRC = ROOT / "01_source"
HARM = ROOT / "03_harmonized"
NORM = HARM / "normalized_sources"
AIM1 = ROOT / "04_aim1"
LOGS = ROOT / "logs"

for p in [META, AIM1, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

SEED = 20261001
RNG = np.random.default_rng(SEED)

PATHS = {
    "cal0c": HARM / "synthetic_context_calibration.tsv",
    "bid_supp": SRC / "BID" / "BID_supplementary_tables.xlsx",
    "bacs_fig1": SRC / "BACS" / "BACS_source_data_fig1.xlsx",
    "elap_source": SRC / "ELAP" / "ELAP_source_data.xlsx",
    "bid_hek": NORM / "BID_HEK293T.tsv",
    "bid_hela": NORM / "BID_HeLa.tsv",
    "elap_hek": NORM / "ELAP_HEK293T.tsv",
    "elap_hela": NORM / "ELAP_HeLa.tsv",
    "praise": NORM / "PRAISE_HEK293T.tsv",
}

PUBLISHED_CALIBRATION = {
    "BID": {
        "primary_response": "fit_R",
        "reason": "R is the BID-seq induced deletion ratio; A is dropout and B is background deletion.",
    },
    "BACS": {
        "primary_response": "mean Psi conversion rate across synthetic NNΨNN contexts",
        "reason": "Fig.1c reports motif-dependent Ψ conversion and U false-positive rates; Ψ conversion is the intrinsic assay response.",
        "published_over_85": 230,
        "published_mean": 0.876,
    },
    "ELAP": {
        "primary_response": "replicate-combined enrichment of synthetic NNΨNN oligos",
        "reason": "Fig.2c-d quantify sequence-context-dependent enrichment across 256 NNΨNN motifs.",
        "positive_control": "+1 U (5'-ΨU-3') preferential enrichment",
    },
}

# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------

def s(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    return str(x).strip()

def to_num(x: Any) -> Optional[float]:
    if x is None:
        return None
    if isinstance(x, (int, float, np.integer, np.floating)):
        x = float(x)
        return None if not np.isfinite(x) else x
    z = s(x).replace(",", "").replace("%", "")
    try:
        v = float(z)
        return v if np.isfinite(v) else None
    except Exception:
        return None

def motif5_token(x: Any, allow_4mer: bool = True) -> str:
    """
    Convert source context labels into RNA 5-mer with central U.
    Accepts:
      AAΨGC -> AAUGC
      AATGC -> AAUGC if center is T
      AAGC  -> AAUGC when central Ψ is implicit (4 flanking bases)
    """
    z = s(x).upper()
    z = z.replace("PSEUDOURIDINE", "Ψ")
    z = z.replace("PSI", "Ψ")
    z = re.sub(r"[\s'\"5′3′5'3'_\-–—/\\]+", "", z)
    z = z.replace("Ψ", "U").replace("T", "U")
    z = re.sub(r"[^ACGU]", "", z)
    if len(z) == 5 and z[2] == "U":
        return z
    if allow_4mer and len(z) == 4:
        return z[:2] + "U" + z[2:]
    return ""

def dinuc_token(x: Any) -> str:
    z = s(x).upper().replace("T", "U")
    z = re.sub(r"[^ACGU]", "", z)
    return z if len(z) == 2 else ""

def pct_rank(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    return x.rank(method="average", pct=True)

def auc_score(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    ok = np.isfinite(score)
    y, score = y[ok], score[ok]
    n1 = int((y == 1).sum())
    n0 = int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = pd.Series(score).rank(method="average").to_numpy()
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

def logistic_cluster(X: np.ndarray, y: np.ndarray, cluster: Sequence[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Logistic IRLS + cluster-robust sandwich covariance.
    Clusters are motif contexts, avoiding pseudo-replication from repeated 5-mers.
    """
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    beta = np.zeros(X.shape[1])
    for _ in range(100):
        eta = np.clip(X @ beta, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-8, None)
        H = X.T @ (W[:, None] * X)
        score = X.T @ (y - p)
        step = np.linalg.pinv(H) @ score
        beta_new = beta + step
        if np.max(np.abs(step)) < 1e-9:
            beta = beta_new
            break
        beta = beta_new

    eta = np.clip(X @ beta, -30, 30)
    p = 1 / (1 + np.exp(-eta))
    W = np.clip(p * (1 - p), 1e-8, None)
    bread = np.linalg.pinv(X.T @ (W[:, None] * X))

    groups = defaultdict(list)
    for i, g in enumerate(cluster):
        groups[str(g)].append(i)

    meat = np.zeros((X.shape[1], X.shape[1]))
    for ids in groups.values():
        Xi = X[ids]
        ui = (y[ids] - p[ids])
        sg = Xi.T @ ui
        meat += np.outer(sg, sg)

    G = len(groups)
    N, K = X.shape
    correction = (G / max(G - 1, 1)) * ((N - 1) / max(N - K, 1))
    cov = bread @ meat @ bread * correction
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    z = beta / np.where(se > 0, se, np.nan)
    return beta, se, z

def p_from_z(z: float) -> float:
    if not np.isfinite(z):
        return np.nan
    return math.erfc(abs(z) / math.sqrt(2))

def cluster_boot_auc(df: pd.DataFrame, n_boot: int = 2000) -> Tuple[float, float]:
    motifs = df["motif"].dropna().unique()
    vals = []
    if len(motifs) < 5:
        return np.nan, np.nan
    grouped = {m: df[df["motif"] == m] for m in motifs}
    for _ in range(n_boot):
        sampled = RNG.choice(motifs, size=len(motifs), replace=True)
        boot = pd.concat([grouped[m] for m in sampled], ignore_index=True)
        a = auc_score(boot["y"].to_numpy(), boot["contrast"].to_numpy())
        if np.isfinite(a):
            vals.append(a)
    if not vals:
        return np.nan, np.nan
    return float(np.quantile(vals, .025)), float(np.quantile(vals, .975))

def motif_permutation_p(df: pd.DataFrame, observed_auc: float, n_perm: int = 5000) -> Tuple[float, float]:
    motif_contrast = df[["motif", "contrast"]].drop_duplicates("motif").set_index("motif")["contrast"]
    motifs = motif_contrast.index.to_numpy()
    values = motif_contrast.to_numpy()
    y = df["y"].to_numpy()
    row_motifs = df["motif"].to_numpy()

    more_one = 0
    more_two = 0
    obs_dev = abs(observed_auc - .5)
    for _ in range(n_perm):
        perm = RNG.permutation(values)
        mp = dict(zip(motifs, perm))
        score = np.array([mp[m] for m in row_motifs], float)
        a = auc_score(y, score)
        if a >= observed_auc:
            more_one += 1
        if abs(a - .5) >= obs_dev:
            more_two += 1
    return (more_one + 1) / (n_perm + 1), (more_two + 1) / (n_perm + 1)

# ---------------------------------------------------------------------
# workbook calibration discovery
# ---------------------------------------------------------------------

def workbook_matrices(path: Path) -> Dict[str, List[List[Any]]]:
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    out = {}
    for ws in wb.worksheets:
        out[ws.title] = [list(r) for r in ws.iter_rows(values_only=True)]
    wb.close()
    return out

def nearest_header(matrix: List[List[Any]], col: int, first_row: int) -> str:
    for r in range(first_row - 1, max(-1, first_row - 10), -1):
        if r < 0:
            break
        if col < len(matrix[r]) and s(matrix[r][col]):
            return s(matrix[r][col])
    return f"col_{col+1}"

def rowwise_candidates(path: Path, assay: str) -> List[Dict[str, Any]]:
    candidates = []
    for sheet, mat in workbook_matrices(path).items():
        if not mat:
            continue
        maxc = max(len(r) for r in mat)
        for c in range(maxc):
            rows = []
            motifs = []
            for ri, row in enumerate(mat):
                v = row[c] if c < len(row) else None
                m = motif5_token(v, allow_4mer=True)
                if m:
                    rows.append(ri)
                    motifs.append(m)
            if len(set(motifs)) < 200:
                continue

            first = min(rows)
            numeric_profiles = []
            for j in range(maxc):
                if j == c:
                    continue
                vals = []
                mids = []
                for ri, m in zip(rows, motifs):
                    row = mat[ri]
                    v = to_num(row[j] if j < len(row) else None)
                    if v is not None:
                        vals.append(v)
                        mids.append(m)
                if len(set(mids)) >= 200:
                    hdr = nearest_header(mat, j, first)
                    tmp = pd.DataFrame({"motif": mids, "v": vals}).groupby("motif", as_index=False)["v"].mean()
                    numeric_profiles.append((j, hdr, tmp))

            candidates.append({
                "layout": "rowwise",
                "sheet": sheet,
                "context_col": c,
                "unique_motifs": len(set(motifs)),
                "numeric_profiles": numeric_profiles,
            })
    return candidates

def matrix16_candidates(path: Path) -> List[Dict[str, Any]]:
    """
    Finds 16x16 heatmap-style blocks where row and column labels are dinucleotides.
    Stores both possible upstream/downstream orientations; ELAP orientation is later
    validated against the published +1U preference.
    """
    out = []
    for sheet, mat in workbook_matrices(path).items():
        if not mat:
            continue
        nr = len(mat)
        nc = max(len(r) for r in mat)

        # Candidate header rows containing >=12 dinucleotide labels
        for hr in range(nr):
            col_labels = []
            for c in range(nc):
                v = mat[hr][c] if c < len(mat[hr]) else None
                d = dinuc_token(v)
                if d:
                    col_labels.append((c, d))
            # de-duplicate by label
            seen = {}
            for c, d in col_labels:
                seen.setdefault(d, c)
            col_labels = [(c, d) for d, c in seen.items()]
            if len(col_labels) < 12:
                continue

            for lc in range(nc):
                row_labels = []
                for r in range(hr + 1, min(nr, hr + 40)):
                    v = mat[r][lc] if lc < len(mat[r]) else None
                    d = dinuc_token(v)
                    if d:
                        row_labels.append((r, d))
                seenr = {}
                for r, d in row_labels:
                    seenr.setdefault(d, r)
                row_labels = [(r, d) for d, r in seenr.items()]
                if len(row_labels) < 12:
                    continue

                records = []
                for r, rd in row_labels:
                    for c, cd in col_labels:
                        if r >= len(mat) or c >= len(mat[r]):
                            continue
                        v = to_num(mat[r][c])
                        if v is not None:
                            records.append((rd, cd, v))
                if len(records) >= 200:
                    out.append({
                        "layout": "matrix16",
                        "sheet": sheet,
                        "header_row": hr,
                        "label_col": lc,
                        "records": records,
                    })
    # keep unique signatures
    uniq = []
    signatures = set()
    for x in out:
        sig = (x["sheet"], x["header_row"], x["label_col"], len(x["records"]))
        if sig not in signatures:
            signatures.add(sig)
            uniq.append(x)
    return uniq

def normalize_fraction(values: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    med = x.median()
    if np.isfinite(med) and med > 1.5:
        return x / 100.0
    return x

def header_score_bacs(h: str) -> int:
    x = h.lower()
    score = 0
    if "conversion" in x: score += 20
    if "bacs" in x or "psi" in x or "ψ" in h: score += 10
    if "mean" in x or "avg" in x: score += 8
    if "rep" in x: score += 4
    if any(k in x for k in ["false", "control", "nnunn", "unmodified", "fp", "sd", "stdev", "error"]):
        score -= 30
    return score

def header_score_elap(h: str) -> int:
    x = h.lower()
    score = 0
    if "enrich" in x: score += 25
    if "pull" in x or "ip" in x: score += 8
    if "rep" in x: score += 6
    if "mean" in x or "avg" in x: score += 6
    if any(k in x for k in ["control", "input", "unmodified", "nnunn", "sd", "stdev", "error"]):
        score -= 25
    return score

def plus1u_effect(df: pd.DataFrame, value_col: str) -> float:
    x = df.copy()
    x[value_col] = pd.to_numeric(x[value_col], errors="coerce")
    a = x.loc[x["motif"].str[3] == "U", value_col].dropna()
    b = x.loc[x["motif"].str[3] != "U", value_col].dropna()
    if len(a) < 5 or len(b) < 20:
        return np.nan
    # robust effect on log1p scale when nonnegative
    amin = min(a.min(), b.min())
    shift = -amin + 1e-9 if amin <= 0 else 0
    return float(np.median(np.log1p(a + shift)) - np.median(np.log1p(b + shift)))

def resolve_bacs(path: Path) -> Tuple[Optional[pd.DataFrame], Dict[str, Any]]:
    cands = rowwise_candidates(path, "BACS")
    scored = []
    for c in cands:
        for j, h, df in c["numeric_profiles"]:
            f = normalize_fraction(df["v"])
            n85 = int((f > .85).sum())
            mean = float(f.mean())
            # Published Fig1c anchors: 230/256 >85%, accumulated mean ~87.6%.
            validation = -abs(n85 - 230) / 5 - abs(mean - .876) * 20
            score = header_score_bacs(h) + validation
            scored.append((score, c, h, df.assign(resp=f), n85, mean))
    if scored:
        scored.sort(key=lambda x: x[0], reverse=True)
        score, c, h, df, n85, mean = scored[0]
        if len(df) >= 250 and abs(n85 - 230) <= 20 and 0.75 <= mean <= 0.98:
            out = df.rename(columns={"v": "raw_value"})
            out["response"] = df["resp"]
            return out[["motif", "raw_value", "response"]], {
                "assay": "BACS", "status": "RESOLVED",
                "layout": c["layout"], "sheet": c["sheet"], "response_source": h,
                "n_motifs": int(out["motif"].nunique()),
                "published_check_n_gt_0.85": n85,
                "published_check_mean": mean,
                "notes": "Selected using source semantics + published Fig1c aggregate anchors, never human outcomes.",
            }

    # Matrix fallback: distribution can validate BACS response, but axis orientation
    # cannot be inferred from distribution alone; do not silently assign motifs.
    mats = matrix16_candidates(path)
    if mats:
        return None, {
            "assay": "BACS", "status": "MATRIX_FOUND_ORIENTATION_REVIEW",
            "layout": "matrix16", "sheet": mats[0]["sheet"],
            "response_source": "", "n_motifs": 256,
            "notes": "Numeric 16x16 block found, but upstream/downstream axis orientation is not inferred from distribution alone. BACS is secondary for Aim1A, so this does not block the primary BID-ELAP test.",
        }
    return None, {
        "assay": "BACS", "status": "UNRESOLVED",
        "layout": "", "sheet": "", "response_source": "", "n_motifs": 0,
        "notes": "No robust 256-context representation recovered.",
    }

def resolve_elap(path: Path) -> Tuple[Optional[pd.DataFrame], Dict[str, Any]]:
    rowc = rowwise_candidates(path, "ELAP")
    scored = []
    for c in rowc:
        for j, h, df in c["numeric_profiles"]:
            effect = plus1u_effect(df.rename(columns={"v": "resp"}), "resp")
            score = header_score_elap(h) + (0 if not np.isfinite(effect) else 10 * effect)
            scored.append((score, c, h, df, effect))
    if scored:
        scored.sort(key=lambda x: x[0], reverse=True)
        score, c, h, df, effect = scored[0]
        if len(df) >= 250 and np.isfinite(effect):
            out = df.rename(columns={"v": "response"})
            return out[["motif", "response"]], {
                "assay": "ELAP", "status": "RESOLVED",
                "layout": c["layout"], "sheet": c["sheet"], "response_source": h,
                "n_motifs": int(out["motif"].nunique()),
                "plus1U_log1p_median_effect": effect,
                "notes": "Mapping validated against the published 5'-PsiU-3' enrichment preference; human outcomes were not inspected.",
            }

    mats = matrix16_candidates(path)
    best = None
    for m in mats:
        rec = pd.DataFrame(m["records"], columns=["row2", "col2", "response"])
        for orientation in ["row_upstream", "col_upstream"]:
            if orientation == "row_upstream":
                rec["motif"] = rec["row2"] + "U" + rec["col2"]
            else:
                rec["motif"] = rec["col2"] + "U" + rec["row2"]
            df = rec.groupby("motif", as_index=False)["response"].mean()
            if df["motif"].nunique() < 240:
                continue
            effect = plus1u_effect(df, "response")
            if not np.isfinite(effect):
                continue
            candidate = (effect, m, orientation, df.copy())
            if best is None or effect > best[0]:
                best = candidate

    if best is not None:
        effect, m, orientation, df = best
        # Positive published chemistry control is used only to orient the source-data axes.
        if effect > 0:
            return df[["motif", "response"]], {
                "assay": "ELAP", "status": "RESOLVED",
                "layout": "matrix16", "sheet": m["sheet"],
                "response_source": f"16x16 heatmap; orientation={orientation}",
                "n_motifs": int(df["motif"].nunique()),
                "plus1U_log1p_median_effect": effect,
                "notes": "Axis orientation chosen by the paper's independently reported +1U preference; no human outcome used.",
            }

    return None, {
        "assay": "ELAP", "status": "UNRESOLVED",
        "layout": "", "sheet": "", "response_source": "", "n_motifs": 0,
        "notes": "Could not recover a validated 256-context enrichment table.",
    }

def resolve_bid(cal0c: Path) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    df = pd.read_csv(cal0c, sep="\t")
    x = df[df["assay"] == "BID"].copy()
    if x["motif_RNA5"].nunique() != 256 or "raw__fit_R" not in x.columns:
        raise RuntimeError("BID calibration does not contain 256 motifs with raw__fit_R.")
    x["response"] = pd.to_numeric(x["raw__fit_R"], errors="coerce")
    if x["response"].notna().sum() != 256:
        raise RuntimeError("BID fit_R is incomplete.")
    return x[["motif_RNA5", "response"]].rename(columns={"motif_RNA5": "motif"}), {
        "assay": "BID", "status": "RESOLVED",
        "layout": "rowwise", "sheet": "Supplementary Table 1",
        "response_source": "fit_R", "n_motifs": 256,
        "notes": "Frozen from published BID model: R = BID-seq induced deletion ratio; A=dropout, B=background.",
    }

# ---------------------------------------------------------------------
# PRAISE parser amendment
# ---------------------------------------------------------------------

def audit_praise_interval(path: Path) -> Dict[str, Any]:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    col = "praise_chr_site"
    if col not in df.columns:
        return {"status": "MISSING_COLUMN"}

    vals = df[col].astype(str)
    dash_count = int(vals.str.contains(r"-", regex=True).sum())

    old_strict = re.compile(r"^(chr[^_]+)_(\d+)(?:-(\d+))?$")
    greedy = re.compile(r"^(chr.+)_(\d+)(?:-(\d+))?$")
    strict_n = sum(bool(old_strict.match(v)) and old_strict.match(v).group(3) is not None for v in vals)
    greedy_n = sum(bool(greedy.match(v)) and greedy.match(v).group(3) is not None for v in vals)

    diff = []
    for v in vals:
        if "-" in v and not old_strict.match(v):
            diff.append(v)
    pd.DataFrame({"chr_site": sorted(set(diff))}).to_csv(
        META / "phase1a_praise_parser_review.tsv", sep="\t", index=False
    )

    return {
        "status": "RESOLVED" if greedy_n == dash_count else "REVIEW",
        "dash_interval_count": dash_count,
        "old_strict_interval_count": strict_n,
        "greedy_chr_interval_count": greedy_n,
        "historical_stage0_interval_count": 1357,
        "n_patterns_missed_by_old_regex": len(set(diff)),
        "interpretation": (
            "The old Phase0C regex did not allow chromosome/contig names containing underscores. "
            "Greedy parsing preserves the coordinate suffix while retaining the full chromosome token."
        ),
    }

# ---------------------------------------------------------------------
# Human 5-mer reconstruction + calibration transfer
# ---------------------------------------------------------------------

def normalize_motif_series(x: pd.Series) -> pd.Series:
    return x.fillna("").map(lambda z: motif5_token(z, allow_4mer=False))

def find_elap_sequence_window(df: pd.DataFrame) -> Tuple[List[str], Dict[str, Any]]:
    raw = [c for c in df.columns if c.startswith("raw__")]
    candidates = []
    for i in range(len(raw) - 4):
        cols = raw[i:i+5]
        fracs = []
        for c in cols:
            z = df[c].fillna("").astype(str).str.upper().str.strip()
            frac = z.isin(["A","C","G","T","U"]).mean()
            fracs.append(frac)
        center = df[cols[2]].fillna("").astype(str).str.upper().str.strip()
        center_u = center.isin(["T","U"]).mean()
        score = np.mean(fracs) + center_u
        if min(fracs) > .90 and center_u > .90:
            candidates.append((score, cols, fracs, center_u))
    if not candidates:
        raise RuntimeError("Could not identify five consecutive ELAP sequence-context columns.")
    candidates.sort(key=lambda x: x[0], reverse=True)
    score, cols, fracs, center_u = candidates[0]
    return cols, {
        "columns": cols,
        "base_fraction_each": fracs,
        "center_U_fraction": center_u,
    }

def prepare_human_context(path: Path, assay: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if assay == "BID":
        df["motif"] = normalize_motif_series(df["motif_5mer_reported"])
        info = {"source": "motif_5mer_reported"}
    else:
        cols, info = find_elap_sequence_window(df)
        seq = df[cols].fillna("").astype(str).agg("".join, axis=1)
        df["motif"] = normalize_motif_series(seq)
        info["source"] = "five consecutive raw sequence columns"

    df["single_u"] = (
        df["motif"].str.len().eq(5)
        & df["motif"].str[2].eq("U")
        & ~df["motif"].str[1].eq("U")
        & ~df["motif"].str[3].eq("U")
    )
    df["locus"] = (
        df["chrom"].astype(str) + ":" +
        df["strand"].astype(str) + ":" +
        df["pos1"].astype(str)
    )
    usable = df[df["motif"].str.len().eq(5)].copy()
    info["n_rows"] = len(df)
    info["n_usable_5mer"] = len(usable)
    info["n_single_u"] = int(df["single_u"].sum())
    return df, info

def build_cell_table(cell: str, bid_path: Path, elap_path: Path, contrast_map: Dict[str, float]) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    b, binfo = prepare_human_context(bid_path, "BID")
    e, einfo = prepare_human_context(elap_path, "ELAP")
    b = b[b["single_u"]].copy()
    e = e[e["single_u"]].copy()

    # Primary single-U locus table should be unique within each assay.
    if b["locus"].duplicated().any():
        raise RuntimeError(f"{cell}: BID duplicate exact single-U loci.")
    if e["locus"].duplicated().any():
        raise RuntimeError(f"{cell}: ELAP duplicate exact single-U loci.")

    bd = b.set_index("locus")
    ed = e.set_index("locus")
    loci = sorted(set(bd.index) | set(ed.index))
    rows = []
    shared_agree = 0
    shared_total = 0

    for loc in loci:
        ib = loc in bd.index
        ie = loc in ed.index
        mb = bd.loc[loc, "motif"] if ib else ""
        me = ed.loc[loc, "motif"] if ie else ""

        if ib and ie:
            cls = "shared"
            shared_total += 1
            if mb == me:
                shared_agree += 1
            motif = mb if mb == me else ""
        elif ib:
            cls = "BID_exclusive_reported"
            motif = mb
        else:
            cls = "ELAP_exclusive_reported"
            motif = me

        rows.append({
            "cell_line": cell,
            "locus": loc,
            "class": cls,
            "bid_reported": int(ib),
            "elap_reported": int(ie),
            "motif": motif,
            "bid_motif": mb,
            "elap_motif": me,
            "contrast": contrast_map.get(motif, np.nan),
        })

    out = pd.DataFrame(rows)
    agree = shared_agree / shared_total if shared_total else np.nan
    if shared_total >= 10 and agree < .95:
        raise RuntimeError(
            f"{cell}: exact shared-site motif agreement only {shared_agree}/{shared_total}; "
            "sequence orientation/coordinate semantics must be reviewed before inference."
        )

    info = {
        "cell_line": cell,
        "bid": binfo,
        "elap": einfo,
        "n_bid_single_u": len(b),
        "n_elap_single_u": len(e),
        "n_shared_single_u": shared_total,
        "shared_motif_agreement": agree,
        "class_counts": out["class"].value_counts().to_dict(),
        "n_with_calibration_contrast": int(out["contrast"].notna().sum()),
    }
    return out, info

def cell_analysis(df: pd.DataFrame, cell: str) -> Dict[str, Any]:
    x = df[
        df["class"].isin(["BID_exclusive_reported", "ELAP_exclusive_reported"])
        & df["contrast"].notna()
        & df["motif"].ne("")
    ].copy()
    x["y"] = (x["class"] == "ELAP_exclusive_reported").astype(int)
    obs_auc = auc_score(x["y"].to_numpy(), x["contrast"].to_numpy())
    lo, hi = cluster_boot_auc(x, 2000)
    p1, p2 = motif_permutation_p(x, obs_auc, 5000)

    X = np.column_stack([
        np.ones(len(x)),
        x["contrast"].to_numpy(float) / 0.1,
    ])
    beta, se, z = logistic_cluster(X, x["y"].to_numpy(), x["motif"].tolist())

    return {
        "cell_line": cell,
        "n": len(x),
        "n_BID_exclusive": int((x["y"] == 0).sum()),
        "n_ELAP_exclusive": int((x["y"] == 1).sum()),
        "n_unique_motifs": int(x["motif"].nunique()),
        "AUC_calibration_contrast": obs_auc,
        "AUC_cluster_boot_CI_low": lo,
        "AUC_cluster_boot_CI_high": hi,
        "motif_permutation_p_one_sided": p1,
        "motif_permutation_p_two_sided": p2,
        "logOR_per_0.1_contrast": beta[1],
        "cluster_robust_SE": se[1],
        "OR_per_0.1_contrast": math.exp(beta[1]),
        "OR_CI_low": math.exp(beta[1] - 1.96 * se[1]),
        "OR_CI_high": math.exp(beta[1] + 1.96 * se[1]),
        "wald_p": p_from_z(z[1]),
    }

def pooled_analysis(site: pd.DataFrame) -> pd.DataFrame:
    x = site[
        site["class"].isin(["BID_exclusive_reported", "ELAP_exclusive_reported"])
        & site["contrast"].notna()
        & site["motif"].ne("")
    ].copy()
    x["y"] = (x["class"] == "ELAP_exclusive_reported").astype(int)
    x["hela"] = (x["cell_line"] == "HeLa").astype(int)
    x["c10"] = x["contrast"].astype(float) / .1
    x["interaction"] = x["c10"] * x["hela"]
    X = np.column_stack([np.ones(len(x)), x["c10"], x["hela"], x["interaction"]])
    beta, se, z = logistic_cluster(X, x["y"].to_numpy(), x["motif"].tolist())
    terms = ["intercept", "calibration_contrast_per_0.1", "HeLa_indicator", "contrast_x_HeLa"]
    return pd.DataFrame({
        "term": terms,
        "beta": beta,
        "cluster_robust_SE": se,
        "z": z,
        "p": [p_from_z(v) for v in z],
        "OR": np.exp(beta),
        "OR_CI_low": np.exp(beta - 1.96*se),
        "OR_CI_high": np.exp(beta + 1.96*se),
    })

def continuum_summary(site: pd.DataFrame) -> pd.DataFrame:
    rows = []
    order = ["BID_exclusive_reported", "shared", "ELAP_exclusive_reported"]
    for cell in ["HEK293T", "HeLa"]:
        d = site[(site["cell_line"] == cell) & site["contrast"].notna()]
        for cls in order:
            z = d[d["class"] == cls]["contrast"].astype(float)
            rows.append({
                "cell_line": cell,
                "class": cls,
                "n": len(z),
                "median_contrast": z.median() if len(z) else np.nan,
                "q25": z.quantile(.25) if len(z) else np.nan,
                "q75": z.quantile(.75) if len(z) else np.nan,
            })
    return pd.DataFrame(rows)

# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------

def main() -> int:
    print("="*96)
    print("TRACE PHASE 1A — CALIBRATION RESOLUTION + FIRST SCIENTIFIC TEST")
    print("="*96)

    required = [
        "cal0c","bacs_fig1","elap_source","bid_hek","bid_hela",
        "elap_hek","elap_hela","praise"
    ]
    missing = [k for k in required if not PATHS[k].exists()]
    if missing:
        raise RuntimeError(f"Missing required files: {missing}")

    # Freeze primary calibration definitions BEFORE human outcomes are inspected.
    bid, bid_meta = resolve_bid(PATHS["cal0c"])
    bacs, bacs_meta = resolve_bacs(PATHS["bacs_fig1"])
    elap, elap_meta = resolve_elap(PATHS["elap_source"])

    mapping = [bid_meta, bacs_meta, elap_meta]
    pd.DataFrame(mapping).to_csv(META / "phase1a_calibration_resolution.tsv", sep="\t", index=False)

    praise_audit = audit_praise_interval(PATHS["praise"])
    (META / "phase1a_praise_parser_summary.json").write_text(
        json.dumps(praise_audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # Primary Aim1A needs BID + ELAP only. BACS calibration is a secondary resource
    # and therefore must not block the headline BID-vs-ELAP experiment.
    if elap is None or elap_meta["status"] != "RESOLVED":
        summary = {
            "status": "ELAP_CALIBRATION_MAPPING_REVIEW_REQUIRED",
            "calibration_resolution": mapping,
            "praise_parser": praise_audit,
            "scientific_note": (
                "Site normalization remains valid. Aim1A is intentionally not run until "
                "ELAP's 256-context synthetic enrichment is mapped without ambiguity."
            ),
        }
        (META / "phase1a_summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        print("STATUS: ELAP_CALIBRATION_MAPPING_REVIEW_REQUIRED")
        print(r"Return 00_meta\phase1a_calibration_resolution.tsv and phase1a_summary.json")
        return 0

    # Assemble resolved calibration.
    bid = bid.copy()
    bid["response_percentile"] = pct_rank(bid["response"])
    bid["assay"] = "BID"

    elap = elap.copy()
    elap["response_percentile"] = pct_rank(elap["response"])
    elap["assay"] = "ELAP"

    resolved = pd.concat([bid, elap], ignore_index=True)
    if bacs is not None and bacs_meta["status"] == "RESOLVED":
        bacs = bacs.copy()
        bacs["response_percentile"] = pct_rank(bacs["response"])
        bacs["assay"] = "BACS"
        resolved = pd.concat([resolved, bacs], ignore_index=True)

    cal_out = AIM1 / "calibration_resolved.tsv"
    resolved.to_csv(cal_out, sep="\t", index=False)

    # All 256 motifs must exist in BID and ELAP.
    br = bid.set_index("motif")["response_percentile"]
    er = elap.set_index("motif")["response_percentile"]
    common = sorted(set(br.index) & set(er.index))
    if len(common) != 256:
        raise RuntimeError(f"BID-ELAP calibration common motifs={len(common)}, expected 256.")
    contrast_map = {m: float(er[m] - br[m]) for m in common}

    # Build exact single-U positive-union tables.
    hek, hek_info = build_cell_table("HEK293T", PATHS["bid_hek"], PATHS["elap_hek"], contrast_map)
    hela, hela_info = build_cell_table("HeLa", PATHS["bid_hela"], PATHS["elap_hela"], contrast_map)
    site = pd.concat([hek, hela], ignore_index=True)
    site.to_csv(AIM1 / "aim1a_site_table.tsv", sep="\t", index=False)

    # Primary independent tests in both biological contexts.
    cellres = pd.DataFrame([
        cell_analysis(site[site["cell_line"] == "HEK293T"], "HEK293T"),
        cell_analysis(site[site["cell_line"] == "HeLa"], "HeLa"),
    ])
    cellres.to_csv(AIM1 / "aim1a_cell_results.tsv", sep="\t", index=False)

    pooled = pooled_analysis(site)
    pooled.to_csv(AIM1 / "aim1a_pooled_logit.tsv", sep="\t", index=False)

    continuum = continuum_summary(site)
    continuum.to_csv(AIM1 / "aim1a_continuum.tsv", sep="\t", index=False)

    # Frozen interpretation gate: do not invent a single arbitrary "success score".
    # Report evidence components independently.
    hekr = cellres[cellres["cell_line"] == "HEK293T"].iloc[0]
    helar = cellres[cellres["cell_line"] == "HeLa"].iloc[0]
    direction_consistent = (
        hekr["AUC_calibration_contrast"] > .5
        and helar["AUC_calibration_contrast"] > .5
        and hekr["logOR_per_0.1_contrast"] > 0
        and helar["logOR_per_0.1_contrast"] > 0
    )
    both_ci_above_half = (
        hekr["AUC_cluster_boot_CI_low"] > .5
        and helar["AUC_cluster_boot_CI_low"] > .5
    )
    interaction = pooled.loc[pooled["term"] == "contrast_x_HeLa"].iloc[0]

    status = (
        "AIM1A_STRONG_CALIBRATION_TRANSFER"
        if direction_consistent and both_ci_above_half
        else "AIM1A_RESULT_READY_FOR_SCIENTIFIC_INTERPRETATION"
    )

    contract = {
        "phase": "1A",
        "primary_question": (
            "Does assay-intrinsic synthetic sequence response explain BID-vs-ELAP "
            "reported-map discordance in human cells?"
        ),
        "primary_assays": ["BID", "ELAP"],
        "BID_response": "fit_R (published BID-induced deletion ratio)",
        "ELAP_response": elap_meta,
        "BACS_calibration_role": "secondary resource; not required to declare Aim1A primary result",
        "normalization": "within-assay percentile; raw assay units are never subtracted",
        "human_universe": "union of sites reported by BID or ELAP, restricted to locally single-U exact loci",
        "labels": {
            "0": "BID-exclusive reported",
            "1": "ELAP-exclusive reported",
            "shared": "held out from fitting; continuum check only",
        },
        "inference": {
            "primary_discrimination": "AUC of fixed synthetic calibration contrast",
            "CI": "motif-cluster bootstrap",
            "null_test": "permute calibration contrast across 5-mer motifs",
            "effect": "logistic OR per 0.1 contrast with motif-cluster-robust SE",
            "cross_context": "HEK293T and HeLa analyzed independently; pooled interaction tests heterogeneity",
        },
        "not_done": [
            "no sequence classifier",
            "no feature tuning",
            "no foundation model",
            "no raw sequencing",
            "no missing call interpreted as biological negative",
        ],
    }
    (META / "phase1a_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    summary = {
        "status": status,
        "calibration_resolution": mapping,
        "praise_parser": praise_audit,
        "human_context_QC": [hek_info, hela_info],
        "cell_results": cellres.to_dict(orient="records"),
        "pooled_interaction": interaction.to_dict(),
        "direction_consistent": bool(direction_consistent),
        "both_cluster_bootstrap_CI_above_0.5": bool(both_ci_above_half),
        "next_step": (
            "If Aim1A supports calibration transfer, proceed directly to chemistry-aware "
            "canonical U-run construction and Aim1B HEK<->HeLa sequence-fingerprint transfer. "
            "If not, do not model-hunt; interpret the calibration result before expanding."
        ),
    }
    (META / "phase1a_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print()
    print("Calibration resolution:")
    for m in mapping:
        print(f"  {m['assay']:5s} {m['status']:32s} n={m.get('n_motifs',0)}  {m.get('response_source','')}")
    print()
    print("Aim1A:")
    print(cellres.to_string(index=False))
    print()
    print("Pooled interaction:")
    print(pooled.to_string(index=False))
    print()
    print("PRAISE parser:", praise_audit["status"],
          "intervals:", praise_audit.get("greedy_chr_interval_count"))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase1a_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase1a_calibration_resolution.tsv")
    print(r"  D:\RNA\Trace\00_meta\phase1a_praise_parser_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase1a_analysis_contract.json")
    print(r"  D:\RNA\Trace\04_aim1\calibration_resolved.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1a_cell_results.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1a_pooled_logit.tsv")
    print(r"  D:\RNA\Trace\04_aim1\aim1a_continuum.tsv")
    print()
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase1a_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase1a_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
