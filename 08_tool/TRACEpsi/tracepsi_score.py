#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE-Ψ v1.1 user-facing scorer
================================

Purpose
-------
Prioritize ALREADY-REPORTED pseudouridine evidence for orthogonal
cross-platform validation.

TRACE-Ψ estimates a RELATIVE EVIDENCE-PORTABILITY RANK.
It does NOT estimate P(true Ψ), biological truth, or causal function.

Validated domain
----------------
Primary model validation was performed on locally-single-U 5-mers:
    central base = U
    immediate flanking bases (-1,+1) are not U

Rows outside this domain are retained in the output but are NOT scored.

Required input
--------------
A TSV containing:
  motif
and either:
  source_assay
or:
  source_BID, source_BACS, source_ELAP

source_assay supports BID, BACS, ELAP.

Optional columns such as chrom/pos1/strand/gene are preserved.

Outputs
-------
TRACEpsi_score
TRACEpsi_reference_percentile_same_pattern
TRACEpsi_reference_percentile_same_breadth
TRACEpsi_reference_percentile
TRACEpsi_within_input_percentile
TRACEpsi_evidence_layer
TRACEpsi_domain_status
TRACEpsi_interpretation

Evidence layers
---------------
multi_technology_supported:
    source breadth >= 2. Orthogonal source support itself is strong evidence;
    TRACE-Ψ is supplementary.

single_source_ranked:
    source breadth = 1. This is the main intended use of TRACE-Ψ ranking.

out_of_validated_domain:
    score is not reported.

No fixed "truth-confidence" tiers are produced. Users should rank candidates by
percentile or by their experimental budget (top K).

Usage
-----
python tracepsi_score.py input.tsv output.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any, Dict, Tuple

import joblib
import numpy as np
import pandas as pd


BASES = "ACGU"
OFFSETS = [-2, -1, 1, 2]


def ss(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    return str(x).strip()


def normalize_motif(x: Any) -> str:
    z = ss(x).upper().replace("Ψ", "U").replace("T", "U")
    z = re.sub(r"[^ACGU]", "", z)
    return z


def motif_domain_status(m: str) -> str:
    if len(m) != 5:
        return "out_of_validated_domain_invalid_length"
    if m[2] != "U":
        return "out_of_validated_domain_center_not_U"
    if m[1] == "U" or m[3] == "U":
        return "out_of_validated_domain_consecutive_U"
    return "in_validated_domain"


def add_source_fields(df: pd.DataFrame) -> pd.DataFrame:
    z = df.copy()

    if {"source_BID", "source_BACS", "source_ELAP"}.issubset(z.columns):
        for c in ["source_BID", "source_BACS", "source_ELAP"]:
            vals = pd.to_numeric(z[c], errors="coerce")
            if vals.isna().any():
                raise ValueError(f"{c} contains non-numeric/missing values.")
            vals = vals.astype(int)
            if not vals.isin([0, 1]).all():
                raise ValueError(f"{c} must contain only 0/1.")
            z[c] = vals
    elif "source_assay" in z.columns:
        a = z["source_assay"].astype(str).str.upper().str.strip()
        bad = ~a.isin(["BID", "BACS", "ELAP"])
        if bad.any():
            raise ValueError(
                "source_assay must be BID, BACS, or ELAP. "
                f"Bad values: {sorted(a[bad].unique())[:10]}"
            )
        z["source_BID"] = a.eq("BID").astype(int)
        z["source_BACS"] = a.eq("BACS").astype(int)
        z["source_ELAP"] = a.eq("ELAP").astype(int)
    else:
        raise ValueError(
            "Input must contain source_assay OR all of "
            "source_BID/source_BACS/source_ELAP."
        )

    z["source_breadth"] = (
        z["source_BID"] + z["source_BACS"] + z["source_ELAP"]
    )
    if (z["source_breadth"] < 1).any():
        raise ValueError("Every row must have at least one source technology.")
    if (z["source_breadth"] > 3).any():
        raise ValueError("source_breadth cannot exceed 3.")

    z["source_pattern"] = z.apply(
        lambda r: "+".join(
            a for a in ["BID", "BACS", "ELAP"]
            if int(r[f"source_{a}"]) == 1
        ),
        axis=1,
    )
    return z


def add_sequence_features(
    df: pd.DataFrame,
    fingerprint: pd.DataFrame,
) -> pd.DataFrame:
    z = df.copy()

    if "motif" not in z.columns:
        raise ValueError("Input requires a column named 'motif'.")

    z["motif"] = z["motif"].map(normalize_motif)
    z["TRACEpsi_domain_status"] = z["motif"].map(motif_domain_status)

    fp = fingerprint.set_index("motif")[
        "fingerprint_logratio_ELAP_vs_BID"
    ].to_dict()

    z["fp_crosscell_ELAPvsBID"] = z["motif"].map(fp)
    z["fp_seen_in_HEK"] = z["fp_crosscell_ELAPvsBID"].notna().astype(int)
    z["fp_crosscell_ELAPvsBID"] = (
        z["fp_crosscell_ELAPvsBID"].fillna(0.0)
    )
    z["fp_abs"] = z["fp_crosscell_ELAPvsBID"].abs()
    z["fp_source_alignment"] = (
        z["fp_crosscell_ELAPvsBID"]
        * (z["source_ELAP"].astype(float) - z["source_BID"].astype(float))
    )

    idx = {-2: 0, -1: 1, 1: 3, 2: 4}
    for off in OFFSETS:
        for base in BASES:
            vals = []
            for m in z["motif"]:
                if len(m) == 5:
                    vals.append(int(m[idx[off]] == base))
                else:
                    vals.append(0)
            z[f"motif_{off:+d}_{base}"] = vals

    return z


def decision_function(bundle: Dict, df: pd.DataFrame) -> np.ndarray:
    features = list(bundle["features"])
    missing = sorted(set(features) - set(df.columns))
    if missing:
        raise ValueError(f"Missing frozen model features: {missing}")

    X = df[features].to_numpy(float)
    Xs = (X - np.asarray(bundle["mean"])) / np.asarray(bundle["scale"])
    return bundle["model"].decision_function(Xs)


def empirical_percentile(score: float, values: np.ndarray) -> float:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    lt = np.sum(values < score)
    eq = np.sum(values == score)
    return float((lt + 0.5 * eq) / len(values))


def reference_percentiles(
    row: pd.Series,
    ref: pd.DataFrame,
) -> Tuple[float, float, float, str]:
    score = float(row["TRACEpsi_score"])
    pattern = str(row["source_pattern"])
    breadth = int(row["source_breadth"])

    r_pattern = ref.loc[
        ref["source_pattern"].astype(str).eq(pattern),
        "TRACEpsi_score",
    ].to_numpy(float)

    r_breadth = ref.loc[
        ref["source_breadth"].astype(int).eq(breadth),
        "TRACEpsi_score",
    ].to_numpy(float)

    p_pattern = (
        empirical_percentile(score, r_pattern)
        if len(r_pattern) >= 20 else np.nan
    )
    p_breadth = (
        empirical_percentile(score, r_breadth)
        if len(r_breadth) >= 20 else empirical_percentile(
            score, ref["TRACEpsi_score"].to_numpy(float)
        )
    )

    if np.isfinite(p_pattern):
        return p_pattern, p_breadth, p_pattern, "same_source_pattern"
    return p_pattern, p_breadth, p_breadth, "same_source_breadth"


def main():
    ap = argparse.ArgumentParser(
        description="TRACE-Ψ v1.1 cross-platform evidence-portability scorer"
    )
    ap.add_argument("input_tsv")
    ap.add_argument("output_tsv")
    ap.add_argument(
        "--bundle-dir",
        default=str(Path(__file__).resolve().parent),
        help="Directory containing frozen TRACE-Ψ model assets.",
    )
    args = ap.parse_args()

    bundle = Path(args.bundle_dir)
    model_path = bundle / "TRACEpsi_frozen_model.joblib"
    fp_path = bundle / "TRACEpsi_fingerprint_HEK.tsv"
    ref_path = bundle / "TRACEpsi_reference_scores.tsv.gz"

    for p in [model_path, fp_path, ref_path]:
        if not p.exists():
            raise FileNotFoundError(f"Missing TRACE-Ψ asset: {p}")

    artifact = joblib.load(model_path)
    if artifact.get("selected_complexity_for_lockbox") != "P2_TRACEpsi":
        raise RuntimeError("Bundle is not the frozen P2 TRACE-Ψ release model.")

    model = artifact["fitted_models"]["P2_TRACEpsi"]
    fp = pd.read_csv(fp_path, sep="\t")
    ref = pd.read_csv(ref_path, sep="\t", compression="gzip")

    required_ref = {"TRACEpsi_score", "source_breadth", "source_pattern"}
    if not required_ref.issubset(ref.columns):
        raise RuntimeError(
            f"Reference score file lacks {sorted(required_ref-set(ref.columns))}"
        )

    x = pd.read_csv(args.input_tsv, sep="\t", low_memory=False)
    x = add_source_fields(x)
    x = add_sequence_features(x, fp)

    x["TRACEpsi_score"] = np.nan

    eligible = x["TRACEpsi_domain_status"].eq("in_validated_domain")
    if eligible.any():
        x.loc[eligible, "TRACEpsi_score"] = decision_function(
            model, x.loc[eligible]
        )

    x["TRACEpsi_reference_percentile_same_pattern"] = np.nan
    x["TRACEpsi_reference_percentile_same_breadth"] = np.nan
    x["TRACEpsi_reference_percentile"] = np.nan
    x["TRACEpsi_percentile_reference_basis"] = ""

    for i in x.index[eligible]:
        pp, pb, p, basis = reference_percentiles(x.loc[i], ref)
        x.at[i, "TRACEpsi_reference_percentile_same_pattern"] = pp
        x.at[i, "TRACEpsi_reference_percentile_same_breadth"] = pb
        x.at[i, "TRACEpsi_reference_percentile"] = p
        x.at[i, "TRACEpsi_percentile_reference_basis"] = basis

    # Within-input percentile only among validated-domain rows with same source pattern.
    x["TRACEpsi_within_input_percentile"] = np.nan
    for pattern, idx in x.loc[eligible].groupby("source_pattern").groups.items():
        vals = x.loc[idx, "TRACEpsi_score"]
        x.loc[idx, "TRACEpsi_within_input_percentile"] = (
            vals.rank(method="average", pct=True)
        )

    x["TRACEpsi_evidence_layer"] = np.where(
        ~eligible,
        "out_of_validated_domain",
        np.where(
            x["source_breadth"].ge(2),
            "multi_technology_supported",
            "single_source_ranked",
        ),
    )

    x["TRACEpsi_interpretation"] = np.where(
        ~eligible,
        "not_scored_outside_validated_locally-single-U_domain",
        np.where(
            x["source_breadth"].ge(2),
            "orthogonal_source_support_present; portability score is supplementary",
            "use portability percentile/rank to prioritize orthogonal validation",
        ),
    )

    out = Path(args.output_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    x.to_csv(out, sep="\t", index=False)

    print("TRACE-Ψ v1.1 scoring complete")
    print("Rows:", len(x))
    print("Scored in validated domain:", int(eligible.sum()))
    print("Out of validated domain:", int((~eligible).sum()))
    print("Single-source ranked:", int(
        x["TRACEpsi_evidence_layer"].eq("single_source_ranked").sum()
    ))
    print("Multi-technology supported:", int(
        x["TRACEpsi_evidence_layer"].eq("multi_technology_supported").sum()
    ))
    print("Output:", out)


if __name__ == "__main__":
    main()
