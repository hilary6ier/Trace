#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE-Ψ user-facing scorer
==========================

Ranks already-reported pseudouridine evidence by cross-platform re-observation
propensity using the frozen Phase4B-v2 P2 model.

This score is NOT P(true Ψ), NOT a biological-confidence probability, and NOT a
causal quantity. It is intended to prioritize loci for orthogonal validation.

Required input TSV:
Either
  A) motif + source_assay
     where source_assay is BID, BACS, or ELAP
or
  B) motif + source_BID + source_BACS + source_ELAP
     using 0/1 indicators.

Optional columns such as chrom, pos1, strand, gene, etc. are preserved.

Example
-------
python tracepsi_score.py input.tsv output.tsv

The tool bundle should contain:
  TRACEpsi_frozen_model.joblib
  TRACEpsi_fingerprint_HEK.tsv
  TRACEpsi_reference_scores.tsv.gz
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any, Dict

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


def motif5(x: Any) -> str:
    z = ss(x).upper().replace("Ψ", "U").replace("T", "U")
    z = re.sub(r"[^ACGU]", "", z)
    return z if len(z) == 5 and z[2] == "U" else ""


def add_source_fields(df: pd.DataFrame) -> pd.DataFrame:
    z = df.copy()

    if {"source_BID", "source_BACS", "source_ELAP"}.issubset(z.columns):
        for c in ["source_BID", "source_BACS", "source_ELAP"]:
            z[c] = pd.to_numeric(z[c], errors="coerce").fillna(0).astype(int)
            if not z[c].isin([0, 1]).all():
                raise ValueError(f"{c} must contain only 0/1.")
    elif "source_assay" in z.columns:
        a = z["source_assay"].astype(str).str.upper().str.strip()
        bad = ~a.isin(["BID", "BACS", "ELAP"])
        if bad.any():
            raise ValueError(
                "source_assay must be BID, BACS or ELAP. "
                f"Bad values: {sorted(a[bad].unique())[:10]}"
            )
        z["source_BID"] = a.eq("BID").astype(int)
        z["source_BACS"] = a.eq("BACS").astype(int)
        z["source_ELAP"] = a.eq("ELAP").astype(int)
    else:
        raise ValueError(
            "Input must contain either source_assay or all of "
            "source_BID/source_BACS/source_ELAP."
        )

    z["source_breadth"] = (
        z["source_BID"] + z["source_BACS"] + z["source_ELAP"]
    )

    if (z["source_breadth"] < 1).any():
        raise ValueError("Every row must have >=1 source technology.")

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
        raise ValueError("Input requires a 5-mer motif column named 'motif'.")

    z["motif"] = z["motif"].map(motif5)
    if z["motif"].eq("").any():
        bad = z.index[z["motif"].eq("")].tolist()[:10]
        raise ValueError(
            f"Invalid motif(s) at rows {bad}. Motif must be a 5-mer with central U."
        )

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
            z[f"motif_{off:+d}_{base}"] = (
                z["motif"].str[idx[off]].eq(base).astype(int)
            )

    return z


def decision_function(bundle: Dict, df: pd.DataFrame) -> np.ndarray:
    features = list(bundle["features"])
    missing = sorted(set(features) - set(df.columns))
    if missing:
        raise ValueError(f"Missing frozen model features: {missing}")

    X = df[features].to_numpy(float)
    Xs = (X - np.asarray(bundle["mean"])) / np.asarray(bundle["scale"])
    return bundle["model"].decision_function(Xs)


def reference_percentile(
    score: float,
    breadth: int,
    ref: pd.DataFrame,
) -> float:
    r = ref.loc[ref["source_breadth"].eq(breadth), "TRACEpsi_score"].to_numpy(float)
    if len(r) < 20:
        r = ref["TRACEpsi_score"].to_numpy(float)
    # midpoint empirical percentile for ties
    lt = np.sum(r < score)
    eq = np.sum(r == score)
    return float((lt + 0.5 * eq) / len(r))


def priority_band(p: float) -> str:
    if p >= 0.90:
        return "high_validation_priority"
    if p >= 0.75:
        return "elevated_validation_priority"
    if p >= 0.25:
        return "intermediate_validation_priority"
    return "lower_validation_priority"


def main():
    ap = argparse.ArgumentParser(
        description="TRACE-Ψ cross-platform evidence-portability scorer"
    )
    ap.add_argument("input_tsv")
    ap.add_argument("output_tsv")
    ap.add_argument(
        "--bundle-dir",
        default=str(Path(__file__).resolve().parent),
        help="Directory containing frozen TRACE-Ψ model assets.",
    )
    args = ap.parse_args()

    bundle_dir = Path(args.bundle_dir)
    model_path = bundle_dir / "TRACEpsi_frozen_model.joblib"
    fp_path = bundle_dir / "TRACEpsi_fingerprint_HEK.tsv"
    ref_path = bundle_dir / "TRACEpsi_reference_scores.tsv.gz"

    for p in [model_path, fp_path, ref_path]:
        if not p.exists():
            raise FileNotFoundError(f"Missing TRACE-Ψ asset: {p}")

    artifact = joblib.load(model_path)
    if artifact.get("selected_complexity_for_lockbox") != "P2_TRACEpsi":
        raise RuntimeError("Bundle does not contain the frozen P2 TRACE-Ψ model.")

    model = artifact["fitted_models"]["P2_TRACEpsi"]
    fingerprint = pd.read_csv(fp_path, sep="\t")
    ref = pd.read_csv(ref_path, sep="\t", compression="gzip")

    x = pd.read_csv(args.input_tsv, sep="\t", low_memory=False)
    x = add_source_fields(x)
    x = add_sequence_features(x, fingerprint)

    x["TRACEpsi_score"] = decision_function(model, x)
    x["TRACEpsi_reference_percentile"] = [
        reference_percentile(float(s), int(b), ref)
        for s, b in zip(x["TRACEpsi_score"], x["source_breadth"])
    ]
    x["TRACEpsi_within_input_percentile"] = (
        x["TRACEpsi_score"].rank(method="average", pct=True)
    )
    x["TRACEpsi_validation_priority"] = (
        x["TRACEpsi_reference_percentile"].map(priority_band)
    )
    x["TRACEpsi_score_semantics"] = (
        "cross-platform evidence-portability rank; not P(true_psi)"
    )

    out = Path(args.output_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    x.to_csv(out, sep="\t", index=False)

    print("TRACE-Ψ scoring complete")
    print("Rows:", len(x))
    print("High validation priority:",
          int(x["TRACEpsi_validation_priority"].eq("high_validation_priority").sum()))
    print("Output:", out)


if __name__ == "__main__":
    main()
