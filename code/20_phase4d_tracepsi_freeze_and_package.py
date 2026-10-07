#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE-Ψ Phase 4D — final freeze, robustness audit, and tool packaging
====================================================================

No model fitting or tuning occurs here.

This phase:
1) verifies that Phase4C externally validated the frozen P2 model;
2) performs a no-retraining post-lockbox robustness audit excluding canonical
   PUS7-like motifs to test whether external ranking is solely a known-motif
   effect;
3) freezes the model card and validation summary;
4) creates a lightweight user-facing TRACE-Ψ scorer bundle.

Project root: D:\RNA\Trace
"""

from __future__ import annotations

import json
import re
import shutil
import traceback
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(r"D:\RNA\Trace")
CODE = ROOT / "code"
META = ROOT / "00_meta"
TRACEPSI = ROOT / "07_tracepsi"
TOOL = ROOT / "08_tool" / "TRACEpsi"
LOGS = ROOT / "logs"

for p in [META, TRACEPSI, TOOL, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

S4B = META / "phase4b_v2_summary.json"
C4B = META / "phase4b_v2_analysis_contract.json"
S4C = META / "phase4c_summary.json"
C4C = META / "phase4c_analysis_contract.json"

MODEL = TRACEPSI / "phase4b_v2_TRACEpsi_development_model.joblib"
FP = TRACEPSI / "phase4a_fingerprint_HEK_for_HeLa.tsv"
HELA_REF = TRACEPSI / "phase4b_v2_TRACEpsi_HeLa_scorecard.tsv.gz"
A549_SCORECARD = TRACEPSI / "phase4c_A549_TRACEpsi_scorecard.tsv.gz"

SCORER_SOURCE = CODE / "21_tracepsi_score.py"

N_BOOT = 5000
SEED = 20261004
PUS7_RE = re.compile(r"^U[ACGU]UA[AG]$")


def safe_auc(y, s):
    return float(roc_auc_score(np.asarray(y, int), np.asarray(s, float)))


def bootstrap_auc(y, s):
    y = np.asarray(y, int)
    s = np.asarray(s, float)
    rng = np.random.default_rng(SEED)
    vals = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(y), len(y))
        yy = y[idx]
        if len(np.unique(yy)) != 2:
            continue
        vals.append(safe_auc(yy, s[idx]))
    return {
        "CI_low": float(np.quantile(vals, .025)),
        "CI_high": float(np.quantile(vals, .975)),
        "bootstrap_support_AUC_gt_0.5": float(np.mean(np.asarray(vals) > .5)),
        "n_boot": int(len(vals)),
    }


def main():
    print("=" * 100)
    print("TRACE-Ψ PHASE 4D — FINAL FREEZE + TOOL PACKAGE")
    print("=" * 100)

    for p in [S4B, C4B, S4C, C4C, MODEL, FP, HELA_REF, A549_SCORECARD, SCORER_SOURCE]:
        if not p.exists():
            raise RuntimeError(f"Missing required frozen asset: {p}")

    with open(S4B, encoding="utf-8") as f:
        s4b = json.load(f)
    with open(S4C, encoding="utf-8") as f:
        s4c = json.load(f)

    if s4b.get("status") != "TRACEPSI_P2_READY_FOR_A549_LOCKBOX":
        raise RuntimeError("Phase4B did not freeze P2 for lockbox.")
    if s4c.get("status") != "TRACEPSI_EXTERNALLY_VALIDATED_A549":
        raise RuntimeError(
            f"Phase4C status is {s4c.get('status')}; expected external validation."
        )

    # ------------------------------------------------------------------
    # No-retraining canonical-motif robustness audit
    # ------------------------------------------------------------------
    a = pd.read_csv(A549_SCORECARD, sep="\t", compression="gzip", low_memory=False)
    a["PUS7_like"] = a["motif"].astype(str).map(lambda x: bool(PUS7_RE.match(x)))

    audit_rows = []
    for label, z in [
        ("all_A549", a),
        ("A549_excluding_PUS7_like", a.loc[~a["PUS7_like"]].copy()),
        ("A549_PUS7_like_only", a.loc[a["PUS7_like"]].copy()),
    ]:
        if z["DRS_reported"].nunique() < 2:
            continue
        y = z["DRS_reported"].astype(int).to_numpy()
        s = z["TRACEpsi_score"].to_numpy(float)
        b = bootstrap_auc(y, s)
        audit_rows.append({
            "subset": label,
            "n": int(len(z)),
            "DRS_positive_n": int(y.sum()),
            "prevalence": float(y.mean()),
            "AUROC": safe_auc(y, s),
            "average_precision": float(average_precision_score(y, s)),
            **b,
        })

    audit = pd.DataFrame(audit_rows)
    audit.to_csv(
        TRACEPSI / "phase4d_external_motif_robustness.tsv",
        sep="\t", index=False
    )

    non_pus = audit[audit["subset"].eq("A549_excluding_PUS7_like")]
    if len(non_pus):
        non_pus_supported = bool(
            non_pus.iloc[0]["AUROC"] > .5
            and non_pus.iloc[0]["bootstrap_support_AUC_gt_0.5"] >= .90
        )
    else:
        non_pus_supported = False

    # ------------------------------------------------------------------
    # Freeze reference scores for future user-facing percentiles
    # ------------------------------------------------------------------
    ref = pd.read_csv(HELA_REF, sep="\t", compression="gzip", low_memory=False)
    if "TRACEpsi_score" not in ref.columns:
        raise RuntimeError("HeLa reference scorecard lacks TRACEpsi_score.")
    ref[["TRACEpsi_score", "source_breadth", "source_pattern"]].to_csv(
        TOOL / "TRACEpsi_reference_scores.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    # Reduce/copy frozen assets.
    shutil.copy2(MODEL, TOOL / "TRACEpsi_frozen_model.joblib")
    shutil.copy2(FP, TOOL / "TRACEpsi_fingerprint_HEK.tsv")
    shutil.copy2(SCORER_SOURCE, TOOL / "tracepsi_score.py")

    # ------------------------------------------------------------------
    # Model card
    # ------------------------------------------------------------------
    p2_metric = next(
        x for x in s4c["lockbox_metrics"]
        if x["score"] == "P2_TRACEpsi_frozen"
    )
    context_metric = next(
        x for x in s4c["lockbox_metrics"]
        if x["score"] == "A_context_only_frozen"
    )

    model_card = {
        "tool": "TRACE-Ψ",
        "title": (
            "A technology-aware framework for cross-platform evidence "
            "portability in pseudouridine mapping"
        ),
        "version": "1.0-frozen",
        "intended_use": (
            "Rank already-reported Ψ loci for orthogonal cross-platform validation."
        ),
        "estimand": (
            "Relative propensity of source-reported evidence to be re-observed by "
            "another profiling technology."
        ),
        "not_intended_as": [
            "P(true Ψ)",
            "biological truth score",
            "causal functional score",
            "replacement for orthogonal experimental validation",
        ],
        "frozen_model": "P2_TRACEpsi",
        "development": {
            "cell_context": "HeLa",
            "source/target roles": "BID, BACS, ELAP LOTO",
            "development_external_platform": "HeLa DRS",
            "breadth1_DRS_AUROC": 0.644274331094766,
            "breadth1_DRS_AUROC_CI": [0.6002652651285907, 0.6882343777873113],
            "breadth1_DRS_Lift10": 1.7499008442404667,
        },
        "external_A549_lockbox": {
            "source": "BID",
            "target": "native DRS / Mod-p ID",
            "n": int(p2_metric["n"]),
            "positives": int(p2_metric["positives"]),
            "prevalence": float(p2_metric["prevalence"]),
            "AUROC": float(p2_metric["AUROC"]),
            "AUROC_CI": [
                float(p2_metric["AUROC_CI_low"]),
                float(p2_metric["AUROC_CI_high"]),
            ],
            "average_precision": float(p2_metric["average_precision"]),
            "Lift10_tieaware": float(p2_metric["Lift10_tieaware"]),
            "budget50_expected_precision": float(
                p2_metric["budget50_expected_precision"]
            ),
            "budget100_expected_precision": float(
                p2_metric["budget100_expected_precision"]
            ),
        },
        "external_ablation_note": {
            "context_only_AUROC": float(context_metric["AUROC"]),
            "P2_minus_context_AUROC": float(
                s4c["lockbox_delta_AUC"][0]["delta_AUROC"]
            ),
            "interpretation": (
                "The full P2 model externally validated as a portability ranker, "
                "but its learned assay-fingerprint component did not outperform "
                "generic local-sequence context in A549. Do not claim external "
                "validation of fingerprint-specific gain."
            ),
        },
        "post_lockbox_no_retraining_audit": {
            "canonical_PUS7_exclusion_supported": non_pus_supported,
            "table": audit.to_dict(orient="records"),
            "interpretation": (
                "Used only to characterize robustness of the already-frozen score; "
                "no model or threshold was changed."
            ),
        },
        "score_output": {
            "TRACEpsi_score": "raw frozen logistic ranking score",
            "TRACEpsi_reference_percentile": (
                "empirical percentile within the frozen HeLa reference distribution "
                "at the same source breadth when possible"
            ),
            "TRACEpsi_validation_priority": (
                "relative experimental-validation priority; not a probability"
            ),
        },
        "priority_bands": {
            "high_validation_priority": "reference percentile >=0.90",
            "elevated_validation_priority": "0.75-0.90",
            "intermediate_validation_priority": "0.25-0.75",
            "lower_validation_priority": "<0.25",
        },
        "limitations": [
            "Validated on processed published call tables rather than uniformly reprocessed raw data.",
            "Primary benchmark is restricted to locally-single-U exact loci.",
            "A549 external benchmark tests BID-to-DRS portability; broader source/platform generalization remains less established.",
            "DRS non-reporting is not biological absence.",
            "Scores are ordinal/ranking-oriented and should not be interpreted as calibrated truth probabilities.",
        ],
    }
    (TOOL / "MODEL_CARD.json").write_text(
        json.dumps(model_card, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    readme = """TRACE-Ψ v1.0
============

Purpose
-------
Prioritize already-reported pseudouridine loci for orthogonal cross-platform
validation.

TRACE-Ψ does NOT estimate whether a locus is a biologically true Ψ site.

Input
-----
TSV with:
  motif
and either:
  source_assay
or:
  source_BID, source_BACS, source_ELAP

Optional genomic annotation columns are preserved.

Run
---
python tracepsi_score.py input.tsv output.tsv

Output
------
TRACEpsi_score
TRACEpsi_reference_percentile
TRACEpsi_within_input_percentile
TRACEpsi_validation_priority

Interpretation
--------------
Higher scores indicate greater cross-platform evidence portability according
to the frozen model. Use scores for ranking/prioritization, not as absolute
probabilities.

External validation
-------------------
Frozen model validated on A549 BID -> native DRS / Mod-p ID:
AUROC ~0.763; top-tail enrichment was strong.

See MODEL_CARD.json for full scope and limitations.
"""
    (TOOL / "README.txt").write_text(readme, encoding="utf-8")

    summary = {
        "status": "TRACEPSI_V1_FROZEN_AND_PACKAGED",
        "tool_dir": str(TOOL),
        "external_motif_robustness": audit.to_dict(orient="records"),
        "non_PUS7_external_signal_supported": non_pus_supported,
        "frozen_assets": [
            str(TOOL / "TRACEpsi_frozen_model.joblib"),
            str(TOOL / "TRACEpsi_fingerprint_HEK.tsv"),
            str(TOOL / "TRACEpsi_reference_scores.tsv.gz"),
            str(TOOL / "tracepsi_score.py"),
            str(TOOL / "MODEL_CARD.json"),
            str(TOOL / "README.txt"),
        ],
        "next_step": (
            "Do not tune the model further. Next work should be manuscript figures, "
            "a small example input/output, and optional Inference Audit packaging "
            "using already-completed Aim3 results."
        ),
    }
    (META / "phase4d_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    print()
    print("STATUS: TRACEPSI_V1_FROZEN_AND_PACKAGED")
    print()
    print("FINAL EXTERNAL ROBUSTNESS AUDIT")
    print(audit.to_string(index=False))
    print()
    print("Tool bundle:", TOOL)
    print("Run scorer with:")
    print(r"  python D:\RNA\Trace\08_tool\TRACEpsi\tracepsi_score.py input.tsv output.tsv")
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase4d_summary.json")
    print(r"  D:\RNA\Trace\08_tool\TRACEpsi\MODEL_CARD.json")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase4d_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        print("\nFATAL ERROR\n")
        print(err)
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
