#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TRACE-Ψ Phase 4E — release hardening (NO model change)
======================================================

This phase does NOT fit, tune, recalibrate or replace the frozen P2 model.

It only improves the software/release layer:
1) replaces the v1 scorer with domain-aware v1.1;
2) deprecates arbitrary probability-like priority bands;
3) adds same-source-pattern reference percentiles;
4) distinguishes multi-technology-supported evidence from single-source
   evidence that actually needs TRACE-Ψ prioritization;
5) writes a compact validation table and example input;
6) performs deterministic package QA.

After Phase4E, computational method development is closed.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\RNA\Trace")
CODE = ROOT / "code"
META = ROOT / "00_meta"
TRACEPSI = ROOT / "07_tracepsi"
TOOL = ROOT / "08_tool" / "TRACEpsi"
LOGS = ROOT / "logs"

for p in [META, TOOL, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

S4C = META / "phase4c_summary.json"
S4D = META / "phase4d_summary.json"
MODEL_CARD = TOOL / "MODEL_CARD.json"
MODEL = TOOL / "TRACEpsi_frozen_model.joblib"
FP = TOOL / "TRACEpsi_fingerprint_HEK.tsv"
REF = TOOL / "TRACEpsi_reference_scores.tsv.gz"
NEW_SCORER = CODE / "23_tracepsi_score_v1_1.py"

for p in [S4C, S4D, MODEL_CARD, MODEL, FP, REF, NEW_SCORER]:
    if not p.exists():
        raise RuntimeError(f"Missing required asset: {p}")


def main():
    print("=" * 100)
    print("TRACE-Ψ PHASE 4E — RELEASE HARDENING")
    print("=" * 100)

    with open(S4C, encoding="utf-8") as f:
        s4c = json.load(f)
    with open(S4D, encoding="utf-8") as f:
        s4d = json.load(f)
    with open(MODEL_CARD, encoding="utf-8") as f:
        card = json.load(f)

    if s4c.get("status") != "TRACEPSI_EXTERNALLY_VALIDATED_A549":
        raise RuntimeError("A549 external validation status is not frozen-valid.")
    if s4d.get("status") != "TRACEPSI_V1_FROZEN_AND_PACKAGED":
        raise RuntimeError("Phase4D package is not frozen.")

    # ------------------------------------------------------------------
    # Install domain-aware scorer v1.1. Preserve the original v1 scorer.
    # ------------------------------------------------------------------
    old = TOOL / "tracepsi_score.py"
    if old.exists() and not (TOOL / "tracepsi_score_v1_0.py").exists():
        shutil.copy2(old, TOOL / "tracepsi_score_v1_0.py")

    shutil.copy2(NEW_SCORER, TOOL / "tracepsi_score.py")

    # ------------------------------------------------------------------
    # Update model card without changing scientific model.
    # ------------------------------------------------------------------
    card["version"] = "1.1-frozen"
    card["release_policy"] = {
        "model_changed_from_v1.0": False,
        "weights_or_features_changed": False,
        "purpose": "software/domain hardening only",
    }
    card["validated_input_domain"] = {
        "primary": "locally-single-U 5-mer with central U and -1/+1 != U",
        "out_of_domain_policy": (
            "retain row but return no TRACE-Ψ score for consecutive-U, "
            "non-central-U or malformed motifs"
        ),
    }
    card["recommended_interpretation"] = {
        "source_breadth_ge_2": (
            "label as multi-technology-supported evidence; chemistry breadth "
            "already carries strong portability evidence"
        ),
        "source_breadth_eq_1": (
            "main TRACE-Ψ use case; rank by score/reference percentile to "
            "prioritize orthogonal validation"
        ),
    }
    card["percentile_policy"] = {
        "primary_reference": (
            "same source-pattern reference distribution when >=20 reference loci"
        ),
        "fallback": "same source-breadth reference distribution",
        "not_a_probability": True,
    }
    card["priority_band_policy"] = {
        "v1_fixed_bands": "deprecated",
        "reason": (
            "0.90/0.75/0.25 cutoffs were heuristic and could look probability-like; "
            "v1.1 reports continuous percentiles/ranks instead"
        ),
        "recommended_use": "choose top K according to experimental budget",
    }
    card["claim_hierarchy"] = {
        "supported": [
            "TRACE-Ψ is an externally validated evidence-portability ranker.",
            "High-scoring A549 BID calls were enriched for independent DRS re-observation.",
            "External signal persisted after excluding canonical PUS7-like motifs.",
        ],
        "qualified": [
            "The full P2 model did not outperform context-only in A549; "
            "do not claim externally validated fingerprint-specific gain.",
            "Broader source/platform generalization beyond BID->DRS remains less established.",
        ],
    }

    with open(MODEL_CARD, "w", encoding="utf-8") as f:
        json.dump(card, f, ensure_ascii=False, indent=2)

    # ------------------------------------------------------------------
    # Compact validation summary
    # ------------------------------------------------------------------
    p2 = next(
        r for r in s4c["lockbox_metrics"]
        if r["score"] == "P2_TRACEpsi_frozen"
    )
    context = next(
        r for r in s4c["lockbox_metrics"]
        if r["score"] == "A_context_only_frozen"
    )
    nonpus = next(
        r for r in s4d["external_motif_robustness"]
        if r["subset"] == "A549_excluding_PUS7_like"
    )

    val = pd.DataFrame([
        {
            "benchmark": "HeLa development DRS breadth=1",
            "n": 1103,
            "positives": 159,
            "AUROC": card["development"]["breadth1_DRS_AUROC"],
            "AUROC_CI_low": card["development"]["breadth1_DRS_AUROC_CI"][0],
            "AUROC_CI_high": card["development"]["breadth1_DRS_AUROC_CI"][1],
            "note": "development benchmark",
        },
        {
            "benchmark": "A549 BID->DRS external lockbox",
            "n": int(p2["n"]),
            "positives": int(p2["positives"]),
            "AUROC": float(p2["AUROC"]),
            "AUROC_CI_low": float(p2["AUROC_CI_low"]),
            "AUROC_CI_high": float(p2["AUROC_CI_high"]),
            "note": "frozen external validation",
        },
        {
            "benchmark": "A549 external excluding PUS7-like motifs",
            "n": int(nonpus["n"]),
            "positives": int(nonpus["DRS_positive_n"]),
            "AUROC": float(nonpus["AUROC"]),
            "AUROC_CI_low": float(nonpus["CI_low"]),
            "AUROC_CI_high": float(nonpus["CI_high"]),
            "note": "post-lockbox no-retraining robustness",
        },
        {
            "benchmark": "A549 context-only ablation",
            "n": int(context["n"]),
            "positives": int(context["positives"]),
            "AUROC": float(context["AUROC"]),
            "AUROC_CI_low": float(context["AUROC_CI_low"]),
            "AUROC_CI_high": float(context["AUROC_CI_high"]),
            "note": "secondary ablation; slightly higher AUC than P2",
        },
    ])
    val.to_csv(
        TOOL / "VALIDATION_SUMMARY.tsv",
        sep="\t",
        index=False,
    )

    # ------------------------------------------------------------------
    # Example input. Use motifs observed in the frozen reference.
    # ------------------------------------------------------------------
    ref = pd.read_csv(REF, sep="\t", compression="gzip", low_memory=False)

    # The compact reference file does not contain motifs; provide safe,
    # syntactically valid illustrative motifs instead. They are examples only.
    example = pd.DataFrame([
        {"site_id": "example_1", "motif": "AGUAC", "source_assay": "BID"},
        {"site_id": "example_2", "motif": "GAUCG", "source_assay": "BACS"},
        {"site_id": "example_3", "motif": "ACUGA", "source_assay": "ELAP"},
        # Consecutive-U example intentionally demonstrates domain guard.
        {"site_id": "example_out_of_domain", "motif": "AUUAG", "source_assay": "BID"},
    ])
    example_path = TOOL / "EXAMPLE_INPUT.tsv"
    example.to_csv(example_path, sep="\t", index=False)

    # ------------------------------------------------------------------
    # Deterministic self-test.
    # ------------------------------------------------------------------
    out_path = TOOL / "EXAMPLE_OUTPUT.tsv"
    cmd = [
        sys.executable,
        str(TOOL / "tracepsi_score.py"),
        str(example_path),
        str(out_path),
        "--bundle-dir",
        str(TOOL),
    ]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "TRACE-Ψ v1.1 self-test failed.\nSTDOUT:\n"
            + proc.stdout
            + "\nSTDERR:\n"
            + proc.stderr
        )

    exout = pd.read_csv(out_path, sep="\t")
    if len(exout) != 4:
        raise RuntimeError("Self-test output row count mismatch.")

    valid = exout["TRACEpsi_domain_status"].eq("in_validated_domain")
    if int(valid.sum()) != 3:
        raise RuntimeError(
            "Self-test expected 3 in-domain and 1 out-of-domain rows."
        )
    if exout.loc[valid, "TRACEpsi_score"].isna().any():
        raise RuntimeError("Self-test in-domain rows were not scored.")
    if exout.loc[~valid, "TRACEpsi_score"].notna().any():
        raise RuntimeError("Self-test out-of-domain row received a score.")

    # ------------------------------------------------------------------
    # README v1.1
    # ------------------------------------------------------------------
    readme = """TRACE-Ψ v1.1
=============

Scientific purpose
------------------
TRACE-Ψ ranks already-reported pseudouridine evidence by its relative
cross-platform portability: how strongly a locus should be prioritized for
orthogonal validation by another profiling technology.

TRACE-Ψ is NOT:
- P(true Ψ)
- a biological truth score
- a causal functional score
- a replacement for experimental validation

Validated domain
----------------
Primary scoring is restricted to locally-single-U 5-mers:
central U, with immediate -1/+1 bases not U.

Consecutive-U or malformed motifs are retained but not scored.

Recommended use
---------------
1. If source_breadth >= 2:
   treat the locus as multi-technology-supported evidence.
   TRACE-Ψ ranking is supplementary.

2. If source_breadth == 1:
   this is the main intended use.
   rank candidates by TRACEpsi_reference_percentile or TRACEpsi_score and
   select top K according to the available orthogonal-validation budget.

Run
---
python tracepsi_score.py input.tsv output.tsv

Input
-----
TSV containing `motif` plus either:
- source_assay (BID/BACS/ELAP)
or:
- source_BID, source_BACS, source_ELAP (0/1)

External validation
-------------------
Frozen P2 model:
- HeLa DRS breadth=1 development AUROC ~0.644
- A549 BID -> native DRS external AUROC ~0.763
- external signal persisted after excluding canonical PUS7-like motifs

Important qualification
-----------------------
In A549, a generic local-sequence context ablation achieved a similar AUROC.
Therefore TRACE-Ψ should be presented as a technology-aware portability
framework, not as proof that the learned assay-fingerprint component is
universally superior to generic local sequence.

See MODEL_CARD.json and VALIDATION_SUMMARY.tsv.
"""
    (TOOL / "README.txt").write_text(readme, encoding="utf-8")

    summary = {
        "status": "TRACEPSI_V1_1_RELEASE_HARDENED",
        "model_changed": False,
        "scientific_result_changed": False,
        "improvements": [
            "validated-domain guard for locally-single-U motifs",
            "same-source-pattern reference percentile",
            "explicit evidence layer: multi-technology-supported vs single-source-ranked",
            "deprecated arbitrary fixed priority bands",
            "example input/output and deterministic self-test",
            "compact validation summary",
            "model card claim hierarchy tightened",
        ],
        "tool_dir": str(TOOL),
        "final_method_status": (
            "Model development closed. TRACE-Ψ v1.1 is the release candidate "
            "for manuscript/tool dissemination."
        ),
        "next_step": (
            "No further model fitting. Build manuscript figures/results and, if desired, "
            "a lightweight UI around the same frozen scorer."
        ),
    }
    (META / "phase4e_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print()
    print("STATUS: TRACEPSI_V1_1_RELEASE_HARDENED")
    print("Model changed: NO")
    print("External validation changed: NO")
    print("Self-test: PASS")
    print("Tool directory:", TOOL)
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase4e_summary.json")
    print(r"  D:\RNA\Trace\08_tool\TRACEpsi\MODEL_CARD.json")
    print(r"  D:\RNA\Trace\08_tool\TRACEpsi\VALIDATION_SUMMARY.tsv")
    print(r"  D:\RNA\Trace\08_tool\TRACEpsi\EXAMPLE_OUTPUT.tsv")
    return 0


if __name__ == "__main__":
    try:
        main()
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase4e_FATAL_ERROR.txt").write_text(
            err, encoding="utf-8"
        )
        print("\nFATAL ERROR\n")
        print(err)
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
