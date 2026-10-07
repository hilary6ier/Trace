#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Trace Phase 0A — local bootstrap and asset audit

Project root:
    D:\RNA\Trace

This script is deliberately LOCAL-ONLY:
- no internet access
- no third-party Python packages
- no FASTQ/BAM work
- never deletes or moves old files
- never overwrites a conflicting destination
- missing assets are reported, not treated as a script failure
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import traceback
from pathlib import Path
from typing import Dict, List

ROOT = Path(r"D:\RNA\Trace")
SEARCH_ROOT = Path(r"D:\RNA")

DIRS = [
    "00_meta",
    "00_meta/inherited_stage0",
    "01_source/BID",
    "01_source/ELAP",
    "01_source/BACS",
    "01_source/PRAISE",
    "01_source/DRS_2026",
    "01_source/code_snapshots",
    "02_reference",
    "03_harmonized",
    "04_aim1",
    "05_aim2",
    "06_aim3",
    "07_figures",
    "08_reports",
    "code",
    "logs",
]

PROJECT_LOCK = r"""# Trace project lock

Root: D:\RNA\Trace

Scientific story:
Calibration -> Measurement -> Evidence -> Inference

Frozen rules:
1. Processed/source-data-first. No FASTQ/BAM reconstruction for the current estimand.
2. A missing call in another assay means "not reported under that published pipeline",
   never "unmodified" or a biological negative.
3. Do not estimate a hidden "true Psi probability" without a valid denominator.
4. Preserve both exact base-level locus identity and chemistry-aware U-run equivalence.
5. For U-runs >1, do not assign arbitrary positional sequence features around an
   uncertain reported nucleotide; analyse them as a localization-bias stratum.
6. Aim 1 separates stable assay-associated effects from biological-context effects.
7. Synthetic assay calibration is a mechanistic anchor; ML is a measurement tool.
8. External cell/platform test sets are never used for feature/hyperparameter tuning.
9. DRS is independent-platform support, not ground truth.
10. Primary thresholds/features/metrics are frozen before external-test inspection.
11. Null primary results are retained; no post-hoc subgroup fishing.
12. Stage0/phase0 files are provenance archives: copy only, never move/delete.
"""

ASSETS = [
    {
        "asset_id": "PRAISE_supplementary",
        "assay": "PRAISE",
        "role": "primary+writer_validation",
        "filenames": ["41589_2023_1304_MOESM3_ESM.xlsx"],
        "dest_rel": r"01_source\PRAISE\41589_2023_1304_MOESM3_ESM.xlsx",
        "expected_sha256": "283c24e702c358b749ebbde024fb9c4edc7a9b134ca198678f170330d4bafcb5",
        "required_now": True,
    },
    {
        "asset_id": "BID_HEK293T_WT",
        "assay": "BID-seq",
        "role": "Aim1_primary",
        "filenames": ["GSE179798_HEK293T_mRNA_WT_BID-seq.xlsx"],
        "dest_rel": r"01_source\BID\GSE179798_HEK293T_mRNA_WT_BID-seq.xlsx",
        "expected_sha256": "47d85533487ce925de014a9ecbda13624925c1204db5b8c453884cf2594de6d1",
        "required_now": True,
    },
    {
        "asset_id": "BID_HeLa_WT",
        "assay": "BID-seq",
        "role": "Aim1_primary",
        "filenames": ["GSE179798_HeLa_mRNA_WT_BID-seq.xlsx"],
        "dest_rel": r"01_source\BID\GSE179798_HeLa_mRNA_WT_BID-seq.xlsx",
        "expected_sha256": "fc5eb85336811b9000d3f758403eb10df7688c311abe6176f3d8264f27175b94",
        "required_now": True,
    },
    {
        "asset_id": "BID_HeLa_shControl",
        "assay": "BID-seq",
        "role": "Aim2_writer_validation",
        "filenames": ["GSE179798_HeLa_mRNA_shControl_BID-seq.xlsx"],
        "dest_rel": r"01_source\BID\GSE179798_HeLa_mRNA_shControl_BID-seq.xlsx",
        "expected_sha256": "17a2bb53cfc517cd6200adcab203dbfafa94bd6d26459ee6907b579cc497f763",
        "required_now": False,
    },
    {
        "asset_id": "BID_A549_WT",
        "assay": "BID-seq",
        "role": "Aim1_replication",
        "filenames": ["GSE179798_A549_mRNA_WT_BID-seq.xlsx"],
        "dest_rel": r"01_source\BID\GSE179798_A549_mRNA_WT_BID-seq.xlsx",
        "expected_sha256": "",
        "required_now": False,
    },
    {
        "asset_id": "ELAP_HEK293T",
        "assay": "ELAP-seq",
        "role": "Aim1_primary",
        "filenames": ["GSE236530_ELAP-HEK-all-update-1.xlsx"],
        "dest_rel": r"01_source\ELAP\GSE236530_ELAP-HEK-all-update-1.xlsx",
        "expected_sha256": "be52b9c12c1bc7e79b9c0fbb2e37103bec5e90fb41754ca3362abc1c9b48fd19",
        "required_now": True,
    },
    {
        "asset_id": "ELAP_HeLa",
        "assay": "ELAP-seq",
        "role": "Aim1_primary",
        "filenames": ["GSE236530_ELAP-HeLa-all-update-1.xlsx"],
        "dest_rel": r"01_source\ELAP\GSE236530_ELAP-HeLa-all-update-1.xlsx",
        "expected_sha256": "730a4bc6dcb7ec06df32975ce7722085ea4e5a418e37c5f91d7838a252aa42b8",
        "required_now": True,
    },
    {
        "asset_id": "ELAP_HEK293T_DKC1",
        "assay": "ELAP-seq",
        "role": "Aim2_writer_validation",
        "filenames": ["GSE236530_DKC1-knockdown-enrichment-change.xlsx"],
        "dest_rel": r"01_source\ELAP\GSE236530_DKC1-knockdown-enrichment-change.xlsx",
        "expected_sha256": "b451609a7ffa0893da64993513169818ba30b65c8d0f424809db90b64020f673",
        "required_now": False,
    },
    {
        "asset_id": "ELAP_supplementary_data",
        "assay": "ELAP-seq",
        "role": "paper_supplement",
        "filenames": ["ELAP_supplementary_data.xlsx", "41467_2026_70597_MOESM3_ESM.xlsx"],
        "dest_rel": r"01_source\ELAP\ELAP_supplementary_data.xlsx",
        "expected_sha256": "c799d971d8bfab758d6028706e03ab6eb5dabc6f4a3f10cfbe63465d029730fc",
        "required_now": True,
    },
    {
        "asset_id": "ELAP_source_data",
        "assay": "ELAP-seq",
        "role": "synthetic_calibration",
        "filenames": ["ELAP_source_data.xlsx", "41467_2026_70597_MOESM6_ESM.xlsx"],
        "dest_rel": r"01_source\ELAP\ELAP_source_data.xlsx",
        "expected_sha256": "a956b8ae54064b8f9b4107c45c814551a2dc13aa808d9f5e229d0634ff14f33e",
        "required_now": True,
    },
    {
        "asset_id": "BACS_HeLa",
        "assay": "BACS",
        "role": "Aim2_primary+localization",
        "filenames": ["GSE241849_Supplementary_Table_3_10.xlsx"],
        "dest_rel": r"01_source\BACS\GSE241849_Supplementary_Table_3_10.xlsx",
        "expected_sha256": "eed449b5481c88a16b0e52ee2bd13caa56a190fd05d69458450bf7128ae6786a",
        "required_now": True,
    },
    {
        "asset_id": "ELAP_code_snapshot",
        "assay": "ELAP-seq",
        "role": "provenance_code",
        "filenames": ["ELAP-seq-v4.zip"],
        "dest_rel": r"01_source\code_snapshots\ELAP-seq-v4.zip",
        "expected_sha256": "9a34b3e8ab644848af002a50c70bb470e2098e7c8905f687491c4e4848be0f78",
        "required_now": False,
    },
    {
        "asset_id": "BACS_code_snapshot",
        "assay": "BACS",
        "role": "provenance_code",
        "filenames": ["bacs.zip"],
        "dest_rel": r"01_source\code_snapshots\bacs.zip",
        "expected_sha256": "9ffd41175b219994b4d6cec4d6d03e67d6b9ff6c4cb15e29053ee6053e9b464f",
        "required_now": False,
    },
]

REPORT_FILENAMES = [
    "01_existing_asset_audit.txt",
    "01_existing_asset_audit.json",
    "02_external_processed_manifest.json",
    "04_processed_evidence_gate.txt",
    "05_coordinate_convention_audit.txt",
    "05_coordinate_audit.json",
    "05b_locus_semantics_audit.txt",
    "05b_locus_semantics_summary.json",
    "06a_recover_published_comparison_evidence.txt",
]

REF_SUFFIXES = (
    ".fa", ".fasta", ".fa.gz", ".fasta.gz", ".gtf", ".gtf.gz",
    ".gff3", ".gff3.gz", ".fai"
)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def write_tsv(path: Path, rows: List[Dict], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def index_files(root: Path):
    idx = {}
    if not root.exists():
        return idx
    trace_resolved = ROOT.resolve()
    for p in root.rglob("*"):
        try:
            if not p.is_file():
                continue
            try:
                if trace_resolved in p.resolve().parents:
                    continue
            except Exception:
                pass
            idx.setdefault(p.name.lower(), []).append(p)
        except (PermissionError, OSError):
            continue
    return idx

def choose_source(asset, file_index):
    candidates = []
    for name in asset["filenames"]:
        candidates.extend(file_index.get(name.lower(), []))
    if not candidates:
        return None

    expected = asset.get("expected_sha256", "")
    if expected:
        for p in candidates:
            try:
                if sha256(p) == expected:
                    return p
            except OSError:
                pass
    return candidates[0]

def copy_asset(asset, src: Path, rows: List[Dict]):
    dst = ROOT / asset["dest_rel"]
    dst.parent.mkdir(parents=True, exist_ok=True)

    src_hash = sha256(src)
    expected = asset.get("expected_sha256", "")
    hash_ok = (not expected) or (src_hash == expected)

    if dst.exists():
        dst_hash = sha256(dst)
        if dst_hash == src_hash:
            status = "PRESENT_VERIFIED"
        else:
            rows.append({
                "asset_id": asset["asset_id"],
                "assay": asset["assay"],
                "role": asset["role"],
                "required_now": asset["required_now"],
                "status": "DESTINATION_CONFLICT",
                "source_path": str(src),
                "trace_path": str(dst),
                "size_bytes": dst.stat().st_size,
                "sha256": dst_hash,
                "expected_sha256": expected,
                "hash_match_expected": False,
                "notes": "Destination differs from discovered source; left untouched.",
            })
            return
    else:
        shutil.copy2(src, dst)
        status = "COPIED_VERIFIED" if hash_ok else "COPIED_HASH_UNCONFIRMED"

    rows.append({
        "asset_id": asset["asset_id"],
        "assay": asset["assay"],
        "role": asset["role"],
        "required_now": asset["required_now"],
        "status": status,
        "source_path": str(src),
        "trace_path": str(dst),
        "size_bytes": dst.stat().st_size,
        "sha256": sha256(dst),
        "expected_sha256": expected,
        "hash_match_expected": hash_ok,
        "notes": "" if hash_ok else "File found by name but did not match the previously recorded hash.",
    })

def main():
    print("=" * 78)
    print("TRACE PHASE 0A — LOCAL BOOTSTRAP")
    print(r"Project root: D:\RNA\Trace")
    print("=" * 78)

    for rel in DIRS:
        (ROOT / rel).mkdir(parents=True, exist_ok=True)

    (ROOT / "00_meta" / "PROJECT_LOCK.md").write_text(PROJECT_LOCK, encoding="utf-8")

    print(r"[1/4] Indexing local files under D:\RNA ...")
    file_index = index_files(SEARCH_ROOT)

    print("[2/4] Discovering and copying known research assets ...")
    rows = []
    for asset in ASSETS:
        src = choose_source(asset, file_index)
        if src is None:
            rows.append({
                "asset_id": asset["asset_id"],
                "assay": asset["assay"],
                "role": asset["role"],
                "required_now": asset["required_now"],
                "status": "NOT_FOUND_LOCAL",
                "source_path": "",
                "trace_path": str(ROOT / asset["dest_rel"]),
                "size_bytes": "",
                "sha256": "",
                "expected_sha256": asset.get("expected_sha256", ""),
                "hash_match_expected": "",
                "notes": "Not a failure. If needed, this will be acquired in the next download step.",
            })
            print(f"  - {asset['asset_id']}: not found locally")
            continue
        copy_asset(asset, src, rows)
        print(f"  + {asset['asset_id']}: {src}")

    print("[3/4] Inheriting Stage0 reports and scanning reference candidates ...")
    inherited_dir = ROOT / "00_meta" / "inherited_stage0"
    inherited_dir.mkdir(parents=True, exist_ok=True)

    report_rows = []
    for name in REPORT_FILENAMES:
        matches = file_index.get(name.lower(), [])
        if not matches:
            report_rows.append({"filename": name, "status": "NOT_FOUND", "source_path": "", "trace_path": ""})
            continue
        src = matches[0]
        dst = inherited_dir / name
        if not dst.exists():
            shutil.copy2(src, dst)
        report_rows.append({"filename": name, "status": "COPIED_OR_PRESENT", "source_path": str(src), "trace_path": str(dst)})

    ref_rows = []
    seen = set()
    for paths in file_index.values():
        for p in paths:
            n = p.name.lower()
            if not n.endswith(REF_SUFFIXES):
                continue
            if not any(k in n for k in ("hg38", "grch38", "gencode", "mane")):
                continue
            key = str(p).lower()
            if key in seen:
                continue
            seen.add(key)
            try:
                size = p.stat().st_size
            except OSError:
                size = ""
            ref_rows.append({
                "path": str(p),
                "filename": p.name,
                "size_bytes": size,
                "reference_hint": (
                    "MANE" if "mane" in n else
                    "GENCODE" if "gencode" in n else
                    "GRCh38/hg38"
                ),
            })

    print("[4/4] Writing manifests ...")
    asset_fields = [
        "asset_id", "assay", "role", "required_now", "status",
        "source_path", "trace_path", "size_bytes", "sha256",
        "expected_sha256", "hash_match_expected", "notes",
    ]
    write_tsv(ROOT / "00_meta" / "local_asset_inventory.tsv", rows, asset_fields)
    write_tsv(ROOT / "00_meta" / "inherited_reports.tsv", report_rows,
              ["filename", "status", "source_path", "trace_path"])
    write_tsv(ROOT / "00_meta" / "reference_candidates.tsv", ref_rows,
              ["path", "filename", "size_bytes", "reference_hint"])

    missing = [r for r in rows if r["status"] == "NOT_FOUND_LOCAL"]
    required_missing = [r for r in missing if str(r["required_now"]).lower() == "true"]
    write_tsv(ROOT / "00_meta" / "missing_assets.tsv", missing, asset_fields)

    summary = {
        "project_root": str(ROOT),
        "n_assets_expected": len(rows),
        "n_assets_found_or_copied": len(rows) - len(missing),
        "n_assets_not_found_local": len(missing),
        "n_required_now_missing": len(required_missing),
        "n_reference_candidates": len(ref_rows),
        "required_missing_asset_ids": [r["asset_id"] for r in required_missing],
        "all_missing_asset_ids": [r["asset_id"] for r in missing],
        "status": "LOCAL_AUDIT_COMPLETE",
    }
    (ROOT / "00_meta" / "phase0a_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print()
    print("=" * 78)
    print("LOCAL AUDIT COMPLETE")
    print("=" * 78)
    print(f"Found/copied:          {summary['n_assets_found_or_copied']}/{summary['n_assets_expected']}")
    print(f"Not found locally:     {summary['n_assets_not_found_local']}")
    print(f"Required-now missing:  {summary['n_required_now_missing']}")
    print(f"Reference candidates:  {summary['n_reference_candidates']}")
    print()
    print(r"Outputs: D:\RNA\Trace\00_meta")
    print("NOT_FOUND_LOCAL does not mean failure.")
    return 0

if __name__ == "__main__":
    try:
        main()
    except Exception:
        err = traceback.format_exc()
        try:
            ROOT.mkdir(parents=True, exist_ok=True)
            (ROOT / "logs").mkdir(parents=True, exist_ok=True)
            (ROOT / "logs" / "phase0a_FATAL_ERROR.txt").write_text(err, encoding="utf-8")
        except Exception:
            pass
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase0a_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
