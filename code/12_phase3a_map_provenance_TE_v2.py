#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRACE Phase 3A v2 — map provenance and translation-efficiency inference
===================================================================

Scientific question
-------------------
On the SAME HeLa functional outcomes and the SAME gene cohort, does the
association between pseudouridine burden and translation efficiency change
when the underlying published Ψ map changes?

This is the first Inference-stage analysis:
    Calibration -> Measurement -> Evidence -> Inference

Primary maps
------------
- BID HeLa
- BACS HeLa
- ELAP HeLa
- DRS HeLa (author psi=1 calls)

Primary outcome
---------------
HeLa translation efficiency values from NAR 2026 Supplementary Table S3.
The paper compiles multiple public HeLa Ribo-seq/RNA-seq datasets; these are
treated as repeated outcome definitions, not pooled on their raw scales.

For each HeLa TE dataset:
    y = z-score(log10(TE)), TE>0

Primary map-burden model
------------------------
Stack all resolved HeLa TE datasets after within-dataset outcome
standardization:

    TE_z ~ standardized log1p(map-specific Ψ count)
           + standardized log(MANE transcript length)
           + TE-dataset fixed effects

SEs are clustered by gene.

Every map model uses the exact same stacked gene-dataset rows.

Primary question is NOT "which map is best". We report:
- standardized burden effect for each map;
- gene-cluster bootstrap pairwise differences between map effects.

Technology-aware secondary
--------------------------
Using the already frozen Phase2B locally-single-U HeLa source union:
- naive_union_count: all BID/BACS/ELAP source-positive single-U loci
- orthogonal_count: loci supported by >=2 source assays
- triple_count: loci supported by all 3 assays
- breadth_sum: sum of BID/BACS/ELAP support counts across loci

These are secondary bridge analyses only. The final sequence-independent test
(M0/M1/M2) is deferred to Phase3B.

Interpretation boundaries
-------------------------
- A gene with count=0 means "no Ψ site assigned to this MANE gene in this
  published map", not "unmodified".
- This is association/inference sensitivity, not a causal Ψ effect.
- DRS is one map, not ground truth.
- No raw sequencing.
- No new threshold tuning.
"""

from __future__ import annotations

import gzip
import json
import math
import re
import traceback
from urllib.parse import unquote
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
REF = ROOT / "02_reference"
NORM = ROOT / "03_harmonized" / "normalized_sources"
AIM2 = ROOT / "05_aim2"
AIM3 = ROOT / "06_aim3"
LOGS = ROOT / "logs"

for p in [META, AIM3, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

GTF = REF / "MANE.GRCh38.v1.5.refseq_genomic.gtf.gz"

BID = NORM / "BID_HeLa.tsv"
BACS = NORM / "BACS_HeLa.tsv"
ELAP = NORM / "ELAP_HeLa.tsv"
DRS_S6 = NORM / "DRS_S6_psi_calls.tsv"

# Phase0C normally created this name. A safe discovery fallback is implemented.
S3_EXPECTED = NORM / "DRS_S3_functional.tsv"

PHASE2B_UNION = AIM2 / "aim2b_HeLa_source_union_DRS.tsv"

SEED = 20261003
RNG = np.random.default_rng(SEED)
N_BOOT = 1500

HELA_TE_STUDIES = [
    "GSE21992",
    "GSE79664",
    "GSE143301",
    "GSE188692",
    "SRA099816",
]


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------

def ss(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    return str(x).strip()


def norm_chr(x: Any) -> str:
    """
    Harmonize public-site chromosome names and NCBI RefSeq genomic accessions.

    MANE RefSeq genomic annotations use accessions such as NC_000001.11,
    whereas the Ψ maps use chr1/chr2/... . Primary GRCh38 accessions are
    converted explicitly; patch accessions (NW_/NT_) are preserved.
    """
    z = ss(x)
    if not z:
        return ""

    if z.lower().startswith("chr"):
        tail = z[3:]
        return "chr" + ("M" if tail.upper() == "MT" else tail)

    if re.fullmatch(r"(?:\d+|X|Y|M|MT)", z, flags=re.I):
        return "chr" + ("M" if z.upper() == "MT" else z.upper())

    m = re.match(r"^NC_(\d{6})(?:\.\d+)?$", z, flags=re.I)
    if m:
        n = int(m.group(1))
        if 1 <= n <= 22:
            return f"chr{n}"
        if n == 23:
            return "chrX"
        if n == 24:
            return "chrY"
        if n == 12920:  # NC_012920.1 mitochondrial genome
            return "chrM"

    # MANE v1.5 includes a small number of genes on GRCh38 patch sequences.
    # Preserve those accessions; they are valid MANE records but will not be
    # forced onto a primary chromosome.
    return z


def parse_attrs(text: str) -> Dict[str, str]:
    """
    Parse either GTF-style:
        gene "TP53"; transcript_id "NM_...";
    or GFF3/NCBI-style:
        gene=TP53;transcript_id=NM_...;Dbxref=GeneID:7157,...

    Repeated keys are retained by comma-concatenation.
    """
    raw = ss(text)
    out: Dict[str, List[str]] = defaultdict(list)

    # GFF3-like key=value syntax.
    if "=" in raw:
        for field in raw.split(";"):
            field = field.strip()
            if not field or "=" not in field:
                continue
            k, v = field.split("=", 1)
            k = k.strip()
            v = unquote(v.strip().strip('"'))
            if k:
                out[k].append(v)

    # GTF quoted syntax. Run even for mixed files.
    for k, v in re.findall(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"', raw):
        if v not in out[k]:
            out[k].append(v)

    # Conservative unquoted GTF fallback: key value;
    if not out:
        for field in raw.split(";"):
            field = field.strip()
            m = re.match(r"^([A-Za-z0-9_.:-]+)\s+(.+)$", field)
            if m:
                k, v = m.group(1), m.group(2).strip().strip('"')
                out[k].append(unquote(v))

    return {k: ",".join(vs) for k, vs in out.items()}


def attr_first(a: Dict[str, str], keys: Sequence[str]) -> str:
    for k in keys:
        v = ss(a.get(k, ""))
        if v:
            return v.split(",")[0].strip()
    return ""


def dbxref_values(a: Dict[str, str]) -> Dict[str, List[str]]:
    out = defaultdict(list)
    raw = ss(a.get("Dbxref", "") or a.get("db_xref", ""))
    for token in raw.split(","):
        token = token.strip()
        if ":" not in token:
            continue
        k, v = token.split(":", 1)
        if k and v:
            out[k].append(v)
    return dict(out)


def mane_gene_symbol(a: Dict[str, str]) -> str:
    # NCBI RefSeq GTF/GFF commonly uses `gene`, not `gene_name`.
    return attr_first(a, ["gene_name", "gene", "gene_symbol", "Name"])


def mane_gene_id(a: Dict[str, str]) -> str:
    v = attr_first(a, ["gene_id", "GeneID"])
    if v:
        return re.sub(r"^gene-", "", v)
    x = dbxref_values(a)
    if x.get("GeneID"):
        return x["GeneID"][0]
    parent = attr_first(a, ["gene_parent"])
    return re.sub(r"^gene-", "", parent)


def mane_transcript_id(a: Dict[str, str]) -> str:
    v = attr_first(a, ["transcript_id", "transcript"])
    if v:
        return re.sub(r"^rna-", "", v)

    name = attr_first(a, ["Name"])
    if re.match(r"^(?:N[MR]|X[MR]|ENST)_?\d|^ENST\d", name, flags=re.I):
        return re.sub(r"^rna-", "", name)

    parent = attr_first(a, ["Parent"])
    if parent:
        for token in parent.split(","):
            token = token.strip()
            if token.startswith("rna-"):
                return token[4:]
            if re.match(r"^(?:NM_|NR_|XM_|XR_|ENST)", token):
                return token

    ident = attr_first(a, ["ID"])
    if ident.startswith("rna-"):
        return ident[4:]
    return ""


def mane_tag(a: Dict[str, str]) -> str:
    return ss(a.get("tag", "")).replace("_", " ").strip()


def strip_version(x: Any) -> str:
    return re.sub(r"\.\d+$", "", ss(x))


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    mu = x.mean()
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index)
    return (x-mu)/sd


def bh_adjust(pvals: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    if len(vals) == 0:
        return out
    order = np.argsort(vals)
    ranked = vals[order]
    q = np.empty(len(vals))
    prev = 1.0
    n = len(vals)
    for i in range(n-1, -1, -1):
        rank = i+1
        v = min(prev, ranked[i]*n/rank)
        q[order[i]] = v
        prev = v
    out[np.where(ok)[0]] = np.minimum(q, 1.0)
    return out


def resolve_s3_path() -> Path:
    if S3_EXPECTED.exists():
        return S3_EXPECTED

    cands = [
        p for p in NORM.glob("*.tsv")
        if "s3" in p.name.lower()
        and any(t in p.name.lower() for t in ["drs", "functional"])
    ]
    if len(cands) == 1:
        return cands[0]

    raise RuntimeError(
        "Could not uniquely resolve normalized DRS Supplementary Table S3. "
        f"Expected {S3_EXPECTED}; candidates={[p.name for p in cands]}"
    )


# -----------------------------------------------------------------------------
# MANE annotation
# -----------------------------------------------------------------------------

def load_mane():
    """
    Parse NCBI/RefSeq MANE v1.5 robustly.

    Scientific choice:
      use MANE Select exons only (one representative transcript per locus),
      rather than mixing MANE Plus Clinical alternatives into the transcript-
      length baseline.

    Returns:
      gene_table
      intervals_strand
      intervals_chrom
      gene_name_map
      gene_id_map
      tx_id_map
      qc
    """
    tx_exon_len = defaultdict(int)
    tx_meta = {}
    intervals_strand = defaultdict(list)
    intervals_chrom = defaultdict(list)

    feature_counts = defaultdict(int)
    exon_attribute_key_counts = defaultdict(int)
    chrom_raw_counts = defaultdict(int)
    chrom_norm_counts = defaultdict(int)

    exon_rows_total = 0
    exon_rows_select = 0
    exon_rows_missing_gene = 0
    exon_rows_missing_tx = 0
    exon_rows_bad_coord = 0
    sample_rows = []

    opener = gzip.open if str(GTF).endswith(".gz") else open

    with opener(GTF, "rt", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, start=1):
            if not line or line.startswith("#"):
                continue

            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                continue

            raw_chrom, source, feature, start, end, score, strand, frame, attrs = p
            feature_counts[feature] += 1

            if feature.lower() != "exon":
                continue

            exon_rows_total += 1
            a = parse_attrs(attrs)

            if len(sample_rows) < 25:
                sample_rows.append({
                    "line_no": line_no,
                    "raw_chrom": raw_chrom,
                    "normalized_chrom": norm_chr(raw_chrom),
                    "source": source,
                    "feature": feature,
                    "strand": strand,
                    "attribute_keys": "|".join(sorted(a.keys())),
                    "gene_name": ss(a.get("gene_name","")),
                    "gene": ss(a.get("gene","")),
                    "gene_id": ss(a.get("gene_id","")),
                    "transcript_id": ss(a.get("transcript_id","")),
                    "tag": ss(a.get("tag","")),
                    "Dbxref": ss(a.get("Dbxref",""))[:500],
                    "raw_attributes": attrs[:1500],
                })

            for k in a:
                exon_attribute_key_counts[k] += 1

            tag = mane_tag(a)
            # NCBI states exon rows associated with MANE carry the MANE tag.
            # Keep MANE Select; explicitly exclude Plus Clinical alternatives.
            is_select = "mane select" in tag.lower()
            if not is_select:
                continue
            exon_rows_select += 1

            gene_name = mane_gene_symbol(a)
            gene_id = mane_gene_id(a)
            tx = mane_transcript_id(a)

            if not gene_name:
                exon_rows_missing_gene += 1
                continue
            if not tx:
                exon_rows_missing_tx += 1
                continue

            try:
                st, en = int(start), int(end)
            except Exception:
                exon_rows_bad_coord += 1
                continue

            chrom = norm_chr(raw_chrom)
            chrom_raw_counts[raw_chrom] += 1
            chrom_norm_counts[chrom] += 1

            # Transcript identity is RefSeq accession; MANE Select should be
            # unique, but include gene symbol defensively against malformed IDs.
            tx_key = (gene_name, tx)

            tx_exon_len[tx_key] += en-st+1
            tx_meta[tx_key] = {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "transcript_id": tx,
                "chrom": chrom,
                "raw_chrom": raw_chrom,
                "strand": strand,
            }

            rec = (st, en, gene_name)
            intervals_strand[(chrom, strand)].append(rec)
            intervals_chrom[chrom].append(rec)

    sample_df = pd.DataFrame(sample_rows)
    sample_df.to_csv(
        META / "phase3a_v2_MANE_attribute_examples.tsv",
        sep="\t", index=False
    )

    # Pre-QC is always written, including on failure.
    pre_qc = {
        "reference_file": str(GTF),
        "feature_counts_top": dict(
            sorted(feature_counts.items(), key=lambda kv: kv[1], reverse=True)[:20]
        ),
        "exon_rows_total": int(exon_rows_total),
        "exon_rows_with_MANE_Select_tag": int(exon_rows_select),
        "exon_rows_missing_gene_after_select": int(exon_rows_missing_gene),
        "exon_rows_missing_transcript_after_select": int(exon_rows_missing_tx),
        "exon_rows_bad_coordinate_after_select": int(exon_rows_bad_coord),
        "top_exon_attribute_keys": dict(
            sorted(
                exon_attribute_key_counts.items(),
                key=lambda kv: kv[1],
                reverse=True
            )[:30]
        ),
    }
    (META / "phase3a_v2_MANE_preQC.json").write_text(
        json.dumps(pre_qc, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    if exon_rows_total == 0:
        raise RuntimeError(
            "MANE reference contains no exon feature rows. "
            r"See 00_meta\phase3a_v2_MANE_preQC.json"
        )

    if exon_rows_select < 100000:
        raise RuntimeError(
            f"Only {exon_rows_select} MANE-Select exon rows were parsed from "
            f"{exon_rows_total} total exon rows. This is inconsistent with the "
            "NCBI MANE v1.5 reference and suggests an attribute/tag parsing problem. "
            r"See 00_meta\phase3a_v2_MANE_attribute_examples.tsv and "
            r"00_meta\phase3a_v2_MANE_preQC.json"
        )

    if not tx_meta:
        raise RuntimeError(
            "MANE Select exon rows were found, but no transcript records could "
            "be constructed. See the MANE attribute audit outputs."
        )

    tx_rows = []
    for tx_key, meta in tx_meta.items():
        tx_rows.append({
            **meta,
            "transcript_length": tx_exon_len[tx_key],
        })
    txdf = pd.DataFrame(tx_rows)

    # MANE Select is intended to provide one representative transcript per locus.
    # If duplicate gene symbols occur, preserve one deterministic representative
    # only after auditing them.
    duplicate_gene_symbols = (
        txdf.groupby("gene_name")["transcript_id"].nunique()
        .loc[lambda s: s > 1]
        .sort_values(ascending=False)
    )
    if len(duplicate_gene_symbols):
        duplicate_gene_symbols.rename("n_select_transcripts").to_csv(
            META / "phase3a_v2_MANE_duplicate_gene_symbols.tsv",
            sep="\t"
        )

    txdf = txdf.sort_values(
        ["gene_name","transcript_length","transcript_id"],
        ascending=[True,False,True]
    )
    gene_table = txdf.drop_duplicates("gene_name", keep="first").copy()

    n_tx = int(len(txdf))
    n_gene = int(gene_table["gene_name"].nunique())

    # NCBI MANE v1.5 reports 19,437 MANE Select transcripts overall.
    # Use a range rather than an exact number to remain robust to harmless
    # parsing differences while still catching reference/schema failures.
    if not (19000 <= n_tx <= 20050):
        raise RuntimeError(
            f"Parsed MANE Select transcripts={n_tx}, outside expected v1.5 "
            "range [19000,20050]. "
            r"See 00_meta\phase3a_v2_MANE_preQC.json"
        )
    if not (19000 <= n_gene <= 20050):
        raise RuntimeError(
            f"Parsed MANE Select gene symbols={n_gene}, outside expected "
            "range [19000,20050]."
        )

    for key in intervals_strand:
        intervals_strand[key].sort(key=lambda z: z[0])
    for key in intervals_chrom:
        intervals_chrom[key].sort(key=lambda z: z[0])

    # Identifier maps for downstream S3 schema discovery.
    gene_name_map = {}
    gene_id_map = {}
    tx_id_map = {}

    for r in txdf.itertuples():
        gene_name_map[r.gene_name.upper()] = r.gene_name

        if r.gene_id:
            gene_id_map[strip_version(r.gene_id)] = r.gene_name

        # RefSeq transcript itself.
        tx_id_map[strip_version(r.transcript_id)] = r.gene_name

    primary_chr_tx = int(
        txdf["chrom"].astype(str).str.match(r"^chr(?:[1-9]|1\d|2[0-2]|X|Y|M)$").sum()
    )
    patch_tx = n_tx-primary_chr_tx

    qc = {
        **pre_qc,
        "MANE_Select_transcripts": n_tx,
        "MANE_Select_genes": n_gene,
        "duplicate_gene_symbols_with_multiple_Select_transcripts": int(
            len(duplicate_gene_symbols)
        ),
        "primary_chromosome_transcripts": primary_chr_tx,
        "patch_or_nonprimary_transcripts": patch_tx,
        "primary_chromosome_fraction": float(primary_chr_tx/n_tx),
        "normalized_chromosome_counts_top": dict(
            sorted(
                chrom_norm_counts.items(),
                key=lambda kv: kv[1],
                reverse=True
            )[:30]
        ),
        "median_transcript_length": float(
            gene_table["transcript_length"].median()
        ),
    }

    (META / "phase3a_v2_MANE_QC.json").write_text(
        json.dumps(qc, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    # Reference should overwhelmingly map to GRCh38 primary chromosomes;
    # MANE v1.5 legitimately retains a small patch subset.
    if qc["primary_chromosome_fraction"] < 0.99:
        raise RuntimeError(
            f"Only {qc['primary_chromosome_fraction']:.3%} of MANE Select "
            "transcripts mapped to primary chr names. Accession harmonization "
            "requires review."
        )

    return (
        gene_table,
        dict(intervals_strand),
        dict(intervals_chrom),
        gene_name_map,
        gene_id_map,
        tx_id_map,
        qc,
    )


def annotate_positions(
    q: pd.DataFrame,
    intervals: Dict,
    stranded: bool,
) -> pd.DataFrame:
    """
    q requires row_id, chrom, pos1 and strand if stranded.
    Returns unique-gene annotation status for each query row.
    """
    out_rows = []

    group_cols = ["chrom","strand"] if stranded else ["chrom"]

    for key_vals, sub in q.groupby(group_cols, sort=False, dropna=False):
        if stranded:
            key = key_vals if isinstance(key_vals, tuple) else (key_vals,)
        else:
            key = key_vals[0] if isinstance(key_vals, tuple) else key_vals

        ints = intervals.get(key, [])
        queries = sorted(
            [(int(r.pos1), int(r.row_id)) for r in sub.itertuples()],
            key=lambda z: z[0]
        )

        active = []
        j = 0
        for pos, row_id in queries:
            while j < len(ints) and ints[j][0] <= pos:
                active.append(ints[j])
                j += 1
            active = [iv for iv in active if iv[1] >= pos]
            genes = sorted({
                iv[2] for iv in active if iv[0] <= pos <= iv[1]
            })

            if len(genes) == 1:
                status = "unique_gene"
                gene = genes[0]
            elif len(genes) == 0:
                status = "no_MANE_exon"
                gene = ""
            else:
                status = "multi_gene_ambiguous"
                gene = ""

            out_rows.append({
                "row_id": row_id,
                "mane_gene": gene,
                "annotation_status": status,
                "n_overlapping_genes": len(genes),
            })

    return pd.DataFrame(out_rows)


# -----------------------------------------------------------------------------
# Published map -> MANE gene burden
# -----------------------------------------------------------------------------

def source_exact_map(path: Path, assay: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    req = {"chrom","pos1","strand"}
    if not req.issubset(df.columns):
        raise RuntimeError(f"{assay}: missing {sorted(req-set(df.columns))}")

    x = pd.DataFrame({
        "chrom": df["chrom"].map(norm_chr),
        "pos1": pd.to_numeric(df["pos1"], errors="coerce"),
        "strand": df["strand"].astype(str).str.strip(),
    })
    x = x[
        x["chrom"].ne("")
        & x["pos1"].notna()
        & x["strand"].isin(["+","-"])
    ].copy()
    x["pos1"] = x["pos1"].astype(int)
    x = x.drop_duplicates(["chrom","strand","pos1"]).reset_index(drop=True)
    x["row_id"] = np.arange(len(x))
    x["assay"] = assay
    return x


def semantic_col(df: pd.DataFrame, candidates: Sequence[str]) -> str:
    canon = {}
    for c in df.columns:
        z = c.lower()
        if z.startswith("raw__"):
            z = z[5:]
        z = re.sub(r"[^a-z0-9]+", "", z)
        canon[c] = z

    for t0 in candidates:
        t = re.sub(r"[^a-z0-9]+", "", t0.lower())
        exact = [c for c,z in canon.items() if z == t]
        if len(exact) == 1:
            return exact[0]

    for t0 in candidates:
        t = re.sub(r"[^a-z0-9]+", "", t0.lower())
        partial = [c for c,z in canon.items() if t in z]
        if len(partial) == 1:
            return partial[0]

    raise RuntimeError(
        f"Cannot resolve {list(candidates)} from {list(df.columns)}"
    )


def drs_hela_map() -> pd.DataFrame:
    df = pd.read_csv(DRS_S6, sep="\t", low_memory=False, dtype=str)
    if "source_sheet" not in df.columns:
        raise RuntimeError("DRS S6 source_sheet missing.")
    h = df[df["source_sheet"].astype(str).str.lower().eq("hela")].copy()

    chr_col = semantic_col(h, ["chr","chromosome"])
    pos_col = semantic_col(h, ["position","pos"])
    psi_col = semantic_col(h, ["psi"])

    x = pd.DataFrame({
        "chrom": h[chr_col].map(norm_chr),
        "pos1": pd.to_numeric(h[pos_col], errors="coerce"),
        "psi": pd.to_numeric(h[psi_col], errors="coerce"),
    })
    x = x[
        x["chrom"].ne("")
        & x["pos1"].notna()
        & x["psi"].eq(1)
    ].copy()
    x["pos1"] = x["pos1"].astype(int)
    x = x.drop_duplicates(["chrom","pos1"]).reset_index(drop=True)
    x["row_id"] = np.arange(len(x))
    x["assay"] = "DRS"
    return x


def annotate_map_to_genes(
    sites: pd.DataFrame,
    intervals_strand: Dict,
    intervals_chrom: Dict,
    stranded: bool,
) -> Tuple[pd.DataFrame, Dict]:
    ann = annotate_positions(
        sites[["row_id","chrom","pos1"] + (["strand"] if stranded else [])],
        intervals_strand if stranded else intervals_chrom,
        stranded=stranded,
    )
    z = sites.merge(ann, on="row_id", how="left", validate="one_to_one")

    qc = {
        "input_unique_sites": int(len(z)),
        "unique_gene_sites": int(z["annotation_status"].eq("unique_gene").sum()),
        "no_MANE_exon": int(z["annotation_status"].eq("no_MANE_exon").sum()),
        "multi_gene_ambiguous": int(z["annotation_status"].eq("multi_gene_ambiguous").sum()),
        "unique_gene_mapping_fraction": float(
            z["annotation_status"].eq("unique_gene").mean()
        ),
    }

    if qc["unique_gene_mapping_fraction"] < 0.50:
        observed = sorted(set(sites["chrom"].astype(str)))[:50]
        reference_keys = sorted(
            set(k[0] for k in intervals_strand.keys())
            if stranded else set(intervals_chrom.keys())
        )[:50]
        raise RuntimeError(
            f"Map annotation fraction unexpectedly low: {qc}. "
            f"Observed map chromosomes={observed}; "
            f"reference chromosomes={reference_keys}. "
            "This usually indicates a genome-build/chromosome naming mismatch."
        )

    return z, qc


def gene_burden_table(
    gene_table: pd.DataFrame,
    intervals_strand: Dict,
    intervals_chrom: Dict,
):
    burden = gene_table[[
        "gene_name","gene_id","transcript_id","transcript_length"
    ]].copy()

    qcs = {}
    site_tables = {}

    for assay, path in [
        ("BID", BID),
        ("BACS", BACS),
        ("ELAP", ELAP),
    ]:
        s0 = source_exact_map(path, assay)
        ann, qc = annotate_map_to_genes(
            s0, intervals_strand, intervals_chrom, stranded=True
        )
        qcs[assay] = qc
        site_tables[assay] = ann
        counts = (
            ann[ann["annotation_status"] == "unique_gene"]
            .groupby("mane_gene").size()
        )
        burden[f"{assay}_count"] = (
            burden["gene_name"].map(counts).fillna(0).astype(int)
        )

    d0 = drs_hela_map()
    dann, dqc = annotate_map_to_genes(
        d0, intervals_strand, intervals_chrom, stranded=False
    )
    qcs["DRS"] = dqc
    site_tables["DRS"] = dann
    dcounts = (
        dann[dann["annotation_status"] == "unique_gene"]
        .groupby("mane_gene").size()
    )
    burden["DRS_count"] = (
        burden["gene_name"].map(dcounts).fillna(0).astype(int)
    )

    for assay in ["BID","BACS","ELAP","DRS"]:
        burden[f"{assay}_density_per_kb"] = (
            1000.0 * burden[f"{assay}_count"]
            / burden["transcript_length"].clip(lower=1)
        )

    return burden, qcs, site_tables


# -----------------------------------------------------------------------------
# Technology-aware Phase2B burden
# -----------------------------------------------------------------------------

def technology_aware_burden(
    gene_table: pd.DataFrame,
    intervals_strand: Dict,
) -> Tuple[pd.DataFrame, Dict]:
    if not PHASE2B_UNION.exists():
        raise RuntimeError(f"Missing frozen Phase2B union: {PHASE2B_UNION}")

    x = pd.read_csv(PHASE2B_UNION, sep="\t", low_memory=False)
    req = {"chrom","pos1","strand","chemistry_breadth"}
    if not req.issubset(x.columns):
        raise RuntimeError(
            f"Phase2B union missing {sorted(req-set(x.columns))}"
        )

    q = pd.DataFrame({
        "chrom": x["chrom"].map(norm_chr),
        "pos1": pd.to_numeric(x["pos1"], errors="coerce"),
        "strand": x["strand"].astype(str).str.strip(),
        "chemistry_breadth": pd.to_numeric(
            x["chemistry_breadth"], errors="coerce"
        ),
    })
    q = q[
        q["chrom"].ne("")
        & q["pos1"].notna()
        & q["strand"].isin(["+","-"])
        & q["chemistry_breadth"].isin([1,2,3])
    ].copy()
    q["pos1"] = q["pos1"].astype(int)
    q = q.drop_duplicates(
        ["chrom","strand","pos1"]
    ).reset_index(drop=True)
    q["row_id"] = np.arange(len(q))

    ann = annotate_positions(
        q[["row_id","chrom","pos1","strand"]],
        intervals_strand,
        stranded=True
    )
    q = q.merge(ann, on="row_id", how="left", validate="one_to_one")
    u = q[q["annotation_status"] == "unique_gene"].copy()

    agg = u.groupby("mane_gene").agg(
        naive_singleU_union_count=("chemistry_breadth","size"),
        orthogonal_singleU_count=(
            "chemistry_breadth", lambda v: int((v>=2).sum())
        ),
        triple_singleU_count=(
            "chemistry_breadth", lambda v: int((v==3).sum())
        ),
        breadth_sum_singleU=("chemistry_breadth","sum"),
    )

    out = gene_table[["gene_name"]].copy()
    for c in [
        "naive_singleU_union_count",
        "orthogonal_singleU_count",
        "triple_singleU_count",
        "breadth_sum_singleU",
    ]:
        out[c] = out["gene_name"].map(agg[c]).fillna(0).astype(int)

    qc = {
        "phase2b_union_rows": int(len(q)),
        "unique_gene_mapped": int(len(u)),
        "mapping_fraction": float(len(u)/len(q)) if len(q) else np.nan,
        "genes_with_any_union_site": int(
            (out["naive_singleU_union_count"] > 0).sum()
        ),
        "genes_with_orthogonal_site": int(
            (out["orthogonal_singleU_count"] > 0).sum()
        ),
    }
    return out, qc


# -----------------------------------------------------------------------------
# S3 functional outcome resolution
# -----------------------------------------------------------------------------

def best_id_column(
    df: pd.DataFrame,
    gene_name_map: Dict[str,str],
    gene_id_map: Dict[str,str],
    tx_id_map: Dict[str,str],
):
    rows = []

    for c in df.columns:
        if c == "source_sheet":
            continue
        vals = df[c].dropna().astype(str)
        vals = vals[vals.str.strip().ne("")]
        if len(vals) == 0:
            continue
        sample = vals.iloc[: min(2500, len(vals))]

        n = len(sample)
        gene_name_hits = sum(
            v.strip().upper() in gene_name_map for v in sample
        )
        gene_id_hits = sum(
            strip_version(v) in gene_id_map for v in sample
        )
        tx_hits = sum(
            strip_version(v) in tx_id_map for v in sample
        )

        best_type, best_hits = max(
            [
                ("gene_name",gene_name_hits),
                ("gene_id",gene_id_hits),
                ("transcript_id",tx_hits),
            ],
            key=lambda z: z[1]
        )
        rows.append({
            "column": c,
            "n_sample": n,
            "best_type": best_type,
            "matches": best_hits,
            "match_fraction": best_hits/n,
        })

    audit = pd.DataFrame(rows).sort_values(
        ["match_fraction","matches"], ascending=False
    )
    if len(audit) == 0 or audit.iloc[0]["match_fraction"] < 0.30:
        raise RuntimeError(
            "Could not identify an S3 MANE-compatible identifier column."
        )

    top = audit.iloc[0]
    return str(top["column"]), str(top["best_type"]), audit


def map_identifier_to_gene(
    value: Any,
    id_type: str,
    gene_name_map: Dict[str,str],
    gene_id_map: Dict[str,str],
    tx_id_map: Dict[str,str],
) -> str:
    z = ss(value)
    if id_type == "gene_name":
        return gene_name_map.get(z.upper(), "")
    if id_type == "gene_id":
        return gene_id_map.get(strip_version(z), "")
    if id_type == "transcript_id":
        return tx_id_map.get(strip_version(z), "")
    return ""


def resolve_hela_te(
    s3_path: Path,
    gene_name_map: Dict[str,str],
    gene_id_map: Dict[str,str],
    tx_id_map: Dict[str,str],
):
    df = pd.read_csv(s3_path, sep="\t", low_memory=False)
    if "source_sheet" not in df.columns:
        raise RuntimeError("Normalized S3 table lacks source_sheet.")

    schema_rows = []
    for sheet, z in df.groupby("source_sheet", dropna=False):
        for c in df.columns:
            if c == "source_sheet":
                continue
            n = int(z[c].notna().sum())
            if n:
                numfrac = float(
                    pd.to_numeric(z[c], errors="coerce").notna().mean()
                )
                schema_rows.append({
                    "source_sheet": sheet,
                    "column": c,
                    "n_nonmissing": n,
                    "numeric_fraction": numfrac,
                })
    schema = pd.DataFrame(schema_rows)
    schema.to_csv(
        META / "phase3a_S3_schema.tsv", sep="\t", index=False
    )

    sheet_values = [str(v) for v in df["source_sheet"].dropna().unique()]
    te_candidates = [
        v for v in sheet_values
        if "te" in v.lower()
        or "translation" in v.lower()
    ]
    if len(te_candidates) != 1:
        raise RuntimeError(
            f"Could not uniquely resolve S3 TE sheet from {sheet_values}; "
            rf"see {META / 'phase3a_S3_schema.tsv'}"
        )
    te_sheet = te_candidates[0]
    te = df[df["source_sheet"].astype(str).eq(te_sheet)].copy()

    id_col, id_type, id_audit = best_id_column(
        te, gene_name_map, gene_id_map, tx_id_map
    )
    id_audit.to_csv(
        META / "phase3a_S3_identifier_audit.tsv",
        sep="\t", index=False
    )
    te["gene_name"] = te[id_col].map(
        lambda v: map_identifier_to_gene(
            v, id_type, gene_name_map, gene_id_map, tx_id_map
        )
    )

    # The paper specifies these five public HeLa TE sources.
    tokens = [t.lower() for t in HELA_TE_STUDIES] + ["hela"]
    te_cols = []
    for c in te.columns:
        if c in {"source_sheet", id_col, "gene_name"}:
            continue
        cl = c.lower()
        if not any(t in cl for t in tokens):
            continue
        n_numeric = int(
            pd.to_numeric(te[c], errors="coerce").notna().sum()
        )
        if n_numeric >= 100:
            te_cols.append(c)

    # If generic "HeLa" produced duplicated aliases, preserve only unique
    # numeric vectors by exact equality.
    unique_cols = []
    signatures = []
    for c in te_cols:
        v = pd.to_numeric(te[c], errors="coerce")
        sig = (
            int(v.notna().sum()),
            float(v.mean()) if v.notna().any() else np.nan,
            float(v.std()) if v.notna().sum() > 1 else np.nan,
        )
        duplicate = False
        for c0, sig0 in zip(unique_cols, signatures):
            v0 = pd.to_numeric(te[c0], errors="coerce")
            if v.equals(v0):
                duplicate = True
                break
        if not duplicate:
            unique_cols.append(c)
            signatures.append(sig)
    te_cols = unique_cols

    if len(te_cols) < 2:
        raise RuntimeError(
            f"Resolved only {len(te_cols)} HeLa TE columns: {te_cols}. "
            rf"Review {META / 'phase3a_S3_schema.tsv'}"
        )

    long_rows = []
    dataset_qc = []
    for c in te_cols:
        v = pd.to_numeric(te[c], errors="coerce")
        ok = te["gene_name"].ne("") & v.notna() & (v > 0)
        z = te.loc[ok, ["gene_name"]].copy()
        z["TE_raw"] = v.loc[ok].to_numpy()
        z["TE_log10"] = np.log10(z["TE_raw"])
        z = z.groupby("gene_name", as_index=False).agg(
            TE_raw=("TE_raw","median"),
            TE_log10=("TE_log10","median"),
        )
        z["TE_z"] = zscore(z["TE_log10"])
        z["TE_dataset"] = c
        long_rows.append(z)
        dataset_qc.append({
            "column": c,
            "positive_MANE_genes": int(len(z)),
            "median_TE": float(z["TE_raw"].median()),
        })

    long = pd.concat(long_rows, ignore_index=True)
    qc = {
        "S3_path": str(s3_path),
        "TE_sheet": te_sheet,
        "identifier_column": id_col,
        "identifier_type": id_type,
        "identifier_match_fraction": float(id_audit.iloc[0]["match_fraction"]),
        "HeLa_TE_columns": te_cols,
        "datasets": dataset_qc,
        "stacked_rows_before_gene_covariates": int(len(long)),
        "unique_genes_before_gene_covariates": int(long["gene_name"].nunique()),
    }
    return long, qc


# -----------------------------------------------------------------------------
# Regression
# -----------------------------------------------------------------------------

def build_design(
    df: pd.DataFrame,
    burden_col: str,
):
    d = df[[
        "gene_name","TE_dataset","TE_z",
        burden_col,"transcript_length"
    ]].dropna().copy()

    # gene-level standardization, then copied to repeated TE rows.
    gene = d[["gene_name",burden_col,"transcript_length"]].drop_duplicates(
        "gene_name"
    ).copy()
    gene["burden_z"] = zscore(np.log1p(gene[burden_col]))
    gene["length_z"] = zscore(np.log(gene["transcript_length"].clip(lower=1)))

    d = d.drop(columns=[burden_col,"transcript_length"]).merge(
        gene[["gene_name","burden_z","length_z"]],
        on="gene_name", how="inner", validate="many_to_one"
    )
    d = d.dropna()

    datasets = sorted(d["TE_dataset"].unique())
    base = datasets[0]
    cols = [
        np.ones(len(d)),
        d["burden_z"].to_numpy(float),
        d["length_z"].to_numpy(float),
    ]
    names = ["intercept","burden_z","length_z"]

    for ds in datasets[1:]:
        cols.append((d["TE_dataset"] == ds).astype(float).to_numpy())
        names.append(f"dataset::{ds}")

    X = np.column_stack(cols)
    y = d["TE_z"].to_numpy(float)
    clusters = d["gene_name"].astype(str).to_numpy()

    return d, X, y, names, clusters, base


def cluster_ols(
    df: pd.DataFrame,
    burden_col: str,
) -> Tuple[pd.DataFrame, Dict, np.ndarray, np.ndarray, List[str], np.ndarray]:
    d, X, y, names, clusters, base = build_design(df, burden_col)

    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ (X.T @ y)
    resid = y - X @ beta

    groups = defaultdict(list)
    for i,g in enumerate(clusters):
        groups[g].append(i)

    meat = np.zeros((X.shape[1],X.shape[1]))
    for ids in groups.values():
        sg = X[ids].T @ resid[ids]
        meat += np.outer(sg,sg)

    G = len(groups)
    N,K = X.shape
    correction = (G/max(G-1,1))*((N-1)/max(N-K,1))
    cov = xtx_inv @ meat @ xtx_inv * correction
    se = np.sqrt(np.clip(np.diag(cov),0,None))

    rows = []
    for name,b,s0 in zip(names,beta,se):
        z = b/s0 if s0 > 0 else np.nan
        p = math.erfc(abs(z)/math.sqrt(2)) if np.isfinite(z) else np.nan
        rows.append({
            "term": name,
            "beta": float(b),
            "SE_gene_cluster": float(s0),
            "z": float(z) if np.isfinite(z) else np.nan,
            "p": float(p) if np.isfinite(p) else np.nan,
            "CI_low": float(b-1.96*s0),
            "CI_high": float(b+1.96*s0),
            "n_stacked_rows": N,
            "n_genes": G,
            "n_TE_datasets": int(d["TE_dataset"].nunique()),
        })

    meta = {
        "base_TE_dataset": base,
        "n_stacked_rows": N,
        "n_genes": G,
        "n_TE_datasets": int(d["TE_dataset"].nunique()),
    }
    return pd.DataFrame(rows), meta, X, y, names, clusters


def per_dataset_effects(
    df: pd.DataFrame,
    burden_col: str,
) -> pd.DataFrame:
    rows = []
    for ds, z0 in df.groupby("TE_dataset"):
        z = z0[[
            "gene_name","TE_z",burden_col,"transcript_length"
        ]].dropna().copy()
        if len(z) < 100:
            continue
        z["burden_z"] = zscore(np.log1p(z[burden_col]))
        z["length_z"] = zscore(
            np.log(z["transcript_length"].clip(lower=1))
        )
        z = z.dropna()

        X = np.column_stack([
            np.ones(len(z)),
            z["burden_z"].to_numpy(float),
            z["length_z"].to_numpy(float),
        ])
        y = z["TE_z"].to_numpy(float)
        inv = np.linalg.pinv(X.T @ X)
        beta = inv @ X.T @ y
        resid = y-X@beta

        # HC3
        h = np.sum((X@inv)*X, axis=1)
        adj = resid / np.clip(1-h, 1e-6, None)
        meat = X.T @ ((adj**2)[:,None]*X)
        cov = inv @ meat @ inv
        se = np.sqrt(np.clip(np.diag(cov),0,None))

        b,s0 = beta[1],se[1]
        p = math.erfc(abs(b/s0)/math.sqrt(2)) if s0 > 0 else np.nan
        rows.append({
            "TE_dataset": ds,
            "burden": burden_col,
            "beta_standardized": b,
            "SE_HC3": s0,
            "CI_low": b-1.96*s0,
            "CI_high": b+1.96*s0,
            "p": p,
            "n_genes": len(z),
        })
    return pd.DataFrame(rows)


def bootstrap_map_differences(
    functional: pd.DataFrame,
    burdens: List[str],
    n_boot: int = N_BOOT,
):
    """
    Gene-cluster bootstrap. Each sampled gene receives an integer bootstrap
    weight; all TE datasets for that gene inherit the same weight.
    """
    fitted = {}
    common_genes = sorted(functional["gene_name"].unique())
    gene_to_idx = {g:i for i,g in enumerate(common_genes)}

    # All map models must share identical rows by construction.
    base_d = None
    for bcol in burdens:
        d,X,y,names,clusters,base = build_design(functional, bcol)
        if base_d is None:
            base_d = d[["gene_name","TE_dataset","TE_z"]].reset_index(drop=True)
            y_ref = y
        else:
            chk = d[["gene_name","TE_dataset","TE_z"]].reset_index(drop=True)
            if not chk.equals(base_d):
                raise RuntimeError(
                    "Map models do not share identical stacked rows."
                )
        fitted[bcol] = (X,y,names,clusters)

    row_gene_index = np.array(
        [gene_to_idx[g] for g in fitted[burdens[0]][3]],
        dtype=int
    )
    n_genes = len(common_genes)

    boot_beta = {b: np.empty(n_boot) for b in burdens}

    for ib in range(n_boot):
        sampled = RNG.integers(0, n_genes, size=n_genes)
        counts = np.bincount(sampled, minlength=n_genes).astype(float)
        w = counts[row_gene_index]

        for bcol in burdens:
            X,y,_,_ = fitted[bcol]
            use = w > 0
            X0 = X[use]
            y0 = y[use]
            w0 = w[use]
            xtwx = X0.T @ (w0[:,None]*X0)
            xtwy = X0.T @ (w0*y0)
            beta = np.linalg.pinv(xtwx) @ xtwy
            boot_beta[bcol][ib] = beta[1]

    rows = []
    for i,a in enumerate(burdens):
        for b in burdens[i+1:]:
            diff = boot_beta[a]-boot_beta[b]
            obs_a = cluster_ols(functional,a)[0]
            obs_b = cluster_ols(functional,b)[0]
            ba = float(obs_a.loc[obs_a["term"]=="burden_z","beta"].iloc[0])
            bb = float(obs_b.loc[obs_b["term"]=="burden_z","beta"].iloc[0])
            obs = ba-bb

            lo,hi = np.quantile(diff,[.025,.975])
            p_lo = (np.sum(diff <= 0)+1)/(len(diff)+1)
            p_hi = (np.sum(diff >= 0)+1)/(len(diff)+1)
            p = min(1.0,2*min(p_lo,p_hi))

            rows.append({
                "map_A": a,
                "map_B": b,
                "beta_A_minus_B": obs,
                "bootstrap_CI_low": float(lo),
                "bootstrap_CI_high": float(hi),
                "bootstrap_p_two_sided": float(p),
                "n_boot": n_boot,
            })

    out = pd.DataFrame(rows)
    out["q_BH_pairwise"] = bh_adjust(out["bootstrap_p_two_sided"])
    return out


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    print("="*100)
    print("TRACE PHASE 3A v2 — MAP PROVENANCE AND TRANSLATION-EFFICIENCY INFERENCE")
    print("="*100)

    for p in [GTF,BID,BACS,ELAP,DRS_S6,PHASE2B_UNION]:
        if not p.exists():
            raise RuntimeError(f"Missing required input: {p}")

    s3_path = resolve_s3_path()

    (
        gene_table,
        intervals_strand,
        intervals_chrom,
        gene_name_map,
        gene_id_map,
        tx_id_map,
        mane_qc,
    ) = load_mane()

    burdens, map_qc, site_tables = gene_burden_table(
        gene_table, intervals_strand, intervals_chrom
    )
    tech, tech_qc = technology_aware_burden(
        gene_table, intervals_strand
    )
    burdens = burdens.merge(
        tech, on="gene_name", how="left", validate="one_to_one"
    )

    burdens.to_csv(
        AIM3 / "aim3a_MANE_gene_burdens.tsv",
        sep="\t", index=False
    )

    te_long, s3_qc = resolve_hela_te(
        s3_path, gene_name_map, gene_id_map, tx_id_map
    )

    functional = te_long.merge(
        burdens,
        on="gene_name",
        how="inner",
        validate="many_to_one"
    )
    functional = functional[
        functional["transcript_length"].notna()
        & functional["TE_z"].notna()
    ].copy()

    # Freeze identical functional row set for every map.
    if functional["gene_name"].nunique() < 500:
        raise RuntimeError(
            f"Only {functional['gene_name'].nunique()} functional MANE genes resolved."
        )

    functional.to_csv(
        AIM3 / "aim3a_functional_cohort.tsv",
        sep="\t", index=False
    )

    map_burdens = [
        "BID_count","BACS_count","ELAP_count","DRS_count"
    ]

    stacked_rows = []
    per_ds_rows = []
    for bcol in map_burdens:
        tab, meta, *_ = cluster_ols(functional, bcol)
        eff = tab[tab["term"]=="burden_z"].copy()
        eff["map"] = bcol.replace("_count","")
        eff["burden_column"] = bcol
        eff["model"] = (
            "TE_z ~ standardized_log1p_map_count + log_transcript_length "
            "+ TE_dataset_fixed_effects"
        )
        stacked_rows.append(eff)

        pds = per_dataset_effects(functional,bcol)
        pds["map"] = bcol.replace("_count","")
        per_ds_rows.append(pds)

    stacked = pd.concat(stacked_rows, ignore_index=True)
    stacked["q_BH_across_maps"] = bh_adjust(stacked["p"])
    stacked.to_csv(
        AIM3 / "aim3a_TE_stacked_map_effects.tsv",
        sep="\t", index=False
    )

    per_ds = pd.concat(per_ds_rows, ignore_index=True)
    per_ds.to_csv(
        AIM3 / "aim3a_TE_per_dataset_effects.tsv",
        sep="\t", index=False
    )

    pairwise = bootstrap_map_differences(
        functional, map_burdens, N_BOOT
    )
    pairwise.to_csv(
        AIM3 / "aim3a_map_pairwise_bootstrap.tsv",
        sep="\t", index=False
    )

    # Technology-aware secondary bridge analysis.
    tech_burdens = [
        "naive_singleU_union_count",
        "orthogonal_singleU_count",
        "triple_singleU_count",
        "breadth_sum_singleU",
    ]
    tech_rows = []
    for bcol in tech_burdens:
        tab, meta, *_ = cluster_ols(functional,bcol)
        eff = tab[tab["term"]=="burden_z"].copy()
        eff["feature"] = bcol
        tech_rows.append(eff)

    tech_eff = pd.concat(tech_rows, ignore_index=True)
    tech_eff["q_BH_secondary"] = bh_adjust(tech_eff["p"])
    tech_eff.to_csv(
        AIM3 / "aim3a_technology_aware_secondary.tsv",
        sep="\t", index=False
    )

    # Simple consistency summaries; not a ranking.
    significant_pair_difference = bool(
        (pairwise["q_BH_pairwise"] < .05).any()
    )
    signs = np.sign(stacked["beta"].to_numpy(float))
    sign_heterogeneity = bool(len(set(signs.tolist())) > 1)

    if significant_pair_difference or sign_heterogeneity:
        status = "AIM3A_MAP_PROVENANCE_CHANGES_TE_INFERENCE"
    else:
        status = "AIM3A_NO_CLEAR_MAP_EFFECT_HETEROGENEITY"

    contract = {
        "phase": "3A_v2",
        "scientific_story": "Calibration -> Measurement -> Evidence -> Inference",
        "primary_question": (
            "Does the Ψ–translation-efficiency association change when only "
            "the underlying published Ψ map is changed?"
        ),
        "cell_line": "HeLa",
        "functional_source": (
            "NAR 2026 Supplementary Table S3 TE values; multiple public HeLa "
            "Ribo-seq/RNA-seq datasets are analyzed as repeated outcome definitions"
        ),
        "maps": ["BID","BACS","ELAP","DRS"],
        "same_cohort_policy": (
            "All four map models use identical stacked gene-by-TE-dataset rows "
            "and the same MANE transcript-length covariate."
        ),
        "primary_model": (
            "within-dataset z(log10 TE) ~ standardized log1p(map-specific Ψ count) "
            "+ standardized log(MANE transcript length) + dataset fixed effects; "
            "gene-cluster robust SE"
        ),
        "map_comparison": (
            "gene-cluster bootstrap pairwise differences of standardized map-burden "
            "effects; BH correction across six map pairs"
        ),
        "secondary_technology_aware_features": tech_burdens,
        "interpretation_boundary": (
            "0 map count means no site assigned to the MANE gene in that published "
            "map, not biological absence. Associations are not causal."
        ),
        "next_step": (
            "Phase3B final sequence-baseline test M0/M1/M2. Phase3A determines "
            "whether map provenance visibly changes the functional association; "
            "Phase3B asks whether technology-aware Ψ adds information beyond sequence."
        ),
        "no_raw_sequencing": True,
    }
    (META / "phase3a_v2_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    summary = {
        "status": status,
        "MANE_QC": mane_qc,
        "map_annotation_QC": map_qc,
        "technology_aware_QC": tech_qc,
        "S3_QC": s3_qc,
        "functional_stacked_rows": int(len(functional)),
        "functional_unique_genes": int(functional["gene_name"].nunique()),
        "n_TE_datasets": int(functional["TE_dataset"].nunique()),
        "map_effects": stacked.to_dict(orient="records"),
        "pairwise_map_effect_differences": pairwise.to_dict(orient="records"),
        "technology_aware_secondary": tech_eff.to_dict(orient="records"),
        "significant_pairwise_map_difference_after_BH": significant_pair_difference,
        "map_effect_sign_heterogeneity": sign_heterogeneity,
        "next_gate": (
            "Freeze Aim3A regardless of direction. Proceed once to Phase3B "
            "M0/M1/M2 with an explicit sequence baseline; do not add additional "
            "functional outcomes unless TE schema itself fails."
        ),
    }
    (META / "phase3a_v2_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print("Functional genes:", functional["gene_name"].nunique())
    print("HeLa TE datasets:", functional["TE_dataset"].nunique())
    print()
    print("MAP EFFECTS")
    print(stacked.to_string(index=False))
    print()
    print("PAIRWISE MAP EFFECT DIFFERENCES")
    print(pairwise.to_string(index=False))
    print()
    print("TECHNOLOGY-AWARE SECONDARY")
    print(tech_eff.to_string(index=False))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase3a_v2_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase3a_v2_analysis_contract.json")
    print(r"  D:\RNA\Trace\00_meta\phase3a_S3_schema.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3a_TE_stacked_map_effects.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3a_TE_per_dataset_effects.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3a_map_pairwise_bootstrap.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3a_technology_aware_secondary.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase3a_v2_FATAL_ERROR.txt").write_text(
            err, encoding="utf-8"
        )
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase3a_v2_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
