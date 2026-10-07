#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRACE Phase 3B — final sequence-baseline test
==============================================

Final scientific question
-------------------------
Within each HeLa translation-efficiency (TE) dataset separately:

    M0: sequence only
    M1: sequence + naive Ψ burden
    M2: sequence + naive Ψ burden + provenance-aware confirmation

Does Ψ add out-of-fold predictive information beyond a biologically strong
sequence baseline, and does cross-technology confirmation add information
beyond total Ψ burden?

Why this design
---------------
Phase 3A closure showed that the five public HeLa TE datasets are not one
homogeneous outcome. Therefore this script NEVER fits a universal pooled TE
model. Each TE dataset is analyzed independently with identical gene folds,
then predictive increments are summarized across datasets.

M0 sequence baseline
--------------------
Features are extracted from the MANE v1.5 RefSeq transcript sequence and its
annotated coding span:
- transcript / 5'UTR / CDS / 3'UTR lengths
- regional GC and mono-nucleotide composition
- 5'UTR and 3'UTR dinucleotide frequencies
- CDS sense-codon composition
- upstream AUG / uORF counts
- Kozak -3 purine and +4 G indicators
- local start- and stop-codon positional one-hot windows

This is intentionally a classical, regularized sequence model rather than a
foundation model. It follows current TE literature showing that UTR sequence,
local initiation context and codon composition carry substantial TE signal.

M1
--
M0 + log1p(number of unique locally-single-U Ψ loci reported by at least one
of BID/BACS/ELAP).

M2
--
M1 + log1p(number of ADDITIONAL assay confirmations beyond the unique-site
count):

    confirmation_excess = breadth_sum_singleU - naive_singleU_union_count

Thus M2 is nested in M1 and asks whether repeated cross-technology support
adds predictive information after total Ψ burden is already known.

Important boundaries
--------------------
- DRS is NOT included in M1/M2 because it served as the held-out platform in
  Aim 2B.
- 0 Ψ burden means no site reported in the frozen BID/BACS/ELAP source union,
  not biological absence.
- No raw sequencing.
- No threshold tuning from final results.
- Same global gene folds are used for M0/M1/M2 and across TE datasets.
- Hyperparameter selection is nested inside the outer CV.
- Primary metric: out-of-fold R^2.
- Primary increments:
      ΔR2_10 = R2(M1) - R2(M0)
      ΔR2_21 = R2(M2) - R2(M1)
- TE-study heterogeneity is preserved; macro summaries do not replace
  per-dataset results.

Project root
------------
D:\\RNA\\Trace
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
import traceback
import urllib.request
from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from urllib.parse import unquote

import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr
    from sklearn.linear_model import Ridge
    from sklearn.metrics import mean_squared_error, r2_score
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
except Exception as e:
    raise RuntimeError(
        "Phase3B requires scipy and scikit-learn."
    ) from e


ROOT = Path(r"D:\RNA\Trace")
META = ROOT / "00_meta"
REF = ROOT / "02_reference"
AIM3 = ROOT / "06_aim3"
LOGS = ROOT / "logs"

for p in [META, REF, AIM3, LOGS]:
    p.mkdir(parents=True, exist_ok=True)

GTF = REF / "MANE.GRCh38.v1.5.refseq_genomic.gtf.gz"
FASTA = REF / "MANE.GRCh38.v1.5.refseq_rna.fna.gz"
FASTA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
    "MANE.GRCh38.v1.5.refseq_rna.fna.gz"
)

FUNCTIONAL = AIM3 / "aim3a_functional_cohort.tsv"

SEED = 20261004
N_OUTER = 5
N_BOOT = 3000
ALPHAS = np.array([1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0], float)

DNA = "ACGT"
DINUCS = ["".join(x) for x in product(DNA, repeat=2)]
STOP_CODONS = {"TAA", "TAG", "TGA"}
CODONS = ["".join(x) for x in product(DNA, repeat=3)]
SENSE_CODONS = [c for c in CODONS if c not in STOP_CODONS]

# RiboNN-style local positional context, kept modest for ~800 genes.
START_OFFSETS = list(range(-12, 16))   # -12 .. +15 relative to start A
STOP_OFFSETS = list(range(-15, 13))    # -15 .. +12 relative to stop first base


# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------

def ss(x: Any) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and np.isnan(x):
        return ""
    return str(x).strip()


def strip_version(x: Any) -> str:
    return re.sub(r"\.\d+$", "", ss(x))


def parse_attrs(text: str) -> Dict[str, str]:
    raw = ss(text)
    out: Dict[str, List[str]] = defaultdict(list)

    if "=" in raw:
        for field in raw.split(";"):
            field = field.strip()
            if not field or "=" not in field:
                continue
            k, v = field.split("=", 1)
            out[k.strip()].append(unquote(v.strip().strip('"')))

    for k, v in re.findall(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"', raw):
        if v not in out[k]:
            out[k].append(v)

    return {k: ",".join(vs) for k, vs in out.items()}


def attr_first(a: Dict[str, str], keys: Sequence[str]) -> str:
    for k in keys:
        v = ss(a.get(k, ""))
        if v:
            return v.split(",")[0].strip()
    return ""


def transcript_id_from_attrs(a: Dict[str, str]) -> str:
    v = attr_first(a, ["transcript_id", "transcript"])
    if v:
        return re.sub(r"^rna-", "", v)
    parent = attr_first(a, ["Parent"])
    if parent.startswith("rna-"):
        return parent[4:]
    ident = attr_first(a, ["ID"])
    if ident.startswith("rna-"):
        return ident[4:]
    return ""


def is_mane_select(a: Dict[str, str]) -> bool:
    tag = ss(a.get("tag", "")).replace("_", " ").lower()
    return "mane select" in tag


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def kmer_freq(seq: str, kmers: Sequence[str]) -> Dict[str, float]:
    if not seq:
        return {k: 0.0 for k in kmers}
    k = len(kmers[0])
    denom = max(len(seq) - k + 1, 0)
    if denom <= 0:
        return {x: 0.0 for x in kmers}

    counts = defaultdict(int)
    for i in range(denom):
        token = seq[i:i+k]
        if set(token) <= set(DNA):
            counts[token] += 1

    return {x: counts[x] / denom for x in kmers}


def mono_freq(seq: str) -> Dict[str, float]:
    if not seq:
        return {b: 0.0 for b in DNA}
    n = len(seq)
    return {b: seq.count(b) / n for b in DNA}


def uorf_count(utr5: str) -> int:
    """
    Count AUG-initiated upstream ORFs with an in-frame stop before the main AUG.
    Overlapping candidate AUGs are counted independently.
    """
    n = 0
    for i in range(max(len(utr5) - 2, 0)):
        if utr5[i:i+3] != "ATG":
            continue
        found = False
        for j in range(i+3, len(utr5)-2, 3):
            if utr5[j:j+3] in STOP_CODONS:
                found = True
                break
        if found:
            n += 1
    return n


def positional_onehot(
    seq: str,
    anchor: int,
    offsets: Sequence[int],
    prefix: str,
) -> Dict[str, float]:
    out = {}
    for off in offsets:
        idx = anchor + off
        base = seq[idx] if 0 <= idx < len(seq) else "N"
        for b in DNA:
            out[f"{prefix}_{off:+d}_{b}"] = float(base == b)
    return out


# -----------------------------------------------------------------------------
# Reference FASTA
# -----------------------------------------------------------------------------

def ensure_fasta() -> Path:
    if FASTA.exists():
        try:
            with gzip.open(FASTA, "rt", encoding="utf-8") as f:
                first = f.readline()
            if first.startswith(">"):
                return FASTA
        except Exception:
            pass

    tmp = FASTA.with_suffix(FASTA.suffix + ".tmp")
    req = urllib.request.Request(
        FASTA_URL,
        headers={"User-Agent": "Mozilla/5.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as r, tmp.open("wb") as w:
            while True:
                block = r.read(1024 * 1024)
                if not block:
                    break
                w.write(block)

        with gzip.open(tmp, "rt", encoding="utf-8") as f:
            first = f.readline()
            if not first.startswith(">"):
                raise RuntimeError("Downloaded MANE RNA FASTA is invalid.")

        tmp.replace(FASTA)
        return FASTA
    except Exception as e:
        try:
            if tmp.exists():
                tmp.unlink()
        except Exception:
            pass
        note = META / "phase3b_MANERNA_download_required.txt"
        note.write_text(
            f"Automatic download failed: {repr(e)}\n\n"
            f"Official file:\n{FASTA_URL}\n\n"
            f"Save it as:\n{FASTA}\n",
            encoding="utf-8",
        )
        raise RuntimeError(
            "MANE v1.5 RefSeq RNA FASTA is required. "
            f"See {note}"
        ) from e


def load_target_fasta(
    path: Path,
    target_ids: Sequence[str],
) -> Dict[str, str]:
    full_targets = set(map(ss, target_ids))
    stripped_targets = {strip_version(x) for x in target_ids}

    seqs: Dict[str, str] = {}
    header = None
    chunks: List[str] = []

    def flush():
        nonlocal header, chunks
        if header is None:
            return
        acc = header.split()[0]
        key_stripped = strip_version(acc)
        if acc in full_targets or key_stripped in stripped_targets:
            seq = "".join(chunks).upper().replace("U", "T")
            seqs[acc] = seq
            seqs[key_stripped] = seq

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        flush()

    return seqs


# -----------------------------------------------------------------------------
# MANE transcript structure
# -----------------------------------------------------------------------------

def parse_target_transcripts(
    target_ids: Sequence[str],
) -> Dict[str, Dict]:
    targets_full = set(map(ss, target_ids))
    targets_strip = {strip_version(x) for x in target_ids}

    recs: Dict[str, Dict] = {}

    with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attrs = p
            a = parse_attrs(attrs)
            if not is_mane_select(a):
                continue

            tx = transcript_id_from_attrs(a)
            if not tx:
                continue
            txs = strip_version(tx)
            if tx not in targets_full and txs not in targets_strip:
                continue

            key = tx if tx in targets_full else txs
            r = recs.setdefault(
                key,
                {
                    "transcript_id": tx,
                    "strand": strand,
                    "exons": [],
                    "CDS": [],
                    "start_codon": [],
                    "stop_codon": [],
                },
            )

            try:
                st, en = int(start), int(end)
            except Exception:
                continue

            if feature == "exon":
                r["exons"].append((st, en))
            elif feature == "CDS":
                r["CDS"].append((st, en))
            elif feature == "start_codon":
                r["start_codon"].append((st, en))
            elif feature == "stop_codon":
                r["stop_codon"].append((st, en))

    return recs


def project_interval_to_tx(
    interval: Tuple[int, int],
    exon_map: Sequence[Tuple[int, int, int]],
    strand: str,
) -> List[Tuple[int, int]]:
    """
    exon_map entries: genomic_start, genomic_end, transcript_start.
    """
    st, en = interval
    out = []

    for ex_st, ex_en, tx_st in exon_map:
        ov_st = max(st, ex_st)
        ov_en = min(en, ex_en)
        if ov_st > ov_en:
            continue

        if strand == "+":
            a = tx_st + (ov_st - ex_st)
            b = tx_st + (ov_en - ex_st)
        else:
            a = tx_st + (ex_en - ov_en)
            b = tx_st + (ex_en - ov_st)

        out.append((min(a, b), max(a, b)))

    return out


def transcript_regions(
    rec: Dict,
    seq_len: int,
) -> Optional[Tuple[int, int]]:
    """
    Return coding_start_tx, coding_end_tx inclusive.
    Prefer start/stop codons; fall back to CDS span.
    """
    exons = sorted(
        rec["exons"],
        key=lambda z: z[0],
        reverse=(rec["strand"] == "-"),
    )
    if not exons:
        return None

    exon_map = []
    cursor = 0
    for st, en in exons:
        exon_map.append((st, en, cursor))
        cursor += en - st + 1

    if abs(cursor - seq_len) > 3:
        return None

    start_parts = []
    for iv in rec["start_codon"]:
        start_parts.extend(project_interval_to_tx(iv, exon_map, rec["strand"]))

    stop_parts = []
    for iv in rec["stop_codon"]:
        stop_parts.extend(project_interval_to_tx(iv, exon_map, rec["strand"]))

    if start_parts and stop_parts:
        coding_start = min(a for a, b in start_parts)
        coding_end = max(b for a, b in stop_parts)
    else:
        cds_parts = []
        for iv in rec["CDS"]:
            cds_parts.extend(project_interval_to_tx(iv, exon_map, rec["strand"]))
        if not cds_parts:
            return None
        coding_start = min(a for a, b in cds_parts)
        coding_end = max(b for a, b in cds_parts)

    if not (0 <= coding_start <= coding_end < seq_len):
        return None

    return int(coding_start), int(coding_end)


# -----------------------------------------------------------------------------
# Sequence features
# -----------------------------------------------------------------------------

def extract_sequence_features(
    genes: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict]:
    fasta_path = ensure_fasta()

    target_ids = genes["transcript_id"].astype(str).tolist()
    seqs = load_target_fasta(fasta_path, target_ids)
    structs = parse_target_transcripts(target_ids)

    rows = []
    failures = []

    for r in genes.itertuples():
        tx = ss(r.transcript_id)
        txs = strip_version(tx)

        seq = seqs.get(tx) or seqs.get(txs)
        rec = structs.get(tx) or structs.get(txs)

        if not seq or rec is None:
            failures.append({
                "gene_name": r.gene_name,
                "transcript_id": tx,
                "reason": "missing_fasta_or_gtf_record",
            })
            continue

        region = transcript_regions(rec, len(seq))
        if region is None:
            failures.append({
                "gene_name": r.gene_name,
                "transcript_id": tx,
                "reason": "coding_span_projection_failed",
            })
            continue

        coding_start, coding_end = region
        utr5 = seq[:coding_start]
        coding = seq[coding_start:coding_end+1]
        utr3 = seq[coding_end+1:]

        unknown_fraction = (
            sum(b not in DNA for b in seq) / max(len(seq), 1)
        )
        if unknown_fraction > 0.01:
            failures.append({
                "gene_name": r.gene_name,
                "transcript_id": tx,
                "reason": f"ambiguous_sequence_fraction={unknown_fraction:.4f}",
            })
            continue

        # Require a plausible annotated start codon for the final baseline.
        start_triplet = seq[coding_start:coding_start+3]
        if start_triplet != "ATG":
            failures.append({
                "gene_name": r.gene_name,
                "transcript_id": tx,
                "reason": f"noncanonical_or_misaligned_start={start_triplet}",
            })
            continue

        feat: Dict[str, float] = {
            "gene_name": r.gene_name,
            "transcript_id": tx,
            "seq_total_len": float(len(seq)),
            "seq_utr5_len": float(len(utr5)),
            "seq_cds_span_len": float(len(coding)),
            "seq_utr3_len": float(len(utr3)),
            "seq_log_total_len": math.log1p(len(seq)),
            "seq_log_utr5_len": math.log1p(len(utr5)),
            "seq_log_cds_len": math.log1p(len(coding)),
            "seq_log_utr3_len": math.log1p(len(utr3)),
            "seq_gc_total": gc_fraction(seq),
            "seq_gc_utr5": gc_fraction(utr5),
            "seq_gc_cds": gc_fraction(coding),
            "seq_gc_utr3": gc_fraction(utr3),
            "seq_uAUG_count": float(
                sum(utr5[i:i+3] == "ATG" for i in range(max(len(utr5)-2, 0)))
            ),
            "seq_uORF_count": float(uorf_count(utr5)),
            "seq_kozak_minus3_purine": float(
                coding_start >= 3 and seq[coding_start-3] in {"A", "G"}
            ),
            "seq_kozak_plus4_G": float(
                coding_start + 3 < len(seq) and seq[coding_start+3] == "G"
            ),
        }

        for region_name, region_seq in [
            ("utr5", utr5),
            ("cds", coding),
            ("utr3", utr3),
        ]:
            mf = mono_freq(region_seq)
            for b in DNA:
                feat[f"seq_{region_name}_mono_{b}"] = mf[b]

        for region_name, region_seq in [("utr5", utr5), ("utr3", utr3)]:
            dfreq = kmer_freq(region_seq, DINUCS)
            for k in DINUCS:
                feat[f"seq_{region_name}_dinuc_{k}"] = dfreq[k]

        # Sense-codon composition from the annotated coding span.
        codon_seq = coding
        if len(codon_seq) >= 3 and codon_seq[-3:] in STOP_CODONS:
            sense_seq = codon_seq[:-3]
        else:
            sense_seq = codon_seq

        codon_counts = defaultdict(int)
        n_codons = 0
        for i in range(0, len(sense_seq)-2, 3):
            c = sense_seq[i:i+3]
            if c in SENSE_CODONS:
                codon_counts[c] += 1
                n_codons += 1

        for c in SENSE_CODONS:
            feat[f"seq_codon_{c}"] = (
                codon_counts[c] / n_codons if n_codons else 0.0
            )

        feat.update(
            positional_onehot(
                seq, coding_start, START_OFFSETS, "seq_startpos"
            )
        )
        stop_anchor = max(coding_end - 2, coding_start)
        feat.update(
            positional_onehot(
                seq, stop_anchor, STOP_OFFSETS, "seq_stoppos"
            )
        )

        rows.append(feat)

    failures_df = pd.DataFrame(failures)
    failures_df.to_csv(
        META / "phase3b_sequence_feature_failures.tsv",
        sep="\t", index=False
    )

    out = pd.DataFrame(rows)
    mapping_fraction = len(out) / max(len(genes), 1)

    if mapping_fraction < 0.95:
        raise RuntimeError(
            f"Only {mapping_fraction:.2%} of Phase3B genes received valid MANE "
            "sequence features. See phase3b_sequence_feature_failures.tsv"
        )

    # MANE transcript length should match the previously frozen transcript length.
    check = out.merge(
        genes[["gene_name", "transcript_length"]],
        on="gene_name",
        how="left",
        validate="one_to_one",
    )
    diff = (
        check["seq_total_len"]
        - pd.to_numeric(check["transcript_length"], errors="coerce")
    ).abs()

    exact_fraction = float((diff <= 3).mean())
    if exact_fraction < 0.95:
        raise RuntimeError(
            f"Only {exact_fraction:.2%} of sequence lengths agree with frozen "
            "MANE transcript length within 3 nt."
        )

    qc = {
        "input_genes": int(len(genes)),
        "sequence_feature_genes": int(len(out)),
        "mapping_fraction": mapping_fraction,
        "length_agreement_within_3nt_fraction": exact_fraction,
        "n_failures": int(len(failures_df)),
        "MANE_RNA_FASTA": str(fasta_path),
        "n_raw_sequence_features": int(
            len([c for c in out.columns if c.startswith("seq_")])
        ),
    }

    return out, qc


# -----------------------------------------------------------------------------
# Global gene folds
# -----------------------------------------------------------------------------

def assign_global_folds(genes: Sequence[str]) -> Dict[str, int]:
    """
    Deterministic shuffled 5-fold assignment shared by every TE dataset/model.
    """
    genes = sorted(set(map(str, genes)))
    rng = np.random.default_rng(SEED)
    perm = np.array(genes, dtype=object)
    rng.shuffle(perm)

    fold_map = {}
    for i, g in enumerate(perm):
        fold_map[str(g)] = int(i % N_OUTER)
    return fold_map


# -----------------------------------------------------------------------------
# Nested Ridge CV
# -----------------------------------------------------------------------------

def make_model(alpha: float) -> Pipeline:
    return Pipeline([
        ("scale", StandardScaler()),
        ("ridge", Ridge(alpha=float(alpha), fit_intercept=True)),
    ])


def tune_alpha(
    d: pd.DataFrame,
    features: Sequence[str],
    outer_fold: int,
) -> Tuple[float, pd.DataFrame]:
    train = d[d["outer_fold"] != outer_fold].copy()
    inner_folds = sorted(train["outer_fold"].unique())

    rows = []

    for alpha in ALPHAS:
        mses = []

        for inner_val in inner_folds:
            tr = train[train["outer_fold"] != inner_val]
            va = train[train["outer_fold"] == inner_val]

            if len(tr) < 50 or len(va) < 20:
                continue

            model = make_model(alpha)
            model.fit(tr[list(features)], tr["TE_z"])
            pred = model.predict(va[list(features)])
            mses.append(mean_squared_error(va["TE_z"], pred))

        rows.append({
            "outer_fold": outer_fold,
            "alpha": float(alpha),
            "inner_mean_MSE": float(np.mean(mses)) if mses else np.nan,
            "inner_folds_used": len(mses),
        })

    tab = pd.DataFrame(rows)
    valid = tab.dropna(subset=["inner_mean_MSE"]).sort_values(
        ["inner_mean_MSE", "alpha"]
    )

    if len(valid) == 0:
        raise RuntimeError(
            f"No valid inner-CV alpha for outer fold {outer_fold}."
        )

    return float(valid.iloc[0]["alpha"]), tab


def nested_oof(
    d: pd.DataFrame,
    features: Sequence[str],
    model_name: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    preds = []
    tune_rows = []

    for outer in range(N_OUTER):
        tr = d[d["outer_fold"] != outer].copy()
        te = d[d["outer_fold"] == outer].copy()

        if len(tr) < 100 or len(te) < 30:
            raise RuntimeError(
                f"{model_name}: outer fold {outer} has train={len(tr)}, test={len(te)}"
            )

        alpha, tuning = tune_alpha(d, features, outer)
        tuning["model"] = model_name
        tune_rows.append(tuning)

        model = make_model(alpha)
        model.fit(tr[list(features)], tr["TE_z"])
        pred = model.predict(te[list(features)])

        q = te[["gene_name", "TE_dataset", "TE_z", "outer_fold"]].copy()
        q["model"] = model_name
        q["prediction"] = pred
        q["selected_alpha"] = alpha
        preds.append(q)

    return (
        pd.concat(preds, ignore_index=True),
        pd.concat(tune_rows, ignore_index=True),
    )


# -----------------------------------------------------------------------------
# Performance / bootstrap
# -----------------------------------------------------------------------------

def metrics(y: np.ndarray, pred: np.ndarray) -> Dict[str, float]:
    rho = spearmanr(y, pred).statistic if len(y) >= 3 else np.nan
    return {
        "R2": float(r2_score(y, pred)),
        "RMSE": float(math.sqrt(mean_squared_error(y, pred))),
        "Spearman": float(rho),
    }


def dataset_performance(pred_long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (ds, model), z in pred_long.groupby(["TE_dataset", "model"]):
        m = metrics(
            z["TE_z"].to_numpy(float),
            z["prediction"].to_numpy(float),
        )
        rows.append({
            "TE_dataset": ds,
            "model": model,
            "n_genes": len(z),
            **m,
        })
    return pd.DataFrame(rows)


def paired_prediction_table(pred_long: pd.DataFrame, ds: str) -> pd.DataFrame:
    z = pred_long[pred_long["TE_dataset"] == ds].copy()
    y = (
        z[["gene_name", "TE_z", "outer_fold"]]
        .drop_duplicates("gene_name")
    )
    p = z.pivot(
        index="gene_name",
        columns="model",
        values="prediction",
    ).reset_index()
    out = y.merge(p, on="gene_name", how="inner", validate="one_to_one")
    required = {"M0_sequence", "M1_sequence_plus_naivePsi", "M2_plus_confirmation"}
    if not required.issubset(out.columns):
        raise RuntimeError(f"{ds}: missing model predictions {required-set(out.columns)}")
    return out


def compute_deltas(z: pd.DataFrame) -> Dict[str, float]:
    y = z["TE_z"].to_numpy(float)
    r0 = r2_score(y, z["M0_sequence"])
    r1 = r2_score(y, z["M1_sequence_plus_naivePsi"])
    r2 = r2_score(y, z["M2_plus_confirmation"])
    return {
        "R2_M0": float(r0),
        "R2_M1": float(r1),
        "R2_M2": float(r2),
        "delta_R2_M1_minus_M0": float(r1-r0),
        "delta_R2_M2_minus_M1": float(r2-r1),
    }


def bootstrap_dataset_deltas(
    z: pd.DataFrame,
    n_boot: int = N_BOOT,
) -> Dict[str, float]:
    n = len(z)
    rng = np.random.default_rng(SEED + n)
    vals10 = []
    vals21 = []

    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        b = z.iloc[idx]

        if np.var(b["TE_z"].to_numpy(float)) < 1e-10:
            continue

        d = compute_deltas(b)
        vals10.append(d["delta_R2_M1_minus_M0"])
        vals21.append(d["delta_R2_M2_minus_M1"])

    def ci(vals):
        if not vals:
            return np.nan, np.nan
        return (
            float(np.quantile(vals, .025)),
            float(np.quantile(vals, .975)),
        )

    lo10, hi10 = ci(vals10)
    lo21, hi21 = ci(vals21)

    return {
        "delta10_CI_low": lo10,
        "delta10_CI_high": hi10,
        "delta21_CI_low": lo21,
        "delta21_CI_high": hi21,
        "bootstrap_replicates": len(vals10),
    }


def hierarchical_macro_bootstrap(
    tables: Dict[str, pd.DataFrame],
    n_boot: int = N_BOOT,
) -> Dict[str, float]:
    datasets = sorted(tables)
    rng = np.random.default_rng(SEED + 999)

    macro10 = []
    macro21 = []

    for _ in range(n_boot):
        chosen = rng.choice(datasets, size=len(datasets), replace=True)
        d10 = []
        d21 = []

        for ds in chosen:
            z = tables[ds]
            n = len(z)
            idx = rng.integers(0, n, size=n)
            b = z.iloc[idx]

            if np.var(b["TE_z"].to_numpy(float)) < 1e-10:
                continue

            d = compute_deltas(b)
            d10.append(d["delta_R2_M1_minus_M0"])
            d21.append(d["delta_R2_M2_minus_M1"])

        if d10:
            macro10.append(float(np.mean(d10)))
            macro21.append(float(np.mean(d21)))

    return {
        "macro_delta10_CI_low": float(np.quantile(macro10, .025)),
        "macro_delta10_CI_high": float(np.quantile(macro10, .975)),
        "macro_delta21_CI_low": float(np.quantile(macro21, .025)),
        "macro_delta21_CI_high": float(np.quantile(macro21, .975)),
        "bootstrap_replicates": len(macro10),
    }


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    print("=" * 104)
    print("TRACE PHASE 3B — FINAL SEQUENCE BASELINE M0/M1/M2")
    print("=" * 104)

    for p in [GTF, FUNCTIONAL]:
        if not p.exists():
            raise RuntimeError(f"Missing required input: {p}")

    functional = pd.read_csv(FUNCTIONAL, sep="\t", low_memory=False)

    required = {
        "gene_name", "transcript_id", "transcript_length",
        "TE_dataset", "TE_z",
        "naive_singleU_union_count", "breadth_sum_singleU",
    }
    missing = sorted(required - set(functional.columns))
    if missing:
        raise RuntimeError(
            f"Frozen Phase3A functional cohort missing: {missing}"
        )

    if functional["TE_dataset"].nunique() != 5:
        raise RuntimeError(
            f"Expected 5 HeLa TE datasets, found {functional['TE_dataset'].nunique()}."
        )

    gene_base = (
        functional[
            [
                "gene_name", "transcript_id", "transcript_length",
                "naive_singleU_union_count", "breadth_sum_singleU",
            ]
        ]
        .drop_duplicates("gene_name")
        .copy()
    )

    # Provenance is nested on top of unique burden.
    gene_base["confirmation_excess"] = (
        pd.to_numeric(gene_base["breadth_sum_singleU"], errors="coerce")
        - pd.to_numeric(gene_base["naive_singleU_union_count"], errors="coerce")
    )
    if (gene_base["confirmation_excess"] < 0).any():
        raise RuntimeError(
            "breadth_sum_singleU < naive_singleU_union_count for one or more genes."
        )

    gene_base["psi_naive_log"] = np.log1p(
        pd.to_numeric(gene_base["naive_singleU_union_count"], errors="coerce")
    )
    gene_base["psi_confirmation_log"] = np.log1p(
        pd.to_numeric(gene_base["confirmation_excess"], errors="coerce")
    )

    seq, seq_qc = extract_sequence_features(gene_base)

    gene = gene_base.merge(
        seq,
        on=["gene_name", "transcript_id"],
        how="inner",
        validate="one_to_one",
    )

    # Same sequence-eligible gene universe for all models.
    analysis = functional.merge(
        gene.drop(columns=["transcript_length"], errors="ignore"),
        on=["gene_name", "transcript_id"],
        how="inner",
        validate="many_to_one",
        suffixes=("", "_gene"),
    )

    analysis = analysis[
        analysis["TE_z"].notna()
        & analysis["psi_naive_log"].notna()
        & analysis["psi_confirmation_log"].notna()
    ].copy()

    if analysis["gene_name"].nunique() < 700:
        raise RuntimeError(
            f"Only {analysis['gene_name'].nunique()} sequence-eligible functional genes."
        )

    # Global fold map shared across every dataset/model.
    fold_map = assign_global_folds(analysis["gene_name"].unique())
    analysis["outer_fold"] = analysis["gene_name"].map(fold_map).astype(int)

    fold_table = (
        analysis[["gene_name", "outer_fold"]]
        .drop_duplicates("gene_name")
        .sort_values("gene_name")
    )
    fold_table.to_csv(
        AIM3 / "aim3b_global_gene_folds.tsv",
        sep="\t", index=False
    )

    # Frozen sequence feature list.
    seq_features = sorted([
        c for c in analysis.columns if c.startswith("seq_")
    ])

    # Remove globally constant sequence columns; this is outcome-independent.
    variable_seq = []
    constant_seq = []
    gene_feature_frame = (
        analysis[["gene_name"] + seq_features]
        .drop_duplicates("gene_name")
    )

    for c in seq_features:
        if pd.to_numeric(
            gene_feature_frame[c], errors="coerce"
        ).nunique(dropna=True) > 1:
            variable_seq.append(c)
        else:
            constant_seq.append(c)

    seq_features = variable_seq

    feature_contract = {
        "M0_sequence": seq_features,
        "M1_sequence_plus_naivePsi": seq_features + ["psi_naive_log"],
        "M2_plus_confirmation": (
            seq_features + ["psi_naive_log", "psi_confirmation_log"]
        ),
    }

    (META / "phase3b_feature_contract.json").write_text(
        json.dumps(
            {
                "sequence_features_used": seq_features,
                "constant_sequence_features_removed": constant_seq,
                "model_features": feature_contract,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    pred_all = []
    tuning_all = []

    for ds, d0 in analysis.groupby("TE_dataset", sort=True):
        d = d0.copy()

        # Verify common fold coverage.
        fold_counts = d["outer_fold"].value_counts().sort_index().to_dict()
        if len(fold_counts) != N_OUTER or min(fold_counts.values()) < 25:
            raise RuntimeError(
                f"{ds}: inadequate shared fold distribution {fold_counts}"
            )

        for model_name, features in feature_contract.items():
            pred, tuning = nested_oof(d, features, model_name)
            pred_all.append(pred)
            tuning["TE_dataset"] = ds
            tuning_all.append(tuning)

    pred_long = pd.concat(pred_all, ignore_index=True)
    tuning = pd.concat(tuning_all, ignore_index=True)

    pred_long.to_csv(
        AIM3 / "aim3b_oof_predictions.tsv",
        sep="\t", index=False
    )
    tuning.to_csv(
        AIM3 / "aim3b_nested_cv_tuning.tsv",
        sep="\t", index=False
    )

    perf = dataset_performance(pred_long)
    perf.to_csv(
        AIM3 / "aim3b_model_performance.tsv",
        sep="\t", index=False
    )

    delta_rows = []
    paired_tables = {}

    for ds in sorted(analysis["TE_dataset"].unique()):
        z = paired_prediction_table(pred_long, ds)
        paired_tables[ds] = z

        d = compute_deltas(z)
        b = bootstrap_dataset_deltas(z)

        delta_rows.append({
            "TE_dataset": ds,
            "n_genes": len(z),
            **d,
            **b,
        })

    deltas = pd.DataFrame(delta_rows)
    deltas.to_csv(
        AIM3 / "aim3b_dataset_deltas.tsv",
        sep="\t", index=False
    )

    macro = {
        "macro_R2_M0": float(deltas["R2_M0"].mean()),
        "macro_R2_M1": float(deltas["R2_M1"].mean()),
        "macro_R2_M2": float(deltas["R2_M2"].mean()),
        "macro_delta_R2_M1_minus_M0": float(
            deltas["delta_R2_M1_minus_M0"].mean()
        ),
        "macro_delta_R2_M2_minus_M1": float(
            deltas["delta_R2_M2_minus_M1"].mean()
        ),
        "datasets_delta10_positive": int(
            (deltas["delta_R2_M1_minus_M0"] > 0).sum()
        ),
        "datasets_delta21_positive": int(
            (deltas["delta_R2_M2_minus_M1"] > 0).sum()
        ),
    }
    macro.update(hierarchical_macro_bootstrap(paired_tables))

    psi_adds = bool(
        macro["macro_delta_R2_M1_minus_M0"] > 0
        and macro["macro_delta10_CI_low"] > 0
        and macro["datasets_delta10_positive"] >= 3
    )
    provenance_adds = bool(
        macro["macro_delta_R2_M2_minus_M1"] > 0
        and macro["macro_delta21_CI_low"] > 0
        and macro["datasets_delta21_positive"] >= 3
    )

    if psi_adds and provenance_adds:
        status = "PHASE3B_PSI_AND_PROVENANCE_ADD_SEQUENCE_INDEPENDENT_TE_INFORMATION"
    elif psi_adds:
        status = "PHASE3B_PSI_ADDS_INFO_PROVENANCE_DOES_NOT_ROBUSTLY_ADD"
    else:
        status = "PHASE3B_NO_ROBUST_SEQUENCE_INDEPENDENT_PSI_INCREMENT"

    contract = {
        "phase": "3B_final",
        "scientific_story": "Calibration -> Measurement -> Evidence -> Inference",
        "primary_question": (
            "Does Ψ burden add held-out TE predictive information beyond sequence, "
            "and does cross-technology confirmation add further information beyond "
            "total Ψ burden?"
        ),
        "TE_policy": (
            "Five HeLa TE datasets analyzed independently because Phase3A showed "
            "strong study-context heterogeneity for several maps."
        ),
        "fold_policy": (
            "One deterministic five-fold gene assignment shared across every "
            "TE dataset and M0/M1/M2."
        ),
        "model": (
            "Ridge regression with StandardScaler; alpha selected by inner CV "
            "using only outer-training genes."
        ),
        "M0": "MANE sequence-derived features only",
        "M1": "M0 + log1p(unique BID/BACS/ELAP single-U union burden)",
        "M2": (
            "M1 + log1p(extra assay confirmations = breadth_sum - unique burden)"
        ),
        "primary_metric": "out-of-fold R2",
        "primary_increments": [
            "delta_R2_M1_minus_M0",
            "delta_R2_M2_minus_M1",
        ],
        "macro_summary": (
            "unweighted mean across the five TE datasets with hierarchical "
            "study+gene bootstrap CI; per-dataset results remain primary context."
        ),
        "positive_increment_rule": (
            "macro delta > 0, hierarchical 95% CI lower bound > 0, "
            "and positive delta in at least 3/5 TE datasets"
        ),
        "interpretation_boundary": (
            "Predictive increment is sequence-independent information conditional "
            "on this classical sequence baseline; it is not causal evidence that "
            "pseudouridylation changes TE."
        ),
        "no_DRS_in_M1_M2": True,
        "no_raw_sequencing": True,
        "final_primary_analysis": True,
    }

    (META / "phase3b_analysis_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2),
        encoding="utf-8"
    )

    summary = {
        "status": status,
        "sequence_QC": seq_qc,
        "sequence_features_used_n": len(seq_features),
        "sequence_constant_features_removed_n": len(constant_seq),
        "functional_genes_sequence_eligible": int(
            analysis["gene_name"].nunique()
        ),
        "TE_datasets": sorted(analysis["TE_dataset"].unique().tolist()),
        "model_performance": perf.to_dict(orient="records"),
        "dataset_deltas": deltas.to_dict(orient="records"),
        "macro_summary": macro,
        "psi_beyond_sequence_supported": psi_adds,
        "provenance_beyond_naive_Psi_supported": provenance_adds,
        "final_interpretation_gate": (
            "This is the final primary analysis. Regardless of direction, "
            "freeze the main study after interpreting M0/M1/M2 together with "
            "Aim1, Aim2 and Phase3A."
        ),
    }

    (META / "phase3b_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8"
    )

    print()
    print("STATUS:", status)
    print("Sequence-eligible genes:", analysis["gene_name"].nunique())
    print("Sequence features:", len(seq_features))
    print()
    print("DATASET DELTAS")
    print(deltas.to_string(index=False))
    print()
    print("MACRO SUMMARY")
    print(json.dumps(macro, indent=2))
    print()
    print(r"Return/upload:")
    print(r"  D:\RNA\Trace\00_meta\phase3b_summary.json")
    print(r"  D:\RNA\Trace\00_meta\phase3b_analysis_contract.json")
    print(r"  D:\RNA\Trace\00_meta\phase3b_feature_contract.json")
    print(r"  D:\RNA\Trace\06_aim3\aim3b_model_performance.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3b_dataset_deltas.tsv")
    print(r"  D:\RNA\Trace\06_aim3\aim3b_oof_predictions.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        err = traceback.format_exc()
        (LOGS / "phase3b_FATAL_ERROR.txt").write_text(
            err, encoding="utf-8"
        )
        print("\nFATAL ERROR\n")
        print(err)
        print(r"Send D:\RNA\Trace\logs\phase3b_FATAL_ERROR.txt back to ChatGPT.")
        try:
            input("\nPress Enter to close...")
        except Exception:
            pass
        raise
