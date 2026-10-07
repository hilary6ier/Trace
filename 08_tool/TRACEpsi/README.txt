TRACE-Ψ v1.1
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
