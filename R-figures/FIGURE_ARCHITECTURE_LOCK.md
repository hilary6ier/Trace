# TRACE-PSI main-figure architecture lock

Status: **REVISED AND LOCKED before final rendering**

Governing workflow: repo-local `nature-figure` v2.8.0

Backend: **R only** (`ggplot2`, `patchwork`, `svglite`, `cairo_pdf`, `ragg`)

Canonical output directory: `R-figures/`
Final width: 183 mm; white page; vector PDF/SVG masters; TIFF at 600 dpi.

## Scientific boundary

TRACE-PSI estimates the portability of positive evidence across profiling
technologies. It does not estimate biological truth or `P(true pseudouridine)`.
DRS re-observation is independent-platform support, not ground truth. A source
non-call means absent from the published positive-call table, not unmodified.

The four main figures form one escalating argument:

`measurement structure -> reproducibility -> single-source prioritization -> frozen external utility`

Aim 3B generic TE prediction is excluded. Aim 3A is excluded from the Main
Figures and may only be used in Supplementary material. No old figure image or
old plot-source export is an evidentiary input; every displayed value is rebuilt
from the frozen files listed below.

## Shared visual and statistical contract

- Typography: Liberation Sans (Arial-compatible open fallback), 6.2--8 pt at
  final size; bold lowercase panel letters; every rendered PDF glyph >=5 pt.
- Background: white; top/right spines removed; no decorative panel boxes or
  default grid lines.
- Model colours: P0 breadth `#9AA0A6`; P1 provenance `#8F6BAE`;
  Context-only `#2A9D8F`; TRACE-PSI `#245B8A`.
- Assay colours: BID orange `#E68632`; BACS blue `#2E86AB`; ELAP violet
  `#7B61A8`; held-out DRS teal `#2A9D8F`.
- Cell/direction colours: HEK293T blue `#277DA1`; HeLa rose `#D25577`.
- Ordered breadth colours: breadth 1 `#86BFD1`; breadth 2 `#2A9D8F`;
  breadth 3 `#D96C4B`.
- Reference/null geometry: `#7A7F85` and `#D6D9DC`, visually subordinate.
- Proportions use points with Wilson 95% CIs; model/effect estimates use points
  with their frozen CIs; paired comparisons use delta estimates with paired
  bootstrap CIs; ROC and precision-recall use locus-level frozen scores;
  ranking utility is derived from the frozen lockbox scorecard with explicit
  tie handling.
- No bar chart is used for a summary estimate. No rainbow, radar, 3D or
  workflow-box panel is permitted.

## Figure 1 — transferable technology-associated sequence structure

**Results-level question.** Does the BID-versus-ELAP sequence structure recur
across cells, and can it be reduced to synthetic calibration or obvious
confounders?

**Figure-level claim.** Technology-associated 5-mer structure transfers across
cellular contexts, remains after removing the synthetic-calibration-aligned
component, and survives gene, confidence and PUS7-like-context attacks.

**Archetype.** Validation envelope. Three-panel asymmetric layout
`AABBB / CCCCC`: two scatter-based structure tests above one full-width
interval envelope.

| Panel | Evidence role | Decisive comparison | Frozen source |
|---|---|---|---|
| a | Primary quantitative evidence | HEK293T versus HeLa motif log-ratios among 105 motifs observed in both cells; Pearson and Spearman correlations reported | `04_aim1/aim1b_motif_effects.tsv`; `00_meta/phase1b_summary.json` |
| b | Decomposition | Synthetic calibration contrast versus human-map motif log-ratio, separately by training cell; point area reflects motif count and weighted fits are descriptive | `04_aim1/aim1b_motif_effects.tsv`; `04_aim1/aim1b_calibration_decomposition.tsv` |
| c | Transfer and robustness envelope | Full fingerprint, calibration residual, gene-disjoint, corrected confidence-restricted and PUS7-like-excluded AUROCs in both cross-cell directions, all in one point/CI grammar | `04_aim1/aim1b_transfer_results.tsv`; `aim1c_gene_disjoint.tsv`; `aim1d_corrected_confidence_sensitivity.tsv`; `aim1c_writer_motif_sensitivity.tsv` |

Panel a is the primary structure evidence; panel c is the full-width validation
envelope. Panel c uses a 0.5 reference line and never implies causality. The
invalid Phase 1C rank-based confidence analysis is excluded; only the corrected
author-defined ELAP confidence result is used.

## Figure 2 — cross-technology evidence predicts DRS re-observation

**Results-level question.** Among source-positive HeLa loci, does support across
distinct chemistries stratify re-observation by held-out DRS?

**Figure-level claim.** DRS re-observation increases strongly and monotonically
with chemistry breadth, is visible within exact source patterns, and survives
motif-conditioned label randomization.

**Archetype.** Evidence landscape with control. Three-panel asymmetric layout
`AAB / AAC`, with the pattern-level table-forest spanning both rows.

| Panel | Evidence role | Decisive comparison | Frozen source |
|---|---|---|---|
| a | Hero evidence landscape | Source-pattern assay membership, DRS+/n and Wilson 95% CI in one table-forest geometry | `05_aim2/aim2b_DRS_pattern_rates.tsv` |
| b | Pooled stratification and effect estimate | Breadth 1, 2 and 3 DRS re-observation proportions with Wilson 95% CIs, plus the motif-cluster-robust OR per additional chemistry | `05_aim2/aim2b_DRS_breadth_rates.tsv`; `05_aim2/aim2b_DRS_breadth_model.tsv` |
| c | Motif-conditioned falsification | Observed mean-breadth contrast versus a 20,000-permutation motif-stratified null; frozen P value retained | `05_aim2/aim2b_DRS_motif_stratified_permutation.tsv`; labels/strata from `aim2b_HeLa_source_union_DRS.tsv` |

The breadth anchors are 14.42%, 46.93% and 73.47%; the frozen OR is 4.614
(95% CI 3.546--6.004). Panel c uses a compact null strip/interval, not a sparse
density plot.

## Figure 3 — resolving breadth=1 prioritization

**Results-level question.** Can TRACE-PSI rank portability when evidence breadth
is fixed at one and therefore contains no ranking information?

**Figure-level claim.** Benchmark evidence shows that TRACE-PSI adds ranking
information within breadth=1 calls and improves over both provenance and
context-only comparators.

**Archetype.** Capability ladder focused on benchmark design. Asymmetric layout
`AAB / AAC`, with panel a spanning both rows.

| Panel | Evidence role | Decisive comparison | Frozen source |
|---|---|---|---|
| a | Benchmark landscape | HeLa leave-one-technology-out AUROC matrix for held-out BID, BACS and ELAP across P0, P1, Context-only and TRACE-PSI | `07_tracepsi/phase4b_v2_internal_metrics.tsv` |
| b | Primary challenge | HeLa DRS breadth=1 AUROC and bootstrap 95% CI for the four-model hierarchy | `07_tracepsi/phase4b_v2_HeLa_DRS_metrics.tsv` |
| c | Paired inference | Frozen paired delta AUROC comparisons, emphasizing TRACE-PSI minus Context-only | `07_tracepsi/phase4b_v2_HeLa_DRS_breadth1_deltaAUC.tsv` |

Panel b anchors: P0 0.500, P1 0.448, Context-only 0.606 and TRACE-PSI
0.644. Panel c anchors TRACE-PSI minus Context-only at +0.0380 (95% CI
+0.0084 to +0.0680). No coefficient panel is allowed.

## Figure 4 — untouched A549 lockbox utility

**Results-level question.** Does the frozen model transfer to a new cell and
native RNA platform, and what validation yield does its ranking provide?

**Figure-level claim.** Without retraining, recalibration or feature changes,
TRACE-PSI ranks A549 BID calls for DRS re-observation and concentrates expected
validation yield at realistic budgets.

**Archetype.** Frozen external validation plus decision utility. Asymmetric
layout `AAB / CCD`, with the recovery panel given the broad lower-left region.

| Panel | Evidence role | Decisive comparison | Frozen source |
|---|---|---|---|
| a | Discrimination | ROC from individual A549 lockbox scores, with diagonal reference and frozen AUROC/CI | `07_tracepsi/phase4c_A549_TRACEpsi_scorecard.tsv.gz`; `phase4c_A549_lockbox_metrics.tsv` |
| b | Rare-positive retrieval | Precision-recall curve from individual scores with the 5.52% prevalence baseline and frozen AP | same files as panel a |
| c | Ranking gain | Tie-aware expected cumulative DRS-positive recovery versus screened fraction, with random-ranking reference and top-10% anchor | frozen scorecard |
| d | Practical budget | Expected hits at K=50, 100 and 200 versus random expectation; K=50 and K=100 are emphasized | frozen scorecard; `phase4c_A549_lockbox_metrics.tsv` |

Lockbox anchors: n=761, DRS+=42, prevalence=5.519%, AUROC=0.7635, AP=0.2898,
top-10% precision=27.85%, Lift@10%=5.046, expected hits=18.324 at K=50 and
25.0 at K=100. Context-only and PUS7 robustness remain outside the Main Figure;
their frozen provenance is documented but not plotted.

## Mandatory render and QA gate

For every figure and after every layout-affecting revision:

1. run `validate_figure.py` on the R source;
2. render SVG, vector PDF, 600-dpi TIFF and an R-generated PNG preview at 183 mm;
3. preserve the patchwork alignment manifest/report and enforce the 1.5-pt gate;
4. run `audit_pdf_text.py --min-pt 5`;
5. run `audit_figure_collisions.py`, inspect every WARN/FAIL and retain JSON;
6. inspect each R-generated preview panel at final physical size for clipping,
   glyph integrity, balance, hidden encoding and statistical geometry;
7. regenerate all panel source TSVs from the frozen inputs and record checksums
   and derivations in `FIGURE_QC_PROVENANCE.md`.

No figure is final until these gates pass or a narrow, documented visual false
positive is justified after direct inspection.
