# TRACE-Ψ main-figure architecture lock

Date locked: 2026-10-07

## Scientific boundary

TRACE-Ψ estimates the portability of cross-platform evidence. It does not estimate biological truth or the probability that a site is truly pseudouridylated. DRS re-observation is treated as independent-platform support, and DRS non-report is not treated as biological absence.

The main-figure evidence chain is fixed as:

1. technology-associated sequence structure persists across cellular contexts;
2. support across distinct technologies stratifies independent-platform re-observation;
3. TRACE-Ψ resolves ranking within breadth=1 single-source evidence;
4. the frozen model transfers to an untouched A549 BID-to-DRS lockbox and improves practical validation yield.

Aim 3B generic TE prediction is excluded. Aim 3A is reserved for Supplementary material.

## Figure 1 — transferable measurement structure

- **A, cross-cell motif concordance:** observed five-mer effects learned independently in HEK293T and HeLa; only motifs observed in both training sets are displayed. Point area reflects the smaller motif count across cells. Pearson and Spearman correlations are descriptive.
- **B, synthetic calibration decomposition:** human-map motif log odds against the frozen synthetic calibration contrast, faceted by training cell. Point area reflects motif count; the line is the frozen weighted decomposition fit. This panel shows that calibration is an anchor but explains only a small component of the human-map structure.
- **C, cross-cell transfer:** full motif score and calibration-residual score AUROC with frozen motif-bootstrap 95% CIs for HEK293T-to-HeLa and HeLa-to-HEK293T transfer.
- **D, robustness:** primary, gene-disjoint, ELAP higher/highest confidence, and PUS7-like-excluded AUROC estimates with frozen 95% CIs in both transfer directions.

Primary sources: `04_aim1/aim1b_motif_effects.tsv`, `aim1b_transfer_results.tsv`, `aim1b_calibration_decomposition.tsv`, `aim1c_gene_disjoint.tsv`, `aim1d_corrected_confidence_sensitivity.tsv`, and `aim1c_writer_motif_sensitivity.tsv`.

## Figure 2 — independent-platform re-observation

- **A, evidence landscape:** each source support pattern is shown as BID/BACS/ELAP membership, DRS-positive count over total, and DRS re-observation proportion with Wilson 95% CI.
- **B, pooled breadth:** breadth 1/2/3 DRS re-observation proportions with Wilson 95% CIs.
- **C, breadth effect:** motif-cluster-robust odds ratio per additional source chemistry with 95% CI.
- **D, motif-conditioned null:** observed difference in mean breadth between DRS-positive and DRS-nonreported loci against a deterministic 20,000-draw within-five-mer permutation null. The displayed P value is the frozen result.

Primary sources: `05_aim2/aim2b_DRS_pattern_rates.tsv`, `aim2b_DRS_breadth_rates.tsv`, `aim2b_DRS_breadth_model.tsv`, `aim2b_DRS_motif_stratified_permutation.tsv`, and `aim2b_HeLa_source_union_DRS.tsv`.

## Figure 3 — ranking within single-source evidence

- **A, held-assay transfer:** HeLa leave-one-technology-out AUROC matrix for P0 breadth, P1 provenance, context-only, and TRACE-Ψ across held BACS, BID, and ELAP targets.
- **B, breadth=1 DRS ranking:** frozen HeLa breadth=1 AUROC and bootstrap 95% CI for the four benchmark models.
- **C, paired improvement:** paired bootstrap delta AUROC with 95% CI for the predeclared breadth=1 comparisons, emphasizing TRACE-Ψ versus context-only.

Primary sources: `07_tracepsi/phase4b_v2_internal_metrics.tsv`, `phase4b_v2_HeLa_DRS_metrics.tsv`, and `phase4b_v2_HeLa_DRS_breadth1_deltaAUC.tsv`. HEK cross-cell metrics remain part of provenance/QC but are not repeated in the main panels.

## Figure 4 — frozen external utility

- **A, ROC:** empirical, tie-aware ROC from individual frozen TRACE-Ψ scores in the untouched A549 BID-to-DRS lockbox.
- **B, precision-recall:** empirical, tie-aware PR curve with the 5.52% lockbox prevalence baseline.
- **C, cumulative recovery:** fraction of all DRS-positive loci recovered against fraction of loci screened, with random-ranking diagonal.
- **D, validation budget:** tie-aware expected DRS-positive hits for exact budgets from 1 to all 761 loci, with random-ranking baseline and frozen K=50/100/200 summaries.

Primary sources: `07_tracepsi/phase4c_A549_TRACEpsi_scorecard.tsv.gz` and `phase4c_A549_lockbox_metrics.tsv`. Context-only and external PUS7 robustness are kept out of the main figure.

## Shared visual and statistical grammar

- Vector master width: 183 mm; white background; sans-serif font; panel letters A–D.
- Restrained, colour-blind-aware palette shared across figures. TRACE-Ψ is blue; context-only is amber; provenance is mauve; breadth is grey; BID/BACS/ELAP use fixed blue/orange/green assay colours.
- Proportions use points with Wilson intervals. Model/effect comparisons use estimates with CIs. Paired comparisons use delta estimates with CIs. ROC/PR and ranking utility derive from individual scorecard rows.
- No summary bars, 3D, radar, rainbow scales, or workflow-box panels. Reference lines are neutral grey and explicitly interpretable.
- All plotting tables are exported to `R-figures/plot_source/`; every plotted number derives from frozen repository outputs or a documented deterministic transformation of them.

