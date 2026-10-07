# Figure QC and provenance

## Frozen analysis state

- Repository commit audited: `005eca95501a1ababf7b5d4c780707e158bbb35f`.
- Scientific constraints were taken from `00_meta/PROJECT_LOCK.md` and the Phase 1B/1C/1D, 2B, 4B-v2, 4C, and 4D contracts/summaries.
- No source TSV, frozen model, scorecard, analysis script, dependency declaration, or lockfile was modified.
- Aim 3B generic TE prediction was not used. Aim 3A was not promoted to a main figure.

## Figure-to-source mapping

| Figure/panel | Frozen source | Plotting transformation |
|---|---|---|
| 1A | `04_aim1/aim1b_motif_effects.tsv` | Inner join of motifs with `n_total > 0` in both cell-specific training maps; point area is the smaller count. |
| 1B | `aim1b_motif_effects.tsv`, `aim1b_calibration_decomposition.tsv`, `aim1b_transfer_results.tsv` | Observed motifs only; frozen weighted-decomposition intercept and slope overlaid without refitting. |
| 1C | `aim1b_transfer_results.tsv` | Frozen full and calibration-residual AUROCs and motif-bootstrap CIs. |
| 1D | `aim1b_transfer_results.tsv`, `aim1c_gene_disjoint.tsv`, `aim1d_corrected_confidence_sensitivity.tsv`, `aim1c_writer_motif_sensitivity.tsv` | Selected predeclared primary, gene-disjoint, higher/highest-confidence, and PUS7-like-exclusion rows. |
| 2A | `05_aim2/aim2b_DRS_pattern_rates.tsv` | Assay membership parsed from support-pattern labels; frozen proportions and Wilson CIs used directly. |
| 2B | `aim2b_DRS_breadth_rates.tsv` | Frozen pooled breadth proportions and Wilson CIs used directly. |
| 2C | `aim2b_DRS_breadth_model.tsv` | Chemistry-breadth OR and motif-cluster-robust CI used directly on a log x-axis. |
| 2D | `aim2b_HeLa_source_union_DRS.tsv`, `aim2b_DRS_motif_stratified_permutation.tsv` | Observed statistic recomputed from source rows; a deterministic 20,000-draw within-motif null is generated for geometry. The displayed P value is the frozen result. |
| 3A | `07_tracepsi/phase4b_v2_internal_metrics.tsv` | Held-target AUROC matrix for all four frozen benchmark models. |
| 3B | `phase4b_v2_HeLa_DRS_metrics.tsv` | `subset == breadth1`; frozen AUROCs and bootstrap CIs. |
| 3C | `phase4b_v2_HeLa_DRS_breadth1_deltaAUC.tsv` | Frozen paired delta AUROCs and bootstrap CIs. |
| 4A–D | `phase4c_A549_TRACEpsi_scorecard.tsv.gz`, `phase4c_A549_lockbox_metrics.tsv` | ROC, PR, tie-aware cumulative recovery, and exact-budget expected hits derived from individual frozen scores and labels. |

Panel-level plotting tables are exported in `R-figures/plot_source/`. The external PUS7 robustness table is copied there as Supplementary provenance but is not plotted in the main Figure 4.

## Independent numerical checks

- Figure 1A retained 105 shared observed motifs. Recalculated Pearson r = 0.50945; the frozen Spearman rho = 0.43012.
- Figure 2 breadth rates were verified from their numerators and denominators: 159/1,103 = 14.415%, 84/179 = 46.927%, and 36/49 = 73.469%. The source-level breadth statistic recalculated to 0.4441208. The R display null had median 0.12215 and central 95% interval 0.06320–0.19017; the observed value lies above all displayed-null support.
- Figure 3 breadth=1 AUROCs were verified as 0.50000, 0.44797, 0.60627, and 0.64427 for P0, P1, context-only, and TRACE-Ψ. The TRACE-Ψ minus context-only delta was 0.037999 with 95% CI 0.008388–0.067977.
- Figure 4 was recomputed from all 761 scorecard rows using an independent Python check: 42 positives, prevalence 0.0551905, AUROC 0.7634943, and average precision 0.2897929. Tie-aware expected hits exactly matched the frozen metrics at K = 50 (18.3243), 100 (25.0000), and 200 (25.4774).

## Rendering and file QC

- R runtime: 4.4.3 in an isolated Pixi/conda-forge environment. Plotting uses `data.table`, `ggplot2`, `patchwork`, `scales`, and `ragg`.
- PDF masters use Cairo vector output. All PDFs are one page and 518 pt wide (approximately 182.7 mm). Embedded subset fonts resolve to Liberation Sans/Liberation Sans Bold.
- TIFF masters are LZW-compressed at 600 × 600 pixels/inch and 4,322 px wide. Heights are 4,200 px (Figure 1), 3,600 px (Figure 2), 2,880 px (Figure 3), and 3,840 px (Figure 4).
- PNG previews were rendered at 240 dpi and inspected as actual images after each material revision.
- Visual inspection confirmed complete panel labels, unclipped titles and axes, aligned text, visible CI endpoints, white backgrounds, consistent assay/model colours, and absence of broken symbols. Long panel titles were shortened after the first render; factor-colour mappings and panel-tag positions were corrected; Figure 2 was compressed to remove excess row spacing.
- Geometry checks: proportions use Wilson intervals; effect estimates and paired differences use CIs; no summary bars are used; ROC/PR curves use score-level rows; recovery/budget curves account for score ties; the permutation panel directly contrasts the observed statistic with the motif-conditioned null.

## Reproduction

Run each script from the repository root with an R installation containing the declared packages, for example:

```bash
Rscript R-figures/Figure1_main.R
Rscript R-figures/Figure2_main.R
Rscript R-figures/Figure3_main.R
Rscript R-figures/Figure4_main.R
```

Each script rewrites its plot-source TSVs and all three final formats (`.pdf`, `_600dpi.tiff`, and `_preview.png`) under `R-figures/`.

