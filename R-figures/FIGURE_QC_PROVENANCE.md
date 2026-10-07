# TRACE-Ψ main figures — QC and provenance

## Governing workflow and scientific boundary

The figures were rebuilt under the repo-local `nature-figure` v2.8.0 workflow.
The complete skill router, core contract, stance, figure contract, multi-panel
evidence architecture, R workflow, design theory, legend conventions and QA
contract were read before implementation. R was the exclusive plotting,
rendering and export backend.

The project lock was enforced: TRACE-Ψ ranks cross-platform re-observation
portability and does not estimate biological truth; DRS is independent-platform
support rather than ground truth; a non-call is not a biological negative.
Aim 3B generic TE prediction and Aim 3A were excluded from all Main Figures.
Existing figure products were not used as data sources. Every plot-source file
was regenerated from the frozen project outputs by the corresponding R script.

## Runtime

- R 4.4.3
- ggplot2 4.0.3
- patchwork 1.3.2
- data.table 1.18.6.1
- svglite 2.2.2
- ragg 1.5.2
- jsonlite 2.0.0
- cowplot 1.2.0
- Font family: generic sans grobs, mapped to Liberation Sans for SVG; vector PDF
  exported with Cairo.
- Canonical width: 183 mm. Heights: Figure 1, 137 mm; Figures 2 and 4, 126 mm;
  Figure 3, 112 mm.

Reproduction from the repository root:

```bash
Rscript R-figures/Figure1_main.R
Rscript R-figures/Figure2_main.R
Rscript R-figures/Figure3_main.R
Rscript R-figures/Figure4_main.R
```

## Frozen-source checksums (SHA-256)

| Frozen input | SHA-256 |
|---|---|
| `00_meta/PROJECT_LOCK.md` | `7a08edb36a40800f93529627f61fb236b15c85857e4ffa5f53a78063a13c0036` |
| `04_aim1/aim1b_transfer_results.tsv` | `009b13cc6b58b5d7d85cc3f966e23e2b9a87e7de7b828825e2873e660703877c` |
| `04_aim1/aim1b_motif_effects.tsv` | `3a3481f295c46df99ca5608c3d52940ca8baaa2571aa49a21175bc7a367c4d33` |
| `04_aim1/aim1b_calibration_decomposition.tsv` | `108d7e4e4186a7a2195c04c5ff49e8b2a090aea2baa1133198b47ea98cc42013` |
| `04_aim1/aim1c_gene_disjoint.tsv` | `25a92760db58aced49fefd2412f94e20d7e7788cf1441966df6b78d4e4d6a0c8` |
| `04_aim1/aim1d_corrected_confidence_sensitivity.tsv` | `4b344c747377d018a24938af402acf3c3e793a16c66063ed5afe2fb1271ac0f8` |
| `04_aim1/aim1c_writer_motif_sensitivity.tsv` | `6ca27a08f862b27fcd2f945d9f3ee147f3b9b6724697d82cd29eb6cd87062fcc` |
| `05_aim2/aim2b_DRS_pattern_rates.tsv` | `76e0e5fe9c22fa5ec9a7ee2ce7e19c91624f5403970b8bccab06d0c7f425152c` |
| `05_aim2/aim2b_DRS_breadth_rates.tsv` | `fbf1f47fee68ef30ae8f665714232f2762ef3317fec94d0d607c99f60f123cb7` |
| `05_aim2/aim2b_DRS_breadth_model.tsv` | `88a4a9f62e7da5260b9dd06a6f80464bd676f01ec2ddb1f6e466e31d262e6b66` |
| `05_aim2/aim2b_DRS_motif_stratified_permutation.tsv` | `373f2281219b21abc95b40700095a6a19fbdbe8f55120e6cb14c4e0f5e2e79e4` |
| `05_aim2/aim2b_HeLa_source_union_DRS.tsv` | `b74a4a927703503888021dec3a33bfca5d138d9b7e6c3c1eae721e79059bd3c1` |
| `07_tracepsi/phase4b_v2_internal_metrics.tsv` | `7d2a1416d66016c761e6cd1324a0743d7b651d2b2c8fcf03ff8a931926d4d2bf` |
| `07_tracepsi/phase4b_v2_HeLa_DRS_metrics.tsv` | `4a92905e4f8759469c7a64816b51fba9b5bdf373c45d03c62a853dc3ccb23629` |
| `07_tracepsi/phase4b_v2_HeLa_DRS_breadth1_deltaAUC.tsv` | `34af70ba26eb8182bec95bf26feb1d5f7d6ef5c5094e0dce6b3d82586d2bf629` |
| `07_tracepsi/phase4c_A549_lockbox_metrics.tsv` | `50230634ed1736f53db229bad36b6ced0d0d88a825ff831b12301af5b9ef1cac` |
| `07_tracepsi/phase4c_A549_TRACEpsi_scorecard.tsv.gz` | `d7e53a5793782796ce1d34e16e40efb9a6f7930339818ac2112f46c602dbd823` |
| `07_tracepsi/phase4d_external_motif_robustness.tsv` | `155a303e7f3c901eff4858387d33addc98fdc33c5039a26f59328d6a237d5ce1` |

## Panel-level provenance and derivation

| Figure/panel | Derivation | Statistical encoding |
|---|---|---|
| 1a | Inner join of motifs with `n_total > 0` in both cell-specific motif tables; 105 motifs retained | Motif log-ratio scatter; least-squares fit; Pearson and Spearman correlations |
| 1b | Observed motifs from both cell-specific motif tables; count-weighted linear fits checked against frozen weighted R² | Scatter with motif-count area and cell colour |
| 1c | Full and calibration-residual rows reshaped directly from frozen transfer results | AUROC point and motif-cluster bootstrap 95% CI |
| 1d | Primary, gene-disjoint, corrected author-confidence and PUS7-like-exclusion rows concatenated without refitting | AUROC point and motif-cluster bootstrap 95% CI |
| 2a | Seven frozen source patterns plus membership parsed from the pattern string | Proportion point and Wilson 95% CI; DRS+/n printed explicitly |
| 2b | Frozen breadth 1/2/3 summary | Proportion point and Wilson 95% CI |
| 2c | Frozen motif-cluster-robust chemistry-breadth coefficient | OR point and 95% CI on an asserted-positive log axis |
| 2d | Frozen observed statistic and P value; 20,000 R permutations reproduce the specified within-5-mer label-randomization scheme from the frozen union table | Full/central-null interval, null median, deterministic 500-point display subset and observed diamond |
| 3a | Twelve rows from frozen internal metrics reshaped to a 4 × 3 matrix | AUROC heatmap with printed values |
| 3b | Frozen HeLa DRS rows restricted to `subset == breadth1` | AUROC point and 3,000-resample bootstrap 95% CI |
| 3c | Frozen paired delta-AUROC table | Delta point and paired-bootstrap 95% CI |
| 4a | Tied score groups accumulated from all 761 frozen A549 scorecard rows; trapezoidal AUROC asserted equal to the frozen metric | Individual-score ROC with random diagonal |
| 4b | Precision and recall accumulated at each unique frozen score; rebuilt average precision asserted equal to the frozen metric | Individual-score step PR curve and prevalence baseline |
| 4c | Expected positives computed for every exact K; a boundary tie contributes its positive fraction | Expected cumulative recovery and random diagonal |
| 4d | Tie-aware expected hits rebuilt at K=50/100/200 and asserted equal to frozen metrics | Point-line utility comparison with random expectation |

## Key numerical assertions

- Figure 1a: 105 common observed motifs; Pearson 0.5094503; Spearman 0.4301159.
- Figure 2b: 0.1441523, 0.4692737 and 0.7346939 for breadth 1/2/3.
- Figure 2c: OR 4.6140668, 95% CI 3.5459269–6.0039625.
- Figure 3b: breadth=1 AUROC P0 0.5000000, P1 0.4479666,
  Context-only 0.6062753 and TRACE-Ψ 0.6442743.
- Figure 3c: TRACE-Ψ minus Context-only delta AUROC 0.0379990, 95% CI
  0.0083880–0.0679765.
- Figure 4: n=761, DRS+=42, prevalence 0.0551905, rebuilt AUROC 0.7634943,
  rebuilt AP 0.2897929, top-10% precision 0.2784810, Lift@10% 5.0458107,
  expected hits 18.3243 at K=50 and 25.0000 at K=100.

## Final render QA

| Figure | Alignment gate | PDF glyph minimum | Collision audit | Physical output | Visual review |
|---|---|---:|---|---|---|
| 1 | PASS; 4 comparisons, 0 fail/warn | 5.199 pt, PASS | PASS; 0 fail/warn | 183 × 137 mm; TIFF 4,322 × 3,236 px | PASS after moving weighted-R² text out of the data field |
| 2 | PASS; 4 comparisons, 0 fail/warn | 5.263 pt, PASS | PASS; 0 fail/warn | 183 × 126 mm; TIFF 4,322 × 2,976 px | PASS after expanding breadth and OR label clearance |
| 3 | PASS; 3 comparisons, 0 fail/warn | 5.405 pt, PASS | 0 fail, 1 reviewed WARN | 183 × 112 mm; TIFF 4,322 × 2,645 px | PASS; WARN is the `0.5` colourbar tick touching its filled colourbar edge, an intentional readable colourbar relationship |
| 4 | PASS; 4 comparisons, 0 fail/warn | 5.300 pt, PASS | 0 fail, 5 reviewed WARN | 183 × 126 mm; TIFF 4,322 × 2,976 px | PASS; WARNs are axis/legend labels adjacent to their intended filled point or endpoint geometry; no glyph is obscured |

All PDFs are one-page vector masters. SVG text remains editable. Preview PNGs
were generated by R at 300 dpi and inspected at native resolution. Checks
confirmed no clipping, broken glyphs, unexplained visual encoding, excessive
blank space, inconsistent model/assay colour, hidden denominator change or
misleading uncertainty geometry.

Static preflight warnings were reviewed: R delimiter checks were independently
confirmed with `Rscript parse()` for all four scripts; Figure 2's log-scale
warning is covered by explicit positivity assertions for OR and both CI bounds.

## Repo-local skill QA compatibility fixes

Two narrow fixes were required for the installed skill under the current
ggplot2 4.0.3 / patchwork 1.3.2 / Cairo-PDF runtime:

1. `panel_alignment.R` now measures actual rendered panel viewports in device
   points when unresolved gtable `null` units otherwise produce zero-area
   rectangles. The 1.5-pt tolerance and strict blocking behaviour are unchanged.
2. `audit_figure_collisions.py` uses PyMuPDF's per-line span boxes instead of
   unioning broad Cairo text traces, and ignores page-spanning off-page Cairo
   construction paths. Text-text, on-page text-stroke, fill-edge and page-clipping
   checks remain active. The script's dependency-free self-test passes.

The measured alignment manifests/reports, diagnostic overlays and collision
JSON/PDF reports are stored in `R-figures/qa/`.

