# TRACE-Ψ main-figure legends — draft

## Fig. 1 | Technology-associated sequence structure transfers across cellular contexts

**a,** Concordance of BID-versus-ELAP 5-mer log-ratios learned independently in
HEK293T and HeLa. Only motifs observed in both cells are shown (n = 105 motifs);
point area is proportional to the pooled number of source loci. The grey dashed
line denotes identity and the black line is the least-squares fit. Pearson
*r* = 0.509 and Spearman rho = 0.430. **b,** Relationship between the measured
synthetic-calibration contrast and the human-map motif log-ratio. Colours denote
the cell in which the human-map fingerprint was learned, point area denotes the
motif-specific number of loci, and lines are count-weighted fits. Synthetic
calibration explains only 4.9% and 8.1% of weighted motif-score variation in
HEK293T and HeLa, respectively. **c,** Cross-cell AUROC for the complete motif
fingerprint and the calibration-residual fingerprint. Horizontal lines are 95%
motif-cluster bootstrap confidence intervals. Test sets contain 492 HeLa loci
for HEK293T-to-HeLa transfer and 623 HEK293T loci for HeLa-to-HEK293T transfer.
**d,** Cross-cell AUROC under the primary analysis and gene-disjoint,
author-defined confidence-restricted and PUS7-like-motif-exclusion analyses.
Horizontal lines are 95% motif-cluster bootstrap confidence intervals; the
dashed vertical line in **c,d** denotes AUROC = 0.5. These results identify a
transferable technology-associated sequence structure but do not establish that
the structure is entirely assay-caused. Source data are provided in the
`R-figures/plot_source` directory.

## Fig. 2 | Chemistry-diverse support stratifies independent-platform re-observation

**a,** Evidence landscape for the 1,331-locus HeLa source-positive union. Filled
circles indicate membership in BID, BACS or ELAP; fractions give the number of
loci re-observed by held-out DRS over the number in each source pattern. Teal
points and horizontal lines show DRS re-observation proportions and Wilson 95%
confidence intervals. **b,** Pooled DRS re-observation by chemistry breadth:
159/1,103 (14.42%) for breadth 1, 84/179 (46.93%) for breadth 2 and 36/49
(73.47%) for breadth 3. Error bars are Wilson 95% confidence intervals.
**c,** Motif-cluster-robust logistic estimate for each additional source
chemistry (odds ratio 4.61, 95% CI 3.55–6.00; n = 1,331 loci in 137 5-mer
clusters; two-sided Wald *P* = 5.13 × 10^-30). **d,** Observed difference in
mean chemistry breadth between DRS-reobserved and DRS-not-reported loci
(diamond) compared with the motif-conditioned null. The thin and thick null
segments denote the full and central 95% permutation ranges, respectively; the
open circle denotes the null median and pale points show a deterministic subset
of permutations for visual context. Labels were permuted within 5-mer motif
strata (20,000 permutations; frozen two-sided *P* = 1.00 × 10^-4). DRS
non-report is not interpreted as biological absence. Source data are provided
in the `R-figures/plot_source` directory.

## Fig. 3 | TRACE-Ψ resolves portability within breadth-one evidence

**a,** HeLa leave-one-technology-out benchmark matrix. Each cell reports AUROC
when BID, BACS or ELAP is held out and predicted from the remaining positive
source evidence. Test sets contain 1,084 BID, 942 BACS and 864 ELAP target rows.
P0 uses evidence breadth, P1 adds source provenance, Context-only adds local
sequence context without the cross-cell fingerprint terms, and TRACE-Ψ (P2)
combines provenance, context and the frozen cross-cell fingerprint.
**b,** Performance on the HeLa DRS breadth=1 challenge (n = 1,103 loci,
159 DRS-reobserved). Points are AUROCs and horizontal lines are 95% bootstrap
confidence intervals from 3,000 locus-level resamples. Breadth alone has no
within-stratum ranking information (AUROC = 0.500), whereas TRACE-Ψ reaches
AUROC = 0.644 (95% CI 0.600–0.688). **c,** Paired differences in AUROC on the
same breadth=1 loci. Horizontal lines are paired-bootstrap 95% confidence
intervals from 3,000 resamples. TRACE-Ψ improves over Context-only by 0.0380
(95% CI 0.0084–0.0680; bootstrap support for delta > 0, 0.994). The benchmark
tests cross-platform evidence portability, not biological truth. Source data
are provided in the `R-figures/plot_source` directory.

## Fig. 4 | Frozen TRACE-Ψ transfers to an untouched A549 BID-to-DRS lockbox

**a,** Receiver-operating-characteristic curve calculated from individual
frozen TRACE-Ψ scores for 761 A549 BID loci, of which 42 were re-observed by
DRS. AUROC is 0.763 (95% bootstrap CI 0.677–0.846; 5,000 locus-level
resamples); the dashed line denotes random ranking. **b,** Precision-recall
curve from the same scores. Average precision is 0.290 and the dashed line is
the lockbox prevalence of 5.52%. **c,** Tie-aware expected cumulative recovery
of the 42 DRS-reobserved loci as progressively more A549 candidates are screened.
The dashed diagonal is random ranking. The nominal top 10% cutoff includes 79
loci after retaining the complete boundary tie and recovers 22/42 DRS-positive
loci (precision 27.8%; lift 5.05-fold). **d,** Expected DRS-positive hits at
fixed experimental budgets, treating loci tied at the boundary as exchangeable.
TRACE-Ψ yields 18.3, 25.0 and 25.5 expected hits at K = 50, 100 and 200,
respectively; grey points show random-ranking expectations at the observed
prevalence. The model, features and hyperparameters were frozen before the A549
DRS lockbox was opened; no retraining, recalibration or feature change was
performed. TRACE-Ψ estimates re-observation portability rather than
`P(true pseudouridine)`. Source data are provided in the
`R-figures/plot_source` directory.

