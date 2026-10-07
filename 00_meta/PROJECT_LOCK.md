# Trace project lock

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
