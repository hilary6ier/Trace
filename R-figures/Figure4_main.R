#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", argv[grepl("^--file=", argv)])
OUT_DIR <- dirname(normalizePath(script_arg))
ROOT <- normalizePath(file.path(OUT_DIR, ".."))
source(file.path(OUT_DIR, "_figure_common.R"))

metrics <- read_tsv(file.path(ROOT, "07_tracepsi", "phase4c_A549_lockbox_metrics.tsv"))
scorecard <- read_tsv(file.path(ROOT, "07_tracepsi", "phase4c_A549_TRACEpsi_scorecard.tsv.gz"))
robustness <- read_tsv(file.path(ROOT, "07_tracepsi", "phase4d_external_motif_robustness.tsv"))
primary <- metrics[score == "P2_TRACEpsi_frozen"]

y <- scorecard$DRS_reported
s <- scorecard$TRACEpsi_score
roc <- roc_points(y, s)
pr <- pr_points(y, s)
budget <- tie_budget_curve(y, s)
prevalence <- mean(y)

write_source(roc, "Figure4_panelA_A549_ROC.tsv")
write_source(pr, "Figure4_panelB_A549_precision_recall.tsv")
write_source(budget[, .(k, expected_positives, fraction_screened, fraction_recovered,
                        expected_precision, expected_lift, boundary_score)],
             "Figure4_panelC-D_A549_ranking_utility.tsv")
write_source(primary, "Figure4_frozen_lockbox_metrics.tsv")
write_source(robustness, "Figure4_supplementary_motif_robustness_provenance.tsv")

# Panel A: individual-score ROC.
pA <- ggplot(roc, aes(FPR, TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = COL["muted"], linewidth = 0.45) +
  geom_step(colour = COL["TRACE"], linewidth = 1.0, direction = "hv") +
  annotate("label", x = 0.96, y = 0.08,
           label = sprintf("AUROC %.3f\n95%% CI %.3f-%.3f", primary$AUROC,
                           primary$AUROC_CI_low, primary$AUROC_CI_high),
           hjust = 1, vjust = 0, family = FONT_FAMILY, size = 2.55,
           linewidth = 0, fill = alpha("white", 0.9), colour = COL["ink"]) +
  coord_equal() +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Frozen A549 lockbox discrimination", subtitle = "761 BID loci; 42 independently re-observed",
       x = "False-positive rate", y = "True-positive rate", tag = "A") +
  theme_manuscript() +
  theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

# Panel B: individual-score precision-recall.
pB <- ggplot(pr, aes(recall, precision)) +
  geom_hline(yintercept = prevalence, linetype = "dashed", colour = COL["muted"], linewidth = 0.5) +
  geom_step(colour = COL["residual"], linewidth = 1.0, direction = "hv") +
  annotate("text", x = 0.98, y = 0.12,
           label = sprintf("Prevalence %.1f%%", 100 * prevalence), hjust = 1,
           family = FONT_FAMILY, size = 2.35, colour = COL["muted"]) +
  annotate("label", x = 0.96, y = 0.82,
           label = sprintf("Average precision %.3f", primary$average_precision),
           hjust = 1, family = FONT_FAMILY, size = 2.5,
           linewidth = 0, fill = alpha("white", 0.9), colour = COL["ink"]) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Precision-recall above prevalence", x = "Recall", y = "Precision", tag = "B") +
  theme_manuscript() +
  theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

# Panel C: cumulative recovery/gain from tie-aware score groups.
pC <- ggplot(budget, aes(fraction_screened, fraction_recovered)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = COL["muted"], linewidth = 0.45) +
  geom_line(colour = COL["TRACE"], linewidth = 1.0) +
  geom_point(data = budget[k %in% c(50, 100, 200)], shape = 21, fill = "white",
             colour = COL["TRACE"], stroke = 0.7, size = 2.4) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Cumulative recovery is front-loaded", x = "Fraction of loci screened",
       y = "Fraction of DRS+ loci recovered", tag = "C") +
  theme_manuscript() +
  theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

# Panel D: exact-budget expected hits with tie handling.
budget_limit <- 250
pd <- budget[k <= budget_limit]
pd[, random_expected := k * prevalence]
frozen_budget <- data.table(
  k = c(50L, 100L, 200L),
  expected_positives = c(primary$budget50_expected_positives,
                         primary$budget100_expected_positives,
                         primary$budget200_expected_positives)
)

pD <- ggplot(pd, aes(k, expected_positives)) +
  geom_line(aes(y = random_expected), linetype = "dashed", colour = COL["muted"], linewidth = 0.55) +
  geom_line(colour = COL["TRACE"], linewidth = 1.0) +
  geom_point(data = frozen_budget, shape = 21, fill = COL["TRACE"], colour = "white",
             stroke = 0.3, size = 2.8) +
  geom_text(data = frozen_budget,
            aes(label = sprintf("K=%d: %.1f hits", k, expected_positives)),
            nudge_y = c(1.7, 1.7, 1.7), family = FONT_FAMILY, size = 2.35, colour = COL["ink"]) +
  annotate("text", x = 244, y = 244 * prevalence + 0.8, label = "Random ranking",
           hjust = 1, family = FONT_FAMILY, size = 2.3, colour = COL["muted"]) +
  scale_x_continuous(limits = c(0, budget_limit), breaks = seq(0, 250, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 31), breaks = seq(0, 30, 5), expand = c(0, 0)) +
  labs(title = "Practical validation-budget yield", x = "Loci selected for validation",
       y = "Expected DRS+ hits", tag = "D") +
  theme_manuscript() +
  theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

design <- "
AAABB
AAACC
DDDDD
"
figure <- wrap_plots(A = pA, B = pB, C = pC, D = pD, design = design) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_figure(figure, "Figure4_main", height_in = 6.4)
