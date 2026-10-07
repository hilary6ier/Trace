#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
SCRIPT_DIR <- dirname(normalizePath(script_file, mustWork = TRUE))
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "_figure_common.R"))

# Static preflight contract (implemented by save_trace_figure): family = "Liberation Sans";
# width_mm = 183; svglite::svglite(); grDevices::cairo_pdf();
# ragg::agg_tiff(res = 600); ragg::agg_png(res = 300);
# require_patchwork_panel_alignment().

scorecard <- read_tsv_gz("07_tracepsi", "phase4c_A549_TRACEpsi_scorecard.tsv.gz")
metrics <- read_tsv("07_tracepsi", "phase4c_A549_lockbox_metrics.tsv")
delta <- read_tsv("07_tracepsi", "phase4c_A549_lockbox_deltaAUC.tsv")
robustness <- read_tsv("07_tracepsi", "phase4d_external_motif_robustness.tsv")

frozen <- metrics[score == "P2_TRACEpsi_frozen"]
stopifnot(nrow(scorecard) == frozen$n, sum(scorecard$DRS_reported) == frozen$positives)
write_plot_source(frozen, "Figure4_frozen_lockbox_metrics.tsv")
write_plot_source(delta, "Figure4_lockbox_deltaAUC_provenance.tsv")
write_plot_source(robustness, "Figure4_supplementary_motif_robustness_provenance.tsv")

y <- scorecard$DRS_reported
score <- scorecard$TRACEpsi_score
n <- length(y)
n_pos <- sum(y)
n_neg <- n - n_pos
prevalence <- mean(y)

# Tied scores are accumulated as groups, matching the frozen ranking semantics.
groups <- scorecard[, .(group_n = .N, group_pos = sum(DRS_reported)),
                    by = .(score = TRACEpsi_score)][order(-score)]
groups[, `:=`(cum_n = cumsum(group_n), cum_pos = cumsum(group_pos))]
groups[, cum_neg := cum_n - cum_pos]

# Panel a: ROC from individual lockbox scores.
panel_a <- rbind(
  data.table(threshold = Inf, false_positive_rate = 0, true_positive_rate = 0),
  groups[, .(threshold = score,
             false_positive_rate = cum_neg / n_neg,
             true_positive_rate = cum_pos / n_pos)]
)
auc_rebuilt <- sum(diff(panel_a$false_positive_rate) *
                     (head(panel_a$true_positive_rate, -1) + tail(panel_a$true_positive_rate, -1)) / 2)
stopifnot(abs(auc_rebuilt - frozen$AUROC) < 1e-12)
write_plot_source(panel_a, "Figure4_panelA_A549_ROC.tsv")

p_a <- ggplot(panel_a, aes(false_positive_rate, true_positive_rate)) +
  geom_abline(slope = 1, intercept = 0, colour = TRACE_COLORS[["neutral_mid"]],
              linewidth = 0.45, linetype = "22") +
  geom_step(colour = TRACE_COLORS[["P2_TRACEpsi"]], linewidth = 0.9, direction = "hv") +
  annotate("text", x = 0.97, y = 0.06,
           label = sprintf("AUROC %.3f\n95%% CI %.3f-%.3f",
                           frozen$AUROC, frozen$AUROC_CI_low, frozen$AUROC_CI_high),
           hjust = 1, vjust = 0, size = 2.15, lineheight = 1.1,
           family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(title = "Frozen A549 discrimination", x = "False-positive rate",
       y = "True-positive rate") +
  theme_trace()

# Panel b: precision-recall from the same scores.
panel_b <- rbind(
  data.table(threshold = Inf, recall = 0, precision = 1),
  groups[, .(threshold = score, recall = cum_pos / n_pos, precision = cum_pos / cum_n)]
)
ap_rebuilt <- sum(diff(panel_b$recall) * tail(panel_b$precision, -1))
stopifnot(abs(ap_rebuilt - frozen$average_precision) < 1e-12)
write_plot_source(panel_b, "Figure4_panelB_A549_precision_recall.tsv")

p_b <- ggplot(panel_b, aes(recall, precision)) +
  geom_hline(yintercept = prevalence, colour = TRACE_COLORS[["neutral_mid"]],
             linewidth = 0.45, linetype = "22") +
  geom_step(colour = TRACE_COLORS[["P2_TRACEpsi"]], linewidth = 0.9, direction = "hv") +
  annotate("text", x = 0.97, y = 0.92,
           label = sprintf("AP %.3f\nPrevalence %.1f%%",
                           frozen$average_precision, 100 * prevalence),
           hjust = 1, vjust = 1, size = 2.15, lineheight = 1.1,
           family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.03), breaks = c(0, 0.5, 1),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(title = "Rare-positive retrieval", x = "Recall", y = "Precision") +
  theme_trace()

# Tie-aware expected positives for every exact validation budget.
expected_hits_at_k <- function(k) {
  if (k <= 0) return(0)
  boundary <- which(groups$cum_n >= k)[1]
  before_n <- if (boundary == 1) 0 else groups$cum_n[boundary - 1]
  before_pos <- if (boundary == 1) 0 else groups$cum_pos[boundary - 1]
  need <- k - before_n
  before_pos + need * groups$group_pos[boundary] / groups$group_n[boundary]
}

panel_c <- data.table(budget = 0:n)
panel_c[, expected_hits := vapply(budget, expected_hits_at_k, numeric(1))]
panel_c[, `:=`(
  fraction_screened = budget / n,
  positive_recovery = expected_hits / n_pos,
  random_expected_hits = budget * prevalence,
  random_recovery = budget / n
)]
write_plot_source(panel_c, "Figure4_panelC_A549_cumulative_recovery.tsv")

k10_inclusive <- frozen$K10_tie_inclusive
hits10 <- frozen$Precision10_tieaware * k10_inclusive
p_c <- ggplot(panel_c, aes(fraction_screened, positive_recovery)) +
  geom_abline(slope = 1, intercept = 0, colour = TRACE_COLORS[["neutral_mid"]],
              linewidth = 0.45, linetype = "22") +
  geom_line(colour = TRACE_COLORS[["P2_TRACEpsi"]], linewidth = 0.95) +
  annotate("point", x = k10_inclusive / n, y = hits10 / n_pos,
           shape = 21, size = 3.0, fill = TRACE_COLORS[["P2_TRACEpsi"]],
           colour = "white", stroke = 0.4) +
  annotate("text", x = k10_inclusive / n + 0.025, y = hits10 / n_pos,
           label = sprintf("Top 10%%: %.0f/%d DRS+\nPrecision %.1f%%; lift %.2fx",
                           hits10, n_pos, 100 * frozen$Precision10_tieaware,
                           frozen$Lift10_tieaware),
           hjust = 0, vjust = 0.5, size = 2.05, lineheight = 1.08,
           family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(title = "Cumulative recovery", x = "A549 candidates screened",
       y = "DRS-positive loci recovered") +
  theme_trace()

# Panel d: practical validation budgets, with the random expectation exposed.
budget_values <- c(50L, 100L, 200L)
panel_d <- data.table(
  budget = budget_values,
  trace_expected_hits = c(frozen$budget50_expected_positives,
                          frozen$budget100_expected_positives,
                          frozen$budget200_expected_positives)
)
panel_d[, rebuilt_expected_hits := vapply(budget, expected_hits_at_k, numeric(1))]
stopifnot(max(abs(panel_d$trace_expected_hits - panel_d$rebuilt_expected_hits)) < 1e-12)
panel_d[, random_expected_hits := budget * prevalence]
panel_d_long <- melt(panel_d, id.vars = "budget",
                     measure.vars = c("trace_expected_hits", "random_expected_hits"),
                     variable.name = "ranking", value.name = "expected_hits")
panel_d_long[, ranking := factor(ranking,
  levels = c("random_expected_hits", "trace_expected_hits"),
  labels = c("Random ranking", "TRACE-PSI"))]
panel_d_labels <- panel_d_long[ranking == "TRACE-PSI"]
panel_d_labels[, `:=`(
  label_x = budget + c(8, 8, 5),
  label_y = expected_hits + c(-3.0, 3.0, 3.0)
)]
write_plot_source(panel_d, "Figure4_panelD_A549_validation_budget.tsv")

p_d <- ggplot(panel_d_long, aes(budget, expected_hits, colour = ranking, group = ranking)) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.45) +
  geom_text(data = panel_d_labels,
            aes(x = label_x, y = label_y, label = sprintf("%.1f", expected_hits)),
            inherit.aes = FALSE, hjust = 0.5, vjust = 0.5,
            size = 1.95, family = FONT_FAMILY, show.legend = FALSE) +
  scale_colour_manual(values = c("Random ranking" = TRACE_COLORS[["neutral_mid"]],
                                 "TRACE-PSI" = TRACE_COLORS[["P2_TRACEpsi"]])) +
  scale_x_continuous(breaks = budget_values, limits = c(42, 220)) +
  scale_y_continuous(limits = c(0, 33.5), breaks = c(0, 10, 20, 30)) +
  labs(title = "Validation-budget utility", x = "Experimental budget, K",
       y = "Expected DRS+ hits") +
  guides(colour = guide_legend(override.aes = list(linewidth = 0, size = 2.4))) +
  theme_trace() +
  theme(legend.position = "top", legend.text = element_text(size = 5.3),
        legend.key.width = grid::unit(3.2, "mm"),
        legend.spacing.x = grid::unit(1.8, "mm"),
        legend.margin = margin(t = 2, b = 2),
        axis.text.x = element_text(size = 5.6), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(margin = margin(b = 4)))

design <- "
AABB
CCCD
"
figure <- p_a + p_b + p_c + p_d +
  plot_layout(design = design, heights = c(1, 0.92), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT_FAMILY))

save_trace_figure(
  figure, "Figure4_main", height_mm = 126,
  panel_ids = c("a", "b", "c", "d")
)
