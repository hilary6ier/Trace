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

pattern <- read_tsv("05_aim2", "aim2b_DRS_pattern_rates.tsv")
breadth <- read_tsv("05_aim2", "aim2b_DRS_breadth_rates.tsv")
model <- read_tsv("05_aim2", "aim2b_DRS_breadth_model.tsv")
perm_frozen <- read_tsv("05_aim2", "aim2b_DRS_motif_stratified_permutation.tsv")
union <- read_tsv("05_aim2", "aim2b_HeLa_source_union_DRS.tsv")

stopifnot(
  all.equal(breadth$DRS_support_rate,
            c(0.1441523118766999, 0.4692737430167598, 0.7346938775510204),
            tolerance = 1e-12) == TRUE,
  abs(model[term == "chemistry_breadth", OR] - 4.614066769509874) < 1e-12
)

# Panel a: one combined assay-membership table and Wilson-CI forest.
pattern_order <- c("BID", "BACS", "ELAP", "BID+BACS", "BID+ELAP",
                   "BACS+ELAP", "BID+BACS+ELAP")
pattern[, y := rev(seq_along(pattern_order))[match(support_pattern, pattern_order)]]
pattern[, fraction_label := sprintf("%d/%d", DRS_supported_n, n)]
assay_x <- c(BID = -0.52, BACS = -0.43, ELAP = -0.34)
membership <- rbindlist(lapply(names(assay_x), function(a) {
  pattern[, .(support_pattern, y, assay = a, x = assay_x[[a]],
              present = grepl(a, support_pattern, fixed = TRUE))]
}))
write_plot_source(pattern, "Figure2_panelA_pattern_rates.tsv")
write_plot_source(membership, "Figure2_panelA_assay_membership.tsv")

header <- data.table(x = unname(assay_x), y = 7.75, label = names(assay_x))
p_a <- ggplot(pattern, aes(DRS_support_rate, y)) +
  geom_segment(aes(x = wilson95_low, xend = wilson95_high, yend = y),
               colour = TRACE_COLORS[["ink"]], linewidth = 0.72) +
  geom_point(shape = 21, size = 2.55, fill = TRACE_COLORS[["DRS"]],
             colour = "white", stroke = 0.35) +
  geom_text(aes(x = -0.22, label = fraction_label), hjust = 0.5,
            size = 2.0, family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]) +
  geom_point(data = membership,
             aes(x = x, y = y, fill = assay, alpha = present),
             inherit.aes = FALSE, shape = 21, size = 2.6,
             colour = TRACE_COLORS[["neutral_mid"]], stroke = 0.35) +
  geom_text(data = header, aes(x, y, label = label), inherit.aes = FALSE,
            size = 1.85, fontface = "bold", family = FONT_FAMILY,
            colour = TRACE_COLORS[["neutral_dark"]]) +
  annotate("text", x = -0.22, y = 7.75, label = "DRS+/n", size = 1.85,
           fontface = "bold", family = FONT_FAMILY,
           colour = TRACE_COLORS[["neutral_dark"]]) +
  scale_fill_manual(values = TRACE_COLORS[c("BID", "BACS", "ELAP")], guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
  scale_y_continuous(breaks = rev(seq_along(pattern_order)), labels = pattern_order,
                     limits = c(0.55, 7.95), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(-0.57, 0.90), breaks = c(0, 0.25, 0.50, 0.75),
    labels = percent_format(accuracy = 1)
  ) +
  labs(title = "Source-pattern evidence landscape",
       x = "DRS re-observation (Wilson 95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_trace() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5.8))

# Panel b: pooled breadth proportions.
panel_b <- copy(breadth)
panel_b[, breadth_label := paste0("Breadth ", chemistry_breadth)]
panel_b[, percent_label := percent1(DRS_support_rate)]
write_plot_source(panel_b, "Figure2_panelB_breadth_rates.tsv")

p_b <- ggplot(panel_b, aes(chemistry_breadth, DRS_support_rate)) +
  geom_line(colour = TRACE_COLORS[["neutral_mid"]], linewidth = 0.55) +
  geom_errorbar(aes(ymin = wilson95_low, ymax = wilson95_high), width = 0.10,
                linewidth = 0.65, colour = TRACE_COLORS[["ink"]]) +
  geom_point(aes(fill = factor(chemistry_breadth)), shape = 21, size = 3.0,
             colour = "white", stroke = 0.35) +
  geom_text(aes(label = percent_label, y = wilson95_high + 0.075), size = 2.05,
            family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]) +
  scale_fill_manual(values = c(`1` = "#B8CDD9", `2` = "#6095AD", `3` = "#245B8A"),
                    guide = "none") +
  scale_x_continuous(limits = c(0.72, 3.18), breaks = 1:3,
                     labels = paste("Breadth", 1:3)) +
  scale_y_continuous(limits = c(0, 0.94), breaks = c(0, 0.25, 0.50, 0.75),
                     labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(title = "Pooled chemistry breadth", x = NULL, y = "DRS re-observation") +
  theme_trace() +
  theme(axis.text.x = element_text(size = 5.45))

# Panel c: motif-cluster robust odds ratio.
panel_c <- model[term == "chemistry_breadth",
                 .(term, OR, OR_CI_low, OR_CI_high, p, n, n_motif_clusters)]
stopifnot(panel_c$OR_CI_low > 0, panel_c$OR > 0, panel_c$OR_CI_high > 0)
panel_c[, label := sprintf("OR %.2f; 95%% CI %.2f-%.2f", OR, OR_CI_low, OR_CI_high)]
write_plot_source(panel_c, "Figure2_panelC_breadth_OR.tsv")

p_c <- ggplot(panel_c, aes(OR, 1)) +
  geom_vline(xintercept = 1, colour = TRACE_COLORS[["neutral_mid"]],
             linetype = "22", linewidth = 0.45) +
  geom_segment(aes(x = OR_CI_low, xend = OR_CI_high, yend = 1),
               linewidth = 0.9, colour = TRACE_COLORS[["ink"]]) +
  geom_point(shape = 21, size = 3.2, fill = TRACE_COLORS[["DRS"]],
             colour = "white", stroke = 0.4) +
  geom_text(aes(x = OR_CI_high * 1.03, label = sub("; ", "\n", label)),
            hjust = 0, vjust = 0.5,
            size = 2.1, lineheight = 1.1, family = FONT_FAMILY) +
  scale_x_log10(limits = c(0.8, 14), breaks = c(1, 2, 4, 8)) +
  scale_y_continuous(limits = c(0.72, 1.28), breaks = NULL) +
  labs(title = "Per-chemistry breadth effect",
       subtitle = "Motif-cluster robust; n = 1,331",
       x = "Odds ratio", y = NULL) +
  theme_trace()

# Panel d: simulate the frozen motif-conditioned randomization scheme in R.
# The frozen observed statistic and P value remain authoritative.
set.seed(20261003)
n_perm <- as.integer(perm_frozen$permutations[1])
y <- union$DRS_reported_support
b <- union$chemistry_breadth
n1 <- sum(y == 1)
n0 <- sum(y == 0)
total_b <- sum(b)
positive_b_sum <- numeric(n_perm)
for (idx in split(seq_len(nrow(union)), union$motif)) {
  k <- sum(y[idx] == 1)
  if (k == 0) next
  if (k == length(idx)) {
    positive_b_sum <- positive_b_sum + sum(b[idx])
  } else {
    draws <- replicate(n_perm, sum(sample(b[idx], size = k, replace = FALSE)))
    positive_b_sum <- positive_b_sum + draws
  }
}
null_stat <- positive_b_sum / n1 - (total_b - positive_b_sum) / n0
panel_d_null <- data.table(permutation_id = seq_len(n_perm), null_statistic = null_stat)
panel_d_frozen <- copy(perm_frozen)
q <- quantile(null_stat, c(0, 0.025, 0.5, 0.975, 1))
panel_d_summary <- data.table(
  null_min = q[[1]], null_q025 = q[[2]], null_median = q[[3]],
  null_q975 = q[[4]], null_max = q[[5]],
  observed = perm_frozen$observed[1], p_two_sided = perm_frozen$p_two_sided[1],
  permutations = n_perm
)
write_plot_source(panel_d_null, "Figure2_panelD_motif_conditioned_null.tsv")
write_plot_source(panel_d_frozen, "Figure2_panelD_frozen_permutation_result.tsv")
write_plot_source(panel_d_summary, "Figure2_panelD_null_summary.tsv")

show_idx <- unique(round(seq(1, n_perm, length.out = 500)))
show_null <- panel_d_null[show_idx]
show_null[, y := 0.84 + 0.035 * sin(seq_len(.N) * 1.7)]
p_d <- ggplot() +
  geom_point(data = show_null, aes(null_statistic, y),
             colour = TRACE_COLORS[["neutral_mid"]], alpha = 0.24, size = 0.55) +
  geom_segment(aes(x = q[[1]], xend = q[[5]], y = 1, yend = 1),
               colour = TRACE_COLORS[["neutral_mid"]], linewidth = 0.45) +
  geom_segment(aes(x = q[[2]], xend = q[[4]], y = 1, yend = 1),
               colour = TRACE_COLORS[["ink"]], linewidth = 2.1, lineend = "round") +
  geom_point(aes(x = q[[3]], y = 1), shape = 21, size = 2.1,
             fill = "white", colour = TRACE_COLORS[["ink"]], stroke = 0.55) +
  geom_point(aes(x = perm_frozen$observed[1], y = 1), shape = 23, size = 3.1,
             fill = TRACE_COLORS[["DRS"]], colour = "white", stroke = 0.45) +
  annotate("text", x = perm_frozen$observed[1], y = 1.14,
           label = sprintf("Observed %.3f\nP = %.1e", perm_frozen$observed[1],
                           perm_frozen$p_two_sided[1]),
           hjust = 1, vjust = 0, size = 2.0, lineheight = 1.05,
           family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]) +
  annotate("text", x = q[[3]], y = 1.09, label = "Motif-conditioned null",
           hjust = 0.5, vjust = 0, size = 1.9, family = FONT_FAMILY,
           colour = TRACE_COLORS[["neutral_dark"]]) +
  scale_x_continuous(limits = c(min(q[[1]] - 0.025, -0.05), perm_frozen$observed[1] + 0.035)) +
  scale_y_continuous(limits = c(0.76, 1.30), breaks = NULL) +
  labs(title = "Motif-conditioned permutation",
       x = "Mean breadth: DRS+ minus DRS not reported", y = NULL) +
  theme_trace() +
  theme(plot.title = element_text(margin = margin(l = 8, b = 2.5)))

design <- "
AAAB
CCDD
"
figure <- p_a + p_b + p_c + p_d +
  plot_layout(design = design, heights = c(1.45, 0.78), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT_FAMILY))

save_trace_figure(
  figure, "Figure2_main", height_mm = 126,
  panel_ids = c("a", "b", "c", "d")
)
