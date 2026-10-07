#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
SCRIPT_DIR <- dirname(normalizePath(script_file, mustWork = TRUE))
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "_figure_common.R"))

# Static preflight contract (implemented by save_trace_figure): family =
# "Liberation Sans"; width_mm = 183; svglite::svglite();
# grDevices::cairo_pdf(); ragg::agg_tiff(res = 600);
# ragg::agg_png(res = 300); require_patchwork_panel_alignment().

pattern <- read_tsv("05_aim2", "aim2b_DRS_pattern_rates.tsv")
breadth <- read_tsv("05_aim2", "aim2b_DRS_breadth_rates.tsv")
model <- read_tsv("05_aim2", "aim2b_DRS_breadth_model.tsv")
perm_frozen <- read_tsv("05_aim2", "aim2b_DRS_motif_stratified_permutation.tsv")
union <- read_tsv("05_aim2", "aim2b_HeLa_source_union_DRS.tsv")

stopifnot(
  all.equal(
    breadth$DRS_support_rate,
    c(0.1441523118766999, 0.4692737430167598, 0.7346938775510204),
    tolerance = 1e-12
  ) == TRUE,
  abs(model[term == "chemistry_breadth", OR] - 4.614066769509874) < 1e-12
)

ASSAY_COLORS <- c(BID = "#E68632", BACS = "#2E86AB", ELAP = "#7B61A8")
BREADTH_COLORS <- c(`1` = "#86BFD1", `2` = "#2A9D8F", `3` = "#D96C4B")

# a | Compact membership table and pattern-specific Wilson intervals.
pattern_order <- c(
  "BID", "BACS", "ELAP", "BID+BACS", "BID+ELAP",
  "BACS+ELAP", "BID+BACS+ELAP"
)
pattern[, y := 8 - match(support_pattern, pattern_order)]
pattern[, breadth := lengths(strsplit(support_pattern, "+", fixed = TRUE))]
pattern[, fraction_label := sprintf("%d/%d", DRS_supported_n, n)]

assay_x <- c(BID = 0.075, BACS = 0.135, ELAP = 0.195)
fraction_x <- 0.295
forest_start <- 0.425
forest_end <- 0.965
forest_max <- 0.85
map_rate <- function(x) forest_start + (x / forest_max) * (forest_end - forest_start)

membership <- rbindlist(lapply(names(assay_x), function(a) {
  pattern[, .(
    support_pattern, y, assay = a, x = assay_x[[a]],
    present = grepl(a, support_pattern, fixed = TRUE)
  )]
}))
pattern[, `:=`(
  forest_x = map_rate(DRS_support_rate),
  forest_low = map_rate(wilson95_low),
  forest_high = map_rate(wilson95_high)
)]
write_plot_source(pattern, "Figure2_panelA_pattern_rates.tsv")
write_plot_source(membership, "Figure2_panelA_assay_membership.tsv")

header <- data.table(
  x = c(unname(assay_x), fraction_x, 0.695),
  y = 7.72,
  label = c(names(assay_x), "DRS+/n", "DRS re-observation")
)

p_a <- ggplot(pattern, aes(y = y)) +
  geom_hline(
    yintercept = seq(1.5, 6.5, by = 1),
    colour = "#EEF1F3", linewidth = 0.30
  ) +
  geom_point(
    data = membership,
    aes(x = x, y = y, fill = assay, alpha = present),
    inherit.aes = FALSE, shape = 21, size = 2.8,
    colour = "#AAB2B8", stroke = 0.32
  ) +
  geom_text(
    aes(x = fraction_x, label = fraction_label),
    size = 2.05, family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]
  ) +
  geom_segment(
    aes(x = forest_low, xend = forest_high, yend = y,
        colour = factor(breadth)),
    linewidth = 0.82
  ) +
  geom_point(
    aes(x = forest_x, fill = factor(breadth)),
    shape = 21, size = 2.85, colour = "white", stroke = 0.34
  ) +
  geom_text(
    data = header, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, size = 1.88, fontface = "bold",
    family = FONT_FAMILY, colour = TRACE_COLORS[["neutral_dark"]]
  ) +
  scale_fill_manual(values = c(ASSAY_COLORS, BREADTH_COLORS), guide = "none") +
  scale_colour_manual(values = BREADTH_COLORS, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
  scale_y_continuous(
    breaks = 7:1, labels = pattern_order,
    limits = c(0.55, 7.92), expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(0.02, 0.99),
    breaks = map_rate(c(0, 0.25, 0.50, 0.75)),
    labels = c("0%", "25%", "50%", "75%"),
    expand = c(0, 0)
  ) +
  labs(title = "Source patterns", tag = "a",
       x = "DRS re-observation (Wilson 95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_trace() +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.8),
    axis.ticks.x = element_line(colour = TRACE_COLORS[["ink"]]),
    plot.margin = margin(5, 6, 5, 5)
  )

# b | Breadth proportions and the adjusted effect estimate share one panel.
panel_b <- copy(breadth)
panel_b[, breadth_label := factor(chemistry_breadth, levels = c(1, 2, 3))]
panel_b[, percent_label := percent1(DRS_support_rate)]
panel_b[, colour_key := factor(chemistry_breadth)]
panel_b_or <- model[term == "chemistry_breadth", .(
  term, OR, OR_CI_low, OR_CI_high, p, n, n_motif_clusters
)]
stopifnot(panel_b_or$OR_CI_low > 0, panel_b_or$OR > 0, panel_b_or$OR_CI_high > 0)
write_plot_source(panel_b, "Figure2_panelB_breadth_rates.tsv")
write_plot_source(panel_b_or, "Figure2_panelB_breadth_OR.tsv")

or_label <- sprintf(
  "OR / chemistry = %.2f  [%.2f, %.2f]",
  panel_b_or$OR, panel_b_or$OR_CI_low, panel_b_or$OR_CI_high
)

p_b <- ggplot(panel_b, aes(DRS_support_rate, breadth_label)) +
  geom_segment(
    aes(x = wilson95_low, xend = wilson95_high, yend = breadth_label,
        colour = colour_key),
    linewidth = 0.92
  ) +
  geom_point(
    aes(fill = colour_key), shape = 21, size = 3.35,
    colour = "white", stroke = 0.38
  ) +
  geom_text(
    aes(x = wilson95_high + 0.035, label = percent_label),
    hjust = 0, size = 2.12, family = FONT_FAMILY,
    colour = TRACE_COLORS[["ink"]]
  ) +
  annotate(
    "text", x = 0.02, y = 3.42, label = or_label,
    hjust = 0, vjust = 0, size = 2.05, family = FONT_FAMILY,
    colour = TRACE_COLORS[["ink"]]
  ) +
  scale_colour_manual(values = BREADTH_COLORS, guide = "none") +
  scale_fill_manual(values = BREADTH_COLORS, guide = "none") +
  scale_x_continuous(
    limits = c(0, 1.06), breaks = c(0, 0.25, 0.50, 0.75),
    labels = percent_format(accuracy = 1), expand = c(0, 0)
  ) +
  scale_y_discrete(expand = expansion(add = c(0.30, 0.66))) +
  labs(title = "Evidence breadth", tag = "b",
       x = "DRS re-observation", y = "Chemistries") +
  theme_trace() +
  theme(
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 8, 3, 5)
  )

# c | Observed statistic against the motif-conditioned permutation null.
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
panel_c_null <- data.table(permutation_id = seq_len(n_perm), null_statistic = null_stat)
panel_c_frozen <- copy(perm_frozen)
q <- quantile(null_stat, c(0, 0.025, 0.5, 0.975, 1))
panel_c_summary <- data.table(
  null_min = q[[1]], null_q025 = q[[2]], null_median = q[[3]],
  null_q975 = q[[4]], null_max = q[[5]],
  observed = perm_frozen$observed[1], p_two_sided = perm_frozen$p_two_sided[1],
  permutations = n_perm
)
write_plot_source(panel_c_null, "Figure2_panelC_motif_conditioned_null.tsv")
write_plot_source(panel_c_frozen, "Figure2_panelC_frozen_permutation_result.tsv")
write_plot_source(panel_c_summary, "Figure2_panelC_null_summary.tsv")

show_idx <- unique(round(seq(1, n_perm, length.out = 360)))
show_null <- panel_c_null[show_idx]
show_null[, y := 0.83 + 0.028 * sin(seq_len(.N) * 1.7)]

p_c <- ggplot() +
  annotate(
    "rect", xmin = q[[2]], xmax = q[[4]], ymin = 0.93, ymax = 1.07,
    fill = "#DCE3E8", colour = NA
  ) +
  geom_point(
    data = show_null, aes(null_statistic, y),
    colour = "#AEB7BE", alpha = 0.28, size = 0.52
  ) +
  geom_segment(
    aes(x = q[[1]], xend = q[[5]], y = 1, yend = 1),
    colour = TRACE_COLORS[["neutral_mid"]], linewidth = 0.55
  ) +
  geom_segment(
    aes(x = q[[2]], xend = q[[4]], y = 1, yend = 1),
    colour = "#46515A", linewidth = 2.25, lineend = "round"
  ) +
  geom_point(
    aes(x = q[[3]], y = 1), shape = 21, size = 2.15,
    fill = "white", colour = "#46515A", stroke = 0.58
  ) +
  geom_point(
    aes(x = perm_frozen$observed[1], y = 1), shape = 23, size = 3.25,
    fill = "#D96C4B", colour = "white", stroke = 0.45
  ) +
  annotate(
    "text", x = q[[3]], y = 1.13, label = "Null 95%",
    hjust = 0.5, vjust = 0, size = 1.9, family = FONT_FAMILY,
    colour = TRACE_COLORS[["neutral_dark"]]
  ) +
  annotate(
    "text", x = perm_frozen$observed[1], y = 1.13,
    label = sprintf("Observed = %.3f\nP = %.1e", perm_frozen$observed[1],
                    perm_frozen$p_two_sided[1]),
    hjust = 1, vjust = 0, size = 2.0, lineheight = 1.05,
    family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]
  ) +
  scale_x_continuous(
    limits = c(min(q[[1]] - 0.025, -0.05), perm_frozen$observed[1] + 0.035),
    breaks = seq(0, 0.4, by = 0.2)
  ) +
  scale_y_continuous(limits = c(0.75, 1.31), breaks = NULL) +
  labs(title = "Motif-conditioned null", tag = "c",
       x = "Mean breadth difference", y = NULL) +
  theme_trace() +
  theme(plot.tag.position = c(-0.07, 1.04),
        plot.margin = margin(3, 5, 5, 5))

design <- "
AAB
AAC
"
figure <- p_a + p_b + p_c +
  plot_layout(design = design, heights = c(1.0, 1.0), guides = "keep") +
  plot_annotation() &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT_FAMILY))

save_trace_figure(
  figure, "Figure2_main", height_mm = 116,
  panel_ids = c("a", "b", "c")
)
