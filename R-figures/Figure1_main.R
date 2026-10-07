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

motif <- read_tsv("04_aim1", "aim1b_motif_effects.tsv")
transfer <- read_tsv("04_aim1", "aim1b_transfer_results.tsv")
decomp <- read_tsv("04_aim1", "aim1b_calibration_decomposition.tsv")
gene <- read_tsv("04_aim1", "aim1c_gene_disjoint.tsv")
confidence <- read_tsv("04_aim1", "aim1d_corrected_confidence_sensitivity.tsv")
writer <- read_tsv("04_aim1", "aim1c_writer_motif_sensitivity.tsv")

CELL_COLORS <- c(HEK293T = "#277DA1", HeLa = "#D25577")
PAIR_COLOR <- "#6C5AA7"

# a | Cross-cell concordance is the hero evidence.
hek <- motif[train_cell == "HEK293T" & n_total > 0,
             .(motif, HEK293T_log_ratio = human_motif_log_ratio,
               HEK293T_n = n_total)]
hela <- motif[train_cell == "HeLa" & n_total > 0,
              .(motif, HeLa_log_ratio = human_motif_log_ratio,
                HeLa_n = n_total)]
panel_a <- merge(hek, hela, by = "motif")
panel_a[, pooled_n := HEK293T_n + HeLa_n]
stopifnot(nrow(panel_a) == 105L)
pearson_r <- cor(panel_a$HEK293T_log_ratio, panel_a$HeLa_log_ratio)
spearman_rho <- cor(panel_a$HEK293T_log_ratio, panel_a$HeLa_log_ratio,
                    method = "spearman")
write_plot_source(panel_a, "Figure1_panelA_motif_concordance.tsv")

lim_a <- range(c(panel_a$HEK293T_log_ratio, panel_a$HeLa_log_ratio))
lim_a <- c(floor(lim_a[1] * 2) / 2, ceiling(lim_a[2] * 2) / 2)

p_a <- ggplot(panel_a, aes(HEK293T_log_ratio, HeLa_log_ratio)) +
  geom_abline(slope = 1, intercept = 0, colour = "#CBD2D8",
              linewidth = 0.45, linetype = "22") +
  geom_smooth(method = "lm", se = FALSE, colour = PAIR_COLOR,
              linewidth = 0.85, formula = y ~ x) +
  geom_point(shape = 21, size = 2.15, fill = "#4C9BC1", colour = "white",
             stroke = 0.28, alpha = 0.84) +
  annotate(
    "text", x = lim_a[1] + 0.10, y = lim_a[2] - 0.12,
    label = sprintf(
      "n == 105 ~~ italic(r) == %.2f ~~ rho == %.2f",
      pearson_r, spearman_rho
    ),
    parse = TRUE,
    hjust = 0, vjust = 1, size = 2.25, family = FONT_FAMILY,
    colour = TRACE_COLORS[["ink"]]
  ) +
  scale_x_continuous(limits = lim_a, breaks = seq(-4, 3, by = 2),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(limits = lim_a, breaks = seq(-4, 3, by = 2),
                     expand = expansion(mult = 0)) +
  coord_equal(clip = "off") +
  labs(title = "Motif concordance", x = "HEK293T log-ratio", y = "HeLa log-ratio") +
  theme_trace() +
  theme(plot.margin = margin(5, 7, 5, 5))

# b | Synthetic calibration explains little of the human-map motif structure.
panel_b <- motif[n_total > 0,
  .(motif, train_cell, test_cell, human_motif_log_ratio,
    calibration_contrast, n_total)]
fits <- panel_b[, {
  fit <- lm(human_motif_log_ratio ~ calibration_contrast, weights = n_total)
  .(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]),
    weighted_R2 = summary(fit)$r.squared, n_motifs = .N)
}, by = .(train_cell, test_cell)]
fits <- merge(
  fits,
  decomp[, .(train_cell, frozen_weighted_R2 = weighted_R2_training_motif_score)],
  by = "train_cell", all.x = TRUE
)
stopifnot(max(abs(fits$weighted_R2 - fits$frozen_weighted_R2)) < 1e-10)
fits[, r2_label := sprintf("%s  |  weighted R2 = %.2f", train_cell, weighted_R2)]
facet_labels <- setNames(fits$r2_label, fits$train_cell)
write_plot_source(panel_b, "Figure1_panelB_calibration_structure.tsv")
write_plot_source(fits, "Figure1_panelB_weighted_fits.tsv")

p_b <- ggplot(
  panel_b,
  aes(calibration_contrast, human_motif_log_ratio,
      fill = train_cell, colour = train_cell)
) +
  geom_hline(yintercept = 0, colour = "#D9DEE2", linewidth = 0.32) +
  geom_vline(xintercept = 0, colour = "#D9DEE2", linewidth = 0.32) +
  geom_smooth(aes(weight = n_total), method = "lm", se = FALSE,
              linewidth = 0.72, formula = y ~ x) +
  geom_point(aes(size = n_total), shape = 21, alpha = 0.60,
             colour = "white", stroke = 0.20) +
  facet_wrap(
    ~train_cell, nrow = 1,
    labeller = as_labeller(facet_labels)
  ) +
  scale_fill_manual(values = CELL_COLORS, guide = "none") +
  scale_colour_manual(values = CELL_COLORS, guide = "none") +
  scale_size_continuous(range = c(0.8, 2.8), guide = "none") +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_y_continuous(limits = c(-4.2, 3.55), breaks = c(-4, -2, 0, 2)) +
  labs(title = "Calibration decomposition", x = "Synthetic contrast",
       y = "Human-map log-ratio") +
  theme_trace() +
  theme(
    strip.text = element_text(size = 6.2, face = "bold", margin = margin(b = 1.5)),
    panel.spacing.x = grid::unit(4.5, "mm"),
    plot.margin = margin(5, 5, 3, 5)
  )

# c | One interval grammar combines the primary, residual and robustness tests.
full <- transfer[, .(
  train_cell, test_cell, condition = "Full",
  AUROC = AUC_cross_cell_5mer,
  CI_low = AUC_CI_low_motif_boot,
  CI_high = AUC_CI_high_motif_boot
)]
residual <- transfer[, .(
  train_cell, test_cell, condition = "Calibration-residual",
  AUROC = residual_fingerprint_AUC_test,
  CI_low = residual_AUC_CI_low_motif_boot,
  CI_high = residual_AUC_CI_high_motif_boot
)]
gene_d <- gene[, .(
  train_cell, test_cell, condition = "Gene-disjoint",
  AUROC = AUC, CI_low, CI_high
)]
confidence_d <- confidence[analysis == "ELAP_higher_or_highest_only", .(
  train_cell, test_cell, condition = "High-confidence",
  AUROC = AUC, CI_low, CI_high
)]
writer_d <- writer[analysis == "remove_PUS7_like_UNUAR", .(
  train_cell, test_cell, condition = "PUS7-like excluded",
  AUROC = AUC, CI_low, CI_high
)]
panel_c <- rbindlist(
  list(full, residual, gene_d, confidence_d, writer_d), use.names = TRUE
)
condition_levels <- c(
  "Full", "Calibration-residual", "Gene-disjoint",
  "High-confidence", "PUS7-like excluded"
)
panel_c[, direction := fifelse(
  train_cell == "HEK293T", "HEK293T to HeLa", "HeLa to HEK293T"
)]
panel_c[, condition_index := 6 - match(condition, condition_levels)]
panel_c[, y := condition_index + fifelse(direction == "HEK293T to HeLa", 0.12, -0.12)]
write_plot_source(panel_c, "Figure1_panelC_transfer_envelope.tsv")

p_c <- ggplot(panel_c, aes(AUROC, y, colour = direction)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.52, ymax = 5.48,
           fill = "#F3F0FA", colour = NA) +
  geom_vline(xintercept = 0.5, colour = TRACE_COLORS[["neutral_mid"]],
             linetype = "22", linewidth = 0.42) +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = y), linewidth = 0.78) +
  geom_point(shape = 21, aes(fill = direction), size = 2.55,
             colour = "white", stroke = 0.32) +
  scale_colour_manual(values = c(
    "HEK293T to HeLa" = CELL_COLORS[["HeLa"]],
    "HeLa to HEK293T" = CELL_COLORS[["HEK293T"]]
  )) +
  scale_fill_manual(values = c(
    "HEK293T to HeLa" = CELL_COLORS[["HeLa"]],
    "HeLa to HEK293T" = CELL_COLORS[["HEK293T"]]
  ), guide = "none") +
  scale_y_continuous(
    breaks = 5:1, labels = condition_levels,
    limits = c(0.52, 5.48), expand = c(0, 0)
  ) +
  scale_x_continuous(limits = c(0.48, 0.94), breaks = c(0.5, 0.7, 0.9)) +
  labs(title = "Transfer envelope", x = "AUROC (95% CI)", y = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.0, size = 2.4))) +
  theme_trace() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.text = element_text(size = 5.45),
    legend.key.width = grid::unit(4.0, "mm"),
    axis.text.y = element_text(size = 5.75),
    plot.margin = margin(3, 5, 5, 5)
  )

design <- "
AABBB
CCCCC
"
p_b_wrapped <- wrap_elements(full = p_b)
figure <- p_a + p_b_wrapped + p_c +
  plot_layout(design = design, heights = c(1.12, 0.88), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT_FAMILY))

save_trace_figure(
  figure, "Figure1_main", height_mm = 120,
  panel_ids = c("a", "b", "c")
)
