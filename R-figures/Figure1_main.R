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

motif <- read_tsv("04_aim1", "aim1b_motif_effects.tsv")
transfer <- read_tsv("04_aim1", "aim1b_transfer_results.tsv")
decomp <- read_tsv("04_aim1", "aim1b_calibration_decomposition.tsv")
gene <- read_tsv("04_aim1", "aim1c_gene_disjoint.tsv")
confidence <- read_tsv("04_aim1", "aim1d_corrected_confidence_sensitivity.tsv")
writer <- read_tsv("04_aim1", "aim1c_writer_motif_sensitivity.tsv")

# Panel a: only motifs observed in both cells are eligible for concordance.
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
spearman_rho <- cor(panel_a$HEK293T_log_ratio, panel_a$HeLa_log_ratio, method = "spearman")
write_plot_source(panel_a, "Figure1_panelA_motif_concordance.tsv")

p_a <- ggplot(panel_a, aes(HEK293T_log_ratio, HeLa_log_ratio)) +
  geom_abline(slope = 1, intercept = 0, colour = TRACE_COLORS[["neutral_light"]],
              linewidth = 0.45, linetype = "22") +
  geom_smooth(method = "lm", se = FALSE, colour = TRACE_COLORS[["ink"]],
              linewidth = 0.65, formula = y ~ x) +
  geom_point(aes(size = pooled_n), shape = 21, fill = TRACE_COLORS[["HEK293T"]],
             colour = "white", stroke = 0.25, alpha = 0.78) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("n = 105 motifs\nPearson r = %.2f\nSpearman rho = %.2f",
                           pearson_r, spearman_rho),
           hjust = -0.05, vjust = 1.12, size = 2.15, family = FONT_FAMILY,
           colour = TRACE_COLORS[["ink"]], lineheight = 1.12) +
  scale_size_continuous(range = c(1.4, 4.2), guide = "none") +
  labs(title = "Cross-cell motif concordance",
       x = "HEK293T motif log-ratio", y = "HeLa motif log-ratio") +
  coord_cartesian(clip = "off") +
  theme_trace()

# Panel b: synthetic calibration explains only a small weighted fraction.
panel_b <- motif[n_total > 0,
  .(motif, train_cell, test_cell, human_motif_log_ratio,
    calibration_contrast, n_total)]
fits <- panel_b[, {
  fit <- lm(human_motif_log_ratio ~ calibration_contrast, weights = n_total)
  .(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]),
    weighted_R2 = summary(fit)$r.squared, n_motifs = .N)
}, by = .(train_cell, test_cell)]
fits <- merge(fits, decomp[, .(train_cell, frozen_weighted_R2 = weighted_R2_training_motif_score)],
              by = "train_cell", all.x = TRUE)
stopifnot(max(abs(fits$weighted_R2 - fits$frozen_weighted_R2)) < 1e-10)
write_plot_source(panel_b, "Figure1_panelB_calibration_structure.tsv")
write_plot_source(fits, "Figure1_panelB_weighted_fits.tsv")

ann_b_label <- sprintf(
  "Weighted R2: HEK293T %.2f; HeLa %.2f",
  fits[train_cell == "HEK293T", weighted_R2],
  fits[train_cell == "HeLa", weighted_R2]
)
p_b <- ggplot(panel_b,
              aes(calibration_contrast, human_motif_log_ratio,
                  colour = train_cell)) +
  geom_hline(yintercept = 0, colour = TRACE_COLORS[["neutral_light"]], linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = TRACE_COLORS[["neutral_light"]], linewidth = 0.35) +
  geom_smooth(aes(weight = n_total), method = "lm", se = FALSE,
              linewidth = 0.65, formula = y ~ x) +
  geom_point(aes(size = n_total), alpha = 0.62) +
  scale_colour_manual(values = TRACE_COLORS[c("HEK293T", "HeLa")]) +
  scale_size_continuous(range = c(0.8, 3.3), guide = "none") +
  labs(title = "Synthetic calibration decomposition",
       subtitle = ann_b_label,
       x = "Synthetic calibration contrast", y = "Human-map motif log-ratio") +
  theme_trace() +
  theme(legend.position = "top")

# Panel c: full and calibration-residual transfer in both directions.
panel_c <- rbindlist(list(
  transfer[, .(train_cell, test_cell, component = "Full fingerprint",
               AUROC = AUC_cross_cell_5mer,
               CI_low = AUC_CI_low_motif_boot,
               CI_high = AUC_CI_high_motif_boot)],
  transfer[, .(train_cell, test_cell, component = "Calibration residual",
               AUROC = residual_fingerprint_AUC_test,
               CI_low = residual_AUC_CI_low_motif_boot,
               CI_high = residual_AUC_CI_high_motif_boot)]
))
panel_c[, direction := paste0(train_cell, " to ", test_cell)]
panel_c[, row_label := paste(direction, component, sep = "  |  ")]
panel_c[, row_label := factor(row_label, levels = rev(c(
  "HEK293T to HeLa  |  Full fingerprint",
  "HEK293T to HeLa  |  Calibration residual",
  "HeLa to HEK293T  |  Full fingerprint",
  "HeLa to HEK293T  |  Calibration residual"
)))]
write_plot_source(panel_c, "Figure1_panelC_crosscell_AUROC.tsv")

p_c <- ggplot(panel_c, aes(AUROC, row_label, colour = test_cell)) +
  geom_vline(xintercept = 0.5, colour = TRACE_COLORS[["neutral_mid"]],
             linetype = "22", linewidth = 0.42) +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = row_label), linewidth = 0.75) +
  geom_point(size = 2.35) +
  scale_colour_manual(values = TRACE_COLORS[c("HEK293T", "HeLa")], guide = "none") +
  scale_x_continuous(limits = c(0.48, 0.94), breaks = c(0.5, 0.7, 0.9)) +
  labs(title = "Cross-cell transfer", x = "AUROC (95% CI)", y = NULL) +
  theme_trace() +
  theme(axis.text.y = element_text(size = 5.55))

# Panel d: compact robustness envelope using only valid frozen analyses.
primary <- transfer[, .(train_cell, test_cell, condition = "Primary",
                        AUROC = AUC_cross_cell_5mer,
                        CI_low = AUC_CI_low_motif_boot,
                        CI_high = AUC_CI_high_motif_boot)]
gene_d <- gene[, .(train_cell, test_cell, condition = "Gene-disjoint",
                   AUROC = AUC, CI_low, CI_high)]
confidence_d <- confidence[analysis == "ELAP_higher_or_highest_only",
  .(train_cell, test_cell, condition = "Confidence-restricted",
    AUROC = AUC, CI_low, CI_high)]
writer_d <- writer[analysis == "remove_PUS7_like_UNUAR",
  .(train_cell, test_cell, condition = "PUS7-like excluded",
    AUROC = AUC, CI_low, CI_high)]
panel_d <- rbindlist(list(primary, gene_d, confidence_d, writer_d), use.names = TRUE)
condition_levels <- c("Primary", "Gene-disjoint", "Confidence-restricted", "PUS7-like excluded")
panel_d[, condition_index := 5 - match(condition, condition_levels)]
panel_d[, direction := paste0(train_cell, " to ", test_cell)]
panel_d[, y := condition_index + ifelse(test_cell == "HeLa", 0.11, -0.11)]
write_plot_source(panel_d, "Figure1_panelD_robustness.tsv")

p_d <- ggplot(panel_d, aes(AUROC, y, colour = test_cell)) +
  geom_vline(xintercept = 0.5, colour = TRACE_COLORS[["neutral_mid"]],
             linetype = "22", linewidth = 0.42) +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = y), linewidth = 0.72) +
  geom_point(aes(shape = test_cell), size = 2.25, stroke = 0.7) +
  scale_colour_manual(values = TRACE_COLORS[c("HEK293T", "HeLa")],
                      labels = c(HEK293T = "HeLa to HEK293T", HeLa = "HEK293T to HeLa")) +
  scale_shape_manual(values = c(HEK293T = 17, HeLa = 16),
                     labels = c(HEK293T = "HeLa to HEK293T", HeLa = "HEK293T to HeLa")) +
  scale_y_continuous(breaks = 4:1, labels = condition_levels,
                     limits = c(0.55, 4.45)) +
  scale_x_continuous(limits = c(0.48, 0.95), breaks = c(0.5, 0.7, 0.9)) +
  labs(title = "Robustness envelope", x = "AUROC (95% CI)", y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.2)), shape = "none") +
  theme_trace() +
  theme(legend.position = "top", legend.text = element_text(size = 5.2),
        axis.text.y = element_text(size = 5.55))

design <- "
AAC
BBD
"
figure <- p_a + p_b + p_c + p_d +
  plot_layout(design = design, heights = c(1.05, 1), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT_FAMILY))

save_trace_figure(
  figure, "Figure1_main", height_mm = 137,
  panel_ids = c("a", "c", "b", "d")
)
