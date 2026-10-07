#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", argv[grepl("^--file=", argv)])
OUT_DIR <- dirname(normalizePath(script_arg))
ROOT <- normalizePath(file.path(OUT_DIR, ".."))
source(file.path(OUT_DIR, "_figure_common.R"))

motif <- read_tsv(file.path(ROOT, "04_aim1", "aim1b_motif_effects.tsv"))
transfer <- read_tsv(file.path(ROOT, "04_aim1", "aim1b_transfer_results.tsv"))
decomp <- read_tsv(file.path(ROOT, "04_aim1", "aim1b_calibration_decomposition.tsv"))
gene <- read_tsv(file.path(ROOT, "04_aim1", "aim1c_gene_disjoint.tsv"))
confidence <- read_tsv(file.path(ROOT, "04_aim1", "aim1d_corrected_confidence_sensitivity.tsv"))
writer <- read_tsv(file.path(ROOT, "04_aim1", "aim1c_writer_motif_sensitivity.tsv"))

direction_label <- function(train_cell, test_cell) {
  ifelse(train_cell == "HEK293T", "HEK293T to HeLa", "HeLa to HEK293T")
}

# Panel A: motif effects observed in both cell-specific training maps.
hek <- motif[train_cell == "HEK293T" & n_total > 0,
             .(motif, HEK293T_log_odds = human_motif_log_ratio, n_HEK293T = n_total)]
hela <- motif[train_cell == "HeLa" & n_total > 0,
              .(motif, HeLa_log_odds = human_motif_log_ratio, n_HeLa = n_total)]
pa_data <- merge(hek, hela, by = "motif")
pa_data[, shared_count := pmin(n_HEK293T, n_HeLa)]
pearson_r <- cor(pa_data$HEK293T_log_odds, pa_data$HeLa_log_odds, method = "pearson")
spearman_r <- cor(pa_data$HEK293T_log_odds, pa_data$HeLa_log_odds, method = "spearman")
write_source(pa_data, "Figure1_panelA_motif_concordance.tsv")

pA <- ggplot(pa_data, aes(HEK293T_log_odds, HeLa_log_odds)) +
  geom_hline(yintercept = 0, colour = COL["grid"], linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = COL["grid"], linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = COL["muted"], linewidth = 0.4) +
  geom_point(aes(size = shared_count), shape = 21, fill = COL["TRACE"], colour = "white",
             stroke = 0.25, alpha = 0.78) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("n = %d motifs\nPearson r = %.2f\nSpearman rho = %.2f",
                           nrow(pa_data), pearson_r, spearman_r),
           hjust = -0.08, vjust = 1.12, family = FONT_FAMILY, size = 2.45,
           linewidth = 0, fill = alpha("white", 0.88), colour = COL["ink"]) +
  scale_size_area(max_size = 4.0, breaks = c(1, 5, 15), name = "Min. calls") +
  coord_equal() +
  labs(title = "Cross-cell motif concordance",
       x = "HEK293T motif effect\nlog odds (ELAP:BID)",
       y = "HeLa motif effect\nlog odds (ELAP:BID)") +
  theme_manuscript() +
  theme(legend.position = c(0.78, 0.17), legend.background = element_blank())

# Panel B: human motif structure versus the synthetic calibration anchor.
pb_data <- motif[n_total > 0,
                 .(motif, train_cell, test_cell, n_BID, n_ELAP, n_total,
                   human_motif_log_ratio, calibration_contrast)]
pb_data[, training_map := factor(train_cell, levels = c("HEK293T", "HeLa"),
                                 labels = c("HEK293T map", "HeLa map"))]
pb_fit <- merge(
  decomp[, .(train_cell, intercept = transfer[match(train_cell, transfer$train_cell), intercept],
             slope = calibration_slope, weighted_R2 = weighted_R2_training_motif_score)],
  unique(pb_data[, .(train_cell, training_map)]), by = "train_cell"
)
write_source(pb_data, "Figure1_panelB_calibration_structure.tsv")
write_source(pb_fit, "Figure1_panelB_weighted_fits.tsv")

pB <- ggplot(pb_data, aes(calibration_contrast, human_motif_log_ratio)) +
  geom_hline(yintercept = 0, colour = COL["grid"], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = COL["grid"], linewidth = 0.3) +
  geom_point(aes(size = n_total), shape = 21, fill = COL["BACS"], colour = "white",
             stroke = 0.2, alpha = 0.72) +
  geom_abline(data = pb_fit, aes(intercept = intercept, slope = slope),
              inherit.aes = FALSE, colour = COL["ink"], linewidth = 0.65) +
  geom_text(data = pb_fit,
            aes(x = -Inf, y = Inf, label = sprintf("weighted R2 = %.2f", weighted_R2)),
            inherit.aes = FALSE, hjust = -0.08, vjust = 1.25, family = FONT_FAMILY,
            size = 2.35, colour = COL["ink"]) +
  facet_wrap(~ training_map, nrow = 1) +
  scale_size_area(max_size = 3.2, guide = "none") +
  labs(title = "Calibration is a partial anchor",
       x = "Synthetic calibration contrast", y = "Human-map motif effect") +
  theme_manuscript() +
  theme(panel.spacing = unit(3, "mm"))

# Panel C: full and calibration-residual cross-cell discrimination.
pc_full <- transfer[, .(
  direction = direction_label(train_cell, test_cell),
  component = "Full motif score", estimate = AUC_cross_cell_5mer,
  low = AUC_CI_low_motif_boot, high = AUC_CI_high_motif_boot, n_test
)]
pc_resid <- transfer[, .(
  direction = direction_label(train_cell, test_cell),
  component = "Calibration residual", estimate = residual_fingerprint_AUC_test,
  low = residual_AUC_CI_low_motif_boot, high = residual_AUC_CI_high_motif_boot, n_test
)]
pc_data <- rbind(pc_full, pc_resid)
pc_data[, direction := factor(direction, levels = c("HeLa to HEK293T", "HEK293T to HeLa"))]
pc_data[, component := factor(component, levels = c("Full motif score", "Calibration residual"))]
write_source(pc_data, "Figure1_panelC_crosscell_AUROC.tsv")

pC <- ggplot(pc_data, aes(estimate, direction, colour = component)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = COL["muted"], linewidth = 0.4) +
  geom_errorbar(aes(xmin = low, xmax = high), position = position_dodge(width = 0.42),
                orientation = "y", width = 0, linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.42), size = 2.15) +
  scale_colour_manual(values = setNames(c(unname(COL["TRACE"]), unname(COL["residual"])),
                                        c("Full motif score", "Calibration residual"))) +
  scale_x_continuous(limits = c(0.48, 0.94), breaks = c(0.5, 0.7, 0.9)) +
  labs(title = "Residual structure transfers", x = "AUROC (95% CI)", y = NULL, colour = NULL) +
  theme_manuscript() +
  theme(legend.position = "top", legend.justification = "left")

# Panel D: compact robustness forest.
pd_primary <- transfer[, .(direction = direction_label(train_cell, test_cell),
                           analysis = "Primary", estimate = AUC_cross_cell_5mer,
                           low = AUC_CI_low_motif_boot, high = AUC_CI_high_motif_boot, n_test)]
pd_gene <- gene[, .(direction = direction_label(train_cell, test_cell),
                   analysis = "Gene-disjoint", estimate = AUC, low = CI_low, high = CI_high, n_test)]
pd_conf <- confidence[analysis == "ELAP_higher_or_highest_only",
                      .(direction = direction_label(train_cell, test_cell),
                        analysis = "Confidence-restricted", estimate = AUC,
                        low = CI_low, high = CI_high, n_test)]
pd_writer <- writer[analysis == "remove_PUS7_like_UNUAR",
                    .(direction = direction_label(train_cell, test_cell),
                      analysis = "Exclude PUS7-like", estimate = AUC,
                      low = CI_low, high = CI_high, n_test)]
pd_data <- rbindlist(list(pd_primary, pd_gene, pd_conf, pd_writer), use.names = TRUE)
pd_data[, analysis := factor(analysis,
  levels = rev(c("Primary", "Gene-disjoint", "Confidence-restricted", "Exclude PUS7-like")))]
pd_data[, direction := factor(direction, levels = c("HEK293T to HeLa", "HeLa to HEK293T"))]
write_source(pd_data, "Figure1_panelD_robustness.tsv")

pD <- ggplot(pd_data, aes(estimate, analysis, colour = direction)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = COL["muted"], linewidth = 0.4) +
  geom_errorbar(aes(xmin = low, xmax = high), position = position_dodge(width = 0.5),
                orientation = "y", width = 0, linewidth = 0.65) +
  geom_point(position = position_dodge(width = 0.5), size = 2.05) +
  scale_colour_manual(values = setNames(c(unname(COL["BID"]), unname(COL["ELAP"])),
                                        c("HEK293T to HeLa", "HeLa to HEK293T"))) +
  scale_x_continuous(limits = c(0.48, 0.96), breaks = c(0.5, 0.7, 0.9)) +
  labs(title = "Robust to locus, confidence and motif restrictions",
       x = "AUROC (95% CI)", y = NULL, colour = NULL) +
  theme_manuscript() +
  theme(legend.position = "top", legend.justification = "left")

design <- "
AAABBB
AAACCC
DDDDDD
"
figure <- wrap_plots(A = pA, B = pB, C = pC, D = pD, design = design,
                     heights = c(0.8, 0.8, 1.0)) +
  plot_annotation(tag_levels = "A", theme = panel_tag_theme) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_figure(figure, "Figure1_main", height_in = 7.0)
