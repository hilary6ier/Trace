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

internal <- read_tsv("07_tracepsi", "phase4b_v2_internal_metrics.tsv")
drs <- read_tsv("07_tracepsi", "phase4b_v2_HeLa_DRS_metrics.tsv")
delta <- read_tsv("07_tracepsi", "phase4b_v2_HeLa_DRS_breadth1_deltaAUC.tsv")

model_order <- c("P0_breadth", "P1_provenance", "A_context_only", "P2_TRACEpsi")
target_order <- c("BID", "BACS", "ELAP")

# Panel a: leave-one-technology-out performance matrix.
panel_a <- internal[, .(held_target, model, n, positives, prevalence,
                        AUROC, AUROC_CI_low, AUROC_CI_high, average_precision)]
panel_a[, model := factor(model, levels = rev(model_order))]
panel_a[, held_target := factor(held_target, levels = target_order)]
panel_a[, value_label := sprintf("%.3f", AUROC)]
write_plot_source(panel_a, "Figure3_panelA_HeLa_LOTO_matrix.tsv")

p_a <- ggplot(panel_a, aes(held_target, model, fill = AUROC)) +
  geom_tile(colour = "white", linewidth = 1.15) +
  geom_text(aes(label = value_label,
                colour = AUROC >= 0.68), size = 2.45, fontface = "bold",
            family = FONT_FAMILY) +
  scale_fill_gradientn(
    colours = c("#F3F5F6", "#C9DCE5", "#72A3B9", "#245B8A"),
    limits = c(0.50, 0.76), breaks = c(0.50, 0.60, 0.70),
    name = "AUROC"
  ) +
  scale_colour_manual(values = c(`FALSE` = TRACE_COLORS[["ink"]], `TRUE` = "white"),
                      guide = "none") +
  scale_y_discrete(labels = MODEL_LABELS[rev(model_order)]) +
  labs(title = "Held-technology benchmark",
       subtitle = "HeLa leave-one-technology-out AUROC",
       x = "Held-out assay", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_trace() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.key.width = grid::unit(14, "mm"),
        plot.margin = margin(5, 8, 5, 5))

# Panel b: breadth=1 DRS ranking challenge.
panel_b <- drs[subset == "breadth1",
               .(model, n, positives, prevalence, AUROC,
                 AUROC_CI_low, AUROC_CI_high, average_precision, Lift10, Precision10)]
panel_b[, model := factor(model, levels = rev(model_order))]
panel_b[, model_label := factor(MODEL_LABELS[as.character(model)],
                                levels = MODEL_LABELS[rev(model_order)])]
write_plot_source(panel_b, "Figure3_panelB_HeLa_breadth1_AUROC.tsv")

p_b <- ggplot(panel_b, aes(AUROC, model_label, colour = model)) +
  geom_vline(xintercept = 0.5, colour = TRACE_COLORS[["neutral_mid"]],
             linetype = "22", linewidth = 0.42) +
  geom_segment(aes(x = AUROC_CI_low, xend = AUROC_CI_high, yend = model_label),
               linewidth = 0.78) +
  geom_point(size = 2.55) +
  geom_text(aes(label = sprintf("%.3f", AUROC),
                x = ifelse(AUROC < 0.5, AUROC_CI_low - 0.008, AUROC_CI_high + 0.008),
                hjust = ifelse(AUROC < 0.5, 1, 0)),
            size = 1.9, family = FONT_FAMILY,
            colour = TRACE_COLORS[["ink"]]) +
  scale_colour_manual(values = TRACE_COLORS[model_order], guide = "none") +
  scale_x_continuous(limits = c(0.39, 0.72), breaks = c(0.4, 0.5, 0.6, 0.7)) +
  labs(title = "Breadth=1 DRS challenge",
       subtitle = "n = 1,103; DRS+ = 159",
       x = "AUROC (bootstrap 95% CI)", y = NULL) +
  theme_trace() +
  theme(axis.text.y = element_text(size = 5.55))

# Panel c: paired bootstrap differences on the same loci.
label_map <- c(
  "P1_provenance-P0_breadth" = "P1 - P0",
  "A_context_only-P1_provenance" = "Context - P1",
  "P2_TRACEpsi-P1_provenance" = "TRACE-PSI - P1",
  "P2_TRACEpsi-A_context_only" = "TRACE-PSI - Context"
)
panel_c <- copy(delta)
panel_c[, comparison_key := paste(model_A, model_B, sep = "-")]
panel_c[, comparison := label_map[comparison_key]]
comparison_order <- c("P1 - P0", "Context - P1", "TRACE-PSI - P1", "TRACE-PSI - Context")
panel_c[, comparison := factor(comparison, levels = rev(comparison_order))]
panel_c[, comparison_colour := fifelse(model_A == "P2_TRACEpsi", "P2_TRACEpsi",
                                        fifelse(model_A == "A_context_only", "A_context_only",
                                                "P1_provenance"))]
write_plot_source(panel_c, "Figure3_panelC_paired_delta_AUROC.tsv")

p_c <- ggplot(panel_c, aes(delta_AUROC, comparison, colour = comparison_colour)) +
  geom_vline(xintercept = 0, colour = TRACE_COLORS[["neutral_mid"]],
             linetype = "22", linewidth = 0.42) +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = comparison), linewidth = 0.78) +
  geom_point(size = 2.55) +
  geom_text(aes(label = sprintf("%+.3f", delta_AUROC),
                x = ifelse(delta_AUROC >= 0, CI_high + 0.008, CI_low - 0.008),
                hjust = ifelse(delta_AUROC >= 0, 0, 1)),
            size = 1.9, family = FONT_FAMILY, colour = TRACE_COLORS[["ink"]]) +
  scale_colour_manual(values = TRACE_COLORS[c("P1_provenance", "A_context_only", "P2_TRACEpsi")],
                      guide = "none") +
  scale_x_continuous(limits = c(-0.155, 0.295), breaks = c(-0.1, 0, 0.1, 0.2)) +
  labs(title = "Paired model gains",
       subtitle = "Same breadth=1 loci; paired bootstrap",
       x = "Delta AUROC (95% CI)", y = NULL) +
  theme_trace() +
  theme(axis.text.y = element_text(size = 5.45))

design <- "
AAB
AAC
"
figure <- p_a + p_b + p_c +
  plot_layout(design = design, widths = c(1, 1, 1.15), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT_FAMILY))

save_trace_figure(
  figure, "Figure3_main", height_mm = 112,
  panel_ids = c("a", "b", "c")
)
