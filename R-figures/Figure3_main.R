#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", argv[grepl("^--file=", argv)])
OUT_DIR <- dirname(normalizePath(script_arg))
ROOT <- normalizePath(file.path(OUT_DIR, ".."))
source(file.path(OUT_DIR, "_figure_common.R"))

internal <- read_tsv(file.path(ROOT, "07_tracepsi", "phase4b_v2_internal_metrics.tsv"))
drs <- read_tsv(file.path(ROOT, "07_tracepsi", "phase4b_v2_HeLa_DRS_metrics.tsv"))
delta <- read_tsv(file.path(ROOT, "07_tracepsi", "phase4b_v2_HeLa_DRS_breadth1_deltaAUC.tsv"))

model_levels <- c("P0_breadth", "P1_provenance", "A_context_only", "P2_TRACEpsi")
model_labels <- c("P0 breadth", "P1 provenance", "Context only", "TRACE-Psi")
model_colours <- c("P0_breadth" = COL["breadth"], "P1_provenance" = COL["provenance"],
                   "A_context_only" = COL["context"], "P2_TRACEpsi" = COL["TRACE"])

# Panel A: held-assay transfer matrix.
pa_data <- internal[, .(held_target, model, n, positives, AUROC, AUROC_CI_low, AUROC_CI_high)]
pa_data[, model := factor(model, levels = rev(model_levels), labels = rev(model_labels))]
pa_data[, held_target := factor(held_target, levels = c("BACS", "BID", "ELAP"))]
write_source(pa_data, "Figure3_panelA_HeLa_LOTO_matrix.tsv")

pA <- ggplot(pa_data, aes(held_target, model, fill = AUROC)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.3f", AUROC)), family = FONT_FAMILY, size = 3.0,
            fontface = "bold", colour = COL["ink"]) +
  scale_fill_gradient2(low = "#F3F4F4", mid = "#B9D6E8", high = COL["TRACE"], midpoint = 0.60,
                       limits = c(0.50, 0.76), guide = "none") +
  coord_equal() +
  labs(title = "Held-assay transfer in HeLa", subtitle = "Leave-one-technology-out AUROC; darker is higher",
       x = "Held target assay", y = NULL, tag = "A") +
  theme_manuscript() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

# Panel B: breadth=1 DRS ranking.
pb_data <- drs[subset == "breadth1",
               .(model, n, positives, prevalence, AUROC, AUROC_CI_low, AUROC_CI_high,
                 average_precision, Lift10, Precision10)]
pb_data[, model := factor(model, levels = rev(model_levels), labels = rev(model_labels))]
write_source(pb_data, "Figure3_panelB_HeLa_breadth1_AUROC.tsv")

pB <- ggplot(pb_data, aes(AUROC, model, colour = model)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = COL["muted"], linewidth = 0.4) +
  geom_errorbar(aes(xmin = AUROC_CI_low, xmax = AUROC_CI_high), orientation = "y", width = 0,
                linewidth = 0.75) +
  geom_point(size = 2.7) +
  geom_text(aes(label = sprintf("%.3f", AUROC)), nudge_y = 0.23, family = FONT_FAMILY,
            size = 2.35, colour = COL["ink"]) +
  scale_colour_manual(values = rev(unname(model_colours)), guide = "none") +
  scale_x_continuous(limits = c(0.38, 0.72), breaks = c(0.4, 0.5, 0.6, 0.7)) +
  labs(title = "Ranking within breadth=1", subtitle = "HeLa source union to independent DRS",
       x = "AUROC (bootstrap 95% CI)", y = NULL, tag = "B") +
  theme_manuscript() +
  theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

# Panel C: paired bootstrap differences.
delta[, comparison := paste(model_A, "minus", model_B)]
comparison_labels <- c(
  "P1_provenance minus P0_breadth" = "P1 - P0",
  "A_context_only minus P1_provenance" = "Context - P1",
  "P2_TRACEpsi minus P1_provenance" = "TRACE-Psi - P1",
  "P2_TRACEpsi minus A_context_only" = "TRACE-Psi - Context"
)
delta[, label := unname(comparison_labels[comparison])]
delta[, label := factor(label, levels = rev(unname(comparison_labels)))]
delta[, emphasis := ifelse(comparison == "P2_TRACEpsi minus A_context_only", "Primary increment", "Benchmark contrast")]
write_source(delta, "Figure3_panelC_paired_delta_AUROC.tsv")

pC <- ggplot(delta, aes(delta_AUROC, label, colour = emphasis)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL["muted"], linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0,
                linewidth = 0.75) +
  geom_point(aes(shape = emphasis), size = 2.7) +
  geom_text(aes(label = sprintf("%+.3f", delta_AUROC)), nudge_y = 0.22,
            family = FONT_FAMILY, size = 2.3, colour = COL["ink"]) +
  scale_colour_manual(values = setNames(c(unname(COL["TRACE"]), unname(COL["provenance"])),
                                        c("Primary increment", "Benchmark contrast")),
                      guide = "none") +
  scale_shape_manual(values = c("Primary increment" = 18, "Benchmark contrast" = 16), guide = "none") +
  scale_x_continuous(limits = c(-0.115, 0.26), breaks = c(-0.1, 0, 0.1, 0.2)) +
  labs(title = "Paired improvement in breadth=1", x = "Delta AUROC (bootstrap 95% CI)", y = NULL, tag = "C") +
  theme_manuscript() +
  theme(plot.margin = margin(12, 6, 5, 5),
        plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

design <- "
AABBB
AACCC
"
figure <- wrap_plots(A = pA, B = pB, C = pC, design = design) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_figure(figure, "Figure3_main", height_in = 4.8)
