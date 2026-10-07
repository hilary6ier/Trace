#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", argv[grepl("^--file=", argv)])
OUT_DIR <- dirname(normalizePath(script_arg))
ROOT <- normalizePath(file.path(OUT_DIR, ".."))
source(file.path(OUT_DIR, "_figure_common.R"))

patterns <- read_tsv(file.path(ROOT, "05_aim2", "aim2b_DRS_pattern_rates.tsv"))
breadth <- read_tsv(file.path(ROOT, "05_aim2", "aim2b_DRS_breadth_rates.tsv"))
model <- read_tsv(file.path(ROOT, "05_aim2", "aim2b_DRS_breadth_model.tsv"))
perm_frozen <- read_tsv(file.path(ROOT, "05_aim2", "aim2b_DRS_motif_stratified_permutation.tsv"))
union <- read_tsv(file.path(ROOT, "05_aim2", "aim2b_HeLa_source_union_DRS.tsv"))

# Panel A: assay membership plus pattern-specific DRS re-observation.
pattern_order <- c("BID+BACS+ELAP", "BID+ELAP", "BID+BACS", "BACS+ELAP", "ELAP", "BACS", "BID")
patterns[, support_pattern := factor(support_pattern, levels = rev(pattern_order))]
patterns[, count_label := sprintf("%d/%d", DRS_supported_n, n)]
membership <- CJ(support_pattern = levels(patterns$support_pattern), assay = c("BID", "BACS", "ELAP"))
membership[, present := mapply(function(p, a) a %in% strsplit(as.character(p), "+", fixed = TRUE)[[1]],
                               support_pattern, assay)]
membership[, assay := factor(assay, levels = c("BID", "BACS", "ELAP"))]
membership[, support_pattern := factor(support_pattern, levels = levels(patterns$support_pattern))]
write_source(patterns, "Figure2_panelA_pattern_rates.tsv")
write_source(membership, "Figure2_panelA_assay_membership.tsv")

pA_table <- ggplot(membership, aes(as.numeric(assay), support_pattern)) +
  geom_point(shape = 22, size = 3.2, stroke = 0.4, fill = COL["light"], colour = COL["grid"]) +
  geom_point(data = membership[present == TRUE], aes(fill = assay), shape = 22, size = 3.2,
             stroke = 0.35, colour = "white") +
  geom_text(data = patterns, aes(x = 4.25, y = support_pattern, label = count_label),
            inherit.aes = FALSE, family = FONT_FAMILY, size = 2.55, colour = COL["ink"])+
  scale_fill_manual(values = unname(c(COL["BID"], COL["BACS"], COL["ELAP"])), guide = "none") +
  scale_x_continuous(breaks = c(1, 2, 3, 4.25), labels = c("BID", "BACS", "ELAP", "DRS+/n"),
                     limits = c(0.55, 4.7), position = "top") +
  labs(title = "Source-pattern evidence landscape", x = NULL, y = NULL, tag = "A") +
  theme_manuscript() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(face = "bold", size = 6.6),
        plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

pA_forest <- ggplot(patterns, aes(DRS_support_rate, support_pattern)) +
  geom_errorbar(aes(xmin = wilson95_low, xmax = wilson95_high), orientation = "y", width = 0,
                linewidth = 0.7, colour = COL["ink"]) +
  geom_point(aes(fill = as.numeric(support_pattern)), shape = 21, size = 2.6,
             colour = "white", stroke = 0.25) +
  scale_fill_gradient(low = COL["BID"], high = COL["observed"], guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.90),
                     breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(x = "DRS re-observation (Wilson 95% CI)", y = NULL) +
  theme_manuscript() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(25, 5, 5, 2))

pA <- pA_table + pA_forest + plot_layout(widths = c(1.18, 1))

# Panel B: pooled breadth rates.
breadth[, breadth_label := factor(chemistry_breadth, levels = 1:3,
                                  labels = c("1 chemistry", "2 chemistries", "3 chemistries"))]
breadth[, label := sprintf("%.1f%%  (%d/%d)", 100 * DRS_support_rate, DRS_supported_n, n)]
breadth[, `:=`(label_x = chemistry_breadth + c(0.12, 0.12, -0.12), label_hjust = c(0, 0, 1))]
write_source(breadth, "Figure2_panelB_breadth_rates.tsv")

pB <- ggplot(breadth, aes(chemistry_breadth, DRS_support_rate)) +
  geom_line(colour = COL["TRACE"], linewidth = 0.75) +
  geom_errorbar(aes(ymin = wilson95_low, ymax = wilson95_high), width = 0.08,
                linewidth = 0.65, colour = COL["TRACE"]) +
  geom_point(shape = 21, fill = COL["TRACE"], colour = "white", stroke = 0.3, size = 3) +
  geom_text(aes(x = label_x, label = label, hjust = label_hjust), family = FONT_FAMILY,
            size = 2.45, colour = COL["ink"]) +
  scale_x_continuous(breaks = 1:3, limits = c(0.85, 3.72)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.9),
                     breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(title = "Pooled breadth gradient", x = "Distinct source chemistries", y = "DRS re-observation", tag = "B") +
  theme_manuscript() +
  theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

# Panel C: motif-cluster-robust breadth effect.
pc_data <- model[term == "chemistry_breadth"]
write_source(pc_data, "Figure2_panelC_breadth_OR.tsv")

pC <- ggplot(pc_data, aes(OR, 1)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = COL["muted"], linewidth = 0.4) +
  geom_errorbar(aes(xmin = OR_CI_low, xmax = OR_CI_high), orientation = "y", width = 0,
                linewidth = 0.8, colour = COL["ink"]) +
  geom_point(shape = 21, fill = COL["observed"], colour = "white", stroke = 0.3, size = 3.2) +
  annotate("text", x = pc_data$OR, y = 1.16,
           label = sprintf("OR %.2f (%.2f-%.2f)", pc_data$OR, pc_data$OR_CI_low, pc_data$OR_CI_high),
           family = FONT_FAMILY, size = 2.55, colour = COL["ink"]) +
  scale_x_log10(limits = c(0.8, 7.2), breaks = c(1, 2, 4, 6)) +
  scale_y_continuous(NULL, breaks = NULL, limits = c(0.82, 1.28)) +
  labs(title = "Per-chemistry effect", x = "Odds ratio (95% CI, log scale)", tag = "C") +
  theme_manuscript() +
  theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

# Panel D: deterministic re-expression of the frozen motif-stratified null.
set.seed(20261003)
perm_data <- union[, .(motif, chemistry_breadth, DRS_reported_support)]
y <- perm_data$DRS_reported_support
b <- perm_data$chemistry_breadth
groups <- split(seq_len(nrow(perm_data)), perm_data$motif)
n1 <- sum(y == 1)
n0 <- sum(y == 0)
null_stat <- numeric(20000)
for (i in seq_along(null_stat)) {
  yp <- y
  for (ids in groups) yp[ids] <- sample(yp[ids], length(ids), replace = FALSE)
  null_stat[i] <- mean(b[yp == 1]) - mean(b[yp == 0])
}
pd_null <- data.table(permutation = seq_along(null_stat), statistic = null_stat)
pd_summary <- pd_null[, .(estimate = median(statistic), low = quantile(statistic, 0.025),
                          high = quantile(statistic, 0.975))]
pd_observed <- data.table(estimate = perm_frozen$observed, low = perm_frozen$observed,
                          high = perm_frozen$observed)
pd_display <- pd_null[seq(1, .N, length.out = 500)]
write_source(pd_null, "Figure2_panelD_motif_conditioned_null.tsv")
write_source(perm_frozen, "Figure2_panelD_frozen_permutation_result.tsv")

pD <- ggplot() +
  geom_jitter(data = pd_display, aes(statistic, 1), height = 0.055, width = 0,
              colour = COL["breadth"], alpha = 0.22, size = 0.65) +
  geom_errorbar(data = pd_summary, aes(y = 1, xmin = low, xmax = high), orientation = "y",
                width = 0, linewidth = 1.0, colour = COL["ink"]) +
  geom_point(data = pd_summary, aes(estimate, 1), shape = 21, fill = "white",
             colour = COL["ink"], stroke = 0.8, size = 2.6) +
  geom_point(data = pd_observed, aes(estimate, 2), shape = 23, fill = COL["observed"],
             colour = "white", stroke = 0.3, size = 3.2) +
  annotate("text", x = perm_frozen$observed, y = 2.23,
           label = sprintf("Delta mean breadth = %.3f\nP = %.4f", perm_frozen$observed,
                           perm_frozen$p_two_sided),
           hjust = 1, family = FONT_FAMILY, size = 2.35, colour = COL["ink"]) +
  scale_y_continuous(NULL, breaks = c(1, 2), labels = c("Motif-conditioned null", "Observed"),
                     limits = c(0.78, 2.48)) +
  scale_x_continuous(limits = range(c(pd_null$statistic, perm_frozen$observed)) + c(-0.03, 0.03)) +
  labs(title = "Motif-conditioned permutation", x = "Mean breadth: DRS+ minus DRS-nonreported", tag = "D") +
  theme_manuscript() +
  theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0, 1))

design <- "
AAB
AAC
DDD
"
figure <- wrap_plots(A = pA, B = pB, C = pC, D = pD, design = design) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_figure(figure, "Figure2_main", height_in = 6.0)
