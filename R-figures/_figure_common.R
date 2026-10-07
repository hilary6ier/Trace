suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
})

options(stringsAsFactors = FALSE)

FIG_WIDTH_IN <- 183 / 25.4
FONT_FAMILY <- "Arial"

COL <- c(
  ink = "#172B3A",
  muted = "#607381",
  grid = "#D9E1E5",
  light = "#EEF2F4",
  BID = "#2468A2",
  BACS = "#D97824",
  ELAP = "#2A8F6A",
  TRACE = "#2468A2",
  context = "#D99A2B",
  provenance = "#8A6B8F",
  breadth = "#8A989F",
  residual = "#C45D4A",
  observed = "#B33D4A"
)

theme_manuscript <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(size = base_size + 1.2, face = "bold", colour = COL["ink"], hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.4, colour = COL["muted"], margin = margin(b = 4)),
      axis.title = element_text(size = base_size, colour = COL["ink"]),
      axis.text = element_text(size = base_size - 0.5, colour = COL["ink"]),
      axis.line = element_line(linewidth = 0.35, colour = COL["ink"]),
      axis.ticks = element_line(linewidth = 0.3, colour = COL["ink"]),
      axis.ticks.length = unit(1.6, "mm"),
      strip.background = element_rect(fill = COL["light"], colour = NA),
      strip.text = element_text(size = base_size, face = "bold", colour = COL["ink"], margin = margin(3, 3, 3, 3)),
      legend.title = element_text(size = base_size - 0.3, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = unit(3.4, "mm"),
      legend.key.width = unit(4.2, "mm"),
      plot.margin = margin(5, 6, 5, 5)
    )
}

panel_tag_theme <- theme(
  plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = 12, colour = COL["ink"]),
  plot.tag.position = c(0, 1)
)

read_tsv <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(as.data.table(read.delim(gzfile(path), sep = "\t", header = TRUE,
                                    na.strings = c("NA", "NaN", ""), check.names = FALSE)))
  }
  fread(path, sep = "\t", na.strings = c("NA", "NaN", ""), check.names = FALSE)
}

write_source <- function(x, filename) {
  dir.create(file.path(OUT_DIR, "plot_source"), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), file.path(OUT_DIR, "plot_source", filename), sep = "\t", na = "NA", quote = FALSE)
}

save_figure <- function(plot, stem, height_in) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  tif_file <- file.path(OUT_DIR, paste0(stem, "_600dpi.tiff"))
  png_file <- file.path(OUT_DIR, paste0(stem, "_preview.png"))

  cairo_pdf(pdf_file, width = FIG_WIDTH_IN, height = height_in, family = FONT_FAMILY, onefile = TRUE)
  print(plot)
  invisible(dev.off())

  agg_tiff(tif_file, width = FIG_WIDTH_IN, height = height_in, units = "in", res = 600,
           compression = "lzw", background = "white")
  print(plot)
  invisible(dev.off())

  agg_png(png_file, width = FIG_WIDTH_IN, height = height_in, units = "in", res = 240,
          background = "white")
  print(plot)
  invisible(dev.off())
}

roc_points <- function(y, score) {
  d <- data.table(y = as.integer(y), score = as.numeric(score))
  g <- d[, .(tp_add = sum(y), fp_add = .N - sum(y)), by = score][order(-score)]
  g[, `:=`(tp = cumsum(tp_add), fp = cumsum(fp_add))]
  out <- rbind(data.table(score = Inf, tp_add = 0L, fp_add = 0L, tp = 0L, fp = 0L), g, fill = TRUE)
  out[, `:=`(TPR = tp / sum(d$y), FPR = fp / sum(1L - d$y))]
  out
}

pr_points <- function(y, score) {
  d <- data.table(y = as.integer(y), score = as.numeric(score))
  g <- d[, .(tp_add = sum(y), n_add = .N), by = score][order(-score)]
  g[, `:=`(tp = cumsum(tp_add), n = cumsum(n_add))]
  g[, `:=`(recall = tp / sum(d$y), precision = tp / n)]
  rbind(data.table(score = Inf, tp_add = 0L, n_add = 0L, tp = 0L, n = 0L, recall = 0, precision = 1), g, fill = TRUE)
}

tie_budget_curve <- function(y, score) {
  d <- data.table(y = as.integer(y), score = as.numeric(score))
  groups <- d[, .(group_n = .N, group_pos = sum(y)), by = score][order(-score)]
  rows <- vector("list", nrow(groups))
  n_before <- 0L
  pos_before <- 0
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    ks <- seq.int(n_before + 1L, n_before + g$group_n)
    expected <- pos_before + (ks - n_before) * g$group_pos / g$group_n
    rows[[i]] <- data.table(k = ks, expected_positives = expected, boundary_score = g$score)
    n_before <- n_before + g$group_n
    pos_before <- pos_before + g$group_pos
  }
  out <- rbindlist(rows)
  out[, `:=`(
    fraction_screened = k / .N,
    fraction_recovered = expected_positives / sum(d$y),
    expected_precision = expected_positives / k,
    expected_lift = (expected_positives / k) / mean(d$y)
  )]
  out
}
