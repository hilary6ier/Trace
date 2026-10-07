suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

if (!exists("SCRIPT_DIR")) {
  stop("SCRIPT_DIR must be defined before sourcing _figure_common.R")
}
if (!exists("REPO_ROOT")) {
  REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
}

OUT_DIR <- file.path(REPO_ROOT, "R-figures")
SOURCE_DIR <- file.path(OUT_DIR, "plot_source")
QA_DIR <- file.path(OUT_DIR, "qa")
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QA_DIR, recursive = TRUE, showWarnings = FALSE)

# Fontconfig needs a writable cache in the restricted cloud runtime.
FONT_CACHE <- file.path(tempdir(), "trace-font-cache")
dir.create(FONT_CACHE, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = FONT_CACHE)

# Generic sans is used in grobs so the alignment probe's base-PDF device remains
# portable; svglite explicitly maps it to Liberation Sans below.
FONT_FAMILY <- "sans"

TRACE_COLORS <- c(
  ink = "#252A30",
  neutral_dark = "#676D73",
  neutral_mid = "#9AA0A6",
  neutral_light = "#D6D9DC",
  neutral_pale = "#F1F3F5",
  HEK293T = "#31688E",
  HeLa = "#C75B7A",
  BID = "#D07A2D",
  BACS = "#3F7EA6",
  ELAP = "#7A5AA6",
  DRS = "#2A9D8F",
  P0_breadth = "#9AA0A6",
  P1_provenance = "#8F6BAE",
  A_context_only = "#2A9D8F",
  P2_TRACEpsi = "#245B8A"
)

MODEL_LABELS <- c(
  P0_breadth = "P0: Breadth",
  P1_provenance = "P1: Provenance",
  A_context_only = "Context-only",
  P2_TRACEpsi = "TRACE-PSI"
)

theme_trace <- function(base_size = 6.6) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = TRACE_COLORS[["ink"]]),
      axis.ticks = element_line(linewidth = 0.30, colour = TRACE_COLORS[["ink"]]),
      axis.ticks.length = grid::unit(1.3, "mm"),
      axis.title = element_text(size = base_size, colour = TRACE_COLORS[["ink"]]),
      axis.text = element_text(size = base_size - 0.4, colour = TRACE_COLORS[["ink"]]),
      plot.title = element_text(size = 7.2, face = "bold", hjust = 0,
                                margin = margin(b = 2.5)),
      plot.subtitle = element_text(size = 6.1, colour = TRACE_COLORS[["neutral_dark"]],
                                   margin = margin(b = 3.0)),
      plot.tag = element_text(size = 8.0, face = "bold", colour = TRACE_COLORS[["ink"]]),
      plot.tag.position = c(0, 1),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = 5.9),
      legend.key.height = grid::unit(3.0, "mm"),
      legend.key.width = grid::unit(4.2, "mm"),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.2, face = "bold", margin = margin(b = 1.5)),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

theme_set(theme_trace())

write_plot_source <- function(x, filename) {
  data.table::fwrite(as.data.table(x), file.path(SOURCE_DIR, filename), sep = "\t", na = "NA")
}

read_tsv <- function(...) {
  data.table::fread(file.path(REPO_ROOT, ...), sep = "\t", na.strings = c("NA", ""))
}

read_tsv_gz <- function(...) {
  path <- file.path(REPO_ROOT, ...)
  data.table::fread(cmd = paste("gzip -cd", shQuote(path)), sep = "\t", na.strings = c("NA", ""))
}

save_trace_figure <- function(plot, stem, height_mm, panel_ids,
                              row_groups = NULL, column_groups = NULL,
                              exemptions = list(), width_mm = 183) {
  stopifnot(inherits(plot, "patchwork"))
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  alignment_helper <- file.path(
    REPO_ROOT, ".agents", "skills", "nature-figure", "scripts", "panel_alignment.R"
  )
  source(alignment_helper, local = environment())
  measured_rows <- .nature_alignment_panel_rows(patchwork::patchworkGrob(plot))
  message("Alignment plot areas: ", nrow(measured_rows), " [",
          paste(measured_rows$name, collapse = ", "), "]")
  all_panel_names <- patchwork::patchworkGrob(plot)$layout$name
  message("All panel-like grobs: ", paste(all_panel_names[grepl("panel", all_panel_names)], collapse = ", "))
  require_patchwork_panel_alignment(
    plot = plot,
    manifest_path = file.path(QA_DIR, paste0(stem, ".alignment-layout.json")),
    report_path = file.path(QA_DIR, paste0(stem, ".alignment.json")),
    overlay_svg = file.path(QA_DIR, paste0(stem, ".alignment.svg")),
    width_in = width_in,
    height_in = height_in,
    panel_ids = panel_ids,
    row_groups = row_groups,
    column_groups = column_groups,
    exemptions = exemptions,
    audit_script = file.path(
      REPO_ROOT, ".agents", "skills", "nature-figure", "scripts", "audit_panel_alignment.py"
    ),
    python = Sys.which("python"),
    tolerance_pt = 1.5,
    gutter_tolerance_pt = 1.5,
    strict = TRUE
  )

  svg_file <- file.path(OUT_DIR, paste0(stem, ".svg"))
  pdf_file <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  tiff_file <- file.path(OUT_DIR, paste0(stem, "_600dpi.tiff"))
  png_file <- file.path(OUT_DIR, paste0(stem, "_preview.png"))

  svglite::svglite(
    svg_file, width = width_in, height = height_in,
    system_fonts = list(sans = "Liberation Sans")
  )
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    pdf_file, width = width_in, height = height_in,
    family = FONT_FAMILY, onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    tiff_file, width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw", background = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    png_file, width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  grDevices::dev.off()

  invisible(c(svg = svg_file, pdf = pdf_file, tiff = tiff_file, png = png_file))
}

percent1 <- function(x) sprintf("%.1f%%", 100 * x)
number2 <- function(x) sprintf("%.2f", x)
