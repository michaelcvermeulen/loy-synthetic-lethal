# Figure style shared by every plotting stage: the thesis house style already
# used in this collection (Open Sans; Y-retaining blue, Y-depleted vermilion).
# Source after common.R.

style_dir <- file.path(
  root_dir, "results", "rna_y_vulnerabilities_v1", "style_snapshots", "thesis-figure-style-v0.3"
)

# One fontconfig file makes Open Sans visible to Cairo (PDF), ragg (PNG) and
# svglite alike. It has to be in place before any graphics package loads.
setup_fonts <- function() {
  font_dir <- file.path(results_dir, ".fontconfig")
  dir.create(file.path(font_dir, "cache"), showWarnings = FALSE, recursive = TRUE)
  conf <- file.path(font_dir, "fonts.conf")
  writeLines(c(
    "<?xml version=\"1.0\"?>",
    "<!DOCTYPE fontconfig SYSTEM \"fonts.dtd\">",
    "<fontconfig>",
    sprintf(
      "  <include ignore_missing=\"yes\">%s</include>",
      normalizePath(file.path(R.home(), "..", "..", "etc", "fonts", "fonts.conf"), mustWork = FALSE)
    ),
    sprintf("  <dir>%s</dir>", file.path(style_dir, "assets", "fonts")),
    sprintf("  <cachedir>%s</cachedir>", file.path(font_dir, "cache")),
    "</fontconfig>"
  ), conf)
  Sys.setenv(FONTCONFIG_FILE = conf)
}
setup_fonts()

suppressPackageStartupMessages(library(ggplot2))
source(file.path(style_dir, "scripts", "thesis_style.R"))

house_palette <- thesis_palette(style_dir)
y_status_colours <- c(
  "Y-positive" = unname(house_palette["Y-retaining"]),
  "Y-negative" = unname(house_palette["Y-depleted"])
)
# Accent colours for classes that are not a Y-status group. chrY is the
# perturbation itself, so it never shares a colour with a result class.
accent_colours <- c(
  chrY = "#C44E52", x_partner = "#00A087", par = "#8172B3",
  hit = "#222222", other = "#B8B8B8", female = "#7E6148"
)

theme_house <- function(base_size = 9) {
  theme_thesis(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size, face = "bold"),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = base_size - 1, colour = "#444444"),
      plot.caption = element_text(size = 6, colour = "#666666", hjust = 0),
      plot.caption.position = "plot",
      axis.text = element_text(size = base_size - 2, colour = "#222222"),
      strip.text = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size - 2),
      legend.key.size = grid::unit(0.7, "lines"),
      legend.margin = margin(0, 0, 0, 0)
    )
}

# Titles and captions wrapped to the width of the panel they sit on, so that
# neighbouring panels of a multi-panel figure never run into each other.
# `panel_in` is the panel width in inches.
wrap_text <- function(text, width) paste(strwrap(text, width = width), collapse = "\n")
panel_labs <- function(panel_in, title, ..., caption = NULL) {
  labs(
    title = wrap_text(title, floor(panel_in / 0.082)),
    caption = if (is.null(caption)) NULL else wrap_text(caption, floor(panel_in / 0.047)),
    ...
  )
}

scale_colour_y_status <- function(...) scale_colour_manual(values = y_status_colours, ...)
scale_fill_y_status <- function(...) scale_fill_manual(values = y_status_colours, ...)

# Save one figure as PDF (for the thesis), SVG and 300 dpi PNG, with identical dimensions.
save_figure <- function(plot, name, width, height, dir = file.path(results_dir, "figures")) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths <- file.path(dir, paste0(name, c(".pdf", ".svg", ".png")))
  ggsave(paths[1], plot, width = width, height = height, device = grDevices::cairo_pdf)
  ggsave(paths[2], plot, width = width, height = height, device = svglite::svglite)
  ggsave(paths[3], plot, width = width, height = height, dpi = 300, device = ragg::agg_png)
  invisible(paths)
}

# Write the data behind a figure next to it.
save_source_data <- function(dt, name) {
  dir <- file.path(results_dir, "figure_source_data")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(dt, file.path(dir, paste0(name, ".csv")))
  invisible(dt)
}
