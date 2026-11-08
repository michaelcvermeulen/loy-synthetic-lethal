# Reproduce frozen Figure 10 without on-canvas titles or explanatory captions.
# No models or simulations are rerun. The manuscript caption remains authoritative.
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
root <- dirname(dirname(dirname(dirname(script))))
fig_dir <- file.path(root, "figures/chapter3")
data_dir <- file.path(fig_dir, "source_data/depmap_benchmark_detection")
export <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07")
script_dir <- file.path(export, "scripts/loy_reanalysis_v1")
derived_dir <- file.path(export, "derived/loy_reanalysis_v1")
style_dir <- file.path(root, "skills/thesis-figure-style")
stem <- "ch03_supp_depmap_benchmark_detection"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Make the same bundled Open Sans fonts available to all three export devices.
cache_dir <- file.path(data_dir, ".fontconfig_cache")
dir.create(cache_dir, showWarnings = FALSE)
fontconfig <- file.path(data_dir, "fonts.conf")
writeLines(c('<?xml version="1.0"?>', '<!DOCTYPE fontconfig SYSTEM "fonts.dtd">',
             '<fontconfig>', sprintf('  <dir>%s/assets/fonts</dir>', style_dir),
             sprintf('  <cachedir>%s</cachedir>', cache_dir), '</fontconfig>'), fontconfig)
Sys.setenv(FONTCONFIG_FILE = fontconfig)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})
source(file.path(style_dir, "scripts/thesis_style.R"))
register_thesis_fonts(style_dir)

# Reuse the original style definitions without its HPC paths or write operations.
style_expr <- parse(file.path(script_dir, "plot_style.R"))
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
for (name in c("theme_house", "wrap_text")) {
    match <- Filter(function(x) is_assignment(x, name), as.list(style_expr))
    stopifnot(length(match) == 1)
    eval(match[[1]])
}
panel_labs <- function(panel_in, title, ..., caption = NULL) labs(title = NULL, caption = NULL, ...)
read_result <- function(name) as.data.table(readRDS(file.path(derived_dir, paste0(name, ".rds"))))
save_source_data <- function(dt, name) {
    original <- fread(file.path(export, "results/loy_reanalysis_v1/figure_source_data", paste0(name, ".csv")), encoding = "UTF-8")
    # Compare the exact plotted values against the frozen CSV before writing local copies.
    current <- copy(dt)
    for (j in names(current)) if (is.factor(current[[j]])) set(current, j = j, value = as.character(current[[j]]))
    comparison <- all.equal(as.data.frame(current), as.data.frame(original), tolerance = 1e-12, check.attributes = FALSE)
    if (!isTRUE(comparison)) {
        print(comparison)
        print(head(current))
        print(head(original))
        stop("Plot source differs from the frozen CSV: ", name)
    }
    fwrite(dt, file.path(data_dir, paste0(name, ".csv")))
}
params <- jsonlite::fromJSON(file.path(script_dir, "parameters.json"))
full <- 9.2

# Reuse only the common saved-data setup and Figure 10 plot construction.
original_lines <- readLines(file.path(script_dir, "21_figures.R"), warn = FALSE, encoding = "UTF-8")
first <- which(original_lines == 'power_types <- read_result("power_cancer_types")')
last <- grep('^type_rows <- tested_types', original_lines)
stopifnot(length(first) == 1, length(last) == 1, last > first)
eval(parse(text = original_lines[first:last]))
first <- grep('^class_levels <- unique\\(power_bench\\$class\\)', original_lines)
last <- grep('^save_figure\\(tagged\\(p10a', original_lines) - 1L
stopifnot(length(first) == 1, length(last) == 1, last > first)
eval(parse(text = original_lines[first:last]))
stopifnot(nrow(bench_heat) == 750, nrow(positive_heat) == 120, recovered == 16)
for (name in c("p10a", "p10b")) {
    assign(name, get(name) + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                                  plot.caption = element_blank(), plot.margin = margin(3, 3, 3, 3)))
}
plot <- (p10a + p10b + plot_layout(widths = c(nrow(bench_columns), 4.6), guides = "collect") +
         plot_annotation(tag_levels = "a")) &
        theme(legend.position = "bottom", plot.tag = element_text(size = 10, face = "bold"))
width <- 9.2
height <- 8.8
ggsave(file.path(fig_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, device = ragg::agg_png)
node <- file.path(Sys.getenv("USERPROFILE"), ".cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node.exe")
status <- system2(node, shQuote(c(file.path(fig_dir, "code/depmap_svg_to_pdf.mjs"),
    file.path(fig_dir, paste0(stem, ".svg")), file.path(fig_dir, paste0(stem, ".pdf")),
    file.path(style_dir, "assets/fonts"))))
stopifnot(status == 0)
jsonlite::write_json(list(source_export = "2026-10-07", source_figure = "fig10_power_benchmarks",
    changes = "Removed on-canvas titles and explanatory captions; compact 8.8-inch height.",
    width_inches = width, height_inches = height, font = "Open Sans", models_refitted = FALSE,
    saved_panel_data_match = TRUE, simulated_benchmark_cells = nrow(bench_heat),
    observed_gene_cells = nrow(positive_heat), observed_recovered = recovered),
    file.path(data_dir, "provenance.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Rendered S03 without on-canvas titles/captions; both panel data tables match the frozen export.\n")

