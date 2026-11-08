# Frozen Figure 9 with former S07 boxplots added as panels d and e.
# No models or simulations are rerun. The manuscript caption remains authoritative.
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), winslash = "/")
root <- dirname(dirname(dirname(dirname(script))))
fig_dir <- file.path(root, "figures/chapter3")
data_dir <- file.path(fig_dir, "source_data/depmap_cancer_type_detection_limits")
export <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07")
script_dir <- file.path(export, "scripts/loy_reanalysis_v1")
derived_dir <- file.path(export, "derived/loy_reanalysis_v1")
style_dir <- file.path(root, "skills/thesis-figure-style")
stem <- "ch03_main_depmap_cancer_type_detection_limits"
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

# Evaluate only the saved-data plotting block for Figure 9, not the analysis pipeline.
original_lines <- readLines(file.path(script_dir, "21_figures.R"), warn = FALSE, encoding = "UTF-8")
first <- which(original_lines == 'power_types <- read_result("power_cancer_types")')
last <- grep('^save_figure\\(tagged\\(\\(p9a', original_lines) - 1L
stopifnot(length(first) == 1, length(last) == 1, last > first)
eval(parse(text = original_lines[first:last]))
stopifnot(identical(as.integer(enough), c(14L, 11L, 7L)), nrow(tested_types) == 30, half_unit == 6)
# Label every eligible group and the largest below-threshold groups in each facet.
# These are display choices only: the underlying points and eligibility stay fixed.
label_limits <- c(lineage = 4L, disease = 6L, subtype = 8L)
below_labels <- both_groups[testable == FALSE][order(-(n_negative + n_positive), label),
    head(.SD, label_limits[as.character(.BY$level)]), by = level]
labelled <- rbind(both_groups[testable == TRUE], below_labels, fill = TRUE)
abbreviations <- fread(file.path(data_dir, "panel_b_abbreviations.tsv"), encoding = "UTF-8")
stopifnot(!anyDuplicated(abbreviations$label))
labelled[, short_label := abbreviations$short_label[match(label, abbreviations$label)]]
if (anyNA(labelled$short_label)) print(labelled[is.na(short_label), .(level, label)])
stopifnot(!anyNA(labelled$short_label))
label_data <- copy(both_groups)
label_data[, short_label := ""]
label_data[labelled, on = .(level, label), short_label := i.short_label]
fwrite(labelled[, .(level, label, short_label, n_negative, n_positive, testable)],
       file.path(data_dir, "panel_b_displayed_labels.csv"))
p9b$layers <- Filter(function(layer) !inherits(layer$geom, "GeomTextRepel"), p9b$layers)
p9b <- p9b + ggrepel::geom_text_repel(data = label_data, aes(label = short_label),
    family = "Open Sans", size = 1.9, lineheight = .85, min.segment.length = 0,
    segment.size = .2, box.padding = .22, point.padding = .15, max.overlaps = Inf,
    max.iter = 30000, max.time = 3, force = 1.5, seed = 1)
for (name in c("p9a", "p9b", "p9c")) {
    assign(name, get(name) + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                                  plot.caption = element_blank(), plot.margin = margin(3, 3, 3, 3)))
}
# The two former S07 panels use their saved points and tests, without refitting.
candidate_dir <- file.path(fig_dir, "source_data/depmap_driver_nominal_candidates")
candidate_points <- fread(file.path(candidate_dir, "plotted_points.csv"), encoding = "UTF-8")
candidate_tests <- fread(file.path(candidate_dir, "saved_tests.csv"), encoding = "UTF-8")
stopifnot(nrow(candidate_points) == 106, setequal(candidate_tests$symbol, c("INSM1", "ATRX")))
fwrite(candidate_points, file.path(data_dir, "panels_de_points.csv"))
fwrite(candidate_tests, file.path(data_dir, "panels_de_saved_tests.csv"))
candidate_plot <- function(target_gene, driver, cancer) {
    pts <- copy(candidate_points[symbol == target_gene])
    hit <- candidate_tests[symbol == target_gene]
    stopifnot(nrow(hit) == 1,
        pts[ctx == 1 & y_status == "Y-positive", .N] == hit$n_pos_in,
        pts[ctx == 1 & y_status == "LOY", .N] == hit$n_neg_in,
        pts[ctx == 0 & y_status == "Y-positive", .N] == hit$n_pos_out,
        pts[ctx == 0 & y_status == "LOY", .N] == hit$n_neg_out)
    pts[, x := fifelse(ctx == 1, 0, 2.6) + as.integer(y_status == "LOY")]
    pts[, group := interaction(ctx, y_status)]
    low <- min(pts$chronos, 0); high <- max(pts$chronos, 0)
    spread <- high - low
    limits <- c(low - .08 * spread, high + .49 * spread)
    yy <- function(fraction) limits[1] + fraction * diff(limits)
    counts <- c(hit$n_pos_in, hit$n_neg_in, hit$n_pos_out, hit$n_neg_out)
    plot <- ggplot(pts, aes(x, chronos)) +
        geom_hline(yintercept = 0, colour = "#999999", linetype = "dashed", linewidth = .3) +
        geom_boxplot(aes(group = group, colour = y_status, fill = y_status),
                     width = .48, alpha = .2, outlier.shape = NA, linewidth = .35) +
        geom_point(aes(fill = y_status), shape = 21, colour = "white", stroke = .2,
                   size = 1.8, alpha = .85, position = position_jitter(width = .13, height = 0, seed = 20261008)) +
        scale_colour_manual(values = c(`Y-positive` = "#3C5488", LOY = "#E64B35"), guide = "none") +
        scale_fill_manual(values = c(`Y-positive` = "#3C5488", LOY = "#E64B35"), guide = "none") +
        scale_x_continuous(breaks = c(0, 1, 2.6, 3.6),
            labels = sprintf("%s\nn = %d", rep(c("Y+", "LOY"), 2), counts),
            limits = c(-.58, 4.18), expand = expansion(mult = 0)) +
        scale_y_continuous(limits = limits, expand = expansion(mult = 0), n.breaks = 4) +
        labs(x = NULL, y = paste(target_gene, "gene effect (Chronos)")) +
        coord_cartesian(clip = "off") + theme_house(8) +
        theme(axis.text = element_text(size = 7), axis.title.y = element_text(size = 8),
              panel.grid = element_blank(), plot.margin = margin(3, 4, 33, 3))
    format_p <- function(p) {
        if (p >= .001) return(sprintf("P = %.3f", p))
        sprintf("P = %.2e", p)
    }
    add_comparison <- function(plot, left, right, effect, p, fdr = NA_real_) {
        centre <- (left + right) / 2
        plot <- plot + annotate("path", x = c(left, left, right, right),
            y = yy(c(.745, .765, .765, .745)), linewidth = .35) +
            annotate("text", x = centre, y = yy(.985), label = sprintf("\u0394 = %+.3f", effect),
                family = "Open Sans", size = 2.65, vjust = 1) +
            annotate("text", x = centre, y = yy(.905), label = format_p(p),
                family = "Open Sans", size = 2.65, vjust = 1)
        if (is.finite(fdr)) plot <- plot + annotate("text", x = centre, y = yy(.825),
            label = sprintf("FDR = %.3f", fdr), family = "Open Sans", size = 2.65, vjust = 1)
        plot
    }
    plot <- add_comparison(plot, 0, 1, hit$inside, hit$inside_p, hit$inside_fdr)
    if (target_gene == "INSM1") plot <- add_comparison(plot, 2.6, 3.6, hit$outside, hit$outside_p)
    else {
        stopifnot(hit$n_neg_out == 0)
        plot <- plot + annotate("text", x = 3.6, y = yy(.39), label = "No models",
                               family = "Open Sans", colour = "#777777", size = 2.5)
    }
    # Place the driver labels below the two Y-group labels without widening the axes.
    plot + annotation_custom(grid::textGrob(paste0(driver, "-altered"),
        x = grid::unit((.5 + .58) / 4.76, "npc"), y = grid::unit(-24, "pt"),
        gp = grid::gpar(fontfamily = "Open Sans", fontsize = 7))) +
        annotation_custom(grid::textGrob(paste0(driver, "-unaltered"),
        x = grid::unit((3.1 + .58) / 4.76, "npc"), y = grid::unit(-24, "pt"),
        gp = grid::gpar(fontfamily = "Open Sans", fontsize = 7))) +
        annotation_custom(grid::textGrob(cancer, x = .5, y = grid::unit(-36, "pt"),
        gp = grid::gpar(fontfamily = "Open Sans", fontsize = 7)))
}
p9d <- candidate_plot("INSM1", "MYC", "Lung cancer")
p9e <- candidate_plot("ATRX", "CDKN2A", "Melanoma")
p9c <- p9c + labs(x = "Shift in LOY lines (Chronos units)") +
    theme(axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))
plot <- ((p9a + p9b + plot_layout(widths = c(1, 2.1))) /
         (free(p9c) + (p9d / p9e) + plot_layout(widths = c(1.6, 1))) +
         plot_layout(heights = c(1, 2.4)) + plot_annotation(tag_levels = "a")) &
        theme(plot.tag = element_text(size = 10, face = "bold"))
width <- 9.2
height <- 9.0
ggsave(file.path(fig_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite::svglite)
ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, device = ragg::agg_png)
node <- file.path(Sys.getenv("USERPROFILE"), ".cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node.exe")
status <- system2(node, shQuote(c(file.path(fig_dir, "code/depmap_svg_to_pdf.mjs"),
    file.path(fig_dir, paste0(stem, ".svg")), file.path(fig_dir, paste0(stem, ".pdf")),
    file.path(style_dir, "assets/fonts"))))
stopifnot(status == 0)
jsonlite::write_json(list(source_export = "2026-10-07", source_figure = "fig9_power_clean_shift",
    changes = "Former S07 INSM1 and ATRX boxplots added as panels d/e; interaction annotations remain omitted. Compact five-panel layout; abbreviated cancer labels in b.",
    width_inches = width, height_inches = height, font = "Open Sans", models_refitted = FALSE,
    saved_panel_data_match = TRUE, eligible_counts = as.integer(enough), simulated_groups = nrow(tested_types),
    groups_with_80_percent_detection_at_half_unit = half_unit,
    panel_b_label_rule = "All groups with at least five lines per Y group; largest 4, 6 and 8 below-threshold groups at lineage, disease and subtype levels.",
    panel_b_labels = as.list(table(labelled$level)),
    panels_de_observations = nrow(candidate_points), panels_de_models_refitted = FALSE),
    file.path(data_dir, "provenance.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Rendered five-panel F03; panels a-c match the frozen export and d-e use the 106 saved S07 observations.\n")
