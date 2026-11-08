options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
outdir <- "integration/feature_scatters"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
links <- fread("integration/feature_overlaps/peak_gene_link_audit.csv",
               colClasses = list(character = "entrez_id"))[mapping_status == "eligible"]
universe <- fread("integration/feature_overlaps/feature_gene_universes.csv",
                  colClasses = list(character = "entrez_id"))
expected <- fread("integration/feature_overlaps/overlap_counts.csv")
effects <- links[, .(atac_median_logFC = median(atac_logFC),
  verified_minimum_FDR = min(atac_FDR), verified_peak_count = uniqueN(peak_id)),
  by = .(feature, entrez_id)]
base <- merge(universe, effects, by = c("feature", "entrez_id"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(base) == nrow(universe), !anyDuplicated(base[, .(feature, entrez_id)]),
  all(is.finite(base$atac_median_logFC)), all(is.finite(base$rna_logFC)),
  all(base$minimum_peak_FDR == base$verified_minimum_FDR),
  all(base$tested_peaks == base$verified_peak_count))
base[, c("verified_minimum_FDR", "verified_peak_count") := NULL]
fwrite(base, file.path(outdir, "feature_gene_effects.csv"))
correlations <- base[, .(eligible_genes = .N,
  spearman_all_genes = cor(atac_median_logFC, rna_logFC, method = "spearman")), by = feature]
fwrite(correlations, file.path(outdir, "descriptive_correlations.csv"))
panels <- list(); summary <- list()
for (i in seq_len(nrow(expected))) {
  e <- expected[i]
  panel <- copy(base[feature == e$feature])
  panel[, `:=`(scheme = e$scheme, ATAC_FDR_cutoff = e$ATAC_FDR, RNA_FDR_cutoff = e$RNA_FDR)]
  panel[, membership := fcase(minimum_peak_FDR < e$ATAC_FDR & rna_FDR < e$RNA_FDR, "Both",
    minimum_peak_FDR < e$ATAC_FDR, "ATAC only", rna_FDR < e$RNA_FDR, "RNA only", default = "Neither")]
  observed <- panel[, .(universe = .N, ATAC_only = sum(membership == "ATAC only"),
    overlap = sum(membership == "Both"), RNA_only = sum(membership == "RNA only"), neither = sum(membership == "Neither"))]
  for (column in names(observed)) stopifnot(observed[[column]] == e[[column]])
  summary[[i]] <- cbind(e[, .(scheme, feature, ATAC_FDR, RNA_FDR)], observed,
    correlations[feature == e$feature, .(spearman_all_genes)])
  panels[[i]] <- panel
}
plot_data <- rbindlist(panels)
stats <- rbindlist(summary)
fwrite(stats, file.path(outdir, "panel_summary.csv"))
fwrite(plot_data[, .(scheme, feature, ATAC_FDR_cutoff, RNA_FDR_cutoff, entrez_id, membership)],
       file.path(outdir, "panel_gene_membership.csv"))
feature_labels <- c(promoter = "Promoter\nTSS overlap", gene_body = "Gene body\nGene-body overlap",
                    intergenic = "Intergenic\nNearest TSS (exploratory)")
palette <- c("Neither" = "#C7CDD1", "ATAC only" = "#0072B2", "RNA only" = "#D55E00", "Both" = "#7B3294")
x_range <- range(base$atac_median_logFC)
y_range <- range(base$rna_logFC)
# The same axes show every eligible gene across all feature and threshold panels.
x_limits <- x_range + c(-1, 1)*diff(x_range)*.04
y_limits <- y_range + c(-1, 1)*diff(y_range)*.12
make_plot <- function(df, panel_stats, heading) {
  df <- copy(df); panel_stats <- copy(panel_stats)
  df[, feature_label := factor(feature, levels = names(feature_labels), labels = feature_labels)]
  panel_stats[, feature_label := factor(feature, levels = names(feature_labels), labels = feature_labels)]
  threshold_levels <- unique(sprintf("ATAC <%.2f | RNA <%.2f", panel_stats$ATAC_FDR, panel_stats$RNA_FDR))
  df[, threshold := factor(sprintf("ATAC <%.2f | RNA <%.2f", ATAC_FDR_cutoff, RNA_FDR_cutoff), levels = threshold_levels)]
  panel_stats[, threshold := factor(sprintf("ATAC <%.2f | RNA <%.2f", ATAC_FDR, RNA_FDR), levels = threshold_levels)]
  df[, membership := factor(membership, levels = names(palette))]
  panel_stats[, label := sprintf("n = %s | Both = %s\nAll-gene Spearman rho = %.3f",
    format(universe, big.mark = ",", trim = TRUE), format(overlap, big.mark = ",", trim = TRUE), spearman_all_genes)]
  ggplot(df, aes(atac_median_logFC, rna_logFC)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = .35) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = .35) +
    geom_point(data = df[membership == "Neither"], aes(colour = membership), size = .55, alpha = .22) +
    geom_point(data = df[membership %in% c("ATAC only", "RNA only")], aes(colour = membership), size = .75, alpha = .55) +
    geom_point(data = df[membership == "Both"], aes(colour = membership), size = 1.35, alpha = .9) +
    geom_text(data = panel_stats, aes(x = x_limits[1], y = y_limits[2], label = label),
      inherit.aes = FALSE, hjust = -.04, vjust = 1.15, size = 3.2, lineheight = 1) +
    scale_colour_manual(values = palette, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    facet_grid(feature_label ~ threshold) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.title = element_blank(),
      panel.grid.minor = element_blank(), strip.text = element_text(size = 10),
      strip.background = element_rect(fill = "#edf3f5", colour = NA),
      plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 10, hjust = 0), panel.spacing = grid::unit(.9, "lines")) +
    labs(title = heading, subtitle = "One point per eligible gene; colors match the corresponding overlap sets",
      x = "Median ATAC log2 fold change across linked peaks in this feature (LOY / CTRL)",
      y = "RNA log2 fold change (LOY / CTRL)",
      caption = paste("ATAC membership: any linked peak passes FDR. Coordinates use all tested linked peaks and stay fixed across thresholds.",
        "ATAC FDR is peak-level; median effects and correlations are descriptive. Intergenic links are nearest-TSS hypotheses.", sep = "\n"))
}
for (scheme_name in unique(stats$scheme)) {
  df <- plot_data[scheme == scheme_name]
  ps <- stats[scheme == scheme_name]
  heading <- if (scheme_name == "matched_FDR") "ATAC versus RNA: matched FDR thresholds" else
    "ATAC versus RNA: RNA FDR fixed at 0.05"
  p <- make_plot(df, ps, heading)
  stem <- paste0(scheme_name, "_scatter_grid")
  ggsave(file.path(outdir, paste0(stem, ".png")), p, width = 14, height = 11, dpi = 220, bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = 14, height = 11, bg = "white", useDingbats = FALSE)
  for (feature_name in names(feature_labels)) {
    p <- make_plot(df[feature == feature_name], ps[feature == feature_name], heading)
    stem <- paste0(scheme_name, "_", feature_name, "_scatters")
    ggsave(file.path(outdir, paste0(stem, ".png")), p, width = 14, height = 5.5, dpi = 220, bg = "white")
    ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = 14, height = 5.5, bg = "white", useDingbats = FALSE)
  }
}
writeLines(c("Scatter plots use exactly the feature-specific eligible gene universes and links used by code/15_feature_overlap_plots.R.",
  "One point per distinct Entrez gene per feature. X: median log2FC across all tested linked ATAC peaks in that feature. Y: existing batch+treatment RNA log2FC.",
  "The same gene coordinates and full axis limits are used at every FDR threshold; no point is dropped for display.",
  "Colors encode Neither, ATAC only, RNA only, or Both. ATAC membership uses any linked peak FDR < cutoff; RNA uses its original gene FDR < cutoff.",
  "Two schemes: matched FDR 0.05/0.10/0.20, and RNA fixed at 0.05 with those ATAC cutoffs. All 18 panel counts exactly reproduce the overlap figures.",
  "Positive effects mean higher in LOY for both assays. No effect-size or direction filter is applied.",
  "A gene's median ATAC effect can differ in sign from an individual significant peak; median effects have no gene-level P-value or FDR.",
  "Spearman correlations use all eligible genes in each feature, not only significant genes; they are descriptive group-effect comparisons, not sample-wise associations.",
  "Genes can appear in multiple feature classes. Intergenic nearest-TSS links remain positional hypotheses.",
  "Exports: two grid PNG/PDF pairs, six individual-feature PNG/PDF pairs, gene effect coordinates, panel membership, summary and verification."),
  file.path(outdir, "METHODS.txt"))
writeLines(c("PASS: unique feature/gene coordinates, finite effects, and identical eligible universes",
  "PASS: ATAC minimum FDR and linked-peak counts reproduce source link data",
  "PASS: all 18 scatter membership partitions exactly match prior overlap counts",
  "PASS: shared axis limits include all gene coordinates; coordinates do not change with FDR"),
  file.path(outdir, "verification.txt"))
print(correlations)
print(stats[, .(scheme, feature, ATAC_FDR, RNA_FDR, universe, overlap)])
