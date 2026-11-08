# Stage 21: figures. House style from plot_style.R; every figure is written as
# PDF, SVG and 300 dpi PNG, and the data behind each panel goes to
# results/loy_reanalysis_v1/figure_source_data/. Titles state the finding and
# captions state the method, so each title is computed from the numbers it claims.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "plot_style.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(ggrepel)
})
assert_freeze()

core <- params$y_genes_core
calls <- read_result("y_status_calls")
cohort <- read_result("cohort")
status_label <- function(y_negative) factor(fifelse(is.na(y_negative), "set aside", fifelse(y_negative == 1, "Y-negative", "Y-positive")), levels = c("Y-positive", "Y-negative", "set aside"))
status_colours <- c(y_status_colours, "set aside" = unname(accent_colours["other"]))
red <- unname(y_status_colours["Y-negative"])
blue <- unname(y_status_colours["Y-positive"])
grey <- unname(accent_colours["other"])
tagged <- function(plot) plot + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 10, face = "bold"))
fmt_p <- function(p) formatC(p, format = "e", digits = 1)
full <- 9.2
power_multiplier <- qnorm(0.975) + qnorm(0.8)

# --- Figures 11 and 12: organoids, RNAi and primary AML (extension stages 25 to 30) -----
# Defined here so that `21_figures.R extension` draws only these two and leaves the first
# ten untouched; a full run draws them last.
draw_extension_figures <- function() {
  ext_dir <- file.path(derived_dir, "extension")
  read_ext <- function(name) read_result(name, dir = ext_dir)
  fdr_cut <- params$gene_scan$fdr
  teal <- unname(accent_colours["x_partner"])

  # Figure 11a: models per Y group in the two cancer types, by source.
  power <- read_ext("extension_power")
  sizes <- melt(power[analysis != "Pan-cancer", .(analysis, source, `Y-negative` = n_negative, `Y-positive` = n_positive)], id.vars = c("analysis", "source"), variable.name = "status", value.name = "models")
  sizes[, source := factor(source, levels = c("Cell lines", "Organoids", "Combined"))]
  sizes[, analysis := sub(" adenocarcinoma", "", analysis)]
  save_source_data(sizes, "fig11a_models_per_source")
  gain <- power[analysis != "Pan-cancer" & source != "Combined", .(models = sum(n_negative + n_positive)), by = source]
  p11a <- ggplot(sizes, aes(source, models, fill = status)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = models), position = position_dodge(width = 0.75), vjust = -0.3, size = 2.1) +
    facet_wrap(~analysis, nrow = 1) +
    scale_fill_manual(values = y_status_colours, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    panel_labs(2.9, title = sprintf("Organoids add %d models to %d cell lines in the two cancer types", gain[source == "Organoids", models], gain[source == "Cell lines", models]),
      x = NULL, y = "Male models with a Y call",
      caption = "CRISPR-screened models per Y group. Combined: both sources, as used in the within-type combination.") +
    theme_house() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")

  # Figure 11b: organoid Y calls, RNA against DNA.
  organoids <- read_ext("organoid_cohort")[male == TRUE]
  organoids[, status := status_label(y_negative)]
  save_source_data(organoids[, .(ModelID, OncotreeLineage, rna_y_mean, dna_y_ratio, rna_call, dna_call, y_tier, status)], "fig11b_organoid_rna_vs_dna")
  p11b <- ggplot(organoids, aes(dna_y_ratio, rna_y_mean, colour = status)) +
    geom_point(size = 0.9, alpha = 0.8) +
    scale_colour_manual(values = status_colours, name = NULL) +
    panel_labs(2.9, title = sprintf("RNA and DNA agree for %d of %d screened male organoids", organoids[y_tier == "concordant", .N], nrow(organoids)),
      x = "chrY copy number relative to ploidy", y = "Mean log2(TPM+1) of seven Y genes",
      caption = sprintf("%d set aside: %d with intermediate RNA, %d where DNA contradicts RNA. RNA leads where DNA is unclear: %d Y-negative points at the right express no Y gene but keep most of the male-specific region in DNA, and %d Y-positive points near zero copy number express all seven Y genes.",
        organoids[is.na(y_negative), .N], organoids[y_tier == "rna_intermediate", .N], organoids[y_tier == "discordant", .N], organoids[y_negative %in% 1 & dna_y_retained > 0.5, .N],
        organoids[y_tier == "rna_led_dna_ambiguous" & rna_call == "expressed" & dna_y_ratio < params$dna_call$retained_ratio_min, .N])) +
    theme_house() + theme(legend.position = "top")

  # Figure 11c: the scan over all organoids.
  scan <- read_ext("organoid_pan_results")
  scan[, class := fifelse(associated, "hit", fifelse(!is.na(y_partner), "X-Y partner, not a hit", "other"))]
  save_source_data(scan[, .(symbol, chromosome, estimate, estimate_chronos_equivalent, p, fdr, label, class)], "fig11c_organoid_volcano")
  p11c <- ggplot(scan[order(class != "other")], aes(estimate_chronos_equivalent, -log10(p), colour = class)) +
    geom_point(aes(size = class), alpha = 0.8) +
    geom_text_repel(data = scan[associated == TRUE], aes(label = symbol), size = 2.3, min.segment.length = 0, show.legend = FALSE, seed = 11) +
    geom_text_repel(data = scan[symbol %in% c("ZFX", "CORO1A")], aes(label = symbol), colour = "#222222", size = 2.3, min.segment.length = 0, box.padding = 0.5, show.legend = FALSE, seed = 11) +
    scale_colour_manual(values = c(hit = red, "X-Y partner, not a hit" = teal, other = grey), name = NULL) +
    scale_size_manual(values = c(hit = 1.6, "X-Y partner, not a hit" = 1.2, other = 0.3), guide = "none") +
    panel_labs(2.9, title = sprintf("%d genes pass in organoids: %s", sum(scan$associated), paste(scan[associated == TRUE, symbol], collapse = " and ")),
      x = "Effect, Y-negative minus Y-positive\n(Chronos-equivalent units)", y = "-log10 p",
      caption = sprintf("%s genes, %d organoids; limma with a weight per organoid; tumour type, library and MSI adjusted.", format(nrow(scan), big.mark = ","), max(scan$n))) +
    theme_house() + theme(legend.position = "top")

  # Figure 11d: the cell-line hits in every source.
  partners <- read_ext("combined_xy_partners")[symbol %in% c("DDX3X", "EIF1AX", "ZFX", "RPS4X")]
  partners[, source := factor(source, levels = rev(c("Cell lines, CRISPR", "Cell lines, RNAi", "Organoids, CRISPR", "CRISPR combined")))]
  partners[, symbol := factor(symbol, levels = c("DDX3X", "EIF1AX", "ZFX", "RPS4X"))]
  save_source_data(partners[, .(symbol, source, smd, smd_low, smd_high, p, n_negative, n_positive)], "fig11d_hits_across_sources")
  p11d <- ggplot(partners, aes(smd, source, colour = source)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888888") +
    geom_linerange(aes(xmin = smd_low, xmax = smd_high), linewidth = 0.4) +
    geom_point(size = 1.4) +
    facet_wrap(~symbol, nrow = 1) +
    scale_colour_manual(values = c("Cell lines, CRISPR" = "#222222", "Cell lines, RNAi" = "#888888", "Organoids, CRISPR" = teal, "CRISPR combined" = red), guide = "none") +
    panel_labs(5.2, title = "DDX3X and RPS4X hold in every source; ZFX does not in organoids, and EIF1AX is not in their library",
      x = "Standardised effect, Y-negative minus Y-positive (95% interval)", y = NULL,
      caption = "Effect in residual SD units of each source, from its own p-value and design. Combined: cell lines and organoids, weighted by design.") +
    theme_house()

  # Figure 11e: the combination within cancer type.
  within <- read_ext("combined_within_type")
  within[, pair := sub(" adenocarcinoma", "", pair)]
  within[, expected := -log10((frank(combined_p) - 0.5) / .N), by = pair]
  save_source_data(within[, .(pair, symbol, cell_inside_p, organoid_inside_p, combined_p, combined_fdr, evidence)], "fig11e_combined_within_type")
  top <- within[, .SD[order(combined_p)][1:3], by = pair]
  p11e <- ggplot(within, aes(expected, -log10(combined_p))) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#888888") +
    geom_point(aes(colour = combined_hit), size = 0.4, alpha = 0.7) +
    geom_text_repel(data = top, aes(label = symbol), size = 2.2, min.segment.length = 0, seed = 11) +
    facet_wrap(~pair, nrow = 1) +
    scale_colour_manual(values = c(`TRUE` = red, `FALSE` = "#555555"), guide = "none") +
    panel_labs(3.8, title = if (any(within$combined_hit)) {
        sprintf("Within cancer type, the combination finds %s", paste(within[combined_hit == TRUE, sprintf("%s (%s)", symbol, tolower(pair))], collapse = ", "))
      } else {
        sprintf("Within cancer type, the combination finds no gene; %s in %s cancer comes closest (FDR %.2f)", within[which.min(combined_p), symbol], tolower(within[which.min(combined_p), pair]), within[which.min(combined_p), combined_fdr])
      },
      x = "Expected -log10 p", y = "Observed -log10 p, combined",
      caption = sprintf("%s genes per cancer type; rank-based test of the Y effect inside the type in each source, combined by design weight. Red: FDR < %.2f.", format(within[, .N, by = pair]$N[1], big.mark = ","), fdr_cut)) +
    theme_house()
  figure11 <- (p11a + p11b + p11c + plot_layout(widths = c(1, 1, 1))) / (p11d + p11e + plot_layout(widths = c(1.35, 1)))
  save_figure(tagged(figure11), "fig11_more_models", width = full, height = 7.4)

  # Figure 12a: Y-loss fraction of the male specimens.
  aml <- read_ext("beataml_cohort")[in_primary == TRUE]
  aml[, tier := factor(evidence_tier, levels = c("karyotype -Y", "karyotype XY", "no karyotype"), labels = c("-Y in the karyotype", "XY karyotype", "No karyotype"))]
  save_source_data(aml[, .(SpecimenID, f, f_no_floor, tier, karyotype_y_fraction, stage)], "fig12a_y_loss_fraction")
  floor_value <- aml$f_floor[1]
  tier_colours <- c("-Y in the karyotype" = red, "XY karyotype" = blue, "No karyotype" = grey)
  p12a <- ggplot(aml, aes(f_no_floor, fill = tier)) +
    geom_histogram(binwidth = 0.05, boundary = 0, colour = "white", linewidth = 0.1) +
    geom_vline(xintercept = floor_value, linetype = "dashed", linewidth = 0.3) +
    scale_fill_manual(values = tier_colours, name = NULL) +
    scale_y_sqrt(breaks = c(1, 5, 20, 50, 100)) +
    panel_labs(4.3, title = sprintf("%d of %d male patients have a Y-loss fraction of at least 0.5; %d of them show -Y in the karyotype", aml[f >= 0.5, .N], nrow(aml), aml[f >= 0.5 & evidence_tier == "karyotype -Y", .N]),
      x = "Y-loss fraction from RNA (before the floor)", y = "Patients (square-root scale)",
      caption = sprintf("One specimen per patient. Seven Y genes, each relative to its X partner and to the median male specimen. Dashed: floor at %.2f, below which the fraction is set to 0.", floor_value)) +
    theme_house() + theme(legend.position = c(0.6, 0.8))

  # Figure 12b: RNA fraction against karyotype clone size.
  clones <- read_ext("beataml_cohort")[eligible == TRUE & karyotype_mentions_y_loss == TRUE & !is.na(karyotype_y_fraction)]
  save_source_data(clones[, .(SpecimenID, f_no_floor, karyotype_y_fraction, stage)], "fig12b_rna_vs_karyotype")
  p12b <- ggplot(clones, aes(karyotype_y_fraction, f_no_floor)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#888888") +
    geom_point(size = 1.2, colour = red, alpha = 0.8) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    panel_labs(4.3, title = sprintf("The RNA fraction follows the karyotype clone size in %d specimens (rank correlation %.2f)", nrow(clones), cor(clones$f_no_floor, clones$karyotype_y_fraction, method = "spearman")),
      x = "Cells in clones without a Y, karyotype", y = "Y-loss fraction from RNA",
      caption = "Male specimens whose karyotype mentions -Y. Line: equality. The karyotype counts dividing cells only.") +
    theme_house()

  # Figure 12c: the inhibitors.
  drugs <- read_ext("beataml_drug_results")
  save_source_data(drugs[, .(inhibitor, n, high_loss, estimate, slope_sd, p, fdr, label)], "fig12c_inhibitor_volcano")
  p12c <- ggplot(drugs, aes(estimate, -log10(p), colour = associated)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_text_repel(data = drugs[order(p)][1:4], aes(label = sub(" \\(.*$", "", inhibitor)), size = 2.2, min.segment.length = 0, show.legend = FALSE, seed = 11) +
    scale_colour_manual(values = c(`TRUE` = red, `FALSE` = "#555555"), guide = "none") +
    panel_labs(4.3, title = if (sum(drugs$associated) == 1) {
      sprintf("1 of %d inhibitors passes: %s, %s", nrow(drugs), sub(" \\(.*$", "", drugs[associated == TRUE, inhibitor]), drugs[associated == TRUE, direction])
    } else {
      sprintf("%d of %d inhibitors pass; the leading one is %s", sum(drugs$associated), nrow(drugs), sub(" \\(.*$", "", drugs[which.min(p), inhibitor]))
    },
      x = "Slope: rank-normalised AUC per complete Y loss (positive = less sensitive)", y = "-log10 permutation p",
      caption = sprintf("Up to %d patients per inhibitor; specimen type, disease stage, age, core-binding-factor fusion and study wave adjusted; %s Freedman-Lane permutations. Red: FDR < %.2f.", max(drugs$n), format(1e5, big.mark = ",", scientific = FALSE), fdr_cut)) +
    theme_house()

  # Figure 12d: the leading inhibitor, specimen by specimen.
  leading <- read_ext("beataml_leading_inhibitor")
  leading[, ras := factor(fifelse(is.na(ras_pathway), "unknown", fifelse(ras_pathway == 1, "RAS-pathway mutant", "RAS-pathway wild type")), levels = c("RAS-pathway wild type", "RAS-pathway mutant", "unknown"))]
  save_source_data(leading[, .(inhibitor, SpecimenID, f, auc, ras, evidence_tier)], "fig12d_leading_inhibitor")
  name <- sub(" \\(.*$", "", leading$inhibitor[1])
  row <- drugs[inhibitor == leading$inhibitor[1]]
  followup <- read_ext("beataml_hit_followup")
  ras_slope <- followup[variant %like% "plus RAS-pathway mutation"]
  p12d <- ggplot(leading, aes(f, auc, colour = ras)) +
    geom_jitter(width = 0.012, height = 0, size = 0.9, alpha = 0.75) +
    scale_colour_manual(values = c("RAS-pathway wild type" = "#555555", "RAS-pathway mutant" = teal, unknown = grey), name = NULL) +
    panel_labs(4.3, title = if (leading[f > 0 & ras_pathway %in% 1, .N] == 0 && leading[ras_pathway %in% 1, median(auc)] < leading[ras_pathway %in% 0, median(auc)]) {
      sprintf("%s: no specimen with Y loss carries a RAS-pathway mutation, the group that responds best", name)
    } else {
      sprintf("%s response against the Y-loss fraction, by RAS-pathway status", name)
    },
      x = "Y-loss fraction", y = sprintf("%s AUC (lower = more sensitive)", name),
      caption = sprintf("%d patients. Slope %.2f SD per complete Y loss (permutation p %s); %.2f with RAS-pathway mutation in the model (p %.3f), a check added after the result.", nrow(leading), row$estimate, fmt_p(row$p), ras_slope$estimate, ras_slope$p)) +
    theme_house() + theme(legend.position = "top")
  save_figure(tagged((p12a + p12b) / (p12c + p12d)), "fig12_primary_aml", width = full, height = 7.6)
}

# --- Figure 13: DepMap NextGen, a third CRISPR source (extension stages 31 to 33) --------
# `21_figures.R nextgen` draws only this figure, so the twelve before it keep their bytes.
# Label placement is seeded and nothing is jittered.
draw_nextgen_figure <- function() {
  ext_dir <- file.path(derived_dir, "extension")
  read_ext <- function(name) read_result(name, dir = ext_dir)
  read_summary <- function(name) jsonlite::fromJSON(file.path(ext_dir, name))
  fdr_cut <- params$gene_scan$fdr
  teal <- unname(accent_colours["x_partner"])
  purple <- unname(accent_colours["par"])
  cohort <- read_ext("nextgen_cohort")
  primary <- cohort[in_primary == TRUE]
  paired <- read_ext("nextgen_paired_genes")
  prepared <- read_summary("nextgen_prepare_summary.json")
  scanned <- read_summary("nextgen_scan_summary.json")
  listed <- function(x) if (length(x) < 2) paste(x) else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])

  # Figure 13a: models per Y group by cancer type in the three sources.
  counts <- read_ext("nextgen_pooled_types")
  shown <- primary[, .N, by = OncotreeLineage][N >= 5][order(-N), OncotreeLineage]
  sizes <- melt(counts[type %in% shown, .(type, source, `Y-negative` = n_negative, `Y-positive` = n_positive)], id.vars = c("type", "source"), variable.name = "status", value.name = "models")
  sizes[, source := factor(source, levels = c("Cell lines", "Sanger organoids", "NextGen"))]
  sizes[, type := factor(type, levels = shown)]
  save_source_data(sizes, "fig13a_models_by_cancer_type")
  # Two-line strip labels, so the longest lineage name fits its facet.
  levels(sizes$type) <- sub("/", "/\n", levels(sizes$type))
  pooled_types <- counts[pooled == TRUE, unique(type)]
  in_pooled <- primary[OncotreeLineage %in% pooled_types]
  p13a <- ggplot(sizes, aes(source, models, fill = status)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = models), position = position_dodge(width = 0.75), vjust = -0.3, size = 2.1) +
    facet_wrap(~type, nrow = 1) +
    scale_fill_manual(values = y_status_colours, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    panel_labs(5.4, title = sprintf("NextGen adds %d male models (%d Y-negative); %s is the one cancer type that only the sources together can test", nrow(primary), sum(primary$y_negative == 1), tolower(listed(pooled_types))),
      x = NULL, y = "Male models with a Y call",
      caption = sprintf("CRISPR-screened models per source, for the cancer types with at least five NextGen models. A test inside a cancer type needs five models per Y group. The %d Y-positive pancreatic cell lines sit in one screening library with %d of the Y-negative ones, so those %d lines are what the pooled test can use.",
        prepared$pooled_floor$n_positive - sum(in_pooled$y_negative == 0), prepared$pooled_floor$n_negative - sum(in_pooled$y_negative == 1), prepared$pooled_floor$n_positive + prepared$pooled_floor$n_negative - nrow(in_pooled))) +
    theme_house() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top", strip.text = element_text(size = 7.5, face = "bold"))

  # Figure 13b: NextGen Y calls, RNA against DNA.
  males <- cohort[eligible == TRUE]
  males[, status := status_label(y_negative)]
  save_source_data(males[, .(ModelID, name, OncotreeLineage, ScreenType, rna_y_mean, dna_y_ratio, rna_call, dna_call, y_tier, status)], "fig13b_nextgen_rna_vs_dna")
  p13b <- ggplot(males[!is.na(rna_y_mean) & !is.na(dna_y_ratio)], aes(dna_y_ratio, rna_y_mean, colour = status)) +
    geom_point(size = 0.9, alpha = 0.8) +
    scale_colour_manual(values = status_colours, name = NULL) +
    panel_labs(3.3, title = sprintf("RNA and DNA agree for %d of %d eligible male models", males[y_tier == "concordant", .N], nrow(males)),
      x = "chrY copy number relative to ploidy\n(mean of seven Y genes)", y = "Mean log2(TPM+1),\nseven Y genes",
      caption = sprintf("The %d screened males without the non-malignant models and those shared with the Sanger cohort. %d of the %d Y-negative points lie on top of each other at zero. %d set aside with intermediate RNA; RNA leads in %d where DNA is unclear; %d without RNA keeps its DNA call and is not drawn. \"No Y expression\" means at most the female maximum or %.0f%% of the typical male level.",
        cohort[recorded_male == TRUE, .N], males[y_negative %in% 1 & rna_y_mean < 0.05 & dna_y_ratio < 0.01, .N], males[y_negative %in% 1, .N],
        males[y_tier == "rna_intermediate", .N], males[y_tier == "rna_led_dna_ambiguous", .N], males[y_tier == "dna_only", .N], 100 * prepared$rna_thresholds$ceiling_fraction)) +
    theme_house() + theme(legend.position = "top")

  # Figure 13c: the scan over all NextGen models.
  scan <- read_ext("nextgen_pan_results")
  classes <- c("paired gene, named in advance", sprintf("FDR < %.2f, does not clear the bar", fdr_cut), "other")
  scan[, class := factor(fifelse(named, classes[1], fifelse(associated & !discovery, classes[2], fifelse(discovery, "discovery", classes[3]))), levels = c(classes[1:2], "discovery", classes[3]))]
  save_source_data(scan[, .(symbol, chromosome, estimate, p, fdr, named, associated, passes_family_wise_bar, discovery, label)], "fig13c_nextgen_volcano")
  others <- scan[associated == TRUE & named == FALSE]
  zfx <- paired[symbol == "ZFX"]
  # The verdict on the permutation p-values: the model's p-values run small in this scan.
  confirmed <- paired[outcome_permutation == "replicated", symbol]
  unconfirmed <- setdiff(paired$symbol, confirmed)
  p13c <- ggplot(scan[order(class == "other", decreasing = TRUE)], aes(estimate, -log10(p), colour = class)) +
    geom_hline(yintercept = -log10(scanned$null$pan_bar), linetype = "dashed", linewidth = 0.3, colour = "#555555") +
    geom_point(aes(size = class), alpha = 0.8) +
    geom_text_repel(data = scan[named == TRUE | associated == TRUE], aes(label = symbol), size = 2.3, min.segment.length = 0, show.legend = FALSE, seed = 13) +
    scale_colour_manual(values = setNames(c(teal, "#222222", red, grey), levels(scan$class)), name = NULL, breaks = intersect(levels(scan$class), as.character(unique(scan$class)))) +
    scale_size_manual(values = setNames(c(1.6, 1.4, 1.6, 0.3), levels(scan$class)), guide = "none") +
    panel_labs(3.3, title = sprintf("%s; %s", if (length(unconfirmed) == 0) sprintf("All %d paired genes replicate in NextGen", nrow(paired)) else
          sprintf("%s replicate in NextGen and %s %s not confirmed", listed(confirmed), listed(unconfirmed), if (length(unconfirmed) == 1) "is" else "are"),
        if (sum(scan$discovery) > 0) {
          sprintf("%d more clear the bar", sum(scan$discovery))
        } else if (nrow(others) == 0) {
          sprintf("no other gene passes FDR %.2f", fdr_cut)
        } else {
          sprintf("%d more pass FDR %.2f, and %s the bar", nrow(others), fdr_cut, if (nrow(others) == 1) "it does not clear" else if (nrow(others) == 2) "neither clears" else "none clears")
        }),
      x = "Effect, Y-negative minus Y-positive\n(Chronos units)", y = "-log10 p",
      caption = sprintf("%s genes, %d models; limma with a weight per model; one term per lineage and screen format. Green: the four paired genes, named in advance. Black: FDR < %.2f. Dashed: the bar from %d shuffles of Y status (p < %s) that a gene outside the four must clear. The four are judged on permutation p-values (%s shuffles): ZFX has p = %.3f in this model and %.3f by permutation.",
        format(nrow(scan), big.mark = ","), max(scan$n), fdr_cut, scanned$null$runs, fmt_p(scanned$null$pan_bar), format(zfx$permutation_shuffles, big.mark = ",", scientific = FALSE), zfx$p, zfx$permutation_p)) +
    theme_house() + theme(legend.position = "none")

  # Figure 13d: the four paired genes in every source.
  three <- read_ext("nextgen_combined_pan")
  sources <- c("Cell lines, CRISPR", "Cell lines, RNAi", "Sanger organoids, CRISPR", "NextGen, CRISPR", "CRISPR combined, three sources")
  partners <- read_ext("nextgen_xy_partners")[symbol %in% paired$symbol & source %in% sources]
  partners[, source := factor(source, levels = rev(sources))]
  partners[, symbol := factor(symbol, levels = paired$symbol)]
  save_source_data(partners[, .(symbol, source, smd, smd_low, smd_high, p, n_negative, n_positive)], "fig13d_paired_genes_across_sources")
  every <- paired$symbol[paired$symbol %in% three[evidence == "supported by every source", symbol]]
  rest <- setdiff(paired$symbol, every)
  p13d <- ggplot(partners, aes(smd, source, colour = source)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888888") +
    geom_linerange(aes(xmin = smd_low, xmax = smd_high), linewidth = 0.4) +
    geom_point(size = 1.4) +
    facet_wrap(~symbol, nrow = 1) +
    scale_colour_manual(values = setNames(c("#222222", "#888888", teal, purple, red), sources), guide = "none") +
    panel_labs(5.4, title = if (length(rest) == 0) {
        "All four paired genes are supported by every CRISPR source that tests them"
      } else {
        sprintf("%s %s supported by every CRISPR source that tests them; %s by %s", listed(every), if (length(every) == 1) "is" else "are", listed(rest), paste(unique(sub("supported by ", "", three[symbol %in% rest, evidence])), collapse = " and "))
      },
      x = "Standardised effect, Y-negative minus Y-positive (95% interval)", y = NULL,
      caption = "Effect in residual SD units of each source, from its own p-value and design. Combined: cell lines, Sanger organoids and NextGen, weighted by design. A source without a point does not test the gene (EIF1AX is not in the Sanger library).") +
    theme_house()

  # Figure 13e: pancreas, the paired genes with the sources pooled.
  by_source <- read_ext("nextgen_pooled_named")
  by_source[, source := factor(source, levels = rev(c("Cell lines", "NextGen", "Pooled")))]
  by_source[, symbol := factor(symbol, levels = rev(paired$symbol))]
  by_source[, p_label := fifelse(source == "Pooled", sprintf("p = %s", formatC(p, format = "g", digits = 2)), "")]
  save_source_data(by_source[, .(type, source, symbol, estimate_chronos, estimate_rank, p, named_fdr, floor, n_negative, n_positive)], "fig13e_pancreas_paired_genes")
  pooled_rows <- by_source[source == "Pooled"]
  holds <- paired$symbol[paired$symbol %in% as.character(pooled_rows[named_fdr < fdr_cut, symbol])]
  fails <- setdiff(paired$symbol, holds)
  size_of <- function(one) by_source[source == one][1, sprintf("%d against %d", n_negative, n_positive)]
  p13e <- ggplot(by_source, aes(estimate_chronos, symbol, colour = source)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888888") +
    geom_point(aes(size = source), position = position_dodge(width = 0.6)) +
    geom_text(aes(label = p_label), position = position_dodge(width = 0.6), hjust = 1.2, size = 2.1, show.legend = FALSE) +
    scale_colour_manual(values = c("Cell lines" = "#222222", NextGen = purple, Pooled = red), name = NULL, breaks = c("Cell lines", "NextGen", "Pooled")) +
    scale_size_manual(values = c("Cell lines" = 1.2, NextGen = 1.2, Pooled = 2), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.4, 0.08))) +
    panel_labs(4.4, title = sprintf("Inside %s, %s show the Y effect when the sources are pooled%s", tolower(listed(unique(by_source$type))), listed(holds),
        if (length(fails) > 0) sprintf("; %s %s not", listed(fails), if (length(fails) == 1) "does" else "do") else ""),
      x = "Difference in Chronos score within strata,\nY-negative minus Y-positive", y = NULL,
      caption = sprintf("Cell lines of one library (%s models), NextGen (%s) and both (%s). p-values from every arrangement of Y status within strata (%s when pooled), on ranks within source; label: the pooled p-value.",
        size_of("Cell lines"), size_of("NextGen"), size_of("Pooled"), format(pooled_rows$arrangements[1], big.mark = ","))) +
    theme_house() + theme(legend.position = "top")

  # Figure 13f: pancreas, every gene.
  everything <- read_ext("nextgen_pooled_type")
  everything[, expected := -log10((frank(p, ties.method = "first") - 0.5) / .N), by = type]
  combined <- read_summary("nextgen_combined_summary.json")
  save_source_data(everything[, .(type, symbol, named, estimate_rank, estimate_chronos, p, p_family_wise, fdr, floor)], "fig13f_pancreas_all_genes")
  one_type <- combined$pooled[[1]]
  p13f <- ggplot(everything, aes(expected, -log10(p))) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#888888") +
    geom_hline(yintercept = -log10(one_type$floor), linetype = "dashed", linewidth = 0.3, colour = "#555555") +
    geom_hline(yintercept = -log10(one_type$bh_critical_p), linetype = "dotted", linewidth = 0.3, colour = red) +
    geom_point(size = 0.3, alpha = 0.6, colour = "#555555") +
    geom_point(data = everything[named == TRUE], colour = teal, size = 1.4) +
    geom_text_repel(data = everything[named == TRUE], aes(label = symbol), colour = teal, size = 2.2, min.segment.length = 0, seed = 13) +
    panel_labs(4.4, title = sprintf("No single gene can be called genome-wide inside %s: the smallest possible p-value is %s, and one discovery would need %s", tolower(names(combined$pooled)[1]), fmt_p(one_type$floor), fmt_p(one_type$bh_critical_p)),
      x = "Expected -log10 p", y = "Observed -log10 p, pooled",
      caption = sprintf("%s genes with a score in every pooled model. Dashed: the floor of the exact test. Dotted: the p-value that one BH discovery at FDR %.2f would need. Smallest family-wise p over all genes: %.2f. The table is a ranking; no gene is called from it.",
        format(nrow(everything), big.mark = ","), fdr_cut, one_type$min_family_wise_p)) +
    theme_house()

  figure13 <- (p13a + p13b + plot_layout(widths = c(1.6, 1))) / (p13c + p13d + plot_layout(widths = c(1, 1.6))) / (p13e + p13f + plot_layout(widths = c(1, 1)))
  save_figure(tagged(figure13), "fig13_nextgen", width = full, height = 10.8)
}
# --- Figure 14: driver status inside cancer types (extension stages 34 to 36) ---------------
# Drawn alone by `21_figures.R driver_type`; a full run draws it last.
draw_driver_type_figure <- function() {
  ext_dir <- file.path(derived_dir, "extension")
  read_ext <- function(name) read_result(name, dir = ext_dir)
  read_summary <- function(name) jsonlite::fromJSON(file.path(ext_dir, name))
  teal <- unname(accent_colours["x_partner"])
  purple <- unname(accent_colours["par"])
  gate <- params$context$min_cell
  big <- function(x) format(x, big.mark = ",", trim = TRUE)
  listed <- function(x) if (length(x) < 2) paste(x) else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
  splits <- read_ext("driver_type_splits")[fitted == TRUE]
  tested <- read_summary("driver_type_tests_summary.json")
  powered <- read_summary("driver_type_power_summary.json")
  level_names <- c(lineage = "Lineage", disease = "Disease", subtype = "Subtype")

  # Figure 14a: how many gene-by-cancer-type pairs keep enough altered lines as the requirement rises.
  distance <- read_ext("driver_type_gate_distance")[source == "Cell lines"]
  step_lines <- c(both_y_groups = 1, at_least_2 = 2, at_least_3 = 3, at_least_5 = 5, at_least_10 = 10, at_least_20 = 20)
  steps <- melt(distance, id.vars = c("level", "pairs", "all_four_cells"), measure.vars = names(step_lines), variable.name = "step", value.name = "n_pairs")
  steps[, lines := step_lines[as.character(step)]]
  steps[, level_name := factor(level_names[level], levels = level_names)]
  save_source_data(steps[, .(level, pairs_with_an_altered_line = pairs, lines_per_y_group_at_least = lines, pairs_left = n_pairs, pairs_with_five_in_all_four_cells = all_four_cells)], "fig14a_pairs_by_lines")
  lineage <- distance[level == "lineage"]
  p14a <- ggplot(steps[n_pairs > 0], aes(lines, n_pairs, colour = level_name)) +
    geom_vline(xintercept = gate, linetype = "dashed", linewidth = 0.3, colour = "#555555") +
    geom_line(linewidth = 0.4) + geom_point(size = 1.2) +
    geom_text_repel(data = steps[lines == gate], aes(label = n_pairs), direction = "y", nudge_x = 0.2, size = 2.1, segment.size = 0.2, min.segment.length = 0, seed = 14, show.legend = FALSE) +
    scale_x_log10(breaks = unname(step_lines)) + scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000)) +
    scale_colour_manual(values = setNames(c("#222222", teal, purple), level_names), name = NULL) +
    panel_labs(4.1, title = sprintf("Of %s gene-by-lineage pairs with an altered line, %d have five altered lines per Y group", big(lineage$pairs), lineage$at_least_5),
      x = "Altered lines in each Y group, at least", y = "Gene-by-cancer-type pairs",
      caption = sprintf("%d male CRISPR lines. A pair is a gene and a cancer type with at least one altered line (driver mutation, amplification or deep deletion). Dashed: the gate of five lines; the numbers are the pairs that reach it. A test of whether the Y effect differs by driver status also needs five unaltered lines per Y group: %d pairs at lineage level, %d at disease and %d at subtype level.",
        tested$lines, lineage$all_four_cells, distance[level == "disease", all_four_cells], distance[level == "subtype", all_four_cells])) +
    theme_house() + theme(legend.position = "top")

  # Figure 14b: the testable splits, lines per Y group among altered lines.
  cells <- splits[, .(split_id, level, type, gene, neg_alt, pos_alt, neg_wt, pos_wt, label, kept_inside, low_information)]
  cells[, kind := fifelse(label == "specificity testable", "difference by driver status testable", "Y effect among altered lines only")]
  cells[, information := fifelse(low_information, "under 80% of the information kept", "80% or more kept")]
  save_source_data(cells, "fig14b_splits")
  named_splits <- cells[pmin(neg_alt, pos_alt) >= 10][, text := sprintf("%s, %s", gene, sub("Oral Cavity Squamous Cell Carcinoma", "oral cavity SCC", sub("Head and Neck Squamous Cell Carcinoma", "HNSCC", sub("Non-Small Cell Lung Cancer", "NSCLC", type))))]
  p14b <- ggplot(cells, aes(neg_alt, pos_alt, colour = kind, shape = information)) +
    geom_vline(xintercept = 10, linetype = "dotted", linewidth = 0.3, colour = "#888888") + geom_hline(yintercept = 10, linetype = "dotted", linewidth = 0.3, colour = "#888888") +
    geom_point(position = position_jitter(width = 0.025, height = 0.025, seed = 14), size = 1.2, alpha = 0.85) +
    geom_text_repel(data = named_splits, aes(label = text), size = 1.9, min.segment.length = 0, seed = 14, show.legend = FALSE, max.overlaps = 20) +
    scale_x_log10(breaks = c(5, 7, 10, 15, 20, 30, 40)) + scale_y_log10(breaks = c(5, 7, 10, 15, 20, 30)) +
    scale_colour_manual(values = c("difference by driver status testable" = red, "Y effect among altered lines only" = "#555555"), name = NULL) +
    scale_shape_manual(values = c("80% or more kept" = 16, "under 80% of the information kept" = 1), name = NULL) +
    guides(colour = guide_legend(nrow = 2, order = 1), shape = guide_legend(nrow = 2, order = 2)) +
    panel_labs(4.9, title = sprintf("%d of the %d testable splits have fewer than ten lines in the smaller altered cell", cells[, sum(pmin(neg_alt, pos_alt) < 10)], nrow(cells)),
      x = "Y-negative altered lines", y = "Y-positive altered lines",
      caption = sprintf("One point per split: a driver gene in a cancer type with five altered lines per Y group (%d splits in %d cancer types and %d genes; a split with the same lines at two levels is drawn once). Named: the %d splits with ten lines in both cells. Open: screening library takes more than a fifth of the information (%d splits). Dotted: ten lines.",
        nrow(cells), uniqueN(cells$type), uniqueN(cells$gene), nrow(named_splits), sum(cells$low_information))) +
    theme_house() + theme(legend.position = "top")

  # Figure 14c: the smallest shift each split finds 80% of the time, at four bars.
  limits <- read_ext("driver_type_limits")[unit == "sd" & target == 0.8]
  bar_names <- c(candidate = "named gene, among altered lines", named_difference = "named gene, difference by driver status",
    family = "search, the layer's rank-based test", family_raw = "search, least squares on Chronos scores")
  limits <- limits[bar %in% names(bar_names)]
  limits[, bar_label := factor(bar_names[bar], levels = bar_names)]
  limits[, smaller := pmin(neg_alt, pos_alt)]
  none_at <- 16
  limits[, drawn := fifelse(reachable, limit, none_at)]
  largest_shift <- max(read_ext("driver_type_shift_detection")[unit == "sd", shift])
  save_source_data(limits[, .(split_id, level, type, gene, neg_alt, pos_alt, bar, limit_sd = limit, reachable, detection_at_largest_shift = ceiling, normal_theory_sd = best_case)], "fig14c_detection_limits")
  at <- function(bar_name) limits[bar == bar_name]
  unreached <- at("family")[reachable == FALSE]
  p14c <- ggplot(limits, aes(smaller, drawn, colour = bar_label, shape = reachable)) +
    geom_point(position = position_jitter(width = 0.03, height = 0, seed = 14), size = 1, alpha = 0.85) +
    scale_x_log10(breaks = c(5, 7, 10, 15, 20, 26)) +
    scale_y_log10(breaks = c(0.5, 1, 2, 4, 10, none_at), labels = c("0.5", "1", "2", "4", "10", "not at 80%")) +
    scale_colour_manual(values = setNames(c(teal, purple, "#222222", "#999999"), bar_names), name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), guide = "none") +
    guides(colour = guide_legend(nrow = 2)) +
    panel_labs(4.4, title = sprintf("A gene named in advance needs %.1f SD in the median split; a search needs %.1f SD where the rank-based test reaches 80%%, which is %d of %d splits",
        median(at("candidate")$limit, na.rm = TRUE), median(at("family")$limit, na.rm = TRUE), sum(at("family")$reachable), nrow(at("family"))),
      x = "Lines in the smaller altered cell", y = "Shift found 80% of the time (SD)",
      caption = sprintf("Shift in the Y-negative altered lines of a split, in SD of the gene's spread between lines, simulated on %d real genes per shift, each tested as the one gene added to the split's fit. Named gene: p < 0.05. Search: the bar of the whole layer (p < %s). Crosses: not found 80%% of the time at any shift up to %g SD; the rank-based test then finds a %g SD shift in %.0f%% to %.0f%% of simulations, because ranks cap how far a few lines can move. Least squares on Chronos scores has no such cap and is not calibrated in cells this small; it is drawn for comparison. The difference by driver status can be tested in %d splits; in a search it reaches 80%% in none of them.",
        powered$simulation$shift_genes, fmt_p(tested$bars$family_inside), largest_shift, largest_shift, 100 * min(unreached$ceiling), 100 * max(unreached$ceiling), nrow(at("named_difference")))) +
    theme_house() + theme(legend.position = "top")

  # Figure 14d: the multiple of today's lines at which a shift of a given size is found.
  needed <- read_ext("driver_type_lines_needed")
  known <- unique(needed[class == "Y paralog", .(reference, convention, effect_sd)])
  needed <- needed[convention == "size" & bar %in% c("family", "specificity")]
  bucket_names <- c("today's lines", "up to 2 times", "up to 4 times", "up to 8 times", "up to 16 times", "more than 16 times")
  needed[, bucket := factor(fifelse(is.infinite(multiple), bucket_names[6], bucket_names[findInterval(multiple, c(1, 2, 4, 8) + 1e-9) + 1]), levels = bucket_names)]
  needed[, size_label := factor(sprintf("%g SD", effect_sd), levels = sprintf("%g SD", sort(unique(effect_sd))))]
  needed[, bar_label := factor(c(family = "Among altered lines", specificity = "Difference by driver status")[bar], levels = c("Among altered lines", "Difference by driver status"))]
  tally <- needed[, .(splits = .N), by = .(bar_label, size_label, bucket)]
  save_source_data(needed[, .(split_id, level, type, gene, bar, shift_sd = effect_sd, multiple, at_most, multiple_normal_theory)], "fig14d_lines_needed")
  whole <- needed[bar == "family"]
  of_size <- function(size) whole[abs(effect_sd - size) < 1e-9]
  range_of <- function(how) known[convention == how, sprintf("%.1f to %.1f", min(effect_sd), max(effect_sd))]
  copies_at_layer <- read_ext("driver_type_cut_check")[bar == "family" & is.finite(copied_against_real_log2) & at_grid_floor == FALSE, mean(copied_against_real_log2)]
  p14d <- ggplot(tally, aes(splits, size_label, fill = bucket)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
    facet_wrap(~bar_label, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = setNames(c(teal, "#DADADA", "#B5B5B5", "#8C8C8C", "#5E5E5E", red), bucket_names), name = NULL, drop = FALSE, guide = guide_legend(nrow = 2, byrow = TRUE)) +
    panel_labs(4.8, title = sprintf("In a search, a shift of 2 SD needs %.1f times today's lines in the median split and one of 1 SD %.1f times; %d of %d splits %s a 2 SD shift today",
        median(of_size(2)$multiple), median(of_size(1)$multiple), of_size(2)[, sum(multiple == 1)], nrow(of_size(2)), if (of_size(2)[, sum(multiple == 1)] == 1) "finds" else "find"),
      x = "Splits", y = "Size of the shift",
      caption = sprintf("For each split and each size: the multiple of today's lines at which the shift is found 80%% of the time at the bar of the whole layer, with the layer's rank-based test, from simulations on %s copies of the split's design (%d genes per shift; each copy carries the spread of a different real gene). Left: the Y effect among altered lines; right: the rule that also needs the difference from unaltered lines. The bar and the variance prior are today's. Checked on real splits cut down and copied back, the copies' limit at this bar is %.0f%% lower than the real one on average, so the multiples are if anything somewhat low. For scale, the four paired genes' effects are %s SD of their own spread and %s SD of a typical gene's.",
        listed(powered$simulation$multiples[-1]), powered$simulation$enlarged_genes, 100 * (1 - 2^copies_at_layer), range_of("own spread"), range_of("typical gene"))) +
    theme_house() + theme(legend.position = "bottom", strip.text = element_text(size = 7.5, face = "bold"))

  # Figure 14e: the four paired genes inside each split.
  paired <- read_ext("driver_type_paired_genes")
  paired[, smaller := pmin(n_neg_in, n_pos_in)]
  paired[, result := fifelse(shown, "shown (FDR < 0.05 among the four genes in all splits)", "not shown")]
  paired[, symbol := factor(symbol, levels = c("DDX3X", "EIF1AX", "ZFX", "RPS4X"))]
  save_source_data(paired[, .(split_id, symbol, n_neg_in, n_pos_in, estimate_chronos = inside, se = inside_se, p = inside_p, named_fdr, shown)], "fig14e_paired_genes_in_splits")
  share <- function(test) paired[test(smaller), round(100 * mean(shown))]
  p14e <- ggplot(paired, aes(smaller, inside, colour = result)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "#888888") +
    geom_point(position = position_jitter(width = 0.03, height = 0, seed = 14), size = 0.9, alpha = 0.85) +
    facet_wrap(~symbol, nrow = 1, scales = "free_y") +
    scale_x_log10(breaks = c(5, 10, 26)) +
    scale_colour_manual(values = c("shown (FDR < 0.05 among the four genes in all splits)" = teal, "not shown" = "#999999"), name = NULL) +
    panel_labs(5.2, title = sprintf("The four known Y-loss dependencies are seen in %d of %d tests inside splits: %d%% with five lines in the smaller cell, %d%% with ten or more",
        sum(paired$shown), nrow(paired), share(function(n) n == 5), share(function(n) n >= 10)),
      x = "Lines in the smaller altered cell", y = "Y effect among altered lines\n(Chronos units)",
      caption = "The four genes are named in advance, so they are judged as a family of their own (four genes in every split). A point is one gene in one split; the estimate is on the raw Chronos scale, the test on ranks.") +
    theme_house() + theme(legend.position = "top")

  # Figure 14f: all tests against their expected p-values.
  screen <- as.data.table(readRDS(file.path(ext_dir, "driver_type_screen.rds")))
  qq_of <- function(p, symbol, test) {
    out <- data.table(p = p, symbol = symbol)[is.finite(p)][order(p)]
    out[, expected := -log10((seq_len(.N) - 0.5) / .N)]
    out[, test := test]
    # Every test below 0.001 and every 200th of the rest: the curve is the same and the file small.
    out[p < 1e-3 | seq_len(.N) %% 200 == 0]
  }
  qq <- rbind(qq_of(screen$inside_p, screen$symbol, "Y effect among altered lines"), qq_of(screen$interaction_p, screen$symbol, "difference by driver status"))
  four <- c("DDX3X", "EIF1AX", "ZFX", "RPS4X")
  qq[, class := fifelse(test == "difference by driver status", "difference by driver status", fifelse(symbol %in% four, "Y effect among altered lines, a paired gene", "Y effect among altered lines"))]
  save_source_data(qq[, .(test, symbol, p, expected_minus_log10_p = expected)], "fig14f_tests")
  # Two splits at different levels can hold the same altered lines and give the same test; such a test is counted once.
  smallest <- unique(screen[inside_p < 1e-5, .(symbol, altered_set, p = signif(inside_p, 10))])
  null <- tested$null
  p14f <- ggplot(qq[order(class != "Y effect among altered lines")], aes(expected, -log10(p), colour = class)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#888888") +
    geom_hline(yintercept = -log10(tested$bars$family_inside), linetype = "dashed", linewidth = 0.3, colour = "#555555") +
    geom_point(aes(size = class), alpha = 0.8) +
    scale_colour_manual(values = c("Y effect among altered lines" = "#555555", "Y effect among altered lines, a paired gene" = teal, "difference by driver status" = purple), name = NULL) +
    scale_size_manual(values = c("Y effect among altered lines" = 0.3, "Y effect among altered lines, a paired gene" = 1.2, "difference by driver status" = 0.3), guide = "none") +
    guides(colour = guide_legend(nrow = 3, override.aes = list(size = 1.2))) +
    panel_labs(4.0, title = sprintf("No test passes; %d of the %d smallest distinct p-values belong to the four paired genes", smallest[symbol %in% four, .N], nrow(smallest)),
      x = "Expected -log10 p", y = "Observed -log10 p",
      caption = sprintf("%s tests of the Y effect among altered lines and %s of its difference from the unaltered lines. Dashed: the p-value one discovery among altered lines would need (%s). Smallest: p below 1e-5, a test that two splits share counted once. With Y status shuffled within driver status, %d of %d runs gave a hit among altered lines and %d a call.",
        big(tested$inside_tests), big(tested$interaction_tests), fmt_p(tested$bars$family_inside), null$runs_with_inside_only, null$runs, null$runs_with_call)) +
    theme_house() + theme(legend.position = "top")

  figure14 <- (p14a + p14b + plot_layout(widths = c(1, 1.3))) / (p14c + p14d + plot_layout(widths = c(1, 1.1))) / (p14e + p14f + plot_layout(widths = c(1.35, 1)))
  save_figure(tagged(figure14), "fig14_driver_type", width = full, height = 11.4)
}
# An argument draws one group of figures and leaves the others as they are; none draws all.
arguments <- commandArgs(trailingOnly = TRUE)
figure_groups <- c("extension", "nextgen", "driver_type")
if (!all(arguments %in% figure_groups)) stop("unknown argument: ", paste(setdiff(arguments, figure_groups), collapse = " "))
if (length(arguments) > 0) {
  if ("extension" %in% arguments) draw_extension_figures()
  if ("nextgen" %in% arguments) draw_nextgen_figure()
  if ("driver_type" %in% arguments) draw_driver_type_figure()
  log_msg("figures written (%s): %d files in the figure folder", paste(arguments, collapse = ", "), length(list.files(file.path(results_dir, "figures"))))
  quit(save = "no", status = 0)
}

# --- Figure 1: Y-status calls ---------------------------------------------------
rna <- calls[has_rna == TRUE & Sex %in% c("Male", "Female"), .(ModelID, Sex, rna_y_mean)]
save_source_data(rna, "fig1a_rna_score_by_sex")
female_max <- max(rna[Sex == "Female", rna_y_mean])
p1a <- ggplot(rna, aes(rna_y_mean, fill = Sex)) +
  geom_histogram(binwidth = 0.2, boundary = 0, colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = c(params$rna_call$depleted_mean_max, params$rna_call$expressed_mean_min), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~Sex, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(Female = unname(accent_colours["female"]), Male = "#555555"), guide = "none") +
  panel_labs(2.7, title = sprintf("No female line exceeds %.2f; males form two modes", female_max),
    x = "Mean log2(TPM+1) of seven Y genes", y = "Cell lines",
    caption = sprintf("%d female and %d male default RNA profiles. Dashed lines: depleted at 1.0 or below, expressed at 2.5 or above.", sum(rna$Sex == "Female"), sum(rna$Sex == "Male"))) +
  theme_house()
both <- calls[Sex == "Male" & has_rna == TRUE & !is.na(dna_y_ratio), .(ModelID, rna_y_mean, dna_y_ratio, dna_source, status = status_label(fifelse(is.na(y_call), NA_real_, as.numeric(y_call == "Y_negative"))))]
save_source_data(both, "fig1b_rna_vs_dna")
definite <- calls[Sex == "Male" & rna_call %in% c("depleted", "expressed") & dna_call %in% c("loss", "retained")]
agree <- sum((definite$rna_call == "depleted") == (definite$dna_call == "loss"))
p1b <- ggplot(both, aes(dna_y_ratio, rna_y_mean, colour = status)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_colour_manual(values = status_colours, name = NULL) +
  panel_labs(2.9, title = sprintf("RNA and DNA agree in %d of %d definite males", agree, nrow(definite)),
    x = "Y copy ratio (median of seven genes)", y = "Mean log2(TPM+1) of seven Y genes",
    caption = "Males with RNA and DNA. The copy ratio is relative to the sample mean: about 0.5 is one Y in a diploid line.") +
  theme_house() + theme(legend.position = c(0.78, 0.22))
tiers <- calls[Sex == "Male", .N, by = .(y_tier, status = status_label(fifelse(is.na(y_call), NA_real_, as.numeric(y_call == "Y_negative"))))]
tier_names <- c(concordant = "RNA and DNA agree", dna_only = "DNA only", rna_only = "RNA only", rna_led_dna_ambiguous = "RNA, DNA partial", rna_intermediate = "RNA intermediate", dna_ambiguous_no_rna = "DNA partial, no RNA", no_data = "no data", discordant = "RNA and DNA disagree")
tiers[, tier := factor(tier_names[y_tier], levels = rev(tier_names))]
save_source_data(tiers, "fig1c_call_tiers")
totals <- calls[Sex == "Male", .(neg = sum(y_call %in% "Y_negative"), pos = sum(y_call %in% "Y_positive"), aside = sum(is.na(y_call)))]
p1c <- ggplot(tiers, aes(N, tier, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(data = tiers[N >= 40], aes(label = N), position = position_stack(vjust = 0.5), size = 1.9, colour = "white") +
  scale_fill_manual(values = status_colours, name = NULL) +
  panel_labs(3.2, title = sprintf("%d Y-negative, %d Y-positive, %d set aside", totals$neg, totals$pos, totals$aside), x = "Male cell lines", y = NULL,
    caption = sprintf("All %s annotated males in Model.csv, by the evidence behind the call. Counts shown for bars of 40 or more.", format(sum(tiers$N), big.mark = ","))) +
  theme_house() + theme(legend.position = "none")
save_figure(tagged(p1a + p1b + p1c + plot_layout(widths = c(1, 1.1, 1.2))), "fig1_y_status_calls", width = full, height = 3.6)

# --- Figure 2: pan-cancer gene scan ----------------------------------------------
genes <- read_result("gene_pan_results")
genes[, class := fifelse(associated, "hit", fifelse(!is.na(y_partner), "X-Y partner, not a hit", "other"))]
save_source_data(genes[, .(symbol, chromosome, estimate, p, fdr, label, class)], "fig2a_gene_volcano")
hit_symbols <- genes[associated == TRUE, symbol]
p2a <- ggplot(genes[order(class != "other")], aes(estimate, -log10(p), colour = class)) +
  geom_point(aes(size = class), alpha = 0.8) +
  geom_text_repel(data = genes[associated == TRUE | (!is.na(y_partner) & p < 0.01)], aes(label = symbol), size = 2.3, min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = c(hit = red, "X-Y partner, not a hit" = unname(accent_colours["x_partner"]), other = grey), name = NULL) +
  scale_size_manual(values = c(hit = 1.6, "X-Y partner, not a hit" = 1.2, other = 0.3), guide = "none") +
  panel_labs(3.5, title = sprintf("%d genes pass FDR 0.05, all X-linked partners of Y genes", sum(genes$associated)),
    x = "Chronos effect, Y-negative minus Y-positive", y = "-log10 p",
    caption = sprintf("%s genes; limma, lineage and library adjusted, donor blocked. Negative: more dependent in Y-negative lines.", format(nrow(genes), big.mark = ","))) +
  theme_house() + theme(legend.position = c(0.7, 0.85))
effect <- readRDS(file.path(matrix_dir, "crispr_effect.rds"))
primary <- cohort[in_primary == TRUE & ModelID %in% rownames(effect)]
box <- rbindlist(lapply(genes[associated == TRUE, feature], function(f) data.table(symbol = gene_symbol(f), ModelID = primary$ModelID, status = status_label(primary$y_negative), chronos = effect[primary$ModelID, f])))
box[, symbol := factor(symbol, levels = hit_symbols)]
save_source_data(box, "fig2b_hit_boxplots")
strip <- genes[associated == TRUE, .(symbol = factor(symbol, levels = hit_symbols), text = sprintf("effect %.2f\nFDR %s", estimate, fmt_p(fdr)))]
p2b <- ggplot(box, aes(status, chronos, colour = status)) +
  geom_jitter(width = 0.2, size = 0.25, alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#222222", linewidth = 0.3, width = 0.5) +
  geom_text(data = strip, aes(x = 1.5, y = Inf, label = text), inherit.aes = FALSE, vjust = 1.1, size = 2.1) +
  facet_wrap(~symbol, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = status_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  panel_labs(5.3, title = "Each hit is a stronger dependency once the Y is gone", x = NULL, y = "Chronos gene effect",
    caption = sprintf("%d Y-negative and %d Y-positive lines. Box: median and quartiles. Effect and FDR from the primary model.", sum(primary$y_negative == 1), sum(primary$y_negative == 0))) +
  theme_house() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(tagged(p2a + p2b + plot_layout(widths = c(1, 1.5))), "fig2_pan_cancer_genes", width = full, height = 3.9)

# --- Figure 3: robustness of the gene hits ----------------------------------------
robust <- read_result("robustness_gene_hits")[symbol %in% hit_symbols]
definitions <- robust[check %in% c("Y-status definition", "one model per donor", "covariates")]
variant_names <- c(primary = "primary call", concordant_only = "RNA and DNA agree only", rna_only = "RNA only", dna_only = "DNA only", one_copy_controls = "one-copy controls only",
  no_relative_gain = "relative Y gain removed", lenient_rna_rule = "lenient RNA rule", xist_high_removed = "XIST-high lines removed", threshold_stable_only = "threshold-stable calls only")
definitions[, variant := fifelse(variant %in% names(variant_names), variant_names[variant], variant)]
definitions[, variant := factor(variant, levels = rev(unique(variant)))]
definitions[, symbol := factor(symbol, levels = hit_symbols)]
save_source_data(definitions, "fig3a_hit_robustness")
p3a <- ggplot(definitions, aes(estimate, variant, colour = check)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dotted") +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), size = 0.15, linewidth = 0.3) +
  facet_wrap(~symbol, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("Y-status definition" = "#222222", "one model per donor" = unname(accent_colours["par"]), covariates = unname(accent_colours["x_partner"])), name = NULL) +
  panel_labs(full, title = "The four effects hold under every Y-status definition and covariate set", x = "Chronos effect, Y-negative minus Y-positive (95% CI)", y = NULL,
    caption = "One robust (donor-clustered) fit per row. Definitions change who counts as Y-negative or Y-positive; covariates are added to lineage.") +
  theme_house() + theme(legend.position = "bottom")
by_lineage <- robust[check == "within lineage"]
by_lineage[, symbol := factor(symbol, levels = hit_symbols)]
save_source_data(by_lineage, "fig3b_hits_within_lineage")
p3b <- ggplot(by_lineage, aes(estimate, reorder(variant, n))) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dotted") +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se), size = 0.12, linewidth = 0.25, colour = red) +
  facet_wrap(~symbol, nrow = 1, scales = "free_x") +
  panel_labs(full, title = sprintf("The effects point the same way in %d of %d lineage-by-gene estimates", sum(by_lineage$estimate < 0), nrow(by_lineage)),
    x = "Chronos effect within lineage (95% CI)", y = NULL,
    caption = "Lineages with at least 10 lines and at least 3 lines per group, largest at the top.") +
  theme_house()
save_figure(tagged(p3a / p3b + plot_layout(heights = c(1.15, 1))), "fig3_gene_hit_robustness", width = full, height = 7.2)

# --- Figure 4: pan-cancer drug scan -----------------------------------------------
drugs <- read_result("drug_pan_results")
dataset_names <- c(prism_REPPRIMARY = "REP.PRIMARY", prism_REP1M = "REP.1M", prism_REP300 = "REP.300", prism_secondary = "PRISM secondary", gdsc1 = "GDSC1", gdsc2 = "GDSC2", ctrp = "CTRP v2")
drugs[, dataset_label := factor(dataset_names[dataset], levels = dataset_names)]
qq <- drugs[, .(observed = sort(-log10(p)), expected = sort(-log10(ppoints(.N)))), by = dataset_label]
save_source_data(qq, "fig4a_drug_qq")
p4a <- ggplot(qq, aes(expected, observed)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted") +
  geom_point(size = 0.3, colour = "#444444") +
  facet_wrap(~dataset_label, nrow = 2) +
  panel_labs(4.9, title = sprintf("No compound passes FDR 0.05 in any of %d dataset-screens", uniqueN(drugs$dataset)),
    x = "Expected -log10 p", y = "Observed -log10 p",
    caption = sprintf("%s compound tests; response on lineage and Y status, donor-clustered robust errors. Points on the dotted line are what chance gives. REP.* are the PRISM primary screens.", format(nrow(drugs), big.mark = ","))) +
  theme_house() + theme(strip.text = element_text(size = 7))
top_drugs <- drugs[order(p)][1:12]
top_drugs[, name := paste0(tolower(compound_name), " (", dataset_label, ")")]
top_drugs[, name := factor(name, levels = rev(name))]
save_source_data(top_drugs[, .(dataset, compound_name, mechanism, standardized, p, fdr, n_negative, n_positive)], "fig4b_top_compounds")
p4b <- ggplot(top_drugs, aes(standardized, name)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dotted") +
  geom_segment(aes(x = 0, xend = standardized, yend = name), linewidth = 0.3, colour = "#888888") +
  geom_point(aes(colour = standardized < 0), size = 1.5) +
  geom_text(aes(label = sprintf("FDR %.3f", fdr), hjust = fifelse(standardized < 0, 1.25, -0.25)), size = 1.9) +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = blue), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.45)) +
  panel_labs(4.1, title = sprintf("The smallest FDR is %.3f, and that compound does not replicate", min(drugs$fdr)), x = "Effect in response SDs", y = NULL,
    caption = "Twelve smallest p-values across all datasets; FDR is within dataset. Negative (red): Y-negative lines more sensitive. Positive (blue): less sensitive.") +
  theme_house() + theme(axis.text.y = element_text(size = 6))
save_figure(tagged(p4a + p4b + plot_layout(widths = c(1.2, 1))), "fig4_pan_cancer_drugs", width = full, height = 4)

# --- Figure 5: context screen -----------------------------------------------------
calibration <- rbind(read_result("context_gene_calibration")[, outcome := "Gene dependencies"], read_result("context_drug_calibration")[, outcome := "Drug responses"])
yield <- calibration[, .(interaction_hits = sum(interaction_fdr_class_hits), context_specific = sum(context_specific)), by = .(outcome, run)]
yield[, run_label := factor(fifelse(run == "observed", "observed", sub("shuffled_", "shuffle ", run)), levels = c("observed", paste("shuffle", 1:5)))]
yield_long <- melt(yield, id.vars = c("outcome", "run", "run_label"), variable.name = "measure", value.name = "count")
yield_long[, panel := paste0(outcome, "\n", fifelse(measure == "interaction_hits", "interactions at FDR 0.05", "context-specific calls"))]
yield_long[, panel := factor(panel, levels = c("Gene dependencies\ninteractions at FDR 0.05", "Drug responses\ninteractions at FDR 0.05", "Gene dependencies\ncontext-specific calls", "Drug responses\ncontext-specific calls"))]
save_source_data(yield_long, "fig5a_context_yield")
p5a <- ggplot(yield_long, aes(run_label, count, fill = run == "observed")) +
  geom_col(width = 0.7) +
  geom_text(aes(label = count), vjust = -0.3, size = 2) +
  facet_wrap(~panel, ncol = 2) +
  scale_fill_manual(values = c(`TRUE` = "#222222", `FALSE` = grey), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  panel_labs(4.4, title = "The context screens find no more than shuffled Y labels do", x = NULL, y = "Tests passing",
    caption = "The same procedure on the real labels and on five shuffles of Y status within lineage and context level. FDR is within context class.") +
  theme_house() + theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(size = 6.5))
tail <- calibration[, .(below = weighted.mean(below_0.001, tests), tests = sum(tests)), by = .(outcome, smallest_cell, shuffled = run != "observed")]
tail[, labels := fifelse(shuffled, "shuffled (mean of 5)", "observed")]
save_source_data(tail, "fig5b_context_tail")
p5b <- ggplot(tail[!is.na(smallest_cell)], aes(smallest_cell, below, colour = labels, group = labels)) +
  geom_hline(yintercept = 0.001, linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 0.3) + geom_point(size = 1) +
  facet_wrap(~outcome) +
  scale_colour_manual(values = c(observed = "#222222", "shuffled (mean of 5)" = unname(accent_colours["par"])), name = NULL) +
  scale_y_continuous(limits = c(0, 0.003)) +
  panel_labs(4.4, title = "Interaction p-values are calibrated even in five-line cells", x = "Lines in the smallest Y-by-context cell", y = "Fraction of tests with p < 0.001",
    caption = "Rank-normalised outcomes. Dotted line: the 0.001 expected under the null. Continuous contexts fall in the 50+ bin. On the raw scale the smallest cells gave 0.0026 for genes.") +
  theme_house() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))
context_summary <- read_result("context_gene_summary")
class_names <- c(lineage = "lineage", disease = "primary disease", subtype = "subtype", molecular_subtype = "molecular subtype", driver = "driver gene", amplification = "amplification", deletion = "deep deletion", genomic = "genome doubling, MSI", genomic_continuous = "CIN, ploidy, aneuploidy, LOH")
context_summary[, class_label := class_names[class]]
# Detectable effect after correcting for the number of tests in the context's class (BH with one hit is Bonferroni).
class_tests <- context_summary[, .(class_tests = sum(tests)), by = class]
context_summary <- merge(context_summary, class_tests, by = "class")
context_summary[, corrected_mde := median_interaction_mde / power_multiplier * (qnorm(1 - 0.05 / (2 * class_tests)) + qnorm(0.8))]
save_source_data(context_summary, "fig5c_context_mde")
p5c <- ggplot(context_summary, aes(median_interaction_mde, reorder(class_label, median_interaction_mde, median))) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6) +
  geom_jitter(height = 0.15, size = 0.4, alpha = 0.6) +
  panel_labs(full, title = sprintf("By the normal-theory formula, a typical context needed an interaction of %.2f Chronos units at a nominal 0.05, and about %.2f after correction for its class",
    median(context_summary$median_interaction_mde, na.rm = TRUE), median(context_summary$corrected_mde, na.rm = TRUE)),
    x = "Median detectable interaction (80% power, nominal 0.05)", y = NULL,
    caption = sprintf("One point per fitted context (%d). Continuous contexts are per standard deviation. For scale, the DDX3X effect is %.2f and the ZFX effect %.2f. These are best-case limits from raw-scale standard errors; Figure 9 gives what the rank-based test, as run, could find.", nrow(context_summary), genes[symbol == "DDX3X", estimate], genes[symbol == "ZFX", estimate])) +
  theme_house()
# The bottom panel has long axis labels; free() stops them from pushing the top row's panels to the right.
save_figure(tagged((p5a + p5b) / free(p5c) + plot_layout(heights = c(1.5, 1))), "fig5_context_screen", width = full, height = 7.2)

# --- Figure 6: expression differences ----------------------------------------------
de <- read_result("expression_de_results")
de[, dosage_class := fifelse(is_y_gene, "chrY", fifelse(par %in% TRUE, "PAR", fifelse(chromosome %in% "X", "chrX non-PAR", fifelse(!is.na(chromosome), "autosomal", NA_character_))))]
family <- de[placed == TRUE & is_y_gene == FALSE]
family[, class := factor(fifelse(par %in% TRUE, "pseudoautosomal", fifelse(de == TRUE, "FDR < 0.05", "other")), levels = c("other", "FDR < 0.05", "pseudoautosomal"))]
save_source_data(family[, .(symbol, chromosome, log_fc, p, fdr, class, label)], "fig6a_expression_volcano")
label_genes <- rbind(family[class == "FDR < 0.05"][order(p)][1:12], family[class == "pseudoautosomal"][order(p)][1:6])
chr_y_median <- median(de[dosage_class == "chrY" & mean_log_cpm > 1, log_fc])
p6a <- ggplot(family[order(class)], aes(log_fc, -log10(p), colour = class)) +
  geom_point(aes(size = class), alpha = 0.7) +
  geom_text_repel(data = label_genes, aes(label = symbol), size = 1.9, min.segment.length = 0, max.overlaps = 40, segment.size = 0.2, show.legend = FALSE) +
  scale_colour_manual(values = c(other = grey, "FDR < 0.05" = "#222222", pseudoautosomal = unname(accent_colours["par"])), name = NULL) +
  scale_size_manual(values = c(other = 0.2, "FDR < 0.05" = 0.6, pseudoautosomal = 0.9), guide = "none") +
  panel_labs(3.4, title = sprintf("%d genes differ outside chrY; pseudoautosomal genes lead", sum(family$de)),
    x = "log2 fold change, Y-negative over Y-positive", y = "-log10 p",
    caption = sprintf("limma-voom, lineage adjusted; %s genes with a known chromosome. chrY genes define the call and are left out (median log2 fold change %.1f).", format(nrow(family), big.mark = ","), chr_y_median)) +
  theme_house() + theme(legend.position = c(0.75, 0.85))
dosage <- de[dosage_class %in% c("autosomal", "chrX non-PAR", "PAR") & mean_log_cpm > 1]
dosage[, dosage_class := factor(dosage_class, levels = c("autosomal", "chrX non-PAR", "PAR"))]
save_source_data(dosage[, .(symbol, dosage_class, log_fc)], "fig6b_dosage_control")
medians <- dosage[, .(median = median(log_fc), n = .N), by = dosage_class]
beyond_axis <- sum(abs(dosage$log_fc) > 1.6)
p6b <- ggplot(dosage, aes(dosage_class, log_fc)) +
  geom_hline(yintercept = c(0, -1), linewidth = 0.3, linetype = c("solid", "dotted")) +
  geom_boxplot(outlier.size = 0.15, linewidth = 0.3, width = 0.6) +
  geom_text(data = medians, aes(y = 1.35, label = sprintf("n = %s\nmedian %.2f", format(n, big.mark = ","), median)), size = 1.9) +
  coord_cartesian(ylim = c(-1.6, 1.6)) +
  panel_labs(2.4, title = sprintf("Pseudoautosomal genes fall by %.2f, near the one-copy loss expected", abs(medians[dosage_class == "PAR", median])),
    x = NULL, y = "log2 fold change",
    caption = paste0("Expressed genes. Dotted line: -1, the loss of one of two copies.", if (beyond_axis > 0) sprintf(" %d genes lie beyond the axis.", beyond_axis) else "")) +
  theme_house() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
pathways <- read_result("expression_de_pathways")[collection == "hallmark" & fdr < 0.05][order(-abs(NES))][1:min(.N, 12)]
# Sentence case, with the acronyms put back.
pretty_set <- function(x) {
  x <- tolower(gsub("_", " ", sub("^HALLMARK_", "", x)))
  for (acronym in c("MYC", "DNA", "UV", "E2F", "G2M", "KRAS", "TNFA", "NFKB", "IL2", "IL6", "JAK", "STAT3", "STAT5", "TGF", "PI3K", "AKT", "WNT")) {
    x <- gsub(paste0("\\b", tolower(acronym), "\\b"), acronym, x)
  }
  gsub("\\bmtorc1\\b", "mTORC1", gsub("\\bmtor\\b", "mTOR", x))
}
pathways[, name := pretty_set(pathway)]
save_source_data(pathways[, .(pathway, NES, fdr, size)], "fig6c_hallmark")
hallmark_title <- if (all(pathways$NES > 0)) "The top Hallmark programmes are all higher in Y-negative lines" else "Hallmark programmes shifted in Y-negative lines"
p6c <- ggplot(pathways, aes(NES, reorder(name, NES), fill = NES > 0)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = red, `FALSE` = blue), guide = "none") +
  panel_labs(3.1, title = hallmark_title, x = "Normalised enrichment score", y = NULL,
    caption = "fgsea on moderated t, chrY removed; Hallmark sets at FDR < 0.05, largest first. Positive: higher in Y-negative lines.") +
  theme_house() + theme(axis.text.y = element_text(size = 6))
save_figure(tagged(p6a + p6b + p6c + plot_layout(widths = c(1.3, 0.9, 1.2))), "fig6_expression_differences", width = full, height = 4.1)

# --- Figure 7: expression associated with the dependencies --------------------------
expr <- readRDS(file.path(matrix_dir, "expr_tpm_pc.rds"))
lines <- cohort[ModelID %in% rownames(effect) & ModelID %in% rownames(expr)]
pairs <- data.table(x_gene = hit_symbols, y_gene = genes[match(hit_symbols, symbol), y_partner])
scatter <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  data.table(pair = sprintf("%s vs %s", pairs$x_gene[i], pairs$y_gene[i]), ModelID = lines$ModelID, status = status_label(lines$y_negative),
    expression = expr[lines$ModelID, match(pairs$y_gene[i], gene_symbol(colnames(expr)))], chronos = effect[lines$ModelID, match(pairs$x_gene[i], gene_symbol(colnames(effect)))])
}))
scatter[, pair := factor(pair, levels = unique(pair))]
save_source_data(scatter, "fig7a_paired_expression")
p7a <- ggplot(scatter[order(status != "set aside")], aes(expression, chronos, colour = status)) +
  geom_point(size = 0.5, alpha = 0.7) +
  facet_wrap(~pair, nrow = 1, scales = "free") +
  scale_colour_manual(values = status_colours, name = NULL) +
  panel_labs(full, title = "Dependency on each X gene tracks expression of its own Y partner", x = "Expression of the Y partner, log2(TPM+1)", y = "Chronos effect of the X gene",
    caption = sprintf("%d eligible males with CRISPR and RNA, including the %d set aside from the two-group call (mostly partial Y deletions).", nrow(lines), sum(is.na(lines$y_negative)))) +
  theme_house() + theme(legend.position = "bottom")
attribution <- read_result("expression_dependency_y_genes")[dependency %in% hit_symbols]
attribution[, dependency := factor(dependency, levels = rev(hit_symbols))]
attribution[, y_gene := factor(y_gene, levels = core)]
attribution[, strength := -log10(joint_p)]
save_source_data(attribution, "fig7b_y_gene_attribution")
strongest <- attribution[, .(top = as.character(y_gene[which.min(joint_p)]), paired = as.character(y_gene[paired == TRUE][1])), by = dependency]
p7b <- ggplot(attribution, aes(y_gene, dependency, fill = strength)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = attribution[paired == TRUE], fill = NA, colour = red, linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.1f", strength), colour = strength > 10), size = 2) +
  scale_fill_gradient(low = "white", high = "#222222", name = "-log10 p") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#222222"), guide = "none") +
  panel_labs(full, title = sprintf("With all seven Y genes in one model, the paired gene is the strongest term for %d of %d dependencies", sum(strongest$top == strongest$paired), nrow(strongest)),
    x = "Y gene expression (all seven in one model)", y = "Dependency",
    caption = "Robust fit of each dependency on all seven Y genes plus lineage and library, on the same lines as panel a. Outlined: the paired gene. Numbers: -log10 p of each term.") +
  theme_house()
save_figure(tagged(p7a / p7b + plot_layout(heights = c(1.5, 1))), "fig7_expression_and_dependency", width = full, height = 6.2)

# --- Figure 8: forests ---------------------------------------------------------------
accuracy <- read_result("forest_expression_accuracy")
accuracy[, model_short := factor(sub(" .*$", "", model), levels = c("m0", "m1", "m2", "m3"))]
named <- c(hit_symbols, "COMPOSITE")
accuracy[, group := fifelse(dependency %in% named, dependency, "other")]
save_source_data(accuracy, "fig8a_nested_accuracy")
# The title is chosen from the numbers: gain from adding all expression to the Y-gene model, per hit.
gain <- dcast(accuracy[dependency %in% hit_symbols], dependency ~ model_short, value.var = "r2_mean")[, .(dependency, gain = m3 - m2)]
helped <- gain[gain > 0.03, dependency]
accuracy_title <- if (length(helped) == 0) {
  "Y status and Y genes predict the hits; wider expression adds nothing"
} else {
  sprintf("Y status and Y genes predict the hits; wider expression adds only for %s", paste(helped, collapse = ", "))
}
line_colours <- c(setNames(c("#222222", red, blue, unname(accent_colours["x_partner"]), unname(accent_colours["par"]))[seq_along(named)], named), other = grey)
p8a <- ggplot(accuracy, aes(model_short, r2_mean, group = dependency, colour = group)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 0.35) + geom_point(size = 0.8) +
  geom_text_repel(data = accuracy[model_short == "m3" & group != "other"], aes(label = dependency), size = 2, hjust = 0, nudge_x = 0.3, direction = "y", segment.size = 0.2, show.legend = FALSE) +
  scale_colour_manual(values = line_colours, guide = "none") +
  scale_x_discrete(labels = c(m0 = "lineage\n+ library", m1 = "+ Y status", m2 = "+ 7 Y genes", m3 = "+ all\nexpression"), expand = expansion(mult = c(0.1, 0.45))) +
  panel_labs(2.9, title = accuracy_title, x = NULL, y = "Out-of-fold R2",
    caption = "Random forests, 10-fold cross-validation repeated 3 times, donors kept within folds. Grey: the other X partners and the context candidate. COMPOSITE: mean of the four hits.") +
  theme_house()
confounder <- read_result("forest_confounders")
confounder[, success := (success_seed1 + success_seed2) / 2]
save_source_data(confounder, "fig8b_confounder_check")
p8b <- ggplot(confounder[role != "hit"], aes(success)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = grey, colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = params$forest$confounder_survival_fraction, linewidth = 0.3, linetype = "dashed") +
  geom_point(data = confounder[role == "hit"], aes(success, 0), colour = red, size = 1.6) +
  annotate("text", x = 0.98, y = Inf, label = paste(confounder[role == "hit", symbol], collapse = "\n"), size = 2, hjust = 1, vjust = 1.2, colour = red) +
  panel_labs(2.9, title = sprintf("Y status beats confounder-matched shadows for %d of %d hits and %d of %d random genes",
    sum(confounder[role == "hit", survives]), sum(confounder$role == "hit"), sum(confounder[role != "hit", survives]), sum(confounder$role != "hit")),
    x = "Share of subsamples won by real Y status", y = "Random non-hit genes",
    caption = "Forest with Y status, lineage, library, genome doubling, TP53, aneuploidy, LOH, ploidy, MSI, CIN and screen quality. Dashed: the 0.8 bar. Red: the hits.") +
  theme_house()
forest_files <- list.files(file.path(derived_dir, "forest_y_status"), pattern = "^forest_y_status_summary_.*\\.rds$", full.names = TRUE)
forest_summary <- rbindlist(lapply(forest_files, readRDS))
forest_labels <- c(crispr_effect = "CRISPR dependencies", dataset_names)
forest_summary[, dataset_label := factor(forest_labels[dataset], levels = rev(forest_labels))]
forest_long <- melt(forest_summary[, .(dataset_label, `stable in both seeds` = as.numeric(stable), `confirmed with shuffled Y (mean)` = mean_null_confirmed)], id.vars = "dataset_label")
save_source_data(forest_summary, "fig8c_y_status_forest")
p8c <- ggplot(forest_long, aes(value, dataset_label, fill = variable)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = round(value, 1)), position = position_dodge(width = 0.7), hjust = -0.25, size = 1.9) +
  scale_fill_manual(values = c("stable in both seeds" = "#222222", "confirmed with shuffled Y (mean)" = grey), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  guides(fill = guide_legend(ncol = 1)) +
  panel_labs(2.9, title = sprintf("Forests that predict Y status confirm %d dependencies and %d compounds", forest_summary[dataset == "crispr_effect", stable], sum(forest_summary[dataset != "crispr_effect", stable])),
    x = "Features confirmed", y = NULL, caption = "Boruta-style selection, shadows shuffled within lineage, 100 iterations per run.") +
  theme_house() + theme(legend.position = "bottom")
save_figure(tagged(p8a + p8b + p8c + plot_layout(widths = c(1, 1, 1))), "fig8_forests", width = full, height = 4.2)


# --- Figures 9 and 10: what the cancer-type tests could detect -------------------------
power_types <- read_result("power_cancer_types")
power_shift <- read_result("power_shift_detection")
power_limits <- read_result("power_limits")
power_bench <- read_result("power_benchmarks")
power_detect <- read_result("power_benchmark_detection")
power_positive <- read_result("power_positive_controls")
power_params <- jsonlite::fromJSON(file.path(script_dir, "power_parameters.json"))
gate <- params$context$min_cell
level_names <- c(lineage = "Lineage", disease = "Primary disease", subtype = "Subtype")
dark <- "#1B3A57"
probability_fill <- scale_fill_gradient(low = "#F3F3F3", high = dark, limits = c(0, 1), breaks = c(0, 0.5, 1), name = "Probability of detection", guide = guide_colourbar(barwidth = 7, barheight = 0.5))
power_types[, level_name := factor(level_names[level], levels = level_names)]

# Rows of the heatmaps: the tested cancer types, largest first within each level.
tested_types <- unique(power_shift[, .(context_id, level, label, n_negative, n_positive)])
tested_types[, size := 1 / (1 / n_negative + 1 / n_positive)]
setorder(tested_types, -size)
# Non-breaking spaces keep the two counts together when a long name wraps.
tested_types[, row_label := vapply(sprintf("%s (%d\u00a0/\u00a0%d)", label, n_negative, n_positive), wrap_text, character(1), width = 44)]
tested_types[, row_label := factor(row_label, levels = rev(row_label))]
tested_types[, level_name := factor(level_names[level], levels = level_names)]
type_rows <- tested_types[, .(context_id, row_label, level_name)]

# Figure 9a: how many cancer types have enough lines.
census_groups <- c("enough lines", "both groups, too few lines", "one Y group missing")
census_long <- power_types[, .(types = .N), by = .(level_name, group = fifelse(smaller_group >= gate, census_groups[1], fifelse(smaller_group >= 1, census_groups[2], census_groups[3])))]
census_long[, group := factor(group, levels = census_groups)]
save_source_data(census_long, "fig9a_census")
enough <- census_long[group == census_groups[1]][order(level_name), types]
totals <- census_long[, .(types = sum(types)), keyby = level_name]$types
p9a <- ggplot(census_long, aes(types, level_name, fill = group)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.62) +
  geom_text(aes(label = types, colour = group), position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.2, show.legend = FALSE) +
  scale_fill_manual(values = setNames(c(dark, "#9FB3C8", "#E3E3E3"), census_groups), name = NULL) +
  scale_colour_manual(values = setNames(c("white", "#222222", "#222222"), census_groups)) +
  scale_y_discrete(limits = rev) +
  guides(fill = guide_legend(ncol = 1)) +
  panel_labs(2.9, title = sprintf("%d of %d lineages, %d of %d diseases and %d of %d subtypes have enough lines", enough[1], totals[1], enough[2], totals[2], enough[3], totals[3]),
    x = "Cancer types", y = NULL,
    caption = sprintf("Male cancer lines with CRISPR data and a Y call (n = %d). Enough: at least %d Y-negative and %d Y-positive lines, the screen's gate.", sum(power_types[level == "lineage", n_negative + n_positive]), gate, gate)) +
  theme_house() + theme(legend.position = "bottom", legend.justification = "left")

# Figure 9b: lines per cancer type, for every type with both groups.
both_groups <- power_types[smaller_group >= 1]
both_groups[, testable := smaller_group >= gate]
named_types <- both_groups[testable == FALSE][order(-(n_negative + n_positive)), head(.SD, 4), by = level_name]
named_types[, short_label := vapply(label, wrap_text, character(1), width = 22)]
save_source_data(both_groups[, .(level, label, n_negative, n_positive, status)], "fig9b_lines_per_cancer_type")
p9b <- ggplot(both_groups, aes(n_negative, n_positive)) +
  annotate("rect", xmin = gate, xmax = Inf, ymin = gate, ymax = Inf, fill = "#E8EEF4") +
  geom_count(aes(colour = testable), alpha = 0.85) +
  ggrepel::geom_text_repel(data = named_types, aes(label = short_label), size = 1.8, lineheight = 0.85, min.segment.length = 0, segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf, seed = 1) +
  facet_wrap(~level_name, nrow = 1) +
  scale_x_continuous(trans = "sqrt", breaks = c(1, 5, 10, 20, 40), limits = c(0.5, NA)) +
  scale_y_continuous(trans = "sqrt", breaks = c(1, 5, 10, 20, 40), limits = c(0.5, NA)) +
  scale_colour_manual(values = c(`TRUE` = dark, `FALSE` = "#9A9A9A"), guide = "none") +
  scale_size_area(max_size = 3.2, breaks = c(1, 3, 6), name = "Cancer types") +
  panel_labs(6, title = "Most cancer types have only a handful of lines in one Y group",
    x = "Y-negative lines", y = "Y-positive lines",
    caption = "One point per cancer type with both Y groups present; shaded area, at least five lines in each group. The largest untestable types are named. Square-root axes.") +
  theme_house() + theme(legend.position = "bottom")

# Figure 9c: probability of finding a clean shift, cancer type by size of the shift.
ceiling_of <- power_limits[unit == "sd" & bar == "one_type" & target == 0.8, .(context_id, ceiling)]
shift_heat <- power_shift[unit == "chronos", .(context_id, column = format(shift, nsmall = 1, drop0trailing = TRUE), shift, detection = one_type)]
shift_columns <- c(unique(shift_heat[order(shift), column]), "any size")
shift_heat <- rbind(shift_heat, ceiling_of[, .(context_id, column = "any size", shift = Inf, detection = ceiling)])
shift_heat <- merge(shift_heat, type_rows, by = "context_id")
shift_heat[, column := factor(trimws(column), levels = trimws(shift_columns))]
first_80 <- shift_heat[is.finite(shift) & detection >= 0.8][order(shift), .SD[1], by = context_id]
save_source_data(shift_heat[, .(level = level_name, cancer_type = row_label, shift_chronos = column, detection)], "fig9c_clean_shift_detection")
half_unit <- shift_heat[column == "0.5", sum(detection >= 0.8)]
never <- ceiling_of[, sum(ceiling < 0.8)]
p9c <- ggplot(shift_heat, aes(column, row_label, fill = detection)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(data = first_80, shape = 16, size = 0.9, colour = "white") +
  geom_text(data = shift_heat[column == "any size"], aes(label = sprintf("%.2f", detection), colour = detection > 0.55), size = 1.9, show.legend = FALSE) +
  facet_grid(level_name ~ ., scales = "free_y", space = "free_y") +
  probability_fill + scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#222222")) +
  panel_labs(full, title = sprintf("Analysed alone, %d of %d tested cancer types would find a shift of 0.5 Chronos units; in %d, no shift of any size is found 80%% of the time", half_unit, nrow(tested_types), never),
    x = "Shift in the Y-negative lines of one cancer type (Chronos units; common-essential genes sit at -1)", y = NULL,
    caption = sprintf("Cancer type (Y-negative / Y-positive lines). %d real selective genes per cell, shifted and tested with the screen's own model, each as the one gene added to the real ones. Found: p below the critical value of that cancer type's own gene family (FDR 0.05). White dot: smallest shift found 80%% of the time. 'any size': a shift of 10 SD.", power_params$shift_genes)) +
  theme_house() + theme(axis.text.y = element_text(size = 6, lineheight = 0.85), legend.position = "bottom", strip.text.y = element_text(size = 7), panel.grid = element_blank())
save_figure(tagged((p9a + p9b + plot_layout(widths = c(1, 2.1))) / free(p9c) + plot_layout(heights = c(1, 2.9))), "fig9_power_clean_shift", width = full, height = 10.2)

# Figure 10a: known vulnerabilities moved into each cancer type's Y sample sizes.
class_levels <- unique(power_bench$class)
bench_columns <- power_bench[, .(benchmark = label, class = factor(class, levels = class_levels), effect, column = sprintf("%s (%+.2f)", label, effect))]
bench_columns <- bench_columns[order(class, -abs(effect))]
bench_columns[, column := factor(column, levels = column)]
bench_heat <- merge(merge(power_detect[design == "own gene", .(context_id, benchmark, detection = one_type, as_run)], bench_columns, by = "benchmark"), type_rows, by = "context_id")
save_source_data(bench_heat[, .(level = level_name, cancer_type = row_label, class, benchmark, effect, detection, detection_as_run = as_run)], "fig10a_benchmark_detection")
reach <- bench_heat[, .(types = sum(detection >= 0.8)), by = benchmark]
example <- "BRAF in BRAF V600E"
p10a <- ggplot(bench_heat, aes(column, row_label, fill = detection)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(level_name ~ class, scales = "free", space = "free", labeller = labeller(class = label_wrap_gen(13))) +
  probability_fill +
  panel_labs(7.6, title = sprintf("Analysed alone, %d of %d tested cancer types would find an effect like %s; %d of %d benchmarks would be found in none", reach[benchmark == example, types], nrow(tested_types), example, reach[types == 0, .N], nrow(reach)),
    x = NULL, y = NULL,
    caption = sprintf("Each benchmark gene resampled %d times with its own lineage levels and noise, its dependent group moved to the Y-negative lines of the cancer type, and tested the same way. Found: as in Figure 9c. Effect in Chronos units in brackets, estimated in the same cohort. A copy is hard to find where the gene is already a dependency of the lineage (BRAF in skin, KRAS and CTNNB1 in bowel). Panel b: what the screen itself saw for the four established dependencies.", power_params$benchmark_replicates)) +
  theme_house() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6), axis.text.y = element_text(size = 6, lineheight = 0.85), strip.text = element_text(size = 6.5), strip.text.y = element_blank(), legend.position = "bottom", panel.grid = element_blank())

# Figure 10b: what the screen saw for the four established Y-loss dependencies.
result_levels <- c("not nominal", "nominal only", "genome-wide in this cancer type alone", "significant in the screen's family")
result_labels <- c("p of 0.05 or more", "p below 0.05 only", "genome-wide, this cancer type alone", "genome-wide, all types of its level")
positive_heat <- merge(power_positive, type_rows, by = "context_id")
# Columns in the order of the Y paralogs in panel a.
positive_heat[, symbol := factor(symbol, levels = sub(" .*$", "", bench_columns[class == "Y paralog", benchmark]))]
positive_heat[, result := factor(result, levels = result_levels)]
save_source_data(positive_heat[, .(level = level_name, cancer_type = row_label, symbol, inside, inside_p, own_fdr, inside_fdr_class, result, predicted_one_type)], "fig10b_positive_controls")
recovered <- positive_heat[, sum(result %in% result_levels[3:4])]
p10b <- ggplot(positive_heat, aes(symbol, row_label, fill = result)) +
  geom_tile(colour = "white", linewidth = 0.25, show.legend = TRUE) +
  facet_grid(level_name ~ "Observed", scales = "free", space = "free") +
  scale_fill_manual(values = setNames(c("#EEEEEE", "#C3CFDD", "#6F8FB5", dark), result_levels), labels = result_labels, name = "Observed", drop = FALSE, guide = guide_legend(nrow = 2, title.position = "top")) +
  panel_labs(1.5, title = sprintf("Seen in the real data: %d of %d", recovered, nrow(positive_heat)), x = NULL, y = NULL) +
  theme_house() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6), axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text = element_text(size = 6.5), strip.text.y = element_text(size = 7), legend.position = "bottom", panel.grid = element_blank())
save_figure(tagged(p10a + p10b + plot_layout(widths = c(nrow(bench_columns), 4.6), guides = "collect") & theme(legend.position = "bottom")), "fig10_power_benchmarks", width = full, height = 9.6)
draw_extension_figures()
draw_nextgen_figure()
draw_driver_type_figure()
log_msg("figures written: %d files", length(list.files(file.path(results_dir, "figures"))))
