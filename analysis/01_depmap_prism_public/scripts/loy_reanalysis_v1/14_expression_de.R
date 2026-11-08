# Stage 14 (plan Stage 6, part 1): expression differences between Y-negative and
# Y-positive lines. limma-voom on raw counts, ~ lineage + Y_negative.
# chrY genes define the exposure, so they are a positive control only: they are
# reported, kept out of the BH family and removed from pathway ranking.
# Positive log fold change = higher in Y-negative lines.
# The weighted fit is run without donor blocking (it would need a per-gene GLS);
# a one-model-per-donor refit is the check on that.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
suppressPackageStartupMessages(library(edgeR))
assert_freeze()

fdr_cut <- params$gene_scan$fdr
cohort <- read_result("cohort")
counts <- read_default_profiles(input_path("expr_counts"))
frame <- analysis_frame(cohort, counts)
log_msg("expression: %d lines (%d Y-negative, %d Y-positive), %d genes in the count table", nrow(frame), sum(frame$exposed == 1), sum(frame$exposed == 0), ncol(counts))

gene <- fread(input_path("gene"), select = c("symbol", "entrez_id", "ensembl_gene_id", "location", "locus_group"))
gene[, chromosome := fifelse(grepl("^([0-9]{1,2}|X|Y)[pq]", location), sub("^([0-9]{1,2}|X|Y)[pq].*$", "\\1", location), NA_character_)]
gene[, par := grepl("^X", location) & grepl(" and Y", location)]
gene[, arm := fifelse(is.na(chromosome), NA_character_, sub("^([0-9]{1,2}|X|Y)([pq]).*$", "\\1\\2", location))]
gene[, entrez := as.character(as.integer(entrez_id))]
annotation <- data.table(feature = colnames(counts), symbol = gene_symbol(colnames(counts)), entrez = gene_entrez(colnames(counts)))
# Match on Entrez ID, then symbol, then Ensembl ID. A missing ID must not match:
# match(NA, x) returns the first NA in x, which would give every unnamed gene one wrong location.
match_known <- function(key, table) {
  out <- match(key, table)
  out[is.na(key) | key == ""] <- NA_integer_
  out
}
idx <- match_known(annotation$entrez, gene$entrez)
idx <- fifelse(is.na(idx), match_known(annotation$symbol, gene$symbol), idx)
idx <- fifelse(is.na(idx), match_known(annotation$symbol, gene$ensembl_gene_id), idx)
annotation[, c("chromosome", "arm", "par", "locus_group") := gene[idx, .(chromosome, arm, par, locus_group)]]
annotation[, is_y_gene := chromosome %in% "Y" & !(par %in% TRUE)]
# Genes with no HGNC record have no location in this collection, so they cannot be told
# apart from unannotated chrY transcripts. They are fitted and reported but kept out of the test family.
annotation[, placed := !is.na(chromosome)]
log_msg("annotation: %d of %d count-table genes have a location; %d chrY, %d PAR", sum(annotation$placed), nrow(annotation), sum(annotation$is_y_gene), sum(annotation$par %in% TRUE))
annotation[, y_partner := params$xy_paralogs$y_gene[match(symbol, params$xy_paralogs$x_partner)]]

fit_de <- function(sub_frame, covariates = character()) {
  design <- build_design(sub_frame, tested = "exposed", covariates = covariates)
  sub_frame <- sub_frame[design$rows]
  dge <- DGEList(t(counts[sub_frame$ModelID, , drop = FALSE]))
  keep <- filterByExpr(dge, design$design)
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  voomed <- voom(dge, design$design)
  fit <- eBayes(lmFit(voomed, design$design), robust = TRUE)
  tab <- topTable(fit, coef = "exposed", number = Inf, sort.by = "none", confint = TRUE)
  data.table(feature = rownames(tab), log_fc = tab$logFC, ci_low = tab$CI.L, ci_high = tab$CI.R, mean_log_cpm = tab$AveExpr, t = tab$t, p = tab$P.Value, lines = nrow(sub_frame))
}

primary <- fit_de(frame)
adjusted <- fit_de(frame, c("WGD", "tp53_altered", "Aneuploidy"))
one_per_donor <- fit_de(frame[order(ModelID)][!duplicated(PatientID)])
results <- merge(annotation, primary, by = "feature", sort = FALSE)
results[, c("log_fc_adjusted", "p_adjusted") := adjusted[match(results$feature, feature), .(log_fc, p)]]
results[, c("log_fc_one_per_donor", "p_one_per_donor") := one_per_donor[match(results$feature, feature), .(log_fc, p)]]
results[is_y_gene == FALSE & placed == TRUE, fdr := bh(p)]
results[, de := !is_y_gene & placed & !is.na(fdr) & fdr < fdr_cut]
results[de == TRUE, fdr_adjusted := bh(p_adjusted)]
results[, survives_adjustment := de & sign(log_fc_adjusted) == sign(log_fc) & fdr_adjusted < fdr_cut]
results[, label := fifelse(de & survives_adjustment %in% TRUE, "LOY-specific", fifelse(de, "LOY-associated", "none"))]
setorder(results, p)
write_result(results, "expression_de_results")
write_result(results[de == TRUE], "expression_de_hits")
log_msg("DE: %d genes fitted; test family %d (%d chrY and %d without a location set aside), %d at FDR < %.2f (%d up, %d down in Y-negative), %d survive adjustment",
  nrow(results), sum(results$placed & !results$is_y_gene), sum(results$is_y_gene), sum(!results$placed), sum(results$de), fdr_cut, sum(results$de & results$log_fc > 0), sum(results$de & results$log_fc < 0), sum(results$label == "LOY-specific"))
# Unplaced genes that behave like chrY genes, for the record.
unplaced <- results[placed == FALSE & p < 1e-10][order(p)]
write_result(unplaced[, .(feature, symbol, log_fc, mean_log_cpm, p)], "expression_de_unplaced_y_like")
log_msg("unplaced genes with p < 1e-10: %d (%d lower in Y-negative)", nrow(unplaced), sum(unplaced$log_fc < 0))

# Positive control: the genes that define the call must be strongly lower in Y-negative lines.
core <- results[symbol %in% params$y_genes_core]
print(core[, .(symbol, log_fc, p)])
stopifnot(nrow(core) == length(params$y_genes_core), all(core$log_fc < -2), all(core$p < 1e-20))

# Internal dosage control: PAR genes lose one of two copies with the Y, X non-PAR and autosomal genes lose none.
results[, dosage_class := fifelse(is_y_gene, "chrY", fifelse(par %in% TRUE, "PAR", fifelse(chromosome %in% "X", "chrX non-PAR", fifelse(!is.na(chromosome), "autosomal", NA_character_))))]
dosage <- results[!is.na(dosage_class) & mean_log_cpm > 1, .(genes = .N, median_log_fc = median(log_fc), fraction_de = mean(de)), by = dosage_class]
par_test <- wilcox.test(results[dosage_class == "PAR" & mean_log_cpm > 1, log_fc], results[dosage_class == "autosomal" & mean_log_cpm > 1, log_fc])
dosage[, par_vs_autosomal_p := par_test$p.value]
write_result(dosage, "expression_de_dosage_control")
print(dosage)
write_result(results[!is.na(y_partner), .(symbol, y_partner, log_fc, p, fdr, log_fc_adjusted, p_adjusted, label)], "expression_de_xy_partners")

# One-model-per-donor check on the unblocked fit.
shared <- results[is.finite(t) & !is_y_gene & placed == TRUE]
donor_check <- data.table(
  lines_all = primary$lines[1], lines_one_per_donor = one_per_donor$lines[1],
  log_fc_correlation = cor(shared$log_fc, shared$log_fc_one_per_donor, use = "complete.obs"),
  de_same_sign = mean(sign(shared[de == TRUE, log_fc]) == sign(shared[de == TRUE, log_fc_one_per_donor]), na.rm = TRUE),
  de_nominal_in_one_per_donor = mean(shared[de == TRUE, p_one_per_donor] < 0.05, na.rm = TRUE)
)
write_result(donor_check, "expression_de_donor_check")
print(donor_check)

# Pathways on the moderated t, chrY removed; one statistic per symbol.
ranked <- results[is_y_gene == FALSE & placed == TRUE & is.finite(t)][!duplicated(symbol)]
stat <- setNames(ranked$t, ranked$symbol)
pathways <- run_fgsea(stat, collections = c("hallmark", "reactome", "gobp"))
write_result(pathways, "expression_de_pathways")
print(pathways[fdr < 0.05, .(collection, pathway = substr(pathway, 1, 55), NES, fdr)][1:min(.N, 25)])
# Positional sets: enrichment of a whole cytoband points at a co-occurring copy-number change, not at Y loss.
positional <- run_fgsea(stat, collections = "positional", min_size = 10, max_size = 2000)
write_result(positional, "expression_de_positional")
print(positional[fdr < 0.05, .(pathway, size, NES, fdr)][1:min(.N, 15)])
arm_tab <- results[!is.na(arm) & is_y_gene == FALSE, .(genes = .N, de_up = sum(de & log_fc > 0), de_down = sum(de & log_fc < 0), median_t = median(t)), by = arm][genes >= 30][order(median_t)]
write_result(arm_tab, "expression_de_by_arm")

jsonlite::write_json(list(
  lines = nrow(frame), y_negative = sum(frame$exposed == 1), y_positive = sum(frame$exposed == 0),
  genes_fitted = nrow(results), genes_tested = sum(results$placed & !results$is_y_gene), genes_without_location = sum(!results$placed), de = sum(results$de), de_up = sum(results$de & results$log_fc > 0), de_down = sum(results$de & results$log_fc < 0),
  loy_specific = sum(results$label == "LOY-specific"), adjusted_lines = adjusted$lines[1],
  pathways_fdr = pathways[, .(significant = sum(fdr < 0.05)), by = collection],
  positional_fdr = sum(positional$fdr < 0.05)
), file.path(derived_dir, "expression_de_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("expression DE done")
