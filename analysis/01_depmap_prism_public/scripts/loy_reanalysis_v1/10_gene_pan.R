# Stage 10 (plan Stage 3): pan-cancer gene dependencies associated with Y loss.
# Primary model: Chronos ~ lineage + library + Y_negative, donor-blocked limma,
# one BH family over all scanned genes. A hit is "LOY-associated"; it is
# "LOY-specific" only if it also survives adjustment for WGD, TP53 and
# aneuploidy, shows a matching contrast in dependency probability, and is not
# explained by a global shift in common-essential depletion.
# Negative estimates always mean "more dependent in Y-negative lines".

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
assert_freeze()
controls <- jsonlite::fromJSON(file.path(derived_dir, "pipeline_controls.json"))
stopifnot(controls$gene_controls_pass, controls$gene_calibration_pass)

gp <- params$gene_scan
cohort <- read_result("cohort")
gene_ann <- read_result("gene_annotation")
effect <- load_matrix("crispr_effect")
dependency <- load_matrix("crispr_dependency")
scan_genes <- gene_ann[in_scan == TRUE, feature]
base_covariates <- c("LibraryKY", "LibraryHumagneCD")
adjusted_covariates <- c(base_covariates, "WGD", "tp53_altered", "Aneuploidy")

# --- primary scan ------------------------------------------------------------
frame <- analysis_frame(cohort, effect)
design <- build_design(frame, tested = "exposed", covariates = base_covariates)
stopifnot(design$estimable, length(design$rows) == nrow(frame))
mat_t <- t(effect[frame$ModelID, scan_genes])
primary <- limma_scan(mat_t, design$design, "exposed", block = frame$PatientID)
donor_correlation <- attr(primary, "correlation")
primary[, fdr := bh(p)]
primary <- merge(primary, group_summary(mat_t, frame$exposed), by = "feature", sort = FALSE)
log_msg("primary scan: %d lines (%d Y-negative, %d Y-positive), %d genes, donor correlation %.3f, %d at FDR < %.2f",
  nrow(frame), sum(frame$exposed == 1), sum(frame$exposed == 0), nrow(primary), donor_correlation, sum(primary$fdr < gp$fdr), gp$fdr)

# --- adjusted for the genomic states that travel with Y loss -----------------
adjusted_design <- build_design(frame, tested = "exposed", covariates = adjusted_covariates)
adjusted_frame <- frame[adjusted_design$rows]
adjusted <- limma_scan(t(effect[adjusted_frame$ModelID, scan_genes]), adjusted_design$design, "exposed",
  block = adjusted_frame$PatientID, correlation = donor_correlation)
primary[, c("estimate_adjusted", "p_adjusted") := adjusted[, .(estimate, p)]]
log_msg("adjusted model: %d of %d lines have WGD, TP53 and aneuploidy", nrow(adjusted_frame), nrow(frame))

# --- global common-essential shift -------------------------------------------
essential <- gene_ann[in_scan == TRUE & common_essential %in% TRUE, feature]
frame[, essential_mean := rowMeans(effect[ModelID, essential], na.rm = TRUE)]
nuisance <- design$design[, colnames(design$design) != "exposed", drop = FALSE]
shift <- fit_robust(frame$essential_mean, frame$exposed, nuisance, cluster = frame$PatientID)
log_msg("common-essential mean, Y-negative minus Y-positive: %.4f (p = %.3g)", shift$estimate, shift$p)
essential_design <- build_design(frame, tested = "exposed", covariates = c(base_covariates, "essential_mean"))
essential_scan <- limma_scan(mat_t, essential_design$design, "exposed", block = frame$PatientID, correlation = donor_correlation)
primary[, c("estimate_essential_adjusted", "p_essential_adjusted") := essential_scan[, .(estimate, p)]]

# --- dependency probability contrast -----------------------------------------
dependent <- (dependency[frame$ModelID, scan_genes] >= gp$dependency_probability_cut) * 1
dep_scan <- scan_robust(dependent, frame$exposed, nuisance, cluster = frame$PatientID)
primary[, fraction_dependent_negative := colMeans(dependent[frame$exposed == 1, , drop = FALSE], na.rm = TRUE)]
primary[, fraction_dependent_positive := colMeans(dependent[frame$exposed == 0, , drop = FALSE], na.rm = TRUE)]
primary[, c("dependent_fraction_difference", "dependent_fraction_p") := dep_scan[, .(estimate, p)]]

# --- female lines as a second group without a Y -------------------------------
model <- fread(input_path("model"), select = c("ModelID", "Sex", "OncotreePrimaryDisease", "EngineeredModel"))
female_lines <- intersect(model[Sex == "Female" & OncotreePrimaryDisease != "Non-Cancerous" & EngineeredModel == "", ModelID], rownames(effect))
primary[, mean_female := colMeans(effect[female_lines, scan_genes], na.rm = TRUE)]

# --- labels ------------------------------------------------------------------
primary[, associated := fdr < gp$fdr & abs(estimate) >= gp$min_effect]
primary[, direction := fifelse(estimate < 0, "more dependent in Y-negative", "less dependent in Y-negative")]
primary[associated == TRUE, fdr_adjusted := bh(p_adjusted)]
primary[, survives_adjustment := associated & sign(estimate_adjusted) == sign(estimate) & fdr_adjusted < gp$fdr]
primary[, survives_essential_shift := associated & sign(estimate_essential_adjusted) == sign(estimate) & p_essential_adjusted < 0.05]
primary[, dependency_contrast := associated & sign(dependent_fraction_difference) == -sign(estimate) & dependent_fraction_p < 0.05]
primary[, label := fifelse(
  associated & survives_adjustment %in% TRUE & survives_essential_shift %in% TRUE & dependency_contrast %in% TRUE, "LOY-specific",
  fifelse(associated, "LOY-associated", "none")
)]
results <- merge(gene_ann[, .(feature, symbol, entrez, chromosome, par, y_partner, common_essential, strongly_selective, selective)], primary, by = "feature", sort = FALSE)
setorder(results, p)
write_result(results, "gene_pan_results")
hits <- results[associated == TRUE]
write_result(hits, "gene_pan_hits")
print(hits[, .(symbol, chromosome, y_partner, estimate, fdr, estimate_adjusted, fdr_adjusted, fraction_dependent_negative, fraction_dependent_positive, label)][1:min(.N, 40)])

# --- secondary scan: continuous Y-loss score on every male with RNA ----------
score_frame <- cohort[ModelID %in% rownames(effect) & !is.na(rna_y_score)]
score_frame[, y_loss_score := 1 - rna_y_score]
score_design <- build_design(score_frame, tested = "y_loss_score", covariates = base_covariates)
score_frame <- score_frame[score_design$rows]
score_scan <- limma_scan(t(effect[score_frame$ModelID, scan_genes]), score_design$design, "y_loss_score", block = score_frame$PatientID)
score_scan[, fdr := bh(p)]
score_scan <- merge(gene_ann[, .(feature, symbol, chromosome, y_partner)], score_scan, by = "feature", sort = FALSE)
setorder(score_scan, p)
write_result(score_scan, "gene_pan_continuous")
log_msg("continuous Y-loss score: %d lines, %d genes at FDR < %.2f", nrow(score_frame), sum(score_scan$fdr < gp$fdr), gp$fdr)

# --- RNAi replication --------------------------------------------------------
rnai <- load_matrix("rnai")
rnai_frame <- analysis_frame(cohort, rnai)
rnai_design <- build_design(rnai_frame, tested = "exposed")
rnai_frame <- rnai_frame[rnai_design$rows]
rnai_t <- t(rnai[rnai_frame$ModelID, ])
enough <- rowSums(!is.na(rnai_t[, rnai_frame$exposed == 1, drop = FALSE])) >= 20 & rowSums(!is.na(rnai_t[, rnai_frame$exposed == 0, drop = FALSE])) >= 20
rnai_scan <- limma_scan(rnai_t[enough, , drop = FALSE], rnai_design$design, "exposed", block = rnai_frame$PatientID)
rnai_scan[, fdr := bh(p)]
setnames(rnai_scan, "feature", "symbol")
setorder(rnai_scan, p)
write_result(rnai_scan, "gene_pan_rnai")
log_msg("RNAi: %d lines (%d Y-negative, %d Y-positive), %d genes, %d at FDR < %.2f",
  nrow(rnai_frame), sum(rnai_frame$exposed == 1), sum(rnai_frame$exposed == 0), nrow(rnai_scan), sum(rnai_scan$fdr < gp$fdr), gp$fdr)
replication <- merge(hits[, .(feature, symbol, estimate, fdr, label)], rnai_scan[, .(symbol, rnai_estimate = estimate, rnai_p = p, rnai_n = n)], by = "symbol", all.x = TRUE)
replication[, rnai_replicates := sign(rnai_estimate) == sign(estimate) & rnai_p < 0.05]
write_result(replication[order(fdr)], "gene_pan_hits_rnai")

# --- pathways ----------------------------------------------------------------
stat <- setNames(results$t, results$symbol)
pathways <- run_fgsea(stat)
write_result(pathways, "gene_pan_pathways")
print(pathways[fdr < 0.05, .(collection, pathway, NES, fdr)][1:min(.N, 20)])

benchmarks <- results[symbol %in% params$xy_paralogs$x_partner, .(symbol, y_partner, estimate, fdr, estimate_adjusted, p_adjusted, label)]
write_result(benchmarks, "gene_pan_xy_partners")
print(benchmarks)

lambda <- median(qchisq(results$p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
jsonlite::write_json(list(
  lines = nrow(frame), y_negative = sum(frame$exposed == 1), y_positive = sum(frame$exposed == 0),
  donors = uniqueN(frame$PatientID), genes_tested = nrow(results), donor_correlation = donor_correlation,
  fdr_hits = sum(results$fdr < gp$fdr), associated = sum(results$associated),
  loy_specific = sum(results$label == "LOY-specific"),
  more_dependent = sum(results$associated & results$estimate < 0), less_dependent = sum(results$associated & results$estimate > 0),
  adjusted_lines = nrow(adjusted_frame), lambda = lambda,
  essential_shift = list(estimate = shift$estimate, p = shift$p),
  continuous_lines = nrow(score_frame), continuous_fdr_hits = sum(score_scan$fdr < gp$fdr),
  rnai_lines = nrow(rnai_frame), hits_with_rnai = sum(!is.na(replication$rnai_p)), hits_rnai_replicated = sum(replication$rnai_replicates, na.rm = TRUE)
), file.path(derived_dir, "gene_pan_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("gene scan done: %d associated (%d LOY-specific), lambda %.2f", sum(results$associated), sum(results$label == "LOY-specific"), lambda)
