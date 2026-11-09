# Stage 26: Sanger organoids, the Y-loss screens.
# In order: known-biology controls (the stage stops if they fail), shuffled Y
# status, the scan over all organoids, the continuous scan, the screen by cancer
# type on rank-normalised scores, the genotype panel, and the follow-up of every
# hit. Scores are the study's scaled log fold changes; a negative estimate means
# "more needed in Y-negative organoids", as in the cell-line stages.
# Raw-scale fits carry one weight per organoid, estimated without Y status in the
# model, because the screens differ in noise. Rank-based fits carry none: with
# the weights, shuffled Y status gave an inflated test inside a cancer type
# (lambda up to 1.44 over four shuffles, against 0.97 to 1.14 without weights).
# LOY_EXTENSION_TRIAL=1 runs a few hundred genes (use with LOY_EXTENSION_OUTPUT).

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
assert_extension_freeze("organoid")

op <- ext_params$organoid
sp <- ext_params$shared
fdr_cut <- sp$fdr
trial <- Sys.getenv("LOY_EXTENSION_TRIAL") == "1"
cohort <- read_ext("organoid_cohort")
gene_ann <- read_ext("organoid_gene_annotation")
effect <- readRDS(file.path(ext_out_dir, "organoid_effect.rds"))
effect_v11 <- readRDS(file.path(ext_out_dir, "organoid_effect_v11.rds"))
depleted <- readRDS(file.path(ext_out_dir, "organoid_depleted.rds"))
prepared <- jsonlite::fromJSON(file.path(ext_out_dir, "organoid_prepare_summary.json"))
effect_scale <- prepared$effect_scale
base_covariates <- c("LibraryV1", "MSI")
adjusted_covariates <- c(base_covariates, "WGD", "tp53_altered", "PGA")
x_partners <- intersect(params$xy_paralogs$x_partner, gene_ann[in_scan == TRUE, symbol])
control_genes <- op$controls$gene
if (trial) {
  keep <- union(head(gene_ann[in_scan == TRUE, symbol], sp$trial_genes), c(x_partners, control_genes, head(gene_ann[selective == TRUE, symbol], sp$trial_genes)))
  gene_ann <- gene_ann[symbol %in% keep | !in_scan]
  log_msg("trial run on %d genes", sum(gene_ann$in_scan))
}
scan_genes <- gene_ann[in_scan == TRUE, symbol]
selective_genes <- gene_ann[selective == TRUE, symbol]
stopifnot(all(control_genes %in% colnames(effect)))

# --- one weight per organoid, estimated without Y status ------------------------------------------
# All malignant screened organoids of either sex, so every analysis below uses the same weights.
weight_frame <- cohort[malignant == TRUE]
weight_design <- cbind(lineage_design(weight_frame$OncotreeLineage), LibraryV1 = weight_frame$LibraryV1, MSI = weight_frame$MSI)
weight_design <- full_rank(weight_design)
array_weights <- arrayWeights(t(effect[weight_frame$ModelID, gene_ann[in_scan == TRUE, symbol]]), weight_design, method = "genebygene")
weights <- setNames(as.numeric(array_weights), weight_frame$ModelID)
weight_table <- weight_frame[, .(ModelID, OncotreeLineage, Sex, library, TPR, SD_LFC, weight = weights[ModelID])]
write_ext(weight_table, "organoid_weights")
log_msg("organoid weights: range %.2f to %.2f; rank correlation with the screen's true-positive rate %.2f", min(weights), max(weights), cor(weight_table$weight, weight_table$TPR, method = "spearman"))

weighted_ols <- function(y, x, z, ids, weighted = TRUE) {
  root <- if (weighted) sqrt(weights[ids]) else rep(1, length(ids))
  fit_ols(y * root, x * root, z * root)
}
lambda_of <- function(p) median(qchisq(p[is.finite(p)], 1, lower.tail = FALSE)) / qchisq(0.5, 1)

# --- 1. known-biology controls, read before any Y result -------------------------------------------
control_frame <- cohort[malignant == TRUE]
biomarker_of <- function(name) switch(name, KRAS = control_frame$alt_KRAS, TP53 = control_frame$alt_TP53, MSI = control_frame$MSI)
control_nuisance <- full_rank(cbind(lineage_design(control_frame$OncotreeLineage), LibraryV1 = control_frame$LibraryV1))
controls <- rbindlist(lapply(seq_along(control_genes), function(i) {
  biomarker <- biomarker_of(op$controls$biomarker[i])
  fit <- weighted_ols(effect[control_frame$ModelID, control_genes[i]], biomarker, control_nuisance, control_frame$ModelID)
  data.table(
    gene = control_genes[i], biomarker = op$controls$biomarker[i], altered = sum(biomarker == 1, na.rm = TRUE), reference = sum(biomarker == 0, na.rm = TRUE),
    estimate = fit$estimate, se = fit$se, p = fit$p, expected_sign = op$controls$expected_sign[i], p_max = op$controls$p_max[i]
  )
}))
controls[, pass := sign(estimate) == expected_sign & p < p_max]
write_ext(controls, "organoid_controls")
print(controls)
if (!all(controls$pass)) stop("known-biology controls failed: ", paste(controls[pass == FALSE, gene], collapse = ", "))

# --- the fits -------------------------------------------------------------------------------------------
primary <- cohort[in_primary == TRUE]
frame <- analysis_frame(primary, effect)
stopifnot(nrow(frame) == nrow(primary))
rank_rows <- function(mat_t) {
  out <- t(apply(mat_t, 1, rank_normalise))
  dimnames(out) <- dimnames(mat_t)
  out
}

pan_scan <- function(scan_frame, covariates = base_covariates, matrix = effect) {
  design <- build_design(scan_frame, tested = "exposed", covariates = covariates)
  stopifnot(design$estimable)
  scan_frame <- scan_frame[design$rows]
  out <- weighted_limma_scan(t(matrix[scan_frame$ModelID, scan_genes]), design$design, "exposed", weights = weights[scan_frame$ModelID])
  setnames(out, "feature", "symbol")
  out[, fdr := bh(p)]
  list(scan = out, frame = scan_frame, design = design$design)
}

# The screen by cancer type: one interaction per gene and the Y effect inside each type.
type_frame <- frame[OncotreeLineage %in% op$tumour_types]
type_genes <- union(selective_genes, union(gene_ann[cell_line_selective == TRUE & in_scan == TRUE, symbol], x_partners))
type_raw <- t(effect[type_frame$ModelID, type_genes])
type_rank <- rank_rows(type_raw)
type_screen <- function(screen_frame, with_raw = TRUE) {
  rbindlist(lapply(op$tumour_types, function(type) {
    ctx_frame <- with_context(screen_frame, as.numeric(screen_frame$OncotreeLineage == type), "binary")
    ranked <- fit_context_patterns(ctx_frame, type_rank, base_covariates, min_cell = sp$min_cell, min_residual_df = sp$min_residual_df, lineage = FALSE)
    out <- ranked[, .(symbol = feature, type = type, interaction_rank = interaction_estimate, interaction_p, interaction_unit, inside_rank = inside_estimate, inside_p, inside_unit, inside_df, n, n_neg_in, n_pos_in, n_neg_out, n_pos_out)]
    if (with_raw) {
      raw <- fit_context_patterns(ctx_frame, type_raw, base_covariates, min_cell = sp$min_cell, min_residual_df = sp$min_residual_df, weights = weights, lineage = FALSE)
      out[, c("interaction", "interaction_se", "inside", "inside_se") := raw[, .(interaction_estimate, interaction_se, inside_estimate, inside_se)]]
    }
    out
  }))
}
call_type_hits <- function(dt) {
  dt[, selective := symbol %in% selective_genes]
  first_type <- op$tumour_types[1]
  # The interaction is one test per gene (the second type's is its mirror image).
  by_gene <- dt[selective == TRUE & type == first_type, .(symbol, value = bh(interaction_p))]
  dt[, interaction_fdr := NA_real_]
  dt[selective == TRUE, interaction_fdr := by_gene$value[match(symbol, by_gene$symbol)]]
  dt[, inside_fdr := NA_real_]
  dt[selective == TRUE, inside_fdr := bh(inside_p)]
  dt[selective == TRUE, inside_fdr_own_type := bh(inside_p), by = type]
  dt[, type_specific := selective & interaction_fdr < fdr_cut & inside_fdr < fdr_cut]
  dt[, detected_inside_only := selective & inside_fdr < fdr_cut & !(type_specific %in% TRUE)]
  dt[]
}

# --- 2. the same procedures under shuffled Y status ------------------------------------------------------
shuffled_frame <- function(b) {
  set.seed(ext_params$seed + 100 * b)
  out <- copy(frame)
  out[, exposed := permute_within(exposed, OncotreeLineage)]
  out
}
calibration <- rbindlist(lapply(seq_len(sp$null_runs), function(b) {
  null_frame <- shuffled_frame(b)
  pan <- pan_scan(null_frame)$scan
  types <- call_type_hits(type_screen(null_frame[OncotreeLineage %in% op$tumour_types], with_raw = FALSE))
  data.table(
    run = paste0("shuffled_", b), pan_fdr_hits = sum(pan$fdr < fdr_cut), pan_below_0.001 = mean(pan$p < 0.001), pan_lambda = lambda_of(pan$p),
    type_specific = sum(types$type_specific, na.rm = TRUE), inside_fdr_hits = sum(types$inside_fdr < fdr_cut, na.rm = TRUE), inside_below_0.001 = mean(types[selective == TRUE, inside_p] < 0.001, na.rm = TRUE),
    inside_lambda = lambda_of(types[selective == TRUE, inside_p])
  )
}))

# --- 3. all organoids together ---------------------------------------------------------------------------
fit <- pan_scan(frame)
primary_scan <- fit$scan
frame <- fit$frame
nuisance <- fit$design[, colnames(fit$design) != "exposed", drop = FALSE]
mat_t <- t(effect[frame$ModelID, scan_genes])
primary_scan <- merge(primary_scan, group_summary(mat_t, frame$exposed)[, .(symbol = feature, mean_negative, mean_positive, n_negative, n_positive)], by = "symbol", sort = FALSE)
primary_scan[, estimate_chronos_equivalent := estimate / effect_scale]

# Adjusted for the genomic states that travel with Y loss.
adjusted <- pan_scan(frame, covariates = adjusted_covariates)
primary_scan[, c("estimate_adjusted", "p_adjusted") := adjusted$scan[match(primary_scan$symbol, symbol), .(estimate, p)]]
log_msg("adjusted model: %d of %d organoids have WGD, TP53 and fraction of genome altered", nrow(adjusted$frame), nrow(frame))

# A global shift in common-essential depletion would show as many small "hits".
essential <- intersect(gene_ann[common_essential == TRUE & in_scan == TRUE, symbol], colnames(effect))
frame[, essential_mean := rowMeans(effect[ModelID, essential])]
shift <- weighted_ols(frame$essential_mean, frame$exposed, nuisance, frame$ModelID)
quality_shift <- weighted_ols(frame[[op$quality_covariate]], frame$exposed, nuisance, frame$ModelID)
log_msg("common-essential mean, Y-negative minus Y-positive: %.4f (p = %.3g); screen true-positive rate: %.3f (p = %.3g)", shift$estimate, shift$p, quality_shift$estimate, quality_shift$p)
essential_scan <- pan_scan(frame, covariates = c(base_covariates, "essential_mean"))$scan
primary_scan[, c("estimate_essential_adjusted", "p_essential_adjusted") := essential_scan[match(primary_scan$symbol, symbol), .(estimate, p)]]

# Contrast in the study's own depletion calls, with screen quality in the model.
quality_nuisance <- cbind(nuisance, quality = frame[[op$quality_covariate]])
depletion_scan <- scan_robust(depleted[frame$ModelID, scan_genes], frame$exposed, quality_nuisance, cluster = frame$PatientID)
# Calls that the covariates explain in full (the same call in every organoid, or a call in the
# lone organoid of a tumour type) leave no Y contrast. Their estimate and p-value would be
# rounding noise, so the estimate is set to 0 and the p-value to missing.
unexplained <- qr.resid(qr(quality_nuisance), depleted[frame$ModelID, scan_genes] * 1)
no_contrast <- scan_genes[colSums(abs(unexplained) > 1e-8) == 0]
depletion_scan[feature %in% no_contrast, c("estimate", "p") := .(0, NA_real_)]
log_msg("depletion calls: %d of %d genes leave no Y contrast", length(no_contrast), length(scan_genes))
primary_scan[, fraction_depleted_negative := colMeans(depleted[frame[exposed == 1, ModelID], scan_genes])]
primary_scan[, fraction_depleted_positive := colMeans(depleted[frame[exposed == 0, ModelID], scan_genes])]
primary_scan[, c("depleted_fraction_difference", "depleted_fraction_p") := depletion_scan[match(primary_scan$symbol, feature), .(estimate, p)]]

primary_scan[, associated := fdr < fdr_cut & abs(estimate_chronos_equivalent) >= op$min_effect_chronos_equivalent]
primary_scan[, direction := fifelse(estimate < 0, "more needed in Y-negative", "less needed in Y-negative")]
primary_scan[associated == TRUE, fdr_adjusted := bh(p_adjusted)]
primary_scan[, survives_adjustment := associated & sign(estimate_adjusted) == sign(estimate) & fdr_adjusted < fdr_cut]
primary_scan[, survives_essential_shift := associated & sign(estimate_essential_adjusted) == sign(estimate) & p_essential_adjusted < 0.05]
primary_scan[, depletion_contrast := associated & sign(depleted_fraction_difference) == -sign(estimate) & depleted_fraction_p < 0.05]
primary_scan[, label := fifelse(
  associated & survives_adjustment %in% TRUE & survives_essential_shift %in% TRUE & depletion_contrast %in% TRUE, "LOY-specific",
  fifelse(associated, "LOY-associated", "none")
)]
primary_scan <- merge(gene_ann[, .(symbol, feature, chromosome, y_partner, common_essential, selective, cell_line_in_scan, seen_before)], primary_scan, by = "symbol", sort = FALSE)
setorder(primary_scan, p)
write_ext(primary_scan, "organoid_pan_results")
lambda <- lambda_of(primary_scan$p)
log_msg("scan over all organoids: %d organoids (%d Y-negative, %d Y-positive), %d genes, %d at FDR < %.2f, %d associated, %d LOY-specific, lambda %.2f",
  nrow(frame), sum(frame$exposed == 1), sum(frame$exposed == 0), nrow(primary_scan), sum(primary_scan$fdr < fdr_cut), fdr_cut, sum(primary_scan$associated), sum(primary_scan$label == "LOY-specific"), lambda)
print(primary_scan[1:15, .(symbol, chromosome, y_partner, estimate, estimate_chronos_equivalent, p, fdr, label, seen_before)])
partners <- primary_scan[symbol %in% params$xy_paralogs$x_partner, .(symbol, y_partner, estimate, estimate_chronos_equivalent, smd, unit, p, fdr, estimate_adjusted, p_adjusted, label)]
write_ext(partners, "organoid_xy_partners")
print(partners)

# --- 4. continuous scan on every eligible screened male ------------------------------------------------------
continuous_frame <- cohort[eligible == TRUE & !is.na(y_loss_fraction)][order(-qc_order)][!duplicated(PatientID)][order(ModelID)]
continuous_design <- build_design(continuous_frame, tested = "y_loss_fraction", covariates = base_covariates)
continuous_frame <- continuous_frame[continuous_design$rows]
continuous <- weighted_limma_scan(t(effect[continuous_frame$ModelID, scan_genes]), continuous_design$design, "y_loss_fraction", weights = weights[continuous_frame$ModelID])
setnames(continuous, "feature", "symbol")
continuous[, fdr := bh(p)]
setorder(continuous, p)
write_ext(continuous, "organoid_continuous")
log_msg("continuous Y-loss fraction: %d organoids, %d genes at FDR < %.2f", nrow(continuous_frame), sum(continuous$fdr < fdr_cut), fdr_cut)

# --- 5. by cancer type -------------------------------------------------------------------------------------------
types <- call_type_hits(type_screen(type_frame))
types[, focus := symbol %in% x_partners]
types[, inside_chronos_equivalent := inside / effect_scale]
types[, cell_line_selective := symbol %in% gene_ann[cell_line_selective == TRUE, symbol]]
types <- merge(types, gene_ann[, .(symbol, feature, chromosome, y_partner, seen_before)], by = "symbol", sort = FALSE)
write_ext(types[order(inside_p)], "organoid_type_screen")
print(types[selective == TRUE, .(tests = .N, min_interaction_p = min(interaction_p, na.rm = TRUE), min_inside_p = min(inside_p, na.rm = TRUE), inside_fdr_hits = sum(inside_fdr < fdr_cut, na.rm = TRUE), type_specific = sum(type_specific, na.rm = TRUE)), by = .(type, n_neg_in, n_pos_in)])

# --- 6. genotype panel (secondary) ---------------------------------------------------------------------------------
test_counts <- read_ext("organoid_test_counts")
panel <- c(setNames(paste0("alt_", op$genotype_genes), op$genotype_genes), WGD = "WGD", MSI = "MSI", `Metastasis-derived` = "metastasis")
panel <- panel[names(panel) %in% test_counts[class == "genotype" & supported == TRUE, context]]
panel_raw <- t(effect[frame$ModelID, selective_genes])
panel_rank <- rank_rows(panel_raw)
genotype <- rbindlist(lapply(names(panel), function(label) {
  ctx_frame <- with_context(frame, frame[[panel[[label]]]], "binary")
  # MSI is the context itself in the MSI row, so it cannot also be a covariate there.
  covariates <- setdiff(base_covariates, panel[[label]])
  ranked <- fit_context_patterns(ctx_frame, panel_rank, covariates, min_cell = sp$min_cell, min_residual_df = sp$min_residual_df)
  raw <- fit_context_patterns(ctx_frame, panel_raw, covariates, min_cell = sp$min_cell, min_residual_df = sp$min_residual_df, weights = weights)
  ranked[, .(symbol = feature, context = label, interaction_rank = interaction_estimate, interaction_p, inside_rank = inside_estimate, inside_p, interaction = raw$interaction_estimate, inside = raw$inside_estimate, n, n_neg_in, n_pos_in, n_neg_out, n_pos_out)]
}))
if (nrow(genotype) > 0) {
  genotype[, interaction_fdr := bh(interaction_p)]
  genotype[, inside_fdr := bh(inside_p)]
  genotype[, context_specific := interaction_fdr < fdr_cut & inside_fdr < fdr_cut]
  genotype[, detected_inside_only := inside_fdr < fdr_cut & !(context_specific %in% TRUE)]
  setorder(genotype, interaction_p)
}
write_ext(genotype, "organoid_genotype_screen")
log_msg("genotype panel: %d contexts (%s), %d tests, %d context-specific, %d inside-only", length(panel), paste(names(panel), collapse = ", "), nrow(genotype), sum(genotype$context_specific, na.rm = TRUE), sum(genotype$detected_inside_only, na.rm = TRUE))

# --- calibration: observed beside shuffled ----------------------------------------------------------------------------
calibration <- rbind(
  data.table(
    run = "observed", pan_fdr_hits = sum(primary_scan$fdr < fdr_cut), pan_below_0.001 = mean(primary_scan$p < 0.001), pan_lambda = lambda_of(primary_scan$p),
    type_specific = sum(types$type_specific, na.rm = TRUE), inside_fdr_hits = sum(types$inside_fdr < fdr_cut, na.rm = TRUE), inside_below_0.001 = mean(types[selective == TRUE, inside_p] < 0.001, na.rm = TRUE),
    inside_lambda = lambda_of(types[selective == TRUE, inside_p])
  ),
  calibration
)
write_ext(calibration, "organoid_calibration")
print(calibration)

# --- hits ----------------------------------------------------------------------------------------------------------------
hits <- rbind(
  primary_scan[associated == TRUE, .(kind = "all organoids", symbol, context = "All organoids", estimate, p, fdr, label, seen_before)],
  types[type_specific == TRUE, .(kind = "specific to the cancer type", symbol, context = type, estimate = inside, p = inside_p, fdr = inside_fdr, label = "type-specific", seen_before)],
  types[detected_inside_only == TRUE, .(kind = "inside the cancer type", symbol, context = type, estimate = inside, p = inside_p, fdr = inside_fdr, label = "inside only", seen_before)],
  if (nrow(genotype) > 0) genotype[context_specific == TRUE, .(kind = "specific to the genotype", symbol, context, estimate = inside, p = inside_p, fdr = inside_fdr, label = "context-specific", seen_before = symbol %in% op$seen_before)],
  fill = TRUE
)
hits[, class := fifelse(seen_before, "seen before", "discovery")]
write_ext(hits, "organoid_hits")
print(hits)

# --- 7. follow-up of every hit and of the X partners -----------------------------------------------------------------------
# Checks may only downgrade a hit. Alternative Y definitions are applied to these genes only.
definitions <- c(primary = "y_call", rna_only = "y_call_rna", dna_confirmed = "y_call_concordant", without_breadth_rule = "y_call_lenient", calls_on_disk = "y_call_on_disk", genome_male = "y_call_genome_male")
frame_for <- function(exposure = "y_call", all_of_patient = FALSE, drop_msi = FALSE, types_only = NULL) {
  out <- cohort[malignant == TRUE & !is.na(get(exposure)) & (male | exposure == "y_call_genome_male")]
  if (!all_of_patient) out <- out[order(-qc_order)][!duplicated(PatientID)]
  if (drop_msi) out <- out[MSI == 0]
  if (!is.null(types_only)) out <- out[OncotreeLineage %in% types_only]
  out <- out[order(ModelID)]
  out[, exposed := as.numeric(get(exposure) == "Y_negative")]
  out
}
one_fit <- function(values, fit_frame, covariates = base_covariates, lineage = TRUE) {
  design <- build_design(fit_frame, tested = "exposed", covariates = covariates, lineage = lineage)
  if (!design$estimable) return(list(estimate = NA_real_, se = NA_real_, p = NA_real_, n = 0L))
  fit_frame <- fit_frame[design$rows]
  weighted_ols(values[fit_frame$ModelID], fit_frame$exposed, design$design[, colnames(design$design) != "exposed", drop = FALSE], fit_frame$ModelID, weighted = lineage)
}
row_of <- function(fit, ...) data.table(..., estimate = fit$estimate, se = fit$se, p = fit$p, n = fit$n)
follow <- function(symbol, context) {
  inside_type <- if (context %in% op$tumour_types) context else NULL
  lineage <- is.null(inside_type)
  # Type-level tests are judged on ranks over the organoids of both types, as in the
  # screen; the scan over all organoids is judged on the raw scale.
  values_of <- function(matrix, ...) {
    if (lineage) return(matrix[, symbol])
    both_types <- frame_for(..., types_only = op$tumour_types)
    setNames(rank_normalise(matrix[both_types$ModelID, symbol]), both_types$ModelID)
  }
  base <- frame_for(types_only = inside_type)
  values <- values_of(effect)
  design <- build_design(base, tested = "exposed", covariates = base_covariates, lineage = lineage)
  base <- base[design$rows]
  checks <- followup_checks(values[base$ModelID], base$exposed, design$design[, colnames(design$design) != "exposed", drop = FALSE], base$OncotreeLineage,
    weights = if (lineage) weights[base$ModelID] else NULL, n_perm = sp$followup_permutations, seed = ext_params$seed + 26)
  by_definition <- rbindlist(lapply(names(definitions), function(name) {
    alt <- frame_for(definitions[[name]], types_only = inside_type)
    row_of(one_fit(values_of(effect, definitions[[name]]), alt, lineage = lineage), check = "Y-status definition", variant = name)
  }))
  no_msi <- frame_for(drop_msi = TRUE, types_only = inside_type)
  everyone <- frame_for(all_of_patient = TRUE, types_only = inside_type)
  other <- rbind(
    row_of(one_fit(values_of(effect, drop_msi = TRUE), no_msi, covariates = "LibraryV1", lineage = lineage), check = "MSI organoids excluded", variant = "MSS only"),
    row_of(one_fit(values_of(effect_v11), base, lineage = lineage), check = "library choice", variant = "v1_1 for dual-library organoids"),
    row_of(one_fit(values_of(effect, all_of_patient = TRUE), everyone, lineage = lineage), check = "all organoids of a patient", variant = "no one-per-patient rule"),
    row_of(one_fit(values, base, covariates = c(base_covariates, op$quality_covariate), lineage = lineage), check = "covariates", variant = "plus screen true-positive rate"),
    row_of(one_fit(values, base, covariates = adjusted_covariates, lineage = lineage), check = "covariates", variant = "plus WGD, TP53, genome altered")
  )
  summary_row <- data.table(
    check = c("permutation within cancer type", "leave one organoid out", "leave one organoid out"),
    variant = c(sprintf("%d shuffles", sp$followup_permutations), "smallest |effect|", "largest p"),
    estimate = c(checks$estimate, sign(checks$estimate) * checks$loo_min_abs, NA_real_), se = NA_real_, p = c(checks$permutation_p, NA_real_, checks$loo_max_p), n = checks$n
  )
  out <- rbind(by_definition, other, summary_row)
  sign_holds <- function(name) isTRUE(sign(by_definition[variant == name, estimate]) == sign(checks$estimate))
  out[, robust := checks$permutation_p < sp$robust_permutation_p & checks$loo_sign_stable & checks$loo_max_p < sp$robust_leave_one_out_p & sign_holds("rna_only") & sign_holds("dna_confirmed")]
  cbind(data.table(symbol = symbol, context = context, scale = if (lineage) "raw" else "rank"), out)
}
targets <- unique(rbind(
  hits[context %in% c("All organoids", op$tumour_types), .(symbol, context)],
  data.table(symbol = x_partners, context = "All organoids"),
  CJ(symbol = x_partners, context = op$tumour_types)
))
followup <- rbindlist(lapply(seq_len(nrow(targets)), function(i) follow(targets$symbol[i], targets$context[i])))
followup[, is_hit := paste(symbol, context) %in% hits[, paste(symbol, context)]]
write_ext(followup, "organoid_hit_followup")
print(followup[check %in% c("permutation within cancer type", "leave one organoid out") | variant %in% c("rna_only", "dna_confirmed"), .(symbol, context, check, variant, estimate = round(estimate, 3), p = signif(p, 3), robust)][1:min(.N, 40)])

jsonlite::write_json(list(
  organoids = nrow(frame), y_negative = sum(frame$exposed == 1), y_positive = sum(frame$exposed == 0), patients = uniqueN(frame$PatientID),
  genes_tested = nrow(primary_scan), fdr_hits = sum(primary_scan$fdr < fdr_cut), associated = sum(primary_scan$associated), loy_specific = sum(primary_scan$label == "LOY-specific"),
  lambda = lambda, essential_shift = list(estimate = shift$estimate, p = shift$p), quality_shift = list(estimate = quality_shift$estimate, p = quality_shift$p),
  adjusted_organoids = nrow(adjusted$frame), continuous_organoids = nrow(continuous_frame), continuous_fdr_hits = sum(continuous$fdr < fdr_cut),
  type_tests = types[selective == TRUE, .N], type_specific = sum(types$type_specific, na.rm = TRUE), inside_only = sum(types$detected_inside_only, na.rm = TRUE),
  genotype_contexts = names(panel), genotype_tests = nrow(genotype), genotype_specific = sum(genotype$context_specific, na.rm = TRUE),
  shuffled_pan_fdr_hits = calibration[run != "observed", pan_fdr_hits], shuffled_inside_fdr_hits = calibration[run != "observed", inside_fdr_hits],
  shuffled_pan_lambda = calibration[run != "observed", pan_lambda], shuffled_inside_lambda = calibration[run != "observed", inside_lambda],
  inside_lambda = calibration[run == "observed", inside_lambda],
  controls_pass = all(controls$pass), weights_range = range(weights), effect_scale = effect_scale, trial = trial,
  parameters_sha256 = ext_section_hash(c("shared", "organoid"))
), file.path(ext_out_dir, "organoid_scan_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("organoid screens done: %d associated over all organoids, %d type-specific, %d inside only, %d genotype-specific; shuffled runs gave %s scan hits",
  sum(primary_scan$associated), sum(types$type_specific, na.rm = TRUE), sum(types$detected_inside_only, na.rm = TRUE), sum(genotype$context_specific, na.rm = TRUE), paste(calibration[run != "observed", pan_fdr_hits], collapse = ", "))
