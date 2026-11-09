# Stage 32: DepMap NextGen models, the Y-loss screens.
# In order: known-biology controls (the stage stops if they fail), shuffled Y
# status, the four paired genes, the scan over all models, the continuous scan,
# the screen by cancer type on rank-normalised scores, and the follow-up of every
# hit. Scores are Chronos gene effects; a negative estimate means "more needed
# in Y-negative models", as in the cell-line stages.
# The model has one term per stratum (lineage and screen format), because Y status
# is tied to both here. Raw-scale fits carry one weight per model, estimated
# without Y status in the model; rank-based fits carry none, as in stage 26.
# Shuffles of Y status within stratum give the family-wise bar that a discovery
# outside the four paired genes must clear.
# The four paired genes also get a permutation p-value from a million shuffles, and
# the same rule is applied to it, because the weighted model's p-values run small here.
# LOY_EXTENSION_TRIAL=1 runs a few hundred genes and few shuffles on Y status
# shuffled within stratum, so a trial never shows a real result (use with
# LOY_EXTENSION_OUTPUT).

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
assert_extension_freeze("nextgen")

np <- ext_params$nextgen
sp <- ext_params$shared
fdr_cut <- sp$fdr
trial <- Sys.getenv("LOY_EXTENSION_TRIAL") == "1"
cohort <- read_ext("nextgen_cohort")
gene_ann <- read_ext("nextgen_gene_annotation")
effect <- readRDS(file.path(ext_out_dir, "nextgen_effect.rds"))
dependency <- readRDS(file.path(ext_out_dir, "nextgen_dependency.rds"))
test_counts <- read_ext("nextgen_test_counts")
chance <- read_ext("nextgen_replication_chance")
cell_pan <- read_result("gene_pan_results")
adjusted_covariates <- c("WGD", "tp53_altered", "Aneuploidy")
quality <- np$quality_covariate
null_runs <- np$null_runs
paired_permutations <- np$paired_permutations
feature_of <- function(symbols) gene_ann[in_scan == TRUE][match(symbols, symbol), feature]
named_features <- feature_of(np$named_genes)
control_features <- feature_of(np$controls$gene)
x_partner_features <- gene_ann[in_scan == TRUE & !is.na(y_partner), feature]
stopifnot(!anyNA(named_features), !anyNA(control_features))
if (trial) {
  keep <- union(head(gene_ann[in_scan == TRUE, feature], sp$trial_genes), c(named_features, control_features, x_partner_features, head(gene_ann[selective == TRUE, feature], sp$trial_genes)))
  gene_ann <- gene_ann[feature %in% keep | !in_scan]
  null_runs <- np$trial_null_runs
  paired_permutations <- sp$followup_permutations
  # Every Y column of the cohort models moves together, within stratum.
  y_columns <- c("y_call", "y_tier", "rna_call", "dna_call", "y_negative", "y_call_rna", "y_call_concordant", "y_call_lenient", "y_call_no_ceiling", "y_call_x_checked", "possible_female", "y_loss_fraction", "rna_y_mean", "dna_y_ratio", "x_ratio")
  set.seed(ext_params$seed)
  rows <- which(cohort$in_primary)
  cohort[rows, (y_columns) := cohort[permute_within(rows, cohort$stratum[rows]), y_columns, with = FALSE]]
  log_msg("trial run: %d genes, %d shuffles, Y status of the cohort shuffled within stratum", sum(gene_ann$in_scan), null_runs)
}
scan_features <- gene_ann[in_scan == TRUE, feature]
selective_features <- gene_ann[selective == TRUE, feature]
essential <- gene_ann[common_essential %in% TRUE & in_scan == TRUE, feature]
symbol_of <- function(features) gene_ann$symbol[match(features, gene_ann$feature)]
lambda_of <- function(p) median(qchisq(p[is.finite(p)], 1, lower.tail = FALSE)) / qchisq(0.5, 1)
# A shuffle that lands on the observed arrangement counts as a tie whatever the last digit says
# (inside a small cancer type the same arrangement comes up many times in 20,000 shuffles).
tie_tolerance <- 1e-9

# --- one weight per model, estimated without Y status ------------------------------------------------
# All malignant screened models of either sex, so every analysis below uses the same weights.
weight_frame <- cohort[malignant == TRUE]
weight_design <- full_rank(lineage_design(weight_frame$stratum))
array_weights <- arrayWeights(t(effect[weight_frame$ModelID, scan_features]), weight_design, method = "genebygene")
weights <- setNames(as.numeric(array_weights), weight_frame$ModelID)
stopifnot(all(is.finite(weights)), all(weights > 0))
weight_table <- weight_frame[, .(ModelID, name, stratum, Sex, ScreenNNMD, ScreenROCAUC, in_primary, y_call, weight = weights[ModelID])]
weight_table[, alone_in_stratum := .N == 1, by = stratum]
write_ext(weight_table, "nextgen_weights")
log_msg("weights: %d models, range %.2f to %.2f; rank correlation with NNMD %.2f; %d models alone in their stratum (weight %s)", nrow(weight_table), min(weights), max(weights),
  cor(weight_table$weight, weight_table$ScreenNNMD, method = "spearman"), sum(weight_table$alone_in_stratum), paste(round(weight_table[alone_in_stratum == TRUE, weight], 2), collapse = ", "))

weighted_ols <- function(y, x, z, ids, weighted = TRUE) {
  root <- if (weighted) sqrt(weights[ids]) else rep(1, length(ids))
  fit_ols(y * root, x * root, z * root)
}

# --- 1. known-biology controls, read before any Y result ---------------------------------------------
control_frame <- cohort[malignant == TRUE]
biomarker_of <- function(name) switch(name, KRAS = control_frame$kras_hotspot, TP53 = control_frame$tp53_altered)
control_nuisance <- full_rank(lineage_design(control_frame$stratum))
controls <- rbindlist(lapply(seq_along(control_features), function(i) {
  biomarker <- biomarker_of(np$controls$biomarker[i])
  fit <- weighted_ols(effect[control_frame$ModelID, control_features[i]], biomarker, control_nuisance, control_frame$ModelID)
  plain <- weighted_ols(effect[control_frame$ModelID, control_features[i]], biomarker, control_nuisance, control_frame$ModelID, weighted = FALSE)
  data.table(
    gene = np$controls$gene[i], biomarker = np$controls$biomarker[i], altered = sum(biomarker == 1, na.rm = TRUE), reference = sum(biomarker == 0, na.rm = TRUE),
    estimate = fit$estimate, se = fit$se, p = fit$p, estimate_unweighted = plain$estimate, p_unweighted = plain$p, expected_sign = np$controls$expected_sign[i], p_max = np$controls$p_max[i]
  )
}))
controls[, pass := sign(estimate) == expected_sign & p < p_max]
write_ext(controls, "nextgen_controls")
print(controls)
if (!all(controls$pass)) stop("known-biology controls failed: ", paste(controls[pass == FALSE, gene], collapse = ", "))

# --- the fits -------------------------------------------------------------------------------------------
# build_design takes its dummies from OncotreeLineage; here that column holds the stratum.
as_frame <- function(rows, exposure = "y_call") {
  out <- copy(rows)[order(ModelID)]
  out[, exposed := as.numeric(get(exposure) == "Y_negative")]
  out[, lineage := OncotreeLineage]
  out[, OncotreeLineage := stratum]
  out
}
frame <- as_frame(cohort[in_primary == TRUE])
stopifnot(!anyNA(frame$exposed), all(frame$ModelID %in% names(weights)))
rank_rows <- function(mat_t) {
  out <- t(apply(mat_t, 1, rank_normalise))
  dimnames(out) <- dimnames(mat_t)
  out
}
pan_scan <- function(scan_frame, covariates = character(), weighted = TRUE) {
  design <- build_design(scan_frame, tested = "exposed", covariates = covariates)
  stopifnot(design$estimable)
  scan_frame <- scan_frame[design$rows]
  out <- weighted_limma_scan(t(effect[scan_frame$ModelID, scan_features]), design$design, "exposed", weights = if (weighted) weights[scan_frame$ModelID] else NULL)
  out[, fdr := bh(p)]
  list(scan = out, frame = scan_frame, design = design$design)
}

# The screen by cancer type: the interaction and the Y effect inside each type with five models per cell.
types <- test_counts[supported == TRUE, context]
type_features <- union(selective_features, union(gene_ann[cell_line_selective == TRUE & in_scan == TRUE, feature], x_partner_features))
type_raw <- t(effect[frame$ModelID, type_features])
type_rank <- rank_rows(type_raw)
type_screen <- function(screen_frame, with_raw = TRUE) {
  rbindlist(lapply(types, function(type) {
    ctx_frame <- with_context(screen_frame, as.numeric(screen_frame$lineage == type), "binary")
    ranked <- fit_context_patterns(ctx_frame, type_rank, character(), min_cell = sp$min_cell, min_residual_df = sp$min_residual_df)
    out <- ranked[, .(feature, type = type, interaction_rank = interaction_estimate, interaction_p, interaction_unit, inside_rank = inside_estimate, inside_p, inside_unit, inside_df, n, n_neg_in, n_pos_in, n_neg_out, n_pos_out)]
    if (with_raw) {
      raw <- fit_context_patterns(ctx_frame, type_raw, character(), min_cell = sp$min_cell, min_residual_df = sp$min_residual_df, weights = weights)
      out[, c("interaction", "interaction_se", "inside", "inside_se") := raw[, .(interaction_estimate, interaction_se, inside_estimate, inside_se)]]
    }
    out
  }))
}
call_type_hits <- function(dt) {
  dt[, selective := feature %in% selective_features]
  dt[, c("interaction_fdr", "inside_fdr") := NA_real_]
  dt[selective == TRUE, interaction_fdr := bh(interaction_p)]
  dt[selective == TRUE, inside_fdr := bh(inside_p)]
  dt[, type_specific := selective & interaction_fdr < fdr_cut & inside_fdr < fdr_cut]
  dt[, detected_inside_only := selective & inside_fdr < fdr_cut & !(type_specific %in% TRUE)]
  dt[]
}

# --- 2. the same procedures under Y status shuffled within stratum -----------------------------------------
reference_sign <- sign(cell_pan[match(np$named_genes, symbol), estimate])
shuffled_frame <- function(b) {
  set.seed(ext_params$seed + 100 * b)
  out <- copy(frame)
  out[, exposed := permute_within(exposed, OncotreeLineage)]
  out
}
null_run <- function(b) {
  null_frame <- shuffled_frame(b)
  pan <- pan_scan(null_frame)$scan
  named <- pan[match(named_features, feature)]
  needed <- pan[feature %in% essential]
  out <- data.table(
    run = sprintf("shuffled_%03d", b), pan_min_p = min(pan$p), pan_fdr_hits = sum(pan$fdr < fdr_cut), pan_associated = sum(pan$fdr < fdr_cut & abs(pan$estimate) >= np$min_effect_chronos),
    pan_below_0.001 = mean(pan$p < 0.001), pan_below_1e5 = sum(pan$p < 1e-5), pan_lambda = lambda_of(pan$p),
    essential_lambda = lambda_of(needed$p), essential_mean_t = mean(needed$t), essential_positive = mean(needed$estimate > 0),
    named_min_p = min(named$p), named_replicated = sum(bh(named$p) < fdr_cut & sign(named$estimate) == reference_sign)
  )
  if (length(types) > 0) {
    by_type <- call_type_hits(type_screen(null_frame, with_raw = FALSE))[selective == TRUE]
    out <- cbind(out, data.table(
      inside_min_p = min(by_type$inside_p, na.rm = TRUE), interaction_min_p = min(by_type$interaction_p, na.rm = TRUE), inside_fdr_hits = sum(by_type$inside_fdr < fdr_cut, na.rm = TRUE),
      type_specific = sum(by_type$type_specific, na.rm = TRUE), inside_below_0.001 = mean(by_type$inside_p < 0.001, na.rm = TRUE), inside_lambda = lambda_of(by_type$inside_p)
    ))
  }
  out
}
null_parts <- parallel::mclapply(seq_len(null_runs), null_run, mc.cores = n_cores)
failed <- !vapply(null_parts, is.data.table, logical(1))
if (any(failed)) stop("shuffled runs failed: ", paste(utils::head(unlist(lapply(null_parts[failed], as.character)), 3), collapse = " | "))
null_table <- rbindlist(null_parts)
pan_bar <- family_wise_bar(null_table$pan_min_p, np$family_wise_level)
inside_bar <- if (length(types) > 0) family_wise_bar(null_table$inside_min_p, np$family_wise_level) else NA_real_
log_msg("%d shuffles: %d with a gene at FDR < %.2f in the scan (%.1f%%), median lambda %.2f; family-wise bar p < %.3g for the scan%s", null_runs, sum(null_table$pan_fdr_hits > 0), fdr_cut,
  100 * mean(null_table$pan_fdr_hits > 0), median(null_table$pan_lambda), pan_bar, if (length(types) > 0) sprintf(" and p < %.3g inside a cancer type", inside_bar) else "")

# --- 3. all models together -------------------------------------------------------------------------------
fit <- pan_scan(frame)
primary_scan <- fit$scan
stopifnot(nrow(fit$frame) == nrow(frame))
nuisance <- fit$design[, colnames(fit$design) != "exposed", drop = FALSE]
mat_t <- t(effect[frame$ModelID, scan_features])
primary_scan <- merge(primary_scan, group_summary(mat_t, frame$exposed), by = "feature", sort = FALSE)
primary_scan[, se := estimate / t]
unweighted <- pan_scan(frame, weighted = FALSE)$scan
primary_scan[, c("estimate_unweighted", "p_unweighted") := unweighted[match(primary_scan$feature, feature), .(estimate, p)]]

# Adjusted for the genomic states that travel with Y loss.
adjusted <- pan_scan(frame, covariates = adjusted_covariates)
primary_scan[, c("estimate_adjusted", "p_adjusted") := adjusted$scan[match(primary_scan$feature, feature), .(estimate, p)]]
log_msg("adjusted model: %d of %d models have WGD, TP53 and aneuploidy", nrow(adjusted$frame), nrow(frame))

# Screen quality differs between the Y groups here, so it gets its own adjusted scan.
quality_scan <- pan_scan(frame, covariates = quality)$scan
primary_scan[, c("estimate_quality_adjusted", "p_quality_adjusted") := quality_scan[match(primary_scan$feature, feature), .(estimate, p)]]

# A global shift in common-essential depletion would show as many small "hits".
frame[, essential_mean := rowMeans(effect[ModelID, essential])]
# These two per-model summaries are tested without weights and against shuffles within stratum. The
# weights follow screen quality, and with them 9.7% of 3,000 shuffles gave p < 0.05 for the
# common-essential mean (5.1% without; checked on shuffled labels before the first real run).
shift_test <- function(values) followup_checks(values, frame$exposed, nuisance, frame$OncotreeLineage, n_perm = sp$followup_permutations, seed = ext_params$seed + 31, tie_tolerance = tie_tolerance)
shift <- shift_test(frame$essential_mean)
quality_shift <- shift_test(frame[[quality]])
log_msg("common-essential mean, Y-negative minus Y-positive: %.4f (p = %.3g, shuffles %.3g); screen NNMD: %.3f (p = %.3g, shuffles %.3g)", shift$estimate, shift$p, shift$permutation_p, quality_shift$estimate, quality_shift$p, quality_shift$permutation_p)
essential_scan <- pan_scan(frame, covariates = "essential_mean")$scan
primary_scan[, c("estimate_essential_adjusted", "p_essential_adjusted") := essential_scan[match(primary_scan$feature, feature), .(estimate, p)]]

# Contrast in the dependency calls (probability of at least the frozen cut), with screen quality in the model.
dependent <- (dependency[frame$ModelID, scan_features] >= params$gene_scan$dependency_probability_cut) * 1
stopifnot(!anyNA(dependent))
quality_nuisance <- cbind(nuisance, quality = frame[[quality]])
dependency_scan <- scan_robust(dependent, frame$exposed, quality_nuisance, cluster = frame$PatientID)
# Calls that the covariates explain in full leave no Y contrast. Their estimate and p-value would be
# rounding noise, so the estimate is set to 0 and the p-value to missing (as in stage 26).
unexplained <- qr.resid(qr(quality_nuisance), dependent)
no_contrast <- scan_features[colSums(abs(unexplained) > 1e-8) == 0]
dependency_scan[feature %in% no_contrast, c("estimate", "p") := .(0, NA_real_)]
log_msg("dependency calls: %d of %d genes leave no Y contrast", length(no_contrast), length(scan_features))
primary_scan[, fraction_dependent_negative := colMeans(dependent[frame$exposed == 1, , drop = FALSE])]
primary_scan[, fraction_dependent_positive := colMeans(dependent[frame$exposed == 0, , drop = FALSE])]
primary_scan[, c("dependent_fraction_difference", "dependent_fraction_p") := dependency_scan[match(primary_scan$feature, feature), .(estimate, p)]]
# The same contrast without screen quality, which is how the cell-line stage has it; reported beside the other, never used for a label.
dependency_plain <- scan_robust(dependent, frame$exposed, nuisance, cluster = frame$PatientID)
# Without screen quality in the model, calls that vary only inside strata holding one Y group say nothing
# about Y: the estimate is zero up to rounding and its p-value would be noise that differs between processors.
exposure_left <- qr.resid(qr(nuisance), frame$exposed)
no_contrast_plain <- scan_features[colSums(abs(qr.resid(qr(nuisance), dependent) * exposure_left) > 1e-8) == 0]
dependency_plain[feature %in% no_contrast_plain, p := NA_real_]
primary_scan[, dependent_fraction_p_without_quality := dependency_plain$p[match(primary_scan$feature, dependency_plain$feature)]]

primary_scan[, named := feature %in% named_features]
primary_scan[, associated := fdr < fdr_cut & abs(estimate) >= np$min_effect_chronos]
primary_scan[, passes_family_wise_bar := p < pan_bar]
primary_scan[, discovery := associated & passes_family_wise_bar & !named]
primary_scan[, direction := fifelse(estimate < 0, "more needed in Y-negative", "less needed in Y-negative")]
primary_scan[associated == TRUE, fdr_adjusted := bh(p_adjusted)]
primary_scan[, survives_adjustment := associated & sign(estimate_adjusted) == sign(estimate) & fdr_adjusted < fdr_cut]
primary_scan[, survives_quality := associated & sign(estimate_quality_adjusted) == sign(estimate) & p_quality_adjusted < 0.05]
primary_scan[, survives_essential_shift := associated & sign(estimate_essential_adjusted) == sign(estimate) & p_essential_adjusted < 0.05]
primary_scan[, dependency_contrast := associated & sign(dependent_fraction_difference) == -sign(estimate) & dependent_fraction_p < 0.05]
primary_scan[, label := fifelse(
  associated & survives_adjustment %in% TRUE & survives_quality %in% TRUE & survives_essential_shift %in% TRUE & dependency_contrast %in% TRUE, "LOY-specific",
  fifelse(associated, "LOY-associated", "none")
)]
primary_scan <- merge(gene_ann[, .(feature, symbol, chromosome, y_partner, common_essential, selective, cell_line_in_scan, organoid_in_scan)], primary_scan, by = "feature", sort = FALSE)
setorder(primary_scan, p)
write_ext(primary_scan, "nextgen_pan_results")
lambda <- lambda_of(primary_scan$p)
log_msg("scan over all models: %d models (%d Y-negative, %d Y-positive), %d genes, %d at FDR < %.2f, %d associated, %d discoveries, %d LOY-specific, lambda %.2f",
  nrow(frame), sum(frame$exposed == 1), sum(frame$exposed == 0), nrow(primary_scan), sum(primary_scan$fdr < fdr_cut), fdr_cut, sum(primary_scan$associated), sum(primary_scan$discovery), sum(primary_scan$label == "LOY-specific"), lambda)
print(primary_scan[1:15, .(symbol, chromosome, y_partner, estimate, p, fdr, passes_family_wise_bar, label)])

# --- 4. the four paired genes, a family of their own ---------------------------------------------------------
# Replicated: the cell-line sign and BH below the cut among the four. Not replicated: not that, and
# significantly weaker than in the cell lines (both are Chronos scores). Otherwise this cohort cannot tell.
paired <- primary_scan[match(named_features, feature), .(feature, symbol, y_partner, estimate, se, ci_low, ci_high, p, df, n_negative, n_positive, mean_negative, mean_positive, estimate_unweighted, p_unweighted, estimate_adjusted, p_adjusted, estimate_quality_adjusted, p_quality_adjusted, genome_wide_fdr = fdr)]
paired[, family_fdr := bh(p)]
reference <- cell_pan[match(paired$symbol, symbol)]
paired[, cell_line_estimate := reference$estimate]
paired[, c("cell_line_mean_negative", "cell_line_mean_positive") := reference[, .(mean_negative, mean_positive)]]
paired[, cell_line_se := reference[, (ci_high - ci_low) / (2 * qt(0.975, df))]]
paired[, difference := estimate - cell_line_estimate]
paired[, difference_se := sqrt(se^2 + cell_line_se^2)]
paired[, difference_p := 2 * pnorm(-abs(difference / difference_se))]
paired[, outcome := replication_outcome(estimate, family_fdr, cell_line_estimate, difference_p, alpha = fdr_cut)]
paired <- merge(paired, chance[, .(symbol, nextgen_spread, expected_se, chance_at_strictest, chance_at_most_lenient)], by = "symbol", sort = FALSE)
print(paired[, .(symbol, estimate, se, p, family_fdr, cell_line_estimate, difference, difference_p, outcome, chance_at_strictest)])
partners <- primary_scan[feature %in% x_partner_features, .(feature, symbol, y_partner, estimate, se, smd, unit, p, fdr, estimate_adjusted, p_adjusted, label)]
write_ext(partners, "nextgen_xy_partner_scan")

# --- 5. continuous scan on every eligible screened male ---------------------------------------------------------
continuous_rows <- cohort[eligible == TRUE & !is.na(y_loss_fraction)]
continuous_frame <- as_frame(continuous_rows[best_of_group(ModelID, PatientID, ScreenNNMD)])
continuous_design <- build_design(continuous_frame, tested = "y_loss_fraction")
stopifnot(continuous_design$estimable)
continuous_frame <- continuous_frame[continuous_design$rows]
continuous <- weighted_limma_scan(t(effect[continuous_frame$ModelID, scan_features]), continuous_design$design, "y_loss_fraction", weights = weights[continuous_frame$ModelID])
continuous[, symbol := symbol_of(feature)]
continuous[, fdr := bh(p)]
setorder(continuous, p)
write_ext(continuous, "nextgen_continuous")
log_msg("continuous Y-loss fraction: %d models, %d genes at FDR < %.2f", nrow(continuous_frame), sum(continuous$fdr < fdr_cut), fdr_cut)

# --- 6. by cancer type ----------------------------------------------------------------------------------------------
by_type <- if (length(types) > 0) call_type_hits(type_screen(frame)) else data.table()
if (nrow(by_type) > 0) {
  by_type[, symbol := symbol_of(feature)]
  by_type[, named := feature %in% named_features]
  by_type[, focus := feature %in% x_partner_features]
  by_type[, passes_family_wise_bar := inside_p < inside_bar]
  by_type[, cell_line_selective := feature %in% gene_ann[cell_line_selective == TRUE, feature]]
  by_type <- merge(by_type, gene_ann[, .(feature, chromosome, y_partner)], by = "feature", sort = FALSE)
  setorder(by_type, inside_p)
  print(by_type[selective == TRUE, .(tests = .N, min_interaction_p = min(interaction_p, na.rm = TRUE), min_inside_p = min(inside_p, na.rm = TRUE), inside_fdr_hits = sum(inside_fdr < fdr_cut, na.rm = TRUE), type_specific = sum(type_specific, na.rm = TRUE)), by = .(type, n_neg_in, n_pos_in)])
}
write_ext(by_type, "nextgen_type_screen")

# --- calibration: observed beside shuffled ------------------------------------------------------------------------------
observed_named <- primary_scan[match(named_features, feature)]
observed_row <- data.table(
  run = "observed", pan_min_p = min(primary_scan$p), pan_fdr_hits = sum(primary_scan$fdr < fdr_cut), pan_associated = sum(primary_scan$associated),
  pan_below_0.001 = mean(primary_scan$p < 0.001), pan_below_1e5 = sum(primary_scan$p < 1e-5), pan_lambda = lambda,
  essential_lambda = lambda_of(primary_scan[feature %in% essential, p]), essential_mean_t = mean(primary_scan[feature %in% essential, t]), essential_positive = mean(primary_scan[feature %in% essential, estimate] > 0),
  named_min_p = min(observed_named$p), named_replicated = sum(paired$outcome == "replicated")
)
if (nrow(by_type) > 0) {
  tested <- by_type[selective == TRUE]
  observed_row <- cbind(observed_row, data.table(
    inside_min_p = min(tested$inside_p, na.rm = TRUE), interaction_min_p = min(tested$interaction_p, na.rm = TRUE), inside_fdr_hits = sum(tested$inside_fdr < fdr_cut, na.rm = TRUE),
    type_specific = sum(tested$type_specific, na.rm = TRUE), inside_below_0.001 = mean(tested$inside_p < 0.001, na.rm = TRUE), inside_lambda = lambda_of(tested$inside_p)
  ))
}
calibration <- rbind(observed_row, null_table)
write_ext(calibration, "nextgen_calibration")
print(calibration[run != "observed", lapply(.SD, function(v) round(mean(v), 4)), .SDcols = setdiff(names(calibration), "run")])

# --- hits --------------------------------------------------------------------------------------------------------------
hits <- rbind(
  primary_scan[associated == TRUE, .(kind = "all models", feature, symbol, context = "All models", estimate, p, fdr, label, named, passes_family_wise_bar)],
  if (nrow(by_type) > 0) by_type[type_specific == TRUE, .(kind = "specific to the cancer type", feature, symbol, context = type, estimate = inside, p = inside_p, fdr = inside_fdr, label = "type-specific", named, passes_family_wise_bar)],
  if (nrow(by_type) > 0) by_type[detected_inside_only == TRUE, .(kind = "inside the cancer type", feature, symbol, context = type, estimate = inside, p = inside_p, fdr = inside_fdr, label = "inside only", named, passes_family_wise_bar)],
  fill = TRUE
)
hits[, class := fifelse(named, "paired gene, named in advance", fifelse(passes_family_wise_bar, "discovery", "above the family-wise bar"))]
write_ext(hits, "nextgen_hits")
print(hits)

# --- 7. follow-up of every hit and of the X partners -------------------------------------------------------------------------
# Checks may only downgrade a hit. Alternative Y definitions are applied to these genes only.
definitions <- c(primary = "y_call", rna_only = "y_call_rna", dna_confirmed = "y_call_concordant", without_breadth_rule = "y_call_lenient", without_ceiling = "y_call_no_ceiling", x_checked = "y_call_x_checked")
frame_for <- function(exposure = "y_call", all_of_donor = FALSE, dome_only = FALSE, type_only = NULL) {
  out <- cohort[eligible == TRUE & !is.na(get(exposure))]
  if (!all_of_donor) out <- out[best_of_group(ModelID, PatientID, ScreenNNMD)]
  if (dome_only) out <- out[dome == 1]
  if (!is.null(type_only)) out <- out[OncotreeLineage == type_only]
  as_frame(out, exposure)
}
stopifnot(identical(frame_for()$ModelID, frame$ModelID))
one_fit <- function(values, fit_frame, covariates = character(), weighted = TRUE, separate = FALSE) {
  fit_frame <- copy(fit_frame)
  if (separate) {
    # Lineage and screen format as separate terms, in place of one term per stratum.
    formats <- lineage_design(fit_frame$ScreenType)[, -1, drop = FALSE]
    fit_frame[, OncotreeLineage := lineage]
    for (j in seq_len(ncol(formats))) fit_frame[, (paste0("format_", j)) := formats[, j]]
    covariates <- c(covariates, paste0("format_", seq_len(ncol(formats))))
  }
  design <- build_design(fit_frame, tested = "exposed", covariates = covariates)
  if (!design$estimable) return(list(estimate = NA_real_, se = NA_real_, p = NA_real_, n = 0L))
  fit_frame <- fit_frame[design$rows]
  weighted_ols(values[fit_frame$ModelID], fit_frame$exposed, design$design[, colnames(design$design) != "exposed", drop = FALSE], fit_frame$ModelID, weighted = weighted)
}
row_of <- function(fit, ...) data.table(..., estimate = fit$estimate, se = fit$se, p = fit$p, n = fit$n)
organoid_pan <- read_ext("organoid_pan_results", dir = ext_dir)
rnai_pan <- read_result("gene_pan_rnai")
look_up <- function(symbol_name) {
  one <- function(table, label, n_of) if (nrow(table) == 1) data.table(check = "other source", variant = label, estimate = table$estimate, se = NA_real_, p = table$p, n = n_of(table))
  rbind(
    one(cell_pan[symbol == symbol_name], "cell lines, CRISPR (Chronos)", function(d) as.integer(d$n_negative + d$n_positive)),
    one(organoid_pan[symbol == symbol_name], "Sanger organoids, CRISPR (scaled log fold change)", function(d) as.integer(d$n_negative + d$n_positive)),
    one(rnai_pan[symbol == symbol_name], "cell lines, RNAi (DEMETER2)", function(d) as.integer(d$n))
  )
}
follow <- function(feature, context) {
  inside_type <- if (context %in% types) context else NULL
  pan_level <- is.null(inside_type)
  # Type-level tests are judged on ranks over the whole cohort, unweighted, as in the
  # screen; the scan over all models is judged on the raw scale, weighted.
  values_of <- function(...) {
    if (pan_level) return(effect[, feature])
    everyone <- frame_for(...)
    setNames(rank_normalise(effect[everyone$ModelID, feature]), everyone$ModelID)
  }
  base <- frame_for(type_only = inside_type)
  values <- values_of()
  design <- build_design(base, tested = "exposed")
  base <- base[design$rows]
  checks <- followup_checks(values[base$ModelID], base$exposed, design$design[, colnames(design$design) != "exposed", drop = FALSE], base$OncotreeLineage,
    weights = if (pan_level) weights[base$ModelID] else NULL, n_perm = sp$followup_permutations, seed = ext_params$seed + 32, tie_tolerance = tie_tolerance)
  by_definition <- rbindlist(lapply(names(definitions), function(name) {
    row_of(one_fit(values_of(definitions[[name]]), frame_for(definitions[[name]], type_only = inside_type), weighted = pan_level), check = "Y-status definition", variant = name)
  }))
  other <- rbind(
    row_of(one_fit(values_of(dome_only = TRUE), frame_for(dome_only = TRUE, type_only = inside_type), weighted = pan_level), check = "screen format", variant = "dome organoids only"),
    row_of(one_fit(values, base, weighted = pan_level, separate = TRUE), check = "screen format", variant = "lineage and format as separate terms"),
    row_of(one_fit(values, base, weighted = FALSE), check = "weights", variant = "unweighted"),
    row_of(one_fit(values_of(all_of_donor = TRUE), frame_for(all_of_donor = TRUE, type_only = inside_type), weighted = pan_level), check = "all models of a donor", variant = "no one-per-donor rule"),
    row_of(one_fit(values, base, covariates = quality, weighted = pan_level), check = "covariates", variant = "plus screen NNMD"),
    row_of(one_fit(values, base, covariates = adjusted_covariates, weighted = pan_level), check = "covariates", variant = "plus WGD, TP53, aneuploidy")
  )
  summary_row <- data.table(
    check = c("permutation within stratum", "leave one model out", "leave one model out"),
    variant = c(sprintf("%d shuffles", sp$followup_permutations), "smallest |effect|", "largest p"),
    estimate = c(checks$estimate, sign(checks$estimate) * checks$loo_min_abs, NA_real_), se = NA_real_, p = c(checks$permutation_p, NA_real_, checks$loo_max_p), n = checks$n
  )
  out <- rbind(by_definition, other, summary_row, look_up(symbol_of(feature)))
  sign_holds <- function(name) isTRUE(sign(by_definition[variant == name, estimate]) == sign(checks$estimate))
  out[, robust := checks$permutation_p < sp$robust_permutation_p & checks$loo_sign_stable & checks$loo_max_p < sp$robust_leave_one_out_p & sign_holds("rna_only") & sign_holds("dna_confirmed")]
  cbind(data.table(feature = feature, symbol = symbol_of(feature), context = context, scale = if (pan_level) "raw" else "rank"), out)
}
targets <- unique(rbind(
  hits[, .(feature, context)],
  data.table(feature = x_partner_features, context = "All models"),
  if (length(types) > 0) CJ(feature = named_features, context = types)
))
followup <- rbindlist(lapply(seq_len(nrow(targets)), function(i) follow(targets$feature[i], targets$context[i])))
followup[, is_hit := paste(feature, context) %in% hits[, paste(feature, context)]]
write_ext(followup, "nextgen_hit_followup")
# The paired genes with the one check that can downgrade a replication. The model's p-values run small
# here (the shuffles above call a paired gene replicated about twice as often as they should), so each
# gene also gets a permutation p-value, from many more shuffles than the follow-up uses, and the rule
# fixed in advance is applied to those as well. Added after the first full run; it can only downgrade.
shuffled_t <- parallel::mclapply(named_features, function(target) followup_statistics(effect[frame$ModelID, target], frame$exposed, nuisance, frame$OncotreeLineage,
  weights = weights[frame$ModelID], n_perm = paired_permutations, seed = ext_params$seed + 34, loo = FALSE), mc.cores = min(n_cores, length(named_features)))
stopifnot(all(vapply(shuffled_t, is.list, logical(1))), all(vapply(shuffled_t, function(part) length(part$perm_t) == paired_permutations, logical(1))))
paired[, permutation_p := vapply(shuffled_t, function(part) (1 + sum(abs(part$perm_t) >= abs(part$t) * (1 - tie_tolerance), na.rm = TRUE)) / (length(part$perm_t) + 1), numeric(1))]
paired[, permutation_shuffles := paired_permutations]
paired[, permutation_fdr := bh(permutation_p)]
paired[, holds_under_permutation := outcome == "replicated" & permutation_fdr < fdr_cut]
paired[, outcome_permutation := replication_outcome(estimate, permutation_fdr, cell_line_estimate, difference_p, alpha = fdr_cut)]
shuffles <- followup[context == "All models" & check == "permutation within stratum"]
paired[, robust := shuffles$robust[match(feature, shuffles$feature)]]
stopifnot(!anyNA(paired$permutation_p), !anyNA(paired$robust))
write_ext(paired, "nextgen_paired_genes")
print(paired[, .(symbol, outcome, permutation_p, permutation_fdr, holds_under_permutation, outcome_permutation, robust)])
print(followup[feature %in% named_features & (check %in% c("permutation within stratum", "leave one model out") | variant %in% c("rna_only", "dna_confirmed", "dome organoids only", "unweighted")), .(symbol, context, check, variant, estimate = round(estimate, 3), p = signif(p, 3), n, robust)])

shuffled <- calibration[run != "observed"]
jsonlite::write_json(list(
  models = nrow(frame), y_negative = sum(frame$exposed == 1), y_positive = sum(frame$exposed == 0), donors = uniqueN(frame$PatientID),
  genes_tested = nrow(primary_scan), fdr_hits = sum(primary_scan$fdr < fdr_cut), associated = sum(primary_scan$associated), discoveries = sum(primary_scan$discovery),
  loy_specific = sum(primary_scan$label == "LOY-specific"), lambda = lambda, unit = unique(round(primary_scan$unit, 12)),
  paired_outcomes = as.list(setNames(paired$outcome, paired$symbol)), paired_hold_under_permutation = as.list(setNames(paired$holds_under_permutation, paired$symbol)),
  paired_outcomes_permutation = as.list(setNames(paired$outcome_permutation, paired$symbol)), paired_permutations = paired_permutations,
  lambda_quality_adjusted = lambda_of(primary_scan$p_quality_adjusted), lambda_essential_adjusted = lambda_of(primary_scan$p_essential_adjusted), lambda_adjusted = lambda_of(primary_scan$p_adjusted),
  essential = list(genes = length(essential), lambda = observed_row$essential_lambda, mean_t = observed_row$essential_mean_t, positive = observed_row$essential_positive,
    lambda_other_genes = lambda_of(primary_scan[!feature %in% essential, p]), runs_at_or_above_lambda = sum(shuffled$essential_lambda >= observed_row$essential_lambda),
    runs_at_or_above_mean_t = sum(shuffled$essential_mean_t >= observed_row$essential_mean_t)),
  essential_shift = list(estimate = shift$estimate, p = shift$p, permutation_p = shift$permutation_p), quality_shift = list(estimate = quality_shift$estimate, p = quality_shift$p, permutation_p = quality_shift$permutation_p),
  adjusted_models = nrow(adjusted$frame), continuous_models = nrow(continuous_frame), continuous_fdr_hits = sum(continuous$fdr < fdr_cut),
  types = types, type_tests = if (nrow(by_type) > 0) by_type[selective == TRUE, .N] else 0L, type_specific = if (nrow(by_type) > 0) sum(by_type$type_specific, na.rm = TRUE) else 0L,
  inside_only = if (nrow(by_type) > 0) sum(by_type$detected_inside_only, na.rm = TRUE) else 0L,
  null = list(
    runs = nrow(shuffled), pan_bar = pan_bar, inside_bar = inside_bar, runs_with_fdr_hit = sum(shuffled$pan_fdr_hits > 0), runs_with_associated = sum(shuffled$pan_associated > 0),
    runs_below_1e5 = sum(shuffled$pan_below_1e5 > 0), expected_below_1e5 = 1e-5 * nrow(primary_scan), mean_below_1e5 = mean(shuffled$pan_below_1e5),
    lambda_median = median(shuffled$pan_lambda), lambda_range = range(shuffled$pan_lambda), runs_at_or_above_lambda = sum(shuffled$pan_lambda >= lambda), named_replicated_runs = sum(shuffled$named_replicated > 0),
    inside_runs_with_fdr_hit = if (length(types) > 0) sum(shuffled$inside_fdr_hits > 0) else NA, type_specific_runs = if (length(types) > 0) sum(shuffled$type_specific > 0) else NA,
    inside_lambda_median = if (length(types) > 0) median(shuffled$inside_lambda) else NA
  ),
  controls_pass = all(controls$pass), weights_range = range(weights), trial = trial, labels_shuffled = trial,
  parameters_sha256 = ext_section_hash(c("shared", "nextgen"))
), file.path(ext_out_dir, "nextgen_scan_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("NextGen screens done: %d associated over all models (%d discoveries), paired genes %s; %d type-specific, %d inside only; %d of %d shuffles gave a scan hit",
  sum(primary_scan$associated), sum(primary_scan$discovery), paste(sprintf("%s %s (%s with permutation p)", paired$symbol, paired$outcome, paired$outcome_permutation), collapse = ", "),
  if (nrow(by_type) > 0) sum(by_type$type_specific, na.rm = TRUE) else 0L, if (nrow(by_type) > 0) sum(by_type$detected_inside_only, na.rm = TRUE) else 0L, sum(shuffled$pan_fdr_hits > 0), nrow(shuffled))
