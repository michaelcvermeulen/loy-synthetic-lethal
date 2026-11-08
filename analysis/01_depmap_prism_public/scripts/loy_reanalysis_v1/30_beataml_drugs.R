# Stage 30: BeatAML2, drug response against the Y-loss fraction.
# One specimen per male patient. Per inhibitor:
#   rank-normalised AUC ~ f + specimen type + disease stage + age + core-binding-factor fusion + study wave
# AUC is the area under the fitted curve, lower = more sensitive, so a negative
# slope means "Y loss goes with sensitivity", as in the cell-line drug stages.
# The slope is carried by 10 to 20 specimens, where the robust test of the
# cell-line scan is not reliable (see 13_context_drugs.R). The p-value is
# therefore a Freedman-Lane permutation of the residuals under the model without
# f; the classical and the robust p-values are reported beside it.
# In order: the known-biology control (the stage stops if it fails), shuffled
# exposure, the scan, drug families, the comparison with the cell-line screens,
# and the follow-up of every hit and of the leading inhibitor.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))
source(file.path(script_dir, "extension_helpers.R"))
assert_freeze()
assert_extension_freeze("beataml")

bp <- ext_params$beataml
sp <- ext_params$shared
fdr_cut <- sp$fdr
trial <- Sys.getenv("LOY_EXTENSION_TRIAL") == "1"
n_perm <- if (trial) 2000L else bp$permutations
cohort <- read_ext("beataml_cohort")
auc <- readRDS(file.path(ext_out_dir, "beataml_auc.rds"))
test_counts <- read_ext("beataml_test_counts")
families <- read_ext("beataml_families")
matches <- read_ext("beataml_compound_match")
tested <- test_counts[tested == TRUE, inhibitor]

dummies <- function(x) {
  x <- factor(as.character(x))
  if (nlevels(x) < 2) return(NULL)
  model.matrix(~x)[, -1, drop = FALSE]
}
nuisance_of <- function(rows, age_and_fusion = TRUE, extra = NULL) {
  parts <- list(intercept = rep(1, nrow(rows)), dummies(rows$specimen_type), dummies(rows$stage), dummies(rows$wave))
  if (age_and_fusion) parts <- c(parts, list(age = rows$age, cbf_fusion = rows$cbf_fusion))
  for (name in extra) parts[[name]] <- rows[[name]]
  do.call(cbind, parts)
}
rank_columns <- function(mat) {
  out <- apply(mat, 2, rank_normalise)
  dimnames(out) <- dimnames(mat)
  out
}

# --- the frame: one specimen per patient, with every covariate of the main model -------------------
frame <- cohort[in_primary == TRUE & !is.na(age)][order(SpecimenID)]
response <- rank_columns(auc[frame$SpecimenID, tested, drop = FALSE])
raw <- auc[frame$SpecimenID, tested, drop = FALSE]
z_main <- nuisance_of(frame)
log_msg("BeatAML2 drug screen: %d specimens (one per patient; %d of the primary cohort lack age), %d with f >= %.1f, %d with f > 0; %d inhibitors; %d permutations",
  nrow(frame), cohort[in_primary == TRUE & is.na(age), .N], sum(frame$f >= bp$high_loss_fraction), bp$high_loss_fraction, sum(frame$f > 0), length(tested), n_perm)

# --- 1. known-biology control, read before any Y result --------------------------------------------------
flt3_inhibitors <- grep(paste(bp$flt3_control_patterns, collapse = "|"), colnames(auc), value = TRUE)
control_frame <- cohort[in_primary == TRUE & !is.na(age) & !is.na(flt3_itd)][order(SpecimenID)]
control_response <- rank_columns(auc[control_frame$SpecimenID, flt3_inhibitors, drop = FALSE])
controls <- scan_ols(control_response, control_frame$flt3_itd, nuisance_of(control_frame))
setnames(controls, "feature", "inhibitor")
controls[, flt3_itd_specimens := sum(control_frame$flt3_itd == 1)]
controls[, pass := median(t) < 0 & min(p) < bp$flt3_control_p_max]
write_ext(controls, "beataml_controls")
print(controls)
if (!all(controls$pass)) stop("the FLT3 control failed")

# --- the scan --------------------------------------------------------------------------------------------
perms <- permutation_matrix(nrow(frame), n_perm, seed = ext_params$seed + 30)
# One inhibitor: least squares on ranks, Freedman-Lane permutation, robust test, raw-scale slope.
fit_inhibitor <- function(j, exposure, z, permutations, keep_null = FALSE) {
  y <- response[, j]
  ok <- which(is.finite(y) & is.finite(exposure))
  classical <- fit_ols(y[ok], exposure[ok], z[ok, , drop = FALSE])
  out <- data.table(inhibitor = j, estimate = classical$estimate, se = classical$se, t = classical$t, p_classical = classical$p, n = classical$n)
  if (!is.null(permutations)) {
    shuffled <- freedman_lane(y[ok], exposure[ok], z[ok, , drop = FALSE], restrict_permutations(permutations, ok))
    out[, p := shuffled$p]
    if (keep_null) return(list(row = out, null = shuffled$t_perm))
  }
  out
}
run_scan <- function(exposure, permutations, keep_null = FALSE) {
  fits <- parallel::mclapply(tested, fit_inhibitor, exposure = exposure, z = z_main, permutations = permutations, keep_null = keep_null, mc.cores = n_cores)
  failed <- vapply(fits, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("inhibitor fits failed: ", paste(tested[failed], collapse = ", "))
  fits
}

# --- 2. shuffled exposure: what the procedure yields when f carries nothing ---------------------------------
# f is shuffled within specimen type and disease stage; fewer permutations are enough for a count of hits.
null_perms <- perms[, seq_len(min(ncol(perms), 20000)), drop = FALSE]
calibration <- rbindlist(lapply(seq_len(sp$null_runs), function(b) {
  set.seed(ext_params$seed + 300 * b)
  shuffled_f <- permute_within(frame$f, paste(frame$specimen_type, frame$stage))
  null <- rbindlist(run_scan(shuffled_f, null_perms))
  null[, fdr := bh(p)]
  data.table(run = paste0("shuffled_", b), tests = nrow(null), fdr_hits = sum(null$fdr < fdr_cut), min_p = min(null$p), below_0.05 = mean(null$p < 0.05), below_0.001 = mean(null$p < 0.001))
}))

# --- 3. the scan with the real exposure ------------------------------------------------------------------------
fits <- run_scan(frame$f, perms, keep_null = TRUE)
results <- rbindlist(lapply(fits, `[[`, "row"))
null_t <- vapply(fits, `[[`, numeric(n_perm), "null")
colnames(null_t) <- tested
rm(fits)
results[, fdr := bh(p)]
robust <- scan_robust(response, frame$f, z_main, cluster = frame$PatientID)
results[, p_robust := robust$p[match(inhibitor, robust$feature)]]
# Slope on the raw scale: AUC standard deviations per complete Y loss.
raw_fit <- scan_ols(raw, frame$f, z_main)
results[, slope_sd := raw_fit$estimate[match(inhibitor, raw_fit$feature)] / apply(raw, 2, sd, na.rm = TRUE)[inhibitor]]
results[, high_loss := colSums(is.finite(response[, inhibitor, drop = FALSE]) & frame$f >= bp$high_loss_fraction)]
results[, any_loss := colSums(is.finite(response[, inhibitor, drop = FALSE]) & frame$f > 0)]
# Beside the main model: without age and fusion, and without the t(8;21) cases.
plain <- scan_ols(response, frame$f, nuisance_of(frame, age_and_fusion = FALSE))
no_t821 <- scan_ols(response[frame$t8_21 == 0, , drop = FALSE], frame$f[frame$t8_21 == 0], nuisance_of(frame[t8_21 == 0]))
results[, c("estimate_without_age_fusion", "p_without_age_fusion") := plain[match(results$inhibitor, feature), .(estimate, p)]]
results[, c("estimate_without_t8_21", "p_without_t8_21") := no_t821[match(results$inhibitor, feature), .(estimate, p)]]
# Adjusted model for the label: TP53, FLT3-ITD, NPM1 and the specimen's mean response over all inhibitors.
frame[, mean_response := rowMeans(response, na.rm = TRUE)]
adjusted_rows <- which(stats::complete.cases(frame[, .(tp53, flt3_itd, npm1, mean_response)]))
adjusted <- scan_ols(response[adjusted_rows, , drop = FALSE], frame$f[adjusted_rows], nuisance_of(frame[adjusted_rows], extra = c("tp53", "flt3_itd", "npm1", "mean_response")))
results[, c("estimate_adjusted", "p_adjusted", "n_adjusted") := adjusted[match(results$inhibitor, feature), .(estimate, p, n)]]
general <- fit_ols(frame$mean_response, frame$f, z_main)
log_msg("mean response over all inhibitors against f: slope %.3f rank units (p = %.3g)", general$estimate, general$p)
results[, associated := fdr < fdr_cut]
results[associated == TRUE, fdr_adjusted := bh(p_adjusted)]
results[, survives_adjustment := associated & sign(estimate_adjusted) == sign(estimate) & fdr_adjusted < fdr_cut]
results[, label := fifelse(associated & survives_adjustment %in% TRUE, "LOY-specific", fifelse(associated, "LOY-associated", "none"))]
results[, direction := fifelse(estimate < 0, "more sensitive with Y loss", "less sensitive with Y loss")]
results[, matched_cell_line_compound := matches$compound_name[match(inhibitor, matches$inhibitor)]]
setorder(results, p)
write_ext(results, "beataml_drug_results")
lambda <- median(qchisq(results$p_classical, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
log_msg("scan: %d inhibitors, %d at FDR < %.2f, smallest permutation p %.2g (%s), lambda %.2f on the classical p-values",
  nrow(results), sum(results$associated), fdr_cut, min(results$p), results$inhibitor[1], lambda)
print(results[1:12, .(inhibitor, n, high_loss, estimate = round(estimate, 2), slope_sd = round(slope_sd, 2), p = signif(p, 3), p_classical = signif(p_classical, 3), p_robust = signif(p_robust, 3), fdr = round(fdr, 3), label)])
calibration <- rbind(data.table(run = "observed", tests = nrow(results), fdr_hits = sum(results$associated), min_p = min(results$p), below_0.05 = mean(results$p < 0.05), below_0.001 = mean(results$p < 0.001)), calibration)
write_ext(calibration, "beataml_calibration")
print(calibration)

# --- 4. drug families: mean statistic against the same permutations ------------------------------------------------
family_sets <- split(families[inhibitor %in% tested, inhibitor], families[inhibitor %in% tested, family])
family_sets <- family_sets[lengths(family_sets) >= bp$family_min_inhibitors]
observed_t <- setNames(results$t, results$inhibitor)
family_results <- rbindlist(lapply(names(family_sets), function(name) {
  members <- family_sets[[name]]
  statistic <- mean(observed_t[members])
  null <- rowMeans(null_t[, members, drop = FALSE])
  data.table(family = name, inhibitors = length(members), mean_t = statistic, p = (1 + sum(abs(null) >= abs(statistic))) / (length(null) + 1))
}))
family_results[, fdr := bh(p)]
setorder(family_results, p)
write_ext(family_results, "beataml_family_results")
print(family_results[1:min(.N, 8)])

# --- 5. the inhibitors shared with the cell-line screens: do the slopes agree in sign? -------------------------------
# The list was fixed in stage 29 (name matches). One test per cell-line reference: the share of
# inhibitors whose slope has the sign of the cell-line estimate, against the same share under the permutations.
cell_pan <- read_result("drug_pan_results")[, .(cell_line_standardized = mean(standardized, na.rm = TRUE), cell_line_datasets = .N, cell_line_min_p = min(p)), by = .(compound_name = toupper(compound_name))]
myeloid <- readRDS(file.path(derived_dir, "context_drug_screen.rds"))
setDT(myeloid)
myeloid <- myeloid[context_id == "lineage__Myeloid", .(myeloid_inside_standardized = mean(inside_standardized, na.rm = TRUE), myeloid_datasets = .N), by = .(compound_name = toupper(compound_name))]
comparison <- merge(matches[inhibitor %in% tested], results[, .(inhibitor, estimate, t, p, slope_sd)], by = "inhibitor")
comparison <- merge(comparison, cell_pan, by = "compound_name", all.x = TRUE)
comparison <- merge(comparison, myeloid, by = "compound_name", all.x = TRUE)
agreement_test <- function(reference, label) {
  rows <- comparison[is.finite(get(reference)) & get(reference) != 0]
  if (nrow(rows) < 5) return(data.table(reference = label, inhibitors = nrow(rows), sign_agreement = NA_real_, expected_under_permutation = NA_real_, p = NA_real_, rank_correlation = NA_real_))
  observed <- mean(sign(rows$t) == sign(rows[[reference]]))
  null <- colMeans(sign(t(null_t[, rows$inhibitor, drop = FALSE])) == sign(rows[[reference]]))
  data.table(reference = label, inhibitors = nrow(rows), sign_agreement = observed, expected_under_permutation = mean(null),
    p = (1 + sum(abs(null - mean(null)) >= abs(observed - mean(null)))) / (length(null) + 1),
    rank_correlation = cor(rows$t, rows[[reference]], method = "spearman"))
}
comparison_summary <- rbind(agreement_test("cell_line_standardized", "cell lines, pan-cancer"), agreement_test("myeloid_inside_standardized", "cell lines, inside myeloid"))
write_ext(comparison[order(p)], "beataml_cell_line_comparison")
write_ext(comparison_summary, "beataml_cell_line_comparison_summary")
print(comparison_summary)

# --- 6. follow-up of every hit, and of the leading inhibitor whether or not it passes ------------------------------------
# Checks may only downgrade a hit.
everyone <- cohort[eligible == TRUE & has_drug == TRUE & !is.na(f) & !is.na(age)][order(SpecimenID)]
# Added after the first result (selumetinib, a MEK inhibitor): MEK-inhibitor response in AML follows
# RAS-pathway mutations, so their share by Y loss is tabulated and they enter two extra checks.
# Read from the exome calls; a specimen without exome data is unknown. Not a pre-specified check.
ras_genes <- c("NRAS", "KRAS", "PTPN11")
exome_calls <- fread(file.path(raw_dir, "drug_response", "beataml2", "beataml_wes_wv1to4_mutations_dbgap.txt"), select = c("dbgap_sample_id", "symbol"), showProgress = FALSE)
frame[, ras_pathway := fifelse(dbgap_dnaseq_sample %in% exome_calls[symbol %in% ras_genes, dbgap_sample_id], 1, fifelse(is.na(tp53), NA_real_, 0))]
ras_table <- frame[, .(specimens = .N, ras_pathway_mutant = sum(ras_pathway == 1, na.rm = TRUE), ras_pathway_wild_type = sum(ras_pathway == 0, na.rm = TRUE), unknown = sum(is.na(ras_pathway))),
  by = .(y_loss = fifelse(f >= bp$high_loss_fraction, sprintf("f at least %.1f", bp$high_loss_fraction), fifelse(f > 0, "f above the floor", "no Y loss")))]
ras_table[, fisher_p := fisher.test(as.matrix(frame[!is.na(ras_pathway), table(f >= bp$high_loss_fraction, ras_pathway)]))$p.value]
write_ext(ras_table, "beataml_ras_by_y_loss")
print(ras_table)
follow <- function(j) {
  variant <- function(check, name, rows, exposure, extra = NULL, cluster = NULL) {
    y <- rank_normalise(auc[rows$SpecimenID, j])
    z <- nuisance_of(rows, extra = extra)
    fit <- if (is.null(cluster)) fit_ols(y, exposure, z) else fit_robust(y, exposure, z, cluster = cluster)
    data.table(check = check, variant = name, estimate = fit$estimate, p = fit$p, n = fit$n, exposed = sum(is.finite(y) & exposure > 0, na.rm = TRUE))
  }
  y <- response[, j]
  ok <- is.finite(y)
  checks <- followup_checks(y[ok], frame$f[ok], z_main[ok, , drop = FALSE], paste(frame$specimen_type, frame$stage)[ok], n_perm = sp$followup_permutations, seed = ext_params$seed + 31)
  confirmed <- frame[f == 0 | evidence_tier == "karyotype -Y"]
  with_karyotype <- frame[!is.na(karyotype_y_fraction)]
  initial <- frame[stage == "Initial Diagnosis"]
  blasty <- frame[!is.na(blasts) & blasts >= bp$blast_minimum]
  with_fab <- frame[!is.na(monocytic)]
  with_ras <- frame[!is.na(ras_pathway)]
  ras_wild_type <- frame[ras_pathway %in% 0]
  out <- rbind(
    variant("exposure", "main: Y relative to X partner, floored", frame, frame$f),
    variant("exposure", "Y genes alone", frame, frame$f_y_only),
    variant("exposure", "without the floor", frame, frame$f_no_floor),
    variant("exposure", sprintf("f at least %.1f against the rest", bp$high_loss_fraction), frame, as.numeric(frame$f >= bp$high_loss_fraction)),
    variant("exposure", "karyotype fraction of cells without Y", with_karyotype, with_karyotype$karyotype_y_fraction),
    variant("specimens", "Y loss confirmed by karyotype, or none", confirmed, confirmed$f),
    variant("specimens", "initial diagnosis only", initial, initial$f),
    variant("specimens", sprintf("blasts at least %d%%", bp$blast_minimum), blasty, blasty$f),
    variant("specimens", "all specimens, patient-clustered errors", everyone, everyone$f, cluster = everyone$PatientID),
    variant("covariates", "plus monocytic morphology", with_fab, with_fab$f, extra = "monocytic"),
    variant("added after the result", "plus RAS-pathway mutation (NRAS, KRAS, PTPN11)", with_ras, with_ras$f, extra = "ras_pathway"),
    variant("added after the result", "RAS-pathway wild-type specimens only", ras_wild_type, ras_wild_type$f),
    data.table(check = c("permutation within specimen type and stage", "leave one specimen out", "leave one specimen out"),
      variant = c(sprintf("%d shuffles of f", sp$followup_permutations), "smallest |slope|", "largest p"),
      estimate = c(checks$estimate, sign(checks$estimate) * checks$loo_min_abs, NA_real_), p = c(checks$permutation_p, NA_real_, checks$loo_max_p), n = checks$n, exposed = NA_integer_)
  )
  sign_holds <- function(name) isTRUE(sign(out[variant == name, estimate]) == sign(checks$estimate))
  main_p <- results[inhibitor == j, p]
  out[, robust := main_p < sp$robust_permutation_p & checks$loo_sign_stable & checks$loo_max_p < sp$robust_leave_one_out_p & sign_holds("Y genes alone") & sign_holds("Y loss confirmed by karyotype, or none")]
  cbind(data.table(inhibitor = j), out)
}
targets <- union(results[associated == TRUE, inhibitor], results$inhibitor[1])
followup <- rbindlist(lapply(targets, follow))
followup[, is_hit := inhibitor %in% results[associated == TRUE, inhibitor]]
write_ext(followup, "beataml_hit_followup")
print(followup[, .(inhibitor, check, variant, estimate = round(estimate, 3), p = signif(p, 3), n, exposed, robust)])

# The leading inhibitor specimen by specimen, for the figure.
leading <- results$inhibitor[1]
leading_table <- frame[, .(inhibitor = leading, SpecimenID, f, f_no_floor, evidence_tier, ras_pathway, stage, specimen_type, auc = raw[, leading], rank_response = response[, leading])]
write_ext(leading_table[!is.na(auc)], "beataml_leading_inhibitor")

# --- what the screen could detect at best -------------------------------------------------------------------------------------
# Standard error of the slope per unit residual SD, for the median inhibitor.
units <- vapply(tested, function(j) {
  ok <- is.finite(response[, j])
  1 / sqrt(sum(qr.resid(qr(z_main[ok, , drop = FALSE]), frame$f[ok])^2))
}, numeric(1))
unit <- median(units)
df <- as.integer(median(results$n)) - qr(z_main)$rank - 1
critical <- bh_next_critical(results$p, fdr_cut)
power <- data.table(
  analysis = "Primary AML, drug response", source = "BeatAML2", n_negative = sum(frame$f >= bp$high_loss_fraction), n_positive = sum(frame$f == 0),
  unit = unit, df = df, family = nrow(results), named = analytic_mde(unit, df, 0.05, sp$power), family_critical_p = critical, genome_wide = analytic_mde(unit, df, critical, sp$power)
)
write_ext(power, "beataml_power")
print(power)

jsonlite::write_json(list(
  specimens = nrow(frame), high_loss = sum(frame$f >= bp$high_loss_fraction), any_loss = sum(frame$f > 0), inhibitors = nrow(results), permutations = n_perm,
  fdr_hits = sum(results$associated), loy_specific = sum(results$label == "LOY-specific"), min_p = min(results$p), leading_inhibitor = results$inhibitor[1], lambda_classical = lambda,
  shuffled_fdr_hits = calibration[run != "observed", fdr_hits], control_pass = all(controls$pass),
  families = nrow(family_results), family_fdr_hits = sum(family_results$fdr < fdr_cut),
  matched_and_tested = nrow(comparison), general_response = list(estimate = general$estimate, p = general$p),
  unit = unit, trial = trial, parameters_sha256 = ext_section_hash(c("shared", "beataml"))
), file.path(ext_out_dir, "beataml_scan_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("BeatAML2 drug screen done: %d of %d inhibitors at FDR < %.2f; %d of %d families; shuffled exposure gave %s hits",
  sum(results$associated), nrow(results), fdr_cut, sum(family_results$fdr < fdr_cut), nrow(family_results), paste(calibration[run != "observed", fdr_hits], collapse = ", "))
