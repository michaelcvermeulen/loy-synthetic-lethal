# Stage 11 (plan Stage 4): pan-cancer drug sensitivities associated with Y loss.
# Every metric here is "lower = more sensitive" (PRISM log2 fold change; AUC in
# PRISM secondary, GDSC and CTRP), so a negative estimate always means
# "Y-negative lines are more sensitive".
# 1. Per dataset: response ~ lineage (+ screen) + Y_negative, donor-clustered HC3, BH within dataset.
# 2. Across datasets: stacked mixed model per shared compound, cell line as a random effect.
# 3. Mechanism and target classes: mean compound statistic against a within-lineage permutation null.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
assert_freeze()
controls <- jsonlite::fromJSON(file.path(derived_dir, "pipeline_controls.json"))
stopifnot(controls$braf_drug_control_pass, controls$tp53_drug_control_pass, controls$drug_calibration_pass)

dp <- params$drug_scan
cohort <- read_result("cohort")
drug_ann <- read_result("drug_feature_annotation")
datasets <- unique(drug_ann$dataset)
adjusted_covariates <- c("WGD", "tp53_altered", "Aneuploidy")
secondary_screen <- readRDS(file.path(matrix_dir, "prism_secondary_screen.rds"))

# Screen used for each PRISM secondary curve; a covariate only when a compound was read in more than one screen.
screen_columns <- function(feature, lines) {
  screen <- secondary_screen[[feature]][match(lines, secondary_screen$depmap_id)]
  levels <- sort(unique(na.omit(screen)))
  if (length(levels) < 2) return(NULL)
  # Lines without a curve keep NA here; their response is missing too, so the fit drops them.
  vapply(levels[-1], function(level) as.numeric(screen == level), numeric(length(screen)))
}

scan_dataset <- function(name, mat, frame, covariates = character()) {
  design <- build_design(frame, tested = "exposed", covariates = covariates)
  frame <- frame[design$rows]
  nuisance <- design$design[, colnames(design$design) != "exposed", drop = FALSE]
  y <- mat[frame$ModelID, , drop = FALSE]
  if (name != "prism_secondary") return(scan_robust(y, frame$exposed, nuisance, cluster = frame$PatientID))
  fits <- parallel::mclapply(colnames(y), function(feature) {
    extra <- screen_columns(feature, frame$ModelID)
    z <- if (is.null(extra)) nuisance else cbind(nuisance, extra)
    fit_robust(y[, feature], frame$exposed, z, cluster = frame$PatientID)
  }, mc.cores = n_cores)
  out <- rbindlist(lapply(fits, as.data.table))
  out[, feature := colnames(y)]
  setcolorder(out, "feature")
  out[]
}

# --- 1. per-dataset scans ----------------------------------------------------
frames <- list()
matrices <- list()
per_dataset <- rbindlist(lapply(datasets, function(name) {
  mat <- load_matrix(name)
  frame <- analysis_frame(cohort, mat)
  features <- intersect(drug_ann[dataset == name & active == TRUE, feature], colnames(mat))
  mat <- mat[, features, drop = FALSE]
  frames[[name]] <<- frame
  matrices[[name]] <<- mat
  scan <- scan_dataset(name, mat, frame)
  adjusted <- scan_dataset(name, mat, frame, covariates = adjusted_covariates)
  scan[, c("estimate_adjusted", "p_adjusted", "n_adjusted") := adjusted[, .(estimate, p, n)]]
  y <- mat[frame$ModelID, , drop = FALSE]
  scan[, response_sd := apply(y, 2, sd, na.rm = TRUE)]
  scan[, mean_negative := colMeans(y[frame$exposed == 1, , drop = FALSE], na.rm = TRUE)]
  scan[, mean_positive := colMeans(y[frame$exposed == 0, , drop = FALSE], na.rm = TRUE)]
  scan <- scan[n_exposed >= dp$min_per_group & n_unexposed >= dp$min_per_group & is.finite(p)]
  scan[, dataset := name]
  scan[, fdr := bh(p)]
  log_msg("%s: %d lines (%d / %d), %d compounds tested, %d at FDR < %.2f, smallest p %.2g",
    name, nrow(frame), sum(frame$exposed == 1), sum(frame$exposed == 0), nrow(scan), sum(scan$fdr < dp$fdr), dp$fdr, min(scan$p))
  scan
}))
per_dataset[, standardized := estimate / response_sd]
per_dataset[, associated := fdr < dp$fdr]
per_dataset[associated == TRUE, fdr_adjusted := bh(p_adjusted)]
per_dataset[, survives_adjustment := associated & sign(estimate_adjusted) == sign(estimate) & fdr_adjusted < dp$fdr]
per_dataset[, label := fifelse(associated & survives_adjustment %in% TRUE, "LOY-specific", fifelse(associated, "LOY-associated", "none"))]
per_dataset <- merge(
  drug_ann[, .(dataset, feature, compound_id, compound_name, mechanism, targets, screen)], per_dataset,
  by = c("dataset", "feature"), sort = FALSE
)
setnames(per_dataset, c("n_exposed", "n_unexposed"), c("n_negative", "n_positive"))
setorder(per_dataset, p)
write_result(per_dataset, "drug_pan_results")
write_result(per_dataset[associated == TRUE], "drug_pan_hits")
print(per_dataset[, .(tested = .N, fdr_hits = sum(associated), min_p = min(p), min_fdr = min(fdr), lambda = median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)), by = dataset])
print(per_dataset[1:15, .(dataset, compound_name, mechanism = substr(mechanism, 1, 28), n_negative, n_positive, standardized, p, fdr)])

# --- 2. stacked model across datasets ----------------------------------------
suppressPackageStartupMessages(library(lme4))
family_of <- function(dataset) sub("_.*$", "", dataset)
tested <- per_dataset[!is.na(compound_id)]
shared <- tested[, .(families = uniqueN(family_of(dataset)), features = .N), by = .(compound_id, compound_name)][families >= 2]
stack_compound <- function(id) {
  rows <- tested[compound_id == id]
  long <- rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    frame <- frames[[rows$dataset[i]]]
    y <- matrices[[rows$dataset[i]]][frame$ModelID, rows$feature[i]]
    data.table(
      ModelID = frame$ModelID, PatientID = frame$PatientID, lineage = frame$OncotreeLineage, exposed = frame$exposed,
      source = paste(rows$dataset[i], rows$feature[i]), z = (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
    )
  }))
  long[is.finite(z)]
}
fit_stacked <- function(id) {
  long <- stack_compound(id)
  out <- data.table(compound_id = id, lines = uniqueN(long$ModelID), observations = nrow(long), sources = uniqueN(long$source))
  mixed <- tryCatch({
    full <- suppressMessages(suppressWarnings(lmer(z ~ exposed + lineage + source + (1 | ModelID), data = long, REML = FALSE)))
    null <- suppressMessages(suppressWarnings(lmer(z ~ lineage + source + (1 | ModelID), data = long, REML = FALSE)))
    statistic <- as.numeric(2 * (logLik(full) - logLik(null)))
    list(estimate = unname(fixef(full)["exposed"]), se = sqrt(diag(as.matrix(vcov(full))))[["exposed"]], lrt = statistic, p = pchisq(statistic, 1, lower.tail = FALSE))
  }, error = function(e) list(estimate = NA_real_, se = NA_real_, lrt = NA_real_, p = NA_real_))
  # Cross-check with the cell line as the unit: one mean z per line, robust test.
  per_line <- long[, .(z = mean(z), lineage = lineage[1], exposed = exposed[1], PatientID = PatientID[1]), by = ModelID]
  robust <- fit_robust(per_line$z, per_line$exposed, lineage_design(per_line$lineage), cluster = per_line$PatientID)
  cbind(out, as.data.table(mixed), data.table(estimate_line_mean = robust$estimate, p_line_mean = robust$p))
}
consensus <- rbindlist(parallel::mclapply(shared$compound_id, fit_stacked, mc.cores = n_cores))
consensus <- merge(shared, consensus, by = "compound_id")
consensus[, fdr := bh(p)]
consensus <- merge(consensus, unique(tested[, .(compound_id, mechanism, targets)])[!duplicated(compound_id)], by = "compound_id", all.x = TRUE)
setorder(consensus, p)
write_result(consensus, "drug_pan_consensus")
log_msg("stacked model: %d shared compounds, %d at FDR < %.2f, smallest p %.2g", nrow(consensus), sum(consensus$fdr < dp$fdr, na.rm = TRUE), dp$fdr, min(consensus$p, na.rm = TRUE))
print(consensus[1:10, .(compound_name, families, lines, estimate, p, fdr, p_line_mean)])

# --- 3. mechanism and target classes -----------------------------------------
# Lineage-centred t statistic for every compound, with missing responses masked.
prepare_masked <- function(mat, frame) {
  y <- mat[frame$ModelID, , drop = FALSE]
  available <- !is.na(y)
  y[!available] <- 0
  lineage <- model.matrix(~ 0 + factor(frame$OncotreeLineage))
  counts <- crossprod(lineage, available * 1)
  centred <- (y - lineage %*% (crossprod(lineage, y) / pmax(counts, 1))) * available
  list(available = available * 1, lineage = lineage, counts = counts, centred = centred,
    syy = colSums(centred^2), df = colSums(available) - colSums(counts > 0) - 1)
}
masked_t <- function(prep, x) {
  x_mean <- crossprod(prep$lineage, prep$available * x) / pmax(prep$counts, 1)
  x_centred <- (x - prep$lineage %*% x_mean) * prep$available
  sxx <- colSums(x_centred^2)
  estimate <- colSums(x_centred * prep$centred) / sxx
  estimate / sqrt((prep$syy - estimate^2 * sxx) / prep$df / sxx)
}
tested_keys <- per_dataset[, paste(dataset, feature)]
prepared <- lapply(datasets, function(name) {
  features <- per_dataset[dataset == name, feature]
  prepare_masked(matrices[[name]][, features, drop = FALSE], frames[[name]])
})
names(prepared) <- datasets
feature_keys <- unlist(lapply(datasets, function(name) paste(name, colnames(prepared[[name]]$centred))))
all_lines <- unique(rbindlist(lapply(frames, function(f) f[, .(ModelID, OncotreeLineage, exposed)])))
stopifnot(!anyDuplicated(all_lines$ModelID))
all_t <- function(exposure_by_line) {
  unlist(lapply(datasets, function(name) masked_t(prepared[[name]], exposure_by_line[frames[[name]]$ModelID])), use.names = FALSE)
}

split_classes <- function(value, pattern) {
  parts <- strsplit(value, pattern)
  lapply(parts, function(p) unique(trimws(p[!is.na(p) & trimws(p) != ""])))
}
class_membership <- function(values, pattern, type) {
  keyed <- per_dataset[match(feature_keys, tested_keys)]
  compound_key <- fifelse(is.na(keyed$compound_id), keyed$compound_name, keyed$compound_id)
  classes <- split_classes(values[match(feature_keys, tested_keys)], pattern)
  long <- data.table(index = rep(seq_along(classes), lengths(classes)), class = unlist(classes))
  long[, compound := compound_key[index]]
  keep <- long[, .(compounds = uniqueN(compound), features = .N), by = class][compounds >= dp$moa_min_compounds]
  long <- long[class %in% keep$class]
  membership <- Matrix::sparseMatrix(i = long$index, j = match(long$class, keep$class), x = 1, dims = c(length(feature_keys), nrow(keep)))
  list(type = type, classes = keep, membership = membership)
}
class_sets <- list(
  class_membership(per_dataset$mechanism, ",", "mechanism"),
  class_membership(per_dataset$targets, "[;,]", "target")
)
class_means <- function(t_stat, set) {
  t_stat[!is.finite(t_stat)] <- 0
  as.numeric(Matrix::crossprod(set$membership, t_stat)) / set$classes$features
}
observed_exposure <- setNames(all_lines$exposed, all_lines$ModelID)
observed_t <- all_t(observed_exposure)
observed <- lapply(class_sets, function(set) class_means(observed_t, set))
set.seed(params$seed + 2)
exceed <- lapply(class_sets, function(set) numeric(nrow(set$classes)))
for (b in seq_len(dp$moa_permutations)) {
  shuffled <- setNames(permute_within(all_lines$exposed, all_lines$OncotreeLineage), all_lines$ModelID)
  shuffled_t <- all_t(shuffled)
  for (k in seq_along(class_sets)) {
    exceed[[k]] <- exceed[[k]] + (abs(class_means(shuffled_t, class_sets[[k]])) >= abs(observed[[k]]))
  }
}
mechanism <- rbindlist(lapply(seq_along(class_sets), function(k) {
  out <- copy(class_sets[[k]]$classes)
  out[, class_type := class_sets[[k]]$type]
  out[, mean_t := observed[[k]]]
  out[, p := (exceed[[k]] + 1) / (dp$moa_permutations + 1)]
  out[, fdr := bh(p)]
  out[order(p)]
}))
write_result(mechanism, "drug_pan_mechanism")
print(mechanism[, .(classes = .N, fdr_hits = sum(fdr < dp$fdr), min_p = min(p)), by = class_type])
print(mechanism[1:12, .(class_type, class = substr(class, 1, 40), compounds, features, mean_t, p, fdr)])

jsonlite::write_json(list(
  per_dataset = per_dataset[, .(tested = .N, fdr_hits = sum(associated), loy_specific = sum(label == "LOY-specific"), min_p = min(p)), by = dataset],
  consensus = list(compounds = nrow(consensus), fdr_hits = sum(consensus$fdr < dp$fdr, na.rm = TRUE), min_p = min(consensus$p, na.rm = TRUE)),
  classes = mechanism[, .(classes = .N, fdr_hits = sum(fdr < dp$fdr), min_p = min(p)), by = class_type],
  permutations = dp$moa_permutations
), file.path(derived_dir, "drug_pan_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("drug scan done")
