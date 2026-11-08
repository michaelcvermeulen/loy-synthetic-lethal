# Stage 08: prove the scan code before it is pointed at Y status.
# - Orientation and power: known genotype-dependency and genotype-drug pairs
#   are run through the same functions with the genotype as the exposure.
# - Calibration: Y status shuffled within lineage must give flat p-values.
# A wrong sign or a miscalibrated null stops the pipeline.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
assert_freeze()

cohort <- read_result("cohort")
context <- read_result("context_matrix")
gene_ann <- read_result("gene_annotation")
drug_ann <- read_result("drug_feature_annotation")
drug_datasets <- unique(drug_ann$dataset)
min_group <- params$drug_scan$min_per_group
gene_covariates <- c("LibraryKY", "LibraryHumagneCD")

exposures <- list(
  BRAF = context$driver__BRAF[match(cohort$ModelID, context$ModelID)],
  KRAS = context$driver__KRAS[match(cohort$ModelID, context$ModelID)],
  TP53 = cohort$tp53_altered
)
stopifnot(!is.null(exposures$BRAF), !is.null(exposures$KRAS))

# --- gene controls -----------------------------------------------------------
effect <- load_matrix("crispr_effect")
scan_genes <- gene_ann[in_scan == TRUE, feature]
gene_control <- function(exposure_name, target_symbol, expected_sign) {
  frame <- copy(cohort)[, exposed := exposures[[exposure_name]]]
  frame <- frame[ModelID %in% rownames(effect) & !is.na(exposed)]
  design <- build_design(frame, tested = "exposed", covariates = gene_covariates)
  frame <- frame[design$rows]
  mat_t <- t(effect[frame$ModelID, scan_genes])
  scan <- limma_scan(mat_t, design$design, "exposed", block = frame$PatientID)
  scan[, rank := frank(p)]
  row <- scan[gene_symbol(feature) == target_symbol]
  data.table(
    control = sprintf("%s altered -> %s dependency", exposure_name, target_symbol), dataset = "crispr_effect",
    feature = row$feature, n_exposed = sum(frame$exposed == 1), n_unexposed = sum(frame$exposed == 0),
    estimate = row$estimate, p = row$p, rank = row$rank, expected_sign = expected_sign,
    pass = sign(row$estimate) == expected_sign && row$p < 1e-4
  )
}
gene_controls <- rbindlist(list(
  gene_control("BRAF", "BRAF", -1),
  gene_control("KRAS", "KRAS", -1),
  gene_control("TP53", "MDM2", 1)
))
print(gene_controls)

# --- drug controls -----------------------------------------------------------
drug_control <- function(exposure_name, pattern, expected_sign, label) {
  rbindlist(lapply(drug_datasets, function(name) {
    features <- drug_ann[dataset == name & (grepl(pattern, compound_name, ignore.case = TRUE) | grepl(pattern, source_name, ignore.case = TRUE)), feature]
    mat <- load_matrix(name)
    features <- intersect(features, colnames(mat))
    if (length(features) == 0) return(NULL)
    frame <- copy(cohort)[, exposed := exposures[[exposure_name]]]
    frame <- frame[ModelID %in% rownames(mat) & !is.na(exposed)]
    design <- build_design(frame, tested = "exposed")
    nuisance <- design$design[, colnames(design$design) != "exposed", drop = FALSE]
    scan <- scan_robust(mat[frame$ModelID, features, drop = FALSE], frame$exposed, nuisance, cluster = frame$PatientID, cores = 1)
    scan <- scan[n_exposed >= min_group & n_unexposed >= min_group]
    if (nrow(scan) == 0) return(NULL)
    best <- scan[which.min(p)]
    data.table(
      control = label, dataset = name, feature = best$feature, n_exposed = best$n_exposed, n_unexposed = best$n_unexposed,
      estimate = best$estimate, p = best$p, rank = NA_real_, expected_sign = expected_sign,
      pass = sign(best$estimate) == expected_sign && best$p < 0.01
    )
  }))
}
drug_controls <- rbind(
  drug_control("BRAF", "DABRAFENIB|VEMURAFENIB|PLX-?4720|PLX-?4032|ENCORAFENIB", -1, "BRAF altered -> BRAF inhibitor sensitivity"),
  drug_control("TP53", "NUTLIN|IDASANUTLIN|RG-?7112|RG-?7388|AMG-?232|SERDEMETAN|MILADEMETAN|HDM-?201|CGM-?097|SAR-?405838|MI-773", 1, "TP53 altered -> MDM2 inhibitor resistance")
)
print(drug_controls)
controls <- rbind(gene_controls, drug_controls)
write_result(controls, "pipeline_controls")

# --- calibration under shuffled Y status -------------------------------------
frame <- analysis_frame(cohort, effect)
design <- build_design(frame, tested = "exposed", covariates = gene_covariates)
frame <- frame[design$rows]
mat_t <- t(effect[frame$ModelID, scan_genes])
mat_t <- mat_t[rowSums(is.na(mat_t)) == 0, , drop = FALSE]
nuisance <- design$design[, colnames(design$design) != "exposed", drop = FALSE]
df <- ncol(mat_t) - qr(nuisance)$rank - 1
set.seed(params$seed)
calibration <- rbindlist(lapply(seq_len(params$calibration_permutations), function(b) {
  shuffled <- permute_within(frame$exposed, frame$OncotreeLineage)
  p <- 2 * pt(-abs(fast_t(mat_t, shuffled, nuisance)), df)
  data.table(
    permutation = b, fraction_below_0.05 = mean(p < 0.05), fdr_hits = sum(bh(p) < 0.05),
    lambda = median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
  )
}))
write_result(calibration, "calibration_gene_permutations")
gene_calibration <- calibration[, .(
  outcome = "crispr_effect", genes = nrow(mat_t), permutations = .N, mean_fraction_below_0.05 = mean(fraction_below_0.05),
  median_lambda = median(lambda), permutations_with_fdr_hits = sum(fdr_hits > 0), max_fdr_hits = max(fdr_hits)
)]
print(gene_calibration)

# The robust drug test, checked the same way on the densest PRISM screen.
prism <- load_matrix("prism_REP1M")
drug_frame <- analysis_frame(cohort, prism)
drug_design <- build_design(drug_frame, tested = "exposed")
drug_nuisance <- drug_design$design[, colnames(drug_design$design) != "exposed", drop = FALSE]
active <- intersect(drug_ann[dataset == "prism_REP1M" & active == TRUE, feature], colnames(prism))
drug_mat <- prism[drug_frame$ModelID, active]
set.seed(params$seed + 1)
drug_calibration <- rbindlist(lapply(seq_len(50), function(b) {
  shuffled <- permute_within(drug_frame$exposed, drug_frame$OncotreeLineage)
  scan <- scan_robust(drug_mat, shuffled, drug_nuisance, cluster = drug_frame$PatientID)
  scan <- scan[n_exposed >= min_group & n_unexposed >= min_group]
  data.table(permutation = b, fraction_below_0.05 = mean(scan$p < 0.05), fdr_hits = sum(bh(scan$p) < 0.05))
}))
write_result(drug_calibration, "calibration_drug_permutations")
drug_summary <- drug_calibration[, .(
  outcome = "prism_REP1M", features = length(active), permutations = .N, mean_fraction_below_0.05 = mean(fraction_below_0.05),
  permutations_with_fdr_hits = sum(fdr_hits > 0), max_fdr_hits = max(fdr_hits)
)]
print(drug_summary)

checks <- list(
  gene_controls_pass = all(gene_controls$pass),
  braf_drug_control_pass = any(drug_controls[grepl("BRAF", control), pass]) && !any(drug_controls[grepl("BRAF", control), sign(estimate) != expected_sign & p < 0.01]),
  tp53_drug_control_pass = any(drug_controls[grepl("TP53", control), pass]) && !any(drug_controls[grepl("TP53", control), sign(estimate) != expected_sign & p < 0.01]),
  gene_calibration_pass = gene_calibration$mean_fraction_below_0.05 < 0.07 && gene_calibration$permutations_with_fdr_hits <= 0.1 * gene_calibration$permutations,
  drug_calibration_pass = drug_summary$mean_fraction_below_0.05 < 0.07 && drug_summary$permutations_with_fdr_hits <= 0.1 * drug_summary$permutations
)
jsonlite::write_json(c(checks, list(
  gene_calibration = as.list(gene_calibration), drug_calibration = as.list(drug_summary)
)), file.path(derived_dir, "pipeline_controls.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
print(unlist(checks))
if (!all(unlist(checks))) stop("pipeline controls failed: ", paste(names(checks)[!unlist(checks)], collapse = ", "))
log_msg("pipeline controls passed")
