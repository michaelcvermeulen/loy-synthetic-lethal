# Stage 16 (plan Stage 7a): genome-wide Boruta-style selection of the
# dependencies, and of the drug responses, that distinguish Y-negative from
# Y-positive lines. Shadows are shuffled within lineage, and lineage and library
# columns are offered to the forest without shadows, so a feature has to carry
# information beyond lineage to be confirmed.
# Forest importance has no FDR: results are reported with two-seed stability and
# with the number of features confirmed when Y status is shuffled within lineage.
# Usage: 16_forest_y_status.R <dataset>   (crispr_effect or a drug dataset name)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "forest_helpers.R"))
assert_freeze()

name <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(name))
out_dir <- file.path(derived_dir, "forest_y_status")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cohort <- read_result("cohort")
mat <- readRDS(file.path(matrix_dir, paste0(name, ".rds")))
frame <- cohort[in_primary == TRUE & ModelID %in% rownames(mat)]
y <- factor(frame$y_negative, levels = c(0, 1))

if (name == "crispr_effect") {
  gene_ann <- read_result("gene_annotation")
  features <- gene_ann[selective == TRUE, feature]
  x <- mat[frame$ModelID, features, drop = FALSE]
  # ranger takes no missing values; only a handful of genes have any.
  x <- x[, colSums(is.na(x)) == 0, drop = FALSE]
  always <- cbind(lineage_columns(frame$OncotreeLineage), LibraryKY = frame$LibraryKY, LibraryHumagneCD = frame$LibraryHumagneCD)
  reference <- read_result("gene_pan_results")[, .(feature, symbol, regression_estimate = estimate, regression_fdr = fdr, regression_label = label)]
} else {
  pan <- read_result("drug_pan_results")
  features <- pan[dataset == name, feature]
  x <- mat[frame$ModelID, features, drop = FALSE]
  complete_enough <- colMeans(!is.na(x)) >= params$drug_scan$forest_min_completeness
  x <- impute_within(x[, complete_enough, drop = FALSE], frame$OncotreeLineage)
  always <- lineage_columns(frame$OncotreeLineage)
  reference <- pan[dataset == name, .(feature, symbol = compound_name, regression_estimate = standardized, regression_fdr = fdr, regression_label = label)]
}
colnames_original <- colnames(x)
colnames(x) <- paste0("f", seq_len(ncol(x)))
log_msg("%s: %d lines (%d Y-negative, %d Y-positive), %d features", name, nrow(frame), sum(frame$y_negative == 1), sum(frame$y_negative == 0), ncol(x))

# How well the features predict Y status at all, out of bag, against a shuffled-label forest.
plain <- ranger(x = cbind(x, always), y = y, num.trees = 1000, probability = TRUE, num.threads = n_cores, seed = forest_par$seeds[1])
auc <- rank_auc(plain$predictions[, "1"], frame$y_negative)
covariates_only <- ranger(x = always, y = y, num.trees = 1000, probability = TRUE, num.threads = n_cores, seed = forest_par$seeds[1])
auc_covariates <- rank_auc(covariates_only$predictions[, "1"], frame$y_negative)
log_msg("%s: out-of-bag AUC %.3f with features, %.3f from lineage and library alone", name, auc, auc_covariates)

selection <- boruta_two_seeds(x, y, frame$OncotreeLineage, always)
result <- selection$features
result[, feature := colnames_original[match(feature, colnames(x))]]
result <- merge(result, reference, by = "feature", all.x = TRUE, sort = FALSE)
result[, dataset := name]
setorder(result, -hits_seed1, -hits_seed2)
write_result(result, paste0("forest_y_status_", name), dir = out_dir)
summary_row <- data.table(
  dataset = name, lines = nrow(frame), features = ncol(x), oob_auc = auc, oob_auc_lineage_only = auc_covariates,
  confirmed_seed1 = sum(result$confirmed_seed1), confirmed_seed2 = sum(result$confirmed_seed2), stable = sum(result$stable),
  null_confirmed = paste(selection$null_confirmed, collapse = ";"), mean_null_confirmed = mean(selection$null_confirmed),
  stable_and_regression_fdr = sum(result$stable & result$regression_fdr < 0.05, na.rm = TRUE)
)
write_result(summary_row, paste0("forest_y_status_summary_", name), dir = out_dir)
print(summary_row)
print(result[confirmed_seed1 | confirmed_seed2][1:min(.N, 25), .(symbol, hits_seed1, hits_seed2, stable, regression_estimate, regression_fdr)])
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info_main.txt"))
