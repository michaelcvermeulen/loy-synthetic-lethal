# Stage 17 (plan Stage 7b): forests that predict each Y-loss dependency from expression.
# One dependency per array task (the task ID indexes expression_dependency_set);
# "collect" binds the per-dependency outputs.
# 1. Nested models compared by out-of-fold accuracy, donors kept within folds:
#      m0 lineage + library | m1 + Y status | m2 + the seven Y genes | m3 + genome-wide expression
#    The m3 pre-filter (top features by correlation) is refitted inside every training fold.
# 2. Boruta-style selection with stratified shadows, run twice: shadows shuffled
#    within lineage ("given lineage"), and within lineage x Y status with Y status
#    offered to the forest ("beyond Y status").
# 3. Partial dependence of the dependency on the top stable features.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "forest_helpers.R"))
assert_freeze()

out_dir <- file.path(derived_dir, "forest_expression")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dependency_set <- read_result("expression_dependency_set")
arg <- commandArgs(trailingOnly = TRUE)[1]

if (identical(arg, "collect")) {
  bind <- function(prefix) rbindlist(lapply(list.files(out_dir, pattern = paste0("^", prefix, "_[0-9]+\\.rds$"), full.names = TRUE), readRDS), fill = TRUE)
  accuracy <- bind("accuracy")
  features <- bind("features")
  partial <- bind("partial")
  nulls <- bind("null")
  stopifnot(uniqueN(accuracy$dependency) == nrow(dependency_set))
  write_result(accuracy, "forest_expression_accuracy")
  write_result(features, "forest_expression_features")
  write_result(partial, "forest_expression_partial_dependence")
  write_result(nulls, "forest_expression_null")
  wide <- dcast(accuracy, dependency ~ model, value.var = "r2_mean")
  print(wide)
  print(features[stable == TRUE, .(stable_features = .N, top = paste(head(gene_symbol[order(-mean_importance)], 6), collapse = ", ")), by = .(dependency, strata)])
  print(nulls)
  quit(status = 0)
}

k <- as.integer(if (identical(arg, "array")) Sys.getenv("SLURM_ARRAY_TASK_ID") else arg)
stopifnot(!is.na(k), k >= 1, k <= nrow(dependency_set))
dependency <- dependency_set$symbol[k]
core <- params$y_genes_core
cohort <- read_result("cohort")
effect <- readRDS(file.path(matrix_dir, "crispr_effect.rds"))
expr <- readRDS(file.path(matrix_dir, "expr_tpm_pc.rds"))
lines <- cohort[ModelID %in% rownames(effect) & ModelID %in% rownames(expr)]

if (dependency == "COMPOSITE") {
  composite_genes <- read_result("gene_pan_hits")[estimate < 0, feature]
  d <- rowMeans(scale(effect[lines$ModelID, composite_genes, drop = FALSE]), na.rm = TRUE)
} else {
  d <- effect[lines$ModelID, dependency_set$feature[k]]
}
# The expression universe is fixed on all lines before any line is dropped for this dependency.
expr_lines <- expr[lines$ModelID, , drop = FALSE]
expressed <- colMeans(expr_lines) >= 1 & apply(expr_lines, 2, sd) >= 0.25
keep <- is.finite(d)
lines <- lines[keep]
d <- d[keep]
expr_lines <- expr_lines[keep, expressed, drop = FALSE]
symbols <- gene_symbol(colnames(expr_lines))
gene_key <- data.table(code = paste0("g", seq_len(ncol(expr_lines))), gene = colnames(expr_lines), gene_symbol = symbols)
colnames(expr_lines) <- gene_key$code
core_codes <- gene_key[gene_symbol %in% core, code]
stopifnot(length(core_codes) == length(core))
log_msg("%s: %d lines (%d Y-negative, %d Y-positive, %d set aside), %d expressed genes", dependency, nrow(lines),
  sum(lines$y_negative %in% 1), sum(lines$y_negative %in% 0), sum(is.na(lines$y_negative)), ncol(expr_lines))

covariates <- cbind(lineage_columns(lines$OncotreeLineage), LibraryKY = lines$LibraryKY, LibraryHumagneCD = lines$LibraryHumagneCD)
y_status <- cbind(y_negative = as.numeric(lines$y_negative %in% 1), y_set_aside = as.numeric(is.na(lines$y_negative)))
y_genes <- expr_lines[, core_codes, drop = FALSE]
genome <- expr_lines[, setdiff(colnames(expr_lines), core_codes), drop = FALSE]

# --- 1. nested models, out of fold ---------------------------------------------
fit_predict <- function(train_x, train_y, test_x, seed) {
  fit <- ranger(x = train_x, y = train_y, num.trees = forest_par$trees, mtry = max(1, floor(ncol(train_x) / 3)), num.threads = n_cores, seed = seed)
  predict(fit, data = test_x, num.threads = n_cores)$predictions
}
models <- c("m0 lineage+library", "m1 + Y status", "m2 + seven Y genes", "m3 + all expression")
accuracy <- rbindlist(lapply(seq_len(forest_par$cv_repeats), function(repeat_id) {
  set.seed(params$seed + 10 * repeat_id + k)
  folds <- donor_folds(lines$PatientID)
  predictions <- matrix(NA_real_, nrow(lines), length(models), dimnames = list(NULL, models))
  for (fold in sort(unique(folds))) {
    test <- folds == fold
    train <- !test
    # Pre-filter on the training fold only.
    correlation <- abs(cor(genome[train, , drop = FALSE], d[train]))[, 1]
    selected <- names(sort(correlation, decreasing = TRUE))[seq_len(min(forest_par$prefilter_features, length(correlation)))]
    designs <- list(covariates, cbind(covariates, y_status), cbind(covariates, y_status, y_genes), cbind(covariates, y_status, y_genes, genome[, selected, drop = FALSE]))
    for (m in seq_along(models)) {
      predictions[test, m] <- fit_predict(designs[[m]][train, , drop = FALSE], d[train], designs[[m]][test, , drop = FALSE], seed = repeat_id * 1000 + fold * 10 + m)
    }
  }
  data.table(
    dependency = dependency, model = models, repeat_id = repeat_id,
    r = apply(predictions, 2, cor, y = d), r2 = 1 - colSums((predictions - d)^2) / sum((d - mean(d))^2)
  )
}))
accuracy_summary <- accuracy[, .(lines = length(d), r_mean = mean(r), r_sd = sd(r), r2_mean = mean(r2), r2_sd = sd(r2)), by = .(dependency, model)]
print(accuracy_summary)
saveRDS(accuracy_summary, file.path(out_dir, sprintf("accuracy_%d.rds", k)))

# --- 2. Boruta-style selection -------------------------------------------------
variable <- names(sort(apply(genome, 2, sd), decreasing = TRUE))[seq_len(min(forest_par$boruta_expression_features, ncol(genome)))]
partner_codes <- gene_key[gene_symbol %in% params$xy_paralogs$x_partner, code]
candidates <- unique(c(core_codes, intersect(partner_codes, colnames(expr_lines)), variable))
x <- expr_lines[, candidates, drop = FALSE]
status_level <- fifelse(is.na(lines$y_negative), "set_aside", fifelse(lines$y_negative == 1, "negative", "positive"))
runs <- list(
  list(strata = "given lineage", by = lines$OncotreeLineage, always = covariates),
  list(strata = "beyond Y status", by = paste(lines$OncotreeLineage, status_level), always = cbind(covariates, y_status))
)
selected <- lapply(runs, function(run) boruta_two_seeds(x, d, run$by, run$always))
features <- rbindlist(lapply(seq_along(runs), function(i) {
  out <- selected[[i]]$features
  out <- merge(out, gene_key, by.x = "feature", by.y = "code", sort = FALSE)
  out[, c("dependency", "strata", "is_core_y_gene") := list(dependency, runs[[i]]$strata, gene_symbol %in% core)]
  out[order(-hits_seed1, -hits_seed2)]
}))
nulls <- rbindlist(lapply(seq_along(runs), function(i) data.table(
  dependency = dependency, strata = runs[[i]]$strata, candidates = ncol(x),
  stable = sum(selected[[i]]$features$stable), null_confirmed = paste(selected[[i]]$null_confirmed, collapse = ";"),
  mean_null_confirmed = mean(selected[[i]]$null_confirmed)
)))
print(nulls)
print(features[stable == TRUE, .(strata, gene_symbol, hits_seed1, hits_seed2, mean_importance)][1:min(.N, 30)])
saveRDS(features[, .(dependency, strata, gene, gene_symbol, is_core_y_gene, hits_seed1, hits_seed2, mean_importance, confirmed_seed1, confirmed_seed2, stable)], file.path(out_dir, sprintf("features_%d.rds", k)))
saveRDS(nulls, file.path(out_dir, sprintf("null_%d.rds", k)))

# --- 3. partial dependence on the top stable features ---------------------------
top <- features[strata == "given lineage" & stable == TRUE][order(-mean_importance)][seq_len(min(.N, 6)), feature]
paired <- gene_key[gene_symbol == dependency_set$y_partner[k], code]
top <- unique(c(top, intersect(paired, colnames(x))))
partial <- NULL
if (length(top) > 0) {
  design <- cbind(covariates, x[, unique(c(top, core_codes)), drop = FALSE])
  fit <- ranger(x = design, y = d, num.trees = 1000, num.threads = n_cores, seed = forest_par$seeds[1])
  partial <- rbindlist(lapply(top, function(top_code) {
    grid <- unique(quantile(design[, top_code], seq(0.02, 0.98, length.out = 25)))
    rbindlist(lapply(grid, function(value) {
      shifted <- design
      shifted[, top_code] <- value
      data.table(dependency = dependency, gene_symbol = gene_key$gene_symbol[match(top_code, gene_key$code)], expression = value,
        predicted = mean(predict(fit, data = shifted, num.threads = n_cores)$predictions))
    }))
  }))
}
saveRDS(partial, file.path(out_dir, sprintf("partial_%d.rds", k)))
log_msg("%s done", dependency)
