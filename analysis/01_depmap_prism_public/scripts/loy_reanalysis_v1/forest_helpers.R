# Forest helpers for the Stage 7 scripts. Source after common.R.
# Runs in the main environment (R 4.2.2, ranger 0.14.1); no limma here.

suppressPackageStartupMessages(library(ranger))

forest_par <- params$forest

# Shadow copy of every column, shuffled within strata. A feature that only
# marks a stratum (a lineage) keeps that information in its shadow, so it
# cannot beat the shadow on stratum information alone.
make_shadows <- function(x, strata) {
  groups <- split(seq_len(nrow(x)), strata)
  groups <- groups[lengths(groups) > 1]
  shadows <- x
  for (idx in groups) {
    for (j in seq_len(ncol(x))) shadows[idx, j] <- x[idx[sample.int(length(idx))], j]
  }
  colnames(shadows) <- paste0("shadow_", colnames(x))
  shadows
}

# Boruta-style selection with stratified shadows. Every iteration fits one
# forest on the candidate features, their shadows and the `always` columns
# (covariates that get no shadow), and counts a hit for each feature whose
# permutation importance beats the best shadow. A feature is confirmed when its
# hit count is significantly above half the iterations (one-sided binomial,
# Bonferroni over features, 0.01), which is the Boruta decision rule. Unlike
# Boruta, no feature is dropped along the way, so the bar never gets lower.
boruta_stratified <- function(x, y, strata, always = NULL, iterations = forest_par$iterations,
                              trees = forest_par$trees, seed = 1, threads = n_cores) {
  stopifnot(!anyNA(x), !anyNA(y), nrow(x) == length(y), !is.null(colnames(x)))
  set.seed(seed)
  hits <- setNames(numeric(ncol(x)), colnames(x))
  importance_sum <- hits
  for (iteration in seq_len(iterations)) {
    shadows <- make_shadows(x, strata)
    data <- cbind(x, shadows)
    if (!is.null(always)) data <- cbind(data, always)
    fit <- ranger(x = data, y = y, num.trees = trees, importance = "permutation", num.threads = threads, seed = seed * 100000 + iteration)
    importance <- fit$variable.importance
    hits <- hits + (importance[colnames(x)] > max(importance[colnames(shadows)]))
    importance_sum <- importance_sum + importance[colnames(x)]
  }
  p <- pbinom(hits - 1, iterations, 0.5, lower.tail = FALSE)
  data.table(
    feature = colnames(x), hits = hits, iterations = iterations, mean_importance = importance_sum / iterations,
    p = p, confirmed = p.adjust(p, "bonferroni") < 0.01
  )
}

# Two seeds on the real outcome plus null runs with the outcome shuffled within strata.
boruta_two_seeds <- function(x, y, strata, always = NULL, null_runs = forest_par$null_runs, ...) {
  seeds <- forest_par$seeds
  real <- lapply(seeds, function(seed) boruta_stratified(x, y, strata, always, seed = seed, ...))
  out <- data.table(
    feature = real[[1]]$feature, hits_seed1 = real[[1]]$hits, hits_seed2 = real[[2]]$hits,
    mean_importance = (real[[1]]$mean_importance + real[[2]]$mean_importance) / 2,
    confirmed_seed1 = real[[1]]$confirmed, confirmed_seed2 = real[[2]]$confirmed
  )
  out[, stable := confirmed_seed1 & confirmed_seed2]
  null_confirmed <- vapply(seq_len(null_runs), function(run) {
    set.seed(9000 + run)
    shuffled <- y
    for (idx in split(seq_along(y), strata)) if (length(idx) > 1) shuffled[idx] <- y[idx[sample.int(length(idx))]]
    sum(boruta_stratified(x, shuffled, strata, always, seed = 7000 + run, ...)$confirmed)
  }, numeric(1))
  list(features = out, null_confirmed = null_confirmed)
}

# One-hot lineage columns for use as plain forest features.
lineage_columns <- function(lineage) {
  lineage <- factor(as.character(lineage))
  out <- model.matrix(~ 0 + lineage)
  colnames(out) <- paste0("lineage_", gsub("[^A-Za-z0-9]+", "_", levels(lineage)))
  out
}

# Median imputation within strata, falling back to the overall median.
impute_within <- function(x, strata) {
  for (j in seq_len(ncol(x))) {
    missing <- is.na(x[, j])
    if (!any(missing)) next
    by_stratum <- tapply(x[, j], strata, median, na.rm = TRUE)
    fill <- by_stratum[as.character(strata[missing])]
    fill[is.na(fill)] <- median(x[, j], na.rm = TRUE)
    x[missing, j] <- fill
  }
  x
}

# Fold assignment that never splits a donor.
donor_folds <- function(donor, k = forest_par$cv_folds) {
  donors <- unique(donor)
  assignment <- setNames(sample(rep_len(seq_len(k), length(donors))), donors)
  unname(assignment[as.character(donor)])
}

rank_auc <- function(score, label) {
  ranks <- rank(score)
  n1 <- sum(label == 1)
  n0 <- sum(label == 0)
  (sum(ranks[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
