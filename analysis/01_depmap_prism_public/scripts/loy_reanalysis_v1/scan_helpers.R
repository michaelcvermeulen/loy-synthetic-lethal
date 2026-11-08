# Scan helpers shared by the regression stages. Source after common.R.
# Needs limma, so these run in the scRNA environment only.

suppressPackageStartupMessages(library(limma))

load_matrix <- function(name) readRDS(file.path(matrix_dir, paste0(name, ".rds")))

# Cohort rows that have the exposure defined and are present in the matrix.
# `exposed` is 1 for Y-negative and 0 for Y-positive.
analysis_frame <- function(cohort, mat, exposure = "y_call") {
  frame <- cohort[ModelID %in% rownames(mat) & !is.na(get(exposure))]
  frame[, exposed := as.numeric(get(exposure) == "Y_negative")]
  frame
}

# Design matrix with the tested column(s) first, so that aliasing never drops
# them, then lineage dummies and numeric covariates. Rows with a missing
# covariate are removed; `rows` gives the rows of `frame` that were kept.
build_design <- function(frame, tested, covariates = character(), lineage = TRUE) {
  columns <- c(tested, covariates)
  complete <- stats::complete.cases(frame[, columns, with = FALSE])
  frame <- frame[complete]
  parts <- list(intercept = rep(1, nrow(frame)))
  for (col in tested) parts[[col]] <- as.numeric(frame[[col]])
  design <- do.call(cbind, parts)
  if (lineage) {
    lineage_part <- lineage_design(frame$OncotreeLineage)[, -1, drop = FALSE]
    design <- cbind(design, lineage_part)
  }
  for (col in covariates) {
    design <- cbind(design, as.numeric(frame[[col]]))
    colnames(design)[ncol(design)] <- col
  }
  # A tested column is estimable only if the other columns do not reproduce it;
  # full_rank alone would keep it and silently drop the nuisance column instead.
  estimable <- all(vapply(tested, function(col) {
    others <- design[, colnames(design) != col, drop = FALSE]
    sum(qr.resid(qr(others), design[, col])^2) > 1e-8 * nrow(design)
  }, logical(1)))
  design <- full_rank(design)
  list(design = design, rows = which(complete), estimable = estimable && all(tested %in% colnames(design)))
}

# Join two limma fits of the same design row-wise (features of the first, then the second).
combine_fits <- function(fit_a, fit_b) {
  fit <- fit_a
  for (component in c("coefficients", "stdev.unscaled")) fit[[component]] <- rbind(fit_a[[component]], fit_b[[component]])
  for (component in c("sigma", "df.residual", "Amean")) fit[[component]] <- c(fit_a[[component]], fit_b[[component]])
  fit
}

# limma scan of one coefficient over features (rows of `mat_t`), samples in
# columns. With repeated donors, a consensus within-donor correlation is
# estimated on a random subset of complete features and used as a blocking
# correlation. Features with missing values are fitted separately, because one
# missing value otherwise sends every feature down limma's slow per-row path.
limma_scan <- function(mat_t, design, coef, block = NULL, correlation = NULL, trend = TRUE) {
  stopifnot(ncol(mat_t) == nrow(design), coef %in% colnames(design))
  complete <- rowSums(is.na(mat_t)) == 0
  use_block <- !is.null(block) && anyDuplicated(block) > 0
  if (use_block && is.null(correlation)) {
    set.seed(params$seed)
    pool <- if (sum(complete) >= 50) which(complete) else seq_len(nrow(mat_t))
    subset <- pool[sample.int(length(pool), min(length(pool), params$gene_scan$dupcor_features))]
    correlation <- duplicateCorrelation(mat_t[subset, , drop = FALSE], design, block = block)$consensus.correlation
    correlation <- min(max(correlation, 0), 0.99)
  }
  fit_part <- function(part) {
    if (use_block && correlation > 0) {
      lmFit(part, design, block = block, correlation = correlation)
    } else {
      lmFit(part, design)
    }
  }
  if (all(complete) || !any(complete)) {
    fit <- fit_part(mat_t)
  } else {
    fit <- combine_fits(fit_part(mat_t[complete, , drop = FALSE]), fit_part(mat_t[!complete, , drop = FALSE]))
    restore <- order(c(which(complete), which(!complete)))
    for (component in c("coefficients", "stdev.unscaled")) fit[[component]] <- fit[[component]][restore, , drop = FALSE]
    for (component in c("sigma", "df.residual", "Amean")) fit[[component]] <- fit[[component]][restore]
  }
  fit <- eBayes(fit, trend = trend, robust = TRUE)
  tab <- topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = TRUE)
  out <- data.table(
    feature = rownames(mat_t), estimate = tab$logFC, ci_low = tab$CI.L, ci_high = tab$CI.R,
    t = tab$t, p = tab$P.Value, df = fit$df.total, n = rowSums(!is.na(mat_t))
  )
  attr(out, "correlation") <- if (use_block) correlation else NA_real_
  out
}

# Plain least-squares t statistics for one tested column, vectorised over
# features with no missing values. Used for permutation nulls, where the same
# statistic is computed for the observed and the shuffled labels.
fast_t <- function(mat_t, tested, nuisance) {
  decomposition <- qr(nuisance)
  x_res <- qr.resid(decomposition, tested)
  sxx <- sum(x_res^2)
  y_res <- t(qr.resid(decomposition, t(mat_t)))
  estimate <- as.numeric(y_res %*% x_res) / sxx
  rss <- rowSums(y_res^2) - estimate^2 * sxx
  df <- ncol(mat_t) - decomposition$rank - 1
  estimate / sqrt(rss / df / sxx)
}

# Group means and sizes of each feature by a 0/1 exposure.
group_summary <- function(mat_t, exposed) {
  data.table(
    feature = rownames(mat_t),
    mean_negative = rowMeans(mat_t[, exposed == 1, drop = FALSE], na.rm = TRUE),
    mean_positive = rowMeans(mat_t[, exposed == 0, drop = FALSE], na.rm = TRUE),
    n_negative = rowSums(!is.na(mat_t[, exposed == 1, drop = FALSE])),
    n_positive = rowSums(!is.na(mat_t[, exposed == 0, drop = FALSE]))
  )
}

# Pre-ranked gene-set enrichment for each collection in the reference snapshot.
run_fgsea <- function(stat, collections = c("hallmark", "reactome", "gobp"), min_size = 15, max_size = 500) {
  reference <- readRDS(file.path(derived_dir, "reference", "gene_sets.rds"))$sets
  stat <- stat[is.finite(stat) & !duplicated(names(stat))]
  rbindlist(lapply(collections, function(collection) {
    set.seed(params$seed)
    res <- fgsea::fgsea(reference[[collection]], stat, minSize = min_size, maxSize = max_size, nproc = 1)
    res <- as.data.table(res)
    res[, leading_edge := vapply(leadingEdge, function(g) paste(head(g, 15), collapse = ";"), character(1))]
    res[, leadingEdge := NULL]
    res[, collection := collection]
    res[, fdr := bh(pval)]
    res[order(pval)]
  }))
}
