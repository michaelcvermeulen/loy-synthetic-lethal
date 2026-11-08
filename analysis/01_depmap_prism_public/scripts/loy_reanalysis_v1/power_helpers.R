# Helpers for the power stage (24_power_context.R). Source after scan_helpers.R
# and context_helpers.R.
#
# The stage asks what the context screen could have detected. Simulated genes
# are fitted with the same calls as scan_subset() in 12_context_genes.R, next to
# the real genes, so every probability refers to the screen as it was run.

# --- significance bars -----------------------------------------------------------

# Largest p-value at which one more test added to `p_family` would be rejected
# by BH at `alpha`: the new test is rejected exactly when its p is below this.
# With k - 1 family p-values at or below k * alpha / (m + 1), the new one makes k.
bh_next_critical <- function(p_family, alpha = 0.05) {
  p <- sort(p_family[is.finite(p_family)])
  m <- length(p) + 1
  thresholds <- seq_len(m) * alpha / m
  max(thresholds[findInterval(thresholds, p) >= seq_len(m) - 1])
}

# Smallest p-value a distribution-free two-sample test can give (two-sided),
# reached when the two groups do not overlap.
exact_floor <- function(n_negative, n_positive) pmin(1, 2 / choose(n_negative + n_positive, n_negative))

# --- normal-theory limits ----------------------------------------------------------

# Power of a two-sided t test with critical p-value `p_crit` when the true
# effect is `shift` and its standard error is `se` (noncentral t).
analytic_power <- function(shift, se, df, p_crit) {
  ncp <- abs(shift) / se
  if (ncp == 0) return(p_crit)
  crit <- qt(1 - p_crit / 2, df)
  pt(crit, df, ncp, lower.tail = FALSE) + pt(-crit, df, ncp)
}

# Smallest effect detected with the given power. R's noncentral t stops at a
# noncentrality of 37.62; an effect beyond that is returned as NA.
analytic_mde <- function(se, df, p_crit, power = 0.8) {
  gap <- function(shift) analytic_power(shift, se, df, p_crit) - power
  upper <- 37.5 * se
  if (gap(upper) < 0) return(NA_real_)
  uniroot(gap, c(0, upper), tol = 1e-12)$root
}

# --- simulated genes ---------------------------------------------------------------

# Rows `genes` of `mat_raw` with `shift` added in the columns flagged by
# `target` (one shift for all genes, or one per gene).
spike_cell <- function(mat_raw, genes, target, shift) {
  stopifnot(length(target) == ncol(mat_raw), length(shift) %in% c(1L, length(genes)))
  out <- mat_raw[genes, , drop = FALSE]
  out[, target] <- out[, target, drop = FALSE] + shift
  out
}

# One known genotype-dependency pair, `effect ~ lineage + covariates + biomarker`
# by least squares on the lines of `frame` with the biomarker known. Keeps the
# coefficient, the lineage and covariate coefficients, and the residuals of
# altered and reference lines separately.
benchmark_fit <- function(frame, effect, biomarker, gene, covariates = character()) {
  stopifnot(length(biomarker) == nrow(frame), gene %in% colnames(effect))
  data <- copy(frame)[, exposed := as.numeric(biomarker)]
  data[, response := effect[match(ModelID, rownames(effect)), gene]]
  data <- data[!is.na(exposed) & is.finite(response)]
  design <- build_design(data, tested = "exposed", covariates = covariates)
  stopifnot(design$estimable)
  data <- data[design$rows]
  x <- design$design
  fit <- lm.fit(x, data$response)
  df <- nrow(x) - fit$rank
  sigma <- sqrt(sum(fit$residuals^2) / df)
  coefficient <- fit$coefficients[["exposed"]]
  se <- sigma * sqrt(solve(crossprod(x))["exposed", "exposed"])
  altered <- fit$residuals[data$exposed == 1]
  reference <- fit$residuals[data$exposed == 0]
  list(
    feature = gene, coefficient = coefficient, se = se, df = df, p = 2 * pt(-abs(coefficient / se), df),
    n_altered = length(altered), n_reference = length(reference),
    residual_sd = sigma, sd_altered = sd(altered), sd_reference = sd(reference),
    altered = unname(altered), reference = unname(reference),
    nuisance = fit$coefficients[names(fit$coefficients) != "exposed"], covariates = covariates
  )
}

# Fitted lineage and covariate effects of a benchmark for the lines of `frame`,
# without the biomarker term. A lineage the benchmark never saw takes the
# reference lineage's level.
benchmark_nuisance <- function(fit, frame) {
  lineage_term <- fit$nuisance[paste0("lineage", frame$OncotreeLineage)]
  out <- fit$nuisance[["intercept"]] + ifelse(is.na(lineage_term), 0, lineage_term)
  for (col in intersect(fit$covariates, names(fit$nuisance))) out <- out + fit$nuisance[[col]] * as.numeric(frame[[col]])
  unname(out)
}

# `reps` resampled copies of a benchmark gene over the lines of `frame`. Every
# line gets its fitted lineage and covariate effect plus a reference residual;
# lines in `target` get the benchmark coefficient as well and, with
# `with_spread`, an altered-group residual in place of the reference one.
benchmark_replicates <- function(fit, frame, target, reps, with_spread = TRUE) {
  stopifnot(length(target) == nrow(frame))
  draw <- function(pool, size) pool[sample.int(length(pool), size, replace = TRUE)]
  out <- matrix(draw(fit$reference, reps * nrow(frame)), reps, nrow(frame))
  if (with_spread) out[, target] <- draw(fit$altered, reps * sum(target))
  out[, target] <- out[, target] + fit$coefficient
  out <- sweep(out, 2, benchmark_nuisance(fit, frame), "+")
  colnames(out) <- frame$ModelID
  out
}

# --- the screen's fits ---------------------------------------------------------------

# The screen's two fits for one binary context on rank-normalised rows: the
# Y-by-context interaction on every line, and the Y effect on the lines inside
# the context. Same calls as scan_subset() in 12_context_genes.R.
fit_context_tests <- function(ctx_frame, mat, covariates) {
  one_fit <- function(sub_frame, tested, fit_covariates) {
    design <- build_design(sub_frame, tested = tested, covariates = fit_covariates)
    if (!design$estimable || nrow(design$design) - ncol(design$design) < 3) stop("context cannot be fitted: ", tested)
    sub_frame <- sub_frame[design$rows]
    scan <- limma_scan(mat[, sub_frame$ModelID, drop = FALSE], design$design, tested, trend = FALSE)
    scan[, se := abs(estimate / t)]
    scan
  }
  interaction <- one_fit(ctx_frame, "interaction", c("exposed", "ctx", covariates))
  inside <- one_fit(ctx_frame[ctx == 1], "exposed", covariates)
  data.table(
    feature = interaction$feature,
    interaction_p = interaction$p, interaction_estimate = interaction$estimate, interaction_se = interaction$se, interaction_df = interaction$df,
    inside_p = inside$p, inside_estimate = inside$estimate, inside_se = inside$se, inside_df = inside$df
  )
}

# The parts of the screen's two fits that stay the same when one gene is added:
# the real genes' residual variances and degrees of freedom, and the algebra
# that gives a new gene's estimate and residual variance. Real genes are fitted
# as limma_scan() fits them, complete rows together and rows with missing values
# apart.
context_fit_parts <- function(ctx_frame, mat, covariates) {
  one_fit <- function(sub_frame, tested, fit_covariates) {
    design <- build_design(sub_frame, tested = tested, covariates = fit_covariates)
    if (!design$estimable || nrow(design$design) - ncol(design$design) < 3) stop("context cannot be fitted: ", tested)
    sub_frame <- sub_frame[design$rows]
    x <- design$design
    y <- mat[, sub_frame$ModelID, drop = FALSE]
    complete <- rowSums(is.na(y)) == 0
    sigma <- rep(NA_real_, nrow(y))
    df_residual <- rep(NA_real_, nrow(y))
    for (part in list(complete, !complete)) {
      if (!any(part)) next
      fit <- lmFit(y[part, , drop = FALSE], x)
      sigma[part] <- fit$sigma
      df_residual[part] <- fit$df.residual
    }
    cross_inverse <- solve(crossprod(x))
    list(
      lines = sub_frame$ModelID, sigma2 = sigma^2, df_residual = df_residual,
      weights = drop(cross_inverse[tested, , drop = FALSE] %*% t(x)), unscaled = sqrt(cross_inverse[tested, tested]),
      basis = qr.Q(qr(x)), df_new = nrow(x) - ncol(x)
    )
  }
  list(
    interaction = one_fit(ctx_frame, "interaction", c("exposed", "ctx", covariates)),
    inside = one_fit(ctx_frame[ctx == 1], "exposed", covariates),
    lines_all = colnames(mat)
  )
}

# The p-values each simulated gene would get if it alone were added to the
# screen: limma's moderated t (eBayes with robust = TRUE, no trend) on the real
# genes plus that one gene. `simulated_raw` holds raw-scale values for every line
# of the screen and is rank-normalised over all lines first, as in the screen.
one_more_tests <- function(parts, simulated_raw) {
  stopifnot(identical(colnames(simulated_raw), parts$lines_all), !anyNA(simulated_raw))
  simulated <- t(apply(simulated_raw, 1, rank_normalise))
  colnames(simulated) <- parts$lines_all
  one_test <- function(part) {
    y <- simulated[, part$lines, drop = FALSE]
    estimate <- drop(y %*% part$weights)
    sigma2 <- (rowSums(y^2) - rowSums((y %*% part$basis)^2)) / part$df_new
    df_pooled <- sum(part$df_residual, na.rm = TRUE) + part$df_new
    last <- length(part$sigma2) + 1
    vapply(seq_len(nrow(y)), function(i) {
      squeezed <- squeezeVar(c(part$sigma2, sigma2[i]), c(part$df_residual, part$df_new), robust = TRUE)
      df_prior <- if (length(squeezed$df.prior) > 1) squeezed$df.prior[last] else squeezed$df.prior
      2 * pt(-abs(estimate[i] / part$unscaled / sqrt(squeezed$var.post[last])), min(part$df_new + df_prior, df_pooled))
    }, numeric(1))
  }
  data.table(row = seq_len(nrow(simulated)), interaction_p = one_test(parts$interaction), inside_p = one_test(parts$inside))
}

# The reference for one_more_tests(): simulated genes are appended to the real
# matrix and the whole screen is refitted, at most `batch` simulated genes per
# fit. With batch = 1 it is slow and exact, and one_more_tests() is tested
# against it. Larger batches shift limma's variance prior: against `baseline`,
# the fit without simulated rows, `drift` is the largest change in log10 p among
# the real rows and `drift_q99` the 99th percentile of that change.
embed_and_fit <- function(mat_rank, simulated_raw, ctx_frame, covariates, batch, baseline) {
  stopifnot(identical(colnames(simulated_raw), colnames(mat_rank)), nrow(baseline) == nrow(mat_rank))
  simulated <- t(apply(simulated_raw, 1, rank_normalise))
  dimnames(simulated) <- list(paste0("sim__", seq_len(nrow(simulated_raw))), colnames(mat_rank))
  real <- seq_len(nrow(mat_rank))
  moved <- function(new, old) abs(log10(new) - log10(old))
  drift <- 0
  drift_q99 <- 0
  batches <- split(seq_len(nrow(simulated)), ceiling(seq_len(nrow(simulated)) / batch))
  parts <- lapply(batches, function(rows) {
    tests <- fit_context_tests(ctx_frame, rbind(mat_rank, simulated[rows, , drop = FALSE]), covariates)
    change <- c(moved(tests$interaction_p[real], baseline$interaction_p), moved(tests$inside_p[real], baseline$inside_p))
    drift <<- max(drift, change, na.rm = TRUE)
    drift_q99 <<- max(drift_q99, quantile(change, 0.99, na.rm = TRUE, names = FALSE))
    tests[-real, .(row = rows, interaction_p, inside_p)]
  })
  list(results = rbindlist(parts), drift = drift, drift_q99 = drift_q99)
}

# --- detection curves ------------------------------------------------------------------

# Shift at which detection first reaches each target, by linear interpolation
# between grid points, with the ceiling (detection at the largest shift). A
# target the curve never reaches is reported as not reachable.
detection_curve <- function(shift, detected, targets = c(0.5, 0.8)) {
  ord <- order(shift)
  shift <- shift[ord]
  detected <- detected[ord]
  limit <- vapply(targets, function(target) {
    at <- which(detected >= target)[1]
    if (is.na(at)) return(NA_real_)
    if (at == 1) return(shift[1])
    shift[at - 1] + (target - detected[at - 1]) / (detected[at] - detected[at - 1]) * (shift[at] - shift[at - 1])
  }, numeric(1))
  data.table(target = targets, limit = limit, reachable = !is.na(limit), ceiling = detected[length(detected)])
}
