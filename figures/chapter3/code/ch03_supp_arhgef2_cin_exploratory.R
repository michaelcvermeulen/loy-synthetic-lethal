# Reconstruct the saved stage-19 ARHGEF2/CIN fits for an exploratory display.
# Base R only; no screen or multiple-testing family is rerun.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
export <- file.path(root, "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07")
derived <- file.path(export, "derived/loy_reanalysis_v1")
out <- file.path(root, "figures/chapter3/source_data/arhgef2_cin_exploratory")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# Load only the four pure helper definitions, avoiding the pipeline bootstrap.
helpers <- parse(file.path(export, "scripts/loy_reanalysis_v1/common.R"))
for (e in helpers) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("full_rank", "lineage_design", "fit_ols", "rank_normalise")) eval(e)
}
cohort <- read.csv(file.path(derived, "cohort.csv"), check.names = FALSE)
context <- read.csv(file.path(derived, "context_matrix.csv"), check.names = FALSE)
saved <- read.csv(file.path(derived, "robustness_context_gene_hits.csv"))
saved <- saved[saved$symbol == "ARHGEF2" & saved$context == "CIN", ]
stopifnot(nrow(saved) == 1)
all_points <- list(); all_curves <- list(); all_fits <- list(); all_difference <- list()
for (assay in c("CRISPR", "RNAi")) {
  matrix_name <- if (assay == "CRISPR") "crispr_effect" else "rnai"
  cat("Reading", assay, "matrix\n"); flush.console()
  mat <- readRDS(file.path(derived, "matrices", paste0(matrix_name, ".rds")))
  feature <- if (assay == "CRISPR") "ARHGEF2 (9181)" else "ARHGEF2"
  values <- mat[, feature]
  frame <- cohort[which(cohort$in_primary & cohort$ModelID %in% rownames(mat)), ]
  rm(mat); gc()
  frame$exposed <- frame$y_negative
  frame$cin <- context$genomic_continuous__CIN[match(frame$ModelID, context$ModelID)]
  frame <- frame[!is.na(frame$cin), ]
  cin_mean <- mean(frame$cin); cin_sd <- sd(frame$cin)
  frame$ctx <- (frame$cin - cin_mean) / cin_sd
  frame$interaction <- frame$exposed * frame$ctx
  frame$raw <- values[frame$ModelID]
  frame <- frame[is.finite(frame$raw), ]
  covariates <- c("exposed", "ctx", if (assay == "CRISPR") c("LibraryKY", "LibraryHumagneCD"))
  frame <- frame[complete.cases(frame[, c("interaction", covariates)]), ]
  X <- cbind(intercept = 1, interaction = frame$interaction,
             lineage_design(frame$OncotreeLineage)[, -1, drop = FALSE])
  for (name in covariates) X <- cbind(X, setNames(data.frame(as.numeric(frame[[name]])), name))
  X <- full_rank(as.matrix(X))
  y <- rank_normalise(frame$raw)
  fit <- lm.fit(X, y)
  beta <- fit$coefficients
  df <- length(y) - fit$rank
  covariance <- solve(crossprod(X)) * sum(fit$residuals^2) / df
  se <- sqrt(diag(covariance))
  crit <- qt(.975, df)
  statistic <- fit_ols(y, frame$interaction, X[, colnames(X) != "interaction", drop = FALSE])
  expected_beta <- if (assay == "CRISPR") saved$interaction_rank else saved$replication_estimate
  expected_p <- if (assay == "CRISPR") saved$interaction_p else saved$replication_p
  stopifnot(abs(statistic$estimate - expected_beta) < 1e-9,
            abs(log(statistic$p / expected_p)) < 1e-7)
  # Remove fitted lineage/library contributions, centred at the assay-wide
  # covariate distribution. Retain CIN, Y status and their interaction.
  nuisance <- setdiff(colnames(X), c("intercept", "interaction", "exposed", "ctx"))
  adjustment <- if (length(nuisance)) as.numeric(sweep(X[, nuisance, drop = FALSE], 2,
                      colMeans(X[, nuisance, drop = FALSE]), "-") %*% beta[nuisance]) else 0
  all_points[[assay]] <- data.frame(assay, ModelID = frame$ModelID,
      cell_line = frame$StrippedCellLineName, lineage = frame$OncotreeLineage,
      y_status = ifelse(frame$exposed == 1, "Y-negative", "Y-positive"),
      cin = frame$cin, cin_z = frame$ctx, raw_dependency = frame$raw,
      rank_dependency = y, adjusted_dependency = y - adjustment)
  curves <- list()
  for (group in 0:1) {
    grid <- seq(min(frame$cin[frame$exposed == group]), max(frame$cin[frame$exposed == group]), length.out = 201)
    XP <- matrix(rep(colMeans(X), each = length(grid)), nrow = length(grid), dimnames = list(NULL, colnames(X)))
    XP[, "ctx"] <- (grid - cin_mean) / cin_sd
    XP[, "exposed"] <- group
    XP[, "interaction"] <- group * XP[, "ctx"]
    mu <- as.numeric(XP %*% beta)
    uncertainty <- sqrt(rowSums((XP %*% covariance) * XP))
    curves[[group + 1]] <- data.frame(assay, y_status = ifelse(group == 1, "Y-negative", "Y-positive"),
        cin = grid, estimate = mu, lower = mu - crit * uncertainty, upper = mu + crit * uncertainty)
  }
  all_curves[[assay]] <- do.call(rbind, curves)
  # Y-negative minus Y-positive contrast, limited to overlapping CIN support.
  lo <- max(tapply(frame$cin, frame$exposed, min)); hi <- min(tapply(frame$cin, frame$exposed, max))
  grid <- seq(lo, hi, length.out = 201)
  C <- matrix(0, nrow = length(grid), ncol = ncol(X), dimnames = list(NULL, colnames(X)))
  C[, "exposed"] <- 1; C[, "interaction"] <- (grid - cin_mean) / cin_sd
  delta <- as.numeric(C %*% beta); delta_se <- sqrt(rowSums((C %*% covariance) * C))
  all_difference[[assay]] <- data.frame(assay, cin = grid, estimate = delta,
      lower = delta - crit * delta_se, upper = delta + crit * delta_se)
  all_fits[[assay]] <- data.frame(assay, n = nrow(frame), n_positive = sum(frame$exposed == 0),
      n_negative = sum(frame$exposed == 1), interaction = statistic$estimate,
      se = statistic$se, df, lower = statistic$estimate - crit * statistic$se,
      upper = statistic$estimate + crit * statistic$se, p = statistic$p,
      cin_mean, cin_sd, verified_against_saved = TRUE)
  cat(assay, "n =", nrow(frame), "interaction =", statistic$estimate, "P =", statistic$p, "\n"); flush.console()
}
write.csv(do.call(rbind, all_points), file.path(out, "points.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_curves), file.path(out, "curves.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_difference), file.path(out, "y_status_contrast.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_fits), file.path(out, "fits.csv"), row.names = FALSE)
write.csv(saved, file.path(out, "saved_robustness_result.csv"), row.names = FALSE)
