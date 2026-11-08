# Unit checks for power_helpers.R (scRNA environment; needs limma).
# Run: Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_power_helpers.R

script_dir <- dirname(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
))))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
source(file.path(script_dir, "power_helpers.R"))

n_pass <- 0L
n_fail <- 0L
check <- function(name, ok) {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
  } else {
    n_fail <<- n_fail + 1L
    cat("FAIL:", name, "\n")
  }
}

# --- bh_next_critical ----------------------------------------------------------
# The new test is rejected exactly when its p-value is below the critical value.
agrees_with_bh <- function(family, alpha = 0.05) {
  crit <- bh_next_critical(family, alpha)
  probes <- c(crit * c(0.2, 0.9, 0.999999), crit * c(1.000001, 1.1, 5))
  rejected <- vapply(probes, function(p) bh(c(family, p))[length(family) + 1] < alpha, logical(1))
  identical(rejected, probes < crit)
}
set.seed(3)
null_family <- runif(400, 0.01, 1)
check("bh_next_critical is alpha / (m + 1) when nothing else passes", isTRUE(all.equal(bh_next_critical(null_family), 0.05 / 401)))
check("bh_next_critical agrees with bh() on a null family", agrees_with_bh(null_family))
hit_family <- c(runif(6, 0, 1e-6), runif(394, 0.01, 1))
check("bh_next_critical is looser when other tests already pass", bh_next_critical(hit_family) > 0.05 / 401)
check("bh_next_critical agrees with bh() with earlier rejections", agrees_with_bh(hit_family))
# One more small p-value can rescue others: 3 values that just miss on their own.
rescue_family <- c(2.6e-4, 3.6e-4, 4.9e-4, runif(396, 0.2, 1))
check("bh_next_critical agrees with bh() when the new test rescues others", agrees_with_bh(rescue_family))
check("bh_next_critical ignores missing p-values", isTRUE(all.equal(bh_next_critical(c(null_family, NA)), 0.05 / 401)))

# --- analytic power -------------------------------------------------------------
check("analytic_mde tends to the normal multiplier at large df", abs(analytic_mde(1, 1e6, 0.05, 0.8) - (qnorm(0.975) + qnorm(0.8))) < 1e-3)
shift <- analytic_mde(0.3, 20, 1e-5, 0.8)
check("analytic_mde inverts analytic_power", abs(analytic_power(shift, 0.3, 20, 1e-5) - 0.8) < 1e-6)
reference <- power.t.test(n = 12, sd = 1, sig.level = 0.05, power = 0.8, strict = TRUE)$delta
check("analytic_mde matches power.t.test for a balanced design", abs(analytic_mde(sqrt(2 / 12), 22, 0.05, 0.8) - reference) < 1e-3)
check("analytic_power is the test size at zero shift", abs(analytic_power(0, 1, 30, 0.01) - 0.01) < 1e-8)

# --- exact_floor -----------------------------------------------------------------
check("exact_floor is 2 / choose(n, n_negative)", isTRUE(all.equal(exact_floor(5, 20), 2 / choose(25, 5))))
separated <- wilcox.test(1:6, 7:13, exact = TRUE)$p.value
check("exact_floor matches an exact Wilcoxon test on separated groups", isTRUE(all.equal(exact_floor(6, 7), separated)))
check("exact_floor never exceeds one", exact_floor(1, 1) == 1)

# --- spike_cell ------------------------------------------------------------------
set.seed(5)
raw <- matrix(rnorm(8 * 10), 8, 10, dimnames = list(paste0("G", 1:8), paste0("L", 1:10)))
target <- rep(c(TRUE, FALSE), c(3, 7))
spiked <- spike_cell(raw, c("G2", "G5"), target, -0.7)
check("spike_cell leaves lines outside the cell unchanged", identical(spiked[, !target], raw[c("G2", "G5"), !target]))
check("spike_cell adds the shift inside the cell", isTRUE(all.equal(spiked[, target], raw[c("G2", "G5"), target] - 0.7)))
per_gene <- spike_cell(raw, c("G2", "G5"), target, c(-1, -2))
check("spike_cell takes one shift per gene", isTRUE(all.equal(per_gene[, target] - raw[c("G2", "G5"), target], matrix(c(-1, -2), 2, 3), check.attributes = FALSE)))
check("spike_cell keeps gene and line names", identical(dimnames(spiked), list(c("G2", "G5"), colnames(raw))))

# --- detection_curve -------------------------------------------------------------
rising <- detection_curve(c(0, 1, 2, 3), c(0, 0.4, 0.9, 1), c(0.5, 0.8))
check("detection_curve interpolates the 50% point", abs(rising$limit[1] - 1.2) < 1e-9)
check("detection_curve interpolates the 80% point", abs(rising$limit[2] - 1.8) < 1e-9)
check("detection_curve reports the ceiling", rising$ceiling[1] == 1 && all(rising$reachable))
plateau <- detection_curve(c(0, 1, 2, 3), c(0, 0.3, 0.6, 0.65), c(0.5, 0.8))
check("detection_curve reports a plateau as not reachable", is.na(plateau$limit[2]) && !plateau$reachable[2] && plateau$ceiling[2] == 0.65)
check("detection_curve still gives the reachable point of a plateau", abs(plateau$limit[1] - (1 + 2 / 3)) < 1e-9)
unsorted <- detection_curve(c(3, 0, 2, 1), c(1, 0, 0.9, 0.4), 0.5)
check("detection_curve sorts by shift", abs(unsorted$limit - 1.2) < 1e-9)
early <- detection_curve(c(1, 2), c(0.9, 1), 0.8)
check("detection_curve returns the first shift when already above target", early$limit == 1)

# --- benchmark_fit and benchmark_replicates ---------------------------------------
set.seed(17)
n <- 400
frame <- data.table(
  ModelID = sprintf("ACH-%06d", seq_len(n)),
  OncotreeLineage = sample(c("Lung", "Skin", "Bowel", "Kidney"), n, replace = TRUE),
  LibraryKY = rbinom(n, 1, 0.3)
)
biomarker <- rbinom(n, 1, 0.25)
biomarker[1:5] <- NA
lineage_mean <- c(Lung = 0, Skin = -0.5, Bowel = 0.3, Kidney = 0.1)[frame$OncotreeLineage]
altered <- biomarker %in% 1
gene_values <- lineage_mean + 0.2 * frame$LibraryKY + ifelse(altered, -1 + rnorm(n, 0, 0.6), rnorm(n, 0, 0.2))
effect <- cbind(gene_values, rnorm(n))
dimnames(effect) <- list(frame$ModelID, c("GENE1 (1)", "GENE2 (2)"))
fit <- benchmark_fit(frame, effect, biomarker, "GENE1 (1)", "LibraryKY")
check("benchmark_fit recovers the biomarker effect", abs(fit$coefficient + 1) < 0.15)
check("benchmark_fit counts altered and reference lines", fit$n_altered == sum(biomarker %in% 1) && fit$n_reference == sum(biomarker %in% 0))
check("benchmark_fit keeps one residual per line in each pool", length(fit$altered) == fit$n_altered && length(fit$reference) == fit$n_reference)
check("benchmark_fit sees the wider spread of the altered group", fit$sd_altered > 2 * fit$sd_reference)
plain <- coef(summary(lm(gene_values ~ factor(frame$OncotreeLineage) + frame$LibraryKY + biomarker)))["biomarker", ]
check("benchmark_fit agrees with lm", abs(fit$coefficient - plain[["Estimate"]]) < 1e-10 && abs(fit$p / plain[["Pr(>|t|)"]] - 1) < 1e-6)
nuisance <- benchmark_nuisance(fit, frame)
check("benchmark_nuisance follows lineage and library, without the biomarker", abs(mean(nuisance[frame$OncotreeLineage == "Skin" & frame$LibraryKY == 0]) - mean(nuisance[frame$OncotreeLineage == "Lung" & frame$LibraryKY == 0]) + 0.5) < 0.15)

cell <- frame$OncotreeLineage == "Kidney" & seq_len(n) %% 2 == 0
set.seed(99)
rep_a <- benchmark_replicates(fit, frame, cell, 300, with_spread = TRUE)
set.seed(99)
rep_b <- benchmark_replicates(fit, frame, cell, 300, with_spread = TRUE)
check("benchmark_replicates is reproducible under a seed", identical(rep_a, rep_b))
check("benchmark_replicates returns one row per replicate and one column per line", identical(dim(rep_a), c(300L, 400L)) && identical(colnames(rep_a), frame$ModelID))
centred <- sweep(rep_a, 2, nuisance)
check("benchmark_replicates leaves lines outside the cell centred", abs(mean(centred[, !cell])) < 0.01)
check("benchmark_replicates moves the cell by the coefficient", abs(mean(centred[, cell]) - fit$coefficient) < 0.03)
check("benchmark_replicates gives the cell the altered spread", sd(centred[, cell]) > 2 * sd(centred[, !cell]))
set.seed(99)
rep_flat <- benchmark_replicates(fit, frame, cell, 300, with_spread = FALSE)
flat <- sweep(rep_flat, 2, nuisance)
check("benchmark_replicates without spread keeps the reference spread in the cell", abs(sd(flat[, cell]) / sd(flat[, !cell]) - 1) < 0.1 && abs(mean(flat[, cell]) - fit$coefficient) < 0.03)

# --- fit_context_tests and embed_and_fit -------------------------------------------
set.seed(23)
n <- 120
ctx_source <- data.table(
  ModelID = sprintf("ACH-%06d", seq_len(n)),
  OncotreeLineage = sample(c("Lung", "Skin", "Bowel"), n, replace = TRUE),
  exposed = rbinom(n, 1, 0.45),
  LibraryKY = rbinom(n, 1, 0.3)
)
ctx_frame <- with_context(ctx_source, as.numeric(ctx_source$OncotreeLineage == "Skin"), "binary")
genes <- paste0("G", 1:300, " (", 1:300, ")")
mat_raw <- matrix(rnorm(300 * n), 300, n, dimnames = list(genes, ctx_source$ModelID))
mat_rank <- t(apply(mat_raw, 1, rank_normalise))
dimnames(mat_rank) <- dimnames(mat_raw)
tests <- fit_context_tests(ctx_frame, mat_rank, "LibraryKY")
design_int <- build_design(ctx_frame, tested = "interaction", covariates = c("exposed", "ctx", "LibraryKY"))
direct_int <- limma_scan(mat_rank[, ctx_frame$ModelID], design_int$design, "interaction", trend = FALSE)
inside_frame <- ctx_frame[ctx == 1]
design_in <- build_design(inside_frame, tested = "exposed", covariates = "LibraryKY")
direct_in <- limma_scan(mat_rank[, inside_frame$ModelID], design_in$design, "exposed", trend = FALSE)
check("fit_context_tests gives the interaction p-values of the screen's fit", isTRUE(all.equal(tests$interaction_p, direct_int$p)))
check("fit_context_tests gives the inside p-values of the screen's fit", isTRUE(all.equal(tests$inside_p, direct_in$p)))
check("fit_context_tests returns one row per feature, in order", identical(tests$feature, genes))
check("fit_context_tests returns standard errors and degrees of freedom", all(tests$inside_se > 0) && all(tests$inside_df > nrow(inside_frame) - 3))

target <- ctx_frame$exposed == 1 & ctx_frame$ctx == 1
hosts <- genes[1:25]
zero <- spike_cell(mat_raw[, ctx_frame$ModelID], hosts, target, 0)
embedded <- embed_and_fit(mat_rank, zero, ctx_frame, "LibraryKY", batch = 10, baseline = tests)
check("embed_and_fit returns one row per simulated gene, in order", nrow(embedded$results) == 25 && identical(embedded$results$row, seq_len(25)))
# A copy of a real gene sits in the same fit as that gene, so the two p-values are equal
# there; against the fit without simulated rows they differ by no more than the reported drift.
host <- tests[match(hosts, feature)]
gap <- max(abs(log10(embedded$results$inside_p) - log10(host$inside_p)), abs(log10(embedded$results$interaction_p) - log10(host$interaction_p)))
check("embed_and_fit reports how far the real p-values moved", is.finite(embedded$drift) && embedded$drift >= 0 && embedded$drift < 0.2)
# The largest move belongs to the few most extreme genes; the 99th percentile describes the rest.
check("embed_and_fit reports the 99th percentile of the move as well", is.numeric(embedded$drift_q99) && length(embedded$drift_q99) == 1 && embedded$drift_q99 >= 0 && embedded$drift_q99 <= embedded$drift)
check("a zero shift returns the host gene's p-value, up to the reported drift", gap <= embedded$drift + 1e-12)
strong <- spike_cell(mat_raw[, ctx_frame$ModelID], hosts, target, -6)
shifted <- embed_and_fit(mat_rank, strong, ctx_frame, "LibraryKY", batch = 10, baseline = tests)
check("a large shift in the cell gives small inside p-values", median(shifted$results$inside_p) < 1e-4)
check("a large shift in the cell gives small interaction p-values", median(shifted$results$interaction_p) < 1e-4)
check("embed_and_fit refuses simulated rows that do not cover the frame's lines", inherits(try(embed_and_fit(mat_rank, strong[, 1:10], ctx_frame, "LibraryKY", batch = 10, baseline = tests), silent = TRUE), "try-error"))

# --- context_fit_parts and one_more_tests ---------------------------------------------
# The fast path must give what the reference path gives with one simulated gene per fit.
agrees_with_reference <- function(real_rank, rows) {
  parts <- context_fit_parts(ctx_frame, real_rank, "LibraryKY")
  reference_base <- suppressWarnings(fit_context_tests(ctx_frame, real_rank, "LibraryKY"))
  fast <- one_more_tests(parts, rows)
  slow <- suppressWarnings(embed_and_fit(real_rank, rows, ctx_frame, "LibraryKY", batch = 1, baseline = reference_base))$results
  isTRUE(all.equal(fast$inside_p, slow$inside_p, tolerance = 1e-8)) && isTRUE(all.equal(fast$interaction_p, slow$interaction_p, tolerance = 1e-8)) && identical(fast$row, slow$row)
}
check("one_more_tests equals fitting each simulated gene alone beside the real ones", agrees_with_reference(mat_rank, strong[1:8, ]))
check("one_more_tests equals the reference for genes with no shift", agrees_with_reference(mat_rank, zero[1:6, ]))
# A gene far more variable than any real one changes how limma treats outliers.
set.seed(31)
wild <- spike_cell(mat_raw[, ctx_frame$ModelID], hosts[1:6], target, -2)
wild[, target] <- wild[, target] + rnorm(6 * sum(target), 0, 4)
check("one_more_tests equals the reference for a very variable gene", agrees_with_reference(mat_rank, wild))
# Real rows with missing values go down limma's separate path, and their degrees of freedom differ.
holes <- mat_raw
holes[5, 1:9] <- NA
holes[40, c(3, 50, 77)] <- NA
holes[41, which(ctx_frame$ctx == 1)[1:4]] <- NA
holes_rank <- t(apply(holes, 1, rank_normalise))
dimnames(holes_rank) <- dimnames(holes)
check("one_more_tests equals the reference when real genes have missing values", agrees_with_reference(holes_rank, strong[1:6, ]))
check("one_more_tests refuses simulated genes with missing values", inherits(try(one_more_tests(context_fit_parts(ctx_frame, mat_rank, "LibraryKY"), holes[5, , drop = FALSE]), silent = TRUE), "try-error"))

cat(sprintf("%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
