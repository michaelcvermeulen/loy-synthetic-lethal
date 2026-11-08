# Unit checks for scan_helpers.R (scRNA environment; needs limma).
# Run: Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_scan_helpers.R

script_dir <- dirname(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
))))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))

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

set.seed(11)
n <- 90
frame <- data.table(
  ModelID = sprintf("ACH-%06d", seq_len(n)),
  OncotreeLineage = sample(c("Lung", "Skin", "Bowel"), n, replace = TRUE),
  exposed = rbinom(n, 1, 0.45),
  LibraryKY = rbinom(n, 1, 0.3),
  WGD = c(NA, NA, rbinom(n - 2, 1, 0.5))
)
mat_t <- matrix(rnorm(40 * n), 40, n, dimnames = list(paste0("G", 1:40, " (", 1:40, ")"), frame$ModelID))
mat_t[1, ] <- mat_t[1, ] - 1.2 * frame$exposed

# --- build_design ------------------------------------------------------------
base <- build_design(frame, tested = "exposed", covariates = "LibraryKY")
check("build_design puts the tested column second", identical(colnames(base$design)[1:2], c("intercept", "exposed")))
check("build_design keeps all complete rows", length(base$rows) == n && base$estimable)
adjusted <- build_design(frame, tested = "exposed", covariates = c("LibraryKY", "WGD"))
check("build_design drops rows with a missing covariate", identical(adjusted$rows, 3:n) && nrow(adjusted$design) == n - 2)
aliased <- copy(frame)[, exposed := as.numeric(OncotreeLineage == "Lung")]
check("build_design flags an exposure aliased with lineage", !build_design(aliased, tested = "exposed")$estimable)

# --- fast_t ------------------------------------------------------------------
nuisance <- base$design[, colnames(base$design) != "exposed", drop = FALSE]
t_fast <- fast_t(mat_t, frame$exposed, nuisance)
t_lm <- vapply(seq_len(nrow(mat_t)), function(i) {
  summary(lm(mat_t[i, ] ~ 0 + base$design))$coefficients[2, "t value"]
}, numeric(1))
check("fast_t equals the lm t statistic", isTRUE(all.equal(unname(t_fast), t_lm, tolerance = 1e-8)))

# --- limma_scan --------------------------------------------------------------
scan <- limma_scan(mat_t, base$design, "exposed", trend = FALSE)
beta_lm <- vapply(seq_len(nrow(mat_t)), function(i) unname(coef(lm(mat_t[i, ] ~ 0 + base$design))[2]), numeric(1))
check("limma_scan estimate equals lm", isTRUE(all.equal(scan$estimate, beta_lm, tolerance = 1e-8)))
check("limma_scan finds the planted effect", which.min(scan$p) == 1 && scan$estimate[1] < -0.8)
mat_na <- mat_t
mat_na[2, 1:9] <- NA
scan_na <- limma_scan(mat_na, base$design, "exposed", trend = FALSE)
check("limma_scan counts non-missing samples", scan_na$n[2] == n - 9)
beta_na <- unname(coef(lm(mat_na[2, ] ~ 0 + base$design))[2])
check("limma_scan keeps feature order when some rows have missing values", identical(scan_na$feature, rownames(mat_na)) &&
  isTRUE(all.equal(scan_na$estimate[2], beta_na, tolerance = 1e-8)) &&
  isTRUE(all.equal(scan_na$estimate[-2], beta_lm[-2], tolerance = 1e-8)))
donor <- c(rep(1:10, each = 2), 11:(n - 10))
blocked <- limma_scan(mat_t, base$design, "exposed", block = donor, trend = FALSE)
check("limma_scan with donor blocking returns a correlation", is.finite(attr(blocked, "correlation")) && nrow(blocked) == 40)

# --- group_summary -----------------------------------------------------------
gs <- group_summary(mat_na, frame$exposed)
check("group_summary sizes add up", all(gs$n_negative + gs$n_positive == rowSums(!is.na(mat_na))))

cat(sprintf("%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
