# Unit checks for the hand-coded helpers in common.R.
# Run: Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_common.R
# Works in both the scRNA (R 4.5.3) and main (R 4.2.2) environments.

script_dir <- dirname(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
))))
source(file.path(script_dir, "common.R"))

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

# --- identifiers -------------------------------------------------------------
check("gene_symbol strips Entrez suffix", identical(gene_symbol(c("DDX3X (1654)", "XIST")), c("DDX3X", "XIST")))
check("gene_entrez parses id", identical(gene_entrez(c("DDX3X (1654)", "XIST")), c("1654", NA_character_)))

# --- multiple testing --------------------------------------------------------
p <- c(0.001, 0.02, NA, 0.5, 0.03)
expected <- rep(NA_real_, 5)
expected[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
check("bh matches p.adjust and keeps NA", isTRUE(all.equal(bh(p), expected)))

# --- permutation within strata ----------------------------------------------
set.seed(1)
x <- c(1, 0, 0, 1, 1, 0, 1)
strata <- c("a", "a", "a", "b", "b", "b", "c")
perm <- permute_within(x, strata)
check("permute_within keeps per-stratum totals", identical(tapply(perm, strata, sum), tapply(x, strata, sum)))
check("permute_within leaves singleton strata", perm[7] == x[7])

# --- robust single-coefficient fit ------------------------------------------
set.seed(2)
n <- 120
lineage <- factor(sample(c("l1", "l2", "l3", "l4"), n, replace = TRUE))
x <- rbinom(n, 1, 0.4)
tech <- rnorm(n)
y <- 0.6 * x + as.numeric(lineage) * 0.3 + 0.2 * tech + rnorm(n) * (1 + x)
z <- model.matrix(~ lineage + tech)
fit <- fit_robust(y, x, z, cluster = seq_len(n))
ref <- lm(y ~ lineage + tech + x)
check("fit_robust estimate equals lm", isTRUE(all.equal(fit$estimate, unname(coef(ref)["x"]), tolerance = 1e-10)))

brute_se <- function(y, x, z, cluster) {
  full <- cbind(z, x)
  xtx_inv <- solve(crossprod(full))
  hat <- rowSums((full %*% xtx_inv) * full)
  resid <- as.numeric(y - full %*% (xtx_inv %*% crossprod(full, y)))
  score <- full * (resid / (1 - hat))
  score <- rowsum(score, cluster, reorder = FALSE)
  vcov <- xtx_inv %*% crossprod(score) %*% xtx_inv
  sqrt(vcov[ncol(full), ncol(full)])
}
check("fit_robust HC3 equals brute force", isTRUE(all.equal(fit$se, brute_se(y, x, z, seq_len(n)), tolerance = 1e-8)))

cluster <- c(rep(1:20, each = 2), 21:100)
fit_cl <- fit_robust(y, x, z, cluster = cluster)
check("fit_robust clustered equals brute force", isTRUE(all.equal(fit_cl$se, brute_se(y, x, z, cluster), tolerance = 1e-8)))
check("fit_robust reports group sizes", fit$n == n && fit$n_exposed == sum(x) && fit$n_unexposed == sum(1 - x))

y_na <- y
y_na[c(3, 50)] <- NA
fit_na <- fit_robust(y_na, x, z, cluster = seq_len(n))
keep <- !is.na(y_na)
check("fit_robust drops missing responses", fit_na$n == n - 2 && isTRUE(all.equal(
  fit_na$estimate, unname(coef(lm(y_na[keep] ~ lineage[keep] + tech[keep] + x[keep]))[6]), tolerance = 1e-10
)))

# a singleton stratum has leverage one; it must contribute nothing and not break the fit
lineage_s <- factor(c(as.character(lineage[-n]), "solo"))
z_s <- model.matrix(~ lineage_s + tech)
fit_s <- fit_robust(y, x, z_s, cluster = seq_len(n))
check("fit_robust handles a singleton stratum", is.finite(fit_s$se) && is.finite(fit_s$p))

# exposure constant within every stratum: not estimable
x_const <- as.numeric(lineage %in% c("l1", "l2"))
fit_c <- fit_robust(y, x_const, model.matrix(~ lineage), cluster = seq_len(n))
check("fit_robust returns NA when not estimable", is.na(fit_c$estimate))

# calibration: under the null with heteroskedastic noise the test is not anti-conservative
set.seed(3)
p_null <- vapply(seq_len(400), function(i) {
  y0 <- as.numeric(lineage) * 0.3 + rnorm(n) * (1 + 2 * x)
  fit_robust(y0, x, model.matrix(~ lineage), cluster = seq_len(n))$p
}, numeric(1))
check("fit_robust null rejection rate near nominal", mean(p_null < 0.05) < 0.09)

# --- classical OLS fit and rank normalisation --------------------------------
fit_plain <- fit_ols(y, x, z)
ref_row <- summary(ref)$coefficients["x", ]
check("fit_ols equals lm estimate, t and p", isTRUE(all.equal(unname(c(fit_plain$estimate, fit_plain$t, fit_plain$p)),
  unname(ref_row[c("Estimate", "t value", "Pr(>|t|)")]), tolerance = 1e-8)))
check("fit_ols returns NA when not estimable", is.na(fit_ols(y, x_const, model.matrix(~ lineage))$estimate))
v <- c(3.2, NA, -1, 10, 0.5)
rn <- rank_normalise(v)
check("rank_normalise keeps order and missing values", is.na(rn[2]) && identical(order(rn[-2]), order(v[-2])) && abs(mean(rn, na.rm = TRUE)) < 1e-12)

# --- scan over many outcomes -------------------------------------------------
set.seed(4)
y_mat <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, paste0("f", 1:6)))
y_mat[, 2] <- y_mat[, 2] + x
y_mat[sample(n, 10), 4] <- NA
scan <- scan_robust(y_mat, x, z, cluster = seq_len(n), cores = 1)
single <- fit_robust(y_mat[, 4], x, z, cluster = seq_len(n))
check("scan_robust returns one row per feature", nrow(scan) == 6 && identical(scan$feature, colnames(y_mat)))
check("scan_robust matches fit_robust", isTRUE(all.equal(scan$estimate[4], single$estimate)) &&
  isTRUE(all.equal(scan$p[4], single$p)))

# --- Cochran's Q -------------------------------------------------------------
q <- cochran_q(c(0.5, 0.5, 0.5), c(0.1, 0.2, 0.3))
check("cochran_q is zero for identical estimates", abs(q$q) < 1e-12 && q$df == 2)
q2 <- cochran_q(c(0, 1), c(0.1, 0.1))
check("cochran_q detects heterogeneity", isTRUE(all.equal(q2$q, 50)) && q2$p < 1e-10)

# --- design helpers ----------------------------------------------------------
m <- cbind(intercept = 1, x = x, a = as.numeric(lineage == "l1"), b = as.numeric(lineage == "l1"))
fr <- full_rank(m)
check("full_rank drops the aliased column and keeps the exposure", ncol(fr) == 3 && "x" %in% colnames(fr))

cat(sprintf("%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
