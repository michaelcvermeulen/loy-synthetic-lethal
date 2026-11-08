# Unit checks for forest_helpers.R (main environment; needs ranger).
# Run: /global/home/hpc5173/miniconda3/envs/main/bin/Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_forest_helpers.R

script_dir <- dirname(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
))))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "forest_helpers.R"))

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

set.seed(21)
n <- 300
lineage <- sample(c("Lung", "Skin", "Bowel", "Blood"), n, replace = TRUE)
x <- matrix(rnorm(n * 40), n, 40, dimnames = list(NULL, paste0("g", 1:40)))
# g1 carries a real within-lineage signal; g2 only marks a lineage.
x[, 2] <- as.numeric(lineage == "Skin") * 3 + rnorm(n, sd = 0.2)
y <- 1.5 * x[, 1] + as.numeric(lineage == "Skin") * 2 + rnorm(n)

# --- shadows -----------------------------------------------------------------
shadows <- make_shadows(x, lineage)
same_content <- all(vapply(unique(lineage), function(l) {
  isTRUE(all.equal(sort(shadows[lineage == l, 1]), sort(x[lineage == l, 1])))
}, logical(1)))
check("make_shadows permutes within strata only", same_content && !identical(shadows[, 1], x[, 1]))
check("make_shadows keeps a lineage marker a lineage marker", cor(shadows[, 2], x[, 2]) > 0.95)

# --- Boruta-style selection ---------------------------------------------------
always <- lineage_columns(lineage)
res <- boruta_stratified(x, y, lineage, always, iterations = 30, trees = 200, seed = 5, threads = 2)
check("planted within-lineage feature is confirmed", res[feature == "g1", confirmed])
check("pure lineage marker is not confirmed", !res[feature == "g2", confirmed])
check("few noise features are confirmed", sum(res[!feature %in% c("g1", "g2"), confirmed]) <= 1)

set.seed(22)
y_null <- y
for (idx in split(seq_len(n), lineage)) y_null[idx] <- y[idx[sample.int(length(idx))]]
res_null <- boruta_stratified(x, y_null, lineage, always, iterations = 30, trees = 200, seed = 6, threads = 2)
check("shuffled outcome confirms almost nothing", sum(res_null$confirmed) <= 1)

label <- factor(as.integer(x[, 1] + rnorm(n) > 0))
res_class <- boruta_stratified(x, label, lineage, always, iterations = 30, trees = 200, seed = 7, threads = 2)
check("classification outcome works and finds the planted feature", res_class[feature == "g1", confirmed])

# --- folds, imputation, AUC ---------------------------------------------------
donor <- c(rep(paste0("d", 1:20), each = 2), paste0("s", seq_len(n - 40)))
set.seed(23)
folds <- donor_folds(donor, k = 10)
check("donor_folds never splits a donor", all(tapply(folds, donor, function(f) length(unique(f))) == 1) && length(unique(folds)) == 10)
x_na <- x
x_na[1:5, 3] <- NA
filled <- impute_within(x_na, lineage)
check("impute_within fills with the stratum median", !anyNA(filled) && isTRUE(all.equal(unname(filled[1, 3]), median(x_na[lineage == lineage[1], 3], na.rm = TRUE))))
check("rank_auc agrees with a perfect and a random score", rank_auc(c(1, 2, 3, 4), c(0, 0, 1, 1)) == 1 && abs(rank_auc(rnorm(2000), rbinom(2000, 1, 0.5)) - 0.5) < 0.05)

cat(sprintf("%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
