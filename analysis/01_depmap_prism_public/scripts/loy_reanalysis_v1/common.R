# Shared helpers for loy_reanalysis_v1.
# Every stage script defines `script_dir` and sources this file from there, so a
# code snapshot is self-contained. Only packages present in both the scRNA
# (R 4.5.3) and main (R 4.2.2) environments are used here.

suppressPackageStartupMessages(library(data.table))

root_dir <- Sys.getenv(
  "LOY_REANALYSIS_ROOT",
  unset = "/global/home/hpc5173/projects/PRISM/data/CCLE/26Q1"
)
analysis_name <- "loy_reanalysis_v1"
raw_dir <- file.path(root_dir, "raw")
depmap_dir <- file.path(raw_dir, "depmap_public_26q1")
drug_dir <- file.path(raw_dir, "drug_response")
derived_dir <- file.path(root_dir, "derived", analysis_name)
results_dir <- file.path(root_dir, "results", analysis_name)
matrix_dir <- file.path(derived_dir, "matrices")
params <- jsonlite::fromJSON(file.path(script_dir, "parameters.json"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
}

# Every source table this analysis reads, relative to the collection root.
input_paths <- c(
  model = "raw/depmap_public_26q1/metadata/Model.csv",
  model_condition = "raw/depmap_public_26q1/metadata/ModelCondition.csv",
  gene = "raw/depmap_public_26q1/metadata/Gene.csv",
  portal_compounds = "raw/depmap_public_26q1/metadata/PortalCompounds.csv",
  crispr_effect = "raw/depmap_public_26q1/crispr/CRISPRGeneEffect.csv",
  crispr_dependency = "raw/depmap_public_26q1/crispr/CRISPRGeneDependency.csv",
  crispr_confounders = "raw/depmap_public_26q1/crispr/CRISPRConfounders.csv",
  expr_tpm_pc = "raw/depmap_public_26q1/expression/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  expr_tpm_all = "raw/depmap_public_26q1/expression/OmicsExpressionTPMLogp1HumanAllGenes.csv",
  expr_counts = "raw/depmap_public_26q1/expression/OmicsExpressionRawReadCountHumanAllGenes.csv",
  cn_wgs = "raw/depmap_public_26q1/copy_number/OmicsCNGeneWGS.csv",
  cn_wes = "raw/depmap_public_26q1/copy_number/OmicsCNGeneMC_WES.csv",
  cn_segments = "raw/depmap_public_26q1/copy_number/OmicsCNSegmentsWGS.csv",
  signatures = "raw/depmap_public_26q1/signatures/OmicsGlobalSignatures.csv",
  mutations = "raw/depmap_public_26q1/mutations/OmicsSomaticMutations.csv",
  mutation_hotspot = "raw/depmap_public_26q1/mutations/OmicsSomaticMutationsMatrixHotspot.csv",
  molecular_subtypes = "raw/depmap_public_26q1/signatures/OmicsInferredMolecularSubtypes.csv",
  rnai = "raw/rnai/combined_demeter2/RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2)_subsetted.csv",
  prism_primary = "raw/drug_response/prism/24Q2/primary/Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv",
  prism_compounds = "raw/drug_response/prism/24Q2/primary/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv",
  prism_secondary = "raw/drug_response/prism/20Q2/secondary/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv",
  gdsc = "raw/drug_response/gdsc/sanger-dose-response.csv",
  ctrp = "derived/ctrp/CTRP_v2.0_AUC_by_ModelID.csv",
  ctrp_cell_map = "derived/ctrp/CTRP_DepMap_cellline_map.csv",
  ctrp_raw_curves = "raw/drug_response/ctrp/CTRPv2.0/v20.data.curves_post_qc.txt",
  ctrp_raw_experiments = "raw/drug_response/ctrp/CTRPv2.0/v20.meta.per_experiment.txt"
)

input_path <- function(key) file.path(root_dir, input_paths[[key]])

# --- identifiers -------------------------------------------------------------

# DepMap gene headers look like "SYMBOL (EntrezID)"; keep the full header as the
# key and parse symbols only for matching and display.
gene_symbol <- function(x) sub(" \\(.*$", "", x)

gene_entrez <- function(x) {
  out <- sub("^.* \\(([0-9]+)\\)$", "\\1", x)
  out[!grepl(" \\([0-9]+\\)$", x)] <- NA_character_
  out
}

# Log how many IDs are absent from a reference and keep the unmatched ones.
audit_ids <- function(label, ids, reference) {
  ids <- unique(ids[!is.na(ids) & ids != ""])
  missing <- setdiff(ids, reference)
  log_msg("ids %s: %d unique, %d not in reference", label, length(ids), length(missing))
  if (length(missing) > 0) {
    audit_dir <- file.path(derived_dir, "join_audit")
    dir.create(audit_dir, showWarnings = FALSE, recursive = TRUE)
    fwrite(
      data.table(label = label, id = missing),
      file.path(audit_dir, paste0(gsub("[^A-Za-z0-9]+", "_", label), "_unmatched.csv"))
    )
  }
  invisible(missing)
}

assert_unique <- function(x, label) {
  if (anyDuplicated(x)) stop(sprintf("%s: key is not unique (%d duplicates)", label, sum(duplicated(x))))
  invisible(TRUE)
}

# --- readers -----------------------------------------------------------------

read_header <- function(path) names(fread(path, nrows = 0, showProgress = FALSE))

# Model-by-feature table whose first, unnamed column is the model ID
# (CRISPR effect, dependency probability, confounders, RNAi).
read_model_matrix <- function(path, features = NULL) {
  header <- read_header(path)
  keep <- header
  if (!is.null(features)) {
    keep <- c(header[1], header[-1][gene_symbol(header[-1]) %in% features | header[-1] %in% features])
  }
  dt <- fread(path, select = keep, showProgress = FALSE)
  ids <- dt[[1]]
  assert_unique(ids, basename(path))
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- ids
  mat
}

# Profile-level table with DepMap's leading metadata columns. Returns one row
# per model, taken from the profile flagged as that model's default.
read_default_profiles <- function(path, features = NULL) {
  header <- read_header(path)
  meta <- intersect(
    c("SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"),
    header
  )
  feature_cols <- setdiff(header, c(meta, "V1"))
  if (!is.null(features)) {
    feature_cols <- feature_cols[gene_symbol(feature_cols) %in% features | feature_cols %in% features]
  }
  dt <- fread(path, select = c(meta, feature_cols), showProgress = FALSE)
  dt <- dt[IsDefaultEntryForModel == "Yes"]
  assert_unique(dt$ModelID, basename(path))
  mat <- as.matrix(dt[, feature_cols, with = FALSE])
  rownames(mat) <- dt$ModelID
  mat
}

write_result <- function(dt, name, dir = derived_dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(dt, file.path(dir, paste0(name, ".csv")))
  saveRDS(dt, file.path(dir, paste0(name, ".rds")))
  log_msg("wrote %s (%d rows)", name, nrow(dt))
  invisible(dt)
}

# Read a table written by write_result. The .rds keeps NA and column types,
# which the .csv does not (a missing call is written as an empty string).
read_result <- function(name, dir = derived_dir) {
  out <- readRDS(file.path(dir, paste0(name, ".rds")))
  setDT(out)
  out
}

# --- freeze ------------------------------------------------------------------

file_sha256 <- function(path) digest::digest(file = path, algo = "sha256")

freeze_file <- file.path(derived_dir, "freeze.json")

frozen_paths <- function() {
  c(
    parameters = file.path(script_dir, "parameters.json"),
    y_status_calls = file.path(derived_dir, "y_status_calls.csv"),
    cohort = file.path(derived_dir, "cohort.csv"),
    context_matrix = file.path(derived_dir, "context_matrix.csv")
  )
}

# Stop if parameters, calls, cohort or contexts changed after the freeze.
assert_freeze <- function() {
  if (!file.exists(freeze_file)) stop("freeze.json is missing; run 07_freeze.R first")
  frozen <- jsonlite::fromJSON(freeze_file)$sha256
  paths <- frozen_paths()
  for (key in names(paths)) {
    if (!identical(file_sha256(paths[[key]]), frozen[[key]])) {
      stop(sprintf("%s changed after the freeze", key))
    }
  }
  log_msg("freeze verified")
  invisible(TRUE)
}

# --- statistics --------------------------------------------------------------

bh <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

# Shuffle x inside each stratum; strata of size one are left as they are.
permute_within <- function(x, strata) {
  out <- x
  for (idx in split(seq_along(x), strata)) {
    if (length(idx) > 1) out[idx] <- x[idx[sample.int(length(idx))]]
  }
  out
}

# Drop aliased columns, keeping earlier columns in preference to later ones.
full_rank <- function(m) {
  decomposition <- qr(m)
  keep <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  m[, keep, drop = FALSE]
}

lineage_design <- function(lineage) {
  lineage <- factor(as.character(lineage))
  if (nlevels(lineage) < 2) {
    return(matrix(1, length(lineage), 1, dimnames = list(NULL, "intercept")))
  }
  design <- model.matrix(~ lineage)
  colnames(design)[1] <- "intercept"
  design
}

# OLS estimate of one coefficient (x) adjusted for z, with a cluster-robust
# HC3-type standard error. x and y are residualised on z (Frisch-Waugh-Lovell),
# each residual is inflated by 1 / (1 - leverage), and scores are summed within
# cluster before squaring. With one observation per cluster this is exactly HC3.
# Rows with leverage one (a stratum of size one) carry no information and are
# given zero weight.
fit_robust <- function(y, x, z, cluster = seq_along(y)) {
  out <- list(
    estimate = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, df = NA_real_,
    n = 0L, n_exposed = NA_integer_, n_unexposed = NA_integer_
  )
  ok <- is.finite(y) & is.finite(x) & (rowSums(!is.finite(z)) == 0)
  out$n <- sum(ok)
  if (out$n < 4) return(out)
  y <- y[ok]
  x <- x[ok]
  z <- z[ok, , drop = FALSE]
  cluster <- cluster[ok]
  if (all(x %in% c(0, 1))) {
    out$n_exposed <- as.integer(sum(x == 1))
    out$n_unexposed <- as.integer(sum(x == 0))
  }
  decomposition <- qr(z)
  x_res <- qr.resid(decomposition, x)
  sxx <- sum(x_res^2)
  if (sxx < 1e-10 * length(x)) return(out)
  y_res <- qr.resid(decomposition, y)
  estimate <- sum(x_res * y_res) / sxx
  resid <- y_res - estimate * x_res
  basis <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
  leverage <- rowSums(basis^2) + x_res^2 / sxx
  weight <- ifelse(leverage < 1 - 1e-8, 1 / (1 - leverage), 0)
  cluster_score <- rowsum(x_res * resid * weight, cluster, reorder = FALSE)
  se <- sqrt(sum(cluster_score^2)) / sxx
  df <- max(1, nrow(cluster_score) - decomposition$rank - 1)
  out$estimate <- estimate
  out$se <- se
  out$t <- estimate / se
  out$p <- 2 * pt(-abs(out$t), df)
  out$df <- df
  out
}

# OLS estimate of one coefficient (x) adjusted for z with the classical,
# model-based standard error. Used on rank-normalised outcomes, where small
# cells make a robust standard error unreliable.
fit_ols <- function(y, x, z) {
  out <- list(estimate = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, df = NA_real_, n = 0L)
  ok <- is.finite(y) & is.finite(x) & (rowSums(!is.finite(z)) == 0)
  out$n <- sum(ok)
  if (out$n < 4) return(out)
  y <- y[ok]
  x <- x[ok]
  z <- z[ok, , drop = FALSE]
  decomposition <- qr(z)
  x_res <- qr.resid(decomposition, x)
  sxx <- sum(x_res^2)
  df <- length(y) - decomposition$rank - 1
  if (sxx < 1e-10 * length(x) || df < 1) return(out)
  y_res <- qr.resid(decomposition, y)
  estimate <- sum(x_res * y_res) / sxx
  se <- sqrt(sum((y_res - estimate * x_res)^2) / df / sxx)
  out$estimate <- estimate
  out$se <- se
  out$t <- estimate / se
  out$p <- 2 * pt(-abs(out$t), df)
  out$df <- df
  out
}

# Rank-based inverse normal transform; missing values stay missing.
rank_normalise <- function(v) qnorm((rank(v, na.last = "keep") - 0.5) / sum(!is.na(v)))

# fit_robust for every column of y_mat (rows = samples). Missing responses are
# handled per column.
scan_robust <- function(y_mat, x, z, cluster = seq_len(nrow(y_mat)), cores = n_cores) {
  fit_one <- function(j) fit_robust(y_mat[, j], x, z, cluster)
  fits <- if (cores > 1) {
    parallel::mclapply(seq_len(ncol(y_mat)), fit_one, mc.cores = cores)
  } else {
    lapply(seq_len(ncol(y_mat)), fit_one)
  }
  out <- rbindlist(lapply(fits, as.data.table))
  out[, feature := colnames(y_mat)]
  setcolorder(out, "feature")
  out[]
}

# Cochran's Q for heterogeneity of independent estimates.
cochran_q <- function(estimate, se) {
  ok <- is.finite(estimate) & is.finite(se) & se > 0
  estimate <- estimate[ok]
  weight <- 1 / se[ok]^2
  pooled <- sum(weight * estimate) / sum(weight)
  q <- sum(weight * (estimate - pooled)^2)
  df <- length(estimate) - 1
  list(
    q = q, df = df,
    p = if (df > 0) pchisq(q, df, lower.tail = FALSE) else NA_real_,
    i2 = if (df > 0 && q > 0) max(0, (q - df) / q) else 0,
    pooled = pooled
  )
}
