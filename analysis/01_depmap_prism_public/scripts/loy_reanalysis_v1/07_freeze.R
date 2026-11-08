# Stage 07: freeze the parameters, the calls, the cohort and the context panel
# before any outcome is tested against Y status. Every later stage asserts
# these hashes. An existing freeze is never overwritten silently: pass
# "--replace" to archive it and write a new one.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))

replace <- "--replace" %in% commandArgs(trailingOnly = TRUE)
paths <- frozen_paths()
stopifnot(all(file.exists(paths)))
hashes <- as.list(vapply(paths, file_sha256, character(1)))

if (file.exists(freeze_file)) {
  previous <- jsonlite::fromJSON(freeze_file)
  if (identical(previous$sha256, hashes)) {
    log_msg("freeze already in place and unchanged")
    quit(status = 0)
  }
  if (!replace) stop("a different freeze already exists; rerun with --replace to archive it")
  archived <- file.path(derived_dir, sprintf("freeze_archived_%s.json", format(Sys.time(), "%Y%m%dT%H%M%S")))
  file.rename(freeze_file, archived)
  log_msg("previous freeze archived as %s", basename(archived))
}

matrices <- list.files(matrix_dir, pattern = "\\.rds$", full.names = TRUE)
matrix_manifest <- rbindlist(lapply(matrices, function(path) {
  object <- readRDS(path)
  data.table(matrix = sub("\\.rds$", "", basename(path)), rows = nrow(object), columns = ncol(object), sha256 = file_sha256(path))
}))
write_result(matrix_manifest, "matrix_manifest")

# Hypothesis ledger: how many tests each family will hold, counted before fitting.
cohort <- read_result("cohort")
primary <- cohort[in_primary == TRUE]
context <- read_result("context_matrix")
definitions <- read_result("context_definitions")[kept == TRUE]
gene_ann <- read_result("gene_annotation")
drug_ann <- read_result("drug_feature_annotation")
min_group <- params$drug_scan$min_per_group
min_cell <- params$context$min_cell
min_continuous <- params$context$min_group_continuous

count_features <- function(mat, lines, exposed, gate) {
  neg <- colSums(!is.na(mat[lines[exposed == 1], , drop = FALSE]))
  pos <- colSums(!is.na(mat[lines[exposed == 0], , drop = FALSE]))
  sum(neg >= gate & pos >= gate)
}
supported_contexts <- function(lines) {
  values <- context[match(lines$ModelID, ModelID)]
  vapply(seq_len(nrow(definitions)), function(i) {
    v <- values[[definitions$context_id[i]]]
    if (definitions$type[i] == "continuous") {
      return(sum(lines$y_negative == 1 & !is.na(v)) >= min_continuous && sum(lines$y_negative == 0 & !is.na(v)) >= min_continuous)
    }
    min(sum(lines$y_negative == 1 & v %in% 1), sum(lines$y_negative == 0 & v %in% 1), sum(lines$y_negative == 1 & v %in% 0), sum(lines$y_negative == 0 & v %in% 0)) >= min_cell
  }, logical(1))
}
ledger <- rbindlist(lapply(c("crispr_effect", unique(drug_ann$dataset)), function(name) {
  mat <- readRDS(file.path(matrix_dir, paste0(name, ".rds")))
  lines <- primary[ModelID %in% rownames(mat)]
  if (name == "crispr_effect") {
    main <- sum(gene_ann$in_scan)
    universe <- sum(gene_ann$selective)
  } else {
    active <- drug_ann[dataset == name & active == TRUE, feature]
    main <- count_features(mat[, intersect(colnames(mat), active), drop = FALSE], lines$ModelID, lines$y_negative, min_group)
    universe <- main
  }
  n_contexts <- sum(supported_contexts(lines))
  data.table(
    dataset = name, y_negative = sum(lines$y_negative == 1), y_positive = sum(lines$y_negative == 0),
    main_effect_tests = main, context_outcomes = universe, contexts_supported = n_contexts,
    interaction_tests_upper_bound = universe * n_contexts
  )
}))
write_result(ledger, "hypothesis_ledger")
print(ledger)

jsonlite::write_json(list(
  created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
  sha256 = hashes,
  primary_cohort = list(y_negative = sum(primary$y_negative == 1), y_positive = sum(primary$y_negative == 0)),
  note = "Thresholds were set from Y evidence only. During design, provisional calls were checked against the X-Y paralog CRISPR columns to confirm the positive controls are detectable; no threshold was changed afterwards."
), freeze_file, auto_unbox = TRUE, pretty = TRUE)
log_msg("freeze written: %s", freeze_file)
