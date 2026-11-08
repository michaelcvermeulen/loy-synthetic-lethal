# Stage 19 (plan Stage 8): robustness of every hit.
# A. Pan-cancer gene hits: alternative Y-status definitions, one model per donor,
#    leave-one-lineage-out, leave-one-line-out, covariate sets, permutation p.
# B. The pan-cancer scans repeated under each alternative definition, to see
#    whether a different definition changes the hit list (genes and drugs).
# C. Context hits (genes and drugs): permutation p, leave-one-line-out,
#    alternative definitions, co-occurring genomic states, replication in RNAi
#    or in another drug dataset.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
source(file.path(script_dir, "context_helpers.R"))
assert_freeze()

fdr_cut <- params$gene_scan$fdr
cohort <- read_result("cohort")
context <- read_result("context_matrix")
definitions <- read_result("context_definitions")
gene_ann <- read_result("gene_annotation")
effect <- load_matrix("crispr_effect")
rnai <- load_matrix("rnai")
pan_hits <- read_result("gene_pan_hits")
continuous <- read_result("gene_pan_continuous")
drug_pan <- read_result("drug_pan_results")
base_covariates <- c("LibraryKY", "LibraryHumagneCD")
exposures <- c(
  primary = "y_call", concordant_only = "y_call_concordant", rna_only = "y_call_rna", dna_only = "y_call_dna",
  one_copy_controls = "y_call_one_copy", no_relative_gain = "y_call_no_relative_gain", lenient_rna_rule = "y_call_lenient",
  xist_high_removed = "y_call_no_xist", threshold_stable_only = "y_call_stable"
)

frame_for <- function(mat, exposure = "y_call", one_per_donor = FALSE) {
  frame <- analysis_frame(cohort, mat, exposure)
  if (one_per_donor) frame <- frame[order(ModelID)][!duplicated(PatientID)]
  frame
}
single_fit <- function(values, frame, covariates = base_covariates, lineage = TRUE) {
  design <- build_design(frame, tested = "exposed", covariates = covariates, lineage = lineage)
  if (!design$estimable) return(list(estimate = NA_real_, se = NA_real_, p = NA_real_, n = 0L, n_exposed = NA_integer_, n_unexposed = NA_integer_))
  frame <- frame[design$rows]
  nuisance <- design$design[, colnames(design$design) != "exposed", drop = FALSE]
  fit_robust(values[frame$ModelID], frame$exposed, nuisance, cluster = frame$PatientID)
}
row_of <- function(fit, ...) data.table(..., estimate = fit$estimate, se = fit$se, p = fit$p, n = fit$n, n_negative = fit$n_exposed, n_positive = fit$n_unexposed)

# Permutation p for a main effect: Y status shuffled within lineage, plain t statistic.
permutation_main <- function(values, frame, covariates, n_perm = 20000) {
  design <- build_design(frame, tested = "exposed", covariates = covariates)
  frame <- frame[design$rows]
  y <- values[frame$ModelID]
  ok <- is.finite(y)
  frame <- frame[ok]
  nuisance <- design$design[ok, colnames(design$design) != "exposed", drop = FALSE]
  y_row <- matrix(y[ok], 1)
  observed <- fast_t(y_row, frame$exposed, nuisance)
  set.seed(params$seed + 19)
  null <- vapply(seq_len(n_perm), function(b) fast_t(y_row, permute_within(frame$exposed, frame$OncotreeLineage), nuisance), numeric(1))
  (1 + sum(abs(null) >= abs(observed))) / (n_perm + 1)
}

# --- A. pan-cancer gene hits ---------------------------------------------------
hit_features <- union(pan_hits$feature, continuous[fdr < fdr_cut, feature])
primary_frame <- frame_for(effect)
gene_robustness <- rbindlist(lapply(hit_features, function(feature) {
  values <- effect[, feature]
  symbol <- gene_symbol(feature)
  by_exposure <- rbindlist(lapply(names(exposures), function(name) row_of(single_fit(values, frame_for(effect, exposures[[name]])), symbol = symbol, check = "Y-status definition", variant = name)))
  donor <- row_of(single_fit(values, frame_for(effect, one_per_donor = TRUE)), symbol = symbol, check = "one model per donor", variant = "first model per donor")
  covariate_sets <- list(
    "lineage only" = character(), "lineage + library (primary)" = base_covariates,
    "+ WGD, TP53, aneuploidy" = c(base_covariates, "WGD", "tp53_altered", "Aneuploidy"),
    "+ LOH fraction, ploidy" = c(base_covariates, "WGD", "tp53_altered", "Aneuploidy", "LoHFraction", "Ploidy"),
    "+ screen quality" = c(base_covariates, "ScreenNNMD", "ScreenROCAUC", "ScreenMedianEssentialDepletion"),
    "+ Cas9 activity, doubling time" = c(base_covariates, "CasActivity", "ScreenDoublingTime")
  )
  by_covariates <- rbindlist(lapply(names(covariate_sets), function(name) row_of(single_fit(values, primary_frame, covariate_sets[[name]]), symbol = symbol, check = "covariates", variant = name)))
  disease_frame <- copy(primary_frame)[, OncotreeLineage := OncotreePrimaryDisease]
  by_disease <- row_of(single_fit(values, disease_frame), symbol = symbol, check = "covariates", variant = "primary disease instead of lineage")
  lineages <- primary_frame[, .N, by = OncotreeLineage][N >= 10, OncotreeLineage]
  loo_lineage <- rbindlist(lapply(lineages, function(l) row_of(single_fit(values, primary_frame[OncotreeLineage != l]), symbol = symbol, check = "leave one lineage out", variant = l)))
  within_lineage <- rbindlist(lapply(lineages, function(l) {
    sub <- primary_frame[OncotreeLineage == l]
    if (min(table(factor(sub$exposed, levels = 0:1))) < 3) return(NULL)
    row_of(single_fit(values, sub, lineage = FALSE), symbol = symbol, check = "within lineage", variant = l)
  }))
  loo_line <- vapply(seq_len(nrow(primary_frame)), function(i) single_fit(values, primary_frame[-i])$estimate, numeric(1))
  line_summary <- data.table(symbol = symbol, check = "leave one line out", variant = c("smallest |effect|", "largest |effect|"),
    estimate = c(loo_line[which.min(abs(loo_line))], loo_line[which.max(abs(loo_line))]), se = NA_real_, p = NA_real_, n = nrow(primary_frame) - 1L, n_negative = NA_integer_, n_positive = NA_integer_)
  permutation <- data.table(symbol = symbol, check = "permutation within lineage", variant = "20,000 shuffles", estimate = NA_real_, se = NA_real_,
    p = permutation_main(values, primary_frame, base_covariates), n = nrow(primary_frame), n_negative = NA_integer_, n_positive = NA_integer_)
  rbind(by_exposure, donor, by_covariates, by_disease, loo_lineage, within_lineage, line_summary, permutation)
}))
write_result(gene_robustness, "robustness_gene_hits")
print(dcast(gene_robustness[check == "Y-status definition"], variant ~ symbol, value.var = "estimate"))
print(gene_robustness[check %in% c("one model per donor", "covariates", "leave one line out", "permutation within lineage"), .(symbol, check, variant, estimate = round(estimate, 3), p = signif(p, 3), n)])

# --- B. the pan-cancer scans under each definition ------------------------------
scan_genes <- gene_ann[in_scan == TRUE, feature]
gene_alt <- rbindlist(lapply(names(exposures), function(name) {
  frame <- frame_for(effect, exposures[[name]])
  design <- build_design(frame, tested = "exposed", covariates = base_covariates)
  frame <- frame[design$rows]
  scan <- limma_scan(t(effect[frame$ModelID, scan_genes]), design$design, "exposed")
  scan[, fdr := bh(p)]
  hits <- scan[fdr < fdr_cut & abs(estimate) >= params$gene_scan$min_effect][order(p)]
  data.table(definition = name, lines = nrow(frame), y_negative = sum(frame$exposed == 1), y_positive = sum(frame$exposed == 0),
    hits = nrow(hits), genes = paste(gene_symbol(hits$feature), collapse = ";"), fdr_below_0.25 = sum(scan$fdr < 0.25))
}))
write_result(gene_alt, "robustness_gene_scan_by_definition")
print(gene_alt)
drug_alt <- rbindlist(lapply(unique(drug_pan$dataset), function(dataset_name) {
  mat <- load_matrix(dataset_name)
  features <- drug_pan[dataset == dataset_name, feature]
  rbindlist(lapply(names(exposures), function(name) {
    frame <- frame_for(mat, exposures[[name]])
    design <- build_design(frame, tested = "exposed")
    frame <- frame[design$rows]
    nuisance <- design$design[, colnames(design$design) != "exposed", drop = FALSE]
    scan <- scan_robust(mat[frame$ModelID, features, drop = FALSE], frame$exposed, nuisance, cluster = frame$PatientID)
    scan <- scan[n_exposed >= params$drug_scan$min_per_group & n_unexposed >= params$drug_scan$min_per_group & is.finite(p)]
    scan[, fdr := bh(p)]
    hit_names <- drug_pan[dataset == dataset_name][match(scan[fdr < fdr_cut, feature], feature), compound_name]
    data.table(dataset = dataset_name, definition = name, lines = nrow(frame), tested = nrow(scan), hits = sum(scan$fdr < fdr_cut), min_p = min(scan$p), min_fdr = min(scan$fdr),
      compounds = paste(hit_names, collapse = ";"))
  }))
}))
write_result(drug_alt, "robustness_drug_scan_by_definition")
print(dcast(drug_alt, dataset ~ definition, value.var = "hits"))

# --- C. context hits ------------------------------------------------------------
# Interaction test on a rank-normalised outcome, as in the screen, for one outcome.
interaction_fit <- function(values, frame, value, type, covariates) {
  ctx_frame <- with_context(frame, value[match(frame$ModelID, cohort$ModelID)], type)
  ctx_frame <- ctx_frame[is.finite(values[ModelID])]
  design <- build_design(ctx_frame, tested = "interaction", covariates = c("exposed", "ctx", covariates))
  if (!design$estimable) return(NULL)
  ctx_frame <- ctx_frame[design$rows]
  nuisance <- design$design[, colnames(design$design) != "interaction", drop = FALSE]
  y <- values[ctx_frame$ModelID]
  list(frame = ctx_frame, nuisance = nuisance, y = y, ranked = rank_normalise(y),
    fit = fit_ols(rank_normalise(y), ctx_frame$interaction, nuisance), raw = fit_ols(y, ctx_frame$interaction, nuisance))
}
# Y status shuffled within lineage and context level; the statistic is recomputed for every shuffle.
permutation_interaction <- function(prepared, type, n_perm = 20000) {
  frame <- prepared$frame
  base <- prepared$nuisance[, !colnames(prepared$nuisance) %in% "exposed", drop = FALSE]
  decomposition <- qr(base)
  d0 <- qr.resid(decomposition, prepared$ranked)
  strata <- if (type == "binary") paste(frame$OncotreeLineage, frame$ctx) else frame$OncotreeLineage
  statistic <- function(exposed_matrix) {
    e0 <- qr.resid(decomposition, exposed_matrix)
    w0 <- qr.resid(decomposition, exposed_matrix * frame$ctx)
    s_ee <- colSums(e0^2); s_ww <- colSums(w0^2); s_ew <- colSums(e0 * w0)
    s_ed <- as.numeric(crossprod(e0, d0)); s_wd <- as.numeric(crossprod(w0, d0))
    det <- s_ee * s_ww - s_ew^2
    beta_e <- (s_ed * s_ww - s_ew * s_wd) / det
    beta_w <- (s_wd * s_ee - s_ew * s_ed) / det
    df <- length(d0) - decomposition$rank - 2
    sigma2 <- (sum(d0^2) - beta_e * s_ed - beta_w * s_wd) / df
    beta_w / sqrt(sigma2 * s_ee / det)
  }
  observed <- statistic(matrix(frame$exposed, ncol = 1))
  set.seed(params$seed + 23)
  exceed <- 0
  for (chunk in split(seq_len(n_perm), ceiling(seq_len(n_perm) / 2000))) {
    shuffled <- vapply(chunk, function(b) permute_within(frame$exposed, strata), numeric(nrow(frame)))
    exceed <- exceed + sum(abs(statistic(shuffled)) >= abs(observed), na.rm = TRUE)
  }
  (1 + exceed) / (n_perm + 1)
}
check_context_hit <- function(values, mat_lines, ctx_id, replicate_values = NULL) {
  type <- definitions[context_id == ctx_id, type][1]
  value <- context[[ctx_id]][match(cohort$ModelID, context$ModelID)]
  frame <- cohort[in_primary == TRUE & ModelID %in% mat_lines][, exposed := y_negative]
  prepared <- interaction_fit(values, frame, value, type, base_covariates[base_covariates %in% names(frame)])
  if (is.null(prepared)) return(NULL)
  ctx_frame <- prepared$frame
  in_context <- if (type == "binary") which(ctx_frame$ctx == 1) else seq_len(nrow(ctx_frame))
  # Leave one line out: only lines in the context can move a binary interaction much.
  loo <- t(vapply(in_context, function(i) {
    fit <- fit_ols(rank_normalise(prepared$y[-i]), ctx_frame$interaction[-i], prepared$nuisance[-i, , drop = FALSE])
    c(fit$estimate, fit$p)
  }, numeric(2)))
  alt <- vapply(c("y_call_rna", "y_call_dna", "y_call_concordant"), function(exposure) {
    alt_frame <- cohort[ModelID %in% mat_lines & !is.na(get(exposure))][, exposed := as.numeric(get(exposure) == "Y_negative")]
    fit <- interaction_fit(values, alt_frame, value, type, base_covariates[base_covariates %in% names(alt_frame)])
    if (is.null(fit)) return(c(NA_real_, NA_real_))
    c(fit$fit$estimate, fit$fit$p)
  }, numeric(2))
  genomic <- interaction_fit(values, frame, value, type, c(base_covariates, "WGD", "tp53_altered", "Aneuploidy"))
  out <- data.table(
    lines = nrow(ctx_frame), interaction_raw = prepared$raw$estimate, interaction_rank = prepared$fit$estimate, interaction_p = prepared$fit$p,
    permutation_p = permutation_interaction(prepared, type),
    loo_max_p = max(loo[, 2]), loo_sign_stable = all(sign(loo[, 1]) == sign(prepared$fit$estimate)),
    rna_only_estimate = alt[1, "y_call_rna"], rna_only_p = alt[2, "y_call_rna"],
    dna_only_estimate = alt[1, "y_call_dna"], dna_only_p = alt[2, "y_call_dna"],
    concordant_estimate = alt[1, "y_call_concordant"], concordant_p = alt[2, "y_call_concordant"],
    genomic_adjusted_estimate = if (is.null(genomic)) NA_real_ else genomic$fit$estimate,
    genomic_adjusted_p = if (is.null(genomic)) NA_real_ else genomic$fit$p
  )
  if (!is.null(replicate_values)) {
    rep_frame <- cohort[in_primary == TRUE & ModelID %in% names(replicate_values)][, exposed := y_negative]
    rep_fit <- interaction_fit(replicate_values, rep_frame, value, type, character())
    cells <- if (is.null(rep_fit)) rep(NA, 4) else context_cells(rep_fit$frame$exposed, rep_fit$frame$ctx, "binary")
    out[, c("replication_estimate", "replication_p", "replication_min_cell") := list(
      if (is.null(rep_fit)) NA_real_ else rep_fit$fit$estimate, if (is.null(rep_fit)) NA_real_ else rep_fit$fit$p,
      if (type == "binary" && !is.null(rep_fit)) min(cells) else NA_real_
    )]
  }
  out
}

gene_context_hits <- read_result("context_gene_hits")
gene_context <- rbindlist(lapply(seq_len(nrow(gene_context_hits)), function(i) {
  hit <- gene_context_hits[i]
  symbol <- gene_symbol(hit$feature)
  replicate_values <- if (symbol %in% colnames(rnai)) rnai[, symbol] else NULL
  res <- check_context_hit(effect[, hit$feature], rownames(effect), hit$context_id, replicate_values)
  if (is.null(res)) return(NULL)
  cbind(data.table(family = hit$family, symbol = symbol, context = hit$label, class = hit$class, n_neg_in = hit$n_neg_in, n_pos_in = hit$n_pos_in), res)
}), fill = TRUE)
if (nrow(gene_context) > 0) {
  gene_context[, robust := permutation_p < 0.001 & loo_sign_stable & loo_max_p < 0.01 & sign(rna_only_estimate) == sign(interaction_rank) & sign(dna_only_estimate) == sign(interaction_rank)]
}
write_result(gene_context, "robustness_context_gene_hits")
print(gene_context)

drug_context_hits <- read_result("context_drug_hits")
drug_context <- rbindlist(lapply(seq_len(nrow(drug_context_hits)), function(i) {
  hit <- drug_context_hits[i]
  mat <- load_matrix(hit$dataset)
  # Replication: the same compound in a dataset from another provider, if there is one.
  others <- drug_pan[!is.na(compound_id) & compound_id %in% hit$compound_id & sub("_.*$", "", dataset) != sub("_.*$", "", hit$dataset)]
  replicate_values <- NULL
  replicate_dataset <- NA_character_
  if (nrow(others) > 0) {
    replicate_dataset <- others$dataset[1]
    replicate_values <- load_matrix(replicate_dataset)[, others$feature[1]]
  }
  res <- check_context_hit(mat[, hit$feature], rownames(mat), hit$context_id, replicate_values)
  if (is.null(res)) return(NULL)
  cbind(data.table(dataset = hit$dataset, compound = hit$compound_name, context = hit$label, class = hit$class, n_neg_in = hit$n_neg_in, n_pos_in = hit$n_pos_in, replicate_dataset = replicate_dataset), res)
}), fill = TRUE)
if (nrow(drug_context) > 0) {
  drug_context[, robust := permutation_p < 0.001 & loo_sign_stable & loo_max_p < 0.01 & sign(rna_only_estimate) == sign(interaction_rank) & sign(dna_only_estimate) == sign(interaction_rank)]
}
write_result(drug_context, "robustness_context_drug_hits")
print(drug_context)
log_msg("robustness done: %d pan-cancer gene hits, %d context gene hits (%d robust), %d context drug hits (%d robust)",
  length(hit_features), nrow(gene_context), if (nrow(gene_context) > 0) sum(gene_context$robust, na.rm = TRUE) else 0,
  nrow(drug_context), if (nrow(drug_context) > 0) sum(drug_context$robust, na.rm = TRUE) else 0)
