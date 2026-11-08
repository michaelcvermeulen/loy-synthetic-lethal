# Stage 15 (plan Stage 6, part 2): which genes' expression is associated with each
# Y-loss dependency. Dependencies: pan-cancer hits, context hits, the X partners of
# the X-Y pairs, and a composite of the pan-cancer hits.
# For each dependency D and expressed gene g:
#   adjusted      D ~ lineage + library + Y status + expr_g   (removes the circular
#                 "Y genes track a Y-selected dependency" association)
#   within_neg    the same inside Y-negative lines only (modifiers once the Y is gone)
#   within_pos    the same inside Y-positive lines only
#   interaction   D ~ ... + Y status * expr_g
# Plus, for the seven core Y genes, which gene carries the association, using
# every male with RNA so that partial-deletion lines contribute.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "scan_helpers.R"))
assert_freeze()

fdr_cut <- params$gene_scan$fdr
core <- params$y_genes_core
cohort <- read_result("cohort")
gene_ann <- read_result("gene_annotation")
effect <- load_matrix("crispr_effect")
expr <- load_matrix("expr_tpm_pc")
pan <- read_result("gene_pan_results")
pan_hits <- read_result("gene_pan_hits")
context_hits <- read_result("context_gene_hits")
de <- read_result("expression_de_results")
base_covariates <- c("LibraryKY", "LibraryHumagneCD")

# --- dependency set ----------------------------------------------------------
partners <- gene_ann[in_scan == TRUE & !is.na(y_partner), feature]
candidates <- rbind(
  pan_hits[, .(feature, source = "pan-cancer hit", rank_p = p)],
  if (nrow(context_hits) > 0) context_hits[, .(feature, source = paste0("context hit: ", label), rank_p = interaction_p)],
  pan[feature %in% partners, .(feature, source = "X-Y partner", rank_p = p)]
)
candidates <- candidates[order(rank_p)][!duplicated(feature)]
keep <- candidates[feature %in% union(pan_hits$feature, partners) | seq_len(.N) <= params$expression$predictor_max_dependencies]
dependency_set <- head(keep, max(params$expression$predictor_max_dependencies, length(union(pan_hits$feature, partners))))
dependency_set[, symbol := gene_symbol(feature)]
dependency_set[, y_partner := gene_ann$y_partner[match(feature, gene_ann$feature)]]

lines_all <- cohort[ModelID %in% rownames(effect) & ModelID %in% rownames(expr)]
primary <- lines_all[in_primary == TRUE]
primary[, exposed := y_negative]
# Composite: mean z-score of the pan-cancer hits that are stronger dependencies in Y-negative lines.
composite_genes <- pan_hits[estimate < 0, feature]
z <- scale(effect[lines_all$ModelID, composite_genes, drop = FALSE])
dep_values <- cbind(effect[lines_all$ModelID, dependency_set$feature, drop = FALSE], COMPOSITE = rowMeans(z, na.rm = TRUE))
dependency_set <- rbind(dependency_set, data.table(feature = "COMPOSITE", source = sprintf("mean z of %d pan-cancer hits", length(composite_genes)), rank_p = NA_real_, symbol = "COMPOSITE", y_partner = NA_character_))
write_result(dependency_set, "expression_dependency_set")
log_msg("dependencies: %d (%s)", nrow(dependency_set), paste(dependency_set$symbol, collapse = ", "))

# --- expression universe (chosen without reference to Y status or dependency) --
expr_lines <- expr[lines_all$ModelID, , drop = FALSE]
expressed <- colMeans(expr_lines) >= 1 & apply(expr_lines, 2, sd) >= 0.25
expr_t <- t(expr_lines[, expressed, drop = FALSE])
gene_table <- data.table(gene = rownames(expr_t), gene_symbol = gene_symbol(rownames(expr_t)))
gene_table[, is_y_gene := gene %in% de[is_y_gene == TRUE, feature]]
gene_table[, de_label := de$label[match(gene, de$feature)]]
gene_table[, de_log_fc := de$log_fc[match(gene, de$feature)]]
log_msg("expression universe: %d genes of %d (mean log2(TPM+1) >= 1, SD >= 0.25); %d lines, %d in the primary call", nrow(expr_t), ncol(expr), nrow(lines_all), nrow(primary))

partial_r <- function(t, df) t / sqrt(t^2 + df)
plain_scan <- function(frame, d, covariates) {
  ok <- is.finite(d)
  frame <- frame[ok]
  design <- build_design(frame, tested = character(), covariates = covariates)
  frame <- frame[design$rows]
  d <- d[ok][design$rows]
  y <- expr_t[, frame$ModelID, drop = FALSE]
  t <- fast_t(y, d, design$design)
  df <- ncol(y) - qr(design$design)$rank - 1
  list(t = unname(t), r = partial_r(unname(t), df), p = 2 * pt(-abs(unname(t)), df), n = ncol(y))
}

# Interaction of Y status with expression in the forward model, vectorised over genes:
# D ~ Z + exposed + e + exposed:e, test of the last term by ordinary least squares.
interaction_scan <- function(frame, d) {
  ok <- is.finite(d)
  frame <- frame[ok]
  design <- build_design(frame, tested = character(), covariates = c(base_covariates, "exposed"))
  frame <- frame[design$rows]
  d <- d[ok][design$rows]
  decomposition <- qr(design$design)
  e <- expr_t[, frame$ModelID, drop = FALSE]
  e0 <- t(qr.resid(decomposition, t(e)))
  w0 <- t(qr.resid(decomposition, t(sweep(e, 2, frame$exposed, `*`))))
  d0 <- qr.resid(decomposition, d)
  s_ee <- rowSums(e0^2)
  s_ww <- rowSums(w0^2)
  s_ew <- rowSums(e0 * w0)
  s_ed <- as.numeric(e0 %*% d0)
  s_wd <- as.numeric(w0 %*% d0)
  det <- s_ee * s_ww - s_ew^2
  beta_e <- (s_ed * s_ww - s_ew * s_wd) / det
  beta_w <- (s_wd * s_ee - s_ew * s_ed) / det
  df <- length(d) - decomposition$rank - 2
  sigma2 <- (sum(d0^2) - beta_e * s_ed - beta_w * s_wd) / df
  t <- beta_w / sqrt(sigma2 * s_ee / det)
  list(estimate = beta_w, t = t, p = 2 * pt(-abs(t), df))
}

associations <- rbindlist(lapply(seq_len(nrow(dependency_set)), function(k) {
  feature <- dependency_set$feature[k]
  d_primary <- dep_values[primary$ModelID, feature]
  adjusted <- plain_scan(primary, d_primary, c(base_covariates, "exposed"))
  negative <- plain_scan(primary[exposed == 1], d_primary[primary$exposed == 1], base_covariates)
  positive <- plain_scan(primary[exposed == 0], d_primary[primary$exposed == 0], base_covariates)
  inter <- interaction_scan(primary, d_primary)
  out <- data.table(
    dependency = dependency_set$symbol[k], gene_table,
    adjusted_r = adjusted$r, adjusted_p = adjusted$p, n = adjusted$n,
    within_negative_r = negative$r, within_negative_p = negative$p, n_negative = negative$n,
    within_positive_r = positive$r, within_positive_p = positive$p, n_positive = positive$n,
    interaction = inter$estimate, interaction_p = inter$p
  )
  out[, adjusted_fdr := bh(adjusted_p)]
  out[, within_negative_fdr := bh(within_negative_p)]
  out[, within_positive_fdr := bh(within_positive_p)]
  out[, interaction_fdr := bh(interaction_p)]
  out
}))
saveRDS(associations, file.path(derived_dir, "expression_dependency_associations.rds"))
fwrite(associations, file.path(derived_dir, "expression_dependency_associations.csv.gz"))
top <- associations[, .SD[order(adjusted_p)][1:25], by = dependency]
write_result(top, "expression_dependency_top")
summary_tab <- associations[, .(
  genes = .N, adjusted_hits = sum(adjusted_fdr < fdr_cut), adjusted_hits_not_y = sum(adjusted_fdr < fdr_cut & !is_y_gene),
  within_negative_hits = sum(within_negative_fdr < fdr_cut), within_positive_hits = sum(within_positive_fdr < fdr_cut),
  interaction_hits = sum(interaction_fdr < fdr_cut),
  top_adjusted = gene_symbol[which.min(adjusted_p)], top_adjusted_r = adjusted_r[which.min(adjusted_p)],
  top_within_negative = gene_symbol[which.min(within_negative_p)], top_within_negative_r = within_negative_r[which.min(within_negative_p)]
), by = dependency]
write_result(summary_tab, "expression_dependency_summary")
print(summary_tab)

# --- which Y gene carries each dependency ------------------------------------
y_expr <- expr[lines_all$ModelID, match(core, gene_symbol(colnames(expr))), drop = FALSE]
colnames(y_expr) <- core
all_design <- build_design(lines_all, tested = character(), covariates = base_covariates)
stopifnot(length(all_design$rows) == nrow(lines_all))
z_all <- all_design$design
score <- 1 - lines_all$rna_y_score
positive_rna <- lines_all$rna_call %in% "expressed"
z_positive <- build_design(lines_all[positive_rna], tested = character(), covariates = base_covariates)$design
attribution <- rbindlist(lapply(seq_len(nrow(dependency_set)), function(k) {
  d <- dep_values[, dependency_set$feature[k]]
  score_fit <- fit_robust(d, score, z_all, cluster = lines_all$PatientID)
  rbindlist(lapply(core, function(g) {
    others <- y_expr[, setdiff(core, g), drop = FALSE]
    marginal <- fit_robust(d, y_expr[, g], z_all, cluster = lines_all$PatientID)
    joint <- fit_robust(d, y_expr[, g], cbind(z_all, others), cluster = lines_all$PatientID)
    in_positive <- fit_robust(d[positive_rna], y_expr[positive_rna, g], z_positive, cluster = lines_all$PatientID[positive_rna])
    given_gene <- fit_robust(d, score, cbind(z_all, y_expr[, g]), cluster = lines_all$PatientID)
    data.table(
      dependency = dependency_set$symbol[k], y_gene = g, paired = identical(dependency_set$y_partner[k], g),
      marginal_estimate = marginal$estimate, marginal_p = marginal$p,
      joint_estimate = joint$estimate, joint_p = joint$p,
      within_positive_estimate = in_positive$estimate, within_positive_p = in_positive$p,
      score_effect = score_fit$estimate, score_effect_given_gene = given_gene$estimate,
      attenuation = 1 - given_gene$estimate / score_fit$estimate, lines = marginal$n
    )
  }))
}))
write_result(attribution, "expression_dependency_y_genes")
print(attribution[dependency %in% pan_hits$symbol, .(dependency, y_gene, paired, marginal_p = signif(marginal_p, 2), joint_estimate = round(joint_estimate, 3), joint_p = signif(joint_p, 2), within_positive_p = signif(within_positive_p, 2), attenuation = round(attenuation, 2))])

# --- intersection with differential expression -------------------------------
# Competitive test of the DE gene sets against each dependency's association ranking.
de_sets <- list(
  up_in_y_negative = which(gene_table$gene %in% de[de == TRUE & log_fc > 0, feature]),
  down_in_y_negative = which(gene_table$gene %in% de[de == TRUE & log_fc < 0, feature]),
  loy_specific_de = which(gene_table$gene %in% de[label == "LOY-specific", feature])
)
de_sets <- de_sets[lengths(de_sets) >= 5]
overlap <- rbindlist(lapply(split(associations, by = "dependency"), function(part) {
  statistic <- qnorm(pmin(pmax(part$adjusted_p / 2, 1e-300), 0.5), lower.tail = FALSE) * sign(part$adjusted_r)
  res <- limma::cameraPR(statistic, de_sets)
  data.table(dependency = part$dependency[1], set = rownames(res), genes = res$NGenes, direction = res$Direction, p = res$PValue)
}))
overlap[, fdr := bh(p)]
write_result(overlap[order(p)], "expression_dependency_de_overlap")
print(overlap[order(p)][1:min(.N, 12)])

# Hallmark programmes associated with each dependency, beyond Y status.
hallmark <- rbindlist(lapply(split(associations[is_y_gene == FALSE], by = "dependency"), function(part) {
  part <- part[!duplicated(gene_symbol)]
  statistic <- setNames(qnorm(pmin(pmax(part$adjusted_p / 2, 1e-300), 0.5), lower.tail = FALSE) * sign(part$adjusted_r), part$gene_symbol)
  res <- run_fgsea(statistic, collections = "hallmark")
  res[, dependency := part$dependency[1]]
  res
}))
write_result(hallmark[order(pval)], "expression_dependency_hallmark")
log_msg("expression-dependency stage done: %d dependencies, %d genes", nrow(dependency_set), nrow(expr_t))
