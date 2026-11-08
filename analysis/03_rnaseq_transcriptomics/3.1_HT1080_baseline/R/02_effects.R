normalize_counts <- function(counts, background, group, cfg) {
  # Filtering uses original all-feature library sizes. TMM is fitted on retained background.
  k <- edgeR::filterByExpr(counts, group = group, min.count = cfg$filter_min_count,
                         min.total.count = cfg$filter_min_total)
  keep <- k & background
  y <- edgeR::DGEList(counts[keep, , drop = FALSE])
  y <- edgeR::calcNormFactors(y, method = 'TMM')
  effective <- y$samples$lib.size * y$samples$norm.factors
  list(keep = keep, expression_keep = k, dge = y, effective = effective,
       logCPM = edgeR::cpm(counts, lib.size = effective, log = TRUE, prior.count = cfg$prior_count),
       CPM = edgeR::cpm(counts, lib.size = effective, log = FALSE))
}
expression_analysis <- function(x, cfg) {
  a <- x$annotation; m <- x$core; pairs <- x$pairs
  counts <- x$counts[, m$sample_key, drop = FALSE]
  group <- interaction(m$batch, m$treatment, drop = TRUE)
  n <- normalize_counts(counts, a$main_background, group, cfg)
  l <- n$logCPM; d <- matched_differences(l, pairs); pm <- period_average(d, pairs)
  e <- a; e$mean_logCPM <- rowMeans(l); e$mean_CPM <- rowMeans(n$CPM)
  e$expressed <- n$expression_keep; e$primary_eligible <- n$keep
  e$exclusion_reason <- ifelse(!a$main_background, 'separate_Y_PAR_MT_or_other_branch',
                               ifelse(!n$expression_keep, 'below_expression_filter', 'retained_primary'))
  e$AUG_2025 <- pm[, 'AUG_2025']; e$MAY_2026 <- pm[, 'MAY_2026']
  e$mean_difference <- rowMeans(pm)
  e$period_difference_MAY_minus_AUG <- e$MAY_2026 - e$AUG_2025
  e$period_direction_agrees <- sign(e$AUG_2025) == sign(e$MAY_2026)
  e$minimum_match <- apply(d, 1, min); e$maximum_match <- apply(d, 1, max)
  e$observed_match_SD <- apply(d, 1, stats::sd)
  e$matches_agreeing_with_mean <- rowSums(sign(d) == sign(e$mean_difference))
  e$median_match_difference <- apply(d, 1, median)
  # Deleting one timing match leaves each period at weight 1/2 (not 2/5 and 3/5).
  loo <- sapply(seq_len(ncol(d)), function(j) rowMeans(period_average(d[, -j, drop = FALSE], pairs[-j, ])))
  colnames(loo) <- paste0('omit_', pairs$block)
  e$LOO_min <- apply(loo, 1, min); e$LOO_max <- apply(loo, 1, max)
  e$LOO_max_abs_change <- apply(abs(loo - e$mean_difference), 1, max)
  e$LOO_all_same_direction <- rowSums(sign(loo) == sign(e$mean_difference)) == ncol(loo)
  low <- edgeR::cpm(counts, lib.size = n$effective, log = TRUE, prior.count = 0.25)
  e$prior025_mean_difference <- rowMeans(period_average(matched_differences(low, pairs), pairs))
  for (total in c(50, 100)) {
    e[[paste0('eligible_min_total_', total)]] <- a$main_background &
      edgeR::filterByExpr(counts, group = group, min.count = cfg$filter_min_count, min.total.count = total)
  }
  # A second normalization supplies a descriptive scale sensitivity, no DE tests.
  y_rle <- edgeR::calcNormFactors(edgeR::DGEList(counts[n$keep, ]), method = 'RLE')
  lr <- edgeR::cpm(counts, lib.size = y_rle$samples$lib.size * y_rle$samples$norm.factors,
                   log = TRUE, prior.count = cfg$prior_count)
  e$RLE_mean_difference <- rowMeans(period_average(matched_differences(lr, pairs), pairs))
  e$magnitude_category <- ifelse(abs(e$mean_difference) < cfg$effect_threshold, 'Below 0.5',
                                  ifelse(e$mean_difference > 0, 'Higher in LOY', 'Lower in LOY'))
  e$descriptive_candidate <- e$primary_eligible & abs(e$mean_difference) >= cfg$effect_threshold &
    e$period_direction_agrees & e$LOO_all_same_direction
  e <- cbind(e, as.data.frame(d), as.data.frame(loo))
  mapping <- choose_symbols(a, rowMeans(l), n$keep)
  tab(mapping, 'symbol_mapping_all_loci'); reps <- mapping[mapping$representative, ]
  tab(reps, 'symbol_mapping_representatives')
  tab(e, 'HT1080_LOY_vs_sgCTRL_all_gene_effects')
  tab(e[e$primary_eligible, ], 'HT1080_LOY_vs_sgCTRL_primary_effects')
  tab(e[e$descriptive_candidate, ], 'HT1080_LOY_vs_sgCTRL_descriptive_candidates')
  tab(data.frame(sample_key = m$sample_key, n$dge$samples,
                 effective_library_reads = n$effective, all_feature_reads = colSums(counts)), 'normalization_core')
  # Contextual PCA has its own all-15 normalization; it never changes the core effects.
  context <- normalize_counts(x$counts, a$main_background,
                               interaction(x$meta$batch, x$meta$treatment, drop = TRUE), cfg)
  pc <- prcomp(t(l[n$keep, ]), center = TRUE, scale. = FALSE)
  pc_context <- prcomp(t(context$logCPM[n$keep, ]), center = TRUE, scale. = FALSE)
  pca <- cbind(m, as.data.frame(pc$x[, 1:4])); tab(pca, 'PCA_core')
  pca_context <- cbind(x$meta, as.data.frame(pc_context$x[, 1:4])); tab(pca_context, 'PCA_context_all15')
  adjusted <- limma::removeBatchEffect(l[n$keep, ], batch = factor(m$batch),
                                      design = model.matrix(~ factor(m$treatment)))
  pc_adjusted <- prcomp(t(adjusted), center = TRUE, scale. = FALSE)
  tab(cbind(m, as.data.frame(pc_adjusted$x[, 1:4])), 'PCA_period_adjusted_visualization_only')
  corr <- cor(l[n$keep, ], method = 'pearson'); tab(matrix_long(corr, 'correlation', 'sample_1'), 'sample_correlations_core')
  qc <- x$meta; qc$assigned_reads <- colSums(x$counts)
  qc$detected_genes <- colSums(x$counts > 0)
  qc$assigned_fraction <- vapply(qc$sample_key, function(s) {
    z <- x$assignment[x$assignment$sample_key == s, ]; z$reads[z$status == 'Assigned'] / sum(z$reads)
  }, numeric(1))
  qc$mitochondrial_fraction <- colSums(x$counts[a$chrom == 'chrM', , drop = FALSE]) / qc$assigned_reads
  qc$ribosomal_protein_fraction <- colSums(x$counts[a$ribosomal_protein, , drop = FALSE]) / qc$assigned_reads
  qc$top10_fraction <- apply(x$counts, 2, function(v) sum(sort(v, decreasing = TRUE)[1:10])) / qc$assigned_reads
  qc$MSY_fraction <- colSums(x$counts[a$chrom == 'chrY' & !a$PAR, , drop = FALSE]) / qc$assigned_reads
  tab(qc, 'library_QC')
  # Parentals lack confirmed matching: unpaired group means, contextual scale only.
  ci <- x$meta$batch == 'AUG_2025'
  parent <- a
  for (tr in c('CTRL', 'LOY')) parent[[paste0(tr, '_minus_PARENTAL_context')]] <-
    rowMeans(context$logCPM[, ci & x$meta$treatment == tr, drop = FALSE]) -
    rowMeans(context$logCPM[, ci & x$meta$treatment == 'PARENTAL', drop = FALSE])
  tab(parent, 'parental_context_expression_differences')
  list(effects = e, norm = n, logCPM = l, differences = d, period_means = pm, loo = loo,
       mapping = reps, pca = pca, pca_variance = pc$sdev^2 / sum(pc$sdev^2),
       pca_context = pca_context, pca_context_variance = pc_context$sdev^2 / sum(pc_context$sdev^2),
       pca_adjusted = pc_adjusted, correlations = corr, qc = qc)
}
