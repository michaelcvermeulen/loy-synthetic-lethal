# CAMERA is an explicitly exploratory fixed-model analysis of these populations.
# Passage blocking does not estimate or remove serial dependence within a population.
camera_design <- function(meta) {
  pair_indices(meta) # Validate one CTRL and one LOY per within-period timing block.
  assert(setequal(unique(meta$batch), c('AUG_2025', 'MAY_2026')), 'CAMERA requires both periods')
  block <- factor(meta$passage_day_block, levels = sort(unique(meta$passage_day_block)))
  design <- model.matrix(~ 0 + block)
  design <- cbind(design,
    LOY_AUG_2025 = as.numeric(meta$treatment == 'LOY' & meta$batch == 'AUG_2025'),
    LOY_MAY_2026 = as.numeric(meta$treatment == 'LOY' & meta$batch == 'MAY_2026'))
  rownames(design) <- meta$sample_key
  assert(!anyDuplicated(rownames(design)), 'Duplicate CAMERA samples')
  assert(qr(design)$rank == ncol(design) && nrow(design) > ncol(design), 'Non-estimable CAMERA model')
  contrasts <- matrix(0, ncol(design), 3,
    dimnames = list(colnames(design), c('combined', 'AUG_2025', 'MAY_2026')))
  contrasts[c('LOY_AUG_2025', 'LOY_MAY_2026'), 'combined'] <- 0.5
  contrasts['LOY_AUG_2025', 'AUG_2025'] <- 1
  contrasts['LOY_MAY_2026', 'MAY_2026'] <- 1
  list(design = design, contrasts = contrasts, residual_df = nrow(design) - ncol(design))
}

enrichment_adjust <- function(z, p_column, family_columns = 'contrast') {
  z$BH_across_collections <- z$BH_within_collection <- z$BH_within_subcollection <- NA_real_
  adjust <- function(columns) {
    groups <- do.call(paste, c(z[, columns, drop = FALSE], sep = '\r'))
    ans <- rep(NA_real_, nrow(z))
    for (j in split(seq_len(nrow(z)), groups)) ans[j] <- p.adjust(z[[p_column]][j], 'BH', n = length(j))
    ans
  }
  z$BH_across_collections <- adjust(family_columns)
  z$BH_within_collection <- adjust(c(family_columns, 'collection'))
  z$BH_within_subcollection <- adjust(c(family_columns, 'collection', 'subcollection'))
  z
}

camera_analysis <- function(x, e, sets, metadata, cfg) {
  model <- camera_design(x$core)
  counts <- x$counts[e$mapping$gene_id, x$core$sample_key, drop = FALSE]
  rownames(counts) <- e$mapping$gene_name
  # Preserve the normalization exposures estimated on the complete retained background.
  # Do not recalculate TMM or library totals on a gene-set-specific subset.
  effective <- setNames(e$norm$effective, colnames(e$logCPM))[colnames(counts)]
  assert(identical(colnames(counts), rownames(model$design)) && all(is.finite(effective)),
    'CAMERA sample or normalization misalignment')
  voom <- limma::voom(counts, design = model$design, lib.size = effective,
    normalize.method = 'none', plot = FALSE, save.plot = TRUE)
  index <- limma::ids2indices(sets, rownames(voom$E), remove.empty = FALSE)
  assert(all(lengths(index) >= cfg$pathway_min_size & lengths(index) <= cfg$pathway_max_size),
    'CAMERA and fgsea eligibility universes differ')
  tab(data.frame(sample_key = rownames(model$design), model$design, check.names = FALSE), 'camera_design_matrix')
  tab(data.frame(coefficient = rownames(model$contrasts), model$contrasts), 'camera_contrasts')
  tab(data.frame(gene_id = e$mapping$gene_id, symbol = rownames(counts)), 'camera_gene_universe')
  tab(data.frame(sample_key = colnames(counts), effective_library_reads = effective,
    mean_voom_weight = colMeans(voom$weights), min_voom_weight = apply(voom$weights, 2, min),
    max_voom_weight = apply(voom$weights, 2, max)), 'camera_voom_sample_weights')
  tab(data.frame(mean_log2_count = voom$voom.xy$x, sqrt_residual_SD = voom$voom.xy$y), 'camera_voom_mean_variance')
  tab(data.frame(mean_log2_count = voom$voom.line$x, fitted_sqrt_residual_SD = voom$voom.line$y), 'camera_voom_trend')
  out <- list()
  modes <- c('fixed_0.01', if (isTRUE(cfg$camera_estimated_correlation_sensitivity)) 'estimated_nonnegative')
  for (mode in modes) for (contrast in colnames(model$contrasts)) {
    message('CAMERA: ', mode, ' / ', contrast)
    rho <- if (mode == 'fixed_0.01') cfg$camera_inter_gene_cor else NA_real_
    z <- limma::camera(voom, index, design = model$design, contrast = model$contrasts[, contrast],
      inter.gene.cor = rho, allow.neg.cor = FALSE, use.ranks = FALSE, trend.var = FALSE, sort = FALSE)
    z$pathway <- rownames(z); rownames(z) <- NULL; z$FDR <- NULL
    z$contrast <- contrast; z$correlation_mode <- mode
    if (!'Correlation' %in% names(z)) z$Correlation <- rho
    z$effective_correlation <- pmax(0, z$Correlation)
    z$residual_df <- model$residual_df
    z$statistical_scope <- 'exploratory_fixed_passage_model_independent_residual_assumption_not_population_replication'
    out[[paste(mode, contrast)]] <- merge(z, metadata, by = 'pathway', sort = FALSE)
  }
  results <- as.data.frame(data.table::rbindlist(out, use.names = TRUE))
  results <- enrichment_adjust(results, 'PValue', c('correlation_mode', 'contrast'))
  results$test_status <- ifelse(is.finite(results$PValue), 'tested', 'numerical_failure')
  main <- results[results$correlation_mode == 'fixed_0.01', ]
  sensitivity <- results[results$correlation_mode != 'fixed_0.01', ]
  tab(main, 'camera_all_collections'); tab(sensitivity, 'camera_estimated_correlation_sensitivity')
  model$voom <- voom; model$gene_id <- e$mapping$gene_id; model$effective_library_reads <- effective
  saveRDS(model, 'outputs/models/camera_voom_model.rds')
  list(results = main, sensitivity = sensitivity, model = model)
}

enrichment_overview <- function(fg, cam, coverage) {
  keys <- c('pathway', 'collection', 'subcollection', 'contrast')
  z <- merge(fg[, c(keys, 'NES', 'pval', 'BH_across_collections')],
    cam[, c(keys, 'Direction', 'PValue', 'BH_across_collections')], by = keys,
    suffixes = c('_fgsea', '_camera'), all = TRUE)
  z$direction_agrees <- (z$NES > 0) == (z$Direction == 'Up')
  z$both_BH05 <- z$BH_across_collections_fgsea < 0.05 & z$BH_across_collections_camera < 0.05
  tab(z, 'fgsea_camera_comparison')
  summary <- do.call(rbind, lapply(split(z, paste(z$collection, z$contrast)), function(d) {
    cv <- coverage[coverage$collection == d$collection[1], ]
    data.frame(collection = d$collection[1], contrast = d$contrast[1],
      resource_sets = nrow(cv), eligible_sets = sum(cv$eligible),
      fgsea_finite = sum(is.finite(d$pval)), camera_finite = sum(is.finite(d$PValue)),
      fgsea_BH05 = sum(d$BH_across_collections_fgsea < 0.05, na.rm = TRUE),
      camera_BH05 = sum(d$BH_across_collections_camera < 0.05, na.rm = TRUE),
      both_BH05 = sum(d$both_BH05, na.rm = TRUE),
      both_BH05_same_direction = sum(d$both_BH05 & d$direction_agrees, na.rm = TRUE),
      direction_agreement_fraction = mean(d$direction_agrees, na.rm = TRUE))
  }))
  tab(summary, 'enrichment_collection_summary')
  list(comparison = z, summary = summary)
}
