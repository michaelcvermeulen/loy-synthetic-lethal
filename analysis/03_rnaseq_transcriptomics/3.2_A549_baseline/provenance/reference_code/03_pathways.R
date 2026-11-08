rank_effect <- function(v, symbols) {
  z <- setNames(v, symbols)
  z[order(-z, names(z))]
}
pathway_analysis <- function(x, e, cfg) {
  sets <- readRDS('resources/msigdb_all_collections.rds')
  all_sets <- do.call(c, unname(sets))
  metadata <- as.data.frame(data.table::fread('resources/msigdb_all_set_metadata.csv'))
  names(metadata)[match(c('gs_name', 'gs_collection', 'gs_subcollection'), names(metadata))] <-
    c('pathway', 'collection', 'subcollection')
  mapping <- e$mapping; ix <- match(mapping$gene_id, e$effects$gene_id)
  effects <- cbind(combined = e$effects$mean_difference[ix],
                   AUG_2025 = e$effects$AUG_2025[ix], MAY_2026 = e$effects$MAY_2026[ix])
  rownames(effects) <- mapping$gene_name
  tab(cbind(mapping[, c('gene_id', 'gene_name')], effects), 'pathway_TF_effect_ranking')
  coverage_sets <- merge(metadata, data.frame(pathway = names(all_sets),
    resource_size = lengths(all_sets), overlap_size = vapply(all_sets,
      function(g) length(intersect(g, rownames(effects))), integer(1))), by = 'pathway', sort = FALSE)
  coverage_sets$eligible <- coverage_sets$overlap_size >= cfg$pathway_min_size &
    coverage_sets$overlap_size <= cfg$pathway_max_size
  coverage_sets$status <- ifelse(coverage_sets$eligible, 'eligible',
    ifelse(coverage_sets$overlap_size < cfg$pathway_min_size, 'below_minimum_overlap', 'above_maximum_overlap'))
  tab(coverage_sets, 'pathway_resource_coverage')
  eligible <- all_sets[coverage_sets$pathway[coverage_sets$eligible]]
  out <- list()
  for (contrast in colnames(effects)) {
    set.seed(cfg$seed + match(contrast, colnames(effects)))
    message('fgsea: all ', length(sets), ' collections / ', contrast, ' (', length(eligible), ' eligible sets)')
    z <- as.data.frame(fgsea::fgseaMultilevel(pathways = eligible,
      stats = rank_effect(effects[, contrast], rownames(effects)), minSize = cfg$pathway_min_size,
      maxSize = cfg$pathway_max_size, eps = cfg$fgsea_eps, nPermSimple = cfg$fgsea_nPermSimple,
      BPPARAM = BiocParallel::SerialParam()))
    z$contrast <- contrast; z$padj <- NULL
    z$p_at_lower_bound <- !is.na(z$pval) & z$pval <= cfg$fgsea_eps
    z$leadingEdge <- vapply(z$leadingEdge, paste, character(1), collapse = ';')
    z$test_status <- ifelse(is.finite(z$pval), 'tested', 'numerical_failure')
    out[[contrast]] <- merge(z, metadata, by = 'pathway', sort = FALSE)
  }
  fg <- enrichment_adjust(do.call(rbind, out), 'pval'); rownames(fg) <- NULL
  fg$statistical_scope <- 'exploratory_preranked_gene_set_null_not_population_replication'
  tab(fg, 'fgsea_all_collections')
  # Preserve the older filename as an explicitly identified subset of the expanded test family.
  tab(fg[fg$collection == 'H' | fg$subcollection == 'CP:REACTOME', ], 'fgsea_hallmark_reactome_all')
  cam <- camera_analysis(x, e, eligible, metadata, cfg)
  overview <- enrichment_overview(fg, cam$results, coverage_sets)
  # Ribosomal-protein exclusion sensitivity on the same frozen representatives.
  keep_rp <- !mapping$ribosomal_protein
  set.seed(cfg$seed + 20)
  sens <- as.data.frame(fgsea::fgseaMultilevel(sets$H,
    rank_effect(effects[keep_rp, 'combined'], rownames(effects)[keep_rp]),
    minSize = cfg$pathway_min_size, maxSize = cfg$pathway_max_size, eps = cfg$fgsea_eps,
    nPermSimple = 10000, BPPARAM = BiocParallel::SerialParam()))
  sens$leadingEdge <- vapply(sens$leadingEdge, paste, character(1), collapse = ';')
  tab(sens, 'fgsea_hallmark_no_ribosomal_sensitivity')
  lc <- e$logCPM[mapping$gene_id, , drop = FALSE]; rownames(lc) <- mapping$gene_name
  # Within-sample ranks have a frozen gene background and no cross-sample fitting.
  percentiles <- apply(lc, 2, function(v) (rank(v, ties.method = 'average') - 1) / (length(v) - 1))
  hs <- t(vapply(sets$H, function(g) colMeans(percentiles[intersect(g, rownames(lc)), , drop = FALSE]),
                 numeric(ncol(lc))))
  colnames(hs) <- colnames(lc)
  tab(matrix_long(hs, 'mean_within_sample_percentile', 'pathway'), 'hallmark_sample_rank_scores')
  hd <- matched_differences(hs, x$pairs)
  tab(matrix_long(hd, 'LOY_minus_CTRL_rank_score', 'pathway'), 'hallmark_matched_rank_score_differences')
  # Freeze the local CollecTRI network, removing only contradictory duplicate edges.
  net <- as.data.frame(data.table::fread('resources/collectri_network.csv'))
  net <- unique(net[, c('source', 'target', 'mor')])
  keys <- paste(net$source, net$target, sep = '\r')
  conflict <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  tab(net[conflict, ], 'collectri_excluded_conflicting_edges')
  net <- net[!conflict & is.finite(net$mor), ]
  coverage <- aggregate(target ~ source, net, function(v) length(unique(v)))
  names(coverage)[2] <- 'resource_targets'
  cov2 <- aggregate(target ~ source, net[net$target %in% rownames(effects), ], function(v) length(unique(v)))
  names(cov2)[2] <- 'expressed_targets'
  coverage <- merge(coverage, cov2, all.x = TRUE); coverage$expressed_targets[is.na(coverage$expressed_targets)] <- 0
  coverage$eligible <- coverage$expressed_targets >= cfg$tf_min_targets
  tab(coverage, 'collectri_target_coverage')
  tf <- as.data.frame(decoupleR::run_ulm(effects, net, .source = 'source', .target = 'target',
                                         .mor = 'mor', minsize = cfg$tf_min_targets))
  names(tf)[names(tf) == 'condition'] <- 'contrast'
  tf$BH_target_pattern <- NA_real_
  for (contrast in unique(tf$contrast)) {
    j <- tf$contrast == contrast; tf$BH_target_pattern[j] <- p.adjust(tf$p_value[j], 'BH')
  }
  tf <- merge(tf, coverage, by = 'source', all.x = TRUE)
  tf$statistical_scope <- 'conditional_target_pattern_not_population_replication'
  tab(tf, 'collectri_ULM_effect_signatures')
  lc_centered <- lc - rowMeans(lc)
  tf_samples <- as.data.frame(decoupleR::run_ulm(lc_centered, net, .source = 'source', .target = 'target',
                                                 .mor = 'mor', minsize = cfg$tf_min_targets))
  tf_samples <- tf_samples[, c('source', 'condition', 'score')]
  names(tf_samples)[2] <- 'sample_key'; tab(tf_samples, 'collectri_ULM_sample_scores_centered')
  target_evidence <- merge(net, data.frame(target = rownames(effects), gene_id = mapping$gene_id, effects), by = 'target')
  target_evidence$signed_combined_effect <- target_evidence$mor * target_evidence$combined
  tab(target_evidence, 'collectri_target_effect_evidence')
  # PROGENy footprint weights, scaled by their absolute sum for transparent effect units.
  model <- as.matrix(readRDS('resources/progeny_model_top100.rds'))
  common <- intersect(rownames(model), rownames(lc)); weights <- model[common, , drop = FALSE]
  denom <- colSums(abs(weights)); assert(all(denom > 0), 'PROGENy footprint without overlap')
  weights <- sweep(weights, 2, denom, '/')
  prog <- t(lc_centered[common, , drop = FALSE]) %*% weights
  prog_d <- matched_differences(t(prog), x$pairs)
  tab(matrix_long(t(prog), 'weighted_centered_logCPM_score', 'pathway'), 'progeny_sample_footprint_scores')
  tab(matrix_long(prog_d, 'LOY_minus_CTRL_footprint_score', 'pathway'), 'progeny_matched_footprint_differences')
  tab(data.frame(pathway = colnames(weights), covered_targets = colSums(weights != 0),
                  model_targets = colSums(model != 0), absolute_weight_sum_before_scaling = denom), 'progeny_coverage')
  # Direct effect heatmap values for the leading edges; no row z-scores on effect differences.
  main <- fg[fg$contrast == 'combined', ]
  leading <- do.call(rbind, lapply(seq_len(nrow(main)), function(i) {
    g <- strsplit(main$leadingEdge[i], ';', fixed = TRUE)[[1]]
    if (!length(g) || identical(g, '')) return(NULL)
    data.frame(collection = main$collection[i], pathway = main$pathway[i], gene_name = g)
  }))
  leading <- merge(leading, data.frame(gene_name = mapping$gene_name, gene_id = mapping$gene_id, effects), by = 'gene_name')
  tab(leading, 'pathway_leading_edge_effects')
  list(fgsea = fg, tf = tf, tf_samples = tf_samples, targets = target_evidence, coverage = coverage,
       hallmark_scores = hs, hallmark_differences = hd, progeny_differences = prog_d,
       effects_symbol = effects, leading = leading, sets = sets, camera = cam$results,
       camera_sensitivity = cam$sensitivity, enrichment_summary = overview$summary,
       enrichment_comparison = overview$comparison, set_coverage = coverage_sets)
}
