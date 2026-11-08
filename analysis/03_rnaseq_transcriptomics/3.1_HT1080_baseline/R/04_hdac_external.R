hdac_external_source <- function(e) {
  cfg <- yaml::read_yaml('config/hdac12_external.yml'); set.seed(cfg$seed)
  meta <- as.data.frame(data.table::fread(cfg$sample_manifest))
  contrasts <- as.data.frame(data.table::fread(cfg$contrasts))
  raw <- as.data.frame(data.table::fread(cfg$counts))
  annotation <- raw[, 1:6]
  counts <- as.matrix(raw[, meta$sample_key]); rownames(counts) <- annotation$gene_id
  assert(nrow(counts) == 60683L && ncol(counts) == 69L && !anyDuplicated(rownames(counts)), 'External count dimensions/IDs changed')
  assert(!anyNA(counts) && all(counts >= 0 & counts == floor(counts)), 'External raw counts are not nonnegative integers')
  assert(setequal(names(raw)[-(1:6)], meta$sample_key), 'External count/sample mismatch')
  annotation$chrom <- sub(':.*', '', annotation$location)
  annotation$start <- as.integer(sub('^[^:]+:([0-9]+)-.*$', '\\1', annotation$location))
  annotation$end <- as.integer(sub('^[^:]+:[0-9]+-([0-9]+):.*$', '\\1', annotation$location))
  assert(!anyNA(annotation$start) && !anyNA(annotation$end), 'External coordinates failed to parse')
  annotation$PAR <- FALSE
  par <- data.table::fread('config/hg38_PAR.bed', header = FALSE)
  for (i in seq_len(nrow(par))) annotation$PAR <- annotation$PAR |
    (annotation$chrom == par[[1]][i] & annotation$start <= par[[3]][i] & annotation$end > par[[2]][i])
  annotation$normalization_background <- annotation$chrom %in% c(paste0('chr', 1:22), 'chrX') & !annotation$PAR
  annotation$versioned_gene_id <- paste(annotation$gene_id, annotation$gene_version, sep = '.')
  meta$library_reads <- colSums(counts); meta$detected_genes <- colSums(counts > 0)
  meta$mitochondrial_fraction <- colSums(counts[annotation$chrom == 'chrM', , drop = FALSE]) / meta$library_reads
  meta$top10_fraction <- apply(counts, 2, function(x) sum(sort(x, decreasing = TRUE)[1:10])) / meta$library_reads
  tab(meta, 'HDACEXT_library_QC')
  fits <- list(); all_results <- list(); normalization <- list(); pca <- list(); model_audit <- list()
  for (stratum in sort(unique(meta$sample_ID_stratum))) {
    message('External edgeR QL fit: ', stratum)
    m <- meta[meta$sample_ID_stratum == stratum & meta$used_in_model, ]
    m$model_group <- factor(m$model_group, levels = sort(unique(m$model_group)))
    design <- model.matrix(~ 0 + model_group, m); colnames(design) <- levels(m$model_group)
    rownames(design) <- m$sample_key
    assert(qr(design)$rank == ncol(design), 'Aliased external design')
    y <- edgeR::DGEList(counts[, m$sample_key], genes = annotation)
    keep <- edgeR::filterByExpr(y, design = design, min.count = cfg$filter_min_count,
      min.total.count = cfg$filter_min_total) & annotation$normalization_background
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::normLibSizes(y, method = 'TMM')
    y <- edgeR::estimateDisp(y, design = design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design = design, robust = TRUE, legacy = FALSE)
    exposure <- y$samples$lib.size * y$samples$norm.factors
    logcpm <- edgeR::cpm(y, log = TRUE, prior.count = cfg$prior_count)
    pc <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
    pca[[stratum]] <- data.frame(m, pc$x[, 1:3], PC1_variance = pc$sdev[1]^2 / sum(pc$sdev^2),
      PC2_variance = pc$sdev[2]^2 / sum(pc$sdev^2))
    normalization[[stratum]] <- data.frame(sample_key = m$sample_key, stratum = stratum,
      y$samples, effective_library_reads = exposure)
    model_audit[[stratum]] <- data.frame(stratum = stratum, libraries = nrow(m), groups = ncol(design),
      residual_df = nrow(design) - ncol(design), genes = nrow(y), common_dispersion = y$common.dispersion,
      model = 'no_intercept_source_group; robust_edgeR_QL; legacy_FALSE; no_replicate_number_pairing')
    tab(data.frame(sample_key = rownames(design), design, check.names = FALSE), paste0('HDACEXT_design_', stratum))
    cc <- contrasts[contrasts$sample_ID_stratum == stratum, ]
    for (i in seq_len(nrow(cc))) {
      con <- setNames(numeric(ncol(design)), colnames(design))
      con[cc$numerator[i]] <- 1; con[cc$denominator[i]] <- -1
      result <- edgeR::glmQLFTest(fit, contrast = con)$table
      result$external_FDR_nonzero <- p.adjust(result$PValue, method = 'BH')
      result$gene_id <- rownames(result); result$contrast <- cc$contrast[i]; result$stratum <- stratum
      all_results[[cc$contrast[i]]] <- result
    }
    fits[[stratum]] <- list(y = y, fit = fit, design = design, metadata = m, keep = keep,
      logCPM = logcpm, effective_library_reads = exposure)
    annotation[[paste0('expressed_', stratum)]] <- keep
  }
  results <- as.data.frame(data.table::rbindlist(all_results))
  results <- merge(results, annotation[, c('gene_id', 'versioned_gene_id', 'gene_name', 'gene_biotype', 'chrom')], by = 'gene_id')
  data.table::fwrite(results, 'outputs/tables/HDACEXT_external_gene_effects.csv.gz', compress = 'gzip')
  tab(annotation, 'HDACEXT_count_annotation_and_filters')
  normalization <- do.call(rbind, normalization); tab(normalization, 'HDACEXT_normalization')
  pca <- do.call(rbind, pca); tab(pca, 'HDACEXT_PCA_within_source_strata')
  model_audit <- do.call(rbind, model_audit); tab(model_audit, 'HDACEXT_model_audit')
  signatures <- list(); coverage <- list()
  for (name in contrasts$contrast) for (arm in c('UP', 'DN')) {
    z <- results[results$contrast == name & is.finite(results$external_FDR_nonzero) &
      results$external_FDR_nonzero < cfg$signature_FDR_nonzero &
      abs(results$logFC) >= cfg$signature_min_abs_logFC & (results$logFC > 0) == (arm == 'UP'), ]
    z <- z[order(-z$F, z$gene_id), ]; z$source_rank <- seq_len(nrow(z))
    z$selected_top250 <- z$source_rank <= cfg$signature_cap_per_arm
    z$source_arm <- rep(arm, nrow(z)); z$pathway <- rep(paste(name, arm, sep = '__'), nrow(z))
    signatures[[paste(name, arm)]] <- z
    coverage[[paste(name, arm)]] <- data.frame(pathway = paste(name, arm, sep = '__'), external_contrast = name,
      source_arm = arm, full_qualifying_genes = nrow(z), selected_source_genes = sum(z$selected_top250))
  }
  signatures <- as.data.frame(data.table::rbindlist(signatures)); coverage <- do.call(rbind, coverage)
  tab(signatures, 'HDACEXT_external_response_membership')
  # ID mapping is independent of external effect sign and never falls back to a symbol-only join.
  map <- e$mapping[, c('gene_id', 'stable_gene_id', 'gene_name', 'chrom')]
  names(map) <- c('HT1080_gene_id', 'gene_id', 'HT1080_symbol', 'HT1080_chrom')
  map <- merge(annotation, map, by = 'gene_id', all.x = TRUE, sort = FALSE)
  map$matched_HT1080_locus <- !is.na(map$HT1080_gene_id)
  map$chromosome_agrees <- map$matched_HT1080_locus & map$chrom == map$HT1080_chrom
  map$common_eligible <- map$expressed_S47 & map$expressed_S52 & map$chromosome_agrees
  map$symbol_changed <- map$matched_HT1080_locus & map$gene_name != map$HT1080_symbol
  tab(map, 'HDACEXT_cross_study_gene_mapping')
  common <- map[map$common_eligible, ]; common <- common[order(common$gene_id), ]
  assert(!anyDuplicated(common$gene_id) && !anyDuplicated(common$HT1080_symbol), 'Ambiguous cross-study mapping')
  tab(common, 'HDACEXT_common_gene_universe')
  coverage$represented_HT1080 <- vapply(coverage$pathway, function(n)
    sum(signatures$pathway == n & signatures$selected_top250 & signatures$gene_id %in% common$gene_id), integer(1))
  coverage$eligible <- coverage$represented_HT1080 >= cfg$pathway_min_size & coverage$represented_HT1080 <= cfg$pathway_max_size
  coverage$status <- ifelse(coverage$eligible, 'eligible', 'below_minimum_overlap')
  coverage <- merge(coverage, contrasts, by.x = 'external_contrast', by.y = 'contrast', all.x = TRUE)
  tab(coverage, 'HDACEXT_response_set_coverage')
  control_genes <- annotation$gene_id[grepl(cfg$HDAC_family_regex, annotation$gene_name)]
  controls <- results[results$gene_id %in% control_genes, ]; tab(controls, 'HDACEXT_HDAC_family_expression_context')
  # Reprocessing check against the published shrinkage/IHW output, without reusing it for selection.
  published <- as.data.frame(data.table::fread('outputs/tables/HDACEXT_published_DESeq2_results.csv.gz'))
  lookup <- c('24h JQ12' = 'JQ12_S52_vs_WT', '24h MS-275' = 'MS275_MS_S52_vs_WT',
    'WT HDAC1 CI' = 'HDAC1_CI_in_WT_vs_WT', 'KO HDAC1 CI' = 'HDAC1_CI_in_KO_vs_WT',
    'WT HDAC2 CI' = 'HDAC2_CI_in_WT_vs_WT', 'KO HDAC2 CI' = 'HDAC2_CI_in_KO_vs_WT',
    'WT HDAC3 CI' = 'HDAC3_CI_in_WT_vs_WT', 'KO HDAC3 CI' = 'HDAC3_CI_in_KO_vs_WT', 'HDAC3 KO' = 'HDAC3_KO_vs_WT')
  compatibility <- do.call(rbind, lapply(names(lookup), function(sheet) {
    a <- published[published$source_worksheet == sheet, ]; b <- results[results$contrast == lookup[[sheet]], ]
    joined <- merge(a, b, by = 'gene_id'); ok <- is.finite(joined$published_logFC) & is.finite(joined$logFC)
    data.frame(source_worksheet = sheet, contrast = lookup[[sheet]], genes = sum(ok),
      pearson = cor(joined$published_logFC[ok], joined$logFC[ok]),
      spearman = cor(joined$published_logFC[ok], joined$logFC[ok], method = 'spearman'),
      scope = 'compatibility_only; published_control_pooling_and_exact_software_differ; no_gene_independent_P')
  }))
  tab(compatibility, 'HDACEXT_published_effect_compatibility')
  list(cfg = cfg, metadata = meta, contrasts = contrasts, annotation = annotation, counts = counts,
    fits = fits, results = results, signatures = signatures, coverage = coverage, common = common,
    mapping = map, pca = pca, model_audit = model_audit, normalization = normalization, compatibility = compatibility)
}

hdac_external_compare <- function(external, x, e, p) {
  cfg <- external$cfg; common <- external$common; ids <- common$gene_id
  assert(length(ids) >= 1000L, 'Insufficient common expression background')
  ht <- cbind(p$effects_symbol[common$HT1080_symbol, , drop = FALSE],
    e$differences[common$HT1080_gene_id, , drop = FALSE], e$loo[common$HT1080_gene_id, , drop = FALSE])
  rownames(ht) <- ids
  ext <- sapply(external$contrasts$contrast, function(name) {
    z <- external$results[external$results$contrast == name, ]
    z$logFC[match(ids, z$gene_id)]
  })
  rownames(ext) <- ids; assert(!anyNA(ext) && !anyNA(ht), 'Cross-study effect alignment failed')
  selected <- external$signatures[external$signatures$selected_top250, ]
  sets <- lapply(external$coverage$pathway, function(name) intersect(selected$gene_id[selected$pathway == name], ids))
  names(sets) <- external$coverage$pathway
  core <- readRDS('outputs/models/04_HDAC12_gene_drivers.rds')
  legacy_mask <- common$gene_id[common$HT1080_symbol %in% unique(unlist(core$masks))]
  hdac_mask <- common$gene_id[grepl(cfg$HDAC_family_regex, common$gene_name) |
    grepl(cfg$HDAC_family_regex, common$HT1080_symbol)]
  masks <- list(original = character(), omit_IFN_ECM_cycle = legacy_mask, omit_HDAC_family = hdac_mask,
    omit_IFN_ECM_cycle_and_HDAC_family = union(legacy_mask, hdac_mask))
  assert(identical(names(masks), unlist(cfg$sensitivity_scenarios)), 'External sensitivity scenarios differ from configuration')
  mask_table <- do.call(rbind, lapply(names(masks), function(n)
    data.frame(scenario = n, gene_id = ids, HT1080_symbol = common$HT1080_symbol, retained = !ids %in% masks[[n]])))
  tab(mask_table, 'HDACEXT_sensitivity_universes')
  correlations <- list(); arm_means <- list(); enrichment <- list(); diagnostic <- list()
  original_camera <- readRDS('outputs/models/camera_voom_model.rds')
  for (scenario in names(masks)) {
    message('HT1080 comparison of external responses: ', scenario)
    keep_ids <- setdiff(ids, masks[[scenario]])
    retained <- lapply(sets, intersect, y = keep_ids)
    eligible <- lengths(retained) >= cfg$pathway_min_size & lengths(retained) <= cfg$pathway_max_size
    for (name in external$contrasts$contrast) {
      u <- retained[[paste0(name, '__UP')]]; d <- retained[[paste0(name, '__DN')]]
      evidence_genes <- intersect(external$signatures$gene_id[external$signatures$contrast == name], keep_ids)
      for (contrast in c('combined', 'AUG_2025', 'MAY_2026')) {
        a <- ext[keep_ids, name]; b <- ht[keep_ids, contrast]
        strong <- keep_ids %in% evidence_genes
        correlations[[paste(scenario, name, contrast)]] <- data.frame(scenario = scenario,
          external_contrast = name, contrast = contrast, genes = length(a), pearson = cor(a, b),
          spearman = cor(a, b, method = 'spearman'), nonzero_genes = sum(a != 0 & b != 0),
          direction_agreement = mean(sign(a[a != 0 & b != 0]) == sign(b[a != 0 & b != 0])),
          external_qualifying_genes = sum(strong),
          external_qualifying_direction_agreement = if (any(strong)) mean(sign(a[strong]) == sign(b[strong])) else NA_real_,
          external_qualifying_spearman = if (sum(strong) >= 3) cor(a[strong], b[strong], method = 'spearman') else NA_real_)
      }
      arm_mean <- function(g) if (length(g)) colMeans(ht[g, , drop = FALSE]) else rep(NA_real_, ncol(ht))
      um <- arm_mean(u); dm <- arm_mean(d)
      arm_means[[paste(scenario, name)]] <- data.frame(scenario = scenario, external_contrast = name,
        contrast = colnames(ht), UP_genes = length(u), DN_genes = length(d), UP_mean = um,
        DN_mean = dm, signed_mean_difference = um - dm,
        source_like_two_arm_pattern = um > 0 & dm < 0)
    }
    v <- original_camera$voom[match(common$HT1080_symbol[match(keep_ids, ids)], rownames(original_camera$voom$E)), ]
    rownames(v$E) <- rownames(v$weights) <- keep_ids
    indices <- limma::ids2indices(retained[eligible], keep_ids)
    for (ci in seq_along(c('combined', 'AUG_2025', 'MAY_2026'))) {
      contrast <- c('combined', 'AUG_2025', 'MAY_2026')[ci]
      set.seed(cfg$seed + 500 + 10 * match(scenario, names(masks)) + ci)
      stats <- ht[keep_ids, contrast]; names(stats) <- keep_ids
      stats <- stats[order(-stats, names(stats))]
      g <- withCallingHandlers(as.data.frame(fgsea::fgseaMultilevel(pathways = retained[eligible], stats = stats,
        minSize = cfg$pathway_min_size, maxSize = cfg$pathway_max_size,
        nPermSimple = cfg$fgsea_nPermSimple, eps = cfg$fgsea_eps, BPPARAM = BiocParallel::SerialParam())),
        warning = function(w) {
          diagnostic[[length(diagnostic) + 1L]] <<- data.frame(scenario = scenario,
            contrast = contrast, method = 'fgseaMultilevel', warning = conditionMessage(w))
          invokeRestart('muffleWarning')
        })
      g$leadingEdge <- vapply(g$leadingEdge, paste, collapse = ';', character(1))
      g <- g[, c('pathway', 'pval', 'NES', 'ES', 'size', 'log2err', 'leadingEdge')]
      row <- external$coverage
      row$scenario <- scenario; row$contrast <- contrast; row$universe_genes <- length(keep_ids)
      row$retained_genes <- lengths(retained)[row$pathway]
      row$eligible_now <- eligible[row$pathway]
      row$status_now <- ifelse(row$eligible_now, 'tested', 'below_minimum_overlap')
      row <- merge(row, g, by = 'pathway', all.x = TRUE)
      for (mode in c('fixed', 'estimated')) {
        rho <- if (mode == 'fixed') cfg$camera_inter_gene_cor else NA_real_
        cam <- limma::camera(v, indices, design = original_camera$design,
          contrast = original_camera$contrasts[, contrast], inter.gene.cor = rho,
          allow.neg.cor = FALSE, use.ranks = FALSE, trend.var = FALSE, sort = FALSE)
        cam$pathway <- rownames(cam)
        if (!'Correlation' %in% names(cam)) cam$Correlation <- rho
        cam <- cam[, c('pathway', 'NGenes', 'Direction', 'PValue', 'Correlation')]
        names(cam)[-1] <- paste0('camera_', mode, '_', c('genes', 'direction', 'P', 'rho'))
        row <- merge(row, cam, by = 'pathway', all.x = TRUE)
      }
      row$fgsea_secondary_BH <- row$camera_fixed_secondary_BH <- row$camera_estimated_secondary_BH <- NA_real_
      row$family_size <- if (scenario == 'original') sum(eligible) else NA_integer_
      if (scenario == 'original') {
        for (method in c('fgsea', 'camera_fixed', 'camera_estimated')) {
          pcol <- if (method == 'fgsea') 'pval' else paste0(method, '_P')
          j <- which(row$eligible_now)
          row[[paste0(method, '_secondary_BH')]][j] <- p.adjust(row[[pcol]][j], 'BH', n = length(j))
        }
      }
      row$camera_estimated_effective_rho <- pmax(0, row$camera_estimated_rho)
      row$fgsea_precision_floor <- !is.na(row$pval) & row$pval <= cfg$fgsea_eps
      row$statistical_scope <- if (scenario == 'original')
        'secondary_external_response_family; conditional_HT1080_camera_independent_residual_assumption' else
        'gene_exclusion_sensitivity; nominal_P_only; no_confirmatory_FDR'
      enrichment[[paste(scenario, contrast)]] <- row
    }
  }
  correlations <- do.call(rbind, correlations); arm_means <- do.call(rbind, arm_means)
  enrichment <- do.call(rbind, enrichment)
  tab(correlations, 'HDACEXT_cross_study_concordance')
  tab(arm_means, 'HDACEXT_HT1080_source_arm_effects')
  tab(enrichment, 'HDACEXT_fgsea_camera_comparison')
  diagnostic <- if (length(diagnostic)) do.call(rbind, diagnostic) else
    data.frame(scenario = character(), contrast = character(), method = character(), warning = character())
  tab(diagnostic, 'HDACEXT_numerical_warnings')
  # Long evidence preserves every common measurable locus, not just a leading edge or selected set.
  evidence <- do.call(rbind, lapply(external$contrasts$contrast, function(name) {
    response <- external$results[external$results$contrast == name, ]
    data.frame(external_contrast = name, common[, c('gene_id', 'versioned_gene_id', 'gene_name', 'HT1080_gene_id', 'HT1080_symbol', 'chrom')],
      external_logFC = ext[, name], external_FDR_nonzero = response$external_FDR_nonzero[match(ids, response$gene_id)],
      ht[, 1:3, drop = FALSE], source_UP = ids %in% sets[[paste0(name, '__UP')]],
      source_DN = ids %in% sets[[paste0(name, '__DN')]], in_existing_IFN_ECM_cycle = ids %in% legacy_mask,
      HDAC_family = ids %in% hdac_mask)
  }))
  data.table::fwrite(evidence, 'outputs/tables/HDACEXT_common_gene_evidence.csv.gz', compress = 'gzip')
  # Include missing external readouts explicitly; missing coverage is not a zero response.
  readouts <- expand.grid(external_contrast = external$contrasts$contrast, HT1080_symbol = cfg$readout_genes,
    stringsAsFactors = FALSE)
  readouts <- merge(readouts, evidence[evidence$HT1080_symbol %in% cfg$readout_genes, ],
    by = c('external_contrast', 'HT1080_symbol'), all.x = TRUE)
  readouts$external_coverage <- ifelse(is.na(readouts$external_logFC), 'not_in_common_expression_universe', 'represented')
  tab(readouts, 'HDACEXT_eight_readout_comparison')
  list(cfg = cfg, common = common, contrasts = external$contrasts, coverage = external$coverage,
    correlations = correlations, arm_means = arm_means, enrichment = enrichment,
    readouts = readouts, sets = sets, masks = masks, HT1080_effects = ht, external_effects = ext,
    numerical_warnings = diagnostic)
}

# Label audit only: a closer published effect profile does not verify dose/time metadata.
# All 36 S52 libraries enter this diagnostic TMM; none of its values define response sets.
hdac_external_ms_audit <- function(external) {
  m <- external$metadata[external$metadata$sample_ID_stratum == 'S52', ]
  y <- edgeR::DGEList(external$counts[external$fits$S52$keep, m$sample_key, drop = FALSE])
  y <- edgeR::normLibSizes(y, method = 'TMM')
  lc <- edgeR::cpm(y, log = TRUE, prior.count = external$cfg$prior_count)
  published <- as.data.frame(data.table::fread('outputs/tables/HDACEXT_published_DESeq2_results.csv.gz'))
  published <- published[published$source_worksheet == '24h MS-275', ]
  out <- do.call(rbind, lapply(c('HAP_MS', 'HAP_MS1', 'HAP_MS3'), function(group) {
    effect <- rowMeans(lc[, m$source_group == group, drop = FALSE]) -
      rowMeans(lc[, m$source_group == 'HAP', drop = FALSE])
    pub <- published$published_logFC[match(rownames(lc), published$gene_id)]
    ok <- is.finite(effect) & is.finite(pub)
    data.frame(source_group = group, stratum = 'S52', genes = sum(ok),
      pearson = cor(effect[ok], pub[ok]), spearman = cor(effect[ok], pub[ok], method = 'spearman'),
      interpretation = 'numerical_compatibility_only_not_a_verified_treatment_assignment',
      method = 'mean_logCPM_difference; diagnostic_TMM_all_36_S52_libraries; original_S52_expression_filter')
  }))
  tab(out, 'HDACEXT_MS_label_compatibility_audit')
  invisible(out)
}
