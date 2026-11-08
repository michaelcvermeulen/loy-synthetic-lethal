chromatin_metadata <- function(metadata) {
  metadata$evidence_kind <- ifelse(grepl('^REACTOME_', metadata$gs_name), 'pathway_components',
    ifelse(grepl('^(MIKKELSEN_|BENPORATH_)', metadata$gs_name), 'occupancy_or_histone_mark_in_source_cells',
      'perturbation_response'))
  metadata$source_perturbation <- ifelse(grepl('knockdown', metadata$gs_description, ignore.case = TRUE), 'knockdown',
    ifelse(grepl('def[f]?icient|defficient', metadata$gs_description, ignore.case = TRUE), 'genetic_deficiency',
      ifelse(grepl('by expression|upon expression', metadata$gs_description, ignore.case = TRUE), 'ectopic_expression', 'not_a_perturbation')))
  metadata$response_arm <- ifelse(grepl('^Genes up-regulated', metadata$gs_description), 'UP',
    ifelse(grepl('^Genes down-regulated', metadata$gs_description), 'DN', NA_character_))
  metadata$source_pair <- ifelse(!is.na(metadata$response_arm), sub('_(UP|DN)$', '', metadata$gs_name), NA_character_)
  metadata$definition_url <- paste0('https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/', metadata$gs_name, '.html')
  metadata
}

chromatin_pair_means <- function(effects, sets, up_name, down_name) {
  # Contradictory membership is excluded from both arms, not assigned an arbitrary sign.
  up <- intersect(sets[[up_name]], rownames(effects)); down <- intersect(sets[[down_name]], rownames(effects))
  ambiguous <- intersect(up, down); up <- setdiff(up, ambiguous); down <- setdiff(down, ambiguous)
  assert(length(up) > 0 && length(down) > 0, 'Empty signed perturbation arm')
  u <- colMeans(effects[up, , drop = FALSE]); v <- colMeans(effects[down, , drop = FALSE])
  data.frame(contrast = colnames(effects), up_genes = length(up), down_genes = length(down),
    ambiguous_genes_excluded = length(ambiguous), UP_mean_difference = u, DN_mean_difference = v,
    source_response_contrast = u - v,
    two_arm_pattern = ifelse(u > 0 & v < 0, 'source_like_two_arm_pattern',
      ifelse(u < 0 & v > 0, 'opposite_source_two_arm_pattern',
        ifelse(u < 0 & v < 0, 'both_arms_lower', 'both_arms_higher_or_zero'))))
}

chromatin_followup <- function(x, e, p, cfg) {
  ccfg <- yaml::read_yaml('config/chromatin_followup.yml')
  all_sets <- do.call(c, unname(p$sets))
  descriptions <- as.data.frame(data.table::fread('resources/msigdb_chromatin_descriptions.csv'))
  screen <- p$set_coverage[grepl(ccfg$screen_pattern, p$set_coverage$pathway, ignore.case = TRUE), ]
  screen$focus <- grepl(paste(ccfg$focus_patterns, collapse = '|'), screen$pathway)
  tab(screen, 'chromatin_signature_screen_coverage')
  selected <- descriptions[descriptions$gs_name %in% screen$pathway[screen$focus], ]
  selected <- chromatin_metadata(selected)
  tab(selected, 'chromatin_signature_definitions')
  screen_results <- p$fgsea[p$fgsea$pathway %in% screen$pathway, ]
  cr <- p$camera[, c('pathway', 'contrast', 'Direction', 'PValue', 'BH_across_collections')]
  names(cr)[3:5] <- c('camera_direction', 'camera_P', 'camera_global_BH')
  cs <- p$camera_sensitivity[, c('pathway', 'contrast', 'Correlation', 'effective_correlation', 'PValue', 'BH_across_collections')]
  names(cs)[3:6] <- c('camera_estimated_rho', 'camera_effective_rho', 'camera_estimated_P', 'camera_estimated_global_BH')
  screen_results <- merge(merge(screen_results, cr, by = c('pathway', 'contrast')), cs, by = c('pathway', 'contrast'))
  screen_results$focus <- screen_results$pathway %in% selected$gs_name
  tab(screen_results, 'chromatin_signature_screen_results')
  focused <- screen_results[screen_results$focus, ]
  tab(focused, 'chromatin_signature_focused_results')
  universe <- rownames(p$effects_symbol)
  sets <- lapply(all_sets[selected$gs_name], intersect, y = universe)
  programs <- lapply(ccfg$overlap_programs, function(names) {
    assert(all(names %in% names(all_sets)), 'Overlap reference set absent')
    intersect(unique(unlist(all_sets[names])), universe)
  })
  tab(do.call(rbind, lapply(names(programs), function(n) data.frame(program = n, gene_name = programs[[n]]))),
    'chromatin_overlap_reference_membership')
  main <- focused[focused$contrast == 'combined', ]
  evidence <- do.call(rbind, lapply(names(sets), function(pathway) {
    if (!length(sets[[pathway]])) return(NULL)
    g <- sets[[pathway]]
    le <- main$leadingEdge[match(pathway, main$pathway)]
    leading <- if (length(le) && !is.na(le)) strsplit(le, ';', fixed = TRUE)[[1]] else character()
    data.frame(pathway = pathway, gene_name = g, in_leading_edge = g %in% leading,
      in_IFN = g %in% programs$IFN, in_ECM_EMT = g %in% programs$ECM_EMT,
      in_cell_cycle = g %in% programs$cell_cycle,
      gene_id = e$mapping$gene_id[match(g, e$mapping$gene_name)], p$effects_symbol[g, , drop = FALSE])
  }))
  annotation <- e$effects[, c('gene_id', 'chrom', 'arm', 'mean_CPM', 'LOO_min', 'LOO_max')]
  evidence <- merge(evidence, annotation, by = 'gene_id', all.x = TRUE)
  tab(evidence, 'chromatin_signature_gene_evidence')
  overlaps <- do.call(rbind, lapply(split(evidence, evidence$pathway), function(z) {
    do.call(rbind, lapply(names(programs), function(program) {
      j <- z$gene_name %in% programs[[program]]; le <- z$in_leading_edge
      data.frame(pathway = z$pathway[1], program = program, represented_genes = nrow(z),
        overlapping_genes = sum(j), membership_fraction = mean(j), leading_edge_genes = sum(le),
        overlapping_leading_genes = sum(j & le), leading_edge_fraction = if (any(le)) mean(j[le]) else NA_real_,
        leading_absolute_effect_fraction = if (sum(abs(z$combined[le])) > 0)
          sum(abs(z$combined[le & j])) / sum(abs(z$combined[le])) else NA_real_)
    }))
  }))
  tab(overlaps, 'chromatin_program_overlap')
  jaccard <- outer(names(sets), names(sets), Vectorize(function(a, b) {
    length(intersect(sets[[a]], sets[[b]])) / length(union(sets[[a]], sets[[b]]))
  }))
  dimnames(jaccard) <- list(names(sets), names(sets)); tab(matrix_long(jaccard, 'jaccard', 'pathway'), 'chromatin_set_overlap')
  chromosome <- do.call(rbind, lapply(split(evidence[evidence$in_leading_edge, ], evidence$pathway[evidence$in_leading_edge]), function(z) {
    t <- as.data.frame(table(z$chrom)); names(t) <- c('chrom', 'leading_genes'); t$pathway <- z$pathway[1]
    t$leading_fraction <- t$leading_genes / nrow(z); t
  }))
  tab(chromosome, 'chromatin_leading_edge_chromosome_distribution')
  # Inspect both response arms regardless of enrichment significance or set-size eligibility.
  pair_def <- selected[!is.na(selected$source_pair), ]
  matched <- e$differences[e$mapping$gene_id, , drop = FALSE]
  loo <- e$loo[e$mapping$gene_id, , drop = FALSE]
  rownames(matched) <- rownames(loo) <- universe
  pair_effects <- cbind(p$effects_symbol, matched, loo)
  pair_results <- do.call(rbind, lapply(split(pair_def, pair_def$source_pair), function(d) {
    assert(setequal(d$response_arm, c('UP', 'DN')), 'Missing source response arm')
    z <- chromatin_pair_means(pair_effects, sets, d$gs_name[d$response_arm == 'UP'], d$gs_name[d$response_arm == 'DN'])
    z$summary_type <- ifelse(z$contrast %in% colnames(p$effects_symbol), 'period_or_combined',
      ifelse(grepl('^omit_', z$contrast), 'leave_one_match_out', 'individual_passage_match'))
    z$source_pair <- d$source_pair[1]; z$source_perturbation <- d$source_perturbation[1]
    z$source_species <- d$gs_source_species[1]; z$source_pmid <- d$gs_pmid[1]
    z
  }))
  tab(pair_results, 'chromatin_paired_response_arms')
  regulators <- e$effects[e$effects$gene_name %in% ccfg$regulators, ]
  tab(regulators, 'chromatin_regulator_gene_effects')
  reg_match <- as.data.frame(matrix_long(e$differences[regulators$gene_id, , drop = FALSE], 'difference'))
  reg_match <- merge(reg_match, regulators[, c('gene_id', 'gene_name', 'chrom')], by = 'gene_id')
  tab(reg_match, 'chromatin_regulator_matched_differences')
  # Post hoc gene-universe exclusions: original global BH remains authoritative for the original screen.
  masks <- list(original = character(), without_MMP1_MMP3 = c('MMP1', 'MMP3'),
    without_IFN = programs$IFN, without_ECM_EMT = programs$ECM_EMT,
    without_IFN_ECM = union(programs$IFN, programs$ECM_EMT),
    without_IFN_ECM_cycle = unique(unlist(programs)))
  eligible_names <- screen$pathway[screen$focus & screen$eligible]
  model <- readRDS('outputs/models/camera_voom_model.rds')
  assert(identical(rownames(model$voom$E), universe), 'Chromatin CAMERA universe mismatch')
  sensitivity <- list(); masks_audit <- list()
  for (i in seq_along(masks)) {
    scenario <- names(masks)[i]; keep <- !universe %in% masks[[i]]
    represented <- lapply(sets[eligible_names], setdiff, y = masks[[i]])
    eligible_now <- lengths(represented) >= cfg$pathway_min_size & lengths(represented) <= cfg$pathway_max_size
    row <- data.frame(pathway = names(represented), scenario = scenario, universe_genes = sum(keep),
      retained_genes = lengths(represented), removed_genes = lengths(sets[names(represented)]) - lengths(represented),
      eligible = eligible_now, status = ifelse(eligible_now, 'tested', 'below_minimum_represented_genes'))
    row$retained_mean_difference <- vapply(represented, function(g) if (length(g)) mean(p$effects_symbol[g, 'combined']) else NA_real_, numeric(1))
    row$retained_mean_AUG <- vapply(represented, function(g) if (length(g)) mean(p$effects_symbol[g, 'AUG_2025']) else NA_real_, numeric(1))
    row$retained_mean_MAY <- vapply(represented, function(g) if (length(g)) mean(p$effects_symbol[g, 'MAY_2026']) else NA_real_, numeric(1))
    masks_audit[[scenario]] <- data.frame(scenario = scenario, gene_name = universe, retained = keep)
    if (scenario == 'original') {
      z <- main[, c('pathway', 'NES', 'pval', 'camera_direction', 'camera_P')]
    } else {
      message('Chromatin overlap sensitivity: ', scenario)
      set.seed(cfg$seed + 300 + i)
      g <- as.data.frame(fgsea::fgseaMultilevel(pathways = represented[eligible_now],
        stats = rank_effect(p$effects_symbol[keep, 'combined'], universe[keep]),
        minSize = cfg$pathway_min_size, maxSize = cfg$pathway_max_size, nPermSimple = cfg$fgsea_nPermSimple,
        eps = cfg$fgsea_eps, BPPARAM = BiocParallel::SerialParam()))
      v <- model$voom[keep, ]; indices <- limma::ids2indices(represented[eligible_now], rownames(v$E))
      cam <- limma::camera(v, indices, design = model$design, contrast = model$contrasts[, 'combined'],
        inter.gene.cor = cfg$camera_inter_gene_cor, allow.neg.cor = FALSE, trend.var = FALSE, sort = FALSE)
      cam <- data.frame(pathway = rownames(cam), camera_direction = cam$Direction, camera_P = cam$PValue)
      z <- merge(g[, c('pathway', 'NES', 'pval')], cam, by = 'pathway')
    }
    row <- merge(row, z, by = 'pathway', all.x = TRUE)
    row$statistical_scope <- 'post_hoc_overlap_sensitivity_nominal_P_no_global_BH_or_causal_adjustment'
    sensitivity[[scenario]] <- row
  }
  sensitivity <- do.call(rbind, sensitivity)
  tab(sensitivity, 'chromatin_overlap_exclusion_sensitivity')
  tab(do.call(rbind, masks_audit), 'chromatin_sensitivity_gene_universes')
  assert(!any(c('BH', 'padj', 'FDR') %in% names(sensitivity)), 'Focused diagnostic P values must not replace global BH')
  list(screen = screen, definitions = selected, results = focused, evidence = evidence,
    overlap = overlaps, pairs = pair_results, regulators = regulators, regulator_matches = reg_match,
    sensitivity = sensitivity, programs = programs, cfg = ccfg)
}
