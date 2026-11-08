zmynd8_followup <- function(x, e, p, cfg) {
  targets <- p$targets[p$targets$source == 'ZMYND8', ]
  assert(nrow(targets) >= cfg$tf_min_targets, 'Insufficient ZMYND8 regulon coverage')
  targets$paper_2016_focused_cooccupied_locus <- targets$target %in% c('SNAI2', 'CD44', 'VEGFA', 'EGFR')
  targets$combined_direction_supports_network_sign <- targets$signed_combined_effect > 0
  targets$share_absolute_signed_effect <- abs(targets$signed_combined_effect) / sum(abs(targets$signed_combined_effect))
  # This share describes effect magnitude, not an exact decomposition of the ULM t score.
  tab(targets, 'ZMYND8_regulon_target_evidence')
  shared <- p$targets[p$targets$target %in% c('MMP1', 'MMP3'), ]
  shared <- merge(shared, p$tf[, c('source', 'contrast', 'score', 'BH_target_pattern')], by = 'source')
  tab(shared, 'ZMYND8_MMP_shared_regulator_evidence')
  gene_names <- unique(c('ZMYND8', 'KDM5D', 'KDM5C', targets$target))
  genes <- e$effects[e$effects$gene_name %in% gene_names, ]
  tab(genes, 'ZMYND8_related_gene_effects')
  sample <- as.data.frame(matrix_long(e$logCPM[genes$gene_id, , drop = FALSE], 'log2_CPM'))
  sample$gene_name <- genes$gene_name[match(sample$gene_id, genes$gene_id)]
  sample$raw_reads <- as.numeric(x$counts[genes$gene_id, colnames(e$logCPM), drop = FALSE])
  sample <- merge(sample, x$core, by = 'sample_key')
  tab(sample, 'ZMYND8_related_sample_expression')
  network <- targets[, c('source', 'target', 'mor')]
  scenarios <- c(list(all_targets = character()), setNames(as.list(network$target), paste0('omit_', network$target)),
    list(omit_MMP1_MMP3 = c('MMP1', 'MMP3')))
  omission <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    net <- network[!network$target %in% scenarios[[scenario]], ]
    z <- as.data.frame(decoupleR::run_ulm(p$effects_symbol, net, minsize = cfg$tf_min_targets))
    names(z)[names(z) == 'condition'] <- 'contrast'
    z$scenario <- scenario; z$excluded_targets <- paste(scenarios[[scenario]], collapse = ';')
    z$remaining_targets <- nrow(net)
    z$statistical_scope <- 'post_hoc_driver_sensitivity_unadjusted_target_pattern_P_not_confirmatory'
    z
  }))
  tab(omission, 'ZMYND8_target_omission_sensitivity')
  matched <- e$differences[e$mapping$gene_id, , drop = FALSE]; rownames(matched) <- e$mapping$gene_name
  match_scores <- as.data.frame(decoupleR::run_ulm(matched, network, minsize = cfg$tf_min_targets))
  names(match_scores)[names(match_scores) == 'condition'] <- 'block'
  match_scores <- merge(match_scores[, c('source', 'block', 'score')], x$pairs, by = 'block')
  tab(match_scores, 'ZMYND8_matched_signature_scores')
  loo <- e$loo[e$mapping$gene_id, , drop = FALSE]; rownames(loo) <- e$mapping$gene_name
  loo_scores <- as.data.frame(decoupleR::run_ulm(loo, network, minsize = cfg$tf_min_targets))
  loo_scores <- loo_scores[, c('source', 'condition', 'score')]
  tab(loo_scores, 'ZMYND8_leave_one_match_out_scores')
  original <- p$tf[p$tf$source == 'ZMYND8', ]
  tab(original, 'ZMYND8_original_ULM_scores')
  heat <- as.data.frame(matrix_long(e$differences[targets$gene_id, , drop = FALSE], 'difference'))
  heat <- merge(heat, targets[, c('gene_id', 'target', 'mor', 'signed_combined_effect')], by = 'gene_id')
  heat$signed_difference <- heat$difference * heat$mor
  tab(heat, 'ZMYND8_target_matched_differences')
  list(targets = targets, genes = genes, sample = sample, omission = omission,
    matched = match_scores, loo = loo_scores, original = original, heat = heat, shared_regulators = shared)
}
