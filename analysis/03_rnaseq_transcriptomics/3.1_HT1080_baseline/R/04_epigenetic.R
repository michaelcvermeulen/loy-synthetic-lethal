# Post hoc interpretation of frozen enrichment; no new gene-set testing family.
epigenetic_contrasts <- function(matched, pairs) {
  means <- period_average(matched, pairs)
  loo <- sapply(seq_len(nrow(pairs)), function(j)
    rowMeans(period_average(matched[, -j, drop = FALSE], pairs[-j, ])))
  colnames(loo) <- paste0('omit_', pairs$block)
  cbind(combined = rowMeans(means), means, matched, loo)
}

epigenetic_semantics <- function(z) {
  z$evidence_kind <- 'name_screen_only_not_individually_audited'
  response <- grepl('^(PRC2_|BMI1_|MEL18_|PRC1_BMI_|PASINI_|SENESE_|HELLER_|WANG_SMARCE1_|SCIBETTA_KDM5B_|JACKSON_DNMT1_|FOSTER_KDM1A_|KONDO_EZH2_)', z$gs_name)
  z$evidence_kind[response] <- 'source_RNA_perturbation_response'
  occupancy <- grepl('^(MIKKELSEN_|MEISSNER_|BENPORATH_|KONDO_.*H3K)', z$gs_name)
  z$evidence_kind[occupancy] <- 'source_occupancy_or_histone_mark'
  z$evidence_kind[grepl('^ACEVEDO_', z$gs_name)] <- 'source_differential_histone_mark'
  z$evidence_kind[grepl('^(REACTOME_|GOBP_|GOCC_|GOMF_|BIOCARTA_)', z$gs_name)] <- 'annotated_components_or_process'
  z$response_arm <- ifelse(response & grepl('_(UP|DN)$', z$gs_name), sub('.*_(UP|DN)$', '\\1', z$gs_name), NA_character_)
  z$direction_note <- 'See frozen description; name alone is not an activity score'
  z$direction_note[grepl('^PRC2_|^PRC1_BMI_', z$gs_name)] <- 'TIG3 regulator knockdown; EZH2_UP, SUZ12_UP and BMI_UP prefixes do not indicate overexpression'
  z$direction_note[grepl('^ACEVEDO_', z$gs_name)] <- 'UP/DN refer to higher/lower promoter histone methylation in liver tumors, not RNA direction'
  z$direction_note[grepl('^(MIKKELSEN_|MEISSNER_)', z$gs_name)] <- 'Mouse-source promoter states mapped to human symbols; no histone measurement in HT1080'
  z$direction_note[grepl('^HELLER_', z$gs_name)] <- 'Combined decitabine and trichostatin A response; not a selective HDAC or DNA methylation readout'
  z$direction_note[grepl('^(WANG_SMARCE1_|SCIBETTA_KDM5B_)', z$gs_name)] <- 'Ectopic expression response; target does not imply direct chromatin occupancy'
  z$direction_note[grepl('^FOSTER_', z$gs_name)] <- 'MSigDB description says heterozygous versus homozygous loss; source-study direction requires separate audit; not treated as a calibrated loss score'
  z$definition_url <- paste0('https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/', z$gs_name, '.html')
  z
}

epigenetic_deep_dive <- function(x, e, p) {
  cfg <- yaml::read_yaml('config/epigenetic_deep_dive.yml'); set.seed(cfg$seed)
  universe <- e$mapping$gene_name
  assert(identical(universe, rownames(p$effects_symbol)) && length(universe) == 25461L,
    'Epigenetic fixed representative universe changed')
  all_sets <- do.call(c, unname(p$sets))
  descriptions <- epigenetic_semantics(as.data.frame(data.table::fread('resources/epigenetic_deep_dive/descriptions.csv')))
  screen <- p$set_coverage[p$set_coverage$pathway %in% descriptions$gs_name, ]
  screen$focus <- grepl(paste(cfg$focus_patterns, collapse = '|'), screen$pathway)
  screen <- merge(screen, descriptions[, setdiff(names(descriptions), c('gs_id','gs_url','db_version'))],
    by.x = 'pathway', by.y = 'gs_name', all.x = TRUE, sort = FALSE)
  tab(screen, 'EPIGEN_catalogue')
  fg <- p$fgsea[p$fgsea$pathway %in% screen$pathway, ]
  camera <- p$camera[, c('pathway','contrast','Direction','PValue','BH_across_collections')]
  names(camera)[3:5] <- c('camera_direction','camera_P','camera_global_BH')
  sens <- p$camera_sensitivity[, c('pathway','contrast','Correlation','effective_correlation','PValue','BH_across_collections')]
  names(sens)[3:6] <- c('estimated_correlation','effective_correlation','camera_estimated_P','camera_estimated_global_BH')
  results <- merge(merge(fg, camera, by = c('pathway','contrast')), sens, by = c('pathway','contrast'))
  results$focus <- results$pathway %in% screen$pathway[screen$focus]
  results$evidence_kind <- screen$evidence_kind[match(results$pathway, screen$pathway)]
  tab(results, 'EPIGEN_original_enrichment')
  focus <- screen$pathway[screen$focus]
  assert(all(focus %in% names(all_sets)), 'Focused frozen membership absent')
  sets <- lapply(all_sets[focus], intersect, y = universe)
  ann <- e$effects[match(e$mapping$gene_id, e$effects$gene_id), ]
  assert(identical(ann$gene_name, universe), 'Epigenetic locus join mismatch')
  values <- cbind(p$effects_symbol, e$differences[e$mapping$gene_id, ], e$loo[e$mapping$gene_id, ])
  rownames(values) <- universe
  chromatin <- readRDS('outputs/models/04_chromatin_followup.rds')
  programs <- chromatin$programs
  program_union <- unique(unlist(programs, use.names = FALSE))
  assert(cfg$nucleosome_reference %in% names(all_sets), 'Nucleosome reference absent')
  nucleosome <- intersect(all_sets[[cfg$nucleosome_reference]], universe)
  evidence <- do.call(rbind, lapply(names(sets), function(pathway) {
    g <- sets[[pathway]]
    if (!length(g)) return(NULL)
    a <- ann[match(g, universe), c('gene_id','gene_name','chrom','arm','mean_CPM','mean_logCPM',
      'prior025_mean_difference','RLE_mean_difference','LOO_min','LOO_max')]
    a$pathway <- pathway
    for (contrast in c('combined','AUG_2025','MAY_2026')) {
      le <- fg$leadingEdge[fg$pathway == pathway & fg$contrast == contrast]
      a[[paste0('leading_edge_', contrast)]] <- if (length(le) == 1 && !is.na(le))
        g %in% strsplit(le, ';', fixed = TRUE)[[1]] else rep(NA, length(g))
    }
    for (name in names(programs)) a[[paste0('in_', name)]] <- g %in% programs[[name]]
    a$in_nucleosome <- g %in% nucleosome
    cbind(a, values[g, , drop = FALSE])
  }))
  tab(evidence, 'EPIGEN_complete_gene_evidence')
  summary <- function(g) hdac_summary(values, g, x$pairs)
  means <- do.call(rbind, lapply(names(sets), function(pathway) data.frame(pathway = pathway, summary(sets[[pathway]]))))
  tab(means, 'EPIGEN_complete_set_effects')
  normalization <- do.call(rbind, lapply(names(sets), function(pathway) {
    a <- ann[match(sets[[pathway]], universe), ]
    data.frame(pathway = pathway, genes = nrow(a), original = mean(a$mean_difference),
      prior025 = mean(a$prior025_mean_difference), RLE = mean(a$RLE_mean_difference))
  }))
  tab(normalization, 'EPIGEN_normalization_sensitivity')
  pair_names <- c(cfg$prc2_pairs, cfg$paired_comparators)
  assert(all(paste0(rep(pair_names, each = 2), c('_UP','_DN')) %in% names(all_sets)), 'Missing source-response pair')
  pairs <- do.call(rbind, lapply(pair_names, function(name) data.frame(source = name,
    chromatin_pair_means(values, all_sets, paste0(name, '_UP'), paste0(name, '_DN')))))
  pairs$interpretation <- ifelse(grepl('^FOSTER_', pairs$source), 'source_label_direction_requires_audit', 'source_RNA_arms_not_enzyme_activity')
  tab(pairs, 'EPIGEN_paired_arm_effects')
  ranks <- hdac_percentiles(x$counts[e$mapping$gene_id, x$core$sample_key, drop = FALSE]); rownames(ranks) <- universe
  score_list <- lapply(pair_names, function(name) {
    z <- hdac_signed_score(ranks, all_sets[[paste0(name, '_UP')]], all_sets[[paste0(name, '_DN')]])
    rownames(z) <- paste(name, rownames(z), sep = '__'); z
  })
  scores <- do.call(rbind, score_list)
  sample_scores <- merge(as.data.frame(matrix_long(scores, 'percentile_score', 'program')), x$core, by = 'sample_key')
  tab(sample_scores, 'EPIGEN_sample_scores')
  score_values <- epigenetic_contrasts(matched_differences(scores, x$pairs), x$pairs)
  score_contrasts <- as.data.frame(matrix_long(score_values, 'percentile_difference', 'program'))
  names(score_contrasts)[2] <- 'contrast'
  score_contrasts$source <- sub('__.*', '', score_contrasts$program)
  score_contrasts$score_arm <- sub('.*__', '', score_contrasts$program)
  tab(score_contrasts, 'EPIGEN_score_contrasts')
  tab(data.frame(gene_name = universe, gene_id = e$mapping$gene_id), 'EPIGEN_score_background')
  # Overlaps use complete represented sets, with no conditioning on an enrichment result.
  prc <- paste0(cfg$prc2_pairs, '_UP')
  prc_union <- unique(unlist(all_sets[prc], use.names = FALSE))
  prc_members <- data.frame(gene_name = prc_union, represented = prc_union %in% universe)
  for (name in cfg$prc2_pairs) for (arm in c('UP','DN'))
    prc_members[[paste(name, arm, sep = '_')]] <- prc_union %in% all_sets[[paste(name, arm, sep = '_')]]
  prc_members$UP_sources <- rowSums(prc_members[, prc])
  prc_members$any_DN_conflict <- rowSums(prc_members[, paste0(cfg$prc2_pairs, '_DN')]) > 0
  tab(prc_members, 'EPIGEN_PRC2_membership')
  eed <- sets[[cfg$primary_program]]
  mark <- unique(unlist(sets[cfg$mark_comparison], use.names = FALSE))
  partitions <- list(EED_only = setdiff(eed, mark), shared_EED_mark = intersect(eed, mark), mark_only = setdiff(mark, eed))
  part_members <- do.call(rbind, lapply(names(partitions), function(name)
    data.frame(partition = name, gene_name = partitions[[name]])))
  for (pathway in c(cfg$primary_program, cfg$mark_comparison)) part_members[[pathway]] <- part_members$gene_name %in% sets[[pathway]]
  tab(part_members, 'EPIGEN_EED_mark_membership')
  part_effects <- do.call(rbind, lapply(names(partitions), function(name) data.frame(partition = name, summary(partitions[[name]]))))
  tab(part_effects, 'EPIGEN_EED_mark_effects')
  shared_groups <- split(prc_members$gene_name[prc_members$represented],
    paste0('UP_in_', prc_members$UP_sources[prc_members$represented], '_sources'))
  support_effects <- do.call(rbind, lapply(names(shared_groups), function(name) data.frame(partition = name, summary(shared_groups[[name]]))))
  tab(support_effects, 'EPIGEN_PRC2_shared_effects')
  # Pairwise membership and leading-edge overlaps expose redundant evidence.
  overlap_sets <- sets[grepl('^PRC2_|^PASINI_|^BENPORATH_|^MIKKELSEN_(MEF|MCV6)_HCP_WITH_H3K27ME3$', names(sets))]
  overlap <- do.call(rbind, lapply(names(overlap_sets), function(a) do.call(rbind, lapply(names(overlap_sets), function(b) {
    ag <- overlap_sets[[a]]; bg <- overlap_sets[[b]]
    al <- evidence$gene_name[evidence$pathway == a & !is.na(evidence$leading_edge_combined) & evidence$leading_edge_combined]
    bl <- evidence$gene_name[evidence$pathway == b & !is.na(evidence$leading_edge_combined) & evidence$leading_edge_combined]
    data.frame(pathway_a = a, pathway_b = b, a_genes = length(ag), b_genes = length(bg),
      shared = length(intersect(ag, bg)), jaccard = length(intersect(ag, bg)) / length(union(ag, bg)),
      a_leading = length(al), b_leading = length(bl), leading_shared = length(intersect(al, bl)))
  }))))
  tab(overlap, 'EPIGEN_PRC2_mark_overlaps')
  influence <- list(); sensitivities <- list()
  for (pathway in names(sets)) {
    g <- sets[[pathway]]; if (!length(g)) next
    g <- g[order(-abs(values[g, 'combined']), ann$gene_id[match(g, universe)])]
    total <- sum(abs(values[g, 'combined']))
    influence[[pathway]] <- data.frame(pathway = pathway, gene_name = g, absolute_rank = seq_along(g),
      combined = values[g, 'combined'], contribution_to_mean = values[g, 'combined'] / length(g),
      absolute_effect_fraction = if (total > 0) abs(values[g, 'combined']) / total else NA_real_)
    mask <- list(original = character(), omit_IFN_ECM_cycle = program_union, omit_nucleosome = nucleosome,
      omit_top_5 = head(g, 5), omit_top_10 = head(g, 10),
      retain_CPM_at_least_5 = universe[ann$mean_CPM < 5])
    for (name in names(mask)) {
      keep <- setdiff(g, mask[[name]]); removed <- setdiff(g, keep)
      sensitivities[[paste(pathway, name)]] <- data.frame(pathway = pathway, scenario = name,
        original_genes = length(g), removed_genes = length(removed),
        removed_absolute_effect_fraction = if (total > 0) sum(abs(values[removed, 'combined'])) / total else NA_real_, summary(keep))
    }
  }
  influence <- do.call(rbind, influence); sensitivity <- do.call(rbind, sensitivities)
  tab(influence, 'EPIGEN_gene_influence'); tab(sensitivity, 'EPIGEN_descriptive_sensitivity')
  candidates <- cbind(ann, values[, setdiff(colnames(values), names(ann)), drop = FALSE])
  candidates$partition <- part_members$partition[match(candidates$gene_name, part_members$gene_name)]
  candidates <- candidates[!is.na(candidates$partition), ]
  candidates$pass_abundance <- candidates$mean_CPM >= 2
  candidates$pass_effect <- abs(candidates$combined) >= 0.5
  candidates$pass_period <- sign(candidates$combined) == sign(candidates$AUG_2025) & sign(candidates$combined) == sign(candidates$MAY_2026)
  candidates$pass_loo <- sign(candidates$combined) == sign(candidates$LOO_min) & sign(candidates$combined) == sign(candidates$LOO_max)
  candidates$pass_normalization <- sign(candidates$combined) == sign(candidates$prior025_mean_difference) & sign(candidates$combined) == sign(candidates$RLE_mean_difference)
  candidates$readout_eligible <- candidates$partition %in% c('EED_only','mark_only') & candidates$pass_abundance &
    candidates$pass_effect & candidates$pass_period & candidates$pass_loo & candidates$pass_normalization
  candidates$direction <- ifelse(candidates$combined > 0, 'higher', 'lower')
  candidates <- candidates[order(candidates$partition, candidates$direction, -abs(candidates$combined), candidates$gene_id), ]
  pick <- unlist(lapply(split(which(candidates$readout_eligible), paste(candidates$partition[candidates$readout_eligible],
    candidates$direction[candidates$readout_eligible])), head, n = 3), use.names = FALSE)
  candidates$selected_post_hoc <- seq_len(nrow(candidates)) %in% pick
  tab(candidates, 'EPIGEN_readout_candidates'); readouts <- candidates[candidates$selected_post_hoc, ]
  tab(readouts, 'EPIGEN_posthoc_readouts')
  atlas <- screen[grepl('^(MIKKELSEN|MEISSNER)_', screen$pathway), ]
  atlas$study <- paste(ifelse(grepl('^MIKKELSEN_', atlas$pathway), 'Mikkelsen', 'Meissner'),
    ifelse(atlas$gs_pmid == '17603471', '2007', '2008'))
  atlas$cell <- sub('^(MIKKELSEN|MEISSNER)_([^_]+)_.*', '\\2', atlas$pathway)
  atlas$CpG_class <- sub('.*_(HCP|ICP|LCP)_.*', '\\1', atlas$pathway)
  assert(all(atlas$CpG_class %in% c('HCP','ICP','LCP')), 'Unparsed promoter CpG class')
  atlas$mark <- ifelse(grepl('H3K4ME2', atlas$pathway),
    ifelse(grepl('H3K27ME3', atlas$pathway), 'K4me2 + K27me3', 'K4me2'),
    ifelse(grepl('H3K4ME3', atlas$pathway) & grepl('H3K27ME3|H327ME3', atlas$pathway), 'K4me3 + K27me3',
      ifelse(grepl('H3K4ME3', atlas$pathway), 'K4me3', 'K27me3')))
  atlas$row_label <- paste(atlas$study, atlas$cell, atlas$CpG_class, sep = ' / ')
  atlas <- merge(atlas, results[, c('pathway','contrast','NES','BH_across_collections','camera_global_BH','camera_estimated_global_BH')], by = 'pathway', all.x = TRUE)
  tab(atlas, 'EPIGEN_histone_context_atlas')
  regulators <- as.data.frame(data.table::rbindlist(lapply(cfg$regulators, function(gene) {
    z <- e$effects[e$effects$gene_name == gene, ]
    z <- z[order(-z$mean_logCPM, z$gene_id), ]
    if (!nrow(z)) return(data.frame(requested_gene = gene, gene_name = gene, status = 'not_mapped_as_symbol'))
    z <- z[1, c('gene_id','gene_name','chrom','primary_eligible','mean_CPM','mean_difference','AUG_2025','MAY_2026',
      'prior025_mean_difference','RLE_mean_difference','LOO_min','LOO_max')]
    z$requested_gene <- gene; z$status <- ifelse(z$primary_eligible, 'primary_RNA_context', 'nonprimary_RNA_context'); z
  }), fill = TRUE))
  tab(regulators, 'EPIGEN_regulator_RNA')
  list(cfg = cfg, catalogue = screen, results = results, sets = sets, evidence = evidence, means = means,
    pairs = pairs, sample_scores = sample_scores, score_contrasts = score_contrasts, score_values = score_values,
    prc_members = prc_members, partitions = partitions, part_members = part_members, part_effects = part_effects,
    support_effects = support_effects, overlap = overlap, influence = influence, sensitivity = sensitivity,
    candidates = candidates, readouts = readouts, atlas = atlas, regulators = regulators, values = values,
    normalization = normalization,
    programs = programs, nucleosome = nucleosome)
}
