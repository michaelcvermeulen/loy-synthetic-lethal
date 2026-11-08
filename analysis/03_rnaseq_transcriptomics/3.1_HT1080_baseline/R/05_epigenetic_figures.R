epigenetic_heat <- function(z, row, column, value, legend, title, limit = max(abs(z[[value]]), na.rm = TRUE)) {
  z$text <- ifelse(is.na(z[[value]]), 'NA', sprintf('%+.2f', z[[value]]))
  heat_effect(z, row, column, value, limit, legend) +
    geom_text(aes(label = text), family = 'Arial', size = 5.5 / ggplot2::.pt) +
    labs(title = title) + theme(legend.position = 'bottom') +
    guides(fill = guide_colourbar(barwidth = grid::unit(30, 'mm'), barheight = grid::unit(2.5, 'mm')))
}

epigenetic_figures <- function(z, x) {
  periods <- c('AUG_2025','MAY_2026','combined')
  pair_labels <- c('PRC2_EED_DN.V1' = 'EED KD', 'PRC2_EZH2_UP.V1' = 'EZH2 KD', 'PRC2_SUZ12_UP.V1' = 'SUZ12 KD')
  arms <- z$pairs[z$pairs$source %in% names(pair_labels) & z$pairs$contrast %in% periods, ]
  arms <- rbind(transform(arms, source_arm = 'UP', value = UP_mean_difference, genes = up_genes),
    transform(arms, source_arm = 'DN', value = DN_mean_difference, genes = down_genes))
  arm_levels <- unlist(lapply(names(pair_labels), function(n) paste(pair_labels[n], c('UP','DN'))), use.names = FALSE)
  arms$label <- factor(paste(pair_labels[as.character(arms$source)], arms$source_arm), levels = rev(arm_levels))
  arms$period <- factor(period_labels[as.character(arms$contrast)], levels = unname(period_labels[periods]))
  a <- epigenetic_heat(arms, 'label', 'period', 'value', 'Mean log2 difference',
    'Complete source RNA arms (TIG3 knockdowns)')
  score <- z$score_contrasts[z$score_contrasts$source %in% names(pair_labels) &
    z$score_contrasts$score_arm == 'signed' & z$score_contrasts$contrast %in% x$pairs$block, ]
  score$occasion <- match(score$contrast, x$pairs$block)
  score$label <- unname(pair_labels[as.character(score$source)])
  b <- ggplot(score, aes(occasion, percentile_difference, colour = label, shape = label)) +
    geom_hline(yintercept = 0, colour = 'gray55', linewidth = 0.25) +
    geom_vline(xintercept = 3.5, linetype = 'dotted', colour = 'gray55', linewidth = 0.25) +
    geom_point(size = 1.65, position = position_dodge(width = 0.4)) +
    scale_colour_manual(values = c('EED KD' = '#0072B2', 'EZH2 KD' = '#D55E00', 'SUZ12 KD' = '#009E73'), name = NULL) +
    scale_shape_manual(values = c('EED KD' = 16, 'EZH2 KD' = 17, 'SUZ12 KD' = 15), name = NULL) +
    scale_x_continuous(breaks = 1:6, labels = c('A1','A2','A3','M1','M2','M3')) +
    labs(x = 'Matched passage (A: August; M: May)', y = 'LOY - sgCTRL signed percentile score',
      title = 'EED and EZH2 scores decrease in all six matches') +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  part <- z$part_effects[z$part_effects$contrast %in% periods, ]
  partition_labels <- c(EED_only = 'EED only (148)', shared_EED_mark = 'Shared (6)', mark_only = 'Mark only (497)')
  # Counts are derived; a test guards against a stale numerical display label.
  partition_labels <- setNames(paste(c('EED only','Shared','Mark only'),
    paste0('(', vapply(z$partitions, length, integer(1)), ')')), names(z$partitions))
  part$label <- factor(partition_labels[as.character(part$partition)], levels = rev(unname(partition_labels)))
  c <- ggplot(part, aes(mean_difference, label, colour = contrast, shape = contrast)) +
    geom_vline(xintercept = 0, colour = 'gray55', linewidth = 0.25) +
    geom_point(size = 1.8, position = position_dodge(width = 0.45)) + shape_periods() +
    scale_colour_manual(values = c(combined = 'black', AUG_2025 = '#0072B2', MAY_2026 = '#D55E00'), labels = period_labels, name = NULL) +
    labs(x = 'Mean LOY - sgCTRL log2 difference', y = NULL,
      title = 'EED response and source-mark genes mostly differ') +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  readout <- do.call(rbind, lapply(x$pairs$block, function(block) data.frame(
    gene_name = z$readouts$gene_name, partition = z$readouts$partition,
    gene_id = z$readouts$gene_id, block = block, value = z$values[z$readouts$gene_name, block])))
  readout$label <- factor(paste(readout$gene_name, ifelse(readout$partition == 'EED_only','(E)','(M)')),
    levels = rev(paste(z$readouts$gene_name, ifelse(z$readouts$partition == 'EED_only','(E)','(M)'))))
  readout$occasion <- factor(c('A1','A2','A3','M1','M2','M3')[match(readout$block, x$pairs$block)], levels = c('A1','A2','A3','M1','M2','M3'))
  d <- epigenetic_heat(readout, 'label', 'occasion', 'value', 'Log2 difference',
    'Post hoc readouts include opposite responses')
  id <- 'Fig10_PRC2_gene_programs'
  caption <- paste0('**Figure 10 | Distinct PRC2-related expression programs in established HT1080 populations.** ',
    '**a**, Complete represented UP and DN arms from EED, EZH2 and SUZ12 knockdowns in human TIG3 fibroblasts (Bracken et al., 2006). Values are mean matched log2 expression differences, with equal period weights for the combined estimate. ',
    '**b**, Matched changes in mean UP-minus-DN within-sample percentile scores on the unchanged 25,461-gene background. A1-A3 and M1-M3 denote August and May passages. Symbols identify source knockdowns; these are descriptive scores, not measured PRC2 activity. ',
    '**c**, Complete membership partitions of the EED-knockdown UP set and the union of the Mikkelsen MEF and MCV6 high-CpG H3K27me3 promoter sets; counts are shown in parentheses. These marks were measured in mouse source cells. ',
    '**d**, Post hoc examples from EED-only (E) and mark-only (M) partitions, including higher and lower genes. Selection requires mean CPM >= 2, absolute combined difference >= 0.5, concordant periods, all six omission estimates and both normalization sensitivities. Up to three genes per direction and partition are selected by absolute effect; only one higher EED-only gene qualifies. Values are unscaled matched effects. ',
    'Colors are centered on zero with panel-specific scales. All 12 libraries repeatedly sample one established population per condition from shared generation history. No gene-independent intervals, new enrichment tests or independent-population validation are implied.')
  manifests <- list()
  manifests[[id]] <- export_figure((a | b) / (c | d) + plot_layout(heights = c(0.9, 1.1)), id, caption,
    c(a = source_panel(id, 'a', arms), b = source_panel(id, 'b', score),
      c = source_panel(id, 'c', part), d = source_panel(id, 'd', readout)),
    'Three complete PRC2 knockdown pairs; six matches; exhaustive EED/mark partitions; balanced post hoc readout rule in frozen config', height = 200)
  # All available source promoter states, including size-ineligible sets and absent combinations.
  atlas <- z$atlas[z$atlas$contrast == 'combined' & !is.na(z$atlas$contrast) | is.na(z$atlas$contrast), ]
  rows <- sort(unique(atlas$row_label))
  mark_levels <- c('K27me3','K4me2','K4me3','K4me2 + K27me3','K4me3 + K27me3')
  grid <- expand.grid(row_label = rows, mark = mark_levels, stringsAsFactors = FALSE)
  atlas <- merge(grid, atlas, by = c('row_label','mark'), all.x = TRUE)
  assert(!anyDuplicated(paste(atlas$row_label, atlas$mark)), 'Duplicate source promoter-state cell')
  atlas$display_status <- ifelse(is.na(atlas$pathway), 'unavailable_combination', ifelse(atlas$eligible, 'tested', 'size_excluded'))
  atlas$text <- ifelse(atlas$display_status == 'unavailable_combination', '-',
    ifelse(atlas$display_status == 'size_excluded', 'NT', paste0(sprintf('%+.1f', atlas$NES),
      ifelse(atlas$BH_across_collections < 0.05, '*', ''))))
  atlas$row <- factor(atlas$row_label, levels = rev(rows)); atlas$column <- factor(atlas$mark, levels = mark_levels)
  a <- ggplot(atlas, aes(column, row)) +
    geom_tile(fill = 'white', colour = 'gray85', linewidth = 0.15) +
    geom_tile(data = atlas[atlas$display_status == 'size_excluded', ], fill = '#D0D0D0', colour = 'gray85', linewidth = 0.15) +
    geom_tile(data = atlas[atlas$display_status == 'tested', ], aes(fill = NES), colour = 'gray85', linewidth = 0.15) +
    geom_text(aes(label = text), family = 'Arial', size = 5.5 / ggplot2::.pt) +
    scale_fill_gradient2(low = '#0072B2', mid = 'white', high = '#D55E00', midpoint = 0,
      limits = c(-2, 2), name = 'Original combined NES') +
    scale_x_discrete(labels = c('K27me3','K4me2','K4me3','K4me2\n+ K27me3','K4me3\n+ K27me3')) +
    labs(x = NULL, y = NULL, title = 'Promoter-state context in mouse source studies') +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 5.5),
      axis.text.x = element_text(size = 5.5), legend.position = 'bottom') +
    guides(fill = guide_colourbar(barwidth = grid::unit(30, 'mm'), barheight = grid::unit(2.5, 'mm')))
  comparator_labels <- c(PASINI_SUZ12_TARGETS = 'SUZ12 loss / mouse ES', BMI1_DN.V1 = 'BMI1 KD / DAOY',
    MEL18_DN.V1 = 'MEL18 KD / DAOY', BMI1_DN_MEL18_DN.V1 = 'BMI1 + MEL18 KD / DAOY',
    PRC1_BMI_UP.V1 = 'BMI1 KD / TIG3', WANG_SMARCE1_TARGETS = 'SMARCE1 expression / BT549',
    SCIBETTA_KDM5B_TARGETS = 'KDM5B expression / HB2', JACKSON_DNMT1_TARGETS = 'DNMT1 loss / mouse MEF',
    FOSTER_KDM1A_TARGETS = 'KDM1A source arms* / mouse ES', SENESE_HDAC1_TARGETS = 'HDAC1 KD / U2OS',
    SENESE_HDAC2_TARGETS = 'HDAC2 KD / U2OS', SENESE_HDAC1_AND_HDAC2_TARGETS = 'HDAC1 + HDAC2 KD / U2OS',
    HELLER_HDAC_TARGETS_SILENCED_BY_METHYLATION = 'Decitabine + TSA / myeloma')
  comp <- do.call(rbind, lapply(names(comparator_labels), function(name) {
    s <- chromatin_pair_means(z$values, z$sets, paste0(name, '_UP'), paste0(name, '_DN'))
    s <- s[s$contrast == 'combined', ]
    data.frame(source = name, source_arm = c('UP','DN'), value = c(s$UP_mean_difference, s$DN_mean_difference),
      genes = c(s$up_genes, s$down_genes))
  }))
  comp$label <- factor(comparator_labels[as.character(comp$source)], levels = rev(unname(comparator_labels)))
  b <- ggplot(comp, aes(value, label, colour = source_arm, shape = source_arm)) +
    geom_vline(xintercept = 0, colour = 'gray55', linewidth = 0.25) +
    geom_point(size = 1.5, position = position_dodge(width = 0.45)) +
    scale_colour_manual(values = c(UP = '#D55E00', DN = '#0072B2'), name = 'Source RNA arm') +
    scale_shape_manual(values = c(UP = 17, DN = 16), name = 'Source RNA arm') +
    labs(x = 'Mean LOY - sgCTRL log2 difference', y = NULL,
      title = 'Other source-response programs\n*Direction caveat in source audit') + theme(axis.text.y = element_text(size = 5.5))
  diagnostic_names <- c(PRC2_EED_DN.V1_UP = 'EED KD UP', PRC2_EZH2_UP.V1_UP = 'EZH2 KD UP',
    PRC2_SUZ12_UP.V1_UP = 'SUZ12 KD UP', MIKKELSEN_MEF_HCP_WITH_H3K27ME3 = 'MEF HCP K27me3',
    MIKKELSEN_MCV6_HCP_WITH_H3K27ME3 = 'MCV6 HCP K27me3')
  diagnostic <- z$sensitivity[z$sensitivity$pathway %in% names(diagnostic_names) & z$sensitivity$contrast == 'combined', ]
  scenario_labels <- c(original = 'All', omit_IFN_ECM_cycle = 'Omit\nprograms', omit_nucleosome = 'Omit\nnucleosome',
    omit_top_5 = 'Omit\ntop 5', omit_top_10 = 'Omit\ntop 10', retain_CPM_at_least_5 = 'CPM\n>= 5')
  diagnostic$row <- factor(diagnostic_names[as.character(diagnostic$pathway)], levels = rev(unname(diagnostic_names)))
  diagnostic$column <- factor(scenario_labels[as.character(diagnostic$scenario)], levels = unname(scenario_labels))
  c <- epigenetic_heat(diagnostic, 'row', 'column', 'mean_difference', 'Mean log2 difference',
    'Direction survives fixed descriptive exclusions') + theme(axis.text.x = element_text(size = 5.5), axis.text.y = element_text(size = 5.5))
  rgenes <- c('EED','EZH1','EZH2','SUZ12','PHF19','KDM6A','KDM6B','KDM5B','TET2','SETDB2','HDAC1','HDAC2')
  reg <- do.call(rbind, lapply(periods, function(contrast) data.frame(gene_name = rgenes, contrast = contrast,
    value = z$regulators[[ifelse(contrast == 'combined','mean_difference',contrast)]][match(rgenes, z$regulators$gene_name)])))
  reg$row <- factor(reg$gene_name, levels = rev(rgenes))
  d <- ggplot(reg, aes(value, row, colour = contrast, shape = contrast)) +
    geom_vline(xintercept = 0, colour = 'gray55', linewidth = 0.25) +
    geom_point(size = 1.45, position = position_dodge(width = 0.5)) + shape_periods() +
    scale_colour_manual(values = c(combined = 'black', AUG_2025 = '#0072B2', MAY_2026 = '#D55E00'), labels = period_labels, name = NULL) +
    labs(x = 'LOY - sgCTRL log2 difference', y = NULL, title = 'Regulator RNA is not a mark or activity assay') +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  id <- 'Fig11_epigenetic_context'
  caption <- paste0('**Figure 11 | Epigenetic enrichments depend on source context and response direction.** ',
    '**a**, All available Mikkelsen/Meissner promoter-state combinations in the 339-set screen, organized by study, source cell and high/intermediate/low CpG class (HCP/ICP/LCP). Values are original combined fgsea normalized enrichment scores (NES); * marks original global BH < 0.05 across 23,696 sets. Gray NT denotes overlap outside 15-500 genes; a dash denotes an unavailable combination. H3 prefixes are shortened to K. MCV6 is partially reprogrammed; IPS denotes source induced pluripotent cells. ',
    '**b**, Complete UP and DN source RNA-arm means. KD denotes knockdown. The KDM1A publication describes loss responses, whereas its MSigDB comparison wording reverses the genotype order; labels are preserved without activity calibration. ',
    '**c**, Mean combined effects after removing the frozen interferon/matrix/cycle union, nucleosome genes, the five or ten largest absolute effects, or genes with mean CPM < 5. Removals are fixed across periods; these are descriptive means, not recalculated enrichment tests. ',
    '**d**, Selected regulator transcripts across periods and their equally weighted mean; the complete regulator table accompanies the figure. ',
    'Heatmaps use separate zero-centered scales. Source chromatin states and RNA responses do not establish HT1080 H3K27me3, H3K9me3, DNA methylation, complex abundance or enzyme activity. Conditional enrichment statistics remain exploratory because the 12 HT1080 libraries repeatedly sample the same established populations; estimated-correlation CAMERA yields no combined global BH < 0.05 result.')
  manifests[[id]] <- export_figure(((a | b) + plot_layout(widths = c(1.2, 1))) /
    ((c | d) + plot_layout(widths = c(1.2, 1))) + plot_layout(heights = c(1.3, 1)), id, caption,
    c(a = source_panel(id, 'a', atlas), b = source_panel(id, 'b', comp),
      c = source_panel(id, 'c', diagnostic), d = source_panel(id, 'd', reg)),
    'All source promoter states in screen; complete selected comparator arms; five named diagnostic sets with all configured masks; post hoc regulator RNA context', height = 240)
  lapply(manifests, function(m) {
    m$analysis_objects <- paste(m$analysis_objects, 'outputs/models/04_epigenetic_deep_dive.rds', sep = ';')
    m$analysis_config <- paste(m$analysis_config, 'config/epigenetic_deep_dive.yml', sep = ';')
    m$statistical_scope <- 'post hoc descriptive diagnostics; unchanged original global enrichment BH; repeated established populations'
    m
  })
}
