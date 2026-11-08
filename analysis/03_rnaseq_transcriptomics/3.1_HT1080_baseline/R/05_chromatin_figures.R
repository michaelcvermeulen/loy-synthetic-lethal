chromatin_figures <- function(z) {
  labels <- c(SENESE_HDAC2_TARGETS_DN = 'HDAC2 KD: down genes',
    SENESE_HDAC1_AND_HDAC2_TARGETS_DN = 'HDAC1/2 KD: down genes',
    SENESE_HDAC1_TARGETS_DN = 'HDAC1 KD: down genes',
    PRC2_EED_DN.V1_UP = 'EED KD: up genes', PRC2_EZH2_UP.V1_UP = 'EZH2 KD: up genes',
    PRC2_SUZ12_UP.V1_UP = 'SUZ12 KD: up genes', PASINI_SUZ12_TARGETS_DN = 'Suz12 loss: down genes',
    BMI1_DN_MEL18_DN.V1_UP = 'BMI1/PCGF2 KD: up genes',
    WANG_SMARCE1_TARGETS_UP = 'SMARCE1 expression: up genes',
    SCIBETTA_KDM5B_TARGETS_UP = 'KDM5B expression: up genes',
    MIKKELSEN_MCV6_HCP_WITH_H3K27ME3 = 'H3K27me3-marked: MCV6',
    MIKKELSEN_MEF_HCP_WITH_H3K27ME3 = 'H3K27me3-marked: MEF')
  data <- z$results[z$results$pathway %in% names(labels), ]
  data$label <- factor(labels[data$pathway], levels = rev(labels))
  a <- ggplot(data, aes(NES, label, shape = contrast, colour = contrast)) +
    geom_vline(xintercept = 0, colour = 'gray70', linewidth = 0.2) + geom_point(size = 1.4) +
    shape_periods() + scale_colour_manual(values = c(combined = 'black', AUG_2025 = '#0072B2', MAY_2026 = '#D55E00'),
      labels = period_labels, name = NULL) + labs(x = 'Original fgsea NES', y = NULL) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  pair_labels <- c(SENESE_HDAC2_TARGETS = 'HDAC2 knockdown', SENESE_HDAC1_AND_HDAC2_TARGETS = 'HDAC1/2 knockdown',
    SENESE_HDAC1_TARGETS = 'HDAC1 knockdown', PRC2_EED_DN.V1 = 'EED knockdown',
    PRC2_EZH2_UP.V1 = 'EZH2 knockdown', PRC2_SUZ12_UP.V1 = 'SUZ12 knockdown',
    PASINI_SUZ12_TARGETS = 'Suz12 deficiency', BMI1_DN_MEL18_DN.V1 = 'BMI1/PCGF2 knockdown',
    WANG_SMARCE1_TARGETS = 'SMARCE1 expression', SCIBETTA_KDM5B_TARGETS = 'KDM5B expression')
  pr <- z$pairs[z$pairs$contrast == 'combined' & z$pairs$source_pair %in% names(pair_labels), ]
  pr$label <- factor(pair_labels[pr$source_pair], levels = rev(pair_labels))
  arms <- rbind(data.frame(pr, arm = 'UP in source', difference = pr$UP_mean_difference),
    data.frame(pr, arm = 'DN in source', difference = pr$DN_mean_difference))
  b <- ggplot(arms, aes(difference, label, shape = arm, colour = arm)) +
    geom_vline(xintercept = 0, colour = 'gray70', linewidth = 0.2) +
    geom_line(aes(group = source_pair), colour = 'gray75', linewidth = 0.3) + geom_point(size = 1.6) +
    scale_colour_manual(values = c('UP in source' = '#D55E00', 'DN in source' = '#0072B2'), name = NULL) +
    scale_shape_manual(values = c('UP in source' = 17, 'DN in source' = 16), name = NULL) +
    labs(x = 'Mean LOY - sgCTRL log2 difference', y = NULL)
  ov <- z$overlap[z$overlap$pathway %in% names(labels), ]; ov$label <- factor(labels[ov$pathway], levels = rev(labels))
  ov$program <- factor(ov$program, levels = c('IFN', 'ECM_EMT', 'cell_cycle'))
  ov$percent <- 100 * ov$membership_fraction
  c <- ggplot(ov, aes(program, label, fill = percent)) + geom_tile() +
    geom_text(aes(label = sprintf('%.0f', percent)), family = 'Arial', size = 6 / ggplot2::.pt,
      colour = ifelse(ov$percent > 30, 'white', 'black')) +
    scale_fill_gradient(low = 'white', high = '#0072B2', limits = c(0, max(50, max(ov$percent))), name = 'Member\noverlap (%)') +
    scale_x_discrete(labels = c(IFN = 'IFN', ECM_EMT = 'Matrix/EMT', cell_cycle = 'Cell cycle')) +
    labs(x = NULL, y = NULL) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = 'bottom', legend.key.width = grid::unit(12, 'mm'))
  se <- z$sensitivity[z$sensitivity$pathway %in% names(labels) & z$sensitivity$scenario %in%
    c('original', 'without_IFN_ECM', 'without_IFN_ECM_cycle'), ]
  se$label <- factor(labels[se$pathway], levels = rev(labels))
  se$scenario_label <- factor(se$scenario, levels = c('original', 'without_IFN_ECM', 'without_IFN_ECM_cycle'),
    labels = c('Original', 'Remove IFN + matrix', 'Also remove cell cycle'))
  q <- ggplot(se, aes(NES, label, shape = scenario_label, colour = scenario_label)) +
    geom_vline(xintercept = 0, colour = 'gray70', linewidth = 0.2) + geom_point(size = 1.4, na.rm = TRUE) +
    scale_colour_manual(values = c('black', '#0072B2', '#D55E00'), name = NULL) +
    scale_shape_manual(values = c(16, 17, 18), name = NULL) + labs(x = 'fgsea NES after universe exclusion', y = NULL) +
    guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
  id <- 'Fig07_chromatin_signatures'
  cap <- paste0('**Figure 7 | Chromatin-associated expression signatures and their specificity.** ',
    '**a**, Original enrichment for 12 selected signatures; names describe the source experiment, not regulator activity in HT1080. KD denotes knockdown. ',
    '**b**, Mean combined expression differences for both response arms of ten source perturbations. UP and DN refer to the source response; opposite directions support a more coherent perturbation-like pattern than a decrease of both arms. ',
    '**c**, Percentage of represented signature genes shared with fixed interferon, matrix/EMT and cell-cycle reference sets, defined without selecting on LOY effects. ',
    '**d**, Combined fgsea after excluding reference genes from both the signature and ranked universe. Each NES is normalized within its own universe and is not an activity difference across scenarios. ',
    'All 37 focused set definitions, eligible and ineligible sets, original global BH values, CAMERA results and five exclusion sensitivities are exported. Original tests retain BH across all 23,696 eligible MSigDB sets; focused sensitivities report nominal P values without confirmatory FDR claims. ',
    'H3K27me3 sets describe marked genes in mouse source cells mapped to human genes; they do not measure histone marks in HT1080. The 12 libraries repeatedly sample one established population per condition from one generation history. CAMERA assumes independent residuals after passage fixed effects; estimated-correlation sensitivity remains inconclusive. No protein activity or chromatin-occupancy claim is established.')
  manifests <- list()
  manifests[[id]] <- export_figure((a|b)/(c|q), id, cap,
    c(a = source_panel(id, 'a', data), b = source_panel(id, 'b', arms),
      c = source_panel(id, 'c', ov), d = source_panel(id, 'd', se)),
    'Twelve named mechanistic examples, including non-enriched controls; ten complete source-response pairs', height = 165)
  reg <- z$regulator_matches[z$regulator_matches$chrom != 'chrY', ]
  gene_order <- z$cfg$regulators[z$cfg$regulators %in% reg$gene_name]
  reg$gene_name <- factor(reg$gene_name, levels = rev(gene_order))
  reg$match <- factor(sub('__passage_', ' P', sub('AUG_2025', 'Aug', sub('MAY_2026', 'May', reg$sample_key))),
    levels = c(paste0('Aug P', 1:3), paste0('May P', 1:3)))
  a <- heat_effect(reg, 'gene_name', 'match', 'difference', max(abs(reg$difference)), 'Log2\ndifference')
  ys <- z$regulator_matches[z$regulator_matches$chrom == 'chrY', ]; ys$match <- reg$match[match(ys$sample_key, reg$sample_key)]
  b <- ggplot(ys, aes(match, difference, group = gene_name, shape = gene_name)) +
    geom_line(linewidth = 0.3, colour = 'gray50') + geom_point(size = 2) +
    scale_x_discrete(labels = function(s) sub('May P', 'M', sub('Aug P', 'A', s))) +
    labs(x = 'Passage/day match', y = 'LOY - sgCTRL log2 difference', shape = NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  id <- 'SuppFig04_chromatin_regulators'
  cap <- paste0('**Supplementary Figure 4 | RNA expression of chromatin regulators.** ',
    '**a**, All 32 mechanism-selected autosomal or X-linked regulators across six matched passages; colors show original TMM-normalized log2 CPM differences with prior count 2, without row scaling. ',
    '**b**, UTY and KDM5D use a separate quantitative axis because their depletion is much larger than the other RNA differences. ',
    'A1-A3 and M1-M3 denote August and May matches. Regulator selection follows the mechanistic questions, without significance or magnitude filtering. ',
    'Gene-level effects, raw counting identities, observed passage ranges and leave-one-match-out summaries accompany the full tables. ',
    'RNA differences do not measure protein abundance, complex assembly, nuclear localization, histone methylation/acetylation or catalytic activity. Both periods sample the same two established populations from one generation history. ',
    'No gene-level P values or independent-population uncertainty intervals are reported.')
  manifests[[id]] <- export_figure(a|b, id, cap,
    c(a = source_panel(id, 'a', reg), b = source_panel(id, 'b', ys)),
    'All 34 mechanism-selected regulators; Y-linked genes separated to preserve the smaller-effect scale', height = 165)
  lapply(manifests, function(m) {
    m$analysis_objects <- paste(m$analysis_objects, 'outputs/models/04_chromatin_followup.rds', sep = ';')
    m
  })
}
