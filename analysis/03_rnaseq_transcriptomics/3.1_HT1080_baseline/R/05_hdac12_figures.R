hdac12_figures <- function(z) {
  source_labels <- c(HDAC1 = 'HDAC1', HDAC2 = 'HDAC2', HDAC1_2 = 'HDAC1/2', HDAC3 = 'HDAC3')
  source_colors <- c(HDAC1 = '#0072B2', HDAC2 = '#D55E00', HDAC1_2 = '#009E73', HDAC3 = '#666666')
  source_shapes <- c(HDAC1 = 16, HDAC2 = 17, HDAC1_2 = 15, HDAC3 = 18)
  source_scales <- list(scale_colour_manual(values = source_colors, labels = source_labels, name = 'Source KD'),
    scale_shape_manual(values = source_shapes, labels = source_labels, name = 'Source KD'))
  part_labels <- c(HDAC2_only_vs_HDAC1 = 'HDAC2 only', HDAC1_only_vs_HDAC2 = 'HDAC1 only',
    concordant_shared = 'Shared, concordant', double_knockdown_only = 'Double KD only',
    conflicting_direction = 'Conflicting direction')
  ov <- z$coverage[z$coverage$represented == 'TRUE' & z$coverage$partition %in% names(part_labels) & z$coverage$genes > 0, ]
  ov$label <- factor(unname(part_labels[as.character(ov$partition)]), levels = rev(part_labels))
  assert(!anyNA(ov$label), 'HDAC membership figure label mapping failed')
  ov$source_arm <- factor(ov$source_arm, levels = c('UP', 'DN', 'MIXED'))
  a <- ggplot(ov, aes(genes, label, fill = source_arm)) + geom_col(width = 0.6) +
    geom_text(aes(label = genes), position = position_stack(vjust = 0.5), size = 6 / ggplot2::.pt,
      family = 'Arial', colour = 'white', data = ov[ov$genes > 10, ]) +
    geom_text(aes(label = genes), hjust = -0.5, size = 6 / ggplot2::.pt, family = 'Arial',
      data = ov[ov$genes <= 10, ]) +
    scale_fill_manual(values = c(UP = '#D55E00', DN = '#0072B2', MIXED = '#666666'),
      labels = c(UP = 'Source UP', DN = 'Source DN', MIXED = 'Conflict'), name = NULL) +
    labs(x = 'Unique represented genes', y = NULL, title = 'Single- and double-knockdown membership') +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06)))
  matches <- z$score_effects[z$score_effects$score_arm == 'signed' & z$score_effects$contrast %in% z$pairs$block, ]
  matches <- merge(matches, z$pairs[, c('block', 'batch', 'passage')], by.x = 'contrast', by.y = 'block')
  b <- ggplot(matches, aes(passage, percentile_difference, group = source, colour = source, shape = source)) +
    geom_hline(yintercept = 0, colour = 'gray70', linewidth = 0.2) +
    geom_line(linewidth = 0.35) + geom_point(size = 1.6) + source_scales +
    facet_grid(~batch, labeller = as_labeller(period_labels)) + scale_x_continuous(breaks = 1:3) +
    labs(x = 'Matched passage', y = 'LOY - sgCTRL signed score (pp)', title = 'Complete source programs across six matches') +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  genes <- z$shortlist$gene_name
  heat <- do.call(rbind, lapply(seq_len(nrow(z$pairs)), function(i) {
    data.frame(gene_name = genes, source_arm = z$shortlist$source_arm, readout_category = z$shortlist$readout_category,
      block = z$pairs$block[i], match = paste0(ifelse(z$pairs$batch[i] == 'AUG_2025', 'A', 'M'), z$pairs$passage[i]),
      difference = z$shortlist[[z$pairs$block[i]]])
  }))
  heat$gene_name <- factor(heat$gene_name, levels = rev(genes))
  heat$match <- factor(heat$match, levels = c(paste0('A', 1:3), paste0('M', 1:3)))
  c <- heat_effect(heat, 'gene_name', 'match', 'difference', max(abs(heat$difference)), 'Log2\ndifference') +
    facet_grid(source_arm ~ ., scales = 'free_y', space = 'free_y',
      labeller = as_labeller(c(UP = 'Source UP', DN = 'Source DN'))) +
    labs(x = 'Passage/day match', title = 'Readouts retain source-UP nonresponses')
  scenario_labels <- c(original = 'All members', omit_top_1 = 'Omit top 1', omit_top_3 = 'Omit top 3',
    omit_top_5 = 'Omit top 5', omit_HDAC1_HDAC2_HDAC3 = 'Omit HDAC1/2/3',
    omit_IFN_ECM_cycle = 'Omit IFN + matrix + cycle', retain_CPM_at_least_5 = 'Retain mean CPM >= 5')
  inf <- z$influence[z$influence$summary_arm == 'signed' & z$influence$contrast == 'combined' &
    z$influence$scenario %in% names(scenario_labels), ]
  inf$label <- factor(scenario_labels[inf$scenario], levels = rev(scenario_labels))
  d <- ggplot(inf, aes(mean_difference, label, colour = source, shape = source)) +
    geom_vline(xintercept = 0, colour = 'gray70', linewidth = 0.2) +
    geom_point(size = 1.6, position = position_dodge(width = 0.55)) + source_scales +
    labs(x = 'Mean source-UP minus source-DN log2 effect', y = NULL, title = 'Sensitivity to gene contributions') +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  id <- 'Fig08_HDAC12_gene_drivers'
  cap <- paste0('**Figure 8 | Genes underlying the HDAC-associated expression pattern.** ',
    '**a**, Mutually exclusive memberships across the Senese HDAC1, HDAC2 and double-knockdown references, restricted to represented HT1080 genes. ',
    'Single-isoform-only labels mean absence from the other single-knockdown list; double-knockdown overlap is allowed. Any directional conflict takes precedence. HDAC3-only comparison genes are omitted here. ',
    '**b**, All six matched differences in signed source-UP-minus-source-DN percentile scores. Each sample ranks its counts on the same 25,461 expressed representative genes; pp denotes percentile points. ',
    '**c**, Eight post hoc readouts selected by the recorded abundance, prior-count and temporal-stability rules, including source-UP genes with opposite or near-zero responses. A1-A3 and M1-M3 denote August and May passages. These genes do not define the program scores. ',
    '**d**, Combined mean expression contrasts after fixed exclusions. Top contributors are ranked by absolute arm-size-weighted combined effects and held fixed across periods. Retained arm means are recomputed. Complete UP and DN summaries accompany every diagnostic. ',
    'Within-source ambiguous UP/DN genes are excluded from signed scores; this removes seven HDAC3 genes from both arms. Reference differences do not calibrate isoform activity. ',
    'The 12 libraries repeatedly sample one population per condition from one generation history. All summaries are descriptive; no new P values, inferential intervals or independent validation are claimed. Protein activity and direct target regulation remain unmeasured.')
  manifests <- list()
  manifests[[id]] <- export_figure((a | b) / (c | d), id, cap,
    c(a = source_panel(id, 'a', ov), b = source_panel(id, 'b', matches),
      c = source_panel(id, 'c', heat), d = source_panel(id, 'd', inf)),
    'Complete source partitions and scores; fixed influence rules; eight post hoc readouts including source-UP nonresponses', height = 180)
  sample <- z$samples[z$samples$score_arm == 'signed', ]
  sample$source <- factor(sample$source, levels = names(source_labels))
  a <- ggplot(sample, aes(passage_number, percentile_score, colour = treatment, shape = treatment, group = treatment)) +
    geom_line(linewidth = 0.35) + geom_point(size = 1.4) +
    scale_colour_manual(values = group_colors[c('CTRL', 'LOY')], breaks = c('CTRL', 'LOY'),
      labels = c('sgCTRL', 'LOY'), name = NULL) +
    scale_shape_manual(values = c(CTRL = 16, LOY = 17), breaks = c('CTRL', 'LOY'),
      labels = c('sgCTRL', 'LOY'), name = NULL) +
    facet_grid(source ~ batch, scales = 'free_y', labeller = labeller(source = source_labels, batch = period_labels)) +
    scale_x_continuous(breaks = 1:3) +
    labs(x = 'Passage', y = 'Within-sample source-UP minus source-DN percentile score', title = 'All 12 libraries; source-specific vertical scales')
  arms <- z$arm_summary[z$arm_summary$contrast %in% c('combined', 'AUG_2025', 'MAY_2026'), ]
  arms$label <- factor(paste(source_labels[arms$source], arms$source_arm),
    levels = rev(unlist(lapply(source_labels, function(s) paste(s, c('UP', 'DN'))))))
  b <- ggplot(arms, aes(mean_difference, label, shape = contrast, colour = contrast)) +
    geom_vline(xintercept = 0, colour = 'gray70', linewidth = 0.2) + geom_point(size = 1.6) + shape_periods() +
    scale_colour_manual(values = c(combined = 'black', AUG_2025 = '#0072B2', MAY_2026 = '#D55E00'),
      labels = period_labels, name = NULL) + labs(x = 'Mean LOY - sgCTRL log2 difference', y = NULL, title = 'Separate source response arms') +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  labels <- c(cell_cycle_mitosis = 'Cell cycle/mitosis', apoptosis_p53 = 'Apoptosis/p53',
    differentiation_regulation = 'Differentiation regulation', adhesion_matrix = 'Adhesion/matrix',
    immune_regulation_response = 'Immune regulation/response', DNA_damage_repair = 'DNA damage/repair',
    any_process_annotation = 'Any of these annotations')
  pr <- z$process_effects[z$process_effects$source == 'HDAC2' & z$process_effects$source_arm == 'DN' &
    z$process_effects$contrast == 'combined' & z$process_effects$annotation %in% names(labels), ]
  pr$label <- factor(labels[pr$annotation], levels = rev(labels))
  pr$membership <- ifelse(pr$inside_annotation, 'Annotated', 'Outside annotation')
  c <- ggplot(pr, aes(mean_difference, label, colour = membership, shape = membership)) +
    geom_vline(xintercept = 0, colour = 'gray70', linewidth = 0.2) +
    geom_point(size = 1.5, position = position_dodge(width = 0.45)) +
    scale_colour_manual(values = c('Annotated' = '#0072B2', 'Outside annotation' = '#D55E00'), name = NULL) +
    scale_shape_manual(values = c('Annotated' = 16, 'Outside annotation' = 17), name = NULL) +
    labs(x = 'HDAC2 source-DN mean log2 difference', y = NULL, title = 'Fixed process annotations') +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  id <- 'SuppFig05_HDAC12_scores_processes'
  cap <- paste0('**Supplementary Figure 5 | Complete sample scores and HDAC response-arm context.** ',
    '**a**, All 12 libraries scored independently using complete represented source UP and DN memberships on the fixed 25,461-gene count-rank background. ',
    'Lines connect successive passages within each established population and period; no line connects August to May. Vertical scales differ by source and are shared between periods within each row. The score has no calibrated activity zero. ',
    '**b**, Mean LOY-minus-sgCTRL expression differences for both arms, averaged equally across periods for the combined result. Seven contradictory HDAC3 symbols are excluded from both arms. ',
    '**c**, HDAC2 source-DN means inside and outside six declared process annotations and their union. Full memberships, counts and absolute-effect fractions are exported. Definitions use frozen Hallmark, Reactome and GO sets; annotations overlap and are incomplete descriptions of biology. The union leaves 69 of 117 genes outside these annotations. ',
    'Additional process labels are descriptive annotations, not outcome-tuned masks or causal adjustments. Existing interferon, matrix/EMT and cell-cycle mask diagnostics remain separately available. ',
    'Repeated samples of one sgCTRL and one LOY population share a generation history. No gene-independent bootstrap, passage-based inferential uncertainty, new enrichment testing or claim of protein activity is used. All scores and source-arm effects are relative RNA summaries.')
  manifests[[id]] <- export_figure(a / (b | c) + plot_layout(heights = c(1.5, 1)), id, cap,
    c(a = source_panel(id, 'a', sample), b = source_panel(id, 'b', arms), c = source_panel(id, 'c', pr)),
    'All 12 libraries for four complete source programs; both response arms; six declared process annotations and their union', height = 215)
  lapply(manifests, function(m) {
    m$analysis_objects <- paste(m$analysis_objects, 'outputs/models/04_HDAC12_gene_drivers.rds', sep = ';')
    m$analysis_config <- paste(m$analysis_config, 'config/hdac12_followup.yml', sep = ';')
    m$statistical_scope <- 'descriptive post hoc; no new enrichment tests or inferential intervals'
    m
  })
}
