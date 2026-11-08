hdac_external_labels <- function(contrasts) {
  suffix <- c(protein_loss = 'KO / WT', inactive_vs_active_transgene = 'CI / active (KO)',
    inactive_transgene_in_WT = 'CI (WT) / WT', inactive_transgene_in_KO_vs_WT = 'CI (KO) / WT',
    inactive_protein_presence_in_KO = 'CI (KO) / KO', active_transgene_control_vs_WT = 'active (KO) / WT')
  labels <- paste(contrasts$target, suffix[as.character(contrasts$comparison_type)])
  drug <- contrasts$comparison_type == 'inhibitor_exposure'
  labels[drug] <- paste(ifelse(contrasts$target[drug] == 'MS275', 'MS-275 label MS', 'JQ12'),
    contrasts$sample_ID_stratum[drug], '/ WT')
  setNames(unname(labels), contrasts$contrast)
}

hdac_external_figures <- function(z, external, e) {
  labels <- hdac_external_labels(z$contrasts)
  selected <- z$contrasts$contrast[z$contrasts$display_tier == 'primary' |
    z$contrasts$comparison_type == 'active_transgene_control_vs_WT']
  make_labels <- function(ids) factor(unname(labels[as.character(ids)]), levels = rev(unname(labels[selected])))
  cc <- z$correlations[z$correlations$scenario == 'original' & z$correlations$external_contrast %in% selected, ]
  cc$label <- make_labels(cc$external_contrast)
  a <- ggplot(cc, aes(spearman, label, colour = contrast, shape = contrast)) +
    geom_vline(xintercept = 0, colour = 'gray60', linewidth = 0.25) +
    geom_point(size = 1.35, position = position_dodge(width = 0.55)) + shape_periods() +
    scale_colour_manual(values = c(combined = 'black', AUG_2025 = '#0072B2', MAY_2026 = '#D55E00'),
      labels = period_labels, name = NULL) +
    labs(x = 'Spearman correlation (14,881 genes)', y = NULL, title = 'Weak agreement of complete effect rankings') +
    theme(axis.text.y = element_text(size = 6)) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  arms <- z$arm_means[z$arm_means$scenario == 'original' & z$arm_means$contrast == 'combined' &
    z$arm_means$external_contrast %in% selected, ]
  arms <- rbind(transform(arms, source_arm = 'UP', mean_effect = UP_mean, genes = UP_genes),
    transform(arms, source_arm = 'DN', mean_effect = DN_mean, genes = DN_genes))
  arms$label <- make_labels(arms$external_contrast)
  b <- ggplot(arms, aes(mean_effect, label, colour = source_arm, shape = source_arm)) +
    geom_vline(xintercept = 0, colour = 'gray60', linewidth = 0.25) +
    geom_point(size = 1.5, position = position_dodge(width = 0.35)) +
    scale_colour_manual(values = c(UP = '#D55E00', DN = '#0072B2'), name = 'Source arm') +
    scale_shape_manual(values = c(UP = 17, DN = 16), name = 'Source arm') +
    labs(x = 'Mean LOY - sgCTRL log2 difference', y = NULL, title = 'Both source arms often decrease in LOY') +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  en <- z$enrichment[z$enrichment$scenario == 'original' & z$enrichment$contrast == 'combined', ]
  en$fgsea_support <- -log10(en$fgsea_secondary_BH)
  en$camera_support <- -log10(en$camera_fixed_secondary_BH)
  en$annotation <- ''
  named <- en$external_contrast %in% c('HDAC2_KO_vs_WT', 'HDAC2_CI_vs_active_in_KO')
  en$annotation[named] <- paste(ifelse(en$external_contrast[named] == 'HDAC2_KO_vs_WT', '2 KO', '2 CI/active'), en$source_arm[named])
  palette <- c(HDAC1 = '#0072B2', HDAC2 = '#D55E00', HDAC3 = '#009E73', HDAC8 = '#CC79A7', MS275 = '#E69F00', JQ12 = '#555555')
  c <- ggplot(en, aes(fgsea_support, camera_support, colour = target, shape = source_arm)) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', colour = 'gray55', linewidth = 0.2) +
    geom_vline(xintercept = -log10(0.05), linetype = 'dashed', colour = 'gray55', linewidth = 0.2) +
    geom_point(size = 1.6) + scale_colour_manual(values = palette, name = NULL) +
    scale_shape_manual(values = c(UP = 17, DN = 16), name = 'Source arm') +
    ggrepel::geom_text_repel(aes(label = annotation), seed = 42, size = 6 / ggplot2::.pt,
      family = 'Arial', min.segment.length = 0, max.overlaps = Inf, colour = 'black', show.legend = FALSE) +
    labs(x = 'fgsea -log10 secondary BH', y = 'Fixed-correlation CAMERA -log10 BH',
      title = '54 sets; estimated-correlation CAMERA: 0 at BH < 0.05') +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 1))
  heat_contrasts <- c('HDAC2_KO_vs_WT', 'HDAC2_CI_vs_active_in_KO', 'HDAC2_CI_in_WT_vs_WT', 'HDAC2_active_in_KO_vs_WT')
  heat_labels <- c(AUG_2025 = 'LOY Aug', MAY_2026 = 'LOY May', HDAC2_KO_vs_WT = 'HDAC2 KO',
    HDAC2_CI_vs_active_in_KO = 'CI / active (KO)', HDAC2_CI_in_WT_vs_WT = 'CI (WT)',
    HDAC2_active_in_KO_vs_WT = 'Active (KO)')
  genes <- unlist(z$cfg$readout_genes)
  heat <- do.call(rbind, lapply(names(heat_labels), function(contrast) {
    if (contrast %in% heat_contrasts) {
      r <- z$readouts[z$readouts$external_contrast == contrast, ]
      values <- r$external_logFC[match(genes, r$HT1080_symbol)]
      scope <- 'external_HAP1_edgeR_logFC'
    } else {
      loci <- e$mapping$gene_id[match(genes, e$mapping$gene_name)]
      values <- e$effects[[contrast]][match(loci, e$effects$gene_id)]
      scope <- 'HT1080_descriptive_matched_log2_difference'
    }
    data.frame(gene = genes, contrast = contrast, value = values, scope = scope)
  }))
  heat$gene <- factor(heat$gene, levels = rev(genes))
  heat$label <- factor(heat_labels[as.character(heat$contrast)], levels = unname(heat_labels))
  heat$text <- ifelse(is.na(heat$value), 'NA', sprintf('%.1f', ifelse(abs(heat$value) < 0.05, 0, heat$value)))
  limit <- max(abs(heat$value), na.rm = TRUE)
  d <- ggplot(heat, aes(label, gene, fill = value)) + geom_tile() +
    geom_text(aes(label = text), size = 6 / ggplot2::.pt, family = 'Arial', colour = 'black') +
    scale_fill_gradient2(low = '#0072B2', mid = 'white', high = '#D55E00', midpoint = 0,
      limits = c(-limit, limit), na.value = '#D0D0D0', name = 'Log2\neffect') +
    labs(x = NULL, y = NULL, title = 'Earlier readouts also respond to active transgene') +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'right',
      legend.key.height = grid::unit(12, 'mm'), legend.key.width = grid::unit(2.5, 'mm'))
  id <- 'Fig09_HDAC_external_comparison'
  cap <- paste0('**Figure 9 | Public HDAC perturbations do not establish a specific HDAC2 catalytic-loss state in LOY.** ',
    '**a**, Spearman correlations between HAP1 perturbation and HT1080 effect rankings on 14,881 common loci. All 15 primary contrasts and all four active-transgene controls are shown. ',
    'A slash denotes numerator versus denominator; CI denotes catalytically inactive transgene, with its cellular background in parentheses. ',
    '**b**, Combined HT1080 means of separately defined source UP and DN arms, in the same row order. Sets use external nonzero-effect BH < 0.05 and absolute log2 fold change >= 0.5, capped at 250 genes per arm by QL F before overlap filtering. ',
    '**c**, Exploratory support across all 54 sets; dashed lines mark secondary BH = 0.05. Labels identify HDAC2 knockout and CI-versus-active arms. The estimated-correlation sensitivity has no combined BH < 0.05 result. ',
    '**d**, Eight previously selected readouts, with unscaled effects. External columns are versus WT except CI/active. Gray NA indicates absent common-background coverage, not zero expression. ',
    'HAP1 models use 63 libraries in two sample-ID strata, with author-reported culture triplicates of established lines; independent clone replication is unverified. MS-275 label MS has unresolved dose/time metadata. HDAC1/2 transgenes are mouse cDNAs aligned to a human reference. ',
    'HT1080 repeatedly samples one population per condition from a shared generation history. Correlations and arm means are descriptive. CAMERA assumes independent residuals; neither method establishes enzyme activity, independent-population validation or drug sensitivity. Full 27-contrast results accompany the figure.')
  manifests <- list()
  manifests[[id]] <- export_figure(((a | b) + plot_layout(widths = c(1.25, 0.75))) / (c | d) +
    plot_layout(heights = c(1.25, 1)), id, cap,
    c(a = source_panel(id, 'a', cc), b = source_panel(id, 'b', arms),
      c = source_panel(id, 'c', en), d = source_panel(id, 'd', heat)),
    'All 15 primary contrasts plus all four active-transgene controls; all 54 sets in panel c; same eight step-2 readouts', height = 230)

  pc <- external$pca[, c('sample_key', 'source_group', 'sample_ID_stratum', 'model_group',
    'replicate_label', 'PC1', 'PC2', 'PC1_variance', 'PC2_variance')]
  pc$facet <- sprintf('%s: PC1 %.1f%%; PC2 %.1f%%', pc$sample_ID_stratum, 100 * pc$PC1_variance, 100 * pc$PC2_variance)
  pc$annotation <- ''
  for (idx in split(seq_len(nrow(pc)), pc$model_group)) {
    dist <- (pc$PC1[idx] - mean(pc$PC1[idx]))^2 + (pc$PC2[idx] - mean(pc$PC2[idx]))^2
    pc$annotation[idx[which.min(dist)]] <- pc$source_group[idx[1]]
  }
  a <- ggplot(pc, aes(PC1, PC2)) + geom_point(size = 1.25, colour = '#0072B2') +
    ggrepel::geom_text_repel(aes(label = annotation), family = 'Arial',
      size = 6 / ggplot2::.pt, seed = 42, min.segment.length = 0, max.overlaps = Inf,
      box.padding = 0.4, point.padding = 0.5, colour = 'black') +
    facet_wrap(~facet, scales = 'free', nrow = 1) +
    labs(x = 'PC1', y = 'PC2', title = 'External expression QC within each sample-ID stratum')
  compat <- external$compatibility
  compat$label <- ifelse(compat$source_worksheet == '24h MS-275', 'MS-275 label MS / published 24h*', compat$source_worksheet)
  compat$label <- factor(compat$label, levels = rev(compat$label))
  b <- ggplot(compat, aes(spearman, label)) + geom_point(size = 1.6, colour = '#0072B2') +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = 'Spearman correlation with published log2 effects', y = NULL,
      title = 'Published-effect compatibility\n*MS correspondence is unverified') +
    theme(axis.text.y = element_text(size = 6))
  sens <- z$correlations[z$correlations$contrast == 'combined', ]
  mask_labels <- c(original = 'Original', omit_IFN_ECM_cycle = 'Omit programs',
    omit_HDAC_family = 'Omit HDAC family', omit_IFN_ECM_cycle_and_HDAC_family = 'Omit both')
  sens$label <- factor(labels[as.character(sens$external_contrast)], levels = rev(unname(labels)))
  sens$mask <- factor(mask_labels[as.character(sens$scenario)], levels = unname(mask_labels))
  sens$text <- sprintf('%+.2f', sens$spearman)
  c <- heat_effect(sens, 'label', 'mask', 'spearman', max(abs(sens$spearman)), 'Spearman\ncorrelation') +
    geom_text(aes(label = text), family = 'Arial', size = 6 / ggplot2::.pt) +
    labs(title = 'All 27 contrasts and four fixed universes') +
    theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'bottom') +
    guides(fill = guide_colourbar(barwidth = grid::unit(32, 'mm'), barheight = grid::unit(2.5, 'mm')))
  id <- 'SuppFig06_HDAC_external_QC_sensitivity'
  cap <- paste0('**Supplementary Figure 6 | Source-study checks and complete perturbation sensitivity.** ',
    '**a**, Centered, unscaled PCA of separately TMM-normalized HAP1 log2 CPM within S47 and S52 sample-ID strata, using 19,734 and 18,980 expressed autosomal/non-PAR X genes. ',
    'All 63 modeled libraries are shown, with one label per source group. Six MS1/MS3 libraries with unresolved treatment metadata are excluded from models and retained in raw QC. Strata are inferred from sample IDs, not verified batch dates. ',
    '**b**, Gene-wise compatibility of reconstructed edgeR effects with nine published DESeq2 worksheets. The MS-275 label-MS comparison has unverified correspondence to the published 24-hour worksheet and should not be read as a confirmed reproduction. ',
    'Separate numerical diagnostics show stronger correspondence for label MS3 but do not establish its treatment identity; primary membership rules remain frozen. ',
    '**c**, Combined HT1080 versus HAP1 rank correlations for all 27 contrasts, on the original common universe and after removing the previously defined IFN/matrix/cycle union, HDAC1/2/3/8-family genes and pseudogenes, or both. Numbers give correlations rounded to two decimals. Genes are removed from both response sets and the analysis background. ',
    'A slash denotes numerator versus denominator; CI indicates inactive transgene. Source-group identifiers and all contrasts are decoded in the sample and contrast manifests. ',
    'These diagnostics do not supply gene-independent uncertainty or independently generated HT1080 replication. Source clonal histories, transgene context, unresolved drug metadata and repeated HT1080 populations limit mechanistic interpretation.')
  manifests[[id]] <- export_figure(a / ((b | c) + plot_layout(widths = c(0.9, 1.1))) +
    plot_layout(heights = c(0.8, 1.2)), id, cap,
    c(a = source_panel(id, 'a', pc), b = source_panel(id, 'b', compat), c = source_panel(id, 'c', sens)),
    'All modeled external libraries; all nine published compatibility comparisons; all 27 contrasts and four masks', height = 240)
  lapply(manifests, function(m) {
    m$analysis_objects <- paste(m$analysis_objects, 'outputs/models/04_HDACEXT_source_models.rds',
      'outputs/models/04_HDACEXT_comparison.rds', sep = ';')
    m$analysis_config <- paste(m$analysis_config, 'config/hdac12_external.yml', sep = ';')
    m$statistical_scope <- 'external clone-conditional QL; descriptive HT1080 concordance; separate exploratory 54-set BH family; conditional CAMERA'
    m
  })
}
