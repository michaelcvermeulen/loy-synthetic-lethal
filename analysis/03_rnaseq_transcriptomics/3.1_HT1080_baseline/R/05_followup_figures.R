enrichment_figures <- function(p) {
  id <- 'Fig05_expanded_enrichment'
  levels <- c('H', paste0('C', 1:9))
  summary <- p$enrichment_summary[p$enrichment_summary$contrast == 'combined', ]
  summary$collection <- factor(summary$collection, levels = rev(levels))
  sizes <- rbind(data.frame(summary[, c('collection', 'resource_sets')], type = 'Resource'),
    data.frame(collection = summary$collection, resource_sets = summary$eligible_sets, type = 'Eligible'))
  a <- ggplot(sizes, aes(resource_sets, collection, fill = type)) +
    geom_col(position = 'dodge', width = 0.75) +
    scale_fill_manual(values = c(Resource = '#BDBDBD', Eligible = '#333333'), name = NULL) +
    scale_x_continuous(labels = scales::label_comma()) + labs(x = 'Gene sets', y = NULL)
  hits <- rbind(data.frame(collection = summary$collection, sets = summary$fgsea_BH05, method = 'fgsea'),
    data.frame(collection = summary$collection, sets = summary$camera_BH05, method = 'CAMERA'))
  b <- ggplot(hits, aes(sets, collection, fill = method)) + geom_col(position = 'dodge', width = 0.75) +
    scale_fill_manual(values = c(fgsea = '#0072B2', CAMERA = '#D55E00'), name = NULL) +
    labs(x = 'Sets with global BH < 0.05', y = NULL)
  compare <- p$enrichment_comparison[p$enrichment_comparison$contrast == 'combined', ]
  compare$fgsea_signed_log10_BH <- sign(compare$NES) * -log10(pmax(compare$BH_across_collections_fgsea, 1e-12))
  compare$camera_signed_log10_BH <- ifelse(compare$Direction == 'Up', 1, -1) *
    -log10(pmax(compare$BH_across_collections_camera, 1e-12))
  compare$display_group <- ifelse(compare$both_BH05, 'Both BH < 0.05', 'Other')
  c <- ggplot(compare, aes(fgsea_signed_log10_BH, camera_signed_log10_BH, colour = display_group)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = 'gray60') +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = 'gray60') +
    geom_point(size = 0.35, alpha = 0.5) +
    scale_colour_manual(values = c('Both BH < 0.05' = '#0072B2', Other = '#BDBDBD'), name = NULL) +
    labs(x = 'fgsea signed -log10(global BH)', y = 'CAMERA signed -log10(global BH)')
  mv <- as.data.frame(data.table::fread('outputs/tables/camera_voom_mean_variance.csv'))
  trend <- as.data.frame(data.table::fread('outputs/tables/camera_voom_trend.csv'))
  q <- ggplot(mv, aes(mean_log2_count, sqrt_residual_SD)) + geom_point(size = 0.3, alpha = 0.25, colour = '#666666') +
    geom_line(data = trend, aes(y = fitted_sqrt_residual_SD), linewidth = 0.6, colour = '#0072B2') +
    labs(x = 'Mean log2 read count', y = 'Square root of residual SD', title = 'CAMERA voom precision trend')
  cap <- paste0('**Figure 5 | Expanded gene-set enrichment and CAMERA.** ',
    '**a**, Available and eligible MSigDB H and C1-C9 sets; eligibility requires 15-500 represented genes in the common expressed autosomal/non-PAR X symbol universe. ',
    '**b**, Combined-contrast sets with BH-adjusted P <0.05 for each method. BH spans all eligible sets across all ten collections, separately per method and contrast. ',
    '**c**, Comparison of direction-signed -log10 global BH values; positive denotes higher expression association in LOY. Values are capped at 12 for display; full values accompany source data. Overlapping sets are not independent discoveries. ',
    '**d**, Voom mean-variance observations and fitted trend. ',
    'fgsea tests the full equal-period-weight effect ranking. CAMERA uses raw counts, existing TMM exposures, voom weights, six passage/day block intercepts and two period-specific LOY effects; the combined contrast weights these effects equally. ',
    'CAMERA uses two-sided competitive tests, fixed inter-gene correlation 0.01, and four model residual degrees of freedom; estimated-correlation sensitivities are exported separately. ',
    'All 12 libraries sample one established population per condition from one generation history. CAMERA assumes independent residuals after fixed effects and does not model residual serial dependence. Its exploratory P values do not establish population replication, causal LOY effects or pathway activity.')
  manifest <- export_figure((a|b)/(c|q), id, cap,
    c(a = source_panel(id, 'a', sizes), b = source_panel(id, 'b', summary),
      c = source_panel(id, 'c', compare), d = source_panel(id, 'd',
        data.table::rbindlist(list(transform(mv, source = 'observed'),
          data.frame(mean_log2_count = trend$mean_log2_count, sqrt_residual_SD = trend$fitted_sqrt_residual_SD,
            source = 'fitted_trend')), use.names = TRUE))),
    'All collections, all eligible combined-contrast sets; no significance filtering before testing', height = 160)
  setNames(list(manifest), id)
}

zmynd8_figure <- function(z) {
  id <- 'Fig06_ZMYND8_followup'
  sample <- z$sample[z$sample$gene_name %in% c('ZMYND8', 'KDM5C', 'KDM5D'), ]
  sample$gene_name <- factor(sample$gene_name, levels = c('ZMYND8', 'KDM5C', 'KDM5D'))
  a <- ggplot(sample, aes(passage_number, log2_CPM, colour = treatment,
    shape = batch, group = interaction(treatment, batch))) +
    geom_line(linewidth = 0.25, alpha = 0.6) + geom_point(size = 1.4) +
    facet_wrap(~ gene_name, nrow = 1, scales = 'free_y') + colour_groups() + shape_periods() +
    scale_x_continuous(breaks = 1:3) + labs(x = 'Within-period passage', y = 'Log2 CPM (prior count 2)')
  original <- z$original
  b <- ggplot(original, aes(score, factor(contrast, levels = c('MAY_2026', 'AUG_2025', 'combined')), shape = contrast)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = 'gray60') + geom_point(size = 2) +
    shape_periods() + scale_y_discrete(labels = period_labels) +
    labs(x = 'ZMYND8 ULM score', y = NULL, title = 'Eight signed regulon targets') + theme(legend.position = 'none')
  h <- z$heat
  order <- z$targets$target[order(-abs(z$targets$signed_combined_effect))]
  labels <- setNames(paste0(z$targets$target, ifelse(z$targets$mor < 0, ' (-)', ' (+)')), z$targets$target)
  h$target_label <- factor(labels[h$target], levels = rev(labels[order]))
  h$match <- factor(sub('__passage_', ' P', sub('AUG_2025', 'Aug', sub('MAY_2026', 'May', h$sample_key))),
    levels = c(paste0('Aug P', 1:3), paste0('May P', 1:3)))
  c <- heat_effect(h, 'target_label', 'match', 'signed_difference', max(abs(h$signed_difference)), 'Signed log2\ndifference')
  omission <- z$omission[z$omission$scenario %in% c('all_targets', 'omit_MMP1', 'omit_MMP3', 'omit_MMP1_MMP3'), ]
  omission$label <- factor(omission$scenario, levels = rev(c('all_targets', 'omit_MMP1', 'omit_MMP3', 'omit_MMP1_MMP3')))
  q <- ggplot(omission, aes(score, label, shape = contrast, colour = contrast)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = 'gray60') + geom_point(size = 1.6) +
    shape_periods() + scale_colour_manual(values = c(combined = 'black', AUG_2025 = '#0072B2', MAY_2026 = '#D55E00'),
      labels = period_labels, name = NULL) +
    scale_y_discrete(labels = c(all_targets = 'All targets', omit_MMP1 = 'Omit MMP1', omit_MMP3 = 'Omit MMP3',
      omit_MMP1_MMP3 = 'Omit both')) + labs(x = 'ULM score after edge omission', y = NULL)
  cap <- paste0('**Figure 6 | ZMYND8 inference is driven by a small target subset.** ',
    '**a**, ZMYND8, KDM5C and KDM5D expression across all 12 core libraries; y-axis scales differ by gene. Lines connect passage order within population and period, without implying known elapsed time. ',
    '**b**, ZMYND8 ULM scores from the combined, August and May effect signatures using all eight expressed CollecTRI targets. Positive scores indicate agreement with the signed regulon, not measured ZMYND8 abundance or biochemical activity. ',
    '**c**, All eight targets across six timing matches, ordered by absolute combined signed effect. Expression differences are multiplied by activating (+) or repressing (-) edge signs. A1-A3 and M1-M3 denote August and May. ',
    '**d**, Post hoc removal of MMP1, MMP3 or both regulon edges while preserving the complete expression background; the full single-target omission series is exported. These are influence diagnostics, not independently validated reduced regulons. ',
    'Combined-signature ZMYND8 RNA changes little, while MMP1 and MMP3 decreases dominate the positive target pattern. Original target-pattern P values and BH corrections across all regulators are in source data. ',
    'There is one established population per condition and one shared generation origin; repeated passages do not establish biological replication. No claim of increased ZMYND8 protein activity, KDM5C compensation, chromatin binding or metastatic behavior follows from these RNA measurements.')
  manifest <- export_figure((a|b)/(wrap_elements(full = ggplotGrob(c))|wrap_elements(full = ggplotGrob(q))), id, cap,
    c(a = source_panel(id, 'a', sample), b = source_panel(id, 'b', original),
      c = source_panel(id, 'c', h), d = source_panel(id, 'd', omission)),
    'User-requested post hoc ZMYND8 follow-up; all eight represented targets; selected MMP omission diagnostics', height = 155)
  manifest$analysis_objects <- paste(manifest$analysis_objects, 'outputs/models/04_ZMYND8_followup.rds', sep = ';')
  setNames(list(manifest), id)
}
