suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
group_colors <- c(CTRL = '#0072B2', LOY = '#D55E00', PARENTAL = '#666666')
period_shapes <- c(AUG_2025 = 16, MAY_2026 = 17, combined = 18)
period_labels <- c(AUG_2025 = 'Aug 2025', MAY_2026 = 'May 2026', combined = 'Combined')
effect_colors <- c('Higher in LOY' = '#D55E00', 'Lower in LOY' = '#0072B2', 'Below 0.5' = '#C7C7C7')
nature_theme <- function() {
  theme_classic(base_size = 7, base_family = 'Arial') +
    theme(axis.line = element_line(linewidth = 0.2), axis.ticks = element_line(linewidth = 0.2),
      axis.ticks.length = grid::unit(1, 'mm'), axis.text = element_text(size = 6.5, colour = 'black'),
      axis.title = element_text(size = 7), plot.title = element_text(size = 7, face = 'plain'),
      strip.background = element_blank(), strip.text = element_text(size = 6.5, colour = 'black'),
      legend.title = element_text(size = 6.5), legend.text = element_text(size = 6.5),
      legend.key.size = grid::unit(3, 'mm'), legend.spacing.x = grid::unit(1, 'mm'),
      legend.margin = margin(0, 0, 0, 0), legend.position = 'bottom',
      plot.tag = element_text(size = 8, face = 'bold', family = 'Arial'),
      plot.margin = margin(3, 4, 3, 3))
}
colour_groups <- function() scale_colour_manual(values = group_colors, breaks = names(group_colors),
                                                  labels = c('sgCTRL', 'LOY', 'Parental'), name = NULL)
shape_periods <- function() scale_shape_manual(values = period_shapes, labels = period_labels, name = NULL)
source_panel <- function(id, panel, data) {
  path <- file.path('outputs/source_data', paste0(id, '_', panel, '.csv'))
  data.table::fwrite(as.data.frame(data), path, na = 'NA')
  path
}
export_figure <- function(plot, id, caption, sources, selection, height = 160) {
  plot <- plot + plot_annotation(tag_levels = 'a', theme = nature_theme())
  saveRDS(plot, file.path('outputs/models', paste0(id, '_plot.rds')))
  stem <- file.path('outputs/figures', id)
  ggsave(paste0(stem, '.pdf'), plot, device = grDevices::cairo_pdf, width = 183, height = height,
         units = 'mm', family = 'Arial', bg = 'white', limitsize = FALSE)
  ggsave(paste0(stem, '.svg'), plot, device = svglite::svglite, width = 183, height = height,
         units = 'mm', bg = 'white', limitsize = FALSE)
  ggsave(paste0(stem, '.png'), plot, device = ragg::agg_png, width = 183, height = height,
         units = 'mm', dpi = 600, bg = 'white', limitsize = FALSE)
  caption_path <- file.path('outputs/captions', paste0(id, '.caption.md'))
  writeLines(caption, caption_path, useBytes = TRUE)
  assert(length(strsplit(caption, '[[:space:]]+')[[1]]) < 250, paste('Caption exceeds 250 words:', id))
  data.frame(figure = id, panel = names(sources), source_data = unname(sources),
    pdf = paste0(stem, '.pdf'), svg = paste0(stem, '.svg'), png = paste0(stem, '.png'),
    caption = caption_path, width_mm = 183, height_mm = height, png_dpi = 600,
    font = 'Arial', base_font_pt = 7, panel_font_pt = 8, selection_rule = selection,
    analysis_objects = paste(file.path('outputs/models', c('01_ingest.rds','02_expression.rds',
      '03_pathways_regulators.rds','04_dosage_repeatability.rds')), collapse = ';'),
    raw_counts_object = 'outputs/models/raw_counts_SE.rds',
    gene_annotation = 'outputs/tables/gene_annotation_gencode49.csv',
    sample_manifest = 'outputs/tables/sample_manifest.csv', analysis_config = 'config/analysis.yml',
    contrast = 'LOY minus sgCTRL, matched within period; equal period weights',
    statistical_scope = 'descriptive; pathway/regulator panels include explicitly conditional exploratory statistics',
    visual_review = 'pending')
}
short_sample <- function(m) paste(ifelse(m$batch == 'AUG_2025', 'A', 'M'),
                                   ifelse(m$treatment == 'CTRL', 'C', ifelse(m$treatment == 'LOY', 'L', 'P')),
                                   m$passage_number, sep = '')
heat_effect <- function(dat, row, column, value, limit, legend) {
  ggplot(dat, aes(x = .data[[column]], y = .data[[row]], fill = .data[[value]])) +
    geom_tile() + scale_x_discrete(labels = function(z) sub('May P', 'M', sub('Aug P', 'A', z))) +
    scale_fill_gradient2(low = '#0072B2', mid = 'white', high = '#D55E00',
      midpoint = 0, limits = c(-limit, limit), name = legend) +
    labs(x = NULL, y = NULL) + nature_theme() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = 'right',
           legend.key.height = grid::unit(12, 'mm'), legend.key.width = grid::unit(2.5, 'mm'))
}
make_figures <- function(x, e, p, d, cfg) {
  # Prebuilt grobs need a device with real Arial metrics, rather than R's default PostScript device.
  capture <- ragg::agg_capture(width = 183, height = 165, units = 'mm', res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  theme_set(nature_theme()); manifests <- list()
  core <- x$core; core$short <- short_sample(core)
  core$short <- factor(core$short, levels = core$short[order(core$batch, core$treatment, core$passage_number)])
  allmeta <- x$meta; allmeta$short <- short_sample(allmeta)
  # Figure 1: core PCA, assignment QC, sample relationships and actual design.
  pc <- merge(e$pca, core[, c('sample_key', 'short')], by = 'sample_key')
  a <- ggplot(pc, aes(PC1, PC2, colour = treatment, shape = batch)) +
    geom_point(size = 2) + ggrepel::geom_text_repel(aes(label = passage_number), size = 6 / ggplot2::.pt,
      seed = 42, min.segment.length = 0, show.legend = FALSE, box.padding = 0.2) +
    colour_groups() + shape_periods() + labs(x = sprintf('PC1 (%.1f%%)', 100 * e$pca_variance[1]),
      y = sprintf('PC2 (%.1f%%)', 100 * e$pca_variance[2]))
  qc <- merge(e$qc, allmeta[, c('sample_key', 'short')], by = 'sample_key')
  qc$short <- factor(qc$short, levels = allmeta$short)
  b <- ggplot(qc, aes(short, assigned_fraction * 100, fill = treatment)) + geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors, labels = c(CTRL = 'sgCTRL', LOY = 'LOY', PARENTAL = 'Parental'), name = NULL) +
    labs(x = 'Period / population / passage', y = 'Assigned reads (%)') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  cr <- as.data.frame(matrix_long(e$correlations, 'r', 'sample_1'))
  cr$sample_1 <- factor(core$short[match(cr$sample_1, core$sample_key)], levels = levels(core$short))
  cr$sample_2 <- factor(core$short[match(cr$sample_key, core$sample_key)], levels = rev(levels(core$short)))
  c <- ggplot(cr, aes(sample_1, sample_2, fill = r)) + geom_tile() +
    scale_fill_viridis_c(option = 'cividis', name = 'Pearson r', breaks = scales::breaks_pretty(n = 3),
      guide = guide_colourbar(barwidth = grid::unit(30, 'mm'), barheight = grid::unit(2.5, 'mm'))) +
    coord_equal() + labs(x = NULL, y = NULL) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  design <- core; design$occasion <- design$passage_number + ifelse(design$batch == 'MAY_2026', 4, 0)
  design$population <- factor(design$treatment, levels = c('LOY', 'CTRL'))
  q <- ggplot(design, aes(occasion, population, colour = treatment, group = treatment)) +
    geom_line(linewidth = 0.4) + geom_point(aes(shape = batch), size = 2.3) +
    geom_text(aes(label = passage_number), nudge_y = 0.15, family = 'Arial', size = 6 / ggplot2::.pt, show.legend = FALSE) +
    scale_x_continuous(breaks = c(2, 6), labels = c('Aug 2025', 'May 2026'), limits = c(0.5, 7.5)) +
    scale_y_discrete(labels = c(LOY = 'LOY', CTRL = 'sgCTRL'), expand = expansion(add = 0.4)) +
    colour_groups() + shape_periods() + labs(x = 'Passage order within each period', y = 'Established population',
      title = 'One shared nucleofection origin') + theme(legend.position = 'none')
  id <- 'Fig01_design_QC'
  cap <- paste0('**Figure 1 | Sampling design and count-level quality control.** ',
    '**a**, Principal-component analysis of ', format(sum(e$effects$primary_eligible), big.mark = ','),
    ' expressed autosomal and non-PAR X genes, using TMM-normalized log2 CPM (prior count 2), centered by gene without variance scaling. ',
    '**b**, Assigned reads divided by all reads in featureCounts summaries; all 15 libraries are shown. ',
    '**c**, Pearson correlations of core normalized expression. Codes denote August (A) or May (M), sgCTRL (C), LOY (L), or parental (P), and passage. ',
    '**d**, Confirmed sampling structure; spacing indicates order, not elapsed time. ',
    'The core comprises 12 RNA libraries: one established sgCTRL and one LOY population, each measured at three successive passages in both periods, from one shared nucleofection origin. ',
    'Passage/day matches are within period. Numbers label passages; circles denote August and triangles May. Blue denotes sgCTRL, vermilion LOY, and gray parental. ',
    'Three August parental libraries provide context and do not enter core normalization or effects. No samples were excluded. All panels are descriptive; libraries and passage matches are not independent population-generation replicates. Source data accompany panels a-d.')
  manifests[[id]] <- export_figure((a | b) / (c | q), id, cap,
    c(a = source_panel(id, 'a', pc), b = source_panel(id, 'b', qc), c = source_panel(id, 'c', cr), d = source_panel(id, 'd', design)),
    'All core samples; all 15 libraries in assignment panel; no outcome-based exclusion')
  # Figure 2: summaries keep all six differences accessible.
  ef <- e$effects[e$effects$primary_eligible, ]
  candidates <- ef[ef$descriptive_candidate & ef$gene_type == 'protein_coding' & ef$mean_CPM >= 2, ]
  candidates <- candidates[order(-abs(candidates$mean_difference), candidates$gene_id), ]
  selected <- head(candidates, 12)
  labels <- head(selected, 6)
  a <- ggplot(ef, aes(mean_logCPM, mean_difference, colour = magnitude_category)) +
    geom_point(size = 0.3, alpha = 0.45) + geom_hline(yintercept = c(-0.5, 0, 0.5), linewidth = 0.2,
      linetype = c('dashed', 'solid', 'dashed'), colour = 'gray55') +
    ggrepel::geom_text_repel(data = labels, aes(label = gene_name), size = 6 / ggplot2::.pt, seed = 42,
      max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
    scale_colour_manual(values = effect_colors, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1), nrow = 1)) +
    labs(x = 'Mean log2 CPM', y = 'Mean matched log2 difference')
  b <- ggplot(ef, aes(AUG_2025, MAY_2026, colour = magnitude_category)) + geom_point(size = 0.3, alpha = 0.45) +
    geom_abline(slope = 1, intercept = 0, colour = 'gray40', linewidth = 0.25) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = 'gray60') +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = 'gray60') +
    scale_colour_manual(values = effect_colors, name = NULL) + guides(colour = 'none') +
    labs(x = 'August mean log2 difference', y = 'May mean log2 difference',
      title = sprintf('All eligible genes: Pearson r = %.2f', cor(ef$AUG_2025, ef$MAY_2026)))
  trajectory_genes <- head(selected, 4)$gene_id
  tr <- merge(as.data.frame(matrix_long(e$logCPM[trajectory_genes, , drop = FALSE], 'logCPM')), core, by = 'sample_key')
  tr$gene_name <- ef$gene_name[match(tr$gene_id, ef$gene_id)]
  c <- ggplot(tr, aes(passage_number, logCPM, colour = treatment, shape = batch,
                       group = interaction(treatment, batch))) + geom_line(linewidth = 0.3) + geom_point(size = 1.3) +
    facet_wrap(~ gene_name, ncol = 2, scales = 'free_y') + colour_groups() + shape_periods() +
    scale_x_continuous(breaks = 1:3) + labs(x = 'Passage', y = 'Normalized log2 CPM')
  obs <- as.data.frame(matrix_long(e$differences[selected$gene_id, , drop = FALSE], 'difference'))
  obs$batch <- x$pairs$batch[match(obs$sample_key, x$pairs$block)]
  obs$gene_name <- factor(selected$gene_name[match(obs$gene_id, selected$gene_id)], levels = rev(selected$gene_name))
  selected$label <- factor(selected$gene_name, levels = rev(selected$gene_name))
  q <- ggplot(obs, aes(difference, gene_name)) +
    geom_vline(xintercept = 0, colour = 'gray50', linewidth = 0.2) +
    geom_segment(data = selected, aes(x = minimum_match, xend = maximum_match, y = label, yend = label),
      inherit.aes = FALSE, linewidth = 0.25, colour = 'gray55') +
    geom_point(aes(shape = batch), size = 1.3, colour = '#555555', alpha = 0.75) +
    geom_point(data = selected, aes(mean_difference, label), inherit.aes = FALSE, shape = 18, size = 1.8, colour = 'black') +
    shape_periods() + labs(x = 'Matched log2 difference (LOY - sgCTRL)', y = NULL)
  id <- 'Fig02_expression_repeatability'
  cap <- paste0('**Figure 2 | Expression differences and their persistence across sampling periods.** ',
    '**a**, Equal-period-weighted mean of six passage-matched LOY-minus-sgCTRL log2 CPM differences versus mean abundance. ',
    'Dashed lines mark an absolute difference of 0.5, a descriptive magnitude threshold. ',
    '**b**, Mean differences in August and May; the diagonal marks equality. ',
    '**c**, Individual normalized expression trajectories for the four largest absolute candidate effects. Lines connect successive passages within population and period. ',
    '**d**, Six matched differences for each of the top 12 candidates; gray lines span observed minima and maxima, and black diamonds mark equal-period means. ',
    'Displayed candidates are protein-coding genes with mean CPM at least 2, absolute mean difference at least 0.5, concordant period directions and no direction reversal after any single-match deletion, ordered by absolute mean effect. ',
    'Circles denote August and triangles May; expression colors denote sgCTRL (blue) and LOY (vermilion). In a-b, colors indicate effect direction beyond the magnitude threshold; gray indicates smaller effects. ',
    'Counts use TMM exposure offsets and log2 CPM with prior count 2. Twelve libraries represent repeated sampling of one population per condition from one origin. Observed ranges are not confidence intervals. No gene-level significance test is used. Source data accompany panels a-d.')
  manifests[[id]] <- export_figure((a | b) / (c | q), id, cap,
    c(a = source_panel(id, 'a', ef), b = source_panel(id, 'b', ef), c = source_panel(id, 'c', tr),
      d = source_panel(id, 'd', merge(obs, selected[, c('gene_id', 'mean_difference', 'minimum_match', 'maximum_match')], by = 'gene_id'))),
    'a-b all eligible genes; c top4 and d top12 protein-coding stable-direction candidates with mean CPM >=2')
  # Figure 3: exploratory pathway and TF patterns plus their actual gene support.
  hg <- p$fgsea[p$fgsea$collection == 'H' & p$fgsea$contrast == 'combined', ]
  hg <- hg[order(hg$BH_across_collections, -abs(hg$NES), hg$pathway), ]
  pathways <- head(hg$pathway, 10)
  fg <- p$fgsea[p$fgsea$pathway %in% pathways, ]
  clean <- function(s) gsub('_', ' ', sub('^HALLMARK_', '', s))
  fg$label <- factor(clean(fg$pathway), levels = rev(clean(pathways)))
  a <- ggplot(fg, aes(NES, label, shape = contrast)) +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = 'gray60') +
    geom_point(aes(colour = contrast), size = 1.7) + shape_periods() +
    scale_colour_manual(values = c(AUG_2025 = '#0072B2', MAY_2026 = '#D55E00', combined = 'black'),
      labels = period_labels, name = NULL) + labs(x = 'Preranked enrichment (NES)', y = NULL)
  t0 <- p$tf[p$tf$contrast == 'combined', ]; t0 <- t0[order(t0$BH_target_pattern, -abs(t0$score), t0$source), ]
  top_tf <- head(t0$source, 10); tf <- p$tf[p$tf$source %in% top_tf, ]; tf$label <- factor(tf$source, levels = rev(top_tf))
  b <- ggplot(tf, aes(score, label, shape = contrast, colour = contrast)) +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = 'gray60') + geom_point(size = 1.7) + shape_periods() +
    scale_colour_manual(values = c(AUG_2025 = '#0072B2', MAY_2026 = '#D55E00', combined = 'black'),
      labels = period_labels, name = NULL) + labs(x = 'CollecTRI ULM score', y = NULL)
  lead <- p$leading[p$leading$pathway %in% head(pathways, 2), ]
  lead <- lead[order(-abs(lead$combined), lead$gene_id), ]; lead <- head(lead[!duplicated(lead$gene_id), ], 14)
  h <- as.data.frame(matrix_long(e$differences[lead$gene_id, , drop = FALSE], 'difference'))
  h$gene_name <- factor(lead$gene_name[match(h$gene_id, lead$gene_id)], levels = rev(lead$gene_name))
  h$match <- factor(sub('__passage_', ' P', sub('AUG_2025', 'Aug', sub('MAY_2026', 'May', h$sample_key))),
                    levels = c(paste('Aug P',1:3,sep=''),paste('May P',1:3,sep='')))
  target <- p$targets[p$targets$source %in% head(top_tf, 2), ]
  target <- target[order(-abs(target$signed_combined_effect), target$gene_id), ]
  target <- head(target[!duplicated(target$gene_id), ], 14)
  th <- as.data.frame(matrix_long(e$differences[target$gene_id, , drop = FALSE], 'difference'))
  th$target_label <- factor(paste0(target$target, ifelse(target$mor > 0, ' (+)', ' (-)'))[match(th$gene_id, target$gene_id)],
                            levels = rev(paste0(target$target, ifelse(target$mor > 0, ' (+)', ' (-)'))))
  th$regulator <- target$source[match(th$gene_id, target$gene_id)]
  th$mor <- target$mor[match(th$gene_id, target$gene_id)]
  th$signed_effect <- th$difference * th$mor
  th$match <- factor(sub('__passage_', ' P', sub('AUG_2025', 'Aug', sub('MAY_2026', 'May', th$sample_key))), levels = levels(h$match))
  lim <- max(abs(c(h$difference, th$signed_effect)))
  c <- heat_effect(h, 'gene_name', 'match', 'difference', lim, 'Log2\ndifference') + labs(title = 'Leading-edge genes')
  q <- heat_effect(th, 'target_label', 'match', 'signed_effect', lim, 'Signed log2\ndifference') +
    labs(title = paste(head(top_tf, 2), collapse = ' / '))
  id <- 'Fig03_pathways_regulators'
  cap <- paste0('**Figure 3 | Exploratory pathway and regulator patterns.** ',
    '**a**, Ten Hallmark pathways ordered by combined-signature BH-adjusted enrichment P value; symbols show combined, August and May NES. ',
    '**b**, Ten CollecTRI regulators ordered by combined-signature BH-adjusted target-pattern P value. ',
    '**c**, Up to 14 unique leading-edge genes from the first two pathways, ordered by absolute combined effect; colors show all six matched log2 differences. ',
    '**d**, Up to 14 unique targets of the first two regulators, ordered by absolute signed combined effect. Target differences are multiplied by the activating (+) or repressing (-) edge sign; their regulator assignments are in source data. ',
    'Heatmap limits are shared and symmetric; A1-A3 and M1-M3 denote August and May matches. fgsea uses all eligible uniquely mapped genes ranked by effect, sets of 15-500 represented genes, signed enrichment, and minimum reported P value 1e-10. ',
    'Primary BH correction spans all eligible H and C1-C9 sets within each contrast. CollecTRI ULM requires at least five represented targets and uses two-sided target-pattern P values, BH-adjusted across regulators within contrast. Exact statistics, coverage and leading edges accompany source data. ',
    'These conditional tests characterize gene-set or target patterns; they do not supply independent population replication. The 12 libraries repeatedly sample one sgCTRL and one LOY population from one origin. NES values are enrichment summaries, not absolute pathway activities.')
  manifests[[id]] <- export_figure((a | b) / (wrap_elements(full = ggplotGrob(c)) | wrap_elements(full = ggplotGrob(q))), id, cap,
    c(a = source_panel(id, 'a', fg), b = source_panel(id, 'b', tf), c = source_panel(id, 'c', h), d = source_panel(id, 'd', th)),
    'Top10 Hallmark and TFs by combined conditional BH P; top14 effect-ranked supporting genes from top2 sets/regulators', height = 165)
  # Figure 4: residual Y signal and baseline-weighted paralog abundance proxies.
  rd <- d$regional[d$regional$region != 'chrM', ]; rd$plot_value <- log2(rd$sum_CPM + 1)
  a <- ggplot(rd, aes(passage_number, plot_value, colour = treatment, shape = batch,
                       group = interaction(treatment, batch))) + geom_line(linewidth = 0.3) + geom_point(size = 1.4) +
    facet_wrap(~ region, ncol = 2, scales = 'free_y') + colour_groups() + shape_periods() +
    scale_x_continuous(breaks = 1:3) + labs(x = 'Passage', y = 'log2(sum normalized CPM + 1)')
  syms <- c('UTY', 'KDM5D', 'DDX3Y', 'RPS4Y1', 'EIF1AY', 'USP9Y', 'ZFY', 'TBL1Y')
  ya <- x$annotation[x$annotation$gene_name %in% syms & x$annotation$chrom == 'chrY' & !x$annotation$PAR, ]
  ya <- ya[order(match(ya$gene_name, syms)), ]
  yh <- as.data.frame(matrix_long(e$logCPM[ya$gene_id, , drop = FALSE], 'logCPM'))
  yh$gene_name <- factor(ya$gene_name[match(yh$gene_id, ya$gene_id)], levels = rev(syms))
  yh$sample <- factor(core$short[match(yh$sample_key, core$sample_key)], levels = levels(core$short))
  b <- ggplot(yh, aes(sample, gene_name, fill = logCPM)) + geom_tile() +
    scale_fill_viridis_c(option = 'cividis', name = 'log2 CPM', breaks = scales::breaks_pretty(n = 3),
      guide = guide_colourbar(barwidth = grid::unit(30, 'mm'), barheight = grid::unit(2.5, 'mm'))) + labs(x = NULL, y = NULL) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  four <- head(d$pair_definitions$X_gene_id, 4)
  xp <- as.data.frame(matrix_long(e$differences[four, , drop = FALSE], 'difference'))
  xp$gene_name <- x$annotation$gene_name[match(xp$gene_id, x$annotation$gene_id)]
  xp$batch <- x$pairs$batch[match(xp$sample_key, x$pairs$block)]
  c <- ggplot(xp, aes(difference, gene_name, shape = batch)) +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = 'gray50') + geom_point(size = 1.7) + shape_periods() +
    labs(x = 'X paralog matched log2 difference', y = NULL)
  totals <- d$pairs; totals$pair <- factor(totals$pair, levels = rev(d$pair_definitions$pair))
  q <- ggplot(totals, aes(log2_relative_total_proxy, pair, shape = batch)) +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = 'gray50') +
    geom_point(size = 1.4, alpha = 0.8, position = position_jitter(width = 0, height = 0.08, seed = 42)) +
    shape_periods() + labs(x = 'log2 relative X + Y abundance proxy', y = NULL)
  id <- 'Fig04_Y_paralog_dosage'
  cap <- paste0('**Figure 4 | Residual Y-linked RNA and paralog abundance.** ',
    '**a**, Summed normalized CPM for non-PAR Yp, non-PAR Yq, X PAR and Y PAR features, shown as log2(sum CPM + 1). ',
    'Lines connect passages within each population and period; panel scales differ. ',
    '**b**, Individual gene-level log2 CPM for eight prespecified Y genes. Sample codes follow Figure 1. ',
    '**c**, All six matched log2 differences for the X partners KDM6A, KDM5C, DDX3X and EIF1AX. ',
    '**d**, Six passage-matched X+Y relative-abundance proxies for each prespecified pair. ',
    'Each gene abundance is CPM divided by its featureCounts union-exon length in kilobases; pair ratios are (X_LOY + Y_LOY)/(X_sgCTRL + Y_sgCTRL), then log2 transformed. ',
    'The vertical line marks an unchanged total proxy. Small vertical jitter separates measurements and has no biological meaning. ',
    'All branches carry the same core TMM exposure offsets; Y genes are not independently normalized. GRCh38 cytobands define arms and pinned PAR intervals define PAR features. ',
    'Blue denotes sgCTRL, vermilion LOY; circles denote August and triangles May. Twelve libraries repeatedly sample one population per condition from one origin. ',
    'These count-derived proxies depend on annotation, isoform composition and mapping specificity. They do not establish protein compensation, Y retention, residual Y-positive cells or independent biological replication. No inferential error bars or gene-level tests are used. Source data include counts, feature lengths and denominator-support flags.')
  manifests[[id]] <- export_figure((a | b) / (c | q), id, cap,
    c(a = source_panel(id, 'a', rd), b = source_panel(id, 'b', yh), c = source_panel(id, 'c', xp), d = source_panel(id, 'd', totals)),
    'All features in regions; eight prespecified Y genes; first four X/Y pairs for X effects; all eleven pair totals', height = 165)
  # Supplementary QC: contextual PCA, period-adjusted core display, count composition and library depth.
  pc <- e$pca_context
  a <- ggplot(pc, aes(PC1, PC2, colour = treatment, shape = batch)) + geom_point(size = 1.8) +
    colour_groups() + shape_periods() + labs(x = sprintf('PC1 (%.1f%%)',100*e$pca_context_variance[1]),
      y = sprintf('PC2 (%.1f%%)',100*e$pca_context_variance[2]), title = 'All 15 libraries; contextual normalization')
  adj <- cbind(core, e$pca_adjusted$x[, 1:2]); variance <- e$pca_adjusted$sdev^2 / sum(e$pca_adjusted$sdev^2)
  b <- ggplot(adj, aes(PC1, PC2, colour = treatment, shape = batch)) + geom_point(size = 1.8) +
    colour_groups() + shape_periods() + labs(x = sprintf('PC1 (%.1f%%)',100*variance[1]),
      y = sprintf('PC2 (%.1f%%)',100*variance[2]), title = 'Period-adjusted display only')
  comp <- data.table::melt(data.table::as.data.table(e$qc), id.vars = c('sample_key','batch','treatment','passage_number'),
    measure.vars = c('mitochondrial_fraction','ribosomal_protein_fraction','top10_fraction'), variable.name='metric',value.name='fraction')
  c <- ggplot(comp, aes(passage_number, fraction*100, colour = treatment, shape = batch,
                         group = interaction(treatment,batch))) +
    geom_line(data = comp[comp$treatment != 'PARENTAL', ], linewidth=0.3) + geom_point(size=1.2) +
    facet_wrap(~metric,ncol=1,scales='free_y') + colour_groups() + shape_periods() +
    scale_x_continuous(breaks=1:3) + labs(x='Passage',y='Assigned counts (%)')
  q <- ggplot(e$qc,aes(assigned_reads/1e6,detected_genes,colour=treatment,shape=batch)) +
    geom_point(size=1.8) + colour_groups() + shape_periods() + labs(x='Assigned reads (millions)',y='Genes with at least one read')
  id <- 'SuppFig01_additional_QC'
  cap <- paste0('**Supplementary Figure 1 | Contextual samples and count composition.** ',
    '**a**, PCA of all 15 libraries using its own TMM normalization, with three parental libraries in gray. ',
    '**b**, Core PCA after removal of an additive sampling-period term while preserving treatment; these values are used only for visualization. ',
    '**c**, Mitochondrial, ribosomal-protein and top-ten-gene fractions of assigned counts. Lines follow passage order within period and group. ',
    '**d**, Detected genes versus assigned read depth. PCA uses the core eligible gene set, centered log2 CPM with prior count 2 and no variance scaling. ',
    'Blue denotes sgCTRL, vermilion LOY, circles August and triangles May. Twelve core libraries represent one established population per condition sampled three times per period across two periods, from one generation history. ',
    'Parental origins and matching are not established; they do not enter core effects. Period removal can subtract biological drift as well as technical differences and is not used in downstream effect estimates. Count composition is descriptive and supplies no automatic sample-exclusion rule. Source data accompany panels a-d.')
  manifests[[id]] <- export_figure((a|b)/(c|q),id,cap,
    c(a=source_panel(id,'a',pc),b=source_panel(id,'b',adj),c=source_panel(id,'c',comp),d=source_panel(id,'d',e$qc)),
    'All 15 contextual samples or all 12 core samples; no expression-outcome exclusions')
  # Supplementary sensitivities and temporally held-out signature.
  a <- ggplot(ef,aes(mean_difference,RLE_mean_difference)) + geom_point(size=0.3,alpha=0.4,colour='#555555') +
    geom_abline(slope=1,intercept=0,linewidth=0.25) + labs(x='TMM mean log2 difference',y='RLE mean log2 difference')
  b <- ggplot(ef,aes(mean_difference,prior025_mean_difference)) + geom_point(size=0.3,alpha=0.4,colour='#555555') +
    geom_abline(slope=1,intercept=0,linewidth=0.25) + labs(x='Prior count 2: mean difference',y='Prior count 0.25: mean difference')
  c <- ggplot(ef,aes(mean_logCPM,LOO_max_abs_change)) + geom_point(size=0.3,alpha=0.4,colour='#555555') +
    labs(x='Mean log2 CPM',y='Maximum absolute change after\nomitting one passage match')
  sc <- d$signature_scores
  q <- ggplot(sc,aes(passage_number,signed_percentile_score,colour=treatment,shape=batch,group=interaction(treatment,batch))) +
    geom_line(linewidth=0.4) + geom_point(size=1.7) + facet_wrap(~batch,ncol=2,labeller=as_labeller(period_labels)) +
    colour_groups() + shape_periods() + scale_x_continuous(breaks=1:3) +
    labs(x='Passage',y='August-derived signed rank score')
  id <- 'SuppFig02_sensitivity_signature'
  cap <- paste0('**Supplementary Figure 2 | Effect sensitivity and an August-derived signature.** ',
    '**a**, TMM versus RLE normalization on the same expressed background. ',
    '**b**, Log2 CPM prior count 2 versus 0.25 under identical TMM offsets. Diagonals mark equal estimates. ',
    '**c**, Largest absolute change in the equal-period mean after deleting any single passage match; the two periods retain equal weight after deletion. ',
    '**d**, A signature selected using August alone, with ', sum(d$signature$direction>0),' higher and ',sum(d$signature$direction<0),
    ' lower genes. August-specific filtering, normalization and symbol selection precede selection by absolute mean difference (at least 0.5), with all three passage differences in the same direction and a maximum of 50 genes per direction. ',
    'The scoring background is frozen from August. Each sample is scored independently as mean within-sample percentile rank of higher genes minus that of lower genes. May observations never determine signature membership or scoring background. ',
    'Blue denotes sgCTRL, vermilion LOY, circles August and triangles May. Lines follow passage order. Twelve libraries repeatedly sample the same two established populations from one generation history. May therefore assesses temporal persistence, not validation in independently generated populations. All panels are descriptive. Full gene/filter and signature source tables accompany the figure.')
  manifests[[id]] <- export_figure((a|b)/(c|q),id,cap,
    c(a=source_panel(id,'a',ef),b=source_panel(id,'b',ef),c=source_panel(id,'c',ef),d=source_panel(id,'d',sc)),
    'All primary eligible genes; signature uses August-only discovery and a frozen scoring universe')
  # Supplementary score diagnostics, with all matched values.
  diag_sets <- c('HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_INTERFERON_ALPHA_RESPONSE',
                 'HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_HYPOXIA',
                 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_OXIDATIVE_PHOSPHORYLATION')
  hd <- as.data.frame(matrix_long(p$hallmark_differences[diag_sets,,drop=FALSE],'difference','pathway'))
  hd$pathway <- factor(clean(hd$pathway),levels=rev(clean(diag_sets)))
  hd$match <- sub('__passage_',' P',sub('AUG_2025','Aug',sub('MAY_2026','May',hd$sample_key)))
  a <- heat_effect(hd,'pathway','match','difference',max(abs(hd$difference)),'Rank-score\ndifference')
  pr <- as.data.frame(matrix_long(p$progeny_differences,'difference','pathway'))
  pr$match <- sub('__passage_',' P',sub('AUG_2025','Aug',sub('MAY_2026','May',pr$sample_key)))
  b <- heat_effect(pr,'pathway','match','difference',max(abs(pr$difference)),'Weighted log2\ndifference')
  id <- 'SuppFig03_program_scores'
  cap <- paste0('**Supplementary Figure 3 | Matched program-score diagnostics.** ',
    '**a**, LOY-minus-sgCTRL differences in mean within-sample percentile ranks for eight prespecified Hallmark programs. ',
    'The rank background contains all core eligible, uniquely mapped genes; scores are computed independently within each sample. ',
    '**b**, Matched PROGENy footprint differences using the human top-100 model. Represented gene weights are divided by their absolute sum within each footprint, ',
    'then applied to normalized log2 expression differences. Thus these are weighted expression summaries, not standardized PROGENy activity units. ',
    'All six passage/day matches are shown, with separate symmetric color scales for the two scoring methods. ',
    'The 12 libraries represent repeated measurements of one sgCTRL and one LOY population from one shared nucleofection origin. ',
    'No group-level significance tests or inferential intervals are used. Cell-cycle, interferon and stress signatures are biological diagnostics and were not automatically regressed out. ',
    'Exact gene membership, footprint weights, target coverage and matched scores are supplied with the source tables.')
  manifests[[id]] <- export_figure(a/b,id,cap,c(a=source_panel(id,'a',hd),b=source_panel(id,'b',pr)),
    'Eight prespecified Hallmark diagnostic programs; all 14 PROGENy footprints',height=160)
  manifests <- c(manifests, enrichment_figures(p), zmynd8_figure(readRDS('outputs/models/04_ZMYND8_followup.rds')))
  manifests <- c(manifests, chromatin_figures(readRDS('outputs/models/04_chromatin_followup.rds')))
  manifests <- c(manifests, hdac12_figures(readRDS('outputs/models/04_HDAC12_gene_drivers.rds')))
  manifests <- c(manifests, hdac_external_figures(readRDS('outputs/models/04_HDACEXT_comparison.rds'),
    readRDS('outputs/models/04_HDACEXT_source_models.rds'), e))
  manifests <- c(manifests, epigenetic_figures(readRDS('outputs/models/04_epigenetic_deep_dive.rds'), x))
  manifests <- c(manifests, chromatin_gsea_figures(readRDS('outputs/models/05_chromatin_GSEA_curves.rds')))
  manifest <- do.call(rbind,manifests); rownames(manifest)<-NULL
  data.table::fwrite(manifest,'outputs/figures/figure_manifest.csv')
  captions <- vapply(unique(manifest$caption),function(f) paste(readLines(f,warn=FALSE),collapse='\n'),character(1))
  writeLines(paste(captions,collapse='\n\n'),'outputs/captions/figure_captions.md',useBytes=TRUE)
  invisible(manifest)
}
