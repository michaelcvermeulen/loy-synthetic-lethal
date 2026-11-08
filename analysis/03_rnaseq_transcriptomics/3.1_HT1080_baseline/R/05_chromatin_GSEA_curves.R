# Draw the exact saved preranked enrichment statistic; never rerun permutation tests.
chromatin_gsea_curve_data <- function(p, e, cfg) {
  ranked <- rank_effect(p$effects_symbol[, cfg$contrast], rownames(p$effects_symbol))
  all_sets <- do.call(c, unname(p$sets)); n <- length(ranked)
  assert(n == 25461L && !anyDuplicated(names(ranked)), 'Original GSEA universe changed')
  definitions <- do.call(rbind, lapply(names(cfg$figures), function(id) {
    a <- cfg$figures[[id]]
    do.call(rbind, lapply(seq_along(a$sets), function(j)
      data.frame(figure = id, panel = letters[j], as.data.frame(a$sets[[j]]))))
  }))
  curves <- list(); summary <- list()
  for (j in seq_len(nrow(definitions))) {
    def <- definitions[j, ]; name <- def$pathway
    saved <- p$fgsea[p$fgsea$pathway == name & p$fgsea$contrast == cfg$contrast, ]
    assert(nrow(saved) == 1 && is.finite(saved$ES), paste('Saved GSEA result missing:', name))
    hits <- names(ranked) %in% all_sets[[name]]; nh <- sum(hits)
    assert(nh == saved$size && nh > 0 && nh < n, paste('GSEA membership mismatch:', name))
    weights <- abs(ranked)^cfg$gsea_parameter
    increments <- ifelse(hits, weights / sum(weights[hits]), -1 / (n - nh))
    running <- c(0, cumsum(increments))
    peak <- if (saved$ES > 0) which.max(running) - 1L else which.min(running) - 1L
    assert(abs(running[peak + 1L] - saved$ES) < 1e-10, paste('Running ES differs from original:', name))
    assert(abs(tail(running, 1)) < 1e-10, paste('Running score does not end at zero:', name))
    le <- hits & if (saved$ES > 0) seq_len(n) <= peak else seq_len(n) > peak
    original_le <- strsplit(saved$leadingEdge, ';', fixed = TRUE)[[1]]
    assert(setequal(names(ranked)[le], original_le), paste('Leading edge differs from original:', name))
    # Cross-check with the installed fgsea statistic implementation, including leading-edge ranks.
    check <- fgsea::calcGseaStat(ranked, which(hits), gseaParam = cfg$gsea_parameter,
      scoreType = cfg$score_type, returnLeadingEdge = TRUE)
    assert(abs(check$res - saved$ES) < 1e-10, paste('fgsea statistic check failed:', name))
    assert(setequal(names(ranked)[check$leadingEdge], original_le), paste('fgsea leading edge check failed:', name))
    curves[[name]] <- data.frame(pathway = name, rank = 0:n,
      gene_name = c(NA_character_, names(ranked)),
      gene_id = c(NA_character_, e$mapping$gene_id[match(names(ranked), e$mapping$gene_name)]),
      ranking_metric = c(NA_real_, unname(ranked)), hit = c(FALSE, hits),
      leading_edge = c(FALSE, le), increment = c(0, unname(increments)),
      running_ES = running, peak_rank = peak)
    cc <- p$camera[p$camera$pathway == name & p$camera$contrast == cfg$contrast, ]
    cs <- p$camera_sensitivity[p$camera_sensitivity$pathway == name & p$camera_sensitivity$contrast == cfg$contrast, ]
    summary[[name]] <- data.frame(def, contrast = cfg$contrast, rank_genes = n,
      genes = nh, leading_genes = sum(le), ES = saved$ES, NES = saved$NES,
      nominal_P = saved$pval, original_global_BH = saved$BH_across_collections,
      camera_fixed_global_BH = cc$BH_across_collections, camera_estimated_global_BH = cs$BH_across_collections,
      peak_rank = peak, peak_ranking_metric = ranked[max(1, peak)],
      gsea_parameter = cfg$gsea_parameter, ES_agreement_error = abs(check$res - saved$ES),
      leading_edge_exact_agreement = TRUE)
  }
  summary <- do.call(rbind, summary); rownames(summary) <- NULL
  tab(summary, 'GSEA_CURVES_selected_signatures')
  tab(data.frame(rank = seq_len(n), gene_name = names(ranked),
    gene_id = e$mapping$gene_id[match(names(ranked), e$mapping$gene_name)],
    combined_log2_difference = unname(ranked)), 'GSEA_CURVES_frozen_ranking')
  list(cfg = cfg, definitions = definitions, curves = curves, summary = summary, ranked = ranked)
}

chromatin_gsea_panel <- function(curve, result, limits) {
  n <- result$rank_genes; peak <- result$peak_rank
  colour <- if (result$ES > 0) '#D55E00' else '#0072B2'
  region <- if (result$ES > 0) c(0, peak) else c(peak, n)
  common_x <- function() scale_x_continuous(limits = c(0, n), expand = expansion(mult = 0),
    breaks = c(0, 5000, 10000, 15000, 20000, 25000), labels = c('0','5k','10k','15k','20k','25k'))
  bh <- if (result$original_global_BH < 0.001) formatC(result$original_global_BH, format = 'e', digits = 2) else
    formatC(result$original_global_BH, format = 'f', digits = 3)
  label <- sprintf('NES %+.2f   |   global BH %s   |   n = %d', result$NES, bh, result$genes)
  a <- ggplot(curve, aes(rank, running_ES)) +
    annotate('rect', xmin = region[1], xmax = region[2], ymin = -Inf, ymax = Inf, fill = colour, alpha = 0.07) +
    geom_hline(yintercept = 0, colour = 'gray50', linewidth = 0.25) +
    geom_line(colour = colour, linewidth = 0.55) +
    annotate('segment', x = peak, xend = peak, y = 0, yend = result$ES,
      colour = colour, linetype = 'dashed', linewidth = 0.3) +
    annotate('point', x = peak, y = result$ES, colour = colour, size = 1.4) +
    common_x() + scale_y_continuous(limits = limits, breaks = c(-0.6,-0.3,0,0.3), expand = expansion(mult = 0)) +
    labs(x = NULL, y = 'Running ES', title = result$title, subtitle = paste(result$context, label, sep = '\n')) +
    nature_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.line.x = element_blank(), plot.title = element_text(size = 7, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 5.8, lineheight = 1.1, margin = margin(b = 4)),
      plot.margin = margin(3, 4, 0, 3))
  hit <- curve[curve$hit, ]
  b <- ggplot(hit, aes(x = rank, xend = rank)) +
    geom_segment(aes(y = 0, yend = 1), colour = '#B8B8B8', linewidth = 0.15) +
    geom_segment(data = hit[hit$leading_edge, ], aes(y = 0, yend = 1), colour = colour, linewidth = 0.2) +
    common_x() + scale_y_continuous(limits = c(0, 1), breaks = NULL, expand = expansion(mult = 0)) +
    labs(x = NULL, y = 'Hits') + nature_theme() + theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
      axis.line = element_blank(), axis.title.y = element_text(size = 5.5, angle = 0, vjust = 0.5),
      plot.margin = margin(0, 4, 0, 3))
  metric <- curve[!is.na(curve$ranking_metric), ]
  c <- ggplot(metric, aes(rank, ranking_metric)) +
    geom_area(data = metric[metric$ranking_metric >= 0, ], fill = '#D55E00', alpha = 0.65, linewidth = 0) +
    geom_area(data = metric[metric$ranking_metric <= 0, ], fill = '#0072B2', alpha = 0.65, linewidth = 0) +
    geom_hline(yintercept = 0, colour = 'gray50', linewidth = 0.2) + common_x() +
    scale_y_continuous(breaks = c(-5,0,5), expand = expansion(mult = c(0.03,0.03))) +
    labs(x = 'Gene rank: higher in LOY to lower in LOY', y = 'Log2\neffect') + nature_theme() +
    theme(axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 5.5),
      axis.text = element_text(size = 5.5), plot.margin = margin(0, 4, 3, 3))
  a / b / c + plot_layout(heights = c(4.6, 0.45, 1.3))
}

chromatin_gsea_figures <- function(z) {
  cfg <- z$cfg; panels <- list(); individual <- list(); manifests <- list()
  directory <- 'outputs/figures/chromatin_GSEA_individual'
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (j in seq_len(nrow(z$summary))) {
    result <- z$summary[j, ]; name <- result$pathway; curve <- z$curves[[name]]
    plot <- chromatin_gsea_panel(curve, result, unlist(cfg$running_ES_limits))
    panels[[name]] <- plot
    source <- source_panel(result$figure, result$panel, curve)
    stem <- file.path(directory, name)
    ggsave(paste0(stem, '.pdf'), plot, device = grDevices::cairo_pdf,
      width = 89, height = 78, units = 'mm', family = 'Arial', bg = 'white')
    ggsave(paste0(stem, '.svg'), plot, device = svglite::svglite,
      width = 89, height = 78, units = 'mm', bg = 'white')
    ggsave(paste0(stem, '.png'), plot, device = ragg::agg_png,
      width = 89, height = 78, units = 'mm', dpi = 600, bg = 'white')
    individual[[name]] <- data.frame(pathway = name, title = result$title, source_context = result$context,
      figure = result$figure, panel = result$panel, source_data = source,
      pdf = paste0(stem, '.pdf'), svg = paste0(stem, '.svg'), png = paste0(stem, '.png'),
      width_mm = 89, height_mm = 78, png_dpi = 600,
      pdf_sha256 = hash_file(paste0(stem, '.pdf')), svg_sha256 = hash_file(paste0(stem, '.svg')),
      png_sha256 = hash_file(paste0(stem, '.png')))
  }
  for (id in names(cfg$figures)) {
    def <- cfg$figures[[id]]; rows <- z$summary[z$summary$figure == id, ]; sources <- character()
    # Wrap each three-track GSEA plot as one tagged panel while preserving Arial metrics.
    plots <- lapply(rows$pathway, function(name) wrap_elements(full = patchworkGrob(panels[[name]])))
    for (j in seq_len(nrow(rows))) sources[rows$panel[j]] <- file.path('outputs/source_data', paste0(id, '_', rows$panel[j], '.csv'))
    caption <- paste0('**Figure ', def$number, ' | ', def$title, '.** ',
      if (def$number == 12) paste0('**a-d**, Source RNA responses to EED, EZH2 or SUZ12 knockdown in human TIG3 fibroblasts; both EED response arms and the weaker SUZ12 result are retained. ',
        '**e,f**, Genes annotated with H3K27me3 at high-CpG promoters in mouse MEF and partially reprogrammed MCV6 source cells. ') else
        paste0('**a,b**, HDAC2 and HDAC1/2 knockdown-reduced genes. **c**, BMI1/MEL18 knockdown-induced genes. ',
          '**d**, SMARCE1 re-expression-induced genes. **e**, Promoters gaining H3K27me3 in human liver tumors; this UP suffix describes the mark, not RNA induction. ',
          '**f**, Suz12-deficiency-reduced genes in mouse embryonic stem cells. '),
      'Each panel shows the exact weighted running enrichment score (ES; exponent 1), all set-member hit positions and the ranked log2 expression differences. ',
      'The 25,461 representative genes are ordered by combined LOY-minus-sgCTRL effect, from higher on the left to lower on the right. ',
      'The colored curve, shaded interval and colored hits identify enrichment direction and leading-edge region; gray hits mark remaining members. A dot and dashed segment locate the ES extremum. ',
      'All panels share ES and ranked-effect scales. Annotations reproduce original normalized enrichment scores (NES), global BH values across 23,696 tested sets and represented-gene counts. ',
      'The combined effect equally weights August and May; 12 libraries repeatedly sample one established population per condition. Display selection is post hoc; no new enrichment tests were run. ',
      'Estimated-correlation CAMERA yields no combined global-BH discoveries. Source signatures do not measure HT1080 enzyme activity or chromatin marks.')
    manifests[[id]] <- export_figure(wrap_plots(plots, ncol = 2, byrow = TRUE), id, caption, sources,
      cfg$selection, height = 240)
    manifests[[id]]$analysis_objects <- paste(manifests[[id]]$analysis_objects, 'outputs/models/05_chromatin_GSEA_curves.rds', sep = ';')
    manifests[[id]]$analysis_config <- paste(manifests[[id]]$analysis_config, 'config/chromatin_GSEA_curves.yml', sep = ';')
    manifests[[id]]$statistical_scope <- cfg$statistics
  }
  individual <- do.call(rbind, individual); rownames(individual) <- NULL
  individual$caption <- paste0('outputs/captions/', individual$figure, '.caption.md')
  data.table::fwrite(individual, file.path(directory, 'manifest.csv'))
  manifests
}
