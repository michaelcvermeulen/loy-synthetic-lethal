# Exploratory gene intersections from the existing, unmodified primary DE fit.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(edgeR)
})
roots <- c(getwd(), file.path(getwd(), 'THESIS_ADDBACK_2026'))
root <- roots[file.exists(file.path(roots, 'outputs/fit_primary.rds'))][1]
if (is.na(root)) stop('Run from the add-back directory or THESIS_ADDBACK_2026')
setwd(root)
outdir <- 'outputs/responsive_genes'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rd <- function(cn) fread(file.path('outputs/de', paste0('de_', cn, '.csv')), data.table = FALSE)
base <- rd('LOY_v_LXHA')
stopifnot(!anyDuplicated(base$gene_id))
all <- base[, c('gene_id', 'gene_name', 'chrom')]
contrasts <- c(LOY = 'LOY_v_LXHA', UTY = 'UTY_v_LXHA', DUAL = 'DUAL_v_LXHA',
               KDM5D = 'KDM5D_v_LXHA', UTY_vs_LOY = 'UTY_v_LOY',
               DUAL_vs_LOY = 'DUAL_v_LOY', DUAL_vs_UTY = 'DUAL_v_UTY')
fields <- c('logFC', 'CI_L', 'CI_R', 'FDR', 'FDR_treat_0.5', 'selected_0.5', 'treat_0.5')
for (nm in names(contrasts)) {
  x <- rd(contrasts[[nm]])
  stopifnot(!anyDuplicated(x$gene_id), setequal(x$gene_id, all$gene_id))
  x <- x[match(all$gene_id, x$gene_id), fields]
  names(x) <- paste0(nm, '_', fields)
  all <- cbind(all, x)
}
stopifnot(!anyNA(all))
all$baseline_small_nonsignificant <- all$LOY_FDR >= 0.05 & abs(all$LOY_logFC) < 0.5
all$baseline_CI_within_margin <- all$LOY_CI_L >= -0.5 & all$LOY_CI_R <= 0.5
all$response_either <- all$UTY_selected_0.5 | all$DUAL_selected_0.5
all$response_both <- all$UTY_selected_0.5 & all$DUAL_selected_0.5
all$response_class <- ifelse(all$response_both, 'Both',
  ifelse(all$UTY_selected_0.5, 'UTY only', ifelse(all$DUAL_selected_0.5, 'Combined only', 'Neither')))
all$same_direction_in_both <- all$response_both & sign(all$UTY_logFC) == sign(all$DUAL_logFC)
sig <- fread('outputs/signature/frozen_signature_mapping.csv', data.table = FALSE)
all$in_August_signature <- all$gene_id %in% sig$gene_id_stable
# These are alternative baseline filters; CI containment is not required to
# have FDR >= 0.05, because a precisely measured small difference may be nonzero.
filters <- list(small_nonsignificant = all$baseline_small_nonsignificant,
                CI_within_margin = all$baseline_CI_within_margin)
summary <- list()
for (key in names(filters)) {
  eligible <- filters[[key]]
  for (method in c('descriptive', 'TREAT')) {
    u <- if (method == 'descriptive') all$UTY_selected_0.5 else all$UTY_treat_0.5
    d <- if (method == 'descriptive') all$DUAL_selected_0.5 else all$DUAL_treat_0.5
    k <- if (method == 'descriptive') all$KDM5D_selected_0.5 else all$KDM5D_treat_0.5
    summary[[paste(key, method)]] <- data.frame(baseline_filter = key, response_test = method,
      eligible = sum(eligible), UTY = sum(eligible & u), DUAL = sum(eligible & d),
      either = sum(eligible & (u | d)), both = sum(eligible & u & d),
      both_same_direction = sum(eligible & u & d & sign(all$UTY_logFC) == sign(all$DUAL_logFC)),
      KDM5D = sum(eligible & k), August_overlap = sum(eligible & (u | d) & all$in_August_signature))
    for (arm in c('UTY', 'DUAL', 'either', 'both')) {
      hit <- switch(arm, UTY = u, DUAL = d, either = u | d, both = u & d)
      z <- all[eligible & hit, ]
      z <- z[order(pmin(z$UTY_FDR, z$DUAL_FDR), z$gene_id), ]
      fwrite(z, file.path(outdir, paste0(key, '_', method, '_', arm, '.csv')))
    }
  }
}
summary <- do.call(rbind, summary)
rownames(summary) <- NULL
fwrite(summary, file.path(outdir, 'selection_summary.csv'))
fwrite(all, file.path(outdir, 'all_genes_classification.csv'))
stopifnot(nrow(all) == 23992,
  with(summary[summary$response_test == 'descriptive', ],
       all(UTY == c(363, 53)) && all(DUAL == c(552, 76)) && all(either == c(657, 91)) && all(both == c(258, 38))))

theme_set(theme_bw(base_size = 10))
save_plot <- function(p, stem, width, height) {
  for (ext in c('png', 'pdf', 'svg')) {
    args <- list(filename = file.path('figures', paste0(stem, '.', ext)), plot = p,
                 width = width, height = height, units = 'in', bg = 'white')
    if (ext == 'png') args$dpi <- 300
    if (ext == 'pdf') args$device <- cairo_pdf
    do.call(ggsave, args)
  }
}
count_data <- rbindlist(lapply(seq_len(nrow(summary)), function(i) {
  z <- summary[i, ]
  data.frame(baseline_filter = z$baseline_filter, response_test = z$response_test,
    category = c('UTY only', 'Both', 'Combined only'), n = c(z$UTY-z$both, z$both, z$DUAL-z$both))
}))
count_data$baseline_filter <- factor(count_data$baseline_filter,
  levels = names(filters), labels = c('Small, non-significant baseline difference', 'Baseline 95% CI within +/-0.5'))
count_data$response_test <- factor(count_data$response_test, c('descriptive', 'TREAT'))
count_data$category <- factor(count_data$category, c('UTY only', 'Both', 'Combined only'))
p <- ggplot(count_data, aes(category, n, fill = category)) + geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
  facet_grid(response_test ~ baseline_filter, scales = 'free_y') +
  scale_fill_manual(values = c('UTY only' = '#0072B2', Both = '#009E73', 'Combined only' = '#D55E00')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = 'Genes', title = 'Add-back responses with limited LOY-vector differences',
       subtitle = 'Disjoint response categories; descriptive filter and TREAT sensitivity') +
  theme(legend.position = 'none', panel.grid.major.x = element_blank())
save_plot(p, 'student_responsive_counts', 10, 6)
fwrite(count_data, file.path(outdir, 'figure_counts.csv'))

# Plot every CI-contained candidate across all 18 individual libraries.
dge <- readRDS('outputs/dge_primary.rds')
samples <- fread('outputs/sample_manifest.csv', data.table = FALSE)
stopifnot(identical(colnames(dge), samples$sample_id))
expression <- cpm(dge, log = TRUE, prior.count = 2)
candidate <- all[all$baseline_CI_within_margin & all$response_either, ]
candidate <- candidate[order(candidate$response_class, candidate$UTY_logFC + candidate$DUAL_logFC), ]
groups <- c('CTRL', 'LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
labels <- c('sgCTRL', 'LOY', 'LXHA', 'UTY', 'KDM5D', 'UTY+KDM5D')
sample_order <- order(match(samples$final_group, groups), samples$replicate_block)
samples <- samples[sample_order, ]
e <- expression[candidate$gene_id, samples$sample_id, drop = FALSE]
centered <- sweep(e, 1, rowMeans(e[, samples$final_group == 'LOY_LXHA', drop = FALSE]), '-')
heat <- data.frame(gene_id = rep(rownames(e), ncol(e)), sample_id = rep(colnames(e), each = nrow(e)),
                   logCPM = as.vector(e), difference_from_vector_mean = as.vector(centered))
heat$gene_label <- paste(candidate$gene_name, candidate$gene_id)[match(heat$gene_id, candidate$gene_id)]
heat$gene_label <- factor(heat$gene_label, levels = rev(paste(candidate$gene_name, candidate$gene_id)))
heat$group <- factor(samples$final_group[match(heat$sample_id, samples$sample_id)], groups, labels)
heat$occasion <- factor(samples$replicate_block[match(heat$sample_id, samples$sample_id)])
fwrite(heat, file.path(outdir, 'figure_sample_expression.csv'))
p <- ggplot(heat, aes(occasion, gene_label, fill = difference_from_vector_mean)) +
  geom_tile() + facet_grid(. ~ group) +
  scale_fill_gradient2(low = '#2166AC', mid = 'white', high = '#B2182B', midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish, name = 'log2 CPM\nminus LXHA mean') +
  labs(x = 'Matched collection occasion', y = NULL,
       title = 'All 91 candidates with baseline confidence intervals within +/-0.5',
       subtitle = 'Observed normalized expression; each gene centered on its LXHA mean; colors clipped at +/-3') +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 5.5),
        strip.text = element_text(size = 9), legend.position = 'bottom')
save_plot(p, 'student_responsive_expression', 10, 16)

# All candidates, sorted from most positive to most negative UTY effect.
top <- candidate[order(-candidate$UTY_logFC, candidate$gene_id), ]
# Compact panel: retain every candidate, average the three collections per population.
panel_mean <- sapply(groups, function(g) rowMeans(expression[top$gene_id,
  samples$sample_id[samples$final_group == g], drop = FALSE]))
panel_delta <- sweep(panel_mean, 1, panel_mean[, 'LOY_LXHA'], '-')
panel_data <- data.frame(gene_id = rep(top$gene_id, length(groups)),
  gene_name = rep(top$gene_name, length(groups)), UTY_rank = rep(seq_len(nrow(top)), length(groups)),
  group = rep(groups, each = nrow(top)), mean_logCPM = as.vector(panel_mean),
  difference_from_vector_mean = as.vector(panel_delta))
panel_data$gene_id <- factor(panel_data$gene_id, levels = rev(top$gene_id))
panel_data$group <- factor(panel_data$group, groups, c('sgCTRL', 'LOY', 'LXHA', 'UTY', 'KDM5D', 'UTY + KDM5D'))
up_n <- sum(top$UTY_logFC >= 0)
direction_labels <- c(sprintf('UTY higher (%d)', up_n), sprintf('UTY lower (%d)', nrow(top) - up_n))
panel_data$UTY_direction <- factor(ifelse(top$UTY_logFC[panel_data$UTY_rank] >= 0,
  direction_labels[1], direction_labels[2]), levels = direction_labels)
fwrite(panel_data, file.path(outdir, 'figure_compact_group_means.csv'))
panel <- ggplot(panel_data[panel_data$group != 'LXHA', ], aes(group, gene_id, fill = difference_from_vector_mean)) +
  geom_tile() + facet_grid(UTY_direction ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient2(low = '#2166AC', mid = 'white', high = '#B2182B', midpoint = 0,
    limits = c(-2, 2), breaks = c(-2, 0, 2), oob = scales::squish,
    name = 'Mean log2 CPM\nrelative to LXHA') +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = '91 genes ordered by UTY log2 fold change', title = 'Add-back-responsive candidates') +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = '#222222'),
    axis.title.y = element_text(size = 8), plot.title = element_text(size = 10, face = 'bold'),
    strip.text.y = element_text(size = 8), legend.position = 'bottom',
    legend.title = element_text(size = 8), legend.text = element_text(size = 7),
    legend.key.height = grid::unit(3, 'mm'), legend.key.width = grid::unit(10, 'mm'))
save_plot(panel, 'student_responsive_compact_heatmap', 3.4, 4.3)
stopifnot(nrow(panel_data) == 546, max(abs(panel_delta[, 'LOY_LXHA'])) == 0)
effect <- rbindlist(lapply(c('LOY', 'UTY', 'KDM5D', 'DUAL'), function(arm) {
  data.frame(gene_id = top$gene_id, gene_label = paste(top$gene_name, top$gene_id), comparison = arm,
             UTY_rank = seq_len(nrow(top)), page = ceiling(seq_len(nrow(top)) / 24),
             logFC = top[[paste0(arm, '_logFC')]], CI_L = top[[paste0(arm, '_CI_L')]],
             CI_R = top[[paste0(arm, '_CI_R')]])
}))
effect$gene_label <- factor(effect$gene_label, levels = rev(paste(top$gene_name, top$gene_id)))
effect$comparison <- factor(effect$comparison, c('LOY', 'UTY', 'KDM5D', 'DUAL'),
  c('LOY vs LXHA', 'UTY vs LXHA', 'KDM5D vs LXHA', 'Combined vs LXHA'))
fwrite(effect, file.path(outdir, 'figure_all91_effects.csv'))
# Remove the superseded generated source table from the analysis and handoff.
for (old in c(file.path(outdir, 'figure_top12_effects.csv'),
              'student_report/data/outputs__responsive_genes__figure_top12_effects.csv'))
  if (file.exists(old)) unlink(old)
plot_effects <- function(z, title) {
  ggplot(z, aes(logFC, gene_label, color = comparison)) +
    annotate('rect', xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = 'grey93') +
    geom_vline(xintercept = 0, color = 'grey60') +
    geom_errorbar(aes(xmin = CI_L, xmax = CI_R), orientation = 'y', width = 0.4) + geom_point(size = 1.7) +
    facet_grid(. ~ comparison) +
    scale_color_manual(values = c('#555555', '#0072B2', '#009E73', '#D55E00')) +
    coord_cartesian(xlim = range(c(effect$CI_L, effect$CI_R))) +
    labs(x = 'Model log2 fold change and pointwise 95% confidence interval', y = NULL,
         title = title,
         subtitle = 'Ordered by UTY vs LXHA log2 fold change (highest to lowest); shaded region is +/-0.5') +
    theme(legend.position = 'none', axis.text.y = element_text(size = 8))
}
save_plot(plot_effects(effect, 'All 91 candidates with precise, small baseline differences'),
          'student_responsive_effects', 16, 23)
page_stems <- sprintf('student_responsive_effects_page%02d', 1:4)
for (i in 1:4) {
  z <- effect[effect$page == i, ]
  title <- sprintf('Add-back responses: UTY-ranked genes %d-%d of 91', min(z$UTY_rank), max(z$UTY_rank))
  save_plot(plot_effects(z, title), page_stems[i], 16, 8)
}
stopifnot(nrow(effect) == 364, all(diff(top$UTY_logFC) <= 0),
          all(table(effect$gene_id) == 4))

methods <- c(
  'Exploratory secondary analysis of the existing collection-occasion-blocked limma-voom fit; no refitting.',
  'Universe: all 23,992 primary tested genes; joins use unique stable gene IDs.',
  'Small/non-significant baseline: LOY-v-LXHA FDR >=0.05 AND abs(logFC)<0.5.',
  'CI-contained baseline: LOY-v-LXHA pointwise 95% CI entirely within [-0.5,0.5]. This is an alternative filter, not necessarily a subset of the non-significant filter.',
  'The 0.5 margin was chosen for this exploratory analysis to match the existing effect-size threshold; it is not a prospectively validated equivalence margin.',
  'Descriptive response: arm-v-LXHA FDR<0.05 AND abs(logFC)>=0.5. TREAT sensitivity: existing arm-v-LXHA TREAT FDR<0.05 at threshold 0.5.',
  'Original FDR values were calculated across all tested genes within each contrast, before this intersection. No FDR claim applies to the final combined selection rule.',
  'Pointwise CI containment is not multiplicity-adjusted equivalence testing. Nonsignificance alone does not prove similarity.',
  'LOY-v-LXHA and add-back-v-LXHA estimates share the LXHA comparator; selection and response estimates are statistically dependent.',
  'UTY only / combined only means passing the threshold only in that comparison; it does not establish different UTY and combined effects. Direct DUAL-v-UTY estimates are exported.',
  'Both counts are intersections irrespective of sign; both_same_direction separately records concordance.',
  'These are population-associated response candidates, not causal UTY targets or rescued LOY genes. Each condition is one established population sampled three times.',
  'Heatmap: all CI-contained descriptive candidates; TMM log2 CPM (prior.count=2), centered per gene on LXHA sample mean, unadjusted for occasion. Source values are unclipped.',
  'Effect figures: all 91 CI-contained descriptive candidates ordered by decreasing UTY-v-LXHA logFC, ties by stable gene ID; LOY, UTY, KDM5D and combined versus LXHA; original pointwise intervals. Four pages of up to 24 genes use the same x-axis limits, plus a full-length export.',
  'August signature membership is annotated but is not used in selecting candidates.'
)
methods <- c(methods, 'Compact heatmap: all 91 candidates, mean TMM log2 CPM across three collections per population minus the gene-specific LXHA mean. The zero-valued LXHA reference column is omitted from the plot but retained in the source table. Rows ordered by decreasing UTY model logFC and split by its sign; no row scaling or occasion adjustment. Colors clipped at +/-2; source table retains unclipped values. UTY higher/lower labels describe effect sign, not separate significance calls.')
writeLines(methods, file.path(outdir, 'METHODS.txt'))
inputs <- c(file.path('outputs/de', paste0('de_', contrasts, '.csv')), 'outputs/dge_primary.rds',
            'outputs/sample_manifest.csv', 'outputs/signature/frozen_signature_mapping.csv', 'code/run_responsive_genes.R')
fwrite(data.frame(file = inputs, sha256 = vapply(inputs, digest::digest, character(1), file = TRUE, algo = 'sha256')),
       file.path(outdir, 'input_manifest.csv'))
writeLines(c('PASS: one-to-one joins for every contrast; 23,992 genes; no missing fields.',
             'PASS: descriptive counts reproduce 363/552/657/258 and 53/76/91/38.',
             'PASS: expression samples align with the frozen sample manifest.',
             capture.output(summary)), file.path(outdir, 'VERIFICATION.txt'))
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))

stems <- c('student_responsive_counts', 'student_responsive_expression', 'student_responsive_effects', page_stems, 'student_responsive_compact_heatmap')
index <- fread('figures/student_figure_index.csv')
index <- index[!stem %in% stems]
index <- rbind(index, data.table(stem = stems,
  description = c('Exploratory response intersections and TREAT sensitivity',
                  'All CI-contained response candidates across 18 libraries',
                  'All 91 candidates ordered by UTY effect: four contrasts',
                  paste('UTY-ranked candidate effect intervals, page', 1:4),
                  'Compact 91-gene group-mean heatmap, UTY effect order, centered on LXHA'), formats = 'PDF;SVG;PNG'))
fwrite(index, 'figures/student_figure_index.csv')
cat('Responsive-gene exports complete.\n')
