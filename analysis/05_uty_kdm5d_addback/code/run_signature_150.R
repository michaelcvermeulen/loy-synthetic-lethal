# Baseline-reversal figure on the balanced 150-gene LOY signature (75 higher / 75 lower in LOY; frozen August rules
# plus August-only limma-voom FDR < 0.20). Same scoring, blocked score model, fry tests and gene-level trajectories as
# run_analysis.R for the original 100-gene signature. The original outputs and figure are left untouched: results go to
# outputs/signature_150/ and figures/*_sig150.*. Run from the add-back analysis directory or from THESIS_ADDBACK_2026.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(ggplot2); library(patchwork); library(digest)})
roots <- c(getwd(), file.path(getwd(), 'THESIS_ADDBACK_2026'))
root <- roots[file.exists(file.path(roots, 'outputs/dge_primary.rds'))][1]
setwd(root)
STYLE <- 'A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/skills/thesis-figure-style'
source(file.path(STYLE, 'scripts/thesis_style.R'))
register_thesis_fonts(STYLE)
theme_set(theme_thesis(base_size = 9))
update_geom_defaults('text', list(family = 'Open Sans', colour = '#222222'))
dir.create('outputs/signature_150', showWarnings = FALSE)
SIG150 <- 'A:/Dropbox/LOY_paralogs/projects/thesis/Ben_Thesis/Aim2_addback_signature_recovery/outputs/signature_variants/signature_genes_by_variant.csv'
BACKGROUND <- 'A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/tables/August_signature_discovery_background.csv'
src <- c(signature_variants = SIG150, background = BACKGROUND, fit = 'outputs/fit_primary.rds', dge = 'outputs/dge_primary.rds',
         samples = 'outputs/sample_manifest.csv', counts = '../SEPT_2026_REANALYSIS/outputs/03_counts/counts_full.rds')
stopifnot(all(file.exists(src)))
fwrite(data.table(role = names(src), path = normalizePath(src, winslash = '/'), sha256 = vapply(src, digest, character(1), algo = 'sha256', file = TRUE)),
       'outputs/signature_150/input_manifest.csv')

group_levels <- c('CTRL', 'LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
cols <- c(CTRL = '#3C5488', LOY = '#E64B35', LOY_LXHA = '#7E6148', LOY_UTY_AB = '#00A087', LOY_KDM5D_AB = '#F39B7F', LOY_KDM5D_UTY_AB = '#4DBBD5')
labs_group <- c(CTRL = 'sgCTRL', LOY = 'Y-depleted', LOY_LXHA = 'Empty vector', LOY_UTY_AB = 'UTY', LOY_KDM5D_AB = 'KDM5D', LOY_KDM5D_UTY_AB = 'UTY + KDM5D')

# ---- signature and scores ----------------------------------------------------------------------------------------
sig <- fread(SIG150, data.table = FALSE); sig <- sig[sig$variant == 'augFDR20_150', ]
sig$gene_id_stable <- sig$stable_gene_id
stopifnot(nrow(sig) == 150, sum(sig$direction > 0) == 75, sum(sig$direction < 0) == 75)
bg <- fread(BACKGROUND, data.table = FALSE)
counts <- readRDS(src[['counts']]); samples <- fread(src[['samples']], data.table = FALSE)
stopifnot(identical(colnames(counts), samples$sample_id))
bg_ids <- intersect(bg$stable_gene_id[bg$representative], rownames(counts))
sig$represented <- sig$gene_id_stable %in% bg_ids; stopifnot(all(sig$represented))
rank_mat <- apply(counts[bg_ids, , drop = FALSE], 2, function(x) (rank(x, ties.method = 'average') - 1) / (length(x) - 1)); rownames(rank_mat) <- bg_ids
up <- sig$gene_id_stable[sig$direction > 0]; down <- sig$gene_id_stable[sig$direction < 0]
score <- colMeans(rank_mat[up, , drop = FALSE]) - colMeans(rank_mat[down, , drop = FALSE])
score_df <- data.frame(sample_id = names(score), final_group = samples$final_group, collection_occasion = samples$replicate_block, signed_percentile_score = score)
fwrite(sig, 'outputs/signature_150/signature_mapping.csv'); fwrite(score_df, 'outputs/signature_150/signature_scores.csv')

fit <- readRDS(src[['fit']]); design <- fit$v$design
fit_score <- eBayes(contrasts.fit(lmFit(matrix(score, nrow = 1), design), fit$cm), robust = TRUE)
score_contrasts <- do.call(rbind, lapply(colnames(fit$cm), function(k) {
  z <- topTable(fit_score, coef = k, number = Inf, confint = TRUE)
  data.frame(contrast = k, difference = z$logFC, CI_L = z$CI.L, CI_R = z$CI.R, t = z$t, PValue = z$P.Value) }))
score_contrasts$FDR_15 <- p.adjust(score_contrasts$PValue, 'BH')
fwrite(score_contrasts, 'outputs/signature_150/score_contrasts.csv')

fry_names <- c('UTY_v_LXHA', 'KDM5D_v_LXHA', 'DUAL_v_LXHA', 'DUAL_v_UTY', 'DUAL_v_KDM5D', 'LOY_v_LXHA')
idx <- list(August_up = which(rownames(fit$v) %in% up), August_down = which(rownames(fit$v) %in% down))
fr <- do.call(rbind, lapply(fry_names, function(k) { z <- fry(fit$v, index = idx, design = design, contrast = fit$cm[, k]); z$signature_arm <- rownames(z); z$contrast <- k; z }))
fr$FDR_12 <- p.adjust(fr$PValue, 'BH')
fr$expected_reversal <- ifelse(fr$signature_arm == 'August_up', 'Down', 'Up'); fr$direction_matches_reversal <- fr$Direction == fr$expected_reversal
fwrite(fr, 'outputs/signature_150/fry_directional.csv')

dge <- readRDS(src[['dge']]); lcpm <- cpm(dge, log = TRUE, prior.count = 2)
sig_ids <- sig$gene_id_stable[sig$gene_id_stable %in% rownames(lcpm)]
gm <- sapply(group_levels, function(g) rowMeans(lcpm[sig_ids, samples$final_group == g, drop = FALSE]))
movement <- do.call(rbind, lapply(c('LOY', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'), function(arm) {
  vector_effect <- gm[, 'LOY_LXHA'] - gm[, 'CTRL']; arm_effect <- gm[, arm] - gm[, 'CTRL']
  crossed <- sign(arm_effect) != sign(vector_effect) & arm_effect != 0 & vector_effect != 0; closer <- abs(arm_effect) < abs(vector_effect)
  data.frame(gene_id = sig_ids, arm = arm, vector_vs_CTRL = vector_effect, arm_vs_CTRL = arm_effect, closer_to_CTRL = closer, crossed_CTRL = crossed,
             category = ifelse(crossed & closer, 'crossed_control_closer', ifelse(crossed & !closer, 'crossed_control_farther', ifelse(!crossed & closer, 'closer_same_side', 'farther_same_side')))) }))
fwrite(movement, 'outputs/signature_150/distance_to_control_by_gene.csv')
ms <- aggregate(gene_id ~ arm + category, movement, length); names(ms)[3] <- 'n_genes'
ms$fraction <- ave(ms$n_genes, ms$arm, FUN = function(x) x / sum(x)); fwrite(ms, 'outputs/signature_150/distance_to_control_summary.csv')

# ---- panels (as plot_signature_panels.R, with LOY-signature labels and 150-gene arms) --------------------------------
save_multi <- function(p, stem, width, height) {
  ggsave(file.path('figures', paste0(stem, '.png')), p, width = width, height = height, dpi = 300, bg = 'white')
  ggsave(file.path('figures', paste0(stem, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path('figures', paste0(stem, '.svg')), p, width = width, height = height, device = svglite::svglite) }
score_df$final_group <- factor(score_df$final_group, group_levels)
p_score <- ggplot(score_df, aes(final_group, signed_percentile_score, colour = final_group, shape = factor(collection_occasion))) +
  geom_point(size = 2.7, position = position_jitter(width = 0.08, height = 0, seed = 7)) +
  stat_summary(fun = mean, geom = 'crossbar', aes(group = final_group), width = 0.48, colour = '#222222') +
  scale_colour_manual(values = cols, guide = 'none') + scale_shape_manual(values = c(16, 17, 15), name = 'Collection\noccasion') +
  scale_x_discrete(labels = labs_group) + labs(x = NULL, y = 'LOY signature score (150 genes)') +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = 'right')
forest_keep <- c('UTY_v_LXHA', 'KDM5D_v_LXHA', 'DUAL_v_LXHA', 'DUAL_v_UTY', 'DUAL_v_KDM5D')
forest <- score_contrasts[match(forest_keep, score_contrasts$contrast), ]
fl <- c('UTY vs vector', 'KDM5D vs vector', 'UTY + KDM5D vs vector', 'UTY + KDM5D vs UTY', 'UTY + KDM5D vs KDM5D')
forest$label <- factor(fl, levels = rev(fl))
p_forest <- ggplot(forest, aes(difference, label, xmin = CI_L, xmax = CI_R)) +
  geom_vline(xintercept = 0, linetype = 2, colour = 'grey60') +
  geom_errorbar(aes(ymin = label, ymax = label), orientation = 'y', width = 0.18) + geom_point(size = 2.4) +
  geom_text(aes(label = ifelse(PValue < 0.001, sprintf('P = %.1e', PValue), paste0('P = ', formatC(PValue, format = 'f', digits = 2)))), x = Inf, hjust = 1, nudge_x = -0.002, size = 2.7) +
  labs(x = 'Difference in LOY signature score', y = NULL) + scale_x_continuous(expand = expansion(mult = c(0.1, 0.5)))
heat <- as.data.frame(gm - gm[, 'CTRL']); heat$gene_id <- rownames(heat)
heat <- merge(heat, sig[, c('gene_id_stable', 'AUG_difference', 'direction')], by.x = 'gene_id', by.y = 'gene_id_stable')
heat <- heat[order(-heat$direction, -abs(heat$AUG_difference)), ]; heat$order <- seq_len(nrow(heat))
heat_long <- melt(as.data.table(heat), id.vars = c('gene_id', 'AUG_difference', 'direction', 'order'),
                  measure.vars = c('LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'), variable.name = 'group', value.name = 'logFC_vs_CTRL')
heat_long$group <- factor(heat_long$group, c('LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'))
p_heat <- ggplot(heat_long, aes(group, factor(order), fill = pmax(pmin(logFC_vs_CTRL, 4), -4))) +
  geom_tile() + scale_fill_gradient2(low = '#2C6DA4', mid = 'white', high = '#D9822B', midpoint = 0, limits = c(-4, 4), name = 'log2 expression\ndifference vs sgCTRL\n(clipped to [-4, 4])') +
  scale_x_discrete(labels = labs_group) + labs(x = NULL, y = sprintf('LOY signature genes (%d expressed of 150)', length(sig_ids))) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank())
fr_plot <- fr
fr_plot$comparison <- factor(fr_plot$contrast, rev(fry_names), labels = rev(c('UTY vs LXHA', 'KDM5D vs LXHA', 'UTY + KDM5D\nvs LXHA', 'UTY + KDM5D\nvs UTY', 'UTY + KDM5D\nvs KDM5D', 'Y-depleted vs LXHA')))
arm_n <- tapply(fr_plot$NGenes, fr_plot$signature_arm, unique)
fr_plot$arm <- factor(fr_plot$signature_arm, c('August_up', 'August_down'),
  labels = c(sprintf('Higher-in-LOY genes (n = %d)\nReversal: DECREASE', arm_n[['August_up']]), sprintf('Lower-in-LOY genes (n = %d)\nReversal: INCREASE', arm_n[['August_down']])))
fr_plot$direction_label <- ifelse(fr_plot$Direction == 'Down', 'Decrease', 'Increase')
fr_plot$cell_label <- paste0(fr_plot$direction_label, '\nFDR = ', vapply(fr_plot$FDR_12, function(x) if (x < 0.001) formatC(x, format = 'e', digits = 1) else formatC(x, format = 'f', digits = 3), character(1)))
fr_plot$label_face <- ifelse(fr_plot$FDR_12 < 0.05, 'bold', 'plain')
p_fry <- ggplot(fr_plot, aes(arm, comparison)) +
  geom_tile(aes(fill = Direction), colour = 'white', linewidth = 1, width = 0.97, height = 0.95) +
  geom_text(aes(label = cell_label, fontface = label_face), size = 2.8, lineheight = 1.15) +
  scale_fill_manual(values = c(Down = '#D6E5F1', Up = '#F7DEC5'), guide = 'none') +
  scale_x_discrete(position = 'top', expand = expansion(add = 0.02)) + scale_y_discrete(expand = expansion(add = 0.05)) +
  labs(x = NULL, y = NULL, title = 'Do both signature halves shift toward reversal?', subtitle = 'fry gene-set tests; each comparison is first group minus second',
       caption = 'Bold: FDR < 0.05 (BH adjustment across all 12 tests).\nDirection and significance are distinct; color indicates direction only.') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 8, face = 'bold', lineheight = 1.2),
        axis.text.y = element_text(size = 8), plot.title = element_text(size = 9, face = 'bold'), plot.subtitle = element_text(size = 8), plot.caption = element_text(size = 7.5, hjust = 0))
fwrite(fr_plot, 'outputs/signature_150/fry_panel_display.csv')
fig2 <- ((p_score + labs(tag = 'a')) | (p_forest + labs(tag = 'b'))) / ((p_heat + labs(tag = 'c')) | (p_fry + labs(tag = 'd'))) + plot_layout(heights = c(1, 1.25), widths = c(1.05, 1))
stems <- c('ch03_main_ht1080_addback_baseline_reversal_sig150', 'student_panels/panel_f_signature_scores_sig150', 'student_panels/panel_g_signature_score_contrasts_sig150',
           'student_panels/panel_h_signature_gene_heatmap_sig150', 'student_panels/panel_i_directional_fry_tests_sig150')
save_multi(fig2, stems[1], 12, 9); save_multi(p_score, stems[2], 8.5, 5.8); save_multi(p_forest, stems[3], 8.5, 5.8); save_multi(p_heat, stems[4], 7.2, 8.5); save_multi(p_fry, stems[5], 7.2, 5.4)
status <- system2(Sys.which('python'), c('code/export_signature_pdfs.py', stems)); if (status != 0) stop('Signature PDF font export failed')
writeLines(c(sprintf('PASS: 150-gene signature (75 up / 75 down) scored on the 23,903-gene August background; %d genes pass the May expression filter', length(sig_ids)),
             sprintf('PASS: UTY vs vector score difference %.3f (P = %.1e); combined %.3f (P = %.1e); KDM5D %.3f (P = %.2f)',
                     forest$difference[1], forest$PValue[1], forest$difference[3], forest$PValue[3], forest$difference[2], forest$PValue[2]),
             'Original 100-gene outputs and figure untouched; this run writes outputs/signature_150/ and figures/*_sig150.*'), 'outputs/signature_150/VERIFICATION.txt')
writeLines(capture.output(sessionInfo()), 'outputs/signature_150/sessionInfo.txt')
print(forest[, c('contrast', 'difference', 'CI_L', 'CI_R', 'PValue')]); print(fr[, c('signature_arm', 'contrast', 'NGenes', 'Direction', 'FDR_12')]); print(ms)
