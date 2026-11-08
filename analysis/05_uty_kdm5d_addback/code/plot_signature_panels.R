# Shared plots: called by run_analysis.R or restyle_signature_figure.R.
score_df$final_group <- factor(score_df$final_group, group_levels)
p_score <- ggplot(score_df, aes(final_group, signed_percentile_score, colour = final_group, shape = factor(collection_occasion))) +
  geom_point(size = 2.7, position = position_jitter(width = 0.08, height = 0, seed = 7)) +
  stat_summary(fun = mean, geom = 'crossbar', aes(group = final_group), width = 0.48, colour = '#222222') +
  scale_colour_manual(values = cols, guide = 'none') + scale_shape_manual(values = c(16, 17, 15), name = 'Collection\noccasion') +
  scale_x_discrete(labels = labs_group) + labs(x = NULL, y = 'Frozen August signature score') +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = 'right')

forest_keep <- c('UTY_v_LXHA', 'KDM5D_v_LXHA', 'DUAL_v_LXHA', 'DUAL_v_UTY', 'DUAL_v_KDM5D')
forest <- score_contrasts[match(forest_keep, score_contrasts$contrast), ]
forest$label <- factor(c('UTY vs vector', 'KDM5D vs vector', 'UTY + KDM5D vs vector',
                         'UTY + KDM5D vs UTY', 'UTY + KDM5D vs KDM5D'), levels = rev(c('UTY vs vector', 'KDM5D vs vector', 'UTY + KDM5D vs vector', 'UTY + KDM5D vs UTY', 'UTY + KDM5D vs KDM5D')))
p_forest <- ggplot(forest, aes(difference, label, xmin = CI_L, xmax = CI_R)) +
  geom_vline(xintercept = 0, linetype = 2, colour = 'grey60') +
  geom_errorbar(aes(ymin = label, ymax = label), orientation = 'y', width = 0.18) + geom_point(size = 2.4) +
  geom_text(aes(label = ifelse(PValue < 0.001, sprintf('P = %.1e', PValue), paste0('P = ', formatC(PValue, format = 'f', digits = 2)))), x = Inf,
            hjust = 1, nudge_x = -0.002, size = 2.7) +
  labs(x = 'Difference in signature score', y = NULL) + scale_x_continuous(expand = expansion(mult = c(0.1, 0.5)))

heat <- as.data.frame(gm - gm[, 'CTRL'])
heat$gene_id <- rownames(heat)
heat <- merge(heat, sig[, c('gene_id_stable', 'AUG_difference', 'direction')], by.x = 'gene_id', by.y = 'gene_id_stable')
heat <- heat[order(-heat$direction, -abs(heat$AUG_difference)), ]
heat$order <- seq_len(nrow(heat))
heat_long <- melt(as.data.table(heat), id.vars = c('gene_id', 'AUG_difference', 'direction', 'order'),
                  measure.vars = c('LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'),
                  variable.name = 'group', value.name = 'logFC_vs_CTRL')
heat_long$group <- factor(heat_long$group, c('LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'))
p_heat <- ggplot(heat_long, aes(group, factor(order), fill = pmax(pmin(logFC_vs_CTRL, 4), -4))) +
  geom_tile() + scale_fill_gradient2(low = '#2C6DA4', mid = 'white', high = '#D9822B', midpoint = 0,
    limits = c(-4, 4), name = 'log2 expression\ndifference vs sgCTRL\n(clipped to [-4, 4])') +
  scale_x_discrete(labels = labs_group) + labs(x = NULL, y = 'Frozen August signature genes') +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Same retained fry results, explicitly labelled directions and all adjusted P values.
fr_plot <- fr
fr_plot$comparison <- factor(fr_plot$contrast, rev(fry_names),
  labels = rev(c('UTY vs LXHA', 'KDM5D vs LXHA', 'UTY + KDM5D\nvs LXHA',
                 'UTY + KDM5D\nvs UTY', 'UTY + KDM5D\nvs KDM5D', 'Y-depleted vs LXHA')))
arm_n <- tapply(fr_plot$NGenes, fr_plot$signature_arm, unique)
fr_plot$arm <- factor(fr_plot$signature_arm, c('August_up', 'August_down'),
  labels = c(sprintf('August-up genes (n = %d)\nReversal: DECREASE', arm_n[['August_up']]),
             sprintf('August-down genes (n = %d)\nReversal: INCREASE', arm_n[['August_down']])))
fr_plot$direction_label <- ifelse(fr_plot$Direction == 'Down', 'Decrease', 'Increase')
fr_plot$cell_label <- paste0(fr_plot$direction_label, '\nFDR = ',
  vapply(fr_plot$FDR_12, function(x) if (x < 0.001) formatC(x, format = 'e', digits = 1)
         else formatC(x, format = 'f', digits = 3), character(1)))
fr_plot$label_face <- ifelse(fr_plot$FDR_12 < 0.05, 'bold', 'plain')
p_fry <- ggplot(fr_plot, aes(arm, comparison)) +
  geom_tile(aes(fill = Direction), colour = 'white', linewidth = 1, width = 0.97, height = 0.95) +
  geom_text(aes(label = cell_label, fontface = label_face), size = 2.8, lineheight = 1.15) +
  scale_fill_manual(values = c(Down = '#D6E5F1', Up = '#F7DEC5'), guide = 'none') +
  scale_x_discrete(position = 'top', expand = expansion(add = 0.02)) +
  scale_y_discrete(expand = expansion(add = 0.05)) +
  labs(x = NULL, y = NULL, title = 'Do both signature halves shift toward reversal?',
       subtitle = 'fry gene-set tests; each comparison is first group minus second',
       caption = 'Bold: FDR < 0.05 (BH adjustment across all 12 tests).\nDirection and significance are distinct; color indicates direction only.') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, face = 'bold', lineheight = 1.2),
        axis.text.y = element_text(size = 8), plot.title = element_text(size = 9, face = 'bold'),
        plot.subtitle = element_text(size = 8), plot.caption = element_text(size = 7.5, hjust = 0))
fwrite(fr_plot, 'outputs/signature/fry_panel_display.csv')

fig2 <- ((p_score + labs(tag = 'a')) | (p_forest + labs(tag = 'b'))) /
        ((p_heat + labs(tag = 'c')) | (p_fry + labs(tag = 'd'))) +
        plot_layout(heights = c(1, 1.25), widths = c(1.05, 1))
save_multi(fig2, 'ch03_main_ht1080_addback_baseline_reversal', 12, 9)
save_multi(p_score, 'student_panels/panel_f_frozen_signature_scores', 8.5, 5.8)
save_multi(p_forest, 'student_panels/panel_g_signature_score_contrasts', 8.5, 5.8)
save_multi(p_heat, 'student_panels/panel_h_signature_gene_heatmap', 7.2, 8.5)
save_multi(p_fry, 'student_panels/panel_i_directional_fry_tests', 7.2, 5.4)


# Windows Cairo substitutes Arial; SVG-to-PDF export embeds the bundled fonts.
status <- system2(Sys.which('python'), 'code/export_signature_pdfs.py')
if (status != 0) stop('Signature PDF font export failed')
