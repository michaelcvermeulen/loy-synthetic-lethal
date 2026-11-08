# Presentation only: reuse saved CAMERA statistics, no model fitting.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- 'A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026'
OUT <- 'A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3'
STYLE <- 'A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/skills/thesis-figure-style'
source(file.path(STYLE, 'scripts/thesis_style.R'))
register_thesis_fonts(STYLE)
theme_set(theme_thesis(base_size = 9))
update_geom_defaults('text', list(family = 'Open Sans', colour = '#222222'))
sens_counts <- fread(file.path(ROOT,'outputs/camera_sensitivity/camera_counts_by_correlation.csv'))
settings <- fread(file.path(ROOT,'outputs/camera_sensitivity/camera_all_settings.csv'))
spec <- fread(file.path(ROOT,'outputs/camera_sensitivity/specificity_vs_naive_LOY.csv'))
collections <- c('Hallmark','Reactome','GO_BP')
contrast_labels <- c(UTY_v_LXHA='UTY',KDM5D_v_LXHA='KDM5D',DUAL_v_LXHA='Dual',DUAL_v_UTY='Dual vs UTY',DUAL_v_KDM5D='Dual vs KDM5D',LOY_v_LXHA='Naive')

lab_setting <- c('fixed_0.01' = 'Default CAMERA (0.01)')
setting_cols <- c('fixed_0.01' = '#3C5488')
lab_short <- c(UTY_v_LXHA = 'UTY vs vector', KDM5D_v_LXHA = 'KDM5D vs vector', DUAL_v_LXHA = 'Dual vs vector',
               DUAL_v_UTY = 'Dual vs UTY', DUAL_v_KDM5D = 'Dual vs KDM5D', LOY_v_LXHA = 'Naive Y-dep. vs vector')

pa_dat <- copy(sens_counts[setting == 'fixed_0.01'])
pa_dat[, contrast := factor(contrast, names(contrast_labels))]
pa_dat[, collection := factor(collection, collections)]
p_counts <- ggplot(pa_dat, aes(contrast, sets_FDR_lt_0.05, fill = setting)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sets_FDR_lt_0.05), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.2) +
  facet_wrap(~ collection, scales = 'free_y', labeller = as_labeller(c(Hallmark = 'Hallmark', Reactome = 'Reactome', GO_BP = 'GO BP'))) +
  scale_x_discrete(labels = lab_short) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = setting_cols, labels = lab_setting, name = 'Inter-gene correlation') +
  labs(x = NULL, y = 'Gene sets with FDR < 0.05') +
  theme(legend.position = 'top', axis.text.x = element_text(angle = 40, hjust = 1, lineheight = 0.9))

hall <- settings[setting == 'fixed_0.01' & collection == 'Hallmark']
hw <- dcast(hall, gene_set ~ contrast, value.var = c('signed_z', 'FDR'))
hw[, label := gsub('_', ' ', sub('^HALLMARK_', '', gene_set))]
hw[, class := fifelse(FDR_UTY_v_LXHA < 0.05 & FDR_LOY_v_LXHA < 0.05 &
                        sign(signed_z_UTY_v_LXHA) == sign(signed_z_LOY_v_LXHA), 'Both FDR < 0.05, same direction',
                fifelse(FDR_UTY_v_LXHA < 0.05, 'UTY only', fifelse(FDR_LOY_v_LXHA < 0.05, 'Naive only', 'Neither')))]
rho_h <- spec[setting == 'fixed_0.01' & collection == 'Hallmark' & contrast == 'UTY_v_LXHA', spearman_vs_naive]
lim <- max(abs(c(hw$signed_z_UTY_v_LXHA, hw$signed_z_LOY_v_LXHA))) * 1.05
to_label <- hw[order(-abs(signed_z_UTY_v_LXHA))][1:8]
# Labels on the right half are placed to the left of their point to stay in panel.
to_label[, hj := ifelse(signed_z_LOY_v_LXHA > 0, 1.08, -0.08)]
p_scatter <- ggplot(hw, aes(signed_z_LOY_v_LXHA, signed_z_UTY_v_LXHA)) +
  geom_hline(yintercept = 0, colour = 'grey80') + geom_vline(xintercept = 0, colour = 'grey80') +
  geom_abline(slope = 1, intercept = 0, colour = 'grey70', linetype = 2) +
  geom_point(aes(colour = class), size = 1.8) +
  geom_text(data = to_label, aes(label = label, hjust = hj), size = 2.1, vjust = 0.5, check_overlap = TRUE) +
  annotate('text', x = -lim, y = lim, hjust = 0, vjust = 1, size = 2.8,
           label = sprintf('Spearman rho = %.2f\n50 Hallmark sets', rho_h)) +
  scale_colour_manual(values = c('Both FDR < 0.05, same direction' = '#E64B35', 'UTY only' = '#00A087',
                                 'Naive only' = '#7E6148', 'Neither' = 'grey70'), name = NULL) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = 'CAMERA signed z: naive Y-depleted vs vector', y = 'CAMERA signed z: UTY vs vector') +
  theme(legend.position = 'bottom') + guides(colour = guide_legend(nrow = 2))

sp_dat <- spec[setting == 'fixed_0.01']
sp_dat[, contrast := factor(contrast, c('UTY_v_LXHA', 'DUAL_v_LXHA', 'KDM5D_v_LXHA'))]
sp_dat[, collection := factor(collection, collections)]
p_shared <- ggplot(sp_dat, aes(collection, spearman_vs_naive, fill = contrast)) +
  geom_hline(yintercept = 0, colour = 'grey70') +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf('%d/%d', shared_same_direction, sets_FDR_lt_0.05)),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 2.2) +
  scale_fill_manual(values = c(UTY_v_LXHA = '#00A087', DUAL_v_LXHA = '#4DBBD5', KDM5D_v_LXHA = '#F39B7F'),
                    labels = c(UTY_v_LXHA = 'UTY', DUAL_v_LXHA = 'UTY + KDM5D', KDM5D_v_LXHA = 'KDM5D'), name = 'vs vector') +
  scale_y_continuous(limits = c(min(0, min(sp_dat$spearman_vs_naive)), 1.08), breaks = seq(0, 1, 0.25)) +
  scale_x_discrete(labels = c(Hallmark = 'Hallmark', Reactome = 'Reactome', GO_BP = 'GO BP')) +
  labs(x = NULL, y = 'Spearman rho with naive\nY-depleted vs vector') +
  theme(legend.position = 'top')

fig <- (p_counts + labs(tag = 'a')) / ((p_scatter + labs(tag = 'b')) | (p_shared + labs(tag = 'c'))) +
  plot_layout(heights = c(1, 1.35))
save_multi <- function(p, stem, width, height) {
  ggsave(file.path(OUT, paste0(stem, '.png')), p, width = width, height = height, dpi = 300, bg = 'white')
  ggsave(file.path(OUT, paste0(stem, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(stem, '.svg')), p, width = width, height = height, device = svglite::svglite)
}
save_multi(fig, 'ch03_supp_ht1080_addback_pathway_sensitivity', 11, 10)
