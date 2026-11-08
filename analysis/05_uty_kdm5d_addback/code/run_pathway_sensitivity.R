options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

# Sensitivity and specificity of the CAMERA pathway results produced by
# run_camera_report.R. Two questions:
#   1. Do the pathway calls depend on the fixed inter-gene correlation of 0.01?
#   2. Are add-back pathway shifts specific to the construct, or shared with the
#      untransduced Y-depleted population compared with the same empty vector?
# Run run_analysis.R and run_camera_report.R first.

script_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)[1]
script <- normalizePath(sub('^--file=', '', script_arg), winslash = '/', mustWork = TRUE)
ROOT <- dirname(dirname(script))
STYLE <- 'A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/skills/thesis-figure-style'
setwd(ROOT)
set.seed(42)
dir.create('outputs/camera_sensitivity', recursive = TRUE, showWarnings = FALSE)

source(file.path(STYLE, 'scripts/thesis_style.R'))
register_thesis_fonts(STYLE)
theme_set(theme_thesis(base_size = 9))
update_geom_defaults('text', list(family = 'Open Sans', colour = '#222222'))

fb <- readRDS('outputs/fit_primary.rds')
v <- fb$v
cm <- fb$cm
map <- fread('outputs/camera/camera_gene_universe.csv')
mem <- fread('outputs/camera/gene_set_membership_tested.csv')
camera_primary <- fread('outputs/camera/camera_all.csv')
v_camera <- v[map$original_row, ]
stopifnot(nrow(v_camera) >= max(mem$camera_row))

contrast_labels <- c(
  UTY_v_LXHA = 'UTY vs vector',
  KDM5D_v_LXHA = 'KDM5D vs vector',
  DUAL_v_LXHA = 'UTY + KDM5D vs vector',
  DUAL_v_UTY = 'UTY + KDM5D vs UTY',
  DUAL_v_KDM5D = 'UTY + KDM5D vs KDM5D',
  LOY_v_LXHA = 'Naive Y-depleted vs vector'
)
collections <- c('Hallmark', 'Reactome', 'GO_BP')
index_by_collection <- lapply(setNames(collections, collections), function(cl) {
  m <- mem[collection == cl]
  split(m$camera_row, m$gene_set)
})

# ---------------- 1. Inter-gene correlation sensitivity ----------------
run_camera <- function(cor_value) {
  rbindlist(lapply(collections, function(cl) {
    rbindlist(lapply(names(contrast_labels), function(cn) {
      message('CAMERA ', cl, ' / ', cn, ' / inter.gene.cor = ', ifelse(is.na(cor_value), 'estimated', cor_value))
      a <- camera(v_camera, index_by_collection[[cl]], design = v$design, contrast = cm[, cn],
                  inter.gene.cor = cor_value, allow.neg.cor = FALSE, sort = FALSE)
      data.table(collection = cl, contrast = cn, gene_set = rownames(a), n_genes = a$NGenes,
                 direction = a$Direction, p_value = a$PValue,
                 estimated_correlation = if ('Correlation' %in% names(a)) a$Correlation else NA_real_)
    }))
  }))
}

est <- run_camera(NA)
cor_summary <- est[, .(median_estimated_correlation = median(estimated_correlation),
                       q25 = quantile(estimated_correlation, 0.25),
                       q75 = quantile(estimated_correlation, 0.75)), by = .(collection, contrast)]
# The fixed intermediate setting is the median per-set estimate across all
# collections and contrasts, rounded to two decimals.
cor_mid <- round(median(est$estimated_correlation), 2)
message('Median estimated inter-gene correlation: ', cor_mid)
mid <- run_camera(cor_mid)

prim <- camera_primary[, .(collection, contrast, gene_set, n_genes, direction, p_value)]
settings <- rbindlist(list(
  prim[, setting := 'fixed_0.01'],
  mid[, .(collection, contrast, gene_set, n_genes, direction, p_value)][, setting := sprintf('fixed_%.2f', cor_mid)],
  est[, .(collection, contrast, gene_set, n_genes, direction, p_value)][, setting := 'estimated_per_set']
))
settings[, FDR := p.adjust(p_value, 'BH'), by = .(setting, collection, contrast)]
settings[, signed_z := ifelse(direction == 'Up', 1, -1) * qnorm(p_value / 2, lower.tail = FALSE)]
setting_levels <- c('fixed_0.01', sprintf('fixed_%.2f', cor_mid), 'estimated_per_set')
settings[, setting := factor(setting, setting_levels)]

sens_counts <- settings[, .(sets_tested = .N, sets_FDR_lt_0.05 = sum(FDR < 0.05),
                            up = sum(FDR < 0.05 & direction == 'Up'),
                            down = sum(FDR < 0.05 & direction == 'Down')),
                        by = .(setting, collection, contrast)]
setorder(sens_counts, collection, contrast, setting)

# ---------------- 2. Specificity versus naive Y-depleted ----------------
# The untransduced Y-depleted population lacks both constructs. A pathway shift
# that it shares with an add-back population, relative to the same empty-vector
# reference, cannot be attributed to the construct.
wide <- dcast(settings, setting + collection + gene_set ~ contrast, value.var = c('signed_z', 'FDR', 'direction'))
spec <- rbindlist(lapply(setdiff(names(contrast_labels), 'LOY_v_LXHA'), function(cn) {
  if (!grepl('_v_LXHA$', cn)) return(NULL)
  w <- copy(wide)
  w[, `:=`(fdr_c = get(paste0('FDR_', cn)), dir_c = get(paste0('direction_', cn)),
           z_c = get(paste0('signed_z_', cn)),
           fdr_l = FDR_LOY_v_LXHA, dir_l = direction_LOY_v_LXHA, z_l = signed_z_LOY_v_LXHA)]
  w[, .(contrast = cn,
        spearman_vs_naive = cor(z_c, z_l, method = 'spearman'),
        sets_FDR_lt_0.05 = sum(fdr_c < 0.05),
        shared_same_direction = sum(fdr_c < 0.05 & fdr_l < 0.05 & dir_c == dir_l),
        naive_opposite_direction = sum(fdr_c < 0.05 & fdr_l < 0.05 & dir_c != dir_l),
        not_detected_in_naive = sum(fdr_c < 0.05 & fdr_l >= 0.05)),
    by = .(setting, collection)]
}))
spec[, fraction_shared := ifelse(sets_FDR_lt_0.05 > 0, shared_same_direction / sets_FDR_lt_0.05, NA_real_)]
setorder(spec, setting, collection, contrast)

# Candidate construct-associated sets: FDR < 0.05 for UTY vs vector at both
# fixed settings, and not significant in the same direction for naive
# Y-depleted vs vector at the primary setting. Listed for inspection only.
p01 <- settings[setting == 'fixed_0.01']
pm <- settings[setting == sprintf('fixed_%.2f', cor_mid)]
cand <- merge(
  p01[contrast == 'UTY_v_LXHA', .(collection, gene_set, n_genes, direction, FDR_0.01 = FDR, z_uty = signed_z)],
  pm[contrast == 'UTY_v_LXHA', .(collection, gene_set, FDR_mid = FDR)], by = c('collection', 'gene_set'))
cand <- merge(cand, p01[contrast == 'LOY_v_LXHA', .(collection, gene_set, direction_naive = direction,
                                                   FDR_naive = FDR, z_naive = signed_z)],
              by = c('collection', 'gene_set'))
cand <- merge(cand, p01[contrast == 'DUAL_v_KDM5D', .(collection, gene_set, direction_dual_v_kdm5d = direction,
                                                     FDR_dual_v_kdm5d = FDR)],
              by = c('collection', 'gene_set'))
cand <- cand[FDR_0.01 < 0.05 & FDR_mid < 0.05 & !(FDR_naive < 0.05 & direction_naive == direction)]
cand[, replicated_in_dual_v_kdm5d := FDR_dual_v_kdm5d < 0.05 & direction_dual_v_kdm5d == direction]
setorder(cand, collection, FDR_mid)
setnames(cand, 'FDR_mid', sprintf('FDR_%.2f', cor_mid))

fwrite(cor_summary, 'outputs/camera_sensitivity/estimated_intergene_correlation.csv')
fwrite(sens_counts, 'outputs/camera_sensitivity/camera_counts_by_correlation.csv')
fwrite(settings, 'outputs/camera_sensitivity/camera_all_settings.csv')
fwrite(spec, 'outputs/camera_sensitivity/specificity_vs_naive_LOY.csv')
fwrite(cand, 'outputs/camera_sensitivity/uty_sets_not_shared_with_naive_LOY.csv')

# ---------------- Figure ----------------
lab_setting <- c('fixed 0.01 (primary)', sprintf('fixed %.2f (median estimate)', cor_mid), 'estimated per set')
names(lab_setting) <- setting_levels
setting_cols <- setNames(c('#3C5488', '#00A087', '#B09C85'), setting_levels)
lab_short <- c(UTY_v_LXHA = 'UTY vs vector', KDM5D_v_LXHA = 'KDM5D vs vector', DUAL_v_LXHA = 'Dual vs vector',
               DUAL_v_UTY = 'Dual vs UTY', DUAL_v_KDM5D = 'Dual vs KDM5D', LOY_v_LXHA = 'Naive Y-dep. vs vector')

pa_dat <- copy(sens_counts)
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
  ggsave(file.path('figures', paste0(stem, '.png')), p, width = width, height = height, dpi = 300, bg = 'white')
  ggsave(file.path('figures', paste0(stem, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path('figures', paste0(stem, '.svg')), p, width = width, height = height, device = svglite::svglite)
}
save_multi(fig, 'ch03_supp_ht1080_addback_pathway_sensitivity', 11, 10)
# Windows Cairo substitutes Arial; reprint the PDF from SVG with embedded Open Sans.
if (system2(Sys.which('python'), c('code/export_signature_pdfs.py', 'ch03_supp_ht1080_addback_pathway_sensitivity')) != 0)
  stop('Open Sans PDF export failed')

capture.output(sessionInfo(), file = 'outputs/camera_sensitivity/sessionInfo.txt')
verification <- c(
  sprintf('Median estimated inter-gene correlation (all sets, all contrasts): %.3f -> fixed setting %.2f',
          median(est$estimated_correlation), cor_mid),
  sprintf('Rows per setting: %s', paste(settings[, .N, by = setting]$N, collapse = ' / ')),
  sprintf('Primary setting reproduces camera_all.csv FDR_collection_contrast: max abs diff %.2e',
          max(abs(merge(settings[setting == 'fixed_0.01', .(collection, contrast, gene_set, FDR)],
                        camera_primary[, .(collection, contrast, gene_set, FDR_collection_contrast)],
                        by = c('collection', 'contrast', 'gene_set'))[, FDR - FDR_collection_contrast]))),
  sprintf('Candidate UTY sets not shared with naive Y-depleted: %d', nrow(cand))
)
writeLines(verification, 'outputs/camera_sensitivity/verification.txt')
cat(verification, sep = '\n')
