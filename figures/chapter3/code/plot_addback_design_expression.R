# Main figure: HT1080 add-back design, construct expression and global RNA structure.
# a, derivation of the six sequenced populations (lineage schematic; revised 22 September 2026);
# b, Western blots; c, raw construct counts; d, saved PCA; e, residual Y-expression context.
# Panels b-d reproduce RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026/code/run_analysis.R
# (Figures block) from its saved output tables; that analysis folder is not modified.
options(warn = 1)


PROJ <- 'A:/Dropbox/LOY_paralogs/projects'
THESIS <- file.path(PROJ, 'SEPT_11_2026_THESIS')
ANALYSIS <- file.path(PROJ, 'RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026')
STYLE <- file.path(THESIS, 'skills/thesis-figure-style')
OUT <- file.path(THESIS, 'figures/chapter3')
SRC_OUT <- file.path(OUT, 'source_data/addback_design_expression')
STEM <- 'ch03_main_ht1080_addback_design_expression'
dir.create(SRC_OUT, recursive = TRUE, showWarnings = FALSE)

# Cairo needs its own fontconfig registration; ragg/svglite use systemfonts below.
font_cache <- file.path(SRC_OUT, '.fontconfig_cache')
dir.create(font_cache, showWarnings = FALSE)
font_config <- file.path(SRC_OUT, 'fonts.conf')
writeLines(c('<?xml version="1.0"?>', '<!DOCTYPE fontconfig SYSTEM "fonts.dtd">',
             '<fontconfig>', sprintf('  <dir>%s/assets/fonts</dir>', STYLE),
             sprintf('  <cachedir>%s</cachedir>', font_cache), '</fontconfig>'), font_config)
Sys.setenv(FONTCONFIG_FILE = font_config)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source(file.path(STYLE, 'scripts/thesis_style.R'))
register_thesis_fonts(STYLE)
theme_set(theme_thesis(base_size = 9))
update_geom_defaults('text', list(family = 'Open Sans', colour = '#222222'))

src <- list(
  construct = file.path(ANALYSIS, 'outputs/construct_rna_abundance.csv'),
  construct_fold = file.path(ANALYSIS, 'outputs/construct_fold_over_sgCTRL.csv'),
  pca_scores = file.path(ANALYSIS, 'outputs/pca_addback_only_scores.csv'),
  pca_variance = file.path(ANALYSIS, 'outputs/pca_addback_only_variance.csv'),
  western = file.path(OUT, 'western_UTY_KDM5D.png'),
  count_table = file.path(OUT, 'source_data/addback_construct_overexpression/ch03_supp_ht1080_addback_construct_overexpression_per_library.csv'),
  y_counts = file.path(dirname(ANALYSIS), 'SEPT_2026_REANALYSIS/outputs/03_counts/counts_full.rds'),
  y_annotation = file.path(dirname(ANALYSIS), 'SEPT_2026_REANALYSIS/outputs/02_annotation/gene_annotation.csv'),
  primary_dge = file.path(ANALYSIS, 'outputs/dge_primary.rds'),
  samples = file.path(ANALYSIS, 'outputs/sample_manifest.csv')
)
stopifnot(all(file.exists(unlist(src))))
manifest <- data.frame(role = names(src), path = unlist(src))
manifest$sha256 <- vapply(manifest$path, function(p) digest::digest(p, algo = 'sha256', file = TRUE), character(1))
fwrite(manifest, file.path(SRC_OUT, 'input_manifest.csv'))

group_levels <- c('CTRL', 'LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
cols <- c(CTRL = '#3C5488', LOY = '#E64B35', LOY_LXHA = '#7E6148',
          LOY_UTY_AB = '#00A087', LOY_KDM5D_AB = '#F39B7F', LOY_KDM5D_UTY_AB = '#4DBBD5')
labs_group <- c(CTRL = 'sgCTRL', LOY = 'Y-depleted', LOY_LXHA = 'Empty vector',
                LOY_UTY_AB = 'UTY', LOY_KDM5D_AB = 'KDM5D', LOY_KDM5D_UTY_AB = 'UTY + KDM5D')
gene_cols <- c(UTY = '#00A087', KDM5D = '#8491B4')
carries <- data.frame(final_group = group_levels,
                      UTY = group_levels %in% c('LOY_UTY_AB', 'LOY_KDM5D_UTY_AB'),
                      KDM5D = group_levels %in% c('LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'))

# ---- panel a: derivation of the sequenced populations -----------------------------------------
node <- function(x, y, label, w = 1.5, h = 0.42) {
  list(annotate('rect', xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
                fill = 'white', colour = 'grey55', linewidth = 0.35),
       annotate('text', x = x, y = y, label = label, size = 2.9))
}
seg <- function(x, y, xend, yend) annotate('segment', x = x, y = y, xend = xend, yend = yend,
                                           colour = 'grey45', linewidth = 0.4)
tile_y <- 1.5; tile_top <- tile_y + 0.25
tiles <- data.frame(final_group = factor(group_levels, group_levels), x = seq_along(group_levels),
                    label = c('sgCTRL', 'Y-depleted', 'Empty\nvector', 'UTY', 'KDM5D', 'UTY +\nKDM5D'))
tiles$text_col <- ifelse(tiles$final_group %in% c('LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'), '#1a1a1a', 'white')
occ <- expand.grid(x = seq_along(group_levels), occasion = factor(1:3))
occ$xo <- occ$x + c(-0.2, 0, 0.2)[as.integer(occ$occasion)]
left_x <- 0.4

p_design <- ggplot() +
  # Parental HT1080 -> sgCTRL and the Y-depleted pool (Cas9 RNP)
  node(2.5, 5.0, 'Parental HT1080', w = 1.3) +
  seg(2.5, 4.79, 2.5, 4.45) +
  annotate('text', x = 2.6, y = 4.62, label = 'Cas9 RNP', hjust = 0, size = 2.6, fontface = 'italic', colour = 'grey30') +
  seg(1, 4.45, 4, 4.45) + seg(1, 4.45, 1, tile_top) + seg(4, 4.45, 4, 4.08) +
  node(4, 3.87, 'Y-depleted pool', w = 1.3) +
  # Y-depleted pool -> untransduced, empty vector, UTY-only and KDM5D-only
  seg(4, 3.66, 4, 3.35) + seg(2, 3.35, 4.5, 3.35) + seg(2, 3.35, 2, tile_top) +
  seg(4.5, 3.35, 4.5, 2.75) + seg(3, 2.75, 5, 2.75) +
  annotate('text', x = 4.6, y = 3.05, label = 'separate transductions', hjust = 0, size = 2.6,
           fontface = 'italic', colour = 'grey30') +
  lapply(3:5, function(x) seg(x, 2.75, x, tile_top)) +
  # Combined population descends from the established UTY-only population.
  seg(4, 1.25, 4, 0.95) + seg(4, 0.95, 6, 0.95) +
  annotate('segment', x = 6, y = 0.95, xend = 6, yend = 1.25,
           colour = 'grey45', linewidth = 0.4, arrow = grid::arrow(length = unit(1.4, 'mm'))) +
  annotate('text', x = 5, y = 0.76, label = '+ KDM5D: second transduction',
           size = 2.35, colour = 'grey30') +
  # Sequenced populations
  geom_tile(data = tiles, aes(x, tile_y, fill = final_group), width = 0.84, height = 0.5, colour = NA) +
  geom_text(data = tiles, aes(x, tile_y, label = label, colour = text_col), size = 2.8,
            fontface = 'bold', lineheight = 0.85) +
  annotate('text', x = 1:2, y = 1.1, label = c('Y-retaining', 'untransduced'), size = 2.5, colour = 'grey35') +
  # Matched RNA collection occasions
  geom_point(data = occ, aes(xo, 0.30, shape = occasion), size = 1.9, colour = 'grey25') +
  annotate('text', x = left_x, y = 0.30, label = 'Matched RNA\ncollections 1-3', hjust = 1,
           size = 2.6, lineheight = 0.9) +
  scale_fill_manual(values = cols, guide = 'none') +
  scale_colour_identity() +
  scale_shape_manual(values = c(`1` = 16, `2` = 17, `3` = 15), guide = 'none') +
  coord_cartesian(xlim = c(-0.55, 6.5), ylim = c(0.05, 5.25), expand = FALSE, clip = 'off') +
  theme_void(base_family = 'Open Sans')

# ---- panel b: unchanged source blot with white vector masks over old A/B letters only --------
wb <- png::readPNG(src$western)
wb_aspect <- dim(wb)[1] / dim(wb)[2]
p_western <- ggplot() +
  annotation_raster(wb, xmin = 0, xmax = 1, ymin = 0, ymax = wb_aspect, interpolate = TRUE) +
  # Pixel bounds in the 667 x 467 source. Both masks lie entirely in the blank left margin.
  annotate('rect', xmin = 9/667, xmax = 40/667, ymin = (467-33)/667,
           ymax = (467-1)/667, fill = 'white', colour = NA) +
  annotate('rect', xmin = 9/667, xmax = 40/667, ymin = (467-275)/667,
           ymax = (467-243)/667, fill = 'white', colour = NA) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, wb_aspect), expand = FALSE) +
  theme_void(base_family = 'Open Sans')

# ---- panel c: construct RNA abundance (as run_analysis.R) --------------------------------------
construct <- fread(src$construct, data.table = FALSE)
construct$final_group <- factor(construct$final_group, group_levels)
construct$gene <- factor(construct$gene, c('UTY', 'KDM5D'))
construct$occasion <- factor(construct$replicate_block)
construct$detected <- ifelse(construct$count > 0, 'reads', 'none')
construct_labels <- c(CTRL = 'sgCTRL', LOY = 'Y-depleted', LOY_LXHA = 'Vector', LOY_UTY_AB = 'UTY',
                      LOY_KDM5D_AB = 'KDM5D', LOY_KDM5D_UTY_AB = 'UTY +\nKDM5D')
construct_map <- data.frame(final_group = factor(rep(group_levels, 2), group_levels),
                            gene = factor(rep(c('UTY', 'KDM5D'), each = length(group_levels)), c('UTY', 'KDM5D')),
                            present = c(carries$UTY, carries$KDM5D))
dodge <- position_dodge(width = 0.55)
c_mean <- aggregate(logCPM ~ final_group + gene, construct, mean)
ctrl_ref <- setNames(c_mean$logCPM[c_mean$final_group == 'CTRL'], as.character(c_mean$gene[c_mean$final_group == 'CTRL']))
c_fold <- merge(c_mean, construct_map[construct_map$present, ], by = c('final_group', 'gene'))
c_fold$fold <- 2^(c_fold$logCPM - ctrl_ref[as.character(c_fold$gene)])
c_fold$label <- sprintf('%.0f×', c_fold$fold)
c_max <- aggregate(logCPM ~ final_group + gene, construct, max)
c_fold$y <- c_max$logCPM[match(paste(c_fold$final_group, c_fold$gene), paste(c_max$final_group, c_max$gene))] + 0.9
c_fold$x <- as.numeric(c_fold$final_group) + ifelse(c_fold$gene == 'UTY', -0.1375, 0.1375)
saved_fold <- fread(src$construct_fold, data.table = FALSE)
chk <- merge(saved_fold, c_fold, by = c('final_group', 'gene'))
stopifnot(nrow(chk) == nrow(saved_fold), max(abs(chk$fold.x - chk$fold.y)) < 1e-8)
floor_lcpm <- min(construct$logCPM)
ref_df <- data.frame(gene = factor(names(ctrl_ref), c('UTY', 'KDM5D')), y = unname(ctrl_ref))

# Raw-count panel transferred from supplementary panel a, using its saved table.
per_lib <- fread(src$count_table)
stopifnot(nrow(per_lib) == 36, uniqueN(per_lib$sample_id) == 18)
count_check <- merge(as.data.table(construct)[, .(sample_id, gene, original_count = count)],
                     per_lib[, .(sample_id, gene, count)], by = c('sample_id', 'gene'))
stopifnot(nrow(count_check) == 36, all(count_check$original_count == count_check$count))
fwrite(per_lib, file.path(SRC_OUT, 'panel_c_raw_construct_counts.csv'))
fmt_n <- function(x) formatC(x, format = 'd', big.mark = ',')
reads_col <- 'Assigned\nreads (M)'
tab <- rbind(
  per_lib[, .(final_group, replicate_block, column = gene, value = as.numeric(count),
              label = fmt_n(count), outline = carries_construct)],
  unique(per_lib[, .(final_group, replicate_block, column = reads_col, value = NA_real_,
                     label = sprintf('%.1f', assigned_reads / 1e6), outline = FALSE)])
)
tab[, column := factor(column, levels = c('UTY', 'KDM5D', reads_col))]
tab[, group := factor(final_group, levels = group_levels, labels = labs_group[group_levels])]
tab[, row := factor(replicate_block, levels = 3:1)]
tab[, text_col := fifelse(is.na(value), '#222222', fifelse(log10(value + 1) > 3.6, '#111111', 'white'))]
p_counts <- ggplot(tab, aes(column, row)) +
  geom_tile(aes(fill = log10(value + 1)), colour = 'white', linewidth = 0.6) +
  geom_tile(data = tab[outline == TRUE], fill = NA, colour = '#111111', linewidth = 0.5,
            width = 0.93, height = 0.88) +
  geom_text(aes(label = label, colour = text_col), size = 2.7, family = 'Open Sans') +
  scale_colour_identity() +
  scale_fill_viridis_c(option = 'magma', na.value = 'grey96', limits = c(0, 6),
                       breaks = c(0, 2, 4, 6), labels = c('0', '100', '10k', '1M'), name = 'Raw reads') +
  facet_grid(group ~ ., switch = 'y') +
  scale_x_discrete(position = 'top', expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 7, colour = 'grey45'),
        strip.placement = 'outside', strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.background = element_blank(), panel.spacing.y = unit(1.1, 'mm'),
        legend.position = 'bottom', legend.key.width = unit(9, 'mm'), legend.key.height = unit(2.2, 'mm'),
        legend.title = element_text(vjust = 0.9))

# ---- panel d: PCA of the four plasmid-bearing populations, with spread ellipses -------------------
# Scores are the retained plasmid-only PCA (TMM log2 CPM, prior count 2, 23,992 primary genes, centred,
# unscaled). Ellipse construction follows thesis/Ben_Thesis/Aim2_addback_pca/code/plot_pca_populations.R:
# covariance of each population's three libraries scaled to 1.5 SD, minor axis floored at 20% of the major.
plasmid_groups <- c('LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
pca_df <- as.data.table(fread(src$pca_scores))
pve <- fread(src$pca_variance, data.table = FALSE)$percent_variance
stopifnot(nrow(pca_df) == 12, setequal(pca_df$final_group, plasmid_groups))
pca_df[, final_group := factor(final_group, plasmid_groups)]
cov_ellipse <- function(x, y, scale = 1.5, floor_frac = 0.2, n = 120) {
  mu <- c(mean(x), mean(y)); e <- eigen(cov(cbind(x, y)))
  radii <- sqrt(pmax(e$values, 0)); radii[2] <- max(radii[2], floor_frac * radii[1])
  t <- seq(0, 2 * pi, length.out = n)
  pts <- t(mu + t(cbind(cos(t), sin(t)) %*% diag(radii * scale) %*% t(e$vectors)))
  data.table(PC1 = pts[, 1], PC2 = pts[, 2])
}
ell <- pca_df[, cov_ellipse(PC1, PC2), by = final_group]
p_pca <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_polygon(data = ell, aes(fill = final_group, colour = final_group, group = final_group),
               alpha = 0.12, linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(colour = final_group, shape = factor(replicate_block)), size = 3) +
  scale_colour_manual(values = cols[plasmid_groups], labels = labs_group[plasmid_groups], name = NULL) +
  scale_fill_manual(values = cols[plasmid_groups], guide = 'none') +
  scale_shape_manual(values = c(`1` = 16, `2` = 17, `3` = 15), name = 'Occasion') +
  labs(x = sprintf('PC1 (%.1f%%)', pve[1]), y = sprintf('PC2 (%.1f%%)', pve[2])) +
  theme(legend.position = 'bottom', legend.box = 'vertical', legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.5, 'mm')) +
  guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 16, size = 2.4)),
         shape = guide_legend(order = 2, override.aes = list(size = 2.2)))
fwrite(ell, file.path(SRC_OUT, 'pca_plasmid_ellipses.csv'))

# ---- panel e: original residual-Y heatmap, transposed for a compact full-width panel ---------
y_genes <- c('RPS4Y1', 'ZFY', 'USP9Y', 'DDX3Y', 'UTY', 'KDM5D', 'EIF1AY')
y_counts <- readRDS(src$y_counts)
y_annot <- fread(src$y_annotation)
y_samples <- fread(src$samples)
y_dge <- readRDS(src$primary_dge)
stopifnot(identical(colnames(y_counts), y_samples$sample_id),
          identical(colnames(y_dge), y_samples$sample_id),
          nrow(y_samples) == 18, all(table(y_samples$final_group) == 3))
y_selected <- y_annot[chrom == 'chrY' & gene_name %in% y_genes]
stopifnot(nrow(y_selected) == 7, !anyDuplicated(y_selected$gene_name))
y_selected <- y_selected[match(y_genes, gene_name)]
y_eff <- y_dge$samples$lib.size * y_dge$samples$norm.factors
y_lcpm <- edgeR::cpm(y_counts[y_selected$gene_id, , drop = FALSE],
                     lib.size = y_eff, log = TRUE, prior.count = 2)
rownames(y_lcpm) <- y_selected$gene_name
y_sample_long <- melt(as.data.table(y_lcpm, keep.rownames = 'gene'),
                      id.vars = 'gene', variable.name = 'sample_id', value.name = 'logCPM')
y_sample_long[, sample_id := as.character(sample_id)]
y_sample_long <- merge(y_sample_long, y_samples[, .(sample_id, final_group, replicate_block)], by = 'sample_id')
# Confirm the construct genes use exactly the same normalisation as panel c's source.
y_check <- merge(y_sample_long[gene %in% c('UTY', 'KDM5D')],
                 as.data.table(construct)[, .(sample_id, gene, saved_logCPM = logCPM)],
                 by = c('sample_id', 'gene'))
stopifnot(nrow(y_check) == 36, max(abs(y_check$logCPM - y_check$saved_logCPM)) < 1e-8)
y_summary <- y_sample_long[, .(mean_logCPM = mean(logCPM), n_collections = .N), by = .(gene, final_group)]
stopifnot(nrow(y_summary) == 42, all(y_summary$n_collections == 3))
fwrite(y_sample_long, file.path(SRC_OUT, 'panel_e_y_expression_per_library.csv'))
fwrite(y_summary, file.path(SRC_OUT, 'panel_e_y_expression_group_means.csv'))
y_summary[, gene := factor(gene, levels = y_genes)]
y_summary[, group := factor(final_group, levels = rev(group_levels), labels = labs_group[rev(group_levels)])]
y_summary[, text_colour := fifelse(mean_logCPM > 7, '#111111', 'white')]
p_y <- ggplot(y_summary, aes(gene, group)) +
  geom_tile(aes(fill = mean_logCPM), colour = 'white', linewidth = 0.5) +
  geom_text(aes(label = sprintf('%.1f', mean_logCPM), colour = text_colour), size = 2.6) +
  scale_colour_identity() +
  scale_fill_viridis_c(option = 'magma', name = 'Mean log2 CPM') +
  scale_x_discrete(position = 'top', expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_text(colour = '#222222'), legend.position = 'right',
        legend.key.height = unit(6, 'mm'), legend.key.width = unit(3, 'mm'))

# ---- assemble --------------------------------------------------------------------------------
# Compact five-panel layout: design/blots, raw counts/PCA, residual-Y heatmap.
top <- (p_design + labs(tag = 'a')) + (p_western + labs(tag = 'b')) + plot_layout(widths = c(1.85, 1))
bottom <- (p_counts + labs(tag = 'c')) + (p_pca + labs(tag = 'd')) + plot_layout(widths = c(1, 1.25))
fig <- top / bottom / (p_y + labs(tag = 'e')) + plot_layout(heights = c(0.95, 1.45, 0.64))
W <- 10.5; H <- 10.3

ggsave(file.path(OUT, paste0(STEM, '.png')), fig, width = W, height = H, dpi = 300, bg = 'white', device = ragg::agg_png)
ggsave(file.path(OUT, paste0(STEM, '.pdf')), fig, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT, paste0(STEM, '.svg')), fig, width = W, height = H, device = svglite::svglite)
fwrite(c_fold[, c('final_group', 'gene', 'logCPM', 'fold')], file.path(SRC_OUT, 'construct_fold_over_sgCTRL.csv'))

# ---- supplementary: optional refresh of the all-population PCA ----------------------------------
if (!('--main-only' %in% commandArgs(trailingOnly = TRUE))) {
all_scores <- file.path(ANALYSIS, 'outputs/pca_scores.csv')
all_var <- file.path(ANALYSIS, 'outputs/pca_variance.csv')
pca_all <- as.data.table(fread(all_scores))
pve_all <- fread(all_var, data.table = FALSE)$percent_variance
stopifnot(nrow(pca_all) == 18)
pca_all[, final_group := factor(final_group, group_levels)]
ell_all <- pca_all[, cov_ellipse(PC1, PC2), by = final_group]
p_pca_all <- ggplot(pca_all, aes(PC1, PC2)) +
  geom_polygon(data = ell_all, aes(fill = final_group, colour = final_group, group = final_group),
               alpha = 0.12, linewidth = 0.3, show.legend = FALSE) +
  geom_point(aes(colour = final_group, shape = factor(replicate_block)), size = 2.4) +
  scale_colour_manual(values = cols, labels = labs_group, name = NULL) +
  scale_fill_manual(values = cols, guide = 'none') +
  scale_shape_manual(values = c(`1` = 16, `2` = 17, `3` = 15), name = 'Collection occasion') +
  labs(x = sprintf('PC1 (%.1f%%)', pve_all[1]), y = sprintf('PC2 (%.1f%%)', pve_all[2])) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 2.4)),
         shape = guide_legend(order = 2))
SUPP <- 'ch03_supp_ht1080_addback_pca_populations'
ggsave(file.path(OUT, paste0(SUPP, '.png')), p_pca_all, width = 6.2, height = 3.9, dpi = 300, bg = 'white')
ggsave(file.path(OUT, paste0(SUPP, '.pdf')), p_pca_all, width = 6.2, height = 3.9, device = cairo_pdf)
ggsave(file.path(OUT, paste0(SUPP, '.svg')), p_pca_all, width = 6.2, height = 3.9, device = svglite::svglite)
fwrite(ell_all, file.path(SRC_OUT, 'pca_all_populations_ellipses.csv'))
}
writeLines(capture.output(sessionInfo()), file.path(SRC_OUT, 'sessionInfo.txt'))
cat('done\n')
