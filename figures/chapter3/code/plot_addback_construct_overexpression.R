# Supplementary figure: scale of UTY and KDM5D overexpression in the HT1080 add-back populations.
# a, raw read counts per library; b, within-library transcriptome rank of UTY and KDM5D.
# log2 CPM uses the primary TMM effective library sizes and prior count 2, as in
# RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026/code/run_analysis.R (construct abundance block).
options(warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

PROJ <- 'A:/Dropbox/LOY_paralogs/projects'
THESIS <- file.path(PROJ, 'SEPT_11_2026_THESIS')
ADDBACK <- file.path(PROJ, 'RNAseq/analysis/HT1080_ADDBACK')
ANALYSIS <- file.path(ADDBACK, 'THESIS_ADDBACK_2026')
STYLE <- file.path(THESIS, 'skills/thesis-figure-style')
OUT <- file.path(THESIS, 'figures/chapter3')
SRC_OUT <- file.path(OUT, 'source_data/addback_construct_overexpression')
STEM <- 'ch03_supp_ht1080_addback_construct_overexpression'
dir.create(SRC_OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(STYLE, 'scripts/thesis_style.R'))
register_thesis_fonts(STYLE)
theme_set(theme_thesis(base_size = 9))

src <- list(
  counts = file.path(ADDBACK, 'SEPT_2026_REANALYSIS/outputs/03_counts/counts_full.rds'),
  annotation = file.path(ADDBACK, 'SEPT_2026_REANALYSIS/outputs/02_annotation/gene_annotation.csv'),
  dge = file.path(ANALYSIS, 'outputs/dge_primary.rds'),
  samples = file.path(ANALYSIS, 'outputs/sample_manifest.csv'),
  construct_check = file.path(ANALYSIS, 'outputs/construct_rna_abundance.csv')
)
stopifnot(all(file.exists(unlist(src))))
manifest <- data.frame(role = names(src), path = unlist(src))
manifest$sha256 <- vapply(manifest$path, function(p) digest::digest(p, algo = 'sha256', file = TRUE), character(1))
fwrite(manifest, file.path(SRC_OUT, 'input_manifest.csv'))

counts <- readRDS(src$counts)
annot <- fread(src$annotation, data.table = FALSE)
dge <- readRDS(src$dge)
samples <- fread(src$samples, data.table = FALSE)
stopifnot(identical(colnames(counts), samples$sample_id), identical(colnames(dge), samples$sample_id))

construct_ids <- annot$gene_id[annot$gene_name %in% c('UTY', 'KDM5D') & annot$chrom == 'chrY']
stopifnot(length(construct_ids) == 2)
names(construct_ids) <- annot$gene_name[match(construct_ids, annot$gene_id)]

# Gene universe: the primary expressed genes plus UTY and KDM5D.
univ <- c(dge$genes$gene_id, unname(construct_ids))
stopifnot(!anyDuplicated(univ), all(univ %in% rownames(counts)))
lib_eff <- dge$samples$lib.size * dge$samples$norm.factors
lcpm <- cpm(counts[univ, , drop = FALSE], lib.size = lib_eff, log = TRUE, prior.count = 2)
n_univ <- nrow(lcpm)

# Reproduce the saved construct abundance exactly.
chk <- fread(src$construct_check, data.table = FALSE)
ri <- match(construct_ids[chk$gene], rownames(lcpm))
ci <- match(chk$sample_id, colnames(lcpm))
stopifnot(max(abs(lcpm[cbind(ri, ci)] - chk$logCPM)) < 1e-8,
          all(counts[cbind(match(construct_ids[chk$gene], rownames(counts)), ci)] == chk$count))

group_levels <- c('CTRL', 'LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
group_labels <- c(CTRL = 'sgCTRL', LOY = 'Y-\ndepleted', LOY_LXHA = 'Vector', LOY_UTY_AB = 'UTY',
                  LOY_KDM5D_AB = 'KDM5D', LOY_KDM5D_UTY_AB = 'UTY +\nKDM5D')
carries <- list(UTY = c('LOY_UTY_AB', 'LOY_KDM5D_UTY_AB'), KDM5D = c('LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'))
samples$assigned_reads <- colSums(counts)[samples$sample_id]

# ---- per-library counts and ranks ------------------------------------------------------------
per_lib <- rbindlist(lapply(names(construct_ids), function(g) {
  id <- construct_ids[[g]]
  rbindlist(lapply(colnames(lcpm), function(s) {
    v <- lcpm[, s]
    data.table(sample_id = s, gene = g, count = counts[id, s], logCPM = v[[id]],
               rank = sum(v > v[[id]]) + 1L, n_genes = n_univ,
               percent_genes_below = 100 * mean(v[names(v) != id] < v[[id]]))
  }))
}))
per_lib <- merge(per_lib, as.data.table(samples)[, .(sample_id, final_group, replicate_block, assigned_reads)],
                 by = 'sample_id')
per_lib[, carries_construct := mapply(function(g, grp) grp %in% carries[[g]], gene, final_group)]
setorder(per_lib, gene, final_group, replicate_block)
fwrite(per_lib, file.path(SRC_OUT, paste0(STEM, '_per_library.csv')))

summ <- per_lib[, .(count_min = min(count), count_max = max(count), mean_logCPM = mean(logCPM),
                    rank_best = min(rank), rank_worst = max(rank), mean_rank = mean(rank),
                    mean_percent_below = mean(percent_genes_below), carries_construct = carries_construct[1]),
                by = .(gene, final_group)]
fwrite(summ, file.path(SRC_OUT, paste0(STEM, '_summary.csv')))
print(summ)

# ---- panel a: raw count table ------------------------------------------------------------------
fmt_n <- function(x) formatC(x, format = 'd', big.mark = ',')
reads_col <- 'Assigned\nreads (M)'
tab <- rbind(
  per_lib[, .(final_group, replicate_block, column = gene, value = as.numeric(count), label = fmt_n(count),
              outline = carries_construct)],
  unique(per_lib[, .(final_group, replicate_block, column = reads_col, value = NA_real_,
                     label = sprintf('%.1f', assigned_reads / 1e6), outline = FALSE)])
)
tab[, column := factor(column, levels = c('UTY', 'KDM5D', reads_col))]
tab[, group := factor(final_group, levels = group_levels, labels = sub('\n', ' ', sub('-\n', '-', group_labels[group_levels])))]
tab[, row := factor(replicate_block, levels = 3:1)]
tab[, text_col := fifelse(is.na(value), '#222222', fifelse(log10(value + 1) > 3.6, '#111111', 'white'))]

pa <- ggplot(tab, aes(column, row)) +
  geom_tile(aes(fill = log10(value + 1)), colour = 'white', linewidth = 0.6) +
  geom_tile(data = tab[outline == TRUE], fill = NA, colour = '#111111', linewidth = 0.5, width = 0.93, height = 0.88) +
  geom_text(aes(label = label, colour = text_col), size = 2.5, family = 'Open Sans') +
  scale_colour_identity() +
  scale_fill_viridis_c(option = 'magma', na.value = 'grey96', limits = c(0, 6), breaks = c(0, 2, 4, 6),
                       labels = c('0', '100', '10k', '1M'), name = 'Raw reads') +
  facet_grid(group ~ ., switch = 'y') +
  scale_x_discrete(position = 'top', expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, tag = 'a') +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 6.5, colour = 'grey45'),
        strip.placement = 'outside', strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.background = element_blank(), panel.spacing.y = unit(1.2, 'mm'),
        legend.position = 'bottom', legend.key.width = unit(9, 'mm'), legend.key.height = unit(2.2, 'mm'),
        legend.title = element_text(vjust = 0.9))

# ---- panel b: transcriptome rank ---------------------------------------------------------------
other <- setdiff(rownames(lcpm), construct_ids)
dist <- rbindlist(lapply(group_levels, function(g) {
  ids <- samples$sample_id[samples$final_group == g]
  data.table(final_group = g, logCPM = as.vector(lcpm[other, ids]))
}))
dist[, group := factor(final_group, levels = group_levels)]
dodge <- c(UTY = -0.17, KDM5D = 0.17)
pts <- copy(per_lib)[, x := match(final_group, group_levels) + dodge[gene]]
gene_cols <- c(UTY = '#00A087', KDM5D = '#8491B4')

ymax <- per_lib[, .(ymax = max(logCPM)), by = .(gene, final_group)]
lab <- merge(summ[carries_construct == TRUE | final_group == 'CTRL'], ymax, by = c('gene', 'final_group'))
lab[, x := match(final_group, group_levels) + dodge[gene]]
lab[, top := 100 - mean_percent_below]
lab[, text := sprintf('rank %s\ntop %s%%', fmt_n(round(mean_rank)),
                      ifelse(top < 0.1, formatC(top, format = 'f', digits = 2), formatC(top, format = 'f', digits = 1)))]

pb <- ggplot() +
  geom_violin(data = dist, aes(group, logCPM), fill = 'grey88', colour = NA, width = 0.85, scale = 'width') +
  geom_boxplot(data = dist, aes(group, logCPM), width = 0.07, outlier.shape = NA, fill = 'white',
               colour = 'grey55', linewidth = 0.3) +
  geom_point(data = pts, aes(x, logCPM, colour = gene, shape = factor(replicate_block)), size = 1.8, stroke = 0.5) +
  geom_label(data = lab, aes(x, ymax + 0.5, label = text, colour = gene), size = 2.2, family = 'Open Sans',
             lineheight = 0.9, vjust = 0, fill = alpha('white', 0.9), linewidth = 0,
             label.padding = unit(0.6, 'mm'), show.legend = FALSE) +
  scale_colour_manual(values = gene_cols, breaks = c('UTY', 'KDM5D'), name = NULL) +
  scale_shape_manual(values = c(`1` = 16, `2` = 15, `3` = 17), name = 'Occasion') +
  scale_x_discrete(labels = group_labels, expand = expansion(add = 0.75)) +
  scale_y_continuous(limits = c(min(lcpm), 17), breaks = seq(-4, 16, 4)) +
  labs(x = NULL, y = sprintf('log2 CPM (%s genes per library)', fmt_n(n_univ)), tag = 'b') +
  theme(legend.position = 'bottom', panel.grid.major.y = element_line(colour = 'grey93', linewidth = 0.3)) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 2.2)))

p <- pa + pb + plot_layout(widths = c(0.85, 1.6))
W <- 180 / 25.4; H <- 125 / 25.4
ggsave(file.path(OUT, paste0(STEM, '.png')), p, width = W, height = H, dpi = 300, bg = 'white')
ggsave(file.path(OUT, paste0(STEM, '.pdf')), p, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT, paste0(STEM, '.svg')), p, width = W, height = H, device = svglite::svglite)
writeLines(capture.output(sessionInfo()), file.path(SRC_OUT, 'sessionInfo.txt'))
cat('done\n')
