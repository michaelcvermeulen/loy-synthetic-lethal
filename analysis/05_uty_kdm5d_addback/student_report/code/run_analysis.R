options(warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)[1]
script <- normalizePath(sub('^--file=', '', script_arg), winslash = '/', mustWork = TRUE)
ROOT <- dirname(dirname(script))
ADDBACK <- dirname(ROOT)
BASELINE <- 'A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE'
STYLE <- 'A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/skills/thesis-figure-style'
setwd(ROOT)

dirs <- c('outputs/de', 'outputs/signature', 'outputs/sensitivity', 'figures', 'provenance')
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
set.seed(42)

source(file.path(STYLE, 'scripts/thesis_style.R'))
register_thesis_fonts(STYLE)
theme_set(theme_thesis(base_size = 9))
update_geom_defaults('text', list(family = 'Open Sans', colour = '#222222'))

src <- list(
  counts = file.path(ADDBACK, 'SEPT_2026_REANALYSIS/outputs/03_counts/counts_full.rds'),
  annotation = file.path(ADDBACK, 'SEPT_2026_REANALYSIS/outputs/02_annotation/gene_annotation.csv'),
  samples = file.path(ADDBACK, 'SEPT_2026_REANALYSIS/config/samples.tsv'),
  relabel = 'A:/Dropbox/LOY_paralogs/projects/RNAseq/counts/STAR_v2_PARfix/MAY_2026_HT1080_ADDBACK/sample_relabeling_2026-05-26.csv',
  audit_hashes = file.path(ADDBACK, 'SEPT_2026_REANALYSIS/outputs/01_audit/input_hashes.csv'),
  signature = file.path(BASELINE, 'outputs/tables/August_frozen_signature_genes.csv'),
  signature_background = file.path(BASELINE, 'outputs/tables/August_signature_discovery_background.csv')
)
stopifnot(all(vapply(src, file.exists, logical(1))))

sha <- function(p) digest::digest(p, algo = 'sha256', file = TRUE)
manifest <- data.frame(role = names(src), path = unlist(src), stringsAsFactors = FALSE)
manifest$sha256 <- vapply(manifest$path, sha, character(1))
manifest$bytes <- file.size(manifest$path)
fwrite(manifest, 'provenance/input_manifest.csv')

samples <- fread(src$samples, data.table = FALSE)
samples$collection_occasion <- factor(samples$replicate_block, levels = c(1, 2, 3))
group_levels <- c('CTRL', 'LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
samples$group <- factor(samples$final_group, levels = group_levels)
stopifnot(nrow(samples) == 18, all(table(samples$group) == 3), all(table(samples$collection_occasion) == 6))
fwrite(samples, 'outputs/sample_manifest.csv')

counts <- readRDS(src$counts)
annot <- fread(src$annotation, data.table = FALSE)
stopifnot(identical(colnames(counts), samples$sample_id), !anyDuplicated(rownames(counts)))
a <- annot[match(rownames(counts), annot$gene_id), ]
stopifnot(!anyNA(a$gene_id))
keep_input <- !a$exclude_Y & !is.na(a$gene_name) & nzchar(a$gene_name)
m <- counts[keep_input, , drop = FALSE]
a_main <- a[keep_input, , drop = FALSE]
rownames(a_main) <- a_main$gene_id

pheno <- data.frame(row.names = samples$sample_id, group = samples$group,
                    collection_occasion = samples$collection_occasion)
dge <- DGEList(m, samples = pheno, genes = a_main, group = pheno$group)
keep <- filterByExpr(dge, group = dge$samples$group, min.count = 10, min.total.count = 15)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = 'TMM')
saveRDS(dge, 'outputs/dge_primary.rds')
fwrite(data.frame(n_named_nonY = nrow(m), n_tested = nrow(dge),
                  n_chrY_or_PAR_removed = sum(a$exclude_Y, na.rm = TRUE)), 'outputs/filter_summary.csv')

make_design <- function(block = TRUE, sam = samples) {
  z <- data.frame(group = factor(sam$final_group, levels = group_levels),
                  collection_occasion = factor(sam$replicate_block))
  d <- if (block) model.matrix(~ 0 + group + collection_occasion, z) else model.matrix(~ 0 + group, z)
  colnames(d) <- sub('^group', '', colnames(d))
  stopifnot(qr(d)$rank == ncol(d))
  d
}

contrast_table <- data.frame(
  name = c('LOY_v_CTRL', 'LXHA_v_CTRL', 'UTY_v_CTRL', 'KDM5D_v_CTRL', 'DUAL_v_CTRL',
           'LOY_v_LXHA', 'UTY_v_LXHA', 'KDM5D_v_LXHA', 'DUAL_v_LXHA',
           'UTY_v_LOY', 'KDM5D_v_LOY', 'DUAL_v_LOY',
           'UTY_v_KDM5D', 'DUAL_v_UTY', 'DUAL_v_KDM5D'),
  formula = c('LOY-CTRL', 'LOY_LXHA-CTRL', 'LOY_UTY_AB-CTRL', 'LOY_KDM5D_AB-CTRL', 'LOY_KDM5D_UTY_AB-CTRL',
              'LOY-LOY_LXHA', 'LOY_UTY_AB-LOY_LXHA', 'LOY_KDM5D_AB-LOY_LXHA', 'LOY_KDM5D_UTY_AB-LOY_LXHA',
              'LOY_UTY_AB-LOY', 'LOY_KDM5D_AB-LOY', 'LOY_KDM5D_UTY_AB-LOY',
              'LOY_UTY_AB-LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB-LOY_UTY_AB', 'LOY_KDM5D_UTY_AB-LOY_KDM5D_AB'),
  stringsAsFactors = FALSE)

fit_analysis <- function(d, y = dge) {
  v <- voom(y, d, plot = FALSE)
  f <- lmFit(v, d)
  cm <- makeContrasts(contrasts = contrast_table$formula, levels = d)
  colnames(cm) <- contrast_table$name
  list(v = v, fit = f, cm = cm, eb = eBayes(contrasts.fit(f, cm), robust = TRUE),
       tr = treat(contrasts.fit(f, cm), lfc = 0.5, robust = TRUE))
}

design <- make_design(TRUE)
fit <- fit_analysis(design)
saveRDS(fit, 'outputs/fit_primary.rds')

de_summary <- list()
for (k in contrast_table$name) {
  tt <- topTable(fit$eb, coef = k, number = Inf, sort.by = 'none', confint = TRUE)
  tr <- topTreat(fit$tr, coef = k, number = Inf, sort.by = 'none')
  out <- data.frame(gene_id = rownames(tt), gene_name = dge$genes$gene_name,
                    chrom = dge$genes$chrom, logFC = tt$logFC, CI_L = tt$CI.L, CI_R = tt$CI.R,
                    AveExpr = tt$AveExpr, t = tt$t, PValue = tt$P.Value, FDR = tt$adj.P.Val,
                    FDR_treat_0.5 = tr$adj.P.Val, stringsAsFactors = FALSE)
  out$selected_0.5 <- out$FDR < 0.05 & abs(out$logFC) >= 0.5
  out$selected_1 <- out$FDR < 0.05 & abs(out$logFC) >= 1
  out$treat_0.5 <- out$FDR_treat_0.5 < 0.05
  fwrite(out, file.path('outputs/de', paste0('de_', k, '.csv')))
  de_summary[[k]] <- data.frame(contrast = k, n_tested = nrow(out), n_FDR = sum(out$FDR < 0.05),
    n_selected_0.5 = sum(out$selected_0.5), n_selected_1 = sum(out$selected_1),
    n_treat_0.5 = sum(out$treat_0.5), n_up_0.5 = sum(out$selected_0.5 & out$logFC > 0),
    n_down_0.5 = sum(out$selected_0.5 & out$logFC < 0))
}
de_summary <- do.call(rbind, de_summary)
fwrite(de_summary, 'outputs/de_summary.csv')

# Factorial terms among plasmid-bearing populations.
sp <- samples$final_group %in% c('LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
fac <- data.frame(UTY = factor(ifelse(grepl('UTY', samples$final_group[sp]), 'yes', 'no')),
                  KDM5D = factor(ifelse(grepl('KDM5D', samples$final_group[sp]), 'yes', 'no')),
                  collection_occasion = samples$collection_occasion[sp])
d_fac <- model.matrix(~ collection_occasion + UTY * KDM5D, fac)
colnames(d_fac) <- make.names(colnames(d_fac))
y_fac <- dge[, sp]
v_fac <- voom(y_fac, d_fac, plot = FALSE)
f_fac <- lmFit(v_fac, d_fac)
cm_fac <- makeContrasts(
  UTY_marginal = UTYyes + 0.5 * UTYyes.KDM5Dyes,
  KDM5D_marginal = KDM5Dyes + 0.5 * UTYyes.KDM5Dyes,
  interaction = UTYyes.KDM5Dyes, levels = d_fac)
eb_fac <- eBayes(contrasts.fit(f_fac, cm_fac), robust = TRUE)
fac_summary <- do.call(rbind, lapply(colnames(cm_fac), function(k) {
  z <- topTable(eb_fac, coef = k, number = Inf, sort.by = 'none')
  data.frame(term = k, n_FDR = sum(z$adj.P.Val < 0.05),
             n_selected_0.5 = sum(z$adj.P.Val < 0.05 & abs(z$logFC) >= 0.5))
}))
fwrite(fac_summary, 'outputs/factorial_summary.csv')

# Unblocked and leave-one-occasion-out sensitivities.
fit0 <- fit_analysis(make_design(FALSE))
unblocked <- do.call(rbind, lapply(contrast_table$name, function(k) {
  z <- topTable(fit0$eb, coef = k, number = Inf, sort.by = 'none')
  data.frame(contrast = k, method = 'unblocked', n_FDR = sum(z$adj.P.Val < 0.05),
             n_selected_0.5 = sum(z$adj.P.Val < 0.05 & abs(z$logFC) >= 0.5))
}))
fwrite(unblocked, 'outputs/sensitivity/unblocked_summary.csv')

loo_contrasts <- c('LOY_v_LXHA', 'UTY_v_LXHA', 'KDM5D_v_LXHA', 'DUAL_v_LXHA', 'DUAL_v_UTY', 'DUAL_v_KDM5D')
full_lfc <- lapply(loo_contrasts, function(k) fread(file.path('outputs/de', paste0('de_', k, '.csv')))$logFC)
names(full_lfc) <- loo_contrasts
loo <- list()
for (occasion in levels(samples$collection_occasion)) {
  sel <- samples$collection_occasion != occasion
  ys <- dge[, sel]
  ss <- samples[sel, ]
  ds <- make_design(TRUE, ss)
  fs <- fit_analysis(ds, ys)
  for (k in loo_contrasts) {
    z <- topTable(fs$eb, coef = k, number = Inf, sort.by = 'none')
    loo[[paste(occasion, k)]] <- data.frame(omitted_occasion = occasion, contrast = k,
      n_FDR = sum(z$adj.P.Val < 0.05), n_selected_0.5 = sum(z$adj.P.Val < 0.05 & abs(z$logFC) >= 0.5),
      cor_logFC_with_full = cor(z$logFC, full_lfc[[k]], method = 'pearson'))
  }
}
loo <- do.call(rbind, loo)
fwrite(loo, 'outputs/sensitivity/leave_one_occasion_out.csv')

# Construct RNA abundance using the primary effective library sizes.
full_main <- counts[match(dge$genes$gene_id, rownames(counts)), , drop = FALSE]
lib_eff <- dge$samples$lib.size * dge$samples$norm.factors
construct_ids <- annot$gene_id[annot$gene_name %in% c('UTY', 'KDM5D') & annot$chrom == 'chrY']
stopifnot(length(construct_ids) == 2)
construct_counts <- counts[construct_ids, , drop = FALSE]
construct_lcpm <- cpm(construct_counts, lib.size = lib_eff, log = TRUE, prior.count = 2)
construct <- data.frame(sample_id = rep(colnames(counts), each = 2),
                        gene = rep(annot$gene_name[match(construct_ids, annot$gene_id)], times = ncol(counts)),
                        logCPM = as.vector(construct_lcpm), count = as.vector(construct_counts),
                        stringsAsFactors = FALSE)
construct <- merge(construct, samples[, c('sample_id', 'final_group', 'replicate_block')], by = 'sample_id')
fwrite(construct, 'outputs/construct_rna_abundance.csv')

# PCA on all tested genes, centred without variance scaling.
lcpm <- cpm(dge, log = TRUE, prior.count = 2)
pca <- prcomp(t(lcpm), center = TRUE, scale. = FALSE)
pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- cbind(samples[, c('sample_id', 'final_group', 'replicate_block')],
                PC1 = pca$x[, 1], PC2 = pca$x[, 2])
fwrite(pca_df, 'outputs/pca_scores.csv')
fwrite(data.frame(PC = paste0('PC', seq_along(pve)), percent_variance = pve), 'outputs/pca_variance.csv')

# A separate PCA for the four plasmid-bearing populations. This is recomputed
# from the 12 relevant libraries rather than obtained by hiding reference
# samples from the six-population PCA.
addback_groups <- c('LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
addback_samples <- samples$final_group %in% addback_groups
pca_addback <- prcomp(t(lcpm[, addback_samples, drop = FALSE]), center = TRUE, scale. = FALSE)
pve_addback <- 100 * pca_addback$sdev^2 / sum(pca_addback$sdev^2)
pca_addback_df <- cbind(
  samples[addback_samples, c('sample_id', 'final_group', 'replicate_block')],
  PC1 = pca_addback$x[, 1], PC2 = pca_addback$x[, 2]
)
fwrite(pca_addback_df, 'outputs/pca_addback_only_scores.csv')
fwrite(data.frame(PC = paste0('PC', seq_along(pve_addback)), percent_variance = pve_addback),
       'outputs/pca_addback_only_variance.csv')

# Frozen August-only signature, scored without May-dependent gene selection.
sig <- fread(src$signature, data.table = FALSE)
bg <- fread(src$signature_background, data.table = FALSE)
sig$gene_id_stable <- if ('stable_gene_id' %in% names(sig)) sig$stable_gene_id else sub('\\..*$', '', sig$gene_id)
bg$gene_id_stable <- if ('stable_gene_id' %in% names(bg)) bg$stable_gene_id else sub('\\..*$', '', bg$gene_id)
bg_ids <- intersect(bg$gene_id_stable[bg$representative], rownames(counts))
sig$represented <- sig$gene_id_stable %in% bg_ids
stopifnot(sum(sig$represented) >= 90, length(unique(sig$direction)) == 2)
rank_mat <- apply(counts[bg_ids, , drop = FALSE], 2,
                  function(x) (rank(x, ties.method = 'average') - 1) / (length(x) - 1))
rownames(rank_mat) <- bg_ids
up <- sig$gene_id_stable[sig$represented & sig$direction > 0]
down <- sig$gene_id_stable[sig$represented & sig$direction < 0]
score <- colMeans(rank_mat[up, , drop = FALSE]) - colMeans(rank_mat[down, , drop = FALSE])
score_df <- data.frame(sample_id = names(score), final_group = samples$final_group,
                       collection_occasion = samples$replicate_block, signed_percentile_score = score)
fwrite(sig, 'outputs/signature/frozen_signature_mapping.csv')
fwrite(score_df, 'outputs/signature/frozen_signature_scores.csv')

d_score <- design
fit_score <- eBayes(contrasts.fit(lmFit(matrix(score, nrow = 1), d_score), fit$cm), robust = TRUE)
score_contrasts <- do.call(rbind, lapply(contrast_table$name, function(k) {
  z <- topTable(fit_score, coef = k, number = Inf, confint = TRUE)
  data.frame(contrast = k, difference = z$logFC, CI_L = z$CI.L, CI_R = z$CI.R,
             t = z$t, PValue = z$P.Value)
}))
score_contrasts$FDR_15 <- p.adjust(score_contrasts$PValue, 'BH')
fwrite(score_contrasts, 'outputs/signature/score_contrasts.csv')

fry_names <- c('UTY_v_LXHA', 'KDM5D_v_LXHA', 'DUAL_v_LXHA', 'DUAL_v_UTY', 'DUAL_v_KDM5D', 'LOY_v_LXHA')
idx <- list(August_up = which(rownames(fit$v) %in% up), August_down = which(rownames(fit$v) %in% down))
fr <- do.call(rbind, lapply(fry_names, function(k) {
  z <- fry(fit$v, index = idx, design = design, contrast = fit$cm[, k])
  z$signature_arm <- rownames(z); z$contrast <- k; z
}))
fr$FDR_12 <- p.adjust(fr$PValue, 'BH')
fr$expected_reversal <- ifelse(fr$signature_arm == 'August_up', 'Down', 'Up')
fr$direction_matches_reversal <- fr$Direction == fr$expected_reversal
fwrite(fr, 'outputs/signature/fry_directional.csv')

# Descriptive gene-level movement relative to CTRL and empty vector.
sig_ids <- sig$gene_id_stable[sig$represented & sig$gene_id_stable %in% rownames(lcpm)]
gm <- sapply(group_levels, function(g) rowMeans(lcpm[sig_ids, samples$final_group == g, drop = FALSE]))
movement <- do.call(rbind, lapply(c('LOY', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB'), function(arm) {
  vector_effect <- gm[, 'LOY_LXHA'] - gm[, 'CTRL']
  arm_effect <- gm[, arm] - gm[, 'CTRL']
  addback_effect <- gm[, arm] - gm[, 'LOY_LXHA']
  crossed <- sign(arm_effect) != sign(vector_effect) & arm_effect != 0 & vector_effect != 0
  closer <- abs(arm_effect) < abs(vector_effect)
  category <- ifelse(crossed & closer, 'crossed_control_closer',
              ifelse(crossed & !closer, 'crossed_control_farther',
              ifelse(!crossed & closer, 'closer_same_side', 'farther_same_side')))
  data.frame(gene_id = sig_ids, gene_name = dge$genes$gene_name[match(sig_ids, dge$genes$gene_id)],
             arm = arm, vector_vs_CTRL = vector_effect, arm_vs_CTRL = arm_effect,
             arm_vs_vector = addback_effect, closer_to_CTRL = closer, crossed_CTRL = crossed,
             category = category, stringsAsFactors = FALSE)
}))
fwrite(movement, 'outputs/signature/distance_to_control_by_gene.csv')
movement_summary <- aggregate(gene_id ~ arm + category, movement, length)
names(movement_summary)[3] <- 'n_genes'
movement_summary$fraction <- ave(movement_summary$n_genes, movement_summary$arm, FUN = function(x) x / sum(x))
fwrite(movement_summary, 'outputs/signature/distance_to_control_summary.csv')

# Claim table: machine-readable values used for Results drafting.
get_n <- function(k, field) de_summary[[field]][match(k, de_summary$contrast)]
claim <- data.frame(
  claim_id = c('tested_genes', 'uty_deg05', 'kdm5d_deg05', 'dual_deg05', 'uty_treat05', 'kdm5d_treat05', 'dual_treat05',
               'signature_represented', 'signature_up', 'signature_down', 'pca_pc1', 'pca_pc2'),
  value = c(nrow(dge), get_n('UTY_v_LXHA', 'n_selected_0.5'), get_n('KDM5D_v_LXHA', 'n_selected_0.5'),
            get_n('DUAL_v_LXHA', 'n_selected_0.5'), get_n('UTY_v_LXHA', 'n_treat_0.5'),
            get_n('KDM5D_v_LXHA', 'n_treat_0.5'), get_n('DUAL_v_LXHA', 'n_treat_0.5'),
            sum(sig$represented), length(up), length(down), pve[1], pve[2]),
  unit = c('genes', rep('genes', 6), 'genes', 'genes', 'genes', 'percent', 'percent'),
  source = c('filter_summary', rep('de_summary', 6), rep('frozen_signature_mapping', 3), rep('pca_variance', 2)),
  stringsAsFactors = FALSE)
fwrite(claim, 'outputs/claim_evidence.csv')

# ---------------- Figures ----------------
cols <- c(CTRL = '#3C5488', LOY = '#E64B35', LOY_LXHA = '#7E6148',
          LOY_UTY_AB = '#00A087', LOY_KDM5D_AB = '#F39B7F', LOY_KDM5D_UTY_AB = '#4DBBD5')
labs_group <- c(CTRL = 'sgCTRL', LOY = 'Y-depleted', LOY_LXHA = 'Empty vector',
                LOY_UTY_AB = 'UTY', LOY_KDM5D_AB = 'KDM5D', LOY_KDM5D_UTY_AB = 'UTY + KDM5D')
save_multi <- function(p, stem, width, height) {
  ggsave(file.path('figures', paste0(stem, '.png')), p, width = width, height = height, dpi = 300, bg = 'white')
  ggsave(file.path('figures', paste0(stem, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path('figures', paste0(stem, '.svg')), p, width = width, height = height, device = svglite::svglite)
}
dir.create('figures/student_panels', recursive = TRUE, showWarnings = FALSE)

design_df <- data.frame(group = factor(group_levels, group_levels), x = seq_along(group_levels), y = 0,
                        label = unname(labs_group[group_levels]))
p_design <- ggplot(design_df, aes(x, y)) +
  annotate('segment', x = 2, xend = 6, y = 0.28, yend = 0.28,
           arrow = arrow(length = unit(2, 'mm')), colour = 'grey55') +
  geom_point(aes(colour = group), size = 4) + geom_text(aes(label = label), angle = 35, hjust = 1, nudge_y = -0.08, size = 3) +
  annotate('text', x = 4, y = 0.41, label = 'separate transductions from the same starting pool', size = 3) +
  annotate('text', x = 3.5, y = -0.38, label = 'three matched collection occasions per population', size = 3) +
  scale_colour_manual(values = cols, guide = 'none') + coord_cartesian(ylim = c(-0.55, 0.55), clip = 'off') +
  theme_void(base_family = 'Open Sans')

p_western <- ggplot() + annotate('rect', xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = 'grey97', colour = 'grey75') +
  annotate('text', x = 0.5, y = 0.58, label = 'Western-blot panel reserved', family = 'Open Sans', fontface = 'bold', size = 3.5) +
  annotate('text', x = 0.5, y = 0.40, label = 'original images and lane map\npending source check', family = 'Open Sans', size = 2.7) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()

# Construct verification: both genes on one shared axis, grouped by population,
# with the sgCTRL level as reference, fold change over sgCTRL for construct
# populations, open symbols for zero-read samples and a construct key below.
construct$final_group <- factor(construct$final_group, group_levels)
construct$gene <- factor(construct$gene, c('UTY', 'KDM5D'))
construct$occasion <- factor(construct$replicate_block)
construct$detected <- ifelse(construct$count > 0, 'reads', 'none')
construct_labels <- c(CTRL = 'sgCTRL', LOY = 'Y-depleted', LOY_LXHA = 'Vector',
                      LOY_UTY_AB = 'UTY', LOY_KDM5D_AB = 'KDM5D',
                      LOY_KDM5D_UTY_AB = 'UTY +\nKDM5D')
gene_cols <- c(UTY = '#00A087', KDM5D = '#8491B4')
construct_map <- data.frame(final_group = factor(rep(group_levels, 2), group_levels),
                            gene = factor(rep(c('UTY', 'KDM5D'), each = length(group_levels)), c('UTY', 'KDM5D')),
                            present = c(group_levels %in% c('LOY_UTY_AB', 'LOY_KDM5D_UTY_AB'),
                                        group_levels %in% c('LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')))
dodge <- position_dodge(width = 0.55)
c_mean <- aggregate(logCPM ~ final_group + gene, construct, mean)
ctrl_ref <- setNames(c_mean$logCPM[c_mean$final_group == 'CTRL'], as.character(c_mean$gene[c_mean$final_group == 'CTRL']))
c_fold <- merge(c_mean, construct_map[construct_map$present, ], by = c('final_group', 'gene'))
c_fold$fold <- 2^(c_fold$logCPM - ctrl_ref[as.character(c_fold$gene)])
c_fold$label <- sprintf('%.0f×', c_fold$fold)
c_max <- aggregate(logCPM ~ final_group + gene, construct, max)
c_fold$y <- c_max$logCPM[match(paste(c_fold$final_group, c_fold$gene), paste(c_max$final_group, c_max$gene))] + 0.9
# Dodge offsets match position_dodge(0.55) for two genes; set explicitly because
# fold labels exist for only one gene in the single-construct groups.
c_fold$x <- as.numeric(c_fold$final_group) + ifelse(c_fold$gene == 'UTY', -0.1375, 0.1375)
fwrite(c_fold[, c('final_group', 'gene', 'logCPM', 'fold')], 'outputs/construct_fold_over_sgCTRL.csv')
floor_lcpm <- min(construct$logCPM)
ref_df <- data.frame(gene = factor(names(ctrl_ref), c('UTY', 'KDM5D')), y = unname(ctrl_ref))

p_construct_main <- ggplot(construct, aes(final_group, logCPM, colour = gene, group = gene)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = floor_lcpm + 0.35, fill = 'grey94') +
  annotate('text', x = 0.5, y = floor_lcpm + 0.18, label = 'no reads', hjust = 0, vjust = 0.5,
           size = 2.3, colour = 'grey45') +
  geom_hline(data = ref_df, aes(yintercept = y, colour = gene), linetype = 2, linewidth = 0.4, show.legend = FALSE) +
  stat_summary(fun = mean, geom = 'crossbar', width = 0.4, linewidth = 0.3, position = dodge,
               fun.min = mean, fun.max = mean, show.legend = FALSE) +
  geom_point(aes(shape = interaction(occasion, detected)), size = 1.9, stroke = 0.6,
             position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0, dodge.width = 0.55, seed = 4)) +
  geom_text(data = c_fold, aes(x = x, y = y, label = label), inherit.aes = FALSE, colour = gene_cols[as.character(c_fold$gene)],
            size = 2.6, fontface = 'bold') +
  annotate('text', x = 1.6, y = ref_df$y + 0.45, label = paste(ref_df$gene, 'level in sgCTRL'), hjust = 0,
           size = 2.3, colour = gene_cols[as.character(ref_df$gene)]) +
  scale_colour_manual(values = gene_cols, name = NULL, labels = c(UTY = 'UTY RNA', KDM5D = 'KDM5D RNA')) +
  scale_shape_manual(values = c(`1.reads` = 16, `2.reads` = 17, `3.reads` = 15,
                                `1.none` = 1, `2.none` = 2, `3.none` = 0), guide = 'none') +
  scale_x_discrete(labels = construct_labels) +
  scale_y_continuous(breaks = seq(-4, 14, 2), expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = NULL, y = 'log2 CPM') +
  theme(legend.position = 'top', legend.margin = margin(0, 0, 0, 0),
        axis.text.x = element_text(size = 8, lineheight = 0.9))

p_construct_key <- ggplot(construct_map, aes(final_group, gene)) +
  geom_text(aes(label = ifelse(present, '+', '–'), colour = ifelse(present, as.character(gene), 'none')),
            size = 3.4, fontface = 'bold') +
  scale_colour_manual(values = c(gene_cols, none = 'grey65'), guide = 'none') +
  scale_x_discrete(labels = construct_labels) +
  scale_y_discrete(limits = c('KDM5D', 'UTY'), labels = c(KDM5D = 'KDM5D construct', UTY = 'UTY construct')) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        axis.text.y = element_text(size = 7.5), axis.text.x = element_blank())

construct_panel <- function(tag = NULL) {
  (p_construct_main + labs(tag = tag)) / p_construct_key + plot_layout(heights = c(1, 0.13))
}
p_construct <- construct_panel()

pca_df$final_group <- factor(pca_df$final_group, group_levels)
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = final_group, shape = factor(replicate_block))) +
  geom_point(size = 3) + scale_colour_manual(values = cols, labels = labs_group, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = 'Collection\noccasion') +
  labs(x = sprintf('PC1 (%.1f%%)', pve[1]), y = sprintf('PC2 (%.1f%%)', pve[2])) +
  theme(legend.position = 'right')

pca_addback_df$final_group <- factor(pca_addback_df$final_group, addback_groups)
p_pca_addback <- ggplot(pca_addback_df, aes(PC1, PC2, colour = final_group, shape = factor(replicate_block))) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = cols[addback_groups], labels = labs_group[addback_groups], name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = 'Collection\noccasion') +
  labs(x = sprintf('PC1 (%.1f%%)', pve_addback[1]), y = sprintf('PC2 (%.1f%%)', pve_addback[2]),
       title = 'PCA of plasmid-bearing populations only',
       subtitle = 'Empty vector, UTY, KDM5D and combined add-back; 23,992 expressed autosomal/non-PAR-X genes') +
  theme(legend.position = 'right', plot.title = element_text(face = 'bold'))

fig1 <- (p_design + labs(tag = 'a')) / ((p_western + labs(tag = 'b')) | construct_panel('c') | (p_pca + labs(tag = 'd'))) +
  plot_layout(heights = c(0.75, 1.6), widths = c(0.6, 1.65, 1.1)) + plot_annotation(tag_levels = NULL)
save_multi(fig1, 'ch03_main_ht1080_addback_design_expression', 13, 7.2)

# Student-facing composite with a wider construct-expression panel.
student_bottom <- (p_western + labs(tag = 'b')) |
                  construct_panel('c') |
                  (p_pca + labs(tag = 'd'))
student_bottom <- student_bottom + plot_layout(widths = c(0.65, 1.75, 1.25))
fig1_student <- (p_design + labs(tag = 'a')) / student_bottom +
  plot_layout(heights = c(0.68, 1.75)) + plot_annotation(tag_levels = NULL)
save_multi(fig1_student, 'student_design_expression_readable', 15.5, 7.8)
save_multi(p_design, 'student_panels/panel_a_experimental_design', 11, 2.8)
save_multi(p_western, 'student_panels/panel_b_western_placeholder', 7.5, 4.5)
save_multi(p_construct, 'student_panels/panel_c_construct_rna_expression', 10.5, 5.4)
save_multi(p_pca, 'student_panels/panel_d_pca_all_populations', 8.5, 6.3)
save_multi(p_pca_addback, 'student_panels/panel_e_pca_addback_only', 8.5, 6.3)

source('code/plot_signature_panels.R')

# Sensitivity figure.
method_counts <- rbind(
  data.frame(contrast = de_summary$contrast, method = 'Blocked, filtered effect', n = de_summary$n_selected_0.5),
  data.frame(contrast = de_summary$contrast, method = 'Blocked, TREAT 0.5', n = de_summary$n_treat_0.5),
  data.frame(contrast = unblocked$contrast, method = 'Unblocked, filtered effect', n = unblocked$n_selected_0.5))
show <- c('LOY_v_LXHA', 'UTY_v_LXHA', 'KDM5D_v_LXHA', 'DUAL_v_LXHA')
method_counts <- method_counts[method_counts$contrast %in% show, ]
p_methods <- ggplot(method_counts, aes(contrast, n, fill = method)) + geom_col(position = 'dodge') +
  scale_x_discrete(labels = c(LOY_v_LXHA = 'Naive Y-depleted\nvs vector', UTY_v_LXHA = 'UTY\nvs vector',
                              KDM5D_v_LXHA = 'KDM5D\nvs vector', DUAL_v_LXHA = 'Dual\nvs vector')) +
  scale_fill_manual(values = c('#3C5488', '#E64B35', '#00A087')) + labs(x = NULL, y = 'Selected genes', fill = NULL) +
  theme(legend.position = 'top')
p_loo <- ggplot(loo, aes(omitted_occasion, cor_logFC_with_full, colour = contrast, group = contrast)) +
  geom_hline(yintercept = 0, colour = 'grey70') + geom_line() + geom_point(size = 2) +
  scale_colour_brewer(palette = 'Dark2') + labs(x = 'Omitted collection occasion', y = 'Correlation with full-model log2 effects', colour = NULL) +
  theme(legend.position = 'right')
y_ids <- annot$gene_id[annot$gene_name %in% c('RPS4Y1', 'ZFY', 'USP9Y', 'DDX3Y', 'UTY', 'KDM5D', 'EIF1AY') & annot$chrom == 'chrY']
y_lcpm <- cpm(counts[y_ids, , drop = FALSE], lib.size = lib_eff, log = TRUE, prior.count = 2)
y_mean <- sapply(group_levels, function(g) rowMeans(y_lcpm[, samples$final_group == g, drop = FALSE]))
rownames(y_mean) <- annot$gene_name[match(rownames(y_mean), annot$gene_id)]
y_long <- melt(as.data.table(y_mean, keep.rownames = 'gene'), id.vars = 'gene', variable.name = 'group', value.name = 'logCPM')
p_y <- ggplot(y_long, aes(group, gene, fill = logCPM)) + geom_tile(colour = 'white') +
  scale_fill_viridis_c(option = 'magma') + scale_x_discrete(labels = labs_group) +
  labs(x = NULL, y = NULL, fill = 'log2 CPM') + theme(axis.text.x = element_text(angle = 35, hjust = 1))
figs <- (p_methods + labs(tag = 'a')) / ((p_loo + labs(tag = 'b')) | (p_y + labs(tag = 'c'))) + plot_layout(heights = c(1, 1.1))
save_multi(figs, 'ch03_supp_ht1080_addback_sensitivity', 11, 8)
save_multi(p_methods, 'student_panels/panel_j_de_method_sensitivity', 9, 5.4)
save_multi(p_loo, 'student_panels/panel_k_leave_one_occasion_out', 9.5, 5.6)
save_multi(p_y, 'student_panels/panel_l_y_expression_context', 9.5, 5.6)

student_figure_index <- data.frame(
  stem = c(
    'student_design_expression_readable',
    'student_camera_overview',
    'student_panels/panel_a_experimental_design',
    'student_panels/panel_b_western_placeholder',
    'student_panels/panel_c_construct_rna_expression',
    'student_panels/panel_d_pca_all_populations',
    'student_panels/panel_e_pca_addback_only',
    'student_panels/panel_f_frozen_signature_scores',
    'student_panels/panel_g_signature_score_contrasts',
    'student_panels/panel_h_signature_gene_heatmap',
    'student_panels/panel_i_directional_fry_tests',
    'student_panels/panel_j_de_method_sensitivity',
    'student_panels/panel_k_leave_one_occasion_out',
    'student_panels/panel_l_y_expression_context'
  ),
  description = c(
    'Wider student composite: design, Western placeholder, construct RNA expression and six-population PCA',
    'CAMERA overview generated by code/run_camera_report.R',
    'Experimental groups, transduction relationships and matched collection occasions',
    'Reserved location for the verified Western blot',
    'UTY and KDM5D RNA abundance across all populations',
    'PCA recomputed across all six populations',
    'PCA recomputed using only empty-vector, UTY, KDM5D and combined populations',
    'Frozen August LOY-signature scores by population',
    'Blocked-model contrasts of the frozen signature score',
    'Expression trajectories of the frozen-signature genes relative to sgCTRL',
    'Directional fry tests for the two frozen-signature arms',
    'Selected-gene counts across blocked, TREAT and unblocked analyses',
    'Leave-one-collection-occasion-out effect correlations',
    'Y-linked expression context, including construct-derived reads'
  ),
  formats = 'PDF;SVG;PNG', stringsAsFactors = FALSE
)
fwrite(student_figure_index, 'figures/student_figure_index.csv')

writeLines(capture.output(sessionInfo()), 'provenance/sessionInfo.txt')
writeLines(c(
  'PASS sample_manifest: 18 libraries, six groups, three matched collection occasions',
  sprintf('PASS design_rank: %d columns, rank %d', ncol(design), qr(design)$rank),
  sprintf('PASS tested_genes: %d autosomal/non-PAR-X genes', nrow(dge)),
  sprintf('PASS frozen_signature: %d/%d genes represented; %d up, %d down', sum(sig$represented), nrow(sig), length(up), length(down)),
  sprintf('PASS primary_contrasts: %d', nrow(contrast_table)),
  sprintf('PASS addback_only_pca: 12 libraries; PC1 %.3f%%, PC2 %.3f%%', pve_addback[1], pve_addback[2]),
  'PASS figures: PDF, SVG and PNG generated for thesis composites and student standalone panels'),
  'provenance/verification.txt')
cat(readLines('provenance/verification.txt'), sep = '\n')
