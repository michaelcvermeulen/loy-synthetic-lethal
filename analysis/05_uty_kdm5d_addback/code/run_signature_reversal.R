options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

# Percent reversal of the frozen August HT1080 LOY signature by UTY, KDM5D and
# combined add-back. Reuses the scoring rule of run_analysis.R exactly (verified
# below) and adds: score components per signature half, occasion-adjusted
# percent reversal with simulation intervals for every reference/add-back
# combination, gene-level reversal, August-to-May concordance and the figures
# used by report/signature_reversal_report.Rmd. Run after run_analysis.R.

script_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)[1]
script <- normalizePath(sub('^--file=', '', script_arg), winslash = '/', mustWork = TRUE)
ROOT <- dirname(dirname(script))
ADDBACK <- dirname(ROOT)
BASELINE <- 'A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE'
STYLE <- 'A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/skills/thesis-figure-style'
setwd(ROOT)
set.seed(42)
out <- 'outputs/signature_reversal'
figdir <- 'report/figures/signature_reversal'
invisible(lapply(c(out, figdir, 'figures'), dir.create, recursive = TRUE, showWarnings = FALSE))

source(file.path(STYLE, 'scripts/thesis_style.R'))
register_thesis_fonts(STYLE)
theme_set(theme_thesis(base_size = 9))
update_geom_defaults('text', list(family = 'Open Sans', colour = '#222222'))

src <- list(
  counts = file.path(ADDBACK, 'SEPT_2026_REANALYSIS/outputs/03_counts/counts_full.rds'),
  signature = file.path(BASELINE, 'outputs/tables/August_frozen_signature_genes.csv'),
  background = file.path(BASELINE, 'outputs/tables/August_signature_discovery_background.csv'),
  august_scores = file.path(BASELINE, 'outputs/tables/August_signature_sample_scores.csv'),
  baseline_code = file.path(BASELINE, 'R/04_dosage.R'),
  baseline_config = file.path(BASELINE, 'config/analysis.yml'),
  baseline_manifest = file.path(BASELINE, 'config/sample_manifest_draft.csv'),
  samples = 'outputs/sample_manifest.csv',
  dge = 'outputs/dge_primary.rds',
  stored_scores = 'outputs/signature/frozen_signature_scores.csv',
  stored_contrasts = 'outputs/signature/score_contrasts.csv',
  de_loy_ctrl = 'outputs/de/de_LOY_v_CTRL.csv',
  de_lxha_ctrl = 'outputs/de/de_LXHA_v_CTRL.csv'
)
stopifnot(all(vapply(src, file.exists, logical(1))))

group_levels <- c('CTRL', 'LOY', 'LOY_LXHA', 'LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB')
labs_group <- c(CTRL = 'sgCTRL', LOY = 'Y-depleted', LOY_LXHA = 'Empty vector',
                LOY_UTY_AB = 'UTY', LOY_KDM5D_AB = 'KDM5D', LOY_KDM5D_UTY_AB = 'UTY + KDM5D')
cols <- c(CTRL = '#3C5488', LOY = '#E64B35', LOY_LXHA = '#7E6148',
          LOY_UTY_AB = '#00A087', LOY_KDM5D_AB = '#F39B7F', LOY_KDM5D_UTY_AB = '#4DBBD5')

# ---------------- 1. August discovery summary ----------------
sig <- fread(src$signature, data.table = FALSE)
bg <- fread(src$background, data.table = FALSE)
sig$gene_id_stable <- sub('\\..*$', '', sig$gene_id)
bg$gene_id_stable <- sub('\\..*$', '', bg$gene_id)
bg$eligible <- bg$all_three_same_direction & abs(bg$AUG_difference) >= 0.5
bg$selected <- bg$gene_id %in% sig$gene_id
discovery <- data.frame(
  step = c('Representative expressed genes in the August background',
           'Same direction in all three passage-matched pairs',
           'Same direction and |mean LOY - sgCTRL| >= 0.5 log2',
           'Eligible, higher in LOY', 'Eligible, lower in LOY',
           'Selected: 50 largest higher in LOY', 'Selected: 50 largest lower in LOY'),
  n = c(nrow(bg), sum(bg$all_three_same_direction), sum(bg$eligible),
        sum(bg$eligible & bg$AUG_difference > 0), sum(bg$eligible & bg$AUG_difference < 0),
        sum(sig$direction > 0), sum(sig$direction < 0)))
fwrite(discovery, file.path(out, 'august_discovery_funnel.csv'))
sig_out <- sig[order(-sig$direction, -abs(sig$AUG_difference)),
               c('gene_name', 'gene_id_stable', 'chrom', 'arm', 'gene_type', 'mean_logCPM', 'AUG_difference', 'direction')]
fwrite(sig_out, file.path(out, 'august_signature_genes.csv'))
fwrite(as.data.frame(table(chrom = sig$chrom, direction = ifelse(sig$direction > 0, 'up', 'down'))),
       file.path(out, 'august_signature_chromosomes.csv'))

# ---------------- 2. Rescore May libraries and verify ----------------
samples <- fread(src$samples, data.table = FALSE)
counts <- readRDS(src$counts)
stopifnot(identical(colnames(counts), samples$sample_id))
bg_ids <- intersect(bg$gene_id_stable[bg$representative], rownames(counts))
rank_mat <- apply(counts[bg_ids, , drop = FALSE], 2,
                  function(x) (rank(x, ties.method = 'average') - 1) / (length(x) - 1))
rownames(rank_mat) <- bg_ids
up <- sig$gene_id_stable[sig$direction > 0 & sig$gene_id_stable %in% bg_ids]
down <- sig$gene_id_stable[sig$direction < 0 & sig$gene_id_stable %in% bg_ids]
stopifnot(length(up) == 50, length(down) == 50)
comp <- data.frame(sample_id = samples$sample_id, group = samples$final_group,
                   occasion = samples$replicate_block,
                   up_percentile = colMeans(rank_mat[up, ]), down_percentile = colMeans(rank_mat[down, ]))
comp$score <- comp$up_percentile - comp$down_percentile
stored <- fread(src$stored_scores, data.table = FALSE)
score_check <- max(abs(comp$score - stored$signed_percentile_score[match(comp$sample_id, stored$sample_id)]))
stopifnot(score_check < 1e-12)
fwrite(comp, file.path(out, 'sample_score_components.csv'))

# Worked example: every signature gene's rank in one sgCTRL and one empty-vector library.
example_ids <- c('HT1080-sgCTRL-1', 'HT1080-LOY-LXHA-1')
worked <- rbindlist(lapply(example_ids, function(s) data.table(
  sample_id = s, gene_id = c(up, down), gene_name = sig$gene_name[match(c(up, down), sig$gene_id_stable)],
  half = rep(c('August-up', 'August-down'), c(length(up), length(down))),
  raw_count = counts[c(up, down), s], percentile = rank_mat[c(up, down), s])))
fwrite(worked, file.path(out, 'worked_example_ranks.csv'))

# ---------------- 3. Occasion-adjusted group means and percent reversal ----------------
comp$group <- factor(comp$group, group_levels)
comp$occasion <- factor(comp$occasion)
X <- model.matrix(~ 0 + group + occasion, comp)
colnames(X) <- sub('^group', '', colnames(X))
fit_outcome <- function(y) {
  f <- lm.fit(X, y)
  df <- nrow(X) - ncol(X)
  s2 <- sum(f$residuals^2) / df
  V <- s2 * chol2inv(qr.R(f$qr))[order(f$qr$pivot), order(f$qr$pivot)]
  dimnames(V) <- list(colnames(X), colnames(X))
  list(b = setNames(f$coefficients, colnames(X)), V = V, df = df)
}
outcomes <- c(score = 'Whole signature score', up_percentile = 'August-up half (mean percentile)',
              down_percentile = 'August-down half (mean percentile)')
fits <- lapply(names(outcomes), function(k) fit_outcome(comp[[k]]))
names(fits) <- names(outcomes)

# Multivariate-t draws of the group coefficients (residual df = 10).
B <- 200000
draw_coef <- function(fx) {
  L <- chol(fx$V[group_levels, group_levels])
  z <- matrix(rnorm(B * length(group_levels)), B) %*% L
  w <- sqrt(rchisq(B, fx$df) / fx$df)
  sweep(z / w, 2, fx$b[group_levels], '+')
}
draws <- lapply(fits, draw_coef)

# Percent reversal of the reference-to-sgCTRL gap achieved by the comparison
# population: 100 * (reference - comparison) / (reference - sgCTRL).
combos <- rbind(
  data.frame(reference = 'LOY_LXHA', comparison = c('LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB', 'LOY'),
             role = c('Add-back vs matched vector', 'Add-back vs matched vector', 'Add-back vs matched vector',
                      'Construct-free negative control')),
  data.frame(reference = 'LOY', comparison = c('LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB', 'LOY_LXHA'),
             role = c('Add-back vs untransduced LOY', 'Add-back vs untransduced LOY', 'Add-back vs untransduced LOY',
                      'Vector alone vs untransduced LOY')),
  data.frame(reference = 'LOY_UTY_AB', comparison = 'LOY_KDM5D_UTY_AB', role = 'Adding KDM5D to UTY'),
  data.frame(reference = 'LOY_KDM5D_AB', comparison = 'LOY_KDM5D_UTY_AB', role = 'Adding UTY to KDM5D'))
rev_table <- rbindlist(lapply(names(outcomes), function(k) {
  fx <- fits[[k]]; dr <- draws[[k]]
  rbindlist(lapply(seq_len(nrow(combos)), function(i) {
    r <- combos$reference[i]; cmp <- combos$comparison[i]
    gap <- fx$b[r] - fx$b['CTRL']
    moved <- fx$b[r] - fx$b[cmp]
    sim <- 100 * (dr[, r] - dr[, cmp]) / (dr[, r] - dr[, 'CTRL'])
    cvec <- setNames(rep(0, ncol(X)), colnames(X)); cvec[r] <- 1; cvec[cmp] <- -1
    se <- sqrt(drop(t(cvec) %*% fx$V %*% cvec))
    data.table(outcome = k, outcome_label = outcomes[[k]], reference = r, comparison = cmp, role = combos$role[i],
               reference_minus_sgCTRL = unname(gap), reference_minus_comparison = unname(moved),
               percent_reversal = unname(100 * moved / gap),
               CI_low = unname(quantile(sim, 0.025)), CI_high = unname(quantile(sim, 0.975)),
               P_difference = 2 * pt(-abs(unname(moved) / se), fx$df),
               denominator_sign_stable = mean(sign(dr[, r] - dr[, 'CTRL']) == sign(gap)))
  }))
}))
fwrite(rev_table, file.path(out, 'percent_reversal_all_combinations.csv'))
group_means <- rbindlist(lapply(names(outcomes), function(k)
  data.table(outcome = k, group = group_levels, adjusted_mean = fits[[k]]$b[group_levels] +
               mean(c(0, fits[[k]]$b[c('occasion2', 'occasion3')])),
             raw_mean = tapply(comp[[k]], comp$group, mean)[group_levels])))
fwrite(group_means, file.path(out, 'group_means.csv'))

# Consistency with the primary moderated score contrasts.
stored_c <- fread(src$stored_contrasts, data.table = FALSE)
chk <- stored_c$difference[stored_c$contrast == 'UTY_v_LXHA']
stopifnot(abs(chk - (fits$score$b['LOY_UTY_AB'] - fits$score$b['LOY_LXHA'])) < 1e-10)

# ---------------- 4. Gene-level reversal ----------------
dge <- readRDS(src$dge)
lcpm <- cpm(dge, log = TRUE, prior.count = 2)
genes_in <- intersect(c(up, down), rownames(lcpm))
gm <- sapply(group_levels, function(g) rowMeans(lcpm[genes_in, samples$final_group == g, drop = FALSE]))
gene_tab <- rbindlist(lapply(c('LOY_LXHA', 'LOY'), function(r) {
  rbindlist(lapply(setdiff(c('LOY_UTY_AB', 'LOY_KDM5D_AB', 'LOY_KDM5D_UTY_AB', 'LOY', 'LOY_LXHA'), r), function(cmp) {
    dirn <- sig$direction[match(genes_in, sig$gene_id_stable)]
    gap <- gm[, r] - gm[, 'CTRL']
    data.table(reference = r, comparison = cmp, gene_id = genes_in,
               gene_name = sig$gene_name[match(genes_in, sig$gene_id_stable)],
               August_direction = ifelse(dirn > 0, 'up', 'down'),
               reference_minus_sgCTRL = gap, reference_minus_comparison = gm[, r] - gm[, cmp],
               LOY_effect_retained = sign(gap) == dirn & abs(gap) >= 0.5,
               percent_reversal = 100 * (gm[, r] - gm[, cmp]) / gap)
  }))
}))
gene_tab[LOY_effect_retained == FALSE, percent_reversal := NA_real_]
fwrite(gene_tab, file.path(out, 'gene_level_reversal.csv'))
gene_summary <- gene_tab[LOY_effect_retained == TRUE, .(
  genes = .N, median_percent = median(percent_reversal),
  q25 = quantile(percent_reversal, 0.25), q75 = quantile(percent_reversal, 0.75),
  moved_towards_sgCTRL = sum(percent_reversal > 0), overshoot_beyond_sgCTRL = sum(percent_reversal > 100),
  moved_away = sum(percent_reversal < 0)), by = .(reference, comparison)]
fwrite(gene_summary, file.path(out, 'gene_level_reversal_summary.csv'))

# ---------------- 5. August-to-May concordance ----------------
de_lc <- fread(src$de_loy_ctrl); de_vc <- fread(src$de_lxha_ctrl)
conc <- data.table(gene_id = genes_in, gene_name = sig$gene_name[match(genes_in, sig$gene_id_stable)],
                   August_LOY_minus_sgCTRL = sig$AUG_difference[match(genes_in, sig$gene_id_stable)],
                   May_LOY_minus_sgCTRL = de_lc$logFC[match(genes_in, de_lc$gene_id)],
                   May_vector_minus_sgCTRL = de_vc$logFC[match(genes_in, de_vc$gene_id)])
fwrite(conc, file.path(out, 'august_may_concordance.csv'))
conc_stats <- data.table(
  comparison = c('May untransduced LOY - sgCTRL', 'May empty vector - sgCTRL'),
  genes = nrow(conc),
  same_direction = c(sum(sign(conc$May_LOY_minus_sgCTRL) == sign(conc$August_LOY_minus_sgCTRL)),
                     sum(sign(conc$May_vector_minus_sgCTRL) == sign(conc$August_LOY_minus_sgCTRL))),
  pearson_r = c(cor(conc$August_LOY_minus_sgCTRL, conc$May_LOY_minus_sgCTRL),
                cor(conc$August_LOY_minus_sgCTRL, conc$May_vector_minus_sgCTRL)),
  median_ratio_to_August = c(median(conc$May_LOY_minus_sgCTRL / conc$August_LOY_minus_sgCTRL),
                             median(conc$May_vector_minus_sgCTRL / conc$August_LOY_minus_sgCTRL)))
fwrite(conc_stats, file.path(out, 'august_may_concordance_summary.csv'))

# ---------------- Figures ----------------
save_multi <- function(p, path_stem, width, height) {
  ggsave(paste0(path_stem, '.png'), p, width = width, height = height, dpi = 300, bg = 'white')
  ggsave(paste0(path_stem, '.pdf'), p, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(path_stem, '.svg'), p, width = width, height = height, device = svglite::svglite)
}
fmt_pct <- function(x) sprintf('%.0f%%', x)
fmt_p <- function(p) ifelse(p < 0.001, sprintf('P = %.1e', p), sprintf('P = %.2g', p))

# Goal figure ----------------------------------------------------------
# Group means (balanced design: identical to occasion-adjusted marginal means).
sc <- tapply(comp$score, comp$group, mean)
ctrl_m <- unname(sc['CTRL']); vec_m <- unname(sc['LOY_LXHA'])
to_pct <- function(s) 100 * (s - ctrl_m) / (vec_m - ctrl_m)   # % of vector LOY state remaining
comp$group_lab <- factor(labs_group[as.character(comp$group)], rev(labs_group[group_levels]))
gm_score <- data.frame(group = group_levels, mean = tapply(comp$score, comp$group, mean)[group_levels])
gm_score$group_lab <- factor(labs_group[gm_score$group], levels(comp$group_lab))
arrows_df <- data.frame(to = c('LOY_UTY_AB', 'LOY_KDM5D_UTY_AB'))
arrows_df$x <- gm_score$mean[gm_score$group == 'LOY_LXHA']
arrows_df$xend <- gm_score$mean[match(arrows_df$to, gm_score$group)]
arrows_df$y <- as.numeric(gm_score$group_lab[gm_score$group == 'LOY_LXHA'])
arrows_df$yend <- as.numeric(gm_score$group_lab[match(arrows_df$to, gm_score$group)])
lab_rev <- rev_table[outcome == 'score' & reference == 'LOY_LXHA']
arrows_df$label <- sprintf('%s reversed\n(95%% CI %s to %s)',
  fmt_pct(lab_rev$percent_reversal[match(arrows_df$to, lab_rev$comparison)]),
  fmt_pct(lab_rev$CI_low[match(arrows_df$to, lab_rev$comparison)]),
  fmt_pct(lab_rev$CI_high[match(arrows_df$to, lab_rev$comparison)]))

comp$yn <- as.numeric(comp$group_lab)
gm_score$yn <- as.numeric(gm_score$group_lab)
p_a <- ggplot(comp, aes(score, yn)) +
  geom_vline(xintercept = c(ctrl_m, vec_m), linetype = 2, colour = 'grey70', linewidth = 0.35) +
  geom_curve(data = arrows_df, aes(x = x, xend = xend, y = y + 0.22, yend = yend + 0.22), inherit.aes = FALSE,
             curvature = 0.25, arrow = arrow(length = unit(1.6, 'mm'), type = 'closed'),
             colour = cols[arrows_df$to], linewidth = 0.5) +
  geom_segment(data = gm_score, aes(x = mean, xend = mean, y = yn - 0.22, yend = yn + 0.22, colour = group),
               linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(colour = group, shape = occasion), size = 2.2) +
  annotate('text', x = c(ctrl_m, vec_m), y = length(levels(comp$group_lab)) + 0.62,
           label = c('sgCTRL mean = 0%', 'empty-vector mean = 100%'), size = 2.3, colour = 'grey35',
           hjust = c(-0.05, 1.05)) +
  scale_y_continuous(breaks = seq_along(levels(comp$group_lab)), labels = levels(comp$group_lab),
                     expand = expansion(add = c(0.5, 0.85))) +
  geom_text(data = arrows_df, aes(x = xend - 0.012, y = yend - 0.32, label = label), inherit.aes = FALSE,
            hjust = 1, size = 2.4, lineheight = 0.9, colour = cols[arrows_df$to]) +
  scale_colour_manual(values = cols, guide = 'none') +
  scale_shape_manual(values = c(16, 17, 15), name = 'Collection\noccasion') +
  scale_x_continuous(name = 'LOY signature score (higher = more LOY-like)',
                     sec.axis = sec_axis(~ to_pct(.), name = 'LOY signature remaining, % of the vector-to-sgCTRL gap (% reversal = 100 - this value)',
                                         breaks = seq(-100, 200, 50))) +
  labs(y = NULL) + theme(legend.position = 'right')

plot_rev <- rev_table[outcome == 'score' & reference == 'LOY_LXHA']
plot_rev[, lab := labs_group[comparison]]
plot_rev[comparison == 'LOY', lab := 'Y-depleted\n(no construct)']
plot_rev[, lab := factor(lab, rev(c('UTY', 'UTY + KDM5D', 'KDM5D', 'Y-depleted\n(no construct)')))]
p_b <- ggplot(plot_rev, aes(percent_reversal, lab, colour = comparison)) +
  annotate('rect', xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, fill = 'grey96') +
  geom_vline(xintercept = c(0, 100), colour = 'grey60', linewidth = 0.35, linetype = c(1, 2)) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25, linewidth = 0.6) +
  geom_point(size = 2.6) +
  geom_text(aes(x = pmax(CI_high, 0) + 6, label = sprintf('%s\n%s', fmt_pct(percent_reversal), fmt_p(P_difference))),
            hjust = 0, size = 2.4, lineheight = 0.9, colour = '#222222') +
  annotate('text', x = c(0, 100), y = 4.55, label = c('vector\nlevel', 'sgCTRL\nlevel'), size = 2.2,
           colour = 'grey40', vjust = 0, lineheight = 0.9) +
  scale_colour_manual(values = cols, guide = 'none') +
  scale_x_continuous(limits = c(-75, 140), breaks = seq(-50, 100, 25)) +
  coord_cartesian(clip = 'off') +
  labs(x = '% reversal of the vector-to-sgCTRL gap', y = NULL)

arm_rev <- rev_table[reference == 'LOY_LXHA' & comparison %in% c('LOY_UTY_AB', 'LOY_KDM5D_UTY_AB', 'LOY_KDM5D_AB')]
arm_rev[, comp_lab := factor(labs_group[comparison], c('UTY', 'UTY + KDM5D', 'KDM5D'))]
out_levels <- c('Whole\nscore', 'August-up\ngenes\n(should fall)', 'August-down\ngenes\n(should rise)')
arm_rev[, out_lab := factor(c(score = out_levels[1], up_percentile = out_levels[2],
                              down_percentile = out_levels[3])[outcome], out_levels)]
arm_rev[, comparison := factor(comparison, c('LOY_UTY_AB', 'LOY_KDM5D_UTY_AB', 'LOY_KDM5D_AB'))]
p_c <- ggplot(arm_rev, aes(out_lab, percent_reversal, colour = comparison)) +
  geom_hline(yintercept = c(0, 100), colour = 'grey60', linewidth = 0.35, linetype = c(1, 2)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.55,
                position = position_dodge(width = 0.6)) +
  geom_point(size = 2.3, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = cols, labels = labs_group, name = NULL,
                      limits = c('LOY_UTY_AB', 'LOY_KDM5D_UTY_AB', 'LOY_KDM5D_AB')) +
  labs(x = NULL, y = '% reversal vs empty vector') + theme(legend.position = 'top')

gl <- gene_tab[reference == 'LOY_LXHA' & LOY_effect_retained == TRUE & comparison %in% c('LOY_UTY_AB', 'LOY_KDM5D_UTY_AB', 'LOY_KDM5D_AB')]
gl[, comp_lab := factor(labs_group[comparison], c('UTY', 'UTY + KDM5D', 'KDM5D'))]
gl[, shown := pmin(pmax(percent_reversal, -150), 250)]
gs <- gene_summary[reference == 'LOY_LXHA' & comparison %in% gl$comparison]
gs[, comp_lab := factor(labs_group[comparison], levels(gl$comp_lab))]
p_d <- ggplot(gl, aes(comp_lab, shown, colour = comparison)) +
  geom_hline(yintercept = c(0, 100), colour = 'grey60', linewidth = 0.35, linetype = c(1, 2)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, colour = 'grey35', fill = NA, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = 3), size = 1.1, alpha = 0.75) +
  geom_text(data = gs, aes(comp_lab, 262, label = sprintf('median %s\n%d/%d toward sgCTRL', fmt_pct(median_percent),
                                                          moved_towards_sgCTRL, genes)),
            inherit.aes = FALSE, size = 2.3, lineheight = 0.9, vjust = 0) +
  scale_colour_manual(values = cols, guide = 'none') +
  scale_y_continuous(limits = c(-150, 300), breaks = seq(-150, 250, 50)) +
  labs(x = NULL, y = 'Per-gene % reversal vs empty vector\n(display clipped to -150 to 250)')

goal <- (p_a + labs(tag = 'a')) / ((p_b + labs(tag = 'b')) | (p_c + labs(tag = 'c')) | (p_d + labs(tag = 'd'))) +
  plot_layout(heights = c(1, 1.05))
save_multi(goal, 'figures/ch03_main_ht1080_addback_signature_percent_reversal', 13, 8.6)
save_multi(goal, file.path(figdir, 'goal_figure'), 13, 8.6)
# Windows Cairo substitutes Arial; reprint the thesis PDF from SVG with embedded Open Sans.
if (system2(Sys.which('python'), c('code/export_signature_pdfs.py', 'ch03_main_ht1080_addback_signature_percent_reversal')) != 0)
  stop('Open Sans PDF export failed')
file.copy('figures/ch03_main_ht1080_addback_signature_percent_reversal.pdf', file.path(figdir, 'goal_figure.pdf'), overwrite = TRUE)

# Supporting figure 1: August discovery ---------------------------------
bg$class <- ifelse(bg$selected, ifelse(bg$AUG_difference > 0, 'Selected: up in LOY', 'Selected: down in LOY'),
                   ifelse(bg$eligible, 'Eligible, not selected', 'Not eligible'))
bg$class <- factor(bg$class, c('Not eligible', 'Eligible, not selected', 'Selected: up in LOY', 'Selected: down in LOY'))
bg <- bg[order(bg$class), ]
lab_d <- sig[order(-abs(sig$AUG_difference)), ][1:12, ]
p_disc <- ggplot(bg, aes(mean_logCPM, AUG_difference, colour = class)) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 2, colour = 'grey70', linewidth = 0.3) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_text(data = lab_d, aes(mean_logCPM, AUG_difference, label = gene_name), inherit.aes = FALSE,
            size = 2.2, vjust = -0.6, check_overlap = TRUE) +
  scale_colour_manual(values = c('Not eligible' = 'grey82', 'Eligible, not selected' = 'grey55',
                                 'Selected: up in LOY' = '#E64B35', 'Selected: down in LOY' = '#3C5488'), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  labs(x = 'Mean August log2 CPM (all six August libraries)', y = 'August LOY - sgCTRL\n(mean of 3 passage-matched differences, log2)') +
  theme(legend.position = 'top')
save_multi(p_disc, file.path(figdir, 'august_discovery'), 8, 5.2)

# Supporting figure 2: score components ----------------------------------
cl <- melt(as.data.table(comp)[, .(sample_id, group, occasion, up_percentile, down_percentile)],
           id.vars = c('sample_id', 'group', 'occasion'), variable.name = 'half', value.name = 'percentile')
cl[, half := factor(ifelse(half == 'up_percentile', 'August-up genes (50)', 'August-down genes (50)'),
                    c('August-up genes (50)', 'August-down genes (50)'))]
cl[, group_lab := factor(labs_group[as.character(group)], labs_group[group_levels])]
p_comp <- ggplot(cl, aes(group_lab, percentile, colour = group)) +
  stat_summary(fun = mean, geom = 'crossbar', width = 0.5, fun.min = mean, fun.max = mean, linewidth = 0.35) +
  geom_point(aes(shape = occasion), size = 2, position = position_jitter(width = 0.08, height = 0, seed = 2)) +
  facet_wrap(~ half, scales = 'free_y') +
  scale_colour_manual(values = cols, guide = 'none') +
  scale_shape_manual(values = c(16, 17, 15), name = 'Collection\noccasion') +
  scale_y_continuous(labels = function(x) sprintf('%.0f', 100 * x)) +
  labs(x = NULL, y = 'Mean within-sample percentile rank') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_multi(p_comp, file.path(figdir, 'score_components'), 9, 4.4)

# Supporting figure 3: worked example -----------------------------------
wk <- copy(worked)
wk[, sample_lab := ifelse(grepl('sgCTRL', sample_id), 'sgCTRL, occasion 1', 'Empty vector, occasion 1')]
wk[, half := factor(half, c('August-up', 'August-down'))]
wk_mean <- wk[, .(m = mean(percentile)), by = .(sample_lab, half)]
p_work <- ggplot(wk, aes(percentile * 100, half, colour = half)) +
  geom_point(position = position_jitter(height = 0.15, width = 0, seed = 5), size = 1.3, alpha = 0.8) +
  geom_errorbarh(data = wk_mean, aes(xmin = m * 100, xmax = m * 100, x = m * 100, y = half), inherit.aes = FALSE,
                 height = 0.6, linewidth = 0.9, colour = '#222222') +
  geom_text(data = wk_mean, aes(m * 100, half, label = sprintf('mean %.1f', m * 100)), inherit.aes = FALSE,
            vjust = -2.3, size = 2.4) +
  facet_wrap(~ sample_lab, ncol = 1) +
  scale_colour_manual(values = c('August-up' = '#E64B35', 'August-down' = '#3C5488'), guide = 'none') +
  labs(x = 'Within-sample percentile of the gene among the frozen background genes', y = NULL)
save_multi(p_work, file.path(figdir, 'worked_example'), 8, 4.2)

# Supporting figure 4: August to May concordance -------------------------
cc <- melt(conc, id.vars = c('gene_id', 'gene_name', 'August_LOY_minus_sgCTRL'),
           measure.vars = c('May_LOY_minus_sgCTRL', 'May_vector_minus_sgCTRL'), variable.name = 'may', value.name = 'may_lfc')
cc[, may := factor(ifelse(may == 'May_LOY_minus_sgCTRL', 'May untransduced Y-depleted - sgCTRL', 'May empty vector - sgCTRL'),
                   c('May untransduced Y-depleted - sgCTRL', 'May empty vector - sgCTRL'))]
cs <- copy(conc_stats); cs[, may := factor(c('May untransduced Y-depleted - sgCTRL', 'May empty vector - sgCTRL'), levels(cc$may))]
p_conc <- ggplot(cc, aes(August_LOY_minus_sgCTRL, may_lfc)) +
  geom_hline(yintercept = 0, colour = 'grey80') + geom_vline(xintercept = 0, colour = 'grey80') +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey60') +
  geom_point(aes(colour = August_LOY_minus_sgCTRL > 0), size = 1.4) +
  geom_text(data = cs, aes(x = -Inf, y = Inf, label = sprintf('r = %.2f\n%d/%d same direction\nmedian May/August = %.2f',
                                                              pearson_r, same_direction, genes, median_ratio_to_August)),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 2.4, lineheight = 0.9) +
  facet_wrap(~ may) + coord_equal() +
  scale_colour_manual(values = c(`TRUE` = '#E64B35', `FALSE` = '#3C5488'), guide = 'none') +
  labs(x = 'August LOY - sgCTRL (log2, discovery)', y = 'May difference from sgCTRL (log2)')
save_multi(p_conc, file.path(figdir, 'august_may_concordance'), 8.5, 4.6)

# ---------------- Provenance ----------------
inputs <- unlist(src)
fwrite(data.frame(role = names(inputs), path = inputs,
                  sha256 = vapply(inputs, digest::digest, character(1), algo = 'sha256', file = TRUE)),
       file.path(out, 'input_manifest.csv'))
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
v <- c(sprintf('PASS rescored May libraries reproduce frozen_signature_scores.csv (max abs diff %.1e)', score_check),
       sprintf('PASS 50 up + 50 down signature genes ranked among %d background genes', length(bg_ids)),
       'PASS unmoderated lm group differences reproduce stored UTY-vs-vector score difference',
       sprintf('PASS %d simulation draws per outcome; residual df %d', B, fits$score$df),
       sprintf('Genes passing May filter used for gene-level summaries: %d', length(genes_in)))
writeLines(v, file.path(out, 'VERIFICATION.txt'))
cat(v, sep = '\n')
