# Gene-level concordance on the common primary universe.
gene_level <- function(x, cfg) {
  u <- x$genes
  md <- sapply(LINES, function(L) { e <- x$eff[[L]]; e$mean_difference[match(u, e$gene_id)] })
  rownames(md) <- u
  sym <- x$eff$DU145$gene_name[match(u, x$eff$DU145$gene_id)]
  fdr <- sapply(LINES, function(L) { d <- x$de[[L]]; d$FDR[match(u, d$gene_id)] })
  thesis <- sapply(LINES, function(L) { d <- x$de[[L]]; d$thesis[match(u, d$gene_id)] })
  pairs <- combn(LINES, 2, simplify = FALSE)
  cor_tab <- do.call(rbind, lapply(pairs, function(p) data.frame(line_1 = p[1], line_2 = p[2], genes = nrow(md),
    Spearman = cor(md[, p[1]], md[, p[2]], method = 'spearman'), Pearson = cor(md[, p[1]], md[, p[2]]))))
  tab(cor_tab, 'gene_level_correlations')
  # Three-way agreement against a permutation null that shuffles each line's effects independently.
  thr <- cfg$gene_agreement_threshold
  agree <- function(m) sum(apply(abs(m) >= thr, 1, all) & apply(sign(m), 1, function(s) length(unique(s)) == 1))
  big <- function(m) sum(apply(abs(m) >= thr, 1, all))
  obs <- agree(md); obs_big <- big(md)
  set.seed(cfg$seed)
  null <- replicate(cfg$gene_permutations, { m <- apply(md, 2, sample); c(agree(m), big(m)) })
  # Two separate questions:
  # (1) magnitude: do large effects co-occur across lines more than under independent shuffling?
  # (2) direction: among genes large in all three, is sign agreement different from the 1/4 expected
  #     for independent signs? (1) is answered by the permutation, (2) by an exact binomial test.
  bt <- binom.test(obs, obs_big, p = 0.25)
  agr <- data.frame(threshold = thr, permutations = ncol(null),
                    genes_large_in_all_three = obs_big, null_large_mean = mean(null[2, ]),
                    null_large_2.5 = quantile(null[2, ], .025), null_large_97.5 = quantile(null[2, ], .975),
                    empirical_P_large_cooccurrence = (sum(null[2, ] >= obs_big) + 1) / (ncol(null) + 1),
                    same_sign_all_three = obs, null_same_sign_mean = mean(null[1, ]),
                    fraction_same_sign_given_large = obs / obs_big, expected_fraction_independent_signs = 0.25,
                    binomial_P_sign_agreement = bt$p.value,
                    binomial_CI_low = bt$conf.int[1], binomial_CI_high = bt$conf.int[2])
  tab(agr, 'gene_level_three_way_agreement')
  # Shared genes under the agreed rule: FDR < 0.05 in >= 2 lines and the same sign of effect in all three.
  sig <- fdr < cfg$alpha; sig[is.na(sig)] <- FALSE
  same <- apply(sign(md), 1, function(s) length(unique(s)) == 1)
  opp <- sapply(seq_len(nrow(md)), function(i) { s <- sign(md[i, sig[i, ]]); length(unique(s)) > 1 })
  g <- data.frame(gene_id = u, gene_name = sym, md, FDR_HT1080 = fdr[, 'HT1080'], FDR_A549 = fdr[, 'A549'], FDR_DU145 = fdr[, 'DU145'],
                  thesis_HT1080 = thesis[, 'HT1080'], thesis_DU145 = thesis[, 'DU145'],
                  n_significant = rowSums(sig), same_sign_all = same, significant_opposite = opp)
  names(g)[3:5] <- paste0('md_', LINES)
  g$class <- ifelse(g$n_significant >= 2 & g$same_sign_all, 'shared',
             ifelse(g$significant_opposite, 'divergent',
             ifelse(g$n_significant == 1, 'line_specific', 'not_significant')))
  tab(g, 'gene_level_all_common_genes'); tab(g[g$class %in% c('shared', 'divergent'), ], 'gene_level_shared_divergent')
  ov <- do.call(rbind, lapply(pairs, function(p) {
    a <- sig[, p[1]]; b <- sig[, p[2]]
    data.frame(line_1 = p[1], line_2 = p[2], sig_1 = sum(a), sig_2 = sum(b), both = sum(a & b),
               both_same_sign = sum(a & b & sign(md[, p[1]]) == sign(md[, p[2]])))
  }))
  # Sign agreement among genes significant in both lines of a pair, against 1/2.
  ov$sign_agreement_fraction <- ifelse(ov$both > 0, ov$both_same_sign / ov$both, NA)
  ov$binomial_P_vs_half <- ifelse(ov$both > 0, mapply(function(k, n) binom.test(k, n, .5)$p.value, ov$both_same_sign, pmax(ov$both, 1)), NA)
  tab(ov, 'gene_level_significant_overlap')
  tab(data.frame(class = names(table(g$class)), genes = as.integer(table(g$class))), 'gene_level_class_summary')
  focus <- g[g$gene_name %in% cfg$replication_genes, ]
  tab(focus, 'gene_level_focus_genes')
  list(genes = g, correlations = cor_tab, agreement = agr, overlap = ov, null = null)
}
