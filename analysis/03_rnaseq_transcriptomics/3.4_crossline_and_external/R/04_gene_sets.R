# Gene-set concordance: classify every commonly tested set, assign themes, screen position.
is_ribo <- function(g) grepl('^(RPL[0-9]|RPS[0-9]|RPLP[0-9]|RPSA|MRPL[0-9]|MRPS[0-9]|FAU|UBA52)', g)
is_hist <- function(g) grepl('^H(1|2A|2B|3|4)C[0-9]|^H1-[0-9]|^H2A[CZ]|^H3C|^H4C', g)

gene_sets <- function(x, cfg) {
  msig <- readRDS(file.path(REPO, cfg$msigdb_rds))
  H <- msig$H; C2 <- msig$C2
  assert(all(cfg$polycomb_reference_sets %in% names(C2)), 'Polycomb reference sets missing from MSigDB C2')
  nbhd <- c(data.table::fread(file.path(REPO, cfg$grhl2_neighbourhood))$gene, 'GRHL2')
  refs <- list(
    interferon_inflammatory = unique(unlist(H[c('HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
                                               'HALLMARK_INFLAMMATORY_RESPONSE')])),
    mitotic = unique(unlist(H[c('HALLMARK_G2M_CHECKPOINT', 'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_E2F_TARGETS')])),
    GRHL2_epithelial = nbhd,
    polycomb_targets = unique(unlist(C2[cfg$polycomb_reference_sets])))
  s <- x$sets
  col <- function(L, v) x$fg[[L]][[v]][match(s, x$fg[[L]]$pathway)]
  d <- data.frame(pathway = s, collection = col('DU145', 'collection'), subcollection = col('DU145', 'subcollection'))
  for (L in LINES) {
    d[[paste0('NES_', L)]] <- col(L, 'NES'); d[[paste0('BH_', L)]] <- col(L, 'BH_across_collections')
    d[[paste0('LE_', L)]] <- col(L, 'leadingEdge')
  }
  nes <- as.matrix(d[, paste0('NES_', LINES)]); bh <- as.matrix(d[, paste0('BH_', LINES)])
  sig <- bh < cfg$alpha; sig[is.na(sig)] <- FALSE
  colnames(sig) <- colnames(nes) <- LINES
  d$n_significant <- rowSums(sig)
  d$same_sign_all <- apply(sign(nes), 1, function(v) length(unique(v)) == 1)
  d$significant_opposite <- sapply(seq_len(nrow(d)), function(i) length(unique(sign(nes[i, sig[i, ]]))) > 1)
  d$class <- ifelse(d$n_significant >= 2 & d$same_sign_all, 'shared',
             ifelse(d$significant_opposite, 'divergent',
             ifelse(d$n_significant >= 2, 'agree_where_significant_third_opposite',
             ifelse(d$n_significant == 1, 'line_specific', 'not_significant'))))
  keep <- d$class != 'not_significant'
  # Themes from the union of leading edges in the lines where the set is significant.
  le_union <- lapply(which(keep), function(i) unique(unlist(strsplit(unlist(d[i, paste0('LE_', LINES)][sig[i, ]]), ';', fixed = TRUE))))
  frac <- function(test) vapply(le_union, function(g) mean(test(g)), numeric(1))
  th <- data.frame(pathway = d$pathway[keep], frac_ribosomal = frac(is_ribo), frac_histone = frac(is_hist))
  for (nm in names(refs)) th[[paste0('frac_', nm)]] <- frac(function(g) g %in% refs[[nm]])
  bio <- as.matrix(th[, paste0('frac_', names(refs))])
  best <- names(refs)[max.col(bio, ties.method = 'first')]; best_frac <- apply(bio, 1, max)
  th$theme <- ifelse(d$collection[keep] == 'C1', 'positional',
              ifelse(th$frac_ribosomal >= cfg$theme_ribosomal_fraction, 'ribosomal/translation',
              ifelse(th$frac_histone >= cfg$theme_histone_fraction, 'histone cluster (6p22)',
              ifelse(best_frac >= cfg$theme_min_fraction, best, 'other'))))
  th$leading_edge_union_size <- lengths(le_union)
  # Post hoc sensitivity (disclosed): widen the inflammatory reference to the Hallmark cytokine-signalling sets.
  ext <- unique(unlist(H[cfg$inflammatory_extended_sensitivity$hallmark_sets]))
  th$frac_inflammatory_extended <- frac(function(g) g %in% ext)
  bio2 <- cbind(bio, frac_inflammatory_extended = th$frac_inflammatory_extended)
  bio2[, 'frac_interferon_inflammatory'] <- 0   # replaced by the extended group
  best2 <- c(names(refs), 'inflammatory_extended')[max.col(bio2, ties.method = 'first')]
  best2[best2 == 'inflammatory_extended'] <- 'interferon_inflammatory'
  th$theme_sensitivity_extended_inflammatory <- ifelse(d$collection[keep] == 'C1', 'positional',
              ifelse(th$frac_ribosomal >= cfg$theme_ribosomal_fraction, 'ribosomal/translation',
              ifelse(th$frac_histone >= cfg$theme_histone_fraction, 'histone cluster (6p22)',
              ifelse(apply(bio2, 1, max) >= cfg$theme_min_fraction, best2, 'other'))))
  # Positional screen: fraction of each significant line's leading edge on its most common arm.
  arm_of <- lapply(LINES, function(L) { e <- x$eff[[L]]; e <- e[!e$PAR & !is.na(e$arm) & e$gene_name != '', ]; e <- e[!duplicated(e$gene_name), ]; setNames(e$arm, e$gene_name) })
  names(arm_of) <- LINES
  for (L in LINES) {
    th[[paste0('top_arm_fraction_', L)]] <- vapply(which(keep), function(i) {
      if (!sig[i, L]) return(NA_real_)
      g <- strsplit(d[i, paste0('LE_', L)], ';', fixed = TRUE)[[1]]; a <- arm_of[[L]][g]; a <- a[!is.na(a)]
      if (!length(a)) NA_real_ else max(table(a)) / length(a)
    }, numeric(1))
  }
  tf_cols <- paste0('top_arm_fraction_', LINES)
  th$positionally_concentrated <- th$theme == 'positional' | apply(th[, tf_cols], 1, function(v) any(v >= cfg$positional_fraction, na.rm = TRUE))
  d <- merge(d, th, by = 'pathway', all.x = TRUE, sort = FALSE)
  # DU145 copy-number rerun carried over; HT1080 and A549 have no GENCODE-keyed copy-number table.
  du <- x$fg$DU145
  for (v in c('survives_CN_screen', 'NES_CN_stable_only', 'BH_CN_stable_only')) d[[paste0('DU145_', v)]] <- du[[v]][match(d$pathway, du$pathway)]
  # HT1080 period sensitivity.
  per <- x$fg$HT1080_periods
  for (p in unique(per$contrast)) {
    q <- per[per$contrast == p, ]
    d[[paste0('HT1080_NES_', p)]] <- q$NES[match(d$pathway, q$pathway)]
    d[[paste0('HT1080_BH_', p)]] <- q$BH_across_collections[match(d$pathway, q$pathway)]
  }
  pcols <- grep('^HT1080_NES_', names(d), value = TRUE)
  d$HT1080_same_sign_both_periods <- apply(sign(d[, c(pcols, 'NES_HT1080')]), 1, function(v) length(unique(v)) == 1)
  d <- d[order(match(d$class, c('shared', 'divergent', 'agree_where_significant_third_opposite', 'line_specific', 'not_significant')), d$pathway), ]
  tab(d, 'gene_sets_all_common')
  sd <- d[d$class %in% c('shared', 'divergent'), ]
  tab(sd, 'gene_sets_shared_divergent')
  summ <- as.data.frame(table(class = d$class, theme = ifelse(is.na(d$theme), 'none', d$theme)))
  summ <- summ[summ$Freq > 0, ]
  tab(summ, 'gene_sets_class_by_theme')
  cs <- as.data.frame(table(class = d$class)); tab(cs, 'gene_sets_class_summary')
  pair_cor <- do.call(rbind, lapply(combn(LINES, 2, simplify = FALSE), function(p) data.frame(line_1 = p[1], line_2 = p[2], sets = nrow(d),
    NES_Spearman = cor(d[[paste0('NES_', p[1])]], d[[paste0('NES_', p[2])]], method = 'spearman', use = 'complete.obs'))))
  tab(pair_cor, 'gene_sets_NES_correlations')
  list(sets = d, shared_divergent = sd, summary = summ, correlations = pair_cor, refs = refs)
}
