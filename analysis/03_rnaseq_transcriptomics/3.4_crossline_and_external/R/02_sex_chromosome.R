# Y inventory, PAR and X-paralog responses on one scale across the three lines.
Y_GENES <- c('RPS4Y1', 'ZFY', 'PRKY', 'USP9Y', 'DDX3Y', 'UTY', 'TMSB4Y', 'NLGN4Y', 'TXLNGY', 'KDM5D', 'EIF1AY', 'TBL1Y', 'AMELY')
XY_PAIRS <- data.frame(X = c('RPS4X', 'ZFX', 'PRKX', 'USP9X', 'DDX3X', 'KDM6A', 'TMSB4X', 'NLGN4X', 'KDM5C', 'EIF1AX'),
                       Y = c('RPS4Y1', 'ZFY', 'PRKY', 'USP9Y', 'DDX3Y', 'UTY', 'TMSB4Y', 'NLGN4Y', 'KDM5D', 'EIF1AY'))

pick_gene <- function(e, symbol, chrom, par = FALSE) {
  i <- which(e$gene_name == symbol & e$chrom == chrom & e$PAR == par)
  if (!length(i)) NA_integer_ else i[1]
}

sex_chromosome <- function(x, cfg) {
  ytab <- ptab <- xtab <- rtab <- list()
  for (L in LINES) {
    e <- x$eff[[L]]; cp <- x$cp[[L]]
    ids <- intersect(e$gene_id, rownames(cp$cpm))
    cm <- condition_means(cp$cpm, cp$meta, ids)
    ctrl <- setNames(cm$CTRL, ids); loy <- setNames(cm$LOY, ids)
    pct <- function(g) ifelse(ctrl[g] >= cfg$expressed_cpm_for_percent, 100 * loy[g] / ctrl[g], NA_real_)
    # Y inventory: male-specific genes.
    for (s in Y_GENES) {
      i <- pick_gene(e, s, 'chrY')
      g <- if (is.na(i)) NA else e$gene_id[i]
      ytab[[length(ytab) + 1]] <- data.frame(line = L, gene = s, gene_id = g,
        CTRL_CPM = if (is.na(i)) NA else ctrl[g], LOY_CPM = if (is.na(i)) NA else loy[g],
        percent_of_control = if (is.na(i)) NA else pct(g),
        mean_difference = if (is.na(i)) NA else e$mean_difference[i],
        position_Mb = if (is.na(i)) NA else e$start[i] / 1e6)
    }
    # PAR genes (counted once on chrX in all three annotations).
    par <- e[e$PAR & e$chrom == 'chrX' & e$gene_id %in% ids, ]
    par$region <- ifelse(par$start < 5e6, 'PAR1', 'PAR2')
    ptab[[L]] <- data.frame(line = L, gene = par$gene_name, gene_id = par$gene_id, region = par$region,
                            CTRL_CPM = ctrl[par$gene_id], LOY_CPM = loy[par$gene_id],
                            percent_of_control = pct(par$gene_id), mean_difference = par$mean_difference)
    # Region sums per period, then averaged: PAR1, PAR2, all male-specific Y.
    for (reg in c('PAR1', 'PAR2', 'MSY')) {
      g <- if (reg == 'MSY') e$gene_id[e$chrom == 'chrY' & !e$PAR & e$gene_id %in% ids] else par$gene_id[par$region == reg]
      ratios <- sapply(cm$per_period, function(pp) sum(pp$LOY[match(g, ids)]) / sum(pp$CTRL[match(g, ids)]))
      rtab[[length(rtab) + 1]] <- data.frame(line = L, region = reg, genes = length(g),
                                             percent_of_control = 100 * mean(ratios), periods = length(ratios))
    }
    # X paralogs and the line's overall chrX shift (median over primary non-PAR chrX genes).
    xmed <- median(e$mean_difference[e$chrom == 'chrX' & !e$PAR & e$primary_eligible])
    for (j in seq_len(nrow(XY_PAIRS))) {
      ix <- pick_gene(e, XY_PAIRS$X[j], 'chrX'); iy <- pick_gene(e, XY_PAIRS$Y[j], 'chrY')
      xtab[[length(xtab) + 1]] <- data.frame(line = L, X = XY_PAIRS$X[j], Y = XY_PAIRS$Y[j],
        X_mean_difference = if (is.na(ix)) NA else e$mean_difference[ix],
        X_CTRL_CPM = if (is.na(ix)) NA else ctrl[e$gene_id[ix]],
        Y_CTRL_CPM = if (is.na(iy)) NA else ctrl[e$gene_id[iy]],
        Y_percent_of_control = if (is.na(iy)) NA else pct(e$gene_id[iy]),
        chrX_median_shift = xmed)
    }
  }
  y <- do.call(rbind, ytab); p <- do.call(rbind, ptab); xp <- do.call(rbind, xtab); r <- do.call(rbind, rtab)
  rownames(y) <- rownames(p) <- rownames(xp) <- rownames(r) <- NULL
  y$status <- ifelse(is.na(y$gene_id), 'absent_from_annotation',
              ifelse(is.na(y$percent_of_control), 'not_expressed_in_control', 'expressed'))
  tab(y, 'Y_inventory_by_line'); tab(p, 'PAR_genes_by_line'); tab(xp, 'XY_paralogs_by_line'); tab(r, 'Y_PAR_region_sums_by_line')
  list(y = y, par = p, xy = xp, regions = r)
}
