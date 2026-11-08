# Replication of selected shared and divergent programmes in the in-house lines and external LOY datasets.
# Metric: AUC that a member gene's Y-loss-polarity effect exceeds that of a non-member gene
# (AUC < 0.5 = members lower with Y absent), with an empirical two-sided P from size-matched random sets.
auc_of <- function(x, member) {
  r <- rank(x); n1 <- sum(member); n2 <- sum(!member)
  (sum(r[member]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

load_effects <- function(cfg, x) {
  out <- list()
  for (L in LINES) {
    e <- x$eff[[L]]; e <- e[e$primary_eligible & e$gene_name != '', ]
    e <- e[!duplicated(e$gene_name), ]
    out[[L]] <- list(group = 'in_house', description = paste(L, 'thesis reanalysis, primary genes'),
                     sym = e$gene_name, lfc = e$mean_difference, fdr = NA)
  }
  reg <- as.data.frame(data.table::fread(cfg$external_registry))
  for (i in seq_len(nrow(reg))) {
    r <- reg[i, ]
    d <- as.data.frame(data.table::fread(reg_path(r$path)))
    sym <- sub('\\.[0-9]+$', '', as.character(d[[r$symbol_col]]))
    if (isTRUE(as.logical(r$upper_case))) sym <- toupper(sym)
    lfc <- as.numeric(d[[r$lfc_col]]) * as.numeric(r$polarity)
    fdr <- if (!is.na(r$fdr_col) && nzchar(r$fdr_col)) as.numeric(d[[r$fdr_col]]) else rep(NA_real_, nrow(d))
    ok <- !is.na(sym) & nzchar(sym) & is.finite(lfc)
    sym <- sym[ok]; lfc <- lfc[ok]; fdr <- fdr[ok]
    first <- !duplicated(sym)
    out[[r$key]] <- list(group = r$group, description = r$description, sym = sym[first], lfc = lfc[first], fdr = fdr[first])
  }
  out
}

external_replication <- function(x, cfg) {
  msig <- readRDS(file.path(REPO, cfg$msigdb_rds))
  all_sets <- do.call(c, unname(msig))
  assert(all(cfg$replication_sets %in% names(all_sets)), 'Replication set missing from MSigDB resource')
  data <- load_effects(cfg, x)
  set.seed(cfg$seed)
  rows <- list(); focus <- list(); pol <- list()
  for (k in names(data)) {
    D <- data[[k]]
    for (s in cfg$replication_sets) {
      member <- D$sym %in% all_sets[[s]]
      n1 <- sum(member)
      if (n1 < 10) { rows[[length(rows) + 1]] <- data.frame(dataset = k, group = D$group, gene_set = s, n_members_tested = n1,
        AUC = NA, median_member = NA, median_other = NA, empirical_P = NA); next }
      a <- auc_of(D$lfc, member)
      null <- replicate(cfg$external_random_sets, { m <- rep(FALSE, length(D$lfc)); m[sample.int(length(D$lfc), n1)] <- TRUE; auc_of(D$lfc, m) })
      p <- (sum(abs(null - .5) >= abs(a - .5)) + 1) / (length(null) + 1)
      rows[[length(rows) + 1]] <- data.frame(dataset = k, group = D$group, gene_set = s, n_members_tested = n1, AUC = a,
        median_member = median(D$lfc[member]), median_other = median(D$lfc[!member]), empirical_P = p)
    }
    for (g in cfg$replication_genes) {
      i <- match(g, D$sym)
      focus[[length(focus) + 1]] <- data.frame(dataset = k, group = D$group, gene = g,
        effect_Yloss_polarity = if (is.na(i)) NA else D$lfc[i], FDR = if (is.na(i)) NA else D$fdr[i])
    }
    # Polarity: the first male-specific Y gene present must be lower with Y absent. In-house lines are
    # checked on their full effect tables (the primary universe excludes chrY by design). Many external
    # DE tables removed Y genes before testing; their direction then rests on the contrast label.
    if (k %in% LINES) { e <- x$eff[[k]]; ysym <- e$gene_name[e$chrom == 'chrY' & !e$PAR]; yl <- e$mean_difference[e$chrom == 'chrY' & !e$PAR] }
    else { ysym <- D$sym; yl <- D$lfc }
    hit <- match(cfg$y_polarity_genes, ysym); hit <- hit[!is.na(hit)]
    g <- if (length(hit)) ysym[hit[1]] else NA
    v <- if (length(hit)) yl[hit[1]] else NA
    pol[[k]] <- data.frame(dataset = k, group = D$group, check_gene = g, effect_Yloss_polarity = v,
                           polarity_ok = if (is.na(v)) NA else v < 0,
                           status = if (is.na(v)) 'not_checkable_Y_genes_absent_direction_from_contrast_label' else 'verified_by_Y_gene',
                           genes_tested = length(D$sym))
  }
  r <- do.call(rbind, rows); f <- do.call(rbind, focus); pc <- do.call(rbind, pol)
  rownames(r) <- rownames(f) <- rownames(pc) <- NULL
  r$direction <- ifelse(is.na(r$AUC), NA, ifelse(r$AUC < .5, 'lower_with_Y_absent', 'higher_with_Y_absent'))
  tab(r, 'external_replication_sets'); tab(f, 'external_focus_genes'); tab(pc, 'external_polarity_check')
  tab(data.frame(dataset = names(data), group = sapply(data, `[[`, 'group'), description = sapply(data, `[[`, 'description'),
                 genes_tested = sapply(data, function(D) length(D$sym))), 'external_dataset_registry_used')
  list(sets = r, focus = f, polarity = pc)
}
