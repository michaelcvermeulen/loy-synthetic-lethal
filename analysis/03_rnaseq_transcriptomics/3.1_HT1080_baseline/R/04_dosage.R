dosage_analysis <- function(x, e, p, cfg) {
  a <- x$annotation; m <- x$core; pairs <- x$pairs
  regions <- c('chrYp_nonPAR', 'chrYq_nonPAR', 'chrX_PAR', 'chrY_PAR', 'chrM')
  regional <- t(vapply(regions, function(r) colSums(e$norm$CPM[a$region == r, , drop = FALSE]), numeric(nrow(m))))
  rownames(regional) <- regions; colnames(regional) <- m$sample_key
  reg <- merge(as.data.frame(matrix_long(regional, 'sum_CPM', 'region')), m, by = 'sample_key')
  tab(reg, 'Y_PAR_MT_region_abundance')
  y <- e$effects[a$chrom == 'chrY' | a$PAR, ]
  tab(y, 'Y_PAR_all_gene_effects')
  # Prespecified widely studied X/Y pairs; one representative per symbol per chromosome.
  pp <- data.frame(X = c('KDM6A', 'KDM5C', 'DDX3X', 'EIF1AX', 'RPS4X', 'ZFX', 'USP9X', 'TBL1X', 'PRKX', 'TXLNG', 'NLGN4X'),
                    Y = c('UTY', 'KDM5D', 'DDX3Y', 'EIF1AY', 'RPS4Y1', 'ZFY', 'USP9Y', 'TBL1Y', 'PRKY', 'TXLNGY', 'NLGN4Y'))
  find_locus <- function(sym, chr) {
    z <- which(a$gene_name == sym & a$chrom == chr & !a$PAR)
    if (!length(z)) return(NA_character_)
    a$gene_id[z[order(-e$effects$mean_CPM[z], a$gene_id[z])][1]]
  }
  pp$X_gene_id <- vapply(pp$X, find_locus, character(1), chr = 'chrX')
  pp$Y_gene_id <- vapply(pp$Y, find_locus, character(1), chr = 'chrY')
  pp$pair <- paste(pp$X, pp$Y, sep = ' / ')
  pp$status <- ifelse(is.na(pp$X_gene_id) | is.na(pp$Y_gene_id), 'missing_counting_feature', 'available')
  tab(pp, 'paralog_pair_definitions')
  # RPKM-like proxy: common TMM CPM divided by union-exon kilobases. No separate Y normalization.
  abundance <- e$norm$CPM / (a$feature_length_bp / 1000)
  out <- list(); single <- list()
  for (i in which(pp$status == 'available')) {
    ix <- match(pp$X_gene_id[i], a$gene_id); iy <- match(pp$Y_gene_id[i], a$gene_id)
    cx <- abundance[ix, pairs$CTRL]; cy <- abundance[iy, pairs$CTRL]
    lx <- abundance[ix, pairs$LOY]; ly <- abundance[iy, pairs$LOY]
    z <- data.frame(pair = pp$pair[i], X_gene_id = pp$X_gene_id[i], Y_gene_id = pp$Y_gene_id[i],
                    pairs, X_CTRL = cx, Y_CTRL = cy, X_LOY = lx, Y_LOY = ly)
    z$relative_total_proxy <- relative_pair_total(lx, ly, cx, cy)
    z$CTRL_Y_fraction_proxy <- ifelse(cx + cy > 0, cy / (cx + cy), NA_real_)
    z$log2_relative_total_proxy <- log2(z$relative_total_proxy)
    # Flag denominators with weak count support rather than stabilizing silently.
    z$CTRL_pair_reads <- x$counts[ix, pairs$CTRL] + x$counts[iy, pairs$CTRL]
    z$low_denominator_support <- z$CTRL_pair_reads < 10
    out[[i]] <- z
    for (j in c(ix, iy)) {
      single[[paste(i, j)]] <- data.frame(pair = pp$pair[i], gene_id = a$gene_id[j], gene_name = a$gene_name[j],
        chrom = a$chrom[j], m, normalized_length_proxy = abundance[j, ], logCPM = e$logCPM[j, ], raw_count = x$counts[j, m$sample_key])
    }
  }
  pair_data <- do.call(rbind, out); rownames(pair_data) <- NULL
  single_data <- do.call(rbind, single); rownames(single_data) <- NULL
  tab(pair_data, 'paralog_matched_relative_total_proxy'); tab(single_data, 'paralog_sample_abundance')
  # August-only discovery including its own filtering, TMM and symbol selection.
  am <- m[m$batch == 'AUG_2025', ]; ac <- x$counts[, am$sample_key, drop = FALSE]
  an <- normalize_counts(ac, a$main_background, factor(am$treatment), cfg)
  ap <- pair_indices(am); ad <- matched_differences(an$logCPM, ap); amean <- rowMeans(ad)
  amap <- choose_symbols(a, rowMeans(an$logCPM), an$keep)
  amap <- amap[amap$representative, ]; ai <- match(amap$gene_id, a$gene_id)
  amap$AUG_difference <- amean[ai]
  amap$all_three_same_direction <- rowSums(sign(ad[ai, , drop = FALSE]) == sign(amean[ai])) == 3
  eligible <- amap$all_three_same_direction & abs(amap$AUG_difference) >= cfg$effect_threshold
  select_direction <- function(s) {
    z <- amap[eligible & sign(amap$AUG_difference) == s, ]
    head(z[order(-abs(z$AUG_difference), z$gene_id), ], cfg$signature_max_genes_each_direction)
  }
  signature <- rbind(select_direction(1), select_direction(-1))
  signature$direction <- sign(signature$AUG_difference)
  assert(all(table(signature$direction) > 0) && length(unique(signature$direction)) == 2,
         'August signature needs genes in both directions')
  tab(signature, 'August_frozen_signature_genes')
  tab(amap, 'August_signature_discovery_background')
  # Within-sample ranks of raw counts are invariant to multiplying a whole library by an exposure factor.
  rank_matrix <- apply(x$counts[amap$gene_id, m$sample_key, drop = FALSE], 2,
                        function(v) (rank(v, ties.method = 'average') - 1) / (length(v) - 1))
  rownames(rank_matrix) <- amap$gene_id
  up <- signature$gene_id[signature$direction > 0]; down <- signature$gene_id[signature$direction < 0]
  score <- colMeans(rank_matrix[up, , drop = FALSE]) - colMeans(rank_matrix[down, , drop = FALSE])
  signature_scores <- data.frame(m, signed_percentile_score = score)
  tab(signature_scores, 'August_signature_sample_scores')
  frozen <- list(discovery_samples = am$sample_key, background_ids = amap$gene_id,
                  selected = signature, score_method = 'mean within-sample percentile(up) minus mean percentile(down)',
                  cutoff = cfg$effect_threshold, max_each_direction = cfg$signature_max_genes_each_direction)
  saveRDS(frozen, 'outputs/models/August_frozen_signature.rds')
  # Record DNA metadata without assuming D0 specimens correspond to these passages.
  cn <- readRDS(cfg$cn_source)
  cn_status <- unique(as.data.frame(cn[, c('sample', 'label')]))
  cn_status$RNA_sample_matching <- 'unconfirmed'
  cn_status$RNA_integration_status <- 'deferred_pending_specimen_and_assembly_correspondence'
  tab(cn_status, 'DNA_specimen_context_status')
  extension <- data.frame(analysis = c('DNA_RNA_coordinate_integration', 'addback_rescue', 'ATAC_RNA_integration', 'formal_gene_DE', 'msVIPER'),
    status = c('deferred', 'optional_not_run', 'optional_not_run', 'not_applicable_primary_design', 'optional_not_run'),
    reason = c('D0 WGS specimens cannot be assigned to RNA passages; WGS T2T versus RNA GRCh38 requires compatible mapping',
      'Outside confirmed core; intervention replication and matching require separate design verification',
      'Outside confirmed core; sample and contrast compatibility require separate verification',
      'One established population per condition from one generation history; serial passages do not establish population replication',
      'CollecTRI ULM implemented; combinatorial regulator inference is an optional extension'))
  tab(extension, 'analysis_status')
  # An evidence table keeps dimensions separate; unmeasured evidence is explicitly missing.
  evidence <- e$effects[e$effects$descriptive_candidate, ]
  le <- split(p$leading$pathway, p$leading$gene_id)
  tf_targets <- split(p$targets$source, p$targets$gene_id)
  lookup <- function(id, lists) if (is.null(lists[[id]])) '' else paste(sort(unique(lists[[id]])), collapse = ';')
  evidence$leading_edge_pathways <- vapply(evidence$gene_id, lookup, character(1), lists = le)
  evidence$collectri_regulators <- vapply(evidence$gene_id, lookup, character(1), lists = tf_targets)
  evidence$gene_level_inferential_support <- 'not_estimated_for_independent_populations'
  evidence$DNA_dosage_support <- NA_character_; evidence$addback_reversal <- NA_real_; evidence$ATAC_support <- NA_real_
  tab(evidence, 'candidate_evidence_dimensions')
  list(regional = reg, pairs = pair_data, single = single_data, pair_definitions = pp,
       signature = signature, signature_scores = signature_scores, evidence = evidence, extension_status = extension)
}
