hdac_partition <- function(membership, sources) {
  genes <- sort(unique(membership$symbol))
  w <- data.frame(gene_name = genes)
  for (s in names(sources)) for (arm in c('UP', 'DN')) {
    w[[paste(s, arm, sep = '_')]] <- genes %in%
      membership$symbol[membership$pathway == paste(sources[[s]], arm, sep = '_')]
  }
  single1 <- w$HDAC1_UP | w$HDAC1_DN; single2 <- w$HDAC2_UP | w$HDAC2_DN
  any_up <- w$HDAC1_UP | w$HDAC2_UP | w$HDAC1_2_UP
  any_dn <- w$HDAC1_DN | w$HDAC2_DN | w$HDAC1_2_DN
  w$partition <- ifelse(any_up & any_dn, 'conflicting_direction',
    ifelse(single1 & single2, 'concordant_shared',
      ifelse(single2, 'HDAC2_only_vs_HDAC1', ifelse(single1, 'HDAC1_only_vs_HDAC2',
        ifelse(any_up | any_dn, 'double_knockdown_only', 'HDAC3_comparison_only')))))
  w$partition_arm <- ifelse(any_up & any_dn, 'MIXED', ifelse(any_up, 'UP', ifelse(any_dn, 'DN', 'HDAC3_only')))
  w$double_knockdown_member <- w$HDAC1_2_UP | w$HDAC1_2_DN
  w$HDAC3_member <- w$HDAC3_UP | w$HDAC3_DN
  w
}

hdac_percentiles <- function(counts) {
  assert(nrow(counts) > 1 && !anyNA(counts), 'Invalid fixed rank background')
  z <- apply(counts, 2, function(v) 100 * (rank(v, ties.method = 'average') - 1) / (length(v) - 1))
  dimnames(z) <- dimnames(counts)
  z
}

hdac_signed_score <- function(values, up, dn) {
  ambiguous <- intersect(up, dn)
  up <- setdiff(intersect(up, rownames(values)), ambiguous)
  dn <- setdiff(intersect(dn, rownames(values)), ambiguous)
  assert(length(up) > 0 && length(dn) > 0, 'Empty HDAC scoring arm')
  u <- colMeans(values[up, , drop = FALSE]); v <- colMeans(values[dn, , drop = FALSE])
  rbind(UP = u, DN = v, signed = u - v)
}

hdac_summary <- function(values, genes, pairs) {
  genes <- intersect(genes, rownames(values))
  if (!length(genes)) return(data.frame(contrast = colnames(values), genes = 0L,
    mean_difference = NA_real_, median_difference = NA_real_, fraction_lower = NA_real_))
  z <- values[genes, , drop = FALSE]
  data.frame(contrast = colnames(z), genes = length(genes), mean_difference = colMeans(z),
    median_difference = apply(z, 2, median), fraction_lower = colMeans(z < 0))
}

hdac12_followup <- function(x, e, p) {
  cfg <- yaml::read_yaml('config/hdac12_followup.yml'); set.seed(cfg$seed)
  universe <- e$mapping$gene_name
  assert(identical(universe, rownames(p$effects_symbol)) && length(universe) == 25461L,
    'HDAC fixed expressed-symbol universe changed')
  all_sets <- do.call(c, unname(p$sets))
  membership <- as.data.frame(data.table::fread('outputs/tables/HDAC12_frozen_membership.csv'))
  wanted <- unlist(lapply(cfg$sources, function(s) paste(s, c('UP', 'DN'), sep = '_')), use.names = FALSE)
  membership <- membership[membership$pathway %in% wanted, ]
  for (name in wanted) assert(setequal(membership$symbol[membership$pathway == name], all_sets[[name]]),
    paste('Audited and frozen HDAC memberships disagree:', name))
  partitions <- hdac_partition(membership, cfg$sources)
  partitions$represented <- partitions$gene_name %in% universe
  partitions$gene_id <- e$mapping$gene_id[match(partitions$gene_name, universe)]
  tab(partitions, 'HDAC12_membership_partitions')
  tab(data.frame(gene_name = universe, gene_id = e$mapping$gene_id), 'HDAC12_score_background')
  values <- cbind(p$effects_symbol, e$differences[e$mapping$gene_id, ], e$loo[e$mapping$gene_id, ])
  rownames(values) <- universe
  annotation <- e$effects[match(e$mapping$gene_id, e$effects$gene_id), ]
  assert(identical(annotation$gene_name, universe), 'HDAC locus join changed symbol order')
  build_programs <- function(def) lapply(def, function(n) {
    assert(all(n %in% names(all_sets)), paste('Missing process reference:', paste(n, collapse = ';')))
    intersect(unique(unlist(all_sets[n], use.names = FALSE)), universe)
  })
  masks <- build_programs(cfg$existing_masks); processes <- build_programs(cfg$process_annotations)
  legacy <- readRDS('outputs/models/04_chromatin_followup.rds')
  for (name in names(masks)) assert(setequal(masks[[name]], legacy$programs[[name]]), 'Existing chromatin mask changed')
  program_members <- do.call(rbind, lapply(c('existing_masks', 'process_annotations'), function(kind) {
    do.call(rbind, lapply(names(cfg[[kind]]), function(name) {
      do.call(rbind, lapply(cfg[[kind]][[name]], function(ref)
        data.frame(kind = kind, program = name, pathway = ref, gene_name = all_sets[[ref]],
          represented = all_sets[[ref]] %in% universe)))
    }))
  }))
  tab(program_members, 'HDAC12_process_reference_membership')
  annotated <- annotation
  for (name in names(masks)) annotated[[paste0('in_', name)]] <- universe %in% masks[[name]]
  for (name in names(processes)) annotated[[paste0('process_', name)]] <- universe %in% processes[[name]]
  annotated$outside_existing_masks <- !universe %in% unique(unlist(masks))
  annotated$outside_all_process_annotations <- !universe %in% unique(unlist(processes))
  annotated$prior_count_change <- annotated$prior025_mean_difference - annotated$mean_difference
  evidence <- do.call(rbind, lapply(names(cfg$sources), function(source) {
    prefix <- cfg$sources[[source]]
    ambiguous <- intersect(all_sets[[paste0(prefix, '_UP')]], all_sets[[paste0(prefix, '_DN')]])
    do.call(rbind, lapply(c('UP', 'DN'), function(arm) {
      name <- paste(prefix, arm, sep = '_'); genes <- intersect(all_sets[[name]], universe)
      z <- annotated[match(genes, universe), ]
      z$source <- source; z$source_arm <- arm; z$pathway <- name
      z$within_source_ambiguous <- genes %in% ambiguous
      for (contrast in c('combined', 'AUG_2025', 'MAY_2026')) {
        le <- p$fgsea$leadingEdge[p$fgsea$pathway == name & p$fgsea$contrast == contrast]
        z[[paste0('leading_edge_', contrast)]] <- if (length(le) == 1 && !is.na(le))
          genes %in% strsplit(le, ';', fixed = TRUE)[[1]] else rep(NA, length(genes))
      }
      z$perturbed_regulator <- genes %in% cfg$regulators
      z$partition <- partitions$partition[match(genes, partitions$gene_name)]
      z$partition_arm <- partitions$partition_arm[match(genes, partitions$gene_name)]
      z$double_knockdown_member <- partitions$double_knockdown_member[match(genes, partitions$gene_name)]
      z$HDAC3_member <- partitions$HDAC3_member[match(genes, partitions$gene_name)]
      z
    }))
  }))
  tab(evidence, 'HDAC12_gene_evidence')
  coverage <- as.data.frame(table(partitions$partition, partitions$partition_arm, partitions$represented))
  names(coverage) <- c('partition', 'source_arm', 'represented', 'genes')
  tab(coverage, 'HDAC12_partition_coverage')
  represented <- partitions[partitions$represented, ]
  partition_summary <- do.call(rbind, lapply(split(represented,
    paste(represented$partition, represented$partition_arm, sep = '__')), function(z) {
      data.frame(partition = z$partition[1], source_arm = z$partition_arm[1],
        hdac_summary(values, z$gene_name, x$pairs))
    }))
  tab(partition_summary, 'HDAC12_partition_effects')
  # Exact UP/DN memberships are frozen independently of HT1080 effect direction.
  frozen <- evidence[, c('source', 'source_arm', 'pathway', 'gene_name', 'gene_id', 'within_source_ambiguous')]
  frozen$used_in_score <- !frozen$within_source_ambiguous
  tab(frozen, 'HDAC12_frozen_scoring_membership')
  rank_values <- hdac_percentiles(x$counts[e$mapping$gene_id, x$core$sample_key, drop = FALSE])
  rownames(rank_values) <- universe
  rank_matched <- matched_differences(rank_values, x$pairs)
  rank_loo <- sapply(seq_len(nrow(x$pairs)), function(j)
    rowMeans(period_average(rank_matched[, -j, drop = FALSE], x$pairs[-j, ])))
  colnames(rank_loo) <- paste0('omit_', x$pairs$block)
  rank_effects <- cbind(combined = rowMeans(period_average(rank_matched, x$pairs)),
    period_average(rank_matched, x$pairs), rank_matched, rank_loo)
  assert(identical(dimnames(rank_effects), dimnames(values)), 'Rank and expression contrasts are misaligned')
  score_list <- lapply(names(cfg$sources), function(source) {
    z <- frozen[frozen$source == source & frozen$used_in_score, ]
    a <- hdac_signed_score(rank_values, z$gene_name[z$source_arm == 'UP'], z$gene_name[z$source_arm == 'DN'])
    dimnames(a)[[1]] <- paste(source, rownames(a), sep = '__'); a
  })
  scores <- do.call(rbind, score_list)
  samples <- as.data.frame(matrix_long(scores, 'percentile_score', 'program'))
  samples$source <- sub('__.*', '', samples$program); samples$score_arm <- sub('.*__', '', samples$program)
  samples <- merge(samples, x$core, by = 'sample_key', sort = FALSE)
  tab(samples, 'HDAC12_sample_percentile_scores')
  matched <- matched_differences(scores, x$pairs)
  score_contrasts <- cbind(combined = rowMeans(period_average(matched, x$pairs)), period_average(matched, x$pairs), matched)
  loo <- sapply(seq_len(nrow(x$pairs)), function(j)
    rowMeans(period_average(matched[, -j, drop = FALSE], x$pairs[-j, ])))
  colnames(loo) <- paste0('omit_', x$pairs$block)
  score_contrasts <- cbind(score_contrasts, loo)
  score_effects <- as.data.frame(matrix_long(score_contrasts, 'percentile_difference', 'program'))
  names(score_effects)[2] <- 'contrast'
  score_effects$source <- sub('__.*', '', score_effects$program)
  score_effects$score_arm <- sub('.*__', '', score_effects$program)
  tab(score_effects, 'HDAC12_percentile_score_contrasts')
  arm_summary <- list(); influence <- list(); contributions <- list(); chromosome <- list(); process_effects <- list()
  for (source in names(cfg$sources)) {
    z <- evidence[evidence$source == source & !evidence$within_source_ambiguous, ]
    up <- z$gene_name[z$source_arm == 'UP']; dn <- z$gene_name[z$source_arm == 'DN']
    for (arm in c('UP', 'DN', 'signed')) {
      target <- if (arm == 'UP') up else if (arm == 'DN') dn else c(up, dn)
      weight <- if (arm == 'signed') ifelse(target %in% up, 1 / length(up), -1 / length(dn)) else rep(1 / length(target), length(target))
      contrib <- data.frame(source = source, summary_arm = arm, gene_name = target,
        source_arm = ifelse(target %in% up, 'UP', 'DN'), weight = weight,
        combined_difference = values[target, 'combined'], contribution = weight * values[target, 'combined'])
      contrib <- contrib[order(-abs(contrib$contribution), contrib$gene_name), ]
      contrib$absolute_rank <- seq_len(nrow(contrib))
      contrib$absolute_effect_fraction <- abs(contrib$contribution) / sum(abs(contrib$contribution))
      contributions[[paste(source, arm)]] <- contrib
      scenarios <- list(original = character(), omit_HDAC1_HDAC2 = c('HDAC1', 'HDAC2'),
        omit_HDAC1_HDAC2_HDAC3 = cfg$regulators, omit_IFN = masks$IFN, omit_ECM_EMT = masks$ECM_EMT,
        omit_IFN_ECM = unique(unlist(masks[c('IFN', 'ECM_EMT')])), omit_IFN_ECM_cycle = unique(unlist(masks)))
      for (k in cfg$top_contributors) scenarios[[paste0('omit_top_', k)]] <- head(contrib$gene_name, k)
      for (cut in cfg$abundance_thresholds_CPM) scenarios[[paste0('retain_CPM_at_least_', cut)]] <- z$gene_name[z$mean_CPM < cut]
      if (arm != 'signed') {
        chr <- annotation$chrom[match(target, universe)]
        ch <- do.call(rbind, lapply(sort(unique(chr)), function(c) {
          j <- chr == c
          data.frame(source = source, source_arm = arm, chrom = c, genes = sum(j), member_fraction = mean(j),
            absolute_effect_fraction = sum(abs(values[target[j], 'combined'])) / sum(abs(values[target, 'combined'])))
        }))
        ch$flagged <- ch$member_fraction >= cfg$chromosome_concentration_threshold |
          ch$absolute_effect_fraction >= cfg$chromosome_concentration_threshold
        chromosome[[paste(source, arm)]] <- ch
        for (c in ch$chrom[ch$flagged]) scenarios[[paste0('omit_', c)]] <- target[chr == c]
        arm_summary[[paste(source, arm)]] <- data.frame(source = source, source_arm = arm,
          hdac_summary(values, target, x$pairs))
        flags <- c(masks, processes, list(any_existing_mask = unique(unlist(masks)),
          any_process_annotation = unique(unlist(processes))))
        process_effects[[paste(source, arm)]] <- do.call(rbind, lapply(names(flags), function(name) {
          do.call(rbind, lapply(c(TRUE, FALSE), function(inside) {
            genes <- target[(target %in% flags[[name]]) == inside]
            ans <- hdac_summary(values[, 1:3, drop = FALSE], genes, x$pairs)
            data.frame(source = source, source_arm = arm, annotation = name, inside_annotation = inside,
              represented_arm_genes = length(target), absolute_effect_fraction =
                sum(abs(values[genes, 'combined'])) / sum(abs(values[target, 'combined'])), ans)
          }))
        }))
      }
      influence[[paste(source, arm)]] <- do.call(rbind, lapply(names(scenarios), function(scenario) {
        removed <- intersect(target, scenarios[[scenario]]); retained <- setdiff(target, removed)
        u <- setdiff(up, removed); v <- setdiff(dn, removed)
        metric <- function(mat) {
          if (arm == 'signed') {
            if (!length(u) || !length(v)) return(rep(NA_real_, ncol(mat)))
            colMeans(mat[u, , drop = FALSE]) - colMeans(mat[v, , drop = FALSE])
          } else if (length(retained)) colMeans(mat[retained, , drop = FALSE]) else rep(NA_real_, ncol(mat))
        }
        prior <- matrix(annotation$prior025_mean_difference, ncol = 1, dimnames = list(universe, 'prior025'))
        rle <- matrix(annotation$RLE_mean_difference, ncol = 1, dimnames = list(universe, 'RLE'))
        ans <- data.frame(source = source, summary_arm = arm, scenario = scenario,
          original_genes = length(target), retained_genes = length(retained),
          retained_UP = if (arm != 'DN') length(u) else 0L, retained_DN = if (arm != 'UP') length(v) else 0L,
          removed_genes = paste(sort(removed), collapse = ';'),
          removed_absolute_effect_fraction = sum(contrib$absolute_effect_fraction[contrib$gene_name %in% removed]),
          contrast = colnames(values), mean_difference = metric(values),
          prior025_combined = metric(prior), RLE_combined = metric(rle),
          percentile_difference = metric(rank_effects))
        ans
      }))
    }
  }
  arm_summary <- do.call(rbind, arm_summary); influence <- do.call(rbind, influence)
  contributions <- do.call(rbind, contributions); chromosome <- do.call(rbind, chromosome)
  process_effects <- do.call(rbind, process_effects)
  tab(arm_summary, 'HDAC12_source_arm_effects'); tab(influence, 'HDAC12_influence_sensitivity')
  tab(contributions, 'HDAC12_gene_contributions'); tab(chromosome, 'HDAC12_chromosome_concentration')
  tab(process_effects, 'HDAC12_process_effects')
  # Outcome-informed readouts are explicitly separate from the complete frozen scoring sets.
  candidates <- evidence[evidence$source == 'HDAC2' & !evidence$within_source_ambiguous, ]
  rules <- cfg$shortlist
  candidates$eligible_abundance <- candidates$mean_CPM >= rules$min_mean_CPM
  candidates$eligible_prior_stability <- abs(candidates$prior_count_change) <= rules$max_absolute_prior_count_change
  base <- candidates$eligible_abundance & candidates$eligible_prior_stability & candidates$outside_existing_masks
  negative <- base & candidates$maximum_match < 0 & candidates$LOO_max < 0 &
    candidates$prior025_mean_difference < 0 & candidates$RLE_mean_difference < 0
  near <- base & apply(abs(as.matrix(candidates[, c('mean_difference', 'AUG_2025', 'MAY_2026',
    'LOO_min', 'LOO_max', 'prior025_mean_difference', 'RLE_mean_difference')])), 1, max) <= rules$max_nonresponse_effect
  candidates$eligible_DN_readout <- negative & candidates$source_arm == 'DN'
  candidates$eligible_UP_opposite <- negative & candidates$source_arm == 'UP'
  candidates$eligible_UP_near_zero <- near & candidates$source_arm == 'UP'
  selections <- lapply(c('DN_readout', 'UP_opposite', 'UP_near_zero'), function(kind) {
    z <- candidates[candidates[[paste0('eligible_', kind)]], ]
    z <- if (kind == 'UP_near_zero') z[order(-z$mean_CPM, z$gene_name), ] else z[order(z$mean_difference, z$gene_name), ]
    n <- rules[[if (kind == 'DN_readout') 'DN_count' else paste0(kind, '_count')]]
    z <- head(z, n); z$readout_category <- rep(kind, nrow(z)); z
  })
  shortlist <- do.call(rbind, selections)
  candidates$selected_readout <- candidates$gene_id %in% shortlist$gene_id
  tab(candidates, 'HDAC12_readout_eligibility'); tab(shortlist, 'HDAC12_readout_shortlist')
  key <- influence[influence$scenario == 'original' & influence$summary_arm == 'signed' & influence$contrast == 'combined', ]
  previous <- legacy$pairs[legacy$pairs$contrast == 'combined', ]
  assert(max(abs(key$mean_difference - previous$source_response_contrast[
    match(unlist(cfg$sources[key$source], use.names = FALSE), previous$source_pair)])) < 1e-10,
    'New signed summaries do not recover the previously reported source contrasts')
  list(cfg = cfg, partitions = partitions, coverage = coverage, partition_summary = partition_summary,
    evidence = evidence, frozen = frozen, background = e$mapping[, c('gene_id', 'gene_name')],
    samples = samples, score_effects = score_effects, arm_summary = arm_summary,
    influence = influence, contributions = contributions, chromosome = chromosome,
    process_effects = process_effects, shortlist = shortlist, candidates = candidates,
    pairs = x$pairs, programs = processes, masks = masks)
}
