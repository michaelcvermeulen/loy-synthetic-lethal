# Custom measured-reference enrichment; original RNA models and statistics are read-only.
options(warn = 1, stringsAsFactors = FALSE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/')))
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(limma)
  library(jsonlite)
})
setDTthreads(2)
set.seed(42)
out <- 'outputs/tables/chromatin_rna'
prov <- 'provenance/chromatin_rna'
partition_arg <- grep('^--partition=', commandArgs(TRUE), value = TRUE)
partition <- if (length(partition_arg)) as.integer(sub('^--partition=', '', partition_arg)) else 0L
stopifnot(partition %in% 0:4)
base_out <- out
cached_fg <- if (file.exists(file.path(base_out, 'fgsea_results.csv'))) fread(file.path(base_out, 'fgsea_results.csv')) else NULL
cached_camera <- if (file.exists(file.path(base_out, 'camera_results.csv'))) fread(file.path(base_out, 'camera_results.csv')) else NULL
if (partition > 0L) out <- file.path(out, 'partitions', paste0('part_', partition))
if (partition > 0L) {
  for (kind in c('fg', 'camera')) {
    local_path <- file.path(out, if (kind == 'fg') 'fgsea_results.csv' else 'camera_results.csv')
    if (file.exists(local_path)) {
      old <- if (kind == 'fg') cached_fg else cached_camera
      local <- fread(local_path)
      joined <- rbindlist(list(local, old), fill = TRUE)
      keys <- c('scenario', 'contrast', 'pathway', if (kind == 'camera') 'correlation_mode')
      joined <- unique(joined, by = keys)
      if (kind == 'fg') cached_fg <- joined else cached_camera <- joined
    }
  }
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(prov, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(out, paste0(name, '.csv')), na = 'NA')
hash <- function(p) digest::digest(file = p, algo = 'sha256')
cfg <- fromJSON('config/chromatin_rna.json')
inputs <- c('config/chromatin_rna.json', 'config/chromatin_atlas.json',
  'outputs/tables/chromatin_atlas/gene_set_catalogue.csv',
  'outputs/tables/chromatin_atlas/gene_set_membership_RNA_universe.csv',
  'outputs/tables/chromatin_atlas/gene_set_membership_all_scenarios.csv.gz',
  'outputs/tables/chromatin_atlas/HT1080_promoter_marks_RNA_universe.csv',
  'outputs/tables/chromatin_atlas/HT1080_GSE261834_measured_symbols.gmt',
  'outputs/tables/pathway_TF_effect_ranking.csv',
  'outputs/tables/HT1080_LOY_vs_sgCTRL_primary_effects.csv',
  'outputs/models/camera_voom_model.rds', 'outputs/models/01_ingest.rds',
  'outputs/tables/fgsea_all_collections.csv',
  'outputs/tables/camera_all_collections.csv',
  'outputs/tables/camera_estimated_correlation_sensitivity.csv')
manifest <- data.frame(path = inputs, sha256 = vapply(inputs, hash, character(1)))
rownames(manifest) <- NULL
mp <- file.path(prov, 'input_manifest.csv')
if (file.exists(mp)) stopifnot(identical(as.data.frame(fread(mp)), manifest)) else {
  fwrite(manifest, mp)
  write_json(list(frozen_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
    config_sha256 = hash(inputs[1]), decisions = cfg), file.path(prov, 'frozen_decisions.json'), pretty = TRUE)
}
# Freeze before reading RNA effects or running any association tests.
if ('--freeze-only' %in% commandArgs(TRUE)) quit(status = 0)
catalogue <- fread(inputs[3])
members <- fread(inputs[4])
all_members <- fread(inputs[5])
reference <- fread(inputs[6])
ranks <- fread(inputs[8])
effect <- fread(inputs[9])
effect <- effect[match(ranks$gene_id, effect$gene_id)]
stopifnot(identical(effect$gene_id, ranks$gene_id), identical(effect$gene_name, ranks$gene_name),
  nrow(ranks) == 25461L, !anyDuplicated(ranks$gene_id), !anyDuplicated(ranks$gene_name),
  max(abs(ranks$combined - effect$mean_difference)) < 1e-12,
  max(abs(ranks$combined - (ranks$AUG_2025 + ranks$MAY_2026) / 2)) < 1e-12)
model <- readRDS(inputs[10])
samples <- readRDS(inputs[11])$core
stopifnot(identical(rownames(model$voom$E), ranks$gene_name),
  identical(model$gene_id, ranks$gene_id), identical(colnames(model$voom$E), rownames(model$design)),
  identical(samples$sample_key, rownames(model$design)), model$residual_df == 4,
  qr(model$design)$rank == 8)
save_table(data.frame(sample_key = rownames(model$design), model$design), 'camera_design')
save_table(data.frame(coefficient = rownames(model$contrasts), model$contrasts), 'camera_contrasts')
save_table(samples, 'sample_metadata')
save_table(effect, 'gene_effect_evidence')
save_table(reference, 'reference_features')
sets_from <- function(scenario) {
  requested_scenario <- scenario
  m <- all_members[scenario == requested_scenario & gene_id %in% ranks$gene_id]
  # Gene-body definitions do not depend on promoter width or TSS choice.
  if (scenario != cfg$primary_scenario) m <- rbind(m,
    members[grepl('GENE_BODY', set_name)], fill = TRUE)
  z <- split(m$gene_name, m$set_name)
  setNames(lapply(catalogue$set_name, function(n) sort(unique(z[[n]]))), catalogue$set_name)
}
primary <- sets_from(cfg$primary_scenario)
stopifnot(identical(as.integer(lengths(primary)), catalogue$n_RNA_universe))
gmt <- fgsea::gmtPathways(inputs[7])
stopifnot(all(vapply(names(primary), function(n) setequal(primary[[n]], gmt[[paste0('HT1080_NEEDSEQ_', n)]]), logical(1))))
rank_vector <- function(v, names = ranks$gene_name) {
  v <- setNames(v, names)
  stopifnot(all(is.finite(v)), !anyDuplicated(names))
  v[order(-v, names(v))]
}
fg_results <- camera_results <- list()
counter <- 0L
run_fg <- function(sets, values, scenario, contrast) {
  counter <<- counter + 1L
  if (partition > 0L && (counter - 1L) %% 4L != partition - 1L) return(invisible(NULL))
  message('fgsea ', counter, ': ', scenario, ' / ', contrast)
  requested_scenario <- scenario; requested_contrast <- contrast
  cached <- if (!is.null(cached_fg)) cached_fg[scenario == requested_scenario & contrast == requested_contrast] else NULL
  if (!is.null(cached) && nrow(cached) == 20L &&
      all(cached$test_status != 'numerical_failure' | cached$nPermSimple_used >= cfg$fgsea_retry_nPermSimple)) {
    message('Reusing hash-verified completed test')
    cached[, `:=`(initial_seed = cfg$seed + counter,
      seed = cfg$seed + counter + fifelse(nPermSimple_used == cfg$fgsea_retry_nPermSimple, 10000L, 0L))]
    fg_results[[length(fg_results) + 1L]] <<- cached
    save_table(rbindlist(fg_results, fill = TRUE), 'fgsea_results')
    return(invisible(cached))
  }
  eligible <- lengths(lapply(sets, intersect, y = names(values))) >= cfg$minimum_set_size &
    lengths(lapply(sets, intersect, y = names(values))) < length(values)
  execute <- function(p, nperm, seed) {
    set.seed(seed)
    as.data.table(fgseaMultilevel(p, values, minSize = cfg$minimum_set_size,
      maxSize = length(values) - 1L, eps = cfg$fgsea_eps, nPermSimple = nperm,
      gseaParam = cfg$fgsea_gseaParam, scoreType = cfg$fgsea_scoreType,
      BPPARAM = BiocParallel::SerialParam()))
  }
  z <- execute(sets[eligible], cfg$fgsea_nPermSimple, cfg$seed + counter)
  z[, nPermSimple_used := cfg$fgsea_nPermSimple]
  failed <- z[!is.finite(pval), pathway]
  if (length(failed)) {
    message('Prespecified higher-permutation retry: ', paste(failed, collapse = ', '))
    # The public single-pathway API avoids the cumulative-prefix kernel's large-K cost.
    # Same test, set size, 100,000 draws and retry seed; only computation is distributed.
    retry <- data.table::rbindlist(BiocParallel::bplapply(sort(failed),
      function(name, all_sets, values, settings, retry_seed) {
        set.seed(retry_seed)
        as.data.frame(fgsea::fgseaMultilevel(all_sets[name], values,
          minSize = settings$minimum_set_size, maxSize = length(values) - 1L,
          eps = settings$fgsea_eps, nPermSimple = settings$fgsea_retry_nPermSimple,
          gseaParam = settings$fgsea_gseaParam, scoreType = settings$fgsea_scoreType,
          BPPARAM = BiocParallel::SerialParam()))
      }, all_sets = sets, values = values, settings = cfg,
      retry_seed = cfg$seed + counter + 10000L,
      BPPARAM = BiocParallel::SnowParam(workers = min(2L, length(failed)), RNGseed = cfg$seed)))
    retry[, nPermSimple_used := cfg$fgsea_retry_nPermSimple]
    z <- rbind(z[!pathway %in% failed], retry)
  }
  z[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ';')]
  z[, padj := NULL]
  z <- merge(data.table(pathway = names(sets), eligible = eligible,
    overlap_size = lengths(lapply(sets, intersect, y = names(values)))), z, by = 'pathway', all.x = TRUE)
  z[, `:=`(scenario = scenario, contrast = contrast, universe_size = length(values),
    initial_seed = cfg$seed + counter,
    seed = cfg$seed + counter + fifelse(nPermSimple_used == cfg$fgsea_retry_nPermSimple, 10000L, 0L), test_status = fifelse(!eligible, 'ineligible_size',
      fifelse(is.finite(pval), 'tested', 'numerical_failure')),
    p_at_lower_bound = is.finite(pval) & pval <= cfg$fgsea_eps,
    BH_custom_20 = p.adjust(pval, 'BH', n = 20), BY_custom_20 = p.adjust(pval, 'BY', n = 20))]
  fg_results[[length(fg_results) + 1L]] <<- z
  save_table(rbindlist(fg_results, fill = TRUE), 'fgsea_results')
  invisible(z)
}
run_camera <- function(sets, scenario, contrast, keep = rep(TRUE, nrow(ranks))) {
  if (partition > 0L && (counter - 1L) %% 4L != partition - 1L) return(invisible(NULL))
  requested_scenario <- scenario; requested_contrast <- contrast
  cached <- if (!is.null(cached_camera)) cached_camera[scenario == requested_scenario & contrast == requested_contrast] else NULL
  if (!is.null(cached) && nrow(cached) == 40L && !any(cached$test_status == 'numerical_failure')) {
    camera_results[[length(camera_results) + 1L]] <<- cached
    save_table(rbindlist(camera_results, fill = TRUE), 'camera_results')
    return(invisible(NULL))
  }
  v <- model$voom[keep, ]
  idx <- limma::ids2indices(sets, rownames(v$E), remove.empty = FALSE)
  eligible <- lengths(idx) >= cfg$minimum_set_size & lengths(idx) < nrow(v$E)
  for (mode in cfg$camera_modes) {
    rho <- if (mode == 'fixed_0.01') .01 else NA_real_
    z <- as.data.table(limma::camera(v, idx[eligible], design = model$design,
      contrast = model$contrasts[, contrast], inter.gene.cor = rho, allow.neg.cor = FALSE,
      use.ranks = FALSE, trend.var = FALSE, sort = FALSE), keep.rownames = 'pathway')
    z[, FDR := NULL]
    if (!'Correlation' %in% names(z)) z[, Correlation := rho]
    z <- merge(data.table(pathway = names(sets), eligible = eligible,
      overlap_size = lengths(idx)), z, by = 'pathway', all.x = TRUE)
    z[, `:=`(scenario = scenario, contrast = contrast, correlation_mode = mode,
      universe_size = nrow(v$E), residual_df = model$residual_df,
      effective_correlation = pmax(0, Correlation),
      BH_custom_20 = p.adjust(PValue, 'BH', n = 20), BY_custom_20 = p.adjust(PValue, 'BY', n = 20),
      test_status = fifelse(!eligible, 'ineligible_size', fifelse(is.finite(PValue), 'tested', 'numerical_failure')))]
    camera_results[[length(camera_results) + 1L]] <<- z
  }
  save_table(rbindlist(camera_results, fill = TRUE), 'camera_results')
}
# Combined is the primary test; period-specific contrasts are temporal sensitivities.
for (contrast in c('combined', 'AUG_2025', 'MAY_2026')) {
  run_fg(primary, rank_vector(ranks[[contrast]]), cfg$primary_scenario, contrast)
  run_camera(primary, cfg$primary_scenario, contrast)
}
message('Primary and period tests complete')
for (scenario in cfg$promoter_sensitivities) {
  sets <- sets_from(scenario)
  run_fg(sets, rank_vector(ranks$combined), scenario, 'combined')
  run_camera(sets, scenario, 'combined')
}
keep <- ranks$gene_id %in% reference[reference_status == 'evaluable_under_frozen_rules', gene_id]
stopifnot(sum(!keep) == 467L)
run_fg(primary, rank_vector(ranks$combined[keep], ranks$gene_name[keep]), 'evaluable_promoters_only', 'combined')
run_camera(primary, 'evaluable_promoters_only', 'combined', keep)
save_table(data.table(gene_id = ranks$gene_id, gene_name = ranks$gene_name,
  included_evaluable_sensitivity = keep), 'background_sensitivity')
for (column in c('prior025_mean_difference', 'RLE_mean_difference', grep('^omit_', names(effect), value = TRUE))) {
  run_fg(primary, rank_vector(effect[[column]]), column, 'combined')
}
# Descriptive confound sensitivities, with no gene-regression P-values or causal interpretation.
stopifnot(all(effect$feature_length_bp > 0), all(!is.na(effect$arm)))
fit1 <- lm(mean_difference ~ splines::ns(mean_logCPM, df = 4) +
  splines::ns(log10(feature_length_bp), df = 4), data = effect)
fit2 <- update(fit1, . ~ . + factor(arm))
residualized <- data.table(gene_id = ranks$gene_id, gene_name = ranks$gene_name,
  original = ranks$combined, abundance_length = residuals(fit1),
  abundance_length_arm = residuals(fit2))
save_table(residualized, 'residualized_effects')
run_fg(primary, rank_vector(residualized$abundance_length), 'residualized_abundance_length', 'combined')
run_fg(primary, rank_vector(residualized$abundance_length_arm), 'residualized_abundance_length_arm', 'combined')
fg <- rbindlist(fg_results, fill = TRUE)
ca <- rbindlist(camera_results, fill = TRUE)
if (partition > 0L) {
  save_table(fg, 'fgsea_results'); save_table(ca, 'camera_results')
  stopifnot(all(vapply(inputs, hash, character(1)) == manifest$sha256))
  write_json(list(status = 'partition_complete', partition = partition, fgsea_rows = nrow(fg), camera_rows = nrow(ca)),
    file.path(out, 'completion.json'), pretty = TRUE, auto_unbox = TRUE)
  quit(status = 0)
}
for (z in list(fg, ca)) {
  z[, analysis_role := fifelse(scenario == cfg$primary_scenario & contrast == 'combined', 'primary', 'sensitivity')]
  z[, BH_pooled_sensitivities := NA_real_]
}
fg[analysis_role == 'sensitivity', BH_pooled_sensitivities := p.adjust(pval, 'BH', n = .N)]
ca[analysis_role == 'sensitivity', BH_pooled_sensitivities := p.adjust(PValue, 'BH', n = .N), by = correlation_mode]
fg[, statistical_scope := 'conditional_preranked_gene_set_null; no_independent_population_replication']
ca[, statistical_scope := 'conditional_fixed_passage_model_assuming_independent_residuals; no_independent_population_replication']
save_table(fg, 'fgsea_results'); save_table(ca, 'camera_results')
saveRDS(list(fgsea = fg, camera = ca, sets = primary, config = cfg), 'outputs/models/chromatin_rna.rds')
writeLines(capture.output(sessionInfo()), file.path(prov, 'R_sessionInfo.txt'))
stopifnot(all(vapply(inputs, hash, character(1)) == manifest$sha256))
write_json(list(status = 'statistics_complete', fgsea_rows = nrow(fg), camera_rows = nrow(ca),
  fgsea_failed = sum(fg$test_status == 'numerical_failure'), camera_failed = sum(ca$test_status == 'numerical_failure'),
  primary_fgsea_BH05 = sum(fg$analysis_role == 'primary' & fg$BH_custom_20 < .05, na.rm = TRUE),
  completed_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE)),
  file.path(prov, 'statistics_completion.json'), pretty = TRUE, auto_unbox = TRUE)
message('All chromatin RNA tests complete')
