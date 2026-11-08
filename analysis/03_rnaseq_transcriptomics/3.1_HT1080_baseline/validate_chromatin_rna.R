# Independently reconstruct sample blocking and verify saved CAMERA test alignment.
options(stringsAsFactors = FALSE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/')))
library(limma)
library(data.table)
meta <- fread('outputs/tables/chromatin_rna/sample_metadata.csv')
model <- readRDS('outputs/models/camera_voom_model.rds')
results <- fread('outputs/tables/chromatin_rna/camera_results.csv')
if (nrow(results[scenario == 'window_2000']) < 120L) {
  parts <- list.files('outputs/tables/chromatin_rna/partitions', pattern = 'camera_results.csv', recursive = TRUE, full.names = TRUE)
  results <- rbindlist(lapply(parts, fread), fill = TRUE)
}
members <- fread('outputs/tables/chromatin_atlas/gene_set_membership_RNA_universe.csv')
genes <- fread('outputs/tables/pathway_TF_effect_ranking.csv')
stopifnot(identical(meta$sample_key, colnames(model$voom$E)))
blocks <- paste(meta$batch, paste0('passage_', meta$passage_number), sep = '__')
design <- sapply(sort(unique(blocks)), function(b) as.numeric(blocks == b))
design <- cbind(design, as.numeric(meta$treatment == 'LOY' & meta$batch == 'AUG_2025'),
  as.numeric(meta$treatment == 'LOY' & meta$batch == 'MAY_2026'))
stopifnot(isTRUE(all.equal(unname(design), unname(model$design))), nrow(design) - qr(design)$rank == 4)
for (block in unique(blocks)) stopifnot(setequal(meta$treatment[blocks == block], c('CTRL', 'LOY')))
indices <- lapply(split(members$gene_id, members$set_name), function(ids) which(genes$gene_id %in% ids))
checks <- 0L
for (contrast in c('combined', 'AUG_2025', 'MAY_2026')) {
  cv <- c(rep(0, 6), switch(contrast, combined = c(.5, .5), AUG_2025 = c(1, 0), MAY_2026 = c(0, 1)))
  for (mode in c('fixed_0.01', 'estimated_nonnegative')) {
    recomputed <- camera(model$voom$E, indices, design = design, contrast = cv,
      weights = model$voom$weights, inter.gene.cor = if (mode == 'fixed_0.01') .01 else NA_real_,
      allow.neg.cor = FALSE, use.ranks = FALSE, trend.var = FALSE, sort = FALSE)
    desired_contrast <- contrast
    z <- results[scenario == 'window_2000' & contrast == desired_contrast & correlation_mode == mode]
    z <- z[match(rownames(recomputed), pathway)]
    stopifnot(max(abs(recomputed$PValue - z$PValue)) < 1e-12,
      identical(recomputed$Direction, z$Direction), identical(as.integer(recomputed$NGenes), z$NGenes))
    checks <- checks + nrow(z)
  }
}
jsonlite::write_json(list(status = 'passed', independent_design_reconstruction = TRUE,
  CAMERA_gene_ID_index_and_explicit_weight_checks = checks),
  'provenance/chromatin_rna/camera_validation.json', auto_unbox = TRUE, pretty = TRUE)
cat('Independent design and CAMERA validation passed:', checks, 'tests\n')
