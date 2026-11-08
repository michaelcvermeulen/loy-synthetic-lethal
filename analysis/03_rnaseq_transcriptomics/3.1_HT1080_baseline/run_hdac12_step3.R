# Independent-study benchmark; original HT1080 expression and enrichment objects are read-only inputs.
options(warn = 1, stringsAsFactors = FALSE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/', mustWork = TRUE)))
invisible(lapply(sort(list.files('R', pattern = '\\.R$', full.names = TRUE)), source))
data.table::setDTthreads(2)
for (manifest in c('resources/manifest.csv', 'resources/hdac12_external/GSE193885/manifest.csv')) {
  r <- as.data.frame(data.table::fread(manifest))
  assert(all(vapply(r$path, hash_file, character(1)) == r$sha256), paste('Frozen resource mismatch:', manifest))
}
paths <- c('outputs/models/01_ingest.rds', 'outputs/models/02_expression.rds',
  'outputs/models/03_pathways_regulators.rds', 'outputs/models/camera_voom_model.rds',
  'outputs/models/04_HDAC12_gene_drivers.rds', 'config/hdac12_external.yml',
  'outputs/tables/HDACEXT_sample_manifest.csv', 'outputs/tables/HDACEXT_contrast_manifest.csv',
  'resources/hdac12_external/GSE193885/GSE193885_counts_raw.tsv.gz')
deps <- data.frame(path = paths, sha256 = unname(vapply(paths, hash_file, character(1))))
hash_path <- 'provenance/HDACEXT_frozen_input_hashes.csv'
if (file.exists(hash_path)) {
  old <- as.data.frame(data.table::fread(hash_path))
  assert(identical(old$path, deps$path) && identical(old$sha256, deps$sha256), 'Step-3 frozen inputs changed')
} else data.table::fwrite(deps, hash_path)
x <- readRDS(paths[1]); e <- readRDS(paths[2]); p <- readRDS(paths[3])
external <- stage('04_HDACEXT_source_models', function() hdac_external_source(e))
hdac_external_ms_audit(external)
if ('--source-only' %in% commandArgs(TRUE)) quit(status = 0)
comparison <- stage('04_HDACEXT_comparison', function() hdac_external_compare(external, x, e, p))
if ('--analysis-only' %in% commandArgs(TRUE)) quit(status = 0)
existing <- as.data.frame(data.table::fread('outputs/figures/figure_manifest.csv'))
prior_path <- 'provenance/HDACEXT_prior_figure_hashes.csv'
if (!file.exists(prior_path)) {
  prior <- existing[!existing$figure %in% c('Fig09_HDAC_external_comparison', 'SuppFig06_HDAC_external_QC_sensitivity'), ]
  data.table::fwrite(unique(prior[, c('figure','pdf','svg','png','pdf_sha256','svg_sha256','png_sha256')]), prior_path)
}
capture <- ragg::agg_capture(width = 183, height = 230, units = 'mm', res = 150)
theme_set(nature_theme())
new <- hdac_external_figures(comparison, external, e)
grDevices::dev.off()
ids <- vapply(new, function(m) m$figure[1], character(1))
manifest <- as.data.frame(data.table::rbindlist(c(list(existing[!existing$figure %in% ids, ]), new), fill = TRUE))
data.table::fwrite(manifest, 'outputs/figures/figure_manifest.csv')
saveRDS(manifest, 'outputs/models/05_figures.rds')
captions <- vapply(unique(manifest$caption), function(f) paste(readLines(f, warn = FALSE), collapse = '\n'), character(1))
writeLines(paste(captions, collapse = '\n\n'), 'outputs/captions/figure_captions.md', useBytes = TRUE)
d <- readRDS('outputs/models/04_dosage_repeatability.rds'); cfg <- yaml::read_yaml('config/analysis.yml')
finalize_run(x, e, p, d, cfg, normalizePath('../../..', winslash = '/'))
