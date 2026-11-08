# Reuse the completed HT1080 models; regenerate only this authorized extension.
options(warn = 1, stringsAsFactors = FALSE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/')))
invisible(lapply(sort(list.files('R', pattern = '\\.R$', full.names = TRUE)), source))
data.table::setDTthreads(2)
for (manifest in c('resources/manifest.csv', 'resources/epigenetic_deep_dive/manifest.csv')) {
  r <- as.data.frame(data.table::fread(manifest))
  assert(all(vapply(r$path, hash_file, character(1)) == r$sha256), paste('Frozen resource mismatch:', manifest))
}
paths <- c('outputs/models/01_ingest.rds','outputs/models/02_expression.rds',
  'outputs/models/03_pathways_regulators.rds','outputs/models/camera_voom_model.rds',
  'outputs/models/04_chromatin_followup.rds','outputs/models/04_HDACEXT_comparison.rds',
  'config/epigenetic_deep_dive.yml','resources/epigenetic_deep_dive/descriptions.csv')
deps <- data.frame(path = paths, sha256 = unname(vapply(paths, hash_file, character(1))))
hash_path <- 'provenance/EPIGEN_frozen_input_hashes.csv'
if (file.exists(hash_path)) {
  old <- as.data.frame(data.table::fread(hash_path))
  assert(identical(old, deps), 'Epigenetic frozen inputs changed')
} else data.table::fwrite(deps, hash_path)
existing <- as.data.frame(data.table::fread('outputs/figures/figure_manifest.csv'))
prior_path <- 'provenance/EPIGEN_prior_figure_hashes.csv'
if (!file.exists(prior_path)) data.table::fwrite(unique(existing[, c('figure','pdf','svg','png',
  'pdf_sha256','svg_sha256','png_sha256')]), prior_path)
x <- readRDS(paths[1]); e <- readRDS(paths[2]); p <- readRDS(paths[3])
if ('--finalize-only' %in% commandArgs(TRUE)) {
  assert(file.exists('outputs/models/04_epigenetic_deep_dive.rds'), 'Run analysis and figures first')
  finalize_run(x, e, p, readRDS('outputs/models/04_dosage_repeatability.rds'),
    yaml::read_yaml('config/analysis.yml'), normalizePath('../../..', winslash = '/'))
  quit(status = 0)
}
if ('--figures-only' %in% commandArgs(TRUE)) z <- readRDS('outputs/models/04_epigenetic_deep_dive.rds') else
  z <- stage('04_epigenetic_deep_dive', function() epigenetic_deep_dive(x, e, p))
if ('--analysis-only' %in% commandArgs(TRUE)) quit(status = 0)
capture <- ragg::agg_capture(width = 183, height = 240, units = 'mm', res = 150)
theme_set(nature_theme()); new <- epigenetic_figures(z, x); grDevices::dev.off()
ids <- vapply(new, function(m) m$figure[1], character(1))
manifest <- as.data.frame(data.table::rbindlist(c(list(existing[!existing$figure %in% ids, ]), new), fill = TRUE))
data.table::fwrite(manifest, 'outputs/figures/figure_manifest.csv'); saveRDS(manifest, 'outputs/models/05_figures.rds')
captions <- vapply(unique(manifest$caption), function(f) paste(readLines(f, warn = FALSE), collapse = '\n'), character(1))
writeLines(paste(captions, collapse = '\n\n'), 'outputs/captions/figure_captions.md', useBytes = TRUE)
if ('--skip-report' %in% commandArgs(TRUE)) quit(status = 0)
d <- readRDS('outputs/models/04_dosage_repeatability.rds'); cfg <- yaml::read_yaml('config/analysis.yml')
finalize_run(x, e, p, d, cfg, normalizePath('../../..', winslash = '/'))
