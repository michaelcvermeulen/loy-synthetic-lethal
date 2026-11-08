# Controlled follow-up using hash-verified completed inputs; no global enrichment refit.
options(warn = 1, stringsAsFactors = FALSE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/', mustWork = TRUE)))
invisible(lapply(sort(list.files('R', pattern = '\\.R$', full.names = TRUE)), source))
data.table::setDTthreads(2)
resources <- as.data.frame(data.table::fread('resources/manifest.csv'))
assert(all(vapply(resources$path, hash_file, character(1)) == resources$sha256), 'Frozen root resource hash mismatch')
paths <- c('outputs/models/01_ingest.rds', 'outputs/models/02_expression.rds',
  'outputs/models/03_pathways_regulators.rds', 'outputs/models/04_chromatin_followup.rds',
  'outputs/tables/HDAC12_frozen_membership.csv', 'config/hdac12_followup.yml')
dependencies <- data.frame(path = paths, sha256 = vapply(paths, hash_file, character(1)))
old <- 'provenance/HDAC12_step2_input_hashes.csv'
if (file.exists(old)) {
  frozen <- as.data.frame(data.table::fread(old))
  assert(identical(dependencies$path, frozen$path) && identical(dependencies$sha256, frozen$sha256),
    'Step-2 inputs changed after the first execution')
} else data.table::fwrite(dependencies, old)
x <- readRDS(paths[1]); e <- readRDS(paths[2]); p <- readRDS(paths[3])
z <- stage('04_HDAC12_gene_drivers', function() hdac12_followup(x, e, p))
if ('--analysis-only' %in% commandArgs(TRUE)) quit(status = 0)
capture <- ragg::agg_capture(width = 183, height = 180, units = 'mm', res = 150)
theme_set(nature_theme())
new <- hdac12_figures(z)
grDevices::dev.off()
existing <- as.data.frame(data.table::fread('outputs/figures/figure_manifest.csv'))
ids <- vapply(new, function(m) m$figure[1], character(1))
manifest <- as.data.frame(data.table::rbindlist(c(list(existing[!existing$figure %in% ids, ]), new), fill = TRUE))
data.table::fwrite(manifest, 'outputs/figures/figure_manifest.csv')
saveRDS(manifest, 'outputs/models/05_figures.rds')
captions <- vapply(unique(manifest$caption), function(f) paste(readLines(f, warn = FALSE), collapse = '\n'), character(1))
writeLines(paste(captions, collapse = '\n\n'), 'outputs/captions/figure_captions.md', useBytes = TRUE)
d <- readRDS('outputs/models/04_dosage_repeatability.rds'); cfg <- yaml::read_yaml('config/analysis.yml')
finalize_run(x, e, p, d, cfg, normalizePath('../../..', winslash = '/'))
