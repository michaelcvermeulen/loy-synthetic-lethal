# Create running-enrichment figures from completed models; original tests remain unchanged.
options(warn = 1, stringsAsFactors = FALSE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/')))
invisible(lapply(sort(list.files('R', pattern = '\\.R$', full.names = TRUE)), source))
data.table::setDTthreads(2)
cfg <- yaml::read_yaml('config/chromatin_GSEA_curves.yml'); set.seed(cfg$seed)
paths <- c('outputs/models/02_expression.rds', 'outputs/models/03_pathways_regulators.rds',
  'config/chromatin_GSEA_curves.yml')
deps <- data.frame(path = paths, sha256 = unname(vapply(paths, hash_file, character(1))))
dep_path <- 'provenance/GSEA_CURVES_frozen_inputs.csv'
if (file.exists(dep_path)) assert(identical(as.data.frame(data.table::fread(dep_path)), deps),
  'GSEA curve inputs changed') else data.table::fwrite(deps, dep_path)
manifest <- as.data.frame(data.table::fread('outputs/figures/figure_manifest.csv'))
prior_path <- 'provenance/GSEA_CURVES_prior_figure_hashes.csv'
if (!file.exists(prior_path)) {
  prior <- unique(manifest[!manifest$figure %in% names(cfg$figures), c('figure','pdf','svg','png')])
  for (kind in c('pdf','svg','png')) prior[[paste0(kind,'_sha256')]] <- vapply(prior[[kind]],hash_file,character(1))
  data.table::fwrite(prior, prior_path)
}
prior <- as.data.frame(data.table::fread(prior_path))
for (kind in c('pdf','svg','png')) assert(all(vapply(prior[[kind]],hash_file,character(1)) == prior[[paste0(kind,'_sha256')]]),
  'A prior publication figure changed')
e <- readRDS(paths[1]); p <- readRDS(paths[2])
z <- stage('05_chromatin_GSEA_curves', function() chromatin_gsea_curve_data(p, e, cfg))
if ('--data-only' %in% commandArgs(TRUE)) quit(status = 0)
capture <- ragg::agg_capture(width = 183, height = 240, units = 'mm', res = 150)
theme_set(nature_theme()); new <- chromatin_gsea_figures(z); grDevices::dev.off()
manifest <- as.data.frame(data.table::rbindlist(c(list(manifest[!manifest$figure %in% names(new), ]), new), fill = TRUE))
data.table::fwrite(manifest, 'outputs/figures/figure_manifest.csv'); saveRDS(manifest, 'outputs/models/05_figures.rds')
captions <- vapply(unique(manifest$caption), function(f) paste(readLines(f, warn = FALSE), collapse = '\n'), character(1))
writeLines(paste(captions, collapse = '\n\n'), 'outputs/captions/figure_captions.md', useBytes = TRUE)
if ('--skip-report' %in% commandArgs(TRUE)) quit(status = 0)
rmarkdown::pandoc_convert(normalizePath('outputs/captions/figure_captions.md', winslash = '/'), to = 'docx',
  output = file.path(normalizePath('outputs/captions', winslash = '/'), 'figure_captions.docx'))
rmarkdown::render('report/chromatin_GSEA_curves.Rmd', output_dir = normalizePath('outputs/report', winslash = '/'), quiet = TRUE)
rmarkdown::render('report/analysis.Rmd', output_dir = normalizePath('outputs/report', winslash = '/'), quiet = TRUE)
message('GSEA figures and reports complete; visual/export review is recorded separately.')
