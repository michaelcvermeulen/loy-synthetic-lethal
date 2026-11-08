# Run after inspecting the composed/individual exports and color-vision simulations.
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(normalizePath(script, winslash = '/')))
source('R/00_common.R')
manifest <- as.data.frame(data.table::fread('outputs/figures/figure_manifest.csv'))
new_ids <- c('Fig12_GSEA_PRC2_H3K27me3','Fig13_GSEA_other_chromatin')
# Figures 10/11 were visually reviewed in the preceding analysis, but their review fields were not written.
review_ids <- c('Fig10_PRC2_gene_programs','Fig11_epigenetic_context',new_ids)
for (id in review_ids) assert(file.exists(file.path('provenance/figure_QA',paste0(id,'_simulations.png'))),
  paste('Missing visual-review artifact:', id))
selected <- manifest$figure %in% review_ids
manifest$visual_review[selected] <- 'reviewed_color_exports_and_CVD_simulations; source_labels_and_values_verified; signs_and_positions_readable_in_grayscale'
manifest$reviewed_at[selected] <- format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')
for (kind in c('pdf','svg','png')) manifest[[paste0(kind,'_sha256')]][selected] <-
  vapply(manifest[[kind]][selected],hash_file,character(1))
data.table::fwrite(manifest,'outputs/figures/figure_manifest.csv'); saveRDS(manifest,'outputs/models/05_figures.rds')
individual <- as.data.frame(data.table::fread('outputs/figures/chromatin_GSEA_individual/manifest.csv'))
individual$visual_review <- 'reviewed_individual_layout_and_composite_counterpart; shared_signed_axes_and_source_labels_verified'
data.table::fwrite(individual,'outputs/figures/chromatin_GSEA_individual/manifest.csv')
captions <- vapply(unique(manifest$caption),function(f) paste(readLines(f,warn=FALSE),collapse='\n'),character(1))
writeLines(paste(captions,collapse='\n\n'),'outputs/captions/figure_captions.md',useBytes=TRUE)
rmarkdown::pandoc_convert(normalizePath('outputs/captions/figure_captions.md',winslash='/'),to='docx',
  output=file.path(normalizePath('outputs/captions',winslash='/'),'figure_captions.docx'))
rmarkdown::render('report/chromatin_GSEA_curves.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),quiet=TRUE)
rmarkdown::render('report/analysis.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),quiet=TRUE)
# Update code/report provenance while retaining verified hashes of unchanged upstream resources.
provenance <- as.data.frame(data.table::fread('provenance/input_code_resource_hashes.csv'))
paths <- unique(c(list.files('R',full.names=TRUE),list.files('config',full.names=TRUE),list.files('tests',full.names=TRUE),
  list.files('report',pattern='\\.Rmd$',full.names=TRUE),'run.R','render_report.R','README.md','FIGURES.md',
  'provenance/DECISIONS.md','run_chromatin_GSEA_plots.R','finalize_chromatin_GSEA_plots.R'))
updated <- data.frame(path=paths,sha256=unname(vapply(paths,hash_file,character(1))),bytes=file.info(paths)$size,
  recorded_at=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'))
provenance <- rbind(provenance[!provenance$path %in% paths,],updated)
data.table::fwrite(provenance,'provenance/input_code_resource_hashes.csv')
message('Reviewed GSEA exports, captions, reports and code provenance finalized.')
