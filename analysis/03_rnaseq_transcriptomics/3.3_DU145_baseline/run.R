options(warn = 1, stringsAsFactors = FALSE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
ROOT <- dirname(normalizePath(script, winslash = '/', mustWork = TRUE))
setwd(ROOT)
REPO <- normalizePath('../../..', winslash = '/', mustWork = TRUE)   # .../projects/RNAseq
for (d in c('outputs/tables','outputs/models','outputs/figures','outputs/source_data','outputs/captions','outputs/report','provenance')) dir.create(d, recursive=TRUE, showWarnings=FALSE)
invisible(lapply(sort(list.files('R', pattern='\\.R$', full.names=TRUE)), source))
cfg <- jsonlite::fromJSON('config/analysis.json'); set.seed(cfg$seed)
data.table::setDTthreads(2)
suppressPackageStartupMessages(library(SummarizedExperiment))
frozen <- data.table::fread('resources/manifest.csv')
assert(all(vapply(frozen$path,hash_file,character(1))==frozen$sha256),'Frozen resources changed')
args <- commandArgs(trailingOnly=TRUE)
if ('--enrichment-only' %in% args || '--report-only' %in% args) {
  x <- readRDS('outputs/models/01_ingest.rds'); e <- readRDS('outputs/models/02_expression.rds')
  id <- readRDS('outputs/models/03_identity.rds'); y <- readRDS('outputs/models/04_y_par.rds'); cn <- readRDS('outputs/models/05_copy_number.rds')
  audit <- data.table::fread('outputs/tables/counting_audit.csv')
  assert(all(vapply(audit$source,hash_file,character(1))==audit$sha256),'Raw inputs changed since core run')
  inputs <- data.table::fread('outputs/tables/input_hashes.csv')
  assert(all(vapply(inputs$path,hash_file,character(1))==inputs$sha256),'Matrix, annotation, manifest or comparison inputs changed since core run')
  assert(identical(cfg,readRDS('outputs/models/core_config.rds')),'Configuration changed; rerun core before resuming')
} else {
  x <- stage('01_ingest',function() ingest(cfg,REPO))
  e <- stage('02_expression',function() du145_expression(x,cfg))
  saveRDS(cfg,'outputs/models/core_config.rds')
  id <- stage('03_identity',function() du145_identity(x,e,cfg))
  y <- stage('04_y_par',function() du145_y_par(x,e,cfg))
  cn <- stage('05_copy_number',function() du145_copy_number(x,e,cfg))
  stage('07_figures_core',function() core_figures(x,e,id,y,cn,cfg))
}
if ('--core-only' %in% args) {
  p <- NULL
} else {
  p <- if ('--report-only' %in% args) readRDS('outputs/models/06_pathways.rds') else stage('06_pathways',function() du145_pathways(x,e,cn,cfg))
  stage('08_figures_pathways',function() pathway_figures(p))
}
write_report(x,e,id,y,cn,p,cfg)
source('validate.R')
capture.output(sessionInfo(),file='provenance/sessionInfo.txt')
code <- c('run.R','validate.R','prepare_workspace.py',list.files('R',full.names=TRUE),list.files('config',full.names=TRUE))
data.table::fwrite(data.frame(path=code,sha256=vapply(code,hash_file,character(1))),'provenance/code_hashes.csv')
message('DU145 analysis completed: ',normalizePath('outputs/report/analysis.html',winslash='/'))
