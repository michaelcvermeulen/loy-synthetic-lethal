.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(QDNAseq);library(Biobase)})
source('WGS/code/preserve_t2t_bin_mask.R')
out <- 'SEPT_11_2026_THESIS/research/model_development/wgs_filter_rerun_20261027'
bins <- readRDS('WGS/t2t_bins_100kb.rds')
results <- list()
for(cohort in c('HT1080','A549','LNCaP')) {
  display <- if(cohort=='LNCaP') 'DU145' else cohort
  raw <- readRDS(file.path(out,'inputs',paste0(display,'_cached_raw_counts.rds')))
  broken <- applyFilters(raw,residual=FALSE,blacklist=FALSE,mappability=FALSE,bases=FALSE,chromosomes=NA,verbose=FALSE)
  stopifnot(sum(fData(broken)$use)==31185)
  fixed <- preserve_t2t_bin_mask(broken,bins)
  stopifnot(identical(fData(fixed)$use,fData(raw)$use),
            isTRUE(all.equal(pData(fixed),pData(raw),check.attributes=FALSE)))
  text <- readLines(paste0('WGS/code/coverage_report_',cohort,'_QDNAseq.Rmd'))
  code <- character(); inside <- FALSE
  for(line in text) {
    if(grepl('^```\\{r',line)) {inside <- TRUE;next}
    if(grepl('^```',line)) {inside <- FALSE;next}
    if(inside) code <- c(code,line)
  }
  parse(text=code)
  final <- readRDS(file.path(out,display,'corrected_QDNAseq_cn.rds'))
  stopifnot(all(is.na(assayDataElement(final,'copynumber')[!fData(final)$use,,drop=FALSE])),
            all(is.na(assayDataElement(final,'segmented')[!fData(final)$use,,drop=FALSE])))
  results[[display]] <- data.frame(model=display,restores_mask_and_read_metadata=TRUE,notebook_parses=TRUE,excluded_bins_are_NA=TRUE)
}
write.csv(do.call(rbind,results),file.path(out,'fix_validation.csv'),row.names=FALSE)
print(do.call(rbind,results))
