.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(QDNAseq);library(Biobase)})
out <- 'SEPT_11_2026_THESIS/research/model_development/wgs_filter_rerun_20261027'
stage_checks <- list()
showMethods('normalizeSegmentedBins',includeDefs=TRUE)
showMethods('smoothOutlierBins',includeDefs=TRUE)
for(cohort in c('HT1080','A549','LNCaP')) {
  display <- if(cohort=='LNCaP') 'DU145' else cohort
  p <- file.path(out,display,'baseline_reproduction_QDNAseq_cn.rds')
  if(!file.exists(p)) next
  old <- readRDS(paste0('WGS/',cohort,'_QDNAseq_cn.rds'))
  new <- readRDS(p)
  cat('\n',display,'\n')
  for(s in sampleNames(old)) {
    a <- assayDataElement(old,'copynumber')[,s]; b <- assayDataElement(new,'copynumber')[,s]
    keep <- a>0.05 & b>0.05
    cat(s,'maxDiff',max(abs(a-b),na.rm=TRUE),'cor',cor(a,b,use='complete.obs'),'ratio quantiles',quantile(b[keep]/a[keep],na.rm=TRUE),'\n')
  }
  cache <- sub('.rdx$','',list.files(paste0('WGS/code/coverage_report_',cohort,'_QDNAseq_cache/html'),pattern='^correct_segment.*.rdx$',full.names=TRUE))
  e <- new.env();lazyLoad(cache,envir=e)
  for(n in c('rc_f','cn_corr','cn_norm','cn_smooth','cn_seg','cn')) {
    x <- get(n,envir=e)
    cat(n,dim(x),'\n')
    if(n=='rc_f') {
      fit <- readRDS(file.path(out,display,'baseline_reproduction_correction_fit.rds'))
      for(a in c('counts','fit')) cat(a,'diff:',max(abs(assayDataElement(x,a)-assayDataElement(fit,a)),na.rm=TRUE),'\n')
      smooth <- smoothOutlierBins(normalizeBins(correctBins(fit)))
    } else cat('CN diff vs baseline final:',max(abs(assayDataElement(x,'copynumber')-assayDataElement(new,'copynumber')),na.rm=TRUE),'\n')
  }
  for(n in c('counts','fit','smoothed_copynumber')) {
    a <- if(n=='smoothed_copynumber') assayDataElement(e$cn_smooth,'copynumber') else assayDataElement(e$rc_f,n)
    b <- if(n=='smoothed_copynumber') assayDataElement(smooth,'copynumber') else assayDataElement(fit,n)
    error <- max(abs(a-b),na.rm=TRUE)
    stopifnot(error < 1e-10, identical(is.na(a),is.na(b)))
    stage_checks[[length(stage_checks)+1]] <- data.frame(model=display,stage=n,max_absolute_difference=error)
  }
}
write.csv(do.call(rbind,stage_checks),file.path(out,'presegmentation_reproduction_checks.csv'),row.names=FALSE)
