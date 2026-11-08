.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5', .libPaths()))
suppressPackageStartupMessages({library(QDNAseq); library(Biobase)})
showMethods('applyFilters', includeDefs=TRUE)
showMethods('binsToUse<-', includeDefs=TRUE)
for (cohort in c('HT1080','A549','LNCaP')) {
  cat('\nCOHORT:', cohort, '\n')
  p <- sub('.rdx$', '', list.files(paste0('WGS/code/coverage_report_',cohort,'_QDNAseq_cache/html'), pattern='^run_qdnaseq.*.rdx$', full.names=TRUE))
  stopifnot(length(p)==1)
  e <- new.env(); lazyLoad(p,envir=e)
  print(pData(e$rc)); print(table(fData(e$rc)$use))
  print(head(fData(e$rc)))
  print(attributes(e$rc@experimentData))
}
print(sessionInfo())
