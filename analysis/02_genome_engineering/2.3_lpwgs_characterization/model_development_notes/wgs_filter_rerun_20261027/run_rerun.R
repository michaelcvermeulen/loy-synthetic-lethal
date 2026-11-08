# Run from A:/Dropbox/LOY_paralogs/projects. All analysis outputs are staged here.
.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5', .libPaths()))
suppressPackageStartupMessages({library(QDNAseq); library(Biobase)})
future::plan(future::sequential)
options('QDNAseq::verbose'=TRUE)
source('WGS/code/preserve_t2t_bin_mask.R')
out <- 'SEPT_11_2026_THESIS/research/model_development/wgs_filter_rerun_20261027'
dir.create(file.path(out,'inputs'),recursive=TRUE,showWarnings=FALSE)
bins <- readRDS('WGS/t2t_bins_100kb.rds')
annotation <- pData(bins)
mask <- annotation$use
stopifnot(sum(mask)==22682, sum(mask & annotation$chromosome=='Y')==58)
write.csv(data.frame(bin_id=rownames(annotation),annotation),file.path(out,'bin_annotation.csv'),row.names=FALSE)
protected <- c(list.files('SEPT_11_2026_THESIS/draftv2',recursive=TRUE,full.names=TRUE),
  'SEPT_11_2026_THESIS/Chapter3_figure_captions.md',
  'SEPT_11_2026_THESIS/figures/chapter3/selected/selected_figure_index.tsv',
  list.files('SEPT_11_2026_THESIS/figures/chapter3/selected', pattern='regional_y|chrY_atlas',recursive=TRUE,full.names=TRUE),
  paste0('WGS/',c('HT1080','A549','LNCaP'),'_QDNAseq_cn.rds'),'WGS/t2t_bins_100kb.rds')
protected <- unique(protected[!file.info(protected)$isdir])
before_hash <- tools::md5sum(protected)
write.csv(data.frame(path=names(before_hash),md5=unname(before_hash)),file.path(out,'protected_files_before.csv'),row.names=FALSE)
qc <- list(); reproduction <- list(); input_manifest <- list()
checkpoint <- function(x, expected, cohort, mode, stage) {
  stopifnot(identical(fData(x)$use,expected))
  qc[[length(qc)+1]] <<- data.frame(cohort,mode,stage,usable_bins=sum(fData(x)$use),usable_Y_bins=sum(fData(x)$use&fData(x)$chromosome=='Y'))
  cat(format(Sys.time()),cohort,mode,stage,'mask verified\n')
  flush.console()
  x
}
export_cn <- function(x, path) {
  a <- fData(x); cn <- assayDataElement(x,'copynumber'); seg <- assayDataElement(x,'segmented')
  tab <- do.call(rbind,lapply(seq_len(ncol(x)),function(j) data.frame(
    bin_id=rownames(a), sample=sampleNames(x)[j], chromosome=a$chromosome,
    start=a$start,end=a$end,use=a$use,copynumber=cn[,j],segmented=seg[,j])))
  write.csv(tab,path,row.names=FALSE,na='')
}
for (cohort in c('HT1080','A549','LNCaP')) {
  display <- if(cohort=='LNCaP') 'DU145' else cohort
  cat('\nSTART',display,format(Sys.time()),'\n')
  dest <- file.path(out,display); dir.create(dest,showWarnings=FALSE)
  cache <- sub('.rdx$','',list.files(paste0('WGS/code/coverage_report_',cohort,'_QDNAseq_cache/html'),pattern='^run_qdnaseq.*.rdx$',full.names=TRUE))
  stopifnot(length(cache)==1)
  env <- new.env(); lazyLoad(cache,envir=env)
  rc <- env$rc
  stopifnot(identical(fData(rc),annotation))
  saveRDS(rc,file.path(out,'inputs',paste0(display,'_cached_raw_counts.rds')))
  input_files <- c(paste0(cache,c('.rdx','.rdb')),paste0('WGS/',cohort,'_QDNAseq_cn.rds'))
  hh <- tools::md5sum(input_files)
  input_manifest[[display]] <- data.frame(display, path=names(hh),md5=unname(hh))
  sample_qc <- data.frame(sample=sampleNames(rc),pData(rc),excluded_low_depth=rc$total.reads < 5e6)
  write.csv(sample_qc,file.path(dest,'input_sample_qc.csv'),row.names=FALSE)
  rc <- rc[,rc$total.reads >= 5e6]
  old <- readRDS(paste0('WGS/',cohort,'_QDNAseq_cn.rds'))
  stopifnot(identical(sampleNames(rc),sampleNames(old)))
  export_cn(old,file.path(dest,'historical_bins.csv'))
  for (mode in c('corrected','baseline_reproduction')) {
    set.seed(20261027)
    expected <- if(mode=='corrected') mask else rep(TRUE,length(mask))
    x <- if(mode=='corrected') preserve_t2t_bin_mask(rc,bins) else
      applyFilters(rc,residual=FALSE,blacklist=FALSE,mappability=FALSE,bases=FALSE,chromosomes=NA)
    x <- checkpoint(x,expected,display,mode,'filtered')
    write.csv(data.frame(sample=sampleNames(x),pData(x)),file.path(dest,paste0(mode,'_sample_qc.csv')),row.names=FALSE)
    x <- checkpoint(estimateCorrection(x),expected,display,mode,'correction_fit')
    saveRDS(x,file.path(dest,paste0(mode,'_correction_fit.rds')))
    x <- checkpoint(correctBins(x),expected,display,mode,'corrected_counts')
    x <- checkpoint(normalizeBins(x),expected,display,mode,'normalized')
    x <- checkpoint(smoothOutlierBins(x),expected,display,mode,'smoothed')
    x <- checkpoint(segmentBins(x,transformFun='sqrt'),expected,display,mode,'segmented')
    x <- checkpoint(normalizeSegmentedBins(x),expected,display,mode,'segment_normalized')
    call_error <- NULL
    x <- tryCatch(callBins(x),error=function(e){call_error <<- conditionMessage(e); x})
    x <- checkpoint(x,expected,display,mode,'final')
    writeLines(if(is.null(call_error)) 'CGHcall succeeded' else paste('CGHcall optional step failed:',call_error),file.path(dest,paste0(mode,'_call_status.txt')))
    saveRDS(x,file.path(dest,paste0(mode,'_QDNAseq_cn.rds')))
    export_cn(x,file.path(dest,paste0(mode,'_bins.csv')))
    if(mode=='baseline_reproduction') {
      for(assay in c('copynumber','segmented')) {
        a <- assayDataElement(old,assay); b <- assayDataElement(x,assay)
        reproduction[[length(reproduction)+1]] <- data.frame(cohort=display,assay,
          same_NA_pattern=identical(is.na(a),is.na(b)),
          max_absolute_difference=max(abs(a-b),na.rm=TRUE),
          median_absolute_difference=median(abs(a-b),na.rm=TRUE),
          correlation=cor(as.vector(a),as.vector(b),use='complete.obs'))
      }
    }
    write.csv(do.call(rbind,qc),file.path(out,'mask_checkpoints.csv'),row.names=FALSE)
  }
}
write.csv(do.call(rbind,reproduction),file.path(out,'baseline_reproduction_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,input_manifest),file.path(out,'input_manifest.csv'),row.names=FALSE)
after_hash <- tools::md5sum(protected)
stopifnot(identical(before_hash,after_hash))
writeLines(paste('Verified unchanged:',length(protected),'protected files'),file.path(out,'protected_files_check.txt'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('\nALL THREE COHORTS COMPLETE',format(Sys.time()),'\n')
