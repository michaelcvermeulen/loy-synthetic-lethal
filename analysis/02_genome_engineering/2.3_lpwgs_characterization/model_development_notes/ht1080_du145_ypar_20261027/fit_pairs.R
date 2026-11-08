.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(QDNAseq);library(Biobase)})
future::plan(future::sequential)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2,args[1] %in% c('HT1080','DU145'),args[2] %in% c('original','masked'))
cell <- args[1];mode <- args[2]
root <- 'SEPT_11_2026_THESIS/research/model_development/ht1080_du145_ypar_20261027'
out <- file.path(root,cell,mode);dir.create(out,recursive=TRUE,showWarnings=FALSE)
samples <- if(cell=='HT1080') c('HT1080-sgCTRL-D0','HT1080-LOY-D0') else c('LNCaP-sgCTRL','LNCaP-LOY-2X')
bins <- readRDS('WGS/t2t_bins_100kb.rds')
if(mode=='original') {
  rc <- readRDS(file.path('SEPT_11_2026_THESIS/research/model_development/wgs_filter_rerun_20261027/inputs',paste0(cell,'_cached_raw_counts.rds')))
  rc <- rc[,match(samples,sampleNames(rc))]
} else {
  target <- file.path(out,'raw_counts.rds')
  if(file.exists(target)) {
    rc <- readRDS(target)
  } else {
    pilot <- 'WGS/reanalysis/Y_PAR_masked_20261027'
    stopifnot(all(file.exists(file.path(pilot,'qc',paste0(samples,'.REALIGNMENT_COMPLETE.txt')))))
    bam <- file.path(pilot,'alignment',paste0(samples,'.YPARmasked.bam'))
    rc <- binReadCounts(bins,bamfiles=bam,bamnames=samples,isPaired=TRUE,minMapq=37,cache=FALSE)
  }
}
stopifnot(identical(sampleNames(rc),samples),identical(fData(rc),pData(bins)),all(rc$total.reads>=5e6))
saveRDS(rc,file.path(out,'raw_counts.rds'))
annotation <- fData(rc)
par <- (annotation$chromosome=='Y' & (annotation$start<=2458320 | annotation$end>62122809)) |
       (annotation$chromosome=='X' & (annotation$start<=2394410 | annotation$end>153925834))
mask <- annotation$bases>=60 & annotation$mappability>=30 & !par
x <- rc
fData(x)$use <- mask
x$used.reads <- colSums(assayDataElement(x,'counts')[mask,,drop=FALSE])
x$expected.variance <- sum(mask)/x$used.reads
if('paired.ends' %in% colnames(pData(x))) x$expected.variance <- x$expected.variance*ifelse(x$paired.ends,2,1)
checkpoints <- list()
check <- function(z,stage) {
  stopifnot(identical(fData(z)$use,mask))
  checkpoints[[length(checkpoints)+1]] <<- data.frame(cell,mode,stage,usable_bins=sum(mask))
  z
}
set.seed(20261027)
x <- check(x,'filtered')
x <- check(estimateCorrection(x),'correction_fit')
saveRDS(x,file.path(out,'correction_fit.rds'))
x <- check(correctBins(x),'corrected')
x <- check(normalizeBins(x),'normalized')
x <- check(smoothOutlierBins(x),'smoothed')
x <- check(segmentBins(x,transformFun='sqrt'),'segmented')
x <- check(normalizeSegmentedBins(x),'segment_normalized')
saveRDS(x,file.path(out,'QDNAseq_cn.rds'))
aa <- fData(x); cn <- assayDataElement(x,'copynumber');seg <- assayDataElement(x,'segmented')
tab <- do.call(rbind,lapply(seq_len(ncol(x)),function(j) data.frame(
  bin_id=rownames(aa),sample=sampleNames(x)[j],chromosome=aa$chromosome,start=aa$start,end=aa$end,
  use=aa$use,copynumber=cn[,j],segmented=seg[,j])))
write.csv(tab,file.path(out,'bins.csv'),row.names=FALSE,na='')
write.csv(data.frame(bin_id=featureNames(rc),assayDataElement(rc,'counts'),check.names=FALSE),file.path(out,'raw_counts.csv'),row.names=FALSE)
write.csv(data.frame(sample=sampleNames(rc),pData(rc)),file.path(out,'sample_qc.csv'),row.names=FALSE)
write.csv(do.call(rbind,checkpoints),file.path(out,'mask_checkpoints.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat(cell,mode,'pair fit complete\n')
