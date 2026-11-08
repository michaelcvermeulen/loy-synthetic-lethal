.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5', .libPaths()))
suppressPackageStartupMessages({library(QDNAseq); library(Biobase)})
future::plan(future::sequential)
options('QDNAseq::verbose'=TRUE)
source('WGS/code/preserve_t2t_bin_mask.R')
out <- 'SEPT_11_2026_THESIS/research/model_development/a549_ypar_realign_20261027'
bins <- readRDS('WGS/t2t_bins_100kb.rds')
a <- readRDS(file.path(out,'A549-CTRL2_raw_counts.rds'))
b <- readRDS(file.path(out,'A549-LOY_raw_counts.rds'))
rc <- combine(a,b)
stopifnot(identical(sampleNames(rc),c('A549-CTRL2','A549-LOY')))
saveRDS(rc,file.path(out,'A549_masked_raw_counts.rds'))
old <- readRDS('SEPT_11_2026_THESIS/research/model_development/wgs_filter_rerun_20261027/inputs/A549_cached_raw_counts.rds')
stopifnot(identical(featureNames(old),featureNames(rc)),identical(sampleNames(old),sampleNames(rc)))
autosomes <- fData(rc)$use & fData(rc)$chromosome %in% as.character(1:22)
raw_qc <- do.call(rbind,lapply(seq_len(ncol(rc)),function(j) {
  before <- assayDataElement(old,'counts')[autosomes,j]
  after <- assayDataElement(rc,'counts')[autosomes,j]
  data.frame(sample=sampleNames(rc)[j],autosomal_bins=sum(autosomes),
    raw_count_correlation=cor(before,after),bins_identical=sum(before==after),
    median_absolute_change=median(abs(after-before)),max_absolute_change=max(abs(after-before)))
}))
write.csv(raw_qc,file.path(out,'raw_count_alignment_comparison.csv'),row.names=FALSE)
write.csv(data.frame(bin_id=featureNames(rc),assayDataElement(rc,'counts'),check.names=FALSE),
          file.path(out,'masked_raw_counts.csv'),row.names=FALSE)
mask <- pData(bins)$use
# Only Y-PAR reference bases changed; every Y-PAR bin already fails this
# frozen annotation's quality mask. All retained reference bases are unchanged.
ypar <- pData(bins)$chromosome=='Y' &
  (pData(bins)$start<=2458320 | pData(bins)$end>62122809)
stopifnot(!any(mask & ypar), sum(mask)==22682)
checkpoints <- list()
check <- function(x,stage) {
  stopifnot(identical(fData(x)$use,mask))
  checkpoints[[length(checkpoints)+1]] <<- data.frame(stage,usable=sum(mask),usable_Y=sum(mask&fData(x)$chromosome=='Y'))
  x
}
set.seed(20261027)
x <- check(preserve_t2t_bin_mask(rc,bins),'filtered')
x <- check(estimateCorrection(x),'correction_fit')
saveRDS(x,file.path(out,'A549_masked_correction_fit.rds'))
x <- check(correctBins(x),'corrected')
x <- check(normalizeBins(x),'normalized')
x <- check(smoothOutlierBins(x),'smoothed')
x <- check(segmentBins(x,transformFun='sqrt'),'segmented')
x <- check(normalizeSegmentedBins(x),'segment_normalized')
# No copy-number calling is needed for the descriptive regional atlas.
saveRDS(x,file.path(out,'A549_masked_QDNAseq_cn.rds'))
aa <- fData(x); cn <- assayDataElement(x,'copynumber'); seg <- assayDataElement(x,'segmented')
tab <- do.call(rbind,lapply(seq_len(ncol(x)),function(j) data.frame(
  bin_id=rownames(aa),sample=sampleNames(x)[j],chromosome=aa$chromosome,
  start=aa$start,end=aa$end,use=aa$use,copynumber=cn[,j],segmented=seg[,j])))
write.csv(tab,file.path(out,'masked_bins.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,checkpoints),file.path(out,'mask_checkpoints.csv'),row.names=FALSE)
write.csv(data.frame(sample=sampleNames(x),pData(x)),file.path(out,'sample_qc.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'fit_sessionInfo.txt'))
cat('New A549 matched-pair QDNAseq analysis complete\n')
