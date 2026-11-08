.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(QDNAseq);library(Biobase)})
future::plan(future::sequential)
out <- 'SEPT_11_2026_THESIS/research/model_development/a549_dense_bins_20261027'
inputs <- 'WGS/reanalysis/Y_PAR_masked_20261027/a549_50kb_comparison'
a <- read.csv(file.path(inputs,'annotation_50kb.csv'),stringsAsFactors=FALSE)
rownames(a) <- a$bin_id;a$bin_id <- NULL
par <- (a$chromosome=='Y' & (a$start<=2458320 | a$end>62122809)) |
       (a$chromosome=='X' & (a$start<=2394410 | a$end>153925834))
a$residual <- NA_real_
a$use <- a$bases>=60 & a$mappability>=30 & !par
original <- readRDS('WGS/t2t_bins_100kb.rds')
bins <- new('AnnotatedDataFrame',data=a,varMetadata=varMetadata(original))
attr(bins,'QDNAseq') <- attr(original,'QDNAseq')
saveRDS(bins,file.path(out,'bins_50kb.rds'))
target <- file.path(out,'A549_50kb_raw_counts.rds')
if(file.exists(target)) {
  rc <- readRDS(target)
  stopifnot(identical(fData(rc),a))
} else {
  samples <- c('A549-CTRL2','A549-LOY')
  bam <- file.path('WGS/reanalysis/Y_PAR_masked_20261027/alignment',paste0(samples,'.YPARmasked.bam'))
  rc <- binReadCounts(bins,bamfiles=bam,bamnames=samples,isPaired=TRUE,minMapq=37,cache=FALSE)
  saveRDS(rc,target)
}
stopifnot(identical(fData(rc)$use,a$use),all(rc$total.reads>=5e6))
mask <- a$use
x <- rc
x$used.reads <- colSums(assayDataElement(x,'counts')[mask,,drop=FALSE])
x$expected.variance <- sum(mask)/x$used.reads
set.seed(20261027)
x <- estimateCorrection(x)
saveRDS(x,file.path(out,'correction_fit_50kb.rds'))
x <- correctBins(x)
x <- normalizeBins(x)
x <- smoothOutlierBins(x)
x <- segmentBins(x,transformFun='sqrt')
x <- normalizeSegmentedBins(x)
stopifnot(identical(fData(x)$use,mask))
saveRDS(x,file.path(out,'QDNAseq_50kb.rds'))
aa <- fData(x);cn <- assayDataElement(x,'copynumber');seg <- assayDataElement(x,'segmented')
tab <- do.call(rbind,lapply(seq_len(ncol(x)),function(j) data.frame(
  bin_id=rownames(aa),sample=sampleNames(x)[j],chromosome=aa$chromosome,start=aa$start,end=aa$end,
  use=aa$use,copynumber=cn[,j],segmented=seg[,j])))
write.csv(tab,file.path(out,'bins_50kb.csv'),row.names=FALSE,na='')
write.csv(data.frame(bin_id=featureNames(rc),assayDataElement(rc,'counts'),check.names=FALSE),file.path(out,'raw_counts_50kb.csv'),row.names=FALSE)
write.csv(data.frame(sample=sampleNames(rc),pData(rc)),file.path(out,'sample_qc_50kb.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('50-kb matched-pair fit complete; mask preserved\n')
