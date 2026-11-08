.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(QDNAseq);library(Biobase)})
future::plan(future::sequential)
options('QDNAseq::verbose'=TRUE)
out <- 'SEPT_11_2026_THESIS/research/model_development/a549_par_track_relaxed_20261027'
base <- 'SEPT_11_2026_THESIS/research/model_development/a549_ypar_realign_20261027'
rc <- readRDS(file.path(base,'A549_masked_raw_counts.rds'))
annotation <- fData(rc)
par <- (annotation$chromosome=='Y' & (annotation$start<=2458320 | annotation$end>62122809)) |
       (annotation$chromosome=='X' & (annotation$start<=2394410 | annotation$end>153925834))
# Original mappability annotations describe an unmasked reference and are not
# valid for estimating the improved X-PAR mappability after Y-PAR masking.
# X-PAR is evaluated by observed primary-read coverage outside QDNAseq.
stopifnot(!any(annotation$use & par))
specs <- data.frame(name=c('map30_difficulty50','map30_difficulty80','map30_no_difficulty'),
                    min_map=c(30,30,30),max_difficulty=c(50,80,100.01))
checks <- list()
for(i in seq_len(nrow(specs))) {
  spec <- specs[i,]
  target <- file.path(out,spec$name); dir.create(target,showWarnings=FALSE)
  mask <- annotation$bases>=60 & annotation$mappability>=spec$min_map &
          annotation$blacklist<spec$max_difficulty & !par
  x <- rc
  fData(x)$use <- mask
  x$used.reads <- colSums(assayDataElement(x,'counts')[mask,,drop=FALSE])
  x$expected.variance <- sum(mask)/x$used.reads
  if('paired.ends' %in% colnames(pData(x))) x$expected.variance <- x$expected.variance*ifelse(x$paired.ends,2,1)
  cat('START',spec$name,'usable bins',sum(mask),'\n');flush.console()
  set.seed(20261027)
  x <- estimateCorrection(x)
  saveRDS(x,file.path(target,'correction_fit.rds'))
  x <- correctBins(x)
  x <- normalizeBins(x)
  x <- smoothOutlierBins(x)
  x <- segmentBins(x,transformFun='sqrt')
  x <- normalizeSegmentedBins(x)
  stopifnot(identical(fData(x)$use,mask))
  saveRDS(x,file.path(target,'QDNAseq_cn.rds'))
  aa <- fData(x); cn <- assayDataElement(x,'copynumber'); seg <- assayDataElement(x,'segmented')
  tab <- do.call(rbind,lapply(seq_len(ncol(x)),function(j) data.frame(
    bin_id=rownames(aa),sample=sampleNames(x)[j],chromosome=aa$chromosome,start=aa$start,end=aa$end,
    use=aa$use,copynumber=cn[,j],segmented=seg[,j])))
  write.csv(tab,file.path(target,'bins.csv'),row.names=FALSE,na='')
  checks[[i]] <- data.frame(filter=spec$name,nonN_min=60,mappability_min=spec$min_map,
    difficulty_max_exclusive=spec$max_difficulty,usable_bins=sum(mask),mask_preserved=TRUE)
  cat('COMPLETE',spec$name,'\n');flush.console()
}
write.csv(do.call(rbind,checks),file.path(out,'fit_validation.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
