.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5', .libPaths()))
suppressPackageStartupMessages({library(QDNAseq); library(Biobase)})
future::plan(future::sequential)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1, args %in% c('A549-CTRL2','A549-LOY'))
sample <- args[1]
out <- 'SEPT_11_2026_THESIS/research/model_development/a549_ypar_realign_20261027'
bam <- file.path('WGS/reanalysis/Y_PAR_masked_20261027/alignment',paste0(sample,'.YPARmasked.bam'))
target <- file.path(out,paste0(sample,'_raw_counts.rds'))
stopifnot(file.exists(bam), !file.exists(target))
bins <- readRDS('WGS/t2t_bins_100kb.rds')
cat('Counting',sample,'using the historical QDNAseq read-count settings\n')
rc <- binReadCounts(bins, bamfiles=bam, bamnames=sample, isPaired=TRUE,
                   minMapq=37, cache=FALSE, verbose=TRUE)
stopifnot(identical(fData(rc),pData(bins)), rc$total.reads>=5e6)
saveRDS(rc,target)
write.csv(data.frame(sample=sampleNames(rc),pData(rc)),
          file.path(out,paste0(sample,'_count_qc.csv')),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'count_sessionInfo.txt'))
cat('Saved',target,'\n')
