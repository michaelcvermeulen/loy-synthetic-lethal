.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(QDNAseq);library(Rsamtools);library(Biobase);library(data.table)})
root<-'A:/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS'
out<-file.path(root,'integration/RNA_interaction_CN')
e<-new.env();load(file.path(out,'sources/hg38.100kbp.SR50.rda'),envir=e)
print(ls(e));bins<-e[[ls(e)[1]]];print(bins);print(head(pData(bins)))
print(summary(pData(bins)));print(attr(bins,'QDNAseq')[c('build','date','author','md5')])
for(s in c('HT1080-sgCTRL-D0','HT1080-LOY-D0')) {
 path<-file.path('A:/Dropbox/LOY_paralogs/projects/WGS/hg38_bam',paste0(s,'.bam'))
 h<-scanBamHeader(path)[[1]];print(head(h$targets));print(h$text[c('@HD','@PG','@RG')])
 fwrite(data.table(chrom=names(h$targets),length=as.numeric(h$targets)),file.path(out,paste0(s,'_sequence_dictionary.csv')))
}
for(fn in c('binReadCounts','applyFilters','normalizeBins','segmentBins','normalizeSegmentedBins'))writeLines(capture.output(get(fn,asNamespace('QDNAseq'))),file.path(out,paste0('QDNAseq_',fn,'_API.txt')))
