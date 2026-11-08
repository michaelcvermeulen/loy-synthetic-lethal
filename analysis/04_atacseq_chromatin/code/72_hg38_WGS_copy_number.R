options(warn=1)
.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(QDNAseq);library(Biobase);library(Rsamtools);library(data.table);library(jsonlite);library(digest)})
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
ROOT<-normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
OUT<-file.path(ROOT,'integration/RNA_interaction_CN');dir.create(file.path(OUT,'figures'),showWarnings=FALSE)
e<-new.env();annotation_path<-file.path(OUT,'sources/hg38.100kbp.SR50.rda');load(annotation_path,envir=e)
bins<-e[['hg38.100kbp.SR50']]
stopifnot(attr(bins,'QDNAseq')$build=='hg38')
samples<-c('HT1080-sgCTRL-D0','HT1080-LOY-D0')
bams<-file.path('A:/Dropbox/LOY_paralogs/projects/WGS/hg38_bam',paste0(samples,'.bam'))
for(p in bams){h<-scanBamHeader(p)[[1]];sizes<-h$targets;stopifnot(sizes['chr1']==248956422,sizes['chr19']==58617616,sizes['chrX']==156040895)}
rcpath<-file.path(OUT,'WGS_hg38_100kb_readcounts.rds')
if(file.exists(rcpath)){rc<-readRDS(rcpath)}else{
 rc<-binReadCounts(bins,bamfiles=bams,bamnames=samples,chunkSize=1000000L,cache=FALSE,
   isPaired=TRUE,isProperPair=TRUE,isFirstMateRead=TRUE,isSecondaryAlignment=FALSE,
   isDuplicate=FALSE,isNotPassingQualityControls=FALSE,minMapq=37,verbose=FALSE)
 saveRDS(rc,rcpath)
}
stopifnot(all(colSums(assayDataElement(rc,'counts'))>1e6))
# No population reference residual is available for some hg38 bins. The standard
# residual filter drops these along with outlier residuals; report all exclusions.
filtered<-applyFilters(rc,residual=TRUE,blacklist=TRUE,mappability=50,bases=50,chromosomes=c('X','Y'))
stopifnot(sum(fData(filtered)$use)>20000)
corrected<-estimateCorrection(filtered)
cn<-correctBins(corrected);cn<-normalizeBins(cn)
normalized<-cn
cn<-smoothOutlierBins(cn)
set.seed(20260909)
seg<-segmentBins(cn,transformFun='sqrt')
seg<-normalizeSegmentedBins(seg)
saveRDS(list(readcounts=rc,filtered=filtered,corrected=corrected,normalized=normalized,smoothed=cn,segmented=seg),file.path(OUT,'WGS_hg38_QDNAseq_objects.rds'))
tab<-as.data.table(fData(seg),keep.rownames='bin_id')
tab[,`:=`(chrom=paste0('chr',chromosome),start0=start-1L,end0=end)]
for(i in seq_along(samples)){
 label<-c('CTRL','LOY')[i]
 tab[,(paste0(label,'_raw_first_mate_counts')):=as.numeric(assayDataElement(rc,'counts')[,i])]
 tab[,(paste0(label,'_normalized_ratio')):=as.numeric(assayDataElement(normalized,'copynumber')[,i])]
 tab[,(paste0(label,'_smoothed_ratio')):=as.numeric(assayDataElement(cn,'copynumber')[,i])]
 tab[,(paste0(label,'_segmented_ratio')):=as.numeric(assayDataElement(seg,'segmented')[,i])]
}
tab[,`:=`(CTRL_segment_log2=log2(CTRL_segmented_ratio),LOY_segment_log2=log2(LOY_segmented_ratio),
 delta_bin_log2=log2(LOY_normalized_ratio/CTRL_normalized_ratio),delta_segment_log2=log2(LOY_segmented_ratio/CTRL_segmented_ratio))]
tab[use==FALSE,c('CTRL_segment_log2','LOY_segment_log2','delta_bin_log2','delta_segment_log2'):=list(NA_real_,NA_real_,NA_real_,NA_real_)]
fwrite(tab,file.path(OUT,'WGS_hg38_100kb_bins.csv'))
fwrite(as.data.table(pData(rc),keep.rownames='sample'),file.path(OUT,'WGS_readcount_QC.csv'))
segments<-list()
for(label in c('CTRL','LOY')){
 good<-tab[use&is.finite(get(paste0(label,'_segment_log2')))]
 good[,value:=get(paste0(label,'_segment_log2'))]
 # Split exports at excluded-bin gaps rather than interpolating missing loci.
 good[,block:=cumsum(c(TRUE,chrom[-1]!=chrom[-.N]|start0[-1]!=end0[-.N]|value[-1]!=value[-.N]))]
 segments[[label]]<-good[,.(chrom=chrom[1],start0=min(start0),end0=max(end0),bins=.N,segment_log2=value[1]),by=block][,sample:=label]
}
fwrite(rbindlist(segments),file.path(OUT,'WGS_hg38_segments.csv'))
chrom<-tab[use&is.finite(delta_segment_log2),.(bins=.N,median_bin_delta=median(delta_bin_log2,na.rm=TRUE),median_segment_delta=median(delta_segment_log2),fraction_abs_segment_delta_ge020=mean(abs(delta_segment_log2)>=.2)),by=chrom]
fwrite(chrom,file.path(OUT,'WGS_chromosome_summary.csv'))
pdf(file.path(OUT,'figures/WGS_QDNAseq_diagnostics.pdf'),width=12,height=5)
for(i in 1:2){plot(rc[,i],logTransform=FALSE,main=samples[i]);noisePlot(corrected[,i]);plot(seg[,i],main=samples[i])}
dev.off()
methods<-c('D0 line-level copy-number context from hg38/GRCh38 BAMs HT1080-sgCTRL-D0 and HT1080-LOY-D0. One WGS sample per condition; direct matching to ATAC or RNA cultures is unconfirmed. No replicate-based CN hypothesis tests or sample-specific offsets are applied.',
 'Reference dictionaries verify hg38 lengths (chr1 248956422, chr19 58617616, chrX 156040895). Existing T2T CN objects are not overlapped. Frozen QDNAseq.hg38 v1.2.0 100kb SR50 annotations provide GC, mappability, blacklist, bases and reference residuals. Annotation bins are 1-based inclusive; exported start0/end0 are zero-based half-open.',
 'QDNAseq read counting uses mapped, properly paired first-mate reads, MAPQ >=37, excluding secondary alignments, failed QC and reads flagged duplicate. Counts therefore represent accepted first-mate alignments; no new duplicate marking is performed. Annotation mappability is a reference SR50 metric rather than a sample-specific measurement.',
 'Filter using QDNAseq reference-residual and blacklist filters, mappability >=50%, known bases >=50%; exclude X and Y from this autosomal dosage assessment and normalization. Excluded bins remain missing, not neutral. GC/mappability LOESS correction, median normalization, outlier smoothing, CBS segmentation with sqrt transform and segmented normalization use QDNAseq defaults; seed 20260909.',
 'Values are relative copy-number ratios, not integer copy numbers or ploidy estimates. Comparison is log2(LOY segmented ratio / CTRL segmented ratio). Absolute global ploidy differences can be removed by median normalization. For descriptive stratification use |delta| >=0.20, with 0.10/0.30 sensitivity thresholds; these are effect-size flags, not significance calls. Resolution is 100kb and does not exclude focal changes below that scale.',
 'Genes and promoter peaks must overlap usable bins; no propagation across filtered gaps. One D0 sample per condition cannot exclude later-acquired copy-number changes or establish causality.',
 'QDNAseq workflow: https://bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf',
 'hg38 annotation: https://github.com/asntech/QDNAseq.hg38/tree/v1.2.0 ; https://doi.org/10.5281/zenodo.17948766')
writeLines(methods,file.path(OUT,'WGS_CN_METHODS.txt'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'WGS_sessionInfo.txt'))
inputs<-c(bams,paste0(bams,'.bai'),annotation_path,script)
fwrite(data.table(path=inputs,bytes=file.info(inputs)$size,sha256=vapply(inputs,digest,character(1),algo='sha256',file=TRUE)),file.path(OUT,'WGS_source_hashes.csv'))
summary<-list(status='complete',assembly='hg38',bin_size=100000,all_bins=nrow(tab),usable_autosomal_bins=sum(tab$use),samples=samples,assay_matching='D0 line-level; direct culture matching unconfirmed',absolute_copy_number=FALSE)
write_json(summary,file.path(OUT,'WGS_model_summary.json'),pretty=TRUE,auto_unbox=TRUE)
print(summary);print(chrom)
