options(warn=1)
.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(DESeq2);library(data.table);library(QDNAseq);library(Rsamtools);library(Biobase);library(jsonlite)})
root<-'A:/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS';out<-file.path(root,'integration/RNA_interaction_CN')
dds<-readRDS(file.path(out,'RNA_interaction_dds.rds'));r<-fread(file.path(out,'RNA_interaction_all_genes.csv'))
stopifnot(identical(r$gene_id,rownames(dds)),all(mcols(dds)$betaConv))
mm<-model.matrix(design(dds),as.data.frame(colData(dds)))
stopifnot(identical(as.numeric(mm[,4]),as.numeric(dds$batch=='MAY_2026'&dds$treatment=='LOY')))
stopifnot(all(table(dds$batch,dds$treatment)==3),qr(mm)$rank==4)
coefname<-'batchMAY_2026.treatmentLOY'
z<-results(dds,name=coefname,independentFiltering=FALSE)
stopifnot(isTRUE(all.equal(r$interaction_log2FoldChange,as.numeric(z$log2FoldChange),tolerance=1e-10)))
stopifnot(isTRUE(all.equal(r$interaction_padj,p.adjust(r$interaction_pvalue,'BH'),tolerance=1e-10)))
stopifnot(max(abs(r$May_log2FoldChange-r$August_log2FoldChange-r$interaction_log2FoldChange))<1e-8)
finite<-is.finite(r$interaction_pvalue)
stopifnot(max(abs(r$interaction_pvalue[finite]-2*pnorm(-abs(r$interaction_stat[finite]))))<1e-10)
stopifnot(max(abs(r$interaction_stat-r$interaction_log2FoldChange/r$interaction_lfcSE),na.rm=TRUE)<1e-8)
# Reparameterize the already fitted dispersions to verify the contrast sign and SE.
reverse<-dds;reverse$batch<-relevel(reverse$batch,ref='MAY_2026')
reverse<-nbinomWaldTest(reverse,betaPrior=FALSE)
zr<-results(reverse,name='batchAUG_2025.treatmentLOY',independentFiltering=FALSE)
max_delta<-max(abs(z$log2FoldChange+zr$log2FoldChange),na.rm=TRUE)
max_se<-max(abs(z$lfcSE-zr$lfcSE),na.rm=TRUE)
fwrite(data.table(gene_id=rownames(dds),symbol=r$gene_name,baseMean=z$baseMean,forward=z$log2FoldChange,reverse=zr$log2FoldChange,difference=z$log2FoldChange+zr$log2FoldChange,SE_forward=z$lfcSE,SE_reverse=zr$lfcSE,FDR_forward=z$padj,FDR_reverse=zr$padj),file.path(out,'RNA_reparameterization_audit.csv'))
cat('Reparameterization max coefficient difference:',max_delta,'SE difference:',max_se,'\n')
# An all-zero design cell gives a boundary coefficient for LENG1. Flag its
# unstable magnitude explicitly and use an LRT sensitivity to assess testing.
cell<-interaction(dds$batch,dds$treatment)
zero_cell<-apply(sapply(levels(cell),function(g)rowSums(counts(dds)[,cell==g,drop=FALSE])==0),1,any)
unstable<-abs(z$log2FoldChange+zr$log2FoldChange)>=.01
stopifnot(all(zero_cell[unstable]),max_se<.01)
lrt<-nbinomLRT(dds,reduced=~batch+treatment)
lr<-results(lrt,independentFiltering=FALSE,alpha=.05)
fwrite(data.table(gene_id=rownames(dds),symbol=r$gene_name,LRT_stat=lr$stat,LRT_pvalue=lr$pvalue,LRT_FDR=lr$padj,Wald_FDR=z$padj,zero_count_design_cell=zero_cell,unstable_Wald_magnitude=unstable),file.path(out,'RNA_interaction_LRT_sensitivity.csv'))

objs<-readRDS(file.path(out,'WGS_hg38_QDNAseq_objects.rds'));bins<-fread(file.path(out,'WGS_hg38_100kb_bins.csv'))
stopifnot(identical(bins$bin_id,featureNames(objs$segmented)),all(bins$start0==bins$start-1),all(bins$end0==bins$end))
stopifnot(all(is.na(bins[use==FALSE,delta_segment_log2])),!any(bins[use==TRUE,chrom] %in% c('chrX','chrY')))
seg<-assayDataElement(objs$segmented,'segmented')
stopifnot(max(abs(bins[use==TRUE,delta_segment_log2]-log2(seg[bins$use,2]/seg[bins$use,1])),na.rm=TRUE)<1e-10)
# Independently scan one well-mapped 100kb bin per autosome and each BAM. Count
# first-mate start positions to match the binning unit (not read-overlap counts).
sel<-bins[use==TRUE][, .SD[which.max(mappability)],by=chrom]
flags<-scanBamFlag(isPaired=TRUE,isProperPair=TRUE,isUnmappedQuery=FALSE,isFirstMateRead=TRUE,isSecondaryAlignment=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE)
audit<-list()
for(s in c('CTRL','LOY')){
 sample<-if(s=='CTRL')'HT1080-sgCTRL-D0' else 'HT1080-LOY-D0'
 bam<-file.path('A:/Dropbox/LOY_paralogs/projects/WGS/hg38_bam',paste0(sample,'.bam'))
 for(i in seq_len(nrow(sel))){
  b<-sel[i];region<-GRanges(b$chrom,IRanges(b$start,b$end))
  positions<-scanBam(bam,param=ScanBamParam(flag=flags,mapqFilter=37,which=region,what='pos'))[[1]]$pos
  count<-sum(positions>=b$start&positions<=b$end)
  expected<-b[[paste0(s,'_raw_first_mate_counts')]]
  audit[[length(audit)+1]]<-data.table(sample=sample,bin_id=b$bin_id,independent_count=count,export_count=expected)
 }
}
audit<-rbindlist(audit);fwrite(audit,file.path(out,'WGS_independent_bin_count_verification.csv'))
stopifnot(all(audit$independent_count==audit$export_count))
result<-list(status='passed_with_reported_boundary_gene',RNA_genes=nrow(r),Wald_BH_and_contrast_checks=TRUE,reversed_reference_max_logFC_error=max_delta,reversed_reference_max_SE_error=max_se,unstable_Wald_magnitude_genes=r$gene_name[unstable],LRT_FDR005=sum(lr$padj<.05,na.rm=TRUE),Wald_and_LRT_FDR005=sum(lr$padj<.05&z$padj<.05,na.rm=TRUE),WGS_independent_bins_checked=nrow(audit),WGS_counts_exact=TRUE,genome_coordinates_and_missing_bin_checks=TRUE)
write_json(result,file.path(out,'model_verification.json'),pretty=TRUE,auto_unbox=TRUE);print(result)
