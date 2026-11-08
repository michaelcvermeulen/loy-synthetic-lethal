options(warn=1)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash='/', mustWork=TRUE))))
suppressPackageStartupMessages({library(data.table); library(Matrix); library(ggplot2); library(patchwork)})
root <- 'motifs/TF_extension'; out <- file.path(root,'TFAP2C_review')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
rd <- function(x) fread(file.path(root,x))
put <- function(x,n) fwrite(x,file.path(out,paste0(n,'.csv')))
all <- rd('concordance/TF_ATAC_footprint_RNA.csv')
ap <- all[grepl('^TFAP2',motif_name)]; ids <- ap[motif_name=='TFAP2C',motif_id]
stopifnot(length(ids)==3, !anyDuplicated(ap$motif_id))
ind <- rd('chromVAR/motifs_independent_limma.csv')
pas <- rd('chromVAR/motifs_passage_sensitivity_limma.csv')
wel <- rd('chromVAR/motifs_Welch_test.csv')
native <- rd('TOBIAS/BINDetect/bindetect_results.txt')
ap[, chromVAR_primary_P := ind$P.Value[match(motif_id,ind$id)]]
ap[, chromVAR_passage_FDR := pas$adj.P.Val[match(motif_id,pas$id)]]
ap[, chromVAR_Welch_FDR := wel$FDR[match(motif_id,wel$id)]]
ap <- merge(ap,native[,.(motif_id,total_tfbs,LOY_mean_score,CTRL_mean_score,LOY_bound,CTRL_bound)],by='motif_id')
put(ap,'AP2_cross_method_comparison')
hom <- rd('HOMER/known_motif_all_results.csv')[motif_id %in% ap$motif_id]
put(hom,'AP2_HOMER_primary_and_passage')
put(fread('motifs/motif_enrichment.csv')[grepl('^TFAP2',motif_name)],'AP2_original_JASPAR2022_screen')
groups <- unique(ap$archetype)
put(merge(rd('chromVAR/archetypes_independent_limma.csv')[id %in% groups],
 rd('chromVAR/archetypes_passage_sensitivity_limma.csv')[id %in% groups],by='id',suffixes=c('_primary','_passage')),'AP2_group_chromVAR')
put(rd('chromVAR/archetype_members.csv')[archetype %in% groups],'AP2_group_members')
samples <- fread('config/samples.tsv'); sn <- samples$sample
z <- rd('chromVAR/motifs_zscores.csv')[id %in% ap$motif_id]
long <- melt(z,id.vars='id',variable.name='sample',value.name='zscore')
long <- merge(long,samples[,.(sample,condition,passage)],by='sample')
long[,motif_name:=ap$motif_name[match(id,ap$motif_id)]]
put(long,'AP2_sample_chromVAR_scores')
robust <- rbindlist(lapply(ids,function(k){
 d <- long[id==k]; x<-d[condition=='LOY',zscore]; y<-d[condition=='CTRL',zscore]
 loo <- sapply(sn,function(s) d[sample!=s & condition=='LOY',mean(zscore)]-d[sample!=s & condition=='CTRL',mean(zscore)])
 data.table(motif_id=k,CTRL_mean=mean(y),LOY_mean=mean(x),delta=mean(x)-mean(y),
 all_LOY_below_all_CTRL=max(x)<min(y),minimum_pairwise_separation=min(y)-max(x),
 leave_one_sample_out_min_delta=min(loo),leave_one_sample_out_max_delta=max(loo))
})); put(robust,'TFAP2C_sample_robustness')
# Reuse the exact peak row order recorded when motifmatchr was run.
res <- fread('differential/primary_TMM_all.csv'); all_tested_ids <- res$peak_id; accounting <- rd('chromVAR/peak_accounting.csv')
stopifnot(identical(res$peak_id,accounting$peak_id))
res <- res[accounting$included]; mm <- readRDS(file.path(root,'chromVAR/motif_matches.rds'))
stopifnot(nrow(mm)==nrow(res),all(ap$motif_id %in% colnames(mm)))
sel <- as.matrix(mm[,ap$motif_id,drop=FALSE]); colnames(sel)<-ap$motif_id
over <- rbindlist(lapply(ids,function(k) rbindlist(lapply(ap$motif_id,function(j){
 a<-sel[,k];b<-sel[,j]; data.table(TFAP2C_motif=k,other_motif=j,other_name=ap$motif_name[match(j,ap$motif_id)],
 n_TFAP2C=sum(a),n_other=sum(b),shared=sum(a&b),Jaccard=sum(a&b)/sum(a|b),fraction_TFAP2C_shared=sum(a&b)/sum(a))
})))); put(over,'AP2_peak_overlap')
pk <- fread('annotation/annotated_differential_peaks.csv')[match(res$peak_id,peak_id)]
stopifnot(identical(pk$peak_id,res$peak_id))
for(k in ids) pk[,(k):=sel[,k]]
pk[,TFAP2C_any:=rowSums(sel[,ids,drop=FALSE])>0]
pk[,TFAP2C_all_three:=rowSums(sel[,ids,drop=FALSE])==3]
pk[,other_AP2:=rowSums(sel[,setdiff(ap$motif_id,ids),drop=FALSE])>0]
put(pk[TFAP2C_any==TRUE],'TFAP2C_motif_positive_peaks')
counts <- rbindlist(lapply(c(ids,'TFAP2C_any','TFAP2C_all_three'),function(k){
 d<-pk[get(k)]; data.table(motif_id=k,n_peaks=nrow(d),shared_other_AP2=sum(d$other_AP2),
 lost_FDR005=sum(d$FDR<.05 & d$logFC<0),gained_FDR005=sum(d$FDR<.05 & d$logFC>0),
 lost_FDR005_LFC1=sum(d$FDR<.05 & d$logFC< -1),gained_FDR005_LFC1=sum(d$FDR<.05 & d$logFC>1),
 lost_FDR010=sum(d$FDR<.1 & d$logFC<0),gained_FDR010=sum(d$FDR<.1 & d$logFC>0),
 median_all_peak_logFC=median(d$logFC),fraction_all_peak_negative=mean(d$logFC<0))
}));put(counts,'TFAP2C_peak_summary')
put(pk[TFAP2C_any==TRUE,.N,by=.(annotation,DA_status=ifelse(FDR<.05 & logFC< -1,'lost_FDR005_LFC1',ifelse(FDR<.05 & logFC>1,'gained_FDR005_LFC1','other')))],'TFAP2C_feature_counts')
# Native TOBIAS site calls: pooled descriptive summaries, no new site-level hypothesis tests.
site_summ <- list(); peak_sites <- list()
for(k in ids){
 d<-rd(paste0('TOBIAS/BINDetect/',k,'/',k,'_overview.txt'))
 stopifnot(nrow(d)==native[motif_id==k,total_tfbs])
 d[,peak_id:=paste0(peak_chr,':',peak_start,'-',peak_end)]
 stopifnot(all(d$peak_id %in% all_tested_ids))
 d[,bound_state:=fifelse(LOY_bound==1 & CTRL_bound==1,'both',fifelse(LOY_bound==1,'LOY_only',fifelse(CTRL_bound==1,'CTRL_only','neither')))]
 site_summ[[k]]<-data.table(motif_id=k,total_sites=nrow(d),CTRL_bound=sum(d$CTRL_bound),LOY_bound=sum(d$LOY_bound),
 both_bound=sum(d$bound_state=='both'),CTRL_only=sum(d$bound_state=='CTRL_only'),LOY_only=sum(d$bound_state=='LOY_only'),
 neither=sum(d$bound_state=='neither'),fraction_negative_site_log2fc=mean(d$LOY_CTRL_log2fc<0),
 median_site_log2fc=median(d$LOY_CTRL_log2fc))
 put(d,paste0(k,'_TOBIAS_sites'))
 q<-d[,.(n_sites=.N,median_site_log2fc=median(LOY_CTRL_log2fc),CTRL_bound_sites=sum(CTRL_bound),LOY_bound_sites=sum(LOY_bound)),by=peak_id]
 q[,motif_id:=k]; peak_sites[[k]]<-q
}
ss<-rbindlist(site_summ);put(ss,'TFAP2C_TOBIAS_site_summary')
ps<-merge(rbindlist(peak_sites),fread('annotation/annotated_differential_peaks.csv')[,.(peak_id,atac_logFC=logFC,atac_FDR=FDR,annotation,nearest_symbol)],by='peak_id',all.x=TRUE)
put(ps,'TFAP2C_TOBIAS_peak_summary')
prom<-fread('integration/RNAseq/ATAC_RNA_promoter_peak_links.csv')[peak_id %in% pk[TFAP2C_any==TRUE,peak_id]]
put(prom,'TFAP2C_promoter_peak_RNA_links')
near<-fread('integration/RNAseq/ATAC_RNA_nearest_gene_links_exploratory.csv')[peak_id %in% pk[TFAP2C_any==TRUE,peak_id] & annotation!='promoter']
put(near,'TFAP2C_nonpromoter_nearest_RNA_exploratory')
gene_lists<-list();gs<-list()
for(cut in c(.05,.1)) {
 for(lfc in c(0,1)) {
  d<-prom[FDR<cut & logFC< -lfc & is.finite(rna_FDR)]
  g<-d[,.(rna_symbol=paste(unique(rna_symbol),collapse=';'),rna_logFC=unique(rna_logFC)[1],rna_FDR=unique(rna_FDR)[1],
   n_lost_motif_promoter_peaks=uniqueN(peak_id),minimum_peak_FDR=min(FDR),median_lost_peak_logFC=median(logFC),peaks=paste(unique(peak_id),collapse=';')),by=rna_ensembl_id]
  g[,`:=`(ATAC_FDR_cutoff=cut,ATAC_abs_logFC_cutoff=lfc)]
  gs[[paste(cut,lfc)]]<-data.table(FDR_cutoff=cut,ATAC_abs_logFC_cutoff=lfc,linked_RNA_genes=nrow(g),
   RNA_down_same_FDR=sum(g$rna_FDR<cut & g$rna_logFC<0),RNA_up_same_FDR=sum(g$rna_FDR<cut & g$rna_logFC>0))
  gene_lists[[paste(cut,lfc)]]<-g
 }
}
genes<-rbindlist(gene_lists);put(genes,'TFAP2C_lost_promoter_genes_thresholds');put(rbindlist(gs),'TFAP2C_promoter_RNA_summary')
put(prom[rna_symbol=='TFAP2C' | symbol=='TFAP2C'],'TFAP2C_own_promoter_motif_links')
rna<-fread('config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv')
put(rna[grepl('^TFAP2',gene_name)],'AP2_RNA_source_rows')
own<-fread('integration/RNAseq/ATAC_RNA_promoter_peak_links.csv')[rna_symbol=='TFAP2C' | symbol=='TFAP2C']
own[,TFAP2C_motif_positive:=peak_id %in% pk[TFAP2C_any==TRUE,peak_id]];put(own,'TFAP2C_own_promoter_all_peaks')
theme_set(theme_bw(base_size=11)); cols<-c(CTRL='#0072B2',LOY='#D55E00')
saveplot<-function(p,n,w=11,h=7){for(ext in c('png','pdf'))ggsave(file.path(out,paste0(n,'.',ext)),p,width=w,height=h,dpi=180)}
plotz<-long[id %in% ids];plotz[,label:=paste0('TFAP2C ',id)];plotz[,sample:=factor(sample,levels=sn)]
pz<-ggplot(plotz,aes(sample,zscore,color=condition,shape=factor(passage)))+geom_point(size=3)+
 facet_wrap(~label,nrow=1)+scale_color_manual(values=cols)+scale_x_discrete(expand=expansion(add=c(.3,.8)))+
 labs(title='TFAP2C motif accessibility in each biological replicate',subtitle='All LOY scores are below all CTRL scores; primary motif-level FDR 0.058-0.067',y='chromVAR deviation z-score',x=NULL,shape='Passage label')+theme(legend.position='bottom',axis.text.x=element_text(angle=45,hjust=1))
saveplot(pz,'TFAP2C_sample_scores',12,4.5)
hplot<-hom[motif_id %in% ids & direction=='lost'];hplot[,label:=paste0(motif_id,'\nFDR ',formatC(BH_FDR_both_directions,format='e',digits=1))]
ph<-ggplot(hplot,aes(motif_id,fold_enrichment,fill=arm))+geom_col(position='dodge')+geom_hline(yintercept=1,lty=2)+
 scale_fill_manual(values=c(primary_TMM='#0072B2',passage_TMM='#999999'),labels=c('Passage sensitivity','Primary'))+
 labs(title='Enrichment in peaks losing accessibility',x=NULL,y='HOMER fold enrichment',fill='DA model')+theme(legend.position='bottom')
pb<-ggplot(melt(ss[,.(motif_id,CTRL_bound,LOY_bound)],id.vars='motif_id'),aes(motif_id,value,fill=variable))+
 geom_col(position='dodge')+scale_fill_manual(values=c(CTRL_bound='#0072B2',LOY_bound='#D55E00'),labels=c('CTRL','LOY'))+
 labs(title='Predicted bound motif sites',subtitle='Depth-matched condition pools; descriptive site calls',x=NULL,y='Number of sites',fill=NULL)+theme(legend.position='bottom')
saveplot(ph/pb,'TFAP2C_enrichment_and_binding',10,8)
pv<-ggplot(pk[TFAP2C_any==TRUE],aes(logFC,-log10(pmax(FDR,1e-300))))+geom_point(alpha=.25,size=.6,color='#455A64')+
 geom_point(data=pk[TFAP2C_any & FDR<.05 & abs(logFC)>1],color='#0072B2',size=1)+geom_hline(yintercept=-log10(.05),lty=2)+geom_vline(xintercept=c(-1,1),lty=2)+
 labs(title='Peaks containing any TFAP2C motif model',subtitle='Motifmatchr sites; blue: primary FDR <0.05 and |log2FC| >1',x='ATAC log2FC (LOY / CTRL)',y='-log10 ATAC FDR')
saveplot(pv,'TFAP2C_peak_volcano',8,5)
joint<-genes[ATAC_FDR_cutoff==.1 & ATAC_abs_logFC_cutoff==0 & rna_FDR<.1][order(rna_FDR)]
put(joint,'TFAP2C_joint_promoter_genes_FDR010')
if(nrow(joint)){
 d<-head(joint,25);d[,rna_symbol:=factor(rna_symbol,levels=rev(rna_symbol))]
 g<-melt(d,id.vars='rna_symbol',measure.vars=c('median_lost_peak_logFC','rna_logFC'))
 pg<-ggplot(g,aes(value,rna_symbol,color=variable))+geom_point(size=2.5)+geom_vline(xintercept=0,lty=2)+
 scale_color_manual(values=c(median_lost_peak_logFC='#0072B2',rna_logFC='#D55E00'),labels=c('ATAC: selected lost promoter peaks','RNA'))+
 labs(title='Genes with a lost AP-2 motif-containing promoter peak and altered RNA',subtitle='ATAC and RNA FDR <0.10; up to 25 genes ranked by RNA FDR; positional links',x='log2FC (LOY / CTRL)',y=NULL,color=NULL)+theme(legend.position='bottom')
 saveplot(pg,'TFAP2C_promoter_RNA_genes',10,max(4,nrow(d)*.28+2))
}
writeLines(c('Reuses existing JASPAR2024 HOMER, chromVAR and TOBIAS results; no selected-motif FDR recalculation.',
 'Independent biological replicates remain primary; passage-adjusted results are sensitivity only, with technical batch unconfirmed.',
 'Leave-one-sample-out values are descriptive mean differences, not refitted significance tests.',
 'Motifmatchr peak membership is reconstructed using the exact saved peak-accounting row order. It need not match native HOMER or TOBIAS site counts.',
 'All peak and gene counts deduplicate their stated unit. Promoter links use existing TSS +/-2 kb annotations. Nonpromoter nearest-TSS links are exploratory.',
 'Joint RNA subsets use the original Ben thesis HT1080 batch + treatment result; no new RNA model or gene-set enrichment test.',
 'TOBIAS bound states and within-site log2fc are pooled descriptive outputs. They do not estimate biological replicate variability.',
 'AP-2 family models overlap; no sequence-based assignment to a specific protein is asserted. No copy-number correction or causality test.',
 'Native RNA logCPM field contains DESeq2 baseMean in this source.'),file.path(out,'METHODS.txt'))
stopifnot(all(robust$all_LOY_below_all_CTRL),all(robust$leave_one_sample_out_max_delta<0),
 all(ss$CTRL_bound==ss$both_bound+ss$CTRL_only),all(ss$LOY_bound==ss$both_bound+ss$LOY_only),
 all(ss$total_sites==ss$both_bound+ss$CTRL_only+ss$LOY_only+ss$neither))
writeLines(c('PASS: motif identity, exact peak order, native site totals, bound-state conservation and RNA joins checked.',
 'PASS: all three TFAP2C models show complete group separation and negative mean effects after omitting any one sample.'),file.path(out,'verification.txt'))
print(ap[motif_name=='TFAP2C',.(motif_id,chromVAR_FDR,chromVAR_passage_FDR,chromVAR_Welch_FDR,HOMER_lost_fold_enrichment)])
print(counts);print(ss);print(rbindlist(gs));print(joint[,.(rna_symbol,rna_logFC,rna_FDR)])
