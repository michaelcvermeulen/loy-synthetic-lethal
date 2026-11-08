# Align add-back follow-up to the current September HT1080 RNA/ATAC integration.
# Reuses frozen add-back coefficients and voom fit; no DE model is refitted.
suppressPackageStartupMessages({library(data.table); library(limma); library(digest)})
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
fig <- dirname(dirname(normalizePath(script, winslash='/')))
projects <- dirname(dirname(dirname(fig)))
out <- file.path(fig, 'source_data/addback_atac_overlap')
add <- file.path(projects, 'RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026/outputs')
src <- c(integration=file.path(fig,'source_data/ht1080_rna_atac_integration/integration_genes.csv'),
 baseline=file.path(projects,'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1/tables/gene_consistency.csv'),
 fit=file.path(add,'fit_primary.rds'), signature=file.path(add,'signature_150/signature_mapping.csv'))
keys <- c('UTY_v_LXHA','KDM5D_v_LXHA','DUAL_v_LXHA','LOY_v_LXHA')
src <- c(src,setNames(file.path(add,'de',paste0('de_',keys,'.csv')),keys),
 may=file.path(add,'de/de_LOY_v_CTRL.csv'))
stopifnot(all(file.exists(src)))
g <- fread(src[['integration']]); stopifnot(nrow(g)==13300, !anyDuplicated(g$ensembl))
sets <- list(joint_FDR005=g[rna_class=='down' & atac_class=='lost',ensembl],
 joint_ATAC010_RNA005=g[rna_class=='down' & atac_class_010=='lost',ensembl],
 ATAC_promoter_lost=g[atac_class=='lost',ensembl],
 ATAC_promoter_gained=g[atac_class=='gained',ensembl])
stopifnot(length(sets[[1]])==39, length(sets[[2]])==58)
strip <- function(x) sub('[.].*$', '',x)
d <- rbindlist(lapply(names(sets),function(s) data.table(set=s,gene_id=sets[[s]])))
d <- merge(d,g[,.(gene_id=ensembl,gene_name,baseline_logFC=rna_logFC,baseline_FDR=rna_FDR,
 atac_median_logFC,atac_min_peak_FDR,atac_class,atac_class_010)],by='gene_id')
for(k in keys){
 x<-fread(src[[k]])[,.(gene_id=strip(gene_id),logFC,FDR)]
 stopifnot(!anyDuplicated(x$gene_id));setnames(x,c('logFC','FDR'),paste0(k,c('_logFC','_FDR')))
 d<-merge(d,x,by='gene_id',all.x=TRUE)
 d[,(paste0(k,'_opposite_baseline')):=sign(get(paste0(k,'_logFC')))==-sign(baseline_logFC)]
}
fit<-readRDS(src[['fit']]); v<-fit$v; cm<-fit$cm
ids<-strip(rownames(v));stopifnot(!anyDuplicated(ids))
index<-lapply(sets,function(z) which(ids %in% z))
d[,tested_in_addback:=gene_id %in% ids]
counts<-rbindlist(lapply(names(sets),function(s) rbindlist(lapply(keys,function(k){
 x<-d[set==s & tested_in_addback]; opp<-x[[paste0(k,'_opposite_baseline')]];fdr<-x[[paste0(k,'_FDR')]]
 data.table(set=s,contrast=k,genes_in_set=length(sets[[s]]),tested_in_addback=nrow(x),
 opposite_baseline=sum(opp),opposite_baseline_FDR005=sum(opp & fdr<.05),same_baseline_FDR005=sum(!opp & fdr<.05))
}))))
# Joint down sets have a defined corrective direction. All-promoter sets can mix RNA
# directions, so their rotation tests report expression direction without calling it rescue.
fr<-rbindlist(lapply(keys,function(k){
 x<-fry(v,index=index,design=v$design,contrast=cm[,k])
 data.table(set=rownames(x),contrast=k,NGenes=x$NGenes,Direction=x$Direction,PValue=x$PValue)
}))
fr[,FDR_across_16_tests:=p.adjust(PValue,'BH')]
stopifnot(all(fr$NGenes==vapply(fr$set,function(s)length(index[[s]]),integer(1))))
fwrite(d,file.path(out,'gene_table.csv'));fwrite(counts,file.path(out,'gene_level_counts.csv'))
fwrite(fr,file.path(out,'fry_tests.csv'))
# Update the baseline agreement lookup as well; signature membership is not changed.
b<-fread(src[['baseline']]);b[,gene_id:=strip(gene_id)]
sig<-fread(src[['signature']])
may<-fread(src[['may']]);may[,gene_id:=strip(gene_id)]
agreement<-merge(b,may[,.(gene_id,May_logFC=logFC,May_FDR=FDR)],by='gene_id')
fwrite(agreement,file.path(out,'baseline_may_agreement.csv'))
siglookup<-merge(sig,b,by.x='stable_gene_id',by.y='gene_id')
fwrite(siglookup,file.path(out,'signature_baseline_lookup.csv'))
fwrite(data.table(signature_genes=nrow(siglookup),same_direction=sum(siglookup$logFC*siglookup$AUG_difference>0),
 signature_FDR005=sum(siglookup$FDR<.05),shared_genes=nrow(agreement),
 Spearman=cor(agreement$logFC,agreement$May_logFC,method='spearman')),file.path(out,'baseline_agreement_summary.csv'))
fwrite(data.table(role=names(src),path=src,sha256=vapply(src,digest,character(1),file=TRUE,algo='sha256')),file.path(out,'input_manifest.csv'))
writeLines(c('Current source: HT1080 SEPT_11_2026_REANALYZE programme model and thesis_RNA_SEPT11 ATAC integration.',
 'RNA-defined joint sets retain FDR <0.05, absolute log2FC >=1 and period-direction agreement; only ATAC relaxes to FDR <0.10 in the 58-gene set.',
 'All-promoter sets include genes with lost or gained peaks at ATAC FDR <0.05, regardless of RNA significance or direction.',
 'Opposite-baseline counts describe direction, not distance or complete rescue. fry uses the saved add-back fit, with BH correction across 16 tests.',
 'Primary RNA, ATAC and add-back DE fits and the August signature are unchanged.'),file.path(out,'README.md'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(counts);print(fr)
