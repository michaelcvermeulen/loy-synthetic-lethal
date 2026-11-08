options(warn=1)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE))))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
out <- 'motifs/TF_extension/concordance'; dir.create(out,recursive=TRUE,showWarnings=FALSE)
cv <- fread('motifs/TF_extension/chromVAR/motifs_independent_limma.csv')
rna <- fread('config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv')
joined <- rbindlist(lapply(seq_len(nrow(cv)),function(i) {
 terms <- strsplit(sub('\\(.*','',cv$motif_name[i]),'::',fixed=TRUE)[[1]]
 rbindlist(lapply(terms,function(tf) {
  hits <- which(rna$gene_name==tf); method <- 'exact_symbol'
  if(!length(hits)){hits <- which(toupper(rna$gene_name)==toupper(tf));method <- 'case_insensitive_candidate'}
  if(length(hits)!=1) return(data.table(motif_id=cv$id[i],motif_name=cv$motif_name[i],TF_label=tf,
   mapping_status=if(length(hits))'ambiguous_RNA_symbol'else'absent_RNA_symbol',rna_gene_id=NA_character_,
   rna_symbol=NA_character_,rna_log2FC=NA_real_,rna_FDR=NA_real_))
  k <- hits[1]
  data.table(motif_id=cv$id[i],motif_name=cv$motif_name[i],TF_label=tf,mapping_status=method,
   rna_gene_id=rna$gene_id[k],rna_symbol=rna$gene_name[k],rna_log2FC=rna$logFC[k],rna_FDR=rna$FDR[k])
 }))
}))
joined <- merge(joined,cv[,.(motif_id=id,archetype,chromVAR_delta_z=LOY_minus_CTRL_zscore,chromVAR_FDR=adj.P.Val)],by='motif_id',all.x=TRUE)
joined[,rna_chromosome := rna$chr[match(rna_gene_id,rna$gene_id)]]
joined[,rna_Y_linked := fifelse(is.na(rna_chromosome),NA,rna_chromosome=='Y')]
joined[,RNA_statistic_status := fifelse(is.na(rna_gene_id),'not_mapped',
 fifelse(is.finite(rna_log2FC)&is.finite(rna_FDR),'finite','nonfinite'))]
coverage <- rbindlist(lapply(c('ZFY','ZFX','SRY'),function(symbol) {
 rr <- rna[gene_name==symbol]; mm <- cv[toupper(motif_name)==toupper(symbol)]
 data.table(TF=symbol,RNA_status=if(nrow(rr)==1)'present'else if(nrow(rr)==0)'absent_from_selected_table'else'ambiguous',
  RNA_chromosome=if(nrow(rr)==1)rr$chr else NA_character_,
  RNA_log2FC=if(nrow(rr)==1)rr$logFC else NA_real_,RNA_FDR=if(nrow(rr)==1)rr$FDR else NA_real_,
  motif_ids=paste(mm$id,collapse=';'),motif_labels=paste(mm$motif_name,collapse=';'),
  motif_label_status=if(nrow(mm))'label_match_including_case_insensitive_candidates'else'no_matching_motif_label')
}))
fwrite(coverage,file.path(out,'sex_chromosome_TF_coverage.csv'))
setorder(joined,chromVAR_FDR,motif_id)
fwrite(joined,file.path(out,'chromVAR_TF_RNA_all.csv'))
fwrite(joined[chromVAR_FDR<.05],file.path(out,'chromVAR_FDR005_TF_RNA.csv'))
plotdata <- joined[is.finite(rna_log2FC)]
p <- ggplot(plotdata,aes(chromVAR_delta_z,rna_log2FC))+geom_point(aes(colour=chromVAR_FDR<.05),alpha=.5,size=1.5)+
 geom_hline(yintercept=0,linetype=2,colour='grey60')+geom_vline(xintercept=0,linetype=2,colour='grey60')+
 geom_text(data=plotdata[chromVAR_FDR<.05],aes(label=paste(TF_label,motif_id)),check_overlap=TRUE,hjust=0,nudge_x=.2,size=3)+
 scale_colour_manual(values=c('FALSE'='grey65','TRUE'='#0072b2'),name='chromVAR FDR <0.05')+
 theme_bw(base_size=11)+labs(x='chromVAR z-score difference (LOY - CTRL)',y='TF RNA log2 fold change (LOY / CTRL)',
 title='Motif accessibility and TF transcript changes',subtitle='RNA: existing HT1080 DESeq2 batch + treatment; ATAC: independent replicates',
 caption='Each point is a motif-TF link; related motifs are not independent evidence. Case-insensitive symbol links are candidates.\nAssays are compared at group level. RNA and motif accessibility directions need not agree.')
for(ext in c('png','pdf')) ggsave(file.path(out,paste0('chromVAR_TF_RNA.',ext)),p,width=11,height=7,dpi=220)
writeLines(c('RNA source unchanged: HT1080_BATCHES DESeq2 batch + treatment LOY/CTRL.',
 'TF transcript annotation retains Y-linked genes; chromosome and Y-linkage flags are recorded. The original primary gene-overlap universe is unchanged.',
 'Present RNA genes with nonfinite statistics are retained and flagged; they are not mislabeled as absent genes.',
 'The selected RNA table contains decreased ZFY RNA, but the JASPAR2024 set has no motif named ZFY. A motif-label coverage gap is not a result about specific TF occupancy.',
 'TF motif labels are joined to unique RNA gene symbols. Exact and case-insensitive candidate matches are labeled; no orthology database was applied.',
 'Dimer motifs retain one link per member; ambiguous and absent symbols are retained as unresolved.',
 'RNA transcript effects do not establish TF protein activity or require concordant motif accessibility direction.',
 'The RNA expression comparison is across groups, not directly matched samples.'),file.path(out,'METHODS.txt'))
print(joined[chromVAR_FDR<.05,.(motif_id,TF_label,mapping_status,chromVAR_delta_z,chromVAR_FDR,rna_log2FC,rna_FDR)])
