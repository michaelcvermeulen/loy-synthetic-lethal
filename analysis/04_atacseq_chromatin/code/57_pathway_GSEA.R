suppressPackageStartupMessages({library(fgsea);library(msigdbr);library(data.table);library(BiocParallel);library(jsonlite)})
args0 <- commandArgs(FALSE)
script <- sub('^--file=', '', args0[grepl('^--file=',args0)])
ROOT <- normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
OUT <- file.path(ROOT,'integration/pathway_GSEA')
dir.create(file.path(OUT,'sources'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(OUT,'figures'),showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
read <- function(path) {x<-fread(file.path(ROOT,path),na.strings=c('','NA'));for(k in intersect(c('entrez_id','entrezid'),names(x)))set(x,j=k,value=as.character(x[[k]]));x}
combined <- read('integration/feature_scatters/feature_gene_effects.csv')[feature=='promoter']
separate <- read('integration/batch_RNA_overlaps/all_gene_results.csv')[feature=='promoter']
universe <- sort(Reduce(intersect,list(combined$entrez_id,separate[batch=='AUG_2025',entrez_id],separate[batch=='MAY_2026',entrez_id])))
stopifnot(length(universe)==11891L,!anyDuplicated(universe))
sources <- c('integration/feature_scatters/feature_gene_effects.csv','integration/batch_RNA_overlaps/all_gene_results.csv',
             'integration/feature_overlaps/peak_gene_link_audit.csv','differential/primary_TMM_all.csv','differential/passage_TMM_all.csv',
             'config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv',
             'integration/batch_RNA_overlaps/sources/ht1080_AUG_2025_LOY_v_CTRL_DESeq2.csv',
             'integration/batch_RNA_overlaps/sources/ht1080_MAY_2026_LOY_v_CTRL_DESeq2.csv')
fwrite(data.table(path=sources,md5=as.character(tools::md5sum(file.path(ROOT,sources)))),file.path(OUT,'sources/input_checksums.csv'))

# Freeze human-native Entrez gene sets. Existing local msigdbr cache is read without modification.
hallmark <- as.data.table(msigdbr(db_species='HS',species='human',collection='H'))
reactome <- as.data.table(msigdbr(db_species='HS',species='human',collection='C2',subcollection='CP:REACTOME'))
sets <- rbindlist(list(hallmark,reactome),fill=TRUE)
sets <- unique(sets[!is.na(ncbi_gene)&ncbi_gene!='',.(pathway=gs_name,collection=ifelse(gs_collection=='H','Hallmark','Reactome'),entrez_id=as.character(ncbi_gene),gene_symbol,db_version,gs_description,gs_url)])
fwrite(sets,file.path(OUT,'sources/gene_set_membership.csv'))
pathways <- lapply(split(sets$entrez_id,sets$pathway),unique)
annotation <- unique(sets[,.(pathway,collection,db_version,gs_description,gs_url)])
stopifnot(!anyDuplicated(annotation$pathway))
fwrite(annotation,file.path(OUT,'sources/pathway_annotation.csv'))
saveRDS(pathways,file.path(OUT,'sources/pathways.rds'))

links <- unique(read('integration/feature_overlaps/peak_gene_link_audit.csv')[feature=='promoter'&entrez_id %in% universe,.(entrez_id,peak_id,chrom)])
primary <- read('differential/primary_TMM_all.csv')
passage <- read('differential/passage_TMM_all.csv')
linked <- merge(links,primary[,.(peak_id,logFC,logCPM,F,PValue,FDR)],by='peak_id',all.x=TRUE)
stopifnot(all(is.finite(linked$F)),all(linked$F>=0),all(is.finite(linked$logFC)))
linked[,signed_stat:=sign(logFC)*sqrt(F)]
linked <- merge(linked,passage[,.(peak_id,passage_logFC=logFC,passage_F=F)],by='peak_id',all.x=TRUE)
stopifnot(all(is.finite(linked$passage_F)),all(linked$passage_F>=0))
linked[,passage_signed_stat:=sign(passage_logFC)*sqrt(passage_F)]
fwrite(linked,file.path(OUT,'promoter_peak_ranking_audit.csv'))
ranking <- linked[,.(chrom=unique(chrom),promoter_peaks=.N,ATAC=median(signed_stat),ATAC_median_logFC=median(logFC),ATAC_passage=median(passage_signed_stat)),by=entrez_id]
# Representative-peak sensitivity uses abundance, never a condition-effect or significance threshold.
representative <- linked[order(-logCPM,peak_id),.SD[1L],by=entrez_id][,.(entrez_id,representative_peak=peak_id,ATAC_representative=signed_stat)]
ranking <- merge(ranking,representative,by='entrez_id')
ranking <- merge(ranking,combined[,.(entrez_id,symbol=rna_symbol)],by='entrez_id')
rna_paths <- c(RNA_Combined=sources[6],RNA_AUG_2025=sources[7],RNA_MAY_2026=sources[8])
for (name in names(rna_paths)) {
  frame <- read(rna_paths[[name]])
  dup <- frame[, .N, by=entrezid][N>1,entrezid]
  frame <- frame[!entrezid %in% dup & entrezid %in% universe]
  stopifnot(nrow(frame)==length(universe),all(is.finite(frame$t)))
  values <- frame[,.(entrez_id=entrezid,stat=t,logFC,FDR)]
  setnames(values,c('stat','logFC','FDR'),paste0(name,c('','_logFC','_FDR')))
  ranking <- merge(ranking,values,by='entrez_id')
}
stopifnot(nrow(ranking)==length(universe),!anyDuplicated(ranking$entrez_id),!any(ranking$chrom %in% c('chrY','chrM')))
setorder(ranking,entrez_id)
fwrite(ranking,file.path(OUT,'gene_ranks.csv'))

main_names <- c('ATAC','RNA_Combined','RNA_AUG_2025','RNA_MAY_2026')
analyses <- c(main_names,'ATAC_representative','ATAC_median_logFC','ATAC_passage',paste0(main_names,'_no_chr19'),'ATAC_tie_shuffle_1','ATAC_tie_shuffle_2')
all_results <- list(); audit <- list(); warnings_seen <- character()
for (i in seq_along(analyses)) {
  name <- analyses[i]; column <- sub('_tie_shuffle_[12]$','',sub('_no_chr19$','',name))
  selected <- if(grepl('_no_chr19$',name)) ranking[chrom!='chr19'] else ranking
  values <- setNames(selected[[column]],selected$entrez_id)
  values <- values[order(-values,names(values))]
  if(grepl('_tie_shuffle_',name)){set.seed(202609100L+i);values<-values[order(-values,runif(length(values)))]}
  stopifnot(all(is.finite(values)),!anyDuplicated(names(values)),any(values>0),any(values<0))
  eligible_sizes <- vapply(pathways,function(p)length(intersect(p,names(values))),integer(1))
  fwrite(data.table(pathway=names(pathways),eligible_genes=eligible_sizes,tested=eligible_sizes>=15&eligible_sizes<=500),file.path(OUT,paste0('pathway_sizes_',name,'.csv')))
  cat('Running',name,'with',length(values),'genes and',sum(eligible_sizes>=15&eligible_sizes<=500),'pathways\n');flush.console()
  set.seed(if(grepl('_tie_shuffle_',name))20260910L else 20260909L+i)
  res <- withCallingHandlers(fgseaMultilevel(pathways=pathways,stats=values,minSize=15,maxSize=500,
    eps=0,sampleSize=101,nPermSimple=5000,scoreType='std',gseaParam=1,BPPARAM=SerialParam()),
    warning=function(w){warnings_seen<<-c(warnings_seen,paste(name,conditionMessage(w)));invokeRestart('muffleWarning')})
  res <- merge(as.data.table(res),annotation[,.(pathway,collection)],by='pathway')
  res[,padj_global:=p.adjust(pval,'BH')]
  res[,padj:=p.adjust(pval,'BH'),by=collection]
  res[,analysis:=name]
  res[,leading_edge_entrez:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))]
  symbols <- setNames(ranking$symbol,ranking$entrez_id)
  res[,leading_edge_symbols:=vapply(leadingEdge,function(g)paste(symbols[g],collapse=';'),character(1))]
  res[,leadingEdge:=NULL]
  fwrite(res,file.path(OUT,paste0('fgsea_',name,'.csv')))
  all_results[[name]] <- res
  audit[[name]] <- data.table(analysis=name,genes=length(values),positive=sum(values>0),negative=sum(values<0),zero=sum(values==0),tied_nonzero=sum(duplicated(values[values!=0])),pathways=nrow(res),finite_NES=sum(is.finite(res$NES)),finite_p=sum(is.finite(res$pval)))
  cat(name,':',sum(res$padj<.05,na.rm=TRUE),'pathways at collection-wise FDR < 0.05\n');flush.console()
}
results <- rbindlist(all_results,use.names=TRUE)
fwrite(results,file.path(OUT,'fgsea_all_results.csv'))
fwrite(rbindlist(audit),file.path(OUT,'rank_analysis_audit.csv'))
writeLines(warnings_seen,file.path(OUT,'fgsea_warnings.txt'))
summary <- list(status='GSEA_complete_pending_figures',common_genes=nrow(ranking),collections=unique(sets[,.(collection,db_version)]),
  minSize=15,maxSize=500,scoreType='std',gseaParam=1,nPermSimple=5000,sampleSize=101,eps=0,
  primary_ATAC_rank='Median of sign(logFC)*sqrt(edgeR QL F) across all tested linked promoter peaks; not a calibrated gene-level test statistic.',
  RNA_rank='Original signed DESeq2 Wald statistic in source column t; LOY versus CTRL.',
  FDR='BH within collection and analysis; padj_global also supplied across both collections within each analysis.',
  sensitivity_analyses=c('Most abundant promoter peak by primary logCPM, tie-break peak ID','Median promoter logFC','Passage-adjusted median signed sqrt(F)','Remove chr19 from every assay ranking and gene set intersection','Two seeded random orderings within tied ATAC ranks, with the same fgsea random seed as primary ATAC'),
  warnings=warnings_seen)
write_json(summary,file.path(OUT,'analysis_settings.json'),pretty=TRUE,auto_unbox=TRUE)
cat('GSEA tables complete.\n')
