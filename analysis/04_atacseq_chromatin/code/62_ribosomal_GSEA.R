suppressPackageStartupMessages({library(data.table);library(fgsea);library(BiocParallel);library(jsonlite)})
arg<-commandArgs(FALSE);script<-sub('^--file=','',arg[grepl('^--file=',arg)])
ROOT<-normalizePath(file.path(dirname(script),'..'),winslash='/')
BASE<-file.path(ROOT,'integration/pathway_GSEA');OUT<-file.path(BASE,'ribosomal_sensitivity')
ranks<-fread(file.path(BASE,'gene_ranks.csv'),colClasses=c(entrez_id='character'))
exclude<-fread(file.path(OUT,'excluded_genes.csv'),colClasses=c(entrez_id='character'))
pathways<-readRDS(file.path(BASE,'sources/pathways.rds'))
annotation<-unique(fread(file.path(BASE,'sources/pathway_annotation.csv'))[,.(pathway,collection)])
main<-c('ATAC','RNA_Combined','RNA_AUG_2025','RNA_MAY_2026')
outputs<-list();audit<-list();warnings_seen<-character()
for(mode in c('no_cytosolic','no_cytosolic_or_mitochondrial')){
  drop<-if(mode=='no_cytosolic')exclude[class=='cytosolic',entrez_id] else exclude$entrez_id
  selected<-ranks[!entrez_id %in% drop]
  fwrite(selected,file.path(OUT,paste0('gene_ranks_',mode,'.csv')))
  for(i in seq_along(main)){
    name<-main[i];values<-setNames(selected[[name]],selected$entrez_id);values<-values[order(-values,names(values))]
    sizes<-vapply(pathways,function(p)length(intersect(p,names(values))),integer(1))
    fwrite(data.table(pathway=names(pathways),eligible_genes=sizes,tested=sizes>=15&sizes<=500),file.path(OUT,paste0('pathway_sizes_',mode,'_',name,'.csv')))
    # Identical algorithm/settings and the assay-specific seed used by the baseline.
    set.seed(20260909L+i)
    cat('Running',mode,name,'with',nrow(selected),'genes\n');flush.console()
    res<-withCallingHandlers(fgseaMultilevel(pathways,values,minSize=15,maxSize=500,eps=0,sampleSize=101,nPermSimple=5000,scoreType='std',gseaParam=1,BPPARAM=SerialParam()),
      warning=function(w){warnings_seen<<-c(warnings_seen,paste(mode,name,conditionMessage(w)));invokeRestart('muffleWarning')})
    res<-merge(as.data.table(res),annotation,by='pathway')
    res[,padj_global:=p.adjust(pval,'BH')];res[,padj:=p.adjust(pval,'BH'),by=collection]
    res[,`:=`(analysis=name,exclusion=mode)]
    symbols<-setNames(ranks$symbol,ranks$entrez_id)
    res[,leading_edge_entrez:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))]
    res[,leading_edge_symbols:=vapply(leadingEdge,function(g)paste(symbols[g],collapse=';'),character(1))]
    res[,leadingEdge:=NULL]
    fwrite(res,file.path(OUT,paste0('fgsea_',mode,'_',name,'.csv')))
    outputs[[paste(mode,name)]]<-res
    audit[[paste(mode,name)]]<-data.table(exclusion=mode,analysis=name,genes=nrow(selected),excluded_genes=length(drop),pathways=nrow(res),significant=sum(res$padj<.05,na.rm=TRUE))
  }
}
fwrite(rbindlist(outputs),file.path(OUT,'fgsea_exclusion_results.csv'));fwrite(rbindlist(audit),file.path(OUT,'run_audit.csv'))
writeLines(warnings_seen,file.path(OUT,'fgsea_warnings.txt'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
cat('Ribosomal exclusion GSEA complete.\n')
