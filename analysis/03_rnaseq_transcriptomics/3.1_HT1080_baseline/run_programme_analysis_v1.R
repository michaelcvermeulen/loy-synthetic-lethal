options(warn=1,stringsAsFactors=FALSE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(normalizePath(script,winslash='/',mustWork=TRUE)))
suppressPackageStartupMessages({library(limma);library(edgeR);library(data.table)})
source('R/00_common.R');source('R/02_effects.R')
setDTthreads(2);set.seed(42)
out <- 'outputs/programme_analysis_v1'
for(d in c('tables','models','provenance','figures')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
put <- function(z,n) fwrite(as.data.frame(z),file.path(out,'tables',paste0(n,'.csv')),na='NA')
inputs <- c('outputs/models/01_ingest.rds','outputs/models/02_expression.rds',
 'outputs/tables/batchDE_batch_condition__combined.csv','outputs/models/batch_aware_DE.rds',
 'resources/msigdb_all_collections.rds','resources/msigdb_all_set_metadata.csv','resources/collectri_network.csv')
hashes <- vapply(inputs,hash_file,character(1))
fwrite(data.frame(path=inputs,sha256=hashes),file.path(out,'provenance/input_hashes.csv'))
x <- readRDS(inputs[1]);e <- readRDS(inputs[2]);cfg <- yaml::read_yaml('config/analysis.yml')
a <- e$effects[e$effects$primary_eligible,];m <- x$core
effective <- setNames(e$norm$effective,colnames(e$logCPM))[m$sample_key]
m$batch <- factor(m$batch,levels=c('AUG_2025','MAY_2026'))
m$treatment <- factor(m$treatment,levels=c('CTRL','LOY'))
D <- model.matrix(~batch*treatment,m);rownames(D) <- m$sample_key
stopifnot(nrow(a)==25526,qr(D)$rank==4,nrow(D)==12)
v <- voom(x$counts[a$gene_id,m$sample_key],design=D,lib.size=effective,normalize.method='none',plot=FALSE)
C <- matrix(0,ncol(D),3,dimnames=list(colnames(D),c('AUG_2025','MAY_2026','interaction')))
C['treatmentLOY',c('AUG_2025','MAY_2026')] <- 1
int <- grep(':',colnames(D),value=TRUE)
C[int,c('MAY_2026','interaction')] <- 1
f <- eBayes(contrasts.fit(lmFit(v,D),C),robust=TRUE,trend=FALSE)
extract <- function(f,coef,a) {
 z <- topTable(f,coef=coef,number=Inf,sort.by='none',confint=TRUE)
 z <- z[match(a$gene_id,rownames(z)),]
 names(z)[names(z)=='P.Value'] <- 'PValue';names(z)[names(z)=='adj.P.Val'] <- 'FDR'
 stopifnot(all(is.finite(z$FDR)),max(abs(z$FDR-p.adjust(z$PValue,'BH')))<1e-12)
 cbind(a[,c('gene_id','gene_name','gene_type','chrom','start','end')],z)
}
results <- list(overall=as.data.frame(fread(inputs[3])))
for(k in colnames(C)) results[[k]] <- extract(f,k,a)
stopifnot(max(abs(results$interaction$logFC-(results$MAY_2026$logFC-results$AUG_2025$logFC)))<1e-10)
z <- results$overall
z$AUG_logFC <- results$AUG_2025$logFC;z$MAY_logFC <- results$MAY_2026$logFC
z$AUG_FDR <- results$AUG_2025$FDR;z$MAY_FDR <- results$MAY_2026$FDR
z$interaction_logFC <- results$interaction$logFC;z$interaction_FDR <- results$interaction$FDR
z$passes_main_cutoffs <- z$FDR<.05 & abs(z$logFC)>=1
z$same_direction <- sign(z$logFC)==sign(z$AUG_logFC) & sign(z$logFC)==sign(z$MAY_logFC)
z$directionally_consistent <- z$passes_main_cutoffs & z$same_direction
z$both_periods_FDR05 <- z$AUG_FDR<.05 & z$MAY_FDR<.05 & z$same_direction
put(z,'gene_consistency');put(z[z$directionally_consistent,],'directionally_consistent_genes')
put(data.frame(sample_key=m$sample_key,D,check.names=FALSE),'interaction_design')
put(data.frame(coefficient=rownames(C),C),'interaction_contrasts')
# Separately filtered and normalized August parental context; do not assume pairing.
pm <- x$meta[x$meta$batch=='AUG_2025' & x$meta$treatment %in% c('CTRL','PARENTAL'),]
pm$treatment <- factor(pm$treatment,levels=c('PARENTAL','CTRL'))
pn <- normalize_counts(x$counts[,pm$sample_key],x$annotation$main_background,pm$treatment,cfg)
pa <- x$annotation[pn$keep,];PD <- model.matrix(~treatment,pm);rownames(PD) <- pm$sample_key
stopifnot(nrow(pm)==6,all(table(pm$treatment)==3),qr(PD)$rank==2)
pv <- voom(x$counts[pa$gene_id,pm$sample_key],design=PD,lib.size=pn$effective,normalize.method='none',plot=FALSE)
pf <- eBayes(lmFit(pv,PD),robust=TRUE,trend=FALSE)
results$CTRL_vs_parental <- extract(pf,'treatmentCTRL',pa)
put(data.frame(sample_key=pm$sample_key,PD,effective_library_reads=pn$effective,check.names=FALSE),'parental_design_normalization')
for(k in names(results)) put(results[[k]],paste0('DE_',k))
summary <- rbindlist(lapply(names(results),function(k) {
 r <- results[[k]];s <- r$FDR<.05 & abs(r$logFC)>=1
 data.frame(contrast=k,tested=nrow(r),FDR05=sum(r$FDR<.05),FDR05_LFC1=sum(s),higher=sum(s&r$logFC>0),lower=sum(s&r$logFC<0))
}))
put(summary,'DE_summary')
overlap <- merge(z,results$CTRL_vs_parental[,c('gene_id','logFC','FDR')],by='gene_id',suffixes=c('','_CTRL_parental'),all.x=TRUE)
overlap$parental_tested <- !is.na(overlap$FDR_CTRL_parental)
overlap$parental_passes_cutoffs <- overlap$parental_tested & overlap$FDR_CTRL_parental<.05 & abs(overlap$logFC_CTRL_parental)>=1
overlap$same_direction_as_CTRL_parental <- sign(overlap$logFC)==sign(overlap$logFC_CTRL_parental)
put(overlap,'parental_context_overlap')
# Freeze unique-symbol mapping, independent of DE significance or test statistics.
mapping <- e$mapping
parent_mapping <- choose_symbols(x$annotation,rowMeans(pn$logCPM),pn$keep)
parent_mapping <- parent_mapping[parent_mapping$representative,]
rankings <- list();effects <- list()
for(k in names(results)) {
 mp <- if(k=='CTRL_vs_parental') parent_mapping else mapping
 rr <- results[[k]][match(mp$gene_id,results[[k]]$gene_id),]
 stopifnot(!anyNA(rr$gene_id),!anyDuplicated(mp$gene_name),all(is.finite(rr$t)))
 rankings[[k]] <- setNames(rr$t,mp$gene_name)
 effects[[k]] <- setNames(rr$logFC,mp$gene_name)
 put(data.frame(gene_id=mp$gene_id,gene_name=mp$gene_name,moderated_t=rr$t,logFC=rr$logFC),paste0('ranking_',k))
}
saveRDS(list(results=results,consistency=z,interaction=list(voom=v,fit=f,design=D,contrasts=C),
 parental=list(voom=pv,fit=pf,design=PD,norm=pn),rankings=rankings,effects=effects),file.path(out,'models/DE.rds'))
print(summary);cat('Directionally consistent:',sum(z$directionally_consistent),'\n')

# Equal primary status for the three chosen collections; one BH family per contrast.
allsets <- do.call(c,unname(readRDS('resources/msigdb_all_collections.rds')))
meta <- as.data.frame(fread('resources/msigdb_all_set_metadata.csv'))
names(meta)[match(c('gs_name','gs_collection','gs_subcollection'),names(meta))] <- c('pathway','collection','subcollection')
meta <- meta[meta$collection=='H' | meta$subcollection %in% c('CP:REACTOME','GO:BP'),]
meta$analysis_collection <- ifelse(meta$collection=='H','Hallmark',ifelse(meta$subcollection=='CP:REACTOME','Reactome','GO Biological Process'))
enrich <- list();coverage <- list()
for(k in names(rankings)) {
 stats <- rankings[[k]];stats <- stats[order(-stats,names(stats))]
 sets <- allsets[meta$pathway]
 size <- vapply(sets,function(s) length(intersect(s,names(stats))),integer(1))
 ok <- size>=15 & size<=500
 coverage[[k]] <- data.frame(contrast=k,pathway=names(sets),overlap_size=size,eligible=ok)
 message('Primary fgsea: ',k,' / ',sum(ok),' sets')
 set.seed(42+match(k,names(rankings)))
 g <- as.data.frame(fgsea::fgseaMultilevel(pathways=sets[ok],stats=stats,minSize=15,maxSize=500,
       eps=1e-10,nPermSimple=10000,BPPARAM=BiocParallel::SerialParam()))
 g$FDR <- p.adjust(g$pval,'BH',n=nrow(g));g$padj <- NULL;g$contrast <- k
 g$leadingEdge <- vapply(g$leadingEdge,paste,character(1),collapse=';')
 g$test_status <- ifelse(is.finite(g$pval),'tested','numerical_failure')
 enrich[[k]] <- merge(g,meta,by='pathway',sort=FALSE)
 put(enrich[[k]],paste0('GSEA_',k))
}
g <- rbindlist(enrich,fill=TRUE);put(g,'GSEA_primary_all');put(rbindlist(coverage),'GSEA_coverage')
gc <- enrich$overall
for(k in c('AUG_2025','MAY_2026')) {
 r <- enrich[[k]][match(gc$pathway,enrich[[k]]$pathway),]
 gc[[paste0(k,'_NES')]] <- r$NES;gc[[paste0(k,'_FDR')]] <- r$FDR
 gc[[paste0(k,'_leadingEdge')]] <- r$leadingEdge
}
gc$same_direction <- sign(gc$NES)==sign(gc$AUG_2025_NES) & sign(gc$NES)==sign(gc$MAY_2026_NES)
gc$directionally_consistent <- gc$FDR<.05 & gc$same_direction
gc$both_periods_FDR05 <- gc$directionally_consistent & gc$AUG_2025_FDR<.05 & gc$MAY_2026_FDR<.05
gc$period_leading_edge_jaccard <- mapply(function(a,b) {
 aa <- strsplit(a,';',fixed=TRUE)[[1]];bb <- strsplit(b,';',fixed=TRUE)[[1]]
 length(intersect(aa,bb))/length(union(aa,bb))
},gc$AUG_2025_leadingEdge,gc$MAY_2026_leadingEdge)
put(gc,'programme_consistency')

# TF activity from signed target-expression changes, not TF transcript abundance.
net <- unique(as.data.frame(fread('resources/collectri_network.csv'))[,c('source','target','mor')])
key <- paste(net$source,net$target,sep='\r');bad <- duplicated(key)|duplicated(key,fromLast=TRUE)
put(net[bad,],'TF_excluded_conflicting_edges');net <- net[!bad & is.finite(net$mor),]
tf <- list()
for(k in names(effects)) {
 mat <- matrix(effects[[k]],ncol=1,dimnames=list(names(effects[[k]]),k))
 r <- as.data.frame(decoupleR::run_ulm(mat,net,.source='source',.target='target',.mor='mor',minsize=5))
 r$FDR <- p.adjust(r$p_value,'BH');r$contrast <- k
 r$statistical_scope <- 'conditional_target_pattern_not_independent_population_replication'
 tf[[k]] <- r
}
put(rbindlist(tf,fill=TRUE),'TF_activity_all')
tc <- tf$overall
for(k in c('AUG_2025','MAY_2026')) {
 r <- tf[[k]][match(tc$source,tf[[k]]$source),]
 tc[[paste0(k,'score')]] <- r$score;tc[[paste0(k,'FDR')]] <- r$FDR
}
tc$same_direction <- sign(tc$score)==sign(tc$AUG_2025score) & sign(tc$score)==sign(tc$MAY_2026score)
tc$directionally_consistent <- tc$FDR<.05 & tc$same_direction
put(tc,'TF_consistency')
saveRDS(list(enrichment=enrich,programme_consistency=gc,tf=tf,tf_consistency=tc),file.path(out,'models/programmes_TF.rds'))
stopifnot(identical(hashes,vapply(inputs,hash_file,character(1))))
writeLines(c('Validation passed: unchanged original inputs; full-rank models; exact interaction contrast; BH checks; original core gene universe and exposure offsets retained.',capture.output(sessionInfo())),file.path(out,'provenance/R_validation.txt'))
cat('RNA programme analysis completed.\n')
