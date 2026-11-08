options(warn=1,stringsAsFactors=FALSE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(normalizePath(script,winslash='/',mustWork=TRUE)))
suppressPackageStartupMessages({library(data.table);library(fgsea)})
setDTthreads(2)
source('R/00_common.R')
base <- 'outputs/programme_analysis_v1'
out <- file.path(base,'all_collections');dir.create(out,recursive=TRUE,showWarnings=FALSE)
inputs <- c(file.path(base,'models/DE.rds'),'resources/msigdb_all_collections.rds','resources/msigdb_all_set_metadata.csv')
hashes <- data.frame(path=inputs,sha256=vapply(inputs,hash_file,character(1)))
hp <- file.path(out,'input_hashes.csv')
if(file.exists(hp)) {
 old <- fread(hp);stopifnot(identical(old$path,hashes$path),identical(old$sha256,hashes$sha256))
} else fwrite(hashes,hp)
rankings <- readRDS(inputs[1])$rankings
sets <- do.call(c,unname(readRDS(inputs[2])))
meta <- as.data.frame(fread(inputs[3]))
names(meta)[match(c('gs_name','gs_collection','gs_subcollection'),names(meta))] <- c('pathway','collection','subcollection')
stopifnot(!anyDuplicated(meta$pathway),all(meta$pathway %in% names(sets)))
res <- list();coverage <- list()
for(k in names(rankings)) {
 stats <- rankings[[k]];stats <- stats[order(-stats,names(stats))]
 universe <- names(stats)
 size <- vapply(sets,function(s)sum(!is.na(fastmatch::fmatch(unique(s),universe))),integer(1))
 ok <- size>=15 & size<=500
 coverage[[k]] <- data.frame(contrast=k,pathway=names(sets),overlap_size=size,eligible=ok)
 path <- file.path(out,paste0('GSEA_',k,'.csv'))
 if(file.exists(path)) {
  message('Resume completed contrast: ',k)
  res[[k]] <- as.data.frame(fread(path));stopifnot(nrow(res[[k]])==sum(ok))
  next
 }
 message('All collections: ',k,' / ',sum(ok),' eligible sets')
 set.seed(42+match(k,names(rankings)))
 warnings <- character()
 z <- withCallingHandlers(as.data.frame(fgseaMultilevel(pathways=sets[ok],stats=stats,minSize=15,maxSize=500,
  eps=1e-10,nPermSimple=10000,BPPARAM=BiocParallel::SerialParam())),warning=function(w){warnings<<-c(warnings,conditionMessage(w))})
 z$padj <- NULL;z$FDR <- p.adjust(z$pval,'BH',n=nrow(z));z$contrast <- k
 z$leadingEdge <- vapply(z$leadingEdge,paste,character(1),collapse=';')
 z$test_status <- ifelse(is.finite(z$pval),'tested','numerical_failure')
 z$p_at_floor <- is.finite(z$pval)&z$pval<=1e-10
 z <- merge(z,meta,by='pathway',sort=FALSE)
 z$FDR_within_collection <- NA_real_
 for(c in unique(z$collection)) {j <- z$collection==c;z$FDR_within_collection[j] <- p.adjust(z$pval[j],'BH',n=sum(j))}
 fwrite(z,path);writeLines(warnings,file.path(out,paste0('warnings_',k,'.txt')));res[[k]] <- z
}
fwrite(rbindlist(coverage),file.path(out,'coverage.csv'))
fwrite(rbindlist(res,fill=TRUE),file.path(out,'GSEA_all_contrasts.csv'))
gc <- res$overall
for(k in c('AUG_2025','MAY_2026')) {
 r <- res[[k]][match(gc$pathway,res[[k]]$pathway),]
 gc[[paste0(k,'_NES')]] <- r$NES;gc[[paste0(k,'_FDR')]] <- r$FDR
 gc[[paste0(k,'_leadingEdge')]] <- r$leadingEdge
}
gc$same_direction <- sign(gc$NES)==sign(gc$AUG_2025_NES)&sign(gc$NES)==sign(gc$MAY_2026_NES)
gc$directionally_consistent <- !is.na(gc$FDR)&gc$FDR<.05&gc$same_direction
gc$both_periods_FDR05 <- gc$directionally_consistent&gc$AUG_2025_FDR<.05&gc$MAY_2026_FDR<.05
fwrite(gc,file.path(out,'programme_consistency.csv'))
stopifnot(identical(hashes$sha256,unname(vapply(inputs,hash_file,character(1)))))
writeLines(c('Global BH across all eligible collections separately per contrast. Within-collection FDR is secondary only.',
 'Current moderated t rankings; autosomal/non-PAR X background. All available frozen MSigDB collections; 15-500 represented symbols.',
 capture.output(sessionInfo())),file.path(out,'METHODS.txt'))
print(data.frame(tested=nrow(gc),overall_FDR05=sum(gc$FDR<.05,na.rm=TRUE),consistent=sum(gc$directionally_consistent,na.rm=TRUE),both=sum(gc$both_periods_FDR05,na.rm=TRUE)))
