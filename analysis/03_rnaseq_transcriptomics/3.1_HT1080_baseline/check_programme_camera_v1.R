options(warn=1,stringsAsFactors=FALSE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(normalizePath(script,winslash='/',mustWork=TRUE)))
suppressPackageStartupMessages({library(limma);library(data.table)})
setDTthreads(2)
out <- 'outputs/programme_analysis_v1'
b <- readRDS('outputs/models/batch_aware_DE.rds')$fits$batch_condition
e <- readRDS('outputs/models/02_expression.rds')
mp <- e$mapping;v <- b$voom[mp$gene_id,];rownames(v) <- mp$gene_name
sets <- do.call(c,unname(readRDS('resources/msigdb_all_collections.rds')))
coverage <- fread(file.path(out,'tables/GSEA_coverage.csv'))
sets <- sets[coverage$pathway[coverage$contrast=='overall' & coverage$eligible]]
ix <- ids2indices(sets,rownames(v),remove.empty=FALSE)
stopifnot(length(ix)==4868,all(lengths(ix)>=15 & lengths(ix)<=500))
rows <- list()
for(mode in c('fixed_0.01','estimated_nonnegative')) {
 message('CAMERA overall: ',mode)
 rho <- if(mode=='fixed_0.01') .01 else NA_real_
 z <- camera(v,ix,design=b$design,contrast=b$contrasts[,1],inter.gene.cor=rho,
             allow.neg.cor=FALSE,sort=FALSE)
 z$pathway <- rownames(z);z$correlation_mode <- mode
 rows[[mode]] <- z
}
fwrite(rbindlist(rows,fill=TRUE),file.path(out,'tables/CAMERA_overall_sensitivity.csv'))
print(rbindlist(lapply(rows,function(z) data.frame(mode=z$correlation_mode[1],FDR05=sum(z$FDR<.05)))))
