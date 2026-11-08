# Stage 0: freeze TF-target priors for tf_modules_v1 (rules in provenance/DECISIONS.md, 24 Sep 2026).
options(warn=1,stringsAsFactors=FALSE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(file.path(dirname(normalizePath(script,winslash='/',mustWork=TRUE)),'..'))
suppressPackageStartupMessages({library(data.table)})
source('R/00_common.R')
res <- 'resources/tf_modules_v1';dir.create(res,showWarnings=FALSE,recursive=TRUE)
TFAP_COMMIT <- '213e0d73e3'
tfap_url <- function(f) paste0('https://raw.githubusercontent.com/HikaruSugimoto/tfactprofiler/',TFAP_COMMIT,'/example/',utils::URLencode(f))

clean <- function(net,name) {
 net <- unique(as.data.frame(net)[,c('source','target','mor')])
 net <- net[is.finite(net$mor) & net$mor!=0 & nzchar(net$source) & nzchar(net$target),]
 key <- paste(net$source,net$target,sep='\r');bad <- duplicated(key)|duplicated(key,fromLast=TRUE)
 if(any(bad)) fwrite(net[bad,],file.path(res,paste0(name,'_excluded_conflicting_edges.csv')))
 net <- net[!bad,];net$self_edge <- net$source==net$target
 fwrite(net,file.path(res,paste0(name,'.csv')))
 data.frame(prior=name,edges=nrow(net),sources=length(unique(net$source)),targets=length(unique(net$target)),
  conflicting_removed=sum(bad),self_edges=sum(net$self_edge),positive=sum(net$mor>0),negative=sum(net$mor<0))
}
summ <- list()
frozen <- fread('resources/collectri_network.csv')
summ$frozen <- clean(frozen,'collectri_frozen_v1')
cur <- decoupleR::get_collectri(organism='human',split_complexes=FALSE)
summ$current <- clean(cur,'collectri_current')
fk <- with(fread(file.path(res,'collectri_frozen_v1.csv')),paste(source,target,mor))
ck <- with(fread(file.path(res,'collectri_current.csv')),paste(source,target,mor))
edge_diff <- data.frame(only_frozen=length(setdiff(fk,ck)),only_current=length(setdiff(ck,fk)),shared=length(intersect(fk,ck)))
fwrite(edge_diff,file.path(res,'collectri_frozen_vs_current_diff.csv'))
d <- dorothea::dorothea_hs;d <- d[d$confidence %in% c('A','B','C'),]
summ$dorothea <- clean(data.frame(source=d$tf,target=d$target,mor=d$mor),'dorothea_ABC')
for(f in c('all','fibroblast')) {
 file <- if(f=='all') 'TFactprofiler_all.csv' else 'TFactprofiler_organ/TFactprofier_fibroblast.csv'
 raw <- file.path(res,paste0('tfactprofiler_',f,'_raw.csv'))
 if(!file.exists(raw)) download.file(tfap_url(file),raw,mode='wb',quiet=TRUE)
 t <- fread(raw);stopifnot(all(c('source','target','weight') %in% names(t)))
 summ[[f]] <- clean(data.frame(source=t$source,target=t$target,mor=t$weight),paste0('tfactprofiler_',f))
}
s <- rbindlist(summ);fwrite(s,file.path(res,'prior_summary.csv'));print(s);print(edge_diff)
files <- list.files(res,full.names=TRUE,pattern='\\.csv$')
man <- data.frame(file=basename(files),sha256=vapply(files,hash_file,character(1)),row.names=NULL)
fwrite(man,file.path(res,'prior_hashes.csv'))
writeLines(c(paste('retrieved',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')),
 paste('OmnipathR',as.character(packageVersion('OmnipathR'))),paste('decoupleR',as.character(packageVersion('decoupleR'))),
 paste('dorothea',as.character(packageVersion('dorothea'))),paste('tfactprofiler commit',TFAP_COMMIT),
 'CollecTRI retrieved with get_collectri(split_complexes=FALSE); OmnipathR may serve from its local cache.',
 capture.output(sessionInfo())),file.path(res,'prior_provenance.txt'))
