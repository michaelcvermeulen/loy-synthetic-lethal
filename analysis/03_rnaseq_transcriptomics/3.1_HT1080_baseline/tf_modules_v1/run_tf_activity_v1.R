# Stages 1-2: multi-method, multi-prior TF target-pattern scores and the pre-registered robustness rule.
# Rules: provenance/DECISIONS.md (24 Sep 2026). Original outputs are read-only and hash-checked.
options(warn=1,stringsAsFactors=FALSE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(file.path(dirname(normalizePath(script,winslash='/',mustWork=TRUE)),'..'))
suppressPackageStartupMessages({library(limma);library(edgeR);library(data.table);library(decoupleR)})
source('R/00_common.R')
setDTthreads(2);set.seed(42)
out <- 'outputs/tf_modules_v1';res <- 'resources/tf_modules_v1'
put <- function(z,n) fwrite(as.data.frame(z),file.path(out,'tables',paste0(n,'.csv')),na='NA')
MIN <- 5;WSUM_TIMES <- 1000
TEST <- Sys.getenv('TFMOD_TEST')=='1'   # quick end-to-end check; outputs go to a test folder
if(TEST) {WSUM_TIMES <- 5;out <- file.path(out,'_test');for(d in c('tables','models','provenance')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)}
inputs <- c('outputs/models/01_ingest.rds','outputs/models/02_expression.rds','outputs/models/batch_aware_DE.rds',
 'outputs/programme_analysis_v1/models/DE.rds','outputs/programme_analysis_v1/tables/TF_activity_all.csv',
 file.path(res,c('collectri_current.csv','collectri_frozen_v1.csv','dorothea_ABC.csv','tfactprofiler_all.csv','tfactprofiler_fibroblast.csv')))
hashes <- vapply(inputs,hash_file,character(1))
fwrite(data.frame(path=inputs,sha256=hashes),file.path(out,'provenance/stage1_input_hashes.csv'))

x <- readRDS(inputs[1]);e <- readRDS(inputs[2]);b <- readRDS(inputs[3]);de <- readRDS(inputs[4])
readnet <- function(p) {n <- fread(p);as.data.frame(n)}
nets <- list(collectri=readnet(inputs[6]),dorothea_ABC=readnet(inputs[8]),
 tfactprofiler_all=readnet(inputs[9]),tfactprofiler_fibroblast=readnet(inputs[10]))
frozen <- readnet(inputs[7])
if(TEST) for(p in c('tfactprofiler_all','tfactprofiler_fibroblast')) nets[[p]] <- nets[[p]][nets[[p]]$source %in% unique(nets[[p]]$source)[1:40],]
mat1 <- function(v,name) matrix(v,ncol=1,dimnames=list(names(v),name))
score <- function(fun,v,net,...) {
 r <- as.data.frame(fun(mat=mat1(v,'x'),network=net,.source='source',.target='target',.mor='mor',minsize=MIN,...))
 r[,c('statistic','source','score','p_value')]
}

# ---- bridge: frozen prior + logFC + ULM reproduces v1 exactly ---------------------------
old <- fread(inputs[5])
bridge <- rbindlist(lapply(names(de$effects),function(k) {
 r <- score(run_ulm,de$effects[[k]],frozen);o <- old[old$contrast==k,]
 m <- merge(r,o[,c('source','score','p_value')],by='source',suffixes=c('_new','_v1'))
 data.frame(contrast=k,n_new=nrow(r),n_v1=nrow(o),matched=nrow(m),max_abs_score_diff=max(abs(m$score_new-m$score_v1)),
  max_abs_p_diff=max(abs(m$p_value_new-m$p_value_v1)))
}))
put(bridge,'TF_bridge_to_v1');print(bridge)
stopifnot(all(bridge$n_new==bridge$n_v1),all(bridge$max_abs_score_diff<1e-10))

# ---- stage 1: scores ------------------------------------------------------------------------
contrasts <- names(de$rankings)
methods <- list(ulm=run_ulm,mlm=run_mlm,wsum=run_wsum,viper=run_viper)
ckpt <- file.path(out,'models/stage1_scores_checkpoint.rds')
all_scores <- if(file.exists(ckpt)) readRDS(ckpt) else list()
for(k in contrasts) for(p in names(nets)) {
 # TFActProfiler priors enter only R4 (overall) and the period comparison; skip them elsewhere for runtime.
 if(grepl('^tfactprofiler',p) && !k %in% c('overall','AUG_2025','MAY_2026')) next
 use <- if(p %in% c('collectri','dorothea_ABC')) names(methods) else 'ulm'
 use <- use[!paste(k,p,use) %in% names(all_scores)]
 for(mname in use) {
  message(k,' / ',p,' / ',mname)
  args <- switch(mname,wsum=list(times=WSUM_TIMES,seed=42),viper=list(pleiotropy=FALSE,verbose=FALSE),list())
  r <- do.call(score,c(list(methods[[mname]],de$rankings[[k]],nets[[p]]),args))
  if(mname=='wsum') r <- r[r$statistic=='norm_wsum',]
  r$method <- mname;r$prior <- p;r$contrast <- k;r$input <- 't'
  all_scores[[paste(k,p,mname)]] <- r
 }
 saveRDS(all_scores,ckpt)
}
# Secondary logFC input (primary prior, ULM) for the input-choice comparison.
for(k in contrasts) {
 r <- score(run_ulm,de$effects[[k]],nets$collectri);r$method <- 'ulm';r$prior <- 'collectri';r$contrast <- k;r$input <- 'logFC'
 all_scores[[paste(k,'logFC')]] <- r
}
S <- rbindlist(all_scores,fill=TRUE)
S[,FDR:=p.adjust(p_value,'BH'),by=.(contrast,prior,method,input)]
S[,signed_z:=sign(score)*qnorm(1-pmax(p_value,1e-300)/2)]
size <- rbindlist(lapply(names(nets),function(p) {
 n <- as.data.table(nets[[p]])[target %in% names(de$rankings$overall)]
 n[,.(targets_covered=.N),by=source][,prior:=p]
}))
S <- merge(S,size,by=c('prior','source'),all.x=TRUE)
put(S,'TF_activity_matrix')

# ---- fry flag (gene correlation-aware), batch + condition voom --------------------------------
bc <- b$fits$batch_condition
mp <- e$mapping[e$mapping$representative,c('gene_id','gene_name','arm','chrom','start')]
cm <- nets$collectri[nets$collectri$target %in% mp$gene_name,]
cm$gene_id <- mp$gene_id[match(cm$target,mp$gene_name)]
idx <- split(cm[,c('gene_id','mor')],cm$source);idx <- idx[vapply(idx,nrow,1L)>=MIN]
idx <- lapply(idx,function(d) data.frame(gene_id=d$gene_id,weight=d$mor))
fr <- fry(bc$voom,index=idx,design=bc$design,contrast=bc$contrasts[,'combined'],sort='none')
fr <- data.frame(source=rownames(fr),fr,row.names=NULL);put(fr,'TF_fry_collectri_overall')

# ---- leave-one-passage-pair-out refits of the batch + condition model -------------------------
m <- x$core;a <- e$effects[e$effects$primary_eligible,]
effective <- setNames(e$norm$effective,colnames(e$logCPM))[m$sample_key]
counts <- x$counts[a$gene_id,m$sample_key]
fit_bc <- function(keep) {
 mm <- m[keep,];mm$batch <- factor(mm$batch,levels=c('AUG_2025','MAY_2026'));mm$treatment <- factor(mm$treatment,levels=c('CTRL','LOY'))
 D <- model.matrix(~batch+treatment,mm);rownames(D) <- mm$sample_key
 v <- voom(counts[,keep],design=D,lib.size=effective[keep],normalize.method='none',plot=FALSE)
 f <- eBayes(lmFit(v,D),trend=FALSE,robust=TRUE)
 t <- f$t[,'treatmentLOY'];t <- t[mp$gene_id];names(t) <- mp$gene_name;t
}
full <- fit_bc(rep(TRUE,nrow(m)))
chk <- max(abs(full-de$rankings$overall[names(full)]))
stopifnot(chk<1e-8)   # refit reproduces the saved overall moderated t
lopo <- rbindlist(lapply(seq_len(nrow(x$pairs)),function(i) {
 drop <- c(x$pairs$CTRL[i],x$pairs$LOY[i]);r <- score(run_ulm,fit_bc(!m$sample_key %in% drop),nets$collectri)
 r$omitted_block <- x$pairs$block[i];r
}))
put(lopo,'TF_leave_one_pair_out_ulm')

# ---- stage 2: rule ----------------------------------------------------------------------------
g <- function(k,p='collectri',mt='ulm') S[contrast==k & prior==p & method==mt & input=='t',.(source,score,p_value,FDR)]
R <- g('overall');setnames(R,c('score','p_value','FDR'),c('ulm_score','ulm_p','ulm_FDR'))
for(k in c('AUG_2025','MAY_2026','interaction','CTRL_vs_parental')) {
 z <- g(k);setnames(z,c('score','p_value','FDR'),paste0(k,c('_score','_p','_FDR')));R <- merge(R,z,by='source',all.x=TRUE)
}
R$R1 <- R$ulm_FDR<.05
R$R2 <- sign(R$AUG_2025_score)==sign(R$ulm_score) & sign(R$MAY_2026_score)==sign(R$ulm_score)
mt <- dcast(S[contrast=='overall' & prior=='collectri' & input=='t'],source~method,value.var=c('score','p_value'))
R <- merge(R,mt,by='source',all.x=TRUE)
agree <- function(s,p) !is.na(s) & sign(s)==sign(R$ulm_score) & p<.05
R$methods_agree <- agree(R$score_ulm,R$p_value_ulm)+agree(R$score_mlm,R$p_value_mlm)+agree(R$score_wsum,R$p_value_wsum)+agree(R$score_viper,R$p_value_viper)
R$R3 <- R$methods_agree>=3
R$consensus_z <- S[contrast=='overall' & prior=='collectri' & input=='t',.(consensus_z=mean(signed_z)),by=source][match(R$source,source),consensus_z]
for(p in c('dorothea_ABC','tfactprofiler_all','tfactprofiler_fibroblast')) {
 z <- S[contrast=='overall' & prior==p & method=='ulm' & input=='t'];R[[paste0(p,'_score')]] <- z$score[match(R$source,z$source)]
 R[[paste0(p,'_p')]] <- z$p_value[match(R$source,z$source)]
}
alt_ok <- function(s) is.na(s) | sign(s)==sign(R$ulm_score)
R$alt_priors_covering <- (!is.na(R$dorothea_ABC_score))+(!is.na(R$tfactprofiler_all_score))
R$R4 <- alt_ok(R$dorothea_ABC_score) & alt_ok(R$tfactprofiler_all_score)
L <- lopo[,.(lopo_same_sign=sum(sign(score)==sign(R$ulm_score[match(source,R$source)])),lopo_min_abs=min(abs(score))),by=source]
R <- merge(R,L,by='source',all.x=TRUE);R$R5 <- R$lopo_same_sign==6

# R6: remove self/complex-subunit edges and the three largest contributors; rescore that TF alone.
complex_sub <- list(AP1=c('JUN','JUNB','JUND','FOS','FOSB','FOSL1','FOSL2'),NFKB=c('NFKB1','NFKB2','RELA','RELB','REL'))
tt <- de$rankings$overall;net <- nets$collectri[nets$collectri$target %in% names(tt),]
contrib <- rbindlist(lapply(R$source[R$R1 %in% TRUE],function(s) {
 n <- net[net$source==s,];sg <- sign(R$ulm_score[R$source==s])
 data.frame(source=s,target=n$target,mor=n$mor,t=tt[n$target],contribution=n$mor*tt[n$target]*sg,
  self_edge=n$target==s | n$target %in% complex_sub[[s]])
}))
contrib <- contrib[order(contrib$source,-contrib$contribution),]
omit <- rbindlist(lapply(unique(contrib$source),function(s) {
 cc <- contrib[contrib$source==s,];drop <- c(cc$target[cc$self_edge],head(cc$target[!cc$self_edge],3))
 n <- net[net$source==s & !net$target %in% drop,]
 r <- if(nrow(n)>=MIN) score(run_ulm,tt,n) else data.frame(statistic='ulm',source=s,score=NA,p_value=NA)
 data.frame(source=s,omitted=paste(drop,collapse=';'),n_omitted=length(drop),remaining=nrow(n),omit_score=r$score,omit_p=r$p_value)
}))
R <- merge(R,omit,by='source',all.x=TRUE)
R$R6 <- !is.na(R$omit_score) & sign(R$omit_score)==sign(R$ulm_score) & R$omit_p<.05
R$fry_FDR <- fr$FDR[match(R$source,fr$source)];R$fry_direction <- fr$Direction[match(R$source,fr$source)]
R$targets_covered <- size[prior=='collectri'][match(R$source,source),targets_covered]
R$self_edge_present <- R$source %in% contrib$source[contrib$self_edge]
rr <- c('R1','R2','R3','R4','R5','R6')
for(r in rr) R[[r]][is.na(R[[r]])] <- FALSE
R$n_R3_R6 <- R$R3+R$R4+R$R5+R$R6
R$tier <- ifelse(R$R1 & R$R2 & R$R3 & R$R4 & R$R5 & R$R6,'robust',ifelse(R$R1 & R$R2 & R$n_R3_R6>=3,'supported',
 ifelse(R$R1 & R$R2,'consistent_only',ifelse(R$R1,'overall_only','not_significant'))))
R$both_periods_FDR05 <- R$R1 & R$R2 & R$AUG_2025_FDR<.05 & R$MAY_2026_FDR<.05
R$parental_drift_flag <- R$tier %in% c('robust','supported') & sign(R$CTRL_vs_parental_score)==sign(R$ulm_score) & R$CTRL_vs_parental_p<.05
R <- R[order(match(R$tier,c('robust','supported','consistent_only','overall_only','not_significant')),-abs(R$consensus_z)),]
put(R,'TF_robustness');put(contrib,'TF_target_contributions_R1')
print(table(R$tier));print(R[R$tier %in% c('robust','supported'),c('source','ulm_score','ulm_FDR','methods_agree','R4','lopo_same_sign','R6','fry_FDR','tier')])
saveRDS(list(scores=S,robustness=R,lopo=lopo,fry=fr,contrib=contrib,bridge=bridge),file.path(out,'models/tf_activity_v1.rds'))
stopifnot(identical(hashes,vapply(inputs,hash_file,character(1))))
writeLines(capture.output(sessionInfo()),file.path(out,'provenance/stage1_sessionInfo.txt'))
cat('Stage 1-2 completed\n')
