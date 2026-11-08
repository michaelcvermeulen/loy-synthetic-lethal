# Stages 1-2 (A549): multi-method, multi-prior TF target-pattern scores and the robustness rule adapted from
# HT1080 tf_modules_v1. Rules: provenance/DECISIONS.md (27 Sep 2026). Original outputs are read-only and hash-checked.
options(warn=1,stringsAsFactors=FALSE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(file.path(dirname(normalizePath(script,winslash='/',mustWork=TRUE)),'..'))
suppressPackageStartupMessages({library(limma);library(edgeR);library(data.table);library(decoupleR)})
source('R/00_common.R')
setDTthreads(2);set.seed(42)
out <- 'outputs/tf_modules_v1';res <- 'resources/tf_modules_v1'
MIN <- 5;WSUM_TIMES <- 1000
TEST <- Sys.getenv('TFMOD_TEST')=='1'   # quick end-to-end check; outputs go to a test folder
if(TEST) {WSUM_TIMES <- 5;out <- file.path(out,'_test')}
for(d in c('tables','models','provenance','figures')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
put <- function(z,n) fwrite(as.data.frame(z),file.path(out,'tables',paste0(n,'.csv')),na='NA')

# ---- priors: the HT1080 tf_modules_v1 files, copied and verified against its manifest --------------
HT <- '../../HT1080/SEPT_11_2026_REANALYZE/resources/tf_modules_v1'
pf <- c('collectri_current.csv','collectri_frozen_v1.csv','dorothea_ABC.csv','tfactprofiler_all.csv','tfactprofiler_fibroblast.csv')
dir.create(res,showWarnings=FALSE,recursive=TRUE)
ref <- fread(file.path(HT,'prior_hashes.csv'))
for(f in pf) {
 if(!file.exists(file.path(res,f))) stopifnot(file.copy(file.path(HT,f),file.path(res,f)))
 stopifnot(hash_file(file.path(res,f))==ref$sha256[ref$file==f])
}
fwrite(ref[ref$file %in% pf],file.path(res,'prior_hashes.csv'))

inputs <- c('outputs/models/01_ingest.rds','outputs/models/02_expression.rds','outputs/tables/collectri_ULM_effect_signatures.csv',file.path(res,pf))
hashes <- vapply(inputs,hash_file,character(1))
fwrite(data.frame(path=inputs,sha256=hashes),file.path(out,'provenance/stage1_input_hashes.csv'))
x <- readRDS(inputs[1]);e <- readRDS(inputs[2])
readnet <- function(p) as.data.frame(fread(p))
nets <- list(collectri=readnet(inputs[4]),dorothea_ABC=readnet(inputs[6]),
 tfactprofiler_all=readnet(inputs[7]),tfactprofiler_fibroblast=readnet(inputs[8]))
frozen <- readnet(inputs[5])
if(TEST) for(p in c('tfactprofiler_all','tfactprofiler_fibroblast')) nets[[p]] <- nets[[p]][nets[[p]]$source %in% unique(nets[[p]]$source)[1:40],]
mat1 <- function(v,name) matrix(v,ncol=1,dimnames=list(names(v),name))
score <- function(fun,v,net,...) {
 r <- as.data.frame(fun(mat=mat1(v,'x'),network=net,.source='source',.target='target',.mor='mor',minsize=MIN,...))
 r[,c('statistic','source','score','p_value')]
}
ulm_mat <- function(M,net) as.data.frame(run_ulm(mat=M,network=net,.source='source',.target='target',.mor='mor',minsize=MIN))

# ---- inputs on the representative-symbol universe ----------------------------------------------------
mp <- e$mapping[e$mapping$representative,c('gene_id','gene_name','arm','chrom','start','mean_logCPM')]
stopifnot(all(mp$gene_id %in% rownames(e$de$fit)),!anyDuplicated(mp$gene_name))
by_symbol <- function(v,ids) {v <- v[match(mp$gene_id,ids)];names(v) <- mp$gene_name;v}
m <- x$core;effective <- setNames(e$norm$effective,colnames(e$logCPM))[m$sample_key]
stopifnot(identical(colnames(e$logCPM),m$sample_key))
counts <- x$counts[rownames(e$de$voom),m$sample_key]
fit_cond <- function(keep) {
 mm <- m[keep,];D <- model.matrix(~condition,data.frame(condition=factor(mm$treatment,levels=c('CTRL','LOY'))));rownames(D) <- mm$sample_key
 v <- voom(counts[,keep],design=D,lib.size=effective[keep],normalize.method='none',plot=FALSE)
 f <- eBayes(lmFit(v,D),trend=FALSE,robust=TRUE)
 by_symbol(f$t[,'conditionLOY'],rownames(f))
}
overall <- by_symbol(e$de$fit$t[,'conditionLOY'],rownames(e$de$fit))
full <- fit_cond(rep(TRUE,nrow(m)))
stopifnot(max(abs(full-overall))<1e-8)   # refit reproduces the saved moderated t

# CTRL versus parental (drift flag only): all 11 libraries, saved context exposures, primary gene universe.
mc <- x$meta;stopifnot(length(e$context$effective)==nrow(mc),identical(colnames(e$context$logCPM),mc$sample_key))
Dp <- model.matrix(~0+group,data.frame(group=factor(mc$treatment,levels=c('CTRL','LOY','PARENTAL'))));rownames(Dp) <- mc$sample_key
vp <- voom(x$counts[rownames(e$de$voom),mc$sample_key],design=Dp,lib.size=e$context$effective,normalize.method='none',plot=FALSE)
fp <- eBayes(contrasts.fit(lmFit(vp,Dp),makeContrasts(groupCTRL-groupPARENTAL,levels=Dp)),trend=FALSE,robust=TRUE)
rankings <- list(overall=overall,CTRL_vs_parental=by_symbol(fp$t[,1],rownames(fp)))
effects <- list(overall=by_symbol(e$effects$mean_difference,e$effects$gene_id),
 CTRL_vs_parental=by_symbol(e$parent$CTRL_minus_PARENTAL_context,e$parent$gene_id))
stopifnot(all(vapply(c(rankings,effects),function(v) all(is.finite(v)),TRUE)))

# ---- bridge: frozen prior + mean difference + ULM reproduces the saved A549 CollecTRI table ---------------
old <- fread(inputs[3]);eff_all <- setNames(e$effects$mean_difference[match(e$mapping$gene_id,e$effects$gene_id)],e$mapping$gene_name)
r <- score(run_ulm,eff_all,frozen);mm <- merge(r,old[,c('source','score','p_value')],by='source',suffixes=c('_new','_saved'))
bridge <- data.frame(contrast='LOY_vs_sgCTRL',n_new=nrow(r),n_saved=nrow(old),matched=nrow(mm),
 max_abs_score_diff=max(abs(mm$score_new-mm$score_saved)),max_abs_p_diff=max(abs(mm$p_value_new-mm$p_value_saved)))
put(bridge,'TF_bridge_to_saved');print(bridge)
stopifnot(bridge$n_new==bridge$n_saved,bridge$max_abs_score_diff<1e-10)

# ---- stage 1: scores ------------------------------------------------------------------------------------
methods <- list(ulm=run_ulm,mlm=run_mlm,wsum=run_wsum,viper=run_viper)
ckpt <- file.path(out,'models/stage1_scores_checkpoint.rds')
all_scores <- if(file.exists(ckpt)) readRDS(ckpt) else list()
for(k in names(rankings)) for(p in names(nets)) {
 # TFActProfiler priors enter only R4 (overall); skip them elsewhere for runtime.
 if(grepl('^tfactprofiler',p) && k!='overall') next
 use <- if(p %in% c('collectri','dorothea_ABC')) names(methods) else 'ulm'
 use <- use[!paste(k,p,use) %in% names(all_scores)]
 for(mname in use) {
  message(k,' / ',p,' / ',mname)
  args <- switch(mname,wsum=list(times=WSUM_TIMES,seed=42),viper=list(pleiotropy=FALSE,verbose=FALSE),list())
  r <- do.call(score,c(list(methods[[mname]],rankings[[k]],nets[[p]]),args))
  if(mname=='wsum') r <- r[r$statistic=='norm_wsum',]
  r$method <- mname;r$prior <- p;r$contrast <- k;r$input <- 't'
  all_scores[[paste(k,p,mname)]] <- r
 }
 saveRDS(all_scores,ckpt)
}
for(k in names(effects)) {
 r <- score(run_ulm,effects[[k]],nets$collectri);r$method <- 'ulm';r$prior <- 'collectri';r$contrast <- k;r$input <- 'logFC'
 all_scores[[paste(k,'logFC')]] <- r
}
S <- rbindlist(all_scores,fill=TRUE)
S[,FDR:=p.adjust(p_value,'BH'),by=.(contrast,prior,method,input)]
S[,signed_z:=sign(score)*qnorm(1-pmax(p_value,1e-300)/2)]
size <- rbindlist(lapply(names(nets),function(p) {
 n <- as.data.table(nets[[p]])[target %in% names(overall)];n[,.(targets_covered=.N),by=source][,prior:=p]
}))
S <- merge(S,size,by=c('prior','source'),all.x=TRUE)
put(S,'TF_activity_matrix')

# ---- fry flag (gene correlation-aware) on the saved unpaired voom -------------------------------------------
cm <- nets$collectri[nets$collectri$target %in% mp$gene_name,];cm$gene_id <- mp$gene_id[match(cm$target,mp$gene_name)]
idx <- split(cm[,c('gene_id','mor')],cm$source);idx <- idx[vapply(idx,nrow,1L)>=MIN]
idx <- lapply(idx,function(d) data.frame(gene_id=d$gene_id,weight=d$mor))
fr <- fry(e$de$voom,index=idx,design=e$de$design,contrast=which(colnames(e$de$design)=='conditionLOY'),sort='none')
fr <- data.frame(source=rownames(fr),fr,row.names=NULL);put(fr,'TF_fry_collectri_overall')

# ---- R2 units: descriptive mean difference with each library omitted (stable-omission rule) ---------------
omit_cols <- grep('^omit_',names(e$effects),value=TRUE);stopifnot(length(omit_cols)==nrow(m))
U <- sapply(omit_cols,function(k) by_symbol(e$effects[[k]],e$effects$gene_id))
r2 <- ulm_mat(U,nets$collectri);names(r2)[names(r2)=='condition'] <- 'unit'
put(r2,'TF_R2_unit_scores')

# ---- R5: leave-one-sample-out refits of the ~ condition model --------------------------------------------
loso <- rbindlist(lapply(seq_len(nrow(m)),function(i) {
 r <- score(run_ulm,fit_cond(seq_len(nrow(m))!=i),nets$collectri);r$omitted_block <- m$sample_key[i];r
}))
put(loso,'TF_leave_one_sample_out_ulm')

# ---- stage 2: rule ---------------------------------------------------------------------------------------
g <- function(k,p='collectri',mt='ulm') S[contrast==k & prior==p & method==mt & input=='t',.(source,score,p_value,FDR)]
R <- g('overall');setnames(R,c('score','p_value','FDR'),c('ulm_score','ulm_p','ulm_FDR'))
z <- g('CTRL_vs_parental');setnames(z,c('score','p_value','FDR'),paste0('CTRL_vs_parental',c('_score','_p','_FDR')));R <- merge(R,z,by='source',all.x=TRUE)
R$R1 <- R$ulm_FDR<.05
R2t <- as.data.table(r2)[,.(R2_units=.N,R2_units_same_sign=sum(sign(score)==sign(R$ulm_score[match(source,R$source)])),R2_min_abs=min(abs(score))),by=source]
R <- merge(R,R2t,by='source',all.x=TRUE);R$R2 <- R$R2_units_same_sign==nrow(m)
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
L <- loso[,.(loo_same_sign=sum(sign(score)==sign(R$ulm_score[match(source,R$source)])),loo_min_abs=min(abs(score))),by=source]
R <- merge(R,L,by='source',all.x=TRUE);R$R5 <- R$loo_same_sign==nrow(m)

# R6: remove self/complex-subunit edges and the three largest contributors; rescore that TF alone.
complex_sub <- list(AP1=c('JUN','JUNB','JUND','FOS','FOSB','FOSL1','FOSL2'),NFKB=c('NFKB1','NFKB2','RELA','RELB','REL'))
tt <- overall;net <- nets$collectri[nets$collectri$target %in% names(tt),]
contrib <- rbindlist(lapply(R$source[R$R1 %in% TRUE],function(s) {
 n <- net[net$source==s,];sg <- sign(R$ulm_score[R$source==s])
 data.frame(source=s,target=n$target,mor=n$mor,t=tt[n$target],contribution=n$mor*tt[n$target]*sg,
  self_edge=n$target==s | n$target %in% complex_sub[[s]])
}))
if(nrow(contrib)) {
 contrib <- contrib[order(contrib$source,-contrib$contribution),]
 omit <- rbindlist(lapply(unique(contrib$source),function(s) {
  cc <- contrib[contrib$source==s,];drop <- c(cc$target[cc$self_edge],head(cc$target[!cc$self_edge],3))
  n <- net[net$source==s & !net$target %in% drop,]
  r <- if(nrow(n)>=MIN) score(run_ulm,tt,n) else data.frame(statistic='ulm',source=s,score=NA,p_value=NA)
  data.frame(source=s,omitted=paste(drop,collapse=';'),n_omitted=length(drop),remaining=nrow(n),omit_score=r$score,omit_p=r$p_value)
 }))
 R <- merge(R,omit,by='source',all.x=TRUE)
} else R[,`:=`(omitted=NA_character_,n_omitted=NA_integer_,remaining=NA_integer_,omit_score=NA_real_,omit_p=NA_real_)]
R$R6 <- !is.na(R$omit_score) & sign(R$omit_score)==sign(R$ulm_score) & R$omit_p<.05
R$fry_FDR <- fr$FDR[match(R$source,fr$source)];R$fry_direction <- fr$Direction[match(R$source,fr$source)]
R$targets_covered <- size[prior=='collectri'][match(R$source,source),targets_covered]
R$self_edge_present <- R$source %in% contrib$source[contrib$self_edge]
rr <- c('R1','R2','R3','R4','R5','R6')
for(r in rr) R[[r]][is.na(R[[r]])] <- FALSE
R$n_R3_R6 <- R$R3+R$R4+R$R5+R$R6
R$tier <- ifelse(R$R1 & R$R2 & R$R3 & R$R4 & R$R5 & R$R6,'robust',ifelse(R$R1 & R$R2 & R$n_R3_R6>=3,'supported',
 ifelse(R$R1 & R$R2,'consistent_only',ifelse(R$R1,'overall_only','not_significant'))))
R$parental_drift_flag <- R$tier %in% c('robust','supported') & sign(R$CTRL_vs_parental_score)==sign(R$ulm_score) & R$CTRL_vs_parental_p<.05
R$rule_scope <- 'A549 stability tiers: R2 = 8/8 sample-omission effects, R5 = 8/8 leave-one-sample-out refits; not replication'
R <- R[order(match(R$tier,c('robust','supported','consistent_only','overall_only','not_significant')),-abs(R$consensus_z)),]
put(R,'TF_robustness');put(contrib,'TF_target_contributions_R1')
print(table(R$tier));print(R[R$tier %in% c('robust','supported'),c('source','ulm_score','ulm_FDR','methods_agree','R4','loo_same_sign','R6','fry_FDR','tier')])
saveRDS(list(scores=S,robustness=R,rankings=rankings,effects=effects,R2_units=U,R2_rule=nrow(m),loso=loso,fry=fr,contrib=contrib,bridge=bridge),
 file.path(out,'models/tf_activity_v1.rds'))
stopifnot(identical(hashes,vapply(inputs,hash_file,character(1))))
writeLines(capture.output(sessionInfo()),file.path(out,'provenance/stage1_sessionInfo.txt'))
cat('Stage 1-2 completed\n')
