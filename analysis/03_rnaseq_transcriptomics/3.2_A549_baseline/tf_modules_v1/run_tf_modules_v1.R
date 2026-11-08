# Stages 3-5: overlap/redundancy, copy-number flags, shared-target modules, stability, nulls, co-activity,
# annotation and confounder checks, adapted from HT1080 tf_modules_v1. Rules: provenance/DECISIONS.md (27 Sep 2026).
options(warn=1,stringsAsFactors=FALSE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(file.path(dirname(normalizePath(script,winslash='/',mustWork=TRUE)),'..'))
suppressPackageStartupMessages({library(limma);library(data.table);library(decoupleR);library(cluster);library(igraph)})
source('R/00_common.R')
setDTthreads(2);set.seed(42)
out <- 'outputs/tf_modules_v1';res <- 'resources/tf_modules_v1'
put <- function(z,n) fwrite(as.data.frame(z),file.path(out,'tables',paste0(n,'.csv')),na='NA')
MIN <- 5;N_BOOT <- 1000;N_DROP <- 200;N_REWIRE <- 200;N_PERM <- 10000;N_RAND <- 1000
TEST <- Sys.getenv('TFMOD_TEST')=='1'
if(TEST) {N_BOOT <- 10;N_DROP <- 5;N_REWIRE <- 3;N_PERM <- 50;N_RAND <- 10;out <- file.path(out,'_test')}
HT <- '../../HT1080/SEPT_11_2026_REANALYZE'

# ==== line-specific block (A549) ==========================================================================
LINE <- 'A549'
line_inputs <- c(cn='A:/Dropbox/LOY_paralogs/projects/WGS/A549_QDNAseq_cn_df.rds',cyto='resources/hg38_cytoBand.txt')
# Arm-level rule as in HT1080: |median log2_seg LOY - sgCTRL| >= 0.15 over used bins.
cn_changed <- function(mp) {
 cyto <- fread(line_inputs[['cyto']],header=FALSE,col.names=c('chrom','start','end','band','stain'))
 cen <- cyto[stain=='acen',.(cen=min(start)),by=chrom]
 cn <- as.data.table(readRDS(line_inputs[['cn']]))[use==TRUE];cn[,chrom:=as.character(chromosome)]
 cn <- merge(cn,cen,by='chrom');cn[,arm:=paste0(chrom,ifelse(start<cen,'p','q'))]
 armcn <- dcast(cn[,.(seg=median(log2_seg)),by=.(arm,label)],arm~label,value.var='seg');armcn[,diff:=LOY-sgCTRL]
 armcn[,changed:=abs(diff)>=.15]
 list(table=armcn,genes=mp$gene_name[mp$arm %in% armcn$arm[armcn$changed]])
}
resid_design <- function(m) model.matrix(~treatment,m)
consistent_set <- function(e,mp) mp$gene_name[mp$gene_id %in% e$effects$gene_id[e$effects$descriptive_candidate & e$effects$primary_eligible]]
extra_refs <- function(e,mp) list()
ieg_table <- function(e,ids) {
 z <- e$effects[match(ids,e$effects$gene_id),c('gene_id','gene_name','mean_difference','LOO_min','LOO_max','LOO_all_same_direction')];z
}
# ==== end line-specific block =============================================================================

inputs <- c(file.path(out,'models/tf_activity_v1.rds'),'outputs/models/01_ingest.rds','outputs/models/02_expression.rds',
 file.path(res,'collectri_current.csv'),'resources/msigdb_all_collections.rds',line_inputs,
 file.path(HT,'outputs/tf_modules_v1/tables/TF_module_driver_targets.csv'))
hashes <- vapply(inputs,hash_file,character(1))
fwrite(data.frame(path=inputs,sha256=hashes),file.path(out,'provenance/stage3_input_hashes.csv'))
A <- readRDS(inputs[1]);x <- readRDS(inputs[2]);e <- readRDS(inputs[3])
R <- as.data.frame(A$robustness);tt <- A$rankings$overall
net <- as.data.frame(fread(inputs[4]));net <- net[net$target %in% names(tt),c('source','target','mor')]
mp <- e$mapping[e$mapping$representative,];stopifnot(identical(mp$gene_name,names(tt)))
consistent <- consistent_set(e,mp)
ht_drivers <- fread(inputs[length(inputs)])[module==1]$target
complex_sub <- list(AP1=c('JUN','JUNB','JUND','FOS','FOSB','FOSL1','FOSL2'),NFKB=c('NFKB1','NFKB2','RELA','RELB','REL'))
ulm1 <- function(v,n) {
 if(nrow(n)<MIN) return(data.frame(source=unique(n$source),score=NA,p_value=NA))
 r <- as.data.frame(run_ulm(mat=matrix(v,ncol=1,dimnames=list(names(v),'x')),network=n,.source='source',.target='target',.mor='mor',minsize=MIN))
 r[,c('source','score','p_value')]
}
ulm_mat <- function(M,n) as.data.frame(run_ulm(mat=M,network=n,.source='source',.target='target',.mor='mor',minsize=MIN))
UM <- cbind(overall=tt,A$R2_units[names(tt),,drop=FALSE])
finish <- function(obj,msg) {
 saveRDS(obj,file.path(out,'models/tf_modules_v1.rds'))
 stopifnot(identical(hashes,vapply(inputs,hash_file,character(1))))
 writeLines(capture.output(sessionInfo()),file.path(out,'provenance/stage3_sessionInfo.txt'))
 cat(msg,'\n');quit(save='no',status=0)
}

# ---- nodes -------------------------------------------------------------------------------------------
robust <- R$source[R$tier=='robust'];supported <- R$source[R$tier=='supported']
nodes <- if(length(robust)>=8) robust else c(robust,supported)
node_rule <- if(length(robust)>=8) 'robust only' else 'robust plus supported (fewer than 8 robust)'
message(LINE,': ',length(robust),' robust, ',length(supported),' supported; nodes: ',node_rule,' (',length(nodes),')')
sgn <- setNames(sign(R$ulm_score[match(nodes,R$source)]),nodes)
shortfall <- function(reason) put(data.frame(line=LINE,robust=length(robust),supported=length(supported),nodes=length(nodes),
 node_rule=node_rule,reason=reason),'TF_module_shortfall')

# ---- contributions and leading contributors ------------------------------------------------------------
contributions <- function(tvec,network,tfs,signs) {
 rbindlist(lapply(tfs,function(s) {
  n <- network[network$source==s & network$target %in% names(tvec),]
  data.table(source=s,target=n$target,mor=n$mor,t=tvec[n$target],contribution=n$mor*tvec[n$target]*signs[[s]])
 }))
}
leading <- function(C) {
 C <- C[contribution>0][order(source,-contribution)]
 C[,cum:=cumsum(contribution)/sum(contribution),by=source]
 C[,keep:=c(TRUE,head(cum,-1)<.8),by=source]
 C[keep==TRUE]
}

# ---- per-regulator copy-number flags (written for any number of nodes) ------------------------------------
CN <- cn_changed(mp);put(CN$table,'TF_copy_number')
changed_genes <- CN$genes
if(!length(nodes)) {shortfall('no robust or supported regulators');finish(list(nodes=nodes,node_rule=node_rule),'Stages 3-5 skipped: no nodes')}
Cn <- contributions(tt,net,nodes,sgn);Ln <- leading(Cn)
conc <- Cn[contribution>0][order(source,-contribution),.(positive_total=sum(contribution),top3_share=sum(head(contribution,3))/sum(contribution),
 leading_n=NA_integer_,top_targets=paste(head(target,5),collapse=';')),by=source]
conc$leading_n <- Ln[,.N,by=source][match(conc$source,source),N];conc$concentration_flag <- conc$top3_share>.5
cnflag <- rbindlist(lapply(nodes,function(s) {
 l <- Ln[source==s];arms <- mp$arm[match(l$target,mp$gene_name)];top <- sort(table(arms),decreasing=TRUE)
 r <- ulm1(tt,net[net$source==s & !net$target %in% changed_genes,])
 data.frame(source=s,top_arm=names(top)[1],top_arm_share=as.numeric(top[1])/nrow(l),cn_removed_score=r$score[1],cn_removed_p=r$p_value[1])
}))
cnflag$cn_flag <- cnflag$top_arm_share>=.5 | is.na(cnflag$cn_removed_score) | sign(cnflag$cn_removed_score)!=sgn[cnflag$source] | cnflag$cn_removed_p>=.05
if(length(nodes)<3) {
 put(merge(as.data.frame(conc),as.data.frame(cnflag),by='source'),'TF_overlap_redundancy')
 shortfall('fewer than 3 nodes; clustering not run');finish(list(nodes=nodes,node_rule=node_rule,conc=conc,cnflag=cnflag),'Stages 3-5 stopped: fewer than 3 nodes')
}

# ---- similarity and clustering ---------------------------------------------------------------------------
simmat <- function(L,tfs,type=c('weighted','jaccard'),restrict=NULL) {
 type <- match.arg(type);if(!is.null(restrict)) L <- L[target %in% restrict]
 S <- matrix(0,length(tfs),length(tfs),dimnames=list(tfs,tfs));diag(S) <- 1
 w <- lapply(tfs,function(s) setNames(L$contribution[L$source==s],L$target[L$source==s]));names(w) <- tfs
 for(i in seq_along(tfs)) for(j in seq_along(tfs)) if(i<j) {
  a <- w[[i]];bb <- w[[j]];sh <- intersect(names(a),names(bb))
  v <- if(!length(a)||!length(bb)) 0 else if(type=='weighted') sum(pmin(a[sh],bb[sh]))/min(sum(a),sum(bb)) else length(sh)/length(union(names(a),names(bb)))
  S[i,j] <- S[j,i] <- v
 }
 S
}
cluster_S <- function(S,k=NULL) {
 d <- as.dist(1-S);h <- hclust(d,method='average');n <- nrow(S)
 ks <- 2:max(2,floor(n/2));sil <- sapply(ks,function(kk) {cl <- cutree(h,kk);if(length(unique(cl))<2) NA else mean(silhouette(cl,d)[,3])})
 if(is.null(k)) k <- ks[which.max(sil)]
 list(h=h,k=k,cl=cutree(h,k),sil=data.frame(k=ks,mean_silhouette=sil),best_sil=max(sil,na.rm=TRUE))
}
S <- simmat(Ln,nodes,'weighted');CL <- cluster_S(S)
clusters <- CL$cl;size <- table(clusters);clusters[clusters %in% as.integer(names(size)[size==1])] <- NA
relabel <- function(cl) {u <- unique(cl[!is.na(cl)]);ord <- u[order(-sapply(u,function(z) sum(cl==z,na.rm=TRUE)))];setNames(match(cl,ord),names(cl))}
module <- relabel(clusters)
put(data.frame(source=rownames(S),S,check.names=FALSE),'TF_module_similarity_weighted')
put(CL$sil,'TF_module_silhouette')
mods <- sort(unique(module[!is.na(module)]))

# ---- overlap and redundancy among nodes --------------------------------------------------------------------
full_sets <- split(net$target,net$source)
unique_scores <- rbindlist(lapply(nodes,function(s) {
 others <- unique(unlist(full_sets[setdiff(nodes,s)]));n <- net[net$source==s & !net$target %in% others,]
 r <- ulm1(tt,n)
 ov <- sapply(setdiff(nodes,s),function(o) length(intersect(full_sets[[s]],full_sets[[o]]))/length(full_sets[[s]]))
 data.frame(source=s,unique_targets=nrow(n),unique_score=r$score[1],unique_p=r$p_value[1],
  most_overlapping=names(which.max(ov)),overlap_fraction=max(ov))
}))
unique_scores$separable <- ifelse(unique_scores$unique_targets<MIN,NA,sign(unique_scores$unique_score)==sgn[unique_scores$source] & unique_scores$unique_p<.05)
cond <- rbindlist(lapply(nodes,function(s) {
 ov <- sort(sapply(setdiff(nodes,s),function(o) length(intersect(full_sets[[s]],full_sets[[o]]))),decreasing=TRUE)
 nb <- names(ov)[ov>0][seq_len(min(5,sum(ov>0)))];n <- net[net$source %in% c(s,nb),]
 r <- tryCatch(as.data.frame(run_mlm(mat=matrix(tt,ncol=1,dimnames=list(names(tt),'x')),network=n,.source='source',.target='target',.mor='mor',minsize=MIN)),error=function(z) NULL)
 if(is.null(r)||!s %in% r$source) return(data.frame(source=s,neighbours=paste(nb,collapse=';'),conditional_score=NA,conditional_p=NA))
 data.frame(source=s,neighbours=paste(nb,collapse=';'),conditional_score=r$score[r$source==s],conditional_p=r$p_value[r$source==s])
}))
red <- Reduce(function(a,bb) merge(a,bb,by='source'),list(as.data.frame(conc),as.data.frame(unique_scores),as.data.frame(cond),as.data.frame(cnflag)))
red$module <- module[red$source];put(red,'TF_overlap_redundancy')

# ---- sensitivities: Jaccard, line-consistent restriction, Leiden -------------------------------------------------
ari <- function(a,bb) {ok <- !is.na(a)&!is.na(bb);if(sum(ok)<2) NA else mclust::adjustedRandIndex(a[ok],bb[ok])}
sens <- list()
for(nm in c('jaccard','line_consistent_set')) {
 Ss <- if(nm=='jaccard') simmat(Ln,nodes,'jaccard') else simmat(Ln,nodes,'weighted',restrict=consistent)
 cs <- cluster_S(Ss);sens[[nm]] <- data.frame(sensitivity=nm,k=cs$k,ARI_vs_primary=ari(module,cs$cl),best_silhouette=cs$best_sil)
}
G <- graph_from_adjacency_matrix(S-diag(nrow(S)),mode='undirected',weighted=TRUE,diag=FALSE)
for(r in c(.25,.5,.75,1,1.5)) {
 set.seed(42);lc <- cluster_leiden(G,objective_function='modularity',resolution=r,n_iterations=-1)
 sens[[paste0('leiden_',r)]] <- data.frame(sensitivity=paste0('leiden_resolution_',r),k=length(unique(membership(lc))),
  ARI_vs_primary=ari(module,setNames(membership(lc),V(G)$name)[nodes]),best_silhouette=NA)
}
put(rbindlist(sens),'TF_module_sensitivity')

# ---- nulls: rewired network (R1 + R2 count, silhouette) and size-matched random TF sets --------------------------------
rew <- rbindlist(lapply(seq_len(N_REWIRE),function(i) {
 set.seed(1000+i);n <- net;n$target <- sample(n$target);n <- n[!duplicated(paste(n$source,n$target)),]
 sc <- ulm_mat(UM,n);o <- sc[sc$condition=='overall',];o$FDR <- p.adjust(o$p_value,'BH')
 u <- sc[sc$condition!='overall',];same <- tapply(sign(u$score)==sign(o$score[match(u$source,o$source)]),u$source,all)
 hit <- o$source[o$FDR<.05 & same[o$source] %in% TRUE];bs <- NA
 if(length(hit)>=4) {sg <- setNames(sign(o$score[match(hit,o$source)]),hit);Ls <- leading(contributions(tt,n,hit,sg))
  hit2 <- intersect(hit,unique(Ls$source));if(length(hit2)>=4) bs <- cluster_S(simmat(Ls,hit2))$best_sil}
 data.frame(iteration=i,R1_R2_regulators=length(hit),best_silhouette=bs)
}))
put(rew,'TF_module_null_rewired')
cov_tfs <- names(which(table(net$source)>=MIN));rs <- table(net$source)[cov_tfs]
bins <- cut(rank(rs),10,labels=FALSE);names(bins) <- cov_tfs
allsg <- setNames(sign(R$ulm_score[match(cov_tfs,R$source)]),cov_tfs);allsg[is.na(allsg)|allsg==0] <- 1
Lall <- leading(contributions(tt,net,cov_tfs,allsg))
meanS <- function(S) mean(S[upper.tri(S)])
bin_tfs <- split(names(bins),bins);node_need <- table(bins[nodes])
obs_overlap <- meanS(S)
rand <- sapply(seq_len(N_RAND),function(i) {
 pick <- unlist(lapply(names(node_need),function(k) sample(bin_tfs[[k]],node_need[[k]])))
 pick <- intersect(unique(pick),unique(Lall$source));meanS(simmat(Lall,pick))
})
nulls <- data.frame(test=c('rewired_R1_R2_count','rewired_best_silhouette','size_matched_mean_overlap'),
 observed=c(sum(R$R1 & R$R2),CL$best_sil,obs_overlap),
 null_mean=c(mean(rew$R1_R2_regulators),mean(rew$best_silhouette,na.rm=TRUE),mean(rand)),
 null_q95=c(quantile(rew$R1_R2_regulators,.95),quantile(rew$best_silhouette,.95,na.rm=TRUE),quantile(rand,.95)),
 empirical_p=c((1+sum(rew$R1_R2_regulators>=sum(R$R1 & R$R2)))/(1+N_REWIRE),
  if(any(!is.na(rew$best_silhouette))) (1+sum(rew$best_silhouette>=CL$best_sil,na.rm=TRUE))/(1+sum(!is.na(rew$best_silhouette))) else NA,
  (1+sum(rand>=obs_overlap))/(1+N_RAND)),   # silhouette null is NA when no rewired network yields >= 4 hits
 null_iterations=c(N_REWIRE,sum(!is.na(rew$best_silhouette)),N_RAND))
put(nulls,'TF_module_nulls')
if(!length(mods)) {shortfall('clustering produced only singleton clusters; no module reported')
 finish(list(nodes=nodes,node_rule=node_rule,S=S,cluster=CL,module=module,nulls=nulls,rewired=rew,random_overlap=rand,redundancy=red),
  'Stages 3-5 stopped: no multi-regulator module')}

# ---- stability: gene bootstrap and edge drop (fixed k) ----------------------------------------------------------
cocluster <- matrix(0,length(nodes),length(nodes),dimnames=list(nodes,nodes))
jac_best <- function(ref,cl) sapply(mods,function(m) {
 a <- names(ref)[ref %in% m];max(sapply(unique(cl),function(z) {bb <- names(cl)[cl==z];length(intersect(a,bb))/length(union(a,bb))}))
})
boot_j <- list();drop_j <- list()
for(i in seq_len(N_BOOT)) {
 g <- sample(names(tt),replace=TRUE);w <- table(g)
 Cb <- Cn[target %in% names(w)];Cb[,contribution:=contribution*as.numeric(w[target])]
 Sb <- simmat(leading(Cb),nodes);cb <- cutree(hclust(as.dist(1-Sb),'average'),CL$k)
 cocluster <- cocluster+outer(cb,cb,'==');boot_j[[i]] <- jac_best(module,cb)
}
for(i in seq_len(N_DROP)) {
 keep <- sample(c(TRUE,FALSE),nrow(Cn),replace=TRUE,prob=c(.8,.2))
 Sd <- simmat(leading(Cn[keep]),nodes);cd <- cutree(hclust(as.dist(1-Sd),'average'),CL$k);drop_j[[i]] <- jac_best(module,cd)
}
bj <- matrix(unlist(boot_j),ncol=length(mods),byrow=TRUE);dj <- matrix(unlist(drop_j),ncol=length(mods),byrow=TRUE)
stab <- data.frame(module=mods,bootstrap_mean_jaccard=colMeans(bj),edge_drop_mean_jaccard=colMeans(dj))
stab$retained <- stab$bootstrap_mean_jaccard>=.6 & stab$edge_drop_mean_jaccard>=.6
stab$stable <- stab$bootstrap_mean_jaccard>=.75 & stab$edge_drop_mean_jaccard>=.75
stab$members <- sapply(stab$module,function(m) paste(names(module)[module %in% m],collapse=';'))
put(stab,'TF_module_stability')
put(data.frame(source=nodes,cocluster/N_BOOT,check.names=FALSE),'TF_module_cocluster_bootstrap')

# ---- module target unions ------------------------------------------------------------------------------------------
union_net <- function(members,n=net) {
 u <- n[n$source %in% members,];u$w <- u$mor*sgn[u$source]
 u <- aggregate(w~target,u,mean);data.frame(source='module',target=u$target,mor=u$w)
}
drivers <- rbindlist(lapply(mods,function(m) {
 mem <- names(module)[module %in% m];l <- Ln[source %in% mem]
 d <- l[,.(n_tfs=uniqueN(source),summed_contribution=sum(contribution)),by=target][n_tfs>=ceiling(length(mem)/2)][order(-summed_contribution)]
 d$module <- m;d$t <- tt[d$target];d$line_consistent <- d$target %in% consistent;d
}))
put(drivers,'TF_module_driver_targets')
fit_voom <- e$de$voom;fit_design <- e$de$design;coef <- which(colnames(fit_design)=='conditionLOY')
expr <- setNames(mp$mean_logCPM,mp$gene_name)[names(tt)];ebin <- cut(rank(expr),20,labels=FALSE);names(ebin) <- names(tt)
bin_genes <- split(names(ebin),ebin)
modtests <- rbindlist(lapply(mods,function(m) {
 u <- union_net(names(module)[module %in% m]);obs <- mean(sign(u$mor)*tt[u$target])
 need <- table(ebin[u$target]);sgn_u <- sign(u$mor)[order(ebin[u$target])]
 perm <- replicate(N_PERM,{g <- unlist(lapply(names(need),function(k) sample(bin_genes[[k]],need[[k]])))
  mean(sgn_u*tt[g])})
 fr <- fry(fit_voom,index=list(module=data.frame(gene_id=mp$gene_id[match(u$target,mp$gene_name)],weight=sign(u$mor))),design=fit_design,contrast=coef)
 r <- ulm1(tt,u)
 data.frame(module=m,union_targets=nrow(u),mean_signed_t=obs,perm_p=(1+sum(perm>=obs))/(1+N_PERM),fry_direction=fr$Direction,fry_p=fr$PValue,
  union_ulm_score=r$score,union_ulm_p=r$p_value)
}))
put(modtests,'TF_module_target_union_tests')

# ---- co-activity (supporting only) -----------------------------------------------------------------------------------
lcpm <- e$logCPM[mp$gene_id,x$core$sample_key];rownames(lcpm) <- mp$gene_name
Dm <- resid_design(x$core);resid <- lcpm-t(Dm%*%solve(crossprod(Dm),t(Dm)%*%t(lcpm)))
pairs <- t(combn(nodes,2))
coact <- rbindlist(lapply(seq_len(nrow(pairs)),function(i) {
 a <- pairs[i,1];bb <- pairs[i,2];sh <- intersect(full_sets[[a]],full_sets[[bb]])
 na <- net[net$source==a & !net$target %in% sh,];nb <- net[net$source==bb & !net$target %in% sh,]
 if(nrow(na)<MIN||nrow(nb)<MIN) return(data.frame(tf1=a,tf2=bb,r=NA))
 sa <- as.data.frame(run_ulm(resid,na,minsize=MIN));sb <- as.data.frame(run_ulm(resid,nb,minsize=MIN))
 data.frame(tf1=a,tf2=bb,r=cor(sa$score[match(colnames(resid),sa$condition)],sb$score[match(colnames(resid),sb$condition)]))
}))
coact$same_module <- !is.na(module[coact$tf1]) & !is.na(module[coact$tf2]) & module[coact$tf1]==module[coact$tf2]
ca <- coact[!is.na(r)]
if(any(ca$same_module) && any(!ca$same_module)) {
 obs <- mean(ca$r[ca$same_module])-mean(ca$r[!ca$same_module])
 perm <- replicate(N_PERM,{pm <- setNames(sample(module),names(module));sm <- !is.na(pm[ca$tf1]) & !is.na(pm[ca$tf2]) & pm[ca$tf1]==pm[ca$tf2]
  if(!any(sm)||all(sm)) NA else mean(ca$r[sm])-mean(ca$r[!sm])})
 nulls <- rbind(nulls,data.frame(test='coactivity_within_minus_between',observed=obs,null_mean=mean(perm,na.rm=TRUE),
  null_q95=quantile(perm,.95,na.rm=TRUE),empirical_p=(1+sum(perm>=obs,na.rm=TRUE))/(1+sum(!is.na(perm))),null_iterations=sum(!is.na(perm))))
}
put(coact,'TF_module_coactivity');put(nulls,'TF_module_nulls')

# ---- TF -> TF edges among nodes ------------------------------------------------------------------------------------------
gene_of <- function(s) if(s %in% names(complex_sub)) complex_sub[[s]] else s
tftf <- rbindlist(lapply(nodes,function(s) {
 tg <- net[net$source==s,];rbindlist(lapply(setdiff(nodes,s),function(o) {h <- tg[tg$target %in% gene_of(o),]
  if(nrow(h)) data.frame(from=s,to=o,target_gene=h$target,mor=h$mor,target_t=tt[h$target]) else NULL}))
}),fill=TRUE)
put(tftf,'TF_module_tf_to_tf_edges')

# ---- annotation: line-consistent genes, line-specific sets and HT1080 module 1 drivers ---------------------------------------
refsets <- c(list(line_consistent_set=consistent,HT1080_module1_drivers=ht_drivers),extra_refs(e,mp))
bg <- unique(net$target[net$source %in% nodes])
ann <- rbindlist(lapply(mods,function(m) {
 d <- drivers[module==m]$target;rbindlist(lapply(names(refsets),function(s) {
  ref <- intersect(refsets[[s]],bg);k <- length(intersect(d,ref))
  data.frame(module=m,reference=s,drivers=length(d),reference_in_background=length(ref),overlap=k,
   overlap_genes=paste(intersect(d,ref),collapse=';'),p=phyper(k-1,length(ref),length(bg)-length(ref),length(d),lower.tail=FALSE))
 }))
}))
ann[,BH_within_module:=p.adjust(p,'BH'),by=module];put(ann,'TF_module_annotation')

# ---- confounding programmes ------------------------------------------------------------------------------------------------
H <- readRDS(inputs[5])$H
conf <- c(H[c('HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT',
 'HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_P53_PATHWAY','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_HYPOXIA')],
 list(IEG_vanDenBrink2017_Fig1=c('FOS','FOSB','JUN','JUNB','ATF3','EGR1','SOCS3','ZFP36','HSP90AA1','CEBPB','CEBPD','HSPB1','HSPA1A','HSPA1B')))
stopifnot(all(lengths(conf)>0))
confr <- rbindlist(lapply(mods,function(m) {
 mem <- names(module)[module %in% m];u <- union_net(mem);d <- drivers[module==m]$target
 rbindlist(lapply(names(conf),function(p) {r <- ulm1(tt,u[!u$target %in% conf[[p]],])
  data.frame(module=m,programme=p,driver_overlap=length(intersect(d,conf[[p]])),driver_fraction=length(intersect(d,conf[[p]]))/max(1,length(d)),
   union_overlap=sum(u$target %in% conf[[p]]),removed_score=r$score,removed_p=r$p_value)}))
}))
confr$collapses <- is.na(confr$removed_score) | confr$removed_score<=0 | confr$removed_p>=.05
put(confr,'TF_module_confounders')
ieg_ids <- mp$gene_id[match(conf$IEG_vanDenBrink2017_Fig1,mp$gene_name)];ieg_ids <- ieg_ids[!is.na(ieg_ids)]
message('IEG genes expressed: ',length(ieg_ids),' of ',length(conf$IEG_vanDenBrink2017_Fig1))
put(ieg_table(e,ieg_ids),'TF_IEG_effects')

# ---- module summary ----------------------------------------------------------------------------------------------------------
msum <- rbindlist(lapply(mods,function(m) {
 mem <- names(module)[module %in% m];st <- stab[stab$module==m,];ct <- confr[confr$module==m,]
 data.frame(line=LINE,module=m,members=paste(mem,collapse=';'),n=length(mem),direction=paste(unique(ifelse(sgn[mem]>0,'higher','lower')),collapse='/'),
  bootstrap_jaccard=st$bootstrap_mean_jaccard,edge_drop_jaccard=st$edge_drop_mean_jaccard,retained=st$retained,stable=st$stable,
  drivers=nrow(drivers[module==m]),top_drivers=paste(head(drivers[module==m]$target,8),collapse=';'),
  confounder_dependent=paste(ct$programme[ct$collapses],collapse=';'))
}))
put(msum,'TF_modules')
memb <- data.frame(source=nodes,module=module[nodes],tier=R$tier[match(nodes,R$source)],ulm_score=R$ulm_score[match(nodes,R$source)],
 consensus_z=R$consensus_z[match(nodes,R$source)],order=match(nodes,CL$h$labels[CL$h$order]))
put(memb,'TF_module_membership')
print(msum);print(nulls);print(rbindlist(sens))
finish(list(nodes=nodes,node_rule=node_rule,S=S,cluster=CL,module=module,stability=stab,drivers=drivers,nulls=nulls,
 rewired=rew,random_overlap=rand,modtests=modtests,coactivity=coact,annotation=ann,confounders=confr,redundancy=red,summary=msum),
 'Stages 3-5 completed')
