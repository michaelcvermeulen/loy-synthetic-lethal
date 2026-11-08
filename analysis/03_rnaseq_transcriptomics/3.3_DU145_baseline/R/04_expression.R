group_difference <- function(l,m) rowMeans(l[,m$treatment=='LOY',drop=FALSE])-rowMeans(l[,m$treatment=='CTRL',drop=FALSE])

# Paired (~ passage + condition) is primary; unpaired (~ condition) is a sensitivity model.
fit_voom <- function(counts,annotation,effective,meta,paired=TRUE) {
  dd <- data.frame(condition=factor(meta$treatment,levels=c('CTRL','LOY')),
                   passage=factor(meta$passage_block))
  design <- if (paired) model.matrix(~passage+condition,dd) else model.matrix(~condition,dd)
  rownames(design) <- meta$sample_key
  assert(qr(design)$rank==ncol(design),'Non-estimable model')
  v <- limma::voom(counts,design,lib.size=effective,normalize.method='none',plot=FALSE,save.plot=TRUE)
  assert(all(is.finite(v$weights)) && all(v$weights>0),'Invalid voom weights')
  fit <- limma::eBayes(limma::lmFit(v,design),robust=TRUE,trend=FALSE)
  z <- limma::topTable(fit,coef='conditionLOY',number=Inf,sort.by='none',adjust.method='BH',confint=TRUE)
  z <- z[match(annotation$gene_id,rownames(z)),,drop=FALSE]
  names(z)[names(z)=='P.Value'] <- 'PValue'; names(z)[names(z)=='adj.P.Val'] <- 'FDR'
  assert(all(is.finite(z$PValue)) && max(abs(z$FDR-p.adjust(z$PValue,'BH')))<1e-12,'BH mismatch')
  z <- cbind(annotation,z); rownames(z) <- NULL
  z$residual_df <- nrow(design)-ncol(design)
  z$prior_df <- if (length(fit$df.prior)==1) fit$df.prior else median(fit$df.prior)
  z$statistical_scope <- if (paired) 'primary_paired_passage_model_one_population_per_arm_not_editing_replication' else
    'sensitivity_unpaired_model_ignores_passage_matching'
  list(table=z,voom=v,fit=fit,design=design)
}

# Thesis rule shared with HT1080: FDR, effect size, and direction agreement across matched units.
thesis_de_call <- function(z,pair_diff,cfg) {
  s <- sign(z$logFC)
  agree <- rowSums(sign(pair_diff)==s)
  z$pairs_agreeing_with_logFC <- agree
  z$thesis_DE <- z$FDR<cfg$de_fdr & abs(z$logFC)>=cfg$de_abs_log2fc & agree>=cfg$de_pairs_required
  z$thesis_direction <- ifelse(z$thesis_DE,ifelse(z$logFC>0,'higher_in_LOY','lower_in_LOY'),NA_character_)
  z
}

du145_expression <- function(x,cfg) {
  a <- x$annotation; m <- x$core; counts <- x$counts[,m$sample_key,drop=FALSE]
  n <- normalize_counts(counts,a$main_background,m$treatment,cfg)
  l <- n$logCPM
  e <- a; e$mean_logCPM <- rowMeans(l); e$mean_CPM <- rowMeans(n$CPM)
  e$CTRL_mean_CPM <- rowMeans(n$CPM[,m$treatment=='CTRL',drop=FALSE])
  e$LOY_mean_CPM <- rowMeans(n$CPM[,m$treatment=='LOY',drop=FALSE])
  e$LOY_to_CTRL_CPM <- ifelse(e$CTRL_mean_CPM>0,e$LOY_mean_CPM/e$CTRL_mean_CPM,NA_real_)
  e$mean_difference <- group_difference(l,m)
  pd <- matched_differences(l,x$pairs)
  colnames(pd) <- paste0('pair_',colnames(pd))
  e$pairs_same_direction_as_mean <- rowSums(sign(pd)==sign(e$mean_difference))
  e$expressed <- n$expression_keep; e$primary_eligible <- n$keep
  e$exclusion_reason <- ifelse(!a$main_background,'separate_Y_PAR_MT_or_other',ifelse(!n$expression_keep,'below_filter','retained_primary'))
  e$CTRL_sample_SD <- apply(l[,m$treatment=='CTRL',drop=FALSE],1,sd)
  e$LOY_sample_SD <- apply(l[,m$treatment=='LOY',drop=FALSE],1,sd)
  # Leave-one-pair-out: drop both libraries of a passage block, refit filtering and TMM.
  lopo <- sapply(x$pairs$block,function(b) {
    keep_s <- m$passage_block!=b
    nn <- normalize_counts(counts[,keep_s,drop=FALSE],a$main_background,m$treatment[keep_s],cfg)
    group_difference(nn$logCPM,m[keep_s,])
  })
  colnames(lopo) <- paste0('omit_pair_',x$pairs$block)
  e$LOPO_min <- apply(lopo,1,min); e$LOPO_max <- apply(lopo,1,max)
  e$LOPO_max_abs_change <- apply(abs(lopo-e$mean_difference),1,max)
  e$LOPO_all_same_direction <- rowSums(sign(lopo)==sign(e$mean_difference))==ncol(lopo)
  lower <- edgeR::cpm(counts,lib.size=n$effective,log=TRUE,prior.count=0.25)
  e$prior025_mean_difference <- group_difference(lower,m)
  rle <- edgeR::calcNormFactors(edgeR::DGEList(counts[n$keep,,drop=FALSE]),method='RLE')
  lr <- edgeR::cpm(counts,lib.size=rle$samples$lib.size*rle$samples$norm.factors,log=TRUE,prior.count=cfg$prior_count)
  e$RLE_mean_difference <- group_difference(lr,m)
  e <- cbind(e,as.data.frame(pd),as.data.frame(lopo))
  tab(e,'DU145_LOY_vs_sgCTRL_all_gene_effects'); tab(e[n$keep,],'DU145_LOY_vs_sgCTRL_primary_effects')
  tab(data.frame(sample_key=m$sample_key,n$dge$samples,effective_library_reads=n$effective,all_feature_reads=colSums(counts)),'normalization_core')
  tab(matrix_long(l,'log2_CPM'),'core_logCPM_all_genes')
  tab(matrix_long(n$CPM,'CPM'),'core_CPM_all_genes')
  mapping <- choose_symbols(a,rowMeans(l),n$keep)
  tab(mapping,'symbol_mapping_all_loci'); mapping <- mapping[mapping$representative,]
  tab(mapping,'symbol_mapping_representatives')
  # Global structure.
  pc <- prcomp(t(l[n$keep,]),center=TRUE,scale.=FALSE)
  pca <- cbind(m,as.data.frame(pc$x[,1:4]))
  tab(pca,'PCA_core')
  pvar <- data.frame(PC=seq_along(pc$sdev),core_variance=pc$sdev^2/sum(pc$sdev^2))
  tab(pvar,'PCA_variance_core')
  corr <- cor(l[n$keep,]); tab(matrix_long(corr,'correlation','sample_1'),'sample_correlations_core')
  qc <- x$meta; qc$assigned_reads <- colSums(x$counts); qc$detected_genes <- colSums(x$counts>0)
  qc$assigned_fraction <- vapply(qc$sample_key,function(s) {z<-x$assignment[x$assignment$sample_key==s,]; z$reads[z$status=='Assigned']/sum(z$reads)},numeric(1))
  qc$mitochondrial_fraction <- colSums(x$counts[a$chrom=='chrM',,drop=FALSE])/qc$assigned_reads
  qc$MSY_fraction <- colSums(x$counts[a$chrom=='chrY' & !a$PAR,,drop=FALSE])/qc$assigned_reads
  qc$ribosomal_protein_fraction <- colSums(x$counts[a$ribosomal_protein,,drop=FALSE])/qc$assigned_reads
  qc$top10_fraction <- apply(x$counts,2,function(v)sum(sort(v,decreasing=TRUE)[1:10]))/qc$assigned_reads
  tab(qc,'library_QC')
  qcore <- qc[match(m$sample_key,qc$sample_key),]
  pc_assoc <- do.call(rbind,lapply(paste0('PC',1:3),function(k) data.frame(axis=k,
    variance=pvar$core_variance[as.integer(sub('PC','',k))],
    treatment_R2=summary(lm(pca[[k]]~m$treatment))$r.squared,
    passage_R2=summary(lm(pca[[k]]~m$passage_block))$r.squared,
    mitochondrial_Spearman=cor(pca[[k]],qcore$mitochondrial_fraction,method='spearman'))))
  tab(pc_assoc,'PCA_design_associations')
  # Primary paired model and the thesis DE call.
  anno <- e[n$keep,c('gene_id','gene_name','gene_type','chrom','arm','start','end','mean_CPM','mean_logCPM','mean_difference','LOPO_all_same_direction')]
  de <- fit_voom(counts[n$keep,],anno,n$effective,m,paired=TRUE)
  de$table <- thesis_de_call(de$table,pd[n$keep,,drop=FALSE],cfg)
  tab(de$table,'DE_paired_voom_primary')
  tab(de$table[de$table$thesis_DE,],'DE_thesis_set')
  tab(data.frame(sample_key=m$sample_key,de$design,check.names=FALSE),'DE_design_paired')
  tab(data.frame(mean_log2_count=de$voom$voom.xy$x,sqrt_residual_SD=de$voom$voom.xy$y),'voom_mean_variance')
  tab(data.frame(mean_log2_count=de$voom$voom.line$x,fitted_sqrt_residual_SD=de$voom$voom.line$y),'voom_trend')
  unp <- fit_voom(counts[n$keep,],anno,n$effective,m,paired=FALSE)
  tab(unp$table,'DE_unpaired_voom_sensitivity')
  t <- de$table
  de_summary <- data.frame(
    rule=c('FDR<0.05','FDR<0.05 & |log2FC|>=1','thesis: FDR<0.05 & |log2FC|>=1 & 3/3 pairs','unpaired FDR<0.05'),
    genes=c(sum(t$FDR<cfg$de_fdr),sum(t$FDR<cfg$de_fdr & abs(t$logFC)>=cfg$de_abs_log2fc),sum(t$thesis_DE),sum(unp$table$FDR<cfg$de_fdr)),
    higher=c(sum(t$FDR<cfg$de_fdr & t$logFC>0),sum(t$FDR<cfg$de_fdr & abs(t$logFC)>=cfg$de_abs_log2fc & t$logFC>0),sum(t$thesis_DE & t$logFC>0),sum(unp$table$FDR<cfg$de_fdr & unp$table$logFC>0)),
    lower=c(sum(t$FDR<cfg$de_fdr & t$logFC<0),sum(t$FDR<cfg$de_fdr & abs(t$logFC)>=cfg$de_abs_log2fc & t$logFC<0),sum(t$thesis_DE & t$logFC<0),sum(unp$table$FDR<cfg$de_fdr & unp$table$logFC<0)),
    residual_df=c(rep(t$residual_df[1],3),unp$table$residual_df[1]))
  tab(de_summary,'DE_summary')
  # Cross-line comparison: exact shared IDs, independently normalised thesis analyses.
  htpath <- file.path(REPO,'analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/tables/HT1080_LOY_vs_sgCTRL_primary_effects.csv')
  apath <- file.path(REPO,'analysis/A549/SEPT_18_2026_REANALYZE/outputs/tables/A549_LOY_vs_sgCTRL_primary_effects.csv')
  ht <- as.data.frame(data.table::fread(htpath,select=c('gene_id','mean_difference')))
  a5 <- as.data.frame(data.table::fread(apath,select=c('gene_id','mean_difference')))
  names(ht)[2] <- 'mean_difference_HT1080'; names(a5)[2] <- 'mean_difference_A549'
  comp <- e[n$keep,c('gene_id','gene_name','chrom','arm','mean_difference')]
  names(comp)[names(comp)=='mean_difference'] <- 'mean_difference_DU145'
  comp <- merge(merge(comp,ht,by='gene_id',all.x=TRUE),a5,by='gene_id',all.x=TRUE)
  comp$DU145_thesis_DE <- comp$gene_id %in% t$gene_id[t$thesis_DE]
  tab(comp,'DU145_cross_line_effect_comparison')
  conc <- function(xa,ya,lab) {ok<-is.finite(xa)&is.finite(ya); data.frame(comparison=lab,shared_genes=sum(ok),Pearson=cor(xa[ok],ya[ok]),Spearman=cor(xa[ok],ya[ok],method='spearman'))}
  tab(rbind(conc(comp$mean_difference_DU145,comp$mean_difference_HT1080,'DU145_vs_HT1080'),
            conc(comp$mean_difference_DU145,comp$mean_difference_A549,'DU145_vs_A549')),'DU145_cross_line_concordance')
  tab(data.frame(analysis=c('paired_passage_model','temporal_repeatability','independent_derivation_replication','DNA_RNA_integration','ATAC','parental_context'),
    status=c('primary','not_estimable','not_claimed','run_gene_arm_window','not_available','not_available'),
    reason=c('Libraries are passage-matched P12/P13/P14 across arms','Only one sampling window (May 2026)','One edited population per arm; passages are repeated sampling, not independent editing events','Author confirmed RNA LOY population is the 2x population profiled by lpWGS LNCaP-LOY-2X; specimens not matched','No DU145 ATAC-seq','No parental DU145 RNA libraries in this dataset')),'analysis_status')
  inputs <- c(file.path(REPO,cfg$raw_counts),file.path(REPO,cfg$gtf),cfg$sample_manifest,htpath,apath)
  tab(data.frame(path=inputs,sha256=vapply(inputs,hash_file,character(1))),'input_hashes')
  list(effects=e,norm=n,logCPM=l,pair_diff=pd,mapping=mapping,pca=pca,pca_variance=pvar$core_variance,pc_assoc=pc_assoc,
       correlations=corr,qc=qc,de=de,unpaired=unp,de_summary=de_summary,comparison=comp)
}
