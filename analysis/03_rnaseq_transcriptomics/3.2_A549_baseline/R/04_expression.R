group_difference <- function(l,m) rowMeans(l[,m$treatment=='LOY',drop=FALSE])-rowMeans(l[,m$treatment=='CTRL',drop=FALSE])

fit_unpaired <- function(counts,annotation,effective,meta,extra=NULL) {
  dd <- data.frame(condition=factor(meta$treatment,levels=c('CTRL','LOY')))
  if (is.null(extra)) design <- model.matrix(~condition,dd) else {
    dd$mitochondrial_fraction <- extra
    design <- model.matrix(~condition+mitochondrial_fraction,dd)
  }
  rownames(design) <- meta$sample_key
  assert(qr(design)$rank==ncol(design),'Non-estimable unpaired model')
  v <- limma::voom(counts,design,lib.size=effective,normalize.method='none',plot=FALSE,save.plot=TRUE)
  fit <- limma::eBayes(limma::lmFit(v,design),robust=TRUE,trend=FALSE)
  z <- limma::topTable(fit,coef='conditionLOY',number=Inf,sort.by='none',adjust.method='BH',confint=TRUE)
  z <- z[match(annotation$gene_id,rownames(z)),,drop=FALSE]
  names(z)[names(z)=='P.Value'] <- 'PValue'; names(z)[names(z)=='adj.P.Val'] <- 'FDR'
  z <- cbind(annotation,z); rownames(z) <- NULL
  z$residual_df <- nrow(design)-ncol(design)
  z$statistical_scope <- 'exploratory_unpaired_model_independent_residual_assumption_derivation_uncertain'
  list(table=z,voom=v,fit=fit,design=design)
}

a549_expression <- function(x,cfg) {
  a <- x$annotation; m <- x$core; counts <- x$counts[,m$sample_key,drop=FALSE]
  n <- normalize_counts(counts,a$main_background,m$treatment,cfg)
  l <- n$logCPM
  e <- a; e$mean_logCPM <- rowMeans(l); e$mean_CPM <- rowMeans(n$CPM)
  e$CTRL_mean_CPM <- rowMeans(n$CPM[,m$treatment=='CTRL',drop=FALSE])
  e$LOY_mean_CPM <- rowMeans(n$CPM[,m$treatment=='LOY',drop=FALSE])
  e$LOY_to_CTRL_CPM <- ifelse(e$CTRL_mean_CPM>0,e$LOY_mean_CPM/e$CTRL_mean_CPM,NA_real_)
  e$mean_difference <- group_difference(l,m)
  e$expressed <- n$expression_keep; e$primary_eligible <- n$keep
  e$exclusion_reason <- ifelse(!a$main_background,'separate_Y_PAR_MT_or_other',ifelse(!n$expression_keep,'below_filter','retained_primary'))
  e$CTRL_sample_SD <- apply(l[,m$treatment=='CTRL',drop=FALSE],1,sd)
  e$LOY_sample_SD <- apply(l[,m$treatment=='LOY',drop=FALSE],1,sd)
  loo <- sapply(seq_len(nrow(m)),function(j) {
    nn <- normalize_counts(counts[,-j,drop=FALSE],a$main_background,m$treatment[-j],cfg)
    group_difference(nn$logCPM,m[-j,])
  })
  colnames(loo) <- paste0('omit_',m$sample_key)
  e$LOO_min <- apply(loo,1,min); e$LOO_max <- apply(loo,1,max)
  e$LOO_max_abs_change <- apply(abs(loo-e$mean_difference),1,max)
  e$LOO_all_same_direction <- rowSums(sign(loo)==sign(e$mean_difference))==ncol(loo)
  lower <- edgeR::cpm(counts,lib.size=n$effective,log=TRUE,prior.count=0.25)
  e$prior025_mean_difference <- group_difference(lower,m)
  rle <- edgeR::calcNormFactors(edgeR::DGEList(counts[n$keep,,drop=FALSE]),method='RLE')
  lr <- edgeR::cpm(counts,lib.size=rle$samples$lib.size*rle$samples$norm.factors,log=TRUE,prior.count=cfg$prior_count)
  e$RLE_mean_difference <- group_difference(lr,m)
  e$descriptive_candidate <- e$primary_eligible & abs(e$mean_difference)>=cfg$effect_threshold & e$LOO_all_same_direction
  e <- cbind(e,as.data.frame(loo))
  tab(e,'A549_LOY_vs_sgCTRL_all_gene_effects'); tab(e[n$keep,],'A549_LOY_vs_sgCTRL_primary_effects')
  tab(e[e$descriptive_candidate,],'A549_LOY_vs_sgCTRL_descriptive_candidates')
  tab(data.frame(sample_key=m$sample_key,n$dge$samples,effective_library_reads=n$effective,all_feature_reads=colSums(counts)),'normalization_core')
  tab(matrix_long(l,'log2_CPM'),'core_logCPM_all_genes')
  tab(matrix_long(n$CPM,'CPM'),'core_CPM_all_genes')
  mapping <- choose_symbols(a,rowMeans(l),n$keep)
  tab(mapping,'symbol_mapping_all_loci'); mapping <- mapping[mapping$representative,]
  tab(mapping,'symbol_mapping_representatives')
  context <- normalize_counts(x$counts,a$main_background,x$meta$treatment,cfg)
  pc <- prcomp(t(l[n$keep,]),center=TRUE,scale.=FALSE)
  pc_context <- prcomp(t(context$logCPM[n$keep,]),center=TRUE,scale.=FALSE)
  pca <- cbind(m,as.data.frame(pc$x[,1:4])); pctx <- cbind(x$meta,as.data.frame(pc_context$x[,1:4]))
  tab(pca,'PCA_core'); tab(pctx,'PCA_context_all11')
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
  tab(data.frame(axis=c('PC1','PC2'),mitochondrial_Spearman=c(cor(pca$PC1,qcore$mitochondrial_fraction,method='spearman'),cor(pca$PC2,qcore$mitochondrial_fraction,method='spearman'))),'PCA_QC_associations')
  parent <- a
  for (tr in c('CTRL','LOY')) parent[[paste0(tr,'_minus_PARENTAL_context')]] <- rowMeans(context$logCPM[,x$meta$treatment==tr,drop=FALSE])-rowMeans(context$logCPM[,x$meta$treatment=='PARENTAL',drop=FALSE])
  tab(parent,'parental_context_expression_differences')
  anno <- e[n$keep,c('gene_id','gene_name','gene_type','chrom','start','end','mean_CPM','mean_logCPM','mean_difference','LOO_all_same_direction')]
  de <- fit_unpaired(counts[n$keep,],anno,n$effective,m)
  tab(de$table,'DE_unpaired_voom_exploratory')
  tab(data.frame(sample_key=m$sample_key,de$design),'DE_design_unpaired')
  tab(data.frame(mean_log2_count=de$voom$voom.xy$x,sqrt_residual_SD=de$voom$voom.xy$y),'voom_mean_variance')
  tab(data.frame(mean_log2_count=de$voom$voom.line$x,fitted_sqrt_residual_SD=de$voom$voom.line$y),'voom_trend')
  sens <- fit_unpaired(counts[n$keep,],anno,n$effective,m,extra=qcore$mitochondrial_fraction)
  sens$table$statistical_scope <- 'exploratory_mitochondrial_covariate_sensitivity_not_primary_data_derived_covariate'
  tab(sens$table,'DE_mitochondrial_covariate_sensitivity')
  tab(data.frame(sample_key=m$sample_key,sens$design),'DE_design_mitochondrial_sensitivity')
  # Separate chromosome branches retain the SAME primary normalization exposures.
  tab(e[a$chrom=='chrY' | a$PAR,],'Y_PAR_all_gene_effects')
  branches <- c('chrYp_nonPAR','chrYq_nonPAR','chrX_PAR','chrY_PAR','chrM')
  region <- do.call(rbind,lapply(branches,function(reg)data.frame(region=reg,m,sum_CPM=colSums(n$CPM[a$region==reg,,drop=FALSE]))))
  tab(region,'Y_PAR_MT_region_abundance')
  panel <- c('RPS4Y1','ZFY','USP9Y','DDX3Y','UTY','NLGN4Y','KDM5D','EIF1AY','DDX3X','EIF1AX','RPS4X','ZFX','KDM5C','KDM6A','USP9X','NLGN4X')
  selected <- which(a$gene_name %in% panel & !a$PAR)
  abundance <- do.call(rbind,lapply(selected,function(j)data.frame(a[j,c('gene_id','gene_name','chrom')],m,CPM=n$CPM[j,],raw_count=counts[j,])))
  tab(abundance,'XY_paralog_sample_abundance')
  pairs <- data.frame(X=c('KDM6A','KDM5C','DDX3X','EIF1AX','RPS4X','ZFX','USP9X','NLGN4X'),Y=c('UTY','KDM5D','DDX3Y','EIF1AY','RPS4Y1','ZFY','USP9Y','NLGN4Y'))
  find <- function(sym,ch) {ix<-which(a$gene_name==sym & a$chrom==ch & !a$PAR); if(!length(ix))NA_integer_ else ix[which.max(e$mean_CPM[ix])]}
  for(j in seq_len(nrow(pairs))) {
    ix<-find(pairs$X[j],'chrX'); iy<-find(pairs$Y[j],'chrY')
    for(s in c('X','Y')) {k<-if(s=='X')ix else iy; pairs[j,paste0(s,'_gene_id')]<-a$gene_id[k]; pairs[j,paste0(s,'_mean_difference')]<-e$mean_difference[k]; pairs[j,paste0(s,'_CTRL_CPM')]<-e$CTRL_mean_CPM[k]; pairs[j,paste0(s,'_LOY_CPM')]<-e$LOY_mean_CPM[k]}
    denom <- e$CTRL_mean_CPM[ix]/a$feature_length_bp[ix]+e$CTRL_mean_CPM[iy]/a$feature_length_bp[iy]
    pairs$relative_total_length_proxy[j] <- if(is.finite(denom) && denom>0)(e$LOY_mean_CPM[ix]/a$feature_length_bp[ix]+e$LOY_mean_CPM[iy]/a$feature_length_bp[iy])/denom else NA_real_
  }
  tab(pairs,'XY_paralog_effects_and_length_proxy')
  # Cross-line concordance is descriptive and uses independently normalised analyses.
  htpath <- file.path(REPO,'analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/tables/HT1080_LOY_vs_sgCTRL_primary_effects.csv')
  ht <- as.data.frame(data.table::fread(htpath,select=c('gene_id','gene_name','mean_difference')))
  comp <- merge(e[n$keep,c('gene_id','gene_name','mean_difference')],ht,by='gene_id',suffixes=c('_A549','_HT1080'))
  tab(comp,'A549_HT1080_effect_comparison')
  tab(data.frame(shared_genes=nrow(comp),Pearson=cor(comp$mean_difference_A549,comp$mean_difference_HT1080),Spearman=cor(comp$mean_difference_A549,comp$mean_difference_HT1080,method='spearman')),'A549_HT1080_concordance')
  tab(data.frame(analysis=c('paired_models','temporal_repeatability','independent_derivation_replication','DNA_RNA_coordinate_integration','HT1080_specific_HDACS_addbacks'),status=c('not_run','not_estimable','unconfirmed','deferred','not_transferred'),reason=c('Author states pairing is uncertain; numerical labels are not passage pairs','Only one sampling batch available','Author states derivation history is uncertain','T2T DNA versus GRCh38 RNA and specimen correspondence require a verified mapping','Cell-specific follow-ups are not part of the A549 baseline reprocessing')),'analysis_status')
  inputs <- c(file.path(REPO,cfg$raw_counts),file.path(REPO,cfg$gtf),cfg$sample_manifest,htpath)
  tab(data.frame(path=inputs,sha256=vapply(inputs,hash_file,character(1))),'input_hashes')
  list(effects=e,norm=n,logCPM=l,mapping=mapping,pca=pca,pca_context=pctx,pca_variance=pvar$core_variance,pca_context_variance=pc_context$sdev^2/sum(pc_context$sdev^2),correlations=corr,qc=qc,context=context,parent=parent,de=de,sensitivity=sens,abundance=abundance,pairs=pairs,comparison=comp)
}
