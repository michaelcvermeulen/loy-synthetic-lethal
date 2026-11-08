rank_effect <- function(v,symbols) {z<-setNames(v,symbols); z[order(-z,names(z))]}
run_fgsea <- function(sets,ranks,cfg,seed) {
  set.seed(seed)
  fg <- as.data.frame(fgsea::fgseaMultilevel(sets,ranks,minSize=cfg$pathway_min_size,maxSize=cfg$pathway_max_size,
    eps=cfg$fgsea_eps,nPermSimple=cfg$fgsea_nPermSimple,BPPARAM=BiocParallel::SerialParam()))
  fg$leadingEdge <- vapply(fg$leadingEdge,paste,character(1),collapse=';'); fg$padj <- NULL
  fg
}
du145_pathways <- function(x,e,cn,cfg) {
  sets <- readRDS('resources/msigdb_all_collections.rds'); all_sets <- do.call(c,unname(sets))
  metadata <- as.data.frame(data.table::fread('resources/msigdb_all_set_metadata.csv'))
  names(metadata)[match(c('gs_name','gs_collection','gs_subcollection'),names(metadata))] <- c('pathway','collection','subcollection')
  mapping <- e$mapping; ix <- match(mapping$gene_id,e$effects$gene_id)
  # Balanced passage-matched design: the group-mean difference equals the mean of pair differences.
  effects <- matrix(e$effects$mean_difference[ix],ncol=1,dimnames=list(mapping$gene_name,'LOY_vs_sgCTRL'))
  tab(cbind(mapping[,c('gene_id','gene_name')],effects),'pathway_TF_effect_ranking')
  coverage <- merge(metadata,data.frame(pathway=names(all_sets),resource_size=lengths(all_sets),
    overlap_size=vapply(all_sets,function(g)length(intersect(g,rownames(effects))),integer(1))),by='pathway',sort=FALSE)
  coverage$eligible <- coverage$overlap_size>=cfg$pathway_min_size & coverage$overlap_size<=cfg$pathway_max_size
  tab(coverage,'pathway_resource_coverage'); eligible <- all_sets[coverage$pathway[coverage$eligible]]
  message('fgsea: ',length(eligible),' eligible sets')
  fg <- run_fgsea(eligible,rank_effect(effects[,1],rownames(effects)),cfg,cfg$seed)
  fg$contrast <- 'LOY_vs_sgCTRL'; fg$p_at_lower_bound <- !is.na(fg$pval) & fg$pval<=cfg$fgsea_eps
  fg <- enrichment_adjust(merge(fg,metadata,by='pathway',sort=FALSE),'pval')
  fg$statistical_scope <- 'exploratory_preranked_gene_set_null_not_population_replication'
  # Copy-number screen: leading-edge composition and a rerun without CN-changed genes.
  cng <- cn$genes
  cn_by_symbol <- setNames(cng$CN_changed[match(mapping$gene_id,cng$gene_id)],mapping$gene_name)
  arm_by_symbol <- setNames(cng$t2t_arm[match(mapping$gene_id,cng$gene_id)],mapping$gene_name)
  le <- strsplit(fg$leadingEdge,';',fixed=TRUE)
  fg$leading_edge_size <- lengths(le)
  fg$leading_edge_CN_changed_fraction <- vapply(le,function(g) mean(cn_by_symbol[g] %in% TRUE),numeric(1))
  fg$leading_edge_top_arm <- vapply(le,function(g) {t<-table(arm_by_symbol[g]); if(!length(t)) NA_character_ else names(t)[which.max(t)]},character(1))
  fg$leading_edge_top_arm_fraction <- vapply(le,function(g) {t<-table(arm_by_symbol[g]); if(!length(t)) NA_real_ else max(t)/length(g)},numeric(1))
  stable <- !(cn_by_symbol %in% TRUE)
  message('fgsea without CN-changed genes: ',sum(!stable),' symbols removed')
  fgs <- run_fgsea(eligible,rank_effect(effects[stable,1],rownames(effects)[stable]),cfg,cfg$seed+10)
  fgs <- fgs[,c('pathway','NES','pval','size')]; names(fgs) <- c('pathway','NES_CN_stable_only','pval_CN_stable_only','size_CN_stable_only')
  fgs$BH_CN_stable_only <- p.adjust(fgs$pval_CN_stable_only,'BH')
  fg <- merge(fg,fgs,by='pathway',all.x=TRUE,sort=FALSE)
  fg$survives_CN_screen <- fg$BH_across_collections<0.05 & is.finite(fg$BH_CN_stable_only) & fg$BH_CN_stable_only<0.05 & sign(fg$NES)==sign(fg$NES_CN_stable_only)
  # Sub-threshold arm gains (chrX +0.18 in LOY-2X) and gene clusters (6p22 histones) escape the
  # gene-level CN screen; flag sets whose leading edge is concentrated on one arm, and C1 positional sets.
  fg$positionally_concentrated <- fg$collection=='C1' | (is.finite(fg$leading_edge_top_arm_fraction) & fg$leading_edge_top_arm_fraction>=0.5)
  fg$interpretable <- fg$survives_CN_screen & !fg$positionally_concentrated
  tab(fg,'fgsea_all_collections'); tab(fg[fg$collection=='H' | fg$subcollection %in% c('CP:REACTOME','GO:BP'),],'fgsea_hallmark_reactome_gobp')
  # CAMERA with the paired design and the same symbol universe.
  counts <- x$counts[mapping$gene_id,x$core$sample_key,drop=FALSE]; rownames(counts) <- mapping$gene_name
  design <- e$de$design
  v <- limma::voom(counts,design,lib.size=e$norm$effective,normalize.method='none',plot=FALSE)
  index <- limma::ids2indices(eligible,rownames(v$E),remove.empty=FALSE)
  cameras <- list()
  for (mode in c('fixed_0.01','estimated_nonnegative')) {
    rho <- if (mode=='fixed_0.01') cfg$camera_inter_gene_cor else NA_real_
    z <- limma::camera(v,index,design=design,contrast=which(colnames(design)=='conditionLOY'),inter.gene.cor=rho,allow.neg.cor=FALSE,use.ranks=FALSE,sort=FALSE)
    z$pathway <- rownames(z); rownames(z) <- NULL; z$FDR <- NULL; z$contrast <- 'LOY_vs_sgCTRL'; z$correlation_mode <- mode
    if (!'Correlation' %in% names(z)) z$Correlation <- rho
    z$residual_df <- nrow(design)-ncol(design)
    z <- enrichment_adjust(merge(z,metadata,by='pathway',sort=FALSE),'PValue',c('contrast','correlation_mode'))
    z$statistical_scope <- 'exploratory_paired_passage_model'
    cameras[[mode]] <- z
  }
  cam <- cameras$fixed_0.01; cam_sens <- cameras$estimated_nonnegative
  tab(cam,'camera_all_collections'); tab(cam_sens,'camera_estimated_correlation_sensitivity')
  comp <- merge(fg[,c('pathway','collection','subcollection','NES','pval','BH_across_collections','survives_CN_screen','interpretable')],
                cam[,c('pathway','Direction','PValue','BH_across_collections')],by='pathway',suffixes=c('_fgsea','_camera'))
  comp$direction_agrees <- (comp$NES>0)==(comp$Direction=='Up')
  comp$both_BH05 <- comp$BH_across_collections_fgsea<0.05 & comp$BH_across_collections_camera<0.05
  tab(comp,'fgsea_camera_comparison')
  summary <- do.call(rbind,lapply(split(comp,comp$collection),function(z) data.frame(collection=z$collection[1],tested=nrow(z),
    fgsea_BH05=sum(z$BH_across_collections_fgsea<0.05,na.rm=TRUE),fgsea_BH05_survive_CN=sum(z$survives_CN_screen,na.rm=TRUE),fgsea_interpretable=sum(z$interpretable,na.rm=TRUE),
    camera_BH05=sum(z$BH_across_collections_camera<0.05,na.rm=TRUE),both_direction_agree=sum(z$both_BH05 & z$direction_agrees,na.rm=TRUE))))
  tab(summary,'enrichment_collection_summary')
  set.seed(cfg$seed+20); keep_rp <- !mapping$ribosomal_protein
  hsens <- run_fgsea(sets$H,rank_effect(effects[keep_rp,1],rownames(effects)[keep_rp]),cfg,cfg$seed+20)
  tab(hsens,'fgsea_hallmark_no_ribosomal_sensitivity')
  # CollecTRI target patterns and PROGENy footprints, as in HT1080/A549.
  net <- unique(as.data.frame(data.table::fread('resources/collectri_network.csv'))[,c('source','target','mor')])
  key <- paste(net$source,net$target,sep='\r'); conflict <- duplicated(key)|duplicated(key,fromLast=TRUE)
  net <- net[!conflict & is.finite(net$mor),]
  tf <- as.data.frame(decoupleR::run_ulm(effects,net,.source='source',.target='target',.mor='mor',minsize=cfg$tf_min_targets))
  names(tf)[names(tf)=='condition'] <- 'contrast'; tf$BH_target_pattern <- p.adjust(tf$p_value,'BH')
  stable_eff <- effects[stable,,drop=FALSE]
  tfs <- as.data.frame(decoupleR::run_ulm(stable_eff,net,.source='source',.target='target',.mor='mor',minsize=cfg$tf_min_targets))
  tf <- merge(tf,data.frame(source=tfs$source,score_CN_stable_only=tfs$score,p_CN_stable_only=tfs$p_value),by='source',all.x=TRUE)
  tf$regulator_expression_mean_difference <- e$effects$mean_difference[match(tf$source,e$effects$gene_name)]
  tf$statistical_scope <- 'conditional_target_pattern_not_measured_regulator_activity'
  tab(tf,'collectri_ULM_effect_signatures')
  lc <- e$logCPM[mapping$gene_id,,drop=FALSE]; rownames(lc) <- mapping$gene_name
  centered <- lc-rowMeans(lc)
  model <- as.matrix(readRDS('resources/progeny_model_top100.rds')); common <- intersect(rownames(model),rownames(lc))
  weights <- model[common,,drop=FALSE]; weights <- sweep(weights,2,colSums(abs(weights)),'/')
  prog <- t(centered[common,,drop=FALSE])%*%weights
  tab(matrix_long(t(prog),'weighted_centered_logCPM_score','pathway'),'progeny_sample_footprint_scores')
  pd <- matched_differences(t(prog),x$pairs)
  tab(data.frame(pathway=colnames(prog),LOY_minus_CTRL=group_difference(t(prog),x$core),pd,pairs_same_sign=rowSums(sign(pd)==sign(group_difference(t(prog),x$core)))),'progeny_group_footprint_differences')
  list(fgsea=fg,camera=cam,camera_sensitivity=cam_sens,comparison=comp,summary=summary,tf=tf,coverage=coverage)
}
