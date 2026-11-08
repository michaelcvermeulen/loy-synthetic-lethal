rank_effect <- function(v,symbols) {z<-setNames(v,symbols); z[order(-z,names(z))]}
a549_pathways <- function(x,e,cfg) {
  sets<-readRDS('resources/msigdb_all_collections.rds'); all_sets<-do.call(c,unname(sets))
  metadata<-as.data.frame(data.table::fread('resources/msigdb_all_set_metadata.csv'))
  names(metadata)[match(c('gs_name','gs_collection','gs_subcollection'),names(metadata))]<-c('pathway','collection','subcollection')
  mapping<-e$mapping; ix<-match(mapping$gene_id,e$effects$gene_id)
  effects<-matrix(e$effects$mean_difference[ix],ncol=1,dimnames=list(mapping$gene_name,'LOY_vs_sgCTRL'))
  tab(cbind(mapping[,c('gene_id','gene_name')],effects),'pathway_TF_effect_ranking')
  coverage<-merge(metadata,data.frame(pathway=names(all_sets),resource_size=lengths(all_sets),overlap_size=vapply(all_sets,function(g)length(intersect(g,rownames(effects))),integer(1))),by='pathway',sort=FALSE)
  coverage$eligible<-coverage$overlap_size>=cfg$pathway_min_size & coverage$overlap_size<=cfg$pathway_max_size
  coverage$status<-ifelse(coverage$eligible,'eligible',ifelse(coverage$overlap_size<cfg$pathway_min_size,'below_minimum_overlap','above_maximum_overlap'))
  tab(coverage,'pathway_resource_coverage'); eligible<-all_sets[coverage$pathway[coverage$eligible]]
  set.seed(cfg$seed)
  message('fgsea: ',length(eligible),' eligible sets across ',length(sets),' collections')
  fg<-as.data.frame(fgsea::fgseaMultilevel(eligible,rank_effect(effects[,1],rownames(effects)),minSize=cfg$pathway_min_size,maxSize=cfg$pathway_max_size,eps=cfg$fgsea_eps,nPermSimple=cfg$fgsea_nPermSimple,BPPARAM=BiocParallel::SerialParam()))
  fg$leadingEdge<-vapply(fg$leadingEdge,paste,character(1),collapse=';'); fg$padj<-NULL
  fg$contrast<-'LOY_vs_sgCTRL'; fg$p_at_lower_bound<-!is.na(fg$pval)&fg$pval<=cfg$fgsea_eps
  fg$test_status<-ifelse(is.finite(fg$pval),'tested','numerical_failure')
  fg<-enrichment_adjust(merge(fg,metadata,by='pathway',sort=FALSE),'pval')
  fg$statistical_scope<-'exploratory_preranked_gene_set_null_not_population_replication'
  tab(fg,'fgsea_all_collections');tab(fg[fg$collection=='H' | fg$subcollection=='CP:REACTOME',],'fgsea_hallmark_reactome')
  # CAMERA uses the same symbol universe and core normalization exposures.
  counts<-x$counts[mapping$gene_id,x$core$sample_key,drop=FALSE];rownames(counts)<-mapping$gene_name
  design<-e$de$design
  v<-limma::voom(counts,design,lib.size=e$norm$effective,normalize.method='none',plot=FALSE)
  index<-limma::ids2indices(eligible,rownames(v$E),remove.empty=FALSE)
  assert(all(lengths(index)>=cfg$pathway_min_size & lengths(index)<=cfg$pathway_max_size),'CAMERA/fgsea universe mismatch')
  cameras<-list()
  for(mode in c('fixed_0.01','estimated_nonnegative')) {
    message('CAMERA ',mode)
    rho<-if(mode=='fixed_0.01')cfg$camera_inter_gene_cor else NA_real_
    z<-limma::camera(v,index,design=design,contrast=which(colnames(design)=='conditionLOY'),inter.gene.cor=rho,allow.neg.cor=FALSE,use.ranks=FALSE,sort=FALSE)
    z$pathway<-rownames(z);rownames(z)<-NULL;z$FDR<-NULL;z$contrast<-'LOY_vs_sgCTRL';z$correlation_mode<-mode
    if(!'Correlation'%in%names(z))z$Correlation<-rho
    z$residual_df<-nrow(design)-ncol(design)
    z<-enrichment_adjust(merge(z,metadata,by='pathway',sort=FALSE),'PValue',c('contrast','correlation_mode'))
    z$statistical_scope<-'exploratory_unpaired_model_independent_residual_assumption_derivation_uncertain'
    cameras[[mode]]<-z
  }
  cam<-cameras$fixed_0.01;cam_sens<-cameras$estimated_nonnegative
  tab(cam,'camera_all_collections');tab(cam_sens,'camera_estimated_correlation_sensitivity')
  saveRDS(list(voom=v,design=design,universe=mapping), 'outputs/models/camera_voom_model.rds')
  comp<-merge(fg[,c('pathway','collection','subcollection','NES','pval','BH_across_collections')],cam[,c('pathway','Direction','PValue','BH_across_collections')],by='pathway',suffixes=c('_fgsea','_camera'))
  comp$direction_agrees<-(comp$NES>0)==(comp$Direction=='Up')
  comp$both_BH05<-comp$BH_across_collections_fgsea<0.05 & comp$BH_across_collections_camera<0.05
  tab(comp,'fgsea_camera_comparison')
  summary<-do.call(rbind,lapply(split(comp,comp$collection),function(z)data.frame(collection=z$collection[1],tested=nrow(z),fgsea_BH05=sum(z$BH_across_collections_fgsea<0.05,na.rm=TRUE),camera_BH05=sum(z$BH_across_collections_camera<0.05,na.rm=TRUE),both_direction_agree=sum(z$both_BH05 & z$direction_agrees,na.rm=TRUE))))
  tab(summary,'enrichment_collection_summary')
  # Hallmark sensitivity excludes ribosomal proteins without changing original global BH.
  set.seed(cfg$seed+20); keep_rp<-!mapping$ribosomal_protein
  hsens<-as.data.frame(fgsea::fgseaMultilevel(sets$H,rank_effect(effects[keep_rp,1],rownames(effects)[keep_rp]),minSize=cfg$pathway_min_size,maxSize=cfg$pathway_max_size,eps=cfg$fgsea_eps,nPermSimple=cfg$fgsea_nPermSimple,BPPARAM=BiocParallel::SerialParam()))
  hsens$leadingEdge<-vapply(hsens$leadingEdge,paste,character(1),collapse=';')
  tab(hsens,'fgsea_hallmark_no_ribosomal_sensitivity')
  lc<-e$logCPM[mapping$gene_id,,drop=FALSE];rownames(lc)<-mapping$gene_name
  percentiles<-apply(lc,2,function(v)(rank(v,ties.method='average')-1)/(length(v)-1))
  hs<-t(vapply(sets$H,function(g)colMeans(percentiles[intersect(g,rownames(lc)),,drop=FALSE]),numeric(ncol(lc))));colnames(hs)<-colnames(lc)
  tab(matrix_long(hs,'mean_within_sample_percentile','pathway'),'hallmark_sample_rank_scores')
  tab(data.frame(pathway=rownames(hs),LOY_minus_CTRL=group_difference(hs,x$core)),'hallmark_group_rank_score_differences')
  net<-unique(as.data.frame(data.table::fread('resources/collectri_network.csv'))[,c('source','target','mor')])
  key<-paste(net$source,net$target,sep='\r');conflict<-duplicated(key)|duplicated(key,fromLast=TRUE)
  tab(net[conflict,],'collectri_excluded_conflicting_edges');net<-net[!conflict & is.finite(net$mor),]
  tc<-aggregate(target~source,net,function(z)length(unique(z)));names(tc)[2]<-'resource_targets'
  tc2<-aggregate(target~source,net[net$target%in%rownames(effects),],function(z)length(unique(z)));names(tc2)[2]<-'expressed_targets'
  tc<-merge(tc,tc2,all.x=TRUE);tc$expressed_targets[is.na(tc$expressed_targets)]<-0;tc$eligible<-tc$expressed_targets>=cfg$tf_min_targets
  tab(tc,'collectri_target_coverage')
  message('CollecTRI ULM and PROGENy footprints')
  tf<-as.data.frame(decoupleR::run_ulm(effects,net,.source='source',.target='target',.mor='mor',minsize=cfg$tf_min_targets))
  names(tf)[names(tf)=='condition']<-'contrast';tf$BH_target_pattern<-p.adjust(tf$p_value,'BH');tf<-merge(tf,tc,by='source',all.x=TRUE)
  tf$statistical_scope<-'conditional_target_pattern_not_measured_regulator_activity'
  tab(tf,'collectri_ULM_effect_signatures')
  centered<-lc-rowMeans(lc)
  ts<-as.data.frame(decoupleR::run_ulm(centered,net,.source='source',.target='target',.mor='mor',minsize=cfg$tf_min_targets))
  tab(ts,'collectri_ULM_sample_scores_centered')
  target<-merge(net,data.frame(target=rownames(effects),gene_id=mapping$gene_id,mean_difference=effects[,1]),by='target')
  target$signed_effect<-target$mor*target$mean_difference;tab(target,'collectri_target_effect_evidence')
  model<-as.matrix(readRDS('resources/progeny_model_top100.rds'));common<-intersect(rownames(model),rownames(lc))
  weights<-model[common,,drop=FALSE];denom<-colSums(abs(weights));assert(all(denom>0),'PROGENy empty overlap')
  weights<-sweep(weights,2,denom,'/');prog<-t(centered[common,,drop=FALSE])%*%weights
  tab(matrix_long(t(prog),'weighted_centered_logCPM_score','pathway'),'progeny_sample_footprint_scores')
  tab(data.frame(pathway=colnames(prog),LOY_minus_CTRL=group_difference(t(prog),x$core)),'progeny_group_footprint_differences')
  tab(data.frame(pathway=colnames(weights),covered_targets=colSums(weights!=0),model_targets=colSums(model!=0),absolute_weight_sum=denom),'progeny_coverage')
  return(list(fgsea=fg,camera=cam,camera_sensitivity=cam_sens,comparison=comp,summary=summary,tf=tf,hallmark_scores=hs,coverage=coverage))
}
