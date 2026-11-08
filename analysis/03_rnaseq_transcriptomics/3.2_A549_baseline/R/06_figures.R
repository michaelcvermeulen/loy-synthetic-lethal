COLORS<-c(CTRL='#0072B2',LOY='#D55E00',PARENTAL='#777777')
save_figure<-function(p,stem,width=9,height=6,caption='') {
  p<-p+ggplot2::theme(plot.title=ggplot2::element_text(face='bold'),panel.grid.minor=ggplot2::element_blank())
  for(ext in c('pdf','svg','png')) {
    destination<-paste0('outputs/figures/',stem,'.',ext)
    temporary<-tempfile(pattern='figure_',tmpdir='outputs/figures',fileext=paste0('.',ext))
    ggplot2::ggsave(temporary,p,width=width,height=height,dpi=600,bg='white',limitsize=FALSE)
    stopifnot(file.exists(temporary),file.info(temporary)$size>1000)
    if(file.exists(destination))stopifnot(file.remove(destination))
    stopifnot(file.rename(temporary,destination))
  }
  writeLines(c(paste0('**',stem,'**'),'',caption),paste0('outputs/captions/',stem,'.md'))
}
core_figures<-function(x,e,cfg) {
  suppressPackageStartupMessages(library(ggplot2))
  theme_set(theme_bw(base_size=11))
  qc<-e$qc
  q<-data.table::melt(data.table::as.data.table(qc),id.vars=c('sample_key','treatment'),measure.vars=c('assigned_fraction','mitochondrial_fraction','MSY_fraction'),variable.name='metric',value.name='fraction')
  p<-ggplot(q,aes(sample_key,100*fraction,fill=treatment))+geom_col()+facet_wrap(~metric,scales='free_y',ncol=1)+scale_fill_manual(values=COLORS)+labs(x=NULL,y='Percent of reads',title='A549 library composition')+theme(axis.text.x=element_text(angle=45,hjust=1))
  save_figure(p,'Fig01_library_QC',10,8,'Eleven available libraries; no samples excluded. Assignment is relative to all featureCounts assignment categories; MSY and mitochondrial fractions are relative to assigned reads. Composition differences do not by themselves identify RNA degradation or another cause.')
  p<-ggplot(e$pca,aes(PC1,PC2,color=treatment,label=sample_index))+geom_point(size=3)+geom_text(vjust=-1,show.legend=FALSE)+scale_color_manual(values=COLORS)+labs(title='A549 core expression structure',x=sprintf('PC1 (%.1f%%)',100*e$pca_variance[1]),y=sprintf('PC2 (%.1f%%)',100*e$pca_variance[2]))
  save_figure(p,'Fig02_PCA_core',7,5.5,'PCA uses eight core libraries and all retained autosomal/non-PAR X genes, TMM log2 CPM (prior count 2), centred without variance scaling. Labels are sample identifiers, not confirmed passage pairs. Parentals do not enter this normalization or PCA.')
  p<-ggplot(e$pca,aes(PC1,PC3,color=treatment,label=sample_index))+geom_hline(yintercept=0,color='#D9D9D9',linewidth=.35)+geom_vline(xintercept=0,color='#D9D9D9',linewidth=.35)+geom_point(size=3)+geom_text(vjust=-1,show.legend=FALSE)+scale_x_continuous(expand=expansion(mult=c(.06,.08)))+scale_y_continuous(expand=expansion(mult=c(.08,.14)))+scale_color_manual(values=COLORS,labels=c(CTRL='sgCTRL',LOY='LOY'))+labs(title='A549 LOY versus sgCTRL expression structure',x=sprintf('PC1 (%.1f%%)',100*e$pca_variance[1]),y=sprintf('PC3 (%.1f%%)',100*e$pca_variance[3]),color=NULL)
  save_figure(p,'SuppFig04_PCA_PC1_PC3_core',7,5.5,'PC1 versus PC3 for the eight A549 core libraries: four LOY and four sgCTRL. PCA uses all retained autosomal/non-PAR X genes, TMM log2 CPM (prior count 2), centred without variance scaling. Labels 1-4 are sample identifiers; they are not treated as confirmed passage pairs. Parentals are excluded.')
  p<-ggplot(e$pca_context,aes(PC1,PC2,color=treatment,label=sample_index))+geom_point(size=3)+geom_text(vjust=-1,show.legend=FALSE)+scale_color_manual(values=COLORS)+labs(title='A549 parental context',x=sprintf('PC1 (%.1f%%)',100*e$pca_context_variance[1]),y=sprintf('PC2 (%.1f%%)',100*e$pca_context_variance[2]))
  save_figure(p,'SuppFig01_PCA_parental_context',7,5.5,'Contextual PCA has a separate all-11 normalization and the fixed core gene universe. These values do not alter primary effects. Differences between parental and nucleofected populations cannot be uniquely assigned to editing or clonal drift.')
  de<-e$de$table;de$FDR_status<-ifelse(de$FDR<.05,'BH < 0.05','BH >= 0.05')
  p<-ggplot(de,aes(logFC,-log10(pmax(PValue,1e-300)),color=FDR_status))+geom_point(size=.6,alpha=.6)+scale_color_manual(values=c('BH < 0.05'='#D55E00','BH >= 0.05'='#999999'))+geom_vline(xintercept=0,linetype=2)+labs(title='A549 exploratory unpaired expression contrast',x='Voom log2 fold change (LOY minus sgCTRL)',y='-log10 nominal P value')
  save_figure(p,'Fig03_unpaired_DE',8,6,'Voom-limma with condition only, robust empirical Bayes, and core TMM exposures. BH covers all primary genes. Pairing and derivation history are uncertain; inferential results assume independent model residuals and do not establish independent editing replication.')
  eff<-e$effects[e$effects$primary_eligible,];eff$status<-ifelse(eff$descriptive_candidate,'Magnitude >= 0.5; stable omission sign','Other')
  p<-ggplot(eff,aes(mean_logCPM,mean_difference,color=status))+geom_point(size=.5,alpha=.5)+geom_hline(yintercept=c(-.5,0,.5),linetype=c(2,1,2))+scale_color_manual(values=c('#D55E00','#999999'))+labs(title='A549 descriptive expression differences',x='Mean log2 CPM',y='Mean log2 CPM difference (LOY minus sgCTRL)')
  save_figure(p,'Fig04_expression_effects',9,6,'Primary descriptive group-mean logCPM differences; candidates require absolute difference >= 0.5 and the same direction after omission of each of eight samples, with filtering and TMM refitted for each omission. This is robustness within the measured libraries, not independent validation.')
  ys<-e$abundance[e$abundance$chrom=='chrY',]
  p<-ggplot(ys,aes(treatment,log2(CPM+1),color=treatment))+geom_point(position=position_jitter(width=.08,seed=42),size=1.8)+facet_wrap(~gene_name,ncol=4,scales='free_y')+scale_color_manual(values=COLORS)+labs(title='Y-linked transcript abundance',x=NULL,y='log2(core-scaled CPM + 1)')
  save_figure(p,'Fig05_Y_expression',10,6,'Each point is one measured library, with the same core normalization used for Y and non-Y genes. Residual or undetected RNA does not establish a whole-chromosome genotype or a Y-negative cell fraction.')
  pair<-e$pairs;pp<-rbind(data.frame(pair=paste(pair$X,pair$Y,sep=' / '),partner='X',effect=pair$X_mean_difference),data.frame(pair=paste(pair$X,pair$Y,sep=' / '),partner='Y',effect=pair$Y_mean_difference))
  p<-ggplot(pp,aes(effect,reorder(pair,effect),color=partner))+geom_point(size=2.5)+geom_vline(xintercept=0,linetype=2)+scale_color_manual(values=c(X='#0072B2',Y='#D55E00'))+labs(title='X/Y paralog expression responses',x='Mean log2 CPM difference (LOY minus sgCTRL)',y=NULL)
  save_figure(p,'Fig06_XY_paralogs',9,6,'Separate X- and Y-linked RNA effects. Expression changes do not establish protein compensation or functional buffering. The accompanying table reports a length-normalised summed-RNA proxy with the same limitation.')
  cross<-e$comparison
  p<-ggplot(cross,aes(mean_difference_HT1080,mean_difference_A549))+geom_point(size=.4,alpha=.3,color='#0072B2')+geom_hline(yintercept=0,linetype=2)+geom_vline(xintercept=0,linetype=2)+labs(title='A549 and HT1080 descriptive effects',x='HT1080 mean log2 CPM difference',y='A549 mean log2 CPM difference')
  save_figure(p,'SuppFig02_cross_line_effects',7,6,'Exact shared gene IDs in independently normalized primary backgrounds. HT1080 averages two sampling periods; A549 has one. Association across genes does not establish common causation or equal regional Y depletion.')
  sensitivity<-data.frame(primary=eff$mean_difference,RLE=eff$RLE_mean_difference,LOO_change=eff$LOO_max_abs_change)
  p<-ggplot(sensitivity,aes(primary,LOO_change))+geom_point(size=.5,alpha=.3,color='#0072B2')+labs(title='Sensitivity to individual A549 libraries',x='Primary mean log2 CPM difference',y='Maximum absolute change after one-sample omission')
  save_figure(p,'SuppFig03_leave_one_out',7,5.5,'Each omission refits filtering and normalization while retaining the original gene IDs for comparison. No library is excluded from the primary analysis on the basis of these results.')
  invisible(TRUE)
}
pathway_figures<-function(p) {
  suppressPackageStartupMessages(library(ggplot2))
  theme_set(theme_bw(base_size=11))
  h<-p$fgsea[p$fgsea$collection=='H',];h<-h[order(h$NES),]
  h$pathway<-factor(sub('^HALLMARK_','',h$pathway),levels=sub('^HALLMARK_','',h$pathway))
  h$global_BH05<-h$BH_across_collections<.05
  g<-ggplot(h,aes(NES,pathway,color=global_BH05))+geom_point(size=2)+geom_vline(xintercept=0,linetype=2)+scale_color_manual(values=c('FALSE'='#999999','TRUE'='#D55E00'))+labs(title='A549 Hallmark enrichment',x='Normalized enrichment score (effect-ranked fgsea)',y=NULL,color='Global BH < 0.05')
  save_figure(g,'Fig07_Hallmark_enrichment',10,11,'All 50 Hallmark sets are shown. Highlighting uses BH across every eligible MSigDB set in all collections, not a Hallmark-only test family. Rank enrichment is distinct from evidence for independent population replication.')
  tf<-p$tf;tf<-head(tf[order(tf$BH_target_pattern,-abs(tf$score)),],25);tf$source<-factor(tf$source,levels=tf$source[order(tf$score)])
  g<-ggplot(tf,aes(score,source,color=score>0))+geom_point(size=2.5)+geom_vline(xintercept=0,linetype=2)+scale_color_manual(values=c('FALSE'='#0072B2','TRUE'='#D55E00'))+labs(title='A549 CollecTRI target-pattern scores',x='ULM score',y=NULL,color='Positive score')
  save_figure(g,'Fig08_TF_target_patterns',9,8,'The 25 lowest-BH regulator target patterns, with absolute score breaking ties. At least five represented targets are required. These are transcript-based target-pattern estimates, not direct measurements of transcription-factor activity.')
  invisible(TRUE)
}
