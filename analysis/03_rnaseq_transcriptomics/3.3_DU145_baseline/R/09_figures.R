COLORS <- c(CTRL='#0072B2',LOY='#D55E00')
save_figure <- function(p,stem,width=9,height=6,caption='') {
  th <- ggplot2::theme(plot.title=ggplot2::element_text(face='bold'),panel.grid.minor=ggplot2::element_blank())
  p <- if (inherits(p,'patchwork')) patchwork:::`&.gg`(p,th) else p+th
  for (ext in c('pdf','svg','png')) {
    destination <- paste0('outputs/figures/',stem,'.',ext)
    temporary <- tempfile(pattern='figure_',tmpdir='outputs/figures',fileext=paste0('.',ext))
    ggplot2::ggsave(temporary,p,width=width,height=height,dpi=600,bg='white',limitsize=FALSE)
    stopifnot(file.exists(temporary),file.info(temporary)$size>1000)
    if (file.exists(destination)) stopifnot(file.remove(destination))
    stopifnot(file.rename(temporary,destination))
  }
  writeLines(c(paste0('**',stem,'**'),'',caption),paste0('outputs/captions/',stem,'.md'))
}
core_figures <- function(x,e,id,y,cn,cfg) {
  suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
  theme_set(theme_bw(base_size=11))
  lab <- c(CTRL='sgCTRL',LOY='LOY (CP + Q53 + P47 2x)')
  q <- data.table::melt(data.table::as.data.table(e$qc),id.vars=c('sample_key','treatment'),measure.vars=c('assigned_fraction','mitochondrial_fraction','MSY_fraction'),variable.name='metric',value.name='fraction')
  p <- ggplot(q,aes(sample_key,100*fraction,fill=treatment))+geom_col()+facet_wrap(~metric,scales='free_y',ncol=1)+scale_fill_manual(values=COLORS,labels=lab)+labs(x=NULL,y='Percent of reads',fill=NULL,title='DU145 library composition')+theme(axis.text.x=element_text(angle=45,hjust=1))
  save_figure(p,'Fig01_library_QC',9,8,'Six libraries, three sgCTRL and three LOY at matched passages P12-P14; none excluded. Assignment is relative to all featureCounts categories; mitochondrial and MSY fractions are relative to assigned reads. Source labels read LNCaP; the RNA identifies the line as DU145.')
  pca <- e$pca
  p <- ggplot(pca,aes(PC1,PC2,color=treatment,shape=passage_block))+geom_point(size=3.2)+scale_color_manual(values=COLORS,labels=lab)+
    labs(color=NULL,shape='Passage',title='DU145 core expression structure',x=sprintf('PC1 (%.1f%%)',100*e$pca_variance[1]),y=sprintf('PC2 (%.1f%%)',100*e$pca_variance[2]))
  save_figure(p,'Fig02_PCA_core',7,5.5,sprintf('PCA of all %s retained autosomal/non-PAR X genes, TMM log2 CPM (prior count 2), centred without scaling. Passage explains PC1 (R2 = %.2f); condition explains PC2 (R2 = %.2f). Passage-matched pairs are therefore modelled as blocks.',format(sum(e$effects$primary_eligible),big.mark=','),e$pc_assoc$passage_R2[1],e$pc_assoc$treatment_R2[2]))
  r <- id$ranked[id$ranked$profile %in% c('sgCTRL_mean','LOY_mean'),]
  top <- r[r$rank<=10 | r$StrippedCellLineName %in% c('DU145','LNCAPCLONEFGC'),]
  top$label <- paste0(top$StrippedCellLineName,' (',top$OncotreeLineage,')')
  top$highlight <- ifelse(top$StrippedCellLineName=='DU145','DU145',ifelse(top$StrippedCellLineName=='LNCAPCLONEFGC','LNCaP','Other'))
  p1 <- ggplot(top,aes(spearman,reorder(label,spearman),color=highlight))+geom_point(size=2.2)+facet_wrap(~profile,ncol=2,labeller=as_labeller(c(sgCTRL_mean='sgCTRL mean',LOY_mean='LOY mean')))+
    geom_text(aes(label=paste0('#',rank)),hjust=-0.3,size=2.8,show.legend=FALSE)+scale_color_manual(values=c(DU145='#D55E00',LNCaP='#0072B2',Other='#777777'))+
    scale_x_continuous(expand=expansion(mult=c(.05,.15)))+labs(x='Spearman correlation with DepMap 26Q1 profile',y=NULL,color=NULL,title='Transcriptome identity')
  mk <- id$markers; mk$group <- ifelse(grepl('sgCTRL',mk$sample_key),'sgCTRL','LOY')
  p2 <- ggplot(mk,aes(gene_name,log2(TPM+1),color=group))+geom_point(position=position_jitter(width=.12,seed=1),size=1.6)+scale_color_manual(values=c(sgCTRL='#0072B2',LOY='#D55E00'))+labs(x=NULL,y='log2(TPM + 1)',color=NULL,title='Prostate lineage and epithelial markers')+theme(axis.text.x=element_text(angle=45,hjust=1))
  save_figure(p1/p2+plot_layout(heights=c(1.3,1)),'Fig03_identity',10,9,sprintf('Spearman correlation of library log2(TPM + 1) with default DepMap 26Q1 profiles over the %s most variable shared protein-coding genes (hash-verified release files). Top ten matches plus DU145 and LNCaP are shown. AR, KLK3, KLK2, FOLH1 and HOXB13 are shown with KRT/VIM context. This supports DU145 identity from RNA; STR profiling was not performed here.',format(id$genes_used,big.mark=',')))
  # Chromosome Y: RNA (percent of sgCTRL per passage pair) against lpWGS segment ratio, on one axis.
  yt <- y$table; yt <- yt[yt$expressed & is.finite(yt$t2t_mid_Mb) & yt$CTRL_mean_CPM>=1,]
  long <- do.call(rbind,lapply(x$pairs$block,function(b) data.frame(gene_name=yt$gene_name,region=yt$region,Mb=yt$t2t_mid_Mb,pair=b,percent=yt[[paste0('percent_of_sgCTRL_',b)]])))
  mean_pts <- data.frame(gene_name=yt$gene_name,Mb=yt$t2t_mid_Mb,percent=yt$percent_of_sgCTRL,region=yt$region)
  sd_dir <- file.path(dirname(REPO),'SEPT_11_2026_THESIS/figures/chapter3/source_data/du145_lpwgs_confirmation')
  segs <- read.csv(file.path(sd_dir,'PAR1_Yp_saved_segments.csv')); reg <- read.csv(file.path(sd_dir,'regional_medians_verified.csv'))
  L <- segs[segs$condition=='LOY',]; C <- segs[segs$condition=='sgCTRL',]
  seg_w <- merge(data.frame(start_bp=L$start_bp,end_bp=L$end_bp,LOY=L$segment_relative_cn),
                 data.frame(cstart=C$start_bp,cend=C$end_bp,CTRL=C$segment_relative_cn),by=NULL)
  seg_w <- seg_w[seg_w$start_bp<=seg_w$cend & seg_w$end_bp>=seg_w$cstart,]
  seg_w$s <- pmax(seg_w$start_bp,seg_w$cstart)/1e6; seg_w$e2 <- pmin(seg_w$end_bp,seg_w$cend)/1e6
  seg_w$ratio <- 100*pmax(seg_w$LOY,0.001)/pmax(seg_w$CTRL,0.001)
  seg_w <- seg_w[seg_w$CTRL>0.0011,]  # ratio undefined where the control segment is at the numerical floor
  yq <- reg[reg$region=='MSY Yq (euchromatic)',]; yq_ratio <- 100*yq$median_relative_cn[yq$condition=='LOY']/yq$median_relative_cn[yq$condition=='sgCTRL']
  dna <- rbind(data.frame(s=seg_w$s,e=seg_w$e2,ratio=seg_w$ratio),data.frame(s=10.88,e=27.45,ratio=yq_ratio))
  shade <- data.frame(xmin=c(0,2.45832,4.8),xmax=c(2.45832,4.8,27.45),lab=c('PAR1','retained distal Yp','depleted MSY'))
  p <- ggplot()+geom_rect(data=shade,aes(xmin=xmin,xmax=xmax,ymin=0.8,ymax=250,fill=lab),alpha=.12)+
    geom_segment(data=dna,aes(x=s,xend=e,y=ratio,yend=ratio),color='#333333',linewidth=1.1)+
    geom_point(data=long[long$Mb<30,],aes(Mb,pmax(percent,1)),color='#D55E00',alpha=.45,size=1.4)+
    geom_point(data=mean_pts[mean_pts$Mb<30,],aes(Mb,pmax(percent,1)),color='#D55E00',size=2.6,shape=21,fill='white',stroke=1)+
    ggrepel::geom_text_repel(data=mean_pts[mean_pts$Mb<30 & (mean_pts$region=='MSY' | mean_pts$gene_name %in% c('CD99','SLC25A6','P2RY8','CRLF2')),],aes(Mb,pmax(percent,1),label=gene_name),size=3,min.segment.length=0,seed=3,max.overlaps=30)+
    geom_hline(yintercept=100,linetype=2,color='#777777')+scale_y_log10(breaks=c(1,5,10,25,50,100,200),limits=c(0.8,250))+
    scale_fill_manual(values=c(PAR1='#56B4E9','retained distal Yp'='#009E73','depleted MSY'='#999999'))+
    labs(x='CHM13v2.0 chromosome Y position (Mb)',y='LOY as percent of sgCTRL (log scale)',fill=NULL,title='DU145 chromosome Y: expression follows retained and depleted DNA')
  save_figure(p,'Fig04_Y_positional_RNA_DNA',11,6,'Orange: RNA for genes with sgCTRL mean CPM >= 1 (small points: each passage pair, LOY/sgCTRL CPM; open circles: group means). Black segments: lpWGS LNCaP-LOY-2X/sgCTRL ratio of saved QDNAseq segment values (0-10.9 Mb) and the euchromatic Yq regional median, both floored at 0.001 relative copy number. RNA values at or below 1% are drawn at 1%. PAR RNA counts are chrX-assigned and plotted at chrY PAR coordinates; PAR RNA cannot be attributed to X or Y. RNA and lpWGS specimens were not matched; the author identifies the RNA LOY population as the 2x population profiled by lpWGS.')
  pr <- y$pairs[is.finite(y$pairs$Y_CTRL_CPM) & y$pairs$Y_CTRL_CPM>=1,]
  pp <- rbind(data.frame(pair=paste(pr$X,pr$Y,sep=' / '),partner='X-linked',effect=pr$X_mean_difference),data.frame(pair=paste(pr$X,pr$Y,sep=' / '),partner='Y-linked',effect=pr$Y_mean_difference))
  p <- ggplot(pp,aes(effect,reorder(pair,effect),color=partner))+geom_vline(xintercept=0,linetype=2)+geom_point(size=2.8)+scale_color_manual(values=c('X-linked'='#0072B2','Y-linked'='#D55E00'))+labs(x='Mean log2 CPM difference (LOY minus sgCTRL)',y=NULL,color=NULL,title='X-Y paralog responses')
  save_figure(p,'Fig05_XY_paralogs',8,5,'Pairs whose Y member has sgCTRL mean CPM >= 1. X-linked partners are not significant in the inclusive paired fit; their small positive shifts are of the size expected from the chrX gain in LOY-2X lpWGS and are not evidence of compensation. RNA does not measure protein or function.')
  par <- y$per_sample[y$per_sample$region %in% c('PAR1','PAR2'),]
  keep_par <- y$table$gene_name[y$table$region %in% c('PAR1','PAR2') & y$table$CTRL_mean_CPM>=5]
  par <- par[par$gene_name %in% keep_par,]
  p <- ggplot(par,aes(treatment,CPM,color=treatment,group=passage_block))+geom_line(color='#BBBBBB')+geom_point(size=1.8)+facet_wrap(~region+gene_name,scales='free_y',ncol=5)+scale_color_manual(values=COLORS,labels=lab)+labs(x=NULL,y='CPM (core TMM scaling)',color=NULL,title='DU145 pseudoautosomal expression')
  save_figure(p,'Fig06_PAR_expression',12,7,'PAR genes with sgCTRL mean CPM >= 5; grey lines join passage-matched libraries. lpWGS shows PAR1 at about 71% of sgCTRL relative copy number in LOY-2X. PAR RNA cannot be assigned to X or Y.')
  de <- e$de$table; cls <- cn$genes[match(de$gene_id,cn$genes$gene_id),'CN_class']
  de$status <- ifelse(de$thesis_DE & cls=='CN_changed_concordant','Thesis DE, CN-concordant',ifelse(de$thesis_DE,'Thesis DE',ifelse(de$FDR<.05 & cls=='CN_changed_concordant','FDR < 0.05, CN-concordant',ifelse(de$FDR<.05,'FDR < 0.05','Other'))))
  p <- ggplot(de,aes(logFC,-log10(PValue),color=status))+geom_point(data=de[de$status=='Other',],size=.5,alpha=.4)+geom_point(data=de[de$status!='Other',],size=1.4)+
    ggrepel::geom_text_repel(data=de[de$thesis_DE,],aes(label=gene_name),size=2.9,max.overlaps=40,seed=2,show.legend=FALSE)+geom_vline(xintercept=c(-1,1),linetype=2,color='#777777')+
    scale_color_manual(values=c('Thesis DE'='#D55E00','Thesis DE, CN-concordant'='#CC79A7','FDR < 0.05'='#E69F00','FDR < 0.05, CN-concordant'='#56B4E9','Other'='#BBBBBB'))+
    labs(x='log2 fold change (LOY vs sgCTRL, paired voom)',y='-log10 P',color=NULL,title='DU145 paired expression contrast')
  save_figure(p,'Fig07_volcano_paired',9,7,sprintf('Paired voom-limma (~ passage + condition), robust empirical Bayes, BH over %s autosomal/non-PAR X genes. Thesis DE: FDR < 0.05, |log2FC| >= 1 and the same direction in all three passage pairs (%s genes). CN-concordant: |lpWGS LOY-2X minus sgCTRL log2 segment| >= %.1f in the direction of the RNA change.',format(nrow(de),big.mark=','),sum(de$thesis_DE),cfg$cn_abs_delta))
  arm <- cn$arm[!grepl('^chrY',cn$arm$unit),]
  p <- ggplot(arm,aes(DNA_delta_LOY2X,RNA_mean_difference))+geom_hline(yintercept=0,color='#CCCCCC')+geom_vline(xintercept=0,color='#CCCCCC')+geom_smooth(method='lm',se=FALSE,color='#777777',linewidth=.6,formula=y~x)+geom_point(aes(size=genes),color='#0072B2',alpha=.7)+
    ggrepel::geom_text_repel(data=arm[abs(arm$DNA_delta_LOY2X)>0.08 | abs(arm$RNA_mean_difference)>0.06,],aes(label=unit),size=3,seed=4)+
    labs(x='lpWGS arm median log2 difference (LOY-2X minus sgCTRL)',y='RNA arm mean log2 CPM difference',size='Genes',title='Arm-level expression follows copy number')
  cc <- cn$concordance
  save_figure(p,'Fig08_RNA_vs_DNA_arm',8,6.5,sprintf('Autosomal and X chromosome arms with >= 20 mapped genes (CHM13v2.0 arm boundaries). Pearson r = %.2f at arm level, %.2f for 10-Mb windows and %.2f for single genes. The 20p deletion and chrX/15q gains were present in LOY-2X; 20p but not chrX/15q was also seen in DMSO-LOY-1X.',cc$Pearson[cc$level=='arm'],cc$Pearson[cc$level=='window_10Mb'],cc$Pearson[cc$level=='gene']))
  pdiff <- e$pair_diff[de$gene_id[de$thesis_DE],,drop=FALSE]; rownames(pdiff) <- de$gene_name[de$thesis_DE]
  pl <- data.frame(gene=rep(rownames(pdiff),ncol(pdiff)),pair=rep(sub('pair_','',colnames(pdiff)),each=nrow(pdiff)),diff=as.vector(pdiff))
  pl$gene <- factor(pl$gene,levels=rownames(pdiff)[order(rowMeans(pdiff))])
  p <- ggplot(pl,aes(diff,gene,color=pair))+geom_vline(xintercept=0,linetype=2)+geom_point(size=2.2)+scale_color_manual(values=c(P12='#009E73',P13='#CC79A7',P14='#E69F00'))+labs(x='LOY minus sgCTRL log2 CPM within passage',y=NULL,color='Passage',title='Thesis DE genes: passage-matched differences')
  save_figure(p,'SuppFig01_thesis_DE_pair_differences',7,7,'Each point is one passage-matched difference (TMM log2 CPM, prior count 2). All genes shown agree in direction across the three pairs, the DU145 analogue of the HT1080 two-window agreement.')
  eff <- e$effects[e$effects$primary_eligible,]
  p <- ggplot(eff,aes(mean_difference,LOPO_max_abs_change))+geom_point(size=.5,alpha=.3,color='#0072B2')+labs(x='Mean log2 CPM difference',y='Max absolute change after omitting one passage pair',title='Sensitivity to individual passage pairs')
  save_figure(p,'SuppFig02_leave_one_pair_out',7,5.5,'Each omission removes both libraries of one passage block and refits filtering and TMM.')
  cr <- e$comparison
  cl <- rbind(data.frame(other='HT1080',x=cr$mean_difference_HT1080,y=cr$mean_difference_DU145,de=cr$DU145_thesis_DE),data.frame(other='A549',x=cr$mean_difference_A549,y=cr$mean_difference_DU145,de=cr$DU145_thesis_DE))
  cl <- cl[is.finite(cl$x),]
  p <- ggplot(cl,aes(x,y))+geom_point(data=cl[!cl$de,],size=.4,alpha=.25,color='#777777')+geom_point(data=cl[cl$de,],size=1.6,color='#D55E00')+geom_hline(yintercept=0,linetype=2)+geom_vline(xintercept=0,linetype=2)+facet_wrap(~other)+labs(x='Other line: mean log2 CPM difference',y='DU145 mean log2 CPM difference',title='DU145 autosomal effects against HT1080 and A549')
  save_figure(p,'SuppFig03_cross_line_effects',10,5,'Exact shared GENCODE v49 IDs from independently normalised thesis analyses; orange marks DU145 thesis DE genes. Different retained Y inventories, copy-number changes and sampling designs limit mechanistic interpretation.')
  invisible(TRUE)
}
pathway_figures <- function(p) {
  suppressPackageStartupMessages(library(ggplot2)); theme_set(theme_bw(base_size=11))
  h <- p$fgsea[p$fgsea$collection=='H',]; h <- h[order(h$NES),]
  h$pathway <- factor(sub('^HALLMARK_','',h$pathway),levels=sub('^HALLMARK_','',h$pathway))
  h$status <- ifelse(h$interpretable %in% TRUE,'Global BH < 0.05, survives CN screen',ifelse(h$BH_across_collections<.05,'Global BH < 0.05 only','n.s.'))
  g <- ggplot(h,aes(y=pathway))+geom_vline(xintercept=0,linetype=2)+geom_segment(aes(x=NES,xend=NES_CN_stable_only,yend=pathway),color='#BBBBBB')+
    geom_point(aes(x=NES_CN_stable_only),shape=1,size=1.8,color='#555555')+geom_point(aes(x=NES,color=status),size=2.2)+
    scale_color_manual(values=c('Global BH < 0.05, survives CN screen'='#D55E00','Global BH < 0.05 only'='#E69F00','n.s.'='#999999'))+labs(x='NES (filled: all genes; open: CN-changed genes removed)',y=NULL,color=NULL,title='DU145 Hallmark enrichment')
  save_figure(g,'Fig09_Hallmark_enrichment',10,11,'Effect-ranked fgsea on all 50 Hallmark sets; BH is across every eligible MSigDB 2026.1.Hs set. Open circles rerun the ranking without genes whose lpWGS LOY-2X/sgCTRL log2 segment difference is >= 0.2 in magnitude. Surviving the CN screen requires global BH < 0.05, BH < 0.05 in the rerun and the same sign.')
  tf <- p$tf; tf <- head(tf[order(tf$BH_target_pattern,-abs(tf$score)),],25); tf$source <- factor(tf$source,levels=tf$source[order(tf$score)])
  g <- ggplot(tf,aes(score,source,color=score>0))+geom_vline(xintercept=0,linetype=2)+geom_point(size=2.5)+geom_point(aes(x=score_CN_stable_only),shape=1,color='#555555')+scale_color_manual(values=c('FALSE'='#0072B2','TRUE'='#D55E00'),guide='none')+labs(x='ULM score (open: CN-changed genes removed)',y=NULL,title='DU145 CollecTRI target-pattern scores')
  save_figure(g,'Fig10_TF_target_patterns',8,8,'The 25 lowest-BH regulator target patterns (>= 5 represented targets). These are transcript-based target-pattern estimates, not measured regulator activity.')
  invisible(TRUE)
}
