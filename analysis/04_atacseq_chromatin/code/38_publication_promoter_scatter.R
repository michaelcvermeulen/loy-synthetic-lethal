options(warn=1)
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE))))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrepel)})
out<-'integration/feature_scatters/publication';dir.create(out,recursive=TRUE,showWarnings=FALSE)
stem<-'promoter_ATAC010_RNA005_labeled'
d<-fread('integration/feature_scatters/promoter_ATAC010_RNA005_enlarged_all_genes.csv')
stopifnot(nrow(d)==11898,!anyDuplicated(d$entrez_id),sum(d$membership=='Both')==82)
d[,concordance:=fcase(membership!='Both','Not jointly significant',atac_median_logFC>0 & rna_logFC>0,'Up in both',
 atac_median_logFC<0 & rna_logFC<0,'Down in both',default='Opposing directions')]
d[,label_rank_score:=pmin(abs(atac_median_logFC),abs(rna_logFC))]
selected<-rbindlist(lapply(c('Up in both','Down in both'),function(s){
 q<-copy(d[concordance==s]);setorder(q,-label_rank_score,rna_FDR,rna_symbol);head(q,6)
}))
stopifnot(nrow(selected)==12,all(selected$minimum_peak_FDR<.1),all(selected$rna_FDR<.05),
 all(selected$atac_median_logFC*selected$rna_logFC>0))
# Verify each label also has a significant promoter peak in the RNA direction.
links<-fread('integration/feature_overlaps/peak_gene_link_audit.csv')[feature=='promoter' & mapping_status=='eligible']
support<-merge(selected[,.(entrez_id,rna_logFC)],links[,.(entrez_id,peak_id,atac_logFC,atac_FDR)],by='entrez_id',allow.cartesian=TRUE)
support<-support[atac_FDR<.1 & atac_logFC*rna_logFC>0]
stopifnot(all(selected$entrez_id %in% support$entrez_id))
fwrite(support,file.path(out,'labeled_genes_supporting_peaks.csv'))
selected[,label_class:='Concordant overlap']
# Add three extremes per sign for each single-assay membership, ranked on that assay's effect.
outliers<-rbindlist(lapply(c('ATAC only','RNA only'),function(m){
 rbindlist(lapply(c(-1,1),function(sign){
  q<-copy(d[membership==m]);q[,outlier_effect:=if(m=='ATAC only')atac_median_logFC else rna_logFC]
  q<-q[outlier_effect*sign>0];q[,label_rank_score:=abs(outlier_effect)]
  setorder(q,-label_rank_score,rna_symbol);q<-head(q,3)
  q[,label_class:=paste(m,if(sign>0)'positive outlier' else 'negative outlier')];q
 }))
}))
stopifnot(nrow(outliers)==12,all(outliers$membership %in% c('ATAC only','RNA only')),
 !anyDuplicated(outliers$entrez_id),!any(outliers$entrez_id %in% selected$entrez_id))
selected<-rbindlist(list(selected,outliers),fill=TRUE)
d[,label:=ifelse(entrez_id %in% selected$entrez_id,rna_symbol,'')]
fwrite(selected,file.path(out,'labeled_genes.csv'));fwrite(d,file.path(out,'all_plotted_genes.csv'))
counts<-d[,.N,by=membership];quadrants<-d[membership=='Both',.N,by=concordance]
fwrite(quadrants,file.path(out,'overlap_direction_counts.csv'))
pal<-c('Neither'='#C7CDD1','ATAC only'='#0072B2','RNA only'='#D55E00','Both'='#7B3294')
d[,membership:=factor(membership,levels=names(pal))]
label_data<-copy(d[membership=='Both' | entrez_id %in% outliers$entrez_id])
label_data[,label_color:=fcase(membership=='ATAC only','#749BAE',membership=='RNA only','#BF947B',default='#20262B')]
outliers[,label_color:=ifelse(membership=='ATAC only','#749BAE','#BF947B')]
outliers[,point_fill:=ifelse(membership=='ATAC only','#A9C9D9','#E8BEA7')]
xl<-c(floor(min(d$atac_median_logFC)*2)/2-.15,ceiling(max(d$atac_median_logFC)*2)/2+.15)
yl<-c(floor(min(d$rna_logFC))-.25,ceiling(max(d$rna_logFC))+.25)
stopifnot(all(d$atac_median_logFC>xl[1] & d$atac_median_logFC<xl[2]),all(d$rna_logFC>yl[1] & d$rna_logFC<yl[2]))
legend_labels<-setNames(sprintf('%s (%s)',names(pal),format(counts$N[match(names(pal),counts$membership)],big.mark=',',trim=TRUE)),names(pal))
rho<-cor(d$atac_median_logFC,d$rna_logFC,method='spearman')
p<-ggplot(d,aes(atac_median_logFC,rna_logFC))+
 geom_hline(yintercept=0,linetype='dashed',colour='#8E979E',linewidth=.3)+
 geom_vline(xintercept=0,linetype='dashed',colour='#8E979E',linewidth=.3)+
 geom_point(data=d[membership=='Neither'],aes(colour=membership),size=.6,alpha=.18,shape=16)+
 geom_point(data=d[membership %in% c('ATAC only','RNA only')],aes(colour=membership),size=.9,alpha=.5,shape=16)+
 geom_point(data=d[membership=='Both'],aes(colour=membership),size=1.65,alpha=.95,shape=16)+
 geom_point(data=outliers,size=1.9,shape=21,fill=outliers$point_fill,colour=outliers$label_color,stroke=.3,show.legend=FALSE)+
 geom_text_repel(data=label_data,aes(label=label),seed=20260908,
  size=3.0,fontface='italic',family='sans',colour=label_data$label_color,box.padding=.45,point.padding=.2,
  min.segment.length=0,segment.colour=label_data$label_color,segment.size=.25,max.overlaps=Inf,show.legend=FALSE,
  force=4,max.time=10,max.iter=30000,xlim=xl,ylim=yl,
  nudge_x=ifelse(label_data$rna_symbol=='GPR37',.65,ifelse(label_data$atac_median_logFC>0,.25,-.25)),
  nudge_y=ifelse(label_data$rna_symbol=='GPR37',-1.65,ifelse(label_data$rna_logFC>0,.4,-.65)))+
 scale_colour_manual(values=pal,labels=legend_labels,drop=FALSE)+
 guides(colour=guide_legend(override.aes=list(alpha=1,size=2.4),nrow=2,byrow=TRUE))+
 scale_x_continuous(breaks=seq(-3,3,1))+scale_y_continuous(breaks=seq(-6,6,2))+
 coord_cartesian(xlim=xl,ylim=yl,expand=FALSE)+theme_classic(base_size=10,base_family='sans')+
 theme(axis.line=element_line(linewidth=.35,colour='#40484E'),axis.ticks=element_line(linewidth=.3),
  axis.text=element_text(colour='#30373D'),axis.title=element_text(size=11),
  legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=9),
  legend.margin=margin(t=5),legend.key.width=grid::unit(.55,'cm'),
  plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=9.5,lineheight=1.2),
  plot.caption=element_text(size=8,colour='#68747C',hjust=0),plot.margin=margin(12,16,8,10))+
 labs(title='Promoter accessibility and gene expression',
  subtitle=sprintf('ATAC FDR <0.10 | RNA FDR <0.05\n11,898 genes; 82 overlap; all-gene Spearman rho = %.3f',rho),
  x=expression('Promoter ATAC log'[2]*' fold change (LOY / CTRL)'),
  y=expression('RNA log'[2]*' fold change (LOY / CTRL)'),
  caption='Dark labels: concordant overlap. Light blue/orange labels: ATAC-only/RNA-only outliers.')
ggsave(file.path(out,paste0(stem,'.png')),p,width=7.2,height=6.8,dpi=600,bg='white')
ggsave(file.path(out,paste0(stem,'.pdf')),p,width=7.2,height=6.8,device=cairo_pdf,bg='white')
ggsave(file.path(out,paste0(stem,'.svg')),p,width=7.2,height=6.8,device=svglite::svglite,bg='white')
ggsave(file.path(out,paste0(stem,'.tiff')),p,width=7.2,height=6.8,dpi=600,compression='lzw',bg='white')
caption<-paste0('Promoter accessibility and gene expression in HT1080 LOY versus CTRL cells. Each point represents one of 11,898 assay-eligible genes. ',
 'The x-axis shows the median ATAC log2 fold change across all tested promoter-linked peaks (TSS +/-2 kb); the y-axis shows the RNA log2 fold change from the Ben thesis combined-batch DESeq2 model (~ batch + treatment). ',
 'Colors indicate genes with at least one promoter peak at ATAC FDR <0.10, RNA FDR <0.05, both criteria (82 genes), or neither. No fold-change filter was applied. ',
 'Among the 82 overlapping genes, ',quadrants[concordance=='Up in both',N],' have positive effects in both assays, ',quadrants[concordance=='Down in both',N],
 ' have negative effects in both, and ',quadrants[concordance=='Opposing directions',N],' have opposing effects, based on the displayed promoter median and RNA effect. ',
 'Twenty-four genes are labeled. Twelve dark labels show six concordant genes per direction, ranked by min(|median ATAC log2FC|, |RNA log2FC|), with ties resolved by RNA FDR and gene symbol. Each concordant labeled gene also has a significant promoter peak in the RNA direction. ',
 'Twelve lighter labels highlight genes significant in only one assay: three positive and three negative ATAC-only extremes ranked by absolute median ATAC log2FC (light blue), and three positive and three negative RNA-only extremes ranked by absolute RNA log2FC (light orange). Outlier ties are resolved by gene symbol. Their highlighted points and leader lines use matching lighter colors. Outlier is a descriptive label, not a separate statistical test. ',
 'Labels are illustrative selections, not an additional statistical test. All 11,898 coordinates and original membership assignments are retained; axis limits are tightened to the full promoter-specific range without excluding points. ',
 'ATAC significance is peak-level; promoter median effects have no gene-level FDR. Spearman correlation is descriptive and uses all eligible genes. ',
 'ATAC uses three independent replicates per condition; technical batch is unconfirmed and no copy-number adjustment is applied.')
writeLines(caption,file.path(out,'figure_caption.txt'))
writeLines(c('Label rule: six genes per concordant direction, descending minimum absolute fold change across assays; ties: RNA FDR, symbol.',
 'Additional 12 labels: three extremes per sign per single-assay membership, ranked by absolute effect in that assay; ties: symbol. Total 24 labels.',
 'Effect estimates and FDRs are reused exactly from the approved overlap plot; no model refit or shrinkage is performed.',
 'Dimensions: 7.2 x 6.8 inches (182.9 x 172.7 mm). PNG/TIFF 600 dpi; PDF/SVG vector. TIFF LZW compression.',
 'Seed 20260908 for ggrepel. Italic gene labels; original four-category colors retained.',
 'Run code/38_publication_promoter_scatter.R, then code/39_publish_scatter_to_report.R after visual review.'),file.path(out,'METHODS.txt'))
writeLines(c('PASS: original 11,898 unique genes and 82 overlap memberships retained.',
 'PASS: 12 concordant labels pass both FDR thresholds and have a significant promoter peak in the RNA direction; 12 outlier labels are significant in exactly one assay.',
 'PASS: axis limits include every plotted gene; no points removed.'),file.path(out,'verification.txt'))
print(selected[,.(rna_symbol,concordance,label_rank_score)]);print(quadrants)
