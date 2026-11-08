options(warn=1)
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE))))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrepel)})
# All-peaks companion to the promoter scatter (code/38): promoter, gene-body and intergenic links pooled.
atac_cut<-.10;rna_cut<-.10
focus<-c('DCC','PDPN','PHF21B','SFN','WNT7A')
out<-'integration/feature_scatters/all_peaks';dir.create(out,recursive=TRUE,showWarnings=FALSE)
stem<-sprintf('all_peaks_ATAC%03d_RNA%03d_labeled',round(atac_cut*100),round(rna_cut*100))
links<-fread('integration/feature_overlaps/peak_gene_link_audit.csv',colClasses=list(character='entrez_id'))[mapping_status=='eligible']
# A peak can link to one gene through more than one feature row; count it once per gene.
links<-unique(links,by=c('entrez_id','peak_id'))
# X = logFC of the gene's best (lowest-FDR) linked peak. The median across pooled gene-body/intergenic peaks
# dilutes the significant peak and misplaces genes by quadrant, so it is kept only as a table column.
d<-links[order(atac_FDR,-abs(atac_logFC)),.(atac_best_logFC=atac_logFC[1],best_peak_id=peak_id[1],best_peak_feature=feature[1],
 atac_median_logFC=median(atac_logFC),minimum_peak_FDR=min(atac_FDR),n_peaks=.N,
 n_sig_peaks=sum(atac_FDR<atac_cut),features=paste(sort(unique(feature)),collapse='+'),
 rna_symbol=rna_symbol[1],rna_logFC=rna_logFC[1],rna_FDR=rna_FDR[1]),by=entrez_id]
stopifnot(!anyDuplicated(d$entrez_id),all(is.finite(d$rna_FDR)),all(is.finite(d$atac_best_logFC)))
d[,atac_sig:=minimum_peak_FDR<atac_cut];d[,rna_sig:=rna_FDR<rna_cut]
d[,membership:=fcase(atac_sig & rna_sig,'Both',atac_sig,'ATAC only',rna_sig,'RNA only',default='Neither')]
d[,concordance:=fcase(membership!='Both','Not jointly significant',atac_best_logFC>0 & rna_logFC>0,'Up in both',
 atac_best_logFC<0 & rna_logFC<0,'Down in both',default='Opposing directions')]
# Direction of the significant peaks themselves, since the median can disagree with them.
sigdir<-links[atac_FDR<atac_cut,.(sig_peak_pattern=fcase(all(atac_logFC>0),'gained',all(atac_logFC<0),'lost',default='mixed')),by=entrez_id]
d<-merge(d,sigdir,by='entrez_id',all.x=TRUE)
d[,label_rank_score:=pmin(abs(atac_best_logFC),abs(rna_logFC))]
selected<-rbindlist(lapply(c('Up in both','Down in both'),function(s){
 q<-copy(d[concordance==s]);setorder(q,-label_rank_score,rna_FDR,rna_symbol);head(q,8)
}))
missing_focus<-setdiff(focus,d$rna_symbol)
if(length(missing_focus)) message('Focus genes not in eligible universe: ',paste(missing_focus,collapse=', '))
label_data<-unique(rbind(selected,d[rna_symbol %in% focus],fill=TRUE),by='entrez_id')
label_data[,label_color:=fcase(rna_symbol %in% focus & membership!='Both','#8A6D3B',default='#20262B')]
label_data[,label:=ifelse(rna_symbol %in% focus & membership!='Both',paste0(rna_symbol,'*'),rna_symbol)]
counts<-d[,.N,by=membership];quadrants<-d[membership=='Both',.N,by=concordance]
fwrite(d,file.path(out,paste0(stem,'_all_genes.csv')))
fwrite(d[membership=='Both'][order(concordance,rna_FDR)],file.path(out,paste0(stem,'_overlap_genes.csv')))
fwrite(label_data,file.path(out,paste0(stem,'_labeled_genes.csv')))
fwrite(links[entrez_id %in% d[rna_symbol %in% focus,entrez_id]][order(rna_symbol,atac_FDR)],
 file.path(out,paste0(stem,'_focus_gene_peaks.csv')))
pal<-c('Neither'='#C7CDD1','ATAC only'='#0072B2','RNA only'='#D55E00','Both'='#7B3294')
d[,membership:=factor(membership,levels=names(pal))]
xl<-c(floor(min(d$atac_best_logFC)*2)/2-.15,ceiling(max(d$atac_best_logFC)*2)/2+.15)
yl<-c(floor(min(d$rna_logFC))-.25,ceiling(max(d$rna_logFC))+.25)
legend_labels<-setNames(sprintf('%s (%s)',names(pal),format(counts$N[match(names(pal),counts$membership)],big.mark=',',trim=TRUE)),names(pal))
rho<-cor(d$atac_best_logFC,d$rna_logFC,method='spearman')
nq<-function(s){n<-quadrants[concordance==s,N];if(length(n))n else 0L}
p<-ggplot(d,aes(atac_best_logFC,rna_logFC))+
 geom_hline(yintercept=0,linetype='dashed',colour='#8E979E',linewidth=.3)+
 geom_vline(xintercept=0,linetype='dashed',colour='#8E979E',linewidth=.3)+
 geom_point(data=d[membership=='Neither'],aes(colour=membership),size=.6,alpha=.18,shape=16)+
 geom_point(data=d[membership %in% c('ATAC only','RNA only')],aes(colour=membership),size=.9,alpha=.5,shape=16)+
 geom_point(data=d[membership=='Both'],aes(colour=membership),size=1.65,alpha=.95,shape=16)+
 geom_point(data=label_data[rna_symbol %in% focus],shape=21,size=2.6,fill=NA,colour='#20262B',stroke=.5)+
 geom_text_repel(data=label_data,aes(label=label),seed=20260922,
  size=3.0,fontface='italic',family='sans',colour=label_data$label_color,box.padding=.45,point.padding=.2,
  min.segment.length=0,segment.colour=label_data$label_color,segment.size=.25,max.overlaps=Inf,
  force=4,max.time=10,max.iter=30000,xlim=xl,ylim=yl,show.legend=FALSE)+
 scale_colour_manual(values=pal,labels=legend_labels,drop=FALSE)+
 guides(colour=guide_legend(override.aes=list(alpha=1,size=2.4),nrow=2,byrow=TRUE))+
 coord_cartesian(xlim=xl,ylim=yl,expand=FALSE)+theme_classic(base_size=10,base_family='sans')+
 theme(axis.line=element_line(linewidth=.35,colour='#40484E'),axis.ticks=element_line(linewidth=.3),
  axis.text=element_text(colour='#30373D'),axis.title=element_text(size=11),
  legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=9),
  plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=9.5,lineheight=1.2),
  plot.caption=element_text(size=8,colour='#68747C',hjust=0),plot.margin=margin(12,16,8,10))+
 labs(title='Accessibility (all linked peaks) and gene expression',
  subtitle=sprintf('ATAC FDR <%.2f (any linked peak) | RNA FDR <%.2f\n%s genes; %d overlap (%d up/up, %d down/down, %d opposing); Spearman rho = %.3f',
   atac_cut,rna_cut,format(nrow(d),big.mark=','),sum(d$membership=='Both'),nq('Up in both'),nq('Down in both'),nq('Opposing directions'),rho),
  x=expression('ATAC log'[2]*' fold change, best linked peak (LOY / CTRL)'),
  y=expression('RNA log'[2]*' fold change (LOY / CTRL)'),
  caption=paste0('Peaks linked by promoter (TSS +/-2 kb), gene-body overlap or nearest TSS (intergenic).\n',
   'X = lowest-FDR linked peak per gene. Dark labels: top 8 concordant genes per direction.\n',
   'Circled: DCC, PDPN, PHF21B, SFN, WNT7A; * = not jointly significant. RNA: batch + treatment DESeq2.'))
ggsave(file.path(out,paste0(stem,'.png')),p,width=7.2,height=7.3,dpi=300,bg='white')
ggsave(file.path(out,paste0(stem,'.pdf')),p,width=7.2,height=7.3,device=cairo_pdf,bg='white')
writeLines(c(sprintf('All-peaks companion to code/38. Eligible links from integration/feature_overlaps/peak_gene_link_audit.csv pooled across promoter, gene_body and intergenic features; each (gene, peak) pair counted once.'),
 'One point per Entrez gene. X: log2FC of the lowest-FDR linked peak (ties: larger |logFC|). Y: batch+treatment DESeq2 RNA log2FC (unchanged).',
 sprintf('ATAC membership: any linked peak FDR <%.2f. RNA membership: gene FDR <%.2f. No fold-change filter.',atac_cut,rna_cut),
 'Median across all linked peaks is kept as a column only: pooling gene-body/intergenic peaks diluted the significant peak and misassigned quadrants.',
 'Concordance uses the best-peak direction; sig_peak_pattern records whether all significant peaks agree (gained/lost/mixed).',
 'Intergenic nearest-TSS and gene-body links are positional hypotheses. Peak FDRs are not gene-level tests; overlap has no calibrated joint FDR.',
 'Labels: top 8 concordant genes per direction by min(|ATAC best-peak logFC|, |RNA|), ties by RNA FDR and symbol; plus the five focus genes. Seed 20260922.'),
 file.path(out,paste0(stem,'_METHODS.txt')))
print(counts);print(quadrants)
print(d[rna_symbol %in% focus,.(rna_symbol,features,n_peaks,n_sig_peaks,atac_best_logFC=round(atac_best_logFC,2),best_peak_feature,minimum_peak_FDR=signif(minimum_peak_FDR,2),sig_peak_pattern,rna_logFC=round(rna_logFC,2),rna_FDR=signif(rna_FDR,2),membership,concordance)])
print(selected[,.(rna_symbol,concordance,features)])
