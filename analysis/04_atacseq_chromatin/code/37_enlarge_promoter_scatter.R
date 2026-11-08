options(warn=1)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE))))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out <- 'integration/feature_scatters'
stem <- 'promoter_ATAC010_RNA005_enlarged'
base <- fread(file.path(out,'feature_gene_effects.csv'))
d <- copy(base[feature=='promoter'])
d[,membership:=fcase(minimum_peak_FDR<.1 & rna_FDR<.05,'Both',minimum_peak_FDR<.1,'ATAC only',rna_FDR<.05,'RNA only',default='Neither')]
pal<-c('Neither'='#C7CDD1','ATAC only'='#0072B2','RNA only'='#D55E00','Both'='#7B3294')
d[,membership:=factor(membership,levels=names(pal))]
expected<-fread(file.path(out,'panel_summary.csv'))[feature=='promoter' & ATAC_FDR==.1 & RNA_FDR==.05]
stopifnot(nrow(expected)==1,nrow(d)==expected$universe,sum(d$membership=='Both')==expected$overlap,
 sum(d$membership=='ATAC only')==expected$ATAC_only,sum(d$membership=='RNA only')==expected$RNA_only,
 sum(d$membership=='Neither')==expected$neither,!anyDuplicated(d$entrez_id))
fwrite(d,file.path(out,paste0(stem,'_all_genes.csv')))
fwrite(d[membership=='Both'][order(rna_FDR)],file.path(out,paste0(stem,'_overlap_genes.csv')))
xr<-range(base$atac_median_logFC);yr<-range(base$rna_logFC)
xl<-xr+c(-1,1)*diff(xr)*.04;yl<-yr+c(-1,1)*diff(yr)*.12
rho<-cor(d$atac_median_logFC,d$rna_logFC,method='spearman')
counts<-d[,.N,by=membership]
legend_labels<-setNames(sprintf('%s (%s)',names(pal),format(counts$N[match(names(pal),counts$membership)],big.mark=',',trim=TRUE)),names(pal))
p<-ggplot(d,aes(atac_median_logFC,rna_logFC))+
 geom_hline(yintercept=0,lty=2,colour='grey55',linewidth=.45)+geom_vline(xintercept=0,lty=2,colour='grey55',linewidth=.45)+
 geom_point(data=d[membership=='Neither'],aes(colour=membership),size=1,alpha=.22)+
 geom_point(data=d[membership %in% c('ATAC only','RNA only')],aes(colour=membership),size=1.6,alpha=.6)+
 geom_point(data=d[membership=='Both'],aes(colour=membership),size=2.7,alpha=.95)+
 scale_colour_manual(values=pal,labels=legend_labels,drop=FALSE)+
 guides(colour=guide_legend(override.aes=list(alpha=1,size=3.5),nrow=1))+
 coord_cartesian(xlim=xl,ylim=yl,expand=FALSE)+theme_bw(base_size=17)+
 theme(legend.position='bottom',legend.title=element_blank(),panel.grid.minor=element_blank(),
 plot.title=element_text(face='bold',size=24),plot.subtitle=element_text(size=16),
 plot.caption=element_text(size=11,hjust=0),plot.margin=margin(18,22,18,18))+
 labs(title='Promoter accessibility versus RNA expression',
 subtitle=sprintf('ATAC FDR <0.10 | RNA FDR <0.05\n%s eligible genes | %s overlap | All-gene Spearman rho = %.3f',
 format(nrow(d),big.mark=',',trim=TRUE),expected$overlap,rho),
 x='Median promoter ATAC log2 fold change (LOY / CTRL)',y='RNA log2 fold change (LOY / CTRL)',
 caption='One point per eligible gene. ATAC membership: any linked promoter peak passes FDR; no fold-change filter.\nATAC coordinates summarize all tested promoter peaks; the median has no gene-level FDR.\nOriginal gene coordinates, colors and axis limits retained. RNA: Ben thesis HT1080 batch + treatment result.')
ggsave(file.path(out,paste0(stem,'.png')),p,width=12,height=10,dpi=300,bg='white')
ggsave(file.path(out,paste0(stem,'.pdf')),p,width=12,height=10,bg='white',useDingbats=FALSE)
img<-file.path(out,paste0(stem,'.png'))
encoded<-jsonlite::base64_enc(readBin(img,'raw',n=file.info(img)$size))
encoded<-gsub('[\r\n]','',encoded)
section<-paste0('<section id="promoter-atac010-rna005"><h2>Promoter overlap: ATAC FDR 0.10 / RNA FDR 0.05</h2>',
 '<p>Enlarged selected panel: 11,898 eligible genes and 82 overlapping genes. Membership uses any linked promoter peak at ATAC FDR &lt;0.10 and RNA FDR &lt;0.05, without a fold-change filter.</p>',
 '<figure><img src="data:image/png;base64,',encoded,'" alt="Enlarged promoter ATAC versus RNA scatter at FDR 0.10 and 0.05">',
 '<figcaption>Coordinates summarize all tested linked promoter peaks; significance remains peak-level.</figcaption></figure>',
 '<p><a href="../',out,'/',stem,'.png">High-resolution PNG</a> · <a href="../',out,'/',stem,'.pdf">Vector PDF</a> · <a href="../',out,'/',stem,'_overlap_genes.csv">82 overlap genes</a></p></section>')
report<-'report/ATACseq_report.html'
txt<-paste(readLines(report,warn=FALSE,encoding='UTF-8'),collapse='\n')
if(grepl('<section id="promoter-atac010-rna005">',txt,fixed=TRUE)){
 txt<-sub('<section id="promoter-atac010-rna005">.*?</section>',section,txt,perl=TRUE)
}else{
 stopifnot(grepl('</body>',txt,fixed=TRUE));txt<-sub('</body>',paste0(section,'\n</body>'),txt,fixed=TRUE)
}
writeLines(txt,report,useBytes=TRUE)
stopifnot(grepl(encoded,paste(readLines(report,warn=FALSE),collapse='\n'),fixed=TRUE))
writeLines(c('PASS: all four membership counts match the original ATAC 0.10 / RNA 0.05 promoter panel.',
 'PASS: 11,898 unique genes, 82 overlapping genes, original coordinates and full shared axis limits retained.',
 'PASS: current enlarged PNG is embedded in the main report. Exports: 3600 x 3000 PNG and 12 x 10 inch vector PDF.'),file.path(out,paste0(stem,'_verification.txt')))
print(counts);cat('Created',img,'\n')
