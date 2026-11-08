script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE))))
out<-'integration/feature_scatters/publication';stem<-'promoter_ATAC010_RNA005_labeled'
png<-file.path(out,paste0(stem,'.png'))
encoded<-gsub('[\r\n]','',jsonlite::base64_enc(readBin(png,'raw',n=file.info(png)$size)))
selected<-data.table::fread(file.path(out,'labeled_genes.csv'))
svg<-paste(readLines(file.path(out,paste0(stem,'.svg')),warn=FALSE),collapse='\n')
stopifnot(nrow(selected)==24,all(vapply(selected$rna_symbol,function(s)grepl(s,svg,fixed=TRUE),logical(1))))
section<-paste0('<section id="promoter-publication"><h2>Publication figure: promoter accessibility and RNA</h2>',
 '<p>ATAC FDR &lt;0.10 and RNA FDR &lt;0.05: 82 overlapping genes (9 concordant increases, 70 concordant decreases, 3 opposing directions). There are 24 labels: 12 concordant overlap genes in dark text, plus six ATAC-only and six RNA-only outliers in lighter blue/orange. Outliers are three positive and three negative extremes in the corresponding assay.</p>',
 '<figure><img src="data:image/png;base64,',encoded,'" alt="Publication promoter ATAC and RNA scatter with selected concordant gene labels">',
 '<figcaption>All eligible genes and original effects are retained. Lighter labels identify single-assay outliers, not jointly significant genes. Labels are illustrative; the median ATAC effect has no gene-level FDR.</figcaption></figure><p>',
 paste(vapply(c('pdf','svg','png','tiff'),function(ext)paste0('<a href="../',out,'/',stem,'.',ext,'">',toupper(ext),'</a>'),character(1)),collapse=' · '),
 ' · <a href="../',out,'/figure_caption.txt">Figure legend</a> · <a href="../',out,'/labeled_genes.csv">Labeled gene statistics</a></p></section>')
report<-'report/ATACseq_report.html';txt<-paste(readLines(report,warn=FALSE,encoding='UTF-8'),collapse='\n')
if(grepl('<section id="promoter-publication">',txt,fixed=TRUE)) {
 txt<-sub('<section id="promoter-publication">.*?</section>',section,txt,perl=TRUE)
} else {
 txt<-sub('</body>',paste0(section,'\n</body>'),txt,fixed=TRUE)
}
writeLines(txt,report,useBytes=TRUE)
stopifnot(grepl(encoded,paste(readLines(report,warn=FALSE),collapse='\n'),fixed=TRUE))
cat('PASS: all 24 selected gene labels are present in the SVG; current PNG embedded in main report.\n',file=file.path(out,'verification.txt'),append=TRUE)
cat('Publication figure added to main report.\n')
