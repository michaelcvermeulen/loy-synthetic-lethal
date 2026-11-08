script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE))))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out<-'integration/down70_ATAC010_RNA005'
rd<-function(n)fread(file.path(out,paste0(n,'.csv')))
go<-rd('GO_all_tested_terms')
# Also expose the commonly used correction over returned, nonzero-hit terms.
go[,`:=`(BH_FDR_hit_terms_only=NA_real_,conditioned_BH_FDR_hit_terms_only=NA_real_)]
go[query_hits>0,`:=`(BH_FDR_hit_terms_only=p.adjust(PValue,'BH'),conditioned_BH_FDR_hit_terms_only=p.adjust(chromosome_conditioned_PValue,'BH')),by=ontology]
fwrite(go,file.path(out,'GO_all_tested_terms.csv'))
chr<-rd('chromosome_enrichment');fam<-rd('ZNF_ZFP_named_enrichment');cl<-rd('within_chromosome_global_clustering')
q<-rd('exact_70_genes');bg<-rd('assay_eligible_background');acc<-rd('GO_accounting')
blocks<-rd('descriptive_1Mb_gap_blocks');sharing<-rd('shared_significant_promoter_peaks')
links<-rd('significant_down_promoter_links')
# Refresh the chromosome plot with explicit color/legend ordering.
cp<-copy(chr);cp[,chrom:=factor(chrom,levels=intersect(paste0('chr',c(1:22,'X','Y')),chr$chrom))]
cp[,query_genes:=as.numeric(query_genes)]
cm<-melt(cp,id.vars='chrom',measure.vars=c('query_genes','expected'))
p<-ggplot(cm,aes(chrom,value,fill=variable))+geom_col(position='dodge')+
 scale_fill_manual(values=c(query_genes='#7B3294',expected='#BFC8CD'),breaks=c('query_genes','expected'),labels=c('Observed down/down','Expected from eligible genes'))+
 labs(title='Chromosome distribution of the exact 70 down/down genes',subtitle='Expected counts use the 11,898 promoter assay-eligible genes',x=NULL,y='Number of genes',fill=NULL)+
 theme_classic(base_size=11)+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='bottom')
for(ext in c('png','pdf'))ggsave(file.path(out,paste0('chromosome_counts.',ext)),p,width=10,height=5,dpi=240,bg='white')
acc[,`:=`(hit_only_FDR005=vapply(ontology,function(o)sum(go[ontology==o,BH_FDR_hit_terms_only]<.05,na.rm=TRUE),integer(1)),
 minimum_FDR_hit_terms_only=vapply(ontology,function(o)min(go[ontology==o,BH_FDR_hit_terms_only],na.rm=TRUE),numeric(1)))]
fwrite(acc,file.path(out,'GO_accounting.csv'))
esc<-function(s){s<-gsub('&','&amp;',as.character(s),fixed=TRUE);s<-gsub('<','&lt;',s,fixed=TRUE);gsub('>','&gt;',s,fixed=TRUE)}
tab<-function(d){
 for(k in names(d))if(is.numeric(d[[k]]))set(d,j=k,value=format(d[[k]],digits=4,trim=TRUE))
 paste0('<table><thead><tr>',paste0('<th>',esc(names(d)),'</th>',collapse=''),'</tr></thead><tbody>',
 paste(apply(as.matrix(d),1,function(r)paste0('<tr>',paste0('<td>',esc(r),'</td>',collapse=''),'</tr>')),collapse=''),'</tbody></table>')
}
fig<-function(n,caption){p<-file.path(out,paste0(n,'.png'));b<-gsub('[\r\n]','',jsonlite::base64_enc(readBin(p,'raw',n=file.info(p)$size)))
 paste0('<figure><img src="data:image/png;base64,',b,'" alt="',esc(caption),'"><figcaption>',esc(caption),'</figcaption></figure>')}
top<-head(go[ontology=='BP' & query_hits>0][order(PValue)],15)
statement<-paste0('The exact 70 concordant down/down genes show a significant ZNF/ZFP-named gene component and chromosome-19 concentration, but no GO term passes FDR 0.05. ',
 'Chromosome 19 contains 19/70 genes versus 4.84 expected (3.92-fold; BH FDR 3.72e-6). ',
 'ZNF/ZFP-prefixed symbols comprise 18/70 genes versus 381/11,898 eligible genes (8.03-fold; P=4.25e-12); the targeted chromosome-conditioned P is 3.89e-8. ',
 'Genes cluster at <=1 Mb beyond chromosome composition (24 pairs versus 9.25 expected; permutation P=0.00110, BH FDR=0.00330). ',
 'The chromosome-19-specific 1 Mb test has P=0.00090 but FDR=0.0621 after testing all 23 chromosomes at three distances. ',
 'These findings support a family and positional pattern, not a validated cell-state pathway, copy-number explanation or causal LOY mechanism.')
writeLines(statement,file.path(out,'CONCLUSIONS.md'))
html<-c('<!doctype html><html><head><meta charset="utf-8"><title>Exact 70-gene enrichment and chromosome clustering</title>',
 '<style>body{font:16px/1.55 system-ui,sans-serif;max-width:1150px;margin:40px auto;padding:0 24px;color:#18242b}h2{margin-top:36px}table{border-collapse:collapse;display:block;overflow:auto;font-size:12px}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:left}th{background:#edf3f5}img{max-width:100%}.lead{background:#edf3f5;padding:18px}figcaption{font-size:13px;color:#46545c}</style></head><body>',
 '<h1>Exact 70 down/down genes: enrichment and chromosome clustering</h1>',
 '<p>ATAC FDR &lt;0.10; RNA FDR &lt;0.05; negative displayed effects in both assays. Background: the same 11,898 eligible promoter genes used in the scatter plot.</p>',
 paste0('<p class="lead">',statement,'</p>'),
 '<h2>Functional enrichment</h2>',
 '<p>GO biological process is the primary functional analysis; molecular function and cellular component are secondary. Tests use each ontology-annotated subset of the assay-eligible background. The primary BH correction includes every GO term with 10–500 background genes, including terms with zero query hits. A chromosome-conditioned test fixes the observed query count on each chromosome and sums the corresponding hypergeometric distributions.</p>',
 tab(copy(acc)),
 '<p><b>No significant terms were found before or after chromosome conditioning.</b> For comparison with workflows that correct only returned terms, a separate BH correction over terms with at least one query hit is supplied. It also yields no significant terms. These correction families differ from the earlier 76-gene, two-direction analysis; their FDR values should not be compared directly.</p>',
 fig('GO_BP_top_terms','Top nominal GO biological-process results. Displayed themes are exploratory; none pass BH FDR 0.05.'),
 tab(top[,.(Description,query_hits,fold_enrichment,PValue,BH_FDR,BH_FDR_hit_terms_only,chromosome_conditioned_PValue,chromosome_conditioned_BH_FDR,query_symbols)]),
 '<h2>Targeted ZNF/ZFP-name enrichment</h2>',tab(copy(fam)),
 '<p>This targeted follow-up tests a category identified in the earlier inspection: symbols beginning ZNF or ZFP. It is not an exhaustive annotation of zinc-finger proteins, nor a genome-wide corrected family discovery. Its enrichment persists after conditioning on chromosome composition, so the chromosome-19 excess alone does not explain it.</p>',
 '<h2>Chromosome enrichment</h2>',fig('chromosome_counts','Observed gene counts versus expected counts under sampling from all assay-eligible genes.'),tab(copy(chr)),
 '<p>Only chromosome 19 passes chromosome-level BH FDR 0.05. Enrichment refers to the down/down gene set relative to eligible genes, not to DNA copy number.</p>',
 '<h2>Clustering within chromosomes</h2>',
 '<p>Ten thousand random gene sets preserve the observed number of genes on every chromosome and sample without replacement from eligible genes. Pair distances use the genomic starts recorded in the selected RNA source. This tests local clustering beyond the chromosome-count imbalance.</p>',tab(copy(cl)),
 fig('gene_positions','Positions of query genes, with ZNF/ZFP-prefixed symbols highlighted. All query gene positions are present.'),
 '<p>The global 1 Mb test is significant, while 5 Mb and 10 Mb tests are not. Of 24 query pairs within 1 Mb, 23 are on chromosome 19 and one is on chromosome 5. Individual chromosome tests undergo a separate BH correction across 23 chromosomes and three distances; none pass FDR 0.05, with chromosome 19 at 1 Mb closest (FDR 0.0621).</p>',
 '<h3>Descriptive blocks joined by gaps of at most 1 Mb</h3>',tab(copy(blocks)),
 '<p>These blocks are descriptive, not independently tested loci. A chain of gaps below 1 Mb can span more than 1 Mb.</p>',
 '<h2>Shared promoter evidence and interpretation limits</h2>',tab(copy(sharing)),
 paste0('<p>The 70 genes have ',uniqueN(links$peak_id),' distinct significant decreased promoter peaks, and ',uniqueN(links$entrez_id),' genes have such support. Three promoter peaks are linked to multiple selected genes. In particular, ZNF345/ZNF790/ZNF568 share one peak and ZNF28/ZNF468 share another. Thus the 70 genes do not represent 70 independent ATAC events.</p>'),
 '<p>Chromosome conditioning controls chromosome composition, but not local dependence between neighboring genes, shared promoter peaks or differences in promoter-peak number. Unconfirmed ATAC batch structure and possible copy-number effects remain unresolved. lpWGS integration is needed to assess copy-number explanations; these data alone do not establish one.</p>',
 '<h2>Saved data and methods</h2><p><a href="../integration/down70_ATAC010_RNA005/exact_70_genes.csv">Exact gene list</a> · <a href="../integration/down70_ATAC010_RNA005/GO_all_tested_terms.csv">All GO tests</a> · <a href="../integration/down70_ATAC010_RNA005/chromosome_enrichment.csv">Chromosome tests</a> · <a href="../integration/down70_ATAC010_RNA005/METHODS.txt">Methods</a> · <a href="ATACseq_report.html">Main report</a></p>',
 '</body></html>')
report<-'report/down70_enrichment.html';writeLines(html,report,useBytes=TRUE)
snippet<-paste0('<section id="down70-enrichment"><h2>Enrichment of the 70 concordant down/down genes</h2><p>No GO term passes FDR 0.05. Chromosome 19 is enriched (19 genes; 3.92-fold; FDR 3.72e-6), and the ZNF/ZFP-named component remains enriched after chromosome conditioning. Local 1 Mb clustering is significant globally (FDR 0.00330), with shared promoter peaks limiting independence. <a href="down70_enrichment.html">Full functional and chromosome analysis</a>.</p></section>')
main<-'report/ATACseq_report.html';text<-paste(readLines(main,warn=FALSE,encoding='UTF-8'),collapse='\n')
if(grepl('<section id="down70-enrichment">',text,fixed=TRUE)){text<-sub('<section id="down70-enrichment">.*?</section>',snippet,text,perl=TRUE)}else{text<-sub('</body>',paste0(snippet,'\n</body>'),text,fixed=TRUE)}
writeLines(text,main,useBytes=TRUE)
stopifnot(all(acc$hit_only_FDR005==0),all(acc$FDR005==0),all(acc$conditioned_FDR005==0),
 nrow(q)==70,nrow(bg)==11898,identical(sort(q$entrez_id),sort(bg[query==TRUE,entrez_id])))
# Independent enumeration of chromosome-count P-values using Fisher's exact test.
for(i in seq_len(nrow(chr))){k<-chr$query_genes[i];M<-chr$background_genes[i]
 p<-fisher.test(matrix(c(k,M-k,70-k,11898-M-70+k),nrow=2),alternative='greater')$p.value
 stopifnot(isTRUE(all.equal(p,chr$PValue[i],tolerance=1e-10)))
}
cat('PASS: chromosome P-values independently match Fisher exact tests; hit-term-only GO correction also yields no significant terms.\n',file=file.path(out,'verification.txt'),append=TRUE)
cat('Sensitivity: GO_all_tested_terms.csv also reports BH over positive-hit terms only within each ontology; no significant terms.\n',file=file.path(out,'METHODS.txt'),append=TRUE)
jsonlite::write_json(list(status='complete',completed=as.character(Sys.time()),query_genes=70,background_genes=11898,
 report=report,GO_significant=0,permutations=10000,validation='passed'),file.path(out,'ANALYSIS_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
print(acc);cat('Saved',report,'\n')
