options(warn=1)
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE))))
suppressPackageStartupMessages({library(data.table);library(AnnotationDbi);library(org.Hs.eg.db);library(GO.db);library(ggplot2);library(jsonlite)})
out<-'integration/down70_ATAC010_RNA005';dir.create(out,recursive=TRUE,showWarnings=FALSE)
put<-function(x,n)fwrite(x,file.path(out,paste0(n,'.csv')))
d<-fread('integration/feature_scatters/publication/all_plotted_genes.csv',colClasses=c(entrez_id='character'))
stopifnot(nrow(d)==11898,!anyDuplicated(d$entrez_id))
d[,query:=minimum_peak_FDR<.1 & rna_FDR<.05 & atac_median_logFC<0 & rna_logFC<0]
stopifnot(sum(d$query)==70,all(d$query==(d$membership=='Both' & d$concordance=='Down in both')))
rna<-fread('config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv')
stopifnot(!anyDuplicated(rna$gene_id),all(d$rna_ensembl_id %in% rna$gene_id))
ix<-match(d$rna_ensembl_id,rna$gene_id)
d[,`:=`(chrom=paste0('chr',sub('^chr','',rna$chr[ix])),gene_start=rna$gene_start[ix],gene_end=rna$gene_end[ix])]
stopifnot(!anyNA(d$chrom),all(is.finite(d$gene_start)),all(d$gene_end>=d$gene_start))
d[,ZNF_ZFP_named:=grepl('^(ZNF|ZFP)',rna_symbol)]
put(d,'assay_eligible_background');put(d[query==TRUE],'exact_70_genes')
put(as.data.table(AnnotationDbi::metadata(org.Hs.eg.db)),'annotation_provenance')
# Exact upper tail of the sum of chromosome-specific hypergeometric variables.
# Conditions on the observed number of query genes on every chromosome.
conditional_tail<-function(K,N,n,k){
 if(k==0)return(1)
 prob<-1
 for(i in seq_along(N)){
  if(n[i]==0)next
  lo<-max(0,n[i]-(N[i]-K[i]));hi<-min(n[i],K[i]);x<-lo:hi
  pr<-dhyper(x,K[i],N[i]-K[i],n[i]);next_prob<-numeric(length(prob)+hi)
  for(j in seq_along(x))next_prob[seq_along(prob)+x[j]]<-next_prob[seq_along(prob)+x[j]]+prob*pr[j]
  prob<-next_prob
 }
 min(1,max(0,if(k>=length(prob))0 else sum(prob[(k+1):length(prob)])))
}
chr<-d[,.(background_genes=.N,query_genes=sum(query),query_symbols=paste(rna_symbol[query],collapse=';')),by=chrom]
chr[,`:=`(expected=70*background_genes/nrow(d),fold_enrichment=query_genes/(70*background_genes/nrow(d)),
 PValue=phyper(query_genes-1,background_genes,nrow(d)-background_genes,70,lower.tail=FALSE))]
chr[,BH_FDR:=p.adjust(PValue,'BH')];put(chr[order(PValue)],'chromosome_enrichment')
N<-chr$background_genes;n<-chr$query_genes
K<-sapply(chr$chrom,function(ch)sum(d$chrom==ch & d$ZNF_ZFP_named))
k<-sum(d$query & d$ZNF_ZFP_named)
znf<-data.table(category='Symbols beginning ZNF or ZFP (not all zinc-finger proteins)',query_count=k,query_total=70,
 background_count=sum(K),background_total=nrow(d),fold_enrichment=(k/70)/(sum(K)/nrow(d)),
 PValue=phyper(k-1,sum(K),nrow(d)-sum(K),70,lower.tail=FALSE),
 chromosome_conditioned_PValue=conditional_tail(K,N,n,k))
put(znf,'ZNF_ZFP_named_enrichment')
cat('Chromosome and named-family tests complete\n')
# GO annotations include ancestor terms. All eligible terms, including zero-hit terms, enter BH.
ann<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=d$entrez_id,keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL')))
ann<-unique(ann[!is.na(GOALL),.(entrez_id=ENTREZID,GO_ID=GOALL,ontology=ONTOLOGYALL)])
terms<-as.data.table(AnnotationDbi::select(GO.db,keys=unique(ann$GO_ID),keytype='GOID',columns='TERM'))
put(ann,'GO_gene_term_membership')
allgo<-list();account<-list()
for(ont in c('BP','MF','CC')){
 a<-ann[ontology==ont];u<-d[entrez_id %in% a$entrez_id];q<-u[query==TRUE]
 sets<-split(a$entrez_id,a$GO_ID);sets<-sets[lengths(sets)>=10 & lengths(sets)<=500]
 strata<-u[,.(N=.N,n=sum(query)),by=chrom]
 result<-rbindlist(lapply(names(sets),function(id){
  member<-sets[[id]];hits<-q[entrez_id %in% member]
  K<-tabulate(match(u[entrez_id %in% member,chrom],strata$chrom),nbins=nrow(strata))
  k<-nrow(hits);M<-length(member)
  data.table(GO_ID=id,ontology=ont,Description=terms$TERM[match(id,terms$GOID)],query_hits=k,query_annotated=nrow(q),
   background_hits=M,background_annotated=nrow(u),fold_enrichment=(k/nrow(q))/(M/nrow(u)),
   PValue=phyper(k-1,M,nrow(u)-M,nrow(q),lower.tail=FALSE),
   chromosome_conditioned_expected=sum(strata$n*K/strata$N),
   chromosome_conditioned_PValue=conditional_tail(K,strata$N,strata$n,k),
   query_symbols=paste(hits$rna_symbol,collapse=';'),query_entrez=paste(hits$entrez_id,collapse=';'))
 }))
 result[,`:=`(BH_FDR=p.adjust(PValue,'BH'),chromosome_conditioned_BH_FDR=p.adjust(chromosome_conditioned_PValue,'BH'))]
 allgo[[ont]]<-result
 account[[ont]]<-data.table(ontology=ont,query_total=70,query_annotated=nrow(q),background_total=nrow(d),background_annotated=nrow(u),
  tested_terms=nrow(result),terms_with_query_hits=sum(result$query_hits>0),FDR005=sum(result$BH_FDR<.05),conditioned_FDR005=sum(result$chromosome_conditioned_BH_FDR<.05))
 cat(ont,nrow(result),'terms; significant',sum(result$BH_FDR<.05),'conditioned',sum(result$chromosome_conditioned_BH_FDR<.05),'\n')
}
go<-rbindlist(allgo);go[,`:=`(BH_FDR_all_ontologies=p.adjust(PValue,'BH'),conditioned_BH_FDR_all_ontologies=p.adjust(chromosome_conditioned_PValue,'BH'))]
put(go[order(ontology,PValue)],'GO_all_tested_terms');put(go[BH_FDR<.05 | chromosome_conditioned_BH_FDR<.05][order(ontology,PValue)],'GO_significant_terms')
put(rbindlist(account),'GO_accounting')
# Chromosome-preserving permutations test whether genes cluster within chromosomes.
set.seed(20260909);B<-10000L;radii<-c(1e6,5e6,1e7)
pair_counts<-function(pos){if(length(pos)<2)return(rep(0L,3));distances<-as.numeric(dist(pos));vapply(radii,function(r)sum(distances<=r),integer(1))}
observed<-sapply(chr$chrom,function(ch)pair_counts(d[query==TRUE & chrom==ch,gene_start]))
null<-array(0L,dim=c(B,3,nrow(chr)))
for(j in seq_len(nrow(chr))){
 ch<-chr$chrom[j];pool<-d[chrom==ch,gene_start];nj<-chr$query_genes[j]
 if(nj>=2)for(b in seq_len(B))null[b,,j]<-pair_counts(pool[sample.int(length(pool),nj)])
}
global_obs<-rowSums(observed);global_null<-apply(null,c(1,2),sum)
cluster<-data.table(distance_bp=radii,observed_gene_pairs=global_obs,expected_pairs=colMeans(global_null),
 permutation_PValue=(1+colSums(sweep(global_null,2,global_obs,'>=')>0))/(B+1),permutations=B)
cluster[,BH_FDR:=p.adjust(permutation_PValue,'BH')];put(cluster,'within_chromosome_global_clustering')
chrcluster<-rbindlist(lapply(seq_len(nrow(chr)),function(j)data.table(chrom=chr$chrom[j],query_genes=chr$query_genes[j],distance_bp=radii,
 observed_pairs=observed[,j],expected_pairs=colMeans(null[,,j]),PValue=(1+colSums(sweep(null[,,j],2,observed[,j],'>=')>0))/(B+1))))
chrcluster[,BH_FDR_all_chromosomes_distances:=p.adjust(PValue,'BH')];put(chrcluster[order(PValue)],'within_chromosome_local_clustering')
positions<-copy(d[query==TRUE][order(chrom,gene_start)])
positions[,block:=cumsum(c(TRUE,diff(gene_start)>1e6)),by=chrom]
blocks<-positions[,.(genes=.N,start=min(gene_start),end=max(gene_end),symbols=paste(rna_symbol,collapse=';')),by=.(chrom,block)]
put(blocks[genes>=2][order(-genes)],'descriptive_1Mb_gap_blocks')
links<-fread('integration/feature_overlaps/peak_gene_link_audit.csv',colClasses=c(entrez_id='character'))[feature=='promoter' & mapping_status=='eligible' & entrez_id %in% d[query==TRUE,entrez_id] & atac_FDR<.1 & atac_logFC<0]
put(links,'significant_down_promoter_links')
sharing<-links[,.(genes=uniqueN(entrez_id),symbols=paste(unique(d$rna_symbol[match(entrez_id,d$entrez_id)]),collapse=';')),by=peak_id]
put(sharing[genes>1][order(-genes)],'shared_significant_promoter_peaks')
theme_set(theme_classic(base_size=11))
savep<-function(p,name,w=9,h=6){for(ext in c('png','pdf'))ggsave(file.path(out,paste0(name,'.',ext)),p,width=w,height=h,dpi=240,bg='white')}
chr_order<-intersect(paste0('chr',c(1:22,'X','Y')),chr$chrom)
cp<-copy(chr);cp[,chrom:=factor(chrom,levels=chr_order)]
cm<-melt(cp,id.vars='chrom',measure.vars=c('query_genes','expected'))
p<-ggplot(cm,aes(chrom,value,fill=variable))+geom_col(position='dodge')+scale_fill_manual(values=c(query_genes='#7B3294',expected='#BFC8CD'),breaks=c('query_genes','expected'),labels=c('Observed down/down','Expected from eligible genes'))+
 labs(title='Chromosome distribution of the exact 70 down/down genes',subtitle='Expected counts use the 11,898 promoter assay-eligible genes',x=NULL,y='Number of genes',fill=NULL)+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='bottom')
savep(p,'chromosome_counts',10,5)
top<-head(go[ontology=='BP' & query_hits>0][order(PValue)],15);top[,Description:=factor(Description,levels=rev(Description))]
p<-ggplot(top,aes(fold_enrichment,Description,size=query_hits,colour=BH_FDR))+geom_point()+scale_colour_gradient(low='#7B3294',high='#ABB6BC',limits=c(0,1))+
 labs(title='GO biological process: exact 70-gene set',subtitle='Top terms by nominal P; FDR includes all eligible BP terms',x='Fold enrichment',y=NULL,size='Query genes',colour='BH FDR')
savep(p,'GO_BP_top_terms',10,7)
positions[,chrom:=factor(chrom,levels=rev(chr_order))]
p<-ggplot(positions,aes(gene_start/1e6,chrom,colour=ZNF_ZFP_named))+geom_point(size=2,alpha=.8)+
 scale_colour_manual(values=c('FALSE'='#7B3294','TRUE'='#D55E00'),labels=c('Other genes','ZNF/ZFP-named'))+
 labs(title='Genomic positions of the 70 down/down genes',subtitle='Positions are RNA-source gene starts; nearby points can represent shared promoter peaks',x='Gene start (Mb)',y=NULL,colour=NULL)+theme(legend.position='bottom')
savep(p,'gene_positions',10,7)
summary<-list(query_genes=70,background_genes=nrow(d),ATAC_FDR=.1,RNA_FDR=.05,chromosome_significant=chr[BH_FDR<.05,chrom],
 GO_accounting=rbindlist(account),ZNF_ZFP=znf,global_clustering=cluster,significant_down_promoter_peaks=uniqueN(links$peak_id),
 query_genes_with_significant_down_peak=uniqueN(links$entrez_id),shared_promoter_peaks=nrow(sharing[genes>1]))
write_json(summary,file.path(out,'summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(c('Exact displayed down/down set: 70 unique Entrez genes, ATAC minimum linked promoter peak FDR <0.10, RNA FDR <0.05, both displayed effects negative.',
 'Background: all 11,898 eligible promoter genes from the same scatter; no whole-genome background or extra expression filter.',
 'GO: org.Hs.eg.db GOALL ancestor annotations, BP primary; MF and CC secondary. Ontology-annotated genes form each effective background/query. Terms of size 10-500 in that background are tested, including zero-hit terms.',
 'Unstratified enrichment uses upper-tail hypergeometric tests. BH covers all eligible terms within each ontology; correction across all ontologies is also supplied. These correction families differ from the old 0.10/0.10 up/down analysis.',
 'Chromosome-conditioned GO uses an exact convolution of chromosome-specific hypergeometric distributions, holding the observed query count on every chromosome fixed; same BH families.',
 'Chromosome enrichment uses upper-tail hypergeometric tests against eligible gene counts; BH across all represented chromosomes.',
 'ZNF/ZFP-name enrichment is one prespecified descriptive family test, also chromosome-conditioned; it is not an exhaustive zinc-finger protein annotation.',
 'Within-chromosome clustering: pair counts at <=1, <=5 and <=10 Mb between gene starts. 10,000 samples without replacement from eligible genes, holding observed per-chromosome query counts fixed; seed 20260909. Empirical upper-tail P=(1+exceedances)/(B+1). BH across three global tests and separately across all chromosome/distance tests.',
 '1 Mb-gap blocks are descriptive single-linkage groups, not independently significant loci. They may span more than 1 Mb.',
 'Neighboring genes may share promoter peaks; gene-level selection is not fully independent. Chromosome conditioning does not resolve this or promoter-peak-number selection bias.',
 'ATAC technical batch remains unconfirmed; no lpWGS copy-number correction, causal mechanism or cell identity is inferred.'),file.path(out,'METHODS.txt'))
stopifnot(all(go$PValue>=0 & go$PValue<=1),all(go$chromosome_conditioned_PValue>=0 & go$chromosome_conditioned_PValue<=1),sum(chr$query_genes)==70,
 all.equal(conditional_tail(30,100,10,5),phyper(4,30,70,10,lower.tail=FALSE))==TRUE)
writeLines(c('PASS: exact 70 IDs reproduce current plot selection, original RNA joins and full eligible background.',
 'PASS: chromosome counts sum to 70; all gene coordinates present; effective GO universes and all tested terms exported.',
 'PASS: conditioned hypergeometric reduces to the standard test with one stratum; probabilities remain within [0,1].',
 'PASS: permutation queries preserve observed chromosome counts; zero-hit GO terms remain in correction families.'),file.path(out,'verification.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(chr[order(PValue)][1:5]);print(znf);print(cluster);print(head(go[ontology=='BP'][order(PValue),.(Description,query_hits,BH_FDR,chromosome_conditioned_BH_FDR,query_symbols)],10))
