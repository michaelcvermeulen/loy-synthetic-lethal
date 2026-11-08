options(warn=1)
.libPaths(c('C:/Users/micha/AppData/Local/R/win-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(DESeq2);library(data.table);library(jsonlite);library(digest)})
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
ROOT <- normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
OUT <- file.path(ROOT,'integration/RNA_interaction_CN')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
src <- 'A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080_BATCHES/count_matrix'
paths <- file.path(src,c('HT1080_BATCHES_counts_raw.csv','HT1080_BATCHES_gene_annotation.csv','HT1080_BATCHES_sample_metadata.csv'))
raw <- fread(paths[1]);anno <- fread(paths[2],colClasses=list(character='entrezid'),na.strings=c('','NA'));meta <- fread(paths[3])
stopifnot(identical(raw$gene_id,anno$gene_id),!anyDuplicated(raw$gene_id))
meta <- meta[treatment %in% c('CTRL','LOY')][order(batch,treatment,sample_id)]
meta[,batch:=factor(batch,levels=c('AUG_2025','MAY_2026'))]
meta[,treatment:=factor(treatment,levels=c('CTRL','LOY'))]
stopifnot(nrow(meta)==12,all(table(meta$batch,meta$treatment)==3))
sample_names <- meta$colname
mat <- as.matrix(raw[,..sample_names]);rownames(mat)<-raw$gene_id
keep <- !is.na(anno$gene_name)&!is.na(anno$entrezid)&rowSums(mat>=10)>=3
fwrite(data.table(gene_id=raw$gene_id,annotation_eligible=!is.na(anno$gene_name)&!is.na(anno$entrezid),samples_count_ge10=rowSums(mat>=10),included=keep),file.path(OUT,'RNA_filter_audit.csv'))
mat <- mat[keep,,drop=FALSE];storage.mode(mat)<-'integer'
cd <- as.data.frame(meta);rownames(cd)<-cd$colname
design <- model.matrix(~batch*treatment,cd)
stopifnot(qr(design)$rank==4L,nrow(design)==12L)
fwrite(meta,file.path(OUT,'RNA_model_sample_metadata.csv'))
write.csv(design,file.path(OUT,'RNA_design_matrix.csv'))
dds <- DESeqDataSetFromMatrix(mat,cd,~batch*treatment)
dds <- DESeq(dds,betaPrior=FALSE,minReplicatesForReplace=Inf,parallel=FALSE)
rn <- resultsNames(dds);print(rn)
interaction <- grep('batch.*treatment',rn,value=TRUE)
stopifnot(length(interaction)==1L,'treatment_LOY_vs_CTRL' %in% rn)
res <- results(dds,name=interaction,alpha=.05,independentFiltering=FALSE)
aug <- results(dds,name='treatment_LOY_vs_CTRL',alpha=.05,independentFiltering=FALSE)
may <- results(dds,contrast=list(c('treatment_LOY_vs_CTRL',interaction)),alpha=.05,independentFiltering=FALSE)
stopifnot(isTRUE(all.equal(as.numeric(may$log2FoldChange-aug$log2FoldChange),as.numeric(res$log2FoldChange),tolerance=1e-7)))
table <- as.data.table(anno[match(rownames(dds),gene_id)])
table[,baseMean:=res$baseMean]
for(item in list(list('interaction',res),list('August',aug),list('May',may))) {
  prefix<-item[[1]];r<-item[[2]]
  for(k in c('log2FoldChange','lfcSE','stat','pvalue','padj'))table[,(paste0(prefix,'_',k)):=as.numeric(r[[k]])]
}
table[,`:=`(interaction_CI_low=interaction_log2FoldChange-1.96*interaction_lfcSE,interaction_CI_high=interaction_log2FoldChange+1.96*interaction_lfcSE,
  beta_converged=as.logical(mcols(dds)$betaConv),maxCooks=as.numeric(mcols(dds)$maxCooks))]
fwrite(table,file.path(OUT,'RNA_interaction_all_genes.csv'))
fwrite(table[!is.na(interaction_padj)&interaction_padj<.05][order(interaction_padj)],file.path(OUT,'RNA_interaction_FDR005.csv'))
fwrite(data.table(sample=colnames(dds),size_factor=sizeFactors(dds)),file.path(OUT,'RNA_interaction_size_factors.csv'))
fwrite(as.data.table(counts(dds,normalized=TRUE),keep.rownames='gene_id'),file.path(OUT,'RNA_interaction_normalized_counts.csv'))
saveRDS(dds,file.path(OUT,'RNA_interaction_dds.rds'))

# Five previously selected pathway questions, using a complementary sample-level
# model to account for correlation between genes. Separate from DESeq2 gene tests.
suppressPackageStartupMessages({library(edgeR);library(limma)})
y <- DGEList(counts=mat);y <- calcNormFactors(y)
pdf(file.path(OUT,'RNA_voom_diagnostic.pdf'),width=7,height=5)
v <- voom(y,design,plot=TRUE);dev.off()
membership <- fread(file.path(ROOT,'integration/pathway_GSEA/sources/gene_set_membership.csv'),colClasses=list(character='entrez_id'))
pathways <- c('HALLMARK_E2F_TARGETS','HALLMARK_MYC_TARGETS_V1','HALLMARK_MYC_TARGETS_V2','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE')
# Prevent ambiguous Entrez mappings from duplicating a gene in a set.
mapped <- table$entrezid
unique_map <- !duplicated(mapped)&!duplicated(mapped,fromLast=TRUE)
indices <- lapply(pathways,function(p)which(unique_map&mapped %in% membership[pathway==p,entrez_id]));names(indices)<-pathways
stopifnot(all(lengths(indices)>=15))
camera_default <- as.data.table(camera(v,index=indices,design=design,contrast=ncol(design),inter.gene.cor=.01),keep.rownames='pathway')
camera_estimated <- as.data.table(camera(v,index=indices,design=design,contrast=ncol(design),inter.gene.cor=NA,allow.neg.cor=FALSE),keep.rownames='pathway')
fwrite(camera_default,file.path(OUT,'RNA_interaction_CAMERA_fixed_correlation.csv'))
fwrite(camera_estimated,file.path(OUT,'RNA_interaction_CAMERA_estimated_correlation.csv'))
fwrite(rbindlist(lapply(names(indices),function(p)data.table(pathway=p,gene_id=rownames(v)[indices[[p]]],entrez_id=mapped[indices[[p]]]))),file.path(OUT,'RNA_interaction_CAMERA_membership.csv'))
saveRDS(list(voom=v,design=design,indices=indices),file.path(OUT,'RNA_interaction_CAMERA_inputs.rds'))

methods <- c('Primary gene test: DESeq2 ~ batch + treatment + batch:treatment on 12 raw-count samples (three CTRL and three LOY per batch). August 2025 and CTRL are reference levels. Parental samples excluded because absent from May.',
 'Interaction log2FC = (May LOY / CTRL log2FC) - (August LOY / CTRL log2FC). Positive means the LOY response is more positive in May; negative means more negative in May. A nonsignificant interaction is not an equivalence test.',
 'Annotation follows the original source count matrix. Retain annotated genes with Entrez ID and >=10 counts in >=3 of these 12 samples. Estimate DESeq2 median-ratio size factors anew on this matrix. Raw counts are not batch-corrected. Existing combined and separate RNA results remain unchanged.',
 'Wald tests use unshrunk coefficients and default Cook outlier filtering, with no automatic count replacement. BH correction spans all genes with a finite P-value for each contrast separately; independent filtering disabled. Primary interaction cutoff FDR <0.05. 95% Wald intervals are unadjusted. No fold-change threshold imposed.',
 'August and May coefficients exported here are contrasts from this new joint interaction fit, sharing its dispersions and normalization. They need not equal the previous separately fitted batch results.',
 'Samples are treated as independent. Sample-number/passage labels do not establish repeated measurements or clone pairing. CTRL in August and sgCTRL in May differ in experimental context. Interaction cannot isolate technical batch, passage, control manipulation, or biological causes.',
 'Complementary pathway interaction: limma-voom with edgeR TMM normalization on the same raw-count gene matrix and four-column design. CAMERA tests five previously selected Hallmark sets. Default fixed inter-gene correlation 0.01 is primary; estimated nonnegative correlation is sensitivity. BH within these five questions. Up/Down denotes relative enrichment of more positive/more negative May-minus-August responses. Ambiguous Entrez mappings excluded from gene-set membership. Gene-set tests are competitive and do not measure absolute pathway activity.',
 'No copy-number adjustment is applied to the RNA or ATAC model. WGS is assessed separately as line-level context.',
 'DESeq2 interaction documentation: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions',
 'CAMERA: https://doi.org/10.1093/nar/gks461')
writeLines(methods,file.path(OUT,'RNA_INTERACTION_METHODS.txt'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'RNA_interaction_sessionInfo.txt'))
fwrite(data.table(path=c(paths,script),sha256=vapply(c(paths,script),digest,character(1),algo='sha256',file=TRUE)),file.path(OUT,'RNA_interaction_source_hashes.csv'))
summary <- list(status='model_complete',samples=ncol(dds),tested_genes=nrow(dds),finite_interaction_p=sum(is.finite(res$pvalue)),interaction_FDR005=sum(res$padj<.05,na.rm=TRUE),interaction_FDR010=sum(res$padj<.1,na.rm=TRUE),nonconverged=sum(!mcols(dds)$betaConv),design_columns=rn)
write_json(summary,file.path(OUT,'RNA_interaction_model_summary.json'),pretty=TRUE,auto_unbox=TRUE)
print(summary);print(camera_default);print(camera_estimated)
