x <- readRDS('RNAseq/analysis/crossline_SEPT_2026/outputs/models/01_ingest.rds')
out <- 'SEPT_11_2026_THESIS/figures/chapter3/source_data/shared_sex_chromosome'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
for (L in names(x$cp)) {
 e <- x$eff[[L]]; ids <- e$gene_id[e$chrom %in% c('chrX','chrY')]
 cp <- x$cp[[L]]; ids <- intersect(ids, rownames(cp$cpm))
 write.csv(data.frame(gene_id=ids,cp$cpm[ids,,drop=FALSE],check.names=FALSE),file.path(out,paste0(L,'_CPM.csv')),row.names=FALSE)
 write.csv(e[e$gene_id %in% ids,],file.path(out,paste0(L,'_effects.csv')),row.names=FALSE)
 write.csv(cp$meta,file.path(out,paste0(L,'_samples.csv')),row.names=FALSE)
}
