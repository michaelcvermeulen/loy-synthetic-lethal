suppressPackageStartupMessages(library(edgeR))
root <- 'A:/Dropbox/LOY_paralogs/projects'
analysis <- file.path(root,'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE')
fig <- file.path(root,'SEPT_11_2026_THESIS/figures/chapter3')
out <- file.path(fig,'source_data/ht1080_rna_expression_changes')
x <- readRDS(file.path(analysis,'outputs/models/01_ingest.rds'))
npath <- file.path(fig,'source_data/ht1080_recurrent_gene_expression/all15_context_normalization.csv')
n <- read.csv(npath)
m <- x$meta
effective <- n$effective_library_reads[match(m$sample_key,n$sample_key)]
stopifnot(nrow(m)==15,!anyNA(effective),all(effective>0))
genes <- c('UNC5A','ROR2','SDK2','COL1A2','FUCA2','MAGEA1')
de <- read.csv(file.path(analysis,'outputs/programme_analysis_v1/tables/gene_consistency.csv'))
ids <- de$gene_id[match(genes,de$gene_name)]
stopifnot(!anyNA(ids))
values <- cpm(x$counts[,m$sample_key],lib.size=effective,log=TRUE,prior.count=2)
z <- do.call(rbind,lapply(seq_along(genes),function(i)data.frame(gene=genes[i],gene_id=ids[i],sample_key=m$sample_key,treatment=m$treatment,batch=m$batch,passage=m$passage_number,log2_CPM=as.numeric(values[ids[i],]))))
stopifnot(nrow(z)==90,all(is.finite(z$log2_CPM)))
write.csv(z,file.path(out,'six_gene_sample_expression.csv'),row.names=FALSE)
writeLines(c('Saved all-15 contextual TMM effective library sizes reused; no period correction.','edgeR log2 CPM with prior.count=2. No DE refit.','Parental: August only (n=3); sgCTRL and LOY: August and May (n=6 each).',paste('Normalization source:',npath)),file.path(out,'six_gene_expression_notes.txt'))
cat('Saved 90 gene/sample expression values.\n')
