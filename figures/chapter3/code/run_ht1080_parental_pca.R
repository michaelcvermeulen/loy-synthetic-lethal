# Reproduce only the all-15 contextual PCA using the author-selected inputs.
args <- commandArgs(trailingOnly=TRUE)
source_dir <- normalizePath(args[1], winslash="/", mustWork=TRUE)
out_dir <- normalizePath(args[2], winslash="/", mustWork=TRUE)
suppressPackageStartupMessages(library(edgeR))
x <- readRDS(file.path(source_dir,"outputs/models/01_ingest.rds"))
e <- readRDS(file.path(source_dir,"outputs/models/02_expression.rds"))
cfg <- yaml::read_yaml(file.path(source_dir,"config/analysis.yml"))
source(file.path(source_dir,"R/02_effects.R"))
stopifnot(ncol(x$counts)==15, identical(colnames(x$counts), x$meta$sample_key))
context <- normalize_counts(x$counts, x$annotation$main_background,
    interaction(x$meta$batch,x$meta$treatment,drop=TRUE),cfg)
keep <- e$norm$keep
stopifnot(sum(keep)==25526)
pca <- prcomp(t(context$logCPM[keep,]),center=TRUE,scale.=FALSE)
# Principal-component signs are arbitrary; align to the existing saved output.
previous <- read.csv(file.path(source_dir,"outputs/tables/PCA_context_all15.csv"),check.names=FALSE)
previous <- previous[match(rownames(pca$x),previous$sample_key),]
for (j in 1:4) {
  if (cor(pca$x[,j], previous[[paste0("PC",j)]]) < 0) {
    pca$x[,j] <- -pca$x[,j]; pca$rotation[,j] <- -pca$rotation[,j]
  }
}
max_difference <- max(abs(pca$x[,1:4] - as.matrix(previous[,paste0("PC",1:4)])))
stopifnot(max_difference < 1e-6)
variance <- pca$sdev^2 / sum(pca$sdev^2)
write.csv(cbind(x$meta,as.data.frame(pca$x[,1:4])),file.path(out_dir,"PCA_parental_all15.csv"),row.names=FALSE)
write.csv(data.frame(component=paste0("PC",seq_along(variance)),variance_fraction=variance),file.path(out_dir,"PCA_variance.csv"),row.names=FALSE)
write.csv(data.frame(sample_key=x$meta$sample_key,context$dge$samples,effective_library_reads=context$effective),file.path(out_dir,"normalization_all15.csv"),row.names=FALSE)
write.csv(data.frame(gene_id=x$annotation$gene_id,core_PCA_included=keep,context_TMM_included=context$keep),file.path(out_dir,"gene_selection.csv"),row.names=FALSE)
saveRDS(pca,file.path(out_dir,"PCA_parental_all15.rds"))
writeLines(c(paste("Libraries:",ncol(x$counts)),paste("PCA genes:",sum(keep)),paste("Context TMM genes:",sum(context$keep)),paste("Maximum difference from saved PC1-PC4:",max_difference),capture.output(sessionInfo())),file.path(out_dir,"run_validation.txt"))
cat("PC1:",round(100*variance[1],1),"%; PC2:",round(100*variance[2],1),"%; maximum saved-coordinate difference:",max_difference,"\n")
