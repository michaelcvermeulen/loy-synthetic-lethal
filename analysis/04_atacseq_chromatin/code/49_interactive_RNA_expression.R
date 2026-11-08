options(warn = 1)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash='/', mustWork=TRUE))))
suppressPackageStartupMessages({library(DESeq2); library(data.table); library(digest)})
out <- 'integration/interactive_atac_rna'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
src <- 'A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080_BATCHES/count_matrix'
counts <- fread(file.path(src,'HT1080_BATCHES_counts_raw.csv'))
anno <- fread(file.path(src,'HT1080_BATCHES_gene_annotation.csv'), colClasses=list(character='entrezid'), na.strings=c('','NA'))
meta <- fread(file.path(src,'HT1080_BATCHES_sample_metadata.csv'))
stopifnot(identical(counts$gene_id,anno$gene_id),identical(names(counts)[-1],meta$colname))
keep <- !is.na(anno$gene_name) & !is.na(anno$entrezid)
mat <- as.matrix(counts[keep,-1]); rownames(mat) <- counts$gene_id[keep]
mat <- mat[rowSums(mat >= 10) >= 3,,drop=FALSE]
results <- c(Combined='config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv',
  AUG_2025='integration/batch_RNA_overlaps/sources/ht1080_AUG_2025_LOY_v_CTRL_DESeq2.csv',
  MAY_2026='integration/batch_RNA_overlaps/sources/ht1080_MAY_2026_LOY_v_CTRL_DESeq2.csv')
factors <- list(); normalized <- list(); checks <- list()
for (name in 'Combined') {
  m <- if(name=='Combined') meta else meta[batch==name & treatment %in% c('CTRL','LOY')]
  cm <- mat[,m$colname,drop=FALSE]
  sf <- estimateSizeFactorsForMatrix(cm)
  nc <- sweep(cm,2,sf,'/')
  result <- fread(results[[name]],colClasses=list(character=c('entrezid','gene_id_full')))
  stopifnot(setequal(rownames(cm),result$gene_id_full))
  recomputed <- rowMeans(nc)[match(result$gene_id_full,rownames(nc))]
  cat(name,'max relative baseMean difference:',max(abs(recomputed-result$logCPM)/pmax(1,result$logCPM)),'\n')
  stopifnot(isTRUE(all.equal(recomputed,result$logCPM,tolerance=1e-8,check.attributes=FALSE)))
  checks[[name]] <- data.table(RNA_result=name,genes=nrow(cm),samples=ncol(cm),
    max_relative_baseMean_difference=max(abs(recomputed-result$logCPM)/pmax(1,result$logCPM)))
  factors[[name]] <- data.table(RNA_result=name,sample=m$colname,size_factor=sf)
  nc <- as.data.table(nc,keep.rownames='gene_id_full')
  nc[,RNA_result:=name]
  normalized[[name]] <- nc
}
fwrite(rbindlist(factors),'integration/interactive_atac_rna/RNA_size_factors.csv')
fwrite(rbindlist(normalized,fill=TRUE),'integration/interactive_atac_rna/RNA_normalized_counts.csv')
fwrite(rbindlist(checks),'integration/interactive_atac_rna/RNA_normalization_verification.csv')
fwrite(meta,file.path(out,'RNA_sample_metadata.csv'))
fwrite(counts[gene_id %in% rownames(mat)],file.path(out,'RNA_raw_counts.csv'))
paths <- c(file.path(src,c('HT1080_BATCHES_counts_raw.csv','HT1080_BATCHES_gene_annotation.csv','HT1080_BATCHES_sample_metadata.csv')),unname(results))
fwrite(data.table(path=paths,sha256=vapply(paths,digest,character(1),algo='sha256',file=TRUE)),file.path(out,'RNA_source_manifest.csv'))
writeLines(c('DESeq2 median-ratio size factors are reconstructed on exactly the combined model prefiltered count matrix and shared by all expression panels.',
  'No differential-expression model is refitted. Normalized-count row means reproduce every original combined-model DESeq2 baseMean (source column logCPM) within tolerance.',
  'Combined normalization includes all 15 original model samples; expression plots display only 12 CTRL/LOY samples, split by RNA batch.',
  'Selecting a separate RNA batch subsets samples but keeps the combined-model size factors for a consistent expression scale. These are not the separate-fit normalization factors. Normalized counts are not batch-corrected. Original logFC and FDR are unchanged.',
  capture.output(sessionInfo())),file.path(out,'RNA_expression_METHODS.txt'))
print(rbindlist(checks))
