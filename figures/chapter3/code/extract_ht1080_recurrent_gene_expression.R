suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
})

analysis_root <- file.path('RNAseq', 'analysis', 'HT1080', 'SEPT_11_2026_REANALYZE')
figure_root <- file.path('SEPT_11_2026_THESIS', 'figures', 'chapter3')
output <- file.path(figure_root, 'source_data', 'ht1080_recurrent_gene_expression')
dir.create(output, recursive = TRUE, showWarnings = FALSE)

genes <- c('WNT7A', 'PAX2', 'CHRNB2', 'GATA4', 'BMP7', 'IHH',
           'NLGN4X', 'MAPK8IP2', 'NRXN2')

se <- readRDS(file.path(analysis_root, 'outputs', 'models', 'raw_counts_SE.rds'))
counts <- assay(se, 'counts')
annotation <- as.data.frame(rowData(se))
metadata <- as.data.frame(colData(se))

group <- interaction(metadata$batch, metadata$treatment, drop = TRUE)
expressed <- filterByExpr(counts, group = group, min.count = 10, min.total.count = 15)
keep <- expressed & annotation$main_background
y <- calcNormFactors(DGEList(counts[keep, , drop = FALSE]), method = 'TMM')
effective <- y$samples$lib.size * y$samples$norm.factors
log_cpm <- cpm(counts, lib.size = effective, log = TRUE, prior.count = 2)
raw_cpm <- cpm(counts, lib.size = effective, log = FALSE)

de <- read.csv(file.path(analysis_root, 'outputs', 'programme_analysis_v1',
                         'tables', 'gene_consistency.csv'), check.names = FALSE)
selected <- de[match(genes, de$gene_name), ]
stopifnot(!anyNA(selected$gene_id), identical(selected$gene_name, genes))
indices <- match(selected$gene_id, rownames(counts))
stopifnot(!anyNA(indices))

rows <- list()
for (i in seq_along(genes)) {
  for (j in seq_len(ncol(counts))) {
    rows[[length(rows) + 1L]] <- data.frame(
      gene = genes[i],
      gene_id = selected$gene_id[i],
      sample_key = metadata$sample_key[j],
      batch = metadata$batch[j],
      treatment = metadata$treatment[j],
      passage_number = metadata$passage_number[j],
      log2_CPM = log_cpm[indices[i], j],
      CPM = raw_cpm[indices[i], j],
      raw_count = counts[indices[i], j],
      stringsAsFactors = FALSE
    )
  }
}
expression <- do.call(rbind, rows)
expression$condition <- factor(expression$treatment,
  levels = c('PARENTAL', 'CTRL', 'LOY'),
  labels = c('Parental', 'sgCTRL', 'Y-depleted'))
expression$batch_label <- factor(expression$batch,
  levels = c('AUG_2025', 'MAY_2026'), labels = c('August', 'May'))
write.csv(expression, file.path(output, 'sample_expression.csv'), row.names = FALSE)

statistics <- selected[, c('gene_id', 'gene_name', 'logFC', 'FDR',
                            'AUG_logFC', 'AUG_FDR', 'MAY_logFC', 'MAY_FDR')]
names(statistics)[names(statistics) == 'logFC'] <- 'LOY_vs_sgCTRL_log2FC'
names(statistics)[names(statistics) == 'FDR'] <- 'LOY_vs_sgCTRL_FDR'
write.csv(statistics, file.path(output, 'batch_aware_statistics.csv'), row.names = FALSE)

normalization <- data.frame(
  sample_key = metadata$sample_key,
  batch = metadata$batch,
  treatment = metadata$treatment,
  effective_library_reads = effective,
  all_feature_reads = colSums(counts)
)
write.csv(normalization, file.path(output, 'all15_context_normalization.csv'), row.names = FALSE)

cat('Exported', nrow(expression), 'gene-sample observations for', length(genes), 'genes\n')
