# Export original counts and saved normalization for the nicotinic follow-up.
# No normalization, model fitting or enrichment testing is repeated.
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))
root <- dirname(normalizePath(script, winslash = '/', mustWork = TRUE))
set.seed(42)
x <- readRDS(file.path(root, 'outputs/models/01_ingest.rds'))
e <- readRDS(file.path(root, 'outputs/models/02_expression.rds'))
genes <- c(paste0('CHRNA', c(1:7, 9:10)), paste0('CHRNB', 1:4),
  'CHRND', 'CHRNE', 'CHRNG', 'CHRFAM7A', 'CHAT', 'ACHE', 'BCHE',
  'SLC18A3', 'SLC5A7', 'RIC3', 'TMEM35A', 'LYNX1', 'LY6H')
selected <- e$effects[match(genes, e$effects$gene_name), ]
stopifnot(nrow(selected) == length(genes), !anyNA(selected$gene_id),
  !anyDuplicated(selected$gene_id), nrow(x$core) == 12)
ids <- selected$gene_id
stopifnot(all(ids %in% rownames(x$counts)), all(ids %in% rownames(e$logCPM)))
out <- file.path(root, 'outputs/tables/nicotinic_receptors')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sample_rows <- lapply(seq_len(nrow(x$core)), function(j) {
  sample <- x$core$sample_key[j]
  data.frame(gene_id = ids, gene_name = genes, sample_key = sample,
    batch = x$core$batch[j], treatment = x$core$treatment[j],
    passage = x$core$passage_number[j], count_unit = 'assigned_reads',
    raw_count = x$counts[ids, sample], CPM = e$norm$CPM[ids, sample],
    logCPM = e$logCPM[ids, sample], stringsAsFactors = FALSE)
})
write.csv(do.call(rbind, sample_rows), file.path(out, 'sample_expression.csv'), row.names = FALSE)
write.csv(selected, file.path(out, 'gene_effects.csv'), row.names = FALSE)
writeLines(c('Original 01_ingest.rds and 02_expression.rds reused.',
  '26 genes; 12 core libraries; counts are assigned reads, not paired fragments.',
  'TMM CPM and logCPM (prior count 2) copied from the original saved normalization.',
  'No new gene-level P values, normalization or enrichment tests.'),
  file.path(out, 'extraction_notes.txt'))
cat(nrow(selected), 'genes and', nrow(x$core), 'samples exported.\n')
