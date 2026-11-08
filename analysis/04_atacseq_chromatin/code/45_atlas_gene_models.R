options(warn = 1)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = '/', mustWork = TRUE))))
suppressPackageStartupMessages({library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(data.table)})
out <- 'integration/gene_atlas/data'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- fread('config/genes_knownGene.tsv', colClasses = list(character = 'gene_id'))
current <- as.data.table(GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE))
stopifnot(all(genes$gene_id %in% current$gene_id))
j <- match(genes$gene_id, current$gene_id)
stopifnot(all(genes$start == current$start[j]), all(genes$end == current$end[j]), all(genes$chrom == as.character(current$seqnames[j])))
symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes$gene_id, keytype = 'ENTREZID', column = 'SYMBOL', multiVals = 'asNA')
genes[, symbol := unname(symbols[gene_id])]
genes[is.na(symbol), symbol := gene_id]
fwrite(genes, file.path(out, 'gene_models.csv'))
ex <- reduce(exonsBy(txdb, by = 'gene'), ignore.strand = FALSE)
flat <- as.data.table(unlist(ex, use.names = FALSE))
flat[, gene_id := rep(names(ex), lengths(ex))]
flat <- flat[gene_id %in% genes$gene_id, .(gene_id, chrom = as.character(seqnames), start0 = start - 1L, end, strand = as.character(strand))]
fwrite(flat, file.path(out, 'merged_gene_exons.csv'))
fwrite(as.data.table(AnnotationDbi::metadata(txdb)), file.path(out, 'gene_model_metadata.csv'))
writeLines(c('Gene bodies match every coordinate in config/genes_knownGene.tsv exactly.',
  'Exons are unioned across all transcripts of each gene; the diagram is a collapsed gene model, not a single canonical transcript.',
  'Gene bodies in gene_models.csv are 1-based closed. Exon start0/end are 0-based half-open.',
  capture.output(sessionInfo())), file.path(out, 'gene_model_provenance.txt'))
cat('Exported', nrow(genes), 'gene bodies and', nrow(flat), 'merged exon intervals.\n')
