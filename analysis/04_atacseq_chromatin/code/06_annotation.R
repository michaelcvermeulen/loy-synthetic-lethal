options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(AnnotationDbi); library(org.Hs.eg.db)
  library(clusterProfiler); library(ggplot2)
})
genes <- fread("config/genes_knownGene.tsv", data.table = FALSE)
genes$gene_id <- as.character(genes$gene_id)
genes$symbol <- unname(mapIds(org.Hs.eg.db, keys = genes$gene_id,
                             keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"))
result <- fread("differential/primary_TMM_all.csv", data.table = FALSE)
peaks <- GRanges(result$chrom, IRanges(result$start+1L, result$end), peak_id = result$peak_id)
tss <- GRanges(genes$chrom, IRanges(genes$tss, genes$tss), strand = genes$strand,
               gene_id = genes$gene_id)
promoters <- GRanges(genes$chrom, IRanges(pmax(1L, genes$tss-2000L), genes$tss+2000L),
                     strand = genes$strand, gene_id = genes$gene_id)
nearest <- distanceToNearest(peaks, tss, ignore.strand = TRUE)
gene_index <- rep(NA_integer_, length(peaks)); gene_index[queryHits(nearest)] <- subjectHits(nearest)
result$nearest_entrez_id <- genes$gene_id[gene_index]
result$nearest_symbol <- genes$symbol[gene_index]
summit <- result$start+251L
result$nearest_tss_distance_signed <- ifelse(genes$strand[gene_index] == "+",
  summit-genes$tss[gene_index], genes$tss[gene_index]-summit)
overlap <- findOverlaps(peaks, promoters, ignore.strand = TRUE)
result$promoter_overlap <- seq_along(peaks) %in% queryHits(overlap)
gene_body <- GRanges(genes$chrom, IRanges(genes$start, genes$end))
result$annotation <- ifelse(result$promoter_overlap, "promoter",
                            ifelse(overlapsAny(peaks, gene_body, ignore.strand = TRUE), "gene_body", "intergenic"))
write.csv(result, "annotation/annotated_differential_peaks.csv", row.names = FALSE)
links <- cbind(result[queryHits(overlap), c("peak_id", "chrom", "start", "end", "logFC", "FDR", "direction")],
               genes[subjectHits(overlap), c("gene_id", "symbol", "tss", "strand")])
links$link_evidence <- "promoter_overlap_TSS_plus_minus_2kb"
write.csv(links, "annotation/promoter_peak_gene_links.csv", row.names = FALSE)
# One descriptive row per promoter-linked gene for later assay integration.
# These summaries do not create gene-level hypothesis tests from peak P-values.
gene_summary <- as.data.table(links)[, .(
  symbol = paste(unique(na.omit(symbol)), collapse = ";"),
  tested_promoter_peaks = uniqueN(peak_id),
  median_peak_logFC = median(logFC),
  minimum_peak_FDR = min(FDR),
  gained_peaks = sum(direction == "gained"),
  lost_peaks = sum(direction == "lost"),
  peak_ids = paste(unique(peak_id), collapse = ";")
), by = gene_id]
gene_summary[, promoter_pattern := fifelse(gained_peaks > 0 & lost_peaks > 0, "mixed",
  fifelse(gained_peaks > 0, "gained", fifelse(lost_peaks > 0, "lost", "ns")))]
write.csv(gene_summary, "annotation/promoter_gene_summary.csv", row.names = FALSE)
universe <- unique(links$gene_id)
enrichment <- list()
for (direction in c("gained", "lost")) {
  query <- unique(links$gene_id[links$direction == direction])
  if (length(query) < 5) next
  go <- enrichGO(gene = query, universe = universe, OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  table <- as.data.frame(go)
  write.csv(table, paste0("annotation/GO_BP_promoter_", direction, ".csv"), row.names = FALSE)
  if (nrow(table)) {
    show <- head(table[order(table$p.adjust), ], 15)
    show$Description <- factor(show$Description, levels = rev(show$Description))
    p <- ggplot(show, aes(-log10(pmax(p.adjust, .Machine$double.xmin)), Description)) +
      geom_point(aes(size = Count), colour = if (direction == "gained") "#D55E00" else "#0072B2") +
      theme_bw(base_size = 11) + labs(x = "-log10 adjusted P", y = NULL)
    ggsave(paste0("annotation/GO_BP_promoter_", direction, ".png"), p, width = 10, height = 7, dpi = 200)
  }
}
writeLines(c("Gene identifiers: UCSC knownGene / Entrez, pinned in config/genes_knownGene.tsv.",
             "Promoter links preserve every gene overlapping a tested peak; distal nearest genes are hypotheses.",
             "GO-BP foreground: genes with a significant gained/lost promoter peak; universe: genes with any tested promoter peak.",
             "Peak FDRs remain peak-level; they are not gene-level differential tests.",
             "Promoter gene summary: median effect across tested promoter peaks, minimum peak FDR, and gained/lost counts; descriptive, not gene-level FDR.",
             "Genes can have both gained and lost peaks. Promoter restriction reduces but does not remove peak-number bias.",
             "RNA annotation and sample matching are required before final cross-assay linkage."), "annotation/METHODS.txt")
summary <- as.data.frame(table(result$annotation, result$direction))
names(summary) <- c("annotation", "direction", "peaks")
write.csv(summary, "annotation/annotation_summary.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "logs/annotation_session.txt")
