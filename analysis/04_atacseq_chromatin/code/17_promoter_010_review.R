options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(AnnotationDbi)
  library(org.Hs.eg.db); library(clusterProfiler); library(jsonlite)})
outdir <- "integration/promoter_FDR010_review"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
g <- fread("integration/FDR_sensitivity/joint_genes_ATAC_0.10_RNA_0.10.csv", colClasses = list(character = "entrez_id"))
links <- fread("integration/RNAseq/ATAC_RNA_promoter_peak_links.csv", colClasses = list(character = "gene_id"))
universe <- fread("integration/RNAseq/ATAC_RNA_gene_universe.csv", colClasses = list(character = "entrez_id"))
stopifnot(nrow(g) == 94, !anyDuplicated(g$entrez_id), all(g$rna_FDR < .1), all(g$atac_minimum_peak_FDR < .1))
eligible_links <- links[gene_id %in% g$entrez_id]
sig_links <- eligible_links[FDR < .1]
best <- sig_links[order(FDR, peak_id), .SD[1], by = gene_id]
idx <- match(g$entrez_id, best$gene_id)
g[, `:=`(chromosome = best$chrom[idx], best_promoter_peak = best$peak_id[idx],
  best_promoter_peak_logFC = best$logFC[idx], best_promoter_peak_FDR = best$FDR[idx])]
descriptions <- mapIds(org.Hs.eg.db, keys = g$entrez_id, keytype = "ENTREZID", column = "GENENAME", multiVals = "first")
g[, gene_description := unname(descriptions[entrez_id])]
g[, `:=`(RNA_direction = ifelse(rna_logFC > 0, "up", "down"),
  core_FDR005_both = rna_FDR < .05 & atac_minimum_peak_FDR < .05,
  ZNF_named_gene = grepl("^ZNF", rna_symbol))]
g[, status := fcase(category_at_cutoff == "concordant_gain", "Up in both",
  category_at_cutoff == "concordant_loss", "Down in both", category_at_cutoff == "discordant", "Opposing directions",
  default = "Mixed promoter")]
stopifnot(all(sign(g$best_promoter_peak_logFC[g$category_at_cutoff == "concordant_gain"]) == 1),
  all(sign(g$best_promoter_peak_logFC[g$category_at_cutoff == "concordant_loss"]) == -1))
setorder(g, rna_FDR, atac_minimum_peak_FDR)
fwrite(g, file.path(outdir, "all_94_genes.csv"))
for (category in unique(g$category_at_cutoff)) fwrite(g[category_at_cutoff == category], file.path(outdir, paste0(category, ".csv")))
fwrite(g[RNA_direction == "up"], file.path(outdir, "RNA_up_17_genes.csv"))
fwrite(g[RNA_direction == "down"], file.path(outdir, "RNA_down_77_genes.csv"))
fwrite(eligible_links, file.path(outdir, "all_promoter_peaks_for_selected_genes.csv"))
fwrite(sig_links, file.path(outdir, "significant_promoter_peaks_FDR010.csv"))
chromosome_summary <- g[, .(genes = .N, symbols = paste(rna_symbol, collapse = ";")), by = .(status, chromosome)]
fwrite(chromosome_summary, file.path(outdir, "chromosome_summary.csv"))
category_summary <- g[, .(genes = .N, retained_at_FDR005_both = sum(core_FDR005_both),
  RNA_abs_logFC_gt075 = sum(abs(rna_logFC) > .75), ZNF_named_genes = sum(ZNF_named_gene)), by = status]
fwrite(category_summary, file.path(outdir, "category_summary.csv"))
for (category in unique(g$category_at_cutoff)) {
  selected <- copy(g[category_at_cutoff == category])
  selected[, display_order := seq_len(.N)]
  selected[, panel := paste0("Genes ", ((display_order-1L)%/%38L)*38L+1L, "-",
    pmin((((display_order-1L)%/%38L)+1L)*38L, .N))]
  selected[, gene_label := factor(rna_symbol, levels = rev(rna_symbol))]
  long <- melt(selected, id.vars = c("entrez_id", "gene_label", "panel"),
    measure.vars = c("atac_median_promoter_logFC", "rna_logFC"), variable.name = "assay", value.name = "logFC")
  long[, assay := factor(assay, levels = c("atac_median_promoter_logFC", "rna_logFC"), labels = c("ATAC promoter median", "RNA"))]
  p <- ggplot(long, aes(logFC, gene_label)) +
    geom_segment(data = selected, aes(x = atac_median_promoter_logFC, xend = rna_logFC,
      y = gene_label, yend = gene_label), colour = "grey75", linewidth = .5, inherit.aes = FALSE) +
    geom_point(aes(colour = assay, shape = assay), size = 2) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    facet_wrap(~panel, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = c("ATAC promoter median" = "#0072B2", RNA = "#D55E00")) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom", legend.title = element_blank(),
      panel.grid.major.y = element_blank(), plot.caption = element_text(hjust = 0)) +
    labs(title = paste0(selected$status[1], " at promoter ATAC/RNA FDR <0.10 (n = ", nrow(selected), ")"),
      subtitle = "Ordered by RNA FDR; one row per gene", x = "log2 fold change (LOY / CTRL)", y = NULL,
      caption = "ATAC effect is the median across all tested promoter peaks; significance requires at least one peak at FDR <0.10.")
  height <- max(4.5, min(38, nrow(selected))*.23+2.2)
  width <- if (nrow(selected) > 38) 13 else 10
  ggsave(file.path(outdir, paste0(category, "_effects.png")), p, width = width, height = height, dpi = 220)
  ggsave(file.path(outdir, paste0(category, "_effects.pdf")), p, width = width, height = height, useDingbats = FALSE)
}
# Exploratory GO-BP: common assay-eligible universe, BH across both directions together.
go_tables <- list(); go_accounting <- list()
for (category in c("concordant_gain", "concordant_loss")) {
  query <- g[category_at_cutoff == category, entrez_id]
  result <- enrichGO(gene = query, universe = universe$entrez_id, OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1,
    qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  table <- as.data.table(as.data.frame(result))
  go_accounting[[category]] <- data.table(category, query_genes = length(query), eligible_universe = nrow(universe),
    GO_annotated_universe = length(result@universe), returned_terms = nrow(table))
  if (nrow(table)) {
    table[, category := category]
    table[, query_symbols := vapply(strsplit(geneID, "/", fixed = TRUE), function(ids) {
      paste(g$rna_symbol[match(ids, g$entrez_id)], collapse = ";")
    }, character(1))]
    go_tables[[category]] <- table
  }
}
go <- rbindlist(go_tables, fill = TRUE)
if (nrow(go)) {
  go[, BH_FDR_across_up_down := p.adjust(pvalue, method = "BH")]
  setorder(go, BH_FDR_across_up_down, pvalue)
  fwrite(go, file.path(outdir, "GO_BP_all_terms.csv"))
  fwrite(go[BH_FDR_across_up_down < .05], file.path(outdir, "GO_BP_significant_terms.csv"))
}
fwrite(rbindlist(go_accounting), file.path(outdir, "GO_accounting.csv"))
summary <- list(ATAC_FDR = .1, RNA_FDR = .1, total = nrow(g),
  up_both = sum(g$status == "Up in both"), down_both = sum(g$status == "Down in both"),
  opposing = sum(g$status == "Opposing directions"), RNA_up = sum(g$RNA_direction == "up"),
  RNA_down = sum(g$RNA_direction == "down"), GO_terms_FDR005 = if(nrow(go)) sum(go$BH_FDR_across_up_down < .05) else 0,
  GO_up_terms_FDR005 = if(nrow(go)) sum(go$BH_FDR_across_up_down < .05 & go$category == "concordant_gain") else 0,
  GO_down_terms_FDR005 = if(nrow(go)) sum(go$BH_FDR_across_up_down < .05 & go$category == "concordant_loss") else 0)
write_json(summary, file.path(outdir, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(c("Review selects the 94 common-universe promoter genes at ATAC peak FDR <0.10 and RNA gene FDR <0.10.",
  "Groups use the direction of significant promoter peaks and RNA. Both significant peak signs would be marked mixed; none occur here.",
  "RNA-up/down lists include opposing-direction genes; concordant lists require agreement between assays.",
  "Median ATAC effects summarize all tested promoter peaks. The best-FDR significant peak is recorded separately. Peak FDR is not gene-level FDR.",
  "Gene symbols are preserved from the RNA source; descriptive names are appended from installed org.Hs.eg.db.",
  "GO biological-process enrichment is exploratory: universe is 11,898 assay-eligible genes; unannotated genes are excluded by enrichGO.",
  "GO term sizes 10-500; BH correction is recomputed across returned up and down GO tests together. Query membership and annotated universe are saved.",
  "Joint selection can retain promoter-peak-number and gene-annotation bias. Related GO terms are redundant and do not demonstrate pathway activity.",
  "Chromosome clustering is tabulated; lpWGS context is not yet included. These results do not distinguish regulation from copy-number effects."),
  file.path(outdir, "METHODS.txt"))
stopifnot(summary$up_both == 14, summary$down_both == 76, summary$opposing == 4,
  summary$RNA_up == 17, summary$RNA_down == 77, sum(category_summary$genes) == 94)
writeLines(c("PASS: 94 unique genes reproduce the 0.10/0.10 promoter overlap", "PASS: 14 concordant up, 76 concordant down, four opposing; RNA partitions 17/77",
  "PASS: best significant promoter peaks satisfy FDR/direction requirements", "PASS: all original promoter links and source statistics retained"),
  file.path(outdir, "verification.txt"))
writeLines(capture.output(sessionInfo()), "logs/promoter_010_review_session.txt")
print(category_summary)
print(summary)
if (nrow(go)) print(head(go[, .(category, Description, Count, BH_FDR_across_up_down, query_symbols)], 12))
