# Follow-up to 81_quadrant_ORA.R (22 Sep 2026): is PRC2 / H3K27me3-target enrichment a property of promoter
# closing in general, or of genes that close AND lose expression? Tabulates membership of three CGP chromatin sets
# by RNA direction x ATAC direction (|log2FC| > 0.25, no FDR filter) and runs Fisher tests:
#   (1) ATAC-closed vs other genes, all genes; (2) the same within RNA-flat genes (|RNA log2FC| <= 0.25).
# Run: Rscript --vanilla <absolute path to this script>
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(msigdbr)})
out <- "integration/thesis_RNA_SEPT11/quadrant_ORA_rna025_atac025"
g <- fread(file.path(out, "gene_quadrants.csv"), colClasses = list(character = "entrez_id"))
g[, atac := fcase(atac_median_logFC < -0.25, "closed", atac_median_logFC > 0.25, "open", default = "flat")]
g[, rna := fcase(rna_logFC < -0.25, "down", rna_logFC > 0.25, "up", default = "flat")]
cgp <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "C2", subcollection = "CGP"))
rows <- list(); tabs <- list()
for (set in c("BENPORATH_ES_WITH_H3K27ME3", "BENPORATH_PRC2_TARGETS", "MEISSNER_NPC_HCP_WITH_H3K4ME2_AND_H3K27ME3")) {
  g[, member := entrez_id %in% as.character(cgp[gs_name == set, ncbi_gene])]
  tabs[[set]] <- g[, .(gene_set = set, n = .N, members = sum(member), pct = 100 * mean(member)), by = .(rna, atac)]
  f1 <- fisher.test(table(g$atac == "closed", g$member))
  f2 <- with(g[rna == "flat"], fisher.test(table(atac == "closed", member)))
  f3 <- with(g[rna == "flat"], fisher.test(table(atac == "open", member)))
  rows[[set]] <- data.table(gene_set = set, overall_pct = 100 * mean(g$member),
                            closed_vs_rest_OR = f1$estimate, closed_vs_rest_P = f1$p.value,
                            closed_vs_rest_RNAflat_OR = f2$estimate, closed_vs_rest_RNAflat_P = f2$p.value,
                            open_vs_rest_RNAflat_OR = f3$estimate, open_vs_rest_RNAflat_P = f3$p.value)
}
tab <- rbindlist(tabs); res <- rbindlist(rows)
fwrite(tab, file.path(out, "PRC2_membership_by_quadrant.csv")); fwrite(res, file.path(out, "PRC2_fisher_tests.csv"))
for (set in unique(tab$gene_set)) {
  cat("\n", set, "(% members; rows RNA, columns ATAC)\n")
  print(dcast(tab[gene_set == set], rna ~ atac, value.var = "pct")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 1) else v)])
}
print(res[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 3) else v)])
