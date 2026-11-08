# Over-representation analysis of RNA-ATAC agreeing quadrants (thesis RNA branch), exploratory.
#
# Genes: the 13,300-gene promoter integration universe (integration/thesis_RNA_SEPT11, 80_integrate_thesis_RNA.R).
# Effects: RNA log2FC (thesis limma-voom) and gene-level promoter ATAC log2FC (median of tested promoter peaks).
# No FDR filter: quadrant membership uses sign plus a magnitude filter on both assays,
#   UP-OPEN     = RNA > +MIN_RNA and ATAC > +MIN_ATAC
#   DOWN-CLOSED = RNA < -MIN_RNA and ATAC < -MIN_ATAC
#   (UP-CLOSED and DOWN-OPEN, the disagreeing quadrants, are run as a comparison.)
# Backgrounds:
#   primary     = all 13,300 universe genes;
#   conditional = genes with the same RNA direction and |RNA log2FC| > MIN_RNA, regardless of ATAC. This asks
#                 whether promoter agreement adds enrichment beyond the RNA change itself.
# Collections: MSigDB Hallmark, Reactome, GO:BP, and C2:CGP (added 22 Sep 2026) (msigdbr 26.1.1, human Entrez); set size 10-500 within each
# background; clusterProfiler::enricher hypergeometric test, BH within collection x quadrant x background.
# Run: Rscript --vanilla <absolute path to this script> [MIN_RNA MIN_ATAC]   (defaults 0.25 0.25)
options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(msigdbr); library(clusterProfiler); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
MIN_RNA <- if (length(args) >= 1) as.numeric(args[1]) else 0.25
MIN_ATAC <- if (length(args) >= 2) as.numeric(args[2]) else 0.25
# Optional third argument "trended" (23 Sep 2026): gene ATAC effect = median promoter log2FC from the trended (csaw loess)
# model (88_trended_normalisation_sensitivity.R), which removes the abundance-dependent trend of the primary estimates.
TRENDED <- length(args) >= 3 && args[3] == "trended"
out <- sprintf("integration/thesis_RNA_SEPT11/quadrant_ORA_rna%03d_atac%03d%s", round(MIN_RNA * 100), round(MIN_ATAC * 100),
               if (TRENDED) "_trended" else "")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

g <- fread("integration/thesis_RNA_SEPT11/RNA_ATAC_promoter_genes.csv", colClasses = list(character = "entrez_id"))
stopifnot(nrow(g) == 13300, !anyDuplicated(g$entrez_id))
if (TRENDED) {
  trd <- fread("differential/trended_normalisation/trended_all.csv")[, .(peak_id, logFC)]
  gp <- g[, .(peak_id = unlist(strsplit(peak_ids, ";"))), by = entrez_id]
  gm <- merge(gp, trd, by = "peak_id")[, .(atac_trended = median(logFC)), by = entrez_id]
  g <- merge(g, gm, by = "entrez_id", all.x = TRUE)
  stopifnot(!anyNA(g$atac_trended))
  g[, atac_median_logFC_primary := atac_median_logFC][, atac_median_logFC := atac_trended]
}
g[, quadrant := fcase(rna_logFC > MIN_RNA & atac_median_logFC > MIN_ATAC, "UP_OPEN",
                      rna_logFC < -MIN_RNA & atac_median_logFC < -MIN_ATAC, "DOWN_CLOSED",
                      rna_logFC > MIN_RNA & atac_median_logFC < -MIN_ATAC, "UP_CLOSED",
                      rna_logFC < -MIN_RNA & atac_median_logFC > MIN_ATAC, "DOWN_OPEN",
                      default = "none")]
fwrite(g[, .(entrez_id, gene_name, chrom, rna_logFC, rna_FDR, rna_class, atac_median_logFC, atac_min_peak_FDR,
             tested_promoter_peaks, quadrant)], file.path(out, "gene_quadrants.csv"))
counts <- g[, .N, by = quadrant][order(quadrant)]
print(counts)

sets <- rbindlist(list(
  as.data.table(msigdbr(db_species = "HS", species = "human", collection = "H"))[, collection := "Hallmark"],
  as.data.table(msigdbr(db_species = "HS", species = "human", collection = "C2", subcollection = "CP:REACTOME"))[, collection := "Reactome"],
  as.data.table(msigdbr(db_species = "HS", species = "human", collection = "C5", subcollection = "GO:BP"))[, collection := "GO:BP"],
  # Added at author request, 22 Sep 2026: C2 chemical and genetic perturbations (exploratory collection).
  as.data.table(msigdbr(db_species = "HS", species = "human", collection = "C2", subcollection = "CGP"))[, collection := "CGP"]))
sets <- unique(sets[, .(collection, gs_name, entrez = as.character(ncbi_gene))])

run <- function(quad, bg_name, bg) {
  genes <- g[quadrant == quad, entrez_id]
  rbindlist(lapply(unique(sets$collection), function(col) {
    t2g <- sets[collection == col & entrez %in% bg, .(gs_name, entrez)]
    e <- enricher(genes, universe = bg, TERM2GENE = t2g, minGSSize = 10, maxGSSize = 500,
                  pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
    if (is.null(e) || !nrow(e@result)) return(NULL)
    r <- as.data.table(e@result)
    r[, `:=`(collection = col, quadrant = quad, background = bg_name, query_n = length(genes), background_n = length(bg))]
    r[, gene_symbols := vapply(strsplit(geneID, "/"), function(ids) paste(g$gene_name[match(ids, g$entrez_id)], collapse = "/"), "")]
    r
  }))
}
res <- rbindlist(list(
  run("UP_OPEN", "universe", g$entrez_id),
  run("DOWN_CLOSED", "universe", g$entrez_id),
  run("UP_CLOSED", "universe", g$entrez_id),
  run("DOWN_OPEN", "universe", g$entrez_id),
  run("UP_OPEN", "RNA_up_conditional", g[rna_logFC > MIN_RNA, entrez_id]),
  run("DOWN_CLOSED", "RNA_down_conditional", g[rna_logFC < -MIN_RNA, entrez_id])), fill = TRUE)
setnames(res, c("p.adjust", "pvalue"), c("BH_FDR", "PValue"))
res[, c("ratio_query", "ratio_bg") := .(sapply(GeneRatio, function(x) eval(parse(text = x))),
                                       sapply(BgRatio, function(x) eval(parse(text = x))))]
res[, fold_enrichment := ratio_query / ratio_bg]
res <- res[order(background, quadrant, collection, PValue)]
fwrite(res[, .(background, quadrant, collection, gs_name = ID, Count, query_n, GeneRatio, BgRatio, fold_enrichment,
               PValue, BH_FDR, gene_symbols)], file.path(out, "ORA_all.csv"))
sig <- res[BH_FDR < 0.05]
fwrite(sig[, .(background, quadrant, collection, gs_name = ID, Count, GeneRatio, BgRatio, fold_enrichment, BH_FDR,
               gene_symbols)], file.path(out, "ORA_FDR005.csv"))
summ <- res[, .(tested_sets = .N, FDR005 = sum(BH_FDR < .05), FDR010 = sum(BH_FDR < .10)),
            by = .(background, quadrant, collection)][order(background, quadrant, collection)]
fwrite(summ, file.path(out, "ORA_summary.csv"))
writeLines(c(sprintf("Quadrants: sign concordance with |RNA log2FC| > %.2f and |median promoter ATAC log2FC| > %.2f; no FDR filter.", MIN_RNA, MIN_ATAC),
             "Universe: 13,300 genes with RNA expression and >=1 tested promoter ATAC peak (thesis RNA branch; ID + symbol join).",
             "ORA: clusterProfiler::enricher (hypergeometric), MSigDB Hallmark / Reactome / GO:BP / C2:CGP (msigdbr 26.1.1), set size 10-500 within background, BH within collection x quadrant x background.",
             "Conditional background: genes with the same RNA direction and |RNA log2FC| above threshold; tests enrichment attributable to promoter agreement beyond the RNA change.",
             "Exploratory: effect thresholds are not significance tests; gene sets share genes; promoter peaks are group-level; no causal inference."),
           file.path(out, "METHODS.txt"))
write_json(list(MIN_RNA = MIN_RNA, MIN_ATAC = MIN_ATAC, counts = setNames(as.list(counts$N), counts$quadrant)),
           file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
print(summ)
top <- res[BH_FDR < 0.10][order(background, quadrant, BH_FDR)][, head(.SD, 8), by = .(background, quadrant)]
print(top[, .(background, quadrant, collection, ID = substr(ID, 1, 55), Count, fe = round(fold_enrichment, 2),
              FDR = signif(BH_FDR, 2))], nrows = 200)
