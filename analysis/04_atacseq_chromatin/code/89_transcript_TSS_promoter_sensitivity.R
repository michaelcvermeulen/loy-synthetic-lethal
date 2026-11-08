# Sensitivity of promoter-based results to the TSS definition (23 Sep 2026).
#
# Primary annotation (06_annotation.R): promoter = gene-level TSS +/- 2 kb, where the gene-level TSS is the
# genes(TxDb knownGene) span start, i.e. the most upstream transcript start. For genes with downstream
# alternative promoters the active promoter can fall outside that window.
# Sensitivity: promoter = +/- 2 kb around every distinct knownGene transcript start of the gene (same genes,
# same chromosomes, same peaks, same differential results). Gene body and intergenic classes as in 06.
# Recomputed: (1) F06d genomic-feature composition and Fisher tests; (2) the F07 RNA-ATAC integration numbers
# (same thesis RNA source and ID + symbol-fallback join as 80_integrate_thesis_RNA.R; ATAC promoter peaks
# changed at FDR < 0.10 for F07, and at 0.05 for the thesis threshold).
# Run (R 4.5, from this run folder to avoid the project renv): Rscript --vanilla <absolute path to this script>
options(warn = 1)
.libPaths(c("C:/Users/micha/AppData/Local/R/win-library/4.5", .libPaths()))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(AnnotationDbi); library(org.Hs.eg.db); library(jsonlite)})
out <- "integration/thesis_RNA_SEPT11/transcript_TSS_sensitivity"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

genes <- fread("config/genes_knownGene.tsv", colClasses = list(character = "gene_id"))
res <- fread("differential/primary_TMM_all.csv")[primary_test_region == TRUE]
stopifnot(nrow(res) == 84510)
peaks <- GRanges(res$chrom, IRanges(res$start + 1L, res$end))
gene_body <- GRanges(genes$chrom, IRanges(genes$start, genes$end))

# Promoter windows under both definitions (genes restricted to those in genes_knownGene.tsv)
prom_gene <- GRanges(genes$chrom, IRanges(pmax(1L, genes$tss - 2000L), genes$tss + 2000L), gene_id = genes$gene_id)
tx <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
tx <- tx[names(tx) %in% genes$gene_id]
td <- as.data.table(unlist(tx))
td[, gene_id := rep(names(tx), lengths(tx))]
td[, tss := fifelse(strand == "+", start, end)]
td <- unique(td[as.character(seqnames) %in% unique(genes$chrom), .(gene_id, chrom = as.character(seqnames), tss)])
prom_tx <- GRanges(td$chrom, IRanges(pmax(1L, td$tss - 2000L), td$tss + 2000L), gene_id = td$gene_id)
cat(nrow(td), "transcript TSSs for", uniqueN(td$gene_id), "genes\n")

annotate <- function(prom) fifelse(overlapsAny(peaks, prom, ignore.strand = TRUE), "promoter",
                                   fifelse(overlapsAny(peaks, gene_body, ignore.strand = TRUE), "gene_body", "intergenic"))
res[, annotation_gene_TSS := annotate(prom_gene)]
res[, annotation_transcript_TSS := annotate(prom_tx)]
fwrite(res[, .N, by = .(annotation_gene_TSS, annotation_transcript_TSS)][order(-N)], file.path(out, "peak_annotation_crosstab.csv"))

# (1) F06d: each direction vs tested ns peaks, two-sided Fisher, BH across 6 tests per definition
feat <- rbindlist(lapply(c("annotation_gene_TSS", "annotation_transcript_TSS"), function(col) {
  r <- rbindlist(lapply(c("lost", "gained"), function(dir) rbindlist(lapply(c("promoter", "gene_body", "intergenic"), function(f) {
    a <- res[direction == dir & get(col) == f, .N]; b <- res[direction == dir, .N] - a
    c_ <- res[direction == "ns" & get(col) == f, .N]; d <- res[direction == "ns", .N] - c_
    ft <- fisher.test(matrix(c(a, c_, b, d), 2))
    data.table(definition = col, direction = dir, feature = f, foreground_in = a, foreground_total = a + b,
               foreground_pct = 100 * a / (a + b), background_pct = 100 * c_ / (c_ + d),
               odds_ratio = unname(ft$estimate), fisher_p = ft$p.value)
  }))))
  r[, fisher_BH := p.adjust(fisher_p, "BH")]
}))
fwrite(feat, file.path(out, "feature_enrichment_fisher_both_definitions.csv"))

# (2) Gene-level promoter summaries and the thesis RNA join (as 80_integrate_thesis_RNA.R)
gene_summary <- function(prom) {
  ov <- findOverlaps(peaks, prom, ignore.strand = TRUE)
  l <- unique(data.table(peak = queryHits(ov), gene_id = mcols(prom)$gene_id[subjectHits(ov)]))
  l[, `:=`(logFC = res$logFC[peak], FDR = res$FDR[peak], chrom = res$chrom[peak])]
  l[, .(tested_promoter_peaks = .N, median_peak_logFC = median(logFC),
        lost_005 = sum(FDR < .05 & logFC < 0), gained_005 = sum(FDR < .05 & logFC > 0),
        lost_010 = sum(FDR < .10 & logFC < 0), gained_010 = sum(FDR < .10 & logFC > 0),
        atac_chrom = paste(sort(unique(chrom)), collapse = ";")), by = gene_id]
}
rna <- fread(normalizePath("../../RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1/tables/gene_consistency.csv",
                           winslash = "/", mustWork = TRUE))
stopifnot(nrow(rna) == 25526)
rna[, ensembl := sub("\\..*$", "", gene_id)]
rna[, rna_class := fifelse(directionally_consistent & logFC > 0, "up", fifelse(directionally_consistent & logFC < 0, "down", "not_DE"))]
genes[, symbol := unname(mapIds(org.Hs.eg.db, keys = gene_id, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"))]
map <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = unique(rna$ensembl), keytype = "ENSEMBL", columns = "ENTREZID"))[!is.na(ENTREZID)]
map[, n_entrez := uniqueN(ENTREZID), by = ENSEMBL][, n_ensembl := uniqueN(ENSEMBL), by = ENTREZID]
map1 <- map[n_entrez == 1 & n_ensembl == 1, .(ensembl = ENSEMBL, entrez = ENTREZID)]

join <- function(atac) {
  atac <- merge(atac, genes[, .(gene_id, symbol)], by = "gene_id")
  setnames(atac, "gene_id", "entrez")
  j <- merge(merge(rna, map1, by = "ensembl"), atac, by = "entrez")
  amb <- atac[, .N, by = symbol][N > 1, symbol]
  fb <- rna[!ensembl %in% map1$ensembl & !gene_name %in% rna[duplicated(gene_name), gene_name]]
  fb <- merge(fb, atac[!symbol %in% amb & !entrez %in% j$entrez], by.x = "gene_name", by.y = "symbol")
  j <- rbind(j, fb, use.names = TRUE, fill = TRUE)
  stopifnot(!anyDuplicated(j$entrez))
  j[chrom == atac_chrom]
}
metrics <- function(j, label) {
  cls <- function(l, g) fifelse(l > 0 & g > 0, "mixed", fifelse(g > 0, "gained", fifelse(l > 0, "lost", "ns")))
  j[, `:=`(a010 = cls(lost_010, gained_010), a005 = cls(lost_005, gained_005))]
  conc <- function(a) j[(rna_class == "down" & a == "lost") | (rna_class == "up" & a == "gained"), .N]
  opp <- function(a) j[(rna_class == "down" & a == "gained") | (rna_class == "up" & a == "lost"), .N]
  data.table(definition = label,
    universe = nrow(j), DE_in_universe = j[rna_class != "not_DE", .N],
    down_in_universe = j[rna_class == "down", .N], up_in_universe = j[rna_class == "up", .N],
    down_with_lost_010 = j[rna_class == "down" & lost_010 > 0, .N],
    pct_down_with_lost_010 = 100 * j[rna_class == "down", mean(lost_010 > 0)],
    notDE_with_lost_010 = j[rna_class == "not_DE" & lost_010 > 0, .N],
    pct_notDE_with_lost_010 = 100 * j[rna_class == "not_DE", mean(lost_010 > 0)],
    up_with_gained_010 = j[rna_class == "up" & gained_010 > 0, .N],
    pct_up_with_gained_010 = 100 * j[rna_class == "up", mean(gained_010 > 0)],
    joint_010 = j[rna_class != "not_DE" & a010 != "ns", .N], concordant_010 = conc(j$a010), opposite_010 = opp(j$a010),
    mixed_010 = j[rna_class != "not_DE" & a010 == "mixed", .N],
    joint_010_chr19 = j[rna_class != "not_DE" & a010 != "ns" & chrom == "chr19", .N],
    joint_005 = j[rna_class != "not_DE" & a005 != "ns", .N], concordant_005 = conc(j$a005),
    spearman_all = cor(j$median_peak_logFC, j$logFC, method = "spearman"))
}
jg <- join(gene_summary(prom_gene))
jt <- join(gene_summary(prom_tx))
m <- rbind(metrics(jg, "gene_TSS (primary)"), metrics(jt, "transcript_TSS (sensitivity)"))
fwrite(m, file.path(out, "integration_metrics_both_definitions.csv"))
jointg <- jg[rna_class != "not_DE" & a010 != "ns", .(entrez, gene_name, rna_class, a010)]
jointt <- jt[rna_class != "not_DE" & a010 != "ns", .(entrez, gene_name, rna_class, a010)]
changes <- merge(jointg, jointt, by = c("entrez", "gene_name", "rna_class"), all = TRUE, suffixes = c("_gene_TSS", "_transcript_TSS"))
fwrite(changes[is.na(a010_gene_TSS) | is.na(a010_transcript_TSS) | a010_gene_TSS != a010_transcript_TSS], file.path(out, "joint_gene_changes_FDR010.csv"))
fwrite(jt[, .(entrez_id = entrez, gene_name, chrom, rna_logFC = logFC, rna_class, tested_promoter_peaks, atac_median_logFC = median_peak_logFC,
              lost_005, gained_005, lost_010, gained_010, in_primary_universe = entrez %in% jg$entrez)],
       file.path(out, "RNA_ATAC_promoter_genes_transcript_TSS.csv"))
print(res[, .N, by = .(annotation_gene_TSS, annotation_transcript_TSS)][order(-N)])
print(feat[, .(definition, direction, feature, n = foreground_in, fg = round(foreground_pct, 1), bg = round(background_pct, 1),
               OR = round(odds_ratio, 2), BH = signif(fisher_BH, 2))])
print(t(m))
cat("Joint-gene changes at FDR 0.10:\n"); print(changes[is.na(a010_gene_TSS) | is.na(a010_transcript_TSS) | a010_gene_TSS != a010_transcript_TSS])
