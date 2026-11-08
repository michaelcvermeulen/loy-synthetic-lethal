# Gene-set analyses of HT1080 accessibility: (A) ranked promoter-accessibility GSEA compared with RNA in the same
# gene universe; (B) nearest-gene set enrichment of distal differential regions (22 Sep 2026).
#
# A. Universe: the 13,300-gene thesis RNA-ATAC integration universe (genes with thesis RNA and >= 1 tested promoter
#    peak; SEPT_11_2026_THESIS figure source data integration_genes.csv). Three rankings of the same genes:
#      ATAC_median   median log2FC of the gene's tested promoter peaks (primary; as in 81_quadrant_ORA.R)
#      ATAC_sqrtF    median over promoter peaks of sign(log2FC) * sqrt(F) (edgeR QL F; sensitivity, precision-aware)
#      ATAC_trended  median promoter log2FC after csaw loess (trended) normalisation (88_trended_normalisation_sensitivity.R).
#                    Added because primary log2FC rise with baseline accessibility (weak peaks down, strong peaks up; the
#                    Y-depleted libraries have higher TSS enrichment), which makes highly accessible housekeeping promoters
#                    rank as "opening". This ranking is the one to interpret for promoter programmes.
#      RNA           thesis RNA log2FC (limma-voom ~batch*treatment, combined contrast)
#    fgsea multilevel (eps = 0, seed 20260922 before each run), MSigDB Hallmark, Reactome and GO:BP (msigdbr 26.1.1),
#    15-500 genes in universe; BH within each collection x ranking.
# B. Distal regions = tested regions not overlapping a knownGene TSS +/- 2 kb (annotation != promoter). Each region is
#    assigned its nearest knownGene TSS (annotation/annotated_differential_peaks.csv) within 1 Mb. Targets: lost or
#    gained at primary FDR < 0.10; background: distal regions with FDR >= 0.20. Region-level two-sided Fisher test per
#    set (target regions near set genes vs background regions near set genes), requiring >= 5 distinct target genes;
#    BH within collection x arm. Gene-level sensitivity: genes near >= 1 target region vs genes near background only
#    (Fisher). Arms: lost, gained, gained excluding chr18 (chromosome-wide dosage shift; copy_number_context/).
#    Nearest-gene assignment is positional and does not identify regulatory targets.
#    Primary distal test (added after rGREAT 2.12.2 was installed, 22 Sep 2026): rGREAT::great() local mode, basal plus
#    extension domains (5 kb up, 1 kb down, extension 1 Mb) around TxDb.Hsapiens.UCSC.hg38.knownGene TSSs, background =
#    all tested distal regions (targets included), gap regions excluded; binomial test restricted to the background, BH
#    per collection x arm, requiring >= 5 target genes. The nearest-gene Fisher test is the sensitivity.
# Run: Rscript --vanilla <absolute path to this script>
options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
.libPaths(c("C:/Users/micha/AppData/Local/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(data.table); library(msigdbr); library(fgsea)})
SEED <- 20260922L
genes_path <- normalizePath("../../SEPT_11_2026_THESIS/figures/chapter3/source_data/ht1080_rna_atac_integration/integration_genes.csv")

load_sets <- function() {
  h <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "H"))
  r <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "C2", subcollection = "CP:REACTOME"))
  b <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "C5", subcollection = "GO:BP"))
  rbindlist(list(h[, .(collection = "Hallmark", gs_name, gene = as.character(ncbi_gene))],
                 r[, .(collection = "Reactome", gs_name, gene = as.character(ncbi_gene))],
                 b[, .(collection = "GO:BP", gs_name, gene = as.character(ncbi_gene))]))[, unique(.SD)]
}
sets <- load_sets()

# ---------------- A. promoter GSEA ----------------
outA <- "integration/thesis_RNA_SEPT11/promoter_accessibility_GSEA"; dir.create(outA, showWarnings = FALSE, recursive = TRUE)
g <- fread(genes_path, colClasses = list(character = "entrez_id"))
da <- fread("differential/primary_TMM_all.csv")[, .(peak_id, logFC, F)]
pk <- g[, .(peak_id = unlist(strsplit(peak_ids, ";"))), by = entrez_id]
pk <- merge(pk, da, by = "peak_id")
sq <- pk[, .(ATAC_sqrtF = median(sign(logFC) * sqrt(F))), by = entrez_id]
trd <- fread("differential/trended_normalisation/trended_all.csv")[, .(peak_id, logFC_trended = logFC)]
sq <- merge(sq, merge(pk, trd, by = "peak_id")[, .(ATAC_trended = median(logFC_trended)), by = entrez_id], by = "entrez_id")
g <- merge(g, sq, by = "entrez_id", all.x = TRUE)
stopifnot(nrow(g) == 13300, !anyNA(g$ATAC_sqrtF), !anyNA(g$ATAC_trended), !anyNA(g$atac_median_logFC), !anyNA(g$rna_logFC))
ranks <- list(ATAC_median = setNames(g$atac_median_logFC, g$entrez_id), ATAC_sqrtF = setNames(g$ATAC_sqrtF, g$entrez_id),
              ATAC_trended = setNames(g$ATAC_trended, g$entrez_id),
              RNA = setNames(g$rna_logFC, g$entrez_id))
res <- rbindlist(lapply(unique(sets$collection), function(col) {
  pw <- split(sets[collection == col & gene %in% g$entrez_id, gene], sets[collection == col & gene %in% g$entrez_id, gs_name])
  rbindlist(lapply(names(ranks), function(rk) {
    set.seed(SEED)
    f <- as.data.table(fgseaMultilevel(pw, ranks[[rk]], minSize = 15, maxSize = 500, eps = 0))
    f[, `:=`(collection = col, ranking = rk, padj = p.adjust(pval, "BH"))]
    f[, leadingEdge_n := lengths(leadingEdge)]
    f[, leadingEdge := vapply(leadingEdge, function(x) paste(g$gene_name[match(x, g$entrez_id)], collapse = ";"), "")]
    f
  }))
}))
fwrite(res, file.path(outA, "fgsea_all.csv"))
wide <- dcast(res, collection + pathway + size ~ ranking, value.var = c("NES", "padj"))
fwrite(wide, file.path(outA, "fgsea_wide.csv"))
cors <- wide[, .(n_sets = .N,
                 rho_ATACmedian_vs_RNA = cor(NES_ATAC_median, NES_RNA, method = "spearman", use = "complete.obs"),
                 rho_ATACsqrtF_vs_RNA = cor(NES_ATAC_sqrtF, NES_RNA, method = "spearman", use = "complete.obs"),
                 rho_ATAC_rankings = cor(NES_ATAC_median, NES_ATAC_sqrtF, method = "spearman", use = "complete.obs"),
                 rho_ATACtrended_vs_RNA = cor(NES_ATAC_trended, NES_RNA, method = "spearman", use = "complete.obs"),
                 rho_ATACtrended_vs_median = cor(NES_ATAC_trended, NES_ATAC_median, method = "spearman", use = "complete.obs"),
                 sig_ATAC_trended_neg = sum(padj_ATAC_trended < .05 & NES_ATAC_trended < 0, na.rm = TRUE),
                 sig_ATAC_trended_pos = sum(padj_ATAC_trended < .05 & NES_ATAC_trended > 0, na.rm = TRUE),
                 both_sig_trended_RNA_same = sum(padj_ATAC_trended < .05 & padj_RNA < .05 & sign(NES_ATAC_trended) == sign(NES_RNA), na.rm = TRUE),
                 sig_ATAC_median_neg = sum(padj_ATAC_median < .05 & NES_ATAC_median < 0, na.rm = TRUE),
                 sig_ATAC_median_pos = sum(padj_ATAC_median < .05 & NES_ATAC_median > 0, na.rm = TRUE),
                 sig_RNA_neg = sum(padj_RNA < .05 & NES_RNA < 0, na.rm = TRUE),
                 sig_RNA_pos = sum(padj_RNA < .05 & NES_RNA > 0, na.rm = TRUE),
                 both_sig_same_sign = sum(padj_ATAC_median < .05 & padj_RNA < .05 & sign(NES_ATAC_median) == sign(NES_RNA), na.rm = TRUE),
                 both_sig_opposite = sum(padj_ATAC_median < .05 & padj_RNA < .05 & sign(NES_ATAC_median) != sign(NES_RNA), na.rm = TRUE)),
             by = collection]
fwrite(cors, file.path(outA, "summary_by_collection.csv"))
THEMES <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_KRAS_SIGNALING_UP",
            "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
            "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", "GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION", "GOBP_AXON_GUIDANCE",
            "GOBP_SYNAPSE_ORGANIZATION", "REACTOME_KRAB_ZINC_FINGER_PROTEINS_ZNFS_REGULATE_TRANSCRIPTION")
fwrite(wide[pathway %in% THEMES], file.path(outA, "themes.csv"))
print(cors)
print(wide[pathway %in% THEMES, .(pathway, NES_ATAC_trended = round(NES_ATAC_trended, 2), padj_ATAC_trended = signif(padj_ATAC_trended, 2),
                                  NES_ATAC_median = round(NES_ATAC_median, 2), padj_ATAC_median = signif(padj_ATAC_median, 2),
                                  NES_ATAC_sqrtF = round(NES_ATAC_sqrtF, 2), NES_RNA = round(NES_RNA, 2), padj_RNA = signif(padj_RNA, 2))])
for (col in unique(res$collection)) {
  cat("\n==", col, "top ATAC_median ==\n")
  print(res[collection == col & ranking == "ATAC_median"][order(padj)][1:8, .(pathway, NES = round(NES, 2), padj = signif(padj, 2), size)])
  cat("
==", col, "top ATAC_trended ==
")
  print(res[collection == col & ranking == "ATAC_trended"][order(padj)][1:10, .(pathway, NES = round(NES, 2), padj = signif(padj, 2), size)])
}

# ---------------- B. distal nearest-gene enrichment ----------------
outB <- "integration/thesis_RNA_SEPT11/distal_nearest_gene_sets"; dir.create(outB, showWarnings = FALSE, recursive = TRUE)
an <- fread("annotation/annotated_differential_peaks.csv", colClasses = list(character = "nearest_entrez_id"))
d <- an[annotation != "promoter" & !is.na(nearest_entrez_id) & nearest_entrez_id != "" & abs(nearest_tss_distance_signed) <= 1e6]
d[, class := fifelse(FDR < .10 & logFC < 0, "lost", fifelse(FDR < .10 & logFC > 0, "gained", fifelse(FDR >= .20, "background", "intermediate")))]
fwrite(d[, .N, by = class], file.path(outB, "region_counts.csv"))
bg <- d[class == "background"]
arms <- list(lost = d[class == "lost"], gained = d[class == "gained"], gained_excl_chr18 = d[class == "gained" & chrom != "chr18"])
univ <- unique(c(bg$nearest_entrez_id, d[class %in% c("lost", "gained"), nearest_entrez_id]))
resB <- rbindlist(lapply(names(arms), function(arm) {
  t <- arms[[arm]]
  tg <- unique(t$nearest_entrez_id); bgg <- setdiff(unique(bg$nearest_entrez_id), tg)
  rbindlist(lapply(unique(sets$collection), function(col) {
    s <- sets[collection == col & gene %in% univ]
    s <- s[, if (.N >= 15 && .N <= 500) .SD, by = gs_name]
    out <- s[, {
      inset <- gene
      a <- sum(t$nearest_entrez_id %in% inset); b <- nrow(t) - a
      c_ <- sum(bg$nearest_entrez_id %in% inset); d_ <- nrow(bg) - c_
      ft <- fisher.test(matrix(c(a, b, c_, d_), 2))
      ga <- sum(tg %in% inset); gb <- length(tg) - ga; gc <- sum(bgg %in% inset); gd <- length(bgg) - gc
      gt <- fisher.test(matrix(c(ga, gb, gc, gd), 2))
      .(set_size = .N, target_regions = nrow(t), target_regions_in_set = a, pct_target = 100 * a / nrow(t),
        pct_background = 100 * c_ / nrow(bg), OR = unname(ft$estimate), P = ft$p.value,
        target_genes_in_set = ga, gene_OR = unname(gt$estimate), gene_P = gt$p.value,
        target_genes = paste(head(unique(t$nearest_symbol[t$nearest_entrez_id %in% inset]), 15), collapse = ";"))
    }, by = gs_name]
    out <- out[target_genes_in_set >= 5]
    out[, `:=`(collection = col, arm = arm, BH = p.adjust(P, "BH"), gene_BH = p.adjust(gene_P, "BH"))]
    out
  }))
}))
fwrite(resB, file.path(outB, "distal_set_enrichment.csv"))
fwrite(resB[, .(sig_region = sum(BH < .05 & OR > 1), sig_region_depleted = sum(BH < .05 & OR < 1),
                sig_both = sum(BH < .05 & gene_BH < .05 & OR > 1)), by = .(arm, collection)], file.path(outB, "summary.csv"))

# ---- rGREAT (primary distal test) ----
suppressPackageStartupMessages({library(rGREAT); library(GenomicRanges)})
dist_all <- an[annotation != "promoter"]
to_gr <- function(x) GRanges(x$chrom, IRanges(x$start + 1L, x$end))  # peak coordinates are 0-based half-open
bg_gr <- to_gr(dist_all)
gt_arms <- list(lost = dist_all[FDR < .10 & logFC < 0], gained = dist_all[FDR < .10 & logFC > 0],
                gained_excl_chr18 = dist_all[FDR < .10 & logFC > 0 & chrom != "chr18"])
resG <- rbindlist(lapply(names(gt_arms), function(arm) {
  rbindlist(lapply(unique(sets$collection), function(col) {
    gs <- split(sets[collection == col, gene], sets[collection == col, gs_name])
    gs <- gs[lengths(gs) >= 15 & lengths(gs) <= 500]
    obj <- great(to_gr(gt_arms[[arm]]), gs, "TxDb.Hsapiens.UCSC.hg38.knownGene", background = bg_gr, verbose = FALSE)
    tb <- as.data.table(getEnrichmentTable(obj))
    tb <- tb[observed_gene_hits >= 5]
    tb[, `:=`(arm = arm, collection = col, target_regions = nrow(gt_arms[[arm]]), BH = p.adjust(p_value, "BH"),
              BH_hyper = p.adjust(p_value_hyper, "BH"))]
    tb
  }))
}))
fwrite(resG, file.path(outB, "rGREAT_enrichment.csv"))
fwrite(resG[, .(sig_enriched = sum(BH < .05 & fold_enrichment > 1), sig_depleted = sum(BH < .05 & fold_enrichment < 1),
                sig_enriched_also_hyper = sum(BH < .05 & fold_enrichment > 1 & BH_hyper < .05)), by = .(arm, collection)],
       file.path(outB, "rGREAT_summary.csv"))
for (a_ in names(gt_arms)) for (c_ in unique(sets$collection)) {
  x <- resG[resG$arm == a_ & resG$collection == c_][order(p_value)][1:8]
  cat("\n== rGREAT", a_, c_, "==\n")
  print(x[, .(id, FE = round(fold_enrichment, 2), regions = observed_region_hits, genes = observed_gene_hits, BH = signif(BH, 2),
              BH_hyper = signif(BH_hyper, 2))])
}
for (a_ in names(arms)) for (c_ in unique(sets$collection)) {
  x <- resB[resB$arm == a_ & resB$collection == c_][order(P)][1:5]
  cat("\n== nearest-gene", a_, c_, "==\n")
  print(x[, .(gs_name, OR = round(OR, 2), BH = signif(BH, 2), gene_BH = signif(gene_BH, 2), target_genes_in_set)])
}
print(fread(file.path(outB, "rGREAT_summary.csv")))
writeLines(capture.output(sessionInfo()), file.path(outB, "sessionInfo.txt"))
