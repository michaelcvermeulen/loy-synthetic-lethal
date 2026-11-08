# Polycomb / H3K27me3-target enrichment at changing promoters, against covariate-matched backgrounds (22 Sep 2026).
#
# Motivation: 81_quadrant_ORA.R / 82_prc2_promoter_check.R found PRC2 / H3K27me3 CGP sets enriched among closing
# promoters regardless of RNA direction. Polycomb targets are CpG-island promoters that are often weakly accessible,
# so the enrichment could reflect baseline accessibility or sequence composition rather than Y depletion.
#
# Covariates per tested peak (501 bp): baseline accessibility (primary logCPM, motifs/peak_covariates.csv), GC fraction
# and CpG observed/expected (computed here from BSgenome.Hsapiens.UCSC.hg38, same reference as chromVAR).
# Strata: logCPM deciles x GC quintiles x CpG o/e quintiles (quantiles over all tested peaks).
#
# A. Gene level (thesis RNA integration universe, 13,300 genes; covariates = mean over the gene's tested promoter
#    peaks): ATAC-closed (median promoter log2FC < -T) and ATAC-open (> +T) genes vs matched ATAC-flat genes
#    (|log2FC| <= 0.25), for T = 0.25 and 0.5; all genes and RNA-flat genes (|RNA log2FC| <= 0.25).
# B. Peak level (promoter peaks = overlapping a knownGene TSS +/- 2 kb): significant lost and gained peaks at
#    FDR < 0.05 and < 0.10 vs matched non-differential promoter peaks (FDR >= 0.20). Set membership = any gene
#    whose promoter the peak overlaps.
# Matching: K = 3 controls per target, sampled without replacement within stratum (seed 20260922); targets in strata
# without enough controls are dropped and counted. Test: two-sided Cochran-Mantel-Haenszel over strata (continuity
# corrected), with Mantel-Haenszel common OR and 95% CI; BH across every test in this script.
# Run: Rscript --vanilla <absolute path to this script> [K]   (default K = 3; K = 1 is a sensitivity run)
options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(msigdbr); library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings)})
out <- "integration/thesis_RNA_SEPT11/polycomb_matched"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.integer(args[1]) else 3L  # K = 1 sensitivity retains more low-accessibility targets
SEED <- 20260922L
TAG <- if (K == 3L) "" else sprintf("_K%d", K)
# Optional second argument "trended" (23 Sep 2026): gene- and peak-level ATAC effects and FDRs from the trended (csaw loess)
# model (88_trended_normalisation_sensitivity.R). Matching covariates (baseline logCPM, GC, CpG o/e) are unchanged.
TRENDED <- length(args) >= 2 && args[2] == "trended"
if (TRENDED) TAG <- paste0(TAG, "_trended")
SETS <- c("BENPORATH_ES_WITH_H3K27ME3", "BENPORATH_PRC2_TARGETS", "BENPORATH_SUZ12_TARGETS", "BENPORATH_EED_TARGETS",
          "MEISSNER_NPC_HCP_WITH_H3K4ME2_AND_H3K27ME3", "MIKKELSEN_MEF_HCP_WITH_H3K27ME3", "MIKKELSEN_NPC_HCP_WITH_H3K27ME3")
cgp <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "C2", subcollection = "CGP"))[gs_name %in% SETS]
members <- split(as.character(cgp$ncbi_gene), cgp$gs_name)
stopifnot(all(SETS %in% names(members)))

# ---- peak covariates ----
da <- fread("differential/primary_TMM_all.csv")[, .(peak_id, chrom, start, end, logFC, FDR)]
if (TRENDED) {
  trd <- fread("differential/trended_normalisation/trended_all.csv")[, .(peak_id, logFC_t = logFC, FDR_t = FDR)]
  da <- merge(da, trd, by = "peak_id")[, `:=`(logFC = logFC_t, FDR = FDR_t)][, c("logFC_t", "FDR_t") := NULL]
  stopifnot(nrow(da) == 84510)
}
cov <- fread("motifs/peak_covariates.csv")[, .(peak_id, GC, logCPM)]
pk <- merge(da, cov, by = "peak_id")
cache <- file.path(out, "peak_CpG_oe.csv")
if (file.exists(cache)) {
  cg <- fread(cache)
} else {
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(pk$chrom, IRanges(pk$start + 1L, pk$end)))
  di <- dinucleotideFrequency(seqs)[, "CG"]; mono <- letterFrequency(seqs, c("C", "G"))
  L <- width(seqs)
  cg <- data.table(peak_id = pk$peak_id, CpG_oe = ifelse(mono[, "C"] * mono[, "G"] > 0, di * L / (mono[, "C"] * mono[, "G"]), 0))
  fwrite(cg, cache)
}
pk <- merge(pk, cg, by = "peak_id")
qbin <- function(x, k) cut(x, unique(quantile(x, seq(0, 1, length.out = k + 1), na.rm = TRUE)), include.lowest = TRUE, labels = FALSE)
pk[, `:=`(b_cpm = qbin(logCPM, 10), b_gc = qbin(GC, 5), b_cpg = qbin(CpG_oe, 5))]

links <- fread("annotation/promoter_peak_gene_links.csv", colClasses = list(character = "gene_id"))
peak_genes <- links[, .(genes = list(unique(gene_id))), by = peak_id]

# ---- matching + CMH ----
match_test <- function(dt, is_target, is_control, strata_cols, member_fun, label) {
  set.seed(SEED)
  dt <- copy(dt)[, stratum := do.call(paste, c(.SD, sep = "_")), .SDcols = strata_cols]
  tg <- dt[is_target]; ct <- dt[is_control]
  picked <- list(); dropped <- 0L
  for (s in unique(tg$stratum)) {
    t_s <- tg[stratum == s]; c_s <- ct[stratum == s]
    n <- min(nrow(t_s), floor(nrow(c_s) / K))
    dropped <- dropped + nrow(t_s) - n
    if (n == 0) next
    picked[[s]] <- rbind(t_s[sample.int(nrow(t_s), n)][, grp := "target"], c_s[sample.int(nrow(c_s), n * K)][, grp := "control"])
  }
  m <- rbindlist(picked)
  rbindlist(lapply(SETS, function(set) {
    m[, member := member_fun(.SD, set)]
    tab <- table(factor(m$grp, c("target", "control")), factor(m$member, c(TRUE, FALSE)), m$stratum)
    tab <- tab[, , apply(tab, 3, function(x) sum(x) > 1 && all(rowSums(x) > 0)), drop = FALSE]
    mh <- mantelhaen.test(tab, correct = TRUE)
    data.table(comparison = label, gene_set = set, targets = sum(m$grp == "target"), controls = sum(m$grp == "control"),
               targets_dropped = dropped, target_pct = 100 * mean(m[grp == "target", member]),
               control_pct = 100 * mean(m[grp == "control", member]), MH_OR = unname(mh$estimate),
               CI_low = mh$conf.int[1], CI_high = mh$conf.int[2], P = mh$p.value)
  }))
}
balance <- function(dt, is_target, is_control, label) {
  x <- rbind(dt[is_target, .(grp = "target", logCPM, GC, CpG_oe)], dt[is_control, .(grp = "all controls", logCPM, GC, CpG_oe)])
  x[, .(n = .N, logCPM = mean(logCPM), GC = mean(GC), CpG_oe = mean(CpG_oe)), by = grp][, comparison := label]
}

# ---- A. gene level ----
g <- fread("integration/thesis_RNA_SEPT11/RNA_ATAC_promoter_genes.csv", colClasses = list(character = "entrez_id"))
if (TRENDED) {
  g_pk <- g[, .(peak_id = unlist(strsplit(peak_ids, ";"))), by = entrez_id]
  g_tm <- merge(g_pk, da[, .(peak_id, logFC)], by = "peak_id")[, .(atac_trended = median(logFC)), by = entrez_id]
  g <- merge(g, g_tm, by = "entrez_id", all.x = TRUE); stopifnot(!anyNA(g$atac_trended))
  g[, atac_median_logFC := atac_trended]
}
gp <- merge(links[, .(peak_id, entrez_id = gene_id)], pk[, .(peak_id, logCPM, GC, CpG_oe)], by = "peak_id")
gc_ <- gp[, .(logCPM = mean(logCPM), GC = mean(GC), CpG_oe = mean(CpG_oe)), by = entrez_id]
g <- merge(g, gc_, by = "entrez_id")
ref <- pk  # bin edges from all tested peaks
brk <- function(x, ref_x, k) { b <- unique(quantile(ref_x, seq(0, 1, length.out = k + 1))); b[1] <- -Inf; b[length(b)] <- Inf; cut(x, b, labels = FALSE) }
g[, `:=`(b_cpm = brk(logCPM, ref$logCPM, 10), b_gc = brk(GC, ref$GC, 5), b_cpg = brk(CpG_oe, ref$CpG_oe, 5))]
gene_member <- function(d, set) d$entrez_id %in% members[[set]]
res <- list(); bal <- list()
for (Tt in c(0.25, 0.5)) for (scope in c("all genes", "RNA-flat genes")) {
  d <- if (scope == "all genes") g else g[abs(rna_logFC) <= 0.25]
  flat <- abs(d$atac_median_logFC) <= 0.25
  for (dir in c("closed", "open")) {
    tgt <- if (dir == "closed") d$atac_median_logFC < -Tt else d$atac_median_logFC > Tt
    lab <- sprintf("gene | ATAC %s (|log2FC| > %.2f) vs matched flat | %s", dir, Tt, scope)
    res[[lab]] <- match_test(d, tgt, flat, c("b_cpm", "b_gc", "b_cpg"), gene_member, lab)
    bal[[lab]] <- balance(d, tgt, flat, lab)
  }
}
# Unmatched reference (for comparison with the matched OR)
for (Tt in c(0.25, 0.5)) {
  tgt <- g$atac_median_logFC < -Tt; flat <- abs(g$atac_median_logFC) <= 0.25
  for (set in SETS) {
    mm <- g$entrez_id %in% members[[set]]
    f <- fisher.test(table(factor(tgt[tgt | flat], c(TRUE, FALSE)), factor(mm[tgt | flat], c(TRUE, FALSE))))
    res[[paste("unmatched", Tt, set)]] <- data.table(comparison = sprintf("gene | ATAC closed (|log2FC| > %.2f) vs ALL flat, UNMATCHED | all genes", Tt),
      gene_set = set, targets = sum(tgt), controls = sum(flat), targets_dropped = 0L,
      target_pct = 100 * mean(mm[tgt]), control_pct = 100 * mean(mm[flat]), MH_OR = unname(f$estimate),
      CI_low = f$conf.int[1], CI_high = f$conf.int[2], P = f$p.value)
  }
}

# ---- B. peak level (promoter peaks) ----
pp <- merge(pk, peak_genes, by = "peak_id")
peak_member <- function(d, set) vapply(d$genes, function(x) any(x %in% members[[set]]), logical(1))
ctrl <- pp$FDR >= 0.20
for (thr in c(0.05, 0.10)) for (dir in c("lost", "gained")) {
  tgt <- pp$FDR < thr & (if (dir == "lost") pp$logFC < 0 else pp$logFC > 0)
  lab <- sprintf("peak | %s promoter peaks FDR < %.2f vs matched non-DA promoter peaks", dir, thr)
  res[[lab]] <- match_test(pp, tgt, ctrl, c("b_cpm", "b_gc", "b_cpg"), peak_member, lab)
  bal[[lab]] <- balance(pp, tgt, ctrl, lab)
}

res <- rbindlist(res); res[, BH_FDR := p.adjust(P, "BH")]
fwrite(res, file.path(out, sprintf("polycomb_matched_tests%s.csv", TAG))); fwrite(rbindlist(bal), file.path(out, paste0("covariate_balance_before_matching", if (TRENDED) "_trended" else "", ".csv")))
writeLines(c("See header of code/83_polycomb_matched.R.",
             "Strata: logCPM deciles x GC quintiles x CpG o/e quintiles over all 84,510 tested peaks; K = 3 matched controls; seed 20260922.",
             "CMH two-sided with continuity correction; BH across all rows of polycomb_matched_tests.csv (the unmatched Fisher rows are included for reference).",
             "Peak level uses promoter peaks only; membership = any overlapping-promoter gene in the set.",
             "Accessibility and sequence matching do not measure H3K27me3; set membership comes from ES/NPC/MEF studies, not HT1080."),
           file.path(out, "METHODS.txt"))
print(rbindlist(bal)[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)], nrows = 100)
show <- res[gene_set %in% c("BENPORATH_ES_WITH_H3K27ME3", "BENPORATH_PRC2_TARGETS", "MEISSNER_NPC_HCP_WITH_H3K4ME2_AND_H3K27ME3")]
print(show[, .(comparison = substr(comparison, 1, 78), set = substr(gene_set, 1, 22), targets, dropped = targets_dropped,
               tgt = round(target_pct, 1), ctl = round(control_pct, 1), OR = round(MH_OR, 2),
               CI = sprintf("%.2f-%.2f", CI_low, CI_high), FDR = signif(BH_FDR, 2))], nrows = 200)
