# Sensitivity of HT1080 differential accessibility to an abundance-dependent (trended) bias (22 Sep 2026).
#
# Motivation: primary log2FC depend on baseline accessibility in every region class (weakest logCPM decile median
# about -0.3, strongest about +0.15), under both peak TMM and background-bin TMM. The Y-depleted library of every
# passage also has higher TSS enrichment (+1.2 to +1.6) and fewer called peaks, consistent with a signal-to-noise
# difference between conditions. Such a difference produces exactly this trend, and it is confounded with condition.
#
# Model: identical to the primary arm of 05_differential.R (same filtered peaks, ~ condition, estimateDisp robust,
# glmQLFit robust, QL F-test, TREAT |log2FC| > 0.5), except that library normalisation uses csaw::normOffsets loess
# offsets (non-linear, abundance-dependent) instead of TMM. The loess removes any abundance-dependent difference,
# including a genuine one, so this is a conservative bound rather than a replacement for the primary model.
# Also written: the passage-paired version (~ passage + condition) with the same offsets.
# Run: Rscript --vanilla <absolute path to this script>
options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
.libPaths(c("C:/Users/micha/AppData/Local/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(csaw); library(SummarizedExperiment); library(data.table)})
set.seed(20260908)
out <- "differential/trended_normalisation"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
samples <- fread("config/samples.tsv", data.table = FALSE)
samples$condition <- factor(samples$condition, levels = c("CTRL", "LOY")); samples$passage <- factor(samples$passage)
count_table <- fread("counts/peak_counts.tsv", data.table = FALSE)
keep <- fread("counts/abundance_filter.csv")
primary <- fread("differential/primary_TMM_all.csv")
counts <- as.matrix(count_table[, samples$sample]); rownames(counts) <- count_table$peak_id
counts <- counts[primary$peak_id, , drop = FALSE]
stopifnot(nrow(counts) == 84510, identical(rownames(counts), primary$peak_id))
qc <- fread("qc/qc_metrics_master.csv")
se <- SummarizedExperiment(list(counts = counts), colData = DataFrame(totals = colSums(counts)))
se$totals <- colSums(counts)
se <- normOffsets(se, se.out = TRUE)  # csaw: loess-based trended offsets (normFactors is the TMM path)
offs <- assay(se, "offset")
res_all <- list(); summ <- list()
for (arm in c("trended", "trended_passage")) {
  y <- DGEList(counts, group = samples$condition)
  y <- scaleOffset(y, offs)
  dm <- if (arm == "trended") model.matrix(~ condition, samples) else model.matrix(~ passage + condition, samples)
  ci <- match("conditionLOY", colnames(dm))
  y <- estimateDisp(y, dm, robust = TRUE)
  fit <- glmQLFit(y, dm, robust = TRUE)
  tb <- topTags(glmQLFTest(fit, coef = ci), n = Inf, sort.by = "none")$table
  tr <- topTags(glmTreat(fit, coef = ci, lfc = 0.5, null = "interval"), n = Inf, sort.by = "none")$table
  r <- data.table(peak_id = rownames(tb), logFC = tb$logFC, logCPM = tb$logCPM, F = tb$F, PValue = tb$PValue, FDR = tb$FDR,
                  TREAT_FDR = tr$FDR)
  fwrite(r, file.path(out, paste0(arm, "_all.csv")))
  m <- merge(r, primary[, .(peak_id, logFC_primary = logFC, FDR_primary = FDR, logCPM_primary = logCPM)], by = "peak_id")
  m[, dec := cut(logCPM_primary, quantile(logCPM_primary, 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
  trend <- m[, .(median_logFC_primary = median(logFC_primary), median_logFC_trended = median(logFC)), by = dec][order(dec)]
  fwrite(trend, file.path(out, paste0(arm, "_logFC_by_logCPM_decile.csv")))
  for (fdr in c(.05, .10)) {
    summ[[length(summ) + 1]] <- data.table(
      arm = arm, FDR = fdr, lost = sum(r$FDR < fdr & r$logFC < 0), gained = sum(r$FDR < fdr & r$logFC > 0),
      primary_lost = sum(m$FDR_primary < fdr & m$logFC_primary < 0), primary_gained = sum(m$FDR_primary < fdr & m$logFC_primary > 0),
      primary_lost_retained = sum(m$FDR_primary < fdr & m$logFC_primary < 0 & m$FDR < fdr & m$logFC < 0),
      primary_gained_retained = sum(m$FDR_primary < fdr & m$logFC_primary > 0 & m$FDR < fdr & m$logFC > 0),
      TREAT_lost = sum(r$TREAT_FDR < fdr & r$logFC < 0), TREAT_gained = sum(r$TREAT_FDR < fdr & r$logFC > 0),
      spearman_logFC_vs_primary = cor(m$logFC, m$logFC_primary, method = "spearman"),
      spearman_logFC_vs_logCPM = cor(m$logFC, m$logCPM_primary, method = "spearman"),
      spearman_primary_logFC_vs_logCPM = cor(m$logFC_primary, m$logCPM_primary, method = "spearman"))
  }
  print(trend)
}
summ <- rbindlist(summ); fwrite(summ, file.path(out, "summary.csv")); print(summ)
tss <- merge(qc[, .(sample, condition, passage, tss_enrichment, n_peaks, frip_consensus)], data.table(sample = samples$sample), by = "sample")
fwrite(tss, file.path(out, "library_signal_to_noise.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
