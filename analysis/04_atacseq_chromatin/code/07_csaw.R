options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({
  library(csaw); library(edgeR); library(data.table); library(GenomicRanges); library(rtracklayer)
})
samples <- fread("config/samples.tsv", data.table = FALSE)
samples$condition <- factor(samples$condition, levels = c("CTRL", "LOY"))
bams <- normalizePath(paste0("bam/", samples$sample, ".bam"), winslash = "/")
blacklist <- import("config/hg38-blacklist.v2.bed", format = "BED")
par <- import("config/hg38_PAR.bed", format = "BED")
discard <- reduce(c(blacklist, par))
max_frag <- max(vapply(samples$sample, function(s) {
  max(fread(paste0("qc/", s, "_fragment_lengths.tsv"))$length)
}, numeric(1)))
param <- readParam(restrict = paste0("chr", c(1:22, "X")), pe = "both",
                   max.frag = max_frag, discard = discard)
cache <- "counts/csaw_filtered_windows.rds"
if (!file.exists(cache)) {
  message("Counting 150 bp windows at 50 bp spacing; maximum observed fragment = ", max_frag)
  windows <- windowCounts(bams, width = 150, spacing = 50, filter = 10, param = param)
  message("Counting 10 kb genome-wide background bins")
  bins <- windowCounts(bams, bin = TRUE, width = 10000, param = param)
  filter <- filterWindowsGlobal(windows, bins)
  keep <- filter$filter > log2(3)
  write.csv(data.frame(candidate_windows = length(keep), enriched_windows = sum(keep),
                       width = 150, spacing = 50, enrichment_fold = 3, max_fragment = max_frag),
            "qc/csaw_filter_summary.csv", row.names = FALSE)
  windows <- windows[keep, ]
  colnames(windows) <- samples$sample
  saveRDS(windows, cache)
} else windows <- readRDS(cache)
stopifnot(identical(colnames(windows), samples$sample), nrow(windows) > 1000)
# Use the same effective sizes as the peak-based background arm for a controlled comparison.
nf <- fread("differential/normalization_factors.csv", data.table = FALSE)
nf <- nf[nf$arm == "background_TMM", ]
nf <- nf[match(samples$sample, nf$sample), ]
y <- asDGEList(windows)
y$samples$lib.size <- nf$effective_size
y$samples$norm.factors <- 1
design <- model.matrix(~ condition, samples)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "conditionLOY")
merged <- mergeWindows(rowRanges(windows), tol = 100, max.width = 5000)
combined <- combineTests(merged$ids, test$table)
best <- getBestTest(merged$ids, test$table)
stopifnot(nrow(combined) == length(merged$regions), nrow(best) == nrow(combined))
effect <- if ("rep.logFC" %in% names(best)) best$rep.logFC else best$logFC
stopifnot(length(effect) == nrow(combined))
out <- data.frame(chrom = as.character(seqnames(merged$regions)),
  start = start(merged$regions)-1L, end = end(merged$regions),
  n_windows = combined$num.tests, logFC_best_window = effect,
  PValue = combined$PValue, FDR = combined$FDR)
out$direction <- ifelse(out$FDR < 0.05, ifelse(effect > 0, "gained", "lost"), "ns")
write.csv(out, "differential/csaw_regions.csv", row.names = FALSE)
for (arm in c("primary_TMM", "background_TMM")) {
  p <- fread(paste0("differential/", arm, "_all.csv"), data.table = FALSE)
  pg <- GRanges(p$chrom, IRanges(p$start+1L, p$end))
  sig <- out[out$FDR < 0.05, ]
  cg <- GRanges(sig$chrom, IRanges(sig$start+1L, sig$end))
  hits <- findOverlaps(pg, cg, ignore.strand = TRUE)
  concordant <- queryHits(hits)[sign(p$logFC[queryHits(hits)]) == sign(sig$logFC_best_window[subjectHits(hits)])]
  p$csaw_significant_same_direction <- seq_len(nrow(p)) %in% concordant
  write.csv(p, paste0("differential/", arm, "_csaw_support.csv"), row.names = FALSE)
  write.csv(p[p$FDR < 0.05 & p$csaw_significant_same_direction, ],
            paste0("differential/", arm, "_corroborated.csv"), row.names = FALSE)
}
writeLines(c("csaw 150 bp windows every 50 bp, >=3-fold enrichment over global 10 kb background.",
             "Window clusters: merge tolerance 100 bp, maximum width 5000 bp; combineTests region-level FDR.",
             "Representative window logFC describes direction and is not a pooled region effect estimate.",
             "Normalization matches the peak-based background sensitivity arm; compare both arms when interpreting corroboration.",
             "Shared read data mean method concordance is not independent experimental validation."), "differential/CSAW_METHODS.txt")
writeLines(capture.output(sessionInfo()), "logs/csaw_session.txt")
