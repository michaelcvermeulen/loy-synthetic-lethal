options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({
  library(edgeR); library(data.table); library(jsonlite); library(ggplot2); library(pheatmap)
})
set.seed(20260908)
samples <- fread("config/samples.tsv", data.table = FALSE)
samples$condition <- factor(samples$condition, levels = c("CTRL", "LOY"))
samples$passage <- factor(samples$passage)
peaks <- fread("peaks/consensus_membership.tsv", data.table = FALSE)
count_table <- fread("counts/peak_counts.tsv", data.table = FALSE)
stopifnot(identical(peaks$peak_id, count_table$peak_id),
          identical(names(count_table)[-1], samples$sample))
counts <- as.matrix(count_table[, samples$sample]); rownames(counts) <- count_table$peak_id
stopifnot(all(counts >= 0), all(counts == round(counts)))
primary_region <- peaks$primary_test_region %in% c(TRUE, "True", "TRUE")
stopifnot(!any(peaks$chrom[primary_region] == "chrY"))
design <- model.matrix(~ condition, samples)
design_passage <- model.matrix(~ passage + condition, samples)
stopifnot(qr(design)$rank == ncol(design), qr(design_passage)$rank == ncol(design_passage))
qc <- do.call(rbind, lapply(samples$sample, function(s) {
  x <- fromJSON(paste0("qc/", s, "_preparation.json"))
  data.frame(sample = s, accepted_fragments = x$accepted_fragments,
             nonY_fragments = x$nonY_fragments, chrY_fragments = x$chrY_fragments,
             tss_enrichment = x$tss_enrichment, nfr_fraction = x$nfr_fraction,
             unmatched_primary_records = x$unmatched_primary_records,
             blacklist_pairs = x$blacklist_pairs)
}))
peak_qc <- fread("qc/peak_qc.tsv", data.table = FALSE)
qc <- cbind(samples[, c("sample", "condition", "passage")], qc[, -1],
            peak_qc[match(samples$sample, peak_qc$sample), c("n_peaks", "frip_own", "frip_consensus", "extra_consensus_assignments")])
qc$chrY_per_million_nonY <- qc$chrY_fragments/qc$nonY_fragments*1e6
qc$original_mitochondrial_fraction <- NA_real_
qc$raw_alignment_rate <- NA_real_
qc$source_pre_dedup_duplicate_rate <- vapply(samples$sample, function(s) {
  lines <- readLines(paste0("../filtered/", s, "_dedup_metrics.txt"), warn = FALSE)
  header <- grep("^LIBRARY\\t", lines)[1]
  keys <- strsplit(lines[header], "\t")[[1]]
  values <- strsplit(lines[header+1L], "\t")[[1]]
  as.numeric(values[match("PERCENT_DUPLICATION", keys)])
}, numeric(1))
write.csv(qc, "qc/qc_metrics_master.csv", row.names = FALSE)
write.csv(design, "config/design_primary.csv")
write.csv(design_passage, "config/design_passage_sensitivity.csv")
y0 <- DGEList(counts[primary_region, , drop = FALSE], group = samples$condition)
keep <- filterByExpr(y0, design = design)
stopifnot(sum(keep) > 1000)
filtered_counts <- y0$counts[keep, , drop = FALSE]
tested_peaks <- peaks[match(rownames(filtered_counts), peaks$peak_id), ]
write.csv(data.frame(peak_id = rownames(y0), retained = keep), "counts/abundance_filter.csv", row.names = FALSE)
background <- fread("counts/background_bin_counts.tsv", data.table = FALSE)
bg_counts <- as.matrix(background[, samples$sample])
bg <- DGEList(bg_counts, lib.size = qc$nonY_fragments)
bg <- bg[rowSums(bg_counts) >= 30, , keep.lib.sizes = TRUE]
stopifnot(nrow(bg) > 1000)
bg <- calcNormFactors(bg, method = "TMM")
bg_effective_sizes <- with(bg$samples, lib.size * norm.factors)
normalization_rows <- list(); results_list <- list(); fit_list <- list()
for (arm in c("primary_TMM", "background_TMM", "passage_TMM")) {
  message("Fitting ", arm)
  y <- DGEList(filtered_counts, group = samples$condition)
  if (arm == "background_TMM") {
    y$samples$lib.size <- bg_effective_sizes
    y$samples$norm.factors <- 1
  } else {
    y <- calcNormFactors(y, method = "TMM")
  }
  dm <- if (arm == "passage_TMM") design_passage else design
  coef_index <- match("conditionLOY", colnames(dm))
  y <- estimateDisp(y, dm, robust = TRUE)
  fit <- glmQLFit(y, dm, robust = TRUE)
  test <- glmQLFTest(fit, coef = coef_index)
  table <- topTags(test, n = Inf, sort.by = "none")$table
  treat_test <- glmTreat(fit, coef = coef_index, lfc = 0.5, null = "interval")
  treat_table <- topTags(treat_test, n = Inf, sort.by = "none")$table
  stopifnot(identical(rownames(table), tested_peaks$peak_id))
  out <- cbind(tested_peaks, table)
  out$TREAT_PValue_lfc0.5 <- treat_table$PValue
  out$TREAT_FDR_lfc0.5 <- treat_table$FDR
  out$direction <- ifelse(out$FDR < 0.05, ifelse(out$logFC > 0, "gained", "lost"), "ns")
  write.csv(out, paste0("differential/", arm, "_all.csv"), row.names = FALSE)
  write.csv(out[out$FDR < 0.05, ], paste0("differential/", arm, "_significant.csv"), row.names = FALSE)
  normalization_rows[[arm]] <- data.frame(arm = arm, sample = samples$sample,
    library_size = y$samples$lib.size, factor = y$samples$norm.factors,
    effective_size = y$samples$lib.size*y$samples$norm.factors)
  results_list[[arm]] <- out
  fit_list[[arm]] <- list(y = y, fit = fit, design = dm)
  if (arm == "primary_TMM") {
    lcpm <- cpm(y, log = TRUE, prior.count = 2)
    write.table(data.frame(peak_id = rownames(lcpm), lcpm), "counts/logCPM_primary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    pdf("differential/primary_model_diagnostics.pdf", width = 8, height = 6)
    plotBCV(y); plotQLDisp(fit); plotMD(test, status = decideTests(test), main = "LOY versus CTRL")
    abline(h = 0, col = "grey50"); dev.off()
    pca <- prcomp(t(lcpm), scale. = FALSE)
    variance <- 100*pca$sdev^2/sum(pca$sdev^2)
    pc <- cbind(samples, as.data.frame(pca$x))
    write.csv(pc, "qc/pca_scores.csv", row.names = FALSE)
    p <- ggplot(pc, aes(PC1, PC2, colour = condition, label = sample)) +
      geom_point(size = 3) + geom_text(nudge_y = 2, show.legend = FALSE) +
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.12)) +
      scale_colour_manual(values = c(CTRL = "#0072B2", LOY = "#D55E00")) +
      theme_bw(base_size = 12) + labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]))
    ggsave("qc/pca.png", p, width = 8, height = 6, dpi = 200)
    correlation <- cor(lcpm, method = "spearman")
    write.csv(correlation, "qc/spearman_correlations.csv")
    pheatmap(correlation, display_numbers = TRUE, filename = "qc/correlation_heatmap.png", width = 7, height = 6)
    mds <- plotMDS(y, plot = FALSE)
    mds_frame <- data.frame(sample = samples$sample, condition = samples$condition,
                            dimension1 = mds$x, dimension2 = mds$y)
    write.csv(mds_frame, "qc/mds_scores.csv", row.names = FALSE)
    p <- ggplot(mds_frame, aes(dimension1, dimension2, colour = condition, label = sample)) +
      geom_point(size = 3) + geom_text(vjust = -0.8, show.legend = FALSE) +
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.12)) +
      scale_colour_manual(values = c(CTRL = "#0072B2", LOY = "#D55E00")) +
      theme_bw(base_size = 12) + labs(x = "MDS dimension 1", y = "MDS dimension 2")
    ggsave("qc/mds.png", p, width = 8, height = 6, dpi = 200)
  }
}
write.csv(do.call(rbind, normalization_rows), "differential/normalization_factors.csv", row.names = FALSE)
saveRDS(fit_list, "differential/models.rds")
primary <- results_list$primary_TMM
colours <- c(gained = "#D55E00", lost = "#0072B2", ns = "grey75")
p <- ggplot(primary, aes(logFC, -log10(pmax(FDR, .Machine$double.xmin)), colour = direction)) +
  geom_point(size = 0.6, alpha = 0.5) + scale_colour_manual(values = colours) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  theme_bw(base_size = 12) + labs(x = "log2 fold change (LOY / CTRL)", y = "-log10 FDR")
ggsave("differential/primary_volcano.png", p, width = 8, height = 6, dpi = 200)
summary <- do.call(rbind, lapply(names(results_list), function(arm) {
  x <- results_list[[arm]]
  data.frame(arm = arm, tested = nrow(x), gained = sum(x$direction == "gained"),
             lost = sum(x$direction == "lost"), treat_fdr005_lfc05 = sum(x$TREAT_FDR_lfc0.5 < 0.05),
             spearman_logFC_vs_primary = cor(primary$logFC, x$logFC, method = "spearman"),
             shared_significant = sum(primary$FDR < 0.05 & x$FDR < 0.05),
             shared_significant_same_direction = sum(primary$FDR < 0.05 & x$FDR < 0.05 & sign(primary$logFC) == sign(x$logFC)))
}))
write.csv(summary, "differential/analysis_summary.csv", row.names = FALSE)
print(summary)
profiles <- do.call(rbind, lapply(samples$sample, function(s) {
  p <- fread(paste0("qc/", s, "_tss_profile.tsv"), data.table = FALSE)
  p$sample <- s
  p$enrichment <- p$insertions/mean(c(head(p$insertions, 100), tail(p$insertions, 100)))
  p$bin <- floor((p$offset+2000)/20)
  aggregate(cbind(offset, enrichment) ~ sample + bin, p, mean)
}))
p <- ggplot(profiles, aes(offset, enrichment, colour = sample)) + geom_line(linewidth = 0.6) +
  theme_bw(base_size = 12) + labs(x = "Distance from TSS (bp)", y = "Insertion enrichment over flanks")
ggsave("qc/tss_profiles.png", p, width = 9, height = 6, dpi = 200)
fragments <- do.call(rbind, lapply(samples$sample, function(s) {
  f <- fread(paste0("qc/", s, "_fragment_lengths.tsv"), data.table = FALSE)
  f$sample <- s; f$fraction <- f$fragments/sum(f$fragments); f
}))
p <- ggplot(fragments[fragments$length <= 1000, ], aes(length, fraction, colour = sample)) +
  geom_line(linewidth = 0.5) + theme_bw(base_size = 12) + labs(x = "Fragment length (bp)", y = "Fraction of eligible fragments")
ggsave("qc/fragment_lengths.png", p, width = 9, height = 6, dpi = 200)
p <- ggplot(qc, aes(sample, chrY_per_million_nonY, fill = condition)) +
  geom_col() + scale_fill_manual(values = c(CTRL = "#0072B2", LOY = "#D55E00")) +
  theme_bw(base_size = 12) + labs(x = NULL, y = "chrY fragments per million non-Y fragments")
ggsave("qc/chrY_signal.png", p, width = 8, height = 5, dpi = 200)
writeLines(capture.output(sessionInfo()), "logs/differential_session.txt")
