options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38); library(TFBSTools)
  library(motifmatchr); library(Matrix); library(ggplot2)
})
set.seed(20260908)
result <- fread("differential/primary_TMM_all.csv", data.table = FALSE)
regions <- GRanges(result$chrom, IRanges(result$start+1L, result$end))
genome <- BSgenome.Hsapiens.UCSC.hg38
sequence <- getSeq(genome, regions)
gc <- rowSums(letterFrequency(sequence, letters = c("G", "C"), as.prob = TRUE))
n_fraction <- as.numeric(letterFrequency(sequence, letters = "N", as.prob = TRUE))
result$GC <- gc; result$N_fraction <- n_fraction
gc_bin <- cut(gc, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)
cpm_breaks <- unique(quantile(result$logCPM, probs = seq(0, 1, by = 0.1)))
cpm_bin <- cut(result$logCPM, breaks = cpm_breaks, include.lowest = TRUE)
result$stratum <- interaction(gc_bin, cpm_bin, drop = TRUE)
write.csv(result[, c("peak_id", "GC", "N_fraction", "logCPM", "stratum")], "motifs/peak_covariates.csv", row.names = FALSE)
background_result <- fread("differential/background_TMM_all.csv", data.table = FALSE)
background_result <- background_result[match(result$peak_id, background_result$peak_id), ]
gc_qc <- rbind(data.frame(GC = gc, logFC = result$logFC, arm = "primary_TMM", N_fraction = n_fraction),
               data.frame(GC = gc, logFC = background_result$logFC, arm = "background_TMM", N_fraction = n_fraction))
gc_qc <- gc_qc[gc_qc$N_fraction == 0, ]
gc_qc$gc_bin <- floor(gc_qc$GC/0.025)
gc_curve <- aggregate(cbind(GC, logFC) ~ arm + gc_bin, gc_qc, median)
p <- ggplot(gc_curve, aes(GC, logFC, colour = arm)) + geom_line() + geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") + theme_bw(base_size = 12) +
  labs(x = "Peak GC fraction (binned)", y = "Median log2 fold change (LOY / CTRL)")
ggsave("qc/gc_effects.png", p, width = 8, height = 5, dpi = 200)
gc_statistics <- do.call(rbind, lapply(split(gc_qc, gc_qc$arm), function(x) {
  data.frame(arm = x$arm[1], spearman_GC_logFC = cor(x$GC, x$logFC, method = "spearman"), n_peaks = nrow(x))
}))
write.csv(gc_statistics, "qc/gc_effect_correlations.csv", row.names = FALSE)
panels <- list(); accounting <- list()
for (direction in c("gained", "lost")) {
  query <- which(result$direction == direction & n_fraction == 0)
  background <- which(result$FDR >= 0.2 & n_fraction == 0)
  chosen_query <- chosen_background <- integer()
  for (stratum in levels(result$stratum)) {
    q <- query[result$stratum[query] == stratum]
    b <- background[result$stratum[background] == stratum]
    n <- min(length(q), length(b))
    if (!n) next
    chosen_query <- c(chosen_query, q[sample.int(length(q), n)])
    chosen_background <- c(chosen_background, b[sample.int(length(b), n)])
  }
  accounting[[direction]] <- data.frame(direction = direction, eligible_query = length(query),
    matched_query = length(chosen_query), matched_background = length(chosen_background),
    unmatched_query = length(query)-length(chosen_query))
  if (length(chosen_query) < 20) next
  panels[[direction]] <- data.frame(index = c(chosen_query, chosen_background),
    role = rep(c("query", "background"), each = length(chosen_query)),
    direction = direction)
}
write.csv(do.call(rbind, accounting), "motifs/matching_accounting.csv", row.names = FALSE)
if (!length(panels)) {
  writeLines("No direction had at least 20 matched query peaks; motif enrichment was not estimable.", "motifs/STATUS.txt")
} else {
  selected <- sort(unique(unlist(lapply(panels, `[[`, "index"))))
  pwms <- readRDS("config/JASPAR2022_human_CORE.rds")
  matches <- motifMatches(matchMotifs(pwms, sequence[selected], out = "matches", p.cutoff = 5e-5))
  stopifnot(nrow(matches) == length(selected), ncol(matches) == length(pwms))
  output <- list()
  for (direction in names(panels)) {
    panel <- panels[[direction]]
    panel$peak_id <- result$peak_id[panel$index]
    panel$stratum <- result$stratum[panel$index]
    write.csv(panel, paste0("motifs/", direction, "_matched_peaks.csv"), row.names = FALSE)
    x <- matches[match(panel$index, selected), , drop = FALSE]
    groups <- split(seq_len(nrow(panel)), droplevels(panel$stratum))
    tests <- lapply(seq_len(ncol(x)), function(j) {
      motif <- as.logical(x[, j])
      tables <- array(0, dim = c(2, 2, length(groups)))
      for (k in seq_along(groups)) {
        ids <- groups[[k]]; is_query <- panel$role[ids] == "query"
        tables[, , k] <- matrix(c(sum(motif[ids] & is_query), sum(!motif[ids] & is_query),
                                  sum(motif[ids] & !is_query), sum(!motif[ids] & !is_query)), nrow = 2)
      }
      fit <- tryCatch(mantelhaen.test(tables, alternative = "greater", correct = FALSE), error = function(e) NULL)
      data.frame(motif_id = names(pwms)[j], motif_name = name(pwms[[j]]), direction = direction,
        query_with_motif = sum(motif & panel$role == "query"),
        background_with_motif = sum(motif & panel$role == "background"),
        odds_ratio = if (is.null(fit)) NA_real_ else unname(fit$estimate),
        PValue = if (is.null(fit) || !is.finite(fit$p.value)) 1 else fit$p.value)
    })
    output[[direction]] <- do.call(rbind, tests)
  }
  table <- do.call(rbind, output)
  table$FDR <- p.adjust(table$PValue, method = "BH")
  table <- table[order(table$FDR), ]
  write.csv(table, "motifs/motif_enrichment.csv", row.names = FALSE)
  for (direction in names(panels)) {
    top <- head(table[table$direction == direction, ], 20)
    top$label <- factor(paste(top$motif_name, top$motif_id), levels = rev(paste(top$motif_name, top$motif_id)))
    p <- ggplot(top, aes(-log10(pmax(FDR, .Machine$double.xmin)), label)) +
      geom_point(colour = if (direction == "gained") "#D55E00" else "#0072B2", size = 2) +
      theme_bw(base_size = 11) + labs(x = "-log10 FDR", y = NULL)
    ggsave(paste0("motifs/", direction, "_motifs.png"), p, width = 8, height = 7, dpi = 200)
  }
}
writeLines(c("JASPAR2022 human CORE motifs; installed database version and R packages are recorded.",
             "Fixed 501 bp peaks. Query: primary FDR<0.05 gained/lost; controls: primary FDR>=0.2.",
             "One-to-one sampling within 5 percentage-point GC bins and baseline logCPM deciles; seed 20260908.",
             "Strata lacking enough controls are subsampled on both sides; unmatched query counts are reported.",
             "N-containing sequences excluded. motifmatchr P cutoff=5e-5, subject sequence background.",
             "One-sided Cochran-Mantel-Haenszel enrichment stratified by matching bins; BH across both directions and all motifs.",
             "Sequence enrichment does not demonstrate TF binding or TF activity."), "motifs/METHODS.txt")
writeLines(capture.output(sessionInfo()), "logs/motifs_session.txt")
