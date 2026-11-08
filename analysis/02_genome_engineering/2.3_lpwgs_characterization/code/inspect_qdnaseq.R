suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
})

# Load cached results
cache_dir <- "coverage_report_LNCaP_QDNAseq_cache/html"
# knitr uses lazyLoad caches (.rdb/.rdx pair)
rc_cache_base <- sub("\\.RData$", "",
                     list.files(cache_dir, pattern = "run_qdnaseq.*RData", full.names = TRUE))
seg_cache_base <- sub("\\.RData$", "",
                      list.files(cache_dir, pattern = "correct_segment.*RData", full.names = TRUE))
lazyLoad(rc_cache_base)
lazyLoad(seg_cache_base)
cat("Loaded objects:", paste(ls(), collapse = ", "), "\n\n")

cat("== rc (raw counts) ==\n")
cat("Total bins:", nrow(rc), "samples:", ncol(rc), "\n")
cat("Bins per chromosome:\n")
print(table(fData(rc)$chromosome))

cat("\n== rc_f (after applyFilters) ==\n")
cat("Total bins kept:", sum(fData(rc_f)$use), "of", nrow(rc_f), "\n")
cat("USED bins per chromosome:\n")
print(table(fData(rc_f)$chromosome[fData(rc_f)$use]))

cat("\n== cn (segmented copy numbers) ==\n")
cat("Total bins:", nrow(cn), "samples:", ncol(cn), "\n")

# Check copynumber and segmented matrices
cn_mat  <- assayDataElement(cn, "copynumber")
seg_mat <- assayDataElement(cn, "segmented")
cat("\nNon-NA copynumber values per chromosome:\n")
chroms <- fData(cn)$chromosome
for (chr in unique(chroms)) {
  n_total <- sum(chroms == chr)
  n_with_cov <- sum(!is.na(cn_mat[chroms == chr, 1]))
  n_with_seg <- sum(!is.na(seg_mat[chroms == chr, 1]))
  cat(sprintf("  chr%-3s: %d bins | covered: %d | segmented: %d\n",
              chr, n_total, n_with_cov, n_with_seg))
}

cat("\n== Sample-level summary ==\n")
print(pData(cn))

cat("\n== Per-chromosome median segmented log2 (all 8 samples) ==\n")
seg_log2 <- log2(pmax(seg_mat, 1e-3))
for (s in colnames(seg_log2)) {
  cat(s, ":\n")
  for (chr in c("1", "7", "8", "X", "Y")) {
    vals <- seg_log2[chroms == chr, s]
    cat(sprintf("  chr%-3s median log2_seg = %.3f (n=%d non-NA)\n",
                chr, median(vals, na.rm = TRUE), sum(!is.na(vals))))
  }
}
