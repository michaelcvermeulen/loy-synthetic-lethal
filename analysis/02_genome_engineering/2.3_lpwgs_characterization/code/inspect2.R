suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
})
cache_dir <- "coverage_report_LNCaP_QDNAseq_cache/html"
e1 <- new.env(); e2 <- new.env()
lazyLoad(sub("\\.RData$", "",
             list.files(cache_dir, pattern = "run_qdnaseq.*RData", full.names = TRUE)),
         envir = e1)
lazyLoad(sub("\\.RData$", "",
             list.files(cache_dir, pattern = "correct_segment.*RData", full.names = TRUE)),
         envir = e2)

cat("== run_qdnaseq cache ==\n")
cat("Objects:", paste(ls(e1), collapse = ", "), "\n")
if ("rc" %in% ls(e1)) {
  rc <- get("rc", envir = e1)
  cat("rc dim:", nrow(rc), "x", ncol(rc), "\n")
  cat("pData:\n")
  print(pData(rc))
  cat("fData$use TRUE:", sum(fData(rc)$use), "\n")
  cat("\nReads counted per chromosome (sample 1):\n")
  cnts <- assayDataElement(rc, "counts")
  print(tapply(cnts[, 1], fData(rc)$chromosome, sum))
}

cat("\n\n== correct_segment cache ==\n")
cat("Objects:", paste(ls(e2), collapse = ", "), "\n")
if ("cn" %in% ls(e2)) {
  cn <- get("cn", envir = e2)
  cat("cn dim:", nrow(cn), "x", ncol(cn), "\n")
  segmat <- assayDataElement(cn, "segmented")
  cat("Segmented values non-NA per chromosome:\n")
  print(tapply(!is.na(segmat[, 1]), fData(cn)$chromosome, sum))
}
