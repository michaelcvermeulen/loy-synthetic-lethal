# Run from the projects directory. Preserve upstream saved QDNAseq values.
root <- "SEPT_11_2026_THESIS/figures/chapter3/source_data"
for (cell in c("a549", "du145")) {
  src <- if (cell == "a549") "WGS/A549_QDNAseq_cn_df.rds" else "WGS/LNCaP_QDNAseq_cn_df.rds"
  keep <- if (cell == "a549") c("A549-CTRL2", "A549-LOY") else c("LNCaP-sgCTRL", "LNCaP-LOY-2X")
  x <- readRDS(src)
  x <- x[x$sample %in% keep, ]
  stopifnot(setequal(unique(x$sample), keep))
  write.csv(x, file.path(root, paste0(cell, "_lpwgs_confirmation"), "qdnaseq_bins.csv"), row.names=FALSE)
}
