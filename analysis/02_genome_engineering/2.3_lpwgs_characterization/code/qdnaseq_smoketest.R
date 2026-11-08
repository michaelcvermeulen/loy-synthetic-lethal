suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
})

bins <- readRDS("../t2t_bins_100kb.rds")
cat("bins is class:", class(bins), "\n")
cat("varLabels:", paste(varLabels(bins), collapse = ", "), "\n\n")

bam <- "../t2t_bam/LNCaP-sgCTRL.bam"
cat("Counting reads from", bam, "...\n")
t0 <- Sys.time()
rc <- binReadCounts(bins, bamfiles = bam, isPaired = TRUE)
cat("binReadCounts done in",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")
print(rc)
cat("Total reads counted:", sum(assayDataElement(rc, "counts"), na.rm = TRUE), "\n")
cat("First 5 chrY bins:\n")
chrY_idx <- which(featureData(rc)$chromosome == "Y")[1:5]
print(assayDataElement(rc, "counts")[chrY_idx, , drop = FALSE])
