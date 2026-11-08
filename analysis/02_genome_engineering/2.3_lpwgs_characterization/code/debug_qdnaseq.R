# Reproduce: count reads from one BAM with our bin file, then trace what
# changes after applyFilters
suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
})

bins <- readRDS("../t2t_bins_100kb.rds")
cat("BIN FILE: use=TRUE in", sum(bins@data$use), "bins\n")
cat("    chromY use=TRUE:", sum(bins@data$use[bins@data$chromosome == "Y"]), "\n")
cat("    chrom1 use=TRUE:", sum(bins@data$use[bins@data$chromosome == "1"]), "\n\n")

rc <- binReadCounts(bins, bamfiles = "../t2t_bam/LNCaP-sgCTRL.bam",
                    isPaired = TRUE, verbose = FALSE)
cat("\nAFTER binReadCounts:\n")
cat("  pData(rc):\n")
print(pData(rc))

# Check bins$use AFTER binReadCounts
cat("\n  fData(rc)$use TRUE:", sum(fData(rc)$use), "\n")
cat("  Per chromosome use=TRUE:\n")
print(table(fData(rc)$chromosome[fData(rc)$use]))

# Direct count of reads in use=TRUE bins from the count matrix
counts <- assayDataElement(rc, "counts")
cat("\n  Direct: total reads counted into bins:", sum(counts), "\n")
cat("  Direct: reads in fData$use bins:", sum(counts[fData(rc)$use, ]), "\n")
cat("  Per chromosome reads in counts matrix:\n")
chrom_totals <- tapply(counts[, 1], fData(rc)$chromosome, sum)
print(chrom_totals)
