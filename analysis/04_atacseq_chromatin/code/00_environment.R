options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
dir.create("logs", showWarnings = FALSE)
packages <- c("Rsamtools", "Rsubread", "edgeR", "csaw", "rtracklayer",
              "ATACseqQC", "GenomicRanges", "GenomicAlignments", "GenomicFeatures",
              "TxDb.Hsapiens.UCSC.hg38.knownGene", "BSgenome.Hsapiens.UCSC.hg38",
              "motifmatchr", "JASPAR2022", "clusterProfiler", "ChIPseeker",
              "data.table", "ggplot2", "pheatmap", "jsonlite", "digest")
installed <- installed.packages()
versions <- data.frame(package = packages,
                       version = installed[match(packages, rownames(installed)), "Version"])
write.csv(versions, "logs/r_packages.csv", row.names = FALSE)
print(versions, row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "", .libPaths()), "logs/r_session.txt")
