script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38); library(TFBSTools)
  library(JASPAR2022); library(motifmatchr); library(GenomicRanges)
})
stopifnot(requireNamespace("statmod", quietly = TRUE), requireNamespace("org.Hs.eg.db", quietly = TRUE))
sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges("chr1", IRanges(1000000, 1000500)))
stopifnot(width(sequence) == 501L)
matrices <- getMatrixSet(JASPAR2022, opts = list(species = 9606, collection = "CORE", all_versions = FALSE))
saveRDS(matrices, "config/JASPAR2022_human_CORE.rds")
write.csv(data.frame(motif_id = names(matrices), motif_name = vapply(matrices, name, character(1))),
          "config/JASPAR2022_human_CORE_manifest.csv", row.names = FALSE)
matches <- motifMatches(matchMotifs(matrices[1:2], sequence, out = "matches", p.cutoff = 5e-5))
stopifnot(identical(dim(matches), c(1L, 2L)))
writeLines(c(paste("Human CORE motifs:", length(matrices)),
             "Local reference sequence access and motif matching passed.",
             capture.output(sessionInfo())), "logs/motif_dependency_check.txt")
cat("Local sequence/motif dependency check passed;", length(matrices), "human CORE motifs available.\n")
