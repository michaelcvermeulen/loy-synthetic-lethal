# Export the retained tidy QDNAseq object without recalculating normalization.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript export_ht1080_lpwgs_bins.R input.rds output.csv")
d <- readRDS(args[1])
stopifnot(all(c("bin_id","sample","log2","log2_seg","chromosome","midpoint","use","label","genome_x") %in% names(d)))
dir.create(dirname(args[2]), recursive=TRUE, showWarnings=FALSE)
write.csv(d, args[2], row.names=FALSE)
