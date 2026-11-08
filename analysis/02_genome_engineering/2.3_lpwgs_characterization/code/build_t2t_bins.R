# Build T2T-CHM13v2.0 100kb bin annotations for QDNAseq.
# Inputs:
#   - BSgenome.Hsapiens.UCSC.hs1                  (T2T genome sequence)
#   - ../mappability_T2T.gz                       (uniquely-mappable regions, BED6)
#   - ../difficulty_regions_T2T.gz                (problematic regions, BED3)
# Output:
#   - ../t2t_bins_100kb.rds                       (AnnotatedDataFrame for QDNAseq)
#
# The QDNAseq bin annotation needs columns:
#   chromosome, start, end, bases, gc, mappability, blacklist, residual, use

suppressPackageStartupMessages({
  library(QDNAseq)
  library(BSgenome.Hsapiens.UCSC.hs1)
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

bsg     <- BSgenome.Hsapiens.UCSC.hs1
binSize <- 100  # kb

# Standard chromosomes (autosomes + X + Y; drop chrM and unplaced contigs)
chroms <- paste0("chr", c(1:22, "X", "Y"))
chr_lengths <- seqlengths(bsg)[chroms]
cat("Chromosome lengths (T2T-CHM13v2.0):\n")
print(chr_lengths)

# --- 1. Build the bin grid ----
bin_size_bp <- binSize * 1000L
bin_list <- lapply(chroms, function(chr) {
  L <- chr_lengths[chr]
  starts <- seq(1L, L, by = bin_size_bp)
  ends   <- pmin(starts + bin_size_bp - 1L, L)
  data.frame(chromosome = sub("^chr", "", chr),
             start = starts, end = ends,
             stringsAsFactors = FALSE)
})
bins <- do.call(rbind, bin_list)
cat("Total bins:", nrow(bins), "\n")

# --- 2. GC content + N count per bin ----
cat("Computing per-bin GC content from BSgenome...\n")
bins_gr <- GRanges(seqnames = paste0("chr", bins$chromosome),
                   ranges = IRanges(bins$start, bins$end))
seqs <- getSeq(bsg, bins_gr)
af   <- alphabetFrequency(seqs, baseOnly = TRUE)
non_n_bases <- rowSums(af[, c("A", "C", "G", "T")])
gc_count    <- rowSums(af[, c("G", "C")])
bins$bases <- 100 * non_n_bases / width(bins_gr)  # % of bin that is non-N
bins$gc    <- 100 * gc_count / pmax(non_n_bases, 1)

# --- 3. Mappability per bin (fraction of bin covered by uniquely-mappable regions)
cat("Loading mappability track...\n")
map_df <- read.table("../mappability_T2T.gz", header = FALSE, sep = "\t",
                     col.names = c("chrom", "start", "end", "name", "score", "strand"),
                     comment.char = "#", stringsAsFactors = FALSE)
map_gr <- GRanges(seqnames = map_df$chrom,
                  ranges = IRanges(map_df$start + 1L, map_df$end))  # BED is 0-based
map_gr <- map_gr[seqnames(map_gr) %in% chroms]
seqlevels(map_gr) <- chroms

# Reduce mappable regions and compute coverage per bin
map_red <- reduce(map_gr)
ov_map <- findOverlaps(bins_gr, map_red)
overlaps_per_bin <- tapply(
  width(pintersect(bins_gr[queryHits(ov_map)], map_red[subjectHits(ov_map)])),
  queryHits(ov_map), sum
)
bins$mappability <- 0
bins$mappability[as.integer(names(overlaps_per_bin))] <-
  100 * as.numeric(overlaps_per_bin) / width(bins_gr[as.integer(names(overlaps_per_bin))])

# --- 4. Blacklist (difficulty regions) per bin: % of bin covered by problematic regions
cat("Loading difficulty regions...\n")
diff_df <- read.table("../difficulty_regions_T2T.gz", header = FALSE, sep = "\t",
                      col.names = c("chrom", "start", "end"),
                      comment.char = "#", stringsAsFactors = FALSE)
diff_gr <- GRanges(seqnames = diff_df$chrom,
                   ranges = IRanges(diff_df$start + 1L, diff_df$end))
diff_gr <- diff_gr[seqnames(diff_gr) %in% chroms]
seqlevels(diff_gr) <- chroms

diff_red <- reduce(diff_gr)
ov_bl <- findOverlaps(bins_gr, diff_red)
bl_per_bin <- tapply(
  width(pintersect(bins_gr[queryHits(ov_bl)], diff_red[subjectHits(ov_bl)])),
  queryHits(ov_bl), sum
)
bins$blacklist <- 0
bins$blacklist[as.integer(names(bl_per_bin))] <-
  100 * as.numeric(bl_per_bin) / width(bins_gr[as.integer(names(bl_per_bin))])

# --- 5. residual + use flags (QDNAseq fills residual after fitting; default NA / TRUE)
bins$residual <- NA_real_
bins$use <- bins$bases > 90 & bins$mappability >= 50 & bins$blacklist < 10

cat("Bins flagged 'use':", sum(bins$use), "/", nrow(bins),
    "(", round(100 * sum(bins$use) / nrow(bins), 1), "%)\n")

# --- 6. Save as AnnotatedDataFrame (QDNAseq format)
rownames(bins) <- paste0(bins$chromosome, ":", bins$start, "-", bins$end)
metadata <- data.frame(labelDescription = c(
  "Chromosome name",
  "Bin start position (1-based)",
  "Bin end position",
  "Percentage of non-N bases",
  "Average GC content (% of non-N bases)",
  "Average mappability (% of bin uniquely mappable)",
  "Percent of bin overlapping blacklisted regions",
  "Median residual after GC/mappability fit (filled in by correctBins)",
  "Logical flag — TRUE for bins to use in normalization"
), row.names = colnames(bins))
bin_adf <- new("AnnotatedDataFrame", data = bins, varMetadata = metadata)

# QDNAseqBinAnnotations requires assembly metadata
attr(bin_adf, "QDNAseq") <- list(
  author      = "michael.vermeulen",
  date        = format(Sys.time()),
  organism    = "Hsapiens",
  build       = "hs1",
  version     = packageVersion("QDNAseq"),
  url         = NA_character_,
  md5         = NA_character_,
  sessionInfo = sessionInfo()
)

saveRDS(bin_adf, "../t2t_bins_100kb.rds")
cat("Wrote ../t2t_bins_100kb.rds\n")
