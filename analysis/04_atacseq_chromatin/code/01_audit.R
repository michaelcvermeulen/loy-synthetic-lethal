options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE)))
setwd(root)
suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(digest); library(jsonlite)
})
samples <- expand.grid(passage = c(11L, 12L, 21L), condition = c("CTRL", "LOY"),
                       stringsAsFactors = FALSE)
samples$sample <- paste(samples$condition, samples$passage, sep = "_")
samples$replicate <- rep(1:3, 2)
samples$bam <- paste0("../filtered/", samples$sample, ".bam")
samples$index <- paste0(samples$bam, ".bai")
samples$analysis_unit <- "independent_replicate_per_user"
samples <- samples[, c("sample", "condition", "replicate", "passage", "analysis_unit", "bam", "index")]
stopifnot(all(file.exists(samples$bam)), all(file.exists(samples$index)))
write.table(samples, "config/samples.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
canonical <- paste0("chr", c(1:22, "X", "Y"))
eof_expected <- as.raw(c(31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,27,0,3,0,0,0,0,0,0,0,0,0))
target_reference <- NULL
audit <- list()
for (i in seq_len(nrow(samples))) {
  sample <- samples$sample[i]; bam <- samples$bam[i]
  message("Auditing ", sample)
  header <- scanBamHeader(bam)[[1]]
  target <- header$targets
  stopifnot(all(canonical %in% names(target)))
  stopifnot(target[["chr1"]] == 248956422, target[["chrX"]] == 156040895,
            target[["chrY"]] == 57227415)
  if (is.null(target_reference)) target_reference <- target
  stopifnot(identical(target, target_reference))
  writeLines(capture.output(str(header, give.attr = FALSE)), paste0("logs/", sample, "_header.txt"))
  header_flat <- paste(unlist(header$text), collapse = " ")
  stopifnot(grepl("SO:coordinate", header_flat), grepl("REMOVE_DUPLICATES=true", header_flat),
            grepl("-q 30", header_flat, fixed = TRUE))
  con <- file(bam, "rb"); seek(con, where = file.info(bam)$size - 28, origin = "start")
  eof <- readBin(con, "raw", n = 28); close(con)
  stopifnot(identical(eof, eof_expected))
  idx <- idxstatsBam(bam)
  write.table(idx, paste0("qc/", sample, "_source_idxstats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  bf <- BamFile(bam, yieldSize = 100000L); open(bf)
  chunk <- scanBam(bf, param = ScanBamParam(what = c("flag", "mapq", "isize", "rname", "qname")))[[1]]
  close(bf)
  flag <- chunk$flag
  audit[[sample]] <- data.frame(sample = sample, bam_bytes = file.info(bam)$size,
    sha256 = digest(bam, algo = "sha256", file = TRUE),
    bai_sha256 = digest(samples$index[i], algo = "sha256", file = TRUE),
    bgzf_eof_present = TRUE, source_mapped_reads = sum(idx$mapped),
    canonical_mapped_reads = sum(idx$mapped[idx$seqnames %in% canonical]),
    source_chrM_reads = sum(idx$mapped[idx$seqnames %in% c("chrM", "MT")]),
    sampled_records = length(flag), sampled_min_mapq = min(chunk$mapq),
    sampled_duplicate_reads = sum(bitwAnd(flag, 1024L) != 0L),
    sampled_supplementary_reads = sum(bitwAnd(flag, 2048L) != 0L),
    sampled_not_proper_pair = sum(bitwAnd(flag, 2L) == 0L))
}
write.csv(do.call(rbind, audit), "qc/source_bam_audit.csv", row.names = FALSE)
sizes <- data.frame(chrom = names(target_reference), length = unname(target_reference))
write.table(sizes, "config/bam.chrom.sizes", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene <- genes(txdb)
gene <- gene[as.character(seqnames(gene)) %in% canonical]
stopifnot(all(seqlengths(gene)[canonical] == target_reference[canonical]))
tss <- ifelse(as.character(strand(gene)) == "+", start(gene), end(gene))
gene_table <- data.frame(gene_id = names(gene), chrom = as.character(seqnames(gene)),
  start = start(gene), end = end(gene), strand = as.character(strand(gene)), tss = tss)
write.table(gene_table, "config/genes_knownGene.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
tss_bed <- unique(data.frame(chrom = gene_table$chrom, start = tss - 1L,
                             end = tss, gene_id = gene_table$gene_id, score = 0,
                             strand = gene_table$strand))
write.table(tss_bed, "config/tss_knownGene.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
blacklist_source <- "../AUG_27_reprocess/config/hg38-blacklist.v2.bed.gz"
stopifnot(file.exists(blacklist_source))
blacklist <- read.table(gzfile(blacklist_source), sep = "\t", header = FALSE)
stopifnot(ncol(blacklist) >= 3, all(blacklist[[2]] >= 0), all(blacklist[[3]] > blacklist[[2]]))
file.copy(blacklist_source, "config/hg38-blacklist.v2.bed.gz", overwrite = TRUE)
write.table(blacklist, "config/hg38-blacklist.v2.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
refs <- c("config/bam.chrom.sizes", "config/genes_knownGene.tsv", "config/tss_knownGene.bed", "config/hg38-blacklist.v2.bed.gz")
write.csv(data.frame(path = refs, sha256 = vapply(refs, digest, "", algo = "sha256", file = TRUE)),
          "config/reference_checksums.csv", row.names = FALSE)
write_json(list(genome = "GRCh38", canonical_chromosomes = canonical,
  primary_chromosomes = paste0("chr", c(1:22, "X")), primary_design = "~ condition",
  reference_condition = "CTRL", contrast = "LOY-CTRL", min_mapq_both_mates = 30,
  peak_caller = "MACS3", peak_format = "BAMPE", peak_qvalue = 0.01,
  effective_genome_size = 2913022398, peak_width = 501, min_support_per_condition = 2,
  primary_fdr = 0.05, seed = 20260908, threads = 8,
  tss_annotation = "TxDb.Hsapiens.UCSC.hg38.knownGene 3.22.0",
  tss_score = "mean insertion density +/-50 bp / mean density in outer 100 bp of +/-2000 bp window; duplicate TSS positions collapsed",
  blacklist_source = blacklist_source, external_mitochondrial_fraction = "not yet recovered",
  normalization = "peak TMM primary; background-bin sensitivity",
  rna_wgs_matching = "pending"), "config/params.json", pretty = TRUE, auto_unbox = TRUE)
writeLines(c("Source BAM audit passed: dictionaries, index readability, BGZF EOF, expected processing history.",
             "Header and 100,000-record sample checks are not full-stream integrity or mate validation.",
             "Full-stream read-pair validation will occur during preparation.",
             "Original chrM contamination is not recoverable from these filtered BAMs."), "qc/source_bam_audit_notes.txt")
message("Source audit and reference preparation complete.")
