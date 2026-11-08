# Transcript-level TSSs for the 13,300-gene promoter integration universe (HT1080 ATAC TSS-group figure).
# The pipeline's config/tss_knownGene.bed has one TSS per gene (gene-span start = most upstream transcript start);
# here every distinct knownGene transcript start of each universe gene is kept, so the most accessible
# annotated TSS can be chosen downstream (prepare_ht1080_atac_tss_groups.py). chrY excluded.
# Run (R 4.5, from ATACseq/SEPT_08_2026_REPROCESS to avoid the project renv): Rscript --vanilla <absolute path>
.libPaths(c("C:/Users/micha/AppData/Local/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(data.table); library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg38.knownGene)})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
fig <- dirname(dirname(normalizePath(script, winslash = "/")))
run <- file.path(dirname(dirname(dirname(fig))), "ATACseq/SEPT_08_2026_REPROCESS")
out <- file.path(fig, "source_data/ht1080_atac_accessibility/tss_transcript_sites.csv")
g <- fread(file.path(run, "integration/thesis_RNA_SEPT11/RNA_ATAC_promoter_genes.csv"), colClasses = list(character = "entrez_id"))
tx <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
tx <- tx[names(tx) %in% g$entrez_id]
d <- as.data.table(unlist(tx))
d[, entrez_id := rep(names(tx), lengths(tx))]
d[, tss := fifelse(strand == "+", start, end)]
d <- unique(d[seqnames != "chrY", .(entrez_id, chrom = as.character(seqnames), tss, strand = as.character(strand))])
fwrite(d, out)
cat(nrow(d), "distinct transcript TSSs for", uniqueN(d$entrez_id), "of", nrow(g), "universe genes\n")
