#!/usr/bin/env Rscript
# Arm-loss meta-analysis: gene annotation table.
# ENSG -> symbol, chromosome, coordinates, chromosome ARM (from the UCSC hg38
# cytoBand acen rows, not hand-entered centromere positions).
# Arm assignment is what lets every downstream contrast drop the cis genes on
# the lost segment before ranking for GSEA -- exactly as chrY/chr3p genes were
# excluded in the GSE336056 single-cell DE.
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "/global/home/hpc5173/projects/LOY/RNA/arm_loss_meta"
GTF  <- "/global/home/hpc5173/projects/LOY/scRNA/metadata/references/gencode.v49.annotation.gtf.gz"
CYTO <- file.path(ROOT, "data/cytoBand.txt.gz")
OUT  <- file.path(ROOT, "results/gene_annotation.tsv.gz")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# ---- centromere midpoint per chromosome: end of the p-side acen band -------------
cyto <- fread(cmd = paste("zcat", shQuote(CYTO)), header = FALSE,
              col.names = c("chrom", "start", "end", "band", "stain"))
acen <- cyto[stain == "acen" & grepl("^p", band), .(cen = max(end)), by = chrom]
stopifnot(nrow(acen) >= 24)

# ---- genes from GENCODE ----------------------------------------------------------
g <- fread(cmd = paste("zcat", shQuote(GTF), "| awk -F'\\t' '$3==\"gene\"'"),
           header = FALSE, sep = "\t",
           select = c(1, 4, 5, 9), col.names = c("chrom", "start", "end", "attr"))
g[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attr)]
g[, ensg     := sub("\\..*$", "", gene_id)]
g[, symbol   := sub('.*gene_name "([^"]+)".*', "\\1", attr)]
g[, biotype  := sub('.*gene_type "([^"]+)".*', "\\1", attr)]
g[, attr := NULL]

g <- merge(g, acen, by = "chrom", all.x = TRUE)
g[, arm := fifelse(is.na(cen), NA_character_,
                   fifelse(end < cen, paste0(sub("^chr", "", chrom), "p"),
                           fifelse(start > cen, paste0(sub("^chr", "", chrom), "q"), NA_character_)))]
# genes spanning the centromere (start < cen < end) stay NA -> treated as cis for
# either arm of that chromosome, which is the conservative choice.
g[, spans_cen := !is.na(cen) & start <= cen & end >= cen]

# one row per ENSG (GENCODE PAR_Y duplicates share a symbol but differ in gene_id)
g <- unique(g, by = "ensg")
setcolorder(g, c("ensg", "gene_id", "symbol", "biotype", "chrom", "start", "end", "arm", "spans_cen"))
fwrite(g[, .(ensg, gene_id, symbol, biotype, chrom, start, end, arm, spans_cen)],
       OUT, sep = "\t")

cat(sprintf("genes: %d | with arm: %d | spanning centromere: %d\n",
            nrow(g), sum(!is.na(g$arm)), sum(g$spans_cen)))
print(g[arm %in% c("1q", "3p", "7p", "8q", "19p"), .N, by = arm][order(arm)])
cat(sprintf("whole-chromosome gene counts 10/13/18: %s\n",
            paste(g[sub("^chr", "", chrom) %in% c("10", "13", "18"), .N,
                    by = chrom][order(chrom)][, paste0(chrom, "=", N)], collapse = ", ")))
cat("wrote", OUT, "\n")
