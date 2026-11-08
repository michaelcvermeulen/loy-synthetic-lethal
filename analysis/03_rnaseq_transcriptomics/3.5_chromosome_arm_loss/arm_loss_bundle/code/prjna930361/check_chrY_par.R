#!/usr/bin/env Rscript
# ============================================================================
# PRJNA930361 biological check: is the A7123 clone actually Y-null, and does
# PAR dosage agree?
#
# This is the check that QC cannot do. Two independent readouts:
#
#   1. MSY (male-specific Y) expression. If chrY really was eliminated, every
#      MSY gene goes to ~0 in the Ynull samples and stays normal in parental.
#
#   2. PAR dosage. PAR1/PAR2 are present on BOTH sex chromosomes, so a male
#      carries two copies. Lose chrY and one copy goes with it, so PAR
#      expression should fall to ~0.5 of parental while autosomes stay at 1.0.
#      This is the readout the PAR-masked reference exists to protect: with an
#      unmasked chrY the PAR reads multi-map between chrX and chrY and are
#      discarded at MAPQ 0, which destroys exactly this signal.
#
# PAR is flagged BY COORDINATE, not from a gene list -- PAR_genes.csv in this
# project is missing PAR2 (VAMP7/IL9R/SPRY3). With the masked reference the PAR
# reads all land on the chrX copy, so PAR genes are chrX entries here.
#
# Usage: Rscript check_chrY_par.R [BATCH_DIR]
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })
BATCH <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else
  "/global/home/hpc5173/projects/LOY/RNA/updated/PRJNA930361"
GTF <- "/global/home/hpc5173/data/gencode.v49.primary_assembly.annotation.gtf"
PAR1_END <- 2781479; PAR2_START <- 155701383          # hg38 chrX

mdir <- c(file.path(BATCH, "matrix"), file.path(BATCH, "matrix_s0"))
mdir <- mdir[file.exists(file.path(mdir, "counts_matrix.tsv"))]
stopifnot("no counts_matrix.tsv found" = length(mdir) > 0)
cat("matrix:", file.path(mdir[1], "counts_matrix.tsv"), "\n")
M <- fread(file.path(mdir[1], "counts_matrix.tsv"))
gid <- M[[1]]; M <- as.matrix(M[, -1]); rownames(M) <- sub("\\..*$", "", gid)

# ---- gene coordinates ------------------------------------------------------------
ann <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' %s", shQuote(GTF)),
             header = FALSE, sep = "\t", select = c(1, 4, 5, 9),
             col.names = c("chrom", "start", "end", "attr"))
ann[, ensg   := sub("\\..*$", "", sub('.*gene_id "([^"]+)".*', "\\1", attr))]
ann[, symbol := sub('.*gene_name "([^"]+)".*', "\\1", attr)]
ann[, attr := NULL]
ann <- unique(ann, by = "ensg")
ann[, region := fcase(
  chrom == "chrY", "MSY",
  chrom == "chrX" & (end < PAR1_END | start > PAR2_START), "PAR",
  chrom == "chrX", "chrX",
  default = "autosome")]
a <- ann[match(rownames(M), ensg)]
cat(sprintf("genes: %d | annotated %d | MSY %d | PAR %d\n",
            nrow(M), sum(!is.na(a$chrom)), sum(a$region == "MSY", na.rm = TRUE),
            sum(a$region == "PAR", na.rm = TRUE)))

cpm <- t(t(M) / colSums(M)) * 1e6
grp <- ifelse(grepl("Ynull", colnames(M)), "Ynull", "WT")
cat("\nlibrary sizes (M reads):\n"); print(round(colSums(M) / 1e6, 1))

# ---- 1. MSY -----------------------------------------------------------------------
CORE <- c("RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "NLGN4Y", "TMSB4Y", "TXLNGY")
i <- match(CORE, a$symbol)
msy <- data.table(gene = CORE, as.data.table(round(cpm[i, , drop = FALSE], 2)))
cat("\n=== 1. MSY genes (CPM) ===\n"); print(msy)
tot <- data.table(sample = colnames(M), group = grp,
                  MSY_cpm = round(colSums(cpm[which(a$region == "MSY"), , drop = FALSE]), 1),
                  XIST_cpm = round(cpm[match("XIST", a$symbol), ], 2))
cat("\ntotal MSY + XIST per sample:\n"); print(tot)

# ---- 2. PAR dosage ------------------------------------------------------------------
# Ratio of Ynull to WT, computed per region. Expect PAR ~0.5, autosome ~1.0.
reg_cpm <- sapply(c("MSY", "PAR", "chrX", "autosome"), function(r)
  colSums(cpm[which(a$region == r), , drop = FALSE]))
R <- data.table(region = colnames(reg_cpm),
                WT    = round(colMeans(reg_cpm[grp == "WT", , drop = FALSE]), 1),
                Ynull = round(colMeans(reg_cpm[grp == "Ynull", , drop = FALSE]), 1))
R[, ratio := round(Ynull / WT, 3)]
cat("\n=== 2. Region dosage, Ynull vs parental ===\n"); print(R)

# per-gene PAR ratio, expressed genes only: a median near 0.5 is the real signal;
# near 1.0 means the Y copy was never counted (unmasked-reference failure mode)
pk <- which(a$region == "PAR" & rowMeans(cpm) > 1)
if (length(pk) >= 3) {
  pg <- data.table(gene = a$symbol[pk],
                   WT    = round(rowMeans(cpm[pk, grp == "WT", drop = FALSE]), 2),
                   Ynull = round(rowMeans(cpm[pk, grp == "Ynull", drop = FALSE]), 2))
  pg[, ratio := round(Ynull / WT, 3)]
  setorder(pg, -WT)
  cat(sprintf("\nexpressed PAR genes (n=%d), median Ynull/WT ratio = %.3f\n",
              nrow(pg), median(pg$ratio)))
  print(pg)
}
hk <- c("GAPDH", "ACTB", "B2M", "HPRT1", "TBP", "PPIA")
hi <- match(hk, a$symbol)
hr <- data.table(gene = hk,
                 WT = round(rowMeans(cpm[hi, grp == "WT", drop = FALSE]), 1),
                 Ynull = round(rowMeans(cpm[hi, grp == "Ynull", drop = FALSE]), 1))
hr[, ratio := round(Ynull / WT, 3)]
cat("\nhousekeeping control (expect ratio ~1.0):\n"); print(hr)

# ---- verdict --------------------------------------------------------------------------
msy_wt <- mean(tot[group == "WT", MSY_cpm]); msy_ko <- mean(tot[group == "Ynull", MSY_cpm])
par_ratio <- R[region == "PAR", ratio]
cat("\n=== VERDICT ===\n")
cat(sprintf("MSY: parental %.1f CPM -> Ynull %.1f CPM (%.1f%% remaining)\n",
            msy_wt, msy_ko, 100 * msy_ko / msy_wt))
cat(sprintf("PAR dosage ratio: %.3f  (expect ~0.5 for true chrY loss; ~1.0 means PAR was not measured)\n", par_ratio))
cat(sprintf("autosome ratio:   %.3f  (sanity, expect ~1.0)\n", R[region == "autosome", ratio]))
cat(if (msy_ko / msy_wt < 0.05 && par_ratio < 0.7)
      "=> chrY elimination CONFIRMED on both readouts.\n"
    else if (msy_ko / msy_wt < 0.05)
      "=> MSY gone, but PAR dosage does NOT halve - check the reference actually masked chrY PAR.\n"
    else "=> MSY is NOT depleted; sample labels or the clone are wrong.\n")
fwrite(tot, file.path(BATCH, "chrY_par_check_per_sample.tsv"), sep = "\t")
fwrite(R, file.path(BATCH, "chrY_par_check_regions.tsv"), sep = "\t")
