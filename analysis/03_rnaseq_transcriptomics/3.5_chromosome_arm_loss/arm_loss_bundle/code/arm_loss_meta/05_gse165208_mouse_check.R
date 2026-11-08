#!/usr/bin/env Rscript
# ============================================================================
# GSE165208 (mouse) characterisation: did each guide do what the paper says?
#
# Three deposited count sets, each with its OWN sgScramble control -- the
# cutting-but-neutral control missing from every public human arm-loss dataset:
#   ssty1_ssty2   sgSsty1 x3, sgSsty2 x3, scr x3
#   vitro_sgKdm5d kdm5d_ko x3, scr x3
#   vivo_sgKdm5d  sg_kd x3, scr x3
#
# The two guide sets test DIFFERENT things and must not be conflated:
#   Ssty1/Ssty2  mouse-specific MULTICOPY Y genes -> many simultaneous cuts ->
#                the plausible route to whole chrY loss (mouse analogue of the
#                ~60-cut ARPE-19 model). Real loss drops EVERY chrY gene.
#   Kdm5d        single-copy, conserved (human KDM5D, one of the core-6).
#                A gene knockout drops its target only, not the chromosome.
#
# Readout per arm: chrY total CPM, per-gene chrY CPM, and Xist (sex sanity).
# There is no PAR dosage check here -- mouse PAR is tiny, poorly annotated, and
# not comparable to the human PAR1 panel, so chrY dosage is the only readout.
# ============================================================================
suppressPackageStartupMessages({ library(data.table) })
ROOT <- "/global/home/hpc5173/projects/LOY/RNA/arm_loss_meta"
DATA <- file.path(ROOT, "data/GSE165208")
GTF  <- "/global/home/hpc5173/data/grcm39/gencode.vM36.primary_assembly.annotation.gtf"
OUT  <- file.path(ROOT, "results")

# ---- mouse gene -> chromosome ----------------------------------------------------
g <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' %s", shQuote(GTF)),
           header = FALSE, sep = "\t", select = c(1, 9), col.names = c("chrom", "attr"))
g[, ensm   := sub("\\..*$", "", sub('.*gene_id "([^"]+)".*', "\\1", attr))]
g[, symbol := sub('.*gene_name "([^"]+)".*', "\\1", attr)]
g[, biotype := sub('.*gene_type "([^"]+)".*', "\\1", attr)]
g[, attr := NULL]
g <- unique(g, by = "ensm"); setkey(g, ensm)
cat(sprintf("mouse GTF: %d genes | chrY %d (protein-coding %d)\n", nrow(g),
            g[chrom == "chrY", .N], g[chrom == "chrY" & biotype == "protein_coding", .N]))

files <- c(ssty = "GSE165208_ssty1_ssty2_all_sample_counts.csv.gz",
           vitro_kdm5d = "GSE165208_vitro_sgKdm5d_counts.csv.gz",
           vivo_kdm5d  = "GSE165208_vivo_sgKdm5d_counts.csv.gz")
YCORE <- c("Ddx3y", "Uty", "Kdm5d", "Eif2s3y", "Zfy1", "Zfy2", "Usp9y",
           "Ssty1", "Ssty2", "Rbmy", "Sly", "Srsy")
res <- list()

for (k in names(files)) {
  M <- fread(cmd = paste("zcat", shQuote(file.path(DATA, files[[k]]))))
  ids <- sub("\\..*$", "", M[[1]]); M <- as.matrix(M[, -1]); rownames(M) <- ids
  a <- g[J(ids)]
  cpm <- t(t(M) / colSums(M)) * 1e6
  arm <- sub("_[0-9]+$", "", colnames(M))

  cat(sprintf("\n================ %s (%d genes x %d libs) ================\n", k, nrow(M), ncol(M)))
  cat("library size (M reads):", paste(round(colSums(M) / 1e6, 1), collapse = " "), "\n")
  ychr <- which(a$chrom == "chrY")
  tot <- data.table(sample = colnames(M), arm = arm,
                    chrY_cpm = round(colSums(cpm[ychr, , drop = FALSE]), 2),
                    Xist_cpm = round(cpm[match("Xist", a$symbol), ], 1))
  print(tot)

  ctrl <- grepl("^scr", arm)
  cat("\nchrY total, mean by arm (ratio vs scr):\n")
  byarm <- tot[, .(chrY = mean(chrY_cpm)), by = arm]
  byarm[, ratio := round(chrY / byarm[grepl("^scr", arm), chrY], 3)]
  print(byarm)

  # per-gene chrY: a knockout hits its target only; chromosome loss hits all
  i <- which(a$symbol %in% YCORE & a$chrom == "chrY")
  if (length(i)) {
    pg <- data.table(gene = a$symbol[i], round(as.data.table(cpm[i, , drop = FALSE]), 2))
    cat("\nper-gene chrY (CPM):\n"); print(pg)
    mg <- melt(pg, id.vars = "gene", variable.name = "sample", value.name = "cpm")
    mg[, arm := sub("_[0-9]+$", "", sample)]
    w <- dcast(mg[, .(cpm = mean(cpm)), by = .(gene, arm)], gene ~ arm, value.var = "cpm")
    sc <- grep("^scr", names(w), value = TRUE)[1]
    for (cc in setdiff(names(w), c("gene", sc)))
      w[[paste0("ratio_", cc)]] <- round(w[[cc]] / pmax(w[[sc]], 1e-6), 3)
    cat("\nper-gene ratio vs scr:\n"); print(w)
  }
  res[[k]] <- tot
}
fwrite(rbindlist(res, idcol = "dataset"), file.path(OUT, "GSE165208_chrY_check.tsv"), sep = "\t")
cat("\nwrote", file.path(OUT, "GSE165208_chrY_check.tsv"), "\n")
