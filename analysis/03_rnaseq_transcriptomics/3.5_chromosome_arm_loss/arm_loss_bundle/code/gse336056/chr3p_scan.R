#!/usr/bin/env Rscript
# GSE336056: what happened on chr3p in the sg3p arm? Expression-dosage scan
# (inferCNV-lite) on staged matrices, cells >= 10k UMI, same-day sgNTC as baseline.
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
ROOT <- Sys.getenv("LOY_SCRNA_ROOT", "/global/home/hpc5173/projects/LOY/scRNA")
STAGE <- file.path(ROOT, "work/GSE336056_staged"); OUT <- file.path(ROOT, "analysis/GSE336056")
source(file.path(ROOT, "scripts/resolve_hg38_msy_genes.R"), local = TRUE)
g <- as.data.table(load_gencode_genes(file.path(ROOT, "metadata/references/gencode.v49.annotation.gtf.gz")))
g <- unique(g[, .(gene = gene_name, chrom = seqname, start, end)], by = "gene")
CEN <- c(chr3 = 90.5e6, chr5 = 48.8e6, chr14 = 17.2e6, chr1 = 123.4e6)   # hg38 centromere midpoints (approx)

mats <- list(); md <- list()
for (f in list.files(STAGE, "\\.rds$", full.names = TRUE)) {
  sid <- sub("\\.rds$", "", basename(f)); m <- readRDS(f)
  n <- Matrix::colSums(m); keep <- n >= 10000; m <- m[, keep]; n <- n[keep]
  p <- tstrsplit(sid, "_"); mats[[sid]] <- m
  md[[sid]] <- data.table(cell = colnames(m), sample = sid, arm = p[[1]], day = p[[2]], nUMI = n)
}
M <- do.call(cbind, mats); D <- rbindlist(md)
det <- Matrix::rowMeans(M > 0); genes <- names(det)[det >= 0.20]
gg <- g[gene %in% genes]
seg <- function(chr, side) gg[chrom == chr & (if (side == "p") end < CEN[[chr]] else start > CEN[[chr]]), gene]
SEGS <- list(chr3p = seg("chr3","p"), chr3q = seg("chr3","q"), chr5p = seg("chr5","p"), chr5q = seg("chr5","q"),
             chr14q = seg("chr14","q"), chr1p = seg("chr1","p"), chr1q = seg("chr1","q"))
cat("genes per segment:", paste(names(SEGS), sapply(SEGS, length)), "\n")
per10k <- function(gs) 1e4 * Matrix::colSums(M[gs, , drop = FALSE]) / D$nUMI
for (s in names(SEGS)) D[[s]] <- per10k(SEGS[[s]])
# log2 ratio to same-day sgNTC median
for (s in names(SEGS)) {
  base <- D[arm == "sgNTC", .(b = median(get(s))), by = day]
  D[[paste0(s, "_l2r")]] <- log2(D[[s]] / base$b[match(D$day, base$day)])
}
# 3p-loss call: 3p ratio below log2(0.8) while 3q stays within +-0.25
D[, loss3p := chr3p_l2r < log2(0.8) & abs(chr3q_l2r) < 0.25]
cat("\n=== per-sample summary (cells >= 10k UMI) ===\n")
S <- D[, .(n = .N, median_3p_l2r = round(median(chr3p_l2r), 3), median_3q_l2r = round(median(chr3q_l2r), 3),
           pct_3p_loss = round(100 * mean(loss3p), 1), n_3p_loss = sum(loss3p),
           median_5q_l2r = round(median(chr5q_l2r), 3), median_14q_l2r = round(median(chr14q_l2r), 3)), by = .(day, arm, sample)]
print(S[order(day, arm)], nrows = 50)
cat("\n=== 3p ratio distribution quantiles, by arm x day ===\n")
print(D[, as.list(round(quantile(chr3p_l2r, c(.01,.05,.1,.25,.5,.75,.9)), 2)), by = .(day, arm)][order(day, arm)])

# breakpoint scan along chr3 (5 Mb bins): 3p-loss sg3p cells vs same-day sgNTC
loss <- D[arm == "sg3p" & loss3p == TRUE, cell]; ctrl <- D[arm == "sgNTC", cell]
c3 <- gg[chrom == "chr3"][order(start)]; c3[, bin := floor(start / 5e6) * 5]
bins <- c3[, .(genes = list(gene), n = .N), by = bin][n >= 5]
B <- rbindlist(lapply(seq_len(nrow(bins)), function(i) {
  gs <- bins$genes[[i]]
  x <- 1e4 * Matrix::colSums(M[gs, loss, drop = FALSE]) / D[match(loss, cell), nUMI]
  y <- 1e4 * Matrix::colSums(M[gs, ctrl, drop = FALSE]) / D[match(ctrl, cell), nUMI]
  data.table(bin_mb = bins$bin[i], n_genes = bins$n[i], l2r_loss_vs_ntc = round(log2(median(x) / median(y)), 3))
}))
cat(sprintf("\n=== chr3 5-Mb bins: 3p-loss sg3p cells (n=%d) vs sgNTC (n=%d); centromere ~90.5 Mb ===\n", length(loss), length(ctrl)))
print(B, nrows = 60)

# same scan for the 3p-loss cells split by day, plus secondary aneuploidies 5q / 14q
cat("\n=== 5q and 14q in 3p-loss vs 3p-intact sg3p cells, by day ===\n")
print(D[arm == "sg3p", .(n = .N, med_5q = round(median(chr5q_l2r), 3), med_14q = round(median(chr14q_l2r), 3),
                          med_1p = round(median(chr1p_l2r), 3)), by = .(day, loss3p)][order(day, loss3p)])

# p53 program in 3p-loss vs 3p-intact sg3p cells (same 12 targets as the Yq DE)
P53 <- c("RPS27L","BAX","GDF15","CDKN1A","TRIM22","ZMAT3","MIR34AHG","MDM2","PHLDA3","DDB2","RRM2B","PVT1")
PROL <- c("TUBA1B","HMGB1","HSP90AA1","PTMA","STAG1","HSPD1","HNRNPD","SMC3","HNRNPF","RAN","LRR1","PSMB1")
E <- function(gs, cells) log2(1 + 1e4 * Matrix::colSums(M[intersect(gs, rownames(M)), cells, drop = FALSE]) / D[match(cells, cell), nUMI])
D[, p53_score := E(P53, cell)]; D[, prolif_score := E(PROL, cell)]
cat("\n=== p53 / proliferation scores (log2 per-10k of gene set) ===\n")
Q <- rbind(
  D[arm == "sg3p", .(group = fifelse(loss3p, "sg3p 3p-loss", "sg3p 3p-intact"), day, p53_score, prolif_score)],
  D[arm == "sgNTC", .(group = "sgNTC", day, p53_score, prolif_score)])
print(Q[, .(n = .N, p53 = round(mean(p53_score), 3), prolif = round(mean(prolif_score), 3)), by = .(day, group)][order(day, group)])
for (d in c("D5", "D10")) {
  a <- Q[day == d & group == "sg3p 3p-loss", p53_score]; b <- Q[day == d & group == "sg3p 3p-intact", p53_score]
  cat(sprintf("%s p53: 3p-loss vs 3p-intact  delta = %+.3f  Wilcoxon P = %.2e\n", d, mean(a) - mean(b), wilcox.test(a, b)$p.value))
  a <- Q[day == d & group == "sg3p 3p-loss", prolif_score]; b <- Q[day == d & group == "sg3p 3p-intact", prolif_score]
  cat(sprintf("%s prolif: 3p-loss vs 3p-intact  delta = %+.3f  Wilcoxon P = %.2e\n", d, mean(a) - mean(b), wilcox.test(a, b)$p.value))
}
fwrite(D[, .(cell, sample, arm, day, nUMI, chr3p_l2r, chr3q_l2r, chr5q_l2r, chr14q_l2r, chr1p_l2r, loss3p, p53_score, prolif_score)],
       file.path(OUT, "GSE336056_chr3p_dosage_per_cell.tsv"), sep = "\t")
fwrite(B, file.path(OUT, "GSE336056_chr3_breakpoint_bins.tsv"), sep = "\t")
cat("done\n")
