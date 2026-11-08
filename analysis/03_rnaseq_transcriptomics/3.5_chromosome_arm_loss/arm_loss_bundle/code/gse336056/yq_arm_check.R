#!/usr/bin/env Rscript
# GSE336056: did sgYq (chrYq double-strand break) produce Yq loss?
# Per staged sample, per chrY gene ordered by hg38 coordinate: fraction of
# cells (>= 10k UMI) with ZERO UMI for that gene, and mean UMI per 10k.
# Compare sgYq vs same-lane sgNTC (L1 = D5, L3 = D10). A breakpoint shows as
# Yp genes unchanged and distal Yq genes dropping in sgYq only.
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
ROOT <- Sys.getenv("LOY_SCRNA_ROOT", "/global/home/hpc5173/projects/LOY/scRNA")
STAGE <- file.path(ROOT, "work", "GSE336056_staged"); OUT <- file.path(ROOT, "analysis", "GSE336056")
# hg38 MSY genes, ordered by start (Mb); Yp then Yq11
PANEL <- c(RPS4Y1=2.84, ZFY=2.94, TBL1Y=6.91, PCDH11Y=4.87, AMELY=6.87, TSPY1=9.47,
           USP9Y=12.70, DDX3Y=12.90, UTY=13.25, TMSB4Y=13.70, NLGN4Y=14.52,
           KDM5D=19.70, EIF1AY=20.58, RPS4Y2=20.76)
PANEL <- sort(PANEL)
arm <- ifelse(PANEL < 10.4, "Yp", "Yq")
files <- list.files(STAGE, pattern="\\.rds$", full.names=TRUE)
res <- list(); cellstats <- list()
for (f in files) {
  sid <- sub("\\.rds$", "", basename(f)); m <- readRDS(f)
  n <- Matrix::colSums(m); keep <- n >= 10000; m <- m[, keep, drop=FALSE]; n <- n[keep]
  g <- intersect(names(PANEL), rownames(m))
  sub <- m[g, , drop=FALSE]
  yp <- Matrix::colSums(m[intersect(names(PANEL)[arm=="Yp"], rownames(m)), , drop=FALSE])
  yq <- Matrix::colSums(m[intersect(names(PANEL)[arm=="Yq"], rownames(m)), , drop=FALSE])
  res[[sid]] <- data.table(sample=sid, gene=g, mb=PANEL[g], arm=arm[match(g, names(PANEL))],
    pct_zero=round(100*Matrix::rowMeans(sub==0),1),
    umi_per10k=round(1e4*Matrix::rowSums(sub)/sum(n),3))
  # per-cell arm-level: Yq-null while Yp retained (Yp >= 5 UMI)
  cellstats[[sid]] <- data.table(sample=sid, n_cells=ncol(m), median_umi=as.integer(median(n)),
    pct_Yq0_given_Yp5=round(100*mean(yq[yp>=5]==0),2), n_Yp5=sum(yp>=5),
    pct_Yp0_given_Yq5=round(100*mean(yp[yq>=5]==0),2),
    pct_bothY0=round(100*mean(yp==0 & yq==0),2),
    mean_Yq_per10k=round(1e4*sum(yq)/sum(n),2), mean_Yp_per10k=round(1e4*sum(yp)/sum(n),2))
}
R <- rbindlist(res); C <- rbindlist(cellstats)
C[, c("armgrp","day","hto") := tstrsplit(sample, "_")]
setorder(C, day, armgrp, hto)
fwrite(R, file.path(OUT, "GSE336056_yq_arm_gene_table.tsv"), sep="\t")
fwrite(C, file.path(OUT, "GSE336056_yq_arm_cell_table.tsv"), sep="\t")
cat("\n=== per-cell arm-level (cells >=10k UMI) ===\n"); print(C, nrows=100)
cat("\n=== pct cells with ZERO UMI per gene (hg38 order), by sample ===\n")
W <- dcast(R, gene + mb + arm ~ sample, value.var="pct_zero"); setorder(W, mb)
setcolorder(W, c("gene","mb","arm", sort(names(W)[-(1:3)])))
print(W, nrows=100)
cat("\n=== UMI per 10k per gene, by sample ===\n")
W2 <- dcast(R, gene + mb + arm ~ sample, value.var="umi_per10k"); setorder(W2, mb)
setcolorder(W2, c("gene","mb","arm", sort(names(W2)[-(1:3)]))); print(W2, nrows=100)
