#!/usr/bin/env Rscript
# Same p53 / proliferation set-scores as chr3p_scan.R, applied to Yq-null vs Yq-intact
# within sgYq, so the two engineered losses can be compared on one metric.
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
ROOT <- Sys.getenv("LOY_SCRNA_ROOT", "/global/home/hpc5173/projects/LOY/scRNA")
STAGE <- file.path(ROOT, "work/GSE336056_staged")
SEL_YP <- c("RPS4Y1","ZFY","PCDH11Y","AMELY","TBL1Y","TSPY1"); SEL_YQ <- c("USP9Y","DDX3Y","UTY","TMSB4Y","NLGN4Y","KDM5D","EIF1AY","RPS4Y2")
P53 <- c("RPS27L","BAX","GDF15","CDKN1A","TRIM22","ZMAT3","MIR34AHG","MDM2","PHLDA3","DDB2","RRM2B","PVT1")
PROL <- c("TUBA1B","HMGB1","HSP90AA1","PTMA","STAG1","HSPD1","HNRNPD","SMC3","HNRNPF","RAN","LRR1","PSMB1")
out <- list()
for (f in list.files(STAGE, "^sgYq_.*\\.rds$", full.names = TRUE)) {
  sid <- sub("\\.rds$", "", basename(f)); m <- readRDS(f); n <- Matrix::colSums(m); k <- n >= 10000; m <- m[, k]; n <- n[k]
  yp <- Matrix::colSums(m[SEL_YP[SEL_YP %in% rownames(m)], ]); yq <- Matrix::colSums(m[SEL_YQ[SEL_YQ %in% rownames(m)], ])
  st <- fifelse(yp >= 5 & yq == 0, "Yq_null", fifelse(yp >= 5 & yq >= 2, "Yq_intact", NA_character_))
  sc <- function(gs) log2(1 + 1e4 * Matrix::colSums(m[intersect(gs, rownames(m)), ]) / n)
  out[[sid]] <- data.table(sample = sid, day = tstrsplit(sid, "_")[[2]], status = st, p53 = sc(P53), prolif = sc(PROL))[!is.na(st)]
}
D <- rbindlist(out)
print(D[, .(n = .N, p53 = round(mean(p53), 3), prolif = round(mean(prolif), 3)), by = .(day, status)][order(day, status)])
for (d in c("D5", "D10")) for (v in c("p53", "prolif")) {
  a <- D[day == d & status == "Yq_null", get(v)]; b <- D[day == d & status == "Yq_intact", get(v)]
  cat(sprintf("%s %s: Yq-null vs intact  delta = %+.3f  Wilcoxon P = %.2e\n", d, v, mean(a) - mean(b), wilcox.test(a, b)$p.value))
}
