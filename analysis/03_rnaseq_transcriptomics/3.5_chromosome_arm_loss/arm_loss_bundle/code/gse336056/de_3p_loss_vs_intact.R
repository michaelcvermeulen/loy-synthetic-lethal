#!/usr/bin/env Rscript
# ============================================================================
# GSE336056: MAST DE, 3p-loss vs 3p-intact cells WITHIN the sg3p arm.
# Mirror of de_yq_null_vs_intact.R so the two engineered arm losses can be
# compared pathway-for-pathway.
#
# Cells: sg3p singlets, >= 10k UMI. Status = loss3p call from chr3p_scan.R
# (GSE336056_chr3p_dosage_per_cell.tsv: 3p l2r < log2(0.8) & |3q l2r| < 0.25
# vs same-day sgNTC). Model identical to the Yq DE:
#   log2(1 + 1e4 * counts / nUMI); detection >= 10% in either group;
#   zlm(~ status + cngeneson + sample, bayesglm, ebayes); hurdle LRT on status.
# Fits: pooled (primary) + D5-only + D10-only (sign-consistency check).
# Gene classes as in the Yq DE plus chr3p (chr3, end < 90.5 Mb): 3p genes are
# half-dosage by construction (and define the status call), so they are
# EXCLUDED from GSEA, exactly like chrY genes in the Yq DE.
# Outputs -> analysis/GSE336056/de/
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(data.table); library(MAST); library(SummarizedExperiment) })
ROOT <- Sys.getenv("LOY_SCRNA_ROOT", "/global/home/hpc5173/projects/LOY/scRNA")
STAGE <- file.path(ROOT, "work/GSE336056_staged"); OUT <- file.path(ROOT, "analysis/GSE336056/de")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
logf <- function(fmt, ...) cat(sprintf(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", fmt, "\n"), ...))
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")); options(mc.cores = NC)
MIN_UMI <- 10000; MIN_PCT <- 0.10
SEL_YP <- c("RPS4Y1","ZFY","PCDH11Y","AMELY","TBL1Y","TSPY1")
SEL_YQ <- c("USP9Y","DDX3Y","UTY","TMSB4Y","NLGN4Y","KDM5D","EIF1AY","RPS4Y2")

# ---- gene classes from GENCODE v49 --------------------------------------------------
source(file.path(ROOT, "scripts/resolve_hg38_msy_genes.R"), local = TRUE)
gtf <- file.path(ROOT, "metadata/references/gencode.v49.annotation.gtf.gz")
genes_all <- as.data.table(load_gencode_genes(gtf))
genes_all <- unique(genes_all[, .(gene_name, seqname, start, end)], by = "gene_name")
par_csv <- read.csv(file.path(ROOT, "..", "single_cell", "PAR_genes.csv"))
par_names <- unique(gsub("_", "-", unlist(par_csv[, sapply(par_csv, is.character), drop = FALSE])))
# Symbol match misses genes whose GENCODE symbol is the bare ENSG id (cellranger
# keeps the old clone-based name, e.g. AC011751.1 = ENSG00000286201, a chrY lncRNA
# at 16.6 Mb). Resolve those through the features file's Ensembl ids.
genes_by_id <- as.data.table(load_gencode_genes(gtf))[, .(gene_id = sub("\\..*$", "", gene_id), seqname, start, end)]
feat <- fread(cmd = paste("zcat", shQuote(file.path(ROOT, "GEO/GSE336056/GSE336056_UTSW20_L1_features.tsv.gz"))), header = FALSE)
sym2id <- setNames(feat$V1, gsub("_", "-", feat$V2))
gene_class <- function(g) {
  m <- genes_all[match(g, gene_name)]
  miss <- is.na(m$seqname)
  if (any(miss)) {
    m2 <- genes_by_id[match(sym2id[g[miss]], gene_id)]
    m$seqname[miss] <- m2$seqname; m$start[miss] <- m2$start; m$end[miss] <- m2$end
  }
  cls <- rep("autosomal", length(g))
  cls[m$seqname %in% c("chrX", "X")] <- "chrX"
  cls[m$seqname %in% c("chrX", "X") & (m$start < 2781479 | m$end > 155701383)] <- "PAR"
  cls[m$seqname %in% c("chrY", "Y")] <- "chrY"
  cls[m$seqname %in% c("chr3", "3") & !is.na(m$end) & m$end < 90.5e6] <- "chr3p"
  cls[g %in% par_names] <- "PAR"
  cls[g %in% c("TTTY15")] <- "chrY"
  cls[is.na(m$seqname) & !(g %in% par_names) & g != "TTTY15"] <- "unannotated"
  cls
}

# ---- cells ------------------------------------------------------------------------------
CALLS <- fread(file.path(ROOT, "analysis/GSE336056/GSE336056_chr3p_dosage_per_cell.tsv"))[arm == "sg3p"]
mats <- list(); md <- list()
for (f in list.files(STAGE, "^sg3p_.*\\.rds$", full.names = TRUE)) {
  sid <- sub("\\.rds$", "", basename(f)); m <- readRDS(f)
  n <- Matrix::colSums(m); keep <- n >= MIN_UMI; m <- m[, keep, drop = FALSE]; n <- n[keep]
  cl <- CALLS[sample == sid]
  lane <- sub("_.*$", "", cl$cell[1])                     # per-cell TSV ids carry the lane prefix
  ids <- if (!startsWith(colnames(m)[1], lane)) paste0(lane, "_", colnames(m)) else colnames(m)
  hit <- match(ids, cl$cell)
  stopifnot(mean(!is.na(hit)) > 0.99)
  st <- fifelse(cl$loss3p[hit], "p3_loss", "p3_intact")
  ok <- !is.na(st)
  mats[[sid]] <- m[, ok, drop = FALSE]; colnames(mats[[sid]]) <- ids[ok]
  md[[sid]] <- data.table(cell = ids[ok], sample = sid, day = tstrsplit(sid, "_")[[2]],
                          nUMI = n[ok], status = st[ok], chr3p_l2r = cl$chr3p_l2r[hit[ok]])
}
M <- do.call(cbind, mats); MD <- rbindlist(md); stopifnot(identical(colnames(M), MD$cell))
MD[, status := factor(status, c("p3_intact", "p3_loss"))]
logf("cells: %d (%s)", nrow(MD), paste(sprintf("%s=%d", names(table(MD$status)), table(MD$status)), collapse = ", "))
print(MD[, .N, by = .(day, sample, status)][order(day, sample, status)])

# ---- expression + gene filter ------------------------------------------------------------
det_null <- Matrix::rowMeans(M[, MD$status == "p3_loss"] > 0)
det_int  <- Matrix::rowMeans(M[, MD$status == "p3_intact"] > 0)
keep_g <- det_null >= MIN_PCT | det_int >= MIN_PCT
logf("genes tested: %d of %d (detection >= %.0f%% in either group)", sum(keep_g), nrow(M), 100 * MIN_PCT)
Mk <- M[keep_g, ]
E <- log2(1 + t(t(Mk) / MD$nUMI) * 1e4)          # dense-ish; ~13k x 6.6k fine
E <- as.matrix(E)
MD[, cngeneson := scale(Matrix::colSums(Mk > 0))[, 1]]
# direct per-cell dosage summaries (memory: MAST logFC is not a dosage estimate)
per10k <- function(idx) 1e4 * Matrix::rowSums(Mk[, idx, drop = FALSE]) / sum(MD$nUMI[idx])
direct <- data.table(gene = rownames(Mk),
  mean_per10k_loss = per10k(MD$status == "p3_loss"), mean_per10k_intact = per10k(MD$status == "p3_intact"),
  pct_det_loss = 100 * det_null[keep_g], pct_det_intact = 100 * det_int[keep_g])
direct[, log2_ratio_direct := log2((mean_per10k_loss + 1e-3) / (mean_per10k_intact + 1e-3))]

# ---- MAST fits -----------------------------------------------------------------------------
fit_mast <- function(idx, tag, covars) {
  cd <- as.data.frame(MD[idx, c("cell", "status", "cngeneson", "sample", "day"), with = FALSE])
  rownames(cd) <- cd$cell; cd$wellKey <- cd$cell
  form <- reformulate(c("status", covars))
  sca <- FromMatrix(E[, idx, drop = FALSE], cData = cd, fData = data.frame(primerid = rownames(E), row.names = rownames(E)))
  logf("[%s] zlm %s on %d cells (%d null) ...", tag, deparse(form), length(idx), sum(cd$status == "p3_loss"))
  z <- zlm(form, sca, method = "bayesglm", ebayes = TRUE, parallel = TRUE)
  s <- summary(z, doLRT = "statusp3_loss", logFC = TRUE)$datatable
  h <- s[contrast == "statusp3_loss" & component == "H", .(gene = primerid, p = `Pr(>Chisq)`)]
  l <- s[contrast == "statusp3_loss" & component == "logFC", .(gene = primerid, mast_logFC = coef, ci_lo = ci.lo, ci_hi = ci.hi)]
  r <- merge(l, h, by = "gene"); r[, padj := p.adjust(p, "BH")]
  setnames(r, c("mast_logFC", "ci_lo", "ci_hi", "p", "padj"), paste0(c("logFC", "ci_lo", "ci_hi", "p", "padj"), "_", tag))
  r
}
res <- fit_mast(seq_len(nrow(MD)), "pooled", c("cngeneson", "sample"))
res <- merge(res, fit_mast(which(MD$day == "D5"), "D5", c("cngeneson", "sample")), by = "gene", all.x = TRUE)
res <- merge(res, fit_mast(which(MD$day == "D10"), "D10", c("cngeneson", "sample")), by = "gene", all.x = TRUE)
res <- merge(res, direct, by = "gene")
res[, gene_class := gene_class(gene)]
res <- merge(res, genes_all[, .(gene = gene_name, chrom = seqname, start)], by = "gene", all.x = TRUE)
res[, sign_consistent := sign(logFC_D5) == sign(logFC_D10) & sign(logFC_D5) == sign(logFC_pooled)]
setorder(res, p_pooled)
fwrite(res, file.path(OUT, "GSE336056_sg3p_3pLoss_vs_intact_MAST.tsv"), sep = "\t")
writeLines(c(
  "GSE336056 sg3p arm: 3p-loss vs 3p-intact (MAST hurdle, bayesglm + ebayes)",
  sprintf("Cells: %d 3p-loss vs %d 3p-intact, >= %d UMI, 4 HTO libraries (2 per day)", sum(MD$status == "p3_loss"), sum(MD$status == "p3_intact"), MIN_UMI),
  "Status: loss3p from chr3p_scan.R (3p l2r < log2(0.8) & |3q l2r| < 0.25 vs same-day sgNTC).",
  "Expression: log2(1 + 1e4 * counts / nUMI). Genes: detection >= 10% in either group.",
  "Pooled: ~ status + cngeneson + sample. Per-day: ~ status + cngeneson + sample (2 libs).",
  "Positive logFC = higher in 3p-loss. padj = BH within each fit.",
  "gene_class: chrY, PAR, chrX, chr3p (chr3 end < 90.5 Mb), autosomal.",
  "GSEA ranks: MAST logFC_pooled, autosomal + chrX only (chr3p/chrY/PAR/unannotated excluded)."),
  file.path(OUT, "model_specification_3p.txt"))

# ---- controls ------------------------------------------------------------------------------
show <- function(d) print(d[, .(gene, gene_class, logFC_pooled, padj_pooled, logFC_D5, logFC_D10, log2_ratio_direct, pct_det_loss, pct_det_intact)], nrows = 200)
logf("=== dosage control: direct log2 ratio by gene class (expect ~ -1 for chr3p) ===")
print(res[, .(n = .N, median_direct_l2r = round(median(log2_ratio_direct), 3), median_mast_logFC = round(median(logFC_pooled, na.rm = TRUE), 3)), by = gene_class])
logf("=== chr3p tumour suppressors (VHL, PBRM1, SETD2, BAP1) ==="); show(res[gene %in% c("VHL", "PBRM1", "SETD2", "BAP1")])
logf("=== top 40 autosomal/X by pooled P (|logFC| sorted within) ===")
show(res[gene_class %in% c("autosomal", "chrX")][1:40])
logf("=== counts: padj_pooled < 0.05 by class & direction ===")
print(res[padj_pooled < 0.05, .N, by = .(gene_class, dir = fifelse(logFC_pooled > 0, "up", "down"))][order(gene_class, dir)])

# ---- GSEA ------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(fgsea); library(msigdbr) })
rk <- res[gene_class %in% c("autosomal", "chrX") & is.finite(logFC_pooled)]
ranks <- sort(setNames(rk$logFC_pooled, rk$gene), decreasing = TRUE)
set.seed(1)
gs <- list(H = msigdbr(species = "Homo sapiens", collection = "H"),
           GOBP = msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
gsea <- rbindlist(lapply(names(gs), function(k) {
  sets <- split(gs[[k]]$gene_symbol, gs[[k]]$gs_name)
  r <- as.data.table(fgsea(sets, ranks, minSize = 15, maxSize = 500, eps = 0, nproc = NC))
  r[, collection := k]; r[, leadingEdge := sapply(leadingEdge, paste, collapse = ",")]; r
}))
setorder(gsea, padj, pval)
fwrite(gsea, file.path(OUT, "GSE336056_sg3p_3pLoss_vs_intact_GSEA.tsv"), sep = "\t")
logf("=== GSEA Hallmark (padj < 0.1) ===")
print(gsea[collection == "H" & padj < 0.1, .(pathway, NES = round(NES, 2), padj = signif(padj, 2), size)], nrows = 60)
logf("=== GSEA GO:BP top 25 ===")
print(gsea[collection == "GOBP"][1:25, .(pathway, NES = round(NES, 2), padj = signif(padj, 2), size)])
logf("done -> %s", OUT)
