#!/usr/bin/env Rscript
# ============================================================================
# GSE336056: MAST DE, Yq-null vs Yq-intact cells WITHIN the sgYq arm.
#
# Cells: sgYq singlets, >= 10k UMI. Status from the 14-gene selection panel
# (yq_arm_check.R): Yq-null = Yq panel UMI == 0 & Yp panel UMI >= 5;
# Yq-intact = Yq panel UMI >= 2 & Yp >= 5. Everything else excluded.
# Model (corpus MAST convention, cf. scripts/mast_GSE309595_donor_cluster.R):
#   log2(1 + 1e4 * counts / nUMI); detection >= 10% in either group;
#   zlm(~ status + cngeneson + sample, bayesglm, ebayes); hurdle LRT on status.
#   sample (4 HTO libraries) encodes day (2 per day) so day is not a separate term.
# Fits: pooled (primary) + D5-only + D10-only (sign-consistency check).
# Gene classes flagged by GENCODE v49 coordinate: chrY (MSY), PAR1/PAR2 (chrX
# < 2.78 Mb or > 155.7 Mb) + PAR_genes.csv, TTTY15 explicit (escapes the MSY
# resolver). These are excluded from GSEA. Yq genes = positive controls; Yp,
# PAR = negative controls (Yp retained). X homologs of lost Yq genes and the
# SLIT/ROBO panel are reported explicitly.
# GSEA: fgsea on MAST logFC (ranked by effect, not signed P), autosomal + X
# only; MSigDB Hallmark and GO:BP.
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
  cls[g %in% par_names] <- "PAR"
  cls[g %in% c("TTTY15")] <- "chrY"
  cls[is.na(m$seqname) & !(g %in% par_names) & g != "TTTY15"] <- "unannotated"
  cls
}

# ---- cells ------------------------------------------------------------------------------
mats <- list(); md <- list()
for (f in list.files(STAGE, "^sgYq_.*\\.rds$", full.names = TRUE)) {
  sid <- sub("\\.rds$", "", basename(f)); m <- readRDS(f)
  n <- Matrix::colSums(m); keep <- n >= MIN_UMI; m <- m[, keep, drop = FALSE]; n <- n[keep]
  yp <- Matrix::colSums(m[intersect(SEL_YP, rownames(m)), , drop = FALSE])
  yq <- Matrix::colSums(m[intersect(SEL_YQ, rownames(m)), , drop = FALSE])
  st <- fifelse(yp >= 5 & yq == 0, "Yq_null", fifelse(yp >= 5 & yq >= 2, "Yq_intact", NA_character_))
  ok <- !is.na(st)
  mats[[sid]] <- m[, ok, drop = FALSE]
  md[[sid]] <- data.table(cell = colnames(m)[ok], sample = sid, day = tstrsplit(sid, "_")[[2]],
                          nUMI = n[ok], status = st[ok], yp_sel = yp[ok], yq_sel = yq[ok])
}
M <- do.call(cbind, mats); MD <- rbindlist(md); stopifnot(identical(colnames(M), MD$cell))
MD[, status := factor(status, c("Yq_intact", "Yq_null"))]
logf("cells: %d (%s)", nrow(MD), paste(sprintf("%s=%d", names(table(MD$status)), table(MD$status)), collapse = ", "))
print(MD[, .N, by = .(day, sample, status)][order(day, sample, status)])

# ---- expression + gene filter ------------------------------------------------------------
det_null <- Matrix::rowMeans(M[, MD$status == "Yq_null"] > 0)
det_int  <- Matrix::rowMeans(M[, MD$status == "Yq_intact"] > 0)
keep_g <- det_null >= MIN_PCT | det_int >= MIN_PCT
logf("genes tested: %d of %d (detection >= %.0f%% in either group)", sum(keep_g), nrow(M), 100 * MIN_PCT)
Mk <- M[keep_g, ]
E <- log2(1 + t(t(Mk) / MD$nUMI) * 1e4)          # dense-ish; ~13k x 6.6k fine
E <- as.matrix(E)
MD[, cngeneson := scale(Matrix::colSums(Mk > 0))[, 1]]
# direct per-cell dosage summaries (memory: MAST logFC is not a dosage estimate)
per10k <- function(idx) 1e4 * Matrix::rowSums(Mk[, idx, drop = FALSE]) / sum(MD$nUMI[idx])
direct <- data.table(gene = rownames(Mk),
  mean_per10k_null = per10k(MD$status == "Yq_null"), mean_per10k_intact = per10k(MD$status == "Yq_intact"),
  pct_det_null = 100 * det_null[keep_g], pct_det_intact = 100 * det_int[keep_g])
direct[, log2_ratio_direct := log2((mean_per10k_null + 1e-3) / (mean_per10k_intact + 1e-3))]

# ---- MAST fits -----------------------------------------------------------------------------
fit_mast <- function(idx, tag, covars) {
  cd <- as.data.frame(MD[idx, c("cell", "status", "cngeneson", "sample", "day"), with = FALSE])
  rownames(cd) <- cd$cell; cd$wellKey <- cd$cell
  form <- reformulate(c("status", covars))
  sca <- FromMatrix(E[, idx, drop = FALSE], cData = cd, fData = data.frame(primerid = rownames(E), row.names = rownames(E)))
  logf("[%s] zlm %s on %d cells (%d null) ...", tag, deparse(form), length(idx), sum(cd$status == "Yq_null"))
  z <- zlm(form, sca, method = "bayesglm", ebayes = TRUE, parallel = TRUE)
  s <- summary(z, doLRT = "statusYq_null", logFC = TRUE)$datatable
  h <- s[contrast == "statusYq_null" & component == "H", .(gene = primerid, p = `Pr(>Chisq)`)]
  l <- s[contrast == "statusYq_null" & component == "logFC", .(gene = primerid, mast_logFC = coef, ci_lo = ci.lo, ci_hi = ci.hi)]
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
fwrite(res, file.path(OUT, "GSE336056_sgYq_YqNull_vs_intact_MAST.tsv"), sep = "\t")
writeLines(c(
  "GSE336056 sgYq arm: Yq-null vs Yq-intact (MAST hurdle, bayesglm + ebayes)",
  sprintf("Cells: %d Yq-null vs %d Yq-intact, >= %d UMI, 4 HTO libraries (2 per day)", sum(MD$status == "Yq_null"), sum(MD$status == "Yq_intact"), MIN_UMI),
  "Expression: log2(1 + 1e4 * counts / nUMI). Genes: detection >= 10% in either group.",
  "Pooled: ~ status + cngeneson + sample. Per-day: ~ status + cngeneson + sample (2 libs).",
  "Positive logFC = higher in Yq-null. padj = BH within each fit.",
  "gene_class by GENCODE v49 coordinate: chrY (MSY + TTTY15), PAR (chrX PAR coords + PAR_genes.csv), chrX, autosomal.",
  "GSEA ranks: MAST logFC_pooled, autosomal + chrX only (chrY/PAR/unannotated excluded)."),
  file.path(OUT, "model_specification.txt"))

# ---- controls + watch panels ------------------------------------------------------------
show <- function(d) print(d[, .(gene, gene_class, logFC_pooled, padj_pooled, logFC_D5, logFC_D10, log2_ratio_direct, pct_det_null, pct_det_intact)], nrows = 200)
logf("=== positive controls: Yq genes ==="); show(res[gene %in% c(SEL_YQ, "TTTY14", "TTTY10", "TXLNGY", "PRORY")][order(start)])
logf("=== negative controls: Yp genes + PAR ==="); show(res[gene %in% c(SEL_YP, "PRKY", "LINC00278", "ZFY-AS1") | gene_class == "PAR"][order(gene_class, start)])
XHOM <- c(USP9Y = "USP9X", DDX3Y = "DDX3X", UTY = "KDM6A", KDM5D = "KDM5C", EIF1AY = "EIF1AX", NLGN4Y = "NLGN4X", TMSB4Y = "TMSB4X", ZFY = "ZFX", RPS4Y1 = "RPS4X", PRKY = "PRKX", TBL1Y = "TBL1X")
logf("=== X homologs of Y genes ==="); show(res[gene %in% XHOM][order(gene)])
ROBO <- c(paste0("ROBO", 1:4), paste0("SLIT", 1:3), paste0("SRGAP", 1:3), "CLASP1", "CLASP2", "PAK1", "NCK1", "NCK2")
logf("=== SLIT/ROBO panel ==="); show(res[gene %in% ROBO][order(gene)])
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
fwrite(gsea, file.path(OUT, "GSE336056_sgYq_YqNull_vs_intact_GSEA.tsv"), sep = "\t")
logf("=== GSEA Hallmark (padj < 0.1) ===")
print(gsea[collection == "H" & padj < 0.1, .(pathway, NES = round(NES, 2), padj = signif(padj, 2), size)], nrows = 60)
logf("=== GSEA GO:BP top 25 ===")
print(gsea[collection == "GOBP"][1:25, .(pathway, NES = round(NES, 2), padj = signif(padj, 2), size)])
logf("done -> %s", OUT)
