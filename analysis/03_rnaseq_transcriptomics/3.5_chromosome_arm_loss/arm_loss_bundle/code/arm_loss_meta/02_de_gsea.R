#!/usr/bin/env Rscript
# ============================================================================
# Arm-loss meta-analysis: per-contrast DESeq2 + GSEA.
#
# Replicate unit = independent clone. Technical/biological replicates of the SAME
# clone are summed (DESeq2 collapseReplicates), so a 3-clone arm contributes n=3,
# not n=6. Control libraries are replicates of one parental genotype and are kept.
#
# cis genes -- those on the lost arm/chromosome -- are half-dosage BY CONSTRUCTION
# and define the contrast. They are reported as a QC check (median log2FC should
# sit near -1) and then EXCLUDED from the GSEA ranks, exactly as chrY and chr3p
# genes were excluded in the GSE336056 single-cell DE.
#
# Two contrast types:
#   arm_loss       loss clones vs the matched control of the same background
#   p53_induction  inducible-p53 +dox vs -dox at FIXED karyotype (GSE150686).
#                  This is the reference for what p53 activation alone looks like,
#                  so the arm-loss programs can be read against it.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(fgsea); library(msigdbr)
})
ROOT <- "/global/home/hpc5173/projects/LOY/RNA/arm_loss_meta"
OUT  <- file.path(ROOT, "results"); DE <- file.path(OUT, "de"); GS <- file.path(OUT, "gsea")
for (d in c(DE, GS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
logf <- function(fmt, ...) cat(sprintf(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", fmt, "\n"), ...))

ds  <- readRDS(file.path(OUT, "datasets.rds"))
ann <- fread(cmd = paste("zcat", shQuote(file.path(OUT, "gene_annotation.tsv.gz"))))
setkey(ann, ensg)

# ---- contrast table --------------------------------------------------------------
# cis: "arm:1q" or "chrom:chr10,chr18" or "" (none to exclude)
CT <- fread(text = '
block|contrast|type|case|control|cis
GSE222379_main|A2058_1q|arm_loss|cell_line=A2058,condition=loss,arm_lost=1q|cell_line=A2058,condition=WT|arm:1q
GSE222379_main|A2058_7p|arm_loss|cell_line=A2058,condition=loss,arm_lost=7p|cell_line=A2058,condition=WT|arm:7p
GSE222379_main|A2058_8q|arm_loss|cell_line=A2058,condition=loss,arm_lost=8q|cell_line=A2058,condition=WT|arm:8q
GSE222379_main|AGS_1q_vsTK|arm_loss|cell_line=AGS,condition=loss,arm_lost=1q|cell_line=AGS,condition=TK|arm:1q
GSE222379_main|HCT116_8q|arm_loss|cell_line=HCT116,condition=loss,arm_lost=8q|cell_line=HCT116,condition=WT|arm:8q
GSE222379_main|MCF10A_1q|arm_loss|cell_line=MCF10A,condition=loss,arm_lost=1q|cell_line=MCF10A,condition=WT|arm:1q
GSE222379_main|MCF10Ahr_1q|arm_loss|cell_line=MCF10A_hRas,condition=loss,arm_lost=1q|cell_line=MCF10A_hRas,condition=WT|arm:1q
GSE222379_A2780|A2780_1q_vsWT|arm_loss|cell_line=A2780,condition=loss,arm_lost=1q|cell_line=A2780,condition=WT|arm:1q
GSE222379_A2780|A2780_1q_vsTK|arm_loss|cell_line=A2780,condition=loss,arm_lost=1q|cell_line=A2780,condition=TK|arm:1q
GSE336576|RPTEC_3p_early|arm_loss|condition=loss,passage=Early|condition=WT|arm:3p
GSE336576|RPTEC_3p_late|arm_loss|condition=loss,passage=Late|condition=WT|arm:3p
GSE150686|RPE1_mono1018|arm_loss|base=RM_1018|base=RPE1_p53KO|chrom:chr10,chr18
GSE150686|RPE1_mono13|arm_loss|base=RM_13|base=RPE1_p53KO|chrom:chr13
GSE150686|RPE1_mono19p|arm_loss|base=RM_19p|base=RPE1_p53_KD|arm:19p
PRJNA930361|ARPE19_chrY|arm_loss|condition=loss|condition=WT|chrom:chrY
GSE150686|ip53_diploid|p53_induction|base=RPE1_ip53_+D|base=RPE1_ip53_-D|
GSE150686|ip53_mono1018|p53_induction|base=RM_10_ip53_+D|base=RM_10_ip53_-D|
GSE150686|ip53_mono13|p53_induction|base=RM_13_ip53_+D|base=RM_13_ip53_-D|
', sep = "|", header = TRUE, colClasses = "character")
CT[is.na(cis), cis := ""]

sel <- function(cd, spec) {          # "col=val,col=val" -> logical index into cd
  keep <- rep(TRUE, nrow(cd))
  for (kv in strsplit(spec, ",")[[1]]) {
    p <- strsplit(kv, "=")[[1]]
    keep <- keep & !is.na(cd[[p[1]]]) & cd[[p[1]]] == p[2]
  }
  keep
}
cis_genes <- function(spec) {
  if (!nzchar(spec)) return(character(0))
  p <- strsplit(spec, ":")[[1]]; v <- strsplit(p[2], ",")[[1]]
  if (p[1] == "arm") ann[arm %in% v, ensg] else ann[chrom %in% v, ensg]
}

# ---- gene sets -------------------------------------------------------------------
logf("loading MSigDB ...")
H  <- msigdbr(species = "Homo sapiens", collection = "H")
BP <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
SETS <- list(H = split(H$gene_symbol, H$gs_name), GOBP = split(BP$gene_symbol, BP$gs_name))
# custom: the ribosomal-protein / translation axis that converges across the
# monosomy, 45,X and trisomy-21 literature. Protein-coding RPL*/RPS* only, so
# the many RP pseudogenes do not enter the set.
rp <- ann[biotype == "protein_coding" & grepl("^RP[LS][0-9]", symbol), unique(symbol)]
SETS$CUSTOM <- list(RIBOSOMAL_PROTEIN = rp)
logf("gene sets: H=%d GOBP=%d | ribosomal-protein set n=%d", length(SETS$H), length(SETS$GOBP), length(rp))

# ---- run -------------------------------------------------------------------------
qc <- list()
for (i in seq_len(nrow(CT))) {
  r <- CT[i]; cd <- copy(ds[[r$block]]$coldata); cnt <- ds[[r$block]]$counts
  i_case <- sel(cd, r$case); i_ctrl <- sel(cd, r$control)
  stopifnot(sum(i_case) > 0, sum(i_ctrl) > 0, !any(i_case & i_ctrl))

  sub <- cd[i_case | i_ctrl]
  sub[, grp := factor(fifelse(sel(sub, r$case), "case", "ctrl"), c("ctrl", "case"))]
  m <- cnt[, sub$sample, drop = FALSE]
  storage.mode(m) <- "integer"

  dds <- DESeqDataSetFromMatrix(m, as.data.frame(sub), ~ grp)
  dds <- collapseReplicates(dds, groupby = factor(sub$unit))     # clone = replicate unit
  n_case <- sum(colData(dds)$grp == "case"); n_ctrl <- sum(colData(dds)$grp == "ctrl")
  if (n_case < 2 || n_ctrl < 2) { logf("SKIP %s: n_case=%d n_ctrl=%d", r$contrast, n_case, n_ctrl); next }
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  logf("%s | %s: %d case vs %d ctrl units, %d genes", r$block, r$contrast, n_case, n_ctrl, nrow(dds))

  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.table(results(dds, contrast = c("grp", "case", "ctrl")), keep.rownames = "ensg")
  res <- merge(res, ann[, .(ensg, symbol, biotype, chrom, start, arm)], by = "ensg", all.x = TRUE)
  cg <- cis_genes(r$cis)
  res[, is_cis := ensg %in% cg]
  res[, `:=`(block = r$block, contrast = r$contrast, type = r$type, cis_spec = r$cis)]
  fwrite(res[order(pvalue)], file.path(DE, sprintf("%s_DE.tsv.gz", r$contrast)), sep = "\t")

  # QC: cis dosage. One lost copy of two -> log2FC near -1 before buffering.
  cis_lfc <- res[is_cis == TRUE & !is.na(log2FoldChange), median(log2FoldChange)]
  tr_lfc  <- res[is_cis == FALSE & !is.na(log2FoldChange), median(log2FoldChange)]
  qc[[r$contrast]] <- data.table(block = r$block, contrast = r$contrast, type = r$type,
    n_case = n_case, n_ctrl = n_ctrl, n_genes = nrow(res), n_cis = sum(res$is_cis),
    median_lfc_cis = round(cis_lfc, 3), median_lfc_trans = round(tr_lfc, 3),
    n_padj05 = res[padj < 0.05, .N])
  logf("   cis median log2FC = %s | trans = %s | padj<0.05: %d",
       ifelse(is.na(cis_lfc), "NA", sprintf("%+.3f", cis_lfc)), sprintf("%+.3f", tr_lfc),
       res[padj < 0.05, .N])

  # GSEA on the Wald statistic, trans genes only (cis excluded: they ARE the lesion)
  rk <- res[is_cis == FALSE & !is.na(stat) & !is.na(symbol) & symbol != ""]
  rk <- rk[, .(stat = mean(stat)), by = symbol]
  ranks <- sort(setNames(rk$stat, rk$symbol), decreasing = TRUE)
  set.seed(1)
  g <- rbindlist(lapply(names(SETS), function(k) {
    f <- as.data.table(fgsea(SETS[[k]], ranks, minSize = 10, maxSize = 500, eps = 0, nproc = NC))
    f[, collection := k]; f[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ",")]; f
  }))
  g[, `:=`(block = r$block, contrast = r$contrast, type = r$type)]
  fwrite(g[order(padj, pval)], file.path(GS, sprintf("%s_GSEA.tsv.gz", r$contrast)), sep = "\t")
}

QC <- rbindlist(qc)
fwrite(QC, file.path(OUT, "contrast_qc.tsv"), sep = "\t")
cat("\n================ cis-dosage QC ================\n"); print(QC, nrows = 40)
logf("done")
