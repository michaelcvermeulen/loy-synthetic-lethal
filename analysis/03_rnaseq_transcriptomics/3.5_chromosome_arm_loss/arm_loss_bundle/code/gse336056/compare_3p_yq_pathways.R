#!/usr/bin/env Rscript
# GSE336056: which pathways / genes move in BOTH engineered arm losses?
# Inputs: de/GSE336056_sgYq_YqNull_vs_intact_{MAST,GSEA}.tsv (de_yq_null_vs_intact.R)
#         de/GSE336056_sg3p_3pLoss_vs_intact_{MAST,GSEA}.tsv (de_3p_loss_vs_intact.R)
# Shared = padj < 0.05 in both, same NES sign. Genes: autosomal + chrX only, excluding
# chr3p (half-dosage in the 3p arm) and chrY/PAR (lost/retained in the Yq arm).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("LOY_SCRNA_ROOT", "/global/home/hpc5173/projects/LOY/scRNA")
DE <- file.path(ROOT, "analysis/GSE336056/de"); FDR <- 0.05

gy <- fread(file.path(DE, "GSE336056_sgYq_YqNull_vs_intact_GSEA.tsv"))
g3 <- fread(file.path(DE, "GSE336056_sg3p_3pLoss_vs_intact_GSEA.tsv"))
G <- merge(gy[, .(pathway, collection, NES_Yq = NES, padj_Yq = padj, size_Yq = size)],
           g3[, .(pathway, collection, NES_3p = NES, padj_3p = padj, size_3p = size)],
           by = c("pathway", "collection"))
G[, class := fcase(padj_Yq < FDR & padj_3p < FDR & sign(NES_Yq) == sign(NES_3p), "shared",
                   padj_Yq < FDR & padj_3p < FDR, "opposite",
                   padj_Yq < FDR, "Yq_only", padj_3p < FDR, "3p_only", default = "neither")]
G[, dir := fifelse(NES_3p > 0, "up", "down")]

cat("=== pathway classes (padj <", FDR, ") ===\n")
print(dcast(G, collection ~ class, fun.aggregate = length, value.var = "pathway"))
cat("\n=== NES concordance (Spearman) ===\n")
for (k in c("H", "GOBP")) {
  x <- G[collection == k]
  cat(sprintf("%-5s all: rho = %.2f (n=%d) | sig in either: rho = %.2f (n=%d)\n", k,
              cor(x$NES_Yq, x$NES_3p, method = "spearman"), nrow(x),
              cor(x[class != "neither", NES_Yq], x[class != "neither", NES_3p], method = "spearman"), x[class != "neither", .N]))
}
fmt <- function(d) d[, .(pathway, NES_3p = round(NES_3p, 2), padj_3p = signif(padj_3p, 2), NES_Yq = round(NES_Yq, 2), padj_Yq = signif(padj_Yq, 2))]
cat("\n=== Hallmark: shared ===\n"); print(fmt(G[collection == "H" & class == "shared"][order(-abs(NES_3p))]), nrows = 60)
cat("\n=== Hallmark: opposite ===\n"); print(fmt(G[collection == "H" & class == "opposite"]), nrows = 60)
cat("\n=== Hallmark: one arm only ===\n"); print(fmt(G[collection == "H" & class %in% c("3p_only", "Yq_only")][order(class, -abs(NES_3p))]), nrows = 60)
cat("\n=== GO:BP shared, top 20 per direction by mean |NES| ===\n")
gb <- G[collection == "GOBP" & class == "shared"][, mNES := (NES_Yq + NES_3p) / 2]
print(fmt(gb[mNES > 0][order(-mNES)][1:min(20, .N)])); print(fmt(gb[mNES < 0][order(mNES)][1:min(20, .N)]))
cat("\n=== GO:BP opposite ===\n"); print(fmt(G[collection == "GOBP" & class == "opposite"][order(padj_3p)]), nrows = 40)

# ---- genes ----
my <- fread(file.path(DE, "GSE336056_sgYq_YqNull_vs_intact_MAST.tsv"))
m3 <- fread(file.path(DE, "GSE336056_sg3p_3pLoss_vs_intact_MAST.tsv"))
K <- merge(my[gene_class %in% c("autosomal", "chrX"), .(gene, lfc_Yq = logFC_pooled, padj_Yq = padj_pooled, sc_Yq = sign_consistent)],
           m3[gene_class %in% c("autosomal", "chrX"), .(gene, lfc_3p = logFC_pooled, padj_3p = padj_pooled, sc_3p = sign_consistent)], by = "gene")
K[, class := fcase(padj_Yq < FDR & padj_3p < FDR & sign(lfc_Yq) == sign(lfc_3p), "shared",
                   padj_Yq < FDR & padj_3p < FDR, "opposite", padj_Yq < FDR, "Yq_only", padj_3p < FDR, "3p_only", default = "neither")]
cat("\n=== genes (autosomal + X, excl. chr3p / chrY / PAR) ===\n"); print(K[, .N, by = class][order(class)])
sy <- K[padj_Yq < FDR]; cat(sprintf("of %d Yq-sig genes: %d sig in 3p, %.0f%% same sign in 3p regardless of sig (binomial P = %.1e)\n",
  nrow(sy), sy[padj_3p < FDR, .N], 100 * mean(sign(sy$lfc_Yq) == sign(sy$lfc_3p), na.rm = TRUE),
  binom.test(sum(sign(sy$lfc_Yq) == sign(sy$lfc_3p), na.rm = TRUE), sum(!is.na(sy$lfc_3p)))$p.value))
cat(sprintf("logFC Spearman, genes sig in either: rho = %.2f\n", K[class != "neither", cor(lfc_Yq, lfc_3p, method = "spearman")]))
cat("\n=== shared genes, top 25 up / down by mean logFC ===\n")
sh <- K[class == "shared"][, m := (lfc_Yq + lfc_3p) / 2]
print(sh[order(-m)][1:min(25, .N), .(gene, lfc_3p = round(lfc_3p, 2), lfc_Yq = round(lfc_Yq, 2), sc_3p, sc_Yq)])
print(sh[order(m)][1:min(25, .N), .(gene, lfc_3p = round(lfc_3p, 2), lfc_Yq = round(lfc_Yq, 2), sc_3p, sc_Yq)])
cat("\n=== opposite-direction genes ===\n"); print(K[class == "opposite", .(gene, lfc_3p = round(lfc_3p, 2), lfc_Yq = round(lfc_Yq, 2))], nrows = 40)

fwrite(G[order(collection, class, padj_3p)], file.path(DE, "GSE336056_3p_vs_Yq_pathway_overlap.tsv"), sep = "\t")
fwrite(K[order(class, padj_3p)], file.path(DE, "GSE336056_3p_vs_Yq_gene_overlap.tsv"), sep = "\t")
cat("done\n")
