# Peak-level RNA-ATAC links for the supplementary peak-type scatter (default RNA 0.05 / ATAC 0.10).
# RNA: thesis limma-voom result (gene_consistency.csv); DE at FDR < RNA_FDR with the thesis
#   |log2FC| >= 1 and same-direction (AUG/MAY) rules (RNA held at the thesis FDR 0.05).
# ATAC: primary edgeR peaks; significant at peak FDR < ATAC_FDR (default 0.10).
# Gene matching: one-to-one Ensembl -> Entrez (org.Hs.eg.db); fallback for unmapped RNA genes:
#   unique SYMBOL -> Entrez in org.Hs.eg.db, only if that Entrez is not already claimed.
# Link sets: (1) promoter links (peak overlaps TSS +/- 2 kb; every link retained);
#            (2) every tested peak to its nearest knownGene TSS (exploratory positional link).
# Run: Rscript --vanilla <absolute path to this script>
suppressPackageStartupMessages({library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)
  library(digest)})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
fig <- dirname(dirname(normalizePath(script, winslash = "/")))
projects <- dirname(dirname(dirname(fig)))
# Optional trailing args: RNA_FDR ATAC_FDR (default 0.05 0.10; author decision 22 Sep 2026 - RNA held at the thesis FDR). Non-default thresholds write to a
# suffixed folder, e.g. source_data/ht1080_rna_atac_peak_scatter_rna005_atac010.
args <- commandArgs(trailingOnly = TRUE)
RNA_FDR <- if (length(args) >= 1) as.numeric(args[1]) else 0.05
ATAC_FDR <- if (length(args) >= 2) as.numeric(args[2]) else 0.10
suffix <- if (RNA_FDR == 0.05 && ATAC_FDR == 0.10) "" else
  sprintf("_rna%03d_atac%03d", round(RNA_FDR * 100), round(ATAC_FDR * 100))
out <- file.path(fig, paste0("source_data/ht1080_rna_atac_peak_scatter", suffix))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
atac_run <- file.path(projects, "ATACseq/SEPT_08_2026_REPROCESS")
src <- c(rna = file.path(projects, "RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1/tables/gene_consistency.csv"),
         promoter_links = file.path(atac_run, "annotation/promoter_peak_gene_links.csv"),
         annotated_peaks = file.path(atac_run, "annotation/annotated_differential_peaks.csv"))

rna <- fread(src[["rna"]])[, ensembl := sub("[.].*$", "", gene_id)]
map <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = rna$ensembl, keytype = "ENSEMBL", columns = "ENTREZID"))[!is.na(ENTREZID)]
map[, n_e := uniqueN(ENTREZID), by = ENSEMBL][, n_g := uniqueN(ENSEMBL), by = ENTREZID]
id <- map[n_e == 1 & n_g == 1, .(ensembl = ENSEMBL, entrez = ENTREZID, join_method = "ID")]
left <- rna[!ensembl %in% id$ensembl & !gene_name %in% rna[duplicated(gene_name), gene_name]]
sym <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = unique(left$gene_name), keytype = "SYMBOL", columns = "ENTREZID"))[!is.na(ENTREZID)]
sym <- sym[, .(ENTREZID = if (uniqueN(ENTREZID) == 1) ENTREZID[1] else NA_character_), by = SYMBOL][!is.na(ENTREZID)]
sym <- sym[!ENTREZID %in% id$entrez & !duplicated(ENTREZID) & !duplicated(ENTREZID, fromLast = TRUE)]
fb <- merge(left[, .(ensembl, gene_name)], sym, by.x = "gene_name", by.y = "SYMBOL")[, .(ensembl, entrez = ENTREZID, join_method = "symbol_fallback")]
key <- rbind(id, fb)
rna <- merge(rna, key, by = "ensembl")
stopifnot(!anyDuplicated(rna$entrez))
rna[, rna_class := fifelse(FDR < RNA_FDR & abs(logFC) >= 1 & same_direction, fifelse(logFC > 0, "up", "down"), "not_DE")]
rna_cols <- rna[, .(entrez, gene_name, rna_logFC = logFC, rna_FDR = FDR, rna_class, join_method)]

classify <- function(dt) {
  dt[, peak_class := fifelse(FDR < ATAC_FDR, fifelse(logFC > 0, "gained", "lost"), "ns")]
  dt[, category := fifelse(rna_class == "not_DE" & peak_class == "ns", "neither",
                   fifelse(rna_class == "not_DE", "peak_only", fifelse(peak_class == "ns", "RNA_only",
                   fifelse((rna_class == "up") == (peak_class == "gained"), "concordant", "opposite"))))]
  dt
}
pl <- fread(src[["promoter_links"]], colClasses = list(character = "gene_id"))
pl <- merge(pl[, .(peak_id, chrom, peak_logFC = logFC, FDR, entrez = gene_id)], rna_cols, by = "entrez")
setnames(pl, "peak_logFC", "logFC"); pl <- classify(pl)
ap <- fread(src[["annotated_peaks"]], colClasses = list(character = "nearest_entrez_id"))
nl <- merge(ap[, .(peak_id, chrom, logFC, FDR, entrez = nearest_entrez_id, tss_distance = nearest_tss_distance_signed,
                   annotation)], rna_cols, by = "entrez")
nl <- classify(nl)
fwrite(pl, file.path(out, "promoter_peak_links.csv"))
fwrite(nl, file.path(out, "nearest_gene_peak_links.csv"))

summ <- function(dt, set) {
  x <- dt[, .(links = .N, peaks = uniqueN(peak_id), genes = uniqueN(entrez)), by = category]
  x[, set := set]
  y <- dt[category %in% c("concordant", "opposite"), .(links = .N, peaks = uniqueN(peak_id), genes = uniqueN(entrez)),
          by = .(category, direction = paste(rna_class, peak_class, sep = "/"))][, set := set]
  list(x, y)
}
s1 <- summ(pl, "promoter"); s2 <- summ(nl, "nearest_gene")
fwrite(rbind(s1[[1]], s2[[1]]), file.path(out, "category_counts.csv"))
fwrite(rbind(s1[[2]], s2[[2]]), file.path(out, "joint_direction_counts.csv"))
sp <- function(dt, sub) { d <- if (sub == "all") dt else dt[category %in% c("concordant", "opposite")]
  data.table(subset = sub, n = nrow(d), spearman = if (nrow(d) > 2) cor(d$rna_logFC, d$logFC, method = "spearman") else NA_real_) }
cors <- rbind(cbind(set = "promoter", rbind(sp(pl, "all"), sp(pl, "both significant"))),
              cbind(set = "nearest_gene", rbind(sp(nl, "all"), sp(nl, "both significant"))))
fwrite(cors, file.path(out, "spearman.csv"))
fwrite(data.table(role = names(src), path = src, sha256 = vapply(src, digest, character(1), algo = "sha256", file = TRUE),
                  RNA_FDR = RNA_FDR, ATAC_FDR = ATAC_FDR), file.path(out, "source_manifest.csv"))
cat("RNA genes matched:", nrow(rna), "(ID", sum(rna$join_method == "ID"), "; symbol", sum(rna$join_method == "symbol_fallback"), ")\n")
print(rbind(s1[[2]], s2[[2]])); print(cors)
print(dcast(rbind(s1[[1]], s2[[1]]), category ~ set, value.var = "links"))
