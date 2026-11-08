# Freeze inputs for the main HT1080 RNA-ATAC promoter integration figure.
# Sources: ATACseq/SEPT_08_2026_REPROCESS/integration/thesis_RNA_SEPT11 (code/80_integrate_thesis_RNA.R;
# ID join with symbol fallback), the thesis RNA programme model (voom log-CPM), and ATAC annotation.
# Nothing is refitted. Locus coverage is extracted afterwards under WSL
# (extract_ht1080_atac_locus_coverage.py with this folder's locus_definitions.csv).
# Run: Rscript --vanilla <absolute path to this script>
suppressPackageStartupMessages({library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)
  library(GenomicRanges); library(limma); library(digest); library(jsonlite)})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
fig <- dirname(dirname(normalizePath(script, winslash = "/")))
projects <- dirname(dirname(dirname(fig)))
out <- file.path(fig, "source_data/ht1080_rna_atac_integration")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
atac_run <- file.path(projects, "ATACseq/SEPT_08_2026_REPROCESS")
integ_dir <- file.path(atac_run, "integration/thesis_RNA_SEPT11")
rna_dir <- file.path(projects, "RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1")
src <- c(genes = file.path(integ_dir, "RNA_ATAC_promoter_genes.csv"),
         correlations = file.path(integ_dir, "continuous_agreement.csv"),
         rna_consistency = file.path(rna_dir, "tables/gene_consistency.csv"),
         rna_model = file.path(rna_dir, "models/DE.rds"),
         atac_all = file.path(atac_run, "differential/primary_TMM_all.csv"),
         promoter_links = file.path(atac_run, "annotation/promoter_peak_gene_links.csv"),
         knownGene = file.path(atac_run, "config/genes_knownGene.tsv"))
strip <- function(x) sub("[.].*$", "", x)

genes <- fread(src[["genes"]], colClasses = list(character = "entrez_id"))
stopifnot(nrow(genes) == 13300, genes[joint == "concordant", .N] == 39)
# Figure threshold (author decision, 22 Sep 2026): RNA stays at the thesis criteria; ATAC promoter peaks at FDR < 0.10.
ATAC_FDR <- 0.10
plinks <- fread(src[["promoter_links"]], colClasses = list(character = "gene_id"))
pk <- plinks[, .(gained_010 = sum(FDR < ATAC_FDR & logFC > 0), lost_010 = sum(FDR < ATAC_FDR & logFC < 0)), by = .(entrez_id = gene_id)]
genes <- merge(genes, pk, by = "entrez_id", all.x = TRUE)
stopifnot(!anyNA(genes$gained_010))
genes[, atac_class_010 := fifelse(gained_010 > 0 & lost_010 > 0, "mixed", fifelse(gained_010 > 0, "gained",
                          fifelse(lost_010 > 0, "lost", "ns")))]
genes[, joint_010 := fifelse(rna_class == "not_DE", "RNA_not_DE", fifelse(atac_class_010 == "ns", "RNA_DE_ATAC_ns",
                     fifelse(atac_class_010 == "mixed", "mixed_promoter",
                     fifelse((rna_class == "up") == (atac_class_010 == "gained"), "concordant", "opposite"))))]
print(genes[, .N, by = .(rna_class, atac_class_010)][order(rna_class, atac_class_010)])
print(genes[rna_class != "not_DE" & atac_class_010 != "ns", .N, by = .(joint_010, rna_class)])
fwrite(genes, file.path(out, "integration_genes.csv"))
file.copy(src[["correlations"]], file.path(out, "continuous_agreement.csv"), overwrite = TRUE)

# Where the 662 thesis DE genes go (panel a)
rna <- fread(src[["rna_consistency"]])
de <- rna[directionally_consistent == TRUE][, ensembl := strip(gene_id)]
kg <- fread(src[["knownGene"]], colClasses = list(character = "gene_id"))
atac_all <- fread(src[["atac_all"]])
map <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = de$ensembl, keytype = "ENSEMBL",
                                           columns = "ENTREZID"))[!is.na(ENTREZID)]
de <- merge(de, map[, .(ENTREZID = ENTREZID[1]), by = ENSEMBL], by.x = "ensembl", by.y = "ENSEMBL", all.x = TRUE)
prom <- kg[, .(gene_id, chrom, s = pmax(1L, tss - 2000L), e = tss + 2000L)]
prom[, has_peak := overlapsAny(GRanges(chrom, IRanges(s, e)), GRanges(atac_all$chrom, IRanges(atac_all$start + 1L, atac_all$end)))]
de[, category := fifelse(ensembl %in% genes$ensembl, "tested_promoter",
                 fifelse(is.na(ENTREZID), "no_ID_or_name_match",
                 fifelse(!ENTREZID %in% kg$gene_id, "no_knownGene_TSS",
                 fifelse(!ENTREZID %in% prom[has_peak == TRUE, gene_id], "no_promoter_peak", "ambiguous_mapping"))))]
de[, direction := fifelse(logFC > 0, "up", "down")]
stopifnot(nrow(de) == 662, de[category == "tested_promoter", .N] == genes[rna_class != "not_DE", .N])
fwrite(de[, .(gene_id, gene_name, gene_type, direction, logFC, FDR, mean_logCPM, category)], file.path(out, "DE_gene_attrition.csv"))
print(dcast(de[, .N, by = .(category, direction)], category ~ direction, value.var = "N", fill = 0))

# Representative loci (panel e): lowest RNA FDR per evidence class, excluding loci in the ATAC main figure
already <- c("DCC", "DDIT4-AS1", "INSR", "CNTN4", "SMPD3")
conc <- genes[joint == "concordant" & !gene_name %in% already & !grepl("^(LOC|ENSG)", gene_name)][order(rna_FDR)]
picks <- rbind(conc[chrom != "chr19"][1:2][, class := "concordant loss"],
               conc[chrom == "chr19" & grepl("^Z(NF|FP)", gene_name)][1][, class := "concordant loss, chr19 KZFP"])
# Author selection, 22 Sep 2026: add the two concordant gains at ATAC FDR < 0.10 and ONECUT2; UNC5A removed.
extra <- data.table(gene_name = c("PHF21B", "PDPN", "ONECUT2"),
                    class = c("concordant gain (ATAC FDR < 0.10)", "concordant gain (ATAC FDR < 0.10)",
                              "RNA up, promoter trending open"))
picks <- rbind(picks, merge(genes, extra, by = "gene_name")[match(extra$gene_name, gene_name)], fill = TRUE)
stopifnot(nrow(picks) == 6, !anyNA(picks$entrez_id))
links <- fread(src[["promoter_links"]], colClasses = list(character = "gene_id"))
loci <- rbindlist(lapply(seq_len(nrow(picks)), function(i) {
  p <- picks[i]; l <- links[gene_id == p$entrez_id]
  tss <- l$tss[1]; strand <- l$strand[1]
  data.table(peak_id = paste0("locus_", i), label = p$gene_name, class = p$class, chrom = p$chrom, tss = tss, strand = strand,
             window_start = tss - 5000L, window_end = tss + 5000L, rna_logFC = p$rna_logFC, rna_FDR = p$rna_FDR,
             atac_median_logFC = p$atac_median_logFC, atac_min_peak_FDR = p$atac_min_peak_FDR,
             peaks = paste(l$peak_id, collapse = ";"))
}))
fwrite(loci, file.path(out, "locus_definitions.csv"))
peaks <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  ids <- strsplit(loci$peaks[i], ";")[[1]]
  atac_all[peak_id %in% ids, .(locus = loci$peak_id[i], peak_id, start, end, logFC, FDR, direction)]
}))
fwrite(peaks, file.path(out, "locus_peaks.csv"))
gm <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  kg[chrom == loci$chrom[i] & end > loci$window_start[i] & start < loci$window_end[i]][, locus := loci$peak_id[i]]
}))
gm[, symbol := unname(mapIds(org.Hs.eg.db, gene_id, "SYMBOL", "ENTREZID"))]
fwrite(gm, file.path(out, "locus_genes.csv"))

# Per-sample RNA expression for the loci (voom log-CPM, programme interaction model, 12 libraries)
m <- readRDS(src[["rna_model"]])
v <- m$interaction$voom
ids <- genes[gene_name %in% loci$label, .(gene_name, ensembl)]
rows <- match(ids$ensembl, strip(rownames(v$E)))
stopifnot(!anyNA(rows))
expr <- rbindlist(lapply(seq_len(nrow(ids)), function(i) {
  data.table(gene_name = ids$gene_name[i], sample = colnames(v$E), logCPM = v$E[rows[i], ])
}))
expr[, condition := fifelse(grepl("-LOY-", sample), "LOY", "CTRL")]
expr[, period := sub("^.*__", "", sample)]
expr[, replicate := as.integer(sub("^HT1080-[A-Za-z]+-([0-9]+)__.*$", "\\1", sample))]
fwrite(expr, file.path(out, "locus_RNA_logCPM.csv"))

fwrite(data.table(role = names(src), path = src, sha256 = vapply(src, digest, character(1), algo = "sha256", file = TRUE)),
       file.path(out, "source_manifest.csv"))
print(loci[, .(label, class, chrom, rna_logFC, rna_FDR, atac_median_logFC)])
