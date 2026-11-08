# Reconcile the HT1080 RNA-ATAC promoter integration onto the thesis RNA source.
#
# Thesis RNA (author-selected 18 Sep 2026): RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE,
#   programme_analysis_v1 gene_consistency.csv; limma-voom ~ batch * treatment, combined
#   Y-depleted vs sgCTRL contrast over 12 core libraries; GENCODE v49; 25,526 expressed
#   autosomal/non-PAR X genes. Thesis DE rule = directionally_consistent
#   (FDR < 0.05, |log2FC| >= 1, same direction in AUG_2025 and MAY_2026).
# ATAC: existing primary edgeR promoter links (TSS +/- 2 kb, UCSC knownGene / Entrez).
# Primary join: RNA Ensembl gene ID -> Entrez via org.Hs.eg.db, unique in both directions,
#   with an exact-symbol fallback for genes lacking a one-to-one ID mapping; same chromosome
#   required in both assays. Sensitivity: exact symbol join only (as used by
#   integrate_programmes_ATAC_v1.py). Also compares the thesis RNA branch with the earlier
#   DESeq2 ~ batch + treatment branch used by integration/RNAseq (unversioned Ensembl IDs).
# No model is refitted; overlap counts are descriptive (no overlap P-value is assigned).
# Run: Rscript --vanilla <absolute path to this script>
options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)
  library(digest); library(jsonlite)})
out <- "integration/thesis_RNA_SEPT11"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sources <- c(thesis_RNA = normalizePath(
  "../../RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1/tables/gene_consistency.csv",
  winslash = "/", mustWork = TRUE),
  earlier_RNA_snapshot = "config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv",
  earlier_joint_genes = "integration/RNAseq/ATAC_RNA_joint_significant_genes.csv",
  ATAC_promoter_links = "annotation/promoter_peak_gene_links.csv",
  ATAC_promoter_summary = "annotation/promoter_gene_summary.csv")
fwrite(data.table(role = names(sources), path = sources,
                  sha256 = vapply(sources, digest, character(1), algo = "sha256", file = TRUE)),
       file.path(out, "source_manifest.csv"))

rna <- fread(sources[["thesis_RNA"]])
stopifnot(nrow(rna) == 25526, !anyDuplicated(rna$gene_id), all(rna$model == "batch_condition"))
rna[, ensembl := sub("\\..*$", "", gene_id)]
rna[, rna_class := fifelse(directionally_consistent & logFC > 0, "up",
                    fifelse(directionally_consistent & logFC < 0, "down", "not_DE"))]
rna[, rna_class_FDRonly := fifelse(FDR < .05 & logFC > 0, "up", fifelse(FDR < .05 & logFC < 0, "down", "not_DE"))]
atac <- fread(sources[["ATAC_promoter_summary"]], colClasses = list(character = "gene_id"))
setnames(atac, "gene_id", "entrez")
links <- fread(sources[["ATAC_promoter_links"]], colClasses = list(character = "gene_id"))
atac_chrom <- links[, .(atac_chrom = paste(sort(unique(chrom)), collapse = ";")), by = .(entrez = gene_id)]
atac <- merge(atac, atac_chrom, by = "entrez")
atac[, atac_class := fifelse(gained_peaks > 0 & lost_peaks > 0, "mixed",
                     fifelse(gained_peaks > 0, "gained", fifelse(lost_peaks > 0, "lost", "ns")))]

# ---- Primary join: Ensembl -> Entrez, unique both ways ----
map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(rna$ensembl), keytype = "ENSEMBL", columns = "ENTREZID")
map <- as.data.table(map)[!is.na(ENTREZID)]
map[, n_entrez := uniqueN(ENTREZID), by = ENSEMBL]
map[, n_ensembl := uniqueN(ENSEMBL), by = ENTREZID]
audit <- rbind(
  data.table(step = "RNA expressed genes", n = nrow(rna)),
  data.table(step = "RNA genes with any Entrez mapping", n = uniqueN(map$ENSEMBL)),
  data.table(step = "RNA genes with one Entrez, not shared", n = map[n_entrez == 1 & n_ensembl == 1, uniqueN(ENSEMBL)]))
map1 <- map[n_entrez == 1 & n_ensembl == 1, .(ensembl = ENSEMBL, entrez = ENTREZID)]
j <- merge(rna, map1, by = "ensembl")
j <- merge(j, atac, by = "entrez")
j[, join_method := "ID"]
audit <- rbind(audit, data.table(step = "  ... with tested ATAC promoter peaks (ID join)", n = nrow(j)))
# Symbol fallback (author decision 22 Sep 2026): only for RNA genes without a one-to-one
# Ensembl->Entrez mapping; the ATAC symbol must be unambiguous and its Entrez not already
# claimed by the ID join. Chromosome agreement is enforced below for both methods.
amb_sym <- atac[, .N, by = symbol][N > 1, symbol]
fb <- rna[!ensembl %in% map1$ensembl & !gene_name %in% rna[duplicated(gene_name), gene_name]]
fb <- merge(fb, atac[!symbol %in% amb_sym & !entrez %in% j$entrez], by.x = "gene_name", by.y = "symbol")
fb[, `:=`(symbol = gene_name, join_method = "symbol_fallback")]
j <- rbind(j, fb, use.names = TRUE, fill = TRUE)
stopifnot(!anyDuplicated(j$entrez), !anyDuplicated(j$ensembl))
audit <- rbind(audit, data.table(step = "  ... added by symbol fallback", n = nrow(fb)))
j[, same_chrom := chrom == atac_chrom]
fwrite(j[same_chrom == FALSE, .(entrez, gene_name, symbol, chrom, atac_chrom)], file.path(out, "chromosome_mismatches.csv"))
j <- j[same_chrom == TRUE]
audit <- rbind(audit, data.table(step = "Eligible universe (same chromosome)", n = nrow(j)))
fwrite(audit, file.path(out, "mapping_audit.csv"))

joint_label <- function(r, a) fifelse(r == "not_DE", "RNA_not_DE",
  fifelse(a == "ns", "RNA_DE_ATAC_ns", fifelse(a == "mixed", "mixed_promoter",
  fifelse((r == "up" & a == "gained") | (r == "down" & a == "lost"), "concordant", "opposite"))))
j[, joint := joint_label(rna_class, atac_class)]
j[, joint_FDRonly := joint_label(rna_class_FDRonly, atac_class)]
fwrite(j[, .(entrez_id = entrez, ensembl, gene_name, atac_symbol = symbol, join_method, chrom, rna_logFC = logFC, rna_FDR = FDR,
             AUG_logFC, MAY_logFC, same_direction, rna_class, tested_promoter_peaks, atac_median_logFC = median_peak_logFC,
             atac_min_peak_FDR = minimum_peak_FDR, gained_peaks, lost_peaks, atac_class, joint, rna_class_FDRonly,
             joint_FDRonly, peak_ids)][order(joint, -abs(rna_logFC))],
       file.path(out, "RNA_ATAC_promoter_genes.csv"))

counts <- function(dt, rcol) {
  x <- dt[, .N, by = .(rna = get(rcol), atac = atac_class)]
  dcast(x, rna ~ atac, value.var = "N", fill = 0)[, rule := rcol]
}
fwrite(rbind(counts(j, "rna_class"), counts(j, "rna_class_FDRonly"), fill = TRUE), file.path(out, "joint_category_table.csv"))

sp <- function(x, y) if (length(x) > 2) unname(cor(x, y, method = "spearman")) else NA_real_
cors <- rbind(
  data.table(subset = "all eligible genes", n = nrow(j), spearman = sp(j$median_peak_logFC, j$logFC)),
  data.table(subset = "RNA DE (thesis rule)", n = j[rna_class != "not_DE", .N],
             spearman = j[rna_class != "not_DE", sp(median_peak_logFC, logFC)]),
  data.table(subset = "ATAC significant promoter (any)", n = j[atac_class != "ns", .N],
             spearman = j[atac_class != "ns", sp(median_peak_logFC, logFC)]),
  data.table(subset = "both significant (thesis rule)", n = j[rna_class != "not_DE" & atac_class != "ns", .N],
             spearman = j[rna_class != "not_DE" & atac_class != "ns", sp(median_peak_logFC, logFC)]))
fwrite(cors, file.path(out, "continuous_agreement.csv"))

# ---- Sensitivity: exact symbol join ----
amb <- atac[, .N, by = symbol][N > 1, symbol]
sj <- merge(rna, atac[!symbol %in% amb], by.x = "gene_name", by.y = "symbol")
sj[, joint := joint_label(rna_class, atac_class)]
sym_summary <- data.table(join = c("ID join + symbol fallback (primary)", "exact symbol only (sensitivity)"),
  universe = c(nrow(j), nrow(sj)),
  both_significant = c(j[rna_class != "not_DE" & atac_class != "ns", .N], sj[rna_class != "not_DE" & atac_class != "ns", .N]),
  concordant = c(j[joint == "concordant", .N], sj[joint == "concordant", .N]),
  opposite = c(j[joint == "opposite", .N], sj[joint == "opposite", .N]),
  mixed = c(j[joint == "mixed_promoter", .N], sj[joint == "mixed_promoter", .N]))
sym_summary[, shared_universe_genes := c(NA, length(intersect(j$entrez, sj$entrez)))]
fwrite(sym_summary, file.path(out, "join_sensitivity.csv"))

# ---- Branch comparison: thesis limma-voom vs earlier DESeq2 ----
ben <- fread(sources[["earlier_RNA_snapshot"]])
ben <- ben[!duplicated(gene_id) & !duplicated(gene_id, fromLast = TRUE)]
bc <- merge(rna[, .(ensembl, gene_name, thesis_logFC = logFC, thesis_FDR = FDR, thesis_class = rna_class)],
            ben[, .(ensembl = gene_id, earlier_logFC = logFC, earlier_FDR = FDR)], by = "ensembl")
bc[, earlier_sig := !is.na(earlier_FDR) & earlier_FDR < .05]
branch <- data.table(
  metric = c("genes in both results", "Spearman logFC", "Pearson logFC", "thesis DE (rule) genes in shared set",
             "earlier FDR<0.05 genes in shared set", "thesis DE also earlier FDR<0.05",
             "... of those, same direction", "earlier FDR<0.05 not thesis DE"),
  value = c(nrow(bc), bc[is.finite(earlier_logFC), cor(thesis_logFC, earlier_logFC, method = "spearman")],
            bc[is.finite(earlier_logFC), cor(thesis_logFC, earlier_logFC)], bc[thesis_class != "not_DE", .N],
            bc[earlier_sig == TRUE, .N], bc[thesis_class != "not_DE" & earlier_sig == TRUE, .N],
            bc[thesis_class != "not_DE" & earlier_sig == TRUE & sign(thesis_logFC) == sign(earlier_logFC), .N],
            bc[thesis_class == "not_DE" & earlier_sig == TRUE, .N]))
fwrite(branch, file.path(out, "RNA_branch_comparison.csv"))

old <- fread(sources[["earlier_joint_genes"]], colClasses = list(character = "entrez_id"))
old_status <- merge(old[, .(entrez_id, earlier_symbol = atac_symbol, earlier_rna_logFC = rna_logFC)],
                    j[, .(entrez_id = entrez, gene_name, rna_logFC = logFC, rna_FDR = FDR, rna_class, atac_class, joint)],
                    by = "entrez_id", all.x = TRUE)
old_status[is.na(joint), joint := "not_in_thesis_universe"]
fwrite(old_status, file.path(out, "earlier_45_joint_genes_status.csv"))

summary <- list(
  thesis_RNA_rule = "directionally_consistent: FDR < 0.05, |log2FC| >= 1, same AUG/MAY direction",
  universe = nrow(j), RNA_DE_in_universe = j[rna_class != "not_DE", .N],
  ATAC_significant_in_universe = j[atac_class != "ns", .N],
  joint = as.list(j[, .N, by = joint][, setNames(N, joint)]),
  joint_FDRonly = as.list(j[, .N, by = joint_FDRonly][, setNames(N, joint_FDRonly)]),
  spearman_all = cors$spearman[1], earlier_45_retained_both_significant = old_status[joint %in% c("concordant", "opposite", "mixed_promoter"), .N])
write_json(summary, file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
writeLines(c(
  "Thesis RNA: SEPT_11_2026_REANALYZE programme_analysis_v1 gene_consistency.csv (limma-voom ~ batch*treatment, combined contrast, 12 core libraries, GENCODE v49). No refit.",
  "Thesis RNA DE: directionally_consistent (FDR < 0.05, |log2FC| >= 1, same AUG_2025/MAY_2026 direction). FDR-only classes are a sensitivity.",
  "ATAC: primary edgeR promoter links (TSS +/- 2 kb, knownGene Entrez). Gene pattern: gained/lost/mixed/ns from significant (FDR < 0.05) promoter peaks; median of all tested promoter peak logFC. Peak FDRs are not gene-level tests.",
  "Primary join: unversioned RNA Ensembl ID -> Entrez via org.Hs.eg.db (one-to-one only), plus exact-symbol fallback for RNA genes without a one-to-one ID mapping when the ATAC symbol is unambiguous and its Entrez is not already claimed; same chromosome required. join_method records which applied. Exact-symbol-only join reported as sensitivity.",
  "Concordant = RNA up with only gained promoter peaks, or RNA down with only lost; opposite = the reverse; mixed = RNA DE with both gained and lost promoter peaks.",
  "Counts and Spearman correlations are descriptive group-level comparisons; RNA and ATAC specimens are not established as matched; no overlap P-value is assigned (genes differ in peak number; observations are dependent).",
  "RNA_branch_comparison.csv compares the thesis RNA result with the earlier DESeq2 ~ batch + treatment branch (EnsDb v86) used by integration/RNAseq, joined on unversioned Ensembl IDs."),
  file.path(out, "METHODS.txt"))
print(audit); print(fread(file.path(out, "joint_category_table.csv"))); print(cors); print(sym_summary); print(branch)
cat("Earlier 45 joint genes:\n"); print(old_status[, .N, by = joint])
