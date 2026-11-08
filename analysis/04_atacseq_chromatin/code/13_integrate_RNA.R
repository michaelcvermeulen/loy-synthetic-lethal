options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(jsonlite); library(digest)})
outdir <- "integration/RNAseq"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create("config/RNAseq", showWarnings = FALSE)
source_root <- "A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080_BATCHES"
sources <- c(result = file.path(source_root, "DE/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv"),
  notebook = file.path(source_root, "analysis_HT1080_BATCHES.Rmd"),
  metadata = file.path(source_root, "count_matrix/HT1080_BATCHES_sample_metadata.csv"))
snapshots <- file.path("config/RNAseq", basename(sources))
stopifnot(all(file.exists(sources)), all(file.copy(sources, snapshots, overwrite = TRUE)))
hashes <- vapply(sources, digest, character(1), algo = "sha256", file = TRUE)
stopifnot(identical(unname(hashes), unname(vapply(snapshots, digest, character(1), algo = "sha256", file = TRUE))))
fwrite(data.table(role = names(sources), source = sources, snapshot = snapshots, sha256 = hashes),
       file.path(outdir, "source_manifest.csv"))
notebook <- readLines(snapshots[2], warn = FALSE)
stopifnot(any(grepl("design    = ~ batch + treatment", notebook, fixed = TRUE)),
          any(grepl('contrast = c("treatment", "LOY", "CTRL")', notebook, fixed = TRUE)))
rna <- fread(snapshots[1], colClasses = list(character = c("entrezid", "gene_id", "gene_id_full")))
metadata <- fread(snapshots[3])
stopifnot(!anyDuplicated(rna$gene_id_full), !anyDuplicated(rna$gene_id),
          sum(metadata$treatment == "CTRL") == 6, sum(metadata$treatment == "LOY") == 6,
          sum(metadata$treatment == "PARENTAL") == 3)
rna[, entrezid := trimws(entrezid)]
rna[entrezid %in% c("", "NA"), entrezid := NA_character_]
rna[, rna_chrom := paste0("chr", sub("^chr", "", chr))]
par <- fread("config/hg38_PAR.bed", col.names = c("chrom", "start", "end", "name"))
rna[, rna_PAR := FALSE]
for (i in seq_len(nrow(par))) {
  rna[rna_chrom == par$chrom[i] & gene_start <= par$end[i] & gene_end > par$start[i], rna_PAR := TRUE]
}
rna[, entrez_multiplicity := .N, by = entrezid]
unique_rna <- rna[!is.na(entrezid) & entrez_multiplicity == 1]
ambiguous <- rna[!is.na(entrezid) & entrez_multiplicity > 1]
fwrite(ambiguous, file.path(outdir, "ambiguous_RNA_Entrez_mappings.csv"))
genes <- fread("annotation/promoter_gene_summary.csv", colClasses = list(character = "gene_id"))
setnames(genes, c("gene_id", "symbol", "median_peak_logFC", "minimum_peak_FDR"),
         c("entrez_id", "atac_symbol", "atac_median_promoter_logFC", "atac_minimum_peak_FDR"))
links <- fread("annotation/promoter_peak_gene_links.csv", colClasses = list(character = "gene_id"))
chroms <- unique(links[, .(gene_id, chrom)])
stopifnot(!anyDuplicated(genes$entrez_id), !anyDuplicated(chroms$gene_id))
genes[, atac_chrom := chroms$chrom[match(entrez_id, chroms$gene_id)]]
idx <- match(genes$entrez_id, unique_rna$entrezid)
fields <- c(rna_ensembl_id = "gene_id", rna_gene_id_full = "gene_id_full", rna_symbol = "gene_name",
  rna_logFC = "logFC", rna_FDR = "FDR", rna_PValue = "PValue", rna_lfcSE = "lfcSE",
  rna_baseMean = "logCPM", rna_chrom = "rna_chrom", rna_PAR = "rna_PAR")
for (field in names(fields)) set(genes, j = field, value = unique_rna[[fields[[field]]]][idx])
genes[, mapping_status := fcase(entrez_id %in% ambiguous$entrezid, "ambiguous_RNA_Entrez",
  is.na(rna_ensembl_id), "not_in_RNA_result",
  atac_chrom != rna_chrom, "chromosome_mismatch",
  rna_chrom == "chrY" | rna_PAR, "RNA_chrY_or_PAR",
  !is.finite(rna_FDR) | !is.finite(rna_PValue) | !is.finite(rna_logFC), "RNA_statistics_unavailable",
  default = "eligible")]
genes[, rna_significant_FDR005 := !is.na(rna_FDR) & rna_FDR < .05]
genes[, rna_thesis_effect_subset := rna_significant_FDR005 & abs(rna_logFC) > .75]
genes[, atac_any_significant_promoter_peak := gained_peaks > 0 | lost_peaks > 0]
genes[, joint_category := fcase(mapping_status != "eligible", "ineligible",
  !rna_significant_FDR005 | !atac_any_significant_promoter_peak, "not_jointly_significant",
  promoter_pattern == "mixed", "mixed_promoter",
  promoter_pattern == "gained" & rna_logFC > 0, "concordant_gain",
  promoter_pattern == "lost" & rna_logFC < 0, "concordant_loss",
  default = "discordant")]
fwrite(genes, file.path(outdir, "gene_mapping_audit.csv"))
eligible <- genes[mapping_status == "eligible"]
stopifnot(nrow(eligible) > 1000, !anyDuplicated(eligible$entrez_id),
          all(eligible$rna_FDR >= 0 & eligible$rna_FDR <= 1),
          all(eligible$atac_chrom == eligible$rna_chrom), !any(eligible$rna_PAR))
fwrite(eligible, file.path(outdir, "ATAC_RNA_gene_universe.csv"))
joint <- eligible[atac_any_significant_promoter_peak & rna_significant_FDR005]
setorder(joint, rna_FDR, atac_minimum_peak_FDR, entrez_id)
fwrite(joint, file.path(outdir, "ATAC_RNA_joint_significant_genes.csv"))
fwrite(joint[rna_thesis_effect_subset == TRUE], file.path(outdir, "ATAC_RNA_joint_genes_RNA_abs_logFC_gt075.csv"))
linked <- cbind(links, genes[match(links$gene_id, genes$entrez_id),
  c(names(fields), "mapping_status", "joint_category"), with = FALSE])
stopifnot(nrow(linked) == nrow(links), !anyDuplicated(linked[, .(peak_id, gene_id)]))
fwrite(linked, file.path(outdir, "ATAC_RNA_promoter_peak_links.csv"))
# Distal nearest-gene assignments are retained as a separate exploratory table.
peaks <- fread("annotation/annotated_differential_peaks.csv",
               colClasses = list(character = "nearest_entrez_id"))
nearest <- peaks[, .(peak_id, chrom, start, end, annotation, nearest_entrez_id, nearest_symbol,
  nearest_tss_distance_signed, atac_logFC = logFC, atac_FDR = FDR)]
nidx <- match(nearest$nearest_entrez_id, unique_rna$entrezid)
for (field in c("rna_ensembl_id", "rna_symbol", "rna_logFC", "rna_FDR", "rna_chrom")) {
  set(nearest, j = field, value = unique_rna[[fields[[field]]]][nidx])
}
nearest[, link_evidence := "nearest_TSS_positional_hypothesis"]
nearest[, same_chromosome := !is.na(rna_chrom) & chrom == rna_chrom]
fwrite(nearest, file.path(outdir, "ATAC_RNA_nearest_gene_links_exploratory.csv"))
mapping_summary <- genes[, .(genes = .N), by = mapping_status]
joint_summary <- eligible[, .(genes = .N), by = joint_category]
fwrite(mapping_summary, file.path(outdir, "mapping_summary.csv"))
fwrite(joint_summary, file.path(outdir, "joint_category_counts.csv"))
rho <- cor(eligible$atac_median_promoter_logFC, eligible$rna_logFC, method = "spearman")
palette <- c(not_jointly_significant = "grey78", concordant_gain = "#D55E00",
  concordant_loss = "#0072B2", discordant = "#CC79A7", mixed_promoter = "#009E73")
p <- ggplot(eligible, aes(atac_median_promoter_logFC, rna_logFC)) +
  geom_point(colour = "grey78", alpha = .3, size = .6) +
  geom_point(data = joint, aes(colour = joint_category), size = 1.8, alpha = .85) +
  scale_colour_manual(values = palette) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = 2) + theme_bw(base_size = 12) +
  labs(x = "Median ATAC promoter log2 fold change (LOY / CTRL)",
       y = "RNA log2 fold change (LOY / CTRL), batch + treatment", colour = NULL,
       subtitle = sprintf("%s eligible genes; descriptive Spearman rho = %.3f", format(nrow(eligible), big.mark = ","), rho),
       caption = "Colored genes: RNA FDR <0.05 and at least one ATAC promoter peak FDR <0.05. ATAC significance remains peak-level.")
ggsave(file.path(outdir, "ATAC_RNA_effects.png"), p, width = 11, height = 7, dpi = 200)
ggsave(file.path(outdir, "ATAC_RNA_effects.pdf"), p, width = 11, height = 7)
summary <- list(rna_model = "DESeq2 ~ batch + treatment", contrast = "LOY / CTRL",
  user_selection = "HT1080 combined-batch RNAseq used in Ben's thesis; batch + treatment",
  selection_evidence = "User specification and primary-model declaration in the source notebook; thesis manuscript itself not located",
  rna_source = unname(sources[1]), rna_source_sha256 = unname(hashes[1]),
  rna_rows = nrow(rna), rna_finite_FDR_rows = sum(is.finite(rna$FDR)),
  rna_contrast_samples = list(CTRL = 6, LOY = 6), parental_samples_in_model = 3,
  rna_annotation = "EnsDb.Hsapiens.v86 annotation in source notebook; source counts documented as GENCODE v49",
  rna_source_logCPM_column = "Actually DESeq2 baseMean, renamed by source notebook; preserved here as rna_baseMean",
  tested_ATAC_promoter_genes = nrow(genes), common_eligible_genes = nrow(eligible),
  joint_significant_genes = nrow(joint), concordant_gain = sum(joint$joint_category == "concordant_gain"),
  concordant_loss = sum(joint$joint_category == "concordant_loss"), discordant = sum(joint$joint_category == "discordant"),
  mixed_promoter = sum(joint$joint_category == "mixed_promoter"),
  joint_genes_RNA_abs_logFC_gt075 = sum(joint$rna_thesis_effect_subset), descriptive_spearman = rho,
  sample_matching = "Group-effect integration; direct ATAC/RNA sample matching not established",
  lpWGS_status = "Pending dataset selection and assay matching; no copy-number adjustment applied")
write_json(summary, file.path(outdir, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(c("RNA selection: HT1080 combined-batch DESeq2 ~ batch + treatment, per user. Source notebook verifies model and LOY/CTRL contrast.",
  "The thesis manuscript was not located; its identity is supplied by the user's specification, not independently verified against the manuscript.",
  "Existing RNA effect estimates, standard errors, P-values and FDR are reused unchanged; no RNA model refit.",
  "Six CTRL and six LOY from AUG_2025/MAY_2026; three parental samples also participate in original model estimation.",
  "RNA source annotation is EnsDb.Hsapiens.v86, not the earlier voom v113 candidate. Source logCPM is actually DESeq2 baseMean.",
  "Join uses existing Entrez IDs. Multiple RNA rows per Entrez are excluded from the gene-level comparison and written to an ambiguity table.",
  "Common universe: tested ATAC promoter genes with one RNA match, same chromosome, non-chrY/non-PAR RNA coordinates, and finite RNA statistics.",
  "RNA gene coordinates are 1-based closed; hg38 PAR BED is 0-based half-open. No symbol-only fallback is used.",
  "RNA significance: existing FDR <0.05. An additional table applies the notebook's |log2FC| >0.75 effect subset; this is not the thresholded Wald/TREAT-like result.",
  "ATAC gene summaries use the median of all tested promoter peak effects. Any significant promoter peak flags a gene; minimum peak FDR is not gene-level FDR.",
  "Concordance uses gained/lost significant-promoter patterns and RNA direction; mixed patterns remain separate. All promoter links are retained.",
  "Nearest-TSS peak links are positional hypotheses, saved separately. Shared gene effects and overlap counts are descriptive; no causal or sample-wise correlation claims.",
  "No formal overlap-enrichment P-value is assigned: genes differ in promoter peak count and genomic observations are dependent.",
  "lpWGS matching and copy-number context remain pending."), file.path(outdir, "METHODS.txt"))
writeLines(c("PASS: source snapshots match SHA256 hashes", "PASS: source model/contrast and 6 CTRL + 6 LOY + 3 parental metadata checked",
  "PASS: unique RNA Ensembl rows and unique gene-universe Entrez IDs", "PASS: ambiguous mapping accounting, chromosome agreement and PAR exclusions",
  "PASS: finite RNA statistics and unchanged one-row-per-promoter-link accounting"), file.path(outdir, "verification.txt"))
writeLines(capture.output(sessionInfo()), "logs/RNA_integration_session.txt")
print(summary)
