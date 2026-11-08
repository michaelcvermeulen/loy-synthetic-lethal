script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages(library(data.table))
outdir <- "integration/FDR_sensitivity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
peaks <- fread("differential/primary_TMM_all.csv")
links <- fread("integration/RNAseq/ATAC_RNA_promoter_peak_links.csv", colClasses = list(character = "gene_id"))
universe <- fread("integration/RNAseq/ATAC_RNA_gene_universe.csv", colClasses = list(character = "entrez_id"))
stopifnot(!anyDuplicated(universe$entrez_id))
links <- links[gene_id %in% universe$entrez_id]
cutoffs <- c(.05, .10, .20)
peak_summary <- rbindlist(lapply(cutoffs, function(alpha) {
  data.table(ATAC_FDR = alpha, tested_peaks = nrow(peaks),
    gained = sum(peaks$FDR < alpha & peaks$logFC > 0),
    lost = sum(peaks$FDR < alpha & peaks$logFC < 0),
    total = sum(peaks$FDR < alpha),
    TREAT_FDR_lfc05 = sum(peaks$TREAT_FDR_lfc0.5 < alpha))
}))
rows <- list()
for (atac_cutoff in cutoffs) {
  counts <- links[, .(gained = sum(FDR < atac_cutoff & logFC > 0),
                       lost = sum(FDR < atac_cutoff & logFC < 0)), by = gene_id]
  g <- copy(universe)
  g[, gained_at_cutoff := counts$gained[match(entrez_id, counts$gene_id)]]
  g[, lost_at_cutoff := counts$lost[match(entrez_id, counts$gene_id)]]
  stopifnot(!anyNA(g$gained_at_cutoff), !anyNA(g$lost_at_cutoff))
  for (rna_cutoff in cutoffs) {
    joint <- g[(gained_at_cutoff > 0 | lost_at_cutoff > 0) & rna_FDR < rna_cutoff]
    joint[, category_at_cutoff := fcase(
      gained_at_cutoff > 0 & lost_at_cutoff > 0, "mixed_promoter",
      gained_at_cutoff > 0 & rna_logFC > 0, "concordant_gain",
      lost_at_cutoff > 0 & rna_logFC < 0, "concordant_loss", default = "discordant")]
    # Exclude fixed-0.05 flags from these threshold-specific exports.
    export <- joint[, .(entrez_id, rna_symbol, atac_median_promoter_logFC, atac_minimum_peak_FDR,
      rna_logFC, rna_FDR, gained_at_cutoff, lost_at_cutoff, category_at_cutoff,
      RNA_abs_logFC_gt075 = abs(rna_logFC) > .75, peak_ids)]
    export[, `:=`(ATAC_FDR_cutoff = atac_cutoff, RNA_FDR_cutoff = rna_cutoff)]
    fwrite(export, file.path(outdir, sprintf("joint_genes_ATAC_%.2f_RNA_%.2f.csv", atac_cutoff, rna_cutoff)))
    rows[[length(rows)+1L]] <- data.table(ATAC_FDR = atac_cutoff, RNA_FDR = rna_cutoff,
      eligible_genes = nrow(g), ATAC_promoter_genes = sum(g$gained_at_cutoff > 0 | g$lost_at_cutoff > 0),
      RNA_genes = sum(g$rna_FDR < rna_cutoff), joint_genes = nrow(joint),
      concordant_gain = sum(joint$category_at_cutoff == "concordant_gain"),
      concordant_loss = sum(joint$category_at_cutoff == "concordant_loss"),
      discordant = sum(joint$category_at_cutoff == "discordant"),
      mixed_promoter = sum(joint$category_at_cutoff == "mixed_promoter"),
      joint_RNA_abs_logFC_gt075 = sum(abs(joint$rna_logFC) > .75))
  }
}
joint_summary <- rbindlist(rows)
stopifnot(peak_summary[ATAC_FDR == .05, total] == 758,
          joint_summary[ATAC_FDR == .05 & RNA_FDR == .05, joint_genes] == 45,
          all(diff(peak_summary$total) >= 0),
          all(joint_summary$joint_genes == with(joint_summary,
            concordant_gain + concordant_loss + discordant + mixed_promoter)))
fwrite(peak_summary, file.path(outdir, "ATAC_peak_summary.csv"))
fwrite(joint_summary, file.path(outdir, "ATAC_RNA_threshold_grid.csv"))
writeLines(c("Sensitivity analysis only. Existing primary model, tested universe and adjusted P-values remain unchanged.",
  "All cutoffs use strict FDR < cutoff. Evaluate ATAC 0.05/0.10/0.20 crossed with RNA 0.05/0.10/0.20.",
  "ATAC-promoter gene calls are recomputed from all promoter peak FDRs at each cutoff; they are not gene-level FDR tests.",
  "Concordance uses significant-promoter directions at the selected cutoff and the RNA effect direction; mixed promoter patterns remain separate.",
  "The RNA |log2FC| >0.75 column is an additional descriptive effect-size subset, not a thresholded Wald test.",
  "The FDR threshold is not the error probability of each hit. Separate assay thresholds do not calibrate a joint gene-list FDR.",
  "Checks: 0.05 reproduces 758 ATAC peaks and 45 joint genes; increasing ATAC cutoff expands peak calls; gene-category totals agree."),
  file.path(outdir, "METHODS.txt"))
print(peak_summary)
print(joint_summary)
