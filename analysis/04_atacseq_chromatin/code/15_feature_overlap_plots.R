options(warn = 1)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(ggplot2)})
outdir <- "integration/feature_overlaps"
dir.create(file.path(outdir, "genes"), recursive = TRUE, showWarnings = FALSE)
peaks <- fread("annotation/annotated_differential_peaks.csv", colClasses = list(character = "nearest_entrez_id"))
genes <- fread("config/genes_knownGene.tsv", colClasses = list(character = "gene_id"))
promoter <- fread("annotation/promoter_peak_gene_links.csv", colClasses = list(character = "gene_id"))
stopifnot(!anyDuplicated(peaks$peak_id), !anyDuplicated(genes$gene_id),
          all(peaks$annotation %in% c("promoter", "gene_body", "intergenic")))
promoter_links <- promoter[, .(feature = "promoter", peak_id, chrom, entrez_id = gene_id,
  atac_logFC = logFC, atac_FDR = FDR, link_evidence = "promoter_overlap_TSS_plus_minus_2kb")]
body <- peaks[annotation == "gene_body"]
hits <- findOverlaps(GRanges(body$chrom, IRanges(body$start+1L, body$end)),
  GRanges(genes$chrom, IRanges(genes$start, genes$end)), ignore.strand = TRUE)
body_links <- body[queryHits(hits), .(feature = "gene_body", peak_id, chrom,
  atac_logFC = logFC, atac_FDR = FDR)]
body_links[, `:=`(entrez_id = genes$gene_id[subjectHits(hits)], link_evidence = "gene_body_overlap")]
distal_links <- peaks[annotation == "intergenic", .(feature = "intergenic", peak_id, chrom,
  entrez_id = nearest_entrez_id, atac_logFC = logFC, atac_FDR = FDR,
  link_evidence = "nearest_TSS_positional_hypothesis")]
links <- unique(rbindlist(list(promoter_links, body_links, distal_links), use.names = TRUE))
stopifnot(!anyDuplicated(links[, .(feature, peak_id, entrez_id)]),
          setequal(links$peak_id, peaks$peak_id))
manifest <- fread("integration/RNAseq/source_manifest.csv")
rna <- fread(manifest[role == "result", snapshot], colClasses = list(character = c("entrezid", "gene_id")))
rna[, entrezid := trimws(entrezid)]
rna[entrezid %in% c("", "NA"), entrezid := NA_character_]
rna[, multiplicity := .N, by = entrezid]
ambiguous <- rna[!is.na(entrezid) & multiplicity > 1, unique(entrezid)]
rna <- rna[!is.na(entrezid) & multiplicity == 1]
rna[, rna_chrom := paste0("chr", sub("^chr", "", chr))]
par <- fread("config/hg38_PAR.bed", col.names = c("chrom", "start", "end", "name"))
rna[, rna_PAR := FALSE]
for (i in seq_len(nrow(par))) {
  rna[rna_chrom == par$chrom[i] & gene_start <= par$end[i] & gene_end > par$start[i], rna_PAR := TRUE]
}
idx <- match(links$entrez_id, rna$entrezid)
links[, `:=`(rna_ensembl_id = rna$gene_id[idx], rna_symbol = rna$gene_name[idx],
  rna_logFC = rna$logFC[idx], rna_FDR = rna$FDR[idx], rna_PValue = rna$PValue[idx],
  rna_chrom = rna$rna_chrom[idx], rna_PAR = rna$rna_PAR[idx])]
links[, mapping_status := fcase(entrez_id %in% ambiguous, "ambiguous_RNA_Entrez",
  is.na(rna_ensembl_id), "not_in_RNA_result", chrom != rna_chrom, "chromosome_mismatch",
  rna_chrom == "chrY" | rna_PAR, "RNA_chrY_or_PAR",
  !is.finite(rna_FDR) | !is.finite(rna_PValue) | !is.finite(rna_logFC), "RNA_statistics_unavailable",
  default = "eligible")]
fwrite(links, file.path(outdir, "peak_gene_link_audit.csv"))
fwrite(links[, .(links = .N, distinct_peaks = uniqueN(peak_id), distinct_genes = uniqueN(entrez_id)),
  by = .(feature, mapping_status)], file.path(outdir, "mapping_summary.csv"))
universe <- links[mapping_status == "eligible", .(
  rna_ensembl_id = unique(rna_ensembl_id), rna_symbol = unique(rna_symbol),
  rna_logFC = unique(rna_logFC), rna_FDR = unique(rna_FDR),
  minimum_peak_FDR = min(atac_FDR), tested_peaks = uniqueN(peak_id),
  peak_ids = paste(unique(peak_id), collapse = ";"), link_evidence = unique(link_evidence)
), by = .(feature, entrez_id)]
stopifnot(!anyDuplicated(universe[, .(feature, entrez_id)]))
fwrite(universe, file.path(outdir, "feature_gene_universes.csv"))
rows <- list()
for (scheme in c("matched_FDR", "RNA_fixed_005")) {
  for (feature_name in c("promoter", "gene_body", "intergenic")) {
    g <- universe[feature == feature_name]
    for (alpha in c(.05, .10, .20)) {
      rna_alpha <- if (scheme == "matched_FDR") alpha else .05
      a <- g$minimum_peak_FDR < alpha
      b <- g$rna_FDR < rna_alpha
      joint <- g[a & b]
      setorder(joint, rna_FDR, minimum_peak_FDR, entrez_id)
      fwrite(joint, file.path(outdir, "genes", sprintf("%s_%s_ATAC_%.2f_RNA_%.2f.csv",
        scheme, feature_name, alpha, rna_alpha)))
      row <- data.table(scheme, feature = feature_name, ATAC_FDR = alpha, RNA_FDR = rna_alpha,
        universe = nrow(g), ATAC_genes = sum(a), RNA_genes = sum(b),
        ATAC_only = sum(a & !b), overlap = sum(a & b), RNA_only = sum(!a & b), neither = sum(!a & !b))
      stopifnot(row$ATAC_only + row$overlap == row$ATAC_genes,
        row$RNA_only + row$overlap == row$RNA_genes,
        row$ATAC_only + row$overlap + row$RNA_only + row$neither == row$universe)
      rows[[length(rows)+1L]] <- row
    }
  }
}
summary <- rbindlist(rows)
previous <- fread("integration/FDR_sensitivity/ATAC_RNA_threshold_grid.csv")
for (i in seq_len(nrow(summary[feature == "promoter"]))) {
  x <- summary[feature == "promoter"][i]
  stopifnot(x$overlap == previous[ATAC_FDR == x$ATAC_FDR & RNA_FDR == x$RNA_FDR, joint_genes],
            x$universe == 11898)
}
stopifnot(all(summary[order(ATAC_FDR), all(diff(overlap) >= 0), by = .(scheme, feature)]$V1))
fwrite(summary, file.path(outdir, "overlap_counts.csv"))
feature_labels <- c(promoter = "Promoter\nTSS overlap", gene_body = "Gene body\nGene-body overlap",
                    intergenic = "Intergenic\nNearest TSS (exploratory)")
make_plot <- function(df, heading) {
  df <- copy(df)
  df[, feature_label := factor(feature, levels = names(feature_labels), labels = feature_labels)]
  df[, threshold := factor(sprintf("ATAC <%.2f | RNA <%.2f", ATAC_FDR, RNA_FDR),
                           levels = unique(sprintf("ATAC <%.2f | RNA <%.2f", ATAC_FDR, RNA_FDR)))]
  angle <- seq(0, 2*pi, length.out = 200)
  circles <- rbindlist(lapply(seq_len(nrow(df)), function(i) {
    rbindlist(lapply(c("ATAC", "RNAseq"), function(assay) {
      data.table(feature_label = df$feature_label[i], threshold = df$threshold[i], assay,
        x = (if (assay == "ATAC") 3.7 else 6.3) + 2.25*cos(angle), y = 3.7 + 2.25*sin(angle))
    }))
  }))
  labels <- rbindlist(lapply(seq_len(nrow(df)), function(i) {
    data.table(feature_label = df$feature_label[i], threshold = df$threshold[i],
      x = c(2.9, 5, 7.1), y = 3.7,
      label = format(c(df$ATAC_only[i], df$overlap[i], df$RNA_only[i]), big.mark = ",", trim = TRUE))
  }))
  plot <- ggplot() + geom_polygon(data = circles, aes(x, y, fill = assay, group = assay),
      alpha = .27, colour = "grey45", linewidth = .45) +
    geom_text(data = labels, aes(x, y, label = label), size = 4.5, fontface = "bold") +
    annotate("text", x = 2.6, y = 6.5, label = "ATAC-linked\ngenes", size = 3.4, lineheight = .95) +
    annotate("text", x = 7.4, y = 6.5, label = "RNAseq\ngenes", size = 3.4, lineheight = .95) +
    geom_text(data = df, aes(x = 5, y = .65, label = paste0("Eligible genes: ", format(universe, big.mark = ",", trim = TRUE),
      "  |  Neither: ", format(neither, big.mark = ",", trim = TRUE))), size = 2.9) +
    scale_fill_manual(values = c(ATAC = "#0072B2", RNAseq = "#D55E00")) +
    facet_grid(feature_label ~ threshold) + coord_equal(xlim = c(0, 10), ylim = c(0, 7.3), expand = FALSE) +
    theme_void(base_size = 12) + theme(legend.position = "none", strip.text = element_text(size = 11),
      strip.background = element_rect(fill = "#edf3f5", colour = NA),
      plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 10, hjust = 0), plot.margin = margin(12, 15, 12, 15),
      panel.spacing = grid::unit(.8, "lines")) +
    labs(title = heading, subtitle = "Distinct genes; any qualifying ATAC peak or RNA change, regardless of direction",
      caption = "Circle areas are schematic. RNA is restricted to each feature's eligible gene universe.\nATAC FDR applies to peaks; overlap does not have a joint gene-level FDR. Genes can occur in multiple feature rows.")
  plot
}
for (scheme_name in unique(summary$scheme)) {
  df <- summary[scheme == scheme_name]
  heading <- if (scheme_name == "matched_FDR") "Genomic feature overlaps: matched ATAC and RNA FDR" else
    "Genomic feature overlaps: RNA FDR fixed at 0.05"
  p <- make_plot(df, heading)
  ggsave(file.path(outdir, paste0(scheme_name, "_overlap_grid.png")), p, width = 13, height = 11, dpi = 220, bg = "white")
  ggsave(file.path(outdir, paste0(scheme_name, "_overlap_grid.pdf")), p, width = 13, height = 11, bg = "white")
  for (feature_name in names(feature_labels)) {
    p <- make_plot(df[feature == feature_name], heading)
    stem <- paste0(scheme_name, "_", feature_name, "_overlaps")
    ggsave(file.path(outdir, paste0(stem, ".png")), p, width = 13, height = 4.8, dpi = 220, bg = "white")
    ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = 13, height = 4.8, bg = "white")
  }
}
writeLines(c("Feature definitions preserve the existing exclusive peak annotation: promoter, then gene body, then intergenic.",
  "Promoter links retain all TSS +/-2 kb overlaps; gene-body links retain all overlapping knownGene gene bodies; intergenic links use the existing nearest-TSS assignment as an exploratory hypothesis.",
  "Circles count unique Entrez genes, not peaks versus genes. A gene enters the ATAC set if any linked peak passes ATAC FDR.",
  "Per-feature universe: genes linked to any tested peak of that feature, uniquely mapped to the selected RNA table, with chromosome agreement, finite RNA statistics and no RNA chrY/PAR overlap.",
  "RNA-only counts include only RNA-significant genes in this feature-specific universe. Neither is printed below each diagram.",
  "ATAC peak classes are mutually exclusive; genes may link to peaks in multiple classes. Feature rows cannot be summed as disjoint gene sets.",
  "Main grid: both assays use FDR <0.05, <0.10, <0.20. Supplementary grid: RNA remains <0.05 while ATAC changes.",
  "No logFC cutoff or direction filter is imposed. Original primary model and adjusted P-values are unchanged.",
  "Venn circle areas are schematic, not proportional to set size. Counts satisfy all set-partition identities.",
  "Overlap does not have a calibrated joint gene-level FDR. Intergenic nearest-TSS links do not establish enhancer targets.",
  "All promoter overlaps reproduce the previous threshold-grid results. Full link audit, universes, counts and overlapping gene lists are saved."),
  file.path(outdir, "METHODS.txt"))
writeLines(c("PASS: exclusive peak classes and complete peak-to-gene link accounting",
  "PASS: unique feature/gene universes and set partition identities for all 18 panels",
  "PASS: promoter universe and overlaps reproduce previous FDR sensitivity results",
  "PASS: overlap counts increase monotonically with relaxed thresholds"), file.path(outdir, "verification.txt"))
print(summary)
