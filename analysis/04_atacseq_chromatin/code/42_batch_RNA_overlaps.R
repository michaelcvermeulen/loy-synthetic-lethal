options(warn = 1)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = '/', mustWork = TRUE))))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel); library(digest); library(jsonlite)})
out <- 'integration/batch_RNA_overlaps'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, 'sources'), showWarnings = FALSE)
src <- 'A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080_BATCHES'
batches <- c('AUG_2025', 'MAY_2026')
files <- c(setNames(file.path(src, 'DE', paste0('ht1080_', batches, '_LOY_v_CTRL_DESeq2.csv')), batches),
  notebook = file.path(src, 'analysis_HT1080_BATCHES.Rmd'),
  metadata = file.path(src, 'count_matrix/HT1080_BATCHES_sample_metadata.csv'))
snap <- file.path(out, 'sources', basename(files))
stopifnot(all(file.copy(files, snap, overwrite = TRUE)))
hash <- vapply(files, digest, character(1), algo = 'sha256', file = TRUE)
stopifnot(identical(unname(hash), unname(vapply(snap, digest, character(1), algo = 'sha256', file = TRUE))))
manifest <- data.table(role = names(files), source = files, snapshot = snap, sha256 = hash)
fwrite(manifest, file.path(out, 'source_manifest.csv'))
nb <- readLines(snap[3], warn = FALSE)
start <- grep('run_batch_de <- function', nb, fixed = TRUE)
block <- nb[start:(start + 24)]
stopifnot(any(grepl('design    = ~ treatment', block, fixed = TRUE)),
  any(grepl('c("CTRL", "LOY")', block, fixed = TRUE)),
  any(grepl('contrast = c("treatment", "LOY", "CTRL")', block, fixed = TRUE)))
meta <- fread(snap[4])[treatment %in% c('CTRL', 'LOY')]
stopifnot(all(meta[, .N, by = .(batch, treatment)]$N == 3))
fwrite(meta, file.path(out, 'contrast_samples.csv'))
# Reuse every original ATAC link, including links ineligible for the combined RNA analysis.
links <- fread('integration/feature_overlaps/peak_gene_link_audit.csv', colClasses = list(character = 'entrez_id'))
links <- unique(links[, .(feature, peak_id, chrom, entrez_id, atac_logFC, atac_FDR, link_evidence)])
stopifnot(!anyDuplicated(links[, .(feature, peak_id, entrez_id)]), all(is.finite(links$atac_FDR)))
atac <- links[, .(chrom = unique(chrom), atac_median_logFC = median(atac_logFC),
  minimum_peak_FDR = min(atac_FDR), tested_peaks = uniqueN(peak_id),
  peak_ids = paste(peak_id, collapse = ';'), link_evidence = unique(link_evidence),
  significant_up_peaks = sum(atac_FDR < .1 & atac_logFC > 0),
  significant_down_peaks = sum(atac_FDR < .1 & atac_logFC < 0)), by = .(feature, entrez_id)]
stopifnot(!anyDuplicated(atac[, .(feature, entrez_id)]))
par <- fread('config/hg38_PAR.bed', col.names = c('chrom', 'start', 'end', 'name'))
all <- list(); audits <- list(); rna_audits <- list()
for (b in batches) {
  r <- fread(manifest[role == b, snapshot], colClasses = list(character = c('entrezid', 'gene_id', 'gene_id_full')))
  stopifnot(all(r$batch_fit == b), !anyDuplicated(r$gene_id), !anyDuplicated(r$gene_id_full))
  r[, entrezid := trimws(entrezid)]
  r[entrezid %in% c('', 'NA'), entrezid := NA_character_]
  r[, multiplicity := .N, by = entrezid]
  r[, rna_chrom := paste0('chr', sub('^chr', '', chr))]
  r[, rna_PAR := FALSE]
  for (i in seq_len(nrow(par))) r[rna_chrom == par$chrom[i] & gene_start <= par$end[i] & gene_end > par$start[i], rna_PAR := TRUE]
  rna_audits[[b]] <- r[, .(batch = b, gene_id, entrezid, gene_name, multiplicity,
    missing_entrez = is.na(entrezid), finite_statistics = is.finite(FDR) & is.finite(PValue) & is.finite(logFC), rna_chrom, rna_PAR)]
  amb <- r[!is.na(entrezid) & multiplicity > 1, unique(entrezid)]
  r <- r[!is.na(entrezid) & multiplicity == 1]
  g <- copy(atac); idx <- match(g$entrez_id, r$entrezid)
  g[, `:=`(batch = b, rna_ensembl_id = r$gene_id[idx], rna_symbol = r$gene_name[idx],
    rna_logFC = r$logFC[idx], rna_FDR = r$FDR[idx], rna_PValue = r$PValue[idx],
    rna_chrom = r$rna_chrom[idx], rna_PAR = r$rna_PAR[idx])]
  g[, mapping_status := fcase(entrez_id %in% amb, 'ambiguous_RNA_Entrez',
    is.na(rna_ensembl_id), 'not_in_RNA_result', chrom != rna_chrom, 'chromosome_mismatch',
    rna_chrom == 'chrY' | rna_PAR, 'RNA_chrY_or_PAR',
    !is.finite(rna_FDR) | !is.finite(rna_PValue) | !is.finite(rna_logFC), 'RNA_statistics_unavailable', default = 'eligible')]
  audits[[b]] <- copy(g)
  g <- g[mapping_status == 'eligible']
  g[, membership := fcase(minimum_peak_FDR < .1 & rna_FDR < .05, 'Both', minimum_peak_FDR < .1, 'ATAC only', rna_FDR < .05, 'RNA only', default = 'Neither')]
  g[, direction := fcase(membership != 'Both', 'Not jointly significant',
    atac_median_logFC > 0 & rna_logFC > 0, 'Up/up', atac_median_logFC < 0 & rna_logFC < 0, 'Down/down',
    atac_median_logFC * rna_logFC < 0, 'Opposing', default = 'Zero effect')]
  all[[b]] <- g
}
audit <- rbindlist(audits); g <- rbindlist(all)
fwrite(rbindlist(rna_audits), file.path(out, 'RNA_mapping_audit.csv'))
fwrite(audit, file.path(out, 'gene_mapping_audit.csv'))
fwrite(audit[, .N, by = .(batch, feature, mapping_status)], file.path(out, 'mapping_summary.csv'))
common <- g[, .(n_batches = uniqueN(batch)), by = .(feature, entrez_id)][n_batches == 2, .(feature, entrez_id)]
g[, common_to_both := paste(feature, entrez_id) %in% paste(common$feature, common$entrez_id)]
fwrite(g, file.path(out, 'all_gene_results.csv'))
summarize <- function(d, universe_type) {
  d[, .(universe_type = universe_type, eligible_genes = .N, ATAC_genes = sum(minimum_peak_FDR < .1),
    RNA_genes = sum(rna_FDR < .05), ATAC_only = sum(membership == 'ATAC only'),
    RNA_only = sum(membership == 'RNA only'), Both = sum(membership == 'Both'), Neither = sum(membership == 'Neither'),
    Up_up = sum(direction == 'Up/up'), Down_down = sum(direction == 'Down/down'), Opposing = sum(direction == 'Opposing'),
    Zero_effect = sum(direction == 'Zero effect'), spearman_all_genes = cor(atac_median_logFC, rna_logFC, method = 'spearman')),
    by = .(batch, feature)]
}
stats <- rbind(summarize(g, 'batch-specific'), summarize(g[common_to_both == TRUE], 'common-to-both'))
stopifnot(all(stats$ATAC_only + stats$RNA_only + stats$Both + stats$Neither == stats$eligible_genes),
  all(stats$Up_up + stats$Down_down + stats$Opposing + stats$Zero_effect == stats$Both),
  !anyDuplicated(g[, .(batch, feature, entrez_id)]))
fwrite(stats, file.path(out, 'overlap_summary.csv'))
for (b in batches) for (f in unique(g$feature)) {
  fwrite(g[batch == b & feature == f & membership == 'Both'], file.path(out, paste0(b, '_', f, '_overlap_genes.csv')))
  fwrite(g[batch == b & feature == f & direction == 'Down/down'], file.path(out, paste0(b, '_', f, '_down_down_genes.csv')))
}
# A direct two-batch comparison uses only genes eligible in both separate RNA fits.
wide <- dcast(g[common_to_both == TRUE], feature + entrez_id ~ batch,
  value.var = c('rna_symbol', 'rna_logFC', 'rna_FDR', 'membership', 'direction'))
fwrite(wide, file.path(out, 'common_gene_batch_comparison.csv'))
replication <- wide[, .(common_eligible = .N,
  overlap_in_both_batches = sum(membership_AUG_2025 == 'Both' & membership_MAY_2026 == 'Both'),
  down_down_in_both_batches = sum(direction_AUG_2025 == 'Down/down' & direction_MAY_2026 == 'Down/down'),
  up_up_in_both_batches = sum(direction_AUG_2025 == 'Up/up' & direction_MAY_2026 == 'Up/up')),
  by = feature]
fwrite(replication, file.path(out, 'replication_summary.csv'))
q <- fread('integration/down70_ATAC010_RNA005/exact_70_genes.csv', colClasses = list(character = 'entrez_id'))
stopifnot(nrow(q) == 70, !anyDuplicated(q$entrez_id))
qtrack <- merge(q[, .(entrez_id, combined_symbol = rna_symbol, ZNF_ZFP_named)],
  audit[feature == 'promoter', .(batch, entrez_id, mapping_status, rna_logFC, rna_FDR)],
  by = 'entrez_id', allow.cartesian = TRUE)
qtrack[, retained_down_down := mapping_status == 'eligible' & !is.na(rna_FDR) & rna_FDR < .05 & rna_logFC < 0]
stopifnot(nrow(qtrack) == 140)
fwrite(qtrack, file.path(out, 'combined_down70_batch_support.csv'))
qstats <- qtrack[, .(original_genes = .N, eligible = sum(mapping_status == 'eligible'),
  RNA_down_significant = sum(retained_down_down), RNA_down_any_FDR = sum(rna_logFC < 0, na.rm = TRUE),
  ZNF_ZFP_down_significant = sum(retained_down_down & ZNF_ZFP_named)), by = batch]
fwrite(qstats, file.path(out, 'combined_down70_support_summary.csv'))
palette <- c('Neither' = '#C7CDD1', 'ATAC only' = '#0072B2', 'RNA only' = '#D55E00', 'Both' = '#7B3294')
fl <- c(promoter = 'Promoter: TSS +/-2 kb', gene_body = 'Gene-body overlap', intergenic = 'Intergenic: nearest TSS (exploratory)')
bl <- c(AUG_2025 = 'August 2025 RNA-seq', MAY_2026 = 'May 2026 RNA-seq')
g[, `:=`(batch_label = factor(batch, levels = batches, labels = bl),
  feature_label = factor(feature, levels = names(fl), labels = fl))]
ps <- copy(stats[universe_type == 'batch-specific'])
ps[, `:=`(batch_label = factor(batch, levels = batches, labels = bl),
  feature_label = factor(feature, levels = names(fl), labels = fl))]
saveplot <- function(p, stem, w, h) {
  ggsave(file.path(out, paste0(stem, '.png')), p, width = w, height = h, dpi = 300, bg = 'white')
  ggsave(file.path(out, paste0(stem, '.pdf')), p, width = w, height = h, bg = 'white', useDingbats = FALSE)
}
plot_scatter <- function(d, s, label = FALSE) {
  s <- copy(s); s[, annotation := sprintf('n = %s | Both = %s\nrho = %.3f', format(eligible_genes, big.mark = ',', trim = TRUE), Both, spearman_all_genes)]
  p <- ggplot(d, aes(atac_median_logFC, rna_logFC)) +
    geom_hline(yintercept = 0, colour = 'grey60', linetype = 2, linewidth = .3) +
    geom_vline(xintercept = 0, colour = 'grey60', linetype = 2, linewidth = .3) +
    geom_point(data = d[membership == 'Neither'], colour = palette['Neither'], alpha = .25, size = .55) +
    geom_point(data = d[membership %in% c('ATAC only', 'RNA only')], aes(colour = membership), alpha = .6, size = .8) +
    geom_point(data = d[membership == 'Both'], aes(colour = membership), alpha = .9, size = 1.45) +
    geom_text(data = s, aes(x = -Inf, y = Inf, label = annotation), inherit.aes = FALSE, hjust = -.07, vjust = 1.15, size = 3.1) +
    scale_colour_manual(values = palette, breaks = names(palette)[-1], drop = FALSE) +
    facet_grid(feature_label ~ batch_label) + theme_bw(base_size = 11) +
    theme(legend.position = 'bottom', legend.title = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0, size = 9)) +
    labs(title = 'ATAC accessibility versus separate RNA-seq batches', subtitle = 'ATAC FDR <0.10 | RNA FDR <0.05 | LOY / CTRL',
      x = 'Median ATAC log2 fold change across linked peaks', y = 'RNA log2 fold change',
      caption = 'One point per eligible gene; same axes in both batches. Grey: neither significant.\nAny linked ATAC peak determines membership; median ATAC effects have no gene-level FDR.')
  if (label) {
    lab <- copy(d[direction %in% c('Up/up', 'Down/down')]); lab[, rank := pmin(abs(atac_median_logFC), abs(rna_logFC))]
    setorder(lab, -rank, entrez_id)
    lab <- lab[, head(.SD, 5), by = .(batch, direction)]
    fwrite(lab, file.path(out, 'promoter_selected_labels.csv'))
    p <- p + geom_text_repel(data = lab, aes(label = rna_symbol), seed = 42, size = 3,
      colour = '#20262B', min.segment.length = 0, max.overlaps = Inf, box.padding = .5)
  }
  p
}
saveplot(plot_scatter(g, ps), 'feature_scatter_grid', 12, 12)
for (f in names(fl)) saveplot(plot_scatter(g[feature == f], ps[feature == f], f == 'promoter'), paste0(f, '_scatters'), 12, 6.2)
angle <- seq(0, 2*pi, length.out = 200)
circles <- rbindlist(lapply(seq_len(nrow(ps)), function(i) rbindlist(lapply(c('ATAC', 'RNA'), function(a)
  data.table(batch_label = ps$batch_label[i], feature_label = ps$feature_label[i], assay = a,
    x = (if (a == 'ATAC') 3.7 else 6.3) + 2.25*cos(angle), y = 3.7 + 2.25*sin(angle))))))
counts <- rbindlist(lapply(seq_len(nrow(ps)), function(i) data.table(batch_label = ps$batch_label[i], feature_label = ps$feature_label[i],
  x = c(2.9, 5, 7.1), y = 3.7, label = c(ps$ATAC_only[i], ps$Both[i], ps$RNA_only[i]))))
p <- ggplot() + geom_polygon(data = circles, aes(x, y, fill = assay, group = assay), alpha = .27, colour = 'grey45', linewidth = .4) +
  geom_text(data = counts, aes(x, y, label = label), fontface = 'bold', size = 4.5) +
  annotate('text', x = 2.7, y = 6.5, label = 'ATAC genes', size = 3.5) + annotate('text', x = 7.3, y = 6.5, label = 'RNA genes', size = 3.5) +
  geom_text(data = ps, aes(x = 5, y = .7, label = paste0('Eligible: ', eligible_genes, ' | Neither: ', Neither)), size = 3) +
  facet_grid(feature_label ~ batch_label) + scale_fill_manual(values = c(ATAC = '#0072B2', RNA = '#D55E00')) +
  coord_equal(xlim = c(0, 10), ylim = c(0, 7.3)) + theme_void(base_size = 11) +
  theme(legend.position = 'none', plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0),
    strip.text.y = element_text(angle = -90), strip.background = element_rect(fill = '#edf3f5', colour = NA),
    plot.margin = margin(15, 20, 15, 20)) +
  labs(title = 'ATAC overlap with each RNA-seq batch', subtitle = 'ATAC FDR <0.10 | RNA FDR <0.05 | No fold-change or direction filter',
    caption = 'Distinct genes within each feature-specific eligible background. Circle areas are schematic.\nIntergenic assignments are positional hypotheses; genes may occur in multiple feature rows.')
saveplot(p, 'feature_overlap_grid', 11, 12)
methods <- c('Existing batch-stratified DESeq2 results from the HT1080_BATCHES notebook are reused without refitting or modifying FDR.',
  'Each separate RNA fit uses ~ treatment: three LOY versus three CTRL; parental and add-back samples are excluded. Positive logFC means higher in LOY.',
  'These are batch-specific fits within the combined notebook, not necessarily historical standalone analysis files. Input genes inherit the combined count filter (at least 10 counts in at least 3 samples across the combined input), then each batch is fit separately with its own DESeq2 normalization, dispersion and independent filtering.',
  'RNA annotation and Entrez IDs are carried from the source notebook (EnsDb.Hsapiens.v86; counts documented as GENCODE v49). Source logCPM is DESeq2 baseMean.',
  'All original ATAC links are reused, including those ineligible in the combined RNA comparison. Primary ATAC differential effects are unchanged.',
  'Eligible universe is determined independently for each batch and feature: unique Entrez RNA mapping, matching chromosome, no chrY/PAR RNA overlap, finite logFC/PValue/FDR.',
  'ATAC FDR <0.10 and RNA FDR <0.05; no fold-change filter. Any qualifying linked peak flags an ATAC gene. Minimum peak FDR is not a gene-level FDR.',
  'Scatter x coordinates and concordant direction use median effects over all linked tested peaks. Significant up/down peak counts are separately exported to detect mixed support.',
  'Promoters: TSS +/-2 kb; gene bodies: overlapping genes; intergenic: nearest TSS, exploratory. Genes may occur in multiple features.',
  'Each scatter pair has shared axes and includes every eligible point. Promoter labels show up to five jointly significant genes per concordant direction per batch, ranked by minimum absolute effect across assays.',
  'Main panels use per-batch eligible universes. Additional common-to-both summaries remove differences in gene eligibility when comparing batches.',
  'The original combined-model 70 down/down genes are tracked descriptively; these batches contributed to their selection, so this is not independent validation.',
  'No overlap significance P-value is assigned because shared peaks, chromosome clustering and variable linked-peak counts violate simple equal-probability gene sampling.',
  'This is group-effect integration, not confirmed sample-matched multiomics. No passage adjustment in these RNA fits and no lpWGS/CN adjustment. Differences in significance alone do not demonstrate a batch-by-treatment interaction.')
writeLines(methods, file.path(out, 'METHODS.txt'))
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
writeLines(c('PASS: source SHA256 matches frozen copies; original model/contrast and sample counts checked',
  'PASS: all original ATAC links considered; unique batch/feature/gene rows and explicit eligibility audit',
  'PASS: every set partition and direction partition sums to its universe or overlap count',
  'PASS: exact original 70 genes tracked once per batch; common-universe comparisons saved'), file.path(out, 'verification.txt'))
print(stats); print(replication); print(qstats)
