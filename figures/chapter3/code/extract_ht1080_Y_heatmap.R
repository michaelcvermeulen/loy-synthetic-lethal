# Display extraction only: reuse saved all-15 TMM effective library sizes.
root <- 'A:/Dropbox/LOY_paralogs/projects'
fig <- file.path(root, 'SEPT_11_2026_THESIS/figures/chapter3')
model <- file.path(root, 'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/models/01_ingest.rds')
x <- readRDS(model)
n <- read.csv(file.path(fig, 'source_data/ht1080_recurrent_gene_expression/all15_context_normalization.csv'))
e <- read.csv(file.path(fig, 'source_data/shared_sex_chromosome/HT1080_effects.csv'))
m <- x$meta
effective <- n$effective_library_reads[match(m$sample_key, n$sample_key)]
stopifnot(nrow(m) == 15, !anyNA(effective), all(effective > 0))
e <- e[e$chrom == 'chrY' & !e$PAR & !grepl('^ENSG', e$gene_name), ]
stopifnot(!anyDuplicated(e$gene_id), all(e$gene_id %in% rownames(x$counts)))
counts <- x$counts[e$gene_id, m$sample_key, drop=FALSE]
cp <- sweep(counts, 2, effective, '/') * 1e6
aug <- rowMeans(cp[, m$treatment == 'CTRL' & m$batch == 'AUG_2025', drop=FALSE])
may <- rowMeans(cp[, m$treatment == 'CTRL' & m$batch == 'MAY_2026', drop=FALSE])
keep <- pmax(aug, may) >= 1
e <- e[keep, ]; cp <- cp[keep, , drop=FALSE]; counts <- counts[keep, , drop=FALSE]
aug <- aug[keep]; may <- may[keep]
short <- paste0(ifelse(m$batch == 'AUG_2025', 'A', 'M'),
                c(CTRL='C', LOY='L', PARENTAL='P')[m$treatment], m$passage_number)
z <- do.call(rbind, lapply(seq_len(nrow(e)), function(i) data.frame(
    gene_id=e$gene_id[i], gene_name=e$gene_name[i], start=e$start[i],
    sample_key=m$sample_key, short=short, treatment=m$treatment,
    batch=m$batch, passage=m$passage_number, effective_library_reads=effective,
    count=as.numeric(counts[i, ]), CPM=as.numeric(cp[i, ]),
    log2_CPM_plus1=log2(as.numeric(cp[i, ]) + 1),
    AUG_sgCTRL_mean_CPM=unname(aug[i]), MAY_sgCTRL_mean_CPM=unname(may[i]))))
stopifnot(nrow(z) == nrow(e) * 15, all(is.finite(z$log2_CPM_plus1)))
out <- file.path(fig, 'source_data/ht1080_rna_expression_changes')
write.csv(z, file.path(out, 'panel_b_Y_sample_expression.csv'), row.names=FALSE)
cat(sprintf('Saved %d Y-linked genes across 15 libraries (%d values); no model refitted.\n', nrow(e), nrow(z)))

# Add three leading qualifying hits in each direction, ranked by saved FDR.
de <- read.csv(file.path(root, 'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1/tables/gene_consistency.csv'))
eligible <- de[de$FDR < .05 & abs(de$logFC) >= 1 &
               de$logFC * de$AUG_logFC > 0 & de$logFC * de$MAY_logFC > 0, ]
hits <- do.call(rbind, lapply(c('Higher', 'Lower'), function(direction) {
    q <- eligible[if (direction == 'Higher') eligible$logFC > 0 else eligible$logFC < 0, ]
    q <- q[order(q$FDR, q$gene_id), ][seq_len(3), ]
    q$heatmap_group <- direction
    q$heatmap_rank <- seq_len(3)
    q
}))
stopifnot(nrow(hits) == 6, !anyDuplicated(hits$gene_id))
hcounts <- x$counts[hits$gene_id, m$sample_key, drop=FALSE]
hcp <- sweep(hcounts, 2, effective, '/') * 1e6
hz <- do.call(rbind, lapply(seq_len(nrow(hits)), function(i) data.frame(
    gene_id=hits$gene_id[i], gene_name=hits$gene_name[i], start=hits$start[i],
    heatmap_group=hits$heatmap_group[i], heatmap_rank=hits$heatmap_rank[i],
    logFC=hits$logFC[i], FDR=hits$FDR[i], sample_key=m$sample_key, short=short,
    treatment=m$treatment, batch=m$batch, passage=m$passage_number,
    effective_library_reads=effective, count=as.numeric(hcounts[i, ]),
    CPM=as.numeric(hcp[i, ]), log2_CPM_plus1=log2(as.numeric(hcp[i, ]) + 1))))
write.csv(hz, file.path(out, 'panel_b_hit_sample_expression.csv'), row.names=FALSE)
write.csv(hits[, c('gene_id', 'gene_name', 'heatmap_group', 'heatmap_rank', 'logFC', 'FDR', 'AUG_logFC', 'MAY_logFC')],
          file.path(out, 'panel_b_hit_selection.csv'), row.names=FALSE)
cat('Added six top hits:', paste(hits$gene_name, collapse=', '), '\n')
