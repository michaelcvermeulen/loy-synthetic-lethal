# Expanded exploratory CAMERA screen; original three-collection results are retained.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(limma); library(msigdbr); library(ggplot2)})
roots <- c(getwd(), file.path(getwd(), 'THESIS_ADDBACK_2026'))
root <- roots[file.exists(file.path(roots, 'outputs/fit_primary.rds'))][1]
setwd(root)
out <- 'outputs/camera_full'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fb <- readRDS('outputs/fit_primary.rds')
map <- fread('outputs/camera/camera_gene_universe.csv')
v <- fb$v[map$original_row, ]
stopifnot(nrow(v) == 17724, !anyDuplicated(map$entrez_id),
          identical(as.character(v$genes$gene_id), map$gene_id))
metadata <- as.data.table(msigdbr_collections())
fwrite(metadata, file.path(out, 'available_collections.csv'))
message('Loading all human MSigDB gene sets...')
resource <- as.data.table(msigdbr(species = 'Homo sapiens'))
stopifnot(all(resource$db_version == '2026.1.Hs'))
eid <- intersect(c('ncbi_gene', 'entrez_gene'), names(resource))[1]
resource[, entrez_id := as.character(get(eid))]
sets <- unique(resource[, .(gene_set = gs_name, collection = gs_collection,
  subcollection = gs_subcollection, collection_name = gs_collection_name, db_version)])
stopifnot(!anyDuplicated(sets$gene_set))
membership <- unique(resource[, .(gene_set = gs_name, entrez_id)])
rm(resource); invisible(gc())
membership[, camera_row := match(entrez_id, as.character(map$entrez_id))]
size_all <- membership[, .(resource_genes = .N), by = gene_set]
mapped <- membership[!is.na(camera_row)]
size_mapped <- mapped[, .(mapped_genes = uniqueN(camera_row)), by = gene_set]
sets <- merge(sets, size_all, by = 'gene_set', all.x = TRUE)
sets <- merge(sets, size_mapped, by = 'gene_set', all.x = TRUE)
sets[is.na(mapped_genes), mapped_genes := 0L]
sets[, status := fifelse(mapped_genes < 10, 'Excluded: fewer than 10 mapped genes',
                 fifelse(mapped_genes > 500, 'Excluded: more than 500 mapped genes', 'Tested'))]
fwrite(sets, file.path(out, 'gene_set_inventory.csv'))
tested <- sets[status == 'Tested']
mapped <- mapped[gene_set %in% tested$gene_set]
index <- split(mapped$camera_row, mapped$gene_set)
index <- lapply(index, unique)
stopifnot(all(lengths(index) >= 10 & lengths(index) <= 500), setequal(names(index), tested$gene_set))
mapped[, gene_id := map$gene_id[camera_row]]
mapped[, gene_symbol := map$gene_symbol[camera_row]]
fwrite(mapped, file.path(out, 'gene_set_membership_tested.csv.gz'))
fwrite(map, file.path(out, 'camera_gene_universe.csv'))
summary_resource <- sets[, .(available = .N, tested = sum(status == 'Tested'),
  below_10 = sum(mapped_genes < 10), above_500 = sum(mapped_genes > 500)),
  by = .(collection, subcollection, collection_name, db_version)]
fwrite(summary_resource, file.path(out, 'resource_summary.csv'))
contrast_labels <- c(UTY_v_LXHA = 'UTY vs empty vector', KDM5D_v_LXHA = 'KDM5D vs empty vector',
  DUAL_v_LXHA = 'UTY + KDM5D vs empty vector', DUAL_v_UTY = 'UTY + KDM5D vs UTY',
  DUAL_v_KDM5D = 'UTY + KDM5D vs KDM5D', LOY_v_LXHA = 'Untransduced LOY vs empty vector')
results <- list()
for (cn in names(contrast_labels)) {
  message('CAMERA full: ', cn, ' / ', length(index), ' sets')
  z <- as.data.table(camera(v, index = index, design = v$design, contrast = fb$cm[, cn],
       inter.gene.cor = 0.01, allow.neg.cor = FALSE, sort = FALSE), keep.rownames = 'gene_set')
  z[, c('FDR') := NULL]
  z[, `:=`(contrast = cn, comparison = unname(contrast_labels[cn]))]
  setnames(z, c('NGenes', 'Direction', 'PValue'), c('n_genes', 'direction', 'p_value'))
  z <- merge(z, tested[, .(gene_set, collection, subcollection)], by = 'gene_set', sort = FALSE)
  # All FDR families are explicit; primary is all sets within the contrast.
  z[, FDR_all_sets_contrast := p.adjust(p_value, 'BH')]
  z[, FDR_collection_contrast := p.adjust(p_value, 'BH'), by = collection]
  z[, FDR_subcollection_contrast := p.adjust(p_value, 'BH'), by = .(collection, subcollection)]
  stats <- lapply(z$gene_set, function(gs) {
    ii <- map$original_row[index[[gs]]]
    lfc <- fb$eb$coefficients[ii, cn]
    tt <- fb$eb$t[ii, cn]
    list(mean_logFC = mean(lfc), top_genes = paste(head(map$gene_symbol[index[[gs]]][order(-abs(tt))], 8), collapse = ';'))
  })
  z[, mean_logFC := vapply(stats, `[[`, numeric(1), 'mean_logFC')]
  z[, top_genes_by_abs_t := vapply(stats, `[[`, character(1), 'top_genes')]
  results[[cn]] <- z
  fwrite(z, file.path(out, paste0('camera_', cn, '.csv')))
}
ans <- rbindlist(results)
ans[, FDR_global := p.adjust(p_value, 'BH')]
ans[, inter_gene_correlation := 0.01]
setorder(ans, contrast, FDR_all_sets_contrast, p_value)
fwrite(ans, file.path(out, 'camera_all.csv'))
fwrite(ans[FDR_all_sets_contrast < 0.05], file.path(out, 'camera_significant.csv'))
overview <- ans[, .(sets_tested = .N, significant = sum(FDR_all_sets_contrast < 0.05),
  up = sum(FDR_all_sets_contrast < 0.05 & direction == 'Up'),
  down = sum(FDR_all_sets_contrast < 0.05 & direction == 'Down')),
  by = .(collection, contrast, comparison)]
fwrite(overview, file.path(out, 'summary.csv'))

# Overlapping tests must reproduce their original raw P values and direction.
old <- fread('outputs/camera/camera_all.csv')
check <- merge(old[, .(gene_set, contrast, original_p = p_value, original_direction = direction)],
               ans[, .(gene_set, contrast, p_value, direction)], by = c('gene_set', 'contrast'))
stopifnot(nrow(check) == nrow(old), all(check$original_direction == check$direction),
          isTRUE(all.equal(check$original_p, check$p_value, tolerance = 1e-10)),
          nrow(ans) == length(index) * length(contrast_labels),
          !anyNA(ans$p_value), all(ans$p_value >= 0 & ans$p_value <= 1))
fwrite(data.table(test = 'Overlap with original CAMERA screen', n = nrow(check),
  max_absolute_p_difference = max(abs(check$original_p-check$p_value))), file.path(out, 'overlap_verification.csv'))

STYLE <- 'A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/skills/thesis-figure-style'
source(file.path(STYLE, 'scripts/thesis_style.R')); register_thesis_fonts(STYLE)
overview[, comparison := factor(comparison, unname(contrast_labels))]
p <- ggplot(overview, aes(comparison, collection, fill = significant)) +
  geom_tile(color = 'white') + geom_text(aes(label = significant,
    color = significant > max(significant)/2), family = 'Open Sans', size = 3) +
  scale_color_manual(values = c('TRUE' = 'black', 'FALSE' = 'white'), guide = 'none') +
  scale_fill_viridis_c(option = 'magma', name = 'Significant\nsets') +
  labs(x = NULL, y = 'MSigDB collection', title = 'CAMERA across all human MSigDB collections',
       subtitle = 'Counts at FDR < 0.05, adjusted across all tested sets within each contrast') +
  theme_thesis() + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
ggsave('figures/student_camera_full_overview.png', p, width = 10, height = 6, dpi = 300, device = ragg::agg_png)
ggsave('figures/student_camera_full_overview.svg', p, width = 10, height = 6, device = svglite::svglite)

writeLines(c(
  'Exploratory CAMERA screen using all human MSigDB collections returned by msigdbr, release 2026.1.Hs (H and C1-C9).',
  'The inventory records every set including those excluded. Retained sets contain 10-500 genes mapped to the frozen 17,724-Entrez-gene universe.',
  'Same primary voom expression, weights, matched-collection-occasion design and six contrasts as the original CAMERA analysis. No restriction to the 91 candidates or August signature.',
  'limma::camera: inter.gene.cor=0.01, allow.neg.cor=FALSE, sort=FALSE. No permutation-based P values.',
  'Primary FDR_all_sets_contrast: BH across all tested sets within each contrast. Also report BH within collection/contrast, within subcollection/contrast and globally across all six contrasts.',
  'Per-contrast files contain the three within-contrast FDR variants; combined camera_all.csv additionally contains the global FDR.',
  'Mean logFC and top genes by absolute moderated t are descriptive annotations; they do not define the CAMERA statistic or identify causal drivers.',
  'Overlapping pathway memberships and related signatures are not independent biological discoveries. C1 tests positional sets; many C3-C9 sets are signatures rather than pathways.',
  'Each construct is one established population sampled on three occasions; results are conditional population associations, not independently replicated construct effects.',
  'https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp',
  'https://igordot.github.io/msigdbr/reference/msigdbr.html'), file.path(out, 'METHODS.txt'))
inputs <- c('outputs/fit_primary.rds', 'outputs/camera/camera_gene_universe.csv', 'code/run_camera_full.R',
            file.path(out, c('available_collections.csv', 'gene_set_inventory.csv', 'gene_set_membership_tested.csv.gz')))
fwrite(data.frame(file = inputs, sha256 = vapply(inputs, digest::digest, character(1), file = TRUE, algo = 'sha256')),
       file.path(out, 'input_manifest.csv'))
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
writeLines(c(sprintf('PASS %d resource sets; %d tested sets; %d collections; %d contrasts; %d tests.',
  nrow(sets), length(index), uniqueN(sets$collection), length(contrast_labels), nrow(ans)),
  sprintf('PASS %d overlapping tests reproduce original P values and directions.', nrow(check))), file.path(out, 'VERIFICATION.txt'))
cat(readLines(file.path(out, 'VERIFICATION.txt')), sep = '\n')
