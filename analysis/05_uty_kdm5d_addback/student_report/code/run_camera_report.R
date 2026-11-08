options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(msigdbr)
})

find_root <- function() {
  candidates <- c(
    normalizePath(getwd(), winslash = '/', mustWork = FALSE),
    normalizePath(file.path(getwd(), 'THESIS_ADDBACK_2026'), winslash = '/', mustWork = FALSE),
    normalizePath(file.path(dirname(getwd()), 'THESIS_ADDBACK_2026'), winslash = '/', mustWork = FALSE)
  )
  hit <- candidates[file.exists(file.path(candidates, 'outputs', 'fit_primary.rds'))]
  if (!length(hit)) stop('Could not locate THESIS_ADDBACK_2026/outputs/fit_primary.rds')
  hit[[1]]
}

root <- find_root()
setwd(root)
dir.create('outputs/camera', recursive = TRUE, showWarnings = FALSE)
dir.create('figures', recursive = TRUE, showWarnings = FALSE)
dir.create('student_report/data', recursive = TRUE, showWarnings = FALSE)
dir.create('student_report/figures', recursive = TRUE, showWarnings = FALSE)
dir.create('student_report/code', recursive = TRUE, showWarnings = FALSE)
dir.create('report/figures', recursive = TRUE, showWarnings = FALSE)
dir.create('report/figures/student_panels', recursive = TRUE, showWarnings = FALSE)

fit_bundle <- readRDS('outputs/fit_primary.rds')
stopifnot(all(c('v', 'cm', 'eb') %in% names(fit_bundle)))
v <- fit_bundle$v
design <- v$design
cm <- fit_bundle$cm
eb <- fit_bundle$eb
genes <- as.data.frame(v$genes)

if (!identical(rownames(v$E), rownames(genes))) stop('Expression and annotation rows do not match')
if (!all(c('gene_name', 'entrez_id', 'gene_id') %in% names(genes))) stop('Required gene annotations are absent')

# CAMERA operates on one representative row per Entrez gene so duplicated symbols
# or identifiers do not receive extra weight. The most highly expressed row is kept.
map <- data.frame(
  original_row = seq_len(nrow(genes)),
  gene_id = as.character(genes$gene_id),
  gene_symbol = as.character(genes$gene_name),
  entrez_id = as.character(genes$entrez_id),
  average_logCPM = rowMeans(v$E),
  stringsAsFactors = FALSE
)
map <- map[!is.na(map$entrez_id) & nzchar(map$entrez_id) & map$entrez_id != 'NA', ]
map <- map[order(map$entrez_id, -map$average_logCPM, map$gene_id), ]
map <- map[!duplicated(map$entrez_id), ]
map <- map[order(map$original_row), ]
map$camera_row <- seq_len(nrow(map))
v_camera <- v[map$original_row, ]

collections <- list(
  Hallmark = list(collection = 'H', subcollection = NULL),
  Reactome = list(collection = 'C2', subcollection = 'CP:REACTOME'),
  GO_BP = list(collection = 'C5', subcollection = 'GO:BP')
)

load_collection <- function(x) {
  args <- list(species = 'Homo sapiens', collection = x$collection)
  if (!is.null(x$subcollection)) args$subcollection <- x$subcollection
  do.call(msigdbr::msigdbr, args)
}

msig_tables <- lapply(collections, load_collection)
message('Loaded MSigDB collections through msigdbr.')
entrez_column <- function(tb) {
  z <- intersect(c('ncbi_gene', 'entrez_gene'), names(tb))
  if (!length(z)) stop('No Entrez identifier column found in msigdbr output')
  z[[1]]
}

membership_all <- rbindlist(lapply(names(msig_tables), function(collection_name) {
  tb <- as.data.frame(msig_tables[[collection_name]])
  eid <- entrez_column(tb)
  data.table(
    collection = collection_name,
    gene_set = as.character(tb$gs_name),
    gene_symbol = as.character(tb$gene_symbol),
    entrez_id = as.character(tb[[eid]]),
    db_version = as.character(tb$db_version)
  )
}), use.names = TRUE)
membership_all <- unique(membership_all[!is.na(entrez_id) & nzchar(entrez_id)])

tested_sets <- list()
tested_membership <- list()
resource_summary <- list()
for (collection_name in names(collections)) {
  mem <- membership_all[collection == collection_name]
  split_ids <- split(mem$entrez_id, mem$gene_set)
  indices <- lapply(split_ids, function(ids) {
    unique(map$camera_row[match(unique(ids), map$entrez_id, nomatch = 0L)])
  })
  indices <- lapply(indices, function(z) z[z > 0L])
  original_sizes <- lengths(indices)
  keep <- original_sizes >= 10L & original_sizes <= 500L
  indices <- indices[keep]
  tested_sets[[collection_name]] <- indices
  used <- rbindlist(lapply(names(indices), function(set_name) {
    ii <- indices[[set_name]]
    data.table(
      collection = collection_name,
      gene_set = set_name,
      camera_row = ii,
      gene_id = map$gene_id[ii],
      gene_symbol = map$gene_symbol[ii],
      entrez_id = map$entrez_id[ii]
    )
  }))
  tested_membership[[collection_name]] <- used
  resource_summary[[collection_name]] <- data.frame(
    collection = collection_name,
    msigdb_version = paste(sort(unique(mem$db_version)), collapse = ';'),
    sets_in_resource = length(split_ids),
    sets_tested = length(indices),
    min_size = 10L,
    max_size = 500L,
    camera_gene_universe = nrow(v_camera)
  )
}
tested_membership <- rbindlist(tested_membership)
resource_summary <- rbindlist(resource_summary)

contrast_labels <- c(
  UTY_v_LXHA = 'UTY vs empty vector',
  KDM5D_v_LXHA = 'KDM5D vs empty vector',
  DUAL_v_LXHA = 'UTY + KDM5D vs empty vector',
  DUAL_v_UTY = 'UTY + KDM5D vs UTY',
  DUAL_v_KDM5D = 'UTY + KDM5D vs KDM5D',
  LOY_v_LXHA = 'Untransduced LOY vs empty vector'
)
stopifnot(all(names(contrast_labels) %in% colnames(cm)))

camera_results <- list()
for (collection_name in names(tested_sets)) {
  index <- tested_sets[[collection_name]]
  for (contrast_name in names(contrast_labels)) {
    message('CAMERA: ', collection_name, ' / ', contrast_name)
    ans <- limma::camera(
      v_camera,
      index = index,
      design = design,
      contrast = cm[, contrast_name],
      inter.gene.cor = 0.01,
      allow.neg.cor = FALSE,
      sort = FALSE
    )
    ans$gene_set <- rownames(ans)
    ans$collection <- collection_name
    ans$contrast <- contrast_name
    ans$contrast_label <- unname(contrast_labels[[contrast_name]])
    camera_results[[paste(collection_name, contrast_name, sep = '__')]] <- as.data.frame(ans)
  }
}
camera_all <- rbindlist(camera_results, use.names = TRUE, fill = TRUE)
message('Annotating CAMERA results with gene-level summaries.')
setnames(camera_all, c('NGenes', 'Direction', 'PValue'), c('n_genes', 'direction', 'p_value'))

# Add transparent gene-level summaries. These are descriptive companions to
# CAMERA and are not used to calculate the competitive-test P value.
membership_by_set <- split(tested_membership$camera_row,
                           paste(tested_membership$collection, tested_membership$gene_set, sep = '__'))
q_by_contrast <- lapply(names(contrast_labels), function(contrast_name) {
  p.adjust(eb$p.value[, contrast_name], method = 'BH')
})
names(q_by_contrast) <- names(contrast_labels)
annotate_row <- function(collection_name, set_name, contrast_name, direction) {
  ii <- unique(membership_by_set[[paste(collection_name, set_name, sep = '__')]])
  lfc <- eb$coefficients[map$original_row[ii], contrast_name]
  tt <- eb$t[map$original_row[ii], contrast_name]
  qv <- q_by_contrast[[contrast_name]][map$original_row[ii]]
  selected <- qv < 0.05 & abs(lfc) >= 0.5
  o <- if (direction == 'Up') order(-tt, -abs(lfc), na.last = NA) else order(tt, -abs(lfc), na.last = NA)
  top <- unique(map$gene_symbol[ii][o])
  top <- top[!is.na(top) & nzchar(top)]
  list(
    mean_logFC = mean(lfc, na.rm = TRUE),
    median_logFC = median(lfc, na.rm = TRUE),
    n_selected_genes = sum(selected, na.rm = TRUE),
    top_directional_genes = paste(head(top, 8), collapse = ';')
  )
}

extras <- lapply(seq_len(nrow(camera_all)), function(i) {
  annotate_row(camera_all$collection[i], camera_all$gene_set[i],
               camera_all$contrast[i], camera_all$direction[i])
})
camera_all[, mean_logFC := vapply(extras, `[[`, numeric(1), 'mean_logFC')]
camera_all[, median_logFC := vapply(extras, `[[`, numeric(1), 'median_logFC')]
camera_all[, n_selected_genes := vapply(extras, `[[`, integer(1), 'n_selected_genes')]
camera_all[, top_directional_genes := vapply(extras, `[[`, character(1), 'top_directional_genes')]
camera_all[, FDR_collection_contrast := p.adjust(p_value, method = 'BH'), by = .(collection, contrast)]
camera_all[, FDR_contrast_all_collections := p.adjust(p_value, method = 'BH'), by = contrast]
camera_all[, FDR_global := p.adjust(p_value, method = 'BH')]
camera_all[, inter_gene_correlation := 0.01]
setcolorder(camera_all, c(
  'collection', 'contrast', 'contrast_label', 'gene_set', 'n_genes', 'direction',
  'p_value', 'FDR_collection_contrast', 'FDR_contrast_all_collections', 'FDR_global',
  'mean_logFC', 'median_logFC', 'n_selected_genes', 'top_directional_genes',
  'inter_gene_correlation'
))
setorder(camera_all, contrast, collection, FDR_collection_contrast, p_value)

camera_sig <- camera_all[FDR_collection_contrast < 0.05]
camera_summary <- camera_all[, .(
  sets_tested = .N,
  sets_FDR_lt_0.05 = sum(FDR_collection_contrast < 0.05),
  up_FDR_lt_0.05 = sum(FDR_collection_contrast < 0.05 & direction == 'Up'),
  down_FDR_lt_0.05 = sum(FDR_collection_contrast < 0.05 & direction == 'Down'),
  minimum_FDR = min(FDR_collection_contrast),
  top_gene_set = gene_set[which.min(FDR_collection_contrast)],
  top_direction = direction[which.min(FDR_collection_contrast)]
), by = .(contrast, contrast_label, collection)]
camera_summary <- camera_summary[
  order(match(contrast, names(contrast_labels)), match(collection, names(collections)))
]

fwrite(camera_all, 'outputs/camera/camera_all.csv')
fwrite(camera_sig, 'outputs/camera/camera_significant_FDR05.csv')
fwrite(camera_summary, 'outputs/camera/camera_summary.csv')
fwrite(resource_summary, 'outputs/camera/gene_set_resource_summary.csv')
fwrite(tested_membership, 'outputs/camera/gene_set_membership_tested.csv')
fwrite(map, 'outputs/camera/camera_gene_universe.csv')
message('Wrote CAMERA result and membership tables.')

method_text <- c(
  'Competitive gene-set testing used limma::camera on the frozen voom object and its observation-level weights.',
  'The design was the primary no-intercept six-group model with fixed effects for matched collection occasion.',
  'The fixed inter-gene correlation was 0.01 and negative correlations were not allowed.',
  'The gene universe was restricted to one expressed representative row per mapped Entrez gene; for duplicates, the row with the highest mean logCPM was retained.',
  'MSigDB Hallmark, Reactome and GO Biological Process sets were obtained through the installed msigdbr package and restricted to 10-500 mapped genes.',
  'The reported primary FDR is Benjamini-Hochberg correction within each contrast and collection. Across-collection and global FDR values are also exported.',
  'Direction is expression in the first named group relative to the second. CAMERA tests coordinated differential expression relative to genes outside the set; it does not test over-representation in a preselected DE list.',
  'Because each construct condition represents one established population sampled on three occasions, pathway results are conditional on these populations and do not estimate variability among independently generated transductions.'
)
writeLines(method_text, 'outputs/camera/METHODS.txt', useBytes = TRUE)

clean_set <- function(x) {
  x <- sub('^HALLMARK_', '', x)
  x <- sub('^REACTOME_', '', x)
  x <- sub('^GOBP_', '', x)
  gsub('_', ' ', x)
}

primary_contrasts <- names(contrast_labels)[1:3]
plot_dat <- camera_all[contrast %in% primary_contrasts]
plot_dat <- plot_dat[, head(.SD, 5), by = .(contrast, collection)]
plot_dat[, display_set := clean_set(gene_set)]
plot_dat[, signed_score := ifelse(direction == 'Up', 1, -1) * pmin(-log10(pmax(FDR_collection_contrast, 1e-300)), 20)]
plot_dat[, contrast_label := factor(contrast_label, levels = unname(contrast_labels[primary_contrasts]))]
plot_dat[, collection := factor(collection, levels = names(collections))]
plot_dat[, display_set := factor(display_set, levels = rev(unique(display_set)))]

p_camera <- ggplot(plot_dat, aes(signed_score, display_set, colour = direction, size = -log10(pmax(FDR_collection_contrast, 1e-300)))) +
  geom_vline(xintercept = 0, colour = 'grey75', linewidth = 0.35) +
  geom_point(alpha = 0.9) +
  facet_grid(collection ~ contrast_label, scales = 'free_y', space = 'free_y') +
  scale_colour_manual(values = c(Down = '#386CB0', Up = '#E34A33')) +
  scale_size_continuous(range = c(1.8, 5.5), guide = 'none') +
  labs(
    x = 'Signed -log10 within-collection FDR', y = NULL,
    colour = 'Direction',
    title = 'CAMERA: strongest competitive gene-set results',
    subtitle = 'Five lowest-FDR sets per collection and primary add-back contrast; positive values are higher in the add-back group'
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.text = element_text(face = 'bold', size = 8),
    panel.grid.major.y = element_blank(),
    legend.position = 'top',
    plot.title = element_text(face = 'bold')
  )
ggsave('figures/student_camera_overview.png', p_camera, width = 14, height = 12, dpi = 220, bg = 'white')
ggsave('figures/student_camera_overview.pdf', p_camera, width = 14, height = 12, device = cairo_pdf)
ggsave('figures/student_camera_overview.svg', p_camera, width = 14, height = 12)

# Copy the compact handoff tables, source figures and code before rendering.
capture.output(sessionInfo(), file = 'outputs/camera/sessionInfo_CAMERA.txt')
handoff_files <- c(
  'outputs/sample_manifest.csv',
  'outputs/de_summary.csv',
  'outputs/construct_rna_abundance.csv',
  'outputs/pca_scores.csv',
  'outputs/pca_variance.csv',
  'outputs/pca_addback_only_scores.csv',
  'outputs/pca_addback_only_variance.csv',
  'outputs/claim_evidence.csv',
  'outputs/signature/frozen_signature_mapping.csv',
  'outputs/signature/frozen_signature_scores.csv',
  'outputs/signature/score_contrasts.csv',
  'outputs/signature/fry_directional.csv',
  'outputs/signature/distance_to_control_summary.csv',
  'outputs/sensitivity/unblocked_summary.csv',
  'outputs/sensitivity/leave_one_occasion_out.csv',
  'outputs/camera/camera_all.csv',
  'outputs/camera/camera_significant_FDR05.csv',
  'outputs/camera/camera_summary.csv',
  'outputs/camera/gene_set_resource_summary.csv',
  'outputs/camera/gene_set_membership_tested.csv',
  'outputs/camera/camera_gene_universe.csv',
  'outputs/camera/METHODS.txt',
  'outputs/camera/sessionInfo_CAMERA.txt',
  'provenance/input_manifest.csv',
  'provenance/verification.txt',
  'provenance/sessionInfo.txt'
)
de_files <- file.path('outputs/de', paste0('de_', names(contrast_labels), '.csv'))
handoff_files <- c(handoff_files, de_files)
for (src in handoff_files) {
  if (!file.exists(src)) stop('Missing handoff source: ', src)
  dest <- file.path('student_report/data', gsub('[/\\\\]', '__', src))
  file.copy(src, dest, overwrite = TRUE)
}
embed_figure_files <- c(
  'figures/student_design_expression_readable.png',
  'figures/ch03_main_ht1080_addback_baseline_reversal.png',
  'figures/ch03_supp_ht1080_addback_sensitivity.png',
  'figures/student_camera_overview.png',
  'figures/student_panels/panel_e_pca_addback_only.png'
)
handoff_figure_files <- c(
  list.files('figures', pattern = '^(ch03_main|ch03_supp|student_).+\\.(png|pdf|svg)$',
             full.names = TRUE, recursive = FALSE),
  list.files('figures/student_panels', pattern = '\\.(png|pdf|svg)$',
             full.names = TRUE, recursive = FALSE)
)
file.copy(handoff_figure_files,
          file.path('student_report/figures', basename(handoff_figure_files)), overwrite = TRUE)
file.copy('figures/student_figure_index.csv',
          'student_report/figures/student_figure_index.csv', overwrite = TRUE)
# Pandoc resolves figure paths relative to the Rmd input directory. These
# temporary report-side copies allow --embed-resources to make the HTML truly
# self-contained; the portable handoff also retains explicit PNG copies.
file.copy(embed_figure_files[!grepl('student_panels', embed_figure_files)],
          file.path('report/figures', basename(embed_figure_files[!grepl('student_panels', embed_figure_files)])),
          overwrite = TRUE)
file.copy(embed_figure_files[grepl('student_panels', embed_figure_files)],
          file.path('report/figures/student_panels', basename(embed_figure_files[grepl('student_panels', embed_figure_files)])),
          overwrite = TRUE)
file.copy('code/run_analysis.R', 'student_report/code/run_analysis.R', overwrite = TRUE)
file.copy('code/run_camera_report.R', 'student_report/code/run_camera_report.R', overwrite = TRUE)
file.copy('report/student_report.Rmd', 'student_report/code/student_report.Rmd', overwrite = TRUE)

hash_targets <- c('outputs/dge_primary.rds', 'outputs/fit_primary.rds',
                  'code/run_analysis.R', 'code/run_camera_report.R',
                  'report/student_report.Rmd',
                  'outputs/camera/gene_set_resource_summary.csv',
                  'outputs/camera/gene_set_membership_tested.csv')
hashes <- vapply(hash_targets, digest::digest, character(1), algo = 'sha256', file = TRUE)
manifest <- data.frame(
  file = names(hashes),
  sha256 = unname(hashes),
  bytes = file.info(names(hashes))$size,
  stringsAsFactors = FALSE
)
fwrite(manifest, 'outputs/camera/camera_input_manifest.csv')
file.copy('outputs/camera/camera_input_manifest.csv',
          'student_report/data/outputs__camera__camera_input_manifest.csv', overwrite = TRUE)

if (!requireNamespace('rmarkdown', quietly = TRUE)) stop('rmarkdown is required to render the report')
source('code/build_student_handoff.R')
message('Rendered student HTML report.')

report_readme <- c(
  '# HT1080 add-back student handoff',
  '',
  'Open `HT1080_addback_student_report_2026-09-21.html` in a web browser.',
  'The HTML is self-contained. The `figures/` directory supplies every composite and component panel as PDF, SVG and high-resolution PNG; `data/` supplies complete CAMERA tables, selected DE tables, PCA coordinates, signature results and provenance; `code/` contains the analysis and report scripts.',
  'Exploratory add-back-responsive gene lists and methods are in data/outputs__responsive_genes__*. The report includes two baseline filters, TREAT sensitivity and three standalone figures. Run code/build_student_handoff.R to refresh these outputs and the ZIP using existing CAMERA results.',
  '',
  'Primary interpretation: UTY-containing established populations show partial molecular reversal of a frozen LOY-associated expression signature. KDM5D alone does not show a comparable retained response. Each condition is one established population sampled on three matched occasions, so inference is conditional on these populations.'
)
writeLines(report_readme, 'student_report/README.md', useBytes = TRUE)

verification <- c(
  sprintf('PASS CAMERA gene universe: %d unique mapped Entrez genes', nrow(v_camera)),
  sprintf('PASS CAMERA collections: %s', paste(sprintf('%s=%d sets', resource_summary$collection, resource_summary$sets_tested), collapse = '; ')),
  sprintf('PASS CAMERA contrasts: %d', length(contrast_labels)),
  sprintf('PASS CAMERA result rows: %d', nrow(camera_all)),
  sprintf('PASS report HTML: %s bytes', file.info('student_report/HT1080_addback_student_report_2026-09-21.html')$size),
  sprintf('PASS handoff data files: %d', length(list.files('student_report/data', full.names = TRUE)))
)
writeLines(verification, 'student_report/VERIFICATION.txt', useBytes = TRUE)
cat(paste(verification, collapse = '\n'), '\n')
finalize_handoff()
