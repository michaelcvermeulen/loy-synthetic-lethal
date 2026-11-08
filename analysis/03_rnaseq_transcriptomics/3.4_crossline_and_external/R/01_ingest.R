# Load the three per-line thesis outputs, verify inputs and define the common universes.
# Nothing here recomputes a per-line statistic: effects, FDRs, NES and BH are carried over.
LINES <- c('HT1080', 'A549', 'DU145')

line_path <- function(cfg, line, key) file.path(REPO, cfg$lines[[line]]$root, cfg$lines[[line]][[key]])

read_cpm <- function(cfg, line) {
  man <- as.data.frame(data.table::fread(line_path(cfg, line, 'manifest')))
  if (line == 'HT1080') {
    e <- readRDS(line_path(cfg, line, 'cpm_model'))
    cpm <- e$norm$CPM
    meta <- man[match(colnames(cpm), man$sample_key), c('sample_key', 'treatment', 'batch')]
  } else {
    long <- data.table::fread(line_path(cfg, line, 'cpm_long'))
    cpm <- as.matrix(data.table::dcast(long, gene_id ~ sample_key, value.var = 'CPM'), rownames = 'gene_id')
    meta <- man[match(colnames(cpm), man$sample_key), c('sample_key', 'treatment', 'batch')]
  }
  assert(!anyNA(meta$treatment) && all(meta$treatment %in% c('CTRL', 'LOY')), paste('Unexpected CPM samples for', line))
  list(cpm = cpm, meta = meta)
}

# Condition means per sampling period, then an equal-weight average across periods (HT1080 has two).
condition_means <- function(cpm, meta, ids) {
  periods <- sort(unique(meta$batch))
  per <- lapply(periods, function(b) {
    m <- meta$batch == b
    data.frame(CTRL = rowMeans(cpm[ids, m & meta$treatment == 'CTRL', drop = FALSE]),
               LOY = rowMeans(cpm[ids, m & meta$treatment == 'LOY', drop = FALSE]))
  })
  list(CTRL = Reduce(`+`, lapply(per, `[[`, 'CTRL')) / length(per),
       LOY = Reduce(`+`, lapply(per, `[[`, 'LOY')) / length(per),
       periods = periods, per_period = setNames(per, periods))
}

ingest <- function(cfg) {
  inputs <- character()
  eff <- de <- fg <- tf <- cp <- list()
  for (L in LINES) {
    p <- line_path(cfg, L, 'effects'); inputs <- c(inputs, p)
    eff[[L]] <- as.data.frame(data.table::fread(p, select = c('gene_id', 'gene_name', 'chrom', 'start', 'end', 'PAR', 'arm',
                                                                  'expressed', 'primary_eligible', 'mean_difference')))
    p <- line_path(cfg, L, 'gene_de'); inputs <- c(inputs, p)
    d <- as.data.frame(data.table::fread(p))
    flag <- cfg$lines[[L]]$gene_thesis_flag
    de[[L]] <- data.frame(gene_id = d$gene_id, logFC = d$logFC, FDR = d$FDR,
                          thesis = if (is.null(flag)) FALSE else as.logical(d[[flag]]))
    p <- line_path(cfg, L, 'fgsea'); inputs <- c(inputs, p)
    f <- as.data.frame(data.table::fread(p))
    if (!is.null(cfg$lines[[L]]$fgsea_contrast)) {
      all_contrasts <- f
      f <- f[f$contrast == cfg$lines[[L]]$fgsea_contrast, ]
      fg[[paste0(L, '_periods')]] <- all_contrasts[all_contrasts$contrast != cfg$lines[[L]]$fgsea_contrast, ]
    }
    assert(!anyDuplicated(f$pathway), paste('Duplicate gene sets after contrast filter for', L))
    fg[[L]] <- f
    p <- line_path(cfg, L, 'collectri'); inputs <- c(inputs, p)
    t <- as.data.frame(data.table::fread(p))
    if (!is.null(cfg$lines[[L]]$fgsea_contrast) && 'contrast' %in% names(t)) t <- t[t$contrast == cfg$lines[[L]]$fgsea_contrast, ]
    assert(!anyDuplicated(t$source), paste('Duplicate regulators for', L))
    tf[[L]] <- t
    cp[[L]] <- read_cpm(cfg, L)
    inputs <- c(inputs, if (L == 'HT1080') line_path(cfg, L, 'cpm_model') else line_path(cfg, L, 'cpm_long'), line_path(cfg, L, 'manifest'))
  }
  # Common universes.
  genes <- Reduce(intersect, lapply(eff, function(d) d$gene_id[d$primary_eligible]))
  sets <- Reduce(intersect, lapply(fg[LINES], `[[`, 'pathway'))
  tab(data.frame(line = LINES, primary_genes = sapply(eff, function(d) sum(d$primary_eligible)),
                 significant_genes_FDR05 = sapply(de, function(d) sum(d$FDR < cfg$alpha, na.rm = TRUE)),
                 thesis_genes = sapply(de, function(d) sum(d$thesis, na.rm = TRUE)),
                 gene_sets_tested = sapply(fg[LINES], nrow),
                 gene_sets_BH05 = sapply(fg[LINES], function(d) sum(d$BH_across_collections < cfg$alpha, na.rm = TRUE))),
      'line_inputs_summary')
  tab(data.frame(universe = c('common_primary_genes', 'common_gene_sets'), n = c(length(genes), length(sets))), 'common_universes')
  extras <- c(file.path(REPO, cfg$msigdb_rds), file.path(REPO, cfg$grhl2_neighbourhood), cfg$external_registry)
  extras <- c(extras, reg_path(data.table::fread(cfg$external_registry)$path))   # external DE tables
  inputs <- unique(c(inputs, extras))
  tab(data.frame(path = inputs, sha256 = vapply(inputs, hash_file, character(1))), 'input_hashes')
  list(eff = eff, de = de, fg = fg, tf = tf, cp = cp, genes = genes, sets = sets)
}
