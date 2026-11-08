normalize_counts <- function(counts, background, group, cfg) {
  # Filtering uses original all-feature library sizes. TMM is fitted on retained background.
  k <- edgeR::filterByExpr(counts, group = group, min.count = cfg$filter_min_count,
                         min.total.count = cfg$filter_min_total)
  keep <- k & background
  y <- edgeR::DGEList(counts[keep, , drop = FALSE])
  y <- edgeR::calcNormFactors(y, method = 'TMM')
  effective <- y$samples$lib.size * y$samples$norm.factors
  list(keep = keep, expression_keep = k, dge = y, effective = effective,
       logCPM = edgeR::cpm(counts, lib.size = effective, log = TRUE, prior.count = cfg$prior_count),
       CPM = edgeR::cpm(counts, lib.size = effective, log = FALSE))
}
