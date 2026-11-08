# Preserve the custom T2T mask. applyFilters() reconstructs its mask and an
# all-disabled call does NOT preserve an existing annotation's use column.
preserve_t2t_bin_mask <- function(rc, bins) {
  annotation <- if (is.data.frame(bins)) bins else Biobase::pData(bins)
  idx <- match(Biobase::featureNames(rc), rownames(annotation))
  stopifnot(!anyNA(idx), !anyDuplicated(rownames(annotation)))
  annotation <- annotation[idx, , drop=FALSE]
  for (field in c('chromosome','start','end','gc','bases','mappability','blacklist')) {
    stopifnot(isTRUE(all.equal(Biobase::fData(rc)[[field]], annotation[[field]], check.attributes=FALSE)))
  }
  mask <- annotation$use
  stopifnot(is.logical(mask), !anyNA(mask),
            identical(mask, annotation$bases >= 60 & annotation$mappability >= 55 & annotation$blacklist < 20),
            all(is.finite(annotation$gc[mask])))
  Biobase::fData(rc)$use <- mask
  rc$used.reads <- colSums(Biobase::assayDataElement(rc, 'counts')[mask,,drop=FALSE])
  # Same calculation as QDNAseq 1.46.0's internal expectedVariance().
  rc$expected.variance <- sum(mask)/rc$used.reads
  if ('paired.ends' %in% colnames(Biobase::pData(rc))) {
    rc$expected.variance <- rc$expected.variance * ifelse(rc$paired.ends,2,1)
  }
  stopifnot(identical(Biobase::fData(rc)$use,mask))
  rc
}
