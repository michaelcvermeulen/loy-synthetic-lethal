# All paths are resolved by run.R; original inputs are read only.
tab <- function(x, name) {
  path <- file.path('outputs/tables', paste0(name, '.csv'))
  data.table::fwrite(as.data.frame(x), path, na = 'NA')
  invisible(path)
}
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = 'sha256')
stable_id <- function(x) sub('\\.[0-9]+(_PAR_Y)?$', '\\1', x)
matrix_long <- function(x, value = 'value', id = 'gene_id') {
  d <- data.table::as.data.table(x, keep.rownames = id)
  data.table::melt(d, id.vars = id, variable.name = 'sample_key',
                  value.name = value, variable.factor = FALSE)
}
align_samples <- function(x, meta) {
  assert(!anyDuplicated(meta$sample_key), 'Duplicate sample metadata keys')
  assert(setequal(colnames(x), meta$sample_key), 'Count/metadata sample mismatch')
  x[, meta$sample_key, drop = FALSE]
}
pair_indices <- function(meta) {
  blocks <- sort(unique(meta$passage_day_block))
  out <- lapply(blocks, function(b) {
    m <- meta[meta$passage_day_block == b, ]
    assert(nrow(m) == 2 && setequal(m$treatment, c('CTRL', 'LOY')),
           paste('Incomplete or duplicated passage match:', b))
    assert(length(unique(m$batch)) == 1 && length(unique(m$passage_number)) == 1,
           paste('Inconsistent passage metadata:', b))
    data.frame(block = b, batch = m$batch[1], passage = m$passage_number[1],
               CTRL = m$sample_key[m$treatment == 'CTRL'],
               LOY = m$sample_key[m$treatment == 'LOY'])
  })
  do.call(rbind, out)
}
matched_differences <- function(x, pairs) {
  assert(all(c(pairs$LOY, pairs$CTRL) %in% colnames(x)), 'Missing paired samples')
  z <- x[, pairs$LOY, drop = FALSE] - x[, pairs$CTRL, drop = FALSE]
  colnames(z) <- pairs$block
  z
}
period_average <- function(d, pairs) {
  periods <- sort(unique(pairs$batch))
  ans <- do.call(cbind, lapply(periods, function(b)
    rowMeans(d[, pairs$batch == b, drop = FALSE])))
  colnames(ans) <- periods
  rownames(ans) <- rownames(d)
  invisible(ans)
}
choose_symbols <- function(annotation, mean_expression, eligible) {
  a <- annotation[eligible & !is.na(annotation$gene_name) & nzchar(annotation$gene_name), ]
  a$mean_logCPM <- mean_expression[match(a$gene_id, annotation$gene_id)]
  a <- a[order(a$gene_name, -a$mean_logCPM, a$gene_id), ]
  a$representative <- !duplicated(a$gene_name)
  a
}
relative_pair_total <- function(x_loy, y_loy, x_ctrl, y_ctrl) {
  denom <- x_ctrl + y_ctrl
  ifelse(denom > 0, (x_loy + y_loy) / denom, NA_real_)
}
stage <- function(name, fn) {
  message('\n[', format(Sys.time(), '%H:%M:%S'), '] ', name)
  ans <- fn()
  saveRDS(ans, file.path('outputs/models', paste0(name, '.rds')))
  invisible(ans)
}
