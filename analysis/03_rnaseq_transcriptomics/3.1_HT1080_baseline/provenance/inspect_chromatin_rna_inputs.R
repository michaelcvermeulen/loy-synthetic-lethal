options(width = 110)
for (f in c('camera_voom_model', '01_ingest', '02_expression')) {
  x <- readRDS(paste0('outputs/models/', f, '.rds'))
  cat('\n', f, ': ', paste(names(x), collapse = ', '), '\n', sep = '')
  if (f == 'camera_voom_model') {
    print(dim(x$voom$E)); print(head(rownames(x$voom$E)))
    print(x$design); print(x$contrasts)
  }
  if (f == '01_ingest') print(x$core)
}
for (p in c('fgsea', 'limma', 'data.table', 'jsonlite', 'digest', 'BiocParallel')) {
  cat(p, as.character(packageVersion(p)), '\n')
}
print(args(fgsea::fgseaMultilevel))
