packages <- c('chromVAR', 'JASPAR2024', 'limma', 'motifmatchr', 'TFBSTools', 'BiocManager')
for (p in packages) cat(p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else 'MISSING', '\n')
cat('R libraries:', .libPaths(), '\n')
