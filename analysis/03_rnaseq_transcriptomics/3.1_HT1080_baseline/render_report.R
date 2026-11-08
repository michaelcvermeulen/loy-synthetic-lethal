# Rendering reads completed outputs and cannot refit models or alter result tables.
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
setwd(dirname(normalizePath(script,winslash='/')))
rmarkdown::render('report/analysis.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),quiet=TRUE)
rmarkdown::render('report/HDAC12_gene_drivers.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),quiet=TRUE)
rmarkdown::render('report/HDAC12_external_comparison.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),quiet=TRUE)
rmarkdown::render('report/epigenetic_deep_dive.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),quiet=TRUE)
rmarkdown::render('report/chromatin_GSEA_curves.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),quiet=TRUE)

rmarkdown::render('report/batch_aware_de.Rmd',output_dir=normalizePath('outputs/report',winslash='/'),quiet=TRUE)
