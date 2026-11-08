script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
rmarkdown::render("report/batch_aware_de.Rmd",output_dir=normalizePath("outputs/report",winslash="/"),quiet=TRUE)
rmarkdown::render("report/analysis.Rmd",output_dir=normalizePath("outputs/report",winslash="/"),quiet=TRUE)
