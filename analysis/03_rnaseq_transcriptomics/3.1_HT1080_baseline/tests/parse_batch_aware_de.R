script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script,winslash="/"))))
for(p in c("R/02_batch_aware_de.R","run_batch_aware_de.R","render_batch_aware_de.R","run.R","render_report.R")) invisible(parse(p))
cat("All modified/new R scripts parse successfully.\n")
