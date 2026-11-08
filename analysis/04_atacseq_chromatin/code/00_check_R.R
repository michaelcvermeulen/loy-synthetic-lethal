script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
for (f in list.files("code", pattern = "\\.R$", full.names = TRUE)) {
  parse(f)
  cat("Parsed", f, "\n")
}
suppressPackageStartupMessages({ library(edgeR); library(csaw) })
print(formals(edgeR::glmQLFit.default))
print(formals(csaw::readParam))
