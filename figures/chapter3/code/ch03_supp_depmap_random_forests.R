# Recreate original ggplot2 histogram bin geometry from the frozen S04 table.
# Then run ch03_supp_depmap_random_forests.py to export the full figure.
args <- commandArgs(FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)]), winslash = "/")
root <- dirname(dirname(dirname(dirname(script))))
data_dir <- file.path(root, "figures/chapter3/source_data/depmap_random_forests")
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
confounder <- fread(file.path(data_dir, "fig8b_confounder_check.csv"))
p <- ggplot(confounder[role != "hit"], aes(success)) +
  geom_histogram(binwidth = 0.05, boundary = 0)
hist <- ggplot_build(p)$data[[1]][, c("xmin", "xmax", "count")]
stopifnot(sum(hist$count) == 200)
fwrite(hist, file.path(data_dir, "panel_b_histogram_geometry.csv"))
cat("Saved original histogram bin geometry.\n")
