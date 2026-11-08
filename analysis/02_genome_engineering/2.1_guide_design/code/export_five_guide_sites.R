# Export original saved site records without changing their scores or orientation.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
labels <- read.csv(args[[1]], stringsAsFactors = FALSE)
sites <- readRDS(args[[2]])
selected <- lapply(sites, function(z) {
  z <- as.data.frame(z)
  z[z$main_guide %in% labels$guide, , drop = FALSE]
})
selected <- do.call(rbind, selected)
stopifnot(setequal(unique(selected$main_guide), labels$guide))
selected$guide_id <- labels$guide_id[match(selected$main_guide, labels$guide)]
dir.create(args[[3]], recursive = TRUE, showWarnings = FALSE)
write.csv(selected, file.path(args[[3]], "saved_nonY_sites.csv"), row.names = FALSE)
print(table(selected$guide_id, selected$direction))
print(head(selected[selected$direction != "forward", ], 3))
cat("Biostrings available:", requireNamespace("Biostrings", quietly = TRUE), "\n")
cat("crisprScore available:", requireNamespace("crisprScore", quietly = TRUE), "\n")
