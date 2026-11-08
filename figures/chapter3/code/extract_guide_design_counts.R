# Read saved analysis outputs; does not run a genome search.
# Rscript --vanilla extract_guide_design_counts.R /path/to/LOY_paralogs/data /path/to/output
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply the original data directory and output directory")
data_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
recurrent <- readRDS(file.path(data_dir, "Y_seq_chm13.RDS"))
restricted <- readRDS(file.path(data_dir, "unique_Y_sequences_chm13.RDS"))
stopifnot(all(nchar(recurrent$V1) == 23L),
          all(nchar(restricted$sequence) == 23L),
          !anyDuplicated(recurrent$V1), !anyDuplicated(restricted$sequence),
          all(restricted$sequence %in% recurrent$V1))
high_copy <- restricted[restricted$count > 20, ]
counts <- data.frame(
  stage = c("recurrent", "Y_restricted", "high_copy"),
  count = c(nrow(recurrent), nrow(restricted), nrow(high_copy)),
  unit = "distinct 23-nt protospacer-PAM sequences",
  minimum_stored_Y_count = c(min(recurrent$count), min(restricted$count), min(high_copy$count)),
  criterion = c("Saved recurrent PAM-compatible input", "Zero exact non-Y matches in saved screen", "Y count > 20"),
  original_object = c("Y_seq_chm13.RDS", "unique_Y_sequences_chm13.RDS", "unique_Y_sequences_chm13.RDS")
)
write.csv(counts, file.path(output_dir, "sequence_filter_counts.csv"), row.names = FALSE)
print(counts, row.names = FALSE)
