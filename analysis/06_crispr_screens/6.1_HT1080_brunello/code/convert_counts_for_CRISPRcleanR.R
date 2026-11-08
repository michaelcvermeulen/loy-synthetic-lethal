# Convert filtered counts to CRISPRcleanR format
# Converts sgRNA IDs from GENE_NUMBER_SEQUENCE to GENE_NUMBER format

library(tidyverse)

# Paths
input_file <- "A:/Dropbox/LOY_paralogs/projects/crispr_screen/ht1080_loy_screen_2/results/module6_d0_filtered_counts.csv"
output_file <- "A:/Dropbox/LOY_paralogs/projects/crispr_screen/ht1080_loy_screen_2/results/module6_d0_filtered_counts_CRISPRcleanR.csv"

# Read filtered counts
counts <- read_csv(input_file, show_col_types = FALSE)

cat("Input file:", input_file, "\n")
cat("Total sgRNAs:", nrow(counts), "\n")

# Convert sgRNA IDs: GENE_NUMBER_SEQUENCE -> GENE_NUMBER
# Extract first two underscore-separated parts
counts_cleanr <- counts %>%
  mutate(
    sgRNA = str_replace(sgRNA, "_[ACGT]+$", "")  # Remove trailing _SEQUENCE
  )

# Verify conversion
cat("\nExample conversions:\n")
head(data.frame(
  original = counts$sgRNA[1:5],
  converted = counts_cleanr$sgRNA[1:5]
))

# Check for duplicates (shouldn't be any)
n_dupes <- sum(duplicated(counts_cleanr$sgRNA))
if (n_dupes > 0) {
  warning(sprintf("Found %d duplicate sgRNA IDs after conversion!", n_dupes))
} else {
  cat("\nNo duplicate sgRNA IDs - conversion successful\n")
}

# Save output
write_csv(counts_cleanr, output_file)
cat("\nOutput saved to:", output_file, "\n")
cat("Total rows:", nrow(counts_cleanr), "\n")
