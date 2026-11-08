#!/usr/bin/env Rscript
# Analyze D0 count distribution and identify outliers

library(tidyverse)
library(ggplot2)

# Load count matrix
counts <- read.table("/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2/mageck/count/ht1080_screen2_R1.count.txt",
                     header = TRUE, sep = "\t", row.names = 1)

# Extract D0 and POST samples
d0_cols <- c("D0_CTRL", "D0_LOY")
post_cols <- c("POST_CTRL_REP1", "POST_CTRL_REP2", "POST_LOY_REP1", "POST_LOY_REP2")

d0_counts <- counts[, d0_cols]
post_counts <- counts[, post_cols]

# Calculate CPM for each sample
calculate_cpm <- function(counts) {
  sweep(counts, 2, colSums(counts), "/") * 1e6
}

d0_cpm <- calculate_cpm(d0_counts)
post_cpm <- calculate_cpm(post_counts)

# Calculate log2 CPM
d0_log2cpm <- log2(d0_cpm + 1)
post_log2cpm <- log2(post_cpm + 1)

# Print summary statistics
cat("===== D0 SAMPLE STATISTICS =====\n")
cat("\nD0-CTRL:\n")
cat(sprintf("  Total reads: %s\n", format(sum(d0_counts$D0_CTRL), big.mark = ",")))
cat(sprintf("  Non-zero guides: %d (%.1f%%)\n",
            sum(d0_counts$D0_CTRL > 0),
            100 * sum(d0_counts$D0_CTRL > 0) / nrow(d0_counts)))
cat(sprintf("  Mean count: %.1f\n", mean(d0_counts$D0_CTRL)))
cat(sprintf("  Median count: %.1f\n", median(d0_counts$D0_CTRL)))
cat(sprintf("  Max count: %d\n", max(d0_counts$D0_CTRL)))

cat("\nD0-LOY:\n")
cat(sprintf("  Total reads: %s\n", format(sum(d0_counts$D0_LOY), big.mark = ",")))
cat(sprintf("  Non-zero guides: %d (%.1f%%)\n",
            sum(d0_counts$D0_LOY > 0),
            100 * sum(d0_counts$D0_LOY > 0) / nrow(d0_counts)))
cat(sprintf("  Mean count: %.1f\n", mean(d0_counts$D0_LOY)))
cat(sprintf("  Median count: %.1f\n", median(d0_counts$D0_LOY)))
cat(sprintf("  Max count: %d\n", max(d0_counts$D0_LOY)))

# Identify outliers using z-score
identify_outliers <- function(log2cpm_vec, z_threshold = 3) {
  z_scores <- (log2cpm_vec - mean(log2cpm_vec)) / sd(log2cpm_vec)
  return(abs(z_scores) > z_threshold)
}

# Identify outliers in D0 samples
d0_ctrl_outliers <- identify_outliers(d0_log2cpm$D0_CTRL)
d0_loy_outliers <- identify_outliers(d0_log2cpm$D0_LOY)

cat("\n===== OUTLIER ANALYSIS (Z-score > 3) =====\n")
cat(sprintf("D0-CTRL outliers: %d guides (%.2f%%)\n",
            sum(d0_ctrl_outliers), 100 * sum(d0_ctrl_outliers) / length(d0_ctrl_outliers)))
cat(sprintf("D0-LOY outliers: %d guides (%.2f%%)\n",
            sum(d0_loy_outliers), 100 * sum(d0_loy_outliers) / length(d0_loy_outliers)))

# Check if outliers are consistent between replicates
both_outliers <- d0_ctrl_outliers & d0_loy_outliers
cat(sprintf("Consistent outliers (both samples): %d guides\n", sum(both_outliers)))

# List top 20 outlier guides
cat("\n===== TOP 20 HIGH-ABUNDANCE GUIDES IN D0-CTRL =====\n")
top_guides <- head(order(d0_counts$D0_CTRL, decreasing = TRUE), 20)
for (i in 1:20) {
  idx <- top_guides[i]
  guide_name <- rownames(d0_counts)[idx]
  ctrl_count <- d0_counts$D0_CTRL[idx]
  loy_count <- d0_counts$D0_LOY[idx]
  ctrl_cpm <- d0_cpm$D0_CTRL[idx]

  # Check if also high in POST samples
  post_avg <- mean(as.numeric(post_counts[idx, ]))

  cat(sprintf("%2d. %s: CTRL=%d (%.1f CPM), LOY=%d, POST_avg=%.1f\n",
              i, guide_name, ctrl_count, ctrl_cpm, loy_count, post_avg))
}

# Calculate filtering impact
cat("\n===== FILTERING IMPACT ANALYSIS =====\n")
thresholds <- c(0.5, 1.0, 2.0, 5.0)
for (thresh in thresholds) {
  low_guides <- d0_log2cpm$D0_CTRL < thresh | d0_log2cpm$D0_LOY < thresh
  cat(sprintf("Guides < %.1f log2 CPM: %d (%.1f%%)\n",
              thresh, sum(low_guides), 100 * sum(low_guides) / length(low_guides)))
}

# Check for PCR duplicates/artifacts
cat("\n===== PCR ARTIFACT DETECTION =====\n")
# Calculate coefficient of variation between D0 replicates
d0_cv <- apply(d0_counts, 1, function(x) sd(x) / mean(x))
high_cv_guides <- sum(d0_cv > 2, na.rm = TRUE)
cat(sprintf("Guides with CV > 2 between D0 samples: %d (%.1f%%)\n",
            high_cv_guides, 100 * high_cv_guides / length(d0_cv)))

# Compare D0 vs POST distributions
cat("\n===== D0 vs POST COMPARISON =====\n")
d0_mean_log2cpm <- rowMeans(d0_log2cpm)
post_mean_log2cpm <- rowMeans(post_log2cpm)

# Correlation
cor_d0_post <- cor(d0_mean_log2cpm, post_mean_log2cpm, method = "spearman")
cat(sprintf("Spearman correlation (D0 vs POST): %.3f\n", cor_d0_post))

# Guides depleted in POST
depleted <- d0_mean_log2cpm - post_mean_log2cpm > 2
cat(sprintf("Guides depleted in POST (>2 log2 fold): %d\n", sum(depleted)))

# Save results for plotting
results <- list(
  d0_counts = d0_counts,
  d0_log2cpm = d0_log2cpm,
  post_log2cpm = post_log2cpm,
  outliers = data.frame(
    guide = rownames(d0_counts),
    d0_ctrl_outlier = d0_ctrl_outliers,
    d0_loy_outlier = d0_loy_outliers
  )
)

saveRDS(results, "analysis/ht1080_screen_2/d0_distribution_analysis.rds")