#!/usr/bin/env Rscript

# HT1080 Screen 2 Count Matrix QC
# Purpose: Perform quality control on MAGeCK count matrix

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(pheatmap)
  library(corrplot)
  library(scales)
})

# Set working directory
setwd("/global/home/hpc5173/projects/CRISPR_screen")

# Define paths
count_file <- "analysis/ht1080_screen_2/mageck/count/ht1080_screen2_R1.count.txt"
output_dir <- "analysis/ht1080_screen_2/qc/count_qc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("===============================================\n")
cat("HT1080 Screen 2 - Count Matrix Quality Control\n")
cat("===============================================\n\n")

# Read count matrix
cat("Reading count matrix...\n")
counts <- read.table(count_file, header = TRUE, sep = "\t", row.names = 1)

# Extract sample columns (skip gene column)
sample_cols <- grep("^(D0|POST)", colnames(counts), value = TRUE)
count_mat <- as.matrix(counts[, sample_cols])

cat(paste("Loaded", nrow(count_mat), "sgRNAs across", ncol(count_mat), "samples\n\n"))

# Basic statistics
cat("Basic Statistics:\n")
cat("-----------------\n")

# Calculate library sizes
lib_sizes <- colSums(count_mat)
cat("\nLibrary sizes (total reads per sample):\n")
print(data.frame(Sample = names(lib_sizes),
                  TotalReads = format(lib_sizes, big.mark = ","),
                  MillionReads = round(lib_sizes/1e6, 1)))

# Calculate mapping rates
cat("\nMapping rate summary:\n")
# This would need the original read counts from logs
# For now, we'll calculate based on expected ~21% mapping rate

# sgRNA detection
detected <- colSums(count_mat > 0)
detection_pct <- round(100 * detected / nrow(count_mat), 1)
cat("\nsgRNA detection:\n")
print(data.frame(Sample = names(detected),
                  Detected = detected,
                  Percentage = detection_pct))

# Zero counts analysis
zero_counts <- rowSums(count_mat == 0)
cat("\nsgRNAs with zero counts:\n")
cat(paste("  In all samples:", sum(zero_counts == ncol(count_mat)), "\n"))
cat(paste("  In at least one sample:", sum(zero_counts > 0), "\n"))
cat(paste("  In no samples:", sum(zero_counts == 0), "\n"))

# Create QC plots
cat("\nGenerating QC plots...\n")

# 1. Library size barplot
pdf(file.path(output_dir, "library_sizes.pdf"), width = 10, height = 6)
par(mar = c(8, 5, 3, 2))
barplot(lib_sizes/1e6,
        names.arg = names(lib_sizes),
        las = 2,
        ylab = "Total reads (millions)",
        main = "Library Sizes",
        col = c(rep("skyblue", 2), rep("coral", 4)))
legend("topright", legend = c("Day 0", "Post-treatment"),
       fill = c("skyblue", "coral"))
dev.off()

# 2. sgRNA detection rate
pdf(file.path(output_dir, "sgrna_detection.pdf"), width = 10, height = 6)
par(mar = c(8, 5, 3, 2))
barplot(detection_pct,
        names.arg = names(detection_pct),
        las = 2,
        ylab = "Detection rate (%)",
        main = "sgRNA Detection Rate",
        col = c(rep("skyblue", 2), rep("coral", 4)),
        ylim = c(0, 100))
abline(h = c(90, 95, 99), lty = 2, col = "gray")
dev.off()

# 3. Count distribution boxplot (log scale)
pdf(file.path(output_dir, "count_distribution.pdf"), width = 10, height = 6)
log_counts <- log10(count_mat + 1)
boxplot(log_counts,
        las = 2,
        ylab = "log10(counts + 1)",
        main = "Count Distribution",
        col = c(rep("skyblue", 2), rep("coral", 4)))
dev.off()

# 4. Correlation heatmap
pdf(file.path(output_dir, "sample_correlation.pdf"), width = 8, height = 7)
cor_mat <- cor(log_counts, method = "spearman")
pheatmap(cor_mat,
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "Sample Correlation (Spearman)",
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

# 5. PCA plot
pdf(file.path(output_dir, "pca_plot.pdf"), width = 8, height = 6)
pca_res <- prcomp(t(log_counts), scale = TRUE)
pca_data <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Sample = rownames(pca_res$x),
  Group = ifelse(grepl("^D0", rownames(pca_res$x)), "Day0",
                 ifelse(grepl("CTRL", rownames(pca_res$x)), "POST_CTRL", "POST_LOY"))
)

var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(hjust = -0.2, vjust = 0.5, size = 3) +
  xlab(paste0("PC1 (", var_explained[1], "%)")) +
  ylab(paste0("PC2 (", var_explained[2], "%)")) +
  ggtitle("PCA of sgRNA Counts") +
  theme_classic() +
  scale_color_manual(values = c("Day0" = "blue", "POST_CTRL" = "green", "POST_LOY" = "red"))
ggsave(file.path(output_dir, "pca_plot_ggplot.pdf"), width = 8, height = 6)
dev.off()

# 6. Replicate correlation scatter plots
pdf(file.path(output_dir, "replicate_correlation.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

# CTRL replicates
if ("POST_CTRL_REP1" %in% colnames(count_mat) & "POST_CTRL_REP2" %in% colnames(count_mat)) {
  plot(log10(count_mat[,"POST_CTRL_REP1"] + 1),
       log10(count_mat[,"POST_CTRL_REP2"] + 1),
       pch = 16, cex = 0.3, col = rgb(0,0,1,0.3),
       xlab = "POST_CTRL_REP1 (log10)",
       ylab = "POST_CTRL_REP2 (log10)",
       main = paste("CTRL Replicates (r =",
                    round(cor(log10(count_mat[,"POST_CTRL_REP1"] + 1),
                             log10(count_mat[,"POST_CTRL_REP2"] + 1)), 3), ")"))
  abline(0, 1, col = "red", lty = 2)
}

# LOY replicates
if ("POST_LOY_REP1" %in% colnames(count_mat) & "POST_LOY_REP2" %in% colnames(count_mat)) {
  plot(log10(count_mat[,"POST_LOY_REP1"] + 1),
       log10(count_mat[,"POST_LOY_REP2"] + 1),
       pch = 16, cex = 0.3, col = rgb(1,0,0,0.3),
       xlab = "POST_LOY_REP1 (log10)",
       ylab = "POST_LOY_REP2 (log10)",
       main = paste("LOY Replicates (r =",
                    round(cor(log10(count_mat[,"POST_LOY_REP1"] + 1),
                             log10(count_mat[,"POST_LOY_REP2"] + 1)), 3), ")"))
  abline(0, 1, col = "red", lty = 2)
}
dev.off()

# 7. Gini coefficient (measure of count inequality)
gini_coef <- function(x) {
  x <- sort(x)
  n <- length(x)
  2 * sum(x * 1:n) / (n * sum(x)) - (n + 1) / n
}

gini_values <- apply(count_mat, 2, gini_coef)
cat("\nGini coefficients (count inequality):\n")
print(data.frame(Sample = names(gini_values),
                  Gini = round(gini_values, 3)))

pdf(file.path(output_dir, "gini_coefficient.pdf"), width = 10, height = 6)
par(mar = c(8, 5, 3, 2))
barplot(gini_values,
        names.arg = names(gini_values),
        las = 2,
        ylab = "Gini coefficient",
        main = "Count Distribution Inequality",
        col = c(rep("skyblue", 2), rep("coral", 4)),
        ylim = c(0, 1))
abline(h = 0.5, lty = 2, col = "gray")
dev.off()

# Save summary statistics
summary_stats <- data.frame(
  Sample = colnames(count_mat),
  TotalReads = lib_sizes,
  DetectedSgRNAs = detected,
  DetectionRate = detection_pct,
  MedianCount = apply(count_mat, 2, median),
  MeanCount = round(apply(count_mat, 2, mean), 1),
  GiniCoef = round(gini_values, 3)
)

write.csv(summary_stats,
          file.path(output_dir, "count_qc_summary.csv"),
          row.names = FALSE)

cat("\n===============================================\n")
cat("Count QC Complete!\n")
cat("===============================================\n")
cat(paste("\nResults saved to:", output_dir, "\n"))
cat("\nQC Summary:\n")
cat("-----------\n")

# Overall QC assessment
qc_pass <- TRUE
qc_warnings <- c()

# Check library sizes
if (any(lib_sizes < 10e6)) {
  qc_warnings <- c(qc_warnings, "- Some samples have < 10M mapped reads")
  qc_pass <- FALSE
}

# Check detection rates
if (any(detection_pct < 90)) {
  qc_warnings <- c(qc_warnings, "- Some samples have < 90% sgRNA detection")
  qc_pass <- FALSE
}

# Check replicate correlation
if ("POST_CTRL_REP1" %in% colnames(count_mat) & "POST_CTRL_REP2" %in% colnames(count_mat)) {
  ctrl_cor <- cor(log10(count_mat[,"POST_CTRL_REP1"] + 1),
                  log10(count_mat[,"POST_CTRL_REP2"] + 1))
  if (ctrl_cor < 0.9) {
    qc_warnings <- c(qc_warnings, paste("- CTRL replicate correlation is low:", round(ctrl_cor, 3)))
    qc_pass <- FALSE
  }
}

if ("POST_LOY_REP1" %in% colnames(count_mat) & "POST_LOY_REP2" %in% colnames(count_mat)) {
  loy_cor <- cor(log10(count_mat[,"POST_LOY_REP1"] + 1),
                 log10(count_mat[,"POST_LOY_REP2"] + 1))
  if (loy_cor < 0.9) {
    qc_warnings <- c(qc_warnings, paste("- LOY replicate correlation is low:", round(loy_cor, 3)))
    qc_pass <- FALSE
  }
}

if (qc_pass) {
  cat("✓ All QC metrics PASS\n")
} else {
  cat("⚠ QC warnings:\n")
  for (warning in qc_warnings) {
    cat(paste(warning, "\n"))
  }
}

cat("\nNext steps:\n")
cat("1. Review QC plots and metrics\n")
cat("2. If QC passes, proceed with MAGeCK test\n")
cat("3. Run MAGeCK-MLE for pathway analysis\n")