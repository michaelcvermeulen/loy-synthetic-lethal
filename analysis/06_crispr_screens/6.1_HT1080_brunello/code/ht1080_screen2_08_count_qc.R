#!/usr/bin/env Rscript

# HT1080 Screen 2 - Count Matrix QC Script
# Purpose: Comprehensive quality control of MAGeCK count matrices

# Load required libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

# Set working directory
setwd("/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2")

# Define paths
count_dir <- "mageck/count"
qc_dir <- "qc"

# Create output directory
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

cat("=========================================\n")
cat("Count Matrix QC Analysis - HT1080 Screen 2\n")
cat("Date:", date(), "\n")
cat("=========================================\n\n")

# Function to calculate Gini coefficient
gini <- function(x) {
  x <- sort(x)
  n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

# Function to analyze count matrix
analyze_count_matrix <- function(count_file, output_prefix, title) {
  if (!file.exists(count_file)) {
    cat("File not found:", count_file, "\n")
    return(NULL)
  }

  cat("\nAnalyzing:", count_file, "\n")
  cat(rep("-", 50), "\n", sep = "")

  # Read count data
  counts <- read.table(count_file, header = TRUE, sep = "\t", row.names = 1)
  samples <- colnames(counts)
  n_samples <- ncol(counts)
  n_sgrnas <- nrow(counts)

  # Initialize results data frame
  qc_results <- data.frame(
    Sample = samples,
    Total_Reads = numeric(n_samples),
    Detected_sgRNAs = numeric(n_samples),
    Detection_Rate = numeric(n_samples),
    Zero_Count_sgRNAs = numeric(n_samples),
    Zero_Rate = numeric(n_samples),
    Mean_Count = numeric(n_samples),
    Median_Count = numeric(n_samples),
    Gini_Coefficient = numeric(n_samples),
    CV = numeric(n_samples)
  )

  # Calculate metrics for each sample
  for (i in 1:n_samples) {
    sample_counts <- counts[, i]
    qc_results$Total_Reads[i] <- sum(sample_counts)
    qc_results$Detected_sgRNAs[i] <- sum(sample_counts > 0)
    qc_results$Detection_Rate[i] <- qc_results$Detected_sgRNAs[i] / n_sgrnas * 100
    qc_results$Zero_Count_sgRNAs[i] <- sum(sample_counts == 0)
    qc_results$Zero_Rate[i] <- qc_results$Zero_Count_sgRNAs[i] / n_sgrnas * 100
    qc_results$Mean_Count[i] <- mean(sample_counts)
    qc_results$Median_Count[i] <- median(sample_counts)
    qc_results$Gini_Coefficient[i] <- gini(sample_counts)
    qc_results$CV[i] <- sd(sample_counts) / mean(sample_counts)
  }

  # Print summary
  cat("\nQC Summary:\n")
  print(qc_results)

  # Save QC results
  write.table(qc_results,
              file.path(qc_dir, paste0(output_prefix, "_qc_metrics.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Create QC plots
  pdf(file.path(qc_dir, paste0(output_prefix, "_qc_plots.pdf")), width = 12, height = 10)

  # Plot 1: Library size
  p1 <- ggplot(qc_results, aes(x = Sample, y = Total_Reads / 1e6, fill = Sample)) +
    geom_bar(stat = "identity") +
    labs(title = "Library Size",
         x = "Sample",
         y = "Total Reads (millions)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Plot 2: Detection rate
  p2 <- ggplot(qc_results, aes(x = Sample, y = Detection_Rate, fill = Sample)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 90, linetype = "dashed", color = "red") +
    labs(title = "sgRNA Detection Rate",
         x = "Sample",
         y = "Detection Rate (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Plot 3: Gini coefficient
  p3 <- ggplot(qc_results, aes(x = Sample, y = Gini_Coefficient, fill = Sample)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "orange") +
    labs(title = "Gini Coefficient (Uniformity)",
         x = "Sample",
         y = "Gini Coefficient") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Plot 4: Zero inflation
  p4 <- ggplot(qc_results, aes(x = Sample, y = Zero_Rate, fill = Sample)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    labs(title = "Zero-Count sgRNAs",
         x = "Sample",
         y = "Zero Rate (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Combine plots
  grid.arrange(p1, p2, p3, p4, ncol = 2, top = paste(title, "- QC Metrics"))

  # Plot 5: Count distribution
  counts_long <- counts %>%
    mutate(sgRNA = row.names(.)) %>%
    pivot_longer(cols = -sgRNA, names_to = "Sample", values_to = "Count") %>%
    mutate(log_count = log2(Count + 1))

  p5 <- ggplot(counts_long, aes(x = log_count, color = Sample)) +
    geom_density(size = 1) +
    labs(title = paste(title, "- Count Distribution"),
         x = "Log2(Count + 1)",
         y = "Density") +
    theme_bw() +
    theme(legend.position = "bottom")

  print(p5)

  # Plot 6: Sample correlation heatmap
  cor_matrix <- cor(counts, method = "spearman", use = "complete.obs")

  pheatmap(cor_matrix,
           display_numbers = TRUE,
           number_format = "%.3f",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(min(cor_matrix), 1, length.out = 101),
           main = paste(title, "- Sample Correlation"),
           fontsize = 10,
           fontsize_number = 8)

  # Plot 7: PCA
  if (n_samples > 2) {
    # Normalize and log transform
    counts_norm <- log2(counts + 1)

    # PCA
    pca_result <- prcomp(t(counts_norm), scale. = TRUE)

    # Create PCA data frame
    pca_df <- data.frame(
      Sample = rownames(pca_result$x),
      PC1 = pca_result$x[,1],
      PC2 = pca_result$x[,2]
    )

    # Add condition information
    pca_df$Condition <- ifelse(grepl("LOY", pca_df$Sample), "LOY", "CTRL")
    pca_df$Timepoint <- ifelse(grepl("D0", pca_df$Sample), "D0", "POST")

    # Calculate variance explained
    var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Timepoint)) +
      geom_point(size = 4) +
      geom_text(aes(label = Sample), hjust = 0, vjust = -1, size = 3) +
      labs(title = paste(title, "- PCA Plot"),
           x = paste0("PC1 (", var_explained[1], "%)"),
           y = paste0("PC2 (", var_explained[2], "%)")) +
      theme_bw()

    print(p_pca)
  }

  dev.off()

  cat("QC plots saved to:", file.path(qc_dir, paste0(output_prefix, "_qc_plots.pdf")), "\n")

  return(qc_results)
}

# Function to compare R1 and R2
compare_r1_r2 <- function(r1_file, r2_file) {
  if (!file.exists(r1_file) || !file.exists(r2_file)) {
    cat("R1 or R2 file not found for comparison\n")
    return(NULL)
  }

  cat("\nComparing R1 and R2 counts...\n")
  cat(rep("-", 50), "\n", sep = "")

  # Read count data
  r1_counts <- read.table(r1_file, header = TRUE, sep = "\t", row.names = 1)
  r2_counts <- read.table(r2_file, header = TRUE, sep = "\t", row.names = 1)

  # Ensure same sgRNA order
  common_sgrnas <- intersect(rownames(r1_counts), rownames(r2_counts))
  r1_counts <- r1_counts[common_sgrnas, ]
  r2_counts <- r2_counts[common_sgrnas, ]

  # Calculate correlations for each sample
  correlations <- numeric(ncol(r1_counts))
  names(correlations) <- colnames(r1_counts)

  for (i in 1:ncol(r1_counts)) {
    correlations[i] <- cor(r1_counts[,i], r2_counts[,i],
                           method = "spearman", use = "complete.obs")
  }

  # Create comparison plot
  pdf(file.path(qc_dir, "R1_R2_comparison.pdf"), width = 10, height = 8)

  # Correlation barplot
  cor_df <- data.frame(Sample = names(correlations),
                       Correlation = correlations)

  p_cor <- ggplot(cor_df, aes(x = Sample, y = Correlation, fill = Sample)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    ylim(0, 1) +
    labs(title = "R1 vs R2 Correlation by Sample",
         x = "Sample",
         y = "Spearman Correlation") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  print(p_cor)

  # Scatter plots for each sample
  for (i in 1:min(6, ncol(r1_counts))) {
    sample_name <- colnames(r1_counts)[i]

    df <- data.frame(
      R1 = log2(r1_counts[,i] + 1),
      R2 = log2(r2_counts[,i] + 1)
    )

    p <- ggplot(df, aes(x = R1, y = R2)) +
      geom_point(alpha = 0.3, size = 0.5) +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
      labs(title = paste(sample_name, "- R1 vs R2"),
           x = "R1 Log2(Count + 1)",
           y = "R2 Log2(Count + 1)") +
      annotate("text", x = Inf, y = -Inf,
               label = paste("r =", round(correlations[i], 3)),
               hjust = 1, vjust = -1) +
      theme_bw()

    print(p)
  }

  dev.off()

  cat("R1 vs R2 comparison saved to:", file.path(qc_dir, "R1_R2_comparison.pdf"), "\n")

  # Print correlation summary
  cat("\nR1 vs R2 Correlations:\n")
  print(cor_df)

  return(correlations)
}

# Main QC pipeline
cat("\n1. Analyzing R1 count matrix...\n")
r1_qc <- analyze_count_matrix(file.path(count_dir, "ht1080_screen2_R1.count.txt"),
                              "R1", "R1 Reads")

cat("\n2. Analyzing R2 count matrix...\n")
r2_qc <- analyze_count_matrix(file.path(count_dir, "ht1080_screen2_R2.count.txt"),
                              "R2", "R2 Reads")

cat("\n3. Comparing R1 and R2...\n")
r1_r2_cor <- compare_r1_r2(file.path(count_dir, "ht1080_screen2_R1.count.txt"),
                           file.path(count_dir, "ht1080_screen2_R2.count.txt"))

# Create summary report
cat("\n=========================================\n")
cat("Creating QC Summary Report...\n")
cat("=========================================\n")

summary_file <- file.path(qc_dir, "qc_summary_report.txt")
sink(summary_file)

cat("HT1080 Screen 2 - Count QC Summary\n")
cat("Date:", date(), "\n")
cat(rep("=", 50), "\n\n", sep = "")

cat("OVERALL STATISTICS:\n")
cat(rep("-", 30), "\n", sep = "")
cat("Total sgRNAs in library: 77,441\n")
cat("Number of samples: 6\n\n")

if (!is.null(r1_qc)) {
  cat("R1 SUMMARY:\n")
  cat(rep("-", 30), "\n", sep = "")
  cat("Average detection rate:", round(mean(r1_qc$Detection_Rate), 2), "%\n")
  cat("Average Gini coefficient:", round(mean(r1_qc$Gini_Coefficient), 3), "\n")
  cat("Average zero rate:", round(mean(r1_qc$Zero_Rate), 2), "%\n\n")
}

if (!is.null(r2_qc)) {
  cat("R2 SUMMARY:\n")
  cat(rep("-", 30), "\n", sep = "")
  cat("Average detection rate:", round(mean(r2_qc$Detection_Rate), 2), "%\n")
  cat("Average Gini coefficient:", round(mean(r2_qc$Gini_Coefficient), 3), "\n")
  cat("Average zero rate:", round(mean(r2_qc$Zero_Rate), 2), "%\n\n")
}

if (!is.null(r1_r2_cor)) {
  cat("R1 vs R2 CORRELATION:\n")
  cat(rep("-", 30), "\n", sep = "")
  cat("Mean correlation:", round(mean(r1_r2_cor), 3), "\n")
  cat("Min correlation:", round(min(r1_r2_cor), 3), "\n")
  cat("Max correlation:", round(max(r1_r2_cor), 3), "\n")
}

sink()

cat("\nQC summary report saved to:", summary_file, "\n")

cat("\n=========================================\n")
cat("Count QC Analysis Complete!\n")
cat("All results saved to:", qc_dir, "\n")
cat("=========================================\n")