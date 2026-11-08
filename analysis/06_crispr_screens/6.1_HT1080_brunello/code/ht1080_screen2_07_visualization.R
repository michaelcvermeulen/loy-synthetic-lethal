#!/usr/bin/env Rscript

# HT1080 Screen 2 - Visualization Script
# Purpose: Generate publication-quality plots from MAGeCK results

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

# Set working directory
setwd("/global/home/hpc5173/projects/CRISPR_screen/analysis/ht1080_screen_2")

# Define paths
count_dir <- "mageck/count"
test_dir <- "mageck/test"
mle_dir <- "mageck_mle"
plot_dir <- "visualization"

# Create output directory
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Set theme for all plots
theme_set(theme_bw() +
          theme(panel.grid.minor = element_blank(),
                text = element_element_text(size = 12)))

cat("=========================================\n")
cat("HT1080 Screen 2 - Visualization Pipeline\n")
cat("Date:", date(), "\n")
cat("=========================================\n\n")

# Function to create volcano plot
create_volcano <- function(data_file, output_name, title, fdr_cutoff = 0.05, lfc_cutoff = 1) {
  if (!file.exists(data_file)) {
    cat("File not found:", data_file, "\n")
    return(NULL)
  }

  # Read data
  data <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Add significance categories
  data <- data %>%
    mutate(
      sig_category = case_when(
        neg.fdr < fdr_cutoff & neg.lfc < -lfc_cutoff ~ "Depleted",
        pos.fdr < fdr_cutoff & pos.lfc > lfc_cutoff ~ "Enriched",
        TRUE ~ "Not significant"
      ),
      label = ifelse(sig_category != "Not significant" &
                     (neg.fdr < 0.01 | pos.fdr < 0.01),
                     Gene, "")
    )

  # Create volcano plot
  p <- ggplot(data, aes(x = pos.lfc - neg.lfc,
                        y = -log10(pmin(neg.fdr, pos.fdr)),
                        color = sig_category)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(aes(label = label),
                    size = 3,
                    max.overlaps = 20,
                    segment.size = 0.3) +
    scale_color_manual(values = c("Depleted" = "blue",
                                  "Enriched" = "red",
                                  "Not significant" = "gray60")) +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed",
               color = "gray40",
               size = 0.5) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "gray40",
               size = 0.5) +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-log10(FDR)",
         color = "Category") +
    theme(legend.position = "bottom")

  # Save plot
  output_file <- file.path(plot_dir, paste0(output_name, "_volcano.pdf"))
  ggsave(output_file, p, width = 8, height = 6)
  cat("Created volcano plot:", output_file, "\n")

  return(p)
}

# Function to create ranked plot
create_ranked_plot <- function(data_file, output_name, title, metric = "pos.lfc") {
  if (!file.exists(data_file)) {
    cat("File not found:", data_file, "\n")
    return(NULL)
  }

  # Read data
  data <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Rank genes
  data <- data %>%
    arrange(desc(get(metric))) %>%
    mutate(rank = row_number(),
           label = ifelse(rank <= 10 | rank >= n() - 9, Gene, ""))

  # Create plot
  p <- ggplot(data, aes(x = rank, y = get(metric))) +
    geom_point(alpha = 0.6, size = 1) +
    geom_text_repel(aes(label = label),
                    size = 3,
                    max.overlaps = 20,
                    segment.size = 0.3) +
    labs(title = title,
         x = "Gene Rank",
         y = gsub("\\.", " ", toupper(metric))) +
    theme(panel.grid.major = element_line(color = "gray90"))

  # Save plot
  output_file <- file.path(plot_dir, paste0(output_name, "_ranked.pdf"))
  ggsave(output_file, p, width = 8, height = 6)
  cat("Created ranked plot:", output_file, "\n")

  return(p)
}

# Function to create correlation heatmap
create_correlation_heatmap <- function(count_file, output_name, title) {
  if (!file.exists(count_file)) {
    cat("File not found:", count_file, "\n")
    return(NULL)
  }

  # Read count data
  counts <- read.table(count_file, header = TRUE, sep = "\t", row.names = 1)

  # Calculate correlation
  cor_matrix <- cor(counts, method = "spearman", use = "complete.obs")

  # Create heatmap
  pdf(file.path(plot_dir, paste0(output_name, "_correlation.pdf")),
      width = 8, height = 7)

  pheatmap(cor_matrix,
           display_numbers = TRUE,
           number_format = "%.2f",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(0.8, 1, length.out = 101),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = title,
           fontsize = 10,
           fontsize_number = 8)

  dev.off()
  cat("Created correlation heatmap:",
      file.path(plot_dir, paste0(output_name, "_correlation.pdf")), "\n")
}

# Function to create sgRNA distribution plot
create_sgrna_distribution <- function(count_file, output_name, title) {
  if (!file.exists(count_file)) {
    cat("File not found:", count_file, "\n")
    return(NULL)
  }

  # Read count data
  counts <- read.table(count_file, header = TRUE, sep = "\t", row.names = 1)

  # Convert to long format
  counts_long <- counts %>%
    rownames_to_column("sgRNA") %>%
    pivot_longer(cols = -sgRNA, names_to = "Sample", values_to = "Count") %>%
    mutate(log_count = log2(Count + 1))

  # Create boxplot
  p1 <- ggplot(counts_long, aes(x = Sample, y = log_count, fill = Sample)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
    labs(title = paste(title, "- sgRNA Count Distribution"),
         x = "Sample",
         y = "Log2(Count + 1)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Create density plot
  p2 <- ggplot(counts_long, aes(x = log_count, color = Sample)) +
    geom_density(size = 1) +
    labs(title = paste(title, "- Count Density"),
         x = "Log2(Count + 1)",
         y = "Density") +
    theme(legend.position = "bottom")

  # Combine plots
  combined <- grid.arrange(p1, p2, ncol = 2)

  # Save plot
  output_file <- file.path(plot_dir, paste0(output_name, "_distribution.pdf"))
  ggsave(output_file, combined, width = 12, height = 5)
  cat("Created distribution plot:", output_file, "\n")
}

# Function to create comparison scatter plot
create_comparison_scatter <- function(file1, file2, label1, label2, output_name, title) {
  if (!file.exists(file1) || !file.exists(file2)) {
    cat("Files not found for comparison\n")
    return(NULL)
  }

  # Read data
  data1 <- read.table(file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data2 <- read.table(file2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Merge data
  merged <- inner_join(data1, data2, by = "Gene", suffix = c(".1", ".2"))

  # Create scatter plot
  p <- ggplot(merged, aes(x = pos.lfc.1 - neg.lfc.1,
                          y = pos.lfc.2 - neg.lfc.2)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    labs(title = title,
         x = paste(label1, "Log2 FC"),
         y = paste(label2, "Log2 FC")) +
    annotate("text", x = Inf, y = -Inf,
             label = paste("R =", round(cor(merged$pos.lfc.1 - merged$neg.lfc.1,
                                            merged$pos.lfc.2 - merged$neg.lfc.2,
                                            use = "complete.obs"), 3)),
             hjust = 1, vjust = -1)

  # Save plot
  output_file <- file.path(plot_dir, paste0(output_name, "_scatter.pdf"))
  ggsave(output_file, p, width = 6, height = 6)
  cat("Created comparison scatter plot:", output_file, "\n")

  return(p)
}

# Main visualization pipeline
cat("\n1. Creating volcano plots...\n")
cat("----------------------------\n")

# Create volcano plots for main comparisons
create_volcano(file.path(test_dir, "D0_LOY_vs_CTRL_R1.gene_summary.txt"),
               "D0_LOY_vs_CTRL",
               "D0: LOY vs CTRL")

create_volcano(file.path(test_dir, "POST_LOY_vs_CTRL_R1.gene_summary.txt"),
               "POST_LOY_vs_CTRL",
               "Post-selection: LOY vs CTRL")

cat("\n2. Creating ranked plots...\n")
cat("---------------------------\n")

# Create ranked plots
create_ranked_plot(file.path(test_dir, "POST_LOY_vs_CTRL_R1.gene_summary.txt"),
                   "POST_LOY_vs_CTRL",
                   "Post-selection: Gene Ranking by Enrichment",
                   metric = "pos.lfc")

cat("\n3. Creating correlation heatmaps...\n")
cat("------------------------------------\n")

# Create correlation heatmap
create_correlation_heatmap(file.path(count_dir, "ht1080_screen2_R1.count.txt"),
                           "samples_R1",
                           "Sample Correlation (R1)")

cat("\n4. Creating sgRNA distribution plots...\n")
cat("----------------------------------------\n")

# Create distribution plots
create_sgrna_distribution(file.path(count_dir, "ht1080_screen2_R1.count.txt"),
                         "R1",
                         "R1 Reads")

cat("\n5. Creating comparison plots...\n")
cat("--------------------------------\n")

# Compare R1 vs R2 if both exist
if (file.exists(file.path(test_dir, "POST_LOY_vs_CTRL_R1.gene_summary.txt")) &&
    file.exists(file.path(test_dir, "POST_LOY_vs_CTRL_R2.gene_summary.txt"))) {

  create_comparison_scatter(file.path(test_dir, "POST_LOY_vs_CTRL_R1.gene_summary.txt"),
                           file.path(test_dir, "POST_LOY_vs_CTRL_R2.gene_summary.txt"),
                           "R1", "R2",
                           "R1_vs_R2_comparison",
                           "R1 vs R2 Concordance")
}

# Compare D0 vs POST
create_comparison_scatter(file.path(test_dir, "D0_LOY_vs_CTRL_R1.gene_summary.txt"),
                         file.path(test_dir, "POST_LOY_vs_CTRL_R1.gene_summary.txt"),
                         "D0", "POST",
                         "D0_vs_POST_comparison",
                         "Temporal Comparison: D0 vs Post-selection")

cat("\n=========================================\n")
cat("Visualization Pipeline Complete!\n")
cat("All plots saved to:", plot_dir, "\n")
cat("=========================================\n")

# Create summary of generated plots
cat("\nGenerated plots:\n")
list.files(plot_dir, pattern = "\\.pdf$", full.names = FALSE) %>%
  paste("-", .) %>%
  cat(sep = "\n")

cat("\n\nVisualization complete at:", date(), "\n")