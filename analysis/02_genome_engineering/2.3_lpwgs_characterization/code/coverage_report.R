#!/usr/bin/env Rscript
# Coverage visualization: genome-wide + chrY zoom for A549 samples (T2T-CHM13v2.0)
# Outputs PNG plots + summary table to coverage/

library(tidyverse)
library(cowplot)

covdir <- "/global/home/hpc5173/projects/LOY/WGS/coverage"
samples <- c("A549-CTRL2", "A549-LOY", "A549-BC14")
sample_colors <- c("A549-CTRL2" = "#2166AC", "A549-LOY" = "#B2182B", "A549-BC14" = "#636363")

# ---- Load genome-wide summary ----
genome_summary <- map_dfr(samples, function(s) {
  read_tsv(file.path(covdir, paste0(s, ".genome.mosdepth.summary.txt")),
           show_col_types = FALSE) %>%
    mutate(sample = s)
})

# ---- Summary stats table ----
autosomes <- paste0("chr", 1:22)
stats <- genome_summary %>%
  filter(chrom %in% c(autosomes, "chrY", "total")) %>%
  group_by(sample) %>%
  summarise(
    genome_mean = mean[chrom == "total"],
    autosome_mean = sum(bases[chrom %in% autosomes]) / sum(length[chrom %in% autosomes]),
    chrY_mean = mean[chrom == "chrY"],
    chrY_autosome_ratio = chrY_mean / autosome_mean,
    .groups = "drop"
  )
write_tsv(stats, file.path(covdir, "coverage_summary.tsv"))
cat("\n--- Coverage Summary ---\n")
print(as.data.frame(stats), digits = 3)

# ---- Panel 1: Genome-wide bar chart ----
chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY")

genome_chrom <- genome_summary %>%
  filter(chrom %in% chrom_order) %>%
  mutate(chrom = factor(chrom, levels = chrom_order))

p1 <- ggplot(genome_chrom, aes(x = chrom, y = mean, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = sample_colors) +
  labs(x = NULL, y = "Mean Coverage",
       title = "Mean Coverage per Chromosome (T2T-CHM13v2.0)",
       fill = "Sample") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "top")

# ---- Load chrY 10kb-window coverage ----
chrY_cov <- map_dfr(samples, function(s) {
  read_tsv(file.path(covdir, paste0(s, ".chrY.regions.bed.gz")),
           col_names = c("chrom", "start", "end", "coverage"),
           show_col_types = FALSE) %>%
    mutate(sample = s,
           midpoint = (start + end) / 2 / 1e6)
})

# ---- Panel 2: chrY coverage track ----
p2 <- ggplot(chrY_cov, aes(x = midpoint, y = coverage, color = sample)) +
  geom_line(alpha = 0.7, size = 0.4) +
  scale_color_manual(values = sample_colors) +
  labs(x = "chrY Position (Mb)", y = "Coverage (10kb windows)",
       title = "chrY Coverage Along the Chromosome",
       color = "Sample") +
  theme_cowplot() +
  theme(legend.position = "top")

# ---- Panel 3: Log2 ratio (LOY / CTRL2) ----
ratio_wide <- chrY_cov[, c("start", "end", "midpoint", "sample", "coverage")]
ratio_wide <- as.data.frame(ratio_wide)
ratio_data <- reshape(ratio_wide, idvar = c("start", "end", "midpoint"),
                      timevar = "sample", v.names = "coverage", direction = "wide")
names(ratio_data) <- gsub("coverage\\.", "", names(ratio_data))
ratio_data$log2_ratio <- log2((ratio_data[["A549-LOY"]] + 0.01) /
                               (ratio_data[["A549-CTRL2"]] + 0.01))
ratio_data$log2_ratio <- pmax(pmin(ratio_data$log2_ratio, 5), -5)

p3 <- ggplot(ratio_data, aes(x = midpoint, y = log2_ratio)) +
  geom_line(color = "#B2182B", alpha = 0.6, size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "chrY Position (Mb)", y = "Log2(LOY / CTRL2)",
       title = "Log2 Coverage Ratio Along chrY (LOY vs CTRL2)") +
  theme_cowplot()

# ---- Save plots ----
ggsave(file.path(covdir, "genome_coverage_barplot.png"), p1,
       width = 14, height = 6, dpi = 200, bg = "white")

ggsave(file.path(covdir, "chrY_coverage_track.png"), p2,
       width = 14, height = 5, dpi = 200, bg = "white")

ggsave(file.path(covdir, "chrY_log2_ratio.png"), p3,
       width = 14, height = 4, dpi = 200, bg = "white")

# Combined multi-panel figure
p_combined <- plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1.2, 1, 0.8),
                        labels = c("A", "B", "C"))
ggsave(file.path(covdir, "coverage_report.png"), p_combined,
       width = 14, height = 14, dpi = 200, bg = "white")

cat("\nPlots saved to:", covdir, "\n")
cat("  - genome_coverage_barplot.png\n")
cat("  - chrY_coverage_track.png\n")
cat("  - chrY_log2_ratio.png\n")
cat("  - coverage_report.png (combined)\n")
cat("  - coverage_summary.tsv\n")
