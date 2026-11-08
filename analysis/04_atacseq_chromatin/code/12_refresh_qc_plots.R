script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setwd(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
suppressPackageStartupMessages(library(ggplot2))
pc <- read.csv("qc/pca_scores.csv")
variance <- vapply(pc[grep("^PC[0-9]+$", names(pc))], var, numeric(1))
variance <- 100*variance/sum(variance)
p <- ggplot(pc, aes(PC1, PC2, colour = condition, label = sample)) +
  geom_point(size = 3) + geom_text(nudge_y = 2, show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  scale_colour_manual(values = c(CTRL = "#0072B2", LOY = "#D55E00")) +
  theme_bw(base_size = 12) + labs(x = sprintf("PC1 (%.1f%%)", variance[1]),
                                  y = sprintf("PC2 (%.1f%%)", variance[2]))
ggsave("qc/pca.png", p, width = 8, height = 6, dpi = 200)
mds <- read.csv("qc/mds_scores.csv")
p <- ggplot(mds, aes(dimension1, dimension2, colour = condition, label = sample)) +
  geom_point(size = 3) + geom_text(vjust = -0.8, show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  scale_colour_manual(values = c(CTRL = "#0072B2", LOY = "#D55E00")) +
  theme_bw(base_size = 12) + labs(x = "MDS dimension 1", y = "MDS dimension 2")
ggsave("qc/mds.png", p, width = 8, height = 6, dpi = 200)
