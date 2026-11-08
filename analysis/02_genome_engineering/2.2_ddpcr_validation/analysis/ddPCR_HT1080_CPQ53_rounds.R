# ---------------------------------------------------------------------------
# HT1080 LOY after successive CP+Q53 nucleofection rounds
#
# Figure: % LOY in HT1080 for Parental (unedited), Round 1 and Round 2 of
# CP+Q53 nucleofection, shown independently in both assay chemistries.
#
# Input : projects/ddPCR/data/summary_table_ddPCR.csv  (read-only)
# Output: projects/ddPCR/figures/HT1080_CPQ53_LOY_rounds.{png,pdf}
#         projects/ddPCR/figures/HT1080_CPQ53_LOY_rounds_source.csv
# ---------------------------------------------------------------------------

library(tidyverse)

proj_dir <- "A:/Dropbox/LOY_paralogs/projects/ddPCR"
fig_dir  <- file.path(proj_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Canonical naming
#
# The master table accumulated several naming schemes over two years. These
# helpers collapse them onto one vocabulary so a sample can be grouped
# regardless of which run it came from. Reusable for any HT1080 figure.
# ---------------------------------------------------------------------------

# Version encodes nucleofection round: 0 = unedited, 1 = first, 2 = second.
canon_round <- function(version) {
  factor(dplyr::recode(as.character(version),
                       "0" = "Parental",
                       "1" = "Round 1",
                       "2" = "Round 2"),
         levels = c("Parental", "Round 1", "Round 2"))
}

# Y-target + reference gene -> a single assay label.
canon_assay <- function(target, ref) {
  factor(paste0(target, "/", ref),
         levels = c("DDX3Y/RNaseP", "AMELY/AMELX"))
}

# Strip the descriptive suffixes that crept in across runs so that e.g.
# "CP LOY pool", "Q53+CP" and "CP+Q53" all resolve to one arm.
canon_arm <- function(name) {
  x <- str_trim(name)
  x <- str_replace_all(x, regex("\\s*LOY\\s*(pool|clone)\\s*$", ignore_case = TRUE), "")
  x <- str_replace_all(x, regex("^HT1080[_ ]+", ignore_case = TRUE), "")
  x <- str_trim(x)
  # order-insensitive: Q53+CP and CP+Q53 are the same arm
  ifelse(x %in% c("Q53+CP", "CP+Q53"), "CP+Q53", x)
}

# ---------------------------------------------------------------------------
# 2. Load and subset
# ---------------------------------------------------------------------------

raw <- read_csv(file.path(proj_dir, "data", "summary_table_ddPCR.csv"),
                show_col_types = FALSE)

ht <- raw %>%
  filter(cell_line == "HT1080") %>%
  mutate(
    arm         = canon_arm(name),
    round       = canon_round(Version),
    assay       = canon_assay(TargetName, RefGene),
    crispr_conc = factor(CRISPR_conc),
    well        = Sample,
    run         = ddPCR_date
  )

# Parental baseline: restrict to the Feb 2025 wells (DNA 2025-02-03) that were
# run alongside the Round 2 samples. The Nov 2025 parental is a different DNA
# prep from the passage experiment and would shift the baseline ~30%.
parental <- ht %>%
  filter(arm == "Parental", DNA_date == "2025-02-03")

# The three groups the figure compares.
plot_dat <- bind_rows(
  parental,
  ht %>% filter(arm == "CP+Q53", round %in% c("Round 1", "Round 2"))
) %>%
  filter(!is.na(assay))

# ---------------------------------------------------------------------------
# 3. Normalise to % LOY
#
#   LOY% = (1 - ratio / parental_mean) * 100,  computed within each assay,
#   because the two chemistries have different baselines (~0.5 vs ~1.0).
# ---------------------------------------------------------------------------

baseline <- parental %>%
  group_by(assay) %>%
  summarise(par_mean = mean(ratio), par_n = n(), .groups = "drop")

plot_dat <- plot_dat %>%
  left_join(baseline, by = "assay") %>%
  mutate(loy_pct = (1 - ratio / par_mean) * 100)

# Group means + n for the bars and the axis annotation.
grp <- plot_dat %>%
  group_by(assay, round) %>%
  summarise(mean_loy = mean(loy_pct), n = n(), .groups = "drop")

# Spread of the parental wells = the assay's practical noise floor. Drawn as a
# band so the reader can judge Round 1 against measurement noise rather than
# against a bare zero line.
noise <- plot_dat %>%
  filter(round == "Parental") %>%
  group_by(assay) %>%
  summarise(lo = min(loy_pct), hi = max(loy_pct), .groups = "drop")

# n differs per assay (AMELY has the extra MARCH_09 wells), so n must be
# annotated per facet — an x-axis label lookup would apply one assay's counts
# to both panels.
n_lab <- grp %>% mutate(lab = paste0("n=", n))

# Put each group's % label clear of its highest point.
lab_pos <- plot_dat %>%
  group_by(assay, round) %>%
  summarise(top = max(loy_pct), .groups = "drop") %>%
  left_join(grp %>% select(assay, round, mean_loy), by = c("assay", "round")) %>%
  mutate(y = pmax(top, mean_loy) + 7)

# ---------------------------------------------------------------------------
# 4. Figure
# ---------------------------------------------------------------------------

pal <- c("Parental" = "#00A087", "Round 1" = "#4DBBD5", "Round 2" = "#E64B35")

p <- ggplot(grp, aes(x = round, y = mean_loy, fill = round)) +
  geom_rect(data = noise, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
            fill = "grey70", alpha = 0.25) +
  geom_col(width = 0.62, alpha = 0.85) +
  geom_point(data = plot_dat,
             aes(x = round, y = loy_pct, shape = crispr_conc),
             inherit.aes = FALSE,
             position = position_jitter(width = 0.13, height = 0, seed = 1),
             size = 2.1, stroke = 0.7, colour = "grey20", fill = "white") +
  geom_text(data = lab_pos, inherit.aes = FALSE,
            aes(x = round, y = y, label = sprintf("%.0f%%", mean_loy)),
            size = 3.5, fontface = "bold", colour = "grey20") +
  geom_text(data = n_lab, inherit.aes = FALSE,
            aes(x = round, y = -22, label = lab),
            size = 3, colour = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ assay) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_shape_manual(values = c("1" = 21, "2" = 24),
                     name = "CRISPR conc.", na.translate = FALSE) +
  scale_y_continuous(limits = c(-25, 105), breaks = seq(-25, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "HT1080 loss of Y after successive CP+Q53 nucleofection rounds",
    subtitle = "LOY normalised to the parental baseline within each assay",
    x        = NULL,
    y        = "% LOY (relative to parental)",
    caption  = paste0(
      "Points are individual wells. Round 1 is a single well (n=1) \u2014 no error bars are shown.\n",
      "Grey band = spread of the parental wells, i.e. the practical noise floor of each assay; ",
      "Round 1 falls within it.\n",
      "Baseline: parental wells from Feb 2025 (DNA 2025-02-03). ",
      "Source: summary_table_ddPCR.csv"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12.5),
    plot.subtitle   = element_text(colour = "grey35", size = 10),
    plot.caption    = element_text(colour = "grey45", size = 7.5, hjust = 0),
    strip.text      = element_text(face = "bold", size = 10.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "HT1080_CPQ53_LOY_rounds.png"), p,
       width = 8, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "HT1080_CPQ53_LOY_rounds.pdf"), p,
       width = 8, height = 5, bg = "white")

# ---------------------------------------------------------------------------
# 5. Traceable source values
# ---------------------------------------------------------------------------

src <- plot_dat %>%
  select(assay, round, arm, well, run, crispr_conc,
         DNA_date, ratio, par_mean, loy_pct) %>%
  arrange(assay, round, well)

src_path <- file.path(fig_dir, "HT1080_CPQ53_LOY_rounds_source.csv")

# The CSV is often open in Excel, which locks it on Windows. Don't let that
# abort the run after the figures have already been written.
tryCatch(
  write_csv(src, src_path),
  error = function(e) {
    alt <- sub("\\.csv$", "_new.csv", src_path)
    write_csv(src, alt)
    warning("Could not overwrite ", basename(src_path),
            " (likely open in Excel). Wrote ", basename(alt), " instead.",
            call. = FALSE)
  }
)

cat("\n--- group means ---\n")
print(as.data.frame(grp %>% select(assay, round, n, mean_loy)))
cat("\n--- parental baselines ---\n")
print(as.data.frame(baseline))
cat("\nWrote figures to:", fig_dir, "\n")
