#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Publication figures: Y+ vs Y- low-pass WGS, per cell line.
#
# Emits five figures for each cell line configured in CELL_LINES below:
#   fig1  genome-wide segmented copy number, Y+ vs Y-
#   fig2  chrY zoom with T2T-CHM13v2.0 ideogram
#   fig3  per-chromosome relative copy number
#   fig4  per-bin difference track, Y- minus Y+
#   fig5  chrY region heatmap (PAR1 / MSY Yp / cen / MSY Yq / Yq12 / PAR2)
#
# Input : ../<LINE>_QDNAseq_cn_df.rds - tidy 100 kb-bin table exported by
#         coverage_report_<LINE>_QDNAseq.Rmd (GC + mappability corrected,
#         median-normalized, CBS-segmented, T2T-CHM13v2.0).
#         ../coverage/<sample>.genome.mosdepth.summary.txt - raw depth.
# Output: ../figures/<LINE>_YposYneg/fig{1..5}_*.{pdf,png} + source stats TSVs
#
# Run:    Rscript figures_YposYneg.R          (from WGS/code/)
#
# Scale choice: Figs 1-2 use log2 ratio (standard CN-track convention). Figs
# 3-5 use LINEAR relative copy number, because chrY in a Y- line has zero
# corrected signal and log2(0) is undefined - QDNAseq floors it at log2(1e-3)
# = -9.97, which is a placeholder, not a measurement. Linear space represents
# zero honestly and makes the Y+ vs Y- difference directly interpretable as a
# fraction of ploidy, which also keeps lines of different ploidy comparable
# (one Y is ~0.5 in pseudodiploid HT1080 but ~0.33 in hypotriploid A549).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(cowplot); library(ragg)
})

# --- Cell lines to build ----------------------------------------------------
# `ctrl`/`loy` are the mosdepth sample names; `ctrl_tag`/`loy_tag` are what
# appears in the figures. The RDS `label` column is "sgCTRL"/"LOY" in both
# files regardless of the underlying sample name.
CELL_LINES <- list(
  list(line = "HT1080",
       rds  = "../HT1080_QDNAseq_cn_df.rds",
       ctrl = "HT1080-sgCTRL-D0", ctrl_tag = "sgCTRL",
       loy  = "HT1080-LOY-D0",    loy_tag  = "LOY"),
  list(line = "A549",
       rds  = "../A549_QDNAseq_cn_df.rds",
       ctrl = "A549-CTRL2",       ctrl_tag = "CTRL2",
       loy  = "A549-LOY",         loy_tag  = "LOY")
)

# --- Constants shared across cell lines --------------------------------------
# Detection limit. At ~0.5x genome depth a 100 kb bin holds ~300 reads at
# ploidy; below log2 = -3 (1/8 of ploidy) lpWGS cannot separate true
# homozygous loss from unmappable sequence, and QDNAseq's zero-floor
# placeholder (log2(1e-3) = -9.97) sits here too. Segments below this are
# drawn pinned to the axis floor rather than plotted as real values.
DETECT_LIMIT <- -3
below_detect <- function(x) !is.na(x) & x < DETECT_LIMIT

# --- chrY landmarks, read from the T2T-CHM13v2.0 reference files in the
# project root rather than hardcoded, so the figures cannot drift from the
# annotation the BAMs were aligned against. All coordinates are Mb.
PAR_BED  <- file.path("..", "chm13v2.0_PAR.bed")               # 3-col, no header
CYTO_BED <- file.path("..", "chm13v2.0_cytobands_allchrs.bed") # 5-col, no header
GENE_BED <- file.path("..", "T2T_genes.gz")                    # bigGenePred + header
stopifnot(file.exists(PAR_BED), file.exists(CYTO_BED), file.exists(GENE_BED))

cytobands_chrY <- read_tsv(CYTO_BED, col_types = "ciicc", progress = FALSE,
                           col_names = c("chrom", "start", "end", "band", "stain")) %>%
  filter(chrom == "chrY") %>%
  transmute(start = start / 1e6, end = end / 1e6, band, stain) %>%
  arrange(start)
stopifnot(nrow(cytobands_chrY) > 0)

CHRY_LEN <- max(cytobands_chrY$end)

par_chrY <- read_tsv(PAR_BED, col_types = "cii", progress = FALSE,
                     col_names = c("chrom", "start", "end")) %>%
  filter(chrom == "chrY") %>% arrange(start)
stopifnot(nrow(par_chrY) == 2)
PAR1 <- c(par_chrY$start[1], par_chrY$end[1]) / 1e6
PAR2 <- c(par_chrY$start[2], par_chrY$end[2]) / 1e6

# Centromere = the span of the acrocentric-stain bands (p11.1 + q11.1).
acen <- cytobands_chrY %>% filter(stain == "acen")
stopifnot(nrow(acen) >= 1)
CEN <- c(min(acen$start), max(acen$end))

# Yq12 heterochromatin. The cytoband file's q12 runs to the end of chrY and so
# contains PAR2; truncate it at PAR2's start so the two regions do not
# double-count the same bins.
q12_band <- cytobands_chrY %>% filter(band == "q12")
stopifnot(nrow(q12_band) == 1)
Q12 <- c(q12_band$start, min(q12_band$end, PAR2[1]))

# Gene ticks for fig 2: union of the PROTEIN-CODING transcripts of each gene.
# Restricting to protein_coding excludes long readthrough/processed isoforms
# (e.g. RPS4Y1-204 extends to 2.61 Mb) that would drag the tick off the locus.
CHRY_GENES <- c("RPS4Y1", "ZFY", "AMELY", "UTY", "KDM5D", "EIF1AY")

gene_coords_chrY <- read_tsv(GENE_BED, col_select = c(1, 2, 3, 13, 18),
                             col_types = cols(.default = col_character()),
                             progress = FALSE) %>%
  rlang::set_names(c("chrom", "start", "end", "gene", "type")) %>%
  filter(chrom == "chrY", type == "protein_coding", gene %in% CHRY_GENES) %>%
  group_by(gene) %>%
  summarise(start_Mb = min(as.numeric(start)) / 1e6,
            end_Mb   = max(as.numeric(end))   / 1e6, .groups = "drop") %>%
  arrange(start_Mb)
stopifnot(setequal(gene_coords_chrY$gene, CHRY_GENES))

cat(sprintf(paste0("chrY reference: %.6f Mb, %d cytobands, ",
                   "PAR1 %.6f-%.6f, PAR2 %.6f-%.6f, cen %.6f-%.6f, ",
                   "Yq12 %.6f-%.6f, %d genes\n"),
            CHRY_LEN, nrow(cytobands_chrY), PAR1[1], PAR1[2], PAR2[1], PAR2[2],
            CEN[1], CEN[2], Q12[1], Q12[2], nrow(gene_coords_chrY)))

BIN_HALF <- 0.05    # half a 100 kb bin, in Mb

CTRL_COL <- "#2166AC"
LOY_COL  <- "#B2182B"

stain_colors <- c(gneg = "white", gpos25 = "grey80", gpos50 = "grey55",
                  gpos75 = "grey30", gpos100 = "black",
                  acen = "#B2182B", gvar = "#4393C3", stalk = "#FED976")

CEN_REGION <- "centromere"
# Every region a bin can be assigned to. The centromere is computed and kept in
# the stats TSV, but is not drawn: it is 3 bins of unmappable sequence, zero in
# both lines, and at 0.5% of chrY it is a hairline in the atlas.
region_all    <- c("PAR1", "MSY Yp", CEN_REGION, "MSY Yq (euchromatic)",
                   "Yq12 heterochromatin", "PAR2")
region_levels <- setdiff(region_all, CEN_REGION)          # drawn regions
MAPPABLE_REGIONS <- c("PAR1", "MSY Yp", "MSY Yq (euchromatic)", "PAR2")
SUMMARY_ROWS  <- c("chrY, mappable regions", "chrY, all bins")
all_rows      <- c(region_levels, SUMMARY_ROWS)

assign_region <- function(m) {
  case_when(m < PAR1[2] ~ "PAR1",
            m < CEN[1]  ~ "MSY Yp",
            m < CEN[2]  ~ "centromere",
            m < Q12[1]  ~ "MSY Yq (euchromatic)",
            m < PAR2[1] ~ "Yq12 heterochromatin",
            TRUE        ~ "PAR2")
}

theme_pub <- function(base_size = 9) {
  theme_cowplot(font_size = base_size, font_family = "sans") +
    theme(
      plot.title    = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle = element_text(size = base_size - 2, colour = "grey30"),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 1),
      legend.title  = element_text(size = base_size),
      legend.text   = element_text(size = base_size - 1),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text    = element_text(size = base_size, face = "bold")
    )
}

mm <- function(x) x / 25.4      # mm -> inches

# Wrap caption text to a character width rather than hand-placing line breaks,
# which silently overflow the panel once a conditional clause is appended.
wrap_sub <- function(txt, width) paste(strwrap(txt, width = width), collapse = "\n")

# Push 1D labels apart until none overlap, staying as close to their true
# positions as possible. Left-to-right pass shoves each label clear of its
# predecessor; right-to-left pass pulls the run back inside [lo, hi] if the
# last label overflowed. Deterministic, unlike a repel solver - the same input
# always gives the same figure.
spread_labels <- function(x, half_w, lo, hi, gap) {
  o <- order(x)
  xs <- x[o]; w <- half_w[o]; n <- length(xs)
  xs[1] <- max(xs[1], lo + w[1])
  if (n > 1) for (i in 2:n) xs[i] <- max(xs[i], xs[i - 1] + w[i - 1] + w[i] + gap)
  xs[n] <- min(xs[n], hi - w[n])
  if (n > 1) for (i in seq(n - 1, 1)) {
    xs[i] <- min(xs[i], xs[i + 1] - w[i + 1] - w[i] - gap)
  }
  out <- numeric(n); out[o] <- xs
  out
}

# Approximate width of one bold character of the gene labels, in Mb, at the
# panel width these figures use. Only affects label spacing, not positions.
GENE_CHAR_MB <- 0.60
GENE_LAB_SIZE <- 2.4

# ===========================================================================
# Build all five figures for one cell line
# ===========================================================================
build_figures <- function(cfg) {

  outdir <- file.path("..", "figures", paste0(cfg$line, "_YposYneg"))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  save_fig <- function(plot, name, width_mm, height_mm) {
    ggsave(file.path(outdir, paste0(name, ".pdf")), plot, device = cairo_pdf,
           width = mm(width_mm), height = mm(height_mm), units = "in")
    agg_png(file.path(outdir, paste0(name, ".png")),
            width = mm(width_mm), height = mm(height_mm),
            units = "in", res = 300, background = "white")
    print(plot); invisible(dev.off())
    cat(sprintf("  %-32s %3g x %3g mm\n", paste0(name, ".{pdf,png}"),
                width_mm, height_mm))
  }

  # Figures carry no title or subtitle: the explanatory text is collected here
  # and written to figure_captions.md so the plots drop straight into a
  # manuscript or slide without duplicated headings, and nothing is lost.
  caps <- list()
  add_cap <- function(file, title, body) {
    caps[[length(caps) + 1]] <<- list(file = file, title = title,
                                      body = wrap_sub(body, 100))
    invisible(NULL)
  }

  lab_levels <- c(sprintf("Y+ (%s)", cfg$ctrl_tag),
                  sprintf("Y− (%s)", cfg$loy_tag))
  lab_map    <- c(sgCTRL = lab_levels[1], LOY = lab_levels[2])
  lab_colors <- setNames(c(CTRL_COL, LOY_COL), lab_levels)

  # --- Load & prepare -------------------------------------------------------
  cn <- readRDS(cfg$rds) %>%
    filter(use, !is.na(log2_seg)) %>%
    mutate(label = factor(lab_map[as.character(label)], levels = lab_levels),
           chromosome = factor(chromosome,
                               levels = paste0("chr", c(1:22, "X", "Y"))),
           clipped = below_detect(log2_seg),
           rel_cn  = 2^log2_seg,
           # Per-bin (unsegmented) corrected coverage. Needed for the region
           # heatmap: PAR1, PAR2 and the centromere are smaller than the local
           # CBS segment, so their segmented value is inherited from a
           # neighbouring region rather than measured.
           clipped_raw = below_detect(log2),
           rel_raw = 2^log2)

  stopifnot(nlevels(droplevels(cn$label)) == 2, !any(is.na(cn$label)))
  cat(sprintf("\n=== %s ===\n%d usable bins x 2 samples (%d segments below detection)\n",
              cfg$line, n_distinct(cn$bin_id), sum(cn$clipped)))

  chrom_axis <- cn %>%
    group_by(chromosome) %>%
    summarise(mid = median(genome_x), xmax = max(genome_x), .groups = "drop") %>%
    mutate(idx  = row_number(),
           xmin = lag(xmax, default = 0),
           fill = ifelse(idx %% 2 == 0, "grey94", "white"),
           short = sub("^chr", "", chromosome))

  chrom_band <- function() {
    list(
      geom_rect(data = chrom_axis,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                    fill = fill), inherit.aes = FALSE),
      scale_fill_identity(),
      scale_x_continuous(breaks = chrom_axis$mid, labels = chrom_axis$short,
                         expand = c(0.005, 0))
    )
  }

  # =========================================================================
  # Figure 1 - genome-wide segmented copy number
  # =========================================================================
  # Floor the panel exactly at the detection limit, so every segment that IS
  # quantifiable fits inside the axis and the line never plunges off-scale.
  YLIM1 <- c(DETECT_LIMIT, 1.5)

  f1_line <- cn %>%
    mutate(seg_y = ifelse(clipped | log2_seg > YLIM1[2], NA_real_, log2_seg))
  f1_pts  <- cn %>% filter(log2 >= YLIM1[1], log2 <= YLIM1[2])
  f1_clip <- cn %>% filter(clipped)
  f1_over <- cn %>% filter(!clipped, log2_seg > YLIM1[2])
  n_over  <- nrow(f1_over)

  # Segments above the axis top get the same treatment as those below it:
  # marked at the edge rather than compressing the whole y-axis and flattening
  # the chrY signal this figure exists to show.
  over_txt <- if (n_over > 0) {
    chrs <- f1_over %>% distinct(chromosome) %>% pull(chromosome) %>%
      as.character() %>% sub("^chr", "", .)
    sprintf(paste0(" Marks on the ceiling are %d off-scale segments ",
                   "(chr%s) reaching log2 %.1f."),
            n_over, paste(chrs, collapse = ", "), max(f1_over$log2_seg))
  } else ""

  fig1 <- ggplot() +
    chrom_band() +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_point(data = f1_pts, aes(genome_x, log2),
               colour = "grey60", size = 0.1, alpha = 0.25) +
    geom_segment(data = f1_clip,
                 aes(x = genome_x - BIN_HALF, xend = genome_x + BIN_HALF,
                     y = YLIM1[1], yend = YLIM1[1], colour = label),
                 linewidth = 1.1, alpha = 0.9) +
    geom_segment(data = f1_over,
                 aes(x = genome_x - BIN_HALF, xend = genome_x + BIN_HALF,
                     y = YLIM1[2], yend = YLIM1[2], colour = label),
                 linewidth = 1.1, alpha = 0.9) +
    geom_line(data = f1_line,
              aes(genome_x, seg_y, group = interaction(label, chromosome),
                  colour = label), linewidth = 0.45, na.rm = TRUE) +
    facet_wrap(~ label, ncol = 1, strip.position = "right") +
    scale_colour_manual(values = lab_colors, guide = "none") +
    coord_cartesian(ylim = YLIM1, expand = FALSE) +
    labs(x = "Chromosome", y = expression(log[2]~"copy-number ratio")) +
    theme_pub() +
    theme(strip.text.y.right = element_text(angle = 0),
          panel.spacing = unit(2, "pt"),
          axis.text.x = element_text(size = 6.5))

  add_cap("fig1_genomewide_cn",
          sprintf("%s low-pass WGS — genome-wide copy number", cfg$line),
          sprintf(paste0(
            "Grey points, 100 kb bins; coloured line, CBS segments. Marks on the ",
            "axis floor are segments below the detection limit (log2 < %g) — ",
            "homozygous deletions and unmappable regions, which lpWGS at this ",
            "depth cannot tell apart.%s"), DETECT_LIMIT, over_txt))

  # =========================================================================
  # Figure 2 - chrY zoom with ideogram
  # =========================================================================
  YLIM2 <- c(-4, 1.4)
  chrY_dat <- cn %>% filter(chromosome == "chrY")

  # Quantify over the mappable male-specific region only: distal to PAR1,
  # proximal to the q12 heterochromatic block.
  msy_stats <- chrY_dat %>%
    filter(midpoint > PAR1[2], midpoint < Q12[1]) %>%
    group_by(label) %>%
    summarise(n_bins = n(), n_below = sum(clipped),
              median_log2 = median(log2_seg),
              median_rel_cn = median(rel_cn), .groups = "drop")

  chrY_panel <- function(lab) {
    d <- chrY_dat %>% filter(label == lab)
    s <- msy_stats %>% filter(label == lab)
    note <- if (s$n_below > s$n_bins / 2)
      sprintf("MSY: %d/%d bins below detection", s$n_below, s$n_bins)
    else
      sprintf("MSY median log2 = %.2f  (%.2f × ploidy)",
              s$median_log2, s$median_rel_cn)

    p <- ggplot() +
      annotate("rect", xmin = Q12[1], xmax = PAR2[1], ymin = -Inf, ymax = Inf,
               fill = "grey85", alpha = 0.55) +
      annotate("rect", xmin = PAR1[1], xmax = PAR1[2], ymin = -Inf, ymax = Inf,
               fill = "#41B6C4", alpha = 0.13) +
      annotate("rect", xmin = PAR2[1], xmax = PAR2[2], ymin = -Inf, ymax = Inf,
               fill = "#41B6C4", alpha = 0.13) +
      geom_vline(data = gene_coords_chrY,
                 aes(xintercept = (start_Mb + end_Mb) / 2),
                 colour = "#1A3E72", linewidth = 0.25, alpha = 0.3) +
      geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
      geom_point(data = filter(d, log2 >= YLIM2[1], log2 <= YLIM2[2]),
                 aes(midpoint, log2), colour = "grey55", size = 0.4,
                 alpha = 0.55) +
      geom_line(data = mutate(d, seg_y = ifelse(clipped, NA_real_, log2_seg)),
                aes(midpoint, seg_y), colour = lab_colors[[lab]],
                linewidth = 0.7, na.rm = TRUE) +
      annotate("text", x = 27.2, y = YLIM2[2] - 0.3, label = note,
               hjust = 1, size = 2.5, colour = "grey25", fontface = "italic") +
      coord_cartesian(xlim = c(0, CHRY_LEN), ylim = YLIM2, expand = FALSE) +
      labs(x = NULL, y = expression(log[2]~"ratio"), title = lab) +
      theme_pub() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            plot.margin = margin(t = 4, r = 8, b = 0, l = 8))

    fl <- d %>% filter(clipped)
    if (nrow(fl))
      p <- p + geom_segment(data = fl,
                            aes(x = midpoint - BIN_HALF,
                                xend = midpoint + BIN_HALF,
                                y = YLIM2[1], yend = YLIM2[1]),
                            colour = lab_colors[[lab]], linewidth = 1.4)
    p
  }

  chrY_ideogram <- function() {
    # Single row of gene labels, spread apart so none overlap, each joined to
    # its true position by an elbow leader. The old two-row stagger placed
    # labels correctly but ran the leader for the lower row straight through
    # the text of the upper row (RPS4Y1/ZFY, KDM5D/EIF1AY are ~0.1 and ~0.8 Mb
    # apart, far closer than a label is wide).
    gene_pos <- gene_coords_chrY %>%
      mutate(gene_x = (start_Mb + end_Mb) / 2,
             half_w = nchar(gene) * GENE_CHAR_MB / 2)
    gene_pos$label_x <- spread_labels(gene_pos$gene_x, gene_pos$half_w,
                                      0, CHRY_LEN, gap = 0.35)
    cen_mid <- mean(CEN)

    ggplot() +
      annotate("rect", xmin = PAR1[1], xmax = PAR1[2], ymin = 1.08, ymax = 1.22,
               fill = "#41B6C4", alpha = 0.9, colour = "grey20",
               linewidth = 0.2) +
      annotate("text", x = PAR1[2] + 0.4, y = 1.15, label = "PAR1", hjust = 0,
               size = 2.5, fontface = "bold", colour = "#0C2C84") +
      annotate("rect", xmin = PAR2[1], xmax = PAR2[2], ymin = 1.08, ymax = 1.22,
               fill = "#41B6C4", alpha = 0.9, colour = "grey20",
               linewidth = 0.2) +
      annotate("text", x = PAR2[1] - 0.4, y = 1.15, label = "PAR2", hjust = 1,
               size = 2.5, fontface = "bold", colour = "#0C2C84") +
      geom_rect(data = cytobands_chrY,
                aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = stain),
                colour = "grey30", linewidth = 0.2) +
      scale_fill_manual(values = stain_colors, guide = "none") +
      annotate("text", x = mean(c(Q12[1], PAR2[1])), y = 0.5,
               label = "q12 heterochromatin — not mappable",
               size = 2.4, colour = "white", fontface = "bold") +
      # cen marker sits above the bands, clear of the gene-label lane below.
      annotate("point", x = cen_mid, y = 1.10, shape = 25, size = 1.5,
               colour = "grey15", fill = "grey15") +
      annotate("text", x = cen_mid, y = 1.36, label = "cen",
               size = 2.4, fontface = "italic") +
      # Elbow leader: straight down out of the band, then across to the label.
      geom_segment(data = gene_pos,
                   aes(x = gene_x, xend = gene_x, y = -0.02, yend = -0.30),
                   colour = "#1A3E72", linewidth = 0.35) +
      geom_segment(data = gene_pos,
                   aes(x = gene_x, xend = label_x, y = -0.30, yend = -0.46),
                   colour = "#1A3E72", linewidth = 0.35) +
      geom_text(data = gene_pos, aes(label_x, -0.50, label = gene),
                colour = "#1A3E72", size = GENE_LAB_SIZE, vjust = 1,
                fontface = "bold") +
      coord_cartesian(xlim = c(0, CHRY_LEN), ylim = c(-0.85, 1.55),
                      expand = FALSE) +
      labs(x = "chrY position (Mb)") +
      theme_pub() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.line.y = element_blank(), axis.title.y = element_blank(),
            plot.margin = margin(t = 0, r = 8, b = 4, l = 8))
  }

  fig2 <- plot_grid(chrY_panel(lab_levels[1]), chrY_panel(lab_levels[2]),
                    chrY_ideogram(), ncol = 1, align = "v", axis = "lr",
                    rel_heights = c(2.3, 2.3, 1.9))

  # =========================================================================
  # Figure 3 - per-chromosome relative copy number
  # =========================================================================
  depth_ratio <- function(sample) {
    d <- read_tsv(file.path("..", "coverage",
                            paste0(sample, ".genome.mosdepth.summary.txt")),
                  show_col_types = FALSE) %>%
      filter(chrom %in% paste0("chr", c(1:22, "X", "Y")))
    auto <- d %>% filter(chrom %in% paste0("chr", 1:22))
    tibble(sample = sample,
           autosome_depth = sum(auto$bases) / sum(auto$length),
           chrY_depth = d %>% filter(chrom == "chrY") %>% with(bases / length)) %>%
      mutate(ratio = chrY_depth / autosome_depth)
  }

  ratios <- bind_rows(depth_ratio(cfg$ctrl), depth_ratio(cfg$loy)) %>%
    mutate(label = factor(lab_levels, levels = lab_levels))

  chrom_med <- cn %>%
    group_by(label, chromosome) %>%
    summarise(n_bins = n(), n_below_detect = sum(clipped),
              median_log2_seg = median(log2_seg),
              median_rel_cn = median(rel_cn), .groups = "drop") %>%
    mutate(short = sub("^chr", "", chromosome))

  # A bar at ~0.001 is invisible and reads as missing data, so label it as the
  # measured zero it is.
  zero_bars <- chrom_med %>% filter(median_rel_cn < 0.02)
  y3_top <- max(1.5, max(chrom_med$median_rel_cn) * 1.06)

  # A whole-chromosome median for chrY is dominated by its unmappable
  # centromere + Yq12 bins (>55% of chrY). Where that drags the Y+ bar well
  # below the line's actual MSY level, say so rather than letting the bar
  # imply the control has no Y.
  chrY_ctrl <- chrom_med %>%
    filter(chromosome == "chrY", label == lab_levels[1]) %>% pull(median_rel_cn)
  msy_ctrl  <- msy_stats$median_rel_cn[msy_stats$label == lab_levels[1]]
  chrY_note <- if (chrY_ctrl < msy_ctrl / 2)
    paste0("\n", paste(strwrap(sprintf(paste0(
      "chrY's whole-chromosome median (%.3f) is dragged down by unmappable ",
      "centromere/Yq12 bins; over its mappable regions the Y+ line reads ",
      "%.3f — see fig 5."), chrY_ctrl, msy_ctrl), width = 92),
      collapse = "\n")) else ""

  fig3 <- ggplot(chrom_med, aes(short, median_rel_cn, fill = label)) +
    geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.3,
               linetype = "22") +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             colour = "grey25", linewidth = 0.15) +
    geom_text(data = zero_bars, aes(short, y3_top * 0.04, label = "0"),
              position = position_dodge(width = 0.75),
              size = 2.6, fontface = "bold", colour = LOY_COL) +
    scale_fill_manual(values = lab_colors, name = NULL) +
    scale_x_discrete(limits = sub("^chr", "", levels(cn$chromosome))) +
    scale_y_continuous(breaks = seq(0, 2, 0.25), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, y3_top)) +
    labs(x = "Chromosome", y = "median relative copy number") +
    theme_pub() +
    theme(legend.position = "top", legend.justification = "left",
          axis.text.x = element_text(size = 6.5))

  add_cap("fig3_per_chromosome",
          sprintf("%s — per-chromosome copy number, Y+ vs Y−", cfg$line),
          paste0("Raw chrY:autosome sequencing depth — ",
                 paste(sprintf("%s %.3f", ratios$label, ratios$ratio),
                       collapse = " | "),
                 sprintf(" (%.1f-fold).", ratios$ratio[1] / ratios$ratio[2]),
                 gsub("\n", " ", chrY_note)))

  # =========================================================================
  # Figure 4 - per-bin difference track, Y- minus Y+
  # =========================================================================
  diff_all <- cn %>%
    select(bin_id, chromosome, genome_x, label, rel_cn, clipped) %>%
    mutate(k = ifelse(label == lab_levels[1], "ypos", "yneg")) %>%
    select(-label) %>%
    pivot_wider(names_from = k, values_from = c(rel_cn, clipped)) %>%
    rename(ypos = rel_cn_ypos, yneg = rel_cn_yneg,
           c_pos = clipped_ypos, c_neg = clipped_yneg) %>%
    filter(!is.na(ypos), !is.na(yneg)) %>%
    mutate(delta = yneg - ypos,
           is_chrY = chromosome == "chrY",
           # Bins below detection in BOTH lines (unmappable sequence, or a
           # deletion shared by both) carry no information about the
           # difference: their delta is a forced ~0 that would dilute every
           # summary.
           uninformative = c_pos & c_neg)

  n_uninf <- sum(diff_all$uninformative)
  diff_df <- diff_all %>% filter(!uninformative)

  noise <- diff_df %>% filter(!chromosome %in% c("chrX", "chrY")) %>%
    pull(delta) %>% mad(na.rm = TRUE)

  chrY_delta <- median(diff_df$delta[diff_df$is_chrY])
  ypos_chrY  <- msy_stats$median_rel_cn[msy_stats$label == lab_levels[1]]

  YLIM4 <- c(-1.1, 1.1)
  f4 <- diff_df %>%
    mutate(off = delta < YLIM4[1] | delta > YLIM4[2],
           y   = pmax(pmin(delta, YLIM4[2]), YLIM4[1]))
  n_off <- sum(f4$off)

  fig4 <- ggplot(filter(f4, !off)) +
    chrom_band() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -noise, ymax = noise,
             fill = "grey65", alpha = 0.5) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.3) +
    geom_hline(yintercept = -ypos_chrY, colour = LOY_COL, linewidth = 0.3,
               linetype = "22") +
    geom_point(aes(genome_x, y, colour = is_chrY, size = is_chrY,
                   alpha = is_chrY)) +
    geom_point(data = filter(f4, off),
               aes(genome_x, y, shape = delta > 0), colour = "grey20",
               size = 0.7, stroke = 0.2, fill = "grey20") +
    scale_shape_manual(values = c(`FALSE` = 25, `TRUE` = 24), guide = "none") +
    scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = LOY_COL),
                        guide = "none") +
    scale_size_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.6), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.55, `TRUE` = 0.9), guide = "none") +
    coord_cartesian(ylim = YLIM4, expand = FALSE) +
    labs(x = "Chromosome", y = "Δ relative copy number  (Y− minus Y+)") +
    theme_pub() +
    theme(axis.text.x = element_text(size = 6.5))

  add_cap("fig4_diff_track",
          sprintf("%s — copy-number difference, Y− minus Y+", cfg$line),
          sprintf(paste0(
            "chrY (red) median Δ = %.2f, against the Y+ chrY level of %.2f ",
            "(dashed line). Grey band = ±1 MAD of the autosomal difference ",
            "(%.3f); autosomal segments outside it are clonal drift between the ",
            "two lines, not Y-linked. %s bins below detection in both lines ",
            "excluded; %d off-scale bins shown as triangles."),
            chrY_delta, ypos_chrY, noise, format(n_uninf, big.mark = ","),
            n_off))

  # =========================================================================
  # Figure 5 - chrY region heatmap
  # =========================================================================
  chrY_reg <- chrY_dat %>% mutate(region = assign_region(midpoint))

  # Medians of per-bin coverage, NOT of CBS segment values - see rel_raw above.
  # median_seg is carried through to the TSV for reference only.
  summarise_region <- function(d) {
    summarise(d, n_bins = n(), n_below = sum(clipped_raw),
              median_rel_cn = median(rel_raw), median_log2 = median(log2),
              median_seg = median(rel_cn), .groups = "drop")
  }

  reg_stats <- bind_rows(
    chrY_reg %>% group_by(region, label) %>% summarise_region(),
    chrY_reg %>% filter(region %in% MAPPABLE_REGIONS) %>% group_by(label) %>%
      summarise_region() %>% mutate(region = SUMMARY_ROWS[1]),
    chrY_reg %>% group_by(label) %>% summarise_region() %>%
      mutate(region = SUMMARY_ROWS[2])
  ) %>%
    # mostly_below still drives which caveats appear in the caption text; it is
    # no longer marked on the cells themselves. Per-region below-detection
    # counts are in the n_below column of figure_chrY_region_stats.tsv.
    mutate(mostly_below = n_below > n_bins / 2,
           cell = sprintf("%.3f", median_rel_cn))

  # Drawn subset: everything except the centromere (kept in reg_stats/TSV).
  reg_plot <- reg_stats %>%
    filter(region != CEN_REGION) %>%
    mutate(row = factor(region, levels = rev(all_rows)))

  row_labs <- reg_plot %>% distinct(region, n_bins) %>%
    mutate(lab = sprintf("%s  (%d)", region, n_bins)) %>%
    with(setNames(lab, region))

  reg_delta <- reg_plot %>%
    mutate(k = ifelse(label == lab_levels[1], "ypos", "yneg")) %>%
    select(row, k, median_rel_cn) %>%
    pivot_wider(names_from = k, values_from = median_rel_cn) %>%
    mutate(delta = yneg - ypos, x = "Δ", cell = sprintf("%+.3f", delta))

  CN_MAX  <- max(0.6, max(reg_plot$median_rel_cn) * 1.05)
  DEL_LIM <- max(0.55, max(abs(reg_delta$delta)) * 1.05)
  seq_pal <- c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B")
  div_pal <- c("#67001F", "#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7",
               "#D1E5F0", "#67A9CF", "#2166AC")

  tile_theme <- function() {
    theme_pub(8) +
      theme(axis.line = element_blank(), axis.ticks = element_blank(),
            axis.title = element_blank(), legend.position = "bottom",
            legend.key.height = unit(3.5, "pt"),
            legend.key.width  = unit(26, "pt"),
            legend.title = element_text(size = 6.5, vjust = 1),
            legend.text  = element_text(size = 6),
            plot.margin = margin(t = 2, r = 2, b = 2, l = 2))
  }

  p_main <- ggplot(reg_plot, aes(label, row, fill = median_rel_cn)) +
    geom_tile(colour = "white", linewidth = 1.1) +
    geom_text(aes(label = cell, colour = median_rel_cn > CN_MAX * 0.55),
              size = 2.7, fontface = "bold") +
    geom_hline(yintercept = length(SUMMARY_ROWS) + 0.5, colour = "grey35",
               linewidth = 0.4) +
    scale_fill_gradientn(colours = seq_pal, limits = c(0, CN_MAX),
                         oob = scales::squish, name = "relative copy number") +
    scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"),
                        guide = "none") +
    scale_y_discrete(labels = row_labs) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    tile_theme() +
    theme(axis.text.x = element_text(size = 8, face = "bold"))

  p_delta <- ggplot(reg_delta, aes(x, row, fill = delta)) +
    geom_tile(colour = "white", linewidth = 1.1) +
    geom_text(aes(label = cell, colour = abs(delta) > DEL_LIM * 0.6),
              size = 2.7, fontface = "bold") +
    geom_hline(yintercept = length(SUMMARY_ROWS) + 0.5, colour = "grey35",
               linewidth = 0.4) +
    scale_fill_gradientn(colours = div_pal, limits = c(-DEL_LIM, DEL_LIM),
                         oob = scales::squish,
                         name = "Δ  (Y− − Y+)") +
    scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"),
                        guide = "none") +
    scale_y_discrete(labels = NULL, breaks = NULL) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    tile_theme() +
    theme(axis.text.x = element_text(size = 8, face = "bold"))

  # Footnotes are derived from this line's own numbers, not hardcoded: which
  # caveats actually apply differs between cell lines.
  get_reg <- function(reg, lab) reg_stats %>% filter(region == reg, label == lab)
  yp_neg  <- get_reg("MSY Yp", lab_levels[2])
  yq_neg  <- get_reg("MSY Yq (euchromatic)", lab_levels[2])
  q12_pos <- get_reg("Yq12 heterochromatin", lab_levels[1])
  par2_n  <- get_reg("PAR2", lab_levels[1])$n_bins
  par2_pos_below <- get_reg("PAR2", lab_levels[1])$mostly_below

  all_pos <- get_reg(SUMMARY_ROWS[2], lab_levels[1])
  q12_n   <- q12_pos$n_bins
  chrY_n  <- all_pos$n_bins

  foot <- c(
    paste0("Values are medians of per-bin corrected coverage, not of CBS ",
           "segment values: PAR1, PAR2 and the centromere are smaller than ",
           "the local segment and would otherwise inherit a neighbouring ",
           "region's value. Segment medians are in the accompanying TSV."),
    sprintf(paste0("PAR1 and PAR2 are sequence-identical on X and Y, so a Y− ",
                   "cell keeps reads there whatever its Y status; PAR2 rests ",
                   "on %d bins%s."),
            par2_n,
            if (par2_pos_below)
              paste0(" and is itself below detection in the Y+ line, so that ",
                     "row is uninformative for this line")
            else ""),
    sprintf(paste0("Yq12 (%d bins) carries no uniquely mappable sequence in ",
                   "either line, so its near-zero value is not evidence of ",
                   "loss. The centromere (%d bins, likewise unmappable and zero ",
                   "in both lines) is not drawn; it remains in the stats TSV."),
            q12_n, get_reg(CEN_REGION, lab_levels[1])$n_bins),
    # Trigger on the gap to the mappable-region value, not on a bare majority
    # of below-detection bins: A549 sits just under that majority yet its
    # all-bins value is still less than half its mappable value.
    if (all_pos$median_rel_cn <
        get_reg(SUMMARY_ROWS[1], lab_levels[1])$median_rel_cn / 1.5)
      sprintf(paste0("Because %d of %d chrY bins are unmappable centromere or ",
                     "Yq12, the 'all bins' row sits near zero in both lines and ",
                     "is not a usable chrY summary — read the 'mappable regions' ",
                     "row instead."),
              q12_n + get_reg("centromere", lab_levels[1])$n_bins, chrY_n)
    else NULL,
    if (yp_neg$mostly_below && yq_neg$mostly_below)
      sprintf(paste0("Both MSY arms are below detection in Y− (%d/%d Yp and ",
                     "%d/%d Yq bins): the male-specific region is entirely absent."),
              yp_neg$n_below, yp_neg$n_bins, yq_neg$n_below, yq_neg$n_bins)
    else if (yp_neg$median_rel_cn > 3 * yq_neg$median_rel_cn)
      sprintf(paste0("MSY Yp retains ~%.2f× in Y−, well above Yq: X-homologous ",
                     "sequence on Yp and residual Y+ cells would both produce ",
                     "this, and lpWGS alone does not separate them."),
              yp_neg$median_rel_cn)
    else NULL
  )
  # Wrap to the figure width rather than relying on hand-placed line breaks,
  # which silently ran off the edge when a conditional clause was added.
  foot <- vapply(foot, function(s) paste(strwrap(s, width = 116), collapse = "\n"),
                 character(1), USE.NAMES = FALSE)
  n_foot <- sum(lengths(strsplit(foot, "\n")))

  fig5 <- plot_grid(p_main, p_delta, nrow = 1, rel_widths = c(2.5, 1.25),
                    align = "h", axis = "tb")

  add_cap("fig5_chrY_region_heatmap",
          sprintf("%s — chrY by region, Y+ (%s) vs Y− (%s)",
                  cfg$line, cfg$ctrl_tag, cfg$loy_tag),
          paste0("Median relative copy number per 100 kb bin, by chrY region ",
                 "(T2T-CHM13v2.0); bin count in each row label. Cells reading ",
                 "~0.001 are measured zeros, not quantities; the per-region ",
                 "count of bins below the detection limit is the n_below column ",
                 "of figure_chrY_region_stats.tsv. ",
                 paste(gsub("\n", " ", foot), collapse = " ")))

  # =========================================================================
  # Figure 6 - chrY atlas: bin track + region strips + ideogram, all sharing
  # one genomic axis, so region values sit above the data they summarise and
  # each region occupies its true width (Yq12 is 55% of chrY; the centromere
  # and PAR2 are 0.5% each).
  # =========================================================================
  YLIM6 <- c(0, 1.3)

  # Region extents, identical to the assign_region() cut points.
  # The centromere gap (CEN[1]-CEN[2]) is deliberately left undrawn.
  region_bounds <- tibble(
    region = region_levels,
    start = c(PAR1[1], PAR1[2], CEN[2], Q12[1], PAR2[1]),
    end   = c(PAR1[2], CEN[1],  Q12[1], Q12[2], PAR2[2])
  ) %>%
    mutate(mid = (start + end) / 2,
           width_frac = (end - start) / CHRY_LEN,
           short = c("PAR1", "MSY Yp", "MSY Yq", "Yq12", "PAR2"),
           # Keep labels inside the panel at both ends of the chromosome.
           hj = case_when(mid < 0.03 * CHRY_LEN ~ 0,
                          mid > 0.92 * CHRY_LEN ~ 1,
                          TRUE ~ 0.5))

  reg_wide <- reg_plot %>%
    filter(region %in% region_levels) %>%
    mutate(k = ifelse(label == lab_levels[1], "ypos", "yneg")) %>%
    select(region, k, median_rel_cn, mostly_below) %>%
    pivot_wider(names_from = k,
                values_from = c(median_rel_cn, mostly_below)) %>%
    left_join(region_bounds, by = "region")

  # Vertical guides at region boundaries tie the three panels together.
  guides6 <- sort(unique(c(region_bounds$start, region_bounds$end)))
  guide_layer <- function() geom_vline(xintercept = guides6, colour = "grey80",
                                       linewidth = 0.2, linetype = "22")

  n_over6 <- chrY_dat %>% filter(rel_raw > YLIM6[2]) %>% nrow()

  bin_panel <- function(lab) {
    d   <- chrY_dat %>% filter(label == lab)
    med <- reg_plot %>% filter(region %in% region_levels, label == lab) %>%
      left_join(region_bounds, by = "region")
    ggplot() +
      annotate("rect", xmin = Q12[1], xmax = Q12[2], ymin = -Inf, ymax = Inf,
               fill = "grey85", alpha = 0.55) +
      annotate("rect", xmin = PAR1[1], xmax = PAR1[2], ymin = -Inf, ymax = Inf,
               fill = "#41B6C4", alpha = 0.13) +
      annotate("rect", xmin = PAR2[1], xmax = PAR2[2], ymin = -Inf, ymax = Inf,
               fill = "#41B6C4", alpha = 0.13) +
      guide_layer() +
      geom_point(data = filter(d, rel_raw <= YLIM6[2]),
                 aes(midpoint, rel_raw), colour = "grey55", size = 0.45,
                 alpha = 0.6) +
      geom_point(data = filter(d, rel_raw > YLIM6[2]),
                 aes(midpoint, YLIM6[2]), shape = 24, size = 0.9, stroke = 0.2,
                 colour = "grey25", fill = "grey25") +
      geom_segment(data = med,
                   aes(x = start, xend = end,
                       y = median_rel_cn, yend = median_rel_cn),
                   colour = lab_colors[[lab]], linewidth = 1) +
      coord_cartesian(xlim = c(0, CHRY_LEN), ylim = YLIM6, expand = FALSE) +
      labs(x = NULL, y = "rel. CN", title = lab) +
      theme_pub() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            plot.title = element_text(size = 9, face = "bold"),
            plot.margin = margin(t = 3, r = 8, b = 0, l = 8))
  }

  # Strip panel geometry (arbitrary units; the y axis is categorical in spirit)
  YP_LO <- 2.15; YP_HI <- 2.95
  YN_LO <- 1.25; YN_HI <- 2.05
  # Regions thinner than this get callout labels instead of an inline value.
  # 0.06 puts PAR1 (3.9% of chrY) in the callout group with the centromere and
  # PAR2: its rectangle is ~7 mm, narrower than the number it would hold.
  NARROW <- 0.06

  strip_dat <- reg_wide %>%
    mutate(narrow = width_frac < NARROW)

  strip_rects <- bind_rows(
    strip_dat %>% transmute(start, end, mid, narrow, hj,
                            v = median_rel_cn_ypos,
                            below = mostly_below_ypos,
                            lo = YP_LO, hi = YP_HI, row = "ypos"),
    strip_dat %>% transmute(start, end, mid, narrow, hj,
                            v = median_rel_cn_yneg,
                            below = mostly_below_yneg,
                            lo = YN_LO, hi = YN_HI, row = "yneg")
  ) %>%
    mutate(txt = sprintf("%.3f", v),
           ctr = (lo + hi) / 2,
           row_col = ifelse(row == "ypos", CTRL_COL, LOY_COL))

  p_strip <- ggplot() +
    guide_layer() +
    geom_rect(data = strip_rects,
              aes(xmin = start, xmax = end, ymin = lo, ymax = hi, fill = v),
              colour = "white", linewidth = 0.4) +
    scale_fill_gradientn(colours = seq_pal, limits = c(0, CN_MAX),
                         oob = scales::squish, guide = "none") +
    # Region names along the top
    geom_text(data = region_bounds,
              aes(mid, 3.18, label = short, hjust = hj),
              size = 2.5, fontface = "bold", colour = "grey20") +
    # Values inside the strip where the region is wide enough...
    geom_text(data = filter(strip_rects, !narrow),
              aes(mid, ctr, label = txt, colour = v > CN_MAX * 0.55),
              size = 2.5, fontface = "bold") +
    scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"),
                        guide = "none") +
    # ...and as callouts below it where it is not (centromere, PAR2).
    geom_segment(data = filter(strip_rects, narrow, row == "yneg"),
                 aes(x = mid, xend = mid, y = YN_LO, yend = 1.02),
                 colour = "grey45", linewidth = 0.25) +
    geom_text(data = filter(strip_rects, narrow, row == "ypos"),
              aes(mid, 0.80, label = txt, hjust = hj),
              size = 2.5, fontface = "bold", colour = CTRL_COL) +
    geom_text(data = filter(strip_rects, narrow, row == "yneg"),
              aes(mid, 0.36, label = txt, hjust = hj),
              size = 2.5, fontface = "bold", colour = LOY_COL) +
    scale_y_continuous(breaks = c((YP_LO + YP_HI) / 2, (YN_LO + YN_HI) / 2),
                       labels = c("Y+", "Y−")) +
    coord_cartesian(xlim = c(0, CHRY_LEN), ylim = c(0.1, 3.5), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(size = 8, face = "bold"),
          plot.margin = margin(t = 0, r = 8, b = 0, l = 8))

  f6_key <- ggdraw() +
    draw_label(sprintf(paste0("▲ bins above %.1f — repeat pile-ups in ",
                              "unmappable sequence, not copy number"), YLIM6[2]),
               size = 6.0, colour = "grey35", x = 0.03, hjust = 0, y = 0.6)

  fig6 <- plot_grid(
    plot_grid(bin_panel(lab_levels[1]), bin_panel(lab_levels[2]),
              p_strip, chrY_ideogram(),
              ncol = 1, align = "v", axis = "lr",
              rel_heights = c(2.1, 2.1, 1.5, 1.55)),
    f6_key,
    ncol = 1, rel_heights = c(1, 0.045)
  )

  add_cap("fig6_chrY_atlas",
          sprintf("%s — chrY atlas: per-bin coverage, region medians, ideogram",
                  cfg$line),
          sprintf(paste0(
            "Grey points are 100 kb bins (relative copy number); the coloured ",
            "bar in each track is that region's median. Strips below repeat ",
            "those medians on the fig-5 colour scale, each region drawn at its ",
            "true width; the centromere (3 unmappable bins, 0.5%% of chrY) is ",
            "left as a gap between MSY Yp and MSY Yq. ▲ marks %d bins above ",
            "%.1f — repeat pile-ups in unmappable sequence, not copy number. ",
            "Strip values reading ~0.001 are measured zeros; per-region counts ",
            "of bins below the detection limit are in the n_below column of ",
            "figure_chrY_region_stats.tsv."), n_over6, YLIM6[2]))

  add_cap("fig2_chrY_ideogram",
          sprintf("%s — chrY log2 copy number with T2T-CHM13v2.0 ideogram",
                  cfg$line),
          sprintf(paste0(
            "Grey points, 100 kb bins; coloured line, CBS segments; marks on ",
            "the axis floor are segments below the detection limit (log2 < %g). ",
            "Cyan bands are PAR1/PAR2, grey is the Yq12 heterochromatic block. ",
            "Quoted medians cover the mappable male-specific region only ",
            "(%.2f–%.2f Mb)."), DETECT_LIMIT, PAR1[2], Q12[1]))

  # =========================================================================
  # Write
  # =========================================================================
  cat("Writing to", normalizePath(outdir, mustWork = FALSE), "\n")
  save_fig(fig1, "fig1_genomewide_cn",       180, 105)
  save_fig(fig2, "fig2_chrY_ideogram",       180, 145)
  save_fig(fig3, "fig3_per_chromosome",      130,  95)
  save_fig(fig4, "fig4_diff_track",          180,  80)
  save_fig(fig5, "fig5_chrY_region_heatmap", 165,  72)
  save_fig(fig6, "fig6_chrY_atlas",          180, 138)

  cap_md <- c(sprintf("# %s — Y+ (%s) vs Y− (%s): figure captions",
                      cfg$line, cfg$ctrl_tag, cfg$loy_tag),
              "",
              paste0("Generated by `figures_YposYneg.R`. Figures carry no ",
                     "title or subtitle; the text below belongs with them."),
              "",
              unlist(lapply(caps[order(vapply(caps, `[[`, "", "file"))],
                            function(x) c(sprintf("## %s", x$file), "",
                                          sprintf("**%s**", x$title), "",
                                          x$body, ""))))
  writeLines(cap_md, file.path(outdir, "figure_captions.md"))

  write_tsv(chrom_med %>% select(-short) %>% arrange(label, chromosome),
            file.path(outdir, "figure_source_stats.tsv"))
  write_tsv(ratios, file.path(outdir, "figure_raw_depth_ratios.tsv"))
  write_tsv(msy_stats, file.path(outdir, "figure_chrY_MSY_stats.tsv"))
  write_tsv(reg_stats %>%
              select(region, label, n_bins, n_below, median_rel_cn, median_log2,
                     median_seg) %>%
              left_join(reg_delta %>% mutate(region = as.character(row)) %>%
                          select(region, delta), by = "region") %>%
              arrange(match(region, all_rows), label),
            file.path(outdir, "figure_chrY_region_stats.tsv"))

  cat("\n--- Raw mosdepth depth ---\n"); print(as.data.frame(ratios))
  cat("\n--- chrY MSY (2.46-27.45 Mb) ---\n"); print(as.data.frame(msy_stats))
  cat(sprintf("\nAutosomal difference MAD: %.4f | chrY median delta: %.4f ",
              noise, chrY_delta))
  cat(sprintf("(Y+ chrY %.4f)\nShared below-detection bins excluded: %d | ",
              ypos_chrY, n_uninf))
  cat(sprintf("off-scale in fig4: %d\n", n_off))
  cat("\n--- chrY by region ---\n")
  print(reg_stats %>% select(region, label, n_bins, n_below, median_rel_cn) %>%
          arrange(match(region, all_rows), label) %>% as.data.frame(),
        digits = 4)

  invisible(list(line = cfg$line, regions = reg_stats, ratios = ratios))
}

invisible(lapply(CELL_LINES, build_figures))
cat("\nDone.\n")
