# Combined screen pathway-enrichment supplement (author requests, 22 October 2026):
# a, Hallmark GSEA (50 sets) and b, expanded Hallmark/Reactome/GO:BP/GO:CC GSEA
# (5,880 sets, top 12 per direction) for the primary HT1080 endpoint drugZ ranking;
# c, Metascape enrichment of the 48 drugZ depletion candidates and d, of the 45
# IsogenicZ depletion candidates (both at screen FDR <0.25).
# Layout: a/b in the left column, c/d in the right column with the shared NES legend below.
# Reads saved plotted values and saved Metascape representatives; nothing is refitted.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

thesis <- "A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS"
skill_dir <- "A:/Dropbox/LOY_paralogs/.claude/skills/thesis-figure-style"
fig_dir <- file.path(thesis, "figures/chapter4")
src <- list(
  hallmark = file.path(fig_dir, "source_data/ht1080_hallmark_screen_gsea_20261021/figure_values.csv"),
  expanded = file.path(fig_dir, "source_data/screen_expanded_gsea_20261021/figure_values.csv"),
  metascape_drugz = file.path(fig_dir, "source_data/ht1080_drugz_metascape_48/panel_a_cluster_representatives.csv"),
  metascape_isogenicz = file.path(fig_dir, "source_data/ht1080_isogenicz_metascape_45_fdr025/supported_cluster_representatives.csv")
)
out_data <- file.path(fig_dir, "source_data/screen_gsea_combined_20261022")
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)

source(file.path(skill_dir, "scripts/thesis_style.R"))
register_thesis_fonts(skill_dir)

sup_num <- function(e) {
  sup <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴", "5" = "⁵",
           "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹", "-" = "⁻")
  paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")
}
fmt_q <- function(q) vapply(q, function(x) {
  if (x >= 0.001) return(formatC(x, format = "fg", digits = 2))
  e <- floor(log10(x)); paste0(formatC(x / 10^e, format = "f", digits = 1), " × 10", sup_num(e))
}, "")
wrap <- function(x, width) vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), "")
n_lines <- function(x) lengths(regmatches(x, gregexpr("\n", x))) + 1L

panel_theme <- function() theme_thesis(base_size = 6.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.subtitle = element_text(size = 6.8, colour = "#555555", margin = margin(0, 0, 4, 0)),
        plot.title.position = "plot",
        plot.tag = element_text(size = 11, face = "bold", family = "Open Sans"),
        legend.title = element_text(face = "plain", size = 6.8),
        legend.text = element_text(size = 6.5),
        plot.margin = margin(3, 4, 3, 2, "pt"))

# ---- a, b: GSEA ----------------------------------------------------------------
hm <- fread(src$hallmark)
ex <- fread(src$expanded)
stopifnot(nrow(hm) == 6, all(hm$sig), nrow(ex) == 24, all(ex$sig))

set_name <- function(p) gsub("_", " ", sub("^(HALLMARK|REACTOME|GOBP|GOCC)_", "", p))
hm[, `:=`(label = set_name(pathway), FDR = FDR_BH_within_ranking, panel = "a")]
# Expanded set names (with collection) wrap at 46 characters to keep the left column narrow.
ex[, `:=`(label = wrap(paste0(set_name(pathway), " [", collection, "]"), 46),
          FDR = FDR_global, panel = "b")]

lim <- ceiling(max(abs(c(hm$NES, ex$NES))) * 2) / 2
fill_scale <- scale_fill_gradient2(
  low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
  limits = c(-lim, lim), breaks = c(-3, 0, 3),
  labels = c("\u22123\nresistor", "0", "3\nsensitizer"), name = "NES (a, b)",
  guide = guide_colourbar(direction = "horizontal", title.position = "top",
                          barwidth = unit(80, "pt"), barheight = unit(6, "pt"), order = 1))
size_scale <- scale_size_continuous(range = c(1.4, 4.4), limits = c(1, lim),
                                    breaks = c(1.5, 2, 2.5, 3, 3.5), name = "|NES|",
                                    guide = guide_legend(order = 2, nrow = 1, title.position = "top",
                                                         label.position = "bottom"))


# All panel coordinates below are in inches, and panel sizes are fixed in inches, so
# row spacing and text positions are exact regardless of the overall canvas.
ROW_IN <- 0.155    # one-line GSEA row
EXTRA_IN <- 0.095  # each additional wrapped line
row_positions <- function(lines) {
  h <- ROW_IN + EXTRA_IN * (lines - 1)
  top <- cumsum(h)
  list(y = max(top) - (top - h / 2) + ROW_IN / 2, total = max(top) + ROW_IN / 2)
}

GSEA_W <- 1.82  # dot + statistics column
dot_panel <- function(d, subtitle, tag) {
  d <- copy(d)[order(-NES)]
  pos <- row_positions(n_lines(d$label))
  d[, y := pos$y]
  p <- ggplot(d, aes(0.12, y)) +
    geom_hline(aes(yintercept = y), colour = "#EEEEEE", linewidth = 0.2) +
    geom_point(aes(fill = NES, size = abs(NES)), shape = 21, colour = "#333333", stroke = 0.25) +
    geom_text(aes(x = 0.32, label = stat), hjust = 0, size = 6.2 / .pt,
              family = "Open Sans", colour = "#222222") +
    fill_scale + size_scale +
    scale_y_continuous(breaks = d$y, labels = d$label, expand = expansion(0)) +
    scale_x_continuous(expand = expansion(0)) +
    coord_cartesian(xlim = c(0, GSEA_W), ylim = c(0, pos$total), clip = "off") +
    labs(x = NULL, y = NULL, subtitle = subtitle, tag = tag) +
    panel_theme() + theme(axis.text.x = element_blank(),
                          axis.text.y = element_text(size = 6.2, lineheight = 0.9, colour = "#222222",
                                                     margin = margin(0, 2, 0, 0)))
  attr(p, "h") <- pos$total
  p
}

pa <- dot_panel(hm, "Hallmark GSEA · 50 sets", "a")
pb <- dot_panel(ex, "Hallmark, Reactome, GO:BP and GO:CC GSEA · 5,880 sets · top 12 per direction", "b")

# ---- c, d: Metascape ---------------------------------------------------------------
# Display labels follow the adopted IsogenicZ Metascape supplement; statistics are the
# saved full-list values (approximate q from rounded log10 q).
term_labels <- c(
  "CORUM:6664" = "STAGA complex, SPT3-linked",
  "WP2011"     = "SREBF / miR33 cholesterol and lipid homeostasis",
  "M17753"     = "BENPORATH MYC/MAX target-gene signature",
  "GO:2000036" = "Regulation of stem-cell population maintenance",
  "M145"       = "PID P53 downstream pathway",
  "M11896"     = "LINDGREN bladder cancer cluster 2A DN",
  "GO:0120545" = "Nucleic acid conformation isomerase activity")

read_metascape <- function(f, keep_supported) {
  d <- fread(f)
  if (keep_supported) d <- d[reported_q_clearly_lt_005 == TRUE]
  stopifnot(all(d$reported_q_clearly_lt_005), !any(d$q_boundary_rounding_ambiguous),
            all(d$GO %in% names(term_labels)))
  d[, .(GO, label = term_labels[GO], logp = LogP, q = approx_q_value,
        n = `#GeneInGOAndHitList`, genes = gsub("\\|", ", ", Hits))]
}
mz <- read_metascape(src$metascape_drugz, keep_supported = TRUE)
mi <- read_metascape(src$metascape_isogenicz, keep_supported = FALSE)
stopifnot(nrow(mz) == 6, nrow(mi) == 4)

# Term name with its contributing genes (grey) beneath, left of the lollipop; n and q at right.
Q_SCALE <- 0.105                # inches per unit of -log10(q)
X_LABEL <- -0.06; X_MIN <- -1.62; X_Q <- 7.6 * Q_SCALE
X_N <- 0.98; X_QTXT <- 1.33; X_MAX <- 1.62
TERM_IN <- 0.102; GENE_IN <- 0.092; PAD_IN <- 0.075
lolli_panel <- function(d, subtitle, tag) {
  d <- copy(d)[order(q, logp)]
  d[, `:=`(term_txt = wrap(label, 31), gene_txt = wrap(genes, 37), x = -log10(q) * Q_SCALE)]
  d[, h := TERM_IN * n_lines(term_txt) + GENE_IN * n_lines(gene_txt) + PAD_IN]
  total <- sum(d$h)
  d[, top := total - (cumsum(h) - h)]
  d[, `:=`(y = top - h / 2, y_term = top - PAD_IN / 2,
           y_gene = top - PAD_IN / 2 - TERM_IN * n_lines(term_txt))]
  head_y <- total + 0.11
  p <- ggplot(d) +
    annotate("segment", x = -log10(0.05) * Q_SCALE, xend = -log10(0.05) * Q_SCALE, y = 0, yend = total,
             colour = "#999999", linetype = "22", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = x, y = y, yend = y), colour = "#D9CCE3", linewidth = 1) +
    geom_point(aes(x = x, y = y), colour = "#6A3D8F", size = 2) +
    geom_text(aes(x = X_LABEL, y = y_term, label = term_txt), hjust = 1, vjust = 1,
              size = 6.2 / .pt, lineheight = 0.88, family = "Open Sans", colour = "#222222") +
    geom_text(aes(x = X_LABEL, y = y_gene, label = gene_txt), hjust = 1, vjust = 1,
              size = 5.5 / .pt, lineheight = 0.88, family = "Open Sans", colour = "#707070") +
    geom_text(aes(x = X_N, y = y, label = n), size = 6.2 / .pt, family = "Open Sans", colour = "#222222") +
    geom_text(aes(x = X_QTXT, y = y, label = fmt_q(q)), size = 6.2 / .pt, family = "Open Sans", colour = "#222222") +
    annotate("text", x = c(X_N, X_QTXT), y = head_y, label = c("n", "q ≈"),
             size = 6.4 / .pt, family = "Open Sans", colour = "#222222") +
    annotate("segment", x = 0, xend = X_Q, y = 0, yend = 0, colour = "#333333", linewidth = 0.3) +
    scale_x_continuous(breaks = c(0, 2, 4, 6) * Q_SCALE, labels = c(0, 2, 4, 6), expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0)) +
    coord_cartesian(xlim = c(X_MIN, X_MAX), ylim = c(0, head_y + 0.06), clip = "off") +
    labs(x = "−log₁₀(Metascape q)", y = NULL, subtitle = subtitle, tag = tag) +
    panel_theme() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6, margin = margin(1.5, 0, 0, 0)),
          axis.title.x = element_text(hjust = (X_Q / 2 - X_MIN) / (X_MAX - X_MIN), size = 6.4,
                                      margin = margin(2, 0, 0, 0)),
          axis.ticks.x = element_line(colour = "#333333", linewidth = 0.3),
          axis.ticks.length.x = unit(2, "pt"))
  attr(p, "h") <- head_y + 0.06
  p
}

pc <- lolli_panel(mz, "Metascape · 48 drugZ depletion candidates", "c")
pd <- lolli_panel(mi, "Metascape · 45 IsogenicZ depletion candidates", "d")

# ---- layout -----------------------------------------------------------------------
# a/b form the left column; c/d and the shared NES legend (extracted from b) the right.
h_of <- function(p) attr(p, "h")
legend <- cowplot::get_plot_component(
  pb + theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "top",
             legend.spacing.x = unit(14, "pt")), "guide-box-bottom", return_all = TRUE)
if (is.list(legend) && !inherits(legend, "grob")) legend <- legend[[1]]
no_legend <- theme(legend.position = "none")
left <- (pa + no_legend) / (pb + no_legend) +
  plot_layout(heights = unit(c(h_of(pa), h_of(pb)), "in"), widths = unit(GSEA_W, "in"))
right <- (pc / pd / plot_spacer() / wrap_elements(full = legend) / plot_spacer()) +
  plot_layout(heights = unit(c(h_of(pc), h_of(pd), 0.12, 0.6, 1), c("in", "in", "in", "in", "null")),
              widths = unit(X_MAX - X_MIN, "in"))
# Each column is wrapped as its own element: patchwork cannot align two nested columns
# whose rows have different fixed-inch heights.
w <- 7.45
fig <- wrap_elements(full = left) | wrap_elements(full = right)
fig <- fig + plot_layout(widths = c(4.2, w - 4.2))

h <- h_of(pa) + h_of(pb) + 0.85
stem <- file.path(fig_dir, "ch04_supp_screen_gsea_combined")
ggsave(paste0(stem, ".pdf"), fig, width = w, height = h, device = cairo_pdf)
ggsave(paste0(stem, ".svg"), fig, width = w, height = h, device = svglite::svglite)
ggsave(paste0(stem, ".png"), fig, width = w, height = h, dpi = 600, device = ragg::agg_png)

plotted <- rbind(hm[, .(panel, analysis = "GSEA", collection = "Hallmark", term = pathway, label, NES, FDR_or_q = FDR, n = NA_integer_, genes = NA_character_)],
                 ex[, .(panel, analysis = "GSEA", collection, term = pathway, label, NES, FDR_or_q = FDR, n = NA_integer_, genes = NA_character_)],
                 mz[, .(panel = "c", analysis = "Metascape drugZ 48", collection = NA_character_, term = GO, label, NES = NA_real_, FDR_or_q = q, n, genes)],
                 mi[, .(panel = "d", analysis = "Metascape IsogenicZ 45", collection = NA_character_, term = GO, label, NES = NA_real_, FDR_or_q = q, n, genes)])
fwrite(plotted, file.path(out_data, "plotted_values.csv"))
jsonlite::write_json(list(
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  script = "figures/chapter4/code/make_screen_gsea_combined.R",
  sources = lapply(src, function(f) list(path = f, md5 = unname(tools::md5sum(f)))),
  panels = list(a = "Hallmark GSEA, 6 sets at BH FDR <= 0.05 within the 50 Hallmark sets",
                b = "Expanded GSEA, top 12 sets per direction at global BH FDR <= 0.05 across 5,880 sets",
                c = "Metascape, 48 drugZ depletion candidates: the 6 cluster representatives with reported q clearly < 0.05",
                d = "Metascape, 45 IsogenicZ depletion candidates: the 4 cluster representatives with q < 0.05"),
  layout = "a/b left column; c/d right column with shared NES legend below",
  shared_nes_limits = c(-lim, lim), size_inches = c(w, h), refit = FALSE
), file.path(out_data, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
message("done: ", stem, "  ", round(w, 2), " x ", round(h, 2), " in")
