# ---------------------------------------------------------------------------
# DDX3X clone panel: gel calls + ddPCR loss-of-Y
#
# Figure: the full genotype record for the DDX3X-edited clones (DDX3X_LOY1-9)
# and their sibling controls (DDX3X_CTRL1-9).
#
#   A. gel calls (Y targets + DDX3X) and the sheet's own Call, recovered from
#      the CELL FILL COLOURS - those columns hold no text.
#   B. raw ddPCR ratio, against the value expected in a Y-bearing cell:
#        RPS4Y1, KDM5D  -> 0.5  (1 Y copy vs the 2-copy RNaseP reference)
#        AMELY/AMELX    -> 1.0  (1 Y copy vs the single X in a male cell)
#   C. the same data normalised within assay to the control clones -> % LOY.
#
# Input : projects/ddPCR/data/LOY_clones.xlsx, sheet "DDX3X_LOY"  (read-only)
# Output: projects/ddPCR/figures/DDX3X_clones_LOY.{png,pdf}
#         projects/ddPCR/figures/DDX3X_clones_LOY_source.csv
# ---------------------------------------------------------------------------

# Loading, normalisation and the shared colour vocabulary live in the common
# file so this figure and the standalone dosage panel cannot drift apart.
source(file.path("A:/Dropbox/LOY_paralogs/projects/ddPCR",
                 "analysis", "ddPCR_DDX3X_common.R"))

# --- Panel A: gel calls -----------------------------------------------------
#
# Continuous x rather than discrete, so the band-description text can sit in
# free space to the right of the four tile columns on the same axis.

tile_x <- setNames(seq_along(levels(calls$track)), levels(calls$track))
calls  <- calls %>% mutate(x = tile_x[as.character(track)])

band_txt <- dat %>% distinct(group, label, band) %>% filter(!is.na(band))

pA <- ggplot(calls, aes(x = x, y = label, fill = status)) +
  geom_tile(width = 0.86, height = 0.78, colour = "white", linewidth = 0.9) +
  geom_text(aes(label = status_glyph[as.character(status)],
                colour = status),
            size = 3.1, fontface = "bold", show.legend = FALSE) +
  geom_text(data = band_txt, inherit.aes = FALSE,
            aes(x = 4.75, y = label, label = band),
            hjust = 0, size = 2.5, colour = "grey30") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = status_pal, drop = FALSE) +
  scale_colour_manual(values = status_ink, drop = FALSE, guide = "none") +
  scale_x_continuous(breaks = tile_x, labels = names(tile_x),
                     limits = c(0.4, 8.1), expand = expansion(0)) +
  labs(title = "A  Gel calls", x = NULL, y = NULL,
       subtitle = "Y targets, DDX3X, and the sheet's own Call\nDDX3X band pattern at right") +
  theme(panel.grid.major.y = element_blank(), strip.text.y = element_blank())

# --- Panel B: raw ratios ----------------------------------------------------
#
# Plot from the full table, not from the non-missing rows: a clone with no
# ddPCR data (LOY6) would otherwise have its factor level dropped from the
# panel and be re-appended out of order by the annotation layer. geom_point
# drops the NA rows itself.

pB <- ggplot(dat, aes(x = ratio, y = label, colour = assay, shape = assay)) +
  geom_vline(xintercept = c(0.5, 1.0), colour = "grey60",
             linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(position = dodge, size = 2.2, stroke = 0.6, fill = "white",
             na.rm = TRUE) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_shape_manual(values = shp, drop = FALSE) +
  scale_x_continuous(limits = c(-0.06, 1.95), breaks = seq(0, 1.75, 0.25)) +
  labs(title = "B  Raw ddPCR ratio", x = "Y target / reference", y = NULL,
       subtitle = paste0("Dashed = ratio expected with one intact Y:\n",
                         "0.5 for RPS4Y1 & KDM5D (1 Y : 2-copy RNaseP), 1.0 for AMELY : AMELX")) +
  theme(strip.text.y = element_blank())

# --- Panel C: Y dosage per cell, carrying a bounded % LOY axis ---------------
#
# Plotted on dosage rather than on LOY% so the quantity has a floor: 0 Y copies
# per cell is as low as it goes, which is 100% LOY. The % LOY axis rides on top
# as a second labelling of the same axis (a linear re-expression, not a second
# measure), and its breaks stop at 0% - so the region beyond one Y per cell is
# deliberately left unlabelled in LOY terms. That is the honest treatment of
# what used to be drawn as negative LOY: it is over-measured Y, not lost Y.

x_hi <- max(1.6, ceiling(max(dat$dosage_hi, dat$dosage, na.rm = TRUE) * 10) / 10)

pC <- ggplot(dat, aes(x = dosage, y = label, colour = assay, shape = assay)) +
  geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.45) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_linerange(aes(xmin = dosage_lo, xmax = dosage_hi),
                 position = dodge, linewidth = 0.45, na.rm = TRUE,
                 show.legend = FALSE) +
  geom_point(position = dodge, size = 2.2, stroke = 0.6, fill = "white",
             na.rm = TRUE) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_shape_manual(values = shp, drop = FALSE) +
  scale_x_continuous(
    limits = c(-0.03, x_hi), breaks = seq(0, 2.5, 0.5),
    sec.axis = sec_axis(transform = ~ (1 - .) * 100, name = "% LOY",
                        breaks = c(100, 75, 50, 25, 0))) +
  labs(title = "C  Y dosage per cell", x = "Y copies per cell", y = NULL,
       subtitle = "dosage = ratio / expected \u00b7 1.0 = one intact Y, 0 = none\nbars = 95% CI from the ddPCR copy numbers")

# Explain the empty ddPCR rows, in both ddPCR panels.
# geom_label, not geom_text: the annotation sits across the dashed expected-
# ratio lines in B, and an unfilled label leaves them running through the words.
if (nrow(not_run) > 0) {
  pB <- pB + geom_label(data = not_run, inherit.aes = FALSE,
                        aes(x = 1.0, y = label, label = why),
                        size = 2.5, colour = "grey45", fontface = "italic",
                        fill = "white", linewidth = 0, label.padding = unit(1, "pt"))
  pC <- pC + geom_label(data = not_run, inherit.aes = FALSE,
                        aes(x = 0.62, y = label, label = why),
                        size = 2.5, colour = "grey45", fontface = "italic",
                        fill = "white", linewidth = 0, label.padding = unit(1, "pt"))
}

# One selective direct label: the only LOY clone whose two arms disagree.
# Labelling every point would be noise; this is the point the figure exists to
# show. The Yq point needs no label - it sits on the 100% end, which already
# reads as "fully lost".
#
# One label for the clone rather than one per assay: the two Yp assays land
# within 0.03 of each other on both the value and the CI end, so separate
# labels collide however they are anchored. Quoted in % LOY, the project's
# usual language, even though the axis underneath is dosage.
split_clone <- dat %>%
  filter(group == "DDX3X LOY clones", !is.na(loy_pct)) %>%
  group_by(label) %>%
  filter(max(loy_pct) - min(loy_pct) > 15) %>%
  ungroup()

if (nrow(split_clone) > 0) {
  yp_lab <- split_clone %>%
    filter(str_detect(assay, "Yp")) %>%
    group_by(group, label) %>%
    summarise(lo = min(loy_pct), hi = max(loy_pct), .groups = "drop") %>%
    mutate(txt = sprintf("Yp retained — %.0f–%.0f%% LOY", lo, hi))

  pC <- pC +
    geom_text(data = yp_lab, inherit.aes = FALSE,
              aes(x = 0.42, y = label, label = txt),
              hjust = 0, vjust = 0.4, size = 2.6,
              fontface = "italic", colour = "grey30")
}

# ---------------------------------------------------------------------------
# 6. Compose
# ---------------------------------------------------------------------------

# How far the control clones scatter in each assay - the practical noise floor,
# now quoted in dosage rather than in LOY%, because that is the scale on which
# it is bounded and interpretable.
spread_txt <- dat %>%
  filter(group == "Control clones", !is.na(dosage)) %>%
  group_by(assay) %>%
  summarise(lo = min(dosage), hi = max(dosage), .groups = "drop") %>%
  mutate(txt = sprintf("%s %.2f–%.2f", assay, lo, hi)) %>%
  pull(txt) %>%
  paste(collapse = "; ")

# Name the flagged wells and their under-reading references in the caption -
# the whole point of the change is that these are attributable, not noise.
flag_txt <- if (nrow(flagged) == 0) "none" else {
  txt <- flagged %>%
    left_join(ref_typical, by = "assay") %>%
    mutate(txt = sprintf("%s %s (ref %.0f vs %.0f typical)",
                         label, str_remove(assay, " \\(Y[pq]\\)"),
                         ref_copies, ref_median)) %>%
    pull(txt)
  # Wrap: the full list on one line overruns the figure width.
  paste(
    vapply(split(txt, (seq_along(txt) - 1) %/% 2),
           paste, character(1), collapse = "; "),
    collapse = ";\n")
}

# base_theme goes on BEFORE the per-panel overrides, not after: it specifies
# strip.text.y itself, so the blanks set on pA and pB when they were built get
# overwritten and the facet strip prints in both gutters as well as at the far
# right. Only the rightmost panel should carry it.
pA_f <- pA + base_theme + theme(panel.grid.major.y = element_blank(),
                                strip.text.y = element_blank())
pB_f <- pB + base_theme + theme(strip.text.y = element_blank())
pC_f <- pC + base_theme

p <- (pA_f | pB_f | pC_f) +
  plot_layout(widths = c(1.18, 1.08, 1.22), guides = "collect") +
  plot_annotation(
    title = "DDX3X clone panel: gel calls and loss of Y by ddPCR",
    subtitle = paste0(
      "RPS4Y1 and AMELY report the Y short arm; KDM5D reports the long arm. ",
      "Gel and Call columns of the sheet hold no text \u2014 they are read from the cell fill colours."),
    caption = paste0(
      "One well per clone per assay. Bars are 95% intervals from the ddPCR copy numbers behind each ratio, ",
      "NOT replicates \u2014 there are none. They are also conservative:\n",
      "the sheet records copies/\u00b5L rather than droplet counts, so a well whose bar clears 1.0 is genuinely off, ",
      "but a well inside its bar is not thereby shown to be sound.\n",
      "Dosage is the ratio over its expected value (0.5 or 1.0), not over a control median. ",
      "That is what keeps the % LOY axis bounded at 100% and stops controls reading as negative LOY;\n",
      "the % LOY axis is deliberately unlabelled beyond one Y per cell, because more Y than one copy is over-measurement, not negative loss. ",
      "Control dosage spans ", spread_txt, ".\n",
      "Controls whose interval excludes one intact Y \u2014 the rows that used to plot as negative LOY. ",
      "Every one is a reference well reading off, not a Y gain:\n",
      flag_txt, ".\n",
      "LOY9 is the one clone the assays split on: both Yp targets sit near 0.2 Y per cell while Yq is undetectable. ",
      "That is well outside the noise, but n=1 and it has no gel \u2014\n",
      "it is equally consistent with a retained Yp fragment and with a ~20% Y-bearing subpopulation.\n",
      "Source: data/LOY_clones.xlsx, sheet DDX3X_LOY (values, formulas and cell fills)."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13.5),
      plot.subtitle = element_text(colour = "grey35", size = 9.2),
      plot.caption  = element_text(colour = "grey45", size = 7, hjust = 0)
    )
  ) &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(fig_dir, "DDX3X_clones_LOY.png"), p,
       width = 15.5, height = 7, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "DDX3X_clones_LOY.pdf"), p,
       width = 15.5, height = 7, bg = "white")

# ---------------------------------------------------------------------------
# 7. Traceable source values
# ---------------------------------------------------------------------------

call_wide <- calls %>%
  select(clone, track, status) %>%
  mutate(track = paste0("call_", track)) %>%
  pivot_wider(names_from = track, values_from = status)

src <- dat %>%
  select(clone, group, assay,
         target_copies, ref_copies, ratio, exp_ratio,
         dosage, dosage_lo, dosage_hi, loy_pct,
         ctrl_median, ctrl_mean, loy_pct_vs_ctrl_median,
         gel_band = notes) %>%
  left_join(call_wide, by = "clone") %>%
  arrange(group, clone, assay)

src_path <- file.path(fig_dir, "DDX3X_clones_LOY_source.csv")

# The CSVs in this project are routinely open in Excel, which locks them on
# Windows. Don't abort after the figures have already been written.
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

cat("\n--- control dosage: empirical median vs the expected value ---\n")
print(as.data.frame(
  baseline %>%
    left_join(expected, by = "assay") %>%
    mutate(median_dosage = ctrl_median / exp_ratio) %>%
    select(assay, ctrl_n, exp_ratio, ctrl_median, median_dosage,
           ctrl_min, ctrl_max)))

cat("\n--- control wells whose interval excludes one intact Y ---\n")
if (nrow(flagged) == 0) cat("none\n") else
  print(as.data.frame(flagged %>% select(-label)))

cat("\n--- LOY clones: Y dosage and % LOY ---\n")
print(as.data.frame(
  dat %>% filter(group == "DDX3X LOY clones") %>%
    mutate(v = if_else(is.na(dosage), NA_character_,
                       sprintf("%.3f (%.0f%% LOY)", dosage, loy_pct))) %>%
    select(clone, assay, v) %>%
    pivot_wider(names_from = assay, values_from = v)))

cat("\n--- % LOY: new absolute scale vs the old control-median scale ---\n")
print(as.data.frame(
  dat %>% filter(!is.na(loy_pct)) %>%
    transmute(clone, assay,
              vs_expected    = round(loy_pct, 1),
              vs_ctrl_median = round(loy_pct_vs_ctrl_median, 1),
              shift          = round(loy_pct - loy_pct_vs_ctrl_median, 1)) %>%
    arrange(desc(abs(shift))) %>%
    head(8)))

cat("\nWrote figures to:", fig_dir, "\n")
