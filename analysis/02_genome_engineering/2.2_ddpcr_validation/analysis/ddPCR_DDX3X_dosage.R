# ---------------------------------------------------------------------------
# DDX3X clones: Y dosage per cell, with the DDX3X gel call alongside
#
# The dosage panel of the three-panel figure, on its own, with two changes:
#
#   - LOY6 is dropped. It has no ddPCR data (its gel found Y, so it never went
#     to ddPCR), and on a ddPCR-only figure an empty annotated row is just a
#     gap. It stays in ddPCR_DDX3X_clones.R, where panel A carries its gel call
#     and the row therefore says something.
#   - the DDX3X gel call is appended at the right, so Y status and edit status
#     can be read off one row.
#
# Input : projects/ddPCR/data/LOY_clones.xlsx, sheet "DDX3X_LOY"  (read-only)
# Output: projects/ddPCR/figures/DDX3X_clones_dosage.{png,pdf}
#         projects/ddPCR/figures/DDX3X_clones_dosage_source.csv
# ---------------------------------------------------------------------------

source(file.path("A:/Dropbox/LOY_paralogs/projects/ddPCR",
                 "analysis", "ddPCR_DDX3X_common.R"))

# ---------------------------------------------------------------------------
# 1. Drop the clones with no ddPCR data
#
# Data-driven rather than a hard-coded "LOY6": if another clone is ever added
# without ddPCR, it drops out the same way instead of appearing as a blank row.
# ---------------------------------------------------------------------------

dropped <- not_run$label
keep_levels <- setdiff(levels(dat$label), as.character(dropped))

ddat <- dat %>%
  filter(!label %in% dropped) %>%
  mutate(label = factor(as.character(label), levels = keep_levels))

dcalls <- calls %>%
  filter(!label %in% dropped) %>%
  mutate(label = factor(as.character(label), levels = keep_levels))

# ---------------------------------------------------------------------------
# 2. Dosage panel
#
# Same construction as panel C of the three-panel figure: dosage on the x axis
# so the quantity has a floor at 0, with % LOY riding on top as a second
# labelling of the same axis, its breaks stopping at 0% because beyond one Y
# per cell the reading is over-measurement, not negative loss.
# ---------------------------------------------------------------------------

# No interval bars on this figure. They are still computed in common.R and
# still drawn in the three-panel figure, but the only intervals derivable from
# this sheet come from treating copies/uL as Poisson counts, which overstates
# the width by ~4x against real droplet-count precision. Left off here rather
# than shown at the wrong width. The copy numbers are in the source CSV.
#
# Axis top now comes from the points alone, not from a CI end, so the panel
# tightens from 2.5 to the highest actual measurement.
x_hi <- max(1.6, ceiling(max(ddat$dosage, na.rm = TRUE) * 10) / 10)

pDose <- ggplot(ddat, aes(x = dosage, y = label, colour = assay, shape = assay)) +
  geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.45) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(position = dodge, size = 2.4, stroke = 0.7, fill = "white",
             na.rm = TRUE) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_shape_manual(values = shp, drop = FALSE) +
  scale_x_continuous(
    limits = c(-0.03, x_hi), breaks = seq(0, 2.5, 0.5),
    sec.axis = sec_axis(transform = ~ (1 - .) * 100, name = "% LOY",
                        breaks = c(100, 75, 50, 25, 0))) +
  labs(x = "Y copies per cell", y = NULL,
       title = "Y dosage per cell",
       subtitle = "1.0 = one intact Y, 0 = none")

# The one clone whose two arms disagree, labelled in % LOY - the project's
# usual language - even though the axis underneath is dosage.
split_clone <- ddat %>%
  filter(group == "DDX3X LOY clones", !is.na(loy_pct)) %>%
  group_by(label) %>%
  filter(max(loy_pct) - min(loy_pct) > 15) %>%
  ungroup()

# One label for the clone, not one per assay: the two Yp assays land within
# 0.03 of each other on both the value and the CI end, so separate labels
# collide however they are anchored.
if (nrow(split_clone) > 0) {
  yp_lab <- split_clone %>%
    filter(str_detect(assay, "Yp")) %>%
    group_by(group, label) %>%
    summarise(lo = min(loy_pct), hi = max(loy_pct), .groups = "drop") %>%
    mutate(txt = sprintf("Yp retained — %.0f–%.0f%% LOY", lo, hi))

  pDose <- pDose +
    geom_text(data = yp_lab, inherit.aes = FALSE,
              aes(x = 0.42, y = label, label = txt),
              hjust = 0, vjust = 0.4, size = 2.8,
              fontface = "italic", colour = "grey30")
}

# ---------------------------------------------------------------------------
# 3. DDX3X gel call, appended
#
# One tile per clone plus the band pattern it was called from, on a continuous
# x so the text can sit beside the tile on the same axis.
# ---------------------------------------------------------------------------

ddx <- dcalls %>% filter(track == "DDX3X")
band_txt <- ddat %>% distinct(group, label, band) %>% filter(!is.na(band))

pGel <- ggplot(ddx, aes(x = 1, y = label, fill = status)) +
  geom_tile(width = 0.62, height = 0.78, colour = "white", linewidth = 0.9) +
  geom_text(aes(label = status_glyph[as.character(status)], colour = status),
            size = 3.2, fontface = "bold", show.legend = FALSE) +
  geom_text(data = band_txt, inherit.aes = FALSE,
            aes(x = 1.45, y = label, label = band),
            hjust = 0, size = 2.6, colour = "grey30") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = status_pal, drop = FALSE) +
  scale_colour_manual(values = status_ink, drop = FALSE, guide = "none") +
  scale_x_continuous(breaks = NULL, limits = c(0.6, 5.0),
                     expand = expansion(0)) +
  labs(x = NULL, y = NULL, title = "DDX3X gel",
       subtitle = "call · band pattern")
# Clone labels are repeated on this panel rather than blanked: it sits across a
# gap from the dosage panel, so a reader tracking one clone needs the row named
# on both sides.

# ---------------------------------------------------------------------------
# 4. Compose
# ---------------------------------------------------------------------------

# No caption and no subtitle on this figure - it is the clean version, for
# talks and slides. What the caption used to carry (the LOY6 omission and why,
# the reference wells reading off, what dosage is normalised against) is in the
# README and in the console summary below, and the three-panel figure stays the
# annotated record.
#
# base_theme goes on BEFORE the per-panel overrides, not after: it specifies
# strip.text.y itself, so a blank set on the panel earlier would be overwritten
# and the facet strip would print on both panels.
pDose_f <- pDose + base_theme + theme(strip.text.y = element_blank())
pGel_f  <- pGel  + base_theme + theme(panel.grid.major.y = element_blank())

p <- (pDose_f | pGel_f) +
  plot_layout(widths = c(4.1, 1.25), guides = "collect") +
  plot_annotation(
    title = "DDX3X clones: Y dosage per cell, with the DDX3X gel call",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(fig_dir, "DDX3X_clones_dosage.png"), p,
       width = 11, height = 6.4, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "DDX3X_clones_dosage.pdf"), p,
       width = 11, height = 6.4, bg = "white")

# ---------------------------------------------------------------------------
# 5. Traceable source values
# ---------------------------------------------------------------------------

# No dosage_lo/dosage_hi here: the figure no longer draws intervals, and the
# CSV is meant to be a faithful record of what is plotted. The copy numbers the
# intervals were derived from are kept, so they can be recomputed.
src <- ddat %>%
  select(clone, group, assay, target_copies, ref_copies, ratio, exp_ratio,
         dosage, loy_pct, gel_band = notes) %>%
  left_join(ddx %>% select(clone, ddx3x_gel_call = status), by = "clone") %>%
  arrange(group, clone, assay)

src_path <- file.path(fig_dir, "DDX3X_clones_dosage_source.csv")

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

cat("\n--- clones omitted (no ddPCR data) ---\n")
print(as.data.frame(not_run %>% select(label, why)))
cat("\n--- Y dosage and DDX3X gel call ---\n")
print(as.data.frame(
  ddat %>%
    mutate(v = if_else(is.na(dosage), NA_character_, sprintf("%.3f", dosage))) %>%
    select(clone, assay, v) %>%
    pivot_wider(names_from = assay, values_from = v) %>%
    left_join(ddx %>% select(clone, DDX3X_gel = status), by = "clone")))
cat("\nWrote figures to:", fig_dir, "\n")
