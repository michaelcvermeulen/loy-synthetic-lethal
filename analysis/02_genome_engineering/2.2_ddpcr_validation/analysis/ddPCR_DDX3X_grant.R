# ---------------------------------------------------------------------------
# DDX3X clones, grant figure: the clean clones only
#
# Condensed single-panel version for a grant. Restricted to the clones that
# will go into the Y+ and LOY pools - see EXCLUDE below for what is left out
# and why, which is the part a reader is entitled to ask about.
#
# Deliberately carries no caption: a grant figure's caveats belong in the
# figure legend the author writes, not burned into the image. The full
# annotated version is DDX3X_clones_LOY.{png,pdf}.
#
# Input : projects/ddPCR/data/LOY_clones.xlsx, sheet "DDX3X_LOY"  (read-only)
# Output: projects/ddPCR/figures/DDX3X_clones_grant.{png,pdf}
#         projects/ddPCR/figures/DDX3X_clones_grant_source.csv
# ---------------------------------------------------------------------------

source(file.path("A:/Dropbox/LOY_paralogs/projects/ddPCR",
                 "analysis", "ddPCR_DDX3X_common.R"))

# ---------------------------------------------------------------------------
# 1. The exclusion list
#
# Held as a named vector rather than a bare c(...) so the figure cannot end up
# with an unexplained set of clones. The reasons print to the console on every
# run and go into the source CSV's header comment.
# ---------------------------------------------------------------------------

EXCLUDE <- c(
  LOY6  = "Y detected on gel; never taken to ddPCR",
  LOY9  = "partial Yp retention (~0.2 Y per cell) - not cleanly LOY",
  LOY2  = "DDX3X gel call ambiguous (single band, smaller than WT)",
  CTRL2 = "DDX3X gel call ambiguous (blank band); AMELX reference well low",
  CTRL5 = "DDX3X gel band is WT size - not edited at DDX3X",
  CTRL9 = "DDX3X and Y gels not assessed; no RPS4Y1 value"
)

keep_levels <- setdiff(levels(dat$label), names(EXCLUDE))

gdat <- dat %>%
  filter(!label %in% names(EXCLUDE), !is.na(dosage)) %>%
  mutate(label = factor(as.character(label), levels = keep_levels))

# The pools these clones are destined for, in the language of the grant.
n_by_group <- gdat %>%
  group_by(group) %>%
  summarise(n = n_distinct(label), .groups = "drop")

gdat <- gdat %>%
  left_join(n_by_group, by = "group") %>%
  mutate(pool = factor(
    if_else(group == "Control clones",
            sprintf("Y+ clones (n=%d)", n),
            sprintf("LOY clones (n=%d)", n)),
    levels = c(sprintf("LOY clones (n=%d)",
                       n_by_group$n[n_by_group$group == "DDX3X LOY clones"]),
               sprintf("Y+ clones (n=%d)",
                       n_by_group$n[n_by_group$group == "Control clones"]))))

# Every retained clone carries the same DDX3X call, so a per-clone tile column
# would be twelve identical tiles. Stated once instead - but computed, not
# asserted, so it cannot go stale if the exclusion list changes.
ddx_status <- calls %>%
  filter(track == "DDX3X", label %in% keep_levels) %>%
  pull(status) %>% as.character() %>% unique()

ddx_note <- if (length(ddx_status) == 1 && ddx_status == "Absent / disrupted") {
  "All clones DDX3X-disrupted on gel"
} else {
  paste0("DDX3X gel calls: ", paste(ddx_status, collapse = ", "))
}

# ---------------------------------------------------------------------------
# 2. Figure
#
# Dosage on x so the quantity has a floor at 0, with % LOY riding on top as a
# second labelling of the same axis - reviewers read % LOY, but it is the
# bounded scale underneath that keeps controls off the negative numbers.
# ---------------------------------------------------------------------------

# facet_wrap gives every panel the same height, which is only right while the
# two groups hold the same number of clones. Warn rather than silently squash
# the rows if a change to EXCLUDE makes them unequal.
if (length(unique(n_by_group$n)) > 1) {
  warning("Groups are no longer equal in size (",
          paste(sprintf("%s=%d", n_by_group$group, n_by_group$n), collapse = ", "),
          "). facet_wrap will force equal panel heights and distort the row ",
          "spacing - switch to facet_grid(pool ~ ., space = \"free_y\").",
          call. = FALSE)
}

# Sized for a small single-column grant slot. Everything here is scaled for a
# ~3.4in wide figure that may be shrunk again in the layout - line weights and
# grid stay light but the type does not go below ~6pt, which is the floor for
# print legibility.
grant_theme <- theme_minimal(base_size = 7) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.22),
    panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.22),
    panel.spacing.y    = unit(5, "pt"),
    # Horizontal strips above each group, not rotated ones down the side: at
    # this size vertical text is the least legible thing on the figure.
    strip.text         = element_text(face = "bold", size = 7, hjust = 0,
                                      margin = margin(b = 2)),
    axis.text          = element_text(size = 6.2),
    axis.title         = element_text(size = 6.8),
    axis.title.x.top   = element_text(margin = margin(b = 2)),
    axis.title.x.bottom = element_text(margin = margin(t = 2)),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.margin      = margin(t = -4),
    legend.text        = element_text(size = 6.2),
    legend.key.size    = unit(7, "pt"),
    legend.spacing.x   = unit(2, "pt"),
    plot.margin        = margin(4, 6, 3, 3)
  )

p <- ggplot(gdat, aes(x = dosage, y = label, colour = assay, shape = assay)) +
  geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
  geom_point(position = position_dodge(width = 0.6),
             size = 1.15, stroke = 0.42, fill = "white") +
  # facet_wrap rather than facet_grid: the two groups hold the same number of
  # clones, so equal-height panels are correct and facet_wrap can put the strip
  # on top. (If the exclusion list ever makes the groups unequal, this needs
  # facet_grid with space = "free_y" again, and rotated strips with it.)
  facet_wrap(~ pool, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_shape_manual(values = shp, drop = FALSE) +
  scale_x_continuous(
    limits = c(-0.03, 1.55), breaks = seq(0, 1.5, 0.25),
    sec.axis = sec_axis(transform = ~ (1 - .) * 100, name = "% LOY",
                        breaks = c(100, 75, 50, 25, 0))) +
  # No title and no subtitle: at this size they cost a disproportionate share
  # of the box, and both belong in the figure legend the author writes. What
  # the subtitle carried - that every clone shown is DDX3X-disrupted - is
  # printed to the console below so it can go straight into that legend.
  labs(x = "Y copies per cell", y = NULL) +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  grant_theme

ggsave(file.path(fig_dir, "DDX3X_clones_grant.png"), p,
       width = 3.4, height = 2.5, dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "DDX3X_clones_grant.pdf"), p,
       width = 3.4, height = 2.5, bg = "white")

# ---------------------------------------------------------------------------
# 3. Traceable source values
# ---------------------------------------------------------------------------

src <- gdat %>%
  select(clone, pool, assay, target_copies, ref_copies, ratio, exp_ratio,
         dosage, loy_pct) %>%
  arrange(pool, clone, assay)

src_path <- file.path(fig_dir, "DDX3X_clones_grant_source.csv")

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

cat("\n--- FOR THE FIGURE LEGEND ---\n")
cat("  ", ddx_note, "\n", sep = "")
cat("   (the figure no longer states this - the subtitle was removed to save\n",
    "    space, so it has to be written into the legend)\n", sep = "")

cat("\n--- clones EXCLUDED from the grant figure ---\n")
for (nm in names(EXCLUDE)) cat(sprintf("  %-6s %s\n", nm, EXCLUDE[[nm]]))

cat("\n--- clones retained ---\n")
print(as.data.frame(gdat %>% group_by(pool) %>%
  summarise(clones = paste(sort(unique(as.character(label))), collapse = ", "),
            .groups = "drop")))

cat("\n--- Y+ clone dosages (LOY clones are all exactly 0) ---\n")
print(as.data.frame(gdat %>% filter(str_detect(pool, "^Y\\+")) %>%
  select(label, assay, dosage) %>%
  pivot_wider(names_from = assay, values_from = dosage) %>%
  mutate(across(where(is.numeric), ~round(., 3)))))

# The wells that will show as scatter away from 1.0. They survived the
# exclusion list because the CLONE is fine - it is the reference well that read
# off - but a reviewer will see them, so they should not be a surprise.
still_flagged <- flagged %>% filter(!label %in% names(EXCLUDE))
cat("\n--- retained wells whose reference read off (visible scatter) ---\n")
if (nrow(still_flagged) == 0) cat("none\n") else
  print(as.data.frame(still_flagged %>%
    left_join(ref_typical, by = "assay") %>%
    transmute(label, assay, dosage = round(dosage, 2),
              ref_copies, typical = round(ref_median))))

cat("\nWrote figures to:", fig_dir, "\n")
