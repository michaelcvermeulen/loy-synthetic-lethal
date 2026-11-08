# ---------------------------------------------------------------------------
# DDX3X clone panel - shared data preparation and scales
#
# Sourced by both DDX3X figures so the loading, the normalisation and the
# colour vocabulary are defined once:
#
#   ddPCR_DDX3X_clones.R  - the full three-panel figure (gel calls, raw ratio,
#                           dosage)
#   ddPCR_DDX3X_dosage.R  - the standalone dosage panel + DDX3X gel call
#
# Defines: dat, calls, baseline, flagged, ref_typical, not_run, expected,
#          pal / shp (assay scales), status_pal / status_glyph / status_ink,
#          dodge, base_theme, proj_dir, fig_dir.
#
# Input: projects/ddPCR/data/LOY_clones.xlsx, sheet "DDX3X_LOY"  (read-only)
# ---------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(patchwork)

proj_dir <- "A:/Dropbox/LOY_paralogs/projects/ddPCR"
fig_dir  <- file.path(proj_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(proj_dir, "analysis", "read_xlsx_cells.R"))

xlsx  <- file.path(proj_dir, "data", "LOY_clones.xlsx")
sheet <- "DDX3X_LOY"

# ---------------------------------------------------------------------------
# 1. Load values
#
# The sheet carries a two-tier header (row 1 = "Gel"/"ddPCR" block labels,
# row 2 = gene names) and repeats gene names across blocks, so read it
# head-less and name the columns by position instead of letting readxl
# de-duplicate them into RPS4Y1...1.
# ---------------------------------------------------------------------------

col_names <- c("clone", "gel_RPS4Y1", "gel_KDM5D", "blank1",
               "RPS4Y1", "AMELY/AMELX", "KDM5D", "blank2",
               "call", "blank3", "gel_DDX3X", "notes")

HEADER_ROWS <- 2   # data starts at sheet row 3

raw <- read_excel(xlsx, sheet = sheet, skip = HEADER_ROWS,
                  col_names = col_names, col_types = "text") %>%
  mutate(sheet_row = row_number() + HEADER_ROWS) %>%
  filter(!is.na(clone))

# ---------------------------------------------------------------------------
# 2. Load the colour-coded calls
#
# The Gel, DDX3X and Call columns record their result as a cell background
# colour with the cell left empty, so readxl returns nothing for them. Recover
# the fills and map the four colours the sheet uses onto one vocabulary. The
# same four colours are used in the ddPCR block, where their meaning is
# unambiguous (green on ~0.5, red on 0), which is what fixes the mapping.
# ---------------------------------------------------------------------------

status_of <- function(argb) {
  case_when(
    argb == "FF00B050" ~ "Present / intact",
    argb == "FFFF0000" ~ "Absent / disrupted",
    argb == "FFFFFF00" ~ "Partial / ambiguous",
    argb == "FF7F7F7F" ~ "Not assessed",
    TRUE               ~ NA_character_
  )
}

status_levels <- c("Present / intact", "Partial / ambiguous",
                   "Absent / disrupted", "Not assessed")

# Column letters: B/C = gel Y targets, I = Call, K = gel DDX3X.
call_cols <- c("2" = "RPS4Y1", "3" = "KDM5D", "11" = "DDX3X", "9" = "Call")

cells <- read_xlsx_cells(xlsx, sheet = sheet)

fills <- cells %>%
  filter(col %in% as.integer(names(call_cols))) %>%
  transmute(sheet_row = row,
            track  = factor(unname(call_cols[as.character(col)]),
                            levels = unname(call_cols)),
            status = factor(status_of(argb), levels = status_levels)) %>%
  inner_join(raw %>% select(sheet_row, clone), by = "sheet_row")

# ---------------------------------------------------------------------------
# 2b. Recover the ddPCR copy numbers
#
# Each ddPCR cell is a formula "=<target copies>/<reference copies>", so the
# two counts behind every ratio are recoverable. They are what lets the figure
# put an uncertainty on a ratio instead of drawing a bare point - and the
# reference count is what explains the controls that read above one Y per cell.
#
# Columns E/F/G = the ddPCR block.
# ---------------------------------------------------------------------------

ddpcr_cols <- c("5" = "RPS4Y1", "6" = "AMELY/AMELX", "7" = "KDM5D")

counts <- cells %>%
  filter(col %in% as.integer(names(ddpcr_cols)))
counts <- bind_cols(counts, split_ratio_formula(counts$formula)) %>%
  transmute(sheet_row = row,
            assay_raw = unname(ddpcr_cols[as.character(col)]),
            target_copies, ref_copies)

# ---------------------------------------------------------------------------
# 3. Tidy the ddPCR block
# ---------------------------------------------------------------------------

# Y-target + the arm it sits on. RPS4Y1 (Yp11.3) and AMELY (Yp11.2) report the
# short arm; KDM5D (Yq11.22) reports the long arm. Keeping the arm in the label
# is what makes a partial-Y clone legible.
assay_levels <- c("RPS4Y1 (Yp)", "AMELY/AMELX (Yp)", "KDM5D (Yq)")
assay_cols   <- c("RPS4Y1", "AMELY/AMELX", "KDM5D")

# Ratio expected in a cell that still carries one intact Y.
expected <- tibble(
  assay    = factor(assay_levels, levels = assay_levels),
  exp_ratio = c(0.5, 1.0, 0.5)
)

# Gel band descriptions are free text; shorten them for the figure. The full
# strings stay in the source CSV.
abbrev_band <- function(x) {
  x %>%
    str_replace_all(regex("\\s*bands?\\b", ignore_case = TRUE), "") %>%
    str_replace_all(regex("main\\s+smaller than WT", ignore_case = TRUE), "main <WT") %>%
    str_replace_all(regex("smaller than WT", ignore_case = TRUE), "<WT") %>%
    str_replace_all("\\s*[;/]\\s*", ", ") %>%
    str_replace_all(",\\s*,", ",") %>%
    str_squish()
}

dat <- raw %>%
  mutate(across(all_of(assay_cols), as.numeric)) %>%
  pivot_longer(all_of(assay_cols), names_to = "assay_raw", values_to = "ratio") %>%
  left_join(counts, by = c("sheet_row", "assay_raw")) %>%
  rename(assay = assay_raw) %>%
  mutate(
    group = if_else(str_detect(clone, "CTRL"), "Control clones", "DDX3X LOY clones"),
    group = factor(group, levels = c("DDX3X LOY clones", "Control clones")),
    assay = factor(recode(assay,
                          "RPS4Y1"      = "RPS4Y1 (Yp)",
                          "AMELY/AMELX" = "AMELY/AMELX (Yp)",
                          "KDM5D"       = "KDM5D (Yq)"),
                   levels = assay_levels),
    # short axis label; the DDX3X_ prefix is constant across the sheet
    label = str_remove(clone, "^DDX3X_"),
    band  = abbrev_band(notes)
  ) %>%
  left_join(expected, by = "assay")

# Axis order: clone 1 at the top of each facet.
clone_levels <- rev(unique(dat$label))
dat <- dat %>% mutate(label = factor(label, levels = clone_levels))

calls <- fills %>%
  mutate(label = factor(str_remove(clone, "^DDX3X_"), levels = clone_levels),
         group = factor(if_else(str_detect(clone, "CTRL"),
                                "Control clones", "DDX3X LOY clones"),
                        levels = c("DDX3X LOY clones", "Control clones"))) %>%
  filter(!is.na(label))

# ---------------------------------------------------------------------------
# 4. Normalise: Y dosage per cell, not % LOY against a control median
#
#   dosage = ratio / expected_ratio        (1.0 = one intact Y per cell, 0 = none)
#   LOY%   = (1 - dosage) * 100
#
# The earlier version divided by the MEDIAN CONTROL RATIO, which made LOY% a
# quantity relative to whichever control happened to land in the middle. Two
# problems with that, both of which the copy numbers now let us fix:
#
#   1. It produced negative LOY% - down to -79% - which is not a thing. It read
#      as "this clone gained Y", when all it meant was "this control's ratio
#      came out above the median control's".
#   2. It hid WHY. Dividing by the theoretical expectation instead makes the
#      scale absolute and the deviation attributable: a control that lands at
#      1.8 is 1.8 Y copies per cell, which is a measurement to explain, not a
#      baseline to normalise other clones against.
#
# Dosage cannot go below 0, so the LOY% axis it carries is bounded at 100%.
# It CAN exceed 1 (LOY% below zero), but that now has a concrete reading -
# the assay measured more Y than one copy per cell - and the CI says whether
# that is counting noise or a bad well. See the flagged wells below: it is not
# noise, it is three under-reading reference wells.
#
# The old control-median normalisation is still computed, and both go to the
# source CSV, so the two can be compared.
# ---------------------------------------------------------------------------

# Uncertainty on a ratio of two Poisson counts a/b: the binomial proportion
# p = a/(a+b) has an exact (Clopper-Pearson) interval, and ratio = p/(1-p) is
# monotone in p, so the interval transforms straight across.
#
# CAVEAT, carried into the caption: a and b here are copies/uL, not the raw
# positive-droplet counts, which this sheet does not record. Droplet counts are
# larger, so these intervals are CONSERVATIVE - a well flagged as off is
# genuinely off; a well inside its interval is not thereby proven fine.
poisson_ratio_ci <- function(a, b, conf = 0.95) {
  alpha <- 1 - conf
  lo_p <- qbeta(alpha / 2,     a,     b + 1)
  hi_p <- qbeta(1 - alpha / 2, a + 1, b)
  lo_p[is.na(a) | is.na(b)] <- NA_real_
  hi_p[is.na(a) | is.na(b)] <- NA_real_
  list(lo = lo_p / (1 - lo_p), hi = hi_p / (1 - hi_p))
}

ci <- poisson_ratio_ci(dat$target_copies, dat$ref_copies)

dat <- dat %>%
  mutate(
    dosage    = ratio / exp_ratio,
    dosage_lo = ci$lo / exp_ratio,
    dosage_hi = ci$hi / exp_ratio,
    loy_pct   = (1 - dosage) * 100
  )

# Kept only for comparison against the previous version of this figure.
baseline <- dat %>%
  filter(group == "Control clones", !is.na(ratio)) %>%
  group_by(assay) %>%
  summarise(ctrl_median = median(ratio),
            ctrl_mean   = mean(ratio),
            ctrl_min    = min(ratio),
            ctrl_max    = max(ratio),
            ctrl_n      = n(), .groups = "drop")

dat <- dat %>%
  left_join(baseline, by = "assay") %>%
  mutate(loy_pct_vs_ctrl_median = (1 - ratio / ctrl_median) * 100)

# Control wells whose dosage interval excludes one intact Y. These are the
# rows that used to show up as negative LOY%. They are not counting noise -
# every one of them is an under-reading REFERENCE well (a low denominator
# inflates the ratio), which is a well to repeat, not a property of the clone.
flagged <- dat %>%
  filter(group == "Control clones", !is.na(dosage_lo),
         dosage_lo > 1 | dosage_hi < 1) %>%
  select(clone, label, assay, target_copies, ref_copies, dosage,
         dosage_lo, dosage_hi)

ref_typical <- dat %>%
  filter(!is.na(ref_copies)) %>%
  group_by(assay) %>%
  summarise(ref_median = median(ref_copies), .groups = "drop")

# Clones with no ddPCR data are kept as rows and annotated with the reason,
# taken from their gel call - a silently missing clone reads as "not made",
# and "not run" alone doesn't say why.
gel_y <- calls %>%
  filter(track %in% c("RPS4Y1", "KDM5D")) %>%
  group_by(group, label) %>%
  summarise(y_present = any(status == "Present / intact"), .groups = "drop")

not_run <- dat %>%
  group_by(group, label) %>%
  summarise(all_na = all(is.na(ratio)), .groups = "drop") %>%
  filter(all_na) %>%
  left_join(gel_y, by = c("group", "label")) %>%
  mutate(why = if_else(y_present %in% TRUE,
                       "Y detected on gel \u2014 not taken to ddPCR",
                       "ddPCR not run"))

# ---------------------------------------------------------------------------
# 5. Figure
# ---------------------------------------------------------------------------

# Assay hues: categorical slots 1-3 of the validated default palette. Shape
# doubles the hue (secondary encoding), so the aqua slot - which sits below
# 3:1 against a white surface - never carries identity on colour alone.
pal <- c("RPS4Y1 (Yp)"      = "#2a78d6",
         "AMELY/AMELX (Yp)" = "#eb6834",
         "KDM5D (Yq)"       = "#1baf7a")
shp <- c("RPS4Y1 (Yp)" = 21, "AMELY/AMELX (Yp)" = 22, "KDM5D (Yq)" = 24)

# Call tiles keep the workbook's own four colours, so the panel is recognisably
# the same coding as the spreadsheet. Each tile also carries a glyph, so the
# status never rests on colour alone (and survives print / CVD).
status_pal <- c("Present / intact"    = "#00B050",
                "Partial / ambiguous" = "#FFD200",
                "Absent / disrupted"  = "#E02020",
                "Not assessed"        = "#9AA0A6")
status_glyph <- c("Present / intact" = "+", "Partial / ambiguous" = "?",
                  "Absent / disrupted" = "\u2212", "Not assessed" = "")
status_ink <- c("Present / intact" = "white", "Partial / ambiguous" = "grey15",
                "Absent / disrupted" = "white", "Not assessed" = "white")

dodge <- position_dodge(width = 0.72)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 10.5),
    plot.subtitle      = element_text(colour = "grey35", size = 8.2),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    strip.text.y       = element_text(face = "bold", size = 9.5),
    axis.text.y        = element_text(size = 8.5),
    legend.position    = "bottom",
    legend.title       = element_blank()
  )
