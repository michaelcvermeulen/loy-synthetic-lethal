# ddPCR project

Analysis project for the ddPCR LOY experiments. Set up 2026-08-26.

## Layout

| Path | Contents |
|------|----------|
| `analysis/` | R scripts and Rmd reports, copied from `code/` |
| `data/` | Curated summary tables only (see note below) |
| `figures/` | Generated figures + their traceable source values |

## Where the raw data lives

**Raw instrument data was NOT copied** — `A:\Dropbox\LOY_paralogs\data\ddPCR`
is 7.3 GB of `.qlb`/`.qlp`/`.qlr` binaries and per-well amplitude CSVs.
Scripts reference it in place at that path.

`data/` holds only the small curated tables: `summary_table_ddPCR.csv`
(master summary), `summary_all_INFO_TABLES.csv`, the two Jan 2025 extracts,
and the XY/XO dilution ratios.

## Figures

### `HT1080_CPQ53_LOY_rounds.{png,pdf}`

HT1080 LOY across successive CP+Q53 nucleofection rounds, from
`analysis/ddPCR_HT1080_CPQ53_rounds.R`.

Rounds come from the `Version` column of the master summary:
`0` = Parental, `1` = first nucleofection, `2` = second round.

| Round | DDX3Y/RNaseP | AMELY/AMELX |
|-------|--------------|-------------|
| Parental | 0% (baseline, n=2) | 0% (baseline, n=2) |
| Round 1 | 1.2% (n=1) | 0.3% (n=1) |
| Round 2 | 88.9% (n=4) | 79.3% (n=7) |

Round 1 produced no detectable LOY; the effect appears only after round 2.
Both assay chemistries agree independently.

**Caveats carried in the figure:**

- Round 1 is a **single well** (C5, Jan 26 2025). No replicates exist anywhere
  in the repo. The "round 1 achieved no LOY" claim rests on one measurement
  and wants a repeat before it goes into a paper.
- The parental baseline uses only the Feb 2025 wells (DNA 2025-02-03) that were
  run alongside the round 2 samples. The Nov 2025 parental is a different DNA
  prep and would shift the baseline ~30%.
- DDX3Y parental wells span -17% to +17%, which is the assay's practical noise
  floor. It is drawn as a grey band; round 1 falls inside it.
- The two assays disagree on absolute magnitude (89% vs 79%). AMELY/AMELX
  consistently understates LOY across this dataset, likely AMELY probe
  cross-reactivity with the AMELX paralog. Prefer DDX3Y/RNaseP or KDM5D/RNaseP
  for absolute LOY.

### `DDX3X_clones_LOY.{png,pdf}`

DDX3X clone panel — 9 edited clones (`DDX3X_LOY1-9`) and 9 sibling controls
(`DDX3X_CTRL1-9`), from `analysis/ddPCR_DDX3X_clones.R`. Input is
`data/LOY_clones.xlsx`, sheet `DDX3X_LOY` — the only figure here that does
**not** come from the master summary table.

Three panels: **A** gel calls, **B** raw ddPCR ratios against their expected
values, **C** ratios normalised to the control clones.

#### The Gel / DDX3X / Call columns are cell fill colours, not text

Those four columns of the sheet are **empty** — the result is recorded as the
cell background colour. `readxl` returns nothing for them, so
`analysis/read_xlsx_cells.R` unzips the `.xlsx` and reads the fills out of the
XML directly (`tidyxl` would do this in one call but is not installed).

| Fill | ARGB | Meaning |
|------|------|---------|
| green | `FF00B050` | present / intact |
| yellow | `FFFFFF00` | partial / ambiguous |
| red | `FFFF0000` | absent / disrupted |
| grey | `FF7F7F7F` | not assessed |

The mapping is fixed by the ddPCR block, where the same four colours sit on
top of numbers whose meaning is unambiguous (green on ~0.5, red on 0).

#### The ddPCR cells are formulas — use the copy numbers

Every ddPCR cell is `=<target copies>/<reference copies>`, e.g. `=147.9/307.5`.
`readxl` returns only the computed ratio and throws the two counts away, so
`read_xlsx_cells()` recovers them from the formula. They matter twice over:
they put an interval on each ratio, and the *reference* count is what explains
every control that reads above one Y per cell.

#### Normalisation: Y dosage per cell, not % LOY against a control median

**This replaced an earlier control-median normalisation that produced negative
% LOY, down to -79%.** Negative loss is not a thing; it meant only "this
control's ratio came out above the median control's".

The figure now divides by the **theoretical expected ratio**, not by a control:

```
dosage = ratio / expected      # 1.0 = one intact Y per cell, 0 = none
LOY%   = (1 - dosage) * 100
```

| Assay | Expected | Why |
|-------|----------|-----|
| RPS4Y1, KDM5D | **0.5** | 1 Y copy vs the 2-copy autosomal RNaseP reference |
| AMELY/AMELX | **1.0** | 1 Y copy vs the single X of a male cell |

Why this is the better denominator:

- **Dosage has a floor at 0**, so the % LOY axis it carries is bounded at 100%.
- **It is absolute.** The scale no longer depends on which control happened to
  land in the middle, so it is comparable across experiments.
- **It attributes the deviation.** A control at 1.8 is *1.8 Y copies per cell* —
  a measurement to explain, not a baseline to normalise other clones against.
- Dosage can still exceed 1, but that now reads as "the assay measured more Y
  than one copy per cell", i.e. over-measurement rather than negative loss.
  Panel C's % LOY axis is deliberately **unlabelled beyond dosage 1.0** to say
  exactly that.

The switch barely moves the LOY clones — the control medians land at dosage
0.98 / 1.02 / 0.97, so the two denominators agree to within ~4% — which is
itself the check that the assays are calibrated. LOY9 goes 80.0→80.3% and
83.0→82.7%. Both scales are in the source CSV (`loy_pct` vs
`loy_pct_vs_ctrl_median`).

#### Confidence intervals

Bars in panel C are 95% intervals on each ratio, from the recovered copy
numbers: for counts `a/b` the binomial proportion `p = a/(a+b)` has an exact
Clopper-Pearson interval, and `ratio = p/(1-p)` is monotone in `p`, so the
interval transforms straight across.

**They are not replicates** — there are none. They are also **conservative**:
the sheet records copies/µL, not raw positive-droplet counts, and droplet
counts are larger. So a well whose bar clears 1.0 is genuinely off; a well
inside its bar is *not* thereby shown to be sound.

#### The negative values were bad reference wells, not noise

Four control measurements have intervals excluding one intact Y — these are
exactly the rows that used to plot as negative % LOY. None is a Y gain; each is
a reference well reading away from its assay's typical level:

| Well | Assay | Target | Reference | Typical ref | Dosage (95% CI) |
|------|-------|--------|-----------|-------------|-----------------|
| CTRL2 | AMELY/AMELX | 129.4 | **71.2** | 146 | 1.82 (1.35–2.46) |
| CTRL7 | AMELY/AMELX | 132 | **89** | 146 | 1.48 (1.13–1.96) |
| CTRL6 | RPS4Y1 | 108 | **153** | 308 | 1.41 (1.09–1.82) |
| CTRL8 | AMELY/AMELX | 223 | **333** | 146 | 0.67 (0.56–0.80) |

A low denominator inflates the ratio. These are wells to repeat, not properties
of the clones — and it is why AMELY/AMELX looked like "the noisiest assay" under
the old normalisation: three of its nine reference wells are off.

#### Results


% LOY, on the absolute (dosage) scale:

| Clone | RPS4Y1 (Yp) | AMELY/AMELX (Yp) | KDM5D (Yq) | Gel |
|-------|-------------|------------------|------------|-----|
| LOY1-5, 7, 8 | 100% | 100% | 100% | Y absent |
| LOY6 | — | — | — | **Y present** |
| LOY9 | 80.3% | 82.7% | 99.9% | not assessed |

- **LOY6 is not an LOY clone.** Its gel shows Y present for both targets and
  the sheet's Call is green; it was dropped before ddPCR. The figure prints
  the reason on the empty row rather than leaving a silent gap.
- **LOY9 is the one clone the assays split on.** Both Yp targets sit near
  0.2 Y copies per cell while Yq is undetectable. That is well outside the
  noise, but it is n=1 with no gel, and it reads equally as a retained
  Yp fragment or a ~20% Y-bearing subpopulation. Needs a repeat plus a gel.
- Every other clone with data is complete LOY in all three assays.

#### Two things in the sheet that need resolving

- **`Call` for CTRL5 is red while all three of its ddPCR assays say Y is
  present** (0.52 / 1.08 / 0.51, dead on the expected values). CTRL5 is also
  the only clone whose DDX3X gel is green, and its band note is
  `single band, WT size` — i.e. unedited at DDX3X. So the red Call is most
  likely "reject this clone, DDX3X is WT" rather than a Y call. That makes
  `Call` mean different things in different rows, since for the LOY clones red
  tracks Y loss. Worth settling before the column is used programmatically.
- **`gel_DDX3X` is red for almost every clone including the controls**, which
  under the colour key above reads as "DDX3X disrupted". That is the intended
  genotype for both arms, so it is probably right — but LOY2 and CTRL1 carry
  the same band note (`single band, smaller than WT`) and are coloured
  differently (yellow vs red), so the column is not fully self-consistent.

**Caveats carried in the figure:**

- One well per clone per assay. **No replicates anywhere.** The bars are
  counting intervals, not replicate error — see above.
- CTRL9 has no RPS4Y1 value at all.
- Most LOY clones have a literal `0` typed in rather than a formula, so there
  are no copy numbers behind them and no interval can be drawn. Their dosage
  is 0 either way; only the upper bound is unknown.
- Control dosage spans RPS4Y1 0.96-1.41, AMELY/AMELX 0.67-1.82, KDM5D
  0.90-1.07. KDM5D is by far the tightest assay of the three — and once the
  three bad AMELX reference wells are set aside, AMELY/AMELX is not
  intrinsically noisy, which is a change from what the old normalisation
  implied.
- Do not compare raw ratios across assays; that is what panel B's two
  different expected lines are there to make obvious. Panel C is the
  comparable one.

### `DDX3X_clones_dosage.{png,pdf}`

Panel C on its own, from `analysis/ddPCR_DDX3X_dosage.R`, with two changes:

- **LOY6 is dropped.** It has no ddPCR data, so on a ddPCR-only figure its row
  is just a gap. It stays in the three-panel figure, where panel A carries its
  gel call and the row therefore says something. The omission is data-driven
  (any clone with no ddPCR data drops out), not a hard-coded `LOY6`, and the
  caption names what was dropped and why.
- **The DDX3X gel call is appended** at the right — call tile plus the band
  pattern it was read from — so Y status and edit status sit on one row.

Clone labels are repeated on the gel panel rather than blanked: it sits across
a gap from the dosage panel, so a reader tracking one clone needs the row named
on both sides.

**This is the clean version — no interval bars, no caption, no subtitle.** Use
`DDX3X_clones_LOY.{png,pdf}` when the caveats need to travel with the figure;
this one is for talks and slides, where they don't fit on the slide anyway.

Two consequences to be aware of when using it:

- **The LOY6 omission is no longer stated on the figure.** It is in this README,
  in the console output of the script, and implicitly in the source CSV. If the
  figure is shown standalone to someone who doesn't know the panel, say it.
- **The interval bars are gone deliberately, not just for tidiness.** The only
  intervals derivable from this sheet come from treating copies/µL as Poisson
  counts, which overstates their width by roughly 4× against real droplet-count
  precision (see the three-panel figure's notes). Drawing them at the wrong
  width invited the reader to conclude "everything else is consistent with one
  intact Y", which is exactly what they do *not* establish. The copy numbers
  stay in the source CSV so proper intervals can be computed if the droplet
  counts are ever recovered from the QuantaSoft exports.

### `DDX3X_clones_grant.{png,pdf}`

Condensed single-panel version for a grant, from
`analysis/ddPCR_DDX3X_grant.R`. **3.4 × 2.5 in at 600 dpi** — sized for a small
single-column slot, no title, no subtitle, no caption. Use the PDF: it is
vector, so it stays sharp if the layout rescales it. Restricted to the clones
that will go into the **Y+** and **LOY** pools.

Type is 6.2–7 pt, which is the floor for print legibility — if the figure gets
shrunk further in the layout, scale the box in the script and re-render rather
than shrinking the placed image, or the labels go under 6 pt.

**Retained (12):** LOY1, LOY3, LOY4, LOY5, LOY7, LOY8 · CTRL1, CTRL3, CTRL4,
CTRL6, CTRL7, CTRL8

**Excluded (6)** — held in the script as a named vector, `EXCLUDE`, so the
figure can never end up with an unexplained clone set; the reasons print on
every run:

| Clone | Reason |
|-------|--------|
| LOY6 | Y detected on gel; never taken to ddPCR |
| LOY9 | partial Yp retention (~0.2 Y per cell) — not cleanly LOY |
| LOY2 | DDX3X gel call ambiguous (single band, smaller than WT) |
| CTRL2 | DDX3X gel call ambiguous (blank band); AMELX reference well low |
| CTRL5 | DDX3X gel band is WT size — not edited at DDX3X |
| CTRL9 | DDX3X and Y gels not assessed; no RPS4Y1 value |

The exclusion list is internally consistent: the four clones dropped for gel
reasons are exactly the ones whose DDX3X call was yellow, green (WT) or grey.
**Every retained clone is DDX3X-disrupted** — which is why there is no DDX3X
tile column, since twelve identical tiles carry one fact. That fact is
**computed, not asserted** (change `EXCLUDE` and it is recomputed) and printed
to the console on every run, because with the subtitle gone the figure no
longer states it. **It needs to go in the figure legend.**

**Three retained wells still scatter away from 1.0**, and a reviewer will see
them. All three are the reference well reading off, not the clone:

| Well | Assay | Dosage | Reference | Typical |
|------|-------|--------|-----------|---------|
| CTRL6 | RPS4Y1 | 1.41 | 153 | 308 |
| CTRL7 | AMELY/AMELX | 1.48 | 89 | 146 |
| CTRL8 | AMELY/AMELX | 0.67 | 333 | 146 |

Each is one assay out of three; the other two agree with 1.0 in every case, so
no clone's Y status is in doubt. If the figure needs to be tidier than this,
the fix is repeating those wells, not dropping the clones.

Note `facet_wrap` is used so the group labels sit horizontally on top — legible
at grant scale where rotated strips are not. That forces **equal panel
heights**, which is only correct while both groups hold the same number of
clones. The script warns if a change to `EXCLUDE` breaks that, and the fix is
`facet_grid(pool ~ ., space = "free_y")`.

### Shared code

`analysis/ddPCR_DDX3X_common.R` holds the loading, the normalisation and the
colour vocabulary for both DDX3X figures, so they cannot drift apart. Each
figure script sources it and does nothing but compose.

`analysis/read_xlsx_cells.R` is the low-level xlsx reader (fills + formulas)
that `common.R` builds on.

**patchwork gotcha, hit twice:** `base_theme` specifies `strip.text.y` and
`axis.text.y`, so a per-panel `theme(strip.text.y = element_blank())` set when
the panel is built gets **overwritten** if `base_theme` is added afterwards in
the compose step — and the facet strip prints in every gutter instead of only
at the far right. Add `base_theme` first, then the per-panel overrides.

## Naming conventions

`analysis/ddPCR_HT1080_CPQ53_rounds.R` defines three reusable canonicalisers,
because the master table accumulated several naming schemes over two years:

- `canon_round()` — `Version` to an ordered Parental/Round 1/Round 2 factor
- `canon_assay()` — target + reference gene to one assay label
- `canon_arm()` — collapses `CP LOY pool`→`CP`, `Q53+CP`→`CP+Q53`,
  `HT1080_P3_REP1_LOY`→`P3_REP1_LOY`, `C4 LOY clone`→`C4`

Ratios are **not comparable across assays**: controls sit near 1.0 for
AMELY/AMELX but near 0.5 for DDX3Y/KDM5D vs RNaseP. Normalise within assay
before pooling.

## Note

The summary CSVs are often open in Excel, which locks them on Windows. The
figure script catches this and writes a `*_new.csv` alongside rather than
aborting after the figures are already written.
