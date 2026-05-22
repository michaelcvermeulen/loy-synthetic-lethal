# 08_functional_screens/

> **Flagged for review.** This section captures two bodies of existing work
> from `LOY_paralogs/projects/` that did not appear in the 1–7 numbering you
> sketched: the CRISPR synthetic-lethal screen and the drug screens. They are
> parked here so they aren't lost in the restructure — please relocate, rename,
> split, or drop as appropriate.

## Planned analyses

| #   | File                                | Description |
|-----|-------------------------------------|-------------|
| 8.1 | `8.1_crispr_synthetic_lethal.Rmd`   | Genome-wide Brunello screen in LOY vs WT isogenic lines; MAGeCK essentiality, BAGEL, hit annotation. TODO. |
| 8.2 | `8.2_drug_screen.Rmd`               | Compound screen (likely SAGA/STAGA-pathway focused), normalization, hit calling, dose-response. TODO. |

## Existing upstream work (not migrated)

- **8.1 — CRISPR screen**: `a:/Dropbox/LOY_paralogs/projects/crispr_screen/`
  - Two screens: `a549_loy_screen1/` (MAGeCK + BAGEL) and `ht1080_loy_screen_2/` (Brunello library, metascape enrichment, passage-tracking metadata).
  - 8 `.Rmd` files total.

- **8.2 — Drug screens**: `a:/Dropbox/LOY_paralogs/projects/drug_screen/`
  - 14 `.Rmd` files: optimization (OPT3–5), normalization, set comparisons (SET2/3/4), power simulation, figures.

## Open questions

- Is the CRISPR functional screen actually part of section 1 (genome
  engineering → screen on engineered lines) rather than its own section?
- Are the drug screens scoped into the manuscript / thesis chapter, or
  separate? If in-scope, do they belong adjacent to ChIP/ATAC mechanism
  validation or under DepMap/dependency reasoning?
