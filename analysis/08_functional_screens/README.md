# 08_functional_screens/

CRISPR synthetic-lethal and small-molecule screens in the engineered LOY vs
WT isogenic lines, used to identify and prioritize genes and compounds that
selectively impair LOY+ cells.

> **Open question (CRISPR screen placement).** The genome-wide CRISPR
> synthetic-lethal screen (8.1) could alternatively live under section 1
> (engineering → screen on engineered lines). Parked here for now; relocate if
> preferred.

## Planned analyses

| #   | File                                | Description |
|-----|-------------------------------------|-------------|
| 8.1 | `8.1_crispr_synthetic_lethal.Rmd`   | Genome-wide Brunello screen in LOY vs WT isogenic lines; MAGeCK essentiality, BAGEL, hit annotation. TODO. |
| 8.2 | `8.2_drug_screen.Rmd`               | Small-molecule screen in LOY vs WT lines; normalization, hit calling, dose-response. **Wet-lab and primary analysis by Ben Snider (Craig Lab); re-analysis and integration for this manuscript by Mike Vermeulen.** TODO. |

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
