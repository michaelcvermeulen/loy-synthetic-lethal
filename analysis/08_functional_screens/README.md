# 08_functional_screens/

CRISPR synthetic-lethal and small-molecule screens in the engineered LOY vs
WT isogenic lines, used to identify and prioritize genes and compounds that
selectively impair LOY+ cells. **A549 and HT1080 are treated as independent
efforts** with their own screen designs, libraries / compound sets, and
analysis pipelines; cross-line concordance is examined in the integration
notebook (8.5).

> **Open question (CRISPR screen placement).** The CRISPR synthetic-lethal
> screens (8.1, 8.2) could alternatively live under section 1 (engineering →
> screen on engineered lines). Parked here for now; relocate if preferred.

## Planned analyses

| #   | File                                       | Description |
|-----|--------------------------------------------|-------------|
| 8.1 | `8.1_crispr_synthetic_lethal_A549.Rmd`     | A549 LOY vs WT CRISPR screen; MAGeCK essentiality, BAGEL, hit annotation. TODO. |
| 8.2 | `8.2_crispr_synthetic_lethal_HT1080.Rmd`   | HT1080 LOY vs WT CRISPR screen (Brunello library); MAGeCK + metascape enrichment, hit annotation. TODO. |
| 8.3 | `8.3_drug_screen_A549.Rmd`                 | A549 small-molecule screen; normalization, hit calling, dose-response. **Wet-lab and primary analysis by Ben Snider (Craig Lab); re-analysis and integration for this manuscript by Mike Vermeulen.** TODO. |
| 8.4 | `8.4_drug_screen_HT1080.Rmd`               | HT1080 small-molecule screen; normalization, hit calling, dose-response. **Wet-lab and primary analysis by Ben Snider (Craig Lab); re-analysis and integration for this manuscript by Mike Vermeulen.** TODO. Confirm whether HT1080 drug screen was performed; drop this row if not. |
| 8.5 | `8.5_screen_integration.Rmd`               | Cross-line concordance: CRISPR vs drug hits across A549 and HT1080; intersection with DepMap dependencies ([`../05_depmap_pancancer/`](../05_depmap_pancancer/)). TODO. |

## Existing upstream work (not migrated)

- **8.1 — A549 CRISPR screen**: `a:/Dropbox/LOY_paralogs/projects/crispr_screen/a549_loy_screen1/` — MAGeCK + BAGEL pipeline.
- **8.2 — HT1080 CRISPR screen**: `a:/Dropbox/LOY_paralogs/projects/crispr_screen/ht1080_loy_screen_2/` — Brunello library, metascape enrichment, passage-tracking metadata.
  - Combined CRISPR-screen `.Rmd` count across both: 8 files.

- **8.3 / 8.4 — Drug screens**: `a:/Dropbox/LOY_paralogs/projects/drug_screen/` — 14 `.Rmd` files: optimization (OPT3–5), normalization, set comparisons (SET2/3/4), power simulation, figures. Cell-line breakdown across A549 vs HT1080 to be confirmed at migration time.

## Open questions

- Is the CRISPR functional screen actually part of section 1 (genome
  engineering → screen on engineered lines) rather than its own section?
- Are the drug screens scoped into the manuscript / thesis chapter, or
  separate? If in-scope, do they belong adjacent to ChIP/ATAC mechanism
  validation or under DepMap/dependency reasoning?
