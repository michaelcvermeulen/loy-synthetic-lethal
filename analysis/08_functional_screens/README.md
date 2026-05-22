# 08_functional_screens/

CRISPR synthetic-lethal and small-molecule screens in the engineered LOY vs
WT isogenic lines, used to identify and prioritize genes and compounds that
selectively impair LOY+ cells. **CRISPR screens were run in both A549 and
HT1080 as independent efforts**; the **drug screen was performed in HT1080
only**. Top CRISPR hits are validated in arrayed follow-up assays (8.4),
then cross-line CRISPR concordance and CRISPR-vs-drug overlap (HT1080) are
examined in the integration notebook (8.5).

## Planned analyses

| #   | File                                       | Description |
|-----|--------------------------------------------|-------------|
| 8.1 | `8.1_crispr_synthetic_lethal_A549.Rmd`     | A549 LOY vs WT CRISPR screen; MAGeCK essentiality, BAGEL, hit annotation. TODO. |
| 8.2 | `8.2_crispr_synthetic_lethal_HT1080.Rmd`   | HT1080 LOY vs WT CRISPR screen (Brunello library); MAGeCK + metascape enrichment, hit annotation. TODO. |
| 8.3 | `8.3_drug_screen_HT1080.Rmd`               | HT1080 small-molecule screen (HT1080 only — no A549 drug screen was performed); normalization, hit calling, dose-response. **Wet-lab and primary analysis by Ben Snider (Craig Lab); re-analysis and integration for this manuscript by Mike Vermeulen.** TODO. |
| 8.4 | `8.4_crispr_hit_validation.Rmd`            | Arrayed validation of top CRISPR hits in LOY vs WT lines: individual sgRNAs / shRNAs / inhibitors, CTG viability + growth curves, dose-response. Confirms screen hits are not artifacts of pooled screen dynamics. TODO. |
| 8.5 | `8.5_screen_integration.Rmd`               | Cross-line CRISPR concordance (A549 vs HT1080), CRISPR-vs-drug overlap in HT1080, validated-hit core set, intersection with DepMap dependencies ([`../05_depmap_pancancer/`](../05_depmap_pancancer/)). TODO. |

## Existing upstream work (not migrated)

- **8.1 — A549 CRISPR screen**: `a:/Dropbox/LOY_paralogs/projects/crispr_screen/a549_loy_screen1/` — MAGeCK + BAGEL pipeline.
- **8.2 — HT1080 CRISPR screen**: `a:/Dropbox/LOY_paralogs/projects/crispr_screen/ht1080_loy_screen_2/` — Brunello library, metascape enrichment, passage-tracking metadata.
  - Combined CRISPR-screen `.Rmd` count across both: 8 files.

- **8.3 — HT1080 drug screen**: `a:/Dropbox/LOY_paralogs/projects/drug_screen/` — 14 `.Rmd` files: optimization (OPT3–5), normalization, set comparisons (SET2/3/4), power simulation, figures. HT1080 only — no A549 drug screen was performed.

- **8.4 — CRISPR hit validation**: `a:/Dropbox/LOY_paralogs/projects/validation/CTG/` — CTG viability data (HT1080 WT vs KO, SREBF1 knockdowns) + `CTG_visualization.Rmd` + rendered figures. Likely seed for the arrayed-validation notebook; expand as additional hits go through validation.

## Open questions

- Is the CRISPR functional screen actually part of section 1 (genome
  engineering → screen on engineered lines) rather than its own section?
- Are the drug screens scoped into the manuscript / thesis chapter, or
  separate? If in-scope, do they belong adjacent to ChIP/ATAC mechanism
  validation or under DepMap/dependency reasoning?
