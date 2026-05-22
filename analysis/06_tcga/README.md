# 06_tcga/

TCGA pan-cancer analysis of LOY: prevalence by tumor type, association with
mutation burden / clinical features, and concordance of LOY-calling methods.

## Planned analyses

| File                              | Description |
|-----------------------------------|-------------|
| `6_tcga.Rmd`                      | LOY calling + QC, differential expression (LOY vs Y-retained), pathway enrichment, mutation enrichment, method concordance. TODO. |

## Existing upstream work (not migrated)

- `a:/Dropbox/LOY_paralogs/projects/depmap_tcga_loy/TCGA/` — 7 `.Rmd` files: LOY calling QC, Boruta feature selection (Qi-based calls), DE, pathway enrichment, method concordance, mutation enrichment. Note potential duplication with the DepMap-level Boruta notebook in [`../05_depmap_pancancer/`](../05_depmap_pancancer/) — worth reconciling at migration time.
