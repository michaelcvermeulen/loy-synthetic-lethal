# 05_depmap_pancancer/

DepMap pan-cancer dependency analysis: which gene knockouts selectively impair
LOY+ cell lines, supporting the SAGA/STAGA addiction model.

## Planned analyses

| File                              | Description |
|-----------------------------------|-------------|
| `5_depmap_pancancer.Rmd`          | LOY-status classification of DepMap lines, Boruta feature selection, ssGSEA pathway analysis, SAGA-complex-focused dependency profile, celligner integration. TODO. |

## Existing upstream work (not migrated)

- `a:/Dropbox/LOY_paralogs/projects/depmap_tcga_loy/` — substantial existing pipeline:
  - `01_Boruta_DepMap.Rmd` (feature selection),
  - `04_SAGA_Complex_Analysis.Rmd` (the closest existing mechanism writeup),
  - ssGSEA / DE / fGSEA notebooks,
  - celligner integration.
  - 10 `.Rmd` files total at the DepMap level (TCGA-specific files live under `TCGA/`, see [`../06_tcga/`](../06_tcga/)).
