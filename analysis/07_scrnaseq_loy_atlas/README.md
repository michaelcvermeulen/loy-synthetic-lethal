# 07_scrnaseq_loy_atlas/

Single-cell LOY atlas across cancer types. Starts with prostate; additional
cancer types to follow. Per-dataset LOY status calling at the cell level, then
cross-sample / cross-cell-type characterization of LOY-associated programs.

## Planned analyses

| File                                | Description |
|-------------------------------------|-------------|
| `7_scrnaseq_loy_atlas.Rmd`          | Combined atlas: LOY calling, integration across datasets, cell-type composition, LOY-associated programs, comparison to bulk RNA-seq signatures from [`../02_rnaseq/`](../02_rnaseq/). TODO. |

Likely to split into per-cancer notebooks (`7.1_prostate.Rmd`, `7.2_<next>.Rmd`,
...) as additional tumor types come in; revisit numbering then.

## Existing upstream work (not migrated)

- `a:/Dropbox/LOY_paralogs/projects/cancer_atlas/` — multi-project Seurat objects (e.g., ccRCC Krishna et al., several GSE datasets), LOY status calling + UMAPs per sample / cell-type / LOY status, per-project processed `.rds` objects, master metadata TSV.
