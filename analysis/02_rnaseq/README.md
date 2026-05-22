# 02_rnaseq/

Bulk RNA-seq across the three engineered isogenic systems (A549, HT1080,
LNCaP), plus a cross-cell-line integrated analysis. Per-line analyses use a
standardized STAR two-pass + featureCounts pipeline with `DESeq2` for
differential expression and `VIPER` for transcription-factor activity inference.

## Planned analyses

| #   | File                          | Description |
|-----|-------------------------------|-------------|
| 2.1 | `2.1_rnaseq_A549.Rmd`         | A549 LOY vs CTRL DE + TF activity. TODO. |
| 2.2 | `2.2_rnaseq_HT1080.Rmd`       | HT1080 LOY vs CTRL DE + TF activity. TODO. |
| 2.3 | `2.3_rnaseq_LNCaP.Rmd`        | LNCaP LOY vs CTRL DE + TF activity. TODO. |
| 2.4 | `2.4_rnaseq_combined.Rmd`     | Cross-cell-line integration; convergent vs cell-context-specific LOY signatures. TODO. |

## Existing upstream work (not migrated)

- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/A549/` — 5 `.Rmd` (analysis, extended, viper).
- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080/` — 5 `.Rmd` (same standardized pipeline).
- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/LNCAP/` — 5 `.Rmd` (same standardized pipeline).
- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/combined/` — `analysis_combined.Rmd`.
- `a:/Dropbox/LOY_paralogs/projects/RNAseq_LOY/` — older/exploratory outputs and DE tables; audit for any analyses to bring forward vs. supersede.
