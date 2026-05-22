# 01_genome_engineering/

Generation and validation of isogenic LOY cell lines via CRISPR-Cas9. This
section documents the engineering arc from sgRNA design through clone QC, so
that all downstream analyses (RNA-seq, ATAC-seq, ChIP-seq, screens) rest on a
validated set of LOY+ vs. WT lines.

## Planned analyses

| #   | File                              | Description |
|-----|-----------------------------------|-------------|
| 1.1 | `1.1_sgRNA_guide_creation.Rmd`    | sgRNA design targeting chrY for clean loss; off-target scoring, library curation. TODO. |
| 1.2 | `1.2_ddPCR_validation.Rmd`        | Droplet digital PCR for chrY copy number in candidate clones. TODO. |
| 1.3 | `1.3_lpWGS_validation.Rmd`        | Low-pass WGS to confirm chrY loss genome-wide and screen for off-target copy-number changes. TODO. |

## Existing upstream work (not migrated)

- **1.1**: `a:/Dropbox/LOY_paralogs/projects/sgRNA_guide_development/` — currently empty; start fresh.
- **1.2**: `a:/Dropbox/LOY_paralogs/projects/validation/CTG/` has CTG viability work; ddPCR-specific work TBD.
- **1.3**: `a:/Dropbox/LOY_paralogs/projects/WGS/` — T2T-aligned BAMs for A549/HT1080/LNCaP, mosdepth coverage, chrY/autosome log2 ratio reports, 11 prior `coverage_report*.Rmd` files. The richest existing pipeline of the three.
