# 04_atacseq_multiomic/

Bulk ATAC-seq in isogenic LOY lines, plus integration with the matched RNA-seq
from [`../02_rnaseq/`](../02_rnaseq/) to identify chromatin-accessibility
changes that coincide with LOY-driven expression changes.

## Planned analyses

| File                                | Description |
|-------------------------------------|-------------|
| `4_atacseq_multiomic.Rmd`           | Differential accessibility (LOY vs CTRL), peak annotation, motif enrichment, ATAC/RNA integration. TODO. |

## Existing upstream work (not migrated)

- `a:/Dropbox/LOY_paralogs/projects/ATACseq/` — aligned BAMs, MACS2 narrow/broad peaks, FastQC, ataqv QC, bigWig tracks for LOY-11/12/21 vs CTRL replicates. No prior ATAC/RNA integration `.Rmd` found.
