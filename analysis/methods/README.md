# methods/

Per-assay / per-technology methods writeups. Each file documents the wet-lab
or computational protocol used so the corresponding analysis `.Rmd` can stay
focused on results and interpretation. These are the documents reviewers will
look at first when asking "how was this generated?".

Survey of existing repo content turned up no dedicated methods files — methods
are currently embedded inline in analysis `.Rmd`s. These pull them out into
canonical, citable references.

## Planned methods writeups

| File                                  | Topic |
|---------------------------------------|-------|
| `m01_genome_engineering.Rmd`          | CRISPR-Cas9 generation of isogenic LOY lines (transfection, selection, clone picking). TODO. |
| `m02_sgRNA_design.Rmd`                | sgRNA design pipeline for chrY targeting; off-target scoring. TODO. |
| `m03_ddPCR.Rmd`                       | ddPCR protocol for chrY copy-number validation. TODO. |
| `m04_lpWGS.Rmd`                       | Low-pass WGS library prep, T2T alignment, mosdepth coverage analysis. TODO. |
| `m05_bulk_rnaseq.Rmd`                 | Bulk RNA-seq library prep + STAR two-pass / featureCounts / DESeq2 + VIPER. TODO. |
| `m06_chipseq.Rmd`                     | ChIP-seq sample prep + computational pipeline. TODO. |
| `m07_atacseq.Rmd`                     | ATAC-seq sample prep + MACS2 / ataqv QC pipeline. TODO. |
| `m08_scrnaseq.Rmd`                    | scRNA-seq processing pipeline (Seurat, integration, LOY calling at cell resolution). TODO. |
| `m09_crispr_functional_screen.Rmd`    | Brunello screen design, transduction, MAGeCK + BAGEL analysis. TODO. |
| `m10_drug_screen.Rmd`                 | Compound screen design, normalization, hit calling. TODO. |
| `m11_depmap_processing.Rmd`           | DepMap data acquisition, LOY classification, Boruta feature selection. TODO. |
| `m12_tcga_processing.Rmd`             | TCGA cohort assembly, LOY calling methods + concordance. TODO. |
