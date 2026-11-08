# `R/` — Shared R Helper Functions

Shared R utility scripts used across the repository's `.Rmd` and `.R` pipelines. Module-specific helper libraries are also colocated inside their respective module directories:

- **Public DepMap & PRISM helpers:** [`../analysis/01_depmap_prism_public/scripts/common.R`](../analysis/01_depmap_prism_public/scripts/common.R), [`context_helpers.R`](../analysis/01_depmap_prism_public/scripts/context_helpers.R), [`extension_helpers.R`](../analysis/01_depmap_prism_public/scripts/extension_helpers.R), [`power_helpers.R`](../analysis/01_depmap_prism_public/scripts/power_helpers.R), [`plot_style.R`](../analysis/01_depmap_prism_public/scripts/plot_style.R).
- **Bulk RNA-seq & Cross-Line helpers:** `00_common.R` inside [`../analysis/03_rnaseq_transcriptomics/3.1_HT1080_baseline/R/`](../analysis/03_rnaseq_transcriptomics/3.1_HT1080_baseline/R/), [`3.2_A549_baseline/R/`](../analysis/03_rnaseq_transcriptomics/3.2_A549_baseline/R/), [`3.3_DU145_baseline/R/`](../analysis/03_rnaseq_transcriptomics/3.3_DU145_baseline/R/), and [`3.4_crossline_and_external/R/`](../analysis/03_rnaseq_transcriptomics/3.4_crossline_and_external/R/).
- **ATAC-seq & Multi-Omic helpers:** [`../analysis/04_atacseq_chromatin/code/00_environment.R`](../analysis/04_atacseq_chromatin/code/00_environment.R).
