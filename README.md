# Molecular Consequences and Synthetic Lethal Dependencies of Loss of the Y Chromosome in Cancer

**Author:** Michael Cory Vermeulen  
**Supervisors:** Dr. Tomas Babak & Dr. Andrew W. Craig  
**Institution:** Queen's University, Kingston, Ontario, Canada (2026)  
**Chapter 2 Publication:** Vermeulen, M., Craig, A. W. & Babak, T. *npj Systems Biology and Applications* **11**, 143 (2025). [https://doi.org/10.1038/s41540-025-00618-7](https://doi.org/10.1038/s41540-025-00618-7)

## Repository Structure

```text
loy-synthetic-lethality/
├── R/                                         # Shared R utility functions
├── data/                                      # Local raw/processed data paths (gitignored)
├── figures/
│   ├── chapter3/code/                         # Chapter 3 figure-generation scripts
│   └── chapter4/code/                         # Chapter 4 figure-generation scripts
└── analysis/
    ├── 01_depmap_prism_public/                # DepMap CRISPR & PRISM drug reanalysis + power simulations
    ├── 02_genome_engineering/                 # T2T-Y sgRNA design, ddPCR & Y-PAR-masked lpWGS
    │   ├── 2.1_guide_design/
    │   ├── 2.2_ddpcr_validation/
    │   └── 2.3_lpwgs_characterization/
    ├── 03_rnaseq_transcriptomics/             # Bulk RNA-seq & cross-dataset comparisons
    │   ├── 3.1_HT1080_baseline/
    │   ├── 3.2_A549_baseline/
    │   ├── 3.3_DU145_baseline/
    │   ├── 3.4_cross_line_synthesis/
    │   └── 3.5_metascape_enrichment/
    ├── 04_atacseq_chromatin/                  # HT1080 ATAC-seq, motif/footprinting & RNA/ATAC integration
    ├── 05_uty_kdm5d_addback/                  # HT1080 lentiviral UTY & KDM5D rescue RNA-seq
    ├── 06_crispr_screens/                     # Genome-wide Brunello CRISPR screens
    │   ├── 6.1_HT1080_brunello/
    │   ├── 6.2_A549_brunello/
    │   └── 6.3_external_isogenic_benchmarks/
    └── 07_competitive_growth_assays/          # Two-colour spectral flow-cytometry competition assays
```
