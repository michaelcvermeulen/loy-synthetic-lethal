# `04_atacseq_chromatin/` — HT1080 Chromatin Accessibility & Multi-Omic Integration

**Thesis Reference:** Chapter 3, §3.2.3.3, §3.3.3.1, and Supplementary Methods S3.3.1 (`Fig 3.9–3.10`, `Supp Fig 3.3`, `Supp Fig 3.21–3.28`).

This module contains the complete reproducible pipeline (`SEPT_08_2026_REPROCESS`) for bulk ATAC-seq in Y-depleted (`LOY`) versus control (`sgCTRL`) HT1080 fibrosarcoma cells, transcription-factor motif and footprinting analyses (`chromVAR`, `HOMER`, `TOBIAS`), and multi-omic integration with the matched HT1080 baseline RNA-seq (`SEPT_11_2026_REANALYZE`), `T2T-CHM13v2.0` / `GRCh38` lpWGS copy number, and ENCODE chromatin reference catalogues.

---

## Directory Structure

```text
04_atacseq_chromatin/
├── README.md
├── PROCESSING_PLAN.md                       # End-to-end processing architecture & QC thresholds
├── TF_EXTENSION_PLAN.md                     # chromVAR, HOMER & TOBIAS footprinting specification
├── RUN_STATUS.md                            # Pipeline execution verification log
├── code/                                    # Numbered R & Python scripts (00_environment.R – 89_transcript_TSS_promoter_sensitivity.R)
├── config/                                  # Sample manifest, primary & passage-adjusted designs, PAR coordinates
├── qc/                                      # Fragment size distributions, TSS enrichment, FRiP & correlation tables
├── motifs/                                  # chromVAR deviations, HOMER de novo/known motifs & TOBIAS summaries
└── integration/                             # RNA/ATAC promoter concordance, quadrant ORA & matched PRC2/H3K27me3 tests
```

---

## Experimental Design & Pipeline Summary

- **Sampling Design (`Supp Fig 3.3`):** Six libraries (`CTRL_11`, `CTRL_12`, `CTRL_21` and `LOY_11`, `LOY_12`, `LOY_21`) collected at passages 11, 12, and 21 from the established `sgCTRL` and `LOY` (`CP + Q53 2x`) HT1080 populations. These represent repeated passage samples from one established population per condition.
- **Peak Calling & Normalisation:** Adapter-trimmed, deduplicated non-mitochondrial alignments (`MAPQ >= 30`, ENCODE v2 blacklist filtered) were called with `MACS3` (`--nomodel --shift -100 --extsize 200`) and merged into **105,647 consensus autosomal and non-PAR X peaks** (with `chrY` and `PAR` analysed in a dedicated branch; `Supp Fig 3.21`). Primary differential accessibility uses background-bin `TMM` normalisation and `limma-voom` (`~ condition`, with `~ passage + condition`, `csaw` 150-bp sliding windows, and trended/loess normalisation as sensitivity analyses; `code/05_differential.R`, `code/07_csaw.R`, `code/88_trended_normalisation_sensitivity.R`).

---

## Headline Findings

1. **Widespread Chromatin Compaction Independent of Autosomal Copy Number (`Fig 3.9`, `Supp Fig 3.21–3.23`):**
   - At $\mathrm{FDR} < 0.05$, Y depletion drives **758 differentially accessible regions (DARs)** across autosomal/non-PAR X peaks, with an **~7:1 skew toward chromatin closing** (**663 closed vs. 95 opened peaks**; `Fig 3.9a,b`).
   - Across `chrY`, `MSY` peaks are completely lost in `LOY` while `PAR1`/`PAR2` accessibility drops in proportion to copy number (`Supp Fig 3.21`).
   - Autosomal accessibility changes show near-zero correlation with 100-kb `lpWGS` DNA copy-number log$_2$ ratios ($\rho = -0.05$), and excluding subclonal copy-number-altered regions (`chr13q`, `chr4`) leaves the closing skew intact (`Supp Fig 3.23`).
2. **Transcription-Factor Motif & Footprinting Remodelling (`Supp Fig 3.24–3.27`):**
   - `chromVAR` (`JASPAR 2022/2024`), covariate-matched `HOMER` motif enrichment (`code/08_motifs.R`, `code/84_homer_all_significant.py`), and `TOBIAS` footprinting (`code/23_run_tobias.py`–`29_footprint_aggregate.py`) demonstrate reduced accessibility and footprint depth at **AP-1 (`FOS`/`JUN`/`BATF`), `CTF`/`NFIC`, and `IRF`** motifs within closing distal enhancers and promoters, contrasted with relative enrichment of **GC-rich `KLF`/`SP`, `ETS`, `TFAP2`, and `Forkhead` (`FOXO`/`FOXA`)** motifs in opening or promoter-proximal regions.
3. **Concordant Promoter Closing, Gene Repression & Polycomb Target Enrichment (`Fig 3.10`, `Supp Fig 3.28`):**
   - Integrating promoter ATAC-seq with the two-window HT1080 baseline RNA-seq (`code/80_integrate_thesis_RNA.R`) shows strong directional concordance between promoter accessibility and transcript abundance (`Fig 3.10a–c`).
   - High-confidence downregulated DEGs are **>30-fold more likely to have a significantly closed promoter peak** than non-DE genes, and aggregate TSS accessibility drops by **21%** even across downregulated DEGs without a discrete $\mathrm{FDR} < 0.05$ peak call (`Fig 3.10d`).
   - Quadrant over-representation (`code/81_quadrant_ORA.R`; `Supp Fig 3.28`) and covariate-matched permutation tests (`code/82_prc2_promoter_check.R`, `code/83_polycomb_matched.R`; matched on baseline accessibility, GC content, and CpG density) demonstrate that **jointly closed and transcriptionally repressed promoters are significantly enriched for canonical Polycomb (`PRC2`: `EZH2`, `SUZ12`) and `H3K27me3` targets** (`Fig 3.10e,f`).
