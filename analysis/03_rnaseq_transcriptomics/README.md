# `03_rnaseq_transcriptomics/` — Bulk RNA-Seq & Comparative Transcriptomics

**Thesis Reference:** Chapter 3, §3.2.3, §3.3.3, and Supplementary Methods S3.2–S3.3 (`Fig 3.7–3.8`, `Fig 3.11–3.13`, `Supp Fig 3.1–3.2`, `Supp Fig 3.15–3.20`, `Supp Fig 3.29–3.41`).

This module contains the finalized thesis bulk RNA-seq pipelines across the three engineered isogenic LOY models (**HT1080**, **A549**, and **DU145**), the cross-cell-line and 18-external-dataset comparative synthesis, and the 16-contrast chromosome/arm-loss benchmark.

---

## Unified Pipeline Standards

- **Alignment & Quantification:** All 44 in-house libraries (Illumina NovaSeq 2×150 bp; [`molecular_samples_inventory/`](molecular_samples_inventory/)) were aligned with `STAR v2` two-pass mode against a **Y-PAR-masked `GRCh38` reference** and quantified with `featureCounts v2.1.1` against **`GENCODE v49`** so that `PAR1`/`PAR2` reads are retained on `chrX` PAR features without multimapper loss (`Supp Fig 3.1–3.2`).
- **Normalisation & Differential Expression:** `TMM` library-size normalisation (`edgeR`) on expressed autosomal and non-PAR X genes, followed by `limma-voom` linear modelling. `MSY` and `PAR` regional retention are evaluated in dedicated sex-chromosome branches outside the primary autosomal/non-PAR X multiple-testing family.
- **Pathway & Network Analysis:** Effect-size/moderated-$t$ `fgsea` and fixed-correlation **`CAMERA`** (default inter-gene correlation $\rho = 0.01$; estimated-correlation `CAMERA` retained only as historical sensitivity) against the frozen **`MSigDB v2026.1.Hs`** collections (primary family: `Hallmark`, `Reactome`, `GO:BP`, 4,868 sets), complemented by `STRING v12` physical/functional PPI clustering and `Metascape` term networks.

---

## Sub-Modules

```text
03_rnaseq_transcriptomics/
├── README.md
├── molecular_samples_inventory/             # Manifest & QC metrics across all 44 in-house RNA-seq libraries
├── 3.1_HT1080_baseline/                     # Ch 3 §3.3.3.1: HT1080 two-window reanalysis (SEPT_11_2026_REANALYZE)
├── 3.2_A549_baseline/                       # Ch 3 §3.3.3.2: A549 RINe-adjusted reanalysis (SEPT_18_2026_REANALYZE)
├── 3.3_DU145_baseline/                      # Ch 3 §3.3.3.3: DU145 passage-paired reanalysis (SEPT_23_2026_REANALYZE)
├── 3.4_crossline_and_external/              # Ch 3 §3.3.3.4: Cross-line synthesis & 18 external Y-loss/mosaic datasets
└── 3.5_chromosome_arm_loss/                 # Ch 3 §3.3.3.5: Comparison with 16 chromosome/arm-loss contrasts
```

---

## Summary of Findings by Sub-Module

### 3.1 HT1080 Baseline RNA-Seq ([`3.1_HT1080_baseline/`](3.1_HT1080_baseline/))

- **Design (`Fig 3.7`, `Supp Fig 3.15`):** 12 core libraries from the established `sgCTRL` and `LOY` (`CP + Q53 2x`) populations sampled across two windows (**August 2025** and **May 2026**; $n=3$ per condition per window) plus 3 August `Parental` libraries (25,526 expressed autosomal/non-PAR X genes).
- **Sex-chromosome dosage (`Fig 3.7b,f`):** Euchromatic `Yq` transcripts (*USP9Y*–*EIF1AY*) are reduced by >99%; distal `Yp` transcripts (*RPS4Y1*, *ZFY*) retain ~20–36%; `PAR1` retains **68.8% / 66.8%** and `PAR2` retains **60.2% / 63.6%** of `sgCTRL` expression in August / May.
- **Differential expression (`Fig 3.7d,e,g–l`):** The primary `~ period + condition` model identifies **1,003 high-confidence non-Y/non-PAR DEGs** ($\mathrm{FDR} < 0.05$, $|\log_2\mathrm{FC}| \ge 1$, matching sign in both windows: **551 higher**, **452 lower**). Top upregulated genes include *UNC5A*, *SDK2*, and *ROR2*; top downregulated genes include *MARCKS*, *P2RX5*, and *CDH12*.
- **Pathways & STRING modules (`Fig 3.8`, `Supp Fig 3.16–3.20`):** `fgsea`, fixed-correlation `CAMERA` (17 discoveries at $\mathrm{FDR} < 0.05$ across 4,868 primary sets), and 42 `STRING` modules reveal strong downregulation of **interferon-$\alpha$/$\gamma$, cytokine/chemokine, complement, and contractile/ECM programmes** alongside upregulation of **axon guidance / synaptic adhesion (`UNC5`–`FLRT`–latrophilin module) and sterol/cholesterol biosynthesis**.

### 3.2 A549 Baseline RNA-Seq ([`3.2_A549_baseline/`](3.2_A549_baseline/))

- **Design (`Supp Fig 3.29–3.31`):** 8 core libraries (`4 sgCTRL`, `4 LOY`) plus 3 `Parental` libraries. Because `PC1` associates with sample RNA integrity (`RINe`), the primary nuclear-gene model adjusts for `RINe` (`~ RINe + condition`) across **22,843 expressed nuclear genes** (including `chrY` and `PAR`), with the unadjusted condition-only model retained as sensitivity (`primary_nuclear_RINe/`).
- **Differential expression:** Identifies **27 $\mathrm{FDR} < 0.05$ genes**: **17 MSY genes** (uniformly depleted) and **10 non-Y/non-PAR genes** (8 higher, 2 lower), with pathway rankings showing trend-level reductions in interferon/inflammatory and ribosomal gene sets (`Supp Fig 3.29c`).

### 3.3 DU145 Baseline RNA-Seq ([`3.3_DU145_baseline/`](3.3_DU145_baseline/))

- **Lineage verification & provenance note:** Raw files for this experiment carry the historical prefix `LNCAP_`. Both STR profiling and transcriptomic identity checks (`05_identity.R`: zero reads for *AR*, *KLK3*, *KLK2*, *FOLH1*, *HOXB13*; highest genome-wide DepMap expression correlation with `DU145`) confirmed that the line is **DU145** androgen-independent prostate carcinoma.
- **Design & DE (`Fig 3.11`, `Supp Fig 3.32–3.33`):** 6 passage-matched libraries (`P12`, `P13`, `P14`, paired `~ passage + treatment`). Identifies **136 non-Y/non-PAR DEGs at $\mathrm{FDR} < 0.05$** and **114 highlighted passage-consistent genes** ($\mathrm{FDR} < 0.05$, $|\log_2\mathrm{FC}| \ge 0.5$, concordant direction in all 3 passage pairs: **83 lower**, **31 higher**).
- **Regional Y retention & biological programmes (`Fig 3.6`, `Fig 3.11–3.12`, `Supp Fig 3.34–3.39`):** While `Yq` (*USP9Y*–*EIF1AY*) is completely lost, distal `Yp` (*RPS4Y1*, *ZFY*) is retained at ~46–54%. Autosomal rewiring (verified after excluding genes on copy-number-altered `chr20p` and `chrX`/`15q`) features strong downregulation of the epithelial master regulator ***GRHL2***, luminal epithelial markers, and **interferon/inflammatory signalling**, paired with upregulation of **basement-membrane collagens (*PXDN*, *COL4A5/6*), ECM organization, and replication-dependent histones** (`chr6p` without local DNA gain; `Supp Fig 3.39`).

### 3.4 Cross-Cell-Line & 18 External Y-Loss Comparisons ([`3.4_crossline_and_external/`](3.4_crossline_and_external/))

- **Universal absence of X-paralog compensation (`Fig 3.13a,b`):** Across `HT1080`, `A549`, and `DU145`, `PAR1` expression drops by 23–33%, `Yq` transcripts are lost, and **none of the X-linked gametologs (*DDX3X*, *EIF1AX*, *RPS4X*, *ZFX*, *KDM6A*, *KDM5C*) show compensatory transcriptional upregulation**.
- **Gene-level divergence vs. pathway-level interferon suppression (`Fig 3.13c,d`, `Supp Fig 3.40`):** Genome-wide autosomal fold-change correlations between the three lines are near zero ($\rho = -0.06$ to $0.11$). However, effect-ranked `fgsea` across our three lines and **18 external Y-loss, naturally Y-mosaic single-cell (`GSE313572`, Tahoe-100M DMSO lines `A549`, `J82`, `SK-MEL-2`, `HCT15`), MB49, 45,X, and *UTY*-KO contrasts** reveals recurrent suppression of **Hallmark Interferon Alpha/Gamma Response and Inflammatory Response** alongside context-dependent shifts in Polycomb (`BMI1`/`MEL18` knockdown reference sets) and ECM programmes.

### 3.5 Comparison with 16 Chromosome- and Arm-Loss Models ([`3.5_chromosome_arm_loss/`](3.5_chromosome_arm_loss/))

- **Separating Y-loss from general monosomy stress (`Supp Fig 3.41`):** Comparing our isogenic LOY models against **16 public chromosome/arm-loss contrasts** (`GSE130306` RPTEC acute CRISPR `Yq` vs. `3p` loss across days 5, 11, 18, and 90; `GSE227499` ReDACT trisomy corrections; `GSE190497` RPE1 monosomies) demonstrates that **cytoplasmic ribosomal-protein suppression** is a shared hallmark of autosomal arm/chromosome loss, whereas acute `Yq` and `3p` loss share early `p53` activation and `E2F` suppression that invert by day 90 of culture.
