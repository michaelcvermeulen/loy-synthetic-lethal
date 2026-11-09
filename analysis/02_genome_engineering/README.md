# `02_genome_engineering/` — Repeat-Targeting Guide Design & Isogenic Model Characterization

**Thesis Reference:** Chapter 3, Section 3.2.1–3.2.2, Section 3.3.1–3.3.2, and Supplementary Methods S3.1 (`Fig 3.1–3.6`, `Supp Fig 3.4–3.14`).

This module covers the *in silico* discovery of chromosome-Y repeat-targeting sgRNAs on the `T2T-CHM13v2.0` reference, transient Cas9 ribonucleoprotein (RNP) engineering of pooled non-clonal Y-depleted populations across human cancer lineages, multiplex droplet digital PCR (`ddPCR`) validation, pharmacological modulation of Y depletion, and `T2T-CHM13v2.0` (Y-PAR-masked) low-pass whole-genome sequencing (`lpWGS`).

---

## Sub-Modules

```text
02_genome_engineering/
├── README.md
├── 2.1_guide_design/                        # Ch 3 Section 3.3.1: T2T-Y sliding-window scan, off-target filtering & 5-guide annotation
├── 2.2_ddpcr_validation/                    # Ch 3 Section 3.3.2: Mixing linearity, HT1080/A549/DU145/HepG2 ddPCR & panobinostat/KU-55933
└── 2.3_lpwgs_characterization/              # Ch 3 Section 3.3.2: T2T-Y Y-PAR-masked lpWGS, 100-kb QDNAseq/mosdepth & regional Y atlases
```

---

## Headline Findings by Sub-Module

### 2.1 *In Silico* T2T-Y Repeat Guide Discovery ([`2.1_guide_design/`](2.1_guide_design/))

- **Workflow (`Supp Fig 3.4`, `Fig 3.1`):**
  1. Exhaustive sliding-window scan of `T2T-CHM13v2.0` (HG002) `chrY` for 20-nt protospacers adjacent to an `NGG` PAM yielded a **14,723-guide candidate table**.
  2. Strict specificity filtering requiring zero exact matches and zero 1-mismatch off-targets on autosomes or `chrX` retained **6,901 Y-specific guides**, ranked into a **215-entry high-copy top-hits table**.
  3. Five lead guides spanning three major Y repeat classes (`CenSat v2.1`) were selected for experimental testing (`Fig 3.1`, all 737 exact `T2T-Y` target coordinates verified):
     - **`CP`** (`5'-TGGABE...` / centromeric `DYZ3` $\alpha$-satellite): 35 exact `T2T-Y` sites.
     - **`CN`** (centromeric `DYZ3` $\alpha$-satellite): 35 exact `T2T-Y` sites.
     - **`Q53`** (`Yq12` `DYZ1` / `HSat3` satellite array): 259 exact `T2T-Y` sites.
     - **`P47`** (`Yp` *TSPY* tandem repeat array): 46 exact `T2T-Y` sites.
     - **`P51`** (`Yp` *TSPY* tandem repeat array): 34 exact `T2T-Y` sites.

### 2.2 Multiplex ddPCR Validation & Pharmacological Modulation ([`2.2_ddpcr_validation/`](2.2_ddpcr_validation/))

- **Cell-mixture linearity check (`Supp Fig 3.7`):** Titration mixtures of Y-positive `H2452` and Y-negative cells across 50 ng and 100 ng total reaction DNA inputs confirm linear recovery of `AMELY`/`AMELX` concentration ratios across the 0–100% range.
- **Reference-copy scaling (`Fig 3.2`):** Main ddPCR figures report estimated Y-linked target copies per cell scaled by the baseline reference locus copy number (`HT1080`: 1 `AMELX`, 2 `RNaseP`; `A549`: 2 `AMELX`, 3 `RNaseP`).
- **Guide screening & successive RNP nucleofections (`Fig 3.2`, `Supp Fig 3.11`):**
  - In `HT1080`, initial RT-qPCR screening (`Fig 3.2a`) and multiplex ddPCR (`Fig 3.2b,c`) show that `CP` and `Q53` strongly outperform `CN`, `P47`, and `P51`—demonstrating that reference repeat copy number alone does not predict cleavage efficacy (`CP` vs. `CN`, both 35 sites in `DYZ3`).
  - Delivering combined **`CP + Q53` Cas9 RNPs across two successive nucleofections (`2x`)** reduces Y target copies by **88–99%** in bulk pooled `HT1080` and `A549` populations without single-cell cloning (`Fig 3.2b–e`).
  - `HepG2` hepatocellular carcinoma populations assessed by `AMELY`/`AMELX` and `DDX3Y`/`RNaseP` ddPCR also show substantial Y copy-number reduction after RNP targeting (Section 3.3.2.4).
- **Pharmacological modulation in DU145 (`Supp Fig 3.13`):** Pre-treating `DU145` cells for 24 h prior to a single `CP + Q53 + P47` nucleofection with the pan-HDAC inhibitor **panobinostat** (30 nM or 75 nM [raw provenance label `70`]) increases Y depletion across `AMELY`, `DDX3Y`, and `KDM5D` compared with vehicle (`DMSO`), whereas post-nucleofection ATM inhibition with **KU-55933** (16 h) does not improve depletion and partially blunts the panobinostat effect.

### 2.3 Low-Pass WGS & Regional Y-Chromosome Atlases ([`2.3_lpwgs_characterization/`](2.3_lpwgs_characterization/))

- **Y-PAR-masked `T2T-CHM13v2.0` alignment & 100-kb bin profiling (`Fig 3.3–3.6`, `Supp Fig 3.5–3.6`, `Supp Fig 3.8–3.10`, `Supp Fig 3.12`, `Supp Fig 3.14`):**
  - Masking `PAR1` and `PAR2` on `chrY` during `BWA-MEM` alignment routes all pseudoautosomal reads to `chrX` PAR bins while preserving unique `MSY` coverage across 100-kb windows (`regional_y_atlases_ypar_20261028/`).
  - **HT1080 (`Fig 3.3`, `Fig 3.4`):** `HT1080-LOY-D0` (`CP + Q53 2x`) exhibits complete loss of euchromatic `Yq` (<1% RNA retention across *USP9Y*–*EIF1AY*) and strong proximal `Yp` reduction, with ~20–36% residual expression of distal `Yp` genes (*RPS4Y1*, *ZFY*) and ~60–69% `PAR1`/`PAR2` expression retention. Autosomal copy-number comparison identifies focal subclonal differences on `chr13q` and `chr4` that are explicitly tracked in downstream RNA-seq and ATAC-seq sensitivity analyses (`Supp Fig 3.8`).
  - **A549 (`Fig 3.5`, `Supp Fig 3.10`):** `A549-LOY` (`CP + Q53 2x`) shows uniform whole-chromosome Y loss across both `Yp` and `Yq` (~10–12% median DNA retention across `Yp` and `Yq` and near-complete loss of all MSY transcripts).
  - **DU145 (`Fig 3.6`, `Supp Fig 3.11–3.12`, `Supp Fig 3.14`):** In `DU145` (`CP + Q53 + P47 2x`), `lpWGS` and RNA-seq reveal complete elimination of euchromatic `Yq` (*USP9Y*–*EIF1AY*) and 96–98% depletion of proximal `Yp` including *AMELY* (~6.5 Mb), but **46–54% DNA and RNA retention of a distal male-specific `Yp` block (~2.46–4.80 Mb)** containing *RPS4Y1* (~2.54 Mb) and *ZFY* (~2.93 Mb). This finding motivated the addition of the *RPS4Y1* ddPCR assay and proves that standard *AMELY*/*AMELX* assays can misclassify partial Y deletions as complete LOY.
