# `06_crispr_screens/` — Isogenic Genome-Wide Brunello CRISPR–Cas9 Screens & Benchmarks

**Thesis Reference:** Chapter 4, Section 4.2.2–4.2.4, Section 4.3.1, Section 4.3.5, Section 4.3.6, and Supplementary Methods S4.2–S4.4 (`Fig 4.2–4.4`, `Fig 4.7–4.8`, `Supp Fig 4.5–4.11`, `Supp Fig 4.13–4.16`).

This module contains the processing, quality-control audits, differential-dependency scoring (`drugZ`, `MAGeCK` RRA/MLE, and baseline-adjusted `IsogenicZ`), network/multi-omic integrations, and 27-screen external isogenic benchmarks for the genome-wide Human CRISPR Knockout Pooled Library (**Brunello**, 76,441 targeting guides across 19,114 genes + 1,000 non-targeting controls) in matched Y-depleted (`LOY`) and control (`sgCTRL`) **HT1080** and **A549** populations.

---

## Sub-Modules

```text
06_crispr_screens/
├── README.md
├── 6.1_HT1080_brunello/                     # Ch 4 Section 4.3.1: HT1080 screen (D0, D12, D21), SAGA/MED12/SREBF, multi-omic integration
├── 6.2_A549_brunello/                       # Ch 4 Section 4.3.5: A549 screen (P3/P4 to P18), EIF1AX/L3MBTL2/UBR5, trimmed-P3 sensitivity
└── 6.3_external_isogenic_benchmarks/        # Ch 4 Section 4.3.6 & Section 4.2.3: 27 external isogenic screens & cross-screen guide QC audits
```

---

## Headline Findings by Sub-Module

### 6.1 HT1080 Genome-Wide Brunello Screen & Multi-Omic Context ([`6.1_HT1080_brunello/`](6.1_HT1080_brunello/))

- **Screen Architecture & Quality Control (`Fig 4.2`, `Supp Fig 4.5`):**
  - Matched `sgCTRL` and `LOY` (`CP + Q53 2x`) HT1080 populations were transduced at low MOI, selected with puromycin to `D0`, and passaged in duplicate (`R1`, `R2`) through `D12` (~8–9 population doublings) and `D21` (~15–16 doublings), with `ddPCR` confirming stable maintenance of Y depletion throughout outgrowth (`Supp Fig 4.5e`).
  - Core essential genes strongly deplete from `D0` to `D12`/`D21` in both arms ($\mathrm{AUROC} > 0.94$; `Fig 4.2c`).
  - **D0 PCR-jackpot diagnosis & scoring hierarchy (`Supp Fig 4.6–4.7`):** Because `D0` baseline libraries contained 173 recurring PCR-jackpotted guides (shared with A549 `P3` at $\rho = 0.997$), primary candidate discovery uses unpaired endpoint **`drugZ`** (`D12 + D21` `LOY` vs. `sgCTRL`), corroborated by **`MAGeCK` RRA** and outlier-masked, baseline-adjusted **`IsogenicZ`** (adapted from Schukken et al., 2026).
- **Recovery of *DDX3X* & Non-Gametolog Candidate Vulnerabilities (`Fig 4.3–4.4`, `Supp Fig 4.7–4.10`):**
  - ***DDX3X*** is the #1 genome-wide synthetic lethal sensitizer in Y-depleted HT1080 cells ($\mathrm{normZ} = -11.45$, $\mathrm{FDR} = 2.13 \times 10^{-26}$; `Fig 4.3a,b`, `Fig 4.4a`), whereas Y-linked guides behave as neutral controls in `LOY` (`Fig 4.3c`).
  - Coordinated non-gametolog sensitizers emerge across four functional modules (`Fig 4.3–4.4`, `Supp Fig 4.8`, `Supp Fig 4.10`):
    1. **SAGA Chromatin Coactivator Complex (6 of 20 subunits):** ***TAF6L*** ($\mathrm{normZ} = -4.87$, $\mathrm{FDR} = 0.0067$), ***TADA2B*** ($\mathrm{normZ} = -4.29$, $\mathrm{FDR} = 0.046$), ***SUPT20H*** ($\mathrm{normZ} = -4.19$, $\mathrm{FDR} = 0.056$), ***ATXN7L3*** ($\mathrm{normZ} = -3.91$, $\mathrm{FDR} = 0.125$), ***TADA1*** ($\mathrm{normZ} = -3.77$, $\mathrm{FDR} = 0.165$), and exploratory ***TAF5L*** (`Fig 4.4b–d`).
    2. **Mediator & Chromatin Regulators:** ***MED12*** ($\mathrm{normZ} = -5.06$, $\mathrm{FDR} = 0.0040$) and X-linked ***BCORL1***.
    3. **SREBP Lipid & Selenoprotein Metabolism:** ***SREBF2***, ***SREBF1***, ***SCLY***, and ***ACSL4***.
    4. **Mitochondrial Import & ECM Adhesion:** X-linked ***TIMM8A*** ($\mathrm{normZ} = -5.73$, $\mathrm{FDR} = 2.02 \times 10^{-4}$), ***LLGL1***, and ***ITGB5***.
- **Multi-Omic Exclusion of Self-Silencing (`molecular_integration/`, `Supp Fig 4.11`):**
  - Integrating Chapter 3 HT1080 RNA-seq, promoter ATAC-seq, and *UTY*/*KDM5D* rescue profiles proves that differential CRISPR dependency does **not** correlate with baseline transcriptional change ($\rho = -0.005$) or promoter closing of the candidate genes themselves (`Supp Fig 4.11a–c`). Instead, SAGA/Mediator/SREBF candidates remain stably expressed and promoter-open while their downstream target programmes (`MYC` targets, sterol metabolism, ECM) are transcriptionally remodelled.

### 6.2 A549 Genome-Wide Brunello Screen & Rules of Gametolog Recovery ([`6.2_A549_brunello/`](6.2_A549_brunello/))

- **Screen Design & Sensitivity Branches (`Fig 4.7–4.8`, `Supp Fig 4.13–4.14`):**
  - Evaluates `P18` endpoint libraries (`LOY` vs. `sgCTRL`) alongside reconstructed `P4` and outlier-trimmed `P3` baseline-adjusted `IsogenicZ` sensitivity models (accounting for a verified `CTRL_P3_R1` sample swap that was genotyped as `LOY` by Y-guide/ddPCR forensics; `qc_docs/`).
- **Candidate Discovery (`Fig 4.8`):**
  - Recovers the X-linked translation-initiation gametolog ***EIF1AX*** (rank 8, $\mathrm{normZ} = -4.09$) alongside the PRC1.6 subunit ***L3MBTL2*** and the HECT E3 ligase ***UBR5*** (which physically regulates SAGA/chromatin complexes).
- **Rules Governing Which X–Y Gametolog Scores (`Supp Fig 4.7`, `Supp Fig 4.9`):**
  - Parental `HT1080` is unusually tolerant of *DDX3X* knockout in DepMap (`Chronos +0.18`), giving a massive dynamic range when *DDX3Y* is lost (`normZ = -11.45`), whereas pan-essential *EIF1AX* (`Chronos -1.48`) and *RPS4X* (`Chronos -1.89`, further buffered by residual distal `Yp` *RPS4Y1* in HT1080) deplete to the sequencing count floor during pre-`D0` puromycin selection (`Supp Fig 4.7d`).
  - Conversely, parental `A549` already has low baseline Y expression and strong baseline *DDX3X* dependency (`Chronos -0.98`), compressing the *DDX3X* differential window while allowing *EIF1AX* sensitization to score at `P18`.

### 6.3 External Isogenic Screen Benchmarks & Technical QC Audits ([`6.3_external_isogenic_benchmarks/`](6.3_external_isogenic_benchmarks/))

- **Specificity Across 27 External Isogenic CRISPR Screens (`isogenicz_27_screens/`, `Supp Fig 4.15`):**
  - Uniformly reprocessing **27 published isogenic CRISPR contrasts** (17 from Desjardins et al., 2026 including `chr13q` and `chr18q` arm losses; 10 from Feng et al., 2022) with both `drugZ` and `IsogenicZ` proves that ***DDX3X* sensitization is strictly unique to LOY (0/27 external contrasts)**, whereas **SAGA complex sensitization recurs in specific non-LOY contexts** including `chr13q` arm loss in `FT282` cells and *KEAP1* knockout in `A549` cells (`Supp Fig 4.15`).
- **Cross-Screen Guide Non-Detection, PCR Jackpots & Baseline Adjustment (`guide_artefact_audit/`, `Supp Fig 4.6`, `Supp Fig 4.16`):**
  - Identifies 618 Brunello guides undetected at `D0`/`P3` across both in-house screens and external Brunello datasets, driven by `5'-TCTC` and `5'-AGACG` vector-overlap sequence prefixes that get stripped during 5′ adapter trimming (`Supp Fig 4.6a–c`).
  - Demonstrates across 29 isogenic screens (`2` in-house + `27` external) that `T0`/`D0` baseline arm divergence quantitatively predicts the rank divergence between endpoint-only `drugZ` and baseline-adjusted `IsogenicZ` ($\rho = 0.61$; `Supp Fig 4.16`).
