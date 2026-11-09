# `data/` — Local Input Data Contract & Public Accessions

Raw sequencing reads, large alignment files (`*.bam`, `*.bw`), flow cytometry `.fcs` files, and large serialized R objects (`*.rds`) are not tracked in git. Smaller derived summary tables, sample manifests, and figure source data are included directly within [`../analysis/`](../analysis/).

To rerun the upstream pipelines from raw/intermediate inputs, place or mount the project data directories under `DATA_ROOT` (default: `A:/Dropbox/LOY_paralogs/projects`):

```text
data/ (or upstream DATA_ROOT)
├── depmap_prism/          # DepMap 24Q4/26Q1, Sanger Score, NextGen 3D, PRISM, GDSC1/2, CTD2, BeatAML
├── sgrna_t2t_y/           # T2T-CHM13v2.0 chrY reference, CenSat v2.1 annotations, candidate guide tables
├── ddpcr/                 # QuantaSoft / QX200 droplet counts & concentration ratio spreadsheets
├── lpwgs/                 # 15 T2T-CHM13v2.0 (Y-PAR-masked) BAMs, 100-kb mosdepth & QDNAseq tracks
├── rnaseq/                # 44 in-house STAR v2 (Y-PAR-masked, GENCODE v49) count matrices & external sets
├── atacseq/               # 6 HT1080 BAMs, MACS3 narrowPeak consensus counts, TOBIAS/HOMER tracks
├── crispr_screens/        # HT1080 & A549 Brunello sgRNA count matrices + 27 external isogenic screens
└── comp_assay/            # BD FACSDiscover A8 & CytoFLEX S spectral flow cytometry .fcs files
```

---

## 1. Public Dependency & Pharmacogenomic Datasets (`01_depmap_prism_public/`)

- **DepMap 2D CRISPR & Multi-Omics:** DepMap Public (`24Q4` and `26Q1` releases): `CRISPRGeneEffect.csv` (Chronos), `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, `OmicsCNGene.csv`, `OmicsSomaticMutations.csv`, `Model.csv` (568 classified male models: 404 Y-retaining, 164 LOY).
- **3D Organoid & Neurosphere Screens:** Sanger Project Score and DepMap NextGen 3D CRISPR viability cohorts.
- **Pharmacological Screens:** PRISM Repurposing (`Primary`, `Secondary`, and `PRISM REP.1M`), GDSC1, GDSC2, CTD2, and BeatAML (`Tyner et al.`).
- **TCGA Pan-Cancer Cohort:** Genomic Data Commons (`GDC`) copy-number (`mLRRY`), WES/WGS somatic mutations, and RNA-seq expression across male tumours.

## 2. In-House Low-Pass Whole-Genome Sequencing (`02_genome_engineering/2.3_lpwgs_characterization/`)

- **Libraries (15 total):**
  - **HT1080 (2):** `HT1080-sgCTRL-D0`, `HT1080-LOY-D0` (`CP + Q53 2x`).
  - **A549 (3):** `A549-CTRL2` (`sgCTRL`), `A549-LOY` (`CP + Q53 2x`), `A549-BC14` (clonal reference).
  - **DU145 (10):** `sgCTRL`, `DMSO-LOY-1X` (`CP + Q53 + P47 1x`), `LOY-2X` (`CP + Q53 + P47 2x`), `30` (panobinostat 30 nM), `70` (panobinostat 75 nM; raw label `70`), `KU` (KU-55933), `KU-30`, `KU-70`, plus technical replicates.
- **Reference:** `T2T-CHM13v2.0` with `chrY:1–2,771,479` (`PAR1`) and `chrY:62,122,809–62,460,029` (`PAR2`) hard-masked so all PAR reads map uniquely to `chrX` PAR intervals.

## 3. In-House Bulk RNA-Seq (`03_rnaseq_transcriptomics/` & `05_uty_kdm5d_addback/`)

- **44 Unique In-House Libraries (Illumina NovaSeq 2×150 bp, `STAR v2` Y-PAR-masked alignment to `GRCh38` / `GENCODE v49`, `featureCounts v2.1.1`):**
  - **HT1080 Baseline (`15` libraries across 2 sampling windows):** August 2025 (`3 Parental`, `3 sgCTRL`, `3 LOY`) and May 2026 (`3 sgCTRL`, `3 LOY`).
  - **A549 Baseline (`11` libraries):** `3 Parental`, `4 sgCTRL`, `4 LOY` (with Agilent TapeStation `RINe` covariates).
  - **DU145 Baseline (`6` libraries):** Passage-matched `P12`, `P13`, `P14` pairs of `sgCTRL` and `LOY` (`CP + Q53 + P47 2x`).
  - **HT1080 *UTY* & *KDM5D* Rescue (`18` libraries total; `12` new transduced libraries + `6` May 2026 baseline libraries):** `LOY + Empty Vector` ($n=3$), `LOY + UTY` ($n=3$), `LOY + KDM5D` ($n=3$), and descendant `LOY + UTY + KDM5D` ($n=3$).
- **Frozen Gene-Set Resource:** `MSigDB v2026.1.Hs` (Hallmark, Reactome, GO:BP, KEGG, C2:CGP, C3, C6).

## 4. External Y-Loss, Mosaic & Chromosome-Loss Transcriptomic Datasets (`03_rnaseq_transcriptomics/`)

- **18 External Y-Loss / Mosaic / 45,X / Y-Paralog KO Contrasts (`3.4_cross_line_synthesis/`):**
  - Constitutional karyotype & iPSC/fibroblast/LCL series: `GSE126712` (45,X vs. 46,XY), Page Lab sex-chromosome dosage series (LCL & fibroblast Y dosage), `GSE264742` (*UTY*-KO hESCs).
  - Engineered & spontaneous cancer LOY models: `GSE301297` (A549 clonal LOY RNA-seq & multiome pseudobulk), `GSE229233` (MB49 mouse bladder LOY), `GSE152001`, `GSE190507` (*DDX3X*/*DDX3Y* KO).
  - Naturally Y-mosaic single-cell cohorts: `GSE313572` (HT1080 3D spheroids control sample) and **Tahoe-100M** DMSO-treated mosaic cell lines passing the inclusion rule (`MOSAIC` call and median PAR $\log_2\mathrm{FC} \le -0.5$: `A549`, `J82`, `SK-MEL-2`, `HCT15`).
- **16 Chromosome- and Arm-Loss Contrasts (`3.4_cross_line_synthesis/`):**
  - `GSE130306` (RPTEC acute CRISPR Yq and 3p arm loss across days 5, 11, 18, and 90), `GSE227499` (ReDACT trisomy-corrected disomic vs. trisomic lines), and `GSE190497` (RPE1 *TP53*-null autosomal monosomies).

## 5. In-House ATAC-Seq (`04_atacseq_chromatin/`)

- **HT1080 ATAC-Seq (`6` libraries):** `CTRL_11`, `CTRL_12`, `CTRL_21` and `LOY_11`, `LOY_12`, `LOY_21` sampled at passages 11, 12, and 21 from the established HT1080 `sgCTRL` and `LOY` populations (`105,647` consensus autosomal/non-PAR X peaks).
- **Chromatin Reference Overlaps:** ENCODE ChIP-seq peak catalogues (`EZH2`, `SUZ12`, `H3K27me3`, `H3K4me3`) and `JASPAR 2022/2024` vertebrate CORE motifs.

## 6. In-House & External Isogenic CRISPR Screens (`06_crispr_screens/`)

- **HT1080 Genome-Wide Brunello Screen (`ht1080_loy_screen_2`):** `76,441` targeting sgRNAs + `1,000` NTCs across `D0` (`sgCTRL`, `LOY`), `D12` (`sgCTRL` R1/R2, `LOY` R1/R2), and `D21` (`sgCTRL` R1/R2, `LOY` R1/R2).
- **A549 Genome-Wide Brunello Screen (`a549_loy_screen1`):** Baseline (`P3` / reconstructed `P4`) and endpoint (`P18`) libraries in matched `sgCTRL` and `LOY` A549 populations.
- **27 External Isogenic CRISPR Screens (`6.3_external_isogenic_benchmarks/`):**
  - **Desjardins et al. (2026) (`Cell Reports` / `bioRxiv`):** 17 Brunello isogenic screens (`Chr13q_FT282`, `Chr18q_RPE1`, `KEAP1_A549`, `ARID1A_BEAS2B`, `KMT2D_BEAS2B`, `STK11_RPE1`, `FBXW7_RPE1`, `CDK12_RPE1`, `IDH1_RPE1`, `TET2_RPE1`, `CCNE1_FT282`, `DNMT3A_K562`, `SF3B1_K562`, `SRSF2_K562`, `U2AF1_K562`, `PKMYT1i_COLhTERT`, `PKMYT1i_DLD1_FNLS_FBXW7`).
  - **Feng et al. (2022) (`Sci Adv`, Figshare `19398332`):** 10 isogenic tumour-suppressor knockout screens (`PTEN`, `RB1`, `CDKN2A`, `TP53`, `KEAP1`, `NF1`, `NF2`, `BAP1`, `ARID1A`, `SMAD4`).
  - **BioGRID ORCS (`v2.0.18`):** Homo sapiens CRISPR screen compendium.

## 7. Spectral Flow-Cytometry Competitive Growth Assays (`07_competitive_growth_assays/`)

- **Calibration & Gating Runs:** `CytoFLEX S` pilot runs (`AUG_13_2026`, `QC2_20260821`) and `BD FACSDiscover A8` spectral unmixing runs (`20260901_Mike_Run1`, `20260904_Mike_Run2`, 9-point reciprocal `eGFP`/`mCherry` mixing series).
- **Targeted RNP Knockout Competition Runs:** `D0` (`20261015`), `D4` (`20261020`), and follow-up timepoints across `NTC`, `OR2J2` cutting control, `DDX3X`, `TADA2B`, `TAF6L`, and `SREBF1` in both reciprocal colour pairings (`LOY-eGFP + CTRL-mCherry` and `LOY-mCherry + CTRL-eGFP`) plus same-genotype controls.
