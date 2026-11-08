# Molecular Consequences and Synthetic Lethal Dependencies of Loss of the Y Chromosome in Cancer

This repository contains the computational pipelines, statistical analyses, multi-omic integrations, and figure-generation code supporting the PhD thesis **_Molecular Consequences and Synthetic Lethal Dependencies of Loss of the Y Chromosome in Cancer_** (Queen's University, December 2026) and its constituent manuscripts:

- **Chapter 2 (Published):** Vermeulen, M., Craig, A. W. & Babak, T. Challenges and opportunities for oncology drug repurposing informed by synthetic lethality. *npj Systems Biology and Applications* **11**, 143 (2025). https://doi.org/10.1038/s41540-025-00618-7
- **Chapter 3 (Manuscript in preparation):** *Generation and Molecular Characterization of Isogenic Loss of the Y Chromosome Cancer Models* (Vermeulen, M., Truesdell, P., Conseil, G., Snider, B., Young, S., Craig, A. W. & Babak, T.)
- **Chapter 4 (Manuscript in preparation):** *Functional Genomics and Synthetic Lethal Vulnerabilities in Loss of the Y Chromosome Cancer Models* (Vermeulen, M., Truesdell, P., Conseil, G., Snider, B., Young, S., Craig, A. W. & Babak, T.)

**Read first** ([`docs/`](docs/)):
[Background](docs/background.html) ·
[Condensed summary](docs/summary.html) ·
[Project & thesis overview](docs/project_plan.html)

---

## Scientific Overview

Somatic loss of the Y chromosome (LOY) is one of the most frequent whole-chromosome alterations in male cancers and is associated with inferior clinical outcomes and altered anti-tumour immunity. However, the direct cell-intrinsic transcriptional, chromatin, and synthetic lethal consequences of losing chromosome Y—separated from inter-tumoural genetic heterogeneity and single-cell clonal bottlenecks—have remained poorly resolved. Across Chapters 2–4, this work establishes:

1. **Pan-cancer synthetic lethality and public cohort resolution limits ([`01_depmap_prism_public/`](analysis/01_depmap_prism_public/)):** Integrating DepMap CRISPR screens (568 classified male models), TCGA mutual exclusivity, Sanger/DepMap NextGen 3D organoid and neurosphere screens, and seven pharmacogenomic drug datasets (6,413 pan-cancer and 495,079 subset-stratified tests) confirms robust pan-cancer LOY synthetic lethality for the X-linked gametologs *DDX3X*, *EIF1AX*, *RPS4X*, and *ZFX* (with *DDX3X*, *EIF1AX*, and *RPS4X* replicating in 3D NextGen models). Empirical power simulations across 25 benchmark synthetic lethalities demonstrate that public panels are statistically underpowered once stratified by lineage or driver alteration, motivating matched isogenic LOY models.
2. **T2T-Y-informed, clone-free pooled LOY engineering and partial Yp retention ([`02_genome_engineering/`](analysis/02_genome_engineering/)):** Exhaustive scanning of the complete `T2T-CHM13v2.0` (HG002) Y assembly identified high-copy repeat guides targeting centromeric `DYZ3` $\alpha$-satellite (`CP`, `CN`), `Yq12` `DYZ1`/`HSat3` (`Q53`), and `Yp` *TSPY* (`P47`, `P51`). Transient Cas9 RNP delivery (`CP + Q53`) across two successive nucleofections—or after transient HDAC inhibition with panobinostat—achieves 88–99% Yq depletion directly in bulk pooled populations across HT1080 fibrosarcoma, A549 lung adenocarcinoma, and DU145 prostate carcinoma cells without single-cell cloning. Multiplex ddPCR and Y-PAR-masked T2T-Y low-pass WGS (`lpWGS`) reveal that while A549 loses the entire Y chromosome, HT1080 and DU145 eliminate euchromatic Yq (*USP9Y*–*EIF1AY*) and proximal Yp (*AMELY*) while partially retaining a distal male-specific Yp segment (~2.46–4.80 Mb) encoding *RPS4Y1* and *ZFY*—exposing a critical blind spot of single-locus *AMELY*/*AMELX* LOY assays.
3. **Lineage-specific autosomal reprogramming vs. conserved cell-intrinsic interferon suppression ([`03_rnaseq_transcriptomics/`](analysis/03_rnaseq_transcriptomics/)):** Across all three engineered models, Y depletion reduces summed `PAR1` expression by 23–33% without compensatory transcriptional upregulation of X-linked gametologs (*DDX3X*, *EIF1AX*, *RPS4X*, *ZFX*, *KDM6A*, *KDM5C*). Autosomal differential expression varies strongly by lineage (1,003 high-confidence non-Y/non-PAR DEGs in HT1080; 136 FDR < 0.05 and 114 passage-consistent DEGs in DU145 marked by *GRHL2* suppression and basement-membrane ECM induction; 10 non-Y/non-PAR DEGs in A549), with near-zero pairwise gene-level fold-change correlations ($\rho = -0.06$ to $0.11$). Despite gene-level divergence, pathway analyses across our models and 18 external Y-loss, naturally Y-mosaic single-cell, and 45,X contrasts reveal recurrent, cell-intrinsic suppression of type I/II interferon and inflammatory signalling, while comparisons across 16 chromosome/arm-loss models separate Y-loss responses from generic monosomy stress.
4. **Chromatin compaction at Polycomb targets and causal *UTY*-dependent reversibility ([`04_atacseq_chromatin/`](analysis/04_atacseq_chromatin/) & [`05_uty_kdm5d_addback/`](analysis/05_uty_kdm5d_addback/)):** In HT1080, ATAC-seq reveals an ~7:1 skew toward chromatin closing (663 closed vs. 95 opened peaks at $\mathrm{FDR} < 0.05$), independent of autosomal copy number ($\rho = -0.05$). Downregulated genes are >30-fold more likely to lose promoter accessibility than non-DE genes, and jointly closed/repressed promoters are significantly enriched for PRC2 (*EZH2*/*SUZ12*) and H3K27me3 targets after matching for baseline accessibility, GC content, and CpG density. Lentiviral re-expression of *UTY* (alone or with *KDM5D*), but not *KDM5D* alone, reverses 51% (and 35%) of the 150-gene HT1080 LOY transcriptional signature and increases expression at 36 of 39 jointly closed/repressed genes.
5. **Isogenic genome-wide CRISPR screens and reciprocal spectral competition validation ([`06_crispr_screens/`](analysis/06_crispr_screens/) & [`07_competitive_growth_assays/`](analysis/07_competitive_growth_assays/)):** Genome-wide Brunello CRISPR–Cas9 knockout screens in matched Y-depleted and control HT1080 and A549 populations establish practical rules governing X–Y gametolog recovery (*DDX3X* as the top genome-wide sensitizer in HT1080 at $\mathrm{normZ} = -11.45$, $\mathrm{FDR} = 2.13 \times 10^{-26}$; *EIF1AX* in A549 at rank 8) and nominate context-specific non-gametolog vulnerabilities including six of 20 subunits of the SAGA coactivator complex (*TAF6L*, *TADA2B*, *SUPT20H*, *ATXN7L3*, *TADA1*, *TAF5L*), Mediator (*MED12*), sterol/selenoprotein metabolism (*SREBF2*, *SREBF1*, *SCLY*, *ACSL4*), *TIMM8A*, *L3MBTL2*, and *UBR5*. Multi-omic integration shows these dependencies are not caused by baseline transcriptional silencing or promoter closure of the candidate genes themselves ($\rho = -0.005$). Benchmarking against 27 external isogenic CRISPR screens shows that *DDX3X* sensitization is unique to LOY (0/27 external contrasts), whereas SAGA dependency recurs in select aneuploid and oncogenic contexts (e.g., 13q arm loss, *KEAP1* mutation). Finally, a reciprocal two-colour `eGFP`/`mCherry` competitive growth assay on the `BD FACSDiscover A8` spectral flow cytometer confirms early LOY-selective fitness penalties after *DDX3X* or *TADA2B* knockout across both reporter orientations.

---

## Repository Structure

```text
loy-synthetic-lethality/
├── README.md
├── LICENSE                                    # MIT
├── .gitignore
├── renv.lock                                  # Pinned R package versions (R 4.5.2)
├── .Rprofile                                  # Activates renv on session start
├── analysis/
│   ├── README.md                              # Index of narrative & numbered analysis modules
│   ├── background.Rmd                         # Thesis background, X–Y gametologs & hypotheses
│   ├── summary.Rmd                            # Condensed headline findings across Chapters 3 & 4
│   ├── project_plan.Rmd                       # Thesis structure, chapter roadmap & completion status
│   ├── 01_depmap_prism_public/                # Ch 4 §4.3.0: DepMap CRISPR, NextGen 3D, PRISM/GDSC/CTD2, power sims
│   ├── 02_genome_engineering/                 # Ch 3 §3.3.1–3.3.2: T2T-Y guide design, multiplex ddPCR, Y-PAR-masked lpWGS
│   ├── 03_rnaseq_transcriptomics/             # Ch 3 §3.3.3: HT1080, A549, DU145 RNA-seq, 18 external LOY & 16 arm-loss sets
│   ├── 04_atacseq_chromatin/                  # Ch 3 §3.3.3.1: HT1080 ATAC-seq, chromVAR/HOMER/TOBIAS motifs, RNA/ATAC integration
│   ├── 05_uty_kdm5d_addback/                  # Ch 3 §3.3.4: HT1080 lentiviral UTY & KDM5D rescue RNA-seq & promoter recovery
│   ├── 06_crispr_screens/                     # Ch 4 §4.3.1, 4.3.5, 4.3.6: HT1080 & A549 Brunello screens, 27-screen benchmark
│   ├── 07_competitive_growth_assays/          # Ch 4 §4.3.2–4.3.4: BD FACSDiscover A8 spectral flow competition assays
│   ├── methods/                               # Per-assay wet-lab & computational methods (Ch 3 §3.2/S3 & Ch 4 §4.2/S4)
│   └── mechanism/                             # Mechanistic synthesis & evidence boundaries across Chapters 3 & 4
├── figures/
│   ├── README.md                              # Figure conventions & numbering overview
│   ├── chapter3/                              # Ch 3 figure-generation scripts (code/), captions & numbering register
│   └── chapter4/                              # Ch 4 figure-generation scripts (code/) & numbering register
├── R/                                         # Shared R helper utilities
├── data/                                      # Local-only input data contract (see data/README.md)
└── docs/                                      # Knitted HTML reports & GitHub Pages landing page
```

---

## Numbered Analysis Modules

| # | Module | Thesis Section | Sub-analyses & Contents |
|---|--------|----------------|-------------------------|
| 1 | [DepMap & PRISM Public Screens](analysis/01_depmap_prism_public/) | Ch 4 §4.3.0 (links Ch 2) | **1.1** Y-status classification (568 male models), RF/Boruta, 3D NextGen organoid/neurosphere replication · **1.2** Context-stratified genetic searches (86 driver contexts) & 25-benchmark power simulations · **1.3** Pan-cancer (6,413) & subset-stratified (495,079) drug screens (PRISM, GDSC, CTD2) |
| 2 | [Genome Engineering & Structural Genomics](analysis/02_genome_engineering/) | Ch 3 §3.3.1–3.3.2 | **2.1** T2T-Y repeat guide discovery (`CP`, `CN`, `Q53`, `P47`, `P51`) · **2.2** Multiplex ddPCR (*AMELY*/*AMELX*, *DDX3Y*/*RNaseP*, *KDM5D*/*RNaseP*, *RPS4Y1*), mixing linearity & panobinostat/KU-55933 modulation · **2.3** T2T-Y Y-PAR-masked lpWGS & regional Y atlases (HT1080, A549, DU145) |
| 3 | [Bulk RNA-Seq & Comparative Transcriptomics](analysis/03_rnaseq_transcriptomics/) | Ch 3 §3.3.3 | **3.1** HT1080 baseline (`SEPT_11_2026_REANALYZE`, 2 windows, 1,003 DEGs, STRING modules) · **3.2** A549 baseline (`SEPT_18_2026_REANALYZE`, RINe-adjusted nuclear model) · **3.3** DU145 baseline (`SEPT_23_2026_REANALYZE`, passage-paired P12–P14, *GRHL2*/ECM) · **3.4** Cross-line & 18 external Y-loss/mosaic datasets · **3.5** 16 chromosome/arm-loss contrasts |
| 4 | [ATAC-Seq & Multi-Omic Integration](analysis/04_atacseq_chromatin/) | Ch 3 §3.3.3.1 | **4.1** HT1080 ATAC-seq QC & differential accessibility (`SEPT_08_2026_REPROCESS`, 663 closed vs. 95 opened peaks) · **4.2** `chromVAR`, HOMER & TOBIAS footprinting · **4.3** Promoter RNA/ATAC integration & matched PRC2/H3K27me3 enrichment |
| 5 | [*UTY* & *KDM5D* Add-Back Rescue](analysis/05_uty_kdm5d_addback/) | Ch 3 §3.3.4 | **5.1** 18-library rescue design & construct verification (`THESIS_ADDBACK_2026`) · **5.2** Partial reversal of 150-gene (51% / 35%) and frozen-100-gene LOY signatures by *UTY* (`fry` directional tests) · **5.3** Recovery stratified by baseline ATAC-seq promoter class (36/39 joint loss genes) |
| 6 | [Isogenic Brunello CRISPR Screens](analysis/06_crispr_screens/) | Ch 4 §4.3.1, 4.3.5, 4.3.6 | **6.1** HT1080 Brunello screen (`drugZ`, `MAGeCK`, `IsogenicZ`; *DDX3X*, SAGA, *MED12*, *SREBF1/2*, *TIMM8A*) · **6.2** A549 Brunello screen (*EIF1AX*, *L3MBTL2*, *UBR5*) · **6.3** Multi-omic & STRING/DepMap network integration · **6.4** 27 external isogenic screens & guide-artefact (`TCTC`/`AGACG`, PCR jackpot) audits |
| 7 | [Competitive Growth Validation](analysis/07_competitive_growth_assays/) | Ch 4 §4.3.2–4.3.4 | **7.1** Reciprocal `eGFP`/`mCherry` spectral unmixing (`BD FACSDiscover A8`), autofluorescence gating & 9-point mixing calibration · **7.2** Cas9 RNP knockout validation of *DDX3X*, *TADA2B*, *TAF6L*, and *SREBF1* (D0-to-D4 reciprocal log-odds fitness shifts) |

---

## Thesis Figure-to-File Map

Complete numbering registers (including all 44 Chapter 3 supplementary figures `Supp Fig 3.1–3.44` and all 16 Chapter 4 supplementary figures `Supp Fig 4.1–4.16`) are maintained in:
- [`figures/chapter3/chapter3_figure_numbering.md`](figures/chapter3/chapter3_figure_numbering.md) and [`figures/chapter3/Chapter3_figure_captions.md`](figures/chapter3/Chapter3_figure_captions.md)
- [`figures/chapter4/chapter4_figure_numbering.md`](figures/chapter4/chapter4_figure_numbering.md)

### Chapter 3 Main Figures (`Fig 3.1–3.15`)

| Thesis Figure | Section | Title | Renderer / Script | Analysis Module |
|---|---|---|---|---|
| **Fig 3.1** | §3.3.1 | Five Y-targeting guides and their repeat regions | [`make_five_guide_target_map.py`](figures/chapter3/code/make_five_guide_target_map.py) | [`02_genome_engineering/2.1_guide_design`](analysis/02_genome_engineering/2.1_guide_design/) |
| **Fig 3.2** | §3.3.2.1 | Guide selection and Y-linked DNA depletion in HT1080 and A549 | [`make_ht1080_guide_screening_ddpcr.py`](figures/chapter3/code/) | [`02_genome_engineering/2.2_ddpcr_validation`](analysis/02_genome_engineering/2.2_ddpcr_validation/) |
| **Fig 3.3** | §3.3.2.1 | Y-targeting reduces Y-chromosome DNA and expression in HT1080 | [`make_ht1080_lpwgs_genomewide.py`](figures/chapter3/code/) | [`02_genome_engineering/2.3_lpwgs_characterization`](analysis/02_genome_engineering/2.3_lpwgs_characterization/) |
| **Fig 3.4** | §3.3.2.1 | Residual Yp DNA and expression persist in Y-depleted HT1080 cells | [`make_regional_y_atlas.py`](figures/chapter3/code/make_regional_y_atlas.py) | [`02_genome_engineering/2.3_lpwgs_characterization`](analysis/02_genome_engineering/2.3_lpwgs_characterization/) |
| **Fig 3.5** | §3.3.2.2 | Regional Y-chromosome depletion and expression in A549 | [`make_regional_y_atlas.py`](figures/chapter3/code/make_regional_y_atlas.py) | [`02_genome_engineering/2.3_lpwgs_characterization`](analysis/02_genome_engineering/2.3_lpwgs_characterization/) |
| **Fig 3.6** | §3.3.2.3 | Regional Y-chromosome depletion and expression in DU145 | [`make_regional_y_atlas.py`](figures/chapter3/code/make_regional_y_atlas.py) | [`02_genome_engineering/2.3_lpwgs_characterization`](analysis/02_genome_engineering/2.3_lpwgs_characterization/) |
| **Fig 3.7** | §3.3.3.1 | HT1080 sampling design, global expression structure and expression differences | [`make_ht1080_rna_combined.py`](figures/chapter3/code/make_ht1080_rna_combined.py) | [`03_rnaseq_transcriptomics/3.1_HT1080_baseline`](analysis/03_rnaseq_transcriptomics/3.1_HT1080_baseline/) |
| **Fig 3.8** | §3.3.3.1 | Transcriptional programmes associated with Y depletion in HT1080 | [`make_ht1080_combined_gsea.py`](figures/chapter3/code/make_ht1080_combined_gsea.py) | [`03_rnaseq_transcriptomics/3.1_HT1080_baseline`](analysis/03_rnaseq_transcriptomics/3.1_HT1080_baseline/) |
| **Fig 3.9** | §3.3.3.1 | Y depletion is associated with widespread reductions in chromatin accessibility in HT1080 cells | [`atac_volcano.py`](figures/chapter3/code/atac_volcano.py) | [`04_atacseq_chromatin`](analysis/04_atacseq_chromatin/) |
| **Fig 3.10** | §3.3.3.1 | Reduced gene expression accompanies promoter closing in Y-depleted HT1080 cells | [`make_ht1080_rna_atac_integration.py`](figures/chapter3/code/) | [`04_atacseq_chromatin`](analysis/04_atacseq_chromatin/) |
| **Fig 3.11** | §3.3.3.3 | Transcriptional variation and expression differences in DU145 cells | [`make_du145_global_expression_combined.py`](figures/chapter3/code/make_du145_global_expression_combined.py) | [`03_rnaseq_transcriptomics/3.3_DU145_baseline`](analysis/03_rnaseq_transcriptomics/3.3_DU145_baseline/) |
| **Fig 3.12** | §3.3.3.3 | Transcriptional programs associated with Y loss in DU145 cells | [`du145_programmes_matched_main.py`](figures/chapter3/code/du145_programmes_matched_main.py) | [`03_rnaseq_transcriptomics/3.3_DU145_baseline`](analysis/03_rnaseq_transcriptomics/3.3_DU145_baseline/) |
| **Fig 3.13** | §3.3.3.4 | Shared and background-dependent expression patterns after Y depletion | [`make_crossline_concise.py`](figures/chapter3/code/make_crossline_concise.py) | [`03_rnaseq_transcriptomics/3.4_crossline_and_external`](analysis/03_rnaseq_transcriptomics/3.4_crossline_and_external/) |
| **Fig 3.14** | §3.3.4 | *UTY* and *KDM5D* rescue in Y-depleted HT1080 cells | [`export_addback_design_pdf.py`](figures/chapter3/code/export_addback_design_pdf.py) | [`05_uty_kdm5d_addback`](analysis/05_uty_kdm5d_addback/) |
| **Fig 3.15** | §3.3.4 | *UTY*-containing populations show partial recovery of the HT1080 LOY expression signature | [`render_addback_signature_final.py`](figures/chapter3/code/render_addback_signature_final.py) | [`05_uty_kdm5d_addback`](analysis/05_uty_kdm5d_addback/) |

### Chapter 4 Main Figures (`Fig 4.1–4.8`)

| Thesis Figure | Section | Title | Renderer / Script | Analysis Module |
|---|---|---|---|---|
| **Fig 4.1** | §4.3.0.1 | DepMap reanalysis recovers four established LOY-associated dependencies | [`ch03_main_depmap_y_status_dependencies.py`](figures/chapter3/code/ch03_main_depmap_y_status_dependencies.py) | [`01_depmap_prism_public`](analysis/01_depmap_prism_public/) |
| **Fig 4.2** | §4.3.1.1 | Genome-wide CRISPR screening of Y-depleted and control HT1080 populations | [`make_ht1080_overview_qc.py`](figures/chapter4/code/make_ht1080_overview_qc.py) | [`06_crispr_screens/6.1_HT1080_brunello`](analysis/06_crispr_screens/6.1_HT1080_brunello/) |
| **Fig 4.3** | §4.3.1.2 | The HT1080 screen recovers *DDX3X* and nominates additional differential dependencies | [`make_ht1080_candidates.py`](figures/chapter4/code/make_ht1080_candidates.py) | [`06_crispr_screens/6.1_HT1080_brunello`](analysis/06_crispr_screens/6.1_HT1080_brunello/) |
| **Fig 4.4** | §4.3.1.2 | Endpoint guide support and SAGA complex organization in the HT1080 screen | [`make_ht1080_support_compact.py`](figures/chapter4/code/make_ht1080_support_compact.py) | [`06_crispr_screens/6.1_HT1080_brunello`](analysis/06_crispr_screens/6.1_HT1080_brunello/) |
| **Fig 4.5** | §4.3.2 | Competition-assay workflow, measurement QC and population characterization | [`make_competition_workflow_qc.py`](figures/chapter4/code/) | [`07_competitive_growth_assays`](analysis/07_competitive_growth_assays/) |
| **Fig 4.6** | §4.3.3–4.3.4 | Competitive growth of LOY and CTRL cells after knockout of *DDX3X*, *SREBF1* or *TADA2B* | [`make_competition_hits_draft.py`](figures/chapter4/code/) | [`07_competitive_growth_assays`](analysis/07_competitive_growth_assays/) |
| **Fig 4.7** | §4.3.5 | A549 screen design and quality control | [`make_a549_matched_qc.py`](figures/chapter4/code/) | [`06_crispr_screens/6.2_A549_brunello`](analysis/06_crispr_screens/6.2_A549_brunello/) |
| **Fig 4.8** | §4.3.5 | The A549 screen recovers *EIF1AX* and nominates additional differential dependencies | [`make_a549_candidate_discovery.py`](figures/chapter4/code/) | [`06_crispr_screens/6.2_A549_brunello`](analysis/06_crispr_screens/6.2_A549_brunello/) |

---

## Software Environment

R package versions are pinned with [renv](https://rstudio.github.io/renv/) (`R 4.5.2`). After cloning:

```r
# From the repository root in R:
renv::restore()
```

Python figure-rendering and bioinformatic pipelines use Python 3.11+ (`numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `pyBigWig`, `pysam`, `flowkit`).

---

## Data Availability

Raw sequencing reads, large alignment files (`*.bam`, `*.bw`), flow cytometry `.fcs` files, and large serialized objects (`*.rds`) are not tracked in git. The local directory contract and public dataset accessions are documented in [`data/README.md`](data/README.md):

- **In-house Bulk RNA-seq (44 libraries across HT1080 baseline, A549 baseline, DU145 baseline, and HT1080 *UTY*/*KDM5D* rescue):** GEO accession pending manuscript submission (`STAR v2` Y-PAR-masked alignments against `GENCODE v49`).
- **In-house ATAC-seq (6 HT1080 libraries across passages 11, 12, 21):** GEO accession pending manuscript submission.
- **In-house Low-Pass WGS (15 libraries across HT1080, A549, and DU145):** Aligned to `T2T-CHM13v2.0` (with Y-PAR masked); summary bin tables and regional Y atlases included in [`analysis/02_genome_engineering/2.3_lpwgs_characterization/`](analysis/02_genome_engineering/2.3_lpwgs_characterization/).
- **In-house Genome-Wide Brunello CRISPR–Cas9 Screens (HT1080 & A549):** Guide and gene-level summary tables (`drugZ`, `MAGeCK`, `IsogenicZ`) included in [`analysis/06_crispr_screens/`](analysis/06_crispr_screens/).
- **Public Pharmacogenomic & Dependency Cohorts:** DepMap Public (`24Q4` / `26Q1`), Sanger Project Score, DepMap NextGen 3D organoid/neurosphere screens, PRISM Repurposing (`Primary`, `Secondary`, `REP.1M`), GDSC1/2, CTD2, BeatAML, and TCGA (`GDC`).
- **External Transcriptomic & Isogenic Screen Benchmarks:** 18 external Y-loss/mosaic/45,X transcriptomic contrasts (`GSE126712`, `GSE264742`, `GSE301297`, `GSE229233`, `GSE152001`, `GSE313572`, Tahoe-100M, Page Lab dataset), 16 chromosome/arm-loss contrasts (`GSE130306` RPTEC Yq/3p, `GSE227499` ReDACT, `GSE190497` RPE1 monosomies), and 27 external isogenic CRISPR screens (Desjardins et al., 2026; Feng et al., 2022; BioGRID ORCS `v2.0.18`).

---

## Contributors & Co-Authorship

- **Michael Cory Vermeulen** — Study conception and design, *in silico* T2T-Y guide discovery, CRISPR–Cas9 RNP Y-chromosome engineering, ddPCR and lpWGS characterization, genome-wide Brunello CRISPR–Cas9 screen execution and library preparation, spectral flow cytometry competitive growth assays, all bioinformatic and statistical analyses, figure generation, and manuscript/thesis writing.
- **Peter Truesdell** (Craig Lab) — Lentiviral *UTY* and *KDM5D* rescue transductions and cell line generation, `eGFP`/`mCherry` fluorescent reporter cell line establishment for competitive growth assays, spheroid and colony growth assays, and assistance with cell culture, CRISPR–Cas9 RNP nucleofections, genome-wide CRISPR screen execution, and CellTiter-Glo viability assays.
- **Gwenaëlle Conseil** (Craig Lab) — Assistance with RT-qPCR and multiplex ddPCR assays.
- **Benjamin Snider** (Craig Lab) — Western blot validations of *UTY* and *KDM5D* rescue lines, *DDX3X* inhibitor testing, HT1080 small-molecule screening, and assistance with tissue culture.
- **Stephanie Young** (Craig Lab) — Assistance with *UTY* and *KDM5D* add-back experiments, STR profiling preparation, and laboratory logistics.
- **Dr. Tomas Babak & Dr. Andrew W. Craig** — Co-supervisors; study supervision, funding acquisition, and editorial guidance.

---

## Funding

Supported by a **Canadian Institutes of Health Research (CIHR) Project Grant (178214)**, an **Ontario Graduate Scholarship (OGS)**, and the **E. G. Bauman Fellowship**.

---

## Citation

- **Thesis:** Vermeulen, M. C. *Molecular Consequences and Synthetic Lethal Dependencies of Loss of the Y Chromosome in Cancer*. PhD thesis, Department of Biology, Queen's University, Kingston, Ontario, Canada (2026).
- **Chapter 2 Publication:** Vermeulen, M., Craig, A. W. & Babak, T. Challenges and opportunities for oncology drug repurposing informed by synthetic lethality. *npj Systems Biology and Applications* **11**, 143 (2025). https://doi.org/10.1038/s41540-025-00618-7

---

## Contact

**Michael Cory Vermeulen**
Babak Lab & Craig Lab · Department of Biology & Sinclair Cancer Research Institute, Queen's University
Kingston, Ontario, Canada
