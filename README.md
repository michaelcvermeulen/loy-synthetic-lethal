# Molecular Consequences and Synthetic Lethal Dependencies of Loss of the Y Chromosome in Cancer

**Author:** Michael Cory Vermeulen  
**Supervisors:** Dr. Tomas Babak & Dr. Andrew W. Craig  
**Institution:** Department of Biology / Queen's Cancer Research Institute, Queen's University, Kingston, Ontario, Canada (December 2026)  
**Published Chapter 2 Reference:** Vermeulen, M., Craig, A. W. & Babak, T. Challenges and opportunities for oncology drug repurposing informed by synthetic lethality. *npj Systems Biology and Applications* **11**, 143 (2025). [https://doi.org/10.1038/s41540-025-00618-7](https://doi.org/10.1038/s41540-025-00618-7)

---

## Thesis Abstract

Loss of tumour suppressor genes through point mutation, deletion, or whole-chromosome loss occurs frequently across cancers but cannot be directly drugged, leaving many patients without targeted therapy options. Synthetic lethality overcomes this barrier by targeting interacting genes whose function becomes essential only when a specific tumour suppressor or chromosome is lost. Despite large-scale screening efforts, prioritizing clinically actionable synthetic lethal interactions remains challenging, and whole-chromosome losses such as somatic loss of the Y chromosome (LOY) in male cancers lack systematic genome-wide molecular and functional profiling. This thesis addresses both gaps by integrating functional and clinical datasets to identify synthetic lethal drug-repurposing opportunities across cancer drivers, and then applying CRISPR-engineered isogenic cell models and genome-wide screens to map the molecular consequences and vulnerabilities of LOY.

First, we established a pan-cancer discovery framework integrating Cancer Dependency Map (DepMap) CRISPR screens, somatic mutual exclusivity in The Cancer Genome Atlas, and PRISM pharmacological screens across cancer driver alterations (**Chapter 2**). While genetic dependencies were enriched for mutual exclusivity in patient tumours and nominated repurposing opportunities such as BRD2 inhibition in *ARID1A*-mutant bladder cancer, small-molecule profiles rarely correlated with CRISPR knockout of their nominal targets. Moreover, applying such pan-cancer analyses to LOY recovered only four broad X-linked paralog dependencies (*DDX3X*, *EIF1AX*, *RPS4X*, and *ZFX*) and lacked power within individual lineages, motivating the generation of matched isogenic models to directly study Y chromosome loss.

Next, we used transient Cas9 ribonucleoprotein delivery targeting Y-chromosome repeats to engineer pooled, non-clonal Y-depleted models in HT1080 fibrosarcoma, A549 lung adenocarcinoma, and DU145 prostate carcinoma cells (**Chapter 3**). Low-pass whole-genome sequencing and transcriptomic profiling revealed that Y depletion did not trigger compensatory X-paralog upregulation and drove strongly lineage-dependent autosomal transcriptional reprogramming. Despite minimal gene-level overlap across lines, interferon and inflammatory signalling was recurrently suppressed. In HT1080 cells, Y depletion caused widespread promoter chromatin closing and Polycomb-target repression, and lentiviral re-expression of *UTY*, but not *KDM5D*, partially reversed the LOY transcriptional signature.

Finally, we performed genome-wide Brunello CRISPR screens in the matched Y-depleted and control HT1080 and A549 populations to identify context-specific genetic dependencies of LOY (**Chapter 4**). These screens recovered *DDX3X* and *EIF1AX*, nominated context-specific vulnerabilities including the SAGA coactivator complex, Mediator, and sterol/selenoprotein metabolism, and targeted competitive growth assays were established to validate candidate dependencies including *DDX3X* and *TADA2B*, informing therapeutic strategies for LOY tumours (**Chapter 5**).

---

## Lay Abstract

Cancer cells divide uncontrollably because of genetic changes that either switch on growth-promoting genes or destroy protective tumour suppressor genes that normally keep cell growth in check. While overactive cancer proteins can often be blocked with targeted drugs, deleted genes or chromosomes cannot be directly targeted because they are no longer present in the cell, leaving many patients without tailored treatments. One way around this challenge is synthetic lethality, an approach that targets a second gene that cancer cells become dependent on to survive after losing a protective gene or chromosome.

The Y chromosome is best known for determining male sex and carries only a small number of genes compared with other chromosomes. In male cancer patients, however, it is one of the most frequently lost chromosomes in tumour cells. Despite how common loss of the Y chromosome (LOY) is across male cancers, how losing it changes cancer cell behaviour and which remaining genes LOY tumours need to survive remain poorly understood.

In this thesis, we first combined public genetic tests in lab-grown cancer cells with tumour records from thousands of patients to find synthetic lethal gene pairs that matter in human cancers, and tested whether existing drugs could be reused to target them. We found that many drugs did not match the effects of removing the genes they were designed to block, showing the need for direct genetic testing. Across these public datasets, we also identified four genes on the X chromosome that normally share essential jobs with matching genes on the Y chromosome; when the Y chromosome is lost, cancer cells must rely on the remaining X-chromosome copy to survive.

However, because public databases mix together hundreds of very different cancers, we wanted to isolate the effect of losing the Y chromosome on its own to see if additional treatable weaknesses emerge. Using CRISPR gene editing, we removed the Y chromosome from human cancer cells in the lab so we could directly compare identical cells with and without a Y chromosome. Studying these cells showed that losing the Y chromosome changes which genes are turned on or off depending on the cancer type, dampens signals that help the immune system recognize tumours, and more tightly packs DNA so certain genes become harder to activate. Finally, by switching off genes one by one across the entire genome in these matched cells, we uncovered new weaknesses in cancer cells that have lost the Y chromosome, pointing toward future strategies to selectively treat LOY tumours.

---

## Overarching Hypothesis & Thesis Organization

> **Overarching Hypothesis (Section 1.5):** We hypothesize that loss of the Y chromosome alters transcription and chromatin accessibility in a strongly lineage-dependent manner, in part through loss of *UTY* and *KDM5D*, and creates both conserved X-paralog dependencies and additional, context-specific genetic vulnerabilities that can be uncovered by combining pan-cancer pharmacogenomic analyses with matched, non-clonal isogenic LOY models.

- **Chapter 1:** General Introduction covering synthetic lethality in precision oncology, the complete 62.5 Mb `T2T-CHM13v2.0` Y-chromosome sequence architecture (`PAR1`, `PAR2`, `XTR`, `XDR`, ampliconic palindromes, `DYZ3` centromere, `DYZ1`/`DYZ2` Yq12 heterochromatin), the 7 broadly expressed ancestral X–Y gametolog pairs (*RPS4Y1*/*RPS4X*, *ZFY*/*ZFX*, *USP9Y*/*USP9X*, *DDX3Y*/*DDX3X*, *UTY*/*KDM6A*, *KDM5D*/*KDM5C*, *EIF1AY*/*EIF1AX*; **Table 1.1**), mosaic LOY in ageing/immunity, and somatic LOY across 5,014 male TCGA tumours (`Fig 1.1–1.4`).
- **Chapter 2 ([Vermeulen et al., 2025](https://doi.org/10.1038/s41540-025-00618-7)):** Pan-cancer synthetic lethal discovery pipeline integrating DepMap CRISPR screens, TCGA somatic mutual exclusivity, and PRISM drug-sensitivity profiles.
- **Chapter 3 (`analysis/02_`–`05_`, `figures/chapter3/code/`):** *Generation and Molecular Characterization of Isogenic Loss of the Y Chromosome Cancer Models* (`Fig 3.1–3.15`, `Supp Fig 3.1–3.44`).
  - Design of 5 repeat-targeting sgRNAs on `T2T-CHM13v2.0` (`Q53` targeting `DYZ1`/`HSat3 A6` in Yq12 [532 matches]; `CN`/`CP` targeting centromeric `DYZ3` $\alpha$-satellite [53/54 matches]; `P47`/`P51` targeting `TSPY`-associated Yp repeats [47/51 matches]) and transient Cas9 RNP delivery to generate pooled, non-clonal Y-depleted models in **HT1080** fibrosarcoma, **A549** lung adenocarcinoma, and **DU145** prostate carcinoma cells.
  - Multiplexed `AMELX`/`AMELY` + `RPP30` ddPCR, Y-PAR-masked low-pass WGS (`~1–2x`), 44-library RNA-seq, HT1080 omni-ATAC-seq, cross-line meta-analysis across 18 external patient/experimental cohorts and 16 chromosome/arm-loss contrasts, and lentiviral cDNA rescue of *UTY* and *KDM5D* in HT1080.
- **Chapter 4 (`analysis/01_`, `06_`–`07_`, `figures/chapter4/code/`):** *Functional Genomics and Synthetic Lethal Vulnerabilities in Loss of the Y Chromosome Cancer Models* (`Fig 4.1–4.8`, `Supp Fig 4.1–4.16`).
  - Re-evaluation of DepMap Chronos CRISPR ($n = 483$ male lines) and PRISM drug screens ($n = 398$ male lines; 6,337 compounds) paired with empirical statistical power simulations across cancer-type and driver-defined subsets (`Fig 4.1`, `Supp Fig 4.1–4.4`).
  - Genome-wide Brunello CRISPR–Cas9 knockout screens (76,441 sgRNAs; 19,114 genes) across Day 0, Day 12, and Day 21 in matched Y-depleted and control **HT1080** (`Fig 4.2–4.4`) and **A549** (`Fig 4.7–4.8`) populations, benchmarked against 8 external isogenic CRISPR screen pairs (`Supp Fig 4.15–4.16`).
  - Reciprocal two-colour fluorescent competitive growth assays quantified on the BD FACSDiscover A8 spectral flow cytometer (`Fig 4.5–4.6`, `Supp Fig 4.12`), confirming early LOY-selective fitness penalties after disruption of *DDX3X*, *SREBF1*, and *TADA2B*.
- **Chapter 5:** General Discussion synthesizing a two-tier model of LOY synthetic lethality (Tier 1 universal X–Y dosage-buffered paralog dependencies vs. Tier 2 lineage- and state-conditioned vulnerabilities in SAGA/Mediator coactivation and sterol/selenoprotein metabolism) and evaluating therapeutic translation barriers (DDX3X/DDX3Y catalytic-pocket conservation, intratumoural karyotypic heterogeneity, and age-related hematopoietic mLOY).

---

## Repository Architecture

```text
loy-synthetic-lethality/
├── README.md                                  # Thesis overview, abstracts, and repository map
├── loy-synthetic-lethality.Rproj              # RStudio project root
├── renv.lock / renv/                          # Reproducible R package lockfile
├── R/                                         # Shared project utilities
├── data/                                      # Local raw/processed data root (gitignored; see data/README.md)
├── figures/                                   # Thesis figure-generation code
│   ├── chapter3/code/                         # R and Python scripts for Fig 3.1–3.15 & Supp Fig 3.1–3.44
│   └── chapter4/code/                         # R and Python scripts for Fig 4.1–4.8 & Supp Fig 4.1–4.16
└── analysis/                                  # Modular analysis pipelines and summary tables aligned with Chapters 3 & 4
    ├── 01_depmap_prism_public/                # DepMap CRISPR, 3D NextGen, PRISM/GDSC/CTD2 drug reanalysis & power simulations
    ├── 02_genome_engineering/                 # T2T-Y sgRNA design, ddPCR screening, and Y-PAR-masked lpWGS
    │   ├── 2.1_guide_design/
    │   ├── 2.2_ddpcr_validation/
    │   └── 2.3_lpwgs_characterization/
    ├── 03_rnaseq_transcriptomics/             # Isogenic RNA-seq (44 libraries), 18 external cohorts, 16 arm-loss contrasts
    │   ├── 3.1_HT1080_baseline/
    │   ├── 3.2_A549_baseline/
    │   ├── 3.3_DU145_baseline/
    │   ├── 3.4_cross_line_synthesis/
    │   └── 3.5_metascape_enrichment/
    ├── 04_atacseq_chromatin/                  # HT1080 omni-ATAC-seq, DiffBind, HOMER/JASPAR motifs, TOBIAS, PRC2 overlap
    ├── 05_uty_kdm5d_addback/                  # HT1080 lentiviral UTY & KDM5D rescue RNA-seq (18 libraries)
    ├── 06_crispr_screens/                     # Genome-wide Brunello CRISPR screens (HT1080, A549, external benchmarks)
    │   ├── 6.1_HT1080_brunello/
    │   ├── 6.2_A549_brunello/
    │   └── 6.3_external_isogenic_benchmarks/
    └── 07_competitive_growth_assays/          # Reciprocal two-colour BD FACSDiscover A8 competitive growth assays
```

> **Note on Cell Line Identity (`DU145` vs. legacy `LNCaP` labels):** During initial culture, one prostate cancer model was provisionally labelled `LNCaP` in raw sequencing filenames and early script headers. Short tandem repeat (STR) profiling and transcriptome-wide correlation against CCLE/DepMap reference panels (`analysis/03_rnaseq_transcriptomics/3.3_DU145_baseline/R/05_identity.R`) definitively authenticated this line as **DU145**. All thesis text, figures, and documentation report its verified identity (**DU145**), while raw filename provenance is preserved for reproducibility.

---

## Co-Authorship & Funding

- **Chapter 2:** Co-authored by Dr. Andrew W. Craig and Dr. Tomas Babak ([Vermeulen et al., *npj Syst. Biol. Appl.* **11**, 143 (2025)](https://doi.org/10.1038/s41540-025-00618-7)).
- **Chapters 3 & 4:** Co-authored by Peter Truesdell and Gwenaëlle Conseil, with experimental contributions from Benjamin Snider and Stephanie Young, under the supervision of Dr. Tomas Babak and Dr. Andrew W. Craig.
- **Funding:** Supported by a Canadian Institutes of Health Research (CIHR) Project Grant (`178214`), an Ontario Graduate Scholarship (OGS), and the E. G. Bauman Fellowship.
