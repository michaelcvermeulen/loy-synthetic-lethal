# `01_depmap_prism_public/` — Public DepMap CRISPR & PRISM Pharmacogenomic Analyses

**Thesis Reference:** Chapter 4, §4.2.1, §4.3.0, and Supplementary Methods S4.1 (`Fig 4.1`, `Supp Fig 4.1–4.4`; builds on the pan-cancer synthetic lethality framework of Chapter 2, *npj Systems Biology and Applications* **11**, 143, 2025).

This module evaluates LOY-associated genetic dependencies and pharmacological sensitivities across public 2D cell-line panels, 3D organoid/neurosphere cohorts, and high-throughput drug screens, and quantifies the statistical power limits of stratified public datasets using empirical benchmark simulations.

---

## Sub-Analyses & Pipeline Architecture

```text
01_depmap_prism_public/
├── README.md
├── EXPORT_README.md                        # HPC export provenance & verification checks
├── EXPORT_CHECKS.txt                       # Automated validation check results
├── scripts/                                # Modular R pipeline (00_smoke_test.R – 36_driver_type_power.R)
│   ├── 01_inputs_check.R – 07_freeze.R     # Y-evidence extraction, consensus Y-status calling & cohort freeze
│   ├── 08_pipeline_controls.R              # Benchmark positive-control SL recovery
│   ├── 10_gene_pan.R – 13_context_drugs.R  # Pan-cancer & context-stratified genetic and drug models
│   ├── 14_expression_de.R – 18_forest_*.R  # Transcriptomic DE, Random Forest & Boruta feature selection
│   ├── 19_hit_robustness.R                 # Alternative Y calls, lineage adjustment & leave-one-lineage-out
│   ├── 24_power_context.R                  # Empirical power simulations across 25 benchmark SLs
│   ├── 25_organoid_prepare.R – 33_combine_three.R # Sanger Project Score, RNAi & DepMap NextGen 3D replication
│   └── 34_driver_type_prepare.R – 36_driver_type_power.R # Cancer-type × driver subset analyses
├── docs/
│   ├── DATA_GUIDE.md                       # Schema & column dictionary for all derived tables
│   └── LOY_REANALYSIS_OUTPUTS.md           # Comprehensive numerical summary of all screen families
├── results/                                # Top-level summary tables across genetic & drug screens
└── tables/
    └── depmap_drug_within_subsets/         # Within-subset drug screen summary tables & SR-1664 records
```

---

## Key Findings

1. **Consensus Y-Status Classification & Core Gametolog Recovery (`Fig 4.1a–e`, `Supp Fig 4.1`):**
   - Integrating chromosome Y copy number, MSY gene expression, and mutation calls across DepMap (`24Q4`/`26Q1`) classifies **568 male CRISPR-screened models** into **404 Y-retaining (`Y+`)** and **164 `LOY`** models.
   - Pan-cancer linear modelling (`limma` and rank-normalised OLS) and Random Forest / Boruta feature selection independently identify four X-linked gametologs as top LOY synthetic lethalities: ***DDX3X***, ***EIF1AX***, ***RPS4X***, and ***ZFX***.
   - Sensitivity analyses confirm that all four dependencies remain robust across alternative Y-status definitions, lineage covariate adjustments, and leave-one-lineage-out jackknifing (`Supp Fig 4.1b`), and correlate inversely with expression of their Y-linked partners (*DDX3Y*, *EIF1AY*, *RPS4Y1*, *ZFY*).
2. **Replication in 3D Malignant Organoid & Neurosphere Models (`Fig 4.1f`):**
   - Extending Y-status classification to Sanger Project Score, DEMETER2 RNAi, and **DepMap NextGen 3D patient-derived organoid and neurosphere CRISPR screens** confirms permutation-significant replication of ***DDX3X***, ***EIF1AX***, and ***RPS4X*** in 3D models, whereas *ZFX* does not replicate in NextGen.
3. **Context-Stratified Genetic Searches & Empirical Power Limits (`Supp Fig 4.2–4.3`):**
   - Testing LOY-versus-`Y+` dependency differences within individual cancer lineages/subtypes, across **86 driver/genomic alteration interaction contexts** (523,654 tests), and within **24 cancer-type $\times$ driver combinations** (414,052 within-subset tests; 146,136 driver-interaction tests) yields no additional non-paralog discoveries at $\mathrm{FDR} < 0.05$ (post hoc nominal signals include weaker dependency on *INSM1* in *MYC*-altered lung models [$\mathrm{FDR} = 0.161$] and *ATRX* in *CDKN2A*-altered melanoma [$\mathrm{FDR} = 0.196$]; `Supp Fig 4.2d,e`).
   - Empirical power simulations transferring **25 established benchmark synthetic lethalities** (e.g., *WRN*–MSI, *PRMT5*–*MTAP*, *ARHGEF2*–CIN, *CIB1*–bladder) onto actual LOY sample splits (`Supp Fig 4.3`) demonstrate that once public panels are subdivided by tumour type or driver mutation, limited sample sizes ($n_{\mathrm{LOY}}$ often $<10$) and inter-line variance prevent detection of even benchmark-strength vulnerabilities—establishing the rationale for matched isogenic LOY screens.
4. **Pan-Cancer & Context-Stratified Drug Screens (`Supp Fig 4.4`):**
   - Across **6,413 pan-cancer drug comparisons** spanning seven screens (PRISM Repurposing Primary & Secondary, PRISM REP.1M, GDSC1, GDSC2, CTD2, plus BeatAML) and **495,079 finite within-subset tests** across 121 binary contexts (13 lineages, 12 primary diseases, 8 subtypes), positive controls (`MDM2`i–*TP53*, `BRAF`i–*BRAF*, `EGFR`i–*EGFR*, `PARP`i–*BRCA1/2*) behave as expected, but no broadly shared LOY drug sensitivity emerges (the sole within-subset $\mathrm{FDR} < 0.05$ association is reduced sensitivity to `SR-1664` in LOY mature B-cell neoplasms in PRISM REP.1M; 8 LOY vs. 16 `Y+`, $\Delta = +0.667$, $\mathrm{FDR} = 0.031$; `Supp Fig 4.4f`).

---

## Thesis Figures Generated from This Module

- **Fig 4.1:** DepMap reanalysis recovers four established LOY-associated dependencies ([`ch03_main_depmap_y_status_dependencies.py`](../../figures/chapter3/code/ch03_main_depmap_y_status_dependencies.py)).
- **Supp Fig 4.1:** Y-status classification, dependency robustness and genetic context analyses ([`ch03_supp_depmap_combined_checks.py`](../../figures/chapter3/code/ch03_supp_depmap_combined_checks.py)).
- **Supp Fig 4.2:** Coverage, detection limits and candidate dependencies in cancer subsets ([`ch03_main_depmap_cancer_type_detection_limits.R`](../../figures/chapter3/code/ch03_main_depmap_cancer_type_detection_limits.R)).
- **Supp Fig 4.3:** Detection of known vulnerabilities moved into one cancer type ([`ch03_supp_depmap_benchmark_detection.R`](../../figures/chapter3/code/ch03_supp_depmap_benchmark_detection.R)).
- **Supp Fig 4.4:** Drug screens recover expected controls but identify no broadly shared LOY sensitivity ([`ch03_main_depmap_drug_response.py`](../../figures/chapter3/code/ch03_main_depmap_drug_response.py)).
