# Module 01: Public DepMap CRISPR and PRISM Drug-Sensitivity Analyses (`loy_reanalysis_v1`)

**Thesis Reference:** Chapter 4, Section 4.3.0 (*Public DepMap and PRISM analyses of LOY-associated vulnerabilities*; Methods 4.2.1; formerly drafted in Chapter 3 Section 3.3.1/3.4.1 and moved to Chapter 4 on 26 October 2026)  
**Thesis Figures:** `Fig 4.1` and `Supp Fig 4.1–4.4`

---

## Overview & Key Thesis Results

This module contains the `loy_reanalysis_v1` pipeline that reclassifies functional Y-chromosome status across the DepMap 26Q1 release and evaluates pan-cancer, cancer-type-specific, and driver-stratified genetic and pharmacological vulnerabilities:

1. **Y-Status Reclassification (`02_extract_y_evidence.R`–`05_build_cohort.R`; `Fig 4.1a–c`, `Supp Fig 4.1a`):**
   - Reassigns functional Y status across 1,143 annotated male models using $\log_2(\mathrm{TPM} + 1)$ expression and DNA copy number of the 7 ubiquitously expressed X-degenerate Y genes (*DDX3Y*, *EIF1AY*, *KDM5D*, *RPS4Y1*, *USP9Y*, *UTY*, *ZFY*) against 711 female reference cell lines: **383 Y-negative (LOY; 39.4%)**, **588 Y-positive (60.6%)**, and 172 unresolved (99.5% agreement across 756 of 760 models with definite calls from both RNA and DNA).
2. **Pan-Cancer Genetic Dependencies & Organoid Replication (`10_gene_pan.R`, `16_forest_y_status.R`, `18_forest_confounders.R`, `19_hit_robustness.R`, `25`–`33`; `Fig 4.1d–g`, `Supp Fig 4.1b–d`):**
   - Across **568 eligible CRISPR-screened male models** (230 Y-negative, 338 Y-positive; 17,205 genes), both lineage/library-adjusted `limma` regression and Boruta random-forest selection (over 6,037 selective dependencies) identify the same four X-linked paralogs as the only genome-wide hits ($\mathrm{FDR} < 0.05$):
     - *EIF1AX* ($\Delta\text{Chronos} = -1.23$)
     - *DDX3X* ($\Delta\text{Chronos} = -0.75$)
     - *RPS4X* ($\Delta\text{Chronos} = -0.60$)
     - *ZFX* ($\Delta\text{Chronos} = -0.27$)
   - Replicates *DDX3X* and *RPS4X* in **95 male Sanger malignant organoids** (`25_organoid_prepare.R`, `26_organoid_scan.R`) and *DDX3X*, *EIF1AX*, and *RPS4X* in **65 male DepMap NextGen 3D organoid/neurosphere models** (`31_nextgen_prepare.R`, `32_nextgen_scan.R`, `33_combine_three.R`).
3. **Cancer-Type & Driver-Context Searches and Empirical Power Simulations (`12_context_genes.R`, `24_power_context.R`, `34`–`36`; `Supp Fig 4.1e,f`, `Supp Fig 4.2`, `Supp Fig 4.3`):**
   - Evaluates 30 eligible Oncotree cancer subsets ($\ge 5$ models per Y group; 182,670 tests), 47 pan-cancer driver-mutation contexts, 86 binary/continuous genomic contexts (523,654 interaction tests), and 68 cancer-by-driver strata (414,052 within-altered tests; 146,136 interaction tests).
   - Empirical power simulations (`24_power_context.R`, `36_driver_type_power.R`) transferring 25 benchmark synthetic lethal effects (*WRN*–MSI, *PRMT5*–*MTAP*, etc.) show that 15 of 25 benchmarks fail to reach 80% detection power in any single cancer type (`Supp Fig 4.2c`, `Supp Fig 4.3a`).
4. **Pharmacological Screens (`11_drug_pan.R`, `13_context_drugs.R`; `Supp Fig 4.4`):**
   - Evaluates **6,413 pan-cancer compound–dataset comparisons** across 7 drug-response screens (`CTRPv2`, `GDSC1`, `GDSC2`, and `PRISM` Primary, `REP.1M`, `REP.300`, and Secondary) and **495,079 context-stratified compound tests** across 121 subsets, recovering BRAF and MDM2 positive controls (`Supp Fig 4.4e`) but no broadly shared LOY small-molecule vulnerability ($\mathrm{FDR} < 0.05$).

> **Note:** Non-thesis exploratory scripts from `loy_reanalysis_v1` (`14_expression_de.R`, `15_expression_dependency.R`, `17_forest_expression.R`, and BeatAML `29`–`30`) are stored locally in `_local_non_thesis/analysis/01_depmap_prism_public/` and excluded from GitHub.
