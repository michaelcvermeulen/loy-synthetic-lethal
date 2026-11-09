# Molecular Consequences and Synthetic Lethal Dependencies of Loss of the Y Chromosome in Cancer

**Michael Cory Vermeulen** — PhD Thesis Companion Repository (Department of Biology, Queen's University, December 2026)  
*Supervisors: Dr. Tomas Babak & Dr. Andrew W. Craig*

---

## Thesis Abstract

Loss of tumour suppressor genes through point mutation, deletion, or whole-chromosome loss occurs frequently across cancers but cannot be directly drugged, leaving many patients without targeted therapy options. Synthetic lethality overcomes this barrier by targeting interacting genes whose function becomes essential only when a specific tumour suppressor or chromosome is lost. Despite large-scale screening efforts, prioritizing clinically actionable synthetic lethal interactions remains challenging, and whole-chromosome losses such as somatic loss of the Y chromosome (LOY) in male cancers lack systematic genome-wide molecular and functional profiling. This thesis addresses both gaps by integrating functional and clinical datasets to identify synthetic lethal drug-repurposing opportunities across cancer drivers, and then applying CRISPR-engineered isogenic cell models and genome-wide screens to map the molecular consequences and vulnerabilities of LOY.

First, we established a pan-cancer discovery framework integrating Cancer Dependency Map (DepMap) CRISPR screens, somatic mutual exclusivity in The Cancer Genome Atlas, and PRISM pharmacological screens across cancer driver alterations (**Chapter 2**). While genetic dependencies were enriched for mutual exclusivity in patient tumours and nominated repurposing opportunities such as BRD2 inhibition in *ARID1A*-mutant bladder cancer, small-molecule profiles rarely correlated with CRISPR knockout of their nominal targets. Moreover, applying such pan-cancer analyses to LOY recovered only four broad X-linked paralog dependencies (*DDX3X*, *EIF1AX*, *RPS4X*, and *ZFX*) and lacked power within individual lineages, motivating the generation of matched isogenic models to directly study Y chromosome loss.

Next, we used transient Cas9 ribonucleoprotein delivery targeting Y-chromosome repeats to engineer pooled, non-clonal Y-depleted models in HT1080 fibrosarcoma, A549 lung adenocarcinoma, and DU145 prostate carcinoma cells (**Chapter 3**). Low-pass whole-genome sequencing and transcriptomic profiling revealed that Y depletion did not trigger compensatory X-paralog upregulation and drove strongly lineage-dependent autosomal transcriptional reprogramming. Despite minimal gene-level overlap across lines, interferon and inflammatory signalling was recurrently suppressed. In HT1080 cells, Y depletion caused widespread promoter chromatin closing and Polycomb-target repression, and lentiviral re-expression of *UTY*, but not *KDM5D*, partially reversed the LOY transcriptional signature.

Finally, we performed genome-wide Brunello CRISPR screens in the matched Y-depleted and control HT1080 and A549 populations to identify context-specific genetic dependencies of LOY (**Chapter 4**). These screens recovered *DDX3X* and *EIF1AX*, nominated context-specific vulnerabilities including the SAGA coactivator complex, Mediator, and sterol/selenoprotein metabolism, and targeted competitive growth assays were established to validate candidate dependencies including *DDX3X* and *TADA2B*, informing therapeutic strategies for LOY tumours (**Chapter 5**).

---

## Front-Matter Reports

- **[Scientific Background (`background.html`)](background.html)** — Chapter 1 General Introduction, `T2T-CHM13v2.0` Y-chromosome sequence architecture, `Table 1.1` MSY gene catalog, and Chapter 2 translational synthetic lethality framework.
- **[Quantitative Results Summary (`summary.html`)](summary.html)** — Integrated findings across Chapter 3 (isogenic LOY models, RNA-seq, ATAC-seq, *UTY*/*KDM5D* rescue) and Chapter 4 (DepMap/PRISM power-simulation reanalysis, HT1080 & A549 Brunello CRISPR screens, and BD FACSDiscover A8 competitive growth validation).
- **[Thesis Architecture & Execution Plan (`project_plan.html`)](project_plan.html)** — Study design, sample inventories, and statistical governance across all 7 analysis modules.

---

## Repository Modules (`analysis/`)

1. **[`01_depmap_prism_public`](https://github.com/michaelcvermeulen/loy-synthetic-lethal/tree/main/analysis/01_depmap_prism_public)** — Pan-cancer & subset-level DepMap Chronos CRISPR ($n = 483$ male lines) and PRISM drug sensitivity ($n = 398$ male lines; 6,337 compounds) benchmark and empirical power simulations (`Fig 4.1`, `Supp Fig 4.1–4.4`).
2. **[`02_genome_engineering`](https://github.com/michaelcvermeulen/loy-synthetic-lethal/tree/main/analysis/02_genome_engineering)** — `T2T-CHM13v2.0` repeat-targeting sgRNA design (`Q53`, `CN`, `CP`, `P47`, `P51`), multiplexed `AMELX`/`AMELY` ddPCR screening, and Y-PAR-masked low-pass WGS (`Fig 3.1–3.6`, `Supp Fig 3.4–3.14`).
3. **[`03_rnaseq_transcriptomics`](https://github.com/michaelcvermeulen/loy-synthetic-lethal/tree/main/analysis/03_rnaseq_transcriptomics)** — Isogenic RNA-seq across HT1080, A549, and DU145 (44 in-house libraries), 18 external patient/experimental LOY cohorts, and 16 chromosome/arm-loss contrasts (`Fig 3.7–3.8`, `Fig 3.11–3.13`, `Supp Fig 3.1–3.2`, `3.15–3.20`, `3.29–3.41`).
4. **[`04_atacseq_chromatin`](https://github.com/michaelcvermeulen/loy-synthetic-lethal/tree/main/analysis/04_atacseq_chromatin)** — HT1080 omni-ATAC-seq, DiffBind `DBA_DESeq2`, HOMER/JASPAR motifs, TOBIAS footprinting, ChIP-Atlas PRC2/H3K27me3 enrichment, and promoter-accessibility/RNA integration (`Fig 3.9–3.10`, `Supp Fig 3.3`, `3.21–3.28`).
5. **[`05_uty_kdm5d_addback`](https://github.com/michaelcvermeulen/loy-synthetic-lethal/tree/main/analysis/05_uty_kdm5d_addback)** — HT1080 lentiviral *UTY* and *KDM5D* rescue RNA-seq across 18 libraries (`Fig 3.14–3.15`, `Supp Fig 3.42–3.44`).
6. **[`06_crispr_screens`](https://github.com/michaelcvermeulen/loy-synthetic-lethal/tree/main/analysis/06_crispr_screens)** — Genome-wide Brunello CRISPR knockout screens in HT1080 and A549 and re-analysis of 8 external isogenic CRISPR screen pairs (`Fig 4.2–4.4`, `Fig 4.7–4.8`, `Supp Fig 4.5–4.11`, `4.13–4.16`).
7. **[`07_competitive_growth_assays`](https://github.com/michaelcvermeulen/loy-synthetic-lethal/tree/main/analysis/07_competitive_growth_assays)** — Reciprocal two-colour BD FACSDiscover A8 spectral flow cytometry competitive growth assays (`Fig 4.5–4.6`, `Supp Fig 4.12`).
