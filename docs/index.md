# LOY Synthetic Lethality — Thesis & Analysis Reports

Landing page for the computational pipelines, methods, and narratives supporting the PhD thesis **_Molecular Consequences and Synthetic Lethal Dependencies of Loss of the Y Chromosome in Cancer_** (Michael Cory Vermeulen, Queen's University, December 2026).

---

## Read First

- **[Background](background.html)** — Somatic LOY in male cancers, `T2T-CHM13v2.0` Y-chromosome architecture, the X–Y gametolog buffering problem, and the central hypotheses tested across Chapters 2–4.
- **[Condensed Summary](summary.html)** — Quantitative headline findings across all seven analysis modules (`01`–`07`), mapped to every Chapter 3 (`Fig 3.1–3.15`, `Supp Fig 3.1–3.44`) and Chapter 4 (`Fig 4.1–4.8`, `Supp Fig 4.1–4.16`) figure.
- **[Project & Thesis Overview](project_plan.html)** — Architecture of Chapters 1–5, detailed section-by-section crosswalk to repository directories, and reproducibility standards.

---

## Analysis Modules (`analysis/`)

1. **[`01_depmap_prism_public/`](../analysis/01_depmap_prism_public/)** — Public DepMap CRISPR (568 male models), 3D NextGen organoid/neurosphere replication, PRISM/GDSC/CTD2 drug screens, and 25-benchmark empirical power simulations (Chapter 4, §4.3.0; links Chapter 2).
2. **[`02_genome_engineering/`](../analysis/02_genome_engineering/)** — `T2T-CHM13v2.0` repeat guide discovery (`CP`, `CN`, `Q53`, `P47`, `P51`), pooled Cas9 RNP Y depletion, multiplex ddPCR, panobinostat/KU-55933 modulation, and Y-PAR-masked lpWGS regional Y atlases across HT1080, A549, and DU145 (Chapter 3, §3.3.1–3.3.2).
3. **[`03_rnaseq_transcriptomics/`](../analysis/03_rnaseq_transcriptomics/)** — Bulk RNA-seq (`STAR v2` Y-PAR-masked / `GENCODE v49`) in HT1080, A549, and DU145, cross-line synthesis across 18 external Y-loss/mosaic datasets, and comparison with 16 chromosome/arm-loss models (Chapter 3, §3.3.3).
4. **[`04_atacseq_chromatin/`](../analysis/04_atacseq_chromatin/)** — HT1080 ATAC-seq differential accessibility (663 closed vs. 95 opened peaks), `chromVAR`/HOMER/TOBIAS motif footprinting, and multi-omic RNA/ATAC integration at PRC2/H3K27me3-enriched promoters (Chapter 3, §3.3.3.1).
5. **[`05_uty_kdm5d_addback/`](../analysis/05_uty_kdm5d_addback/)** — Lentiviral *UTY* and *KDM5D* rescue RNA-seq in Y-depleted HT1080 cells, demonstrating 51% (and 35%) reversal of the 150-gene LOY signature by *UTY* and preferential re-induction of closed-promoter genes (Chapter 3, §3.3.4).
6. **[`06_crispr_screens/`](../analysis/06_crispr_screens/)** — Genome-wide Brunello CRISPR–Cas9 knockout screens in matched Y-depleted and control HT1080 and A549 populations (*DDX3X*, *EIF1AX*, SAGA complex, *MED12*, *SREBF1/2*, *TIMM8A*, *L3MBTL2*, *UBR5*), multi-omic integration, and benchmarking across 27 external isogenic CRISPR screens (Chapter 4, §4.3.1, §4.3.5, §4.3.6).
7. **[`07_competitive_growth_assays/`](../analysis/07_competitive_growth_assays/)** — Reciprocal two-colour `eGFP`/`mCherry` spectral flow-cytometry competitive growth assays (`BD FACSDiscover A8`) validating *DDX3X* and *TADA2B* knockout fitness penalties (Chapter 4, §4.3.2–4.3.4).

---

## Methods, Mechanism & Figure Registers

- **[`analysis/methods/`](../analysis/methods/)** — Complete wet-lab and bioinformatic methods writeups aligned with Chapter 3 (§3.2, S3.1–S3.4) and Chapter 4 (§4.2, S4.1–S4.5).
- **[`analysis/mechanism/`](../analysis/mechanism/)** — Mechanistic synthesis documents covering X–Y gametolog buffering rules, lineage specificity vs. conserved interferon suppression, *UTY*-reversible chromatin compaction, and SAGA/metabolic synthetic lethalities.
- **[`figures/chapter3/chapter3_figure_numbering.md`](../figures/chapter3/chapter3_figure_numbering.md)** & **[`figures/chapter4/chapter4_figure_numbering.md`](../figures/chapter4/chapter4_figure_numbering.md)** — Canonical figure registers and captions.
