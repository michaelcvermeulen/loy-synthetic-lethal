# `analysis/`

Computational pipelines, statistical scripts, configuration manifests, and summary output tables supporting Chapters 3 and 4 (and linking Chapter 2) of the PhD thesis **_Molecular Consequences and Synthetic Lethal Dependencies of Loss of the Y Chromosome in Cancer_**.

---

## Numbered Analysis Modules

```text
analysis/
├── 01_depmap_prism_public/           # Ch 4 Section 4.3.0: DepMap CRISPR, 3D NextGen, PRISM/GDSC/CTD2, power sims
├── 02_genome_engineering/            # Ch 3 Section 3.3.1–3.3.2: T2T-Y guide design, multiplex ddPCR, Y-PAR-masked lpWGS
├── 03_rnaseq_transcriptomics/        # Ch 3 Section 3.3.3: HT1080, A549, DU145 RNA-seq, 18 external LOY & 16 arm-loss sets
├── 04_atacseq_chromatin/             # Ch 3 Section 3.3.3.1: HT1080 ATAC-seq, chromVAR/HOMER/TOBIAS, RNA/ATAC integration
├── 05_uty_kdm5d_addback/             # Ch 3 Section 3.3.4: HT1080 lentiviral UTY & KDM5D rescue RNA-seq & promoter recovery
├── 06_crispr_screens/                # Ch 4 Section 4.3.1, 4.3.5, 4.3.6: HT1080 & A549 Brunello screens, 27-screen benchmark
└── 07_competitive_growth_assays/     # Ch 4 Section 4.3.2–4.3.4: BD FACSDiscover A8 spectral flow competition assays
```

### Workflow Sequence

1. **Public pharmacogenomic & CRISPR dependency screens ([`01_depmap_prism_public/`](01_depmap_prism_public/))** establish pan-cancer X–Y gametolog dependencies (*DDX3X*, *EIF1AX*, *RPS4X*, *ZFX*), replicate them in 3D organoid/neurosphere cohorts, and prove via empirical power simulations that lineage- and driver-stratified public subsets are underpowered to detect context-specific vulnerabilities—motivating matched isogenic LOY models.
2. **Isogenic model engineering & structural characterization ([`02_genome_engineering/`](02_genome_engineering/))** designs repeat-targeting Cas9 RNP guides against `T2T-CHM13v2.0` chrY (`CP`, `CN`, `Q53`, `P47`, `P51`), generates non-clonal Y-depleted populations across HT1080, A549, DU145, and HepG2, and combines multiplex ddPCR with Y-PAR-masked T2T-Y lpWGS to resolve whole-Y loss (A549) versus complete Yq loss with partial distal Yp (*RPS4Y1*/*ZFY*) retention (HT1080, DU145).
3. **Comparative transcriptomics & chromatin profiling ([`03_rnaseq_transcriptomics/`](03_rnaseq_transcriptomics/), [`04_atacseq_chromatin/`](04_atacseq_chromatin/), [`05_uty_kdm5d_addback/`](05_uty_kdm5d_addback/))** demonstrate universal absence of X-paralog compensation, strong lineage specificity of individual autosomal DEGs alongside conserved suppression of interferon/inflammatory signalling across 18 external Y-loss/mosaic datasets, widespread promoter chromatin closing (~7:1 skew) enriched at PRC2/H3K27me3 targets in HT1080, and partial transcriptional reversal by lentiviral *UTY* (but not *KDM5D*) rescue.
4. **Isogenic genome-wide CRISPR screens & competitive growth validation ([`06_crispr_screens/`](06_crispr_screens/), [`07_competitive_growth_assays/`](07_competitive_growth_assays/))** uncover *DDX3X* (`normZ = -11.45` in HT1080), *EIF1AX* (in A549), and candidate vulnerabilities in SAGA (*TAF6L*, *TADA2B*, *SUPT20H*, *ATXN7L3*, *TADA1*, *TAF5L*), Mediator (*MED12*), sterol/selenoprotein metabolism (*SREBF2*, *SREBF1*, *SCLY*, *ACSL4*), *TIMM8A*, *L3MBTL2*, and *UBR5*, benchmark their specificity across 27 external isogenic screens, and validate *DDX3X* and *TADA2B* via reciprocal two-colour spectral flow cytometry competition assays.

Shared R helpers live in [`../R/`](../R/), thesis figure renderers and numbering registers live in [`../figures/`](../figures/), and local input paths follow [`../data/README.md`](../data/README.md).
