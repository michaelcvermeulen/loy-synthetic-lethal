# analysis/

R Markdown source for the manuscript and thesis chapter. Organized by topic
section; sub-analyses live in numbered subdirectories. Each `.Rmd` knits to a
self-contained HTML report in [`../docs/`](../docs/).

```
analysis/
├── 01_genome_engineering/        # 1.1 sgRNA design · 1.2 ddPCR · 1.3 lpWGS
├── 02_rnaseq/                    # 2.1 A549 · 2.2 HT1080 · 2.3 LNCaP · 2.4 combined
├── 03_chipseq_public/            # public ChIP-seq (Postovit ChIP not yet generated)
├── 04_atacseq_multiomic/         # ATAC-seq + ATAC/RNA-seq integration
├── 05_depmap_pancancer/          # DepMap dependency analysis
├── 06_tcga/                      # TCGA pan-cancer LOY analysis
├── 07_scrnaseq_loy_atlas/        # single-cell LOY atlas (prostate + others)
├── 08_functional_screens/        # CRISPR synthetic-lethal screen + drug screens
├── methods/                      # per-assay / per-technology methods writeups
└── mechanism/                    # mechanism-of-action theories and supporting evidence
```

Sequence: **methods** documents protocols → **analysis** sections execute and
interpret per-assay results → **mechanism** integrates findings across sections
to argue specific theories (KDM5D/UTY loss, PRC2 redistribution,
SAGA/STAGA addiction).

Helpers shared across analyses go in [`../R/`](../R/); inputs follow
[`../data/README.md`](../data/README.md).
