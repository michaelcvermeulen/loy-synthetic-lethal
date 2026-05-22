# LOY synthetic lethality

This repository contains the analyses supporting **[MANUSCRIPT TITLE PLACEHOLDER]**
and the corresponding chapter of my PhD thesis. Together they characterize loss
of Y chromosome (LOY) as a context-specific vulnerability in cancer. We propose
a mechanistic model in which LOY-driven loss of *KDM5D* and *UTY* leads to
H3K4me3 accumulation and PRC2 redistribution, derepressing Polycomb-silenced
neural and developmental genes, and creating addiction to the SAGA/STAGA
coactivator complex as a synthetic lethal dependency.

Source `.Rmd` files live in [`analysis/`](analysis/); knitted HTML reports are
served from [`docs/`](docs/) via GitHub Pages once the repo is made public at
manuscript submission.

---

## Repository structure

```
loy-synthetic-lethality/
├── README.md
├── LICENSE                          # MIT
├── .gitignore
├── renv.lock                        # pinned R package versions
├── .Rprofile                        # activates renv on session start
├── renv/                            # activate.R tracked; library/ ignored
├── analysis/
│   ├── 01_genome_engineering/       # 1.1 sgRNA · 1.2 ddPCR · 1.3 lpWGS
│   ├── 02_rnaseq/                   # 2.1 A549 · 2.2 HT1080 · 2.3 LNCaP · 2.4 combined
│   ├── 03_chipseq_public/           # public ChIP-seq (Postovit ChIP not yet generated)
│   ├── 04_atacseq_multiomic/        # ATAC-seq + ATAC/RNA-seq integration
│   ├── 05_depmap_pancancer/         # DepMap dependency analysis
│   ├── 06_tcga/                     # TCGA pan-cancer LOY analysis
│   ├── 07_scrnaseq_loy_atlas/       # single-cell LOY atlas (prostate + others)
│   ├── 08_functional_screens/       # CRISPR A549/HT1080 + drug A549/HT1080 + integration
│   ├── methods/                     # per-assay methods writeups
│   └── mechanism/                   # mechanism-of-action theories + supporting evidence
├── R/                               # shared helper functions sourced by .Rmd files
├── data/                            # local-only; see data/README.md for layout
├── figures/                         # publication-ready figure exports
└── docs/                            # knitted HTML reports + GitHub Pages landing page
```

## Analyses

Each numbered section lives in its own subdirectory under
[`analysis/`](analysis/); see the per-section `README.md` for sub-analyses and
links to existing upstream work.

| #   | Section                                                       | Sub-analyses                                                           |
|-----|---------------------------------------------------------------|------------------------------------------------------------------------|
| 1   | [Genome engineering](analysis/01_genome_engineering/)         | 1.1 sgRNA guide creation · 1.2 ddPCR · 1.3 lpWGS                       |
| 2   | [Bulk RNA-seq](analysis/02_rnaseq/)                           | 2.1 A549 · 2.2 HT1080 · 2.3 LNCaP · 2.4 combined cross-cell-line       |
| 3   | [ChIP-seq (public)](analysis/03_chipseq_public/)              | Postovit in-house ChIP pending; public sources only for now            |
| 4   | [ATAC-seq + multiomic](analysis/04_atacseq_multiomic/)        | Differential accessibility + ATAC/RNA integration                      |
| 5   | [DepMap pan-cancer](analysis/05_depmap_pancancer/)            | Dependency profiling, SAGA-focused                                     |
| 6   | [TCGA](analysis/06_tcga/)                                     | LOY calling, DE, pathway + mutation enrichment                         |
| 7   | [scRNA-seq LOY atlas](analysis/07_scrnaseq_loy_atlas/)        | Prostate first; additional cancer types to come                        |
| 8   | [Functional screens](analysis/08_functional_screens/)         | 8.1 CRISPR A549 · 8.2 CRISPR HT1080 · 8.3 Drug A549 · 8.4 Drug HT1080 · 8.5 cross-line integration |

### Methods and mechanism

- [`analysis/methods/`](analysis/methods/) — per-assay / per-technology
  protocol writeups (genome engineering, sgRNA design, ddPCR, lpWGS, bulk
  RNA-seq, ChIP-seq, ATAC-seq, scRNA-seq, screens, DepMap/TCGA processing).
- [`analysis/mechanism/`](analysis/mechanism/) — mechanism-of-action theories
  and how each numbered analysis supports them (KDM5D/UTY → H3K4me3
  accumulation, PRC2 redistribution, Polycomb derepression, SAGA/STAGA
  synthetic lethality).

## Figure-to-file map

| Manuscript Figure | Repo file | Section anchor |
|-------------------|-----------|----------------|
| TODO              | TODO      | TODO           |

## Software environment

Analyses are pinned with [renv](https://rstudio.github.io/renv/). After cloning:

```r
# from the repo root, in R
renv::restore()
```

- R version: R 4.5.2 (pinned in `renv.lock`)
- Full `sessionInfo()` output: `PLACEHOLDER: to be appended after the analyses
  are finalized.`

## Data availability

Raw and processed data are **not** tracked in this repository. The expected
local layout is documented in [`data/README.md`](data/README.md). At manuscript
submission, datasets will be deposited as follows (all marked TBD until then):

- **RNA-seq**: GEO, `TODO: accession`
- **ATAC-seq**: GEO, `TODO: accession`
- **ChIP-seq (in-house)**: GEO, `TODO: accession`
- **CRISPR screen counts**: `TODO: deposition target + accession`
- **Public ChIP-seq**: GSE261834, GSE135580
- **Public scRNA-seq**: `TODO: list GEO accessions`
- **DepMap**: `TODO: pin release (e.g., 24Q4)`
- **ENCODE**: `TODO: list ENCFF/ENCSR portal IDs`

## Contributors

- **Mike Vermeulen** — repository author; primary analyses, integration, and manuscript work.
- **Ben Snider** (Craig Lab) — small-molecule drug screens in LOY vs WT isogenic lines ([§8.2](analysis/08_functional_screens/)); wet-lab execution and primary analysis. Re-analysis for this manuscript by M.V.
- `TODO: list additional collaborators (Postovit Lab ChIP-seq contributors, etc.) as in-house data sets are folded in.`

## Citation

Thesis: `PLACEHOLDER: Vermeulen M. [Thesis title]. PhD thesis, Queen's
University, 2026.`

Publication: `PLACEHOLDER: full citation once accepted.`

## Contact

Mike Vermeulen
Babak Lab + Craig Lab
Department of Biology, Queen's University
Email: `PLACEHOLDER: email`
