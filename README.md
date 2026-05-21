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
├── LICENSE                  # MIT
├── .gitignore
├── renv.lock                # pinned R package versions
├── .Rprofile                # activates renv on session start
├── renv/                    # activate.R tracked; library/ ignored
├── analysis/                # .Rmd source for each numbered analysis
├── R/                       # shared helper functions sourced by .Rmd files
├── data/                    # local-only; see data/README.md for layout
├── figures/                 # publication-ready figure exports
└── docs/                    # knitted HTML reports + GitHub Pages landing page
```

## Analyses

| #  | Analysis                | File                                         | Description |
|----|-------------------------|----------------------------------------------|-------------|
| 01 | CRISPR screen           | `analysis/01_crispr_screen.Rmd`              | TODO        |
| 02 | RNA-seq (HT1080)        | `analysis/02_rnaseq_ht1080.Rmd`              | TODO        |
| 03 | RNA-seq (LNCaP, HepG2)  | `analysis/03_rnaseq_lncap_hepg2.Rmd`         | TODO        |
| 04 | ATAC-seq                | `analysis/04_atacseq.Rmd`                    | TODO        |
| 05 | ChIP-seq (Polycomb)     | `analysis/05_chipseq_polycomb.Rmd`           | TODO        |
| 06 | scRNA-seq (LOY atlas)   | `analysis/06_scrnaseq_loy_atlas.Rmd`         | TODO        |
| 07 | DepMap pan-cancer       | `analysis/07_depmap_pancancer.Rmd`           | TODO        |
| 08 | ENCODE tissue ATAC      | `analysis/08_encode_tissue_atac.Rmd`         | TODO        |

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

## Citation

Thesis: `PLACEHOLDER: Vermeulen M. [Thesis title]. PhD thesis, Queen's
University, 2026.`

Publication: `PLACEHOLDER: full citation once accepted.`

## Contact

Mike Vermeulen
Babak Lab + Craig Lab
Department of Biology, Queen's University
Email: `PLACEHOLDER: email`
