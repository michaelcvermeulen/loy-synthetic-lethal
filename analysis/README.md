# analysis/

R Markdown source for each analysis in the manuscript and thesis chapter.
Files are numbered to match the [analyses table in the repo README](../README.md#analyses):

```
01_crispr_screen.Rmd
02_rnaseq_ht1080.Rmd
03_rnaseq_lncap_hepg2.Rmd
04_atacseq.Rmd
05_chipseq_polycomb.Rmd
06_scrnaseq_loy_atlas.Rmd
07_depmap_pancancer.Rmd
08_encode_tissue_atac.Rmd
```

Each `.Rmd` is knitted to a self-contained HTML report written to
[`../docs/`](../docs/), which is served as GitHub Pages once the repo is public.
Shared helper functions are sourced from [`../R/`](../R/); local data inputs
follow the layout in [`../data/README.md`](../data/README.md).
