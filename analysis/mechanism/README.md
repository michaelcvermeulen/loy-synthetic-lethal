# mechanism/

Mechanism-of-action writeups. Each file lays out a specific theory about how
LOY drives the observed phenotypes and walks through which results from the
numbered analyses support, refine, or qualify that theory. This is where the
narrative argument of the manuscript gets assembled — analyses on their own
report findings; mechanism notebooks make claims and defend them.

Survey of existing repo content turned up no dedicated mechanism writeups
beyond `depmap_tcga_loy/04_SAGA_Complex_Analysis.Rmd`, which is partially
playing this role for the SAGA story. These extract and unify those arguments.

## Planned mechanism writeups

| File                                       | Theory |
|--------------------------------------------|--------|
| `th01_kdm5d_uty_h3k4me3.Rmd`               | Loss of *KDM5D* (H3K4me2/3 demethylase) and *UTY* (H3K27me3 demethylase) on LOY drives global H3K4me3 accumulation at Polycomb-marked loci. Evidence from RNA-seq + public ChIP-seq + ATAC-seq. TODO. |
| `th02_prc2_redistribution.Rmd`             | H3K4me3 accumulation forces PRC2 to redistribute, depleting H3K27me3 from canonical Polycomb targets. Evidence from ChIP-seq + ATAC-seq + scRNA-seq. TODO. |
| `th03_polycomb_derepression.Rmd`           | PRC2 redistribution derepresses Polycomb-silenced neural and developmental genes — the "wrong-program" signature seen in LOY RNA-seq across cell lines. TODO. |
| `th04_saga_synthetic_lethality.Rmd`        | The derepressed program creates addiction to the SAGA/STAGA coactivator complex, exposing a synthetic-lethal vulnerability. Evidence from DepMap + CRISPR screen + drug screen. TODO. |

Add more mechanism notebooks as the story develops (e.g., escape genes,
cell-context dependencies, prostate-specific axes). Numbering is loose —
renumber as needed.
