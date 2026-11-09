# `04_atacseq_chromatin/`

HT1080 bulk ATAC-seq processing, differential accessibility, motif/footprinting (`chromVAR`, `HOMER`, `TOBIAS`), and RNA-seq/ENCODE chromatin integration (`Chapter 3, Section 3.3.3.1`; `Fig 3.9–3.10`, `Supp Fig 3.3`, `Supp Fig 3.21–3.28`).

- [`code/`](code/) — R and Python analysis scripts (`00_environment.R` – `89_transcript_TSS_promoter_sensitivity.R`).
- [`config/`](config/) — Sample manifests, design matrices, and PAR coordinate definitions.
- [`qc/`](qc/) — Alignment, fragment size, TSS enrichment, and FRiP QC tables.
- [`motifs/`](motifs/) — `chromVAR`, `HOMER`, and `TOBIAS` motif/footprinting summary tables.
- [`outputs/`](outputs/) — Differential accessibility, promoter-RNA integration, and PRC2/H3K27me3 overlap tables.
