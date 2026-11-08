# DU145 RNA-seq thesis reanalysis, 23 September 2026

Open [the analysis report](outputs/report/analysis.html). Complete numerical results are in `outputs/tables`, models in `outputs/models`, and working figure exports in `outputs/figures` (PDF, SVG and 600-dpi PNG) with captions in `outputs/captions`. Thesis-styled figures are built separately under `SEPT_11_2026_THESIS/figures/chapter3/`.

This is the numerical source for Chapter 3 Results 3.4.7.3 (author decision, 23 September 2026). It is a new analysis of the existing featureCounts tables. It follows the HT1080 `SEPT_11_2026_REANALYZE` and A549 `SEPT_18_2026_REANALYZE` core: GENCODE v49, filterByExpr (min.count 10, min.total.count 15), TMM on expressed autosomes and non-PAR X, voom-limma with robust empirical Bayes, and frozen MSigDB 2026.1.Hs, CollecTRI and PROGENy resources copied with SHA-256 manifests. It does not repeat alignment or counting. Original count files and earlier branches in `analysis/LNCAP/` are unchanged; they are provenance and sensitivity sources only.

## Identity

The directory, sample IDs and WGS files say LNCaP. The RNA identifies the line as DU145: AR, KLK3, KLK2, FOLH1 and HOXB13 have zero reads in every library, and DU145 is the best DepMap 26Q1 match. Hash-verified release files are listed in `outputs/tables/identity_depmap_input_hashes.csv`. Source labels are preserved. STR profiling was not performed here.

## Design decisions (author, 23 September 2026)

- Six libraries from one sampling window (May 2026): three sgCTRL and three LOY at matched passages P12, P13 and P14. There are no parental libraries.
- The LOY population is the CP + Q53 + P47 2x population, so lpWGS `LNCaP-LOY-2X` is its genomic context. RNA and WGS specimens were not matched. The RNA shares the LOY-2X chrX and 15q gains that are absent from `LNCaP-DMSO-LOY-1X`.
- The primary model is paired, `~ passage + condition` (residual df 2). Passages are repeated sampling of one population per arm, not independent editing events. The unpaired model is a sensitivity check.
- The thesis DE rule is FDR < 0.05, |log2FC| >= 1 and the same direction in all three passage pairs. This is the DU145 analogue of the HT1080 rule requiring agreement across the August and May windows.
- Chromosome Y and PAR genes are excluded from the primary test family. They are analysed in a separate inclusive paired fit that uses the same TMM exposures (`Y_PAR_gene_table`, `DE_inclusive_nuclear_paired_voom`) and are positioned on CHM13v2.0 chrY against saved lpWGS segments.

## Copy number

GENCODE v49 genes are mapped to CHM13v2.0 by stable Ensembl ID from the CAT/Liftoff annotation (`WGS/T2T_genes.gz`), with a unique same-chromosome symbol fallback. Gene-level dosage is the median saved QDNAseq `log2_seg` over usable 100-kb bins overlapping the gene; no resegmentation is done. A gene counts as copy-number-changed when |LOY-2X minus sgCTRL| >= 0.2. Pathway results carry leading-edge arm composition and a rerun with copy-number-changed genes removed (`survives_CN_screen`).

## Reproduce

From this directory in PowerShell:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' --vanilla run.R
```

`--core-only` skips enrichment, `--enrichment-only` resumes from verified core models, and `--report-only` re-renders. The DepMap 26Q1 files are read from `H:/PRISM/data/CCLE/26Q1` and verified against that release's manifest. `provenance/` holds validation, code hashes and session information.
