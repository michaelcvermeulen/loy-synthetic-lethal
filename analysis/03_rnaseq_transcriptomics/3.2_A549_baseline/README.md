# A549 RNA-seq reprocessing, 18 September 2026

Open [the analysis report](outputs/report/analysis.html). Complete numerical results are in `outputs/tables`, models in `outputs/models`, and publication exports in `outputs/figures` (PDF, SVG and 600-dpi PNG). Figure captions and their source-table manifest accompany the exports.

This is a new analysis of the existing featureCounts matrix, following the core HT1080 `SEPT_11_2026_REANALYZE` workflow. It does not repeat FASTQ alignment or counting. Original count files and earlier analyses are preserved.

## Design decision

The author stated on 18 September 2026 that A549 pairing and population-derivation history are uncertain. The four LOY and four sgCTRL libraries are therefore compared **without pairing**. Sample suffixes 1–4 are identifiers, not inferred passage matches. There is one documented sampling batch, so the HT1080 two-period model is not transferred. Independent population derivation is not asserted.

The primary outputs are descriptive group-mean effects with one-sample-omission sensitivity. Condition-only voom-limma results are exploratory: their statistical inference assumes independent residuals, which cannot be established from the available derivation history. Repeated sampling of one population, if later confirmed, would not constitute independent editing replication. The three parental libraries enter a separate contextual normalization and analysis only.

## Implemented analyses

- Exact comparison of all 11 source count tables with the matrix and assignment summaries; exact versioned GENCODE v49 annotation.
- Core filtering and TMM normalization on expressed autosomal/non-PAR X genes, matching HT1080 settings. Original counts are reads because the featureCounts v2.1.1 commands lack `--countReadPairs`.
- Library composition, PCA, sample correlations, descriptive effects, unpaired voom-limma and mitochondrial-covariate sensitivity.
- Refitted leave-one-sample-out normalization, RLE and logCPM-prior sensitivity.
- Separate Y/PAR abundance and X/Y paralog expression summaries using the same core scaling.
- A secondary PAR1/PAR2 follow-up requested after the initial primary results: an inclusive expressed-nuclear voom fit, both inclusive and PAR-family BH, all PAR effect estimates and per-sample plots. This does not alter the primary non-PAR test family.
- Full frozen MSigDB 2026.1.Hs resource coverage, effect-ranked fgsea, fixed- and estimated-correlation CAMERA, global and within-collection BH, and Hallmark ribosomal-protein sensitivity.
- CollecTRI target-pattern estimates, PROGENy footprint summaries and descriptive exact-gene comparison with HT1080.

The 0.5 logCPM candidate cutoff plus stable omission direction is a descriptive screening rule, not an FDR threshold or proof of replication. Gene-set P values address a different null from the gene-level tests. RNA abundance does not establish chromosome structure or functional paralog compensation.

DNA/RNA coordinate integration is deferred pending verified assemblies and specimen correspondence. HT1080-specific rescue, HDAC and other mechanistic extensions are not part of this baseline A549 analysis.

## Reproduce

From this directory in PowerShell:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' --vanilla run.R
```

Use `--core-only` for expression/QC, `--enrichment-only` to resume enrichment from verified core inputs, or `--report-only` to regenerate the report from saved models. Resuming checks raw files, annotations, the sample manifest, comparison inputs, frozen resources and configuration. If those change, rerun the full analysis.

The resource cache is copied from HT1080 with SHA-256 manifests. `prepare_workspace.py` records workspace preparation; it is not needed to rerun an already prepared analysis. R scripts and configuration determine all results. `provenance` holds run logs, scientific validation, code hashes, summary statistics and package/session versions. Input hashes are in `outputs/tables/input_hashes.csv` and the source count audit. The original HT1080 scripts used as references are retained under `provenance/reference_code`.

The Windows environment may report unsupported `C.UTF-8` locale settings on R startup. Completion is determined from the native R exit code and `provenance/validation.json`, not PowerShell's formatting of warnings as error-stream records.
