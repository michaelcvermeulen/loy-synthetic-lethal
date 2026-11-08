# loy_reanalysis_v1 — export of 2026-10-07

A copy of the Y-chromosome-loss analysis of DepMap 26Q1 as it stood on 2026-10-07: gene dependencies, drug sensitivities, the context screens, expression, the random forests, the power analysis (stage 24), the organoid, RNAi and primary-AML stages (25 to 30), DepMap NextGen as a third CRISPR source (31 to 33) and driver status inside cancer types (34 to 36).

It was taken from `/global/home/hpc5173/projects/PRISM/data/CCLE/26Q1` on the cluster. The folders keep the layout they have there, so every path written in the report, the READMEs and the scripts resolves inside this folder.

It replaces the export of 2026-10-05. What changed since then is listed at the end.

## Where to look

| For | Open |
|---|---|
| The results, written up | `results/loy_reanalysis_v1/report.html` or `report.pdf`. Both stand alone. |
| The headline numbers and how to read them | `results/loy_reanalysis_v1/README.md` |
| Figures for the thesis | `results/loy_reanalysis_v1/figures/`: fourteen figures, each as PDF, SVG and 300 dpi PNG |
| The data behind each figure panel | `results/loy_reanalysis_v1/figure_source_data/` |
| Result tables | `derived/loy_reanalysis_v1/final/`, each as CSV and RDS. `docs/LOY_REANALYSIS_OUTPUTS.md` says what each table holds, one row per what. |
| Methods | `scripts/loy_reanalysis_v1/` (`README.md` and the three parameter files), and `results/loy_reanalysis_v1/session_info_*.txt` for R and package versions |
| Decisions taken during the run | `results/loy_reanalysis_v1/progress_ledger.md` (each "Ruling" line) |
| Which job produced what | `results/loy_reanalysis_v1/submission_ledger.tsv`, `code_snapshots/`, `logs/` |
| Data sources and licences | `docs/DATA_GUIDE.md` |

`derived/loy_reanalysis_v1/final/` repeats tables that are also in `derived/loy_reanalysis_v1/` and its `extension/` folder. The copies are identical; the scripts read the originals and `final/` is the set meant for handing on.

`derived/loy_reanalysis_v1/matrices/` holds the outcome matrices the analysis used (models by features, as R matrices): CRISPR gene effect and dependency probability, protein-coding expression, RNAi, and the PRISM, GDSC and CTRP drug responses. They are subsets of the public source data, so cite the sources in `docs/DATA_GUIDE.md` when you use them.

## What is not here

- **The raw source data** (`raw/`, about 19 GB) and the papers and supplements in `references/`. The 26 input tables of DepMap 26Q1 and the drug panels are listed with size and SHA-256 in `derived/loy_reanalysis_v1/input_manifest.csv`. The 12 public inputs of stages 25 to 30 and the 6 of stages 31 to 33 are listed with source and checksum in `scripts/loy_reanalysis_v1/extension_parameters.json` (sections `inputs` and `nextgen_inputs`), and again as checksum lists in `docs/SHA256SUMS_extension_inputs` and `docs/SHA256SUMS_nextgen_inputs`. Stages 25 and 34 also read supplementary tables 2 to 6 of the Sanger organoid study (`raw/studies/sanger_organoid/supplementary_tables/`), which were in the collection before this analysis; they are listed in `docs/SHA256SUMS`, the collection's own checksums of its original 197 files. Stages 34 to 36 read no input beyond these. These checksums record what was on the cluster; they do not verify the files against their publishers.
- **The September analyses** (`loy_dependencies_v1`, `y_dosage_vulnerabilities_v1`, `rna_y_vulnerabilities_v1`) and the earlier sex-chromosome calling work. This analysis compares itself with them in the `comparison_*` tables and in the report. One piece of the earlier work is included because the plotting code reads it: the figure style in `results/rna_y_vulnerabilities_v1/style_snapshots/thesis-figure-style-v0.3/`.
- The font cache `results/loy_reanalysis_v1/.fontconfig/`, which the plotting code rebuilds.

## Using it on another machine

Reading tables needs nothing special. The RDS tables are `data.table` objects, which also behave as data frames:

```r
hits <- readRDS("derived/loy_reanalysis_v1/final/gene_pan_results.rds")
effect <- readRDS("derived/loy_reanalysis_v1/matrices/crispr_effect.rds")  # models by genes
```

To run the scripts against this copy, point them at it. Without this variable they look for the cluster path.

```bash
export LOY_REANALYSIS_ROOT=/path/to/loy_reanalysis_v1_export_2026-10-07
Rscript --vanilla "$LOY_REANALYSIS_ROOT/scripts/loy_reanalysis_v1/21_figures.R"
```

- **Figures (stage 21)** need only what is in this folder. This was tested on the cluster on 2026-10-07: the archive was unpacked to a scratch folder and stage 21 was run against it. The outcome of that test, and the R and package versions it ran with, are in `EXPORT_CHECKS.txt`. It has not been run on any other machine. Without an argument the stage draws all fourteen figures. With `extension`, `nextgen` or `driver_type` it draws Figures 11 and 12, Figure 13 or Figure 14 alone and leaves the others as they are. The stage loads `data.table`, `jsonlite`, `digest`, `ggplot2`, `patchwork`, `ggrepel`, `svglite`, `ragg` and `systemfonts`. It finds Open Sans through fontconfig, which is how Linux and macOS locate fonts; on Windows, expect to install Open Sans or adjust `plot_style.R`.
- **The report (stage 22)** also reads only files in this folder, but rendering it was not tested from this copy. It needs pandoc, and for the PDF a LaTeX installation with XeLaTeX and the DejaVu Serif font. Stage 22 also rewrites `derived/loy_reanalysis_v1/final/`.
- **Every other stage** reads `raw/` and will stop without it.
- `submit_stage.sh` and `job.sbatch` are SLURM wrappers for the cluster. Call `Rscript` directly instead.
- Stage 21 overwrites `results/loy_reanalysis_v1/figures/` and `figure_source_data/` in the folder it is pointed at. Keep an untouched copy if you change the plotting code.

## Checking the copy

`MANIFEST.sha256` lists the SHA-256 of every file here except itself and `EXPORT_CHECKS.txt`. From this folder:

```bash
sha256sum -c --quiet MANIFEST.sha256
```

On Windows, run that in Git Bash or WSL, or compare single files with `Get-FileHash`.

`EXPORT_CHECKS.txt` records what was checked when the archive was built. The analysis's own checks are in `results/loy_reanalysis_v1/validation_summary.json` (115 checks, all passed when last run on 2026-10-06). All of these are checks on computation and file integrity. None of them is biological validation.

## What changed since the export of 2026-10-05

- **One file from the collection was added that the earlier export left out:** `docs/SHA256SUMS`, because it lists inputs that the other checksum lists do not.
- **Nothing that export held was removed, and its result tables are unchanged.** Of its 4,255 files, 4,237 are byte-identical here. That includes every table under `derived/`, the 36 files of Figures 1 to 12 and their 34 source-data files.
- **18 files changed.** The report (`report.html`, `report.pdf`, `report.Rmd`) and `22_report.Rmd` gained two sections. The other 14 are documentation, ledgers, `validation_summary.json`, and the scripts and helpers that the new stages extend (`21_figures.R`, `22_render.R`, `23_validate_outputs.R`, `extension_helpers.R` and its tests, `extension_parameters.json`).
- **New: DepMap NextGen as a third CRISPR source** (stages 31 to 33, added 2026-10-05): the `nextgen_*` tables, Figure 13 and the report section "DepMap NextGen models".
- **New: driver status inside cancer types** (stages 34 to 36, added 2026-10-06): the `driver_type_*` tables, Figure 14 and the report section "Driver status inside cancer types".
- The headline did not change: DDX3X, EIF1AX, ZFX and RPS4X are the only Y-loss dependencies, and no drug, cancer-type or driver-specific effect passes.
