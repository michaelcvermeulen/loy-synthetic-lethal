# loy_reanalysis_v1 — how to run

A fresh Y-chromosome-loss analysis of DepMap 26Q1 male cell lines, built from the raw tables in this collection. Results and the report are in `results/loy_reanalysis_v1/`; tables are in `derived/loy_reanalysis_v1/`.

Run every command from the collection root, `/global/home/hpc5173/projects/PRISM/data/CCLE/26Q1`.

## One stage

```bash
bash scripts/loy_reanalysis_v1/submit_stage.sh <stage_name> <script> [args]
```

This copies the script directory to `results/loy_reanalysis_v1/code_snapshots/<stage>_XXXXXXXX/`, submits the copy to SLURM (`--account=def-hpcg1766_cpu`), and appends the job to `results/loy_reanalysis_v1/submission_ledger.tsv`. Do not edit `submit_stage.sh` while a submission is waiting; bash reads it as it runs.

| Variable | Effect |
|---|---|
| `LR_TIME`, `LR_MEMORY`, `LR_CPUS` | Resources (default 2 h, 16 G, 4 CPUs) |
| `LR_WAIT=1` | Block until the job ends and return its exit status |
| `LR_ARRAY=1-N` | Run as an array; the script reads `SLURM_ARRAY_TASK_ID` |

Scripts named `*_forest_*.R` run in the `main` conda environment (R 4.2.2, ranger); everything else runs in `scRNA` (R 4.5.3, limma, edgeR, fgsea, lme4).

## Order

| Step | Script | Needs | Resources used |
|---|---|---|---|
| Environment check | `00_smoke_test.R` | — | 2 CPU, 8 G |
| Input checksums and dimensions | `01_inputs_check.R` | — | 8 CPU, 16 G |
| Y evidence, calls, validation | `02_extract_y_evidence.R`, `03_call_y_status.R`, `04_validate_calls.R` | — | 4 CPU, 24 G |
| Cohort, contexts, matrices | `05_build_cohort.R`, `06_prepare_outcomes.R` | calls | 4 CPU, 32 G |
| Freeze | `07_freeze.R` (`--replace` to archive and redo) | cohort, matrices | 2 CPU, 16 G |
| Known-biology controls, calibration | `08_pipeline_controls.R` | freeze | 8 CPU, 24 G |
| Pan-cancer genes, drugs | `10_gene_pan.R`, `11_drug_pan.R` | controls | 8 CPU, 32 G |
| Context interactions | `12_context_genes.R`, `13_context_drugs.R` | pan-cancer results | 8 CPU, 32 G |
| Expression | `14_expression_de.R`, then `15_expression_dependency.R` | context gene hits | 4–8 CPU, 48 G |
| Forests | `16_forest_y_status.R <dataset>`; `17_forest_expression.R array` then `collect`; `18_forest_confounders.R` | stage 15 | 8 CPU, 24 G |
| Robustness, comparison | `19_hit_robustness.R`, `20_compare_prior.R` | all scans | 8 CPU, 32 G |
| What the cancer-type tests could detect | `24_power_context.R` (`resume` skips the long simulation and reads it back) | context gene screen, controls | 16 CPU, 48 G, about 80 min; 32 CPU, about 45 min |
| Public inputs for stages 25 to 30 | `bash scripts/loy_reanalysis_v1/fetch_extension_inputs.sh` (not a SLURM stage; needs internet; refuses to overwrite) | — | about 355 MB |
| Organoids | `25_organoid_prepare.R`, then `26_organoid_scan.R` | downloads, freeze | 4 CPU, 24 G; about 1 and 2 min |
| RNAi context screen | `27_rnai_context.R` | freeze, context definitions | 8 CPU, 32 G; about 2 min |
| Cell lines with organoids | `28_combine_sources.R` | stages 26 and 27, stored gene scan and context screen | 4 CPU, 16 G; under 1 min |
| Primary AML | `29_beataml_prepare.R`, then `30_beataml_drugs.R` | downloads | 2 CPU, 16 G; then 8 CPU, 48 G, about 1 min |
| Inputs for stages 31 to 33 | `bash scripts/loy_reanalysis_v1/stage_nextgen_inputs.sh <folder with the downloads>` (not a SLURM stage; copies six browser downloads, verifies each, refuses to overwrite) | — | about 320 MB |
| DepMap NextGen | `31_nextgen_prepare.R`, then `32_nextgen_scan.R` | NextGen inputs, both earlier freezes | 4 CPU, 16 G, under 1 min; 8 CPU, 32 G, about 2 min |
| Three CRISPR sources | `33_combine_three.R` | stages 26, 28 and 32 | 4 CPU, 24 G; about 4 min |
| Driver status inside cancer types | `34_driver_type_prepare.R`, then `35_driver_type_tests.R`, then `36_driver_type_power.R` | the frozen cohort and context panel, stages 12, 24, 25, 26 and 31 | 4 CPU, 24 G, under 1 min; 16 CPU, 48 G, about 4 min; 16 CPU, 64 G, about 1 hour (`LR_TIME=04:00:00`) |
| Figures, report, checks | `21_figures.R` (`extension` draws only Figures 11 and 12, `nextgen` only Figure 13, `driver_type` only Figure 14; any other argument stops the script), `22_render.R`, `23_validate_outputs.R` | everything | 4 CPU, 16 G |

`16_forest_y_status.R` takes one of `crispr_effect`, `prism_REPPRIMARY`, `prism_REP1M`, `prism_REP300`, `prism_secondary`, `gdsc1`, `gdsc2`, `ctrp`. `17_forest_expression.R array` needs `LR_ARRAY=1-<rows of expression_dependency_set.csv>`.

## Shared files

| File | Holds |
|---|---|
| `parameters.json` | Every threshold, gate and seed. Hashed by the freeze; later stages refuse to run if it changes. |
| `common.R` | Paths, readers, robust and classical single-coefficient fits, permutation and FDR helpers |
| `scan_helpers.R` | limma scan with donor blocking, design builder, fast permutation statistic, fgsea wrapper |
| `context_helpers.R` | Cell counts, gates and FDR families for the context stages |
| `power_parameters.json` | Simulation sizes, shift grids, seeds and the benchmark list of stage 24. Not part of the freeze; its hash is written to `power_summary.json`. |
| `power_helpers.R` | Critical p-values, simulated genes, the screen's two fits, detection curves, normal-theory limits |
| `extension_parameters.json` | Seeds, gates and thresholds of stages 25 to 36, the matched cancer types, the counts asserted before any test, the 12 public inputs of stages 25 to 30 with URL and checksum (`inputs`) and the six NextGen files (`nextgen_inputs`). Hashed by section in the four extension freezes; the NextGen settings are the section `nextgen`, those of stages 34 to 36 the section `driver_type`. |
| `extension_helpers.R` | Guarded writers for `derived/loy_reanalysis_v1/extension/`, the second freeze, Y calls on another dataset's scale, the Y-loss fraction, the karyotype reader, weighted limma, context fits with missing lines, the combination across sources, permutation tests, follow-up checks, the report readers stage 23 uses, and the helpers of stages 34 to 36 (gene status, census of splits, the split model, simulated genes on enlarged copies of a design, lines needed) |
| `fetch_extension_inputs.sh` | Downloads the public inputs of stages 25 to 30 and verifies each against the publisher's checksum |
| `stage_nextgen_inputs.sh` | Copies the NextGen inputs of stages 31 to 33 from a folder of browser downloads and verifies each against the checksum recorded in `nextgen_inputs` |
| `forest_helpers.R` | Boruta-style selection with stratified shadows, donor-aware folds |
| `plot_style.R` | House figure style (Open Sans; Y-positive `#3C5488`, Y-negative `#E64B35`) |
| `tests/` | Unit checks for the five helper files; run each with the matching `Rscript --vanilla` |

## Tests

```bash
/global/home/hpc5173/miniconda3/envs/scRNA/bin/Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_common.R
/global/home/hpc5173/miniconda3/envs/scRNA/bin/Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_scan_helpers.R
/global/home/hpc5173/miniconda3/envs/scRNA/bin/Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_power_helpers.R
/global/home/hpc5173/miniconda3/envs/scRNA/bin/Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_extension_helpers.R
/global/home/hpc5173/miniconda3/envs/main/bin/Rscript --vanilla scripts/loy_reanalysis_v1/tests/test_forest_helpers.R
```

## Stage 24 trial runs

`24_power_context.R` takes about 80 minutes on 16 CPUs (about 32 ms per simulated gene, 51,500 genes per cancer type). To try a change first, copy `power_parameters.json`, shrink `shift_genes`, `benchmark_replicates` and `genotype_genes`, add `"only_contexts": ["lineage__Lung", "driver__KRAS"]`, and point the stage at the copy and at a scratch folder:

```bash
LOY_POWER_PARAMETERS=/path/to/trial.json LOY_POWER_OUTPUT=/path/to/scratch \
  /global/home/hpc5173/miniconda3/envs/scRNA/bin/Rscript --vanilla scripts/loy_reanalysis_v1/24_power_context.R
```

Nothing in `derived/` is written in a trial run.

## Stages 25 to 30: organoids, RNAi and primary AML

Added 2026-10-04. These stages screen three public sources for vulnerabilities that go with Y loss. They read the frozen cohort and the stored cell-line results and write only to `derived/loy_reanalysis_v1/extension/`; `write_ext` refuses any other folder and any table name without one of the prefixes `organoid_`, `rnai_`, `combined_`, `beataml_`, `extension_`.

- **Order.** 25 before 26; 26 and 27 before 28; 29 before 30. Stage 21 then draws Figures 11 and 12, and stages 22 and 23 render and check.
- **Second freeze.** Stages 25 and 29 each write a freeze (`extension_freeze_organoid.json`, `extension_freeze_beataml.json`) with the hash of the parameter sections they depend on and of the tables they prepared. Stages 26, 28, 30, 22 and 23 stop if either has changed. `--replace` on a preparation stage archives the old freeze and writes a new one.
- **Counts fixed in advance.** `extension_parameters.json` holds the group counts seen before any test (`expected_counts`); the preparation stages stop if a count differs by more than `count_tolerance`.
- **Parameters.** `parameters.json` and `power_parameters.json` are not touched by these stages. Each stage writes the hash of its own sections to its summary file, and stage 23 compares all six.
- **Rank-based tests carry no weights.** The organoid fits on the raw scale weight each organoid by its noise; the rank-based tests (cancer type, genotype panel, their follow-up) do not, because weights inflated them under shuffled Y status.
- **Trial runs.** Copy the extension folder to a scratch folder and point a stage at it. `LOY_EXTENSION_TRIAL=1` cuts stages 26 and 27 to a few hundred genes and stage 30 to 2,000 permutations; without it the stage runs in full, which is how a rerun is compared with the stored tables. Nothing in `derived/` is written.

```bash
mkdir -p /path/to/scratch && cp -p derived/loy_reanalysis_v1/extension/* /path/to/scratch/
LOY_EXTENSION_OUTPUT=/path/to/scratch LOY_EXTENSION_TRIAL=1 \
  /global/home/hpc5173/miniconda3/envs/scRNA/bin/Rscript --vanilla scripts/loy_reanalysis_v1/26_organoid_scan.R
```

- **Arithmetic differs between nodes.** A fit on a node with AVX-512 and one without differ in their last digits (about 1e-10 of a number's size). Tables are therefore reproducible in every label, count and call, and in their numbers to that precision, but not always byte for byte. Stage 23 keeps its small refit tables as stored when the recomputed values agree to 1e-8 (`write_result_stable`), so the "nothing earlier changed" check does not depend on the node.

## Stages 31 to 33: DepMap NextGen as a third CRISPR source

Added 2026-10-05. Stage 31 prepares the NextGen cohort, stage 32 screens it alone, and stage 33 combines it with the cell lines and the Sanger organoids. They follow the rules of stages 25 to 30 (own prefix `nextgen_`, tables only in `derived/loy_reanalysis_v1/extension/`, nothing earlier rewritten) and add these:

- **Order.** 31 before 32, 32 before 33; stage 33 also needs the stored results of stages 10, 12, 26 and 28. `21_figures.R nextgen` then draws Figure 13, and stages 22 and 23 render and check.
- **Third freeze.** Stage 31 writes `extension_freeze_nextgen.json` over the `shared` and `nextgen` parameter sections, the cohort, the Y calls, both matrices, the gene annotation, the cancer types to pool, the chance of replication and the test counts. Stages 32, 33, 22 and 23 stop if any of them changed. No NextGen gene was compared with Y status before the first freeze existed. That freeze was replaced once (`31_nextgen_prepare.R --replace`), after the review of 2026-10-05, to correct which Sanger organoid one shared model is paired with and to add two settings; the first is archived beside it, and stage 23 proves that nothing else in it changed.
- **The four paired genes are named.** DDX3X, EIF1AX, ZFX and RPS4X form a family of their own in NextGen (`named_genes`). Any other gene needs FDR below 0.05, an effect of at least 0.10 Chronos units and a p-value below a bar from 300 shuffles of Y status within stratum (`null_runs`, `family_wise_level`).
- **Permutation p-values decide.** The weighted model's p-values run small in NextGen, so the four paired genes are judged on permutation p-values from a million shuffles (`paired_permutations`), and each combined hit of stage 33 is followed up with a million joint shuffles (`combined_permutations`). `followup_statistics` draws them 20,000 at a time, so memory stays flat.
- **One stratum per lineage and screen format.** `build_design` takes its dummies from `OncotreeLineage`, so the NextGen frames put the stratum label there and keep the lineage in `lineage`.
- **Trial runs never show a real result.** `LOY_EXTENSION_TRIAL=1` cuts stage 32 to a few hundred genes and 20 shuffles and shuffles Y status within stratum before any fit; stage 33 then reads that trial and shuffles NextGen's Y status in its pooled tests and in the follow-up of combined hits. A trial needs the outputs of stage 31 in the same output folder:

```bash
mkdir -p /path/to/scratch
LOY_EXTENSION_OUTPUT=/path/to/scratch bash scripts/loy_reanalysis_v1/submit_stage.sh nextgen_prepare_dry 31_nextgen_prepare.R
LOY_EXTENSION_OUTPUT=/path/to/scratch LOY_EXTENSION_TRIAL=1 bash scripts/loy_reanalysis_v1/submit_stage.sh nextgen_scan_trial 32_nextgen_scan.R
LOY_EXTENSION_OUTPUT=/path/to/scratch LOY_EXTENSION_TRIAL=1 bash scripts/loy_reanalysis_v1/submit_stage.sh nextgen_combine_trial 33_combine_three.R
```

  The scratch folder must be on a filesystem the compute nodes see (`/global/scratch/...`, not `/tmp`). Without `LOY_EXTENSION_TRIAL` the same three commands rerun the stages in full into the scratch folder, which is how a rerun is compared with the stored tables.
- **Ties in permutation counts.** Inside a small cancer type the same arrangement of Y status comes up many times in 20,000 shuffles, and each time it ties with the observed one. Stages 32 and 33 count a shuffle within 1e-9 of the observed statistic as a tie (`tie_tolerance`), so the count does not depend on the processor. Stages 26 and 30 keep the plain comparison they ran with.
- **Reproducibility.** Stages 32 and 33 were rerun on AMD nodes (`SBATCH_CONSTRAINT=amd`, jobs 12306255 and 12306190) against stored tables from an Intel AVX-512 node: over all 28 NextGen tables every label, count, call and missing value is identical, and no number differs by more than 1.2e-10 of its size. Two things had to be fixed for that to hold. Permutation counts now treat a shuffle within 1e-9 of the observed statistic as a tie (`tie_tolerance`), and a dependency-call contrast carries no p-value where the calls leave no Y contrast.
- **Figure 13 only.** Draw it with `21_figures.R nextgen`. A run without an argument redraws every figure, and five of the first twelve place points or labels without a seed, so their bytes would change.

## Stages 34 to 36: driver status inside cancer types

Added 2026-10-06. The context screen (stage 12) tests each driver alteration across all cancer types at once. These stages ask the finer question inside one cancer type, mainly to show what the collection can and cannot test at that level. Stage 34 builds one status per gene and the list of splits, stage 35 runs the tests, and stage 36 simulates what each split could detect. They follow the rules of stages 25 to 33 (own prefix `driver_type_`, tables only in `derived/loy_reanalysis_v1/extension/`, nothing earlier rewritten) and add these:

- **Order.** 34 before 35, 35 before 36. `21_figures.R driver_type` then draws Figure 14, and stages 22 and 23 render and check.
- **One status per gene.** A gene is altered in a line with any of a driver-flagged mutation, an amplification or a deep deletion, with the definitions, thresholds and gene lists of `05_build_cohort.R` and `parameters.json`. Stage 34 rebuilds the statuses from the raw mutation and copy-number tables and stops unless they equal every alteration context that the frozen context matrix kept.
- **Split.** A gene in a cancer type (lineage, disease or subtype) with five altered lines per Y group. The test of a difference by driver status also needs five unaltered lines per Y group. A split with the same lines at two levels is fitted once.
- **Fourth freeze.** Stage 34 writes `extension_freeze_driver_type.json` over the `shared` and `driver_type` parameter sections, both status tables, the census, the splits and their lines for cell lines and organoids, the pairs to combine and the gene families. Stages 35, 36, 22 and 23 stop if any of them changed. No gene was compared with Y status inside a split before that freeze existed, apart from one slip in a trial run that the progress ledger records.
- **One model per split.** `fit_split` fits rank-normalised score ~ library + driver status + Y among altered lines + Y among unaltered lines on all lines of the cancer type, with limma's moderated t as in the screen. The context screen fits altered lines alone; here that would let screening library take up to three quarters of a split's information. The existing fit helpers (`fit_context_tests`, `context_fit_parts`, `fit_context_patterns`) cannot fit a split without unaltered lines and are left as they are; a unit check ties the new helpers to them on a split with four cells.
- **Controls first.** Stage 35 tests KRAS by KRAS status and MDM2 by TP53 status inside single cancer types before any Y result and stops if they fail.
- **Shuffles.** 200 runs with Y status shuffled within driver status inside each split (`null_runs`), each a full run.
- **Least squares beside the rule.** limma's moderated t takes a gene's variance almost wholly from the other genes at these sizes. `fit_split_ols` gives each test by plain least squares on the same ranks (`inside_p_ols`); those p-values describe and make no call.
- **Power.** Stage 36 reads its grids, host-gene numbers and benchmarks from `power_parameters.json`, so its limits compare with those of stage 24. It reproduces the stored p-values of every split before it simulates, by refitting and by sending real genes through its own simulation path. Its limits are limits of the layer's rank-based test; `family_raw` gives the same bar for least squares on Chronos scores, which has no ceiling and is not calibrated in small cells. "Lines needed" comes from `enlarged_split_tests`: the design copied 2, 4, 8 and 16 times, each copy carrying the spread of a different real gene, with the bar and the variance prior held at today's values. A known effect is sized two ways, over its own gene's spread and over a typical gene's (`convention`). The check on real lines cuts a real split to 1/k, copies it k times and sets it against the real split.
- **Trial runs never show a real result.** `LOY_EXTENSION_TRIAL=1` cuts stage 35 to a few hundred genes and 5 shuffles and shuffles Y status within each split's cells before any fit, in every variant of the follow-up as well; stage 36 then reads that trial and simulates four splits. A trial needs the outputs of stage 34 in the same output folder:

```bash
mkdir -p /path/to/scratch
LOY_EXTENSION_OUTPUT=/path/to/scratch bash scripts/loy_reanalysis_v1/submit_stage.sh driver_prepare_dry 34_driver_type_prepare.R
LOY_EXTENSION_OUTPUT=/path/to/scratch LOY_EXTENSION_TRIAL=1 LR_CPUS=16 LR_MEMORY=48G bash scripts/loy_reanalysis_v1/submit_stage.sh driver_tests_trial 35_driver_type_tests.R
LOY_EXTENSION_OUTPUT=/path/to/scratch LOY_EXTENSION_TRIAL=1 LR_CPUS=16 LR_MEMORY=48G bash scripts/loy_reanalysis_v1/submit_stage.sh driver_power_trial 36_driver_type_power.R
```

  Without `LOY_EXTENSION_TRIAL` the same three commands rerun the stages in full into the scratch folder, which is how a rerun is compared with the stored tables.
- **Reproducibility.** Stages 34 to 36 were rerun on AMD nodes (`SBATCH_CONSTRAINT=amd`, jobs 12315272, 12316474 and 12315290) against stored tables from Intel AVX-512 nodes: over all 31 `driver_type_` tables every label, count, call, flag and missing value is identical and the rows come in the same order, and no number differs by more than 1.5e-9 of its size. The one exception is `path_gap` in `driver_type_power_splits`, a difference between two nearly equal numbers, where two values differ by 7.5e-12. One thing had to be fixed for that to hold. Head and neck at lineage level and head-and-neck squamous cell carcinoma at disease level give the same test among TP53-altered lines, their p-values tie to the last bits, and which row came first depended on the processor. The follow-up takes the seed of its shuffles from the row, so one permutation p-value differed by one shuffle in 20,000. Tied tests now keep the order of the split table (`tied_order`, with a unit check). In the stored tables this moved four rows of `driver_type_paired_genes` and 40 of `driver_type_hit_followup` and changed that one p-value from 2 to 1 in 20,001 (stage 35 rerun as job 12316457); no flag, count or reported number changed. One thing was left as it is: `driver_type_splits.csv` and `driver_type_organoid_splits.csv` differ between processor types in the last digit of the columns that give the information a design keeps, so the hashes in the freeze are those of the Intel run. A rerun of stage 34 on another processor type belongs in a scratch folder, as above; in place it would rewrite those two tables and then stop at the freeze.
- **Figure 14 only.** Draw it with `21_figures.R driver_type`.
