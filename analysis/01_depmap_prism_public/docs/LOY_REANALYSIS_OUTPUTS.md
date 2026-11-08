# loy_reanalysis_v1 — output inventory and provenance

Added 2026-10-02. A fresh Y-chromosome-loss analysis of DepMap 26Q1 male cell lines. Source data in `raw/` was read only; nothing there was edited, moved or deleted. These outputs are later additions and are not part of the original 197-file relocation manifest or `docs/SHA256SUMS`.

## Where things are

| Path | Contents |
|---|---|
| `scripts/loy_reanalysis_v1/` | R code, `parameters.json`, SLURM wrappers, unit tests, run instructions (`README.md`) |
| `derived/loy_reanalysis_v1/` | Calls, cohort, context panel, per-stage tables |
| `derived/loy_reanalysis_v1/final/` | Deliverable tables as CSV and RDS |
| `derived/loy_reanalysis_v1/extension/` | Tables of stages 25 to 30 (organoids, RNAi, combination, BeatAML2) of stages 31 to 33 (DepMap NextGen) and of stages 34 to 36 (driver status inside cancer types), their four freezes and the three baselines; see "Addition 2026-10-04", "Addition 2026-10-05" and "Addition 2026-10-06" |
| `derived/loy_reanalysis_v1/matrices/` | Outcome matrices (models by features) as `.rds` |
| `derived/loy_reanalysis_v1/join_audit/` | Identifiers that did not match the 26Q1 model registry, per input |
| `derived/loy_reanalysis_v1/reference/gene_sets.rds` | Snapshot of the MSigDB collections used |
| `results/loy_reanalysis_v1/` | Report, figures, source data, logs, code snapshots, ledgers |

## Inputs

Twenty-six source tables, listed with path, size, SHA-256 and dimensions in `derived/loy_reanalysis_v1/input_manifest.csv`. All are under `raw/` except the CTRP mapping in `derived/ctrp/`, which was spot-checked against the raw CTRP curves (`ctrp_derived_vs_raw.csv`). Checksums are compared with `docs/SHA256SUMS`, which is a local relocation check and not publisher verification.

| Input | Used for |
|---|---|
| `raw/depmap_public_26q1/metadata/Model.csv`, `ModelCondition.csv`, `Gene.csv`, `PortalCompounds.csv` | Registry, sex, lineage, donors, gene locations, compound key |
| `raw/depmap_public_26q1/expression/` (protein-coding TPM, all-gene TPM, raw counts) | RNA Y call, XIST, differential expression, expression associations |
| `raw/depmap_public_26q1/copy_number/` (WGS genes, WES genes, WGS segments) | DNA Y call, copy-number contexts |
| `raw/depmap_public_26q1/crispr/` (effect, dependency probability, confounders) | Gene dependencies, library and screen covariates |
| `raw/depmap_public_26q1/signatures/` (global signatures, inferred molecular subtypes) | Genome doubling, MSI, CIN, ploidy, aneuploidy, LOH, molecular subtypes |
| `raw/depmap_public_26q1/mutations/` (long table, hotspot matrix) | Driver contexts, TP53 status, set of mutation-profiled models |
| `raw/rnai/combined_demeter2/` | RNAi replication |
| `raw/drug_response/prism/24Q2/primary/`, `prism/20Q2/secondary/`, `gdsc/`, `derived/ctrp/` | Drug responses |

Two references outside `raw/` are read: the figure style in `results/rna_y_vulnerabilities_v1/style_snapshots/thesis-figure-style-v0.3/`, and, in the comparison stage only, the earlier analyses' result tables.

## Key tables in `derived/loy_reanalysis_v1/final/`

| Table | One row per | Notes |
|---|---|---|
| `y_status_calls` | model (2,154) | `y_call`, `y_tier`, RNA and DNA evidence, every alternative definition, sex flags |
| `cohort` | eligible male (1,076) | Calls plus covariates; `in_primary` marks the two-group comparison |
| `context_definitions` | candidate context (378) | Class, source, cell counts, whether kept |
| `gene_pan_results` | gene (17,205) | Effect, FDR, adjusted effect, dependency fractions, label |
| `gene_pan_continuous`, `gene_pan_rnai` | gene | Secondary scan on the continuous Y score; RNAi scan |
| `drug_pan_results` | dataset and compound (6,413) | Effect, robust error, FDR within dataset, adjusted effect |
| `drug_pan_consensus`, `drug_pan_mechanism` | shared compound; class | Cross-dataset mixed model; class-level permutation test |
| `context_gene_screen`, `context_drug_screen` (`.csv.gz`, `.rds`) | outcome and context | Interaction and inside-context effects, rank-based p-values, raw-scale estimates |
| `context_gene_calibration`, `context_drug_calibration` | cell-size bin and run | Observed and shuffled-label yields |
| `expression_de_results` | gene (45,269 fitted) | `placed` and `is_y_gene` mark the test family |
| `expression_dependency_associations` (`.csv.gz`, `.rds`) | dependency and gene | Adjusted, within-group and interaction associations |
| `expression_dependency_y_genes` | dependency and Y gene | Marginal, joint and within-Y-positive effects |
| `forest_y_status_features`, `forest_expression_*`, `forest_confounders` | feature or outcome | Forest results with seed stability and null runs |
| `robustness_*`, `comparison_*` | check | Robustness of hits; comparison with the earlier analyses |

Gene identifiers keep the full DepMap header (`SYMBOL (EntrezID)`) in the `feature` column, with the parsed symbol alongside.

## Provenance and validation

- Every job ran from a copy of the code in `results/loy_reanalysis_v1/code_snapshots/`, recorded in `submission_ledger.tsv` with its SLURM job ID.
- `derived/loy_reanalysis_v1/freeze.json` holds the hashes of the parameters, calls, cohort and context panel taken before any outcome was tested against Y status. One earlier freeze is archived beside it; the change (a run-time setting) is explained in `results/loy_reanalysis_v1/progress_ledger.md`.
- `progress_ledger.md` records every decision that departed from or completed the plan.
- Checks: input checksums and dimensions (`input_manifest.csv`); known-biology and calibration controls (`pipeline_controls.csv`, `pipeline_controls.json`); final independent checks (`results/loy_reanalysis_v1/validation_summary.json`).
- These are computational checks on public cell-line data. They are not biological validation.

## Addition 2026-10-03: what the cancer-type tests could detect (stage 24)

A power analysis of the context screen at three levels of cancer type (lineage, primary disease, subtype). It asks how large a Y-loss effect confined to one cancer type would have had to be for the screen to find it, and whether known vulnerabilities would have been found at the same sample sizes. No source data was added: the stage reads the frozen cohort, the context panel, the CRISPR matrix and the stored screen.

| Path | Contents |
|---|---|
| `scripts/loy_reanalysis_v1/24_power_context.R` | The stage |
| `scripts/loy_reanalysis_v1/power_helpers.R`, `tests/test_power_helpers.R` | Helpers and their unit checks |
| `scripts/loy_reanalysis_v1/power_parameters.json` | Simulation sizes, shift grids, seeds, the benchmark list. Not part of the freeze; its hash is in `power_summary.json`. |
| `results/loy_reanalysis_v1/figures/fig9_power_clean_shift.*`, `fig10_power_benchmarks.*` | The two figures, with source data in `figure_source_data/` |

Tables in `derived/loy_reanalysis_v1/`, copied to `final/`:

| Table | One row per | Notes |
|---|---|---|
| `power_census` | level (3) | Cancer types with lines, with both Y groups, and with at least 3, 5 and 10 lines per group |
| `power_cancer_types` | cancer type (190) | Counts, status (tested, tested under another name, below the gate, one Y group missing), exact-test floor, normal-theory limits |
| `power_shift_detection` | tested cancer type and shift (690) | Probability of finding a clean shift at each bar; shifts in SD units and in Chronos units |
| `power_limits` | tested cancer type, bar, unit and target (480) | Shift found 50% and 80% of the time, the ceiling, the normal-theory limit |
| `power_benchmarks` | benchmark (25) | Effect, spread and group sizes of each known vulnerability in the eligible male lines |
| `power_benchmark_detection` | tested cancer type, benchmark and design (3,000) | Probability of finding a Y-loss effect like the benchmark. Designs: `own gene`, `typical gene` (with and without the altered group's spread), `own gene, no effect` |
| `power_benchmark_summary`, `power_benchmark_null`, `power_benchmark_levels` | benchmark (by design and level; pooled; by cancer type) | Cancer types reaching 50% and 80%; false-positive rates; the benchmark gene's own level in each cancer type |
| `power_positive_controls` | tested cancer type and hit (120) | What the screen saw for DDX3X, EIF1AX, ZFX and RPS4X inside each cancer type, and what the simulation predicts with the effect in every Y-negative line (`predicted_everywhere_*`) or confined to the one cancer type (`predicted_*`) |
| `power_alone_hits` | cancer type and gene (17) | Every gene that passes FDR 0.05 with a cancer type analysed alone; `established` marks the four pan-cancer hits |
| `power_genotype_context_limits` | genotype context, bar and target | Clean-shift limits for the 82 driver, copy-number, molecular-subtype and genomic contexts |
| `power_drug_limits` | drug dataset and lineage | Lines per Y group, normal-theory limits, exact-test floor, BRAF-inhibitor benchmark |
| `power_summary.json` | — | Headline numbers, the stage's own checks, the parameter hash |

How to read them:

- A simulated gene counts as found when its p-value is below the value at which one more test would pass FDR 0.05 in the stored family. `candidate` is p below 0.05 for one gene named in advance; `one_type` uses the cancer type's own gene family ("alone" in the report); `class_inside` is the inside effect within the family of the context class; `as_run` is the screen's rule (`class_inside` and the interaction, the latter also within the class family).
- Each simulated gene is tested as the single gene added to the screen. The stage checks this fast computation against a full refit of the screen with one gene appended.
- Probabilities come from 400 or 500 simulated genes each and carry a Monte Carlo error of about 0.02.
- The simulation tables were written by job 12302772 and the remaining tables by a resumed run of the same stage; a full rerun under the same seeds is recorded in the progress ledger.
- These are limits of the design, computed by simulation on public cell-line data. They are not observed power, and they are not biological validation.

## Addition 2026-10-04: organoids, RNAi and primary AML (stages 25 to 30)

Three public sources screened for vulnerabilities that go with Y loss, because the cell lines are too few per cancer type: the Sanger organoid CRISPR screens (Herranz-Ors et al., Nature 2026), RNAi (DEMETER2) on the cell lines, and BeatAML2 drug response in primary AML. Each source is screened alone; cell lines and organoids are then combined per gene. The stages read the frozen cohort and the stored cell-line results and change none of them: all 459 earlier files (457 tables and matrices, two parameter files) keep the SHA-256 they had before the first new stage ran, which stage 23 checks against `extension/extension_baseline_sha256.txt`.

Source data added (read-only; details, licence and checksums in `docs/DATA_GUIDE.md` and `docs/SHA256SUMS_extension_inputs`):

| Location | Used for |
|---|---|
| `raw/studies/sanger_organoid/figshare/` (5 files, figshare 10.6084/m9.figshare.28339340, CC BY 4.0) | chrY copy number, MSI and WGD, sex and morphology, control genes, gene positions |
| `raw/drug_response/beataml2/` (7 files, `biodev/beataml2.0_data`, CC BY 4.0) | Clinical table with karyotypes, drug response curve fits, normalised RNA-seq, exome mutation calls, drug families, sample mapping, licence |

The organoid scores, expression, QC and annotation already in `raw/studies/sanger_organoid/` and the RNAi matrix in `raw/rnai/combined_demeter2/` are read as they were.

| Path | Contents |
|---|---|
| `scripts/loy_reanalysis_v1/25_organoid_prepare.R`, `26_organoid_scan.R` | Organoid Y calls, cohort and matrices; controls, scans, shuffles, follow-up |
| `scripts/loy_reanalysis_v1/27_rnai_context.R` | RNAi context screen |
| `scripts/loy_reanalysis_v1/28_combine_sources.R` | Cell lines with organoids; CRISPR against RNAi; what each test could detect |
| `scripts/loy_reanalysis_v1/29_beataml_prepare.R`, `30_beataml_drugs.R` | Y-loss fraction, karyotypes and cohort; inhibitor screen, families, follow-up |
| `scripts/loy_reanalysis_v1/extension_helpers.R`, `tests/test_extension_helpers.R` | Helpers and their 92 unit checks |
| `scripts/loy_reanalysis_v1/extension_parameters.json`, `fetch_extension_inputs.sh` | Settings, asserted counts and input checksums; the downloads |
| `results/loy_reanalysis_v1/figures/fig11_more_models.*`, `fig12_primary_aml.*` | The two figures, with source data in `figure_source_data/` (`fig11a` to `fig12d`) |

Tables in `derived/loy_reanalysis_v1/extension/`, copied to `final/` (CSV and RDS; the two largest as `.csv.gz` and RDS):

| Table | One row per | Notes |
|---|---|---|
| `organoid_y_calls` | organoid with data (256) | Sex by record and genome, RNA level of the seven Y genes, chrY copy number relative to ploidy, RNA call, DNA call, tier, the integrated call and five alternative definitions |
| `organoid_cohort` | screened organoid (162) | Screen quality, tumour type, library, MSI, WGD, TP53, fraction of genome altered, patient, `in_primary` (95: 47 Y-negative, 48 Y-positive) |
| `organoid_gene_annotation` | gene (16,562) | Chromosome, Y partner, common-essential, selective, in scan (16,547; the 15 chrY genes are out) |
| `organoid_weights` | malignant screened organoid (161) | The limma weight of each organoid, estimated without Y status |
| `organoid_controls` | control (3) | KRAS in KRAS-mutant, WRN in MSI, MDM2 by TP53 status |
| `organoid_pan_results` | gene (16,547) | Effect, Chronos-equivalent effect, p, FDR, adjusted estimates, depletion-call contrast, label |
| `organoid_continuous` | gene (16,547) | The same scan with Y level as a continuous score |
| `organoid_xy_partners` | X partner in the scan (10) | Their results as a family |
| `organoid_type_screen` | gene and cancer type (13,372) | Interaction with cancer type and the Y effect inside colorectal and oesophageal organoids; `selective` marks the 3,392 genes of the screen's family, the rest are cell-line selective genes kept for the combination |
| `organoid_genotype_screen` | selective gene and genotype group (23,744) | The secondary panel: TP53, APC, KRAS, PIK3CA, SMAD4, WGD, metastasis-derived |
| `organoid_hits`, `organoid_hit_followup` | hit (2); gene, context, check and variant (420) | The hits and their follow-up; `robust` applies the rule of stage 19 |
| `organoid_calibration` | run (observed and 5 shuffled) | Hit counts, tail fractions and lambda with Y status shuffled within cancer type |
| `organoid_test_counts`, `organoid_input_manifest` | context (10); input file (5) | Counts written before any test; inputs with checksums |
| `rnai_context_screen` (`.csv.gz`) | gene and context (240,509) | Interaction and inside-context tests on the lines with data for that gene, with the cell counts |
| `rnai_context_focused` | focus gene and context (691) | The 11 X-linked focus genes as their own family |
| `rnai_context_hits`, `rnai_context_inside_only` | call (0; 7) | Context-specific calls; detections inside a context without an interaction |
| `rnai_context_summary`, `rnai_context_calibration` | context (77); run and smallest cell | Yield per context; yield under five shuffles of Y status |
| `rnai_hit_followup` | detection, check and variant (63) | Follow-up of the seven inside-only detections |
| `combined_pan` | gene in both scans (16,252) | Each source's p and z, the pooled z, p and FDR, the difference between sources with its interval, the evidence class |
| `combined_within_type` | gene and cancer type (11,022) | The same for the Y effect inside colorectal and oesophagogastric models; `cell_inside_p_stored_model` is the p-value of the stored screen, `cell_inside_p` the one used (with MSI for colorectal) |
| `combined_cross_assay` (`.csv.gz`), `combined_cross_assay_summary` | gene and context tested in both assays (109,057); cancer-type flag (2) | CRISPR and RNAi side by side; counts nominal in both against chance |
| `combined_cell_line_leads_in_rnai` | cell-line lead (68) | Each lead of the cell-line screens looked up in RNAi |
| `combined_seen_before`, `combined_xy_partners` | gene (3); partner and source (42) | RPS4X, DDX3X and CORO1A across sources; the X partners across sources |
| `extension_power`, `beataml_power` | test and source (10; 1) | Models per Y group and the best-case effect found 80% of the time |
| `beataml_cohort` | specimen with RNA-seq (671) | Sex, guards, Y-loss fraction with and without the floor, karyotype evidence, covariates, `in_primary` (263 patients) |
| `beataml_karyotypes` | specimen (671) | The karyotype as read: cells, cells without Y, female clone, chimerism, parsed flag |
| `beataml_test_counts`, `beataml_families`, `beataml_compound_match` | inhibitor (161); inhibitor and family (376); matched inhibitor (81) | The gate counted before any fit; drug families; names matched to cell-line compounds |
| `beataml_controls` | FLT3 inhibitor (4) | Response by FLT3-ITD |
| `beataml_drug_results` | inhibitor (95) | Slope per complete Y loss on rank-normalised AUC, permutation p (the primary one), classical and robust p, FDR, adjusted estimates, label |
| `beataml_family_results`, `beataml_calibration` | family (22); run (6) | Family tests on the same permutations; yield with Y loss shuffled |
| `beataml_cell_line_comparison`, `_summary` | matched inhibitor (47); reference (2) | Slopes beside the cell-line estimates; the sign-agreement test |
| `beataml_hit_followup`, `beataml_ras_by_y_loss`, `beataml_leading_inhibitor` | check (15); Y-loss group (3); specimen (226) | Follow-up of the hit, including the two RAS-pathway checks added after it was seen; RAS-pathway mutations by Y loss; the data of Figure 12d |

Not copied to `final/`: the matrices (`organoid_effect.rds`, `organoid_effect_v11.rds`, `organoid_depleted.rds`, `beataml_auc.rds`), six stage summaries (`*_summary.json`), the two freezes, `extension_validation_organoid_refit` (stage 23's refit), and the baseline (`extension_baseline_sha256.txt`, `extension_baseline_report.html`).

How to read them:

- **Kinds of finding.** RPS4X, DDX3X and CORO1A were seen in an undocumented organoid comparison on disk before the screen was run, so in organoids they are replication, not discovery (`seen_before`). RNAi against CRISPR is replication in the same lines; organoids against cell lines, and BeatAML2 against cell-line drug data, are independent.
- **Organoid Y status.** RNA leads and DNA is the check. Ten calls rest on RNA with an unclear DNA result; three of those are Y-negative organoids that express no Y gene but keep most of the male-specific region in DNA. The definition `dna_confirmed` leaves all ten out.
- **Weights.** Fits on the organoid score carry a weight per organoid. Rank-based tests carry none: a first run that weighted them was inflated under shuffled Y status, and its two detections (DDX3X inside colorectal organoids, and the combined colorectal result) were withdrawn.
- **Combination.** Each source contributes a signed z from its own p-value, weighted by its design; effect sizes are not pooled. A combined result inside a cancer type is "detected across sources"; "specific" stays with the interaction rule inside one source. The difference between sources (`smd_difference*`) is in standardised units, which a noisier screen lowers, so a difference there is not by itself a weaker dependency.
- **CRISPR against RNAi.** `inside_expected_by_chance` in `combined_cross_assay_summary` assumes two assays that share nothing. These share 260 of 307 lines, so agreement above that figure is expected with no Y effect and is not evidence of one.
- **BeatAML2.** The exposure is a Y-loss fraction from RNA (seven Y genes, each relative to its X partner), set to 0 below 0.29. The p-value is a Freedman-Lane permutation p (100,000 permutations). `high_loss` in `beataml_drug_results` counts the patients with a fraction of at least 0.5 among those tested with that inhibitor (11 for the one hit, 15 in the cohort). `fisher_p` in `beataml_ras_by_y_loss` compares a fraction of at least 0.5 with all other patients. The two RAS-pathway rows of `beataml_hit_followup` were added after the result and do not enter the robust rule.
- **Reproducibility.** Stage 30 was rerun on the same node type with identical bytes. Stages 26 and 28 were rerun on a node with AVX-512 against stored tables from a node without it: every label, count and call is identical, and no number differs by more than 4.5e-10 of its size. Jobs and code snapshots are in `results/loy_reanalysis_v1/submission_ledger.tsv`.
- **Checks.** `results/loy_reanalysis_v1/validation_summary.json` holds 56 checks, 28 of them for these stages. The downloads are checked against the publishers' checksums, which is a local integrity check. Nothing here is biological validation: these are associations in public data.

## Addition 2026-10-05: DepMap NextGen as a third CRISPR source (stages 31 to 33)

DepMap NextGen (Neiswender et al., Nature 2026) screened 147 organoid and neurosphere models with one Cas12a library; its scores are Chronos gene effects. It is screened alone for vulnerabilities that go with Y loss and then combined with the cell lines and the Sanger organoids. The cell-line side is the frozen 26Q1 analysis; the traditional screens inside the NextGen files are not used.

Source data added (read-only; sources and checksums in `docs/DATA_GUIDE.md` and `docs/SHA256SUMS_nextgen_inputs`):

| Location | Used for |
|---|---|
| `raw/studies/depmap_nextgen/next_gen_expression.csv`, `next_gen_copy_number.csv` | RNA level and copy number of the seven Y genes; X copy number |
| `raw/studies/depmap_nextgen/supplementary_tables/41586_2026_10843_MOESM3_ESM.xlsx` | Recorded sex, model names, source laboratory, medium, treatment status |
| `raw/studies/depmap_nextgen/figshare/model_metadata.csv`, `references/depmap_nextgen/` | The figshare model table (kept for reference), the paper and the figshare archive |

The screens, QC table, mutation tables and signatures already in `raw/studies/depmap_nextgen/` are read as they were and checked against `docs/SHA256SUMS`.

| Path | Contents |
|---|---|
| `scripts/loy_reanalysis_v1/stage_nextgen_inputs.sh` | Copies and verifies the six supplied files |
| `scripts/loy_reanalysis_v1/31_nextgen_prepare.R` | Sex, Y calls, donors, cohort, matrices, gene annotation, the cancer types to pool, the chance of replication; the NextGen freeze |
| `scripts/loy_reanalysis_v1/32_nextgen_scan.R` | Controls, 300 shuffles, the scan, the four paired genes, the continuous scan, the screen by cancer type, follow-up |
| `scripts/loy_reanalysis_v1/33_combine_three.R` | Three sources together, the pooled test, shared donors, follow-up of combined hits, what each test could detect |
| `scripts/loy_reanalysis_v1/extension_helpers.R`, `tests/test_extension_helpers.R` | Helpers and their unit checks (132 in all, 40 of them for these stages) |
| `results/loy_reanalysis_v1/figures/fig13_nextgen.*` | Figure 13, with source data in `figure_source_data/` (`fig13a` to `fig13f`) |

Tables in `derived/loy_reanalysis_v1/extension/`, copied to `final/` (CSV and RDS):

| Table | One row per | Notes |
|---|---|---|
| `nextgen_y_calls` | NextGen model (314) | Recorded sex, RNA level of the seven Y genes, chrY and chrX copy number relative to ploidy, RNA call, DNA call, tier, the integrated call and five alternative definitions |
| `nextgen_cohort` | screened model (147) | Screen format and quality, lineage, stratum, donor, covariates, the flags `eligible`, `called`, `first_of_donor` and `in_primary` (65: 27 Y-negative, 38 Y-positive); `sanger_model` is the Sanger organoid of the same name for the 10 models whose donor Sanger also screened |
| `nextgen_donor_groups` | pair of models of one donor (6) or model shared with the Sanger cohort (10) | What joined them (HCMI donor ID, name root, genotype similarity), which model was kept (for a shared donor, the organoid in the Sanger cohort), and every organoid Sanger screened of that donor |
| `nextgen_gene_annotation` | gene (18,435) | Chromosome, Y partner, common-essential, missing scores, `in_scan` (18,137), `selective` (2,622), `named` |
| `nextgen_test_counts` | lineage in the cohort (7) | Models per Y group inside and outside, whether the gate of five is met, family sizes; written before any fit |
| `nextgen_pooled_types` | cancer type and source (34) | Models per Y group, whether the source can test the type alone, whether the type is pooled (pancreas) |
| `nextgen_replication_chance` | paired gene (4) | The chance of replication at the cell-line effect size, computed before any Y test |
| `nextgen_input_manifest` | input file (13) | The six supplied files and the seven already in the collection, with checksums |
| `nextgen_weights` | malignant screened model (144) | The limma weight of each model, estimated without Y status |
| `nextgen_controls` | control (2) | KRAS in KRAS-hotspot models, MDM2 by TP53 status |
| `nextgen_calibration` | run (observed and 300 shuffled) | Smallest p, FDR hits, tail counts and lambda of the scan, lambda, mean t and share of positive estimates among the common-essential genes, calls of the paired-gene rule, and the same for the test inside oesophagus/stomach |
| `nextgen_pan_results` | gene (18,137) | Effect in Chronos units, p, FDR, `passes_family_wise_bar`, `discovery`, unweighted, adjusted and quality-adjusted estimates, dependency-call contrast with and without screen quality, label |
| `nextgen_paired_genes` | paired gene (4) | Effect, group means here and in the cell lines, BH among the four, the cell-line effect and the p-value of the difference, `outcome` (the rule on the model's p-values), the permutation p from a million shuffles, `outcome_permutation` (the rule on those, which is the verdict reported), `robust`, the chance of replication |
| `nextgen_xy_partner_scan`, `nextgen_continuous` | X partner in the scan (11); gene (18,137) | The X partners as a family; the scan with Y level as a continuous score |
| `nextgen_type_screen` | gene inside oesophagus/stomach (6,324) | Interaction and the Y effect inside the type; `selective` marks the 2,622 genes of the family |
| `nextgen_hits`, `nextgen_hit_followup` | gene at FDR 0.05 (3); gene, context, check and variant (302) | `class` says whether a hit is a paired gene, a discovery or above the bar; the follow-up covers the hits, the X partners and the paired genes inside oesophagus/stomach |
| `nextgen_combined_pan` | gene in the cell-line scan and at least one other (17,159) | Each source's p and z, the pooled z, p and FDR, the number of supporting sources and the label, the two-source result beside it, and the cell-line against NextGen difference in Chronos units |
| `nextgen_combined_followup` | combined hit, check and variant (59) | The follow-up on the combined statistic: a million joint shuffles, each model left out, RNA-only and DNA-confirmed calls, adjusted models; `bh_cut` is the p-value the hit needs to stay one at its rank |
| `nextgen_combined_hits` | combined hit (6) | Rank and `bh_cut`, the permutation p and the largest leave-one-out p against that cut, the FDR with two sources, with screen quality or the genomic covariates in NextGen's model and with every source's adjusted model, group means and the share of models called dependent in the cell lines and in NextGen, the RNAi p |
| `nextgen_combined_within_type` | gene inside oesophagogastric adenocarcinoma (6,064) | The Y effect inside the type in each of the three sources and combined |
| `nextgen_pooled_type` | gene (17,097) | The pooled pancreas test: estimate in rank and Chronos units, exact p, family-wise p, the floor; a ranking, no calls |
| `nextgen_pooled_named` | paired gene and source (12) | The four paired genes pooled and in each source, with BH among the four |
| `nextgen_pooled_sensitivity` | paired gene and lineage (12) | The same test in bowel, oesophagus/stomach and brain; secondary, no calls; `estimate_chronos` is empty where the pool holds Sanger scores, which are not Chronos scores |
| `nextgen_shared_donors` | NextGen model whose donor Sanger also screened (10) | The same-named Sanger organoid, the Y call and Y evidence at each institute, every Sanger organoid of the donor, and the similarity of the two dependency profiles |
| `nextgen_xy_partners` | X partner and source (64) | Standardised effects in every source and combined |
| `nextgen_power` | test and source (7) | Models per Y group and the best-case effect found 80% of the time |

Not copied to `final/`: the matrices (`nextgen_effect.rds`, `nextgen_dependency.rds`), three stage summaries (`nextgen_*_summary.json`), the freeze with the first one archived beside it (`extension_freeze_nextgen_archived_*.json`), and the baseline taken before these stages (`nextgen_baseline_sha256.txt`, `nextgen_baseline_report.html`).

How to read them:

- **Kinds of finding.** DDX3X, EIF1AX, ZFX and RPS4X were named before the screen, so in NextGen they are replication. Any other gene is a discovery and, in NextGen alone, must clear the bar from the shuffles as well as FDR 0.05 and 0.10 Chronos units. POLR1G and OCEL1 pass FDR and not the bar; they are not called.
- **The verdict on the four.** `outcome` is the rule fixed in advance on the model's p-values, and all four pass it. `outcome_permutation` is the same rule on permutation p-values from a million shuffles, added after the first full run because the shuffles showed the model's p-values to run small (the rule fired in 18 of 300 null runs, about twice its rate). On those, DDX3X, EIF1AX and RPS4X replicate and ZFX does not (p 0.051). The report calls ZFX "not confirmed", not "not replicated": the group means show its smaller contrast to come from more dependent Y-positive models, which a Cas12a guide that also cuts the Y copy would give.
- **NextGen Y status.** Sex is the recorded sex; no file here can catch a mislabelled female among Y-negative males (`possible_female` marks the Y-negative models with a female-like X copy number, and the variant `x_checked` leaves them out). Recorded females show almost no Y signal in the expression table, so "no Y expression" is at most the female maximum or 2% of the typical male level. Y status belongs to the culture: HCM-SANG-0299-C15-B is Y-positive at Sanger and Y-negative at the Broad.
- **Strata.** The model has one term per lineage and screen format (`stratum`). Seven cohort models sit in a stratum with one Y group and carry no information.
- **Calibration.** Shuffled Y status gives a gene at FDR 0.05 in 26 of 300 runs. The observed inflation factor (1.16) is reached by 9 of the 300 shuffles and sits in the common-essential genes (1.85, 63% of estimates positive), as noisier Y-negative screens would give; with screen quality in the model it is 1.05. The shuffles do not protect against this bias, because shuffling Y status also breaks its link with screen quality.
- **Combined hits.** `evidence` counts the sources that support a combined hit. POT1 and ARHGAP26 cross FDR 0.05 in the three-source combination with effects below 0.10 Chronos units in the cell lines. POT1 is less needed with Y loss and ARHGAP26 is not a dependency, so neither is a vulnerability of Y loss; neither was a hit with two sources, and neither stays one with genome doubling, TP53 and aneuploidy in NextGen's model. The follow-up in `nextgen_combined_followup` and `nextgen_combined_hits` was designed after these hits were seen.
- **Pooled test.** Ranks are taken within a source over its models of the cancer type; Y groups are compared within strata (library in the cell lines, screen format in NextGen); the p-value counts every arrangement of Y status within strata, or samples 100,000 of them when there are more than two million. `floor` is the smallest p-value the test can give. The cell-line half of the pancreas test uses lines that were part of the scan that named the four genes; the NextGen half does not.
- **Freeze.** `extension_freeze_nextgen.json` replaced the first freeze after the review of 2026-10-05, when the pairing of one shared model with its Sanger organoid was corrected and two settings (the number of shuffles for the paired genes and for combined hits) were added. Stage 23 proves that nothing else changed: every other frozen file has its first hash, the cohort has it again once the pairing column is put back, and the parameters have theirs without the two settings.
- **Reproducibility.** Stages 32 and 33 were rerun on AMD nodes (`SBATCH_CONSTRAINT=amd`, jobs 12306255 and 12306190) against stored tables from an Intel AVX-512 node: over all 28 NextGen tables every label, count, call and missing value is identical, and no number differs by more than 1.2e-10 of its size. Two things had to be fixed for that to hold. Permutation counts now treat a shuffle within 1e-9 of the observed statistic as a tie (`tie_tolerance`), and a dependency-call contrast carries no p-value where the calls leave no Y contrast.
- **Checks.** `results/loy_reanalysis_v1/validation_summary.json` holds 81 checks, 25 of them for these stages. The six supplied files are checked against checksums taken when they were copied (figshare's own MD5 for one of them), which is a local integrity check. Nothing here is biological validation: these are associations in public data.

## Addition 2026-10-06: driver status inside cancer types (stages 34 to 36)

The context screen tests each driver alteration across all cancer types at once, with cancer type as an adjustment. These stages ask the question inside one cancer type: is there a Y effect among the lines that carry a driver alteration, and does it differ from the lines of the same type without it? The purpose is to show what the collection can and cannot test at this level, so the census and the detection limits lead and the tests support. No new source data: the stages read the 26Q1 mutation and copy-number tables, the frozen cohort and context panel, and the Sanger organoid and NextGen tables of the earlier additions.

| Path | Contents |
|---|---|
| `scripts/loy_reanalysis_v1/34_driver_type_prepare.R` | One status per gene (driver mutation, amplification or deep deletion), the census of every gene-by-cancer-type pair, the splits and what each design keeps, the same for the organoids, a census for NextGen, the pairs to combine, the gene families; the freeze |
| `scripts/loy_reanalysis_v1/35_driver_type_tests.R` | Controls, one model per split, the rule, least-squares p-values beside it, the four paired genes inside each split, 200 shuffles, organoids, the combination, follow-up |
| `scripts/loy_reanalysis_v1/36_driver_type_power.R` | Clean shifts at six bars, benchmarks moved in two ways, enlarged copies of each design and lines needed, the check that cuts real splits and copies them back, organoid limits |
| `scripts/loy_reanalysis_v1/extension_helpers.R`, `tests/test_extension_helpers.R` | Helpers and their unit checks (172 in all, 40 of them for these stages) |
| `results/loy_reanalysis_v1/figures/fig14_driver_type.*` | Figure 14, with source data in `figure_source_data/` (`fig14a` to `fig14f`) |

Tables in `derived/loy_reanalysis_v1/extension/`, copied to `final/` (CSV and RDS; the screen as `.csv.gz` and RDS):

| Table | One row per | Notes |
|---|---|---|
| `driver_type_events`, `driver_type_status` | altered line and gene (5,624); line (568, one column per gene, 838 genes) | Which kind of event makes a line altered; status 1, 0 or missing |
| `driver_type_organoid_events`, `driver_type_organoid_status` | altered organoid and gene (853); organoid (95, 261 genes) | From the study's table of driver events; amplification only for the curated genes and at twice ploidy |
| `driver_type_census`, `driver_type_census_summary` | gene and cancer type with an altered line (10,222); source and level (5) | The four cell counts and the label: specificity testable, Y effect among altered only, below the gate, one Y group missing |
| `driver_type_gate_distance` | source and level (5) | Pairs with at least 1, 2, 3, 5, 10, 20 altered lines per Y group, and with five in all four cells |
| `driver_type_splits`, `driver_type_split_lines` | split (75, of which 68 fitted); line of a fitted split (2,548) | Cells, donors per cell, information kept (`kept_inside`, and `kept_inside_alone` for altered lines fitted alone), the smallest p a two-sided exact test could give (`exact_floor`; `robust_possible` compares it with the robustness threshold and is stricter than the follow-up's own floor, see `robust_reachable` in `driver_type_power_splits`), `shared_group`, `same_lines_as`; the frozen lines with their Y call and status |
| `driver_type_organoid_splits`, `driver_type_organoid_split_lines` | organoid split (8); organoid of a split (366) | As above for the organoids |
| `driver_type_pairs`, `driver_type_genes` | split present in both sources (3); gene and source (12,776) | What is combined; the gene families (6,094 for cell lines, 6,682 for organoids) |
| `driver_type_controls` | control (4) | KRAS by KRAS status inside lung and bowel, MDM2 by TP53 status inside CNS/brain and lymphoid |
| `driver_type_screen` (`.csv.gz`, RDS) | split and selective gene (414,052) | `inside` (Y effect among altered lines), `outside` (among unaltered), `interaction` (their difference, only with four cells of five), on the rank scale (`_rank`, with the p-values) and in Chronos units; FDR over the layer and within the split; `inside_p_ols`, `interaction_p_ols` and `inside_fdr_ols`: plain least squares on the same ranks, a description and no call; `inside_p_raw`: limma on the Chronos scale, kept with the estimate and not calibrated in cells this small; `altered_set`; the same alteration in other cancer types (`elsewhere`) |
| `driver_type_split_results` | split (68) | Tests, smallest p, inflation factor, hits, the split's own bar |
| `driver_type_paired_genes` | paired gene and split (272) | The four paired genes as a named family (`named_fdr`, `shown`) |
| `driver_type_hits`, `driver_type_hit_followup` | hit (none); test, check and variant (1,611) | The follow-up covers the 161 tests where a paired gene is seen: permutation within driver status, leave one line out, three Y definitions, adjusted model, one line per donor, altered lines alone, the organoid counterpart |
| `driver_type_calibration` | run and source (402) | Observed and 200 shuffled runs for cell lines and organoids: tail shares, inflation factor, hits, and the same for the least-squares p-values (`ols_`) |
| `driver_type_organoid_screen`, `driver_type_combined` | organoid split and gene (27,136); pair and gene (6,666) | The organoid tests; cell lines and organoids combined in three colorectal splits |
| `driver_type_power_splits` | split (68) | Bars, design unit, the gap to the stored p-values by refitting (zero) and through the simulation path (`path_gap`), detection at the largest shift, the smallest p the follow-up's shuffles can give (`followup_floor`) and `robust_reachable` |
| `driver_type_shift_detection`, `driver_type_limits` | split and shift (1,564); split, bar, unit and target (912) | Share of simulated genes found at each of six bars; the shift found 50% and 80% of the time, `reachable`, `ceiling` (detection at the largest shift), and the normal-theory figure (`best_case`) |
| `driver_type_benchmark_detection` | split, benchmark and design (3,400) | Known effects moved in as a typical gene and as a copy of the gene; `own_driver` marks a benchmark moved into a split of its own driver |
| `driver_type_enlarged_detection`, `driver_type_enlarged_limits` | split, multiple and shift (3,740); split, bar and multiple (1,140) | The same on 1, 2, 4, 8 and 16 copies of the design |
| `driver_type_lines_needed` | split, bar and reference effect | `convention`: a known effect over its own gene's spread (`own spread`), the same Chronos shift over a typical gene's spread (`typical gene`), or a plain size in SD (`size`); `multiple`: the multiple of today's lines at which it is found 80% of the time (infinite: more than 16); `at_most` where it could not be read between two simulated sizes; `multiple_normal_theory` for comparison; `below_grid`, `above_grid` where the effect lies outside the simulated shifts |
| `driver_type_cut_check` | set, fraction and bar | The largest sets cut to a half or a quarter and copied back to full size: the limit of the real split, of the cut split and of the copied split, and the copied limit against the real one (`copied_against_real_log2`; below zero: copies are kinder than lines) |
| `driver_type_organoid_limits` | organoid or combined split and bar (33) | Best cases from the fitted designs |

Not copied to `final/`: three stage summaries (`driver_type_*_summary.json`), the freeze (`extension_freeze_driver_type.json`), and the baseline taken before these stages (`driver_type_baseline_sha256.txt`, `driver_type_baseline_report.html`).

How to read them:

- **Split.** A gene in a cancer type with five altered lines per Y group. 68 are fitted; 24 also have five unaltered lines per Y group and can test a difference by driver status. Lineage, disease and subtype overlap, and many splits share most of their lines (`shared_group`); hits in splits that share lines would be one finding.
- **Wording.** A call means that a gene's Y effect differs by driver status inside a cancer type and is present among the altered lines. It never means "only in": the same alteration in other cancer types is a description (`elsewhere`), not part of the rule.
- **Bars.** `candidate`: one gene named in advance, p below 0.05. `split_alone`: BH within one split; kept for the power statement only, since 191 of 200 shuffled runs have such a test somewhere. `family`: the layer's rule among altered lines (p below 1.2e-7 for one more discovery). `specificity`: the rule that also needs the difference. `named_difference`: the difference for one gene named in advance, p below 0.05. `family_raw`: the `family` bar met by plain least squares on Chronos scores, for comparison.
- **Limits.** Properties of the designs and of the layer's test, not observed power. A bar with `reachable` false is not reached 80% of the time at any shift up to 10 SD; `ceiling` says how often a 10 SD shift is found. That ceiling belongs to the rank-based test, because ranks cap how far a few lines can move; least squares on Chronos scores (`family_raw`) has none, and is not calibrated in cells this small. Normal theory (`best_case`) agrees for a named gene and is too kind at the genome-wide bars.
- **The test borrows its variance.** limma's moderated t takes a gene's variance almost wholly from the other genes here, so it is conservative for a quiet gene and generous for a variable one. `inside_p_ols` lets each gene keep its own: by it, 11 cell-line tests and 2 organoid tests pass FDR 0.05, all DDX3X or EIF1AX, and no difference by driver status does. These p-values were added after the review and make no call.
- **Lines needed.** For the layer's test, from simulations on copies of each split's design, each copy with the spread of a different real gene and ranked over its own 568 lines; the bar and the variance prior are today's. Multiples are interpolated between 1, 2, 4, 8 and 16. How large a known effect counts as depends on the gene it lands on: DDX3X's 0.75 Chronos units are 1.9 SD of its own spread and 4.2 SD of a typical selective gene's, and the table carries both. The check on real lines (a split cut down and copied back against the real split) shows that copies are a little kinder than lines, so the multiples are somewhat low.
- **The same test twice.** Head and neck at lineage level and head-and-neck squamous cell carcinoma at disease level hold the same altered lines for TP53 and give the same p-values, so 6,089 of the 414,052 tests are repeats; the family counts them.
- **Added after the review.** Four things were added after an independent review of the first full run and are marked in the report: the least-squares p-values, the bar `named_difference`, the bar `family_raw`, and the second way of sizing a known effect. The check on real lines was replaced, and `robust_reachable` corrects `robust_possible`.
- **Paired genes.** DDX3X, EIF1AX, ZFX and RPS4X are a family named in advance inside the splits and serve as the positive control: 161 of 272 tests seen, 40% where the smaller altered cell has five lines.
- **Freeze and trial.** No gene was compared with Y status inside a split before `extension_freeze_driver_type.json` existed, with one exception recorded in the progress ledger: a trial run printed three real estimates for RPS4X in one split, because three follow-up variants read real Y calls. Trials now shuffle those as well.
- **Reproducibility.** Stages 34 to 36 were rerun on AMD nodes (`SBATCH_CONSTRAINT=amd`, jobs 12315272, 12316474 and 12315290) against stored tables from Intel AVX-512 nodes: over all 31 `driver_type_` tables every label, count, call, flag and missing value is identical and the rows come in the same order, and no number differs by more than 1.5e-9 of its size. The one exception is `path_gap` in `driver_type_power_splits`, a difference between two nearly equal numbers, where two values differ by 7.5e-12. One thing had to be fixed for that to hold. Head and neck at lineage level and head-and-neck squamous cell carcinoma at disease level give the same test among TP53-altered lines, their p-values tie to the last bits, and which row came first depended on the processor. The follow-up takes the seed of its shuffles from the row, so one permutation p-value differed by one shuffle in 20,000. Tied tests now keep the order of the split table (`tied_order`, with a unit check). In the stored tables this moved four rows of `driver_type_paired_genes` and 40 of `driver_type_hit_followup` and changed that one p-value from 2 to 1 in 20,001 (stage 35 rerun as job 12316457); no flag, count or reported number changed. One thing was left as it is: `driver_type_splits.csv` and `driver_type_organoid_splits.csv` differ between processor types in the last digit of the columns that give the information a design keeps, so the hashes in the freeze are those of the Intel run. A rerun of stage 34 on another processor type belongs in a scratch folder (`LOY_EXTENSION_OUTPUT`, see `scripts/loy_reanalysis_v1/README.md`); in place it would rewrite those two tables and then stop at the freeze.
- **Checks.** `results/loy_reanalysis_v1/validation_summary.json` holds 115 checks, 34 of them for these stages. Nothing here is biological validation: these are associations in public data, and the limits are simulations on it.
