# HT1080 RNA-seq: robust TF regulators and shared-target TF modules (thesis 3.4.7.1)

## Context
The thesis Results currently report 60 "inferred regulators" and say "No regulator figure is currently selected". These come from one `decoupleR::run_ulm` call on the limma **logFC**, using one frozen CollecTRI file of unknown date (`run_programme_analysis_v1.R` L126–147).

That approach has three weaknesses:
- It is exactly what the project's own critique (`github/.../warnings_about_TF_regulons.md`) warns against: one method, logFC input, overlapping regulons left unhandled.
- ULM p-values ignore correlation between genes. CAMERA with estimated correlation found no significant gene sets, so BH on these p-values is optimistic.
- AP-1-family hits (FOS/JUN/JUNB/AP1/SRF) are open to self-edge and immediate-early/handling artefacts.

**Author decisions:**
- Output: thesis main figure + supplement + Results paragraph.
- "TF modules" = robust regulators grouped by shared targets.
- RNA only, no ATAC.
- Compute: local R 4.5.2 plus a conda env or WSL for Python.

**Intended outcome:** a small, defensible set of regulators and 2–5 shared-target modules, each with stability, a null model and confounder checks, and a figure that replaces the unused regulator panel.

**Tool choice** (from a 2024–26 survey):
- **Methods:** decoupler (ULM/MLM/wsum/VIPER), keeping CollecTRI as the primary prior.
- **Alternative priors:** DoRothEA A–C and TFActProfiler (reported *NAR* 2026, ChIP-Atlas-informed signed prior), because the benchmarks show prior choice drives results.
- **Rotation test:** limma `mroast`, which allows for correlation between genes.
- **Rejected:** single-cell GRN tools (SCENIC+, CellOracle, LINGER), data-driven GRN learning (GENIE3/WGCNA/NMF; n=12 is too small) and LLM priors (about 50% accuracy).

## Location and conventions
- **Analysis root:** `RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/`. New code goes in `tf_modules_v1/`, outputs in `outputs/tf_modules_v1/{tables,models,figures,provenance}`. Original outputs are read-only and hash-checked.
- **Planning records:** first write `PLAN_20260924_tf_modules.md` and a dated `provenance/DECISIONS.md` entry that **fixes the robustness and module rules before any run**, because the current 60 hits are already known.
- **R:** run `"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" --vanilla`. The Rscript on PATH is anaconda R 4.1.
- **Installs needed:**
  - Python env `envs/tfmod_py.yml` with decoupler ≥2 and TFActProfiler (Linux deps → WSL Ubuntu-24.04).
  - Nothing new in R. decoupleR 2.16, OmnipathR 3.18.4, dorothea 1.22, viper 1.44 and limma 3.66 are already installed.
- **Reused inputs:**
  - `outputs/programme_analysis_v1/models/DE.rds`: voom EList, fit, design, contrasts, `rankings` (moderated t).
  - `outputs/models/batch_aware_DE.rds`.
  - `outputs/models/02_expression.rds`: logCPM, 6 matched differences, leave-one-out effects.
  - `outputs/programme_analysis_v1/tables/`: `TF_activity_all.csv`, `directionally_consistent_genes.csv` (662), `TF_programme_links.csv`.
  - STRING `modules.csv` and programme communities.
  - The leave-one-out and target-omission pattern in `R/04_zmynd8.R`.

## Stages
**0. Priors** (`prepare_tf_priors_v1.R`, `tfactprofiler_prior.py`). Everything saved to `resources/tf_modules_v1/`:
- **CollecTRI:** a current dated pull (`get_collectri(split_complexes=FALSE)`), recording retrieval date, version, sha256 and an edge diff against the frozen file.
- **Frozen CollecTRI:** kept for the bridge check.
- **DoRothEA A–C.**
- **TFActProfiler:** first verify the repo, licence, size and edge-list availability. If it only outputs its own scores, use them for direction agreement only. If it can't be installed, mark it "not assessable" and don't block.
- **Edge cleanup:** conflicting duplicate edges removed as in v1. Self-edges are flagged, not dropped.

**1. TF scores** (`run_tf_activity_v1.R`):
- **Input:** moderated t for overall, AUG_2025, MAY_2026, interaction and CTRL_vs_parental. logFC is secondary.
- **Methods:** ULM, MLM, wsum (normalised, 1,000 permutations, seed 42) and VIPER (analytic), minsize 5, across 3 priors.
- **Bridge:** frozen prior + logFC + ULM must reproduce `TF_activity_all.csv` (difference below 1e-10).
- **Significance:** ULM(t) with BH across TFs per contrast, the same family as before.
- **Consensus:** mean signed z across methods, used for ranking and display only; p-values are never combined across methods.
- **mroast:** on each regulon (gene weights = mor, 9,999 rotations), reported as a flag that allows for correlation between genes.
- **Leave-one-passage-pair-out:** 6 limma refits.

**2. Robust regulator rule** (fixed in DECISIONS.md before running):
- R1: overall ULM(t) FDR <0.05
- R2: same sign in Aug and May
- R3: at least 3 of 4 methods agree in sign with p <0.05
- R4: same sign under every alternative prior that covers the TF
- R5: sign kept in 6/6 leave-one-pair-out refits
- R6: holds after dropping its top 3 contributing targets and its self-edge

Tiers: **Robust** meets all six. **Supported** meets R1–R2 plus at least 3 of R3–R6. State that these filters reuse the same reads.

**3. Overlap and redundancy** (`run_tf_overlap_v1.R`):
- **Target contributions:** c_g = mor·t·sign(score). Leading contributors are the targets supplying 80% of the score.
- **Concentration flag:** top 3 targets supply more than 50% (the ZMYND8/MMP pattern).
- **Unique vs shared targets:** for overlapping pairs such as AP1/FOS/JUN/JUNB/SRF, score targets unique to one TF and shared targets separately.
- **Conditional MLM:** each TF fitted jointly with its most-overlapping neighbours.
- **Flag:** "not separable from X" when a TF's unique-target score fails.
- **Copy-number flag:** re-score without genes on chromosome 18/13 or lpWGS-changed segments. Flag a regulon when half or more of its leading contributors sit on one arm.

**4. Modules** (`run_tf_modules_v1.R`):
- **Nodes:** robust regulators; add supported ones if there are fewer than 8.
- **Similarity:** weighted overlap coefficient on leading-contributor targets. Sensitivities: Jaccard, and overlap restricted to the 662 consistent genes.
- **Primary clustering:** average-linkage hierarchical, k chosen by the highest mean silhouette. Leiden across resolutions as a sensitivity (report ARI).
- **Stability:** 1,000 gene bootstraps plus 200 runs dropping 20% of edges. A module is kept at co-clustering Jaccard ≥0.6 and called "stable" at ≥0.75.
- **Nulls** (ULM(t) only, to keep runtime down):
  - 200 degree-preserving rewired CollecTRI networks, giving a null count of robust regulators and null silhouette.
  - Size-matched random TF sets for the overlap level.
  - Size- and expression-matched gene-set permutation plus `mroast` on each module's signed target union.
- **Supporting evidence (not used to define modules):**
  - Co-activity: per-library ULM on logCPM with batch and condition removed, using non-shared targets and a label-permutation test. Weak evidence at n=12.
  - CollecTRI TF→TF edges, shown with each TF's own transcript t.

**5. Annotation and confounders** (`annotate_tf_modules_v1.R`):
- **Driver targets:** leading contributors in at least half of a module's TFs.
- **Links to existing results:** hypergeometric tests against STRING modules, GSEA communities and the 662 genes, with the regulon-covered background and BH within each module.
- **Confounding programmes:** IFN-α/γ, E2F/G2M, TNFα/NF-κB, P53, UPR, hypoxia, and the dissociation/immediate-early genes of van den Brink 2017. For each module, report the overlap and re-score with each programme removed. Label a module "confounder-dependent" if it collapses.
- **AP-1/SRF check:** test whether IEG shifts are uniform across the 6 matched pairs and both periods. Check the lab records for harvest or handling order.
- **Drift flag:** a robust TF that also shifts in CTRL_vs_parental.

**6. Thesis outputs:**
- **Main figure** `ch03_main_ht1080_tf_modules`, 183 mm wide, script `figures/chapter3/code/make_ht1080_tf_modules.py` using `thesis_style.apply_style`, with source data in `source_data/ht1080_tf_modules/` plus `provenance.json` and `_notes.md`:
  - **a:** regulator robustness matrix (rows grouped by module; Aug/May/overall t-scores in blue #3C5488 / vermilion #E64B35, method/prior/leave-one-out/omission checks).
  - **b:** module similarity heatmap with stability values.
  - **c:** module × programme dot plot, with confounder columns set apart.
  - **d:** driver-target strips for 2–3 modules.
- **Supplement** `ch03_supp_ht1080_tf_module_robustness`: old-to-new bridge, prior concordance, leave-one-out traces, omission/self-edge/copy-number checks, nulls, co-activity, drift and confounder re-scoring.
- **Captions** go in `Chapter3_figure_captions.md`. Mark `ch03_main_ht1080_rna_programmes_regulators` panel b as superseded (archive first). Add to the selected index only after the author picks.
- **Results:** draft about 150 words to replace the "No regulator figure…" sentence in `Chapter3_results.md` (~L215). Use "target-pattern" wording and make only targeted edits (per the thesis editing feedback).

## Verification
- The bridge reproduces v1 exactly, and the original input hashes are unchanged.
- R decoupleR 2.16 and Python decoupler v2 agree (r >0.999) for ULM/MLM/wsum on the overall t.
- **Planted signal:** shift a random TF's targets in t; the pipeline must recover it and place it in the right module.
- Rewired networks give close to 0 robust regulators.
- Every module member satisfies its tier, and the figure `source_data` equals the analysis tables.
- Visual QA at thesis page size with the thesis-figure-style skill.
- `sessionInfo` and the conda env export are saved to provenance.

## Limitations to state
- One population per condition, so inference is conditional on the model.
- Modules are shared-target footprints; overlapping regulons can't be attributed to one TF, and well-studied TFs have larger regulons.
- A new CollecTRI pull may change the hits (documented in the bridge).
- AP-1/IEG and copy-number or IFN programmes are possible confounders.
- Co-activity evidence is weak at n=12.
- Scores reflect coordinated target expression, not TF abundance, binding or activity. ATAC support is out of scope here.
