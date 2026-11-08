## Expanded all-collection GSEA

[All-collection results and figures](outputs/programme_analysis_v1/all_collections/report.html): 23,696 eligible sets; 3,281 overall global-FDR discoveries. The earlier three-collection analysis is retained.

## Current programme workflow (18 September 2026)

[Agreed analysis plan](PLAN_20260918.md) | [Applied analysis report](outputs/programme_analysis_v1/report.html) | [Execution status and reproduction](outputs/programme_analysis_v1/README.md)

# SEPT_11_2026_REANALYZE

Core analysis of HT1080 LOY versus sgCTRL in August 2025 and May 2026. The same established populations were sampled again: one population per condition, one shared nucleofection origin, six within-period passage/day matches. Primary outputs are descriptive expression differences and temporal repeatability.

Read the [current synthesis of the HT1080 LOY transcriptional phenotype](HT1080_LOY_TRANSCRIPTION_SYNTHESIS.md) for the consolidated biological interpretation, supporting evidence, HDAC external-benchmark update, PRC2/H3K27me3 findings and remaining mechanistic questions (14 September 2026).

The subsequent [protein-association network report](outputs/report/network_modules.html) builds fresh STRING/MCODE modules from the current descriptive candidates. It includes [lower-RNA](outputs/figures/network_modules/HT1080_LOY_network_modules_DOWN.pdf) and [higher-RNA](outputs/figures/network_modules/HT1080_LOY_network_modules_UP.pdf) figures, a [complete six-page atlas](outputs/figures/network_modules/HT1080_LOY_network_module_atlas.pdf), physical-interaction sensitivity, annotations and [captions](outputs/captions/network_modules_captions.md). These additional exports have their own [manifest](outputs/figures/network_modules/manifest.csv). Run `network_modules_analysis.py` and then `plot_network_modules.py` from this analysis directory to reproduce the extension; the original expression/GSEA models are reused.

The [all-module significance plot](outputs/figures/network_modules/HT1080_LOY_all_module_annotation_significance.pdf) compares the existing raw annotation P values and BH q values for all 24 lower-RNA and 18 higher-RNA modules. [Source values and tested terms](outputs/tables/network_modules/annotation_significance_plot_data.csv) and a [caption](outputs/captions/network_module_significance.md) accompany the PDF/SVG/600 dpi PNG exports. Run `plot_network_module_significance.py` to reproduce this plot without rerunning tests. These are functional-annotation statistics; MCODE itself supplies no P value.

The [interactive GSEA pathway network](outputs/report/gsea_pathway_network.html) connects all 808 combined global-BH < 0.05 sets by shared leading-edge genes. The [leading-edge figure](outputs/figures/gsea_network/HT1080_GSEA_leading_edge_network.pdf) and [full-membership comparison](outputs/figures/gsea_network/HT1080_GSEA_full_membership_network.pdf) retain all 32 higher-RNA and 776 lower-RNA sets, including isolates, with matching positions and PDF/SVG/600 dpi PNG exports. [Captions](outputs/captions/gsea_pathway_network.md), complete [node statistics](outputs/tables/gsea_network/network_nodes.csv), [edge similarities](outputs/tables/gsea_network/network_edges.csv), community gene drivers, threshold sensitivity and GraphML networks accompany the figures. Run `build_gsea_pathway_network.py` followed by `plot_gsea_pathway_network.py` to reproduce this extension from frozen inputs. The existing global BH values are preserved; no new enrichment tests are run. Period support reflects repeated sampling of established populations.

The [nicotinic receptor follow-up](NICOTINIC_RECEPTOR_FOLLOWUP.md) examines all receptor subunits, the five UP_01 members and cholinergic context genes using saved counts and normalization. Its [publication figure](outputs/figures/nicotinic_receptors/HT1080_nicotinic_receptor_followup.pdf) shows family-wide effects, all six matched differences and individual sample abundance, with [caption](outputs/captions/nicotinic_receptor_followup.md) and source tables. It highlights stronger read support for CHRNA7/CHRNB2, low-count sensitivity of CHRNA6/CHRNB4, mixed family directions, CHRFAM7A ambiguity and the CHRNB2 link to both Mikkelsen source H3K27me3 leading edges. Run `extract_nicotinic_followup.R` then `nicotinic_followup.py`; original models and enrichment statistics remain unchanged.

Open [the report](outputs/report/analysis.html). Figures are in [outputs/figures](outputs/figures), with editable PDF/SVG and 600 dpi PNG exports. [Captions](outputs/captions/figure_captions.md) also have a [DOCX copy](outputs/captions/figure_captions.docx). Every panel has source data linked through the [figure manifest](outputs/figures/figure_manifest.csv).

The [focused review of the original notebooks](provenance/legacy_code_review.md) records the design, chromosome-arm, compensation, annotation and VIPER issues addressed or avoided in this reanalysis.

The expanded enrichment run tests every eligible human MSigDB H and C1-C9 set and subcollection (2026.1.Hs), with complete coverage/exclusion records. [fgsea results](outputs/tables/fgsea_all_collections.csv) use the full effect rankings; [CAMERA results](outputs/tables/camera_all_collections.csv) use voom and a fixed passage-blocked model with separate period effects. Global BH spans all eligible sets within each method and contrast. CAMERA assumes independent model residuals and remains exploratory for these serially sampled populations; estimated inter-gene correlation is a separate sensitivity. Figure 5 compares methods and documents voom precision estimation.

The user-requested [ZMYND8 follow-up](provenance/ZMYND8_interpretation.md) examines all eight regulon targets, period consistency, RNA abundance and target-omission sensitivity. Figure 6 shows why the positive inferred score is strongly dependent on MMP1/MMP3 and does not establish increased ZMYND8 expression or compensation for KDM5D loss.

The [chromatin follow-up](provenance/chromatin_interpretation.md) examines 37 source-defined signatures (34 eligible for enrichment), paired UP/DN response arms, 34 regulator RNAs and five gene-exclusion sensitivities. Figure 7 and Supplementary Figure 4 distinguish persistent HDAC/Polycomb-associated target patterns from measured chromatin activity. The original global BH results are retained; post hoc exclusion scenarios report nominal P values only. Source descriptions are frozen by `prepare_chromatin_resources.R` from the same MSigDB release, with the selections and reference masks recorded in `config/chromatin_followup.yml`.

The [HDAC1/2 plan](HDAC12_PLAN.md) has completed [step 1: reference audit](outputs/report/HDAC12_reference_audit.html), [step 2: gene drivers](outputs/report/HDAC12_gene_drivers.html) and [step 3: direct external benchmark](outputs/report/HDAC12_external_comparison.html). The audit checks 18 source definitions, mappings and study-design limits. Step 2 partitions eight Senese response arms, examines every represented gene, freezes complete sample scores and records gene influence, abundance, temporal consistency and process annotations. Figure 8 and Supplementary Figure 5 show a persistent HDAC2 source-DOWN pattern with little source-UP response.

The GSE193885 HAP1 benchmark reconstructs 27 contrasts from 63 libraries, compares complete effects on 14,881 common loci and tests 54 externally selected UP/DN sets using fgsea and CAMERA. Figure 9 and Supplementary Figure 6 show that active-transgene controls resemble LOY at least as well as HDAC2 knockout; the inactive-versus-active HDAC2 ranking is weakly inverse. This weakens a specific catalytic-loss interpretation. The source has unresolved MS-275 dose/time labels and unverified independent clone replication. Optional LINCS, complex analysis and addback remain unrun.

To recompute the external benchmark from its frozen resources and original HT1080 inputs, regenerate its two figures, run scientific checks and render the reports:

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla 'analysis/HT1080/SEPT_11_2026_REANALYZE/run_hdac12_step3.R'
```

The source preparation scripts are `prepare_hdac12_external.py --fetch` and `audit_hdac12_external.py`. Raw downloads are cached under `resources/hdac12_external/GSE193885/`; the XLSB reader is pyxlsb 1.0.10 installed only in `tools/hdac12_python`. The step-3 runner verifies input/resource hashes and refits external models; it does not rerun the original global MSigDB analysis. `--source-only` and `--analysis-only` deliberately stop before later stages. Re-rendering saved results is available through `render_report.R`.

To recompute only the HDAC step-2 follow-up from its hash-verified saved inputs, regenerate its figures, run the scientific checks and render both reports:

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla 'analysis/HT1080/SEPT_11_2026_REANALYZE/run_hdac12_step2.R'
```

The step-2 configuration and input hashes were recorded before its first computation. This controlled runner leaves the prior global enrichment objects and original figures intact. The complete runner also includes the new analysis and figures. `audit_hdac12_sources.py` rebuilds the reference audit tables offline from cached source documents; its `--fetch` option is reserved for explicit source retrieval.

From the RNAseq repository root, use the verified R installation:

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla 'analysis/HT1080/SEPT_11_2026_REANALYZE/run.R'
```

The runner checks frozen resource hashes and recomputes every stage. It cannot silently reuse stale `.done` markers. It reads original count files without modifying them. `renv.lock` and `outputs/tables/package_versions.csv` record the environment. `targets` was proposed but is not installed; an explicit full rebuild is the implemented orchestration. No package upgrade is required.

To re-render only the report, without recalculating or overwriting analysis objects:

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla 'analysis/HT1080/SEPT_11_2026_REANALYZE/render_report.R'
```

`prepare_resources.R` freezes the MSigDB sets, CollecTRI network, PROGENy model and UCSC cytobands on initial preparation. The included resources are sufficient for rerunning; no live annotation/network request occurs in `run.R`. Intentionally changing resources requires updating their manifest and documenting the change. The original CollecTRI retrieval date is unknown; its exact file hash is recorded.

The main effect table is `outputs/tables/HT1080_LOY_vs_sgCTRL_primary_effects.csv`. The full retained-feature table also includes low-expression genes, Y, PAR, mitochondrial and other features, each with an eligibility flag. `candidate_evidence_dimensions.csv` keeps expression magnitude, consistency, influence, pathway and TF information separate; unavailable DNA/addback/ATAC evidence is missing rather than fabricated.

DNA-RNA coordinate integration is deferred because WGS D0 specimen correspondence and genome compatibility are unresolved. Addback and ATAC analyses remain extensions; they are not represented as completed core results. No primary gene-level significance, causal LOY, protein compensation or independently replicated validation claim is supported by this design.

See [PLAN.md](PLAN.md) for the agreed scientific design, [provenance/DECISIONS.md](provenance/DECISIONS.md) for implementation choices, and [provenance/scientific_checks.txt](provenance/scientific_checks.txt) for checks of sample pairing, gene identity, arm annotation, signature selection and paralog arithmetic.

## PRC2 and epigenetic deep dive

The September 14 [HT1080 chromatin resource review and integration plan](HT1080_CHROMATIN_RESOURCE_REVIEW_AND_PLAN.md) audits the earlier `A:/multiome/results/ChIPseq` work and public datasets. Its [resource catalogue](provenance/chipseq_resources/chromatin_resource_catalogue.csv) prioritizes the measured GSE261834 histone panel, independent ChIP references and group-level ATAC integration. ENCODE imputations are explicitly distinguished from measurements. The user confirmed that the six ATAC libraries are successive passages and their direct RNA partners are uncertain; older independent-replicate ATAC inference must be revised for this integration.

Step 1 is complete: the [measured HT1080 chromatin atlas](CHROMATIN_REFERENCE_ATLAS.md) contains eight histone marks and two lamin peak tracks, all 25,461 ranked genes, 20 gene-set definitions with GMT exports, exact signal/coverage summaries and TSS/window/threshold sensitivity. The [illustrated report](outputs/report/chromatin_reference_atlas.html) includes a PDF/SVG/600-dpi PNG figure and caption. Of the ranked promoters, 467 remain unknown because of blacklist masking; source replicate pooling and coarse lamin-signal normalization remain qualified. Independent interval and basewise signal checks passed. Run `python run_chromatin_atlas.py` from this directory to rebuild from frozen inputs. Atlas construction is separate from enrichment and ATAC integration, and is not added to the original full-analysis runner.

Step 2 is complete: [measured chromatin/RNA enrichment](CHROMATIN_RNA_ENRICHMENT.md) tests all 20 custom sets with GSEA and both CAMERA correlation settings, retaining large sets in full and a separate custom BH family. Five primary GSEA sets enrich at the lower-RNA end; neither CAMERA setting supports a discovery. H3K4me3/H3K27me3 co-marked promoters remain unenriched, while broad H3K27me3 changes direction between August and May. Abundance/length/arm sensitivity retains three of the five GSEA associations. The [standalone report](outputs/report/chromatin_rna_enrichment.html) contains four publication figures, [captions](outputs/captions/chromatin_rna.md), full gene evidence and numerical limitations. Two large-set estimates remain missing in one passage sensitivity after the prescribed retry. Original global MSigDB statistics are unchanged. Run `python run_chromatin_rna.py` to reuse the frozen inputs and completed partitions, validate the results and regenerate this extension; `--statistics-only` omits figure/report regeneration. This step measures associations with reference chromatin, not LOY-associated histone changes, and performs no ATAC integration.

The [focused report](outputs/report/epigenetic_deep_dive.html) extends interpretation to 339 epigenetic catalogue sets, complete PRC2 response arms, gene overlap/influence, normalization and passage sensitivity, and 41 source promoter-state definitions. Figure 10 and Figure 11 add eight panels, bringing the project to 17 figures and 62 panels. EED-depletion-induced genes mostly decrease in LOY, while two source H3K27me3 promoter sets increase; their small gene overlap argues against merging these results into one PRC2 activity score. The report distinguishes source histone marks, source RNA responses and HT1080 regulator RNA.

Run only this extension from frozen inputs:

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla 'analysis/HT1080/SEPT_11_2026_REANALYZE/run_epigenetic_deep_dive.R'
```

`--analysis-only`, `--figures-only --skip-report` and `--finalize-only` support explicit stages. The complete runner also includes this extension. `prepare_epigenetic_resources.R` extracts same-release cached metadata; `audit_epigenetic_sources.py` caches primary bibliographic records. Original global BH values are reused without new tests. Frozen input and prior figure hashes are checked by `tests/check_epigenetic_deliverables.py`; new scientific contracts run during finalization. No new public PRC2 expression dataset, chromatin measurement or independently generated HT1080 population was analyzed. See [the plan](EPIGENETIC_DEEP_DIVE_PLAN.md).

## Chromatin GSEA publication plots

The [GSEA figure report](outputs/report/chromatin_GSEA_curves.html) contains two six-panel running-enrichment figures, Figure 12 and Figure 13, with individual PDF/SVG/600 dpi PNG downloads and captions. The panels cover PRC2 depletion responses, H3K27me3 source promoter sets, HDAC1/2, BMI1/MEL18, SMARCE1 and Suz12 deficiency. Their full curves and leading edges reproduce the original fgsea results exactly. Original global BH values remain unchanged; the EED DN and weaker SUZ12 UP results are retained as comparisons.

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla 'analysis/HT1080/SEPT_11_2026_REANALYZE/run_chromatin_GSEA_plots.R'
```

This runner reuses completed models and writes only plot data, figures, captions and reports. `--data-only` and `--skip-report` allow explicit stages. The individual export manifest is under `outputs/figures/chromatin_GSEA_individual/`; the composite figures use the main figure manifest. The project now has 19 composite figures and 74 source panels, plus 12 separately exported GSEA panels. Validate with `tests/check_chromatin_GSEA_exports.py` and the existing `tests/check_deliverables.py`. After visual review, `finalize_chromatin_GSEA_plots.R` records review status and renders reports. Historical completion/output-hash manifests describe the state at their respective task completions; original-model and prior-figure preservation snapshots remain checked separately.

## Batch-aware gene-level extension (18 September 2026)

The author requested a gene-level batch-aware analysis in this report. The new [report section](outputs/report/analysis.html#batch-aware-gene-level-differential-expression) and [standalone report](outputs/report/batch_aware_de.html) add limma-voom models without replacing the descriptive core summaries. Original raw counts, 12 core libraries, 25,526 eligible genes and core TMM offsets are reused; the batch-adjusted PCA matrix and parentals are excluded.

The main extension model uses six fixed within-period passage-pair intercepts and two period-specific LOY effects, with an equal-weight combined contrast (four residual df). A batch + condition model is included as a sensitivity (nine residual df). Both use robust empirical Bayes moderation and BH across genes separately per contrast. All P values and confidence intervals are conditional on these fixed models and independent residuals, not independent generation replication. This supersedes the earlier absence of any gene-level significance output only for this explicitly exploratory extension; primary descriptive tables and existing pathway/network outputs remain unchanged.

Reproduce only the extension from this directory:

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla run_batch_aware_de.R
python plot_batch_aware_de.py
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla render_batch_aware_de.R
python tests/check_batch_aware_de.py
```

The complete `run.R` also fits the extension after core expression analysis. Extension outputs use the `batchDE_` table prefix and `batch_aware_DE.rds` model. Its figure has a separate `outputs/figures/batchDE_figure_manifest.csv`, and the thesis-style master exports are under the thesis Chapter 3 figure folder. Existing global figure numbering is unchanged. Input hashes and checks are in `provenance/batchDE_input_hashes.csv` and `batchDE_validation.txt`.

The STRING/MCODE modules still use persistent descriptive candidates, not this new FDR list. No module, pathway or transcription-factor analysis was rerun.
