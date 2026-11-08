# SEPT_11_2026_REANALYZE

Analysis plan, updated 11 September 2026. The confirmed core has been implemented and run. The [analysis report](outputs/report/analysis.html) contains count/annotation audits, descriptive expression effects, sensitivities, pathways/regulators, temporal signature assessment and Y/paralog dosage. Four main figures and three supplementary figures have captions and PDF/SVG/PNG exports. See [README.md](README.md) for execution and [implementation decisions](provenance/DECISIONS.md) for changes from the proposed tools. DNA integration, addback and ATAC remain explicitly deferred or optional.

**Confirmed scope and design:** HT1080, both August 2025 and May 2026, comparing the same established LOY and sgCTRL populations sampled again. Suffixes 1/2/3 identify successive passages, matched by passage/day between LOY and sgCTRL within each sampling period. The populations trace to the same original nucleofected pool; the batches are not independently generated LOY experiments. The core contains six RNA libraries per condition across two periods, representing repeated measurements of one established population per condition. Parental samples provide context; addback rescue and ATAC integration remain proposed extensions.

The aim is to characterize expression differences between the established LOY and sgCTRL populations, assess their stability over matched passages and sampling periods, examine DNA dosage and regulatory context, and prioritize follow-up mechanisms. Within this dataset, LOY status cannot be separated from other differences specific to the two established populations. Methods come from `../reannotate_v105/`, with the code-review corrections incorporated. Existing analyses provide reusable code and comparison results; their conclusions are not acceptance criteria for this run.

The user-supplied [RNA-seq methods review](A:/Dropbox/LOY_paralogs/projects/proposals_summaries/deep_research/RNAseq_methods_review.md) was read in full and incorporated on 11 September 2026. It informs quantification checks, pathway hypotheses, TF interpretation and optional network work. TMM normalization, expression summaries, curated pathways and regulons remain useful. The primary analysis is now descriptive at the established-population level; an ordinary passage-blocked limma model will not supply claims of independently replicated LOY effects. Formal DE, TREAT and model-based pathway testing are conditional extensions rather than mandatory primary outputs.

**Figure requirement:** all plots will use a consistent Nature-style, publication-quality design, with a complete caption for every figure. [FIGURES.md](FIGURES.md) defines the layout, fonts, colours, exports and visual review. [FIGURE_CAPTIONS.md](FIGURE_CAPTIONS.md) provides draft panel-by-panel captions to populate from verified metadata and fitted results.

## 1. Starting material and confirmed design

Paths below are relative to the RNAseq repository unless otherwise specified.

| Input | Proposed use | Verified during planning |
|---|---|---|
| `analysis/HT1080_BATCHES/count_matrix/HT1080_BATCHES_counts_raw.csv` | Starting raw matrix | 78,899 unique gene IDs x 15 samples; all gene IDs and counts exactly match all 15 source featureCounts tables |
| `analysis/HT1080_BATCHES/count_matrix/HT1080_BATCHES_sample_metadata.csv` | Starting sample manifest | Sample keys exactly match matrix column order; CTRL/LOY passage matching and repeated sampling of the same populations within batch are confirmed |
| `counts/STAR_v2_PARfix/AUG_2025_HT1080_LOY/` | Original counts and assignment summaries | 3 CTRL, 3 LOY, 3 PARENTAL |
| `counts/STAR_v2_PARfix/MAY_2026_HT1080_ADDBACK/` | May primary samples and optional addbacks | 3 sgCTRL, 3 LOY, and 12 additional plasmid-arm samples |
| `annotation/gencode.v49.primary_assembly.annotation.gtf.gz` | Canonical gene annotation | Present; inspected count headers identify GENCODE v49 |
| `analysis/HT1080_ADDBACK/SEPT_2026_REANALYSIS/` | Reusable annotation, rescue and integration helpers | Existing pipeline located; reuse follows a focused code review |
| `A:/Dropbox/LOY_paralogs/projects/WGS/HT1080_QDNAseq_cn_df.rds` | DNA copy-number context | Present; biological sample matching and coordinate provenance must be established |
| `A:/Dropbox/LOY_paralogs/projects/ATACseq/SEPT_08_2026_REPROCESS/annotation/promoter_gene_summary.csv` | Optional promoter accessibility integration | Present; sample, contrast and genome-build compatibility must be established |

Core scope is **12 LOY/CTRL libraries**, plus **3 parental libraries** for QC and secondary comparisons. Adding all May plasmid arms gives **27 unique libraries**, not 33: May LOY and sgCTRL already occur in the core matrix. Raw counts are referenced, not overwritten. SHA256 hashes of the three starting files and the supplied methods review are recorded in `input_inventory.csv`.

**User clarification, 11 September 2026:** suffixes 1/2/3 are passage numbers; CTRL-1 and LOY-1 were passaged first and on the same day. Corresponding passage numbers will be represented as passage/day matches within each batch. These matches do not cross batches, and a passage day does not automatically identify the RNA harvest date or sequencing batch.

**Further user clarification:** passages 1-3 come from the same starting populations. Each batch therefore contains one observed CTRL population series and one observed LOY population series, each sampled at three passages. Do not report n = 3 independent populations per condition within a batch. These serial measurements are not technical sequencing replicates and should not be automatically pooled as if they were repeated lanes of one library.

**Final user clarification:** August and May sampled the same established populations again, originating from the same nucleofected pool. Both batches are explicitly included in the primary scope. The core design questions are resolved: one population-generation history, one established population per condition, two sampling periods and three passage-matched measurements per condition per period. The distinct sampling periods support an assessment of repeatability, not independent LOY-generation replication. Population origins and matching for parental and addback groups remain separate extension-specific questions and do not block the core descriptive analysis.

The user clarification resolves conflicting descriptions in `analysis/HT1080_BATCHES/analysis_HT1080_BATCHES.Rmd` (near lines 47-58). The manifest in `config/sample_manifest_draft.csv` records passage matching, serial-measurement status, within-period series IDs, shared starting-population IDs across periods and the common nucleofection origin. Its draft name reflects unconfirmed metadata for contextual/extension groups, not an unresolved core design. IDs are analysis labels rather than asserted laboratory barcodes. Sample order follows the verified raw matrix. [Design record](config/README.md)

## 2. Audit and rebuild the analysis objects

- Freeze a sample manifest with unique sample key, raw filename, batch, treatment, culture/experimental-unit ID, generation-experiment ID, passage, optional matched block, transduction/construct and available preparation/sequencing metadata. Record missing values explicitly.
- Verify nonnegative integer counts, unique exact gene IDs, sample alignment, identical counting features across inputs, library sizes and featureCounts assignment fractions. Join by identifiers, never positional assumptions alone.
- Audit the counting command, STAR reference/PAR masking, strandedness and count units. The inspected featureCounts 2.1.1 headers use `-p` without `--countReadPairs`, and omit `-s`; this indicates read counting and unstranded assignment under the documented defaults. Verify all headers and library preparation records before calling these fragment counts. A matrix cannot retrospectively resolve strandedness or mapping ambiguity. An upstream recount becomes a separate task only if the audit establishes a material incompatibility. [Subread 2.1.1 manual](https://subread.sourceforge.net/SubreadUsersGuide.pdf)
- Derive symbols, gene types, chromosomes and coordinates from the actual GENCODE v49 counting GTF. Keep the exact counting gene ID as the primary key and a separate stable-ID mapping. Retain count features lacking a symbol or Entrez ID. Audit duplicate symbols and mapping failures.
- Do not use the existing combined-matrix annotation as the canonical annotation: its README identifies EnsDb v86. Do not select the newest AnnotationHub record at render time. Freeze any supplementary Entrez, HGNC, pathway and regulon mappings separately.
- Review the documented May addback relabeling and retain original/final labels in the manifest. Source filenames already reflect the relabeling; do not apply the permutation again or relabel samples solely from expected expression.

Create one immutable raw-count object, one annotation table, and explicitly named analysis subsets. A `SummarizedExperiment` can hold counts, annotation and metadata; model-specific `DGEList` objects derive from it.

## 3. QC, filtering and normalization

Inspect library complexity proxies available from counts, assignment summaries, detected genes, count distributions, sample correlations, MDS/PCA, hierarchical clustering and high-abundance gene contributions. Plot batch, treatment, culture/block, passage and transduction labels together. Add mitochondrial, ribosomal, MSY, Yp/Yq and PAR summaries. Cell-cycle and interferon/stress scores are biological diagnostics, not automatic adjustment covariates.

Use expression filtering defined before DE. For the primary 12-library comparison, propose `edgeR::filterByExpr` with treatment-by-batch groups and `min.total.count = 15` for continuity. The three measurements in each group support expression detection across passages; they do not establish three independent replicates. Save the complete keep/drop table and effective thresholds, and revisit the filter if the modeling input is summarized at population level. Parental and addback samples do not determine the primary tested universe. Sensitivities at 50 and 100 total counts are optional diagnostics; compare methods on explicitly shared tested genes. [edgeR user guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

For the main trans-expression analysis, use expressed autosomal and non-PAR X genes. Estimate TMM normalization on this prespecified background. Analyze MSY, PAR and mitochondrial genes in a separate retained-feature branch, carrying the same sample exposure offsets when evaluating dosage. Do not independently rescale the small Y-gene panel. Keep ribosomal protein genes in the primary universe; an RPS/RPL-excluded pathway run can reproduce the older sensitivity.

Produce PCA before and, if informative, after removing a fitted sampling-period contribution while preserving the condition contrast. Corrected values are for visualization; retain the original normalized expression for effect summaries. Sampling period may include biological drift as well as technical batch differences, so do not label every period effect technical. Use leave-one-passage-out diagnostics before making exclusion decisions. Exclusions require a documented sample-level reason and a complete sensitivity comparison.

All normalized abundance results are relative. Without suitable spike-ins or independent cell-normalized measurements, this design cannot establish global changes in RNA output per cell.

## 4. Expression differences and repeatability

**Primary estimand:** the average passage-matched log-expression difference between these established LOY and sgCTRL populations, giving August and May equal weight. For each gene, calculate the difference in normalized log2 expression for LOY versus sgCTRL within each of the six passage/day matches. Summarize the three differences in each period, then average the two period means. Preserve all six observed differences, normalized abundances, within-period ranges and direction agreement. Document the log-CPM prior-count convention and examine low-count sensitivity.

The primary tables contain effect sizes and observed variability, without gene-level significance or FDR claims about independently generated LOY populations. Ranges or standard deviations describe variation across the observed passages; they are not confidence intervals for a general LOY effect. A magnitude cutoff such as absolute log2 difference >= 0.5 can label a descriptive candidate list, with its meaning stated explicitly. This scope follows the experimental structure; a descriptive fold-change analysis is also an established option when independent biological replication is absent. [edgeR user guide, section 2.13](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

Planned comparisons:

- LOY versus sgCTRL, averaged across matched passages and equally across the two sampling periods.
- August-only and May-only differences, displayed together to assess repeatability and drift.
- Changes in the LOY-minus-sgCTRL difference between periods and over passage order, described without an independent-replication heterogeneity test.
- August CTRL versus PARENTAL and LOY versus PARENTAL as contextual expression comparisons; preserve the separate status of unconfirmed parental matching.

Use leave-one-passage-match-out summaries, median-versus-mean effects, expression-filter sensitivity and selected normalization/weighting checks to identify fragile candidates. Resampling six passage matches cannot supply a confidence interval across independent nucleofection experiments. Neither a population fixed effect, a random intercept nor a selected dispersion can separate population history from LOY status in this design.

**Role of the original DE methods:** limma-voom point estimates and quality-weight diagnostics can be retained as sensitivity/comparison outputs. Previously saved limma, edgeR, DESeq2 and TREAT calls are historical model-based results, not independently replicated confirmation. Do not make their P values the primary discovery criterion. Any future within-population inferential model must state its temporal/covariance assumptions and restricted estimand before fitting; independent-population DE would require additional independently generated populations. There is no automatic mixed-model substitution that resolves the present design.

For exploratory gene-set or regulon procedures that return P values, preserve their actual null hypotheses and declared multiplicity family. Expanded fgsea and CAMERA use global BH across all eligible MSigDB collections within each method and contrast; collection/subcollection BH values are secondary labeled outputs. These values must not be relabeled as gene-level biological-replication FDR.

## 5. Pathways and regulators to retain

The core uses descriptive effect rankings and sample pathway/TF profiles. Genes and regulon targets are not experimental replicates. The user subsequently requested all MSigDB collections and CAMERA. CAMERA is now an explicitly exploratory fixed-model extension: six passage/day block intercepts and two period-specific LOY effects, voom weights, and an equal-period combined contrast. Residual independence is an assumption, not a property established by serial passage blocking; no population-level inference is claimed.

| Analysis | Tool and input | Planned improvement |
|---|---|---|
| Ranked enrichment | `fgsea` on complete combined and period-specific effect rankings | All available human H and C1-C9 collections/subcollections; 15-500 represented genes; complete eligibility and leading-edge exports |
| Competitive gene-set test | `limma::camera` on voom expression/weights with matched passage fixed effects | Same set/gene universe; fixed inter-gene correlation 0.01 plus estimated-correlation sensitivity; exploratory independent-residual assumption, four residual degrees of freedom |
| Targeted signature response | ROAST/fry, deferred inferential extension | Requires a defensible model; descriptive signature trajectories proceed now |
| Sample pathway profiles | One prespecified GSVA/ssGSEA or targeted signature scorer; PROGENy for signaling | Show passage trajectories and matched score differences, with the gene universe and scale frozen |
| TF activity | `decoupleR` with frozen CollecTRI; DoRothEA as a resource sensitivity | Contrast-level ULM for prioritization and descriptive sample activities; target-pattern P values are explicitly conditional |
| Master regulators | Corrected msVIPER as an optional secondary module | Descriptive prioritization with audited networks; no passage-bootstrap claims of independent population replication |

The fgsea workflow and CollecTRI/ULM implementation are documented by their maintainers. Record exact gene-set/regulon versions, downloads and hashes; a package version alone does not freeze its data. [fgsea documentation](https://bioconductor.org/packages/release/bioc/html/fgsea.html), [decoupleR TF activity vignette](https://saezlab.github.io/decoupleR/articles/tf_bk.html)

State the question beside each pathway result: fgsea describes enrichment along a fixed ranking; CAMERA asks whether a set is more associated than the background; ROAST asks whether the set contains a response to a modeled contrast. CAMERA is now run with its covariance limitation made explicit; ROAST/fry remains deferred. Inter-gene correlation adjustment does not estimate within-population serial correlation. Never describe preranked P values as sample-label randomization. [limma gene-set documentation](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)

Group overlapping pathways by their leading-edge genes while preserving complete tables. If ORA is added for a targeted list, define the eligible expressed/mapped universe and assess abundance/length selection bias. Record the actual multiplicity procedure; classic GSEA's exploratory FDR convention is not interchangeable with fgsea BH-adjusted p-values. These interpretation rules follow sections 6.2-6.5 of the supplied review.

For leading TFs, show activating/repressing target contributions, shared targets among candidate regulators and sensitivity to uncertain edges. A ULM target-pattern p-value is conditional on the input signature and network; the number of targets is not the number of biological replicates. Keep this separate from testing group differences in sample activity scores. Treat inferred TF responses downstream of UTY/KDM5D as candidate mechanisms, not direct measurements of chromatin-regulator catalytic activity (review section 7).

The requested chromatin follow-up is implemented in `R/04_chromatin.R`. It freezes source descriptions and species, separates perturbation responses from occupancy and component sets, compares both UP/DN response arms across periods, timing matches and leave-one-match-out summaries, and exports 34 regulator profiles. The focused 37-set catalogue includes non-enriched comparisons; 34 sets meet the original size limits. Five post hoc sensitivities remove specified MMP1/MMP3, IFN, matrix/EMT and cell-cycle genes from both sets and their background. Original all-collection BH values remain the discovery reference; the subset sensitivities have nominal P values only and do not estimate causal adjustment or absolute chromatin activity. Figure 7 and Supplementary Figure 4 present the completed analysis.

Keep locus-level DE keyed by gene ID. Where a tool requires unique symbols, create one deterministic mapping before consulting differential statistics, for example a documented GTF-informed representative with highest mean abundance as a tie-breaker. Use that same mapping throughout the run. Do not choose the locus with the largest absolute t statistic separately in each analysis.

For VIPER, verify the installed API and use the supported combinatorial `level` argument; export actual synergy results separately from individual enrichment p-values. Match the sampling unit to the confirmed experimental design and record skipped/failed stages explicitly. Treat context-mismatched ARACNe networks as exploratory, with target coverage shown. [VIPER manual](https://bioconductor.org/packages/release/bioc/manuals/viper/man/viper.pdf)

## 6. Additions ranked by scientific value

**A. Agreement and influence across sampling periods — core.** Build one gene-level evidence table containing combined and period-specific effects, all matched passage differences and influence measures. Evaluate agreement of pathways and TF activities as well as individual genes. Define a signature using August only, freeze it, and assess it in May as repeatability across sampling occasions. The shared populations make this a retrospective temporal assessment, not independent biological validation; using separate RNA matrices does not change that limitation.

Freeze the signature genes, directions, scoring background and any fitted preprocessing in the discovery data. Check whether sample scores change when batches are scored jointly versus separately; some scoring methods depend on cohort composition. Do not treat separately calculated GSVA/ssGSEA scores or enrichment NES values as automatically comparable absolute activity units (review sections 6.5-6.6).

**B. Residual Y expression and paralog dosage — core.** Compare MSY, Yp/Yq and PAR profiles by batch and treatment. Existing HT1080 reports give different interpretations of residual Yp expression; revisit the underlying gene-level evidence without assuming either interpretation. Counts alone cannot distinguish residual Y-positive cells, partial Y retention and X/Y mapping ambiguity. Select any decisive BAM coverage or orthogonal DNA/FISH/ddPCR checks as explicit follow-up inputs.

For X/Y pairs, report each gene's abundance and log2FC, and a baseline-weighted total-abundance estimate where quantification supports it. A proposed relative-total measure is `(A_X_LOY + A_Y_LOY) / (A_X_CTRL + A_Y_CTRL)`, using consistently normalized, length-aware abundance estimates. Raw gene counts have different length/mappability biases; count-only pair totals must be labeled as proxies. Show individual passage measurements and their observed variability; do not bootstrap passages into independent-population confidence intervals. Do not label `log2FC_X + log2FC_Y = 0` as perfect compensation. RNA abundance alone does not demonstrate functional or protein-level compensation.

If locus ambiguity materially changes a conclusion, scope a separate reference-aware quantification comparison using the original reads and a compatible decoy-aware Salmon reference. Gene-level import would use supported tximport/tximeta handling of estimated counts and length offsets, or an appropriate abundance-derived count route without double length correction. This cannot be reconstructed from the current gene-count matrix. Compare locus attribution and uncertainty; inferential draws are not additional biological samples. [tximport documentation](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)

**C. DNA copy-number context — core if samples can be matched.** Map RNA genes and DNA segments using assembly-compatible coordinates. Prefer actual DNA segments; use proper centromere/cytoband boundaries for arm summaries. Never split a chromosome at its midpoint or hard-code an arm as concordant. Display DNA dosage and RNA effects together, and repeat selected enrichment analyses after excluding affected segments with a correspondingly revised background. If DNA has only one specimen per condition, this is descriptive context: do not invent per-RNA-sample CN covariates or claim a replicated CN-adjusted causal effect. Genomic correlation should not be tested as if adjacent genes were independent observations.

**D. UTY/KDM5D addback rescue — highest-priority extension.** Reuse the September addback pipeline's input/provenance work after reviewing its helpers. Compare each addback with LOY+LXHA; compare LXHA with naive LOY to measure the vector/transduction difference. Fit the four plasmid arms as a UTY x KDM5D factorial design, including the interaction, with blocks only if confirmed. Check transgene expression and co-perturbations before mechanistic attribution.

Use an August-derived signature to assess May rescue without reusing May RNA observations for signature selection; the shared population history still prevents calling August an independently generated validation cohort. Report directional reversal and residual differences from CTRL. Inferential intervals, factorial interaction tests and equivalence claims require confirmed replication of the addback interventions and an appropriate model; otherwise show descriptive effects. Returning to nonsignificance is not evidence of restoration. Rescue ratios are unstable for small baseline effects, and shared controls create covariance that must be included in any uncertainty calculation. Distinguish a factorial expression contrast from a demonstrated molecular interaction.

**E. ATAC-RNA-regulator concordance — next extension.** Use the existing September ATAC promoter and TF outputs if biological samples and contrasts match. Test concordant promoter accessibility/expression and agreement among RNA-derived activity, chromVAR and TOBIAS. Document peak-to-gene mapping, multiple-promoter aggregation and the common tested background. Report support for regulatory hypotheses, not causal proof. Distal enhancer links require additional evidence; nearest-gene assignment alone is insufficient.

**F. External perturbation or dependency support — targeted extension.** Compare prioritized programs with available independent Y-gene knockout/dosage datasets and HT1080 screen results after checking metadata, sample overlap and appropriate gene universes. Separate expression response from essentiality. Use other cell lines to assess generalization after their own within-cell-line designs and QC, without directly merging raw counts from incompatible studies.

**G. Focused curated interaction/complex context — optional extension.** Following review section 8, construct a small network around supported rescue targets or regulatory candidates. Select a resource for the actual question: physical interactions, complex membership, directed regulation or functional association. Record edge provenance and distinguish measured input genes from added neighbors. Assess sensitivity to evidence channels and, for network-enrichment claims, an appropriate measurable/degree-aware background. A connected hub is a follow-up candidate; network degree and shared pathway annotations do not establish independent mechanistic validation.

Do not make small-sample de novo coexpression networks, inferred bulk-RNA CN, or prediction classifiers primary additions. Splicing, isoform usage, fusions and allele-specific expression need junction/transcript/alignment inputs beyond this gene-count matrix. They can be scoped later if a specific result motivates them.

## 7. Tools and implementation layout

Use the verified R 4.5.2 installation explicitly: `C:/Program Files/R/R-4.5.2/bin/Rscript.exe`. The current PATH resolves a different R installation (4.1.3). Freeze a compatible environment rather than upgrading everything during implementation.

| Purpose | Available tools observed on 11 September |
|---|---|
| Count models | limma 3.66.0; edgeR 4.8.2; DESeq2 1.50.2 |
| Annotation and data objects | rtracklayer 1.70.1; GenomicRanges 1.62.1; SummarizedExperiment 1.40.0 |
| Pathways | fgsea 1.36.2; msigdbr 26.1.1; GSVA 2.4.9; progeny 1.32.0 |
| Regulators | decoupleR 2.16.0; viper 1.44.0 |
| Conditional repeated-measure analysis | variancePartition 1.40.2 |
| Reproducibility and focused checks | renv 1.2.3; testthat 3.3.2 |
| Pipeline execution | Recommend `targets`; not currently installed, would be a new dependency |

Use `targets` to track inputs, configurations and downstream dependencies, and `renv` to record packages. Ensure files/resources and configuration values are explicit dependencies so changes invalidate affected results. A `.done` marker alone is insufficient. Keep a project-local persistent cache separate from curated deliverables, with its location configurable for Dropbox performance. [targets manual](https://books.ropensci.org/targets/), [renv documentation](https://rstudio.github.io/renv/articles/renv.html)

Proposed implementation structure; the core descriptive design is resolved:

```text
SEPT_11_2026_REANALYZE/
  PLAN.md
  FIGURES.md              # publication style, figure catalogue, export requirements
  FIGURE_CAPTIONS.md      # draft captions; final text generated from verified results
  input_inventory.csv
  config/                 # sample manifest, contrasts, thresholds, resources
  R/                      # named ingest, QC, DE, pathway and integration functions
  _targets.R
  renv.lock
  report/analysis.Rmd      # reads saved results; performs no model refits
  outputs/tables/
  outputs/models/
  outputs/figures/         # PDF + SVG + 600 dpi PNG; source data and figure manifest
  outputs/captions/        # caption per figure; combined Markdown and DOCX
  outputs/report/
  provenance/             # input/resource hashes, session info, decisions
  tests/testthat/         # focused tests for scientific failure modes
```

One computation owns each result file. Sensitivity runs receive distinct method/contrast names. Reuse or extract audited functions instead of copying whole notebooks. Save seed 42 and stage-specific random seeds for reproducible resampling. Generate methods text from the actual configuration and model objects, including which samples contributed to filtering, normalization and fitting.

## 8. Execution order and completion criteria

The user-requested next investigation is detailed in [HDAC12_PLAN.md](HDAC12_PLAN.md). It proposes source-signature auditing, HDAC1-versus-HDAC2 gene specificity, independent genetic/catalytic/inhibitor comparisons, complex context and conditional addback integration. This is an additional-work plan; those new analyses are not represented as completed results.

| Phase | Work | Reviewable deliverable |
|---|---|---|
| 1: establish the core | Audit inputs/annotation using confirmed metadata; QC; matched expression differences and passage sensitivities | Sample/design audit, count/annotation objects, effect tables and QC report |
| 2: interpret the core | Exploratory pathways/TFs, temporal repeatability, Y/paralog dosage, matched CN context | Main biological report and gene/pathway evidence tables |
| 3: test mechanisms | Selected addback, ATAC and external comparisons | Clearly separated integration results and a prioritized validation shortlist |

Meaningful implementation checks should cover sample reordering, duplicate-symbol loci, counting-GTF joins, contrast direction and estimability, correct chromosome arms, baseline-weighted paralog calculations, signature sample independence, and explicit handling of failed optional stages. Test these scientific failure modes with small known examples rather than reproducing implementation code in tests.

Completion means: every result links to its source counts, annotation, model, contrast and tested universe; all figures satisfy FIGURES.md, have source tables and PDF/SVG/PNG exports, and have complete captions with verified sample sizes, statistics and uncertainty definitions; no unexplained sample exclusions or arm overrides remain; primary and sensitivity results are distinguishable; and the report renders from saved objects. Changing a declared input must invalidate affected results, while re-rendering the report must not overwrite models. Compare the new August-only fit with the old report and explain differences caused by annotation, filtering, design or bug fixes. Matching the old number of DEGs is not a success criterion.

The final priority table should show the individual evidence dimensions: batch agreement, effect size, statistical support, unit influence, method sensitivity, DNA dosage, TF/pathway membership, addback reversal and ATAC support. Do not collapse these into an uncalibrated majority-vote significance score.
