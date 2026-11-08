# HDAC1/2 investigation plan

Prepared 12 September; updated 13 September 2026. Steps 1, 2 and the direct public-data benchmark in step 3 are complete: [reference audit](outputs/report/HDAC12_reference_audit.html), [gene drivers](outputs/report/HDAC12_gene_drivers.html), [external perturbation comparison](outputs/report/HDAC12_external_comparison.html). Step 1 retains documented source-recovery limits. Optional LINCS and steps 4-5 remain unrun.

The GSE193885 benchmark narrows the interpretation: HDAC2 knockout weakly agrees with LOY, but active-transgene controls agree more strongly and the inactive-versus-active HDAC2 ranking is weakly inverse. Both response arms frequently decrease. The evidence supports a persistent partial gene pattern without establishing reduced HDAC2 catalytic activity. Drug dose/time labels and independent source-clone replication remain unresolved. Figure 9 and Supplementary Figure 6 report the comparison and its controls; complex analysis is separately scoped step 4.

The objective is to determine which parts of the established HT1080 LOY expression difference resemble HDAC1/2 perturbation, whether the evidence is more specific to HDAC2 or a shared complex, and which observations would distinguish catalytic changes from altered recruitment or a broader cell-state response.

## Starting evidence and questions

The completed analysis shows a stronger HDAC2-knockdown response pattern than HDAC1 alone. For the source UP-minus-DN mean log2 expression contrast, the combined values are +0.349 for HDAC2, +0.232 for HDAC1/2 double knockdown and +0.030 for HDAC1. The HDAC1 contrast changes sign between August and May (-0.042, +0.102); the HDAC2 contrast is positive in both (+0.204, +0.494). The HDAC2 down arm is lower in both periods and persists after the existing IFN, matrix/EMT and cell-cycle exclusions. HDAC1 and HDAC2 RNA are modestly lower, while CHD4 is modestly higher. These observations motivate specificity testing; they do not measure enzyme activity or complex assembly.

All core observations remain repeated samples of one established population per condition from one shared generation history. The six timing matches cannot support a general causal LOY effect or independent-population validation. The new work is exploratory and follows inspection of the original enrichment results.

## 1. Audit the HDAC reference evidence and define the comparison catalogue

**Question:** What perturbation does each signature actually represent, and how much independent evidence is available?

- Audit all Senese HDAC1, HDAC2, HDAC1/2 and HDAC3 UP/DN arms. Recover the source table or underlying data if available; record RNAi timing, controls, probe mapping, selection thresholds, biological versus technical replication and whether a result was validated independently.
- The original paper describes a comparatively limited HDAC2-only response, despite the presence of HDAC2 gene lists in MSigDB. Resolve the relationship between those lists and the original analysis before interpreting a strong HT1080 enrichment as a calibrated HDAC2 loss signature. If the original supplement cannot be recovered, retain the frozen MSigDB definitions and explicitly record that verification limit. [Senese et al., 2007](https://pubmed.ncbi.nlm.nih.gov/17470557/)
- Review the additional HDAC-related sets already in the frozen resource, including ROPERO_HDAC2_TARGETS, HELLER_HDAC_TARGETS_UP/DN, PEART_HDAC_PROLIFERATION_CLUSTER_UP/DN and MARKS_HDAC_TARGETS_UP/DN. Separate genetic perturbation, inhibitor exposure, combined drug exposure, literature summaries, physical occupancy, complex components and observational coexpression.
- Treat multiple gene lists from one experiment as one source of biological evidence. Drug-plus-decitabine sets are not HDAC-only experiments; HDAC3 and unrelated growth-arrest perturbations are specificity comparisons, not guaranteed null controls.
- The C9 GENES_CORRELATED_WITH_HDAC1_DELETION set is derived from expression associations with CRISPR dependency across cell lines. It is not a measured transcriptional response following HDAC1 deletion and must remain in a separate evidence category. [MSigDB definition](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GENES_CORRELATED_WITH_HDAC1_DELETION.html)

**Output:** a frozen reference manifest linking every set to source study, intervention, direction, cell type, time, membership, eligibility and limitations. Expand the existing chromatin resource metadata only for explicitly added records; preserve the prior resource hashes.

## 2. Identify the genes and processes responsible for the pattern

**Question:** Is the signal HDAC2-specific, shared with HDAC1/2, or driven by a small group of downstream genes?

- Partition source memberships into HDAC2-only, HDAC1-only, concordant shared, conflicting-direction and double-knockdown-only groups, retaining the UP and DN arms. These labels describe specificity within the reference catalogue, not proof of isoform-specific regulation in HT1080.
- Examine every represented member, with leading-edge membership as an annotation rather than the sole analysis universe. Record exact gene identity, average abundance, six matched differences, August/May effects, leave-one-match-out ranges and reference-program overlaps.
- Extend the completed passage sensitivity with gene-contribution diagnostics: report the fraction of absolute effect contributed by the largest genes and repeat summaries after removing the top 1, 3 and 5 contributors under a fixed rule. Retain the original results beside each diagnostic.
- Annotate the remaining genes against frozen process definitions: cell cycle/mitosis, apoptosis/p53, differentiation, adhesion, immune response and DNA damage/repair. Display both the annotated overlap and the expression of genes outside those annotations. First reuse existing reference masks; add only justified, declared masks rather than trying exclusions until a preferred result appears.
- Inspect abundance, prior-count sensitivity and chromosomal concentration. Use leave-one-chromosome-out summaries where a cluster dominates. RNA chromosome clustering is not evidence of a DNA copy-number alteration; DNA adjustment remains conditional on compatible, matched specimens.
- Freeze complete source-derived UP/DN programs for sample-level signed percentile scoring on a fixed, expressed-gene background. Plot all 12 libraries and the six within-period differences. Do not select a score's genes using both periods and then call May validation.

**Output:** a complete gene evidence table, membership-overlap figure, passage-score plots and a short list of robust downstream readouts. The list should include source-UP genes that fail to change as expected, as well as the strongest source-DN examples.

## 3. Compare independent genetic and inhibitor perturbations

**Question:** Does LOY resemble loss of HDAC protein, impaired catalytic activity, or a nonspecific treatment response?

**First external candidate: GSE193885.** The public GEO record describes human HAP1 RNA-seq comparing class I HDAC knockouts, active/inactive transgenes, and inhibitor exposure, with three reported replicates per condition. This provides a useful design for separating protein loss from catalytic inhibition. The available count files, exact reference contrasts, clone histories and independence of those replicates still need audit before analysis; metadata discovery is not a completed benchmark. [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193885)

- Reconstruct each valid perturbation-versus-matched-control contrast separately within the external study. Use its verified design, annotation and counting units. Never normalize or batch-correct the public and HT1080 samples into a single inferential cohort.
- Prioritize HDAC1 and HDAC2 loss, their catalytic mutants and corresponding active-protein controls. Include HDAC3/8 and inhibitor contrasts as comparators where the design permits. Do not equate a stable knockout with acute inhibition.
- Compare complete signed effect rankings on their common measurable genes. Report directional concordance and effect/rank correlation descriptively; do not use naive gene-independent correlation P values. Also test eligible externally defined UP/DN sets with fgsea and CAMERA, retaining both response arms and the existing overlap sensitivities.
- Compare August and May separately, as well as the equally weighted combined HT1080 effect. Report source-study and cell-context consistency without counting overlapping signatures as independent replications.
- Use LINCS/CMap as a secondary extension only after the direct RNA-seq comparison: separate genetic perturbations from compounds, track cell line, time, dose and signature quality, and distinguish measured L1000 landmark genes from inferred genes. Include proliferation/cytotoxicity comparators and aggregate at the perturbation/study level. Similarity to an inhibitor profile is a mechanistic clue, not a prediction of drug sensitivity. Download public references for local comparison where practical. [Official LINCS workflow](https://lincsproject.org/LINCS/tools/workflows/query-a-gene-expression-signature-against-the-lincs-11000-data)

**Output:** a perturbation-concordance matrix with rows for source experiment/condition and columns for August, May, combined effects and declared exclusions. Keep genetic loss, catalytic mutants and drug exposure visibly separate.

## 4. Examine HDAC-containing complexes and the ZMYND8 connection

**Question:** Is the expression pattern more consistent with a particular recruitment complex?

- Curate experimentally supported NuRD, SIN3 and CoREST memberships from primary references or a versioned complex resource. Include HDAC1/2 and shared components without counting them repeatedly as independent support.
- Expand regulator expression context to CHD3/4, MTA1/2/3, MBD2/3, GATAD2A/B, RBBP4/7, SIN3A/B, SUDS3, SAP30, RCOR1/2/3 and KDM1A, subject to verified membership and expressed-gene coverage. Display shared and complex-distinguishing subunits separately.
- Include ZMYND8, all four KDM5 paralogs and GATAD1. The primary DNA-repair study identified a KDM5A-ZMYND8-NuRD connection, including HDAC1/2 and CHD4; this motivates a broader comparison and does not demonstrate the same mechanism in HT1080 LOY. [Gong et al., 2017](https://rupress.org/jcb/article/216/7/1959/39050/Histone-demethylase-KDM5A-regulates-the-ZMYND8)
- Compare curated complex-perturbation signatures when available, with source context and UP/DN semantics preserved. Separate component RNA from downstream target responses.
- For candidate downstream loci, annotate public HDAC1/2, CHD4 or ZMYND8 occupancy only after checking cell type, genome build, assay and peak-to-gene mapping. Context-mismatched binding provides supporting annotation rather than direct HT1080 occupancy. Do not infer a complex from a 12-sample coexpression network or assign protein stoichiometry from RNA.

**Output:** a compact source-annotated complex diagram and regulator/target evidence matrix. Explicitly retain missing protein, binding and catalytic evidence.

## 5. Test whether existing addback data can distinguish mechanisms

**Question:** Do KDM5D, UTY or their combination reverse an independently defined HDAC-associated program?

An existing HT1080 addback reanalysis was located. Integration is a conditional extension: audit construct identity and expression, relabel provenance, empty-vector controls, sample overlap, passage matching and independent experimental units before reusing its reported contrasts.

- Freeze HDAC programs from the public sources before examining addback outcomes. If any HT1080-selected genes are used, identify shared samples and avoid claims of independent validation.
- Use the appropriate empty-vector control for each valid construct comparison. Plot source-UP and source-DN arms separately, alongside individual candidate genes; direction toward the relevant CTRL reference is more informative than a single unsigned enrichment statistic.
- Compare KDM5D, UTY and dual addback only where the actual design supports those contrasts. Report partial reversal descriptively; no automatic full-rescue label or percent-rescue ratio with an unstable denominator.
- Reuse compatible ATAC data only after specimen and assembly checks. Accessibility alone does not measure histone acetylation or identify the responsible HDAC complex.

**Output:** a clearly separated addback comparison, or a precise record of the unresolved design dependency. The first four workstreams can finish without this extension.

## Tools, statistical rules and implementation

Use the existing R 4.5.2 environment: edgeR/limma for validated external expression contrasts, fgsea and CAMERA for gene-set comparisons, existing matched-difference/rank-score helpers for HT1080, GenomicRanges for compatible coordinates, and ggplot2/patchwork for figures. Use GEOquery or direct public processed files according to their availability; a new package is not required merely to read a downloaded count table. Use documented LINCS data access only if that extension is pursued.

Preserve the original 23,696-set global BH results. Define and freeze any new secondary testing family before running it; report exploratory BH over the complete declared family per method and contrast, including negative results. Membership partitions smaller than the minimum receive descriptive summaries only. Gene-removal diagnostics retain nominal statistics without a confirmatory FDR claim. Fixed-correlation CAMERA and its estimated-correlation sensitivity remain separate; neither resolves serial dependence in the HT1080 populations. No bootstrap over genes or six timing matches is presented as independent-population uncertainty.

Proposed implementation: `config/hdac12_followup.yml`, source files under `resources/hdac12/`, `R/04_hdac12.R`, `R/05_hdac12_figures.R`, and outputs prefixed `HDAC12_`. Record source versions, accessions, hashes, model formulas, effective universes and skipped comparisons. Keep source-audit, computation and report rendering separate. Original counts remain untouched.

## Figure and completion plan

Completed Figure 8 shows source-signature relationships, matched program differences, source-UP/DN gene evidence and gene influence. Completed Figure 9 shows external perturbation concordance, both response arms, fgsea/CAMERA and the eight existing readouts. Supplementary Figures 5-6 provide complete sample scores, process annotations, source QC and exclusion sensitivity. Complex context and addback require separately completed analyses before receiving figures.

Use the existing Nature-style standard: 183 mm width, Arial, editable PDF/SVG, 600 dpi PNG, source data for every panel and complete captions. Captions must identify repeated population measurements, source-defined gene selection, score direction, exclusions and the inferential limits. Keep complete evidence in source tables even when figures show selected examples.

Run workstreams 1 and 2 first, then 3, with 4 supplying mechanistic context. Pursue 5 after the addback audit. Completion should answer: which reference evidence is reliable; what distinguishes HDAC2 from HDAC1/shared responses; whether independent perturbations agree; which genes carry the result; and what remains unmeasured. A failed external match is a useful result and should narrow the claim to a source-specific expression overlap.

The resulting experimental shortlist should prioritize protein measurements of HDAC1/2 and relevant complex components, quantitative histone acetylation, and selected-locus chromatin profiling. If new perturbations are undertaken, compare HDAC1, HDAC2 and combined perturbation in CTRL and LOY with appropriate controls and independently generated populations. Define early transcriptional readouts and later growth/viability readouts separately. The shortlist would support a subsequent assay design with replication based on independent populations.
