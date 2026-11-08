# HT1080 LOY: current understanding of the transcriptional phenotype

**Updated: 14 September 2026.** Synthesis of the completed `SEPT_11_2026_REANALYZE` analyses, including the expanded GSEA/CAMERA screen, ZMYND8 investigation, HDAC1/2 source audit and external benchmark, and PRC2/epigenetic deep dive. This document consolidates existing results; it introduces no new statistical tests.

## Current biological interpretation

The established HT1080 LOY population has a recurring expression phenotype characterized by **lower interferon/inflammatory and extracellular-matrix-associated RNA programs**, marked depletion of several Y-linked transcripts, and **selective changes in chromatin-associated gene programs**. Important directions recur in August 2025 and May 2026, although the transcriptome also changes substantially between periods.

Chromatin regulation is a plausible contributor, but the signatures do not identify one responsible enzyme or a uniform change in repression. The positive inferred ZMYND8 score depends strongly on MMP1/MMP3; the HDAC2-associated program persists, but the external benchmark weakens a specific catalytic-loss explanation; and the PRC2-associated decreases and source H3K27me3-associated increases involve largely different genes. The best current model is **selective remodeling of expression programs in these established populations**, with the contribution of Y loss, population history and individual chromatin regulators still unresolved.

Here, “transcriptional phenotype” means the pattern of **steady-state, normalized RNA abundance**. These data do not directly measure transcription rates, total RNA per cell, chromatin marks, protein abundance or enzymatic activity.

## Experimental scope and strength of evidence

The primary comparison contains 12 RNA libraries: LOY and sgCTRL at matched passages 1–3 in each of two periods. Matching refers to the passage/day within each period. August and May sampled the **same established populations**, originating from the same nucleofected pool; they are not independently generated LOY experiments. There is one established population per condition. Three August parental libraries provide context only.

Primary effects are LOY-minus-sgCTRL log2-normalized expression differences, averaged within each period and then equally across periods. Positive values mean higher relative RNA in LOY. Across 25,526 eligible autosomal/non-PAR X genes, August–May effects have Pearson r = 0.490 and Spearman rho = 0.377; 64.3% have the same direction. That percentage includes small effects near zero. The analysis identifies 4,892 **descriptive candidates** with an absolute combined difference of at least 0.5, agreement between period directions, and retained direction after each single-match omission. These are not gene-level FDR discoveries.

Temporal recurrence and normalization/omission checks support the description of these populations. They do not supply independent biological replication or separate a causal LOY effect from other differences established during their history. “Persistent” below means recurrence in the sampled periods, without assuming uninterrupted culture or unchanged biology between them. [Design and primary report](outputs/report/analysis.html), [primary gene effects](outputs/tables/HT1080_LOY_vs_sgCTRL_primary_effects.csv), [run summary](provenance/run_summary.yml).

## 1. The broad expression phenotype

The expanded screen used MSigDB 2026.1.Hs, covering H and C1–C9. Of 35,361 resource sets, 23,696 met the original 15–500 represented-gene limits on a 25,461-symbol background. The combined fgsea ranking contains 808 sets with global Benjamini–Hochberg (BH) adjusted P < 0.05: 776 negative and 32 positive. Overlapping sets describe related genes, so these are not 808 independent biological mechanisms.

Representative broad programs are shown below. NES is the normalized enrichment score; negative values place set members toward genes with lower relative RNA in LOY. BH values retain the original correction across all 23,696 eligible sets in the combined contrast.

| Source program | Combined NES | Global BH adjusted P | Direction in August / May |
| --- | ---: | ---: | --- |
| Hallmark interferon alpha response | -1.897 | 0.00352 | Lower / lower |
| Hallmark interferon gamma response | -1.739 | 0.00715 | Lower / lower |
| Hallmark inflammatory response | -1.811 | 0.00177 | Lower / lower |
| Hallmark epithelial–mesenchymal transition | -1.741 | 0.00522 | Lower / lower |
| Reactome extracellular matrix organization | -1.748 | 0.000845 | Lower / lower |
| Hallmark angiogenesis | -1.876 | 0.0256 | Lower / lower |
| Hallmark IL6–JAK–STAT3 signaling | -1.754 | 0.0326 | Lower / lower |

All seven directions recur; not every set passes global BH < 0.05 separately in both periods. MMP1 and MMP3 are concrete matrix-associated examples: combined differences are -4.95 and -2.76 log2 units, respectively, with decreases in all six timing matches. MMP3 is relatively low abundance (mean CPM 2.56), which matters when choosing validation readouts.

Regulator footprints provide compatible, overlapping evidence: combined CollecTRI ULM scores include FOS -5.43, JUN -4.65, AP1 -4.55 and ETS2 -5.41. These are inferred target-expression patterns, not measurements of TF protein activity, and shared targets prevent treating them as independent explanations.

This pattern motivates investigation of altered inflammatory signaling and matrix interactions. It does not establish reduced invasion, metastasis, angiogenesis or immune infiltration. In particular, immune signatures in this cell-line experiment do not measure immune-cell abundance. The predominance of negative enrichment also does not demonstrate globally reduced transcription: normalization measures relative abundance.

A uniform proliferation explanation is unsupported. Hallmark E2F targets switch from positive enrichment in August (NES +1.482) to negative in May (-1.251), with combined NES -0.337 and BH = 1. Period-dependent programs need to remain distinct from recurring changes. [All GSEA results](outputs/tables/fgsea_all_collections.csv), [regulator scores](outputs/tables/collectri_ULM_effect_signatures.csv), [pathway figure](outputs/figures/Fig03_pathways_regulators.pdf).

## 2. Y-linked RNA loss and incomplete preservation of paralog totals

Non-PAR Yq RNA is strongly depleted: LOY/sgCTRL ratios of period-mean CPM are approximately 0.025 in August and 0.017 in May. More Yp signal remains, with corresponding ratios of 0.250 and 0.387. These are regional RNA abundance ratios, not estimates of the fraction of Y-positive cells. Residual signal could reflect retained Y material, population mixture or mapping behavior; RNA alone cannot distinguish them.

| Chromatin-related gene | Combined log2 RNA difference | Interpretation at the RNA level |
| --- | ---: | --- |
| KDM5D | -7.497 | Strong depletion in both periods |
| UTY | -6.435 | Strong depletion in both periods |
| KDM5C | -0.051 | Little change in the X-linked partner |
| KDM6A | +0.155 | Modest increase in the X-linked partner |

The length-adjusted X+Y relative-abundance proxies remain below sgCTRL in every timing match for KDM5C/KDM5D (ratios 0.73–0.79), KDM6A/UTY (0.52–0.76), and DDX3X/DDX3Y (0.40–0.54). Thus, the observed X-partner RNA does not preserve these measured pair totals. These proxies use gene lengths and TMM exposures; they are not transcript-resolved TPM, protein dosage or functional equivalence tests. Other pairs behave differently, so a blanket claim of absent or complete compensation is inappropriate.

Matched DNA confirmation remains important. Available WGS specimen correspondence and genome compatibility were unresolved, so DNA–RNA integration and copy-number adjustment were not performed. [Regional abundance](outputs/tables/Y_PAR_MT_period_summary.csv), [all gene effects](outputs/tables/HT1080_LOY_vs_sgCTRL_all_gene_effects.csv), [matched paralog totals](outputs/tables/paralog_matched_relative_total_proxy.csv), [Y/paralog figure](outputs/figures/Fig04_Y_paralog_dosage.pdf).

## 3. ZMYND8: a positive inferred score with a narrow driver base

ZMYND8 RNA changes little: +0.053 log2 units overall, approximately a 1.04-fold normalized ratio. Its inferred ULM score is positive (+4.457; target-pattern BH = 0.000529 across 687 regulators). The score is positive in both periods and all six timing matches, although the May-only BH is 0.075.

The regulon contains only eight represented targets, six with repressive edges. Lower expression of a repressive target raises the inferred score. Removing MMP1 reduces the combined score to +1.29; removing both MMP1 and MMP3 changes it to -0.70. The pattern therefore provides much stronger evidence for recurring MMP repression than for uniquely increased ZMYND8 activity. FOS, JUN and ETS1 also have network edges to these genes.

The direct KDM5D/JARID1D–ZMYND8 study provides a mechanistic hypothesis, but its focused co-occupied loci do not show coordinated repression here: CD44 and EGFR rise modestly (+0.263 and +0.325), whereas SNAI2 and VEGFA change little. Increased ZMYND8 protein activity or compensation for KDM5D loss has not been demonstrated. The direct reference is [Li et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4975651/); the previously supplied 2017 article concerns [ZMYND8–Drebrin structure/localization](https://pubmed.ncbi.nlm.nih.gov/28966017/).

**Current conclusion:** a recurring, MMP-dominated expression pattern is compatible with several regulatory mechanisms. ZMYND8 remains a candidate to perturb, with protein, occupancy and dependency evidence still missing. [ZMYND8 interpretation](provenance/ZMYND8_interpretation.md), [target omissions](outputs/tables/ZMYND8_target_omission_sensitivity.csv), [Figure 6](outputs/figures/Fig06_ZMYND8_followup.pdf).

## 4. HDAC1/2: persistent source-program similarity, weak catalytic specificity

The Senese HDAC2-knockdown DOWN set is among the strongest chromatin enrichments: NES -2.252, global BH 3.78e-7. The double-HDAC1/2-knockdown DOWN set is also negative (NES -2.099, BH 3.78e-7). These are genes decreased after the source perturbations in U2OS cells.

The complete HDAC2 DOWN arm contains 117 measured genes and averages -0.368 log2 units in LOY; the 113-gene UP arm averages only -0.019. The signed response score is positive in all six timing matches, but largely because the DOWN arm decreases. Direction survives gene-influence, abundance, normalization and predefined process-overlap checks. Removing the five largest absolute effects retains a DOWN-arm mean of -0.272. PLAC8, PPP2R2B, BEX1 and SORL1 provide recurring gene-level readouts. HDAC1 and HDAC2 RNA themselves are modestly lower (-0.300 and -0.205).

**The completed external benchmark changes the mechanistic interpretation.** GSE193885 supplies HAP1 HDAC knockouts, active and catalytically inactive transgenes, and drug conditions. The analysis modeled 63 libraries, reconstructed 27 contrasts and compared complete rankings across 14,881 common loci. For HDAC2:

| External HAP1 contrast | Spearman correlation with combined HT1080 LOY effects |
| --- | ---: |
| HDAC2 knockout versus wild type | +0.068 |
| Active HDAC2 transgene in knockout background versus wild type | +0.106 |
| Catalytically inactive versus active HDAC2 in knockout background | -0.046 |

All correlations are small. The active-transgene control agrees at least as well as knockout, while the inactive-versus-active comparison is weakly inverse. These observations **weaken a specific HDAC2 catalytic-loss explanation**. They do not exclude an HDAC contribution, but source-signature similarity alone cannot identify that contribution. Cross-cell-line context, unresolved source clone independence and some drug dose/time labels further limit calibration.

HDAC steps 1–3 are complete; optional LINCS, further complex analysis and addback remain unrun. The current conclusion incorporates step 3 and supersedes treating the earlier source signature as sufficient evidence for HDAC2 loss of function. [Gene-driver report](outputs/report/HDAC12_gene_drivers.html), [external benchmark report](outputs/report/HDAC12_external_comparison.html), [full-ranking concordance](outputs/tables/HDACEXT_cross_study_concordance.csv), [Figure 9](outputs/figures/Fig09_HDAC_external_comparison.pdf).

## 5. PRC2 and H3K27me3: distinct programs with opposing directions

The deeper epigenetic screen revisited 339 catalogue sets, of which 235 were eligible original tests; 79 eligible sets received focused analysis. Original global BH values were retained. The most informative PRC2/mark comparisons are:

| Source-defined signature | Combined NES | Global BH adjusted P |
| --- | ---: | ---: |
| EED knockdown UP, TIG3 | -2.063 | 3.27e-5 |
| EED knockdown DN, TIG3 | +0.651 | 1 |
| EZH2 knockdown UP, TIG3 | -1.783 | 0.00521 |
| SUZ12 knockdown UP, TIG3 | -1.610 | 0.0602 |
| H3K27me3 at high-CpG promoters, mouse MEF | +1.845 | 8.56e-6 |
| H3K27me3 at high-CpG promoters, mouse MCV6 | +1.818 | 5.93e-5 |

### The EED-depletion-induced program decreases

Of 154 measured genes induced by EED knockdown in the source experiment, 108 (70.1%) are lower in LOY; their mean difference is -0.297 log2 units. The source DN arm is near unchanged (+0.032). This is consistent with suppression of an EED-depletion-induced program, but does not establish greater PRC2 activity.

The EED direction remains after removing predefined interferon/matrix/cell-cycle genes, nucleosome genes or the ten largest absolute effects. Abundance restriction attenuates it: mean CPM >= 5 gives -0.121. SUZ12 UP becomes nearly unchanged under that restriction (-0.011), making EED a stronger descriptive priority. EED, EZH2 and SUZ12 source signatures come from the same study, so their agreement is not three independent source replications.

### Source H3K27me3-marked genes increase, largely elsewhere

The two positive promoter sets describe histone states in **mouse source cells**, mapped to human gene symbols. They do not measure H3K27me3 in HT1080. MCV6 is a partially reprogrammed population; this enrichment does not establish pluripotency in LOY.

Their 503-gene union overlaps the 154-gene EED UP arm by only six genes. The partition comprises 148 EED-only genes, six shared genes and 497 mark-only genes. Combined mean RNA differences are -0.288, -0.516 and +0.180, respectively. Even within the positive mark-only partition, 40.6% of genes decrease. These enrichments therefore describe different distributions of RNA changes, without implying uniform repression or derepression.

Illustrative post hoc readouts include GPR37, ADGRL4 and SERPINB2 among EED-only decreases; WNT7A, SDK2 and ADGRB1 among mark-only increases; and RYR2, CLSTN2 and POU3F3 as opposite-direction mark-only controls. They are candidate readouts, not validated direct PRC2 targets.

PRC2-component RNA is also mixed: EED +0.056, EZH2 -0.188 and SUZ12 -0.139. Neither these transcripts nor source-set enrichments establish globally increased PRC2 abundance or H3K27me3. [Deep-dive report and source definitions](outputs/report/epigenetic_deep_dive.html), [original enrichment evidence](outputs/tables/EPIGEN_original_enrichment.csv), [EED/mark partitions](outputs/tables/EPIGEN_EED_mark_effects.csv), [Figure 10](outputs/figures/Fig10_PRC2_gene_programs.pdf).

## 6. Other epigenetic observations refine the model

| Finding | Current interpretation |
| --- | --- |
| BMI1/MEL18-knockdown UP genes decrease; several DN arms also decrease | A Polycomb-associated program is present, but the paired directions do not support a simple increase in PRC1 activity. |
| SMARCE1 re-expression UP genes decrease; the DN arm also decreases | Source-program overlap does not demonstrate SMARCE1 inactivation. The source perturbation was re-expression. |
| Heller UP and DN sets both decrease | These arose from combined decitabine and trichostatin A treatment; HDAC and DNA-methylation contributions cannot be separated. |
| Liver-tumor H3K27me3 UP set is negative | “UP” denotes source promoters gaining H3K27me3, not RNA induction. This differs from the positive mouse promoter sets and emphasizes context. |
| Other promoter-state sets differ by source cell and CpG class | The H3K27me3 result cannot be generalized to all chromatin states. |
| Focused H3K9me3, DNMT1, KDM1A and KDM5B response sets lack combined global-BH support | No coherent enrichment evidence for a general change in these mechanisms; non-enrichment does not prove no involvement. |
| SETDB2 RNA rises (+0.624); TET2 falls (-0.376) | Additional mechanistic hypotheses, without a measured corresponding histone or DNA modification. |

Taken together, the chromatin findings support investigation of **specific gene programs and loci**. They do not justify assigning the entire phenotype to “HDAC2 off,” “PRC2 on,” increased chromatin repression or KDM5-family compensation. [Epigenetic context figure](outputs/figures/Fig11_epigenetic_context.pdf), [regulator RNA](outputs/tables/EPIGEN_regulator_RNA.csv).

## 7. How much inferential weight should GSEA and CAMERA carry?

GSEA evaluates enrichment in the descriptive gene ranking. Fixed-correlation CAMERA adds a competitive, precision-weighted expression model and finds 733 combined sets with global BH < 0.05. However, **none passes that threshold when CAMERA estimates within-set gene correlation**. The fixed-correlation results are sensitive to this assumption and should not be presented as robust statistical confirmation of GSEA.

This sensitivity does not erase the observed expression directions. It means the available model does not support the same discovery claims under the alternative correlation treatment. CAMERA has only four residual degrees of freedom, making correlation estimation limited. The methods also use different statistics, so disagreement is not a biological absence test.

Gene correlation and sample independence are separate issues. Genes in a set can move together, while passage samples can be dependent because they come from the same established population. The fixed passage-block model addresses timing matches but does not estimate serial population dependence. Neither method converts six timing matches into six independently generated LOY experiments.

Accordingly, prioritize effect direction, complete-set behavior, temporal recurrence, influence checks and source specificity. The focused chromatin selections are post hoc; their displayed BH values remain those of the original full screen. The external HDAC benchmark's 54 response sets form a separate analysis family. [Method comparison](outputs/tables/fgsea_camera_comparison.csv), [estimated-correlation sensitivity](outputs/tables/camera_estimated_correlation_sensitivity.csv), [expanded-enrichment figure](outputs/figures/Fig05_expanded_enrichment.pdf).

## 8. Questions carried forward from the completed analyses

These are proposed follow-ups, not completed experiments or validated mechanisms.

| Priority | Question | Previously proposed discriminating work |
| --- | --- | --- |
| 1 | Which changes follow LOY across independent generation histories? | Generate independent matched sgCTRL/LOY populations, confirm DNA-level Y status and align passage/harvest metadata. Resolve available WGS-to-RNA specimen correspondence. |
| 2 | Do the opposing PRC2-associated RNA programs reflect local chromatin changes? | Measure H3K27me3 localization with matched RNA, alongside a calibrated bulk mark measurement; include H3K27ac. Test EED-only and mark-only loci with opposite-direction and unchanged controls. |
| 3 | Which regulators are required for the observed programs? | Verify protein/target engagement and perturb ZMYND8, HDAC1/2 or PRC2 in the relevant HT1080 backgrounds. Read both response arms. A public PRC2 perturbation benchmark could first improve source specificity. |
| 4 | Which Y-linked losses contribute? | Use appropriately controlled, confirmed KDM5D or UTY addback to test reversal of predefined RNA and chromatin readouts. Distinguish rescue of a Y-gene contribution from a general whole-LOY mechanism. |
| 5 | Do matrix-associated RNA changes alter cell behavior? | Measure relevant protein/secreted outputs and migration or invasion under matched conditions, accounting for growth differences. |

The strongest next mechanistic distinction is between **locus-specific chromatin redistribution**, **a broader cell-state change**, and **regulator-specific dependency**. Existing RNA signatures alone cannot choose between them. Proposed methods and their source references are detailed in the [epigenetic report](outputs/report/epigenetic_deep_dive.html), [ZMYND8 notes](provenance/ZMYND8_interpretation.md) and [HDAC plan](HDAC12_PLAN.md).

## Suggested manuscript wording

> Repeated sampling of established HT1080 LOY and sgCTRL populations identified recurring differences in relative RNA abundance, including lower interferon/inflammatory and extracellular-matrix-associated programs. Chromatin-associated enrichments included reduced expression of HDAC2-depletion-downregulated and EED-depletion-induced genes, together with increased expression of distinct genes annotated with H3K27me3 at promoters in external source models. Gene-level and source-perturbation analyses indicated selective, context-dependent programs rather than a uniform change in chromatin repression. These observations nominate regulatory mechanisms for testing but do not establish altered enzyme activity or a general causal effect of LOY, because the sampled passages derive from the same established populations.

## Figure and evidence guide

| Purpose | Existing deliverable |
| --- | --- |
| Design, period effects and expression recurrence | [Figure 1](outputs/figures/Fig01_design_QC.pdf), [Figure 2](outputs/figures/Fig02_expression_repeatability.pdf), [main report](outputs/report/analysis.html) |
| Broad pathways and Y/paralog context | [Figure 3](outputs/figures/Fig03_pathways_regulators.pdf), [Figure 4](outputs/figures/Fig04_Y_paralog_dosage.pdf) |
| ZMYND8 and HDAC drivers/counterevidence | [Figure 6](outputs/figures/Fig06_ZMYND8_followup.pdf), [Figure 8](outputs/figures/Fig08_HDAC12_gene_drivers.pdf), [Figure 9](outputs/figures/Fig09_HDAC_external_comparison.pdf) |
| PRC2 genes and epigenetic context | [Figure 10](outputs/figures/Fig10_PRC2_gene_programs.pdf), [Figure 11](outputs/figures/Fig11_epigenetic_context.pdf) |
| Publication-quality running GSEA curves | [PRC2/H3K27me3, Figure 12](outputs/figures/Fig12_GSEA_PRC2_H3K27me3.pdf), [other chromatin signatures, Figure 13](outputs/figures/Fig13_GSEA_other_chromatin.pdf), [individual panels and report](outputs/report/chromatin_GSEA_curves.html) |
| Captions and panel source data | [Figure captions](outputs/captions/figure_captions.md), [figure manifest](outputs/figures/figure_manifest.csv) |

All quantitative statements above refer to the saved analysis outputs available on the update date. The underlying counts, expression models, enrichment results and figures were not recomputed for this synthesis. Earlier extension notes describe their completion-time state; the current interpretation includes the completed HDAC external benchmark and subsequent PRC2 work.
