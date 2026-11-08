# Chapter 3 — Figure captions

Latest review update (30 October 2026): current draftv2 Fig 3.11 (A549 global transcriptional variation) is now Supp Fig 3.44. Later main numbers remain stable during review. Panels and caption wording are unchanged.

Current page-reduction update (30 October 2026): former draftv2 main Figs 3.10, 3.11, 3.17 and 3.20 are now Supp Figs 3.40–3.43. Current draftv2 inline captions and its figure-numbering register govern export; other numbers in this historical master may reflect earlier sequences. Panels and caption wording are retained.

This file is the active source for Chapter 3 figure captions. Figure numbers remain provisional; headings use descriptive titles and identify main or supplementary figures.

## supp - Draft H3K27 and H3K4 signatures across all available comparisons

File: `figures/chapter3/ch03_supp_draft_h3k27_h3k4_all_comparisons.pdf`

*Exploratory four-page review draft, 28 September 2026; not selected for the thesis.*

**H3K27- and H3K4-related expression signatures across 46 comparisons.** Pages show H3K27 marks/functions, H3K4 marks/functions, signatures containing both marks, and PRC2/H3K4-regulator signatures. Add-back comparisons are excluded. Of the original 87 sets, 64 are displayed because they were tested in at least 23 of the remaining 46 comparisons. This display filter does not change the original testing families or FDR values. Rows are hierarchically clustered within each page using average linkage and Pearson correlation distance (1 - r), calculated from uncapped NES values across jointly tested comparisons, with optimal leaf ordering. Untested cells are excluded from each pairwise distance; column order is unchanged. Each displayed comparison has equal weight, so the clustering is descriptive and does not account for sample dependence. GSEA colours encode NES with limits of -4 to +4; complete uncapped results remain in the source table. Positive values refer to the higher end of the displayed comparison, including p53 induction versus uninduced controls. Dots identify FDR <0.05. Grey crosses denote sets absent from saved results or outside the eligible size range, rather than nonsignificant findings. Thesis models retain their saved adjustment across all tested sets within contrast. The external, RPTEC and autosomal/p53 comparisons use new exploratory fgsea tests with BH adjustment across eligible sets in the original focused panel within contrast; significance is therefore not directly equivalent across these blocks. New tests use 10–500 mapped genes, or 15–500 for early RPTEC. Existing comparison-specific gene rankings, loss-region exclusions and sample dependence are retained. HT1080 periods, RPTEC days, shared clones and reused controls are not independent experiments. HCP, ICP and LCP denote high-, intermediate- and low-CpG promoter sets. Gene-set enrichment describes transcription associated with the named reference signatures and does not measure histone modifications or establish regulator activity. The H3K27 group includes an acetylation signature; mixed-mark and perturbation signatures are displayed separately. ARPE-19 is pending and is not included.

## supp - Supp Fig 3.42. Transcriptional responses vary across chromosome-loss models and sampling times.

<a id="fig-3-28"></a>

**Figure:** `ch03_main_arm_loss_comparison` | **Results:** 3.3.4.5 | **Chapter number:** Supp Fig 3.42 (current draftv2; formerly main Fig 3.17)

[PDF](figures/chapter3/selected/page_reduction_20261030/ch03_supp_arm_loss_comparison.pdf) | [Editable SVG](figures/chapter3/selected/page_reduction_20261030/ch03_supp_arm_loss_comparison.svg) | [PNG](figures/chapter3/selected/page_reduction_20261030/ch03_supp_arm_loss_comparison_preview_841e717e428f.png)

**Supp Fig 3.42. Transcriptional responses vary across chromosome-loss models and sampling times.**

**a,b,** Hallmark (a) and GO:BP (b) NES from scRNAseq data (GSE336056) in renal proximal tubule epithelial cells (RPTECs) after CRISPR targeting of Yq or 3p. RNA-defined loss and intact groups were compared within each guide treatment, pooling days 5 and 10: 194 Yq-null versus 6,407 intact cells, and 649 3p-loss versus 8,934 remaining sg3p-treated cells. Colours/shapes indicate FDR < 0.05; ρ is Spearman correlation. **c,** Selected programmes across these and 14 bulk contrasts. Black dots mark significant enrichment (FDR < 0.05) for the corresponding programme and comparison. Brackets mark reused samples: [1], A2058 controls; [2], A2780 loss clones; [3], RPTEC clones/controls; [4], RPE1 p53-KO controls. WT/TK distinguish controls. **d,** Hallmark enrichment in the RPTEC CRISPR 3p-loss experiment at days 5, 10, 18 and 90 after targeting. Days 5 and 10 show scRNAseq comparisons of RNA-defined 3p-loss cells against remaining sg3p-treated cells (GSE336056). Days 18 and 90 show bulk RNAseq comparisons of the same three established 3p-loss clones against unedited controls (GSE336576). The divider separates the single-cell and bulk experiments. In c,d, blue/red indicate negative/positive NES and grey untested sets; asterisks in d mark FDR < 0.05. Correction is within collection/contrast. Rankings exclude lost regions, using MAST effects (single-cell) or DESeq2 Wald statistics (bulk). RNA calls do not independently establish DNA loss. [Dataset and analysis details](tables/chapter3/arm_loss_comparison/README.md).

[Source and reproduction notes](figures/chapter3/ch03_main_arm_loss_comparison_notes.md)

## Retired supp - RNA-based classification and regional expression in chromosome-loss models

**Status:** Removed from active Results on 8 October 2026 during the author-requested simplification. Historical exports and caption are retained; regional-expression checks remain in the linked source and design tables.

File: `figures/chapter3/ch03_supp_arm_loss_rna_checks.pdf`

*Draft caption for review, 28 September 2026.*

**Regional RNA depletion in bulk chromosome-loss models.** Median DESeq2 gene log2 fold changes within and outside the designated lost region for 14 bulk contrasts. Lines connect medians and do not indicate uncertainty. Bracketed numbers identify reused samples as in the main figure. ReDACT comparisons remove extra arm copies, whereas RPE1 comparisons involve monosomy. These regional RNA checks do not independently establish DNA copy number or precise deletion breakpoints.


## main - Complete LOY prevalence across TCGA cancer types

**Figure:** `ch03_main_intro_tcga_loy_prevalence` | **Introduction:** 3.2.1 | **Panel:** a (working designation)

[PDF](figures/chapter3/ch03_main_intro_tcga_loy_prevalence.pdf) | [Editable SVG](figures/chapter3/ch03_main_intro_tcga_loy_prevalence.svg) | [PNG](figures/chapter3/ch03_main_intro_tcga_loy_prevalence.png)

**a,** Percentage of male primary tumours with complete genomic loss of chromosome Y in each of the 28 cancer types with male prevalence estimates reported in Table S1C of Qi et al. (2023). Cancer types are ordered by decreasing prevalence. Bar-end labels give percentages rounded to one decimal place; the right-hand column gives the number of complete-LOY tumours divided by the number with a genomic Y-status call. Complete LOY corresponds to the published `Y_status = LOY` category and excludes relative LOY, partial LOY and expression-defined functional LOY. The percentages and denominators were independently checked against tumour-level calls in Table S1A. Data replotted from [Qi et al., Cell (2023), Fig. 2A and Tables S1A/S1C](https://doi.org/10.1016/j.cell.2023.06.006).

The overall cohort contained 1,390 complete-LOY tumours among 5,014 male tumours with genomic Y-status calls (27.7%, reported as 28% in the article). The 28 plotted cancer types contain 5,003 of these tumours. Table S1C does not provide a male breast-cancer prevalence estimate; Table S1A contains 11 called male breast tumours, none with complete LOY, which are included only in the overall denominator. Twenty-seven additional male records in Table S1A lack a genomic Y-status call and do not contribute to these denominators. These are descriptive frequencies in the published TCGA cohort; no confidence intervals or inferential tests were added.

[Source and reproduction notes](figures/chapter3/ch03_main_intro_tcga_loy_prevalence_notes.md)


## main - Genomic Y-chromosome status across TCGA cancer types

**Figure:** `ch03_main_intro_tcga_y_status_distribution` | **Introduction:** 3.2.1 | **Panel:** a (earlier standalone version)

**Status:** Earlier standalone version, superseded for thesis placement by the selected paired prevalence/hazard figure. Retained for reference.

[PDF](figures/chapter3/ch03_main_intro_tcga_y_status_distribution.pdf) | [Editable SVG](figures/chapter3/ch03_main_intro_tcga_y_status_distribution.svg) | [PNG](figures/chapter3/ch03_main_intro_tcga_y_status_distribution.png)

**a,** Distribution of the six published genomic Y-chromosome status categories among male primary tumours: complete LOY, relative LOY (rLOY), partial LOY (pLOY), gain, gain plus loss and WT. Statuses were taken directly from the mutually exclusive `Y_status` calls in Table S1A of [Qi et al., Cell (2023)](https://doi.org/10.1016/j.cell.2023.06.006). No calls were reclassified, and expression-defined functional LOY was not used. The gain-plus-loss label preserves the source category `Gain_Loss` and does not establish the temporal sequence of the alterations.

Each bar sums to 100%; the right-hand column gives its number of tumours with a genomic Y-status call. The 28 cancer types and their order match the complete-LOY-only panel, ranked by decreasing complete-LOY prevalence. Segment labels show rounded percentages for categories occupying at least 8% of a bar; smaller segments are retained without labels. Colours are consistent across cancer types, and gain plus loss is additionally marked by diagonal hatching. Exact counts and unrounded percentages for every category are retained in the source table.

The separate all-tumour bar includes all 5,014 called male tumours: 1,390 complete LOY (27.7%), 114 rLOY (2.3%), 206 pLOY (4.1%), 239 gain (4.8%), 77 gain plus loss (1.5%) and 2,988 WT (59.6%). The cancer-type bars contain 5,003 tumours, following the Table S1C male cohort breakdown used in the first version. Eleven male breast tumours are included in the all-tumour summary but not shown as a separate cancer-type bar; their calls comprise eight WT, one gain, one partial LOY and one gain plus loss. Twenty-seven male records without genomic Y-status calls are excluded. These are descriptive frequencies; no inferential tests were performed.

[Source and reproduction notes](figures/chapter3/ch03_main_intro_tcga_y_status_distribution_notes.md)

## main - Complete LOY and progression-free survival across TCGA cancer types

**Figure:** `ch03_main_intro_tcga_loy_pfs_forest` | **Introduction:** 3.2.1 | **Panel:** b (earlier standalone version)

**Status:** Earlier standalone version, superseded for thesis placement by the selected paired prevalence/hazard figure. Retained for reference.

[PDF](figures/chapter3/ch03_main_intro_tcga_loy_pfs_forest.pdf) | [Editable SVG](figures/chapter3/ch03_main_intro_tcga_loy_pfs_forest.svg) | [PNG](figures/chapter3/ch03_main_intro_tcga_loy_pfs_forest.png)

**b,** Published hazard ratios for progression-free survival (PFS), comparing male tumours with complete LOY against WT tumours. Cancer-specific estimates are from Table S4 of [Qi et al., Cell (2023)](https://doi.org/10.1016/j.cell.2023.06.006) and are adjusted for age at diagnosis and the number of autosomal chromosome arms lost. Points indicate hazard ratios and horizontal lines nominal 95% confidence intervals; the dashed reference line marks HR = 1. Hazard ratios above 1 indicate a higher event hazard in the complete-LOY group. Cancer types follow the order in panel a. LAML and PCPG are retained as labelled rows without estimates because neither is reported in Table S4.

Printed q values reproduce the published Benjamini-Hochberg adjustment across 26 cancer-specific LOY coefficient tests. Bold labels and vermilion q values identify q < 0.05. Table S4 provides full-precision HRs and coefficient Wald P values but not confidence limits. The 95% Wald intervals were recovered algebraically from these values, using the authors' published notebook to verify the coefficient test; no survival models were refitted. These intervals are not adjusted for multiple testing.

The grey diamond in the separate top section shows the published **unadjusted** pan-cancer estimate (HR 1.17, 95% CI 1.06-1.28; 4,271 observations and 1,956 events in the authors' saved notebook output). This pooled association is influenced by cancer-type composition and is not a pooled estimate of the adjusted cancer-specific models. Its HR and interval are reproduced at the two-decimal precision available in that output. The source uses the PFS endpoint in the TCGA-CDR ExtraEndpoints sheet; PFS and the resource's PFI endpoint are distinct, and the published estimates are not relabelled as PFI. The panel depicts associations rather than causal effects.

[Source and reproduction notes](figures/chapter3/ch03_main_intro_tcga_loy_pfs_forest_notes.md)

## main - Fig 3.1. LOY prevalence and progression-free survival across TCGA cancer types.

**Figure:** `ch03_main_intro_tcga_prevalence_hazard_paired` | **Introduction:** draftv1, section 3.1 | **Chapter number:** Fig 3.1 | **Status:** Selected by the author on 17 September 2026; inserted into the shorter Introduction at the author's request on 7 October 2026. This combined figure supplies Introduction panels a and b; the separate panels remain earlier alternatives.

[PDF](figures/chapter3/ch03_main_intro_tcga_prevalence_hazard_paired.pdf) | [Editable SVG](figures/chapter3/ch03_main_intro_tcga_prevalence_hazard_paired.svg) | [PNG](figures/chapter3/ch03_main_intro_tcga_prevalence_hazard_paired.png) | [Interactive version](figures/chapter3/ch03_main_intro_tcga_prevalence_hazard_paired.html)

**Fig 3.1. LOY prevalence and progression-free survival across TCGA cancer types.**

**a,** Genomic Y status in male TCGA tumours, ordered by complete-LOY prevalence. Bars show complete LOY (no Y copies), relative LOY (rLOY; fewer whole-Y copies than expected for tumour ploidy, but at least one retained), partial LOY (pLOY; regional Y loss), gain (Y copy gains without losses), gain+loss (both regional gains and losses), and WT (Y copy number expected for tumour ploidy); *n* denotes tumours with a call. **b,** Progression-free survival hazard ratios comparing complete LOY with unaltered Y status, adjusted within each cancer type for age and autosomal arm losses. Error bars show nominal 95% confidence intervals; *q* values use Benjamini–Hochberg correction. Curves connect cancer types; red curves and bold labels mark *q* < 0.05, and dashed curves mark unavailable estimates. Bottom rows show the overall Y-status distribution (5,014 tumours) and unadjusted pan-cancer estimate (4,271 observations). Adapted from [Qi et al. (2023)](https://doi.org/10.1016/j.cell.2023.06.006), Tables S1A/S4 and the published analysis notebook.

[Source and reproduction notes](figures/chapter3/ch03_main_intro_tcga_prevalence_hazard_paired_notes.md)

## main - Fig 3.2. Complete Y loss and autosomal arm-loss burden across TCGA cancer types.

**Figure:** `ch03_main_intro_tcga_loy_arm_losses` | **Introduction:** draftv1, section 3.1 | **Chapter number:** Fig 3.2 | **Status:** Added at the author's request on 7 October 2026.

[PDF](figures/chapter3/ch03_main_intro_tcga_loy_arm_losses.pdf) | [Editable SVG](figures/chapter3/ch03_main_intro_tcga_loy_arm_losses.svg) | [PNG](figures/chapter3/ch03_main_intro_tcga_loy_arm_losses.png)

**Fig 3.2. Complete Y loss and autosomal arm-loss burden across TCGA cancer types.**

Complete-LOY prevalence versus mean autosomal chromosome-arm losses across 28 TCGA cancer types. LOY proportions use tumours with Y-status calls; arm-loss means use all available male records, including repeat records per case. Red points highlight papillary renal cell carcinoma (KIRP), clear cell renal cell carcinoma (KIRC) and uveal melanoma (UVM). Pearson *r* = 0.333; two-sided *P* = 0.0833; *n* = 28 cancer types. Adapted from [Qi et al. (2023)](https://doi.org/10.1016/j.cell.2023.06.006), Figure 2E and Tables S1A/S1C.

[Source data](figures/chapter3/source_data/intro_tcga_loy_arm_losses/plotted_values.csv) | [Source and reproduction notes](figures/chapter3/ch03_main_intro_tcga_loy_arm_losses_notes.md)

## main - LOY prevalence and chromosome-arm loss

**Figure:** `ch03_main_intro_loy_arm_aneuploidy` | **Introduction:** 3.2.1 | **Status:** Reprocessed candidate figure; numbering remains provisional.

[PDF](figures/chapter3/ch03_main_intro_loy_arm_aneuploidy.pdf) | [Editable SVG](figures/chapter3/ch03_main_intro_loy_arm_aneuploidy.svg) | [PNG](figures/chapter3/ch03_main_intro_loy_arm_aneuploidy.png)

**a,** Complete-LOY prevalence versus the mean number of autosomal chromosome arms lost across 28 TCGA cancer types, reprocessed from the sources used for Figure 2E of [Qi et al. (2023)](https://doi.org/10.1016/j.cell.2023.06.006). Complete-LOY percentages use the published genomic Y-status calls and denominators in Tables S1A/S1C. Arm-loss means use Taylor et al. (2018) Table S2 `AS_del`, with male samples identified through TCGA-CDR. Following the authors' notebook, means include all available male arm-score records for each cancer type and are not restricted to samples with a Qi Y-status call. The two axes therefore need not share sample denominators. Multiple sample records per case are retained as in the source analysis. The plotted cancer types contribute 5,028 arm-score records; per-type record and unique-case counts are provided in the source table. No missing means were imputed. Pearson r = 0.333 and two-sided P = 0.0833 reproduce the published notebook result; n = 28 cancer types, not individual tumours. This measure counts losses only, rather than combined arm gains and losses.

**b,** Yp deletion frequency versus the published GISTIC Yp deletion z-score for the same 28 cancer types, replotted directly from Qi Table S2A (Figure 2F). These scores derive from the paper's arm-level background analysis; no GISTIC model was rerun. Yp deletion frequency is an arm-level measure and is not interchangeable with the complete-LOY prevalence in a. KIRP (papillary renal cell carcinoma), KIRC (clear cell renal cell carcinoma) and UVM (uveal melanoma) are highlighted in both panels because they are discussed in the Introduction, rather than according to a new significance test.

Each point represents one cancer type, labelled by its TCGA code. The original notebooks' robust regression curves are omitted. The descriptive cross-cancer correlation in a is distinct from the background-model analysis in b and does not by itself establish selection for LOY. The arm-by-arm gene-count plots in Qi Figure 2G were not reconstructed. No new test of selection or causal effect was performed.

[Source and reproduction notes](figures/chapter3/ch03_main_intro_loy_arm_aneuploidy_notes.md)

## main - Fig 3.6. Y-chromosome guide filtering and the top-hits table.

<a id="fig-3-6"></a>

**Figure:** `ch03_main_guide_design_edit` · **Working ID:** CH3-S3.4.4-F01 · **Results:** 3.4.4 | **Chapter number:** Fig 3.6

[Selected PDF](figures/chapter3/selected/CH3_S3_4_4_F01_guide_discovery_selection/CH3_S3_4_4_F01_guide_discovery_selection.pdf) · [Selected PNG](figures/chapter3/selected/CH3_S3_4_4_F01_guide_discovery_selection/CH3_S3_4_4_F01_guide_discovery_selection.png)

[PDF](figures/chapter3/ch03_main_guide_design_edit.pdf) · [Editable SVG](figures/chapter3/ch03_main_guide_design_edit.svg) · [PNG](figures/chapter3/ch03_main_guide_design_edit.png)

**Fig 3.6. Y-chromosome guide filtering and the top-hits table.**

**a,** The complete HG002 T2T-Y chromosome sequence (62,460,029 bp; CP086569.2) was scanned using a sliding 23-nt window. Sequences comprising a 20-nt candidate protospacer followed by an NGG PAM were retained, considering both strands. Each tile represents a nucleotide position.

**b,** The candidate table contained 14,723 distinct 20-nt guides. Requiring at least 10 exact Y matches and no exact non-Y matches, combining both orientations, retained 6,901 guides. The shortlist contained 215 guides with no recorded exact non-Y matches or non-Y sites with one or two mismatches.
**c,** All five initially selected guides (Q53, CN, CP, P47 and P51) were represented in the top-hits table. Vertical dividers separate Yq, centromeric and Yp targets. The 20-nt guide sequences are shown below their labels in the 5′-to-3′ direction, excluding the PAM.

[Source and reproduction notes](figures/chapter3/ch03_main_guide_design_edit_caption.md#source-and-reproduction-notes)

## main - Fig 3.7. Five Y-targeting guides and their repeat regions.

<a id="fig-3-7"></a>

**Figure:** `ch03_main_guide_target_map_stacked` · **Working ID:** CH3-S3.4.4-F02 · **Results:** 3.4.4 | **Chapter number:** Fig 3.7

[Selected PDF](figures/chapter3/selected/CH3_S3_4_4_F02_guide_target_map/CH3_S3_4_4_F02_guide_target_map.pdf) · [Selected PNG](figures/chapter3/selected/CH3_S3_4_4_F02_guide_target_map/CH3_S3_4_4_F02_guide_target_map.png)

[PDF](figures/chapter3/ch03_main_guide_target_map_stacked.pdf) · [Editable SVG](figures/chapter3/ch03_main_guide_target_map_stacked.svg) · [PNG](figures/chapter3/ch03_main_guide_target_map_stacked.png)

**Fig 3.7. Five Y-targeting guides and their repeat regions.**

**a,** Exact matches to five guides and annotated repeat regions across HG002 T2T-Y (CP086569.2). Q53 targets DYZ1/HSat3 A6 repeats in Yq12 (blue; 532 matches); CN/CP target centromeric DYZ3 α-satellite (red; 53/54); P47/P51 target TSPY-associated Yp repeats (teal; 47/51). Tracks show all intervals in each repeat class, including those without guide matches. **b,c,** Enlargements of DYZ3 (b) and TSPY-associated repeats (c), retaining all matches for the displayed guides. Points mark exact 20-nucleotide matches with correctly oriented NGG PAMs; vertical offsets separate overlaps. Sequences are 5′ to 3′ without PAMs; coordinates are megabases. Reference: [Rhie et al. (2023)](https://doi.org/10.1038/s41586-023-06457-y).

[Source and reproduction notes](figures/chapter3/ch03_main_guide_target_map_stacked_notes.md)

## supp - Supp Fig 3.7. ddPCR measurements across known cell mixtures.

<a id="supp-fig-3-7"></a>

**Figure:** `ch03_main_ddpcr_mixing_validation` | **Working ID:** CH3-S3.4.5.1-S03 | **Results:** 3.4.5.1 | **Chapter number:** Supp Fig 3.7

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_1_S03_ddpcr_mixing_check/CH3_S3_4_5_1_S03_ddpcr_mixing_check.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_1_S03_ddpcr_mixing_check/CH3_S3_4_5_1_S03_ddpcr_mixing_check.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_1_S03_ddpcr_mixing_check/CH3_S3_4_5_1_S03_ddpcr_mixing_check.png)

[PDF](figures/chapter3/ch03_main_ddpcr_mixing_validation.pdf) · [Editable SVG](figures/chapter3/ch03_main_ddpcr_mixing_validation.svg) · [PNG](figures/chapter3/ch03_main_ddpcr_mixing_validation.png)

**Supp Fig 3.7. ddPCR measurements across known cell mixtures.**

Y-negative KP4 and Y-positive H2452 cells were mixed before DNA extraction at expected H2452 fractions of 0, 0.25, 0.5, 0.75 and 1. **a,b,** Observed AMELY/AMELX ratios measured using the AMELX/Y ddPCR assay with 50 and 100 ng DNA input per reaction, respectively. Points show technical measurements; solid lines show linear fits and dashed lines indicate numerical equality. Pearson correlations and measurement counts are displayed.

[Source and reproduction notes](figures/chapter3/ch03_main_ddpcr_mixing_validation_notes.md)

## main — HT1080 origin, Y-linked expression and baseline dependencies

**Figure:** `ch03_main_ht1080_profile` · **Results:** 3.4.5

[PDF](figures/chapter3/ch03_main_ht1080_profile.pdf) · [Editable SVG](figures/chapter3/ch03_main_ht1080_profile.svg) · [PNG](figures/chapter3/ch03_main_ht1080_profile.png)

**a,** Reported background of HT1080 (RRID:CVCL_0317), derived from a pelvic fibrosarcoma in a 35-year-old man. The approximately 30-h doubling time is a DSMZ culture estimate; the pseudodiploid karyotype (modal chromosome count 46) and genotype annotations are published or repository-reported properties. These annotations do not represent new measurements of the author's stocks. [DSMZ ACC 315](https://www.dsmz.de/collection/catalogue/details/culture/ACC-315); [ECACC 85111505](https://www.culturecollections.org.uk/nop/product/ht-1080); [ATCC CCL-121](https://www.atcc.org/products/ccl-121); [Cellosaurus](https://www.cellosaurus.org/CVCL_0317). TP53 wild type refers to reports from ATCC-derived cultures; TP53-mutant HT1080 cultures have also been described. [Moran and Maki, 2010](https://doi.org/10.1158/1535-7163.MCT-09-1220); [TP53-mutant culture report](https://pubmed.ncbi.nlm.nih.gov/7519049/).

**b,** Reference microscopy image of HT-1080 from the Leibniz Institute DSMZ (ACC 315), reproduced with source attribution. The original image, including its 50-µm scale bar and source labels, is shown in full. This is a repository culture image. [DSMZ image](https://www.dsmz.de/fileadmin/downloads/ACC/cells_in_culture/HT-1080_Bild_2.jpg).

**c,** Expression of nine Y-linked genes in the default DepMap Public 26Q1 RNA profile for HT1080 (ACH-000054), expressed as log₂(TPM + 1). Genes are ordered by their GRCh38 start coordinates from Yp to Yq: RPS4Y1, ZFY, USP9Y, DDX3Y, UTY, TMSB4Y, NLGN4Y, KDM5D and EIF1AY. Connectors link gene midpoints on the full-length reference chromosome diagram to equally spaced expression rows. The diagram uses Ensembl gene coordinates and UCSC hg38 cytobands; it describes reference positions rather than the karyotype of the profiled cells. Labels give expression values to two decimal places.

**d,** Sixteen relatively selective baseline genetic dependencies in HT1080. Blue points show integrated Chronos gene effects; grey marks show the median and interquartile range across the 1,207 other screened models, including both sexes and all lineages. More negative gene effects indicate stronger dependency. Genes were required to have an HT1080 gene effect below −0.5 and an effect below −0.5 in no more than 20% of other models. Of 103 qualifying genes, the 16 with the largest negative differences from the reference median are displayed. This is a descriptive ranking of baseline public-data profiles, not a test of LOY-associated synthetic lethality. Interquartile ranges describe variation between reference models; no inferential P values are reported.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_profile_notes.md)

## supp — Y-linked expression in the four selected cell lines

**Figure:** `ch03_supp_selected_models_y_expression` · **Results:** 3.4.5

[PDF](figures/chapter3/ch03_supp_selected_models_y_expression.pdf) · [Editable SVG](figures/chapter3/ch03_supp_selected_models_y_expression.svg) · [PNG](figures/chapter3/ch03_supp_selected_models_y_expression.png)

Expression of nine Y-linked genes in HT1080, A549, DU145 and NCI-H2452 from DepMap Public 26Q1. Each column represents the default RNA profile for one model; colour and printed values show log₂(TPM + 1), with a common colour scale of 0–10. Values are rounded to two decimal places for display. The genes include the Y paralogs of DDX3X, EIF1AX, RPS4X and ZFX, together with five additional Y-linked genes. PRKY was absent from the protein-coding expression matrix and was omitted without imputation. These baseline RNA profiles describe gene expression and do not independently establish chromosome-Y copy number or complete LOY.

Rows follow GRCh38 genomic order from Yp to Yq, matching panel c of the HT1080 profile. A full-length chromosome diagram with UCSC cytobands is linked to the rows through the midpoint of each Ensembl gene interval. Gene positions are shown on the reference chromosome; expression rows are equally spaced for readability.

[Source and reproduction notes](figures/chapter3/ch03_supp_selected_models_y_expression_notes.md)


## main - Initial assessment of five Y-targeting guides by RT-qPCR in HT1080

**Figure:** `ch03_main_ht1080_initial_guides_rtqpcr` | **Results:** 3.4.5

**Status:** Alternative draft retained for individual-well and run-level inspection; not the selected main figure.

[PDF](figures/chapter3/ch03_main_ht1080_initial_guides_rtqpcr.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_initial_guides_rtqpcr.svg) | [PNG](figures/chapter3/ch03_main_ht1080_initial_guides_rtqpcr.png)

Relative expression of **a,** EIF1AY; **b,** KDM5D; **c,** RPS4Y1; and **d,** UTY in HT1080 WT and populations treated individually with CN, CP, Q53, P47 or P51.

Circles and triangles denote technical wells from Run 1 and Run 2, respectively (source files labelled ATTEMPT1 and ATTEMPT2). The runs are repeat assays of the same edited populations, not independent transfections. Each available target/population/run combination contains three technical wells; CN was available only in Run 1, giving 132 plotted wells across the four targets. Horizontal offsets separate runs and wells and have no quantitative meaning. Short horizontal marks show expression ratios calculated from the reported target mean Ct; they are not error bars. Dashed lines indicate WT = 1.

Expression was calculated as 2^(-delta-delta-Ct), normalising each target to the reported mean GAPDH Ct of the same population and then to the target-minus-GAPDH mean Ct of WT within the same run. Individual points retain target-well variation; reference-well uncertainty is not propagated. Runs were not pooled, and no inferential tests were performed. These measurements describe transcript abundance and do not establish the fraction of cells with chromosome loss. RNA-collection times and whether the repeat assays used the same RNA/cDNA preparation remain to be confirmed.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_initial_guides_rtqpcr_notes.md)

## supp - Additional transcripts measured in the initial HT1080 RT-qPCR assay

**Figure:** `ch03_supp_ht1080_initial_rtqpcr_additional_targets` | **Results:** 3.4.5

[PDF](figures/chapter3/ch03_supp_ht1080_initial_rtqpcr_additional_targets.pdf) | [Editable SVG](figures/chapter3/ch03_supp_ht1080_initial_rtqpcr_additional_targets.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_initial_rtqpcr_additional_targets.png)

**a,** USP9Y and **b,** ZFY expression in the six HT1080 populations in Run 1. Points represent three technical wells per target/population combination (36 wells total); short horizontal marks show ratios calculated from the reported mean Ct. Normalisation and the WT = 1 reference follow the main figure. These targets were absent from the Run 2 table. No uncertainty intervals or inferential tests are shown. NLGN4Y was also assayed in Run 1 but could not be normalised because all three WT wells were undetermined; one P51 well had a recorded Ct and all other NLGN4Y wells were undetermined. Undetermined values were not replaced with zero expression or an arbitrary Ct.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_initial_guides_rtqpcr_notes.md)

## main — Selected cell lines in Y-expression-defined dependency distributions

**Figure:** `ch03_main_selected_models_paralog_dependencies` · **Results:** 3.4.5

[PDF](figures/chapter3/ch03_main_selected_models_paralog_dependencies.pdf) · [Editable SVG](figures/chapter3/ch03_main_selected_models_paralog_dependencies.svg) · [PNG](figures/chapter3/ch03_main_selected_models_paralog_dependencies.png)

Integrated Chronos gene effects for **a,** DDX3X; **b,** EIF1AX; **c,** RPS4X; and **d,** ZFX in DepMap Public 26Q1. The reference distributions contain 352 Y-expressing and 219 Y-depleted male-derived models, retaining one donor representative with an eligible integrated CRISPR profile under the executed RNA-based analysis. Intermediate and unassessable RNA states are excluded. Blue and vermilion denote Y-expressing and Y-depleted groups, respectively.

Violins show kernel-density distributions; boxes indicate the interquartile range, white lines the median, and whiskers the most extreme observations within 1.5 interquartile ranges of the box. HT1080, A549, DU145 and NCI-H2452 are displayed separately on the same gene-effect axis, with their scores printed to two decimal places. These four models are also included in the Y-expressing reference distribution. More negative scores indicate stronger dependency; axis ranges differ between genes.

Displayed mean differences are Y-depleted minus Y-expressing. The reference pools the analysis's discovery and held-out partitions for descriptive comparison, without covariate adjustment or a new inferential test. The group labels describe Y-linked RNA expression and do not independently establish a diploid XY karyotype or physical chromosome loss. The accompanying DNA-supported comparison retains genomic classification and assay-matching distinctions.

[Source and reproduction notes](figures/chapter3/ch03_main_selected_models_paralog_dependencies_notes.md)

## supp — Selected cell lines against the DNA-supported LOY cohort

**Figure:** `ch03_supp_selected_models_dna_paralog_dependencies` · **Results:** 3.4.5

[PDF](figures/chapter3/ch03_supp_selected_models_dna_paralog_dependencies.pdf) · [Editable SVG](figures/chapter3/ch03_supp_selected_models_dna_paralog_dependencies.svg) · [PNG](figures/chapter3/ch03_supp_selected_models_dna_paralog_dependencies.png)

DepMap Public 26Q1 Chronos gene effects for **a,** DDX3X; **b,** EIF1AX; **c,** RPS4X; and **d,** ZFX, using the 151-model primary cohort from the executed DNA-supported dependency analysis: 81 Y-retaining and 70 LOY donor representatives. Violin, box and whisker definitions follow the main comparison. HT1080, A549, DU145 and NCI-H2452 are shown as separately labelled points with their integrated gene-effect scores.

HT1080 belongs to this reference cohort. Grey points mark models outside it: A549 has low-confidence, DNA–RNA-discordant Y evidence; DU145 fails the primary DNA/CRISPR culture-condition matching criterion; and H2452 has unresolved Y-dosage evidence. Their points do not contribute to the reference distributions or the reported statistical tests. All four retain Y-linked expression in the complementary RNA analysis.

Adjusted differences and FDR-adjusted P values are from the saved primary model, which conditions on the recorded lineage/screen strata. Differences are LOY minus Y-retaining. P values derive from two-sided moderated t tests and Benjamini–Hochberg correction across the 17,782 genes assessed in that scan. Group sizes and mean gene effects were matched to the original results before these statistics were added. More negative effects indicate stronger dependency, and axis ranges differ between genes.

[Source and reproduction notes](figures/chapter3/ch03_supp_selected_models_dna_paralog_dependencies_notes.md)


## main - Mean expression reduction across the Y chromosome after individual-guide treatment

**Figure:** `ch03_main_ht1080_rtqpcr_reduction_ideogram` | **Results:** 3.4.5

**Status:** Alternative four-gene draft retained; superseded as the main figure by the selected three-gene version.

[PDF](figures/chapter3/ch03_main_ht1080_rtqpcr_reduction_ideogram.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_rtqpcr_reduction_ideogram.svg) | [PNG](figures/chapter3/ch03_main_ht1080_rtqpcr_reduction_ideogram.png)

Mean relative-expression reductions for RPS4Y1, UTY, KDM5D and EIF1AY after individual treatment of HT1080 with CN, CP, Q53, P47 or P51. The GRCh38 Y-chromosome ideogram is drawn to scale; leader lines connect annotated gene midpoints to equally spaced columns. Gene coordinates are from Ensembl and cytobands from UCSC, using the frozen reference tables retained with the figure data. The ideogram provides reference context and does not represent an observed structural state of the edited populations.

For each gene and guide, expression ratios calculated from reported mean Ct values were first normalised to GAPDH and WT within each run. The arithmetic mean of the two run-specific ratios was then used to calculate reduction as 100 x (1 - mean relative expression). CN uses its single available run. Each run comprised three technical wells per measured target/population combination; runs were repeat assays of the same edited populations, not independent transfections. Printed percentages are rounded to the nearest integer and colour encodes unrounded values on a scale centred at zero. Positive values indicate expression below WT; negative values indicate expression above WT and were not truncated. CD99 is omitted at the author's request.

This descriptive summary does not estimate the percentage of cells that lost chromosome Y and includes no inferential tests. Run-to-run and well-level variation is shown in the companion individual-measurement figure.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_initial_guides_rtqpcr_notes.md)

## main - Fig 3.8. Guide selection and Y-linked DNA depletion in HT1080 and A549.

<a id="fig-3-8"></a>

**Figure:** `ch03_main_ht1080_guide_screening_ddpcr` | **Working ID:** CH3-S3.4.5.1-F01 | **Results:** HT1080 and A549 model development | **Status:** Five-panel combination adopted 27 October 2026 | **Chapter number:** Fig 3.8

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_1_F01_ht1080_guide_screening_ddpcr/CH3_S3_4_5_1_F01_ht1080_guide_screening_ddpcr.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_1_F01_ht1080_guide_screening_ddpcr/CH3_S3_4_5_1_F01_ht1080_guide_screening_ddpcr.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_1_F01_ht1080_guide_screening_ddpcr/CH3_S3_4_5_1_F01_ht1080_guide_screening_ddpcr.png)

[PDF](figures/chapter3/ch03_main_ht1080_guide_screening_ddpcr.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_guide_screening_ddpcr.svg) | [PNG](figures/chapter3/ch03_main_ht1080_guide_screening_ddpcr.png)

**Fig 3.8. Guide selection and Y-linked DNA depletion in HT1080 and A549.**

**a,** HT1080 RT–qPCR measurements of RPS4Y1, KDM5D and EIF1AY after individual guide delivery. Expression was normalised to GAPDH and sgCTRL. Values show percentage reductions relative to sgCTRL, averaged across two assay runs of the same edited populations, except CN (one run), with three technical wells per target and population per run. Negative values indicate expression above sgCTRL. The GRCh38 Y ideogram marks the RNA targets and ddPCR target genes AMELY (Yp; b,d) and DDX3Y (Yq; c,e).

**b–e,** ddPCR estimates of AMELY (b,d) and DDX3Y (c,e) copies per cell in HT1080 (b,c) and A549 (d,e). Bars show medians, whiskers interquartile ranges, and points individual measurements; n includes repeat assays. Labels 1× and 2× denote one and two successive nucleofections; hatching identifies the latter. The shape key applies to HT1080 only; A549 circles represent assay measurements. Estimates assume AMELX/RNaseP reference-copy numbers of 1/2 in HT1080 and 2/3 in A549, unchanged after editing. Dashed lines mark one target copy per cell; NA denotes the unmeasured HT1080 DDX3Y result after two nucleofections.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_guide_screening_ddpcr_notes.md)

## Retired - Complete HT1080 RT-qPCR guide comparison

CH3-S3.4.5.1-S02 was removed from the active supplementary figures by the author. The main guide comparison remains CH3-S3.4.5.1-F01. Historical caption and exports are preserved in `archive/retire_ht1080_all_rtqpcr_20260923_232334`.


## supp - Supp Fig 3.6. Reference genomic and expression profiles of HT1080.

<a id="supp-fig-3-6"></a>

**Figure:** `ch03_main_ht1080_profile_wgs_cn` | **Working ID:** CH3-S3.4.5.1-S01 | **Results:** 3.4.5.1 | **Title:** Reference genomic and expression profiles of HT1080. | **Chapter number:** Supp Fig 3.6

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_1_S01_ht1080_profile/CH3_S3_4_5_1_S01_ht1080_profile.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_1_S01_ht1080_profile/CH3_S3_4_5_1_S01_ht1080_profile.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_1_S01_ht1080_profile/CH3_S3_4_5_1_S01_ht1080_profile.png)

[PDF](figures/chapter3/ch03_main_ht1080_profile_wgs_cn.pdf) · [Editable SVG](figures/chapter3/ch03_main_ht1080_profile_wgs_cn.svg) · [PNG](figures/chapter3/ch03_main_ht1080_profile_wgs_cn.png)

**Supp Fig 3.6. Reference genomic and expression profiles of HT1080.**

**a,b,** Reported cell-line characteristics and reference culture image of HT1080, a male-derived fibrosarcoma line. Image from DSMZ (ACC 315); scale bar, 50 µm.

**c,** Y-linked expression and gene-level WGS copy-number estimates from DepMap Public 26Q1, ordered by GRCh38 position. Bars show log₂(TPM + 1).

**d,** Top selective dependencies in HT1080. Blue points show Chronos gene effects; grey marks show the median and interquartile range across 1,207 other models. Genes had effects below −0.5 in HT1080 and no more than 20% of reference models, ranked by their difference from the reference median. More negative values indicate stronger dependency.

**e,** Genome-wide WGS copy-number estimates in GRCh38 order. Horizontal marks show segments; points mark their midpoints. Blue identifies chromosome Y, and the dashed line marks two copies.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_profile_wgs_cn_notes.md)

## main - HT1080 CP2 RT-qPCR summary

**Figure:** `ch03_main_ht1080_cp2_rtqpcr_three_genes_ideogram` | **Results:** 3.4.5

**Status:** Processed candidate figure; placement in the narrative remains to be selected.

[PDF](figures/chapter3/ch03_main_ht1080_cp2_rtqpcr_three_genes_ideogram.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_cp2_rtqpcr_three_genes_ideogram.svg) | [PNG](figures/chapter3/ch03_main_ht1080_cp2_rtqpcr_three_genes_ideogram.png)

Mean reductions in RPS4Y1, KDM5D and EIF1AY expression in HT1080 CP2 relative to run-matched HT1080 WT. The two TRANSFECTION2 files are treated as repeat assays of the same CP2 population, with CP2 interpreted as a second CP transfection on the basis of the author's working interpretation. Each run contained three technical wells per target/population combination. The figure uses an equal-weight arithmetic mean of the two separately normalised expression ratios before calculating reduction. The source label E1F1AY is harmonised to EIF1AY; the original label is retained in the derived table.

Expression was calculated as 2^(-delta-delta-Ct), using the reported mean GAPDH Ct for each sample and the matched reference sample within each run. Reduction is 100 x (1 - relative expression); positive values indicate lower expression and negative values indicate higher expression than the reference. No undetermined Ct was replaced with 40. Numerical labels are rounded, with values above 99.9% displayed as >99.9%; this is a display convention for calculated estimates, not a validated detection limit. Colours use the same -100 to 100% scale as the initial HT1080 figure. The GRCh38 Y ideogram supplies reference gene positions, not an observed chromosome structure. No inferential tests were performed.

[Source and reproduction notes](figures/chapter3/other_rtqpcr_processing_notes.md)

## supp - All targets in the HT1080 CP2 assays

**Figure:** `ch03_supp_ht1080_cp2_rtqpcr_all_genes_ideogram` | **Results:** 3.4.5

**Status:** Processed candidate figure; placement in the narrative remains to be selected.

[PDF](figures/chapter3/ch03_supp_ht1080_cp2_rtqpcr_all_genes_ideogram.pdf) | [Editable SVG](figures/chapter3/ch03_supp_ht1080_cp2_rtqpcr_all_genes_ideogram.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_cp2_rtqpcr_all_genes_ideogram.png)

All six recorded targets are shown. CP2 uses the arithmetic mean of two separately normalised assay runs; A375 was available only in Run 1 and KP4 only in Run 2. These additional cell lines are contextual comparators, not matched edited controls. Three technical wells were recorded per target/sample/run combination. E1F1AY is labelled EIF1AY; UTY2 retains its source assay label and is positioned at the UTY locus for reference context, pending primer-record confirmation. Asterisks indicate ratios calculated from a partially detected target: A375 KDM5D and EIF1AY were detected in one of three wells, and KP4 RPS4Y1 in two of three wells. ND indicates that all KP4 KDM5D target wells were undetermined, with no numerical reduction assigned. GAPDH is marked ref.; CD99 is an X/Y target. Partial-detection ratios use the reported mean of determined wells and are not fully quantified sample estimates.

Expression was calculated as 2^(-delta-delta-Ct), using the reported mean GAPDH Ct for each sample and the matched reference sample within each run. Reduction is 100 x (1 - relative expression); positive values indicate lower expression and negative values indicate higher expression than the reference. No undetermined Ct was replaced with 40. Numerical labels are rounded, with values above 99.9% displayed as >99.9%; this is a display convention for calculated estimates, not a validated detection limit. Colours use the same -100 to 100% scale as the initial HT1080 figure. The GRCh38 Y ideogram supplies reference gene positions, not an observed chromosome structure. No inferential tests were performed.

[Source and reproduction notes](figures/chapter3/other_rtqpcr_processing_notes.md)

## main - HT1080 C24 clone expression relative to Y-positive C7

**Figure:** `ch03_main_ht1080_clones_rtqpcr_three_genes_ideogram` | **Results:** 3.4.5

**Status:** Processed candidate figure; placement in the narrative remains to be selected.

[PDF](figures/chapter3/ch03_main_ht1080_clones_rtqpcr_three_genes_ideogram.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_clones_rtqpcr_three_genes_ideogram.svg) | [PNG](figures/chapter3/ch03_main_ht1080_clones_rtqpcr_three_genes_ideogram.png)

RPS4Y1, KDM5D and EIF1AY expression in C24 and the separate KP4 comparator relative to C7, which the author confirmed as Y-positive. Four technical wells were recorded per target/sample combination. C24 RPS4Y1 was detected in all four wells (reported mean Ct 32.167), C24 EIF1AY in one of four wells (Ct 33.626), and C24 KDM5D in none. KP4 RPS4Y1 was detected in one of four wells, with EIF1AY and KDM5D undetermined in all wells. Asterisks mark detected-well estimates from partial detection; ND marks all-undetermined targets, not 100% reduction. The RPS4Y1 NTC was detected in two of four wells (reported mean Ct 33.046), and GAPDH NTC in three of four (mean Ct 30.180). These late NTC signals and partial target detection limit interpretation of the very low clone signals; amplification/melt-curve review and a validated quantification range remain outstanding. The figure does not establish complete chromosome loss.

Expression was calculated as 2^(-delta-delta-Ct), using the reported mean GAPDH Ct for each sample and the matched reference sample within each run. Reduction is 100 x (1 - relative expression); positive values indicate lower expression and negative values indicate higher expression than the reference. No undetermined Ct was replaced with 40. Numerical labels are rounded, with values above 99.9% displayed as >99.9%; this is a display convention for calculated estimates, not a validated detection limit. Colours use the same -100 to 100% scale as the initial HT1080 figure. The GRCh38 Y ideogram supplies reference gene positions, not an observed chromosome structure. No inferential tests were performed.

[Source and reproduction notes](figures/chapter3/other_rtqpcr_processing_notes.md)

## supp - All targets in the December 2024 HT1080 clone assay

**Figure:** `ch03_supp_ht1080_clones_rtqpcr_all_genes_ideogram` | **Results:** 3.4.5

**Status:** Processed candidate figure; placement in the narrative remains to be selected.

[PDF](figures/chapter3/ch03_supp_ht1080_clones_rtqpcr_all_genes_ideogram.pdf) | [Editable SVG](figures/chapter3/ch03_supp_ht1080_clones_rtqpcr_all_genes_ideogram.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_clones_rtqpcr_all_genes_ideogram.png)

All four measured targets are shown for C24 and KP4 relative to the author-confirmed Y-positive C7 reference. GAPDH is explicitly marked as the normalising reference rather than shown as an independent reduction endpoint. Four technical wells were recorded per target/sample combination. The partial-detection and NTC qualifications in the three-gene clone figure apply: C24 EIF1AY and KP4 RPS4Y1 each had one determined well, C24 and KP4 KDM5D had none, and KP4 EIF1AY had none. Asterisks identify ratios from the determined-well mean only; ND indicates all-undetermined targets. RPS4Y1 and GAPDH had late NTC amplification, so low detected signals do not provide validated quantitative evidence of complete Y loss.

Expression was calculated as 2^(-delta-delta-Ct), using the reported mean GAPDH Ct for each sample and the matched reference sample within each run. Reduction is 100 x (1 - relative expression); positive values indicate lower expression and negative values indicate higher expression than the reference. No undetermined Ct was replaced with 40. Numerical labels are rounded, with values above 99.9% displayed as >99.9%; this is a display convention for calculated estimates, not a validated detection limit. Colours use the same -100 to 100% scale as the initial HT1080 figure. The GRCh38 Y ideogram supplies reference gene positions, not an observed chromosome structure. No inferential tests were performed.

[Source and reproduction notes](figures/chapter3/other_rtqpcr_processing_notes.md)


## Retired - HT1080 standalone ddPCR figure

CH3-S3.4.5.1-F02 has been incorporated into CH3-S3.4.5.1-F01, panels b and c. Previous captions are preserved in [the pre-merge archive](archive/ht1080_merged_figure_20260923/Chapter3_figure_captions.md).

## main - Fig 3.9. Y-targeting reduces Y-chromosome DNA and expression in HT1080.

<a id="fig-3-9"></a>

**Figure:** `ch03_main_ht1080_lpwgs_genomewide` | **Working ID:** CH3-S3.4.5.1-F03 | **Chapter number:** Fig 3.9

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_1_F03_ht1080_genomewide/CH3_S3_4_5_1_F03_ht1080_genomewide.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_1_F03_ht1080_genomewide/CH3_S3_4_5_1_F03_ht1080_genomewide.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_1_F03_ht1080_genomewide/CH3_S3_4_5_1_F03_ht1080_genomewide.png)

[PDF](figures/chapter3/ch03_main_ht1080_lpwgs_genomewide.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_lpwgs_genomewide.svg) | [PNG](figures/chapter3/ch03_main_ht1080_lpwgs_genomewide.png)

**Fig 3.9. Y-targeting reduces Y-chromosome DNA and expression in HT1080.**

HT1080 Y-depleted (two CP + Q53 nucleofections) and sgCTRL populations. **a,** Y-depleted-minus-sgCTRL lpWGS log₂ relative copy-number differences across 630 5-Mb windows. Inset, Y-chromosome detail before Yq12, showing back-transformed 1-Mb median signals in sgCTRL (blue circles) and Y-depleted cells (red triangles); grey shading marks the excluded centromere. **b,c,** DNA–RNA differences across 1,597 1-Mb windows (b) and 43 chromosome arms (c), with Pearson correlations including/excluding Y. Red triangles denote Y and grey circles autosomes/X; in c, 13q is green and labelled non-Y deviations are blue. Panels a,b use evenly spaced log₂ units. Panel c shares b’s linear RNA axis, while its DNA axis retains symmetric-log spacing, linear within ±0.3. The inset shows unlogged relative copy number. RNA values are TMM-normalised mean log₂ CPM differences, equally weighted across two periods. Non-Y windows/arms require ≥ 5/≥ 20 expressed genes; Y regions require ≥ 1. Coordinates follow CHM13v2.0; DNA uses usable 100-kb bins, Yq12 is excluded (grey in a), and RNA also excludes PAR. Each population has one lpWGS sample; DNA–RNA specimen matching is not established. Floor-level DNA values do not measure precise fold depletion.

[Source and reproduction notes](figures/chapter3/source_data/ht1080_lpwgs_difference_draft/README.md)

## supp - Supp Fig 3.41. Additional copy-number differences in HT1080.

<a id="supp-fig-3-41"></a>

**Figure:** `ch03_supp_ht1080_lpwgs_copy_number_detail` | **Working ID:** CH3-S3.4.5.1-S04 | **Chapter number:** Supp Fig 3.41 (provisional addition during draft 2 review; final first-mention ordering pending)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_1_S04_ht1080_lpwgs_copy_number_detail/CH3_S3_4_5_1_S04_ht1080_lpwgs_copy_number_detail.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_1_S04_ht1080_lpwgs_copy_number_detail/CH3_S3_4_5_1_S04_ht1080_lpwgs_copy_number_detail.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_1_S04_ht1080_lpwgs_copy_number_detail/CH3_S3_4_5_1_S04_ht1080_lpwgs_copy_number_detail.png)

**Supp Fig 3.41. Additional copy-number differences in HT1080.**

Y-depleted-minus-sgCTRL differences after two CP + Q53 nucleofections. **a,** Median paired-bin lpWGS differences across 46 autosomal/X arms; green marks 13q and hatching poorly resolved acrocentric short arms. **b,c,** Chromosomes 13 (b) and 4 (c) in 1-Mb windows, excluding centromere-overlapping bins (grey). The distal-13q segment difference is +0.368 log₂ (1.29-fold; green); chromosome-4 arm medians are −0.281 on 4p and −0.237 on 4q (0.82- and 0.85-fold; blue). Coordinates follow CHM13v2.0; DNA uses usable 100-kb bins. One lpWGS sample represents each population.

[Source and reproduction notes](figures/chapter3/source_data/ht1080_lpwgs_difference_draft/README.md)

## main - Fig 3.10. Residual Yp DNA and expression persist for a subset of genes in Y-depleted HT1080 cells.

<a id="fig-3-10"></a>

**Figure:** `ch03_main_ht1080_lpwgs_chrY_atlas` | **Working ID:** CH3-S3.4.5.1-F04 | **Results:** 3.4.5.1 | **Title:** Residual Yp DNA and expression persist for a subset of genes in Y-depleted HT1080 cells. | **Chapter number:** Fig 3.10 | **Current draft 2 number:** Fig 3.4

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_1_F04_ht1080_regional_y/CH3_S3_4_5_1_F04_ht1080_regional_y.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_1_F04_ht1080_regional_y/CH3_S3_4_5_1_F04_ht1080_regional_y.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_1_F04_ht1080_regional_y/CH3_S3_4_5_1_F04_ht1080_regional_y.png)

[PDF](figures/chapter3/ch03_main_ht1080_lpwgs_chrY_atlas.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_lpwgs_chrY_atlas.svg) | [PNG](figures/chapter3/ch03_main_ht1080_lpwgs_chrY_atlas.png)

**Fig 3.10. Residual Yp DNA and expression persist for a subset of genes in Y-depleted HT1080 cells.**

**a,** lpWGS retention after two CP + Q53 nucleofections, relative to sgCTRL. Grey points, corrected 100-kb Y bins; red lines, regional medians; teal, combined X/Y PAR coverage mapped to X using a Y-PAR-masked reference. Hatched regions are excluded; Yq12 is compressed. One DNA sample per population; <1% Yq signal does not establish complete absence.

**b,** HG002 T2T-Y ideogram aligned to a, with RNA-seq genes above and ddPCR targets below.

**c,** RNA expression relative to mean sgCTRL within each sampling period (blue, sgCTRL; red, Y-depleted). Circles, August 2025; triangles, May 2026; three passage samples per condition per period. Horizontal marks show six-sample means. Teal identifies PAR genes, whose expression is not Y-specific.

Dashed lines mark 100% of sgCTRL. Samples derive from one established population per condition. Inferred XY-equivalent dosage predicts ~50% PAR retention after complete Y-PAR loss with unchanged X dosage; this is not a confirmed karyotype. Filtering details are in **Supplementary Methods S3.1.3**.

[Source and reproduction notes](figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md)

## main - Y-linked and pseudoautosomal expression in HT1080

**Figure:** `ch03_main_ht1080_y_expression_confirmation` | **Results:** 3.4.5.1

[PDF](figures/chapter3/ch03_main_ht1080_y_expression_confirmation.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_y_expression_confirmation.svg) | [PNG](figures/chapter3/ch03_main_ht1080_y_expression_confirmation.png)

**a,** Positions of the selected genes on the HG002 T2T-Y map, with connectors to chromosome-ordered gene columns. **b,** Expression in sgCTRL and Y-depleted HT1080 populations. Bars show mean log2(CPM + 1) across six passage samples per condition; points identify August 2025 (circles) and May 2026 (triangles), with three passage samples per condition in each window. **c,** The same values for all 12 samples on a common magma colour scale, without row standardisation or clustering. Sample indices denote within-window passage positions. These samples derive from one original edited pool and one original control pool and are not independent nucleofection replicates.

The panel retains eight expressed male-specific Y genes and three selected pseudoautosomal genes (SLC25A6, CD99 and VAMP7), using the existing criterion of CPM >= 1 in at least three of the 12 samples in either condition. AMELY, SRY, TBL1Y and TMSB4Y did not meet that filter. CPM is the existing gene count divided by the sum of all assigned gene counts in that library, multiplied by one million. Quantification uses GRCh38 PAR-fix counts; the gene-position illustration uses the T2T-Y annotation. Pseudoautosomal genes are quantified at chrX in this workflow and do not provide Y-specific expression measurements. The display preserves the original Ben_Thesis counts and normalization, without batch subtraction, TMM correction or significance testing. Ratios quoted in the Results use within-window mean CPM, not ratios of log-transformed bar heights. Residual RNA signal does not establish the proportion of Y-retaining cells or the physical location of retained DNA.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_lpwgs_confirmation_notes.md)

## Retired main - Former Fig 3.11: A549 ddPCR panels merged into Fig 3.8

<a id="fig-3-11"></a>

On 27 October 2026 the author adopted the combined figure. Former Fig 3.11a,b are now **Fig 3.8d,e**; use the [combined figure and caption](#fig-3-8). The standalone figure is retired from draft 2. Its exports and numerical sources remain historical records.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_2_F01_a549_ddpcr/CH3_S3_4_5_2_F01_a549_ddpcr.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_2_F01_a549_ddpcr/CH3_S3_4_5_2_F01_a549_ddpcr.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_2_F01_a549_ddpcr/CH3_S3_4_5_2_F01_a549_ddpcr.png)

[Previous caption and chapter state](archive/fig3_8_ht1080_a549_adopted_20261027_6de881a2/Chapter3_figure_captions.md#fig-3-11)

## supp - A549 unscaled ddPCR concentration ratios

**Figure:** `ch03_supp_a549_ddpcr_raw_ratios`

[PDF](figures/chapter3/ch03_supp_a549_ddpcr_raw_ratios.pdf) | [Editable SVG](figures/chapter3/ch03_supp_a549_ddpcr_raw_ratios.svg) | [PNG](figures/chapter3/ch03_supp_a549_ddpcr_raw_ratios.png)

**a,** AMELY/AMELX and **b,** DDX3Y/RNaseP concentration ratios in A549 populations treated with CP, Q53 or CP + Q53. Labels 1× and 2× denote one nucleofection and two successive nucleofections, respectively. Vermilion bars indicate guide-treated populations; blue bars indicate the non-targeting-guide control (sgCTRL). Hatching marks CP + Q53 2×. Bars and printed values show medians, whiskers show interquartile ranges, and points show individual assay measurements. Printed n denotes measurements, including repeated assays of the same DNA specimen, rather than independent nucleofections.

Measurements were obtained on 26 and 27 January, 24 and 27 February, and 9 and 25 March 2025. CP and Q53 measurements derive from DNA specimens E1 and E3, respectively. Single-nucleofection CP + Q53 measurements derive from E5, E6, H4, H7, I1, I2, K6 and L1; assay availability differs between specimens. Double-nucleofection measurements derive from L3 and L6 for AMELY/AMELX and L6 for DDX3Y/RNaseP. Control specimens are E2, E7, H3, H9, I3 and I5. The two I1 AMELY measurements are excluded, retaining the exclusion in the saved A549 analysis; its DDX3Y measurements are retained. The comparison combines measurements collected on different dates and does not establish the number of independently derived populations.

Ratios were calculated from instrument-reported target and reference concentrations, without copy-number scaling or normalisation to a control mean. The y-axis scales differ between assays. No inferential tests were performed. These locus-specific measurements do not establish complete chromosome loss.

[Source and reproduction notes](figures/chapter3/ch03_main_a549_ddpcr_nucleofection_rounds_notes.md)

The corresponding [main figure](figures/chapter3/ch03_main_a549_ddpcr_nucleofection_rounds.pdf) converts these same measurements to estimated target copies per cell using reference-locus copy numbers.

## main - Fig 3.13. DU145 ddPCR comparison of initial and second-nucleofection populations.

<a id="fig-3-13"></a>

**Figure:** `ch03_main_du145_ddpcr_assay_comparison` | **Working ID:** CH3-S3.4.5.3-F01 | **Results:** 3.4.5.3 | **Title:** ddPCR reveals Y chromosome DNA depletion in DU145. | **Chapter number:** Fig 3.13

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_3_F01_du145_ddpcr/CH3_S3_4_5_3_F01_du145_ddpcr.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_3_F01_du145_ddpcr/CH3_S3_4_5_3_F01_du145_ddpcr.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_3_F01_du145_ddpcr/CH3_S3_4_5_3_F01_du145_ddpcr.png)

[PDF](figures/chapter3/ch03_main_du145_ddpcr_assay_comparison.pdf) | [Editable SVG](figures/chapter3/ch03_main_du145_ddpcr_assay_comparison.svg) | [PNG](figures/chapter3/ch03_main_du145_ddpcr_assay_comparison.png)

**Fig 3.13. DU145 ddPCR comparison of initial and second-nucleofection populations.**

**a,** AMELY/AMELX ddPCR ratios following one (1×) or two successive (2×) CP + Q53 + P47 nucleofections, with corresponding sgCTRL populations and parental cells in the 1× comparison.

**b,** DDX3Y/RNaseP ratios after one nucleofection and KDM5D/RNaseP ratios after two. KDM5D was assessed in place of DDX3Y after the second nucleofection.

Bars show medians, whiskers interquartile ranges where available, and points individual measurements; n denotes assay measurements. Ratios are unscaled, and controls are shown separately for each assay and round.

[Source and reproduction notes](figures/chapter3/ch03_main_du145_ddpcr_assay_comparison_notes.md)

## main - A549 Y-linked DNA dosage over time

**Figure:** `ch03_main_a549_ddpcr_timecourse` | **Results:** 3.4.5.2

[PDF](figures/chapter3/ch03_main_a549_ddpcr_timecourse.pdf) | [Editable SVG](figures/chapter3/ch03_main_a549_ddpcr_timecourse.svg) | [PNG](figures/chapter3/ch03_main_a549_ddpcr_timecourse.png)

Estimated Y-linked target dosage in A549 populations after one nucleofection with CP + Q53 (vermilion squares) or a non-targeting guide (sgCTRL; blue circles). **a,** AMELY, measured relative to AMELX. **b,** DDX3Y, measured relative to RNaseP. The x-axis indicates DNA collection date: 24 and 31 January and 3, 5 and 19 February 2025. Lines connect date-specific medians, whiskers indicate interquartile ranges, and hollow points show individual assay observations with small horizontal offsets for visibility. The dashed line indicates one target copy per cell.

AMELY observations per date are 3, 1, 1, 1 and 2 for CP + Q53 and 2, 1, 1, 1 and 5 for sgCTRL. DDX3Y observations are 3, 1, 1, 1 and 3 for CP + Q53 and 2, 1, 1, 2 and 4 for sgCTRL, respectively. These counts include repeat assays of DNA specimens and do not indicate independent transfections. Collection dates are distinct from ddPCR run dates. Lines describe condition summaries across sampled dates; no inferential tests were performed.

Instrument concentration ratios were multiplied by two AMELX copies or three RNaseP copies, assuming stable reference dosage after editing. Both panels share a copies-per-cell scale; values are not normalised to a control mean or interpreted as Y-negative cell fractions. This uses the measurement convention adopted for the thesis A549 guide-comparison figure, replacing the historical time-course script's positive-droplet ratio and half-scaled DDX3Y values. The two I1 AMELY readings remain excluded according to the retained analysis, while I1 DDX3Y readings remain included. Double-nucleofection populations are shown separately in the guide-comparison figure.

[Source and reproduction notes](figures/chapter3/ch03_main_a549_ddpcr_timecourse_notes.md)


## main - Fig 3.12. Regional Y-chromosome depletion and expression in A549.

<a id="fig-3-12"></a>

**Figure:** `ch03_main_a549_lpwgs_chrY_atlas` | **Working ID:** CH3-S3.4.5.2-F02 | **Results:** 3.4.5.2 | **Title:** DNA and RNA profiling supports Y chromosome loss in A549. | **Chapter number:** Fig 3.12 | **Current draft 2 number:** Fig 3.5

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_2_F02_a549_regional_y/CH3_S3_4_5_2_F02_a549_regional_y.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_2_F02_a549_regional_y/CH3_S3_4_5_2_F02_a549_regional_y.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_2_F02_a549_regional_y/CH3_S3_4_5_2_F02_a549_regional_y.png)

[PDF](figures/chapter3/ch03_main_a549_lpwgs_chrY_atlas.pdf) | [Editable SVG](figures/chapter3/ch03_main_a549_lpwgs_chrY_atlas.svg) | [PNG](figures/chapter3/ch03_main_a549_lpwgs_chrY_atlas.png)

**Fig 3.12. Regional Y-chromosome depletion and expression in A549.**

**a,** lpWGS retention after two CP + Q53 nucleofections, relative to sgCTRL. Grey points, corrected 100-kb Y bins; red lines, regional medians; teal, combined X/Y PAR coverage mapped to X using a Y-PAR-masked reference. Hatched regions are excluded; Yq12 is compressed. One DNA sample per population.

**b,** HG002 T2T-Y ideogram aligned to a, with RNA-seq genes above and ddPCR targets below.

**c,** RNA expression relative to mean sgCTRL (blue, sgCTRL; red, Y-depleted). Points show four samples per condition; horizontal marks show group means. Teal identifies PAR genes, whose expression is not Y-specific.

Dashed lines mark 100% of sgCTRL. Samples derive from one established population per condition. Inferred XXY-equivalent dosage predicts ~67% PAR retention after complete Y-PAR loss with unchanged X dosage; this is not a confirmed karyotype. Filtering details are in **Supplementary Methods S3.1.3**.

[Source and reproduction notes](figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md)

## main - Fig 3.14. Regional Y-chromosome depletion and expression in DU145.

<a id="fig-3-14"></a>

**Figure:** `ch03_main_du145_lpwgs_chrY_atlas` | **Working ID:** CH3-S3.4.5.3-F02 | **Results:** 3.4.5.3 | **Title:** Distal Yp DNA and expression persist after Y chromosome targeting in DU145. | **Chapter number:** Fig 3.14 | **Current draft 2 number:** Fig 3.6

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_3_F02_du145_regional_y/CH3_S3_4_5_3_F02_du145_regional_y.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_3_F02_du145_regional_y/CH3_S3_4_5_3_F02_du145_regional_y.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_3_F02_du145_regional_y/CH3_S3_4_5_3_F02_du145_regional_y.png)

[PDF](figures/chapter3/ch03_main_du145_lpwgs_chrY_atlas.pdf) | [Editable SVG](figures/chapter3/ch03_main_du145_lpwgs_chrY_atlas.svg) | [PNG](figures/chapter3/ch03_main_du145_lpwgs_chrY_atlas.png)

**Fig 3.14. Regional Y-chromosome depletion and expression in DU145.**

**a,** lpWGS retention after two CP + Q53 + P47 nucleofections, relative to sgCTRL. Grey points, corrected 100-kb Y bins; red lines, regional medians; teal, combined X/Y PAR coverage mapped to X using a Y-PAR-masked reference. Hatched regions are excluded; Yq12 is compressed. One DNA sample per population. Distal/proximal Yp summaries do not define a precise breakpoint.

**b,** HG002 T2T-Y ideogram aligned to a, with RNA-seq genes above and ddPCR targets below.

**c,** RNA expression relative to mean sgCTRL (blue, sgCTRL; red, Y-depleted). Points show three passage samples per condition; horizontal marks show group means. NLGN4Y is NA because control expression was zero. Teal identifies PAR genes, whose expression is not Y-specific.

Dashed lines mark 100% of sgCTRL. Samples derive from one established population per condition. Using local X:Yq dosage and the higher X dosage in the depleted pool, complete Y-derived PAR loss predicts ~73–74% PAR retention. This conditional dosage model does not establish the structural karyotype. Filtering and dosage details are in **Supplementary Methods S3.1.3**.

[Source and reproduction notes](figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md)

## supp - Supp Fig 3.42. Separate chromosome-Y copy-number profiles in HT1080, A549 and DU145.

<a id="supp-fig-3-42"></a>

**Figure:** `ch03_supp_regional_y_individual_tracks` | **Working ID:** CH3-S3.4.5-S05 | **Results:** 3.3.2 | **Chapter number:** Supp Fig 3.42 (provisional)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_S05_regional_y_individual_tracks/CH3_S3_4_5_S05_regional_y_individual_tracks.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_S05_regional_y_individual_tracks/CH3_S3_4_5_S05_regional_y_individual_tracks.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_S05_regional_y_individual_tracks/CH3_S3_4_5_S05_regional_y_individual_tracks.png)

**Supp Fig 3.42. Separate chromosome-Y copy-number profiles in HT1080, A549 and DU145.**

**a–f,** Individual lpWGS profiles for sgCTRL (blue) and Y-depleted populations (red) in HT1080 (a,b), A549 (c,d) and DU145 (e,f), following two successive CP + Q53 nucleofections in HT1080 and A549 or CP + Q53 + P47 nucleofections in DU145. Points show corrected relative copy number in the earlier unmasked-reference analysis and its original bin selection; coloured lines show regional medians within each population. These historical tracks differ from the updated main atlases, which use Y-PAR-masked realignments, revised filtering and X-mapped PAR coverage. These values retain the original genome-normalised copy-number scale and are not divided by sgCTRL. All panels use the same linear axes and full, uncompressed HG002 T2T-Y coordinates. Triangles mark bins above the displayed range of 1.3; their counts and maximum values are indicated. Teal shading marks pseudoautosomal regions. Hatched shading marks the centromere and Yq12, which are excluded rather than plotted as DNA loss. The numerical floor of 0.001 is not a measured residual copy number or a detection limit. One DNA sample was sequenced per population; bins are not biological replicates.

[Source and reproduction notes](figures/chapter3/source_data/regional_y_atlases_normalised/README.md)

## supp - Supp Fig 3.10. Distal Yp DNA persists after Y-chromosome targeting in DU145.

<a id="supp-fig-3-10"></a>

**Figure:** `ch03_supp_du145_lpwgs_PAR1_Yp_zoom` | **Working ID:** CH3-S3.4.5.3-S01 | **Results:** 3.4.5.3 | **Title:** Distal Yp DNA persists after Y-chromosome targeting in DU145. | **Chapter number:** Supp Fig 3.10

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_3_S01_du145_distal_yp/CH3_S3_4_5_3_S01_du145_distal_yp.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_3_S01_du145_distal_yp/CH3_S3_4_5_3_S01_du145_distal_yp.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_3_S01_du145_distal_yp/CH3_S3_4_5_3_S01_du145_distal_yp.png)

[PDF](figures/chapter3/ch03_supp_du145_lpwgs_PAR1_Yp_zoom.pdf) | [Editable SVG](figures/chapter3/ch03_supp_du145_lpwgs_PAR1_Yp_zoom.svg) | [PNG](figures/chapter3/ch03_supp_du145_lpwgs_PAR1_Yp_zoom.png)

**Supp Fig 3.10. Distal Yp DNA persists after Y-chromosome targeting in DU145.**

**a,b,** lpWGS relative copy-number profiles across Yp in sgCTRL and the CP + Q53 + P47 population after two successive nucleofections. Grey points show corrected 100-kb bins; coloured lines show segmented copy-number estimates. Teal shading marks PAR1 and grey shading marks the centromere. Brackets identify retained-signal intervals in the edited population. One sample was sequenced per condition.

**c,** Aligned cytoband and PAR1 tracks with selected genes, including RPS4Y1, ZFY and AMELY. Coordinates follow HG002 T2T-Y; segment boundaries do not establish precise structural breakpoints.

[Source and reproduction notes](figures/chapter3/ch03_supp_du145_lpwgs_PAR1_Yp_zoom_notes.md)

## main - H2452 initial-guide ddPCR comparison

**Figure:** `ch03_main_h2452_ddpcr_initial_guides` | **Results:** 3.4.5.4

[PDF](figures/chapter3/ch03_main_h2452_ddpcr_initial_guides.pdf) | [Editable SVG](figures/chapter3/ch03_main_h2452_ddpcr_initial_guides.svg) | [PNG](figures/chapter3/ch03_main_h2452_ddpcr_initial_guides.png)

Recorded AMELY/AMELX ratios from the initial H2452 individual-guide experiment, combining Batch 1 and Batch 2. CP, CN and Q53 populations received one nucleofection (1x). Vermilion indicates guide-treated populations, blue the CTRL2 control (displayed as sgCTRL), and grey parental H2452. Bars and printed values show medians of the combined measurements; whiskers show interquartile ranges. Hollow circles identify Batch 1 and hollow triangles Batch 2.

The figure includes 22 entries from `summary_table_ddPCR_edit.csv`: six CP measurements and four each for CN, Q53, sgCTRL and parental cells. CP measurements derive from DNA specimen A7, CN from B1, Q53 from A6 and sgCTRL from A1; parental specimens retain the source label H2452. The batches are repeat assay batches, not independent nucleofections. Printed n counts measurements rather than biological replicates. Each measurement receives equal weight in the median and IQR; CP has four measurements in Batch 1 and two in Batch 2, while every other condition has two per batch.

Values reproduce the summary table's ratio field without conversion to copies per cell or normalisation to a control mean. Historical ratio conventions and target/channel labels remain subject to source reconciliation; pooling follows the author's selected comparison and does not establish an independently validated concentration-ratio analysis. No inferential tests were performed. Other assay dates and the later five-guide experiment are not included in this plot.

[Source and reproduction notes](figures/chapter3/ch03_main_h2452_ddpcr_initial_guides_notes.md)


## main - Fig 3.15. Panobinostat treatment is associated with greater Y-chromosome depletion in DU145.

<a id="fig-3-15"></a>

**Figure:** `ch03_main_du145_drug_treatment_comparison` | **Working ID:** CH3-S3.4.5.5-F01 | **Results:** 3.4.5.5 | **Title:** Panobinostat treatment is associated with greater Y-chromosome depletion in DU145. | **Chapter number:** Fig 3.15

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_5_F01_du145_all_treatments/CH3_S3_4_5_5_F01_du145_all_treatments.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_5_F01_du145_all_treatments/CH3_S3_4_5_5_F01_du145_all_treatments.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_5_F01_du145_all_treatments/CH3_S3_4_5_5_F01_du145_all_treatments.png)

[PDF](figures/chapter3/ch03_main_du145_drug_treatment_comparison.pdf) | [Editable SVG](figures/chapter3/ch03_main_du145_drug_treatment_comparison.svg) | [PNG](figures/chapter3/ch03_main_du145_drug_treatment_comparison.png)

**Fig 3.15. Panobinostat treatment is associated with greater Y-chromosome depletion in DU145.**

All drug conditions received one CP + Q53 + P47 nucleofection. Panobinostat (30 or 75 nM) or DMSO vehicle was applied for 24 h before nucleofection; KU-55933 (KU) was applied for 16 h afterwards.

**a,** Y-linked ddPCR target/reference ratios. Dotted lines mark the vehicle-treated, nucleofected population (1× DMSO).

**b,** Target-dosage reductions relative to sgCTRL (upper) or 1× DMSO (lower), calculated as 100 × (1 − condition/comparator). The separate two-round, no-drug population (2×) uses its own assay-run sgCTRL; DDX3Y was not measured.

**c,** Median lpWGS relative copy number across distal Yp (2.46–4.8 Mb), proximal Yp (6–10 Mb) and euchromatic Yq (10.883–27.450 Mb), in T2T-Y coordinates. Floor denotes the numerical floor of 0.001.

**d,** Depletion relative to sgCTRL in whole-Y/autosome coverage and median regional copy number across male-specific Yp (2.458–10.566 Mb) and euchromatic Yq. Points represent independent nucleofections without (*n* = 2) or with panobinostat (*n* = 4); horizontal lines show medians.

Each condition contributed one lpWGS sample and one ddPCR measurement per target.

[Source and reproduction notes](figures/chapter3/ch03_supp_du145_drug_treatment_comparison_notes.md)

## Retired - DU145 grouped panobinostat comparison

CH3-S3.4.5.5-F02 was removed from the active figure set by the author. The approved treatment comparison is CH3-S3.4.5.5-F01. Historical caption and exports are preserved in `archive/retire_du145_grouped_20260923_231413`.


## Retired - DU145 supplementary pharmacological comparison

CH3-S3.4.5.5-S02 is now the main figure CH3-S3.4.5.5-F01.

## Retired - DU145 whole-Y and regional DNA content

CH3-S3.4.5.5-S01 was removed from the active supplementary figures by the author. The main treatment comparison remains CH3-S3.4.5.5-F01. Historical caption and exports are preserved in `archive/retire_du145_y_content_20260924_085215`.


## supp - Supp Fig 3.8. Reference genomic and expression profiles of A549 and DU145.

<a id="supp-fig-3-8"></a>

**Figure:** `ch03_main_a549_du145_profiles` | **Working ID:** CH3-S3.4.5.2-S01 | **Results:** 3.4.5.2 | **Title:** Reference genomic and expression profiles of A549 and DU145. | **Chapter number:** Supp Fig 3.8

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_2_S01_a549_du145_profiles/CH3_S3_4_5_2_S01_a549_du145_profiles.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_2_S01_a549_du145_profiles/CH3_S3_4_5_2_S01_a549_du145_profiles.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_2_S01_a549_du145_profiles/CH3_S3_4_5_2_S01_a549_du145_profiles.png)

[PDF](figures/chapter3/ch03_main_a549_du145_profiles.pdf) | [Editable SVG](figures/chapter3/ch03_main_a549_du145_profiles.svg) | [PNG](figures/chapter3/ch03_main_a549_du145_profiles.png)

**Supp Fig 3.8. Reference genomic and expression profiles of A549 and DU145.**

**a,b,** Reported cell-line characteristics and DepMap Public 26Q1 profiles of A549 and DU145, respectively. Y-linked genes follow GRCh38 genomic order; blue bars show log₂(TPM + 1), with WGS copy-number estimates alongside.

Top selective dependencies are shown as Chronos gene effects (blue), with the median and interquartile range across 1,207 other models (grey). Genes had effects below −0.5 in the profiled line and no more than 20% of reference models, ranked by their difference from the reference median. More negative values indicate stronger dependency.

Genome-wide WGS profiles follow GRCh38 chromosome order. Blue identifies chromosome Y, and dashed lines mark two copies. Public A549 profiles show near-zero Y copy-number estimates despite detectable Y-linked expression.

[Source and reproduction notes](figures/chapter3/ch03_main_a549_du145_profiles_notes.md)

## main - Chromosome-Y copy number across the five study cell-line backgrounds

**Figure:** `ch03_main_selected_models_Y_copy_number` | **Results:** 3.4.5

[PDF](figures/chapter3/ch03_main_selected_models_Y_copy_number.pdf) | [Editable SVG](figures/chapter3/ch03_main_selected_models_Y_copy_number.svg) | [PNG](figures/chapter3/ch03_main_selected_models_Y_copy_number.png)

Baseline chromosome-Y copy-number profiles from DepMap Public 26Q1 for HT1080, A549, DU145, NCI-H2452 (displayed as H2452) and HepG2. Both panels show released linear relative copy-number ratios on a shared scale from 0 to 1. These are not absolute copies per cell: equal ratios need not imply equal absolute dosage when model ploidy differs. Relative ratios are retained for all five models because the default HepG2 WGS profile lacks a matched ploidy estimate. Each model contributes one default WGS profile. Public profiles describe the corresponding repository models, not the author's experimental stocks or edited populations.

**a,** Saved WGS segments across the reported interval, with the GRCh38 Yp, centromere and proximal Yq ideogram above. Horizontal lines represent released segment estimates; vertical connections mark adjacent segment transitions. Numbers give segment ratios to three decimal places for segments longer than 200 kb; values below 0.001 are labelled <0.001. The narrow A549 segment at 10.120-10.136 Mb (ratio 0.00637) is retained without a printed label. All five released profiles end at 26.637 Mb; the display extends to 27 Mb, and the unreported terminal interval is shaded and hatched. Distal Yq, including most of Yq12 and PAR2, is outside the reported segment coverage and must not be interpreted as zero copy number. H2452 shows a transition at 2.928 Mb from a ratio of 0.00330 to 0.60727. The transition is a public segmentation boundary, not an experimentally validated breakpoint. Light central shading marks the GRCh38 centromeric interval. Teal shading marks PAR1 at the distal Yp end (GRCh38 chrY:10,001-2,781,479, 1-based inclusive), with a dotted boundary on the tracks. PAR boundaries are distinct from cytoband boundaries. Although the released segment coordinates span PAR1, the segment table does not establish an independently measured PAR1 copy number or its upstream masking/mapping treatment. PAR2 is beyond the displayed and reported interval. [GRCh38 PAR coordinates and reference masking](https://gatk.broadinstitute.org/hc/en-us/articles/360041155232-Reference-Genome-Components).

**b,** Released gene-level relative copy-number values for nine Y-linked genes, ordered by GRCh38 start position from Yp to Yq. Colour uses the same 0-1 scale for every model; printed values are rounded to two decimal places, except values below 0.001. Gene-level and segment-level estimates are separate released products and are shown without forcing them to agree. In particular, the DU145 gene-level profile has higher values for RPS4Y1 and ZFY than for the remaining displayed genes despite a single released segment over this interval. H2452 has a near-zero RPS4Y1 estimate with higher values for the other genes.

A549 has near-zero Y-DNA estimates despite retained Y-linked RNA expression in its separately collected public RNA profile; the reviewed call table flags this disagreement. No stock-specific Y-loss conclusion is inferred from these public values. No inferential tests or replicate error bars are shown.

[Source and reproduction notes](figures/chapter3/ch03_main_selected_models_Y_copy_number_notes.md)

## main - HT1080 RNA-seq sampling design and global expression structure

**Historical ID:** `CH3-S3.4.6.1-F01` (superseded global-structure layout) | **Figure:** `ch03_main_ht1080_rna_global_structure` | **Results:** 3.4.6.1 | **Status:** Superseded; panels d, c and a incorporated into the combined expression figure, formerly F02 and now the active F01. Historical caption retained.

[PDF](figures/chapter3/ch03_main_ht1080_rna_global_structure.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_rna_global_structure.svg) | [PNG](figures/chapter3/ch03_main_ht1080_rna_global_structure.png)

**HT1080 RNA-seq sampling design and global expression structure.**

**a,** PCA of 15 libraries across 25,526 expressed autosomal and non-PAR X-linked genes. Expression was adjusted for sampling period while preserving population differences, for visualisation only. Numbers denote passages. **b,** Percentage of reads assigned to genes. **c,** Pearson correlations and hierarchical clustering of the 12 sgCTRL and Y-depleted libraries using TMM-normalised expression without period adjustment. **d,** Sampling design: one established population per condition, sampled at three successive passages in each of two periods. Populations underwent two nucleofections with non-targeting sgRNA or CP + Q53.

Blue denotes sgCTRL, red Y-depleted and grey parental samples; circles and triangles distinguish the two sampling periods. Sample codes indicate period (A, August; M, May), population (C, sgCTRL; L, Y-depleted; P, parental) and passage. Parental samples are included in a and b only. Passage samples are repeated measurements of the same established populations.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_rna_global_structure_notes.md)

## supp - HT1080 parental-inclusive principal component analysis

**Figure:** `ch03_supp_ht1080_parental_pca` | **Results:** 3.4.6.1

[PDF](figures/chapter3/ch03_supp_ht1080_parental_pca.pdf) | [Editable SVG](figures/chapter3/ch03_supp_ht1080_parental_pca.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_parental_pca.png)

Principal component analysis of all 15 HT1080 RNA-seq libraries: six sgCTRL, six Y-depleted and three parental samples. The contextual analysis fits its own TMM normalisation across all 15 libraries, with expression filtering and normalisation based on autosomal and non-pseudoautosomal X-linked genes. PCA uses the same 25,526-gene set retained for the primary 12-library comparison, log2 counts per million with prior count 2, gene centring and no variance scaling. PC1 and PC2 explain 37.7% and 12.4% of the variance, respectively. Blue denotes sgCTRL, vermilion Y-depleted and grey parental samples. Circles denote August 2025 and triangles May 2026; numbers identify passages within each population and sampling period. Parental samples were collected only in August and are contextual; their matching to the nucleofected populations is not established. The sgCTRL and Y-depleted libraries repeatedly sample one established population per condition across the two periods. This descriptive analysis does not alter the normalisation or expression effects of the primary comparison. No inferential tests are shown.

[Source and reproduction notes](figures/chapter3/ch03_supp_ht1080_parental_pca_notes.md)

## supp - Supp Fig 3.12. HT1080 expression structure and sgCTRL–parental comparisons.

<a id="supp-fig-3-12"></a>

**Figure:** `ch03_supp_ht1080_expression_structure_parental` | **Chapter number:** Supp Fig 3.12

**Working ID:** `CH3-S3.3.6.1-S01` | **Figure:** `ch03_supp_ht1080_expression_structure_parental` | **Results:** 3.3.4.1 | **Status:** selected; merged S01/S02 on 9 October 2026

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S01_HT1080_PCA_comparisons/CH3_S3_4_6_1_S01_HT1080_PCA_comparisons.pdf) | [Editable SVG](figures/chapter3/selected/CH3_S3_4_6_1_S01_HT1080_PCA_comparisons/CH3_S3_4_6_1_S01_HT1080_PCA_comparisons.svg) | [PNG](figures/chapter3/selected/CH3_S3_4_6_1_S01_HT1080_PCA_comparisons/CH3_S3_4_6_1_S01_HT1080_PCA_comparisons_correlation_f97de51dcc92.png)

**Supp Fig 3.12. HT1080 expression structure and sgCTRL–parental comparisons.**

**a,b,** PCA of 15 libraries before/after additive period adjustment preserving population differences. **c,d,** Separate unadjusted PCAs for August (nine libraries) and May (six), using the same 25,526 autosomal/non-PAR X genes. Adjustment is for visualisation only. Blue/red/grey denote sgCTRL/Y-depleted/parental; circles/triangles denote August/May; numbers give passages. Parental samples occur only in August. **e,f,** August sgCTRL-minus-parental comparisons, with three passage samples per group. **e,** Limma-voom analysis with BH correction; all genes meeting FDR < 0.05 and |log₂ fold change| ≥ 1 are labelled. Blue/red mark decreases/increases; dashed lines mark thresholds. **f,** NES and FDR for 17 sets selected by FDR (up to three per direction/collection), corrected across all 23,526 sets. H, Hallmark; GO, GO Biological Process; R, Reactome. Blue/red indicate lower/higher sgCTRL expression. **g,** Pearson correlations of unadjusted TMM expression across 12 sgCTRL/Y-depleted libraries, clustered by average linkage (1 − r). Codes give period (A/M), population (C/L) and passage. All samples repeatedly measure established populations.

[Source and reproduction notes](figures/chapter3/ch03_supp_ht1080_expression_structure_parental_notes.md)

## supp - HT1080 batch-aware gene-level differential expression

**Figure:** `ch03_supp_ht1080_batch_aware_de` | **Status:** report extension; not yet incorporated into thesis Results

[PDF](figures/chapter3/ch03_supp_ht1080_batch_aware_de.pdf) | [Editable SVG](figures/chapter3/ch03_supp_ht1080_batch_aware_de.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_batch_aware_de.png)

**supp | HT1080 batch-aware gene-level differential expression.** **a,** Combined passage-aware limma-voom log2 fold change versus BH-adjusted P value for all 25,526 genes. Blue and vermilion identify lower and higher expression at FDR <0.05; grey denotes other genes. The dashed line marks FDR 0.05. **b,** Model coefficients compared with the original descriptive equal-period mean differences. **c,** Coefficients from the passage-aware and simpler batch + condition models. **d,** August and May coefficients from the passage-aware model. Colours in b-d retain the combined-model categories; dashed diagonals indicate equal effects. Pearson correlations describe gene-level effect agreement. **e,** Voom mean-variance data for the passage-aware model and fitted trends from both designs. **f,** Raw combined-contrast P-value histograms for both models, with a logarithmic count axis. Model colours in e/f denote designs rather than effect direction. The passage-aware design uses six within-period passage-pair intercepts and two period-specific population effects, combined with equal weights; it has four residual degrees of freedom. The simpler model has nine. Both use original integer counts and the same core TMM exposure offsets. Robust empirical Bayes tests assess nonzero effects; BH adjustment spans 25,526 genes per model/contrast. Twelve libraries repeatedly sample one sgCTRL and one Y-depleted population. Inference assumes independent model residuals and does not establish independent population-generation replication. All core libraries were retained; existing network and pathway outputs are unchanged.

[Source and reproduction notes](figures/chapter3/ch03_supp_ht1080_batch_aware_de_notes.md)


## supp - HT1080 batch + condition volcano with PAR genes

**Figure:** `ch03_supp_ht1080_batch_condition_volcano` | **Status:** PAR-inclusive companion analysis

[PDF](figures/chapter3/ch03_supp_ht1080_batch_condition_volcano.pdf) | [Editable SVG](figures/chapter3/ch03_supp_ht1080_batch_condition_volcano.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_batch_condition_volcano.png)

**supp | Batch-adjusted expression differences in HT1080, including pseudoautosomal genes.** Estimated log2 fold changes in Y-depleted relative to sgCTRL cells are plotted against negative log10 BH-adjusted P values. Blue and vermilion identify non-PAR genes with lower and higher expression, respectively, satisfying FDR < 0.05 and absolute log2 fold change >= 1; grey denotes other non-PAR genes. Purple diamonds identify all expressed features overlapping the annotated pseudoautosomal regions (PAR), irrespective of statistical significance or effect size. Dashed lines mark the FDR and fold-change thresholds. Overall, 437 features had lower and 232 higher expression at both thresholds, including four PAR features. Labels identify the 18 lowest-FDR and four largest-magnitude qualifying non-PAR features in each direction, plus expressed protein-coding PAR genes and all PAR features passing both cutoffs; overlapping selections are labelled once.

The companion limma-voom model includes batch and condition, using the same 12 core libraries and TMM effective library sizes as the preceding analysis. The existing expression filter retained 25,526 autosomal/non-PAR X features and 27 PAR features (25,553 tested features). Voom weights and robust empirical Bayes moderation were fitted across this expanded set, with BH correction across all 25,553 tests. The original PAR-excluded model is preserved; the count differs from its 662 qualifying genes because four qualifying PAR features were added and three non-PAR features crossed the threshold after refitting. Two-sided moderated tests assess a nonzero condition effect; the fold-change threshold is an additional filter on estimated effects. The PAR designation uses the existing GRCh38 interval-overlap annotation, including boundary-overlapping features; all 27 expressed PAR features are represented by chrX annotation records, and their counts do not establish X- versus Y-derived transcription. Libraries repeatedly sample one established population per condition, so model-based inference does not establish reproducibility across independently generated populations.

[Source and reproduction notes](figures/chapter3/ch03_supp_ht1080_batch_condition_volcano_notes.md)


## retired supp - HT1080 RNA programme analysis overview

**Figure:** `ch03_supp_ht1080_rna_programme_overview` | **Status:** merged into CH3-S3.3.6.1-S01, panels e–f, on 9 October 2026. Standalone exports and the following caption are historical.

[PDF](figures/chapter3/ch03_supp_ht1080_rna_programme_overview.pdf) | [Editable SVG](figures/chapter3/ch03_supp_ht1080_rna_programme_overview.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_rna_programme_overview.png)

**Gene-level and pathway differences between sgCTRL and parental HT1080 cells.** Both panels compare sgCTRL with parental cells sampled in August 2025, with three passage samples per group.

**a,** Differential-expression volcano plot. Blue and red identify genes with lower and higher expression in sgCTRL, respectively, at FDR <0.05 and absolute log2 fold change >=1; grey indicates other genes. All 12 qualifying genes are labelled. Dashed lines mark the FDR and effect-size thresholds. Expression differences were estimated using limma-voom, with Benjamini-Hochberg correction across tested genes.

**b,** Top GSEA hits from Hallmark (H), GO Biological Process (GO) and Reactome (R). Up to three significant sets per direction per collection were selected by FDR, with nominal P value and pathway name resolving ties. Hallmark has only two significant positive sets; both are included. Heatmap colours and cell labels show normalised enrichment scores (NES); the adjacent column gives FDR. Blue denotes enrichment toward lower expression in sgCTRL and red enrichment toward higher expression. White separators distinguish collections. All 17 displayed sets meet FDR <0.05, calculated across all 23,526 sets tested in this contrast rather than separately for the displayed collections. GSEA used saved moderated-t rankings and sets with 15-500 represented genes. Signature suffixes such as DN are part of the set name and do not directly specify pathway activity.

Samples represent repeated measurements of established populations rather than independent editing experiments. Saved differential-expression and enrichment results were used without refitting either analysis.

[Source and reproduction notes](figures/chapter3/ch03_supp_ht1080_rna_programme_overview_notes.md)

## main - Fig 3.16. HT1080 sampling design, global expression structure and expression differences.

<a id="fig-3-16"></a>

**Figure:** `ch03_main_ht1080_rna_expression_changes` | **Chapter number:** Fig 3.16

**Working ID:** `CH3-S3.4.6.1-F01` | **Figure:** `ch03_main_ht1080_rna_expression_changes` | **Status:** selected main figure

[PDF](figures/chapter3/ch03_main_ht1080_rna_expression_changes.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_rna_expression_changes.svg) | [PNG](figures/chapter3/ch03_main_ht1080_rna_expression_changes.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_F01_HT1080_expression_changes/CH3_S3_4_6_1_F01_HT1080_expression_changes.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_1_F01_HT1080_expression_changes/CH3_S3_4_6_1_F01_HT1080_expression_changes.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_1_F01_HT1080_expression_changes/CH3_S3_4_6_1_F01_HT1080_expression_changes_PAR_no_CSF2RA_6042e01ddab5.png)

**Fig 3.16. HT1080 sampling design, global expression structure and expression differences.**

**a,** One sgCTRL and one Y-depleted population, each generated by two nucleofections, sampled at three passages in August 2025 and May 2026 (12 libraries). These are repeated population samples. **b,** Fourteen Y genes and three leading genes per direction across all 15 libraries, including three August parental samples: unadjusted log₂(CPM + 1), without row scaling. Both axes use average-linkage clustering with Pearson correlation distance. Y genes require mean sgCTRL CPM ≥ 1 in either period; higher/lower genes have the lowest FDR among those meeting the criteria in d. Labels are grey/red/blue for Y-linked/higher/lower genes. Strips identify period/population; sample codes give period (A/M), population (P, parental; C, sgCTRL; L, Y-depleted) and passage. **c,** PCA of 25,526 autosomal/non-PAR X genes after removing an additive period effect using limma::removeBatchEffect while preserving population differences, for visualisation only. **d,** Period-adjusted limma-voom analysis, with BH correction across these genes. Blue/red mark lower/higher expression at FDR < 0.05, |log₂ fold change| ≥ 1 and matching period directions; grey marks other genes. Dashed lines mark thresholds. **e,** Period-specific effects, coloured as in d; the diagonal denotes equality. **f,** PAR expression relative to sgCTRL: PAR1 (cyan), PAR2 (blue-grey), and summed PAR1 + PAR2 (grey). Solid/hatched bars denote August/May. Bars use group-mean CPM from three passages per condition; individual genes require sgCTRL CPM ≥ 5. CSF2RA is omitted individually but retained in regional sums. PAR counts combine X/Y transcripts and are descriptive; the dashed line marks 100%. **g–l,** Unadjusted TMM-normalised log₂ CPM (prior count 2) for the indicated genes in all 15 samples; axes differ. Boxes show medians and interquartile ranges (IQR); whiskers reach observations within 1.5 IQR.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_rna_expression_changes_notes.md)

## main - HT1080 pathway and regulator patterns across sampling periods

**Figure:** `ch03_main_ht1080_rna_programmes_regulators` | **Status:** working alternative; not selected. Panel a is superseded by the GSEA-only figure and panel b is not currently selected.

[PDF](figures/chapter3/ch03_main_ht1080_rna_programmes_regulators.pdf) | [Editable SVG](figures/chapter3/ch03_main_ht1080_rna_programmes_regulators.svg) | [PNG](figures/chapter3/ch03_main_ht1080_rna_programmes_regulators.png)

**a,** Normalised enrichment scores for selected directionally consistent pathways in the overall and period-specific comparisons. Displayed sets are the four smallest-overall-FDR discoveries from each of Hallmark, Reactome and GO Biological Process, plus the three smallest-FDR positively enriched sets, with duplicates removed (15 sets). **b,** CollecTRI ULM target-pattern scores for the 12 smallest-overall-FDR directionally consistent regulators. Colours distinguish positive (vermilion) and negative (blue) scores, with separate scales for the two panels. Stars denote FDR <0.05 in the corresponding contrast; the overall FDR is printed beside each row. Complete gene-set and regulator results are retained in the programme-analysis tables.

Pathway enrichment uses signed moderated t rankings of all 25,461 uniquely represented autosomal/non-PAR X symbols; 4,868 gene sets with 15-500 represented members were tested, with BH correction across the three collections per contrast. Regulatory inference uses model log-fold changes and at least five CollecTRI targets, with BH correction across tested regulators per contrast. Display eligibility requires overall FDR <0.05 and matching score directions in August and May. The gene-level absolute log2FC cutoff is not applied to these ranked analyses. Target-pattern scores do not directly measure TF abundance, occupancy or activity, and source pathway names do not establish cell identity. The accompanying CAMERA analysis identifies 17 sets at FDR <0.05 using the default inter-gene correlation of 0.01.

[Source and reproduction notes](figures/chapter3/ch03_main_ht1080_rna_programmes_regulators_notes.md)


## main - HT1080 programme and regulator patterns across all MSigDB collections

**Figure:** `ch03_main_ht1080_rna_programmes_regulators_all_collections` | **Status:** all-collection alternative main draft; three-collection version retained

[PDF](figures/chapter3/ch03_main_ht1080_rna_programmes_regulators_all_collections.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_rna_programmes_regulators_all_collections.svg) | [PNG](figures/chapter3/ch03_main_ht1080_rna_programmes_regulators_all_collections.png)

**a,** Representative GSEA pathways and signatures grouped into eight author-selected biological themes, with nine positively enriched sets (Up in Y-depleted) and twelve negatively enriched sets (Down in Y-depleted). Upregulated themes comprise axon guidance and synaptic organisation, Polycomb-associated developmental signatures, N-glycan processing, and sterol/cholesterol responses. Downregulated themes comprise interferon and cytokine responses, extracellular matrix and adhesion, migration/cytoskeleton/wound response, and growth-factor signalling. These are interpretive groupings rather than newly computed network modules. Displayed sets were explicitly selected for biological relevance within the agreed themes and all pass overall global FDR <0.05 with matching NES directions in August and May; significance in both periods was not required. Colours indicate normalised enrichment scores (NES), using one shared scale for the up and down groups; stars indicate global FDR <0.05 in the indicated contrast. Overall FDR is printed beside each row. H3K27me3-associated and PRC2-target signatures refer to gene-set membership and do not measure histone modification or PRC2 activity in these cells. Neuronal signatures do not establish neuronal differentiation. Exact pathway identifiers, display labels and theme assignments are retained in source tables. **b,** The existing 12 smallest-overall-FDR directionally consistent CollecTRI regulator patterns; this panel is unchanged by the expanded GSEA screen and uses its own TF testing family and score scale.

GSEA uses signed moderated t rankings from the current overall batch + condition model and the August/May effects of the condition-by-period model. The primary ranking comprises 25,461 unique autosomal/non-PAR X symbols. Sets with 15-500 represented symbols were tested; BH correction spans all 23,696 eligible sets across every collection, separately per contrast. There were 3,281 overall discoveries, of which 3,264 had matching NES directions in both periods and 966 were also significant in each period. Positional, perturbation, disease and cell-type signatures describe their source definitions, not demonstrated mechanisms or identities in HT1080. Long labels are abbreviated for display; exact identifiers are retained in source tables. The three-collection CAMERA sensitivity does not constitute a correlation-aware test of this expanded family. GSEA remains exploratory.

Collection codes: H, Hallmark; C1, positional; C2, curated; C3, regulatory targets; C4, computational; C5, ontology; C6, oncogenic; C7, immunologic; C8, cell type; C9, cell perturbation.

## supp - Complete HT1080 all-collection GSEA significant-set atlas

**Figure:** `ch03_supp_ht1080_gsea_all_collections_atlas` | **Status:** 137-page supplementary atlas

[PDF](figures/chapter3/ch03_supp_ht1080_gsea_all_collections_atlas.pdf)

All 3,281 pathways passing overall global FDR <0.05 are shown, ordered by collection and overall FDR. Each row displays overall, August and May NES on a common colour scale across pages. Stars mark global FDR <0.05 in each contrast, and the overall FDR is printed beside the row. All overall discoveries are retained, including those without period-direction agreement; per-period statistics are therefore shown without a persistence filter. Full pathway identifiers, leading-edge genes, collection metadata, and page assignments are in the source tables. All tested pathways, including nonsignificant sets, are available in the searchable all-collection report and full CSV results.


## main - Fig 3.17. Transcriptional programmes associated with Y depletion in HT1080.

<a id="fig-3-17"></a>

**Figure:** `ch03_main_ht1080_combined_gsea` | **Chapter number:** Fig 3.17

**Working ID:** `CH3-S3.3.6.1-F03` | **Figure:** `ch03_main_ht1080_combined_gsea` | **Status:** selected; CAMERA FDR column added to panel a on 9 October 2026

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_F03_HT1080_GSEA_programmes/CH3_S3_4_6_1_F03_HT1080_GSEA_programmes.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_1_F03_HT1080_GSEA_programmes/CH3_S3_4_6_1_F03_HT1080_GSEA_programmes.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_1_F03_HT1080_GSEA_programmes/CH3_S3_4_6_1_F03_HT1080_GSEA_programmes_camera_247ccbcbe974.png)

**Fig 3.17. Transcriptional programmes associated with Y depletion in HT1080.**

**a,** Selected Hallmark, Reactome and GO Biological Process sets, with overall fgsea and CAMERA FDRs at right. **b,** Additional sets from the all-MSigDB analysis, grouped by theme, with overall fgsea FDRs. Colours show fgsea normalised enrichment scores for the overall and two sampling-period comparisons: blue indicates negative and red positive enrichment. Displayed sets met overall fgsea FDR < 0.05 with matching directions in both periods; stars indicate fgsea FDR < 0.05 in each comparison. CAMERA used fixed inter-gene correlation 0.01; bold values indicate CAMERA FDR < 0.05. Benjamini–Hochberg correction was applied separately within each method and comparison across 4,868 sets in a and 23,696 in b.

H, Hallmark; R, Reactome; GO, Gene Ontology Biological Process; C2, curated gene sets; C5, ontology gene sets.

[Source data and provenance](figures/chapter3/source_data/ht1080_combined_gsea/provenance.json)

## Retired - Separate HT1080 all-MSigDB GSEA figure

CH3-S3.4.6.1-S03 has been incorporated into CH3-S3.4.6.1-F03. Original figures and captions are retained in `archive/merge_ht1080_gsea_20260924_110638/`.


## main - HT1080 Metascape up pathway network

**Figure:** `ch03_main_ht1080_metascape_up_network` | **Status:** restyled existing 8 September 2026 combined-run Metascape analysis

[PDF](figures/chapter3/ch03_main_ht1080_metascape_up_network.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_metascape_up_network.svg) | [PNG](figures/chapter3/ch03_main_ht1080_metascape_up_network.png)

Network of enriched terms from the upregulated gene list in the existing combined August/May HT1080 DESeq2 batch + treatment analysis (BH FDR <0.05 and absolute log2FC >0.5; custom background). This analysis is separate from the current limma/GSEA workflow. Nodes represent 114 connected enriched terms and edges reproduce the 537 connections in the Metascape GONetwork export dated 2026-09-08 13:03:39. The existing Louvain community assignments are preserved, with curated descriptive theme names. Up to eight major communities are coloured; smaller communities are grey. Node area scales as 12 + 4 times negative log10 of the raw enrichment P value; table statistics are raw enrichment P values, not FDR. Network export log-q values remain in the original source files.

Layouts use Kamada-Kawai within connected components; separately packed components have no interpretable between-component distance. As in the source figure, 2 isolated terms are omitted. Tables give community term counts, strongest enrichment P value, two representative terms ordered by enrichment significance, and the four most frequently recurring hit genes (ties alphabetical). These are enrichment-hit lists, not GSEA leading-edge genes. Labels describing reference perturbations, chromatin marks or cell states denote source signatures and do not demonstrate the same perturbation or phenotype in HT1080. Full node descriptions, hits, communities and exported edges are saved in source_data/ht1080_metascape_latest.


## main - HT1080 Metascape down pathway network

**Figure:** `ch03_main_ht1080_metascape_down_network` | **Status:** restyled existing 8 September 2026 combined-run Metascape analysis

[PDF](figures/chapter3/ch03_main_ht1080_metascape_down_network.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_metascape_down_network.svg) | [PNG](figures/chapter3/ch03_main_ht1080_metascape_down_network.png)

Network of enriched terms from the downregulated gene list in the existing combined August/May HT1080 DESeq2 batch + treatment analysis (BH FDR <0.05 and absolute log2FC >0.5; custom background). This analysis is separate from the current limma/GSEA workflow. Nodes represent 49 connected enriched terms and edges reproduce the 203 connections in the Metascape GONetwork export dated 2026-09-08 13:07:22. The existing Louvain community assignments are preserved, with curated descriptive theme names. Up to eight major communities are coloured; smaller communities are grey. Node area scales as 12 + 4 times negative log10 of the raw enrichment P value; table statistics are raw enrichment P values, not FDR. Network export log-q values remain in the original source files.

Layouts use Kamada-Kawai within connected components; separately packed components have no interpretable between-component distance. As in the source figure, 10 isolated terms are omitted. Tables give community term counts, strongest enrichment P value, two representative terms ordered by enrichment significance, and the four most frequently recurring hit genes (ties alphabetical). These are enrichment-hit lists, not GSEA leading-edge genes. Labels describing reference perturbations, chromatin marks or cell states denote source signatures and do not demonstrate the same perturbation or phenotype in HT1080. Full node descriptions, hits, communities and exported edges are saved in source_data/ht1080_metascape_latest.


## supp - HT1080 Metascape networks, combined viewing copy

**Figure:** `ch03_supp_ht1080_metascape_networks_combined`

[PDF](figures/chapter3/ch03_supp_ht1080_metascape_networks_combined.pdf). Two-page viewing copy of the main up- and downregulated Metascape networks; captions and methodological definitions are given above.


## supp - Supp Fig 3.14. Pathway associations among genes with higher expression in Y-depleted HT1080 cells.

<a id="supp-fig-3-14"></a>

**Figure:** `ch03_main_ht1080_metascape_current_up_network_v2` | **Status:** aligned 19 September 2026 Metascape rerun | **Chapter number:** Supp Fig 3.14

[PDF](figures/chapter3/ch03_main_ht1080_metascape_current_up_network_v2.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_metascape_current_up_network_v2.svg) | [PNG](figures/chapter3/ch03_main_ht1080_metascape_current_up_network_v2.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_F04_HT1080_Metascape_pathway_networks/CH3_S3_4_6_1_F04A_HT1080_Metascape_higher.pdf)

**Supp Fig 3.14. Pathway associations among genes with higher expression in Y-depleted HT1080 cells.**

Metascape term-similarity network for genes with higher expression (FDR < 0.05, |log₂ fold change| ≥ 1 and matching directions in both sampling periods), using tested genes as background. Nodes represent enriched terms; node area scales with −log₁₀ enrichment *P* value, and connections indicate term similarity. Colours identify Louvain communities. The table gives community size, strongest enrichment, representative terms and recurring contributing genes. Network positions are for display and do not represent biological distances.

## supp - Supp Fig 3.15. Pathway associations among genes with lower expression in Y-depleted HT1080 cells.

<a id="supp-fig-3-15"></a>

**Figure:** `ch03_main_ht1080_metascape_current_down_network_v2` | **Status:** aligned 19 September 2026 Metascape rerun | **Chapter number:** Supp Fig 3.15

[PDF](figures/chapter3/ch03_main_ht1080_metascape_current_down_network_v2.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_metascape_current_down_network_v2.svg) | [PNG](figures/chapter3/ch03_main_ht1080_metascape_current_down_network_v2.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_F04_HT1080_Metascape_pathway_networks/CH3_S3_4_6_1_F04B_HT1080_Metascape_lower.pdf)

**Supp Fig 3.15. Pathway associations among genes with lower expression in Y-depleted HT1080 cells.**

Metascape term-similarity network for genes with lower expression (FDR < 0.05, |log₂ fold change| ≥ 1 and matching directions in both sampling periods), using tested genes as background. Nodes represent enriched terms; node area scales with −log₁₀ enrichment *P* value, and connections indicate term similarity. Colours identify Louvain communities. The table gives community size, strongest enrichment, representative terms and recurring contributing genes. Six isolated terms are omitted from the network and retained in the table. Network positions are for display and do not represent biological distances.

## supp - Current HT1080 Metascape networks, two-page copy

**Figure:** `ch03_supp_ht1080_metascape_current_networks`

[PDF](figures/chapter3/ch03_supp_ht1080_metascape_current_networks.pdf). Two-page viewing copy containing the aligned current UP and DOWN Metascape figures described above.

## Revised display - Current HT1080 Metascape networks (version 2)

**Working IDs:** `CH3-S3.4.6.1-F04A` (higher expression) and `CH3-S3.4.6.1-F04B` (lower expression) | **Status:** selected companion main figures

[Higher-expression network](figures/chapter3/ch03_main_ht1080_metascape_current_up_network_v2.pdf) | [Lower-expression network](figures/chapter3/ch03_main_ht1080_metascape_current_down_network_v2.pdf) | [Two-page supplementary copy](figures/chapter3/ch03_supp_ht1080_metascape_current_networks_v2.pdf)

The gene selection, enrichment results, community membership, colours and node-size encoding are as described above. The revised UP display contains all 124 terms and 772 connections. The revised DOWN display contains 70 terms and all 318 unique connections; six isolated grey reference-signature terms are omitted from the network display and retained in the source tables and summary table. Communities are positioned using fixed display anchors, with seeded spring layouts within each community. These positions improve legibility and do not represent biological coordinates or distances. Larger theme labels are placed outside the node groups. Earlier figure versions are archived under `figures/chapter3/archive/metascape_layout_before_revision_20260919`.

## supp - Supp Fig 3.13. Enriched pathways and recurring contributing genes in HT1080.

<a id="supp-fig-3-13"></a>

**Figure:** `ch03_supp_ht1080_metascape_terms_and_gene_modules` | **Chapter number:** Supp Fig 3.13

**Working ID:** `CH3-S3.4.6.1-S04` | **Figure:** `ch03_supp_ht1080_metascape_terms_and_gene_modules` | **Status:** selected supplementary figure

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S04_HT1080_Metascape_terms_and_gene_modules/CH3_S3_4_6_1_S04_HT1080_Metascape_terms_and_gene_modules.pdf) | [Working SVG](figures/chapter3/ch03_supp_ht1080_metascape_terms_and_gene_modules.svg) | [Working PNG](figures/chapter3/ch03_supp_ht1080_metascape_terms_and_gene_modules.png)

**Supp Fig 3.13. Enriched pathways and recurring contributing genes in HT1080.**

**a,b,** Top 20 Metascape terms for genes with higher (a) or lower (b) expression, ranked by enrichment *P* value. Bars show −log₁₀ *P*; labels give contributing-gene counts and reported *q* values. **c,** Genes recurring most frequently across higher-expression terms. Bar length indicates the number of terms containing each gene; colours identify pathway communities from **Supp Fig** 3.14. Rows are grouped by the community contributing the most memberships; daggers mark ties.

Input genes met FDR < 0.05 and |log₂ fold change| ≥ 1 with matching directions in both sampling periods. Enrichment used the tested-gene background.

[Source data and provenance](figures/chapter3/source_data/ht1080_metascape_terms_and_gene_modules/provenance.json)

## supp - Supp Fig 3.16. Protein-association modules organise HT1080 expression changes.

<a id="supp-fig-3-16"></a>

**Figure:** `ch03_main_ht1080_expanded_string_modules` | **Chapter number:** Supp Fig 3.16

**Working ID:** `CH3-S3.3.6.1-F05` | **Figure:** `ch03_main_ht1080_expanded_string_modules` | **Status:** selected; module labels reviewed 9 October 2026

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_F05_HT1080_STRING_modules/CH3_S3_4_6_1_F05_HT1080_STRING_modules.pdf) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_1_F05_HT1080_STRING_modules/CH3_S3_4_6_1_F05_HT1080_STRING_modules_clean_2f23c0ddeb11.png)

**Supp Fig 3.16. Protein-association modules organise HT1080 expression changes.**

**a,b,** Selected modules among genes with lower (a) or higher (b) expression. Candidates had |mean log₂ expression difference| ≥ 0.5, matching directions between periods and stability across leave-one-match-out checks, without a gene-level significance threshold. STRING associations (score ≥ 0.7) were clustered using MCODE. Node area reflects absolute expression difference and edge width reflects association score. Solid edges have physical-network support ≥ 0.7; dashed edges indicate other STRING associations. Labels give module size and the top annotation’s Benjamini–Hochberg-adjusted *P* value across all module–term tests. These STRING-connected groups may span several functions; labels do not establish a single pathway or an assembled complex in HT1080.

## supp - Supp Fig 3.17. Functional annotations of HT1080 protein-association modules.

<a id="supp-fig-3-17"></a>

**Figure:** `ch03_supp_ht1080_module_annotation_significance` | **Chapter number:** Supp Fig 3.17

**Working ID:** `CH3-S3.3.6.1-S17` | **Figure:** `ch03_supp_ht1080_module_annotation_significance` | **Status:** selected; top terms added 9 October 2026

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S17_HT1080_module_annotations/CH3_S3_4_6_1_S17_HT1080_module_annotations.pdf) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_1_S17_HT1080_module_annotations/CH3_S3_4_6_1_S17_HT1080_module_annotations_clean_4ba47da93531.png)

**Supp Fig 3.17. Functional annotations of HT1080 protein-association modules.**

**a,b,** Top annotations of lower-expression (a) and higher-expression (b) modules. Each descriptive module name is followed by its top GO Biological Process (GO:BP) or Reactome term. Circles show one-sided hypergeometric *P* values and diamonds show Benjamini–Hochberg FDRs across all module–term tests, using expressed STRING-mapped genes as background. Statistics refer to the displayed term. The dashed line marks 0.05; Genes gives module size. Module names match **Supp Fig** 3.16; enrichment tests do not test expression changes.

## main - Fig 3.29. UTY and KDM5D re-expression in Y-depleted HT1080 cells.

<a id="fig-3-29"></a>

**Figure:** `ch03_main_ht1080_addback_design_expression` | **Working ID:** CH3-S3.3.5-F01 | **Results:** 3.3.5 | **Status:** revised 8 October 2026: sequential UTY-to-combined derivation (a), obsolete blot sub-labels covered (b), supplementary raw-count panel incorporated (c), saved plasmid-only PCA retained (d), residual Y-expression heatmap incorporated (e) | **Chapter number:** Fig 3.29

[Selected PDF](figures/chapter3/selected/CH3_S3_3_5_F01_ht1080_addback_design_expression/ch03_main_ht1080_addback_design_expression.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_3_5_F01_ht1080_addback_design_expression/ch03_main_ht1080_addback_design_expression.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_3_5_F01_ht1080_addback_design_expression/ch03_main_ht1080_addback_design_expression.png)

[PDF](figures/chapter3/ch03_main_ht1080_addback_design_expression.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_addback_design_expression.svg) | [PNG](figures/chapter3/ch03_main_ht1080_addback_design_expression.png)

**Fig 3.29. UTY and KDM5D re-expression in Y-depleted HT1080 cells.**

**a,** Derivation of the six populations. The combined population was generated by introducing KDM5D into the established UTY-only population. Symbols indicate three matched RNA collection occasions. **b,** Western blots for KDM5D (top) and UTY (bottom), with α-actinin loading controls. LOY denotes the untransduced Y-depleted population and pLX_HA the empty vector. **c,** Raw reads assigned to UTY and KDM5D, including transgene-derived reads, across collection occasions 1–3. Outlines mark populations carrying the corresponding construct; the right-hand column gives total assigned reads in millions. **d,** PCA of the four plasmid-bearing populations using expressed autosomal and non-PAR X-linked genes. Colours identify populations and shapes identify collection occasions. Ellipses summarise within-population spread and are not confidence regions. **e,** Mean log₂ CPM for seven Y-linked genes across three collection occasions per population, using the primary effective library sizes. Numbers give group means; UTY and KDM5D include transgene-derived reads.

## supp - Scale of UTY and KDM5D overexpression in the HT1080 add-back populations

**Figure:** `ch03_supp_ht1080_addback_construct_overexpression` | **Status:** removed from Results_edit.md at the author’s request on 8 October 2026. Panel a is retained as panel c of the main add-back design figure; this standalone caption and its exports are historical records.

[PDF](figures/chapter3/ch03_supp_ht1080_addback_construct_overexpression.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_construct_overexpression.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_construct_overexpression.png)

**UTY and KDM5D abundance following re-expression.**

**a,** Reads assigned to UTY and KDM5D across three matched collection occasions per population. Outlined tiles identify populations carrying the corresponding construct; the right-hand column shows total assigned reads per library. **b,** Expression distributions across 23,994 genes, pooling three libraries per population. Coloured points show UTY and KDM5D abundance; shapes identify collection occasions. Labels give mean within-library expression ranks and the percentage of genes at or above that abundance. Counts include endogenous and transgene-derived transcripts.


## supp - Supp Fig 3.39. Global expression structure of HT1080 add-back populations.

<a id="supp-fig-3-39"></a>

**Figure:** `ch03_supp_ht1080_addback_pca_populations` | **Working ID:** CH3-S3.3.5-S01 | **Results:** 3.3.5 | **Status:** 22 September 2026; replaces `ch03_supp_ht1080_addback_pca_plasmid_only`, whose content moved to main panel d (archived) | **Chapter number:** Supp Fig 3.39

[Selected PDF](figures/chapter3/selected/CH3_S3_3_5_S01_ht1080_addback_pca_populations/ch03_supp_ht1080_addback_pca_populations.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_3_5_S01_ht1080_addback_pca_populations/ch03_supp_ht1080_addback_pca_populations.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_3_5_S01_ht1080_addback_pca_populations/ch03_supp_ht1080_addback_pca_populations.png)

[PDF](figures/chapter3/ch03_supp_ht1080_addback_pca_populations.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_pca_populations.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_pca_populations.png)

**Supp Fig 3.39. Global expression structure of HT1080 add-back populations.**

PCA of 18 libraries using TMM-normalised log₂ CPM for 23,992 expressed autosomal and non-PAR X-linked genes. Colours identify populations and shapes identify three matched collection occasions per established population. Ellipses summarise within-population spread and are not confidence regions. Axis labels indicate variance explained.

## main - Fig 3.30. UTY-containing populations show partial recovery of the HT1080 LOY expression signature.

<a id="fig-3-30"></a>

**Figure:** `ch03_main_ht1080_addback_signature_recovery` | **Working ID:** CH3-S3.3.5-F02 | **Results:** 3.3.5 | **Status:** Three-panel main figure updated 9 October 2026; panel c is the 150-gene directional-test panel formerly S02b. | **Chapter number:** Fig 3.30

[Selected PDF](figures/chapter3/selected/CH3_S3_3_5_F02_ht1080_addback_signature_recovery/ch03_main_ht1080_addback_signature_recovery.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_3_5_F02_ht1080_addback_signature_recovery/ch03_main_ht1080_addback_signature_recovery.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_3_5_F02_ht1080_addback_signature_recovery/ch03_main_ht1080_addback_signature_recovery.png)

[PDF](figures/chapter3/ch03_main_ht1080_addback_signature_recovery.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_addback_signature_recovery.svg) | [PNG](figures/chapter3/ch03_main_ht1080_addback_signature_recovery.png)

**Fig 3.30. UTY-containing populations show partial recovery of the HT1080 LOY expression signature.**

**a,** Scores for the August-defined 150-gene signature (75 higher- and 75 lower-in-LOY genes). Higher scores indicate a more LOY-like profile. Points show three matched collection occasions per established population; shapes identify occasions and bars show adjusted means. Shading marks the empty-vector-to-sgCTRL gap. **b,** Percentage of this score gap recovered, with vector at 0% and sgCTRL at 100%. Points and bars show estimates and 95% intervals from the collection-adjusted model; *P* values compare each population with vector. **c,** Directional fry tests of the 75 higher- and 75 lower-in-LOY genes across six contrasts. Headers give the expected reversal direction. Cells show observed direction and BH FDR across 12 tests; blue denotes decrease, orange increase and bold text FDR < 0.05.

## supp - Supp Fig 3.40. HT1080 add-back expression responses across signature and promoter-overlap genes.

<a id="supp-fig-3-40"></a>

**Figure:** `ch03_supp_ht1080_addback_baseline_reversal` | **Working ID:** CH3-S3.3.5-S02 | **Results:** 3.3.5 | **Status:** Four-panel supplement updated 9 October 2026: a, 150-gene heatmap; b–c, frozen-100 scores and contrasts; d, former S04 39-gene RNA/ATAC overlap. The 150-gene tests are now F02c. Frozen-100 directional tests remain omitted. | **Chapter number:** Supp Fig 3.40

[Selected PDF](figures/chapter3/selected/CH3_S3_3_5_S02_ht1080_addback_baseline_reversal/ch03_supp_ht1080_addback_baseline_reversal.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_3_5_S02_ht1080_addback_baseline_reversal/ch03_supp_ht1080_addback_baseline_reversal.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_3_5_S02_ht1080_addback_baseline_reversal/ch03_supp_ht1080_addback_baseline_reversal.png)

[PDF](figures/chapter3/ch03_supp_ht1080_addback_baseline_reversal.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_baseline_reversal.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_baseline_reversal.png)

**Supp Fig 3.40. HT1080 add-back expression responses across signature and promoter-overlap genes.**

**a,** Group-mean log₂ differences from sgCTRL for 150 August-defined genes, ordered by direction/effect and clipped to ±4. **b,** Frozen 100-gene signature scores (50 higher/50 lower in LOY); higher scores are more LOY-like. Shapes distinguish three matched collection occasions; bars show means. **c,** Collection-adjusted score differences (first population minus second), with 95% intervals and nominal two-sided moderated-t *P*. **d,** Effects for 39 genes with baseline RNA/promoter loss (both FDR < 0.05; RNA |log₂ fold change| ≥ 1 and matching period directions). Baseline Y-depleted cells are compared with sgCTRL; other populations with empty vector. Red/blue denote higher/lower expression; dots mark within-contrast gene-level FDR < 0.05. Each group is one established population sampled three times; August/May reuse the original baseline populations. Accessibility was not measured after add-back.

## supp - Derivation of the 150-gene HT1080 LOY signature and its score

**Figure:** `ch03_supp_ht1080_addback_signature_schematic` | **Status:** removed from Results_edit.md at the author’s request on 8 October 2026. Caption, exports and source data are retained as historical records; the signature-derivation prose remains in Results.

[PDF](figures/chapter3/ch03_supp_ht1080_addback_signature_schematic.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_signature_schematic.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_signature_schematic.png)

Schematic of signature selection and scoring. Left: gene selection used only the August 2025 sgCTRL and Y-depleted libraries (three passage-matched pairs), with chromosome-Y and pseudoautosomal genes excluded. Of 23,903 background genes, 10,113 had the same sign of matched difference in all three pairs and 3,738 also had a mean absolute difference of at least 0.5 log2 units. Of these, 218 (75 higher and 143 lower in Y-depleted cells) were supported in the same direction by an August-only limma-voom fit at FDR < 0.20. The 75 largest in each direction formed the 150-gene signature. Right: each May library is scored from its own raw counts by converting the 23,903 background genes to within-library percentile ranks and subtracting the mean percentile of the 75 lower-in-LOY genes from that of the 75 higher-in-LOY genes; higher scores are more Y-depleted-like. This follows the single-sample rank-based approach of singscore (Foroutan et al., 2018), computed directly on the frozen background. Percent recovered expresses a population's change from the empty-vector mean as a fraction of the vector-to-sgCTRL gap.

## supp - Sensitivity of LOY-signature recovery to signature size and an FDR rule

**Figure:** `ch03_supp_ht1080_addback_signature_variants` | **Status:** removed from Results_edit.md at the author’s request on 8 October 2026. Caption, exports and source data are retained as historical records; the signature-sensitivity prose remains in Results.

[PDF](figures/chapter3/ch03_supp_ht1080_addback_signature_variants.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_signature_variants.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_signature_variants.png)

Four signatures built from the same August 2025 background with the same concordance and effect-size rules (same sign in all three passage pairs; absolute mean difference ≥ 0.5 log2): the originally frozen 100 genes (50/50, no FDR rule), 200 genes without an FDR rule (100/100), 100 genes with an August-only limma-voom FDR < 0.20 (50/50), and the balanced 150 genes with that FDR rule (75/75; only 75 higher-in-LOY genes pass). **a,** LOY signature score per library under each signature. **b,** Percentage of the vector-to-sgCTRL gap recovered, with 95% intervals from multivariate-t draws of the blocked score model. Across the four signatures, UTY recovered 45–57% of the gap and UTY + KDM5D 35–41%. KDM5D remained between −15% and −13%, and the untransduced Y-depleted population moved further from sgCTRL under every variant. The FDR-supported variants share 42 of the original 100 genes. An August-only FDR < 0.05 leaves only 26 rule-passing genes and cannot support a 100-gene signature.

## supp - Movement of the originally frozen 100-gene HT1080 LOY signature after add-back

**Figure:** `ch03_supp_ht1080_addback_baseline_reversal_frozen100` | **Working ID:** CH3-S3.3.5-S03 | **Results:** 3.3.5 | **Status:** removed from active Results on 9 October 2026. Panels a–b are retained in S02b–c; the second heatmap and frozen-100 directional-test panel are omitted from active display. This caption and all S03 exports remain historical records.

[Selected PDF](figures/chapter3/selected/CH3_S3_3_5_S03_ht1080_addback_baseline_reversal_frozen100/ch03_supp_ht1080_addback_baseline_reversal_frozen100.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_3_5_S03_ht1080_addback_baseline_reversal_frozen100/ch03_supp_ht1080_addback_baseline_reversal_frozen100.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_3_5_S03_ht1080_addback_baseline_reversal_frozen100/ch03_supp_ht1080_addback_baseline_reversal_frozen100.png)

[PDF](figures/chapter3/ch03_supp_ht1080_addback_baseline_reversal_frozen100.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_baseline_reversal_frozen100.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_baseline_reversal_frozen100.png)

**Movement of the originally frozen 100-gene HT1080 LOY signature after add-back.**

**a,** Within-sample scores for the originally frozen August signature (50 higher- and 50 lower-in-LOY genes). Points show three matched collection occasions per population, distinguished by shape; bars show means. **b,** Collection-adjusted score differences and 95% intervals; P values are nominal. **c,** Mean log₂ expression differences from sgCTRL for the 75 signature genes passing the May filter, clipped to −4 to 4. **d,** Directional fry tests for those genes (38 higher, 37 lower in LOY), with BH FDR across 12 tests. Blue denotes decrease, orange increase and bold text FDR <0.05. Headers give the expected reversal directions; LXHA denotes empty vector. Each group represents one established population.

## supp - Sensitivity and residual Y-expression context for HT1080 add-back RNA-seq

**Figure:** `ch03_supp_ht1080_addback_sensitivity` | **Status:** removed from Results_edit.md at the author’s request on 8 October 2026. Panel c is retained as panel e of the main add-back design figure; this standalone caption and exports are historical records.

[PDF](figures/chapter3/ch03_supp_ht1080_addback_sensitivity.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_sensitivity.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_sensitivity.png)

Sensitivity analyses and Y-linked expression context. **a,** Numbers of genes selected using the collection-occasion-blocked nonzero-effect test followed by an absolute log2 fold-change filter of 0.5, blocked TREAT at 0.5 and an unblocked model with the nonzero-effect filter. FDR correction was applied within each gene-level contrast. **b,** Pearson correlations between full-model and leave-one-collection-occasion-out log2 fold changes for six comparisons. The correlations describe effect-estimate agreement; selected-gene counts varied more strongly and are retained in the source table. **c,** Group-mean log2 CPM for selected Y-linked genes, calculated with the primary effective library sizes. UTY and KDM5D include transgene-derived reads in the corresponding construct populations. Regional residual-Y differences provide population context and are not part of the primary autosomal/non-PAR X differential-expression background.

## supp - Correlation sensitivity and construct specificity of add-back pathway tests

**Figure:** `ch03_supp_ht1080_addback_pathway_sensitivity` | **Status:** removed from Results_edit.md with its associated pathway-analysis prose at the author’s request on 8 October 2026. Caption, exports and source data are retained as historical records.

[PDF](figures/chapter3/ch03_supp_ht1080_addback_pathway_sensitivity.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_pathway_sensitivity.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_pathway_sensitivity.png)

Competitive CAMERA gene-set tests for the HT1080 add-back populations. **a,** Numbers of MSigDB Hallmark (50 tested), Reactome (1,340) and GO Biological Process (4,684) sets with FDR <0.05 for six contrasts using the default fixed inter-gene correlation of 0.01. Benjamini-Hochberg correction was applied within each collection and contrast. **b,** Signed CAMERA z statistics for the 50 Hallmark sets, comparing UTY versus empty vector with untransduced Y-depleted versus the same empty vector (fixed correlation 0.01). The dashed line marks identity. Colour shows which comparisons reached FDR <0.05; the eight sets with the largest UTY statistics are labelled. **c,** Spearman correlations of signed CAMERA statistics between each add-back-versus-vector contrast and naive Y-depleted versus vector, by collection (fixed correlation 0.01). Labels give the number of sets significant for the add-back contrast that were also significant, in the same direction, for naive Y-depleted versus vector. None were significant in the opposite direction for UTY or combined re-expression. Each group is one established population sampled on three matched collection occasions; pathway results are conditional on these populations. Complete per-set results are in the source tables.

## supp - Expression responses of genes with baseline RNA and promoter-accessibility loss

**Figure:** `ch03_supp_ht1080_addback_atac_overlap` | **Working ID:** CH3-S3.3.5-S04 | **Results:** 3.3.5 | **Status:** removed from active Results on 9 October 2026 and incorporated into S02d. All 39 genes, five contrasts, original colour limits and FDR markers are retained. This standalone caption and exports remain historical.

[Selected PDF](figures/chapter3/selected/CH3_S3_3_5_S04_ht1080_addback_atac_overlap/ch03_supp_ht1080_addback_atac_overlap.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_3_5_S04_ht1080_addback_atac_overlap/ch03_supp_ht1080_addback_atac_overlap.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_3_5_S04_ht1080_addback_atac_overlap/ch03_supp_ht1080_addback_atac_overlap.png)

[PDF](figures/chapter3/ch03_supp_ht1080_addback_atac_overlap.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_addback_atac_overlap.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_addback_atac_overlap.png)

**Expression responses of genes with baseline RNA and promoter-accessibility loss.**

Expression differences for 39 genes with baseline RNA and promoter-accessibility loss (FDR <0.05 in both assays; RNA |log₂ fold change| ≥1 with matching period directions). Columns compare baseline Y-depleted cells with sgCTRL, and the indicated add-back populations with empty vector. Red and blue indicate higher and lower expression; dots mark gene-level FDR <0.05 within each add-back contrast. Each established population was sampled on three matched occasions. Accessibility was not measured after add-back.

## main - Fig 3.18. Y depletion is associated with widespread reductions in chromatin accessibility in HT1080 cells.

<a id="fig-3-18"></a>

**Figure:** `ch03_main_ht1080_atac_accessibility` | **Status:** draft for review (22 September 2026) | **Chapter number:** Fig 3.18

[PDF](figures/chapter3/ch03_main_ht1080_atac_accessibility.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_atac_accessibility.svg) | [PNG](figures/chapter3/ch03_main_ht1080_atac_accessibility.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_F06_HT1080_ATAC_accessibility/CH3_S3_4_6_1_F06_HT1080_ATAC_accessibility.pdf)

**Fig 3.18. Y depletion is associated with widespread reductions in chromatin accessibility in HT1080 cells.**

One sgCTRL and one Y-depleted population sampled at passages 11, 12 and 21. **a,** Y fragments per million non-Y fragments. **b,** PCA across 84,510 consensus regions excluding Y/PAR; arrows connect matching passages. **c,** Differential accessibility (edgeR, FDR < 0.05; dashed threshold). Open circles mark higher-confidence changes passing edgeR’s glmTreat minimum-effect test (|log₂ fold change| > 0.5; FDR < 0.05); labels identify nearest genes, with promoters in bold. **d,** Proportions of peaks assigned to each genomic feature; promoters are TSS ±2 kb. Comparisons use two-sided Fisher tests with BH correction. **e,** Significant-region fold changes; dashed lines mark twofold changes. **f,** Promoter coverage for genes of interest in 50-bp bins per million non-Y fragments. Thin lines show samples; shading shows condition means. Grey bands mark tested regions; arrows show transcription direction. Blue/red denote sgCTRL/Y-depleted in a,b,f and accessibility losses/gains in c,e.

## supp - Supp Fig 3.4. HT1080 ATAC-seq library quality.

**Figure:** `ch03_methods_ht1080_atac_qc` | **Chapter number:** Supp Fig 3.4 | **Placement:** Methods | **Status:** moved from Results supplement CH3-S3.3.6.1-S05 to Methods section 3.2.4.3 on 9 October 2026; panels unchanged

[PDF](figures/chapter3/ch03_methods_ht1080_atac_qc.pdf) | [SVG](figures/chapter3/ch03_methods_ht1080_atac_qc.svg) | [PNG](figures/chapter3/ch03_methods_ht1080_atac_qc.png)

**Supp Fig 3.4. HT1080 ATAC-seq library quality.**

Six libraries from one sgCTRL and one Y-depleted population, each sampled at passages 11, 12 and 21. **a,** Fragment-length distributions. **b,** Tn5 insertion profiles around transcription start sites (TSSs), normalised to flanks; enrichment compares the central 100 bp with outer flanks. **c,** Fragment counts, fraction overlapping each library's own peaks (FRiP), and peak counts before/after downsampling to 13.6 million fragments. **d,** Spearman correlations across 84,510 regions, excluding Y and PAR. **e,** Coverage around 34,174 TSSs per million non-Y fragments, with shared row order and colour scale capped at the 99th percentile. Blue/red denote sgCTRL (C)/Y-depleted (L); solid/dashed/dotted lines denote passages 11/12/21. Passage samples are repeated measurements of established populations.

[Source and reproduction notes](figures/chapter3/ch03_methods_ht1080_atac_qc_notes.md)

## supp - Fully labelled HT1080 differential-accessibility volcano

**Figure:** `ch03_supp_ht1080_atac_volcano_labelled` | **Status:** draft for review (22 September 2026)

[PDF](figures/chapter3/ch03_supp_ht1080_atac_volcano_labelled.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_atac_volcano_labelled.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_atac_volcano_labelled.png)

The differential-accessibility result of the main figure (panel c), shown with a broader label set that includes LOC-annotated loci. The data, test and display match the main panel: edgeR quasi-likelihood F-test, design ~ condition, 84,510 tested regions. Positive log2 fold changes indicate higher accessibility in Y-depleted cells. The dashed line marks the largest raw P value with Benjamini–Hochberg FDR < 0.05, and open circles mark the 32 regions that also passed TREAT at |log2FC| > 0.5. Labels follow a fixed rule applied to all FDR < 0.05 regions: the 20 lowest P values, the five most extreme fold changes in each direction, and the five loci shown in the main figure. Each label gives the gene with the nearest knownGene TSS. Bold labels mark promoter-overlapping regions; for other regions the gene name is positional only. Several of the strongest changes lie near uncharacterised LOC loci, which is why they are not labelled in the main figure. Labels for less extreme lost regions are placed to the right of the lost cloud for legibility.

## supp - Supp Fig 3.20. Accessibility and motif content at HT1080 differential regions.

<a id="supp-fig-3-20"></a>

**Figure:** `ch03_supp_ht1080_atac_differential_heatmaps` | **Status:** draft for review (22 September 2026) | **Chapter number:** Supp Fig 3.20

[PDF](figures/chapter3/ch03_supp_ht1080_atac_differential_heatmaps.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_atac_differential_heatmaps.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_atac_differential_heatmaps.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S08_HT1080_ATAC_differential_regions/CH3_S3_4_6_1_S08_HT1080_ATAC_differential_regions.pdf)

**Supp Fig 3.20. Accessibility and motif content at HT1080 differential regions.**

Differential regions (ATAC FDR < 0.05) were compared with three non-differential controls per region, matched for GC content and baseline accessibility. **a,** Mean fragment coverage around region centres; shading spans three passage samples per condition. Values give Y-depleted/sgCTRL coverage ratios within the central 200 bp. **b,** Enrichment of selected motifs, shown as Mantel–Haenszel odds ratios with 95% confidence intervals, stratified by matching group. Filled points indicate Benjamini–Hochberg FDR < 0.05. **c,** Individual-sample coverage, log₂ ratios of condition means and motif presence. Percentages indicate regions containing each motif; bold values mark significant target–control differences. Rows share the same order across columns; block heights are not proportional to region counts.

Coverage is normalised per million non-Y fragments. Blue denotes sgCTRL (C) and red Y-depleted (L) in a; blue and red identify lost and gained regions in b.

## supp - Supp Fig 3.25. Reduced promoter accessibility extends beyond significant peak calls.

<a id="supp-fig-3-25"></a>

**Figure:** `ch03_supp_ht1080_atac_tss_groups` | **Status:** draft for review (23 September 2026) | **Chapter number:** Supp Fig 3.25

[PDF](figures/chapter3/ch03_supp_ht1080_atac_tss_groups.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_atac_tss_groups.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_atac_tss_groups.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S16_HT1080_ATAC_TSS_groups/CH3_S3_4_6_1_S16_HT1080_ATAC_TSS_groups.pdf)

**Supp Fig 3.25. Reduced promoter accessibility extends beyond significant peak calls.**

TSS-centred ATAC coverage for the indicated gene groups and matched non-differentially-expressed controls (three per gene). Promoter-loss groups use ATAC FDR < 0.05. Coverage is normalised per million non-Y fragments in 50-bp bins. Mean profiles are shown in **Fig** 3.21d.

**a,** TSS-centred coverage and log₂ Y-depleted/sgCTRL ratios, with genes above their matched controls. Rows are ordered by pooled central coverage.

**b,** Per-gene central log₂ coverage ratios within TSS ±100 bp. Boxes show medians and quartiles; whiskers extend to 1.5 × IQR. Δ summarises the median difference between genes and controls across 200 control draws, with the 2.5th–97.5th percentile range. *P* values are median two-sided Mann–Whitney values across those draws.

**c,** Central coverage ratios for genes (circles) and controls (squares), with one point per passage pair.

## main - Fig 3.21. Reduced gene expression accompanies promoter closing in Y-depleted HT1080 cells.

<a id="fig-3-21"></a>

**Figure:** `ch03_main_ht1080_rna_atac_integration` | **Status:** draft for review (22 September 2026) | **Chapter number:** Fig 3.21

[PDF](figures/chapter3/ch03_main_ht1080_rna_atac_integration.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_rna_atac_integration.svg) | [PNG](figures/chapter3/ch03_main_ht1080_rna_atac_integration.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_F09_HT1080_RNA_ATAC_integration/CH3_S3_4_6_1_F09_HT1080_RNA_ATAC_integration.pdf)

**Fig 3.21. Reduced gene expression accompanies promoter closing in Y-depleted HT1080 cells.**

RNA criteria: FDR < 0.05, |log₂ fold change| ≥ 1 and matching period directions. ATAC FDR < 0.10, except promoter-loss groups in d (< 0.05); promoters are TSS ±2 kb. **a,b,** Tested-promoter coverage among 662 RNA hits (a) and promoter losses/gains by RNA class (b). **c,** Effects across 17,032 promoter peak–gene links. Blue/red mark concordant decreases/increases; diamonds mark opposite changes. Lighter points are significant in ≤ 1 assay; Spearman correlations describe all/jointly significant links. **d,** Mean TSS coverage for indicated genes (solid) and matched non-DE controls (dashed; three per gene). Shading spans three passage samples per condition; values give Y-depleted/sgCTRL coverage within TSS ±100 bp. **e,** Representative loci: ATAC lines show samples and shading condition means; peak bands mark losses/gains/nonsignificant regions (blue/red/grey). RNA circles/triangles denote periods; horizontal marks show means. Coverage uses 50-bp bins per million non-Y fragments. In d,e, population colours are blue (sgCTRL) and red (Y-depleted).

## supp - Supp Fig 3.24. RNA–ATAC concordance differs by peak location.

<a id="supp-fig-3-24"></a>

**Figure:** `ch03_supp_ht1080_rna_atac_peak_scatter` | **Status:** draft for review (22 September 2026) | **Chapter number:** Supp Fig 3.24

[PDF](figures/chapter3/ch03_supp_ht1080_rna_atac_peak_scatter.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_rna_atac_peak_scatter.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_rna_atac_peak_scatter.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S12_HT1080_RNA_ATAC_peak_types/CH3_S3_4_6_1_S12_HT1080_RNA_ATAC_peak_types.pdf)

**Supp Fig 3.24. RNA–ATAC concordance differs by peak location.**

RNA and ATAC log₂ fold changes for **a,** promoter, **b,** gene-body and **c,** intergenic peaks. Promoter peaks overlap TSS ±2 kb and retain every overlapping gene assignment; other peaks are assigned to the nearest annotated TSS. Each point represents a peak–gene link.

RNA differential expression required FDR < 0.05, |log₂ fold change| ≥ 1 and matching directions in both sampling periods; ATAC significance required FDR < 0.10. Blue/red indicates concordant decreases/increases; diamonds indicate opposite changes. Lighter points denote significance in one assay or neither. Panels report link and gene counts and Spearman correlations.

## supp - Supp Fig 3.26. Gene-set enrichment across expression and promoter-accessibility groups.

<a id="supp-fig-3-26"></a>

**Figure:** `ch03_supp_ht1080_quadrant_ORA` | **Status:** draft for review, exploratory (22 September 2026) | **Chapter number:** Supp Fig 3.26

[PDF](figures/chapter3/ch03_supp_ht1080_quadrant_ORA.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_quadrant_ORA.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_quadrant_ORA.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S13_HT1080_RNA_ATAC_quadrant_ORA/CH3_S3_4_6_1_S13_HT1080_RNA_ATAC_quadrant_ORA.pdf)

**Supp Fig 3.26. Gene-set enrichment across expression and promoter-accessibility groups.**

**a,** RNA and median promoter ATAC log₂ fold changes for 13,300 genes. Dashed and dotted lines mark ±0.25 and ±0.5, respectively. Quadrant counts use these thresholds in both assays, with > 0.5 counts in brackets; no significance filter was applied. **b,** Gene-set over-representation among genes with concordant changes, showing up to four significant sets per collection. Dot size indicates gene count and colour indicates FDR. Black outlines identify sets remaining significant against a background matched for RNA direction and effect threshold. **c,** Percentage of embryonic stem-cell H3K27me3 targets within each RNA–ATAC group; labels give members/total genes. Flat denotes |log₂ fold change| ≤ 0.25. **d,** H3K27me3- and PRC2-target enrichment in the indicated promoter groups. Controls were matched for baseline accessibility, GC content and CpG content, except in the unmatched comparison. Points show Mantel–Haenszel odds ratios with 95% confidence intervals; open points indicate FDR ≥ 0.05. Vertical ticks show estimates using one control per target instead of three.

## supp - Supp Fig 3.21. Accessibility, footprinting and expression evidence for selected regulatory motifs.

<a id="supp-fig-3-21"></a>

**Figure:** `ch03_supp_ht1080_atac_motif_open_closed` | **Status:** draft for review (22 September 2026) | **Chapter number:** Supp Fig 3.21

[PDF](figures/chapter3/ch03_supp_ht1080_atac_motif_open_closed.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_atac_motif_open_closed.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_atac_motif_open_closed.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S10_HT1080_ATAC_motif_evidence/CH3_S3_4_6_1_S10_HT1080_ATAC_motif_evidence.pdf)

**Supp Fig 3.21. Accessibility, footprinting and expression evidence for selected regulatory motifs.**

**a,** Selected motifs across HOMER enrichment, chromVAR deviations, TOBIAS footprint changes and TF RNA effects (combined, August and May). HOMER blue/red denotes enrichment in closing/opening regions (log₂ enrichment, 0–1.5; ≤ 0 white, > 1.5 saturated), with BH correction across 879 motifs in both directions. For other assays, blue/red denotes lower/higher values in Y-depleted cells; RNA scales span ±3 log₂ fold change. Black dots mark FDR < 0.05 for HOMER, chromVAR and RNA; rings additionally mark RNA |log₂ fold change| ≥ 1. Grey RNA tiles are unavailable. RNA retains its original testing families; nonsignificant estimates remain coloured. **b,** Per-sample chromVAR scores, centred within motifs. **c,** TOBIAS changes and pooled-site *P* values across all motifs; coloured points identify a. TOBIAS pools conditions, does not estimate between-sample variability and has no significance markers in a. C, sgCTRL; L, Y-depleted.

## exploratory - Motifs enriched at closing and opening chromatin (fold-change-filtered), with sequence logos

**Figure:** `ch03_exploratory_ht1080_atac_motif_logos_fold_change_filtered` | **Status:** exploratory. The author selected the all-significant version (`ch03_supp_ht1080_atac_motif_logos_all_significant`) for the supplement on 22 September 2026.

[PDF](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_fold_change_filtered.pdf) | [SVG](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_fold_change_filtered.svg) | [PNG](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_fold_change_filtered.png)

HOMER known-motif enrichment (JASPAR2024 CORE vertebrates, 879 motifs) in regions that close or open in Y-depleted HT1080 cells. Targets are the primary-analysis regions with FDR < 0.05 and |log2FC| > 1 (609 closing, 47 opening). Each set was compared with GC- and accessibility-matched non-differential regions (FDR ≥ 0.20, 10 per target). The rows are the ten lowest-FDR motif archetypes enriched in closing regions and every motif enriched in opening regions (BH FDR < 0.05 across both directions), one JASPAR motif per archetype (Vierstra v2.1beta). Sequence logos are drawn from the JASPAR count matrices used in the analysis, with letter height showing information content in bits (pseudocount 0.25 per base). Bars show the percentage of target regions (colour) and matched background regions (grey) containing the motif, followed by the fold enrichment. The final column gives the HOMER P value and the BH-adjusted FDR. Motif enrichment does not establish binding by a specific family member.

## supp - Fig 3.20. Transcription-factor motifs enriched in regions with altered accessibility in Y-depleted HT1080 cells.

<a id="fig-3-20"></a>

**Figure:** `ch03_main_ht1080_atac_motif_logos` | **Chapter number:** Fig 3.20

**Working ID:** `CH3-S3.4.6.1-F08` | **Figure:** `ch03_main_ht1080_atac_motif_logos` | **Status:** selected main figure

[Selected PDF](figures/chapter3/selected/page_reduction_20261030/ch03_supp_ht1080_atac_motif_logos.pdf) | [Selected PNG](figures/chapter3/selected/page_reduction_20261030/ch03_supp_ht1080_atac_motif_logos.png)

**Supp Fig 3.41. Transcription-factor motifs enriched in regions with altered accessibility in Y-depleted HT1080 cells.**

**a,b,** Selected motifs enriched in regions with lower (a) or higher (b) accessibility (ATAC FDR < 0.05). HOMER enrichment used GC- and baseline-accessibility-matched controls, with Benjamini–Hochberg correction across both directions. Letter heights indicate sequence information content. Coloured and grey bars show motif frequencies in target and background regions; annotations give fold enrichment, *P* values and FDR.

## exploratory - Motif enrichment restricted to promoter regions, with sequence logos

**Figure:** `ch03_exploratory_ht1080_atac_motif_logos_promoter` | **Status:** exploratory (22 September 2026)

[PDF](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_promoter.pdf) | [SVG](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_promoter.svg) | [PNG](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_promoter.png)

This is the HOMER analysis of `ch03_supp_ht1080_atac_motif_logos_all_significant`, restricted to promoter regions (overlapping a knownGene TSS ± 2 kb): 129 closing and 17 opening regions at FDR < 0.05. Each set was compared with GC- and accessibility-matched non-differential promoter regions (FDR ≥ 0.20, 10 per target), so enrichment is promoter against promoter. No motif reached BH FDR < 0.05 (best: KLF11 and NR4A1, FDR 0.05). The rows show the ten lowest-P archetypes for closing and the four lowest-P archetypes for opening regions, and should not be read as significant enrichment. Motif families enriched in the all-region analysis (TFAP2, ZIC, RUNX, TEAD, FOX, BCL11B, GATA) retain fold enrichments of about 1.4–1.7 in closing promoters, without significance. No opening promoter contains a Sox6 match. The motif signal in the full analysis is therefore carried mainly by non-promoter regions, which make up about 80% of the differential regions, although the smaller promoter sets also limit power.

## exploratory - Motif enrichment in promoter regions at ATAC FDR < 0.10, with sequence logos

**Figure:** `ch03_exploratory_ht1080_atac_motif_logos_promoter_fdr010` | **Status:** exploratory (22 September 2026)

[PDF](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_promoter_fdr010.pdf) | [SVG](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_promoter_fdr010.svg) | [PNG](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_promoter_fdr010.png)

This is the promoter-restricted HOMER analysis of `ch03_exploratory_ht1080_atac_motif_logos_promoter`, repeated at ATAC FDR < 0.10 to match the RNA–ATAC integration figure. It includes 248 closing and 67 opening promoter regions (overlapping a knownGene TSS ± 2 kb). Each set was compared with GC- and accessibility-matched non-differential promoter regions (FDR ≥ 0.20). Limited promoter controls allowed five per closing target and ten per opening target. BH correction was applied across 879 motifs × 2 directions. Only NFYB (the NF-Y CCAAT box) reached FDR < 0.05, in closing promoters (13% vs 6%, 2.3×, FDR 0.008). The remaining rows are the lowest-P archetypes and are not significant. They include FOXC1, KLF13 and RARB for closing promoters, and ZBTB18, SOX15, DBP and MAF for opening promoters. Motif families that were significant across all regions (TFAP2, EBF, ZIC, RUNX, TEAD, BCL11B) show fold enrichments of only about 1.2–1.4 in closing promoters, and Sox6 is depleted in opening promoters.

## exploratory - Motif enrichment in distal regions at ATAC FDR < 0.10, with sequence logos

**Figure:** `ch03_exploratory_ht1080_atac_motif_logos_distal_fdr010` | **Status:** exploratory (22 September 2026)

[PDF](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_distal_fdr010.pdf) | [SVG](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_distal_fdr010.svg) | [PNG](figures/chapter3/ch03_exploratory_ht1080_atac_motif_logos_distal_fdr010.png)

HOMER known-motif enrichment in distal regions (gene body or intergenic, not overlapping a knownGene TSS ± 2 kb) at ATAC FDR < 0.10: 1,407 closing and 331 opening regions. Each set was compared with GC- and accessibility-matched non-differential distal regions (FDR ≥ 0.20). The number of available controls in the sparsest strata limited matching to one control per closing target and six per opening target. BH correction was applied across 879 motifs × 2 directions, and 96 closing and 20 opening motifs reached FDR < 0.05. The rows are the ten lowest-FDR archetypes for closing regions and six for opening regions; opening motifs are dominated by related SOX models. Closing distal regions reproduce the all-region motif families at higher significance: TFAP2 (2.1×, FDR 2e-43), ZIC, BCL11B, FOX, GATA, RUNX, EBF and REL. Opening distal regions show a consistent SOX signal (Sox5 2.9×, Sox6 2.1×, SOX4, SOX8, SOX9) and a nuclear-receptor DR3 motif (RXRA::VDR 2.6×). Together with the promoter-only analyses, which found only NFYB at FDR < 0.05, these results place the motif signal mainly in distal regions.

## supp - Supp Fig 3.22. Motif enrichment at promoter and distal regions.

<a id="supp-fig-3-22"></a>

**Figure:** `ch03_supp_ht1080_atac_promoter_vs_distal` | **Status:** draft for review (22 September 2026) | **Chapter number:** Supp Fig 3.22

[PDF](figures/chapter3/ch03_supp_ht1080_atac_promoter_vs_distal.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_atac_promoter_vs_distal.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_atac_promoter_vs_distal.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S11_HT1080_ATAC_promoter_vs_distal/CH3_S3_4_6_1_S11_HT1080_ATAC_promoter_vs_distal.pdf)

**Supp Fig 3.22. Motif enrichment at promoter and distal regions.**

Differential regions (ATAC FDR < 0.10) were separated into promoters (TSS ±2 kb) and distal regions. HOMER enrichment used GC- and accessibility-matched non-differential controls of the same class, with Benjamini–Hochberg correction across motifs and both directions within each class. **a,** Selected motifs; colour shows log₂ fold enrichment and dot size indicates FDR. Filled points mark FDR < 0.05. **b,** Number of significantly enriched motifs in each class and direction.

## supp - Fig 3.19. Genomic context of chromatin-accessibility changes in Y-depleted HT1080 cells.

<a id="fig-3-19"></a>

**Figure:** `ch03_main_ht1080_atac_copy_number` | **Chapter number:** Fig 3.19

**Working ID:** `CH3-S3.4.6.1-F07` | **Figure:** `ch03_main_ht1080_atac_copy_number` | **Status:** selected main figure (24 September 2026)

[PDF](figures/chapter3/ch03_main_ht1080_atac_copy_number.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_atac_copy_number.svg) | [PNG](figures/chapter3/ch03_main_ht1080_atac_copy_number.png)

[Selected PDF](figures/chapter3/selected/page_reduction_20261030/ch03_supp_ht1080_atac_genomic_context.pdf)

**Supp Fig 3.40. Genomic context of chromatin-accessibility changes in Y-depleted HT1080 cells.**

**a,** Genome-centred chromosome-median differences in ATAC (passage points), RNA (diamonds) and lpWGS (bars); the lower strip gives percentages of regions losing/gaining accessibility at FDR < 0.05. **b,** Pairwise assay comparisons across the same 39 autosomal arms, without significance filtering. ATAC averages three passage-specific arm medians; RNA uses genome-centred arm-median effects for combined, August and May comparisons; lpWGS uses arm-median segment differences. Lines show linear fits; correlations have two-sided, unadjusted *P* values. Arms are genomic summaries, not biological replicates; sampling times were not established as matched across assays. **c,** Accessibility along chromosomes 18 and 13, with 5-Mb medians (black), lpWGS differences (purple) and significant regions (ticks). **d,** Motif enrichment before/after excluding chromosome 13/18 losses or chromosome 18 gains: Mantel–Haenszel odds ratios and 95% confidence intervals, adjusted for GC content and baseline accessibility. **e,** Genes with promoter-accessibility loss (ATAC FDR < 0.10), before/after restricting to small lpWGS differences (|Δlog₂ copy ratio| < 0.2); labels give gene counts. lpWGS came from earlier D0 samples and does not establish copy number in the ATAC-profiled populations.

## supp - Fig 3.31. Expression responses by baseline promoter-accessibility class.

<a id="fig-3-31"></a>

**Figure:** `ch03_main_ht1080_addback_by_atac_class` | **Working ID:** CH3-S3.3.5-F03 | **Results:** 3.3.5 | **Status:** Main figure simplified to two panels on 8 October 2026; former c is retained in historical exports. | **Chapter number:** Fig 3.31

[Selected PDF](figures/chapter3/selected/page_reduction_20261030/ch03_supp_ht1080_addback_by_atac_class.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_3_5_F03_ht1080_addback_by_atac_class/ch03_main_ht1080_addback_by_atac_class.svg) | [Selected PNG](figures/chapter3/selected/page_reduction_20261030/ch03_supp_ht1080_addback_by_atac_class.png)

[PDF](figures/chapter3/ch03_main_ht1080_addback_by_atac_class.pdf) | [SVG](figures/chapter3/ch03_main_ht1080_addback_by_atac_class.svg) | [PNG](figures/chapter3/ch03_main_ht1080_addback_by_atac_class.png)

**Supp Fig 3.43. Expression responses by baseline promoter-accessibility class.**

**a,** Expression responses in HT1080 genes down-regulated following Y loss, comparing 58 genes with reduced promoter accessibility (ATAC FDR < 0.10; blue) with 116 other down-regulated genes matched for baseline RNA effect and mean expression (two controls per promoter-loss gene; grey). Points represent genes; lines show medians and interquartile ranges. Expression differences are the first named population minus the second, so positive values indicate higher expression in the first. Annotations give the adjusted additional response in promoter-loss genes (Δ) and its FDR. **b,** Adjusted differences in expression response between promoter-loss and other down-regulated genes across eight population comparisons. Positive estimates indicate a larger response in promoter-loss genes; shaded rows compare each population with empty vector. The Δ annotations in a and estimates in b use all 225 eligible genes (58 promoter-loss and 167 other down-regulated genes), adjusting for baseline RNA effect and mean expression, with HC3 standard errors. Bars show 95% confidence intervals; filled points indicate BH FDR < 0.05 across eight comparisons. Each population had three matched collection occasions; the combined population descended from the sampled UTY population. Accessibility was measured only before add-back.

## supp - Supp Fig 3.19. Accessibility trends and pathway associations in Y-depleted HT1080 cells.

<a id="supp-fig-3-19"></a>

**Figure:** `ch03_supp_ht1080_atac_gene_set_context` | **Status:** draft for review (22 September 2026) | **Chapter number:** Supp Fig 3.19

[PDF](figures/chapter3/ch03_supp_ht1080_atac_gene_set_context.pdf) | [SVG](figures/chapter3/ch03_supp_ht1080_atac_gene_set_context.svg) | [PNG](figures/chapter3/ch03_supp_ht1080_atac_gene_set_context.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S15_HT1080_ATAC_gene_set_context/CH3_S3_4_6_1_S15_HT1080_ATAC_gene_set_context.pdf)

**Supp Fig 3.19. Accessibility trends and pathway associations in Y-depleted HT1080 cells.**

**a,** Median accessibility changes across ten baseline-accessibility quantiles, before and after loess correction of the relationship between baseline signal and estimated fold change. **b,** Hallmark promoter-accessibility GSEA under the primary and corrected models, showing sets significant in either analysis. **c,** Hallmark normalised enrichment scores (NES) for corrected promoter accessibility and RNA, using the same 13,300 genes. Filled teal points in b and c indicate accessibility FDR < 0.05. **d,** GREAT enrichment for distal regions with lower or higher accessibility (ATAC FDR < 0.10), showing the top 20 closing-associated and all opening-associated sets passing both region- and gene-based FDR < 0.05. Points indicate region-based fold enrichment; columns give both FDRs. **e,** Opening-associated enrichment before and after chromosome 18 exclusion. Open squares indicate sets no longer passing both tests; labels identify sets below the five-gene reporting threshold.

## main - Fig 3.24. Transcriptional variation and expression differences in DU145 cells.

<a id="fig-3-24"></a>

**Figure:** `ch03_main_du145_rna_global_structure` | **Working ID:** CH3-S3.4.6.3-F01 | **Section:** 3.4.6.3 | **Status:** Selected | **Chapter number:** Fig 3.24

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_F01_DU145_RNA_global/CH3_S3_4_6_3_F01_DU145_RNA_global.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_F01_DU145_RNA_global/CH3_S3_4_6_3_F01_DU145_RNA_global.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_F01_DU145_RNA_global/CH3_S3_4_6_3_F01_DU145_RNA_global_genes_6a63603091d2.png)

**Fig 3.24. Transcriptional variation and expression differences in DU145 cells.**

Three passages (P12–P14) of one sgCTRL and one Y-depleted population. **a,** PCA of 18,767 autosomal/non-PAR X genes; shapes identify passages and dashed lines connect matching passages. **b,** Differential expression between Y-depleted (LOY) and sgCTRL cells, using passage-blocked voom–limma analysis (residual df = 2) with BH correction across these genes. Red marks FDR < 0.05, |log₂ fold change| ≥ 0.5 and agreement across all passages; blue/orange mark other significant genes with/without copy-number concordance; grey marks nonsignificant genes. **c,** Minimum versus maximum passage-specific expression differences for each gene, used to assess whether changes retain the same direction across all three passages. Shading marks genes consistently higher or lower in Y-depleted cells. Red denotes the highlighted subset, orange other significant genes. **d,** Mean RNA versus median lpWGS arm differences (≥ 20 genes/arm). Point area reflects gene count; insets give Pearson correlations at gene, window and arm resolution. **e,** Copy-number context of all 147 differentially expressed genes (FDR < 0.05) and the 114-gene subset with |log₂ fold change| ≥ 0.5 and consistent direction across passages. Stacked bars show whether expression changes coincide with DNA dosage changes: stable (grey), |lpWGS difference| < 0.2; concordant (blue) or discordant (orange), |lpWGS difference| ≥ 0.2 with DNA and RNA changes in the same or opposite direction. **f,** GRHL2 CPM and group means. **g,h,** TAF9 (g) and PXDN (h) TMM-normalised log₂ CPM (prior count 2). Points show the three passage samples per condition; boxes show medians and interquartile ranges, with whiskers extending to observations within 1.5 interquartile ranges. FDRs and fold changes in f–h are from the model in b.

[Source and reproduction notes](figures/chapter3/source_data/du145_global_expression_combined/README.md)

## main - Fig 3.25. Y-linked expression tracks regional Y loss in DU145 cells.

<a id="fig-3-25"></a>

**Figure:** `ch03_main_du145_rna_sex_chromosome` | **Working ID:** CH3-S3.4.6.3-F02 | **Section:** 3.4.6.3 | **Status:** Selected | **Chapter number:** Fig 3.25

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_F02_DU145_RNA_sex_chromosomes/CH3_S3_4_6_3_F02_DU145_RNA_sex_chromosomes.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_F02_DU145_RNA_sex_chromosomes/CH3_S3_4_6_3_F02_DU145_RNA_sex_chromosomes.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_F02_DU145_RNA_sex_chromosomes/CH3_S3_4_6_3_F02_DU145_RNA_sex_chromosomes.png)

**Fig 3.25. Y-linked expression tracks regional Y loss in DU145 cells.**

**a,** Y-linked RNA relative to sgCTRL, with lpWGS ratios (black); shading marks PAR1, distal Yp (to 4.8 Mb) and remaining MSY. Values ≤ 1% are displayed at 1%; dashed DNA segments fall below this limit and control copy ratios < 0.05 are omitted. **b,** Summed regional expression. **c,d,** Residual Y (c, log axis) and PAR expression (d, linear), requiring mean sgCTRL CPM ≥ 1/≥ 5. Combined PAR sums all PAR1/PAR2 features. **e,f,** X–Y paralog log₂ CPM differences (e) and summed X + Y expression relative to sgCTRL (f), requiring Y-linked sgCTRL CPM ≥ 5. Blue/red denote sgCTRL/Y-depleted in b,f and X/Y in e.

[Source and reproduction notes](figures/chapter3/source_data/du145_section_revision/README.md)

## main - DU145 autosomal expression differences and their copy-number context (retired standalone figure)

**Figure:** `ch03_main_du145_rna_expression_changes` | **Working ID:** CH3-S3.4.6.3-F03 | **Section:** 3.4.6.3 | **Status:** Retired — merged into F01, panels b–f (9 October 2026)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_F03_DU145_RNA_expression/CH3_S3_4_6_3_F03_DU145_RNA_expression.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_F03_DU145_RNA_expression/CH3_S3_4_6_3_F03_DU145_RNA_expression.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_F03_DU145_RNA_expression/CH3_S3_4_6_3_F03_DU145_RNA_expression.png)

**DU145 expression differences and their copy-number context.** Differences are Y-depleted minus sgCTRL; lpWGS uses the LOY-2X population.

**a,** Paired voom–limma analysis (~ passage + condition; residual df = 2), with Benjamini–Hochberg correction across 18,767 autosomal/non-PAR X genes. Red marks FDR <0.05, |log₂ fold change| ≥0.5 and matching directions across all three passage pairs. Dashed lines mark the effect threshold.

**b,** Minimum and maximum log₂ CPM differences across passage pairs P12–P14. Shading indicates consistent direction. Red denotes genes meeting the criteria in a; orange, other FDR <0.05 genes; grey, nonsignificant genes.

**c,** Mean RNA log₂ CPM differences versus median lpWGS log₂ differences for autosomal/X arms with ≥20 mapped genes. Point area reflects gene count; inset values give descriptive Pearson correlations at gene, 10-Mb-window and arm resolution.

**d,** Copy-number classes among significant genes and the highlighted subset. Stable indicates an absolute lpWGS segment difference <0.2; concordant/discordant changes meet the ≥0.2 threshold and match/oppose RNA direction.

**e,** GRHL2 CPM in three passage samples per condition; horizontal marks show means.

Genes were mapped to CHM13v2.0; gene-level DNA summaries use median saved QDNAseq segment values across overlapping usable 100-kb bins.

[Source and reproduction notes](figures/chapter3/source_data/du145_section_revision/README.md)


## main - Fig 3.26. Transcriptional programmes associated with Y loss in DU145 cells.

<a id="fig-3-26"></a>

**Figure:** `ch03_main_du145_rna_programmes_regulators` | **Working ID:** CH3-S3.4.6.3-F04 | **Section:** 3.4.6.3 | **Status:** Selected | **Chapter number:** Fig 3.26

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_F04_DU145_RNA_programmes/CH3_S3_4_6_3_F04_DU145_RNA_programmes.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_F04_DU145_RNA_programmes/CH3_S3_4_6_3_F04_DU145_RNA_programmes.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_F04_DU145_RNA_programmes/CH3_S3_4_6_3_F04_DU145_RNA_programmes.png)

**Fig 3.26. Transcriptional programmes associated with Y loss in DU145 cells.**

Three passages of one sgCTRL and one Y-depleted population. **a,** Effect-ranked fgsea NES for 19 selected gene sets; blue/red denote negative/positive enrichment and asterisks FDR < 0.05. Adjacent genes are the three leading-edge members with largest absolute mean expression differences, without a gene-level significance filter. **b,** fgsea and passage-blocked CAMERA (fixed correlation 0.01) FDRs, corrected separately across all 23,252 sets; bold marks FDR < 0.05. H, Hallmark; R, Reactome; CGP, perturbation signatures; MF/CC, GO Molecular Function/Cellular Component. CGP up/down suffixes describe reference signatures. **c,** All 154 significant fgsea sets, grouped by leading-edge composition; colour indicates retained significance and direction after excluding genes with |lpWGS segment difference| ≥ 0.2. **d,** Passage-specific Y-depleted-minus-sgCTRL log₂ CPM differences for selected genes (circles/squares/triangles, P12/P13/P14). Filled points indicate primary gene-level FDR < 0.05; grey lines span passage differences, not confidence intervals. The right column gives mean sgCTRL CPM. Gene sets overlap and measure coordinated expression, not pathway or enzyme activity. Complete statistics and sensitivity results are in the [supporting tables](tables/chapter3/du145_section_review/README.md); individual significant sets appear in **Supp Fig** 3.36.

[Source and reproduction notes](figures/chapter3/source_data/du145_rna_programmes_regulators/README.md)

## supp - Supp Fig 3.31. Expression patterns of passage-consistent genes in DU145 cells.

<a id="supp-fig-3-31"></a>

**Figure:** `ch03_supp_du145_rna_gene_detail` | **Working ID:** CH3-S3.4.6.3-S01 | **Section:** 3.4.6.3 | **Status:** Selected | **Chapter number:** Supp Fig 3.31

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_S01_DU145_RNA_gene_detail/CH3_S3_4_6_3_S01_DU145_RNA_gene_detail.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_S01_DU145_RNA_gene_detail/CH3_S3_4_6_3_S01_DU145_RNA_gene_detail.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_S01_DU145_RNA_gene_detail/CH3_S3_4_6_3_S01_DU145_RNA_gene_detail.png)

**Supp Fig 3.31. Expression patterns of passage-consistent genes in DU145 cells.**

**a,** TMM-normalised log₂ CPM for 114 genes with FDR < 0.05, |log₂ fold change| ≥ 0.5 and the same direction of change across all three passage pairs. Values are standardised within each gene across six samples (row z-scores; colour limits ±2). Genes are ordered by model log₂ fold change across three consecutive blocks. Columns are grouped by condition, with sgCTRL (C) followed by Y-depleted (Y) samples and passages 12–14 ordered within each group; the white divider separates conditions. Expression is not adjusted for passage. **b,** Model log₂ fold change against mean CPM for all tested genes. Red marks the 114-gene subset; orange marks other genes with FDR < 0.05. Dashed lines mark log₂ fold changes of ±0.5; the dotted line marks 5 CPM.

For external context, the DepMap 26Q1 DU145 reference profile (ACH-000979) reports GRHL2 at log₂(TPM + 1) = 0.02. This is a separate stock and uses different expression units from the experimental CPM values; the measurements are not directly comparable.

[Source and reproduction notes](figures/chapter3/source_data/du145_section_revision/README.md)

## supp - Supp Fig 3.30. DU145 copy-number context, model sensitivity and sample structure.

<a id="supp-fig-3-30"></a>

**Figure:** `ch03_supp_du145_rna_copy_number` | **Working ID:** CH3-S3.4.6.3-S02 | **Section:** 3.4.6.3 | **Status:** Selected | **Chapter number:** Supp Fig 3.30

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_S02_DU145_RNA_copy_number/CH3_S3_4_6_3_S02_DU145_RNA_copy_number.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_S02_DU145_RNA_copy_number/CH3_S3_4_6_3_S02_DU145_RNA_copy_number.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_S02_DU145_RNA_copy_number/CH3_S3_4_6_3_S02_DU145_RNA_copy_number_combined_d7d755c6f4e5.png)

**Supp Fig 3.30. DU145 copy-number context, model sensitivity and sample structure.**

Y-depleted-minus-sgCTRL differences. **a,b,** RNA versus lpWGS differences at 10-Mb-window (a) and gene (b) resolution, excluding Y. Labels in a summarise extreme windows by chromosome; b states off-scale counts. **c,** Arms containing the most FDR < 0.05 genes. **d,** Two Y-depleted lpWGS populations versus the same control, distinguishing shared from population-specific changes. **e,** Passage-paired versus unpaired model effects for the same libraries; the unpaired model has no FDR < 0.05 genes. **f,** PC1 by passage, with within-condition lines and categorical-passage R². **g,** MSY reads as a percentage of assigned reads; annotation gives the ratio of condition means. **h,** Pearson correlations of TMM log₂ CPM, with average-linkage clustering (1 − r). Blue/red denote sgCTRL/Y-depleted in f–h. C/Y denote populations; numbers denote P12–P14. RNA samples are three passages of one established population per condition. Copy-number-restricted refit results are in the [supplementary DE table](tables/chapter3/ch03_supp_du145_copy_number_stable_DE.csv).

[Source and reproduction notes](figures/chapter3/source_data/du145_global_expression_combined/README.md)

## supp - Supp Fig 3.36. DU145 gene-set enrichment and contributing genes.

<a id="supp-fig-3-36"></a>

**Figure:** `ch03_supp_du145_rna_gene_set_themes` | **Working ID:** CH3-S3.4.6.3-S03 | **Section:** 3.4.6.3 | **Status:** Selected | **Chapter number:** Supp Fig 3.36

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_S03_DU145_RNA_gene_sets/CH3_S3_4_6_3_S03_DU145_RNA_gene_sets.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_S03_DU145_RNA_gene_sets/CH3_S3_4_6_3_S03_DU145_RNA_gene_sets.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_S03_DU145_RNA_gene_sets/CH3_S3_4_6_3_S03_DU145_RNA_gene_sets_two_panels_94f7d4c1ebd5.png)

**Supp Fig 3.36. DU145 gene-set enrichment and contributing genes.**

**a,** All 154 gene sets passing BH correction across 23,252 MSigDB sets, positioned by NES within leading-edge themes. Colour marks retained significance/direction after excluding genes with |lpWGS segment difference| ≥ 0.2; counts give retained sets. **b,** Six recurrent leading-edge genes per theme, selected by membership count with effect-size/symbol tie-breaking and no gene-level significance filter. Tiles show unscaled Y-depleted-minus-sgCTRL log₂ CPM at P12–P14 (blue/red, lower/higher), clipped to ±1.5; printed values retain the full range. Sets gives gene memberships/retained sets in that theme. Asterisks denote passage-blocked gene-level FDR < 0.05 across 18,767 genes. Memberships reflect overlapping annotations. Other and chromosomal-position themes appear only in a. Complete statistics, theme assignments and leading edges accompany the source data.

[Source and reproduction notes](figures/chapter3/source_data/du145_rna_gene_set_themes/README.md)

## main - Fig 3.27. Shared and background-dependent expression patterns after Y depletion.

<a id="fig-3-27"></a>

**Figure:** `ch03_main_crossline_rna_comparison` | **Results:** 3.3.4.4 | **Chapter number:** Fig 3.27

[PDF](figures/chapter3/ch03_main_crossline_rna_comparison.pdf) | [Editable SVG](figures/chapter3/ch03_main_crossline_rna_comparison.svg) | [PNG](figures/chapter3/ch03_main_crossline_rna_comparison.png)

**Fig 3.27. Shared and background-dependent expression patterns after Y depletion.**

**a,** Grouped bars show male-specific Y-gene expression relative to sgCTRL on a logarithmic axis. Numbers above bars give the underlying percentages. Shading marks distal Yp; annotations give summed PAR1 expression. Bars begin at 0.07%; values ≤ 0.1% are displayed at 0.1%. Genes below 1 CPM in sgCTRL are omitted. **b,** Grouped horizontal bars show X-linked paralog expression differences (Y-depleted minus sgCTRL mean log₂ CPM), extending from zero. Dotted lines mark median differences across non-PAR X genes. Colours in a and b identify cell lines. **c,** Pairwise Spearman correlations between cell lines of LOY-associated differential expression (Y-depleted minus sgCTRL) across 16,165 commonly tested genes. **d,** Normalised enrichment scores (NES) for selected immune, rRNA-processing and cell-cycle sets. Blue and red indicate lower and higher ranked expression; asterisks mark FDR < 0.05 across all eligible sets within each background. The comparison uses mean-expression effects, with equal weighting of HT1080 sampling periods and the earlier condition-only A549 analysis. Gene-set scores describe expression patterns rather than direct pathway activity.

[Source and reproduction notes](figures/chapter3/ch03_main_crossline_rna_comparison_notes.md)

## supp - Supp Fig 3.38. Expression signatures across external Y-loss and related datasets.

<a id="supp-fig-3-38"></a>

**Figure:** `ch03_supp_crossline_external_comparison` | **Results:** 3.3.4.4 | **Chapter number:** Supp Fig 3.38

[PDF](figures/chapter3/ch03_supp_crossline_external_comparison.pdf) | [Editable SVG](figures/chapter3/ch03_supp_crossline_external_comparison.svg) | [PNG](figures/chapter3/ch03_supp_crossline_external_comparison.png)

**Supp Fig 3.38. Expression signatures across external Y-loss and related datasets.**

Eight selected signatures across the three thesis models and 18 external contrasts. Colour shows AUC −0.5: blue indicates lower ranked expression and red higher ranked expression in the Y-loss-oriented comparison. Dots mark nominal empirical *P* < 0.05 from 500 size-matched random gene sets, without correction across comparisons. IFN-α, IFN-γ, inflammatory, G2/M and E2F are Hallmark sets. The other columns show the CDH1-target signature and genes increased after BMI1 or MEL18 knockdown in DAOY medulloblastoma cells (BMI1_DN.V1_UP and MEL18_DN.V1_UP). Blue in these two columns therefore indicates the opposite direction to the reference knockdown response. Groups distinguish somatic LOY, natural mosaicism, constitutional 45,X and related Y-dosage comparisons. The asterisk identifies multiome profiles from the same A549 clones used for bulk RNA-seq. Sources and processing are in the [external dataset table](tables/chapter3/Table_CH3_external_expression_datasets.md).

[Source and reproduction notes](figures/chapter3/ch03_supp_crossline_external_comparison_notes.md)

## Retired supp - Shared and divergent gene sets across the three lines

**Status:** Removed from the active Results during the author-requested cross-line simplification on 8 October 2026. Caption, source data and exports are retained as historical records. The complete results remain in the [shared and divergent gene-set table](tables/chapter3/crossline_comparison/shared_and_divergent_gene_sets.tsv).

**Figure:** `ch03_supp_crossline_gene_sets` | **Results:** 3.4.6.4

[PDF](figures/chapter3/ch03_supp_crossline_gene_sets.pdf) | [Editable SVG](figures/chapter3/ch03_supp_crossline_gene_sets.svg) | [PNG](figures/chapter3/ch03_supp_crossline_gene_sets.png)

**Shared and divergent gene-set enrichment across cell backgrounds.**

Heatmap of normalised enrichment scores for 46 shared and 23 divergent gene sets, separated by the horizontal line. Dots mark adjusted P <0.05 across all eligible sets within each background. The two adjacent HT1080 columns indicate enrichment direction in each sampling period. Right-hand annotations identify leading-edge themes, chromosomal concentration and sensitivity to exclusion of copy-number-changed genes in DU145. Full gene-set names, statistics and leading edges are provided in the source tables.

[Source and reproduction notes](figures/chapter3/source_data/crossline_gene_sets/provenance.json)

## Retired supp - PXDN and dataset polarity in the cross-line replication panel

**Status:** Removed from Results_edit.md at the author’s request on 8 October 2026. Caption and exports are retained as historical records.

**Figure:** `ch03_supp_crossline_external_detail` | **Results:** 3.4.6.4

[PDF](figures/chapter3/ch03_supp_crossline_external_detail.pdf) | [Editable SVG](figures/chapter3/ch03_supp_crossline_external_detail.svg) | [PNG](figures/chapter3/ch03_supp_crossline_external_detail.png)

**PXDN expression across Y-loss and related datasets.**

**a,** PXDN expression differences in the thesis models and reanalysed external datasets, oriented so that negative values indicate lower expression in the Y-loss comparison. Thesis-model values are mean log2 CPM differences; external values are log2 fold changes. Labels show estimates and FDR where available. Bars are clipped at ±3, with full values displayed; “not tested” indicates absence from the imported results. Colours identify dataset groups. **b,** Comparison-direction checks, aligned with panel a. Green indicates a male-specific Y-gene estimate consistent with the assigned direction; grey indicates direction taken from the documented contrast.

[Source and reproduction notes](figures/chapter3/source_data/crossline_external_detail/provenance.json)


## Retired supp - Y-linked DNA dosage following individual-guide targeting in H2452

**Status:** Removed from Results_edit.md with the H2452 Y-targeting account at the author’s request on 8 October 2026. Caption and exports are retained as historical records.

**Figure:** `ch03_supp_h2452_ddpcr_verified_initial_guides` | **Working ID:** CH3-S3.4.5.4-S01 | **Results:** 3.4.5.4 | **Title:** Y-linked DNA dosage following individual-guide targeting in H2452.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_4_S01_h2452_initial_guides/CH3_S3_4_5_4_S01_h2452_initial_guides.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_5_4_S01_h2452_initial_guides/CH3_S3_4_5_4_S01_h2452_initial_guides.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_5_4_S01_h2452_initial_guides/CH3_S3_4_5_4_S01_h2452_initial_guides.png)

[PDF](figures/chapter3/ch03_supp_h2452_ddpcr_verified_initial_guides.pdf) | [Editable SVG](figures/chapter3/ch03_supp_h2452_ddpcr_verified_initial_guides.svg) | [PNG](figures/chapter3/ch03_supp_h2452_ddpcr_verified_initial_guides.png)

AMELY/AMELX ratios measured by ddPCR in sgCTRL and CP-, CN- or Q53-targeted populations. Points show individual DNA measurements from two batches, distinguished by circles and triangles. Horizontal lines and labels show medians. Two measurements were available for sgCTRL, CP and CN, and one for Q53.

[Source and reproduction notes](figures/chapter3/ch03_supp_h2452_ddpcr_verified_initial_guides_notes.md)

<a id="supp---gene-level-expression-and-promoter-accessibility-changes-in-ht1080"></a>

## supp - Supp Fig 3.23. Chromosome 19 context of concordant expression and promoter-accessibility changes.

<a id="supp-fig-3-23"></a>

**Figure:** `ch03_supp_ht1080_joint_promoter_genes` | **Selected:** CH3-S3.4.6.1-S18 | **Chapter number:** Supp Fig 3.23

[PDF](figures/chapter3/selected/CH3_S3_4_6_1_S18_HT1080_joint_promoter_genes/CH3_S3_4_6_1_S18_HT1080_joint_promoter_genes.pdf)

**Supp Fig 3.23. Chromosome 19 context of concordant expression and promoter-accessibility changes.**

**a,** RNA and promoter ATAC log₂ fold changes for 63 jointly significant genes. Filled points show RNA effects; open points show median promoter ATAC effects. Blue/red indicates concordant decreases/increases; black indicates opposite changes. Bold names identify chromosome 19 genes. RNA differential expression required FDR < 0.05, |log₂ fold change| ≥ 1 and matching directions in both sampling periods; ATAC significance required FDR < 0.10.

**b,** Chromosome 19 ATAC changes across tested peaks (grey), with significant losses at FDR < 0.05 highlighted in blue. The black line connects 1-Mb medians; purple points show earlier D0 lpWGS differences.

**c,** Enlargement of 50–59 Mb, with ZNF-associated promoter peaks labelled. Shared peak assignments are listed together. Shading marks 52–53 Mb, where 22 of 26 peaks had lower ATAC signal (median log₂ fold change, −0.53). Coordinates are hg38.

## supp - Supp Fig 3.44. Global transcriptional variation in A549 cells.

<a id="fig-3-22"></a>
<a id="supp-fig-3-44"></a>

**Figure:** `ch03_main_a549_rna_global_structure` | **Working ID:** CH3-S3.4.6.2-F01 | **Section:** 3.4.6.2 | **Status:** Selected supplementary | **Chapter number:** Supp Fig 3.44 (current draftv2; former Fig 3.11)

[Selected PDF](figures/chapter3/selected/a549_rna_to_supp_20261030/ch03_supp_a549_rna_global_structure.pdf) | [Selected SVG](figures/chapter3/selected/a549_rna_to_supp_20261030/ch03_supp_a549_rna_global_structure.svg) | [Selected PNG](figures/chapter3/selected/a549_rna_to_supp_20261030/ch03_supp_a549_rna_global_structure_preview_aec24a3459c7.png)

[PDF](figures/chapter3/selected/a549_rna_to_supp_20261030/ch03_supp_a549_rna_global_structure.pdf) | [Editable SVG](figures/chapter3/selected/a549_rna_to_supp_20261030/ch03_supp_a549_rna_global_structure.svg) | [PNG](figures/chapter3/selected/a549_rna_to_supp_20261030/ch03_supp_a549_rna_global_structure_preview_aec24a3459c7.png)

**Supp Fig 3.44. Global transcriptional variation in A549 cells.**

Four sgCTRL and four Y-depleted samples. **a,** Unadjusted PCA of 22,792 autosomal/non-PAR X genes; blue/red denote sgCTRL/Y-depleted and numbers identify samples. **b,** RINe-adjusted, unpaired voom–limma analysis with BH correction across 22,843 nuclear genes. Brown/teal mark non-PAR Y/PAR; blue/red mark significant decreases/increases among other genes; grey marks other nonsignificant genes. Positive effects indicate higher Y-depleted expression. The horizontal line marks the largest nominal *P* among FDR < 0.05 genes; ±1 lines are visual guides. **c,** Selected gene-set enrichment: blue/red tiles show negative/positive NES; columns give fgsea and fixed-correlation CAMERA (0.01) FDRs, corrected across 23,457 sets. Bold values indicate FDR < 0.05. These pathway analyses use unadjusted expression differences, without RINe adjustment. H, Hallmark; R, Reactome; GO, Gene Ontology Biological Process.

**Supplementary asset:** `ch03_supp_a549_rna_global_structure`. Legacy Figure/Working IDs are retained as provenance keys.

[Source and reproduction notes](figures/chapter3/selected/a549_rna_to_supp_20261030/source_notes.md)

## main - Fig 3.23. Regional Y-linked and pseudoautosomal expression in A549 cells.

<a id="fig-3-23"></a>

**Figure:** `ch03_main_a549_rna_sex_chromosomes` | **Working ID:** CH3-S3.4.6.2-F02 | **Section:** 3.4.6.2 | **Status:** Selected | **Chapter number:** Fig 3.23

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_2_F02_A549_RNA_sex_chromosomes/CH3_S3_4_6_2_F02_A549_RNA_sex_chromosomes.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_2_F02_A549_RNA_sex_chromosomes/CH3_S3_4_6_2_F02_A549_RNA_sex_chromosomes.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_2_F02_A549_RNA_sex_chromosomes/CH3_S3_4_6_2_F02_A549_RNA_sex_chromosomes.png)

[PDF](figures/chapter3/ch03_main_a549_rna_sex_chromosomes.pdf) | [Editable SVG](figures/chapter3/ch03_main_a549_rna_sex_chromosomes.svg) | [PNG](figures/chapter3/ch03_main_a549_rna_sex_chromosomes.png)

**Fig 3.23. Regional Y-linked and pseudoautosomal expression in A549 cells.**

**a,** Y-linked RNA and lpWGS ratios relative to sgCTRL. Black lines show DNA; dotted PAR1 uses corrected-bin medians. Shading marks PAR1, distal Yp (to 4.8 Mb) and remaining MSY. Values ≤ 1% are displayed at 1%, with dashed DNA segments below this limit; control copy ratios < 0.05 are omitted. **b,** Summed regional expression. **c,d,** Residual Y (c, log axis) and PAR expression (d, linear), requiring mean sgCTRL CPM ≥ 1/≥ 5. In c, positive values below 1% are retained and open triangles mark zero. Combined PAR sums all PAR features. **e,f,** X–Y paralog log₂ CPM differences (e) and summed X + Y expression relative to sgCTRL (f), requiring Y-linked sgCTRL CPM ≥ 5. Blue/red denote sgCTRL/Y-depleted in b,f and X/Y in e. RNA uses descriptive, unadjusted TMM-normalised means from four samples per condition; DNA was sampled separately. Pooled expression is abundance-weighted and does not establish functional compensation.

[Source and reproduction notes](figures/chapter3/source_data/shared_sex_chromosome/A549/provenance.json)

## main - Transcriptional programme differences in Y-depleted A549 cells

**Figure:** `ch03_main_a549_rna_pathways` | **Working ID:** CH3-S3.4.6.2-F03 | **Section:** 3.4.6.2 | **Status:** Historical — merged into F01 panel c, 9 October 2026

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_2_F03_A549_RNA_pathways/CH3_S3_4_6_2_F03_A549_RNA_pathways.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_2_F03_A549_RNA_pathways/CH3_S3_4_6_2_F03_A549_RNA_pathways.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_2_F03_A549_RNA_pathways/CH3_S3_4_6_2_F03_A549_RNA_pathways.png)

[PDF](figures/chapter3/ch03_main_a549_rna_pathways.pdf) | [Editable SVG](figures/chapter3/ch03_main_a549_rna_pathways.svg) | [PNG](figures/chapter3/ch03_main_a549_rna_pathways.png)

**Transcriptional programme differences in Y-depleted A549 cells.** **a,** The eight significant Hallmark sets (H) and selected Reactome (R) and GO Biological Process (GO) sets, grouped by theme. Rankings use mean log2CPM differences. Heatmap colours show normalised enrichment scores, with blue negative and red positive; asterisks mark fgsea FDR <0.05. **b,** Benjamini–Hochberg-adjusted P values across 23,457 MSigDB sets for fgsea and CAMERA using the default fixed inter-gene correlation of 0.01; bold values meet FDR <0.05. The additional sets illustrate themes rather than an independent test family. Several nominal fgsea P values reached the analysis floor of 10⁻¹⁰.

[Source and reproduction notes](figures/chapter3/a549_rna_thesis_notes.md)


## supp - Supp Fig 3.28. Parental A549 cells show higher interferon and inflammatory expression.

<a id="supp-fig-3-28"></a>

**Figure:** `ch03_supp_a549_parental_comparisons` | **Working ID:** CH3-S3.4.6.2-S03 | **Section:** 3.4.6.2 | **Status:** Selected | **Chapter number:** Supp Fig 3.28

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_2_S03_A549_RNA_parental/CH3_S3_4_6_2_S03_A549_RNA_parental.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_2_S03_A549_RNA_parental/CH3_S3_4_6_2_S03_A549_RNA_parental.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_2_S03_A549_RNA_parental/CH3_S3_4_6_2_S03_A549_RNA_parental.png)

[PDF](figures/chapter3/ch03_supp_a549_parental_comparisons.pdf) | [Editable SVG](figures/chapter3/ch03_supp_a549_parental_comparisons.svg) | [PNG](figures/chapter3/ch03_supp_a549_parental_comparisons.png)

**Supp Fig 3.28. Parental A549 cells show higher interferon and inflammatory expression.**

**a,b,** Unpaired voom–limma differential expression comparisons of three parental samples with four sgCTRL (a) or four Y-depleted samples (b). Positive log₂ fold changes indicate higher parental expression. Genes meeting FDR < 0.05 are dark grey when higher in parental cells, blue when higher in sgCTRL and red when higher in Y-depleted cells; other genes are light grey. Dashed lines mark the comparison-specific nominal P thresholds corresponding to Benjamini–Hochberg FDR < 0.05. **c,** Normalised enrichment scores for Hallmark sets meeting global FDR < 0.05 in either comparison, using genes ranked by mean log₂ CPM differences. Positive scores indicate enrichment towards parental cells. Asterisks mark Benjamini–Hochberg FDR < 0.05 across all tested MSigDB sets within each comparison.

[Source and reproduction notes](figures/chapter3/a549_parental_comparisons_notes.md)

## supp - Supp Fig 3.27. RNA integrity, model sensitivity and Y-linked expression in A549.

<a id="supp-fig-3-27"></a>

**Figure:** `ch03_supp_a549_rna_qc` | **Working ID:** CH3-S3.4.6.2-S01 | **Section:** 3.4.6.2 | **Status:** Selected | **Chapter number:** Supp Fig 3.27

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_2_S01_A549_RNA_QC/CH3_S3_4_6_2_S01_A549_RNA_QC.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_2_S01_A549_RNA_QC/CH3_S3_4_6_2_S01_A549_RNA_QC.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_2_S01_A549_RNA_QC/CH3_S3_4_6_2_S01_A549_RNA_QC.png)

[PDF](figures/chapter3/ch03_supp_a549_rna_qc.pdf) | [Editable SVG](figures/chapter3/ch03_supp_a549_rna_qc.svg) | [PNG](figures/chapter3/ch03_supp_a549_rna_qc.png)

**Supp Fig 3.27. RNA integrity, model sensitivity and Y-linked expression in A549.**

Four sgCTRL (blue) and four Y-depleted samples (red). **a–c,** PC1 versus RINe, exonic fraction and mitochondrial fraction; PCA uses 22,792 autosomal/non-PAR X genes. Numbers identify unpaired samples; ρ denotes Spearman correlation. **d,e,** Condition-only versus RINe-adjusted voom–limma models, sharing TMM normalisation and BH correction across 22,843 nuclear genes. Brown/teal mark non-PAR Y/PAR; blue/red mark significant decreases/increases among other genes; grey marks other nonsignificant genes. Positive effects indicate higher Y-depleted expression. Horizontal lines mark the largest nominal *P* among significant genes; ±1 lines are visual guides. **f,** Unadjusted, unscaled TMM log₂ CPM for 13 Y genes, ordered by nominal *P*. Statistics come from the primary RINe-adjusted model with correction across all 22,843 genes; KDM5D is bold.

[Source and reproduction notes](figures/chapter3/a549_qc_expanded_notes.md)

## supp - Supp Fig 3.29. Hallmark expression patterns across A549 samples.

<a id="supp-fig-3-29"></a>

**Figure:** `ch03_supp_a549_rna_hallmark_ssgsea` | **Working ID:** CH3-S3.4.6.2-S06 | **Section:** 3.4.6.2 | **Status:** Selected | **Chapter number:** Supp Fig 3.29

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_2_S06_A549_RNA_Hallmark_ssGSEA/CH3_S3_4_6_2_S06_A549_RNA_Hallmark_ssGSEA.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_2_S06_A549_RNA_Hallmark_ssGSEA/CH3_S3_4_6_2_S06_A549_RNA_Hallmark_ssGSEA.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_2_S06_A549_RNA_Hallmark_ssGSEA/CH3_S3_4_6_2_S06_A549_RNA_Hallmark_ssGSEA.png)

[PDF](figures/chapter3/ch03_supp_a549_rna_hallmark_ssgsea.pdf) | [Editable SVG](figures/chapter3/ch03_supp_a549_rna_hallmark_ssgsea.svg) | [PNG](figures/chapter3/ch03_supp_a549_rna_hallmark_ssgsea.png)

**Supp Fig 3.29. Hallmark expression patterns across A549 samples.**

ssGSEA scores for 50 Hallmark gene sets in four sgCTRL and four Y-depleted samples. Scores were standardised within each gene set across all eight samples (row z-scores). Red indicates scores above the row mean and blue below; colours do not indicate statistical significance. Numbers identify unpaired samples.

[Source and reproduction notes](figures/chapter3/a549_hallmark_ssgsea_notes.md)

## supp - Supp Fig 3.2. RNA-seq quality control across 44 in-house libraries.

**Figure:** `ch03_methods_rna_library_qc` | **Chapter number:** Supp Fig 3.2 | **Placement:** Methods | **Status:** Draft for author review; compact methods version of `ch03_supp_rna_all_samples_qc_draft`

[PDF](figures/chapter3/ch03_methods_rna_library_qc.pdf) | [Editable SVG](figures/chapter3/ch03_methods_rna_library_qc.svg) | [PNG](figures/chapter3/ch03_methods_rna_library_qc.png)

**Supp Fig 3.2. RNA-seq quality control across 44 in-house libraries.**

Libraries comprise A549 (11), August HT1080 (9), May HT1080 including add-backs (18), and DU145 (6), ordered identically across panels. **a,** Trimmed read pairs entering STAR. **b,** STAR alignment classes as percentages of input pairs. **c,** Qualimap exonic, intronic and intergenic fractions of classified reads. **d,** fastp duplication estimates. **e,** Mitochondrial alignments as a percentage of mapped alignment records (blue), and rRNA-associated bases as a percentage of Picard aligned bases (orange); denominators differ. **f,** Qualimap 5′/3′ coverage bias; the dashed line denotes equal coverage. Condition colours apply in a,d,f. P, parental; C, sgCTRL; Yd, Y-depleted; EV, empty vector; U, UTY; K, KDM5D; UK, both add-backs. Indices identify samples or collection occasions, not independent population derivations. All 44 libraries are shown.

[Source and reproduction notes](figures/chapter3/ch03_methods_rna_library_qc_notes.md)

## supp - Supp Fig 3.3. Read quality, input RNA and gene-body coverage across 44 RNA-seq libraries.

**Figure:** `ch03_methods_rna_input_quality_coverage` | **Chapter number:** Supp Fig 3.3 | **Placement:** Methods | **Status:** Draft for author review; compact methods version of `ch03_supp_rna_quality_coverage_draft`, with panels reordered

[PDF](figures/chapter3/ch03_methods_rna_input_quality_coverage.pdf) | [Editable SVG](figures/chapter3/ch03_methods_rna_input_quality_coverage.svg) | [PNG](figures/chapter3/ch03_methods_rna_input_quality_coverage.png)

**Supp Fig 3.3. Read quality, input RNA and gene-body coverage across 44 RNA-seq libraries.**

All 44 libraries, sample order and condition abbreviations follow **Supp Fig** 3.2. **a,** Bases with Phred scores ≥ 30 before (grey) and after (blue) fastp processing. **b,c,** Adapter-trimmed (b) and retained reads (c), each as a percentage of raw reads; trimmed reads can be retained. **d,** Input RNA integrity: RINe for A549/August HT1080, RQN for May HT1080, and provider-reported RIN for DU145. **e,** DV200 (percentage of RNA >200 nt), except for May HT1080, which shows the fraction within 200–10,000 nt. **f,** Qualimap coverage across 100 relative gene-body positions, divided by each library's mean. Curves represent libraries; dashed lines mark mean-normalised coverage of one.

[Source and reproduction notes](figures/chapter3/ch03_methods_rna_input_quality_coverage_notes.md)

## Retired - RNA-seq QC supplementary drafts

`ch03_supp_rna_all_samples_qc_draft` and `ch03_supp_rna_quality_coverage_draft` were replaced by the two methods figures above at the author's request on 27 September 2026. Their exports and notes remain in `figures/chapter3/`; the historical captions are preserved in `archive/methods_rna_qc_captions_20260927_143008`. In the methods version of the second figure, the draft's panels d, e and f are now f, d and e.

## supp - Sensitivity of A549 expression differences to RNA integrity

**Figure:** `ch03_supp_a549_rna_RINe_sensitivity` | **Working ID:** CH3-S3.4.6.2-S02 | **Section:** 3.4.6.2 | **Status:** Historical — merged into S01 panels d and e, 9 October 2026

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_2_S02_A549_RNA_RINe/CH3_S3_4_6_2_S02_A549_RNA_RINe.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_2_S02_A549_RNA_RINe/CH3_S3_4_6_2_S02_A549_RNA_RINe.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_2_S02_A549_RNA_RINe/CH3_S3_4_6_2_S02_A549_RNA_RINe.png)

[PDF](figures/chapter3/ch03_supp_a549_rna_RINe_sensitivity.pdf) | [Editable SVG](figures/chapter3/ch03_supp_a549_rna_RINe_sensitivity.svg) | [PNG](figures/chapter3/ch03_supp_a549_rna_RINe_sensitivity.png)

**Sensitivity of A549 expression differences to RNA integrity.** **a, b,** Voom–limma results for the same 22,843 expressed nuclear genes and four samples per condition using a condition-only sensitivity model (a) or the primary model including continuous RINe (b). Both models used the same TMM effective library sizes, robust empirical Bayes moderation and Benjamini–Hochberg correction across all 22,843 genes. The condition-only model identified 16 genes at FDR <0.05, all non-PAR Y-linked; the RINe-adjusted model identified 27, comprising 17 non-PAR Y-linked genes and ten non-Y/non-PAR genes. Brown, non-PAR Y; teal, PAR. Among other genes, blue and red identify significant decreases and increases, respectively, and grey indicates FDR ≥0.05. Each horizontal dashed line marks the largest nominal P value among significant genes in that model; vertical dotted lines at log₂ fold changes of ±1 are visual guides. Axes are shared. Complete results are available in the [gene-level results tables](tables/chapter3/a549_primary_nuclear_RINe/README.md).

[Source and reproduction notes](figures/chapter3/a549_RINe_sensitivity_notes.md)

## supp - Y-linked expression changes dominate the A549 volcano

**Figure:** `ch03_supp_a549_rna_Y_inclusive_volcano` | **Working ID:** CH3-S3.4.6.2-S04 | **Section:** 3.4.6.2 | **Status:** Historical — heatmap merged into S01f; duplicate volcano omitted, 9 October 2026

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_2_S04_A549_RNA_Y_inclusive/CH3_S3_4_6_2_S04_A549_RNA_Y_inclusive.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_2_S04_A549_RNA_Y_inclusive/CH3_S3_4_6_2_S04_A549_RNA_Y_inclusive.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_2_S04_A549_RNA_Y_inclusive/CH3_S3_4_6_2_S04_A549_RNA_Y_inclusive.png)

[PDF](figures/chapter3/ch03_supp_a549_rna_Y_inclusive_volcano.pdf) | [Editable SVG](figures/chapter3/ch03_supp_a549_rna_Y_inclusive_volcano.svg) | [PNG](figures/chapter3/ch03_supp_a549_rna_Y_inclusive_volcano.png)

**Y-linked expression depletion in A549 cells.** **a,** RINe-adjusted voom–limma differential expression results for 22,843 expressed nuclear genes in four sgCTRL and four Y-depleted samples, using the primary model described in Figure CH3-S3.4.6.2-F01. Brown, non-PAR Y-linked genes; teal, PAR genes; grey, other genes. Dashed lines mark nominal P = 0.05 and log₂ fold changes of ±1 as visual guides; significance was assessed by Benjamini–Hochberg FDR <0.05 across all 22,843 genes. **b,** Per-sample TMM-normalised log2CPM for selected Y-linked genes, alongside nominal P values and full-family FDR from the RINe-adjusted model. Rows are ordered by nominal P; displayed expression values are neither RINe-adjusted nor row-scaled. KDM5D is bold.

[Source and reproduction notes](figures/chapter3/a549_Y_inclusive_volcano_notes.md)

## main - Regional Y-linked and pseudoautosomal expression in HT1080 cells

**Figure:** `ch03_main_ht1080_rna_sex_chromosomes` | **Former Results placement:** 3.4.6.1 | **Status:** Omitted from current section by author decision, 30 September 2026; regional characterisation is covered in 3.4.5.1. Historical draft caption and assets retained.

[PDF](figures/chapter3/ch03_main_ht1080_rna_sex_chromosomes.pdf) | [PNG](figures/chapter3/ch03_main_ht1080_rna_sex_chromosomes.png)

**Regional Y-linked and pseudoautosomal expression in HT1080 cells.** **a,** Mean Y-linked expression relative to sgCTRL along chromosome Y, with saved lpWGS segment ratios in black. Shading marks PAR1, distal Yp and the remaining male-specific region (MSY); the common distal-Yp boundary is 4.8 Mb. **b,** Mean summed regional CPM; percentages indicate residual expression relative to sgCTRL. **c, d,** Residual expression of individual Y-linked (c) and pseudoautosomal (d) genes. Gene order is shared across cell lines; points require mean sgCTRL CPM ≥1 for Y genes and ≥5 for PAR genes. **e,** Mean log2CPM differences (Y-depleted minus sgCTRL) for X–Y paralog pairs with Y-linked sgCTRL CPM ≥5. Blue and red identify sgCTRL and Y-depleted cells in b, and X- and Y-linked genes in e. Values ≤1% are displayed at 1% in a, c and d; dashed lpWGS segments indicate values below this display limit. Intervals with control relative copy number <0.05 are omitted from the DNA comparison. Circles and triangles identify the two sampling periods, each with three passage samples per condition. Solid and hatched bars identify the same periods; paired percentages follow this order. RNA and lpWGS samples were collected separately.

[Source and reproduction notes](figures/chapter3/source_data/shared_sex_chromosome/HT1080/provenance.json)


## supp - Clustered expression of 114 genes in DU145

**Figure:** `ch03_supp_du145_rna_clustered_heatmap` | **Results:** 3.4.6.3 | **Status:** Alternative for review

[PDF](figures/chapter3/ch03_supp_du145_rna_clustered_heatmap.pdf) | [PNG](figures/chapter3/ch03_supp_du145_rna_clustered_heatmap.png)

**Clustered expression of passage-consistent genes in DU145 cells.** Row-standardised TMM log2CPM for 114 genes meeting FDR <0.05, |log2 fold change| ≥0.5 and consistent direction across all three passage pairs. Genes are clustered by average linkage using correlation distance; samples are clustered within each condition using Euclidean distance on row z-scores. Conditions remain grouped. Colours show row z-scores across all six libraries, limited to ±2. Labels indicate passage. Expression is not adjusted for passage.

[Source and reproduction notes](figures/chapter3/source_data/du145_rna_clustered_heatmap/provenance.json)


## supp - Supp Fig 3.37. Histone-associated expression increases without a detected local DNA copy-number gain in DU145 cells.

<a id="supp-fig-3-37"></a>

**Figure:** `ch03_supp_du145_chr6_dna_rna` | **Working ID:** CH3-S3.4.6.3-S05 | **Section:** 3.4.6.3 | **Status:** Selected | **Chapter number:** Supp Fig 3.37

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_S05_DU145_chr6_DNA_RNA/CH3_S3_4_6_3_S05_DU145_chr6_DNA_RNA.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_S05_DU145_chr6_DNA_RNA/CH3_S3_4_6_3_S05_DU145_chr6_DNA_RNA.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_S05_DU145_chr6_DNA_RNA/CH3_S3_4_6_3_S05_DU145_chr6_DNA_RNA.png)

**Supp Fig 3.37. Histone-associated expression increases without a detected local DNA copy-number gain in DU145 cells.**

**a,b,** lpWGS log₂ copy ratios across chromosome 6 and the histone-cluster interval, respectively. Points show usable 100-kb bins; lines show saved segments for sgCTRL (blue) and Y-depleted cells (red). **c,d,** Corresponding Y-depleted minus sgCTRL DNA differences; dashed lines mark the absolute 0.2 screening threshold. **e,f,** RNA log₂ fold changes from the primary passage-blocked model. Teal identifies histone genes, with black outlines indicating FDR < 0.05; grey points show other tested genes. Shading marks the span of mapped chromosome 6p histone genes. Coordinates use T2T-CHM13v2.0. RNA represents three passage samples per population; lpWGS used one separately collected sample per condition. DNA points outside the displayed axes are omitted visually. Matched axes are used for whole-chromosome and enlarged views.

[Source and reproduction notes](figures/chapter3/source_data/du145_section_revision/README.md)

## Retired main - Inferred regulator patterns vary across cell backgrounds

**Status:** Removed from the active Results during the author-requested cross-line simplification on 8 October 2026. Caption, source data and exports are retained as historical records. The former stability/module Methods are preserved in the [pre-edit archive](archive/crossline_condense_20261008_175215/Chapter3_methods.md).

File: `figures/chapter3/ch03_main_crossline_tf_regulators.pdf`

**Inferred regulator patterns vary across cell backgrounds.**

**a,** CollecTRI ULM scores for regulators passing all six stability checks in at least one background, plus FOXA2 and GRHL2. Blue and red indicate lower and higher inferred activity; filled and open circles mark robust and supported regulators, respectively. **b,** Regulators reaching ULM FDR <0.05, grouped by stability tier. Labels give significant/tested counts. **c,** Module stability across gene bootstraps and random removal of regulatory edges, measured by mean Jaccard similarity to the original module. Shading marks the ≥0.75 threshold required in both tests; asterisks identify modules below it. Adjacent columns show unadjusted fry P values, module members and inferred directions. Stability criteria and line-specific sampling checks are defined in Methods.


## Retired supp - Sensitivity of inferred regulator modules

**Status:** Removed from Results_edit.md at the author’s request on 8 October 2026. Caption and exports are retained as historical records.

File: `figures/chapter3/ch03_supp_crossline_tf_module_sensitivity.pdf`

**Sensitivity of inferred regulator modules.**

**a,** Empirical P values for regulator counts, clustering separation, target overlap and residual coactivity. The dashed line marks P = 0.05. Regulator-count tests assess the first two stability criteria; permutation procedures are described in Methods. **b,** Signed module-target scores after excluding members of the indicated gene sets. Vertical lines mark scores before exclusion; open red circles indicate loss of nominal significance. Positive scores indicate agreement with inferred regulator directions. Annotations give expression-matched permutation and unadjusted fry P values for the complete target sets. Axes are scaled separately; P values are not corrected across modules or sensitivity tests.


## supp - Inferred regulators and shared targets in HT1080

File: `figures/chapter3/ch03_supp_ht1080_tf_regulator_modules.pdf`

**Inferred regulators and shared targets in HT1080.**

**a,** ULM scores, FDR and six stability checks for robust (bold) and supported regulators. Filled and open circles indicate passed and failed checks, respectively; criteria are defined in Methods. **b,** Weighted overlap among contributing target genes. The outlined module contains 27 regulators; NFE2L2 remains unclustered. **c,** Leading module targets, ranked by summed contribution, with their moderated t statistics. Blue indicates lower inferred activity or expression in Y-depleted cells. Module assignments correspond across panels.


## supp - Inferred regulators and shared targets in A549

File: `figures/chapter3/ch03_supp_a549_tf_regulator_modules.pdf`

**Inferred regulators and shared targets in A549.**

**a,** ULM scores, FDR and six stability checks for robust (bold) and supported regulators. Filled and open circles indicate passed and failed checks, respectively; criteria are defined in Methods. **b,** Weighted overlap among contributing target genes. The outlined module comprises HNF1A, HNF1B and supported regulator FOXA2; NCOR1 remains unclustered. **c,** Leading module targets, ranked by summed contribution, with their moderated t statistics. Blue and red indicate lower and higher inferred activity or expression in Y-depleted cells, respectively. Module assignments correspond across panels.


## supp - Inferred regulators and shared targets in DU145

File: `figures/chapter3/ch03_supp_du145_tf_regulator_modules.pdf`

**Inferred regulators and shared targets in DU145.**

**a,** ULM scores, FDR and six stability checks for robust (bold) and supported regulators. Filled and open circles indicate passed and failed checks, respectively; criteria are defined in Methods. **b,** Weighted overlap among contributing target genes, with module outlines. **c,** Leading targets within each module, ranked by summed contribution, with their moderated t statistics. Dashed outlines and asterisks identify modules below the stability threshold. Blue and red indicate lower and higher inferred activity or expression in Y-depleted cells, respectively. Module assignments correspond across panels.


## Retired supp - RNA and copy-number correspondence in A549 and DU145

File: `figures/chapter3/ch03_supp_a549_du145_rna_dna_copy_number.pdf`

**Status:** Removed from Results_edit.md at the author’s request on 8 October 2026. Caption and exports are retained as historical records. Current genome-wide, window and arm-level comparisons are shown in CH3-S3.3.6.2-S07; its Y-region definitions differ from the whole-Y-arm calculation retained below.

**RNA and copy-number correspondence in A549 and DU145.** **a, e,** Genome-wide lpWGS differences (Y-depleted minus sgCTRL) in A549 (a) and DU145 (e), summarised in fixed 5-Mb windows. Corrected unsegmented log2 relative copy-number values were summarised by the median separately in each population across paired usable 100-kb bins, then subtracted. Source bins overlapping Yq12 (chrY:27,449,937-62,460,029; zero-based BED coordinates) were excluded before aggregation, including the proximal boundary-crossing bin. The excluded interval is shaded; boundary windows terminate at the retained source span. Each track contains 630 windows, including six Y windows. Chromosome Y is red; other chromosomes alternate grey. The vertical axis uses symmetric-log spacing, linear within +/-0.3; no resulting window values are clipped.

**b, c, f, g,** DNA differences versus mean RNA log2 CPM differences at gene (b, f) and 10-Mb-window (c, g) resolution. These panels retain the expressed autosomal/non-PAR X gene sets from the current thesis analyses. Genes were mapped to CHM13v2.0; gene-level DNA differences use the difference between population-specific median saved segment values across overlapping bins. Window summaries compare mean gene-expression differences with median paired DNA-bin segment differences, requiring at least 10 mapped expressed genes. Gene-panel DNA axes use symmetric-log spacing, linear within +/-0.25. In g, the four windows furthest from the median DNA/RNA point after standard-deviation scaling, plus the windows with the largest positive DNA and RNA differences, are highlighted and labelled by their full overlapping CHM13 cytoband spans. A label spanning p and q denotes a window crossing the centromere; it is not a focal alteration assigned to a single band.

**d, h,** Arm-level comparisons, now including Yp and Yq in red. Non-Y summaries retain their original values and require at least 20 mapped expressed genes (41 arms). Y arms require at least one mapped expressed non-PAR Y gene: A549 Yp/Yq contain four/12 genes and DU145 Yp/Yq five/nine genes. Saved Y-gene expression effects use the same normalized log2 CPM scale as the non-Y effects. Genes were mapped by stable Ensembl identifier on the same chromosome, with unique same-chromosome symbol matching as fallback. Yq12-overlapping DNA bins and RNA gene spans were excluded from Y-arm summaries. DNA arm values are median paired saved segment log2 differences; RNA values are mean gene-expression differences. Both axes use symmetric-log spacing, linear within +/-0.3. Pearson correlations are shown for all 43 arms and for the 41 non-Y arms separately, because the large Y changes strongly influence the overall correlations. Yp and Yq are labelled in both panels, with selected non-Y arms labelled in DU145.

Corresponding panels share axis limits. Correlations use the original log2 differences, without applying the display-axis transformation. One lpWGS sample represents each population; RNA and DNA samples were not established as matched, and genomic units are not biological replicates. These comparisons are descriptive, not tests of causality. No RNA model or DNA segmentation was refitted. Local differences in repetitive regions are not independently validated structural alterations, and numerical-floor values do not establish precise fold depletion. Original LNCaP sample labels for DU145 are preserved in the source data.


## supp - Supp Fig 3.9. Genome-wide DNA differences and regional RNA correspondence in A549 and DU145.

<a id="supp-fig-3-9"></a>

**Figure:** `ch03_supp_a549_du145_dna_rna_arm_comparison` | **Chapter number:** Supp Fig 3.9

**Working ID:** CH3-S3.4.6.2-S07 | **Role:** supp | **Status:** Selected | **Results:** 3.4.5.2 and 3.4.5.3

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_2_S07_a549_du145_dna_rna_arms/CH3_S3_4_6_2_S07_a549_du145_dna_rna_arms.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_2_S07_a549_du145_dna_rna_arms/CH3_S3_4_6_2_S07_a549_du145_dna_rna_arms.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_2_S07_a549_du145_dna_rna_arms/CH3_S3_4_6_2_S07_a549_du145_dna_rna_arms.png)

[Working PDF](figures/chapter3/ch03_supp_a549_du145_dna_rna_arm_comparison.pdf) | [Source and reproduction notes](figures/chapter3/source_data/a549_du145_rna_dna/README.md)

**Supp Fig 3.9. Genome-wide DNA differences and regional RNA correspondence in A549 and DU145.**

Y-depleted-minus-sgCTRL differences in A549 (a–c,g) and DU145 (d–f,h). **a,d,** lpWGS log₂ relative copy-number differences across 630 5-Mb windows from usable 100-kb bins. Y is red; shading marks excluded Yq12. **b,e,** Median paired-bin DNA segment differences (grey) and mean TMM-normalised RNA log₂ CPM differences (green) across 41 non-Y arms with ≥ 20 genes. **c,f,** DNA–RNA comparisons for these arms plus male-specific Yp and euchromatic Yq (red), excluding PAR/centromere/Yq12 and requiring ≥ 1 gene. Pearson correlations include/exclude Y; labels mark four large non-Y deviations. **g,h,** Median DNA versus mean RNA differences across 10-Mb windows with ≥ 10 autosomal/non-PAR X genes (A549, 299; DU145, 296). Labels mark selected deviations. CHM13v2.0 coordinates. All DNA and RNA axes use linear spacing of the plotted log₂ differences. Panels c,f share DNA limits of −8 to 1 and RNA limits of −8.3 to 2.5; DNA tracks a,d use a shared symmetric range. Each population has one lpWGS sample, without established DNA–RNA matching. Genomic units are not biological replicates; DNA signals are relative and floor-level values limit interpretation.

## supp - Supp Fig 3.18. HT1080 ATAC-seq coverage and differential accessibility of chromosome Y and pseudoautosomal regions.

<a id="supp-fig-3-18"></a>

**Figure:** `ch03_supp_ht1080_atac_Y_PAR` | **Chapter number:** Supp Fig 3.18

Figure ID: CH3-S3.4.6.1-S19

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_1_S19_HT1080_ATAC_Y_PAR/CH3_S3_4_6_1_S19_HT1080_ATAC_Y_PAR.pdf)

**Supp Fig 3.18. HT1080 ATAC-seq coverage and differential accessibility of chromosome Y and pseudoautosomal regions.**

One sgCTRL (blue) and one Y-depleted population (red), each sampled at three passages. **a,** Coverage across the first 27 Mb of hg38 Y in 10-kb bins per million non-Y fragments. Shading marks PAR1; its gold overlay uses chrX alignments and is not Y-specific. **b,** Differential accessibility including MSY (purple circles) and chrX PAR peaks (gold triangles); other regions are grey. Y signal can reflect chromosome dosage as well as accessibility. **c–f,** Coverage at P2RY8, AKAP17A (chrX PAR1), RPS4Y1 and DDX3Y, respectively, in 100-bp bins. Lines show condition means; shading spans the three libraries. Black bars mark genes; insets show 3-kb windows around the strongest mean sgCTRL bin within each gene.

## supp - Supp Fig 3.11. Shared chromosome-arm copy-number differences across Y-targeted DU145 populations.

<a id="supp-fig-3-11"></a>

**Figure:** `ch03_supp_du145_arm_recurrence` | **Chapter number:** Supp Fig 3.11

[Selected PDF](figures/chapter3/selected/CH3_S3_4_5_5_S03_du145_arm_recurrence/CH3_S3_4_5_5_S03_du145_arm_recurrence.pdf)

**Supp Fig 3.11. Shared chromosome-arm copy-number differences across Y-targeted DU145 populations.**

Arm-median lpWGS log₂ segment differences in seven Y-targeted DU145 populations versus the same sgCTRL, using usable 100-kb bins and CHM13v2.0 cytobands. Guides, CP + Q53 + P47; 1×/2×, one/two successive nucleofections; Pano, panobinostat; KU, KU-55933. Blue/red indicate lower/higher relative dosage, saturated at ±0.8; 20p is outlined. Poorly resolved acrocentric short arms (13p, 14p, 15p, 21p and 22p) are masked in light grey with diagonal hatching and excluded from the colour display. Low-depth KU-only is shown but excluded from the six-profile recurrence summary. The bottom row is the guides-2× RNA comparison: mean differences across ≥ 20 genes/arm; plain grey marks unavailable values. This RNA comparison does not describe the other treatments, and DNA/RNA samples were not matched. Y-arm summaries include all usable bins. Signals are relative, and independent editing ancestry is unresolved.

## supp - Supp Fig 3.32. Pathway associations among genes with higher expression in Y-depleted DU145 cells.

<a id="supp-fig-3-32"></a>

**Figure:** `ch03_supp_du145_metascape_up_network` | **Chapter number:** Supp Fig 3.32

**Working ID:** CH3-S3.4.6.3-S06 | **Role:** supp | **Status:** Selected | **Results:** 3.4.6.3

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_S06_DU145_Metascape_UP/CH3_S3_4_6_3_S06_DU145_Metascape_UP.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_S06_DU145_Metascape_UP/CH3_S3_4_6_3_S06_DU145_Metascape_UP.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_S06_DU145_Metascape_UP/CH3_S3_4_6_3_S06_DU145_Metascape_UP.png)

[Working PDF](figures/chapter3/ch03_supp_du145_metascape_up_network.pdf)

**Supp Fig 3.32. Pathway associations among genes with higher expression in Y-depleted DU145 cells.**

Metascape term-similarity network for 31 higher-expression genes (Y-depleted versus sgCTRL; FDR < 0.05, |log₂ fold change| ≥ 0.5 and matching directions across three passages). Nodes show all 102 exported terms and 472 connections; area scales with −log₁₀ enrichment *P* and colour identifies six Louvain communities. The table gives term counts, strongest enrichment, representative terms and recurring genes. Positions have no biological-distance meaning. Statistics use the complete 31-gene list and a 15,651-gene enrichment denominator from the submitted tested-gene background. Reference-signature labels do not establish cell identity or perturbation in DU145.

[Source data, review and reproduction](figures/chapter3/source_data/du145_metascape_current/README.md)

## supp - Supp Fig 3.33. Pathway associations among genes with lower expression in Y-depleted DU145 cells.

<a id="supp-fig-3-33"></a>

**Figure:** `ch03_supp_du145_metascape_down_network` | **Chapter number:** Supp Fig 3.33

**Working ID:** CH3-S3.4.6.3-S07 | **Role:** supp | **Status:** Selected | **Results:** 3.4.6.3

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_S07_DU145_Metascape_DOWN/CH3_S3_4_6_3_S07_DU145_Metascape_DOWN.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_S07_DU145_Metascape_DOWN/CH3_S3_4_6_3_S07_DU145_Metascape_DOWN.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_S07_DU145_Metascape_DOWN/CH3_S3_4_6_3_S07_DU145_Metascape_DOWN.png)

[Working PDF](figures/chapter3/ch03_supp_du145_metascape_down_network.pdf)

**Supp Fig 3.33. Pathway associations among genes with lower expression in Y-depleted DU145 cells.**

Metascape term-similarity network for 83 lower-expression genes (Y-depleted versus sgCTRL; FDR < 0.05, |log₂ fold change| ≥ 0.5 and matching directions across three passages). Nodes show all 105 exported terms and 689 connections, including the isolated 20p13 term; area scales with −log₁₀ enrichment *P* and colour identifies six Louvain communities. The table gives term counts, strongest enrichment, representative terms and recurring genes. Positions have no biological-distance meaning. Statistics use the complete 83-gene list and a 15,651-gene enrichment denominator from the submitted tested-gene background. Reference-signature labels do not establish cell identity or perturbation in DU145.

[Source data, review and reproduction](figures/chapter3/source_data/du145_metascape_current/README.md)

## supp - Supp Fig 3.34. Enriched pathways and recurring contributing genes in DU145.

<a id="supp-fig-3-34"></a>

**Figure:** `ch03_supp_du145_metascape_terms_and_gene_modules` | **Chapter number:** Supp Fig 3.34

**Working ID:** CH3-S3.4.6.3-S08 | **Role:** supp | **Status:** Selected | **Results:** 3.4.6.3

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_S08_DU145_Metascape_terms_genes/CH3_S3_4_6_3_S08_DU145_Metascape_terms_genes.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_S08_DU145_Metascape_terms_genes/CH3_S3_4_6_3_S08_DU145_Metascape_terms_genes.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_S08_DU145_Metascape_terms_genes/CH3_S3_4_6_3_S08_DU145_Metascape_terms_genes.png)

[Working PDF](figures/chapter3/ch03_supp_du145_metascape_terms_and_gene_modules.pdf)

**Supp Fig 3.34. Enriched pathways and recurring contributing genes in DU145.**

**a,b,** Top 20 Metascape terms ranked by nominal *P* for higher- (a) and lower-expression genes (b). Bars show −log₁₀ *P*; labels give contributing-gene counts and approximate BH *q* values reconstructed from rounded exports. Rank selection does not require *q* < 0.05. **c,** Membership of all 31 higher-expression genes in the 102 UP network terms, coloured by community and grouped by the largest contribution. Genes are ordered by total recurrence; daggers mark tied community assignments. Recurrence reflects overlapping annotations. Inputs are the complete 31 higher-/83 lower-expression lists (FDR < 0.05, |log₂ fold change| ≥ 0.5, agreement across three passages). The enrichment denominator is 15,651 genes from the submitted tested-gene background.

[Source data and provenance](figures/chapter3/source_data/du145_metascape_terms_and_gene_modules/provenance.json)

## supp - DU145 Metascape figures, combined review copy

**Figure:** `ch03_supp_du145_metascape_review` | **Status:** Draft viewing copy

[Three-page PDF](figures/chapter3/ch03_supp_du145_metascape_review.pdf). Pages show the higher-expression network, lower-expression network, and top-term/recurring-gene figure, in that order. Captions and statistical definitions are given in the three entries above.


## supp - Supp Fig 3.35. Protein-interaction modules among differentially expressed genes in DU145.

<a id="supp-fig-3-35"></a>

**Figure:** `ch03_supp_du145_metascape_ppi_modules` | **Chapter number:** Supp Fig 3.35

**Working ID:** CH3-S3.4.6.3-S09 | **Role:** supp | **Status:** Selected | **Results:** 3.4.6.3

[Selected PDF](figures/chapter3/selected/CH3_S3_4_6_3_S09_DU145_Metascape_PPI/CH3_S3_4_6_3_S09_DU145_Metascape_PPI.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_6_3_S09_DU145_Metascape_PPI/CH3_S3_4_6_3_S09_DU145_Metascape_PPI.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_6_3_S09_DU145_Metascape_PPI/CH3_S3_4_6_3_S09_DU145_Metascape_PPI.png)

[Working PDF](figures/chapter3/ch03_supp_du145_metascape_ppi_modules.pdf)

**Supp Fig 3.35. Protein-interaction modules among differentially expressed genes in DU145.**

Metascape MCODE modules from 31 higher- and 83 lower-expression genes (Y-depleted versus sgCTRL; FDR < 0.05, |log₂ fold change| ≥ 0.5, agreement across three passages). **a,** Higher-expression collagen/extracellular-matrix module. **b,c,** Lower-expression immune/adhesion (b) and chemokine (c) modules. All 12 proteins and 17 exported database links are shown. Colours identify modules; node size, edge width and position do not encode expression or interaction strength. MCODE identifiers are assigned separately within UP/DOWN analyses. Links include physical, functional and inferred associations, not interactions measured in these samples. CCR9 has very low sgCTRL expression.

[Source data and reproduction notes](figures/chapter3/source_data/du145_metascape_ppi_modules/README.md)

## supp - DU145 original and copy-number-filtered enrichment scores

**Figure:** `ch03_supp_du145_gsea_original_vs_cn_scatter` | **Status:** Draft alternative for S03a; selected S03 unchanged

[PDF](figures/chapter3/ch03_supp_du145_gsea_original_vs_cn_scatter.pdf) | [SVG](figures/chapter3/ch03_supp_du145_gsea_original_vs_cn_scatter.svg) | [PNG](figures/chapter3/ch03_supp_du145_gsea_original_vs_cn_scatter.png)

**Comparison of enrichment scores before and after excluding genes in copy-number-altered intervals in DU145.** Each point represents one gene set significant in the original analysis (Benjamini–Hochberg FDR <0.05 across all tested collections). Original normalised enrichment scores (NES; horizontal axis) are compared with saved NES after removing genes in intervals with an absolute lpWGS segment difference ≥0.2 (vertical axis). Positive NES denotes enrichment towards higher expression in Y-depleted cells relative to sgCTRL; negative NES denotes enrichment towards lower expression. Colours identify the leading-edge themes used in S03. Filled symbols retain significance and matching enrichment direction after copy-number exclusion; open symbols do not. The dashed line indicates unchanged NES. Axes share the same linear scale and no jitter is applied.

Of the 154 originally significant sets, 150 have both NES values and are plotted. Four chromosomal-position sets (chr20p11, chr20p12, chr20p13 and chr11q24) have no NES in the saved filtered results and are listed separately, not assigned zero. Overall, 139 of the original 154 sets retain significance and direction. The chr11q23 set is labelled because it is the sole paired set with an absolute NES change ≥0.2; its filtered result does not meet the retained criterion. All other 149 paired NES changes are smaller than 0.2. The plotted sets were selected on original significance, and loss of significance after filtering does not by itself establish a copy-number-driven pathway effect. No enrichment analysis was rerun.

[Source data and reproduction notes](figures/chapter3/source_data/du145_gsea_original_vs_cn_scatter/README.md)

## main - Fig 3.3. DepMap reanalysis recovers four established LOY-associated dependencies.

<a id="fig-3-3"></a>

**Figure:** `ch03_main_depmap_y_status_dependencies` | **Working ID:** CH3-S3.3.1-F01 | **Earlier source ID:** CH3-S3.4.1-F01 | **Editing copy:** Results 3.3.1 | **Title:** DepMap reanalysis recovers four established LOY-associated dependencies. | **Chapter number:** Fig 3.3

[Selected PDF](figures/chapter3/selected/CH3_S3_4_1_F01_depmap_y_status_dependencies/CH3_S3_4_1_F01_depmap_y_status_dependencies.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_1_F01_depmap_y_status_dependencies/CH3_S3_4_1_F01_depmap_y_status_dependencies.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_1_F01_depmap_y_status_dependencies/CH3_S3_4_1_F01_depmap_y_status_dependencies.png)

[PDF](figures/chapter3/ch03_main_depmap_y_status_dependencies.pdf) | [Editable SVG](figures/chapter3/ch03_main_depmap_y_status_dependencies.svg) | [PNG](figures/chapter3/ch03_main_depmap_y_status_dependencies.png)

**Fig 3.3. DepMap reanalysis recovers four established LOY-associated dependencies.**

**a,** Mean log₂(TPM + 1) for seven Y genes in 711 female and 905 male profiles; dashed lines mark classification thresholds. **b,** Expression score versus normalised Y-gene copy ratio in 882 male models. Blue/red/grey denote Y-positive/Y-negative/unresolved calls. **c,** LOY prevalence in 17 lineages with ≥ 10 classified models (*n* = 875). **d,** Chronos associations across 17,205 genes, adjusted for lineage/library. Red, FDR < 0.05; green, other X–Y paralogs; grey, other genes. **e,** Four dependencies in 338 Y-positive male, 230 Y-negative male and 478 female models (blue/red/brown). Brackets give adjusted male-group effects (Δ, Y-negative minus Y-positive) and FDR; females are descriptive references. Boxes show medians and interquartile ranges (IQR); whiskers reach observations within 1.5 IQR. **f,** Source-scaled associations across CRISPR screening datasets and guide libraries: cell lines, Sanger organoids (47 Y-negative/48 Y-positive; minLib/v1_1 libraries) and NextGen models (27/38), with approximate model-based 95% intervals; filled symbols meet each source's FDR criterion. **g,** Library-adjusted effects across 17 lineages (≥ 10 models; ≥ 3 per Y group), with donor-clustered 95% intervals. Negative effects indicate stronger LOY dependency.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_1_F01_depmap_y_status_dependencies/selection_notes.md)

## supp — Evidence supporting integrated Y-status calls in DepMap

**Figure:** `ch03_supp_depmap_y_call_evidence` | **Working ID:** CH3-S3.4.1-S01 | **Results:** 3.4.1 | **Title:** Evidence supporting integrated Y-status calls in DepMap. | **Status:** Superseded in Results_edit.md by the combined CH3-S3.4.1-S01 supplement on 8 October 2026; this source caption is retained for historical exports.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_1_S01_depmap_y_call_evidence/CH3_S3_4_1_S01_depmap_y_call_evidence.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_1_S01_depmap_y_call_evidence/CH3_S3_4_1_S01_depmap_y_call_evidence.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_1_S01_depmap_y_call_evidence/CH3_S3_4_1_S01_depmap_y_call_evidence.png)

[PDF](figures/chapter3/ch03_supp_depmap_y_call_evidence.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_y_call_evidence.svg) | [PNG](figures/chapter3/ch03_supp_depmap_y_call_evidence.png)

**Evidence supporting integrated Y-status calls in DepMap.** Counts of all 1,143 annotated male models by the evidence supporting their integrated Y-status call. Blue, red and grey denote Y-positive (588), Y-negative (383) and set-aside (172) models, respectively; counts are printed for segments containing at least 40 models. Expression determines the integrated call, DNA supplies a call where RNA is unavailable, and definite disagreement sets a model aside. Intermediate RNA profiles, ambiguous DNA without RNA, and models lacking both measurements are also set aside. DNA copy ratios derive from whole-genome sequencing where available and otherwise from exome sequencing. This inventory precedes the cancer-model and engineering exclusions used for the lineage frequencies in [Figure CH3-S3.4.1-F01, panel c](figures/chapter3/selected/CH3_S3_4_1_F01_depmap_y_status_dependencies/CH3_S3_4_1_F01_depmap_y_status_dependencies.pdf).

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_1_S01_depmap_y_call_evidence/selection_notes.md)

## main — X-paralog dependency and expression of the Y-linked counterpart

**Figure:** `ch03_main_depmap_y_partner_expression` | **Working ID:** CH3-S3.4.2-F02 | **Results:** 3.4.2 | **Title:** X-paralog dependency and expression of the Y-linked counterpart.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_2_F02_depmap_y_partner_expression/CH3_S3_4_2_F02_depmap_y_partner_expression.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_2_F02_depmap_y_partner_expression/CH3_S3_4_2_F02_depmap_y_partner_expression.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_2_F02_depmap_y_partner_expression/CH3_S3_4_2_F02_depmap_y_partner_expression.png)

[PDF](figures/chapter3/ch03_main_depmap_y_partner_expression.pdf) | [Editable SVG](figures/chapter3/ch03_main_depmap_y_partner_expression.svg) | [PNG](figures/chapter3/ch03_main_depmap_y_partner_expression.png)

**X-paralog dependency and expression of the Y-linked counterpart.** **a,** Chronos gene effect of each X-linked gene against expression of its Y-linked paralog, in the 605 eligible male lines with both CRISPR and RNA data. This set includes the 69 lines set aside from the two-group call, most of which carry partial Y deletions and which provide the intermediate expression values. Colours denote the integrated Y-status call. **b,** Significance of each Y-linked gene as a term in a joint model containing all seven Y-linked genes together with lineage and library, fitted on the same lines. Numbers and shading show the negative log10 p of each term; outlined cells identify the paired gene of each dependency. Y-linked genes that are usually lost together carry overlapping signal, so the joint model localises each dependency to a gene or short interval rather than establishing a single causal gene. These are observational associations across unrelated lines.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_2_F02_depmap_y_partner_expression/selection_notes.md)

## main - Fig 3.4. Coverage, detection limits and candidate dependencies in cancer subsets.

<a id="fig-3-4"></a>

**Figure:** `ch03_main_depmap_cancer_type_detection_limits` | **Working ID:** CH3-S3.4.2-F03 | **Results:** 3.4.2 | **Title:** Coverage, detection limits and candidate dependencies in cancer subsets. | **Chapter number:** Fig 3.4

[Selected PDF](figures/chapter3/selected/CH3_S3_4_2_F03_depmap_cancer_type_detection_limits/CH3_S3_4_2_F03_depmap_cancer_type_detection_limits.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_2_F03_depmap_cancer_type_detection_limits/CH3_S3_4_2_F03_depmap_cancer_type_detection_limits.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_2_F03_depmap_cancer_type_detection_limits/CH3_S3_4_2_F03_depmap_cancer_type_detection_limits.png)

[PDF](figures/chapter3/ch03_main_depmap_cancer_type_detection_limits.pdf) | [Editable SVG](figures/chapter3/ch03_main_depmap_cancer_type_detection_limits.svg) | [PNG](figures/chapter3/ch03_main_depmap_cancer_type_detection_limits.png)

**Fig 3.4. Coverage, detection limits and candidate dependencies in cancer subsets.**

**a,b,** Cancer-group coverage among 568 Y-classified male CRISPR models. Eligibility requires ≥ 5 models per Y group; b uses square-root axes. **c,** Detection probability across 30 cancer groups after uniformly lowering LOY Chronos scores. Each cell uses 500 real gene profiles tested at within-group genome-wide FDR < 0.05. White dots mark the smallest shift reaching 80% detection; the final column is a 10-SD shift. Simulations retain between-line variation; Oncotree levels reuse models. **d,** INSM1 dependency in MYC-altered (7/13) and unaltered (21/31) lung models. **e,** ATRX dependency in CDKN2A-altered (20/5) and unaltered (9/0) melanoma models; counts are Y-positive/LOY. Points show raw Chronos scores (blue/red); lower scores indicate stronger dependency. Boxes show medians and interquartile ranges (IQR); whiskers reach observations within 1.5 IQR.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_2_F03_depmap_cancer_type_detection_limits/selection_notes.md)

## main — Replication of the Y-loss dependencies in organoids and other models

**Figure:** `ch03_main_organoid_rnai_replication` | **Working ID:** CH3-S3.4.3-F01 | **Results:** 3.4.3 | **Title:** Replication of the Y-loss dependencies in organoids and other models.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_3_F01_organoid_rnai_replication/CH3_S3_4_3_F01_organoid_rnai_replication.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_3_F01_organoid_rnai_replication/CH3_S3_4_3_F01_organoid_rnai_replication.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_3_F01_organoid_rnai_replication/CH3_S3_4_3_F01_organoid_rnai_replication.png)

[PDF](figures/chapter3/ch03_main_organoid_rnai_replication.pdf) | [Editable SVG](figures/chapter3/ch03_main_organoid_rnai_replication.svg) | [PNG](figures/chapter3/ch03_main_organoid_rnai_replication.png)

**Replication of the Y-loss dependencies in organoids and other models.** **a,** Numbers of Y-negative and Y-positive cell lines, organoids and combined models used in the colorectal and oesophagogastric comparisons. These panels include 92 organoids; the pan-cancer organoid analysis includes three additional models of other tumour types. **b,** Y-expression score against chromosome-Y copy number in 111 screened male organoids, coloured by the integrated call. Thirteen models were set aside because of intermediate or discordant evidence; further eligibility and patient-level selection yielded 95 organoids for the primary analysis. **c,** Association between Y status and gene effect across 16,547 genes in the 95 primary organoids, adjusting for tumour type, library and MSI. Negative effects indicate greater dependency in Y-negative organoids; effects are in Chronos-equivalent units, obtained by dividing the organoid score by 0.73. **d,** Standardised effects and 95% intervals for the four X-linked dependencies across cell-line CRISPR, cell-line RNAi, organoid CRISPR and the combined CRISPR analysis. Effects are expressed in residual-standard-deviation units for each source; EIF1AX was unavailable in RNAi and organoids. This panel displays the four prespecified dependencies, rather than a list restricted to significant combined results. **e,** Quantile–quantile plots of combined within-cancer-type P values across 5,511 genes in each of the colorectal and oesophagogastric comparisons. Neither comparison yielded a gene at FDR <0.05. Combination uses test statistics with design-based weights because the original assay scales differ. Organoid fits account for screen noise; cell-line CRISPR and RNAi cohorts overlap, whereas the organoids provide independent models from one study.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_3_F01_organoid_rnai_replication/selection_notes.md)

## supp — Robustness of the Y-loss dependencies to Y-status definition and covariates

**Figure:** `ch03_supp_depmap_dependency_robustness` | **Working ID:** CH3-S3.4.2-S01 | **Results:** 3.4.2 | **Status:** Superseded in Results_edit.md by the combined CH3-S3.4.1-S01 supplement on 8 October 2026; this source caption is retained for historical exports.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_2_S01_depmap_dependency_robustness/CH3_S3_4_2_S01_depmap_dependency_robustness.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_2_S01_depmap_dependency_robustness/CH3_S3_4_2_S01_depmap_dependency_robustness.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_2_S01_depmap_dependency_robustness/CH3_S3_4_2_S01_depmap_dependency_robustness.png)

[PDF](figures/chapter3/ch03_supp_depmap_dependency_robustness.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_dependency_robustness.svg) | [PNG](figures/chapter3/ch03_supp_depmap_dependency_robustness.png)

**Robustness of the four Y-loss dependencies.** **a,** Effect estimate and 95% confidence interval for each gene under each Y-status definition, each covariate set and a fit retaining one model per donor. One robust, donor-clustered fit per row. Definitions change which lines count as Y-negative or Y-positive; covariate sets are additions to lineage. **b,** Effect estimate and 95% confidence interval within each lineage holding at least 10 lines and at least 3 lines per Y group, ordered by size. Effects are Chronos gene effect, Y-negative minus Y-positive. Confidence intervals within lineages are wide and these within-lineage fits are not independent tests of the pan-cancer result.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_2_S01_depmap_dependency_robustness/selection_notes.md)

## supp — Genetic context analyses and within-cancer dependencies

**Figure:** `ch03_supp_depmap_context_screens` | **Working ID:** CH3-S3.4.2-S02 | **Results:** 3.4.2 | **Status:** Superseded in Results_edit.md by the combined CH3-S3.4.1-S01 supplement on 8 October 2026; this source caption is retained for historical exports.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_2_S02_depmap_context_screens/CH3_S3_4_2_S02_depmap_context_screens.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_2_S02_depmap_context_screens/CH3_S3_4_2_S02_depmap_context_screens.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_2_S02_depmap_context_screens/CH3_S3_4_2_S02_depmap_context_screens.png)

[PDF](figures/chapter3/ch03_supp_depmap_context_screens.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_context_screens.svg) | [PNG](figures/chapter3/ch03_supp_depmap_context_screens.png)

**Genetic context analyses and within-cancer dependencies.** **a,** Median detectable genetic interaction at 80% power and nominal P = 0.05 across 86 driver, copy-number, molecular-subtype and genomic-state contexts. Each point is one context. Effects are in Chronos units, per standard deviation for continuous traits. These normal-theory estimates use raw-scale standard errors and provide best-case approximations for the rank-based tests. **b,** ARHGEF2 dependency in 473 CRISPR-screened models, grouped into CIN thirds for display. Points show rank-normalised gene effects after removal of centred lineage and library contributions from the follow-up model; lower values indicate stronger dependency. Blue and red denote Y-positive and Y-negative models; group counts are shown. The FDR = 0.038 refers to the continuous LOY × CIN interaction, corrected within the continuous-genomic class. **c,** Untransformed CIB1 Chronos scores in ten Y-positive and eight LOY bladder/urinary tract models. Points show individual models; the dashed line marks zero. Lower scores indicate stronger dependency. The FDR = 0.030 comes from the original rank-normalised limma screen, corrected across 6,089 genes within this cancer subset; library covariates were constant. Boxes in a–c show medians and interquartile ranges (IQRs), with whiskers extending to the furthest values within 1.5 IQRs. Y status used the primary RNA-led classification incorporating DNA copy-number evidence.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_2_S02_depmap_context_screens/selection_notes.md)

## supp - Supp Fig 3.5. Detection of known vulnerabilities moved into one cancer type.

<a id="supp-fig-3-5"></a>

**Figure:** `ch03_supp_depmap_benchmark_detection` | **Working ID:** CH3-S3.4.2-S03 | **Results:** 3.4.2 | **Chapter number:** Supp Fig 3.5

[Selected PDF](figures/chapter3/selected/CH3_S3_4_2_S03_depmap_benchmark_detection/CH3_S3_4_2_S03_depmap_benchmark_detection.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_2_S03_depmap_benchmark_detection/CH3_S3_4_2_S03_depmap_benchmark_detection.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_2_S03_depmap_benchmark_detection/CH3_S3_4_2_S03_depmap_benchmark_detection.png)

[PDF](figures/chapter3/ch03_supp_depmap_benchmark_detection.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_benchmark_detection.svg) | [PNG](figures/chapter3/ch03_supp_depmap_benchmark_detection.png)

**Supp Fig 3.5. Detection of known vulnerabilities moved into one cancer type.**

**a,** Probability of detecting each benchmark vulnerability when its effect is applied to the Y-negative lines of a cancer-type group and that group is analysed individually across the genome. The display includes 25 benchmarks and 30 evaluated groups across the nested lineage, primary-disease and subtype classifications. Each benchmark was estimated in the eligible male lines, then resampled 400 times with its own lineage levels and noise. Groups of benchmarks denote oncogene addiction, paralog or synthetic-lethal relationships, collateral dependencies, lost dependencies and the four Y-linked paralog dependencies; brackets give the fitted Chronos effect. **b,** Observed results for the four established dependencies in the same groups, distinguishing nominal significance, significance within each group's gene family and significance across all groups at that classification level. An effect is harder to detect where the lineage already depends strongly on the gene. These simulations retain benchmark-specific variability and do not imply equal detection probability for every gene of the same mean effect size.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_2_S03_depmap_benchmark_detection/selection_notes.md)

## supp — Drug-response associations with Y status across seven datasets

**Figure:** `ch03_supp_depmap_drug_response` | **Working ID:** CH3-S3.4.3-S01 | **Results:** 3.4.3 | **Status:** Removed from Results_edit.md at the author’s request on 8 October 2026; this source caption is retained for historical exports.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_3_S01_depmap_drug_response/CH3_S3_4_3_S01_depmap_drug_response.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_3_S01_depmap_drug_response/CH3_S3_4_3_S01_depmap_drug_response.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_3_S01_depmap_drug_response/CH3_S3_4_3_S01_depmap_drug_response.png)

[PDF](figures/chapter3/ch03_supp_depmap_drug_response.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_drug_response.svg) | [PNG](figures/chapter3/ch03_supp_depmap_drug_response.png)

**Drug-response associations with Y status across seven datasets.** **a,** Quantile-quantile plots of observed against expected negative log10 p across 6,413 compound tests, one panel per dataset-screen. Response was modelled on lineage and Y status with donor-clustered robust standard errors. Points on the dotted line are what chance produces. **b,** The twelve smallest p-values across all datasets, with the effect in response standard deviations and the within-dataset FDR. Negative effects denote greater sensitivity in Y-negative lines and positive effects lower sensitivity. Every response metric used is oriented so that lower values mean greater sensitivity, confirmed with BRAF-inhibitor and MDM2-inhibitor controls in each dataset. The datasets share cell lines, use different assays and doses, and predate current culture stocks. No compound met FDR <0.05 after correction within its dataset.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_3_S01_depmap_drug_response/selection_notes.md)

## supp — Y loss and drug response in primary acute myeloid leukaemia

**Figure:** `ch03_supp_beataml_y_loss_drug_response` | **Working ID:** CH3-S3.4.3-S02 | **Results:** 3.4.3 | **Status:** Removed from Results_edit.md at the author’s request on 8 October 2026; this source caption is retained for historical exports.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_3_S02_beataml_y_loss_drug_response/CH3_S3_4_3_S02_beataml_y_loss_drug_response.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_3_S02_beataml_y_loss_drug_response/CH3_S3_4_3_S02_beataml_y_loss_drug_response.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_3_S02_beataml_y_loss_drug_response/CH3_S3_4_3_S02_beataml_y_loss_drug_response.png)

[PDF](figures/chapter3/ch03_supp_beataml_y_loss_drug_response.pdf) | [Editable SVG](figures/chapter3/ch03_supp_beataml_y_loss_drug_response.svg) | [PNG](figures/chapter3/ch03_supp_beataml_y_loss_drug_response.png)

**Y loss and drug response in primary acute myeloid leukaemia.** **a,** Distribution of the expression-based Y-loss fraction before flooring in 263 male BeatAML2 patients, with one sample per patient. Colours indicate recorded karyotypic Y loss, an XY karyotype or unavailable karyotype data; the dashed line marks 0.29, below which the fraction used for association testing was set to zero. **b,** Expression-based Y-loss fraction against karyotype clone size in 14 samples with both measurements. The line denotes equality; karyotyping assesses dividing cells. **c,** Association of Y-loss fraction with rank-normalised drug-response AUC across 95 inhibitors. Positive effects indicate lower sensitivity with greater Y loss. The y axis shows permutation P values; red denotes FDR <0.05 across the inhibitor family. Sample numbers vary by inhibitor and are smaller than the full RNA cohort in a. **d,** Selumetinib AUC against Y-loss fraction in 226 patients, coloured by RAS-pathway mutation status; 11 patients had a Y-loss fraction of at least 0.5. This genotype comparison was added after the association was observed. Drug response is measured for the whole sample and does not isolate the Y-negative clone.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_3_S02_beataml_y_loss_drug_response/selection_notes.md)

## supp — Expression differences between Y-negative and Y-positive cell lines

**Figure:** `ch03_supp_depmap_expression_differences` | **Working ID:** CH3-S3.4.3-S03 | **Results:** 3.4.3 | **Status:** Removed from Results_edit.md at the author’s request on 8 October 2026; this source caption is retained for historical exports.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_3_S03_depmap_expression_differences/CH3_S3_4_3_S03_depmap_expression_differences.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_3_S03_depmap_expression_differences/CH3_S3_4_3_S03_depmap_expression_differences.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_3_S03_depmap_expression_differences/CH3_S3_4_3_S03_depmap_expression_differences.png)

[PDF](figures/chapter3/ch03_supp_depmap_expression_differences.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_expression_differences.svg) | [PNG](figures/chapter3/ch03_supp_depmap_expression_differences.png)

**Expression differences between Y-negative and Y-positive cell lines.** **a,** Volcano plot of differential expression between 321 Y-negative and 447 Y-positive lines, fitted by limma-voom across the 32,190 genes with a known chromosome outside chrY. Chromosome-Y genes define the Y-status call and serve only as a positive control. **b,** Median log2 fold change by dosage class, comparing pseudoautosomal genes with autosomal and non-pseudoautosomal X-linked genes. Loss of one of two copies predicts a fall of 1.0 log2 units. **c,** Hallmark gene sets shifted between the groups. Pathway shifts should be read with care: 87 positional gene sets also pass, which is the signature of co-occurring copy-number change, so these shifts describe the kind of tumour that loses its Y as much as the consequences of the loss. A further 12,937 genes in the count table have no chromosome assignment in this data collection and were excluded from the test family; several behave like chromosome-Y transcripts.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_3_S03_depmap_expression_differences/selection_notes.md)

## supp — Random-forest analyses of Y status, dependencies and drug response

**Figure:** `ch03_supp_depmap_random_forests` | **Working ID:** CH3-S3.4.3-S04 | **Results:** 3.4.3 | **Status:** Superseded in Results_edit.md by the combined CH3-S3.4.1-S01 supplement on 8 October 2026; this source caption is retained for historical exports.

[Selected PDF](figures/chapter3/selected/CH3_S3_4_3_S04_depmap_random_forests/CH3_S3_4_3_S04_depmap_random_forests.pdf) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_3_S04_depmap_random_forests/CH3_S3_4_3_S04_depmap_random_forests.svg) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_3_S04_depmap_random_forests/CH3_S3_4_3_S04_depmap_random_forests.png)

[PDF](figures/chapter3/ch03_supp_depmap_random_forests.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_random_forests.svg) | [PNG](figures/chapter3/ch03_supp_depmap_random_forests.png)

**Random-forest analyses of Y status, dependencies and drug response.** **a,** Out-of-fold R² for nested random forests predicting gene dependency, adding Y status, seven Y-linked expression measurements and wider gene expression to lineage and library. COMPOSITE is the mean of the four main dependencies; grey lines show the other X-linked partners and the context candidate. These prediction models used 10-fold cross-validation repeated three times, with donors kept within folds. **b,** Fraction of subsamples in which real Y status outperformed confounder-matched shuffled versions of itself. The five candidate genes include KDM5C as well as the four genome-wide regression hits; grey bars show 200 random control genes. The dashed line marks a success fraction of 0.8. The model also includes lineage, library, genome doubling, TP53 status, aneuploidy, loss of heterozygosity, ploidy, MSI, chromosomal instability and screen quality. **c,** Features selected when predicting Y status from CRISPR dependencies or from each of seven drug-response datasets. Boruta-style selection compares feature importance with shuffled shadow features, using within-lineage shuffling and 100 iterations per run. Black bars count features selected in both random-seed runs; grey bars show the mean number selected with shuffled Y-status labels. Variable importance has no FDR of its own. The four stable dependency features are supported separately by genome-wide regression; no compound was selected.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_3_S04_depmap_random_forests/selection_notes.md)


### supp — Exploratory PRISM REP.1M LOY drug-response volcano

**LOY-associated drug-response differences in PRISM REP.1M.** Each point represents one of 1,265 compound tests. The x axis shows the regression estimate for LOY relative to Y-positive lines, divided by the response standard deviation; negative values indicate greater sensitivity in LOY and positive values indicate lower sensitivity. The y axis shows −log10 nominal P from the saved lineage-adjusted analysis with donor-clustered robust standard errors. The dashed horizontal line marks nominal P = 0.05, not an FDR threshold. No compound met FDR <0.05 after Benjamini–Hochberg correction within this dataset (minimum FDR = 0.976). Labels identify the four smallest nominal P values in each effect direction; darker points highlight these illustrative compounds and do not denote significant associations. Smaller two-line labels give the compound name and its recorded target, where available, or its recorded mechanism; mechanisms are not evidence of a specific human protein target. No effect-size filter was applied. Exploratory export: `figures/chapter3/ch03_supp_prism_rep1m_volcano.pdf`; not assigned a coded figure number or inserted into the Results.


### supp — Exploratory PRISM REP.1M compound-response box plots

**PRISM REP.1M responses for the eight compounds labelled in the exploratory volcano.** Compounds were selected as the four smallest nominal P values in each direction of the LOY-associated effect. Points represent individual cell lines, with Y-positive lines in blue and LOY lines in red. Boxes show medians and interquartile ranges; whiskers extend to the most extreme observations within 1.5 interquartile ranges. All observations, including outliers, are shown. Responses are the saved PRISM log2 fold changes, with lower values indicating greater sensitivity. Each panel has an independent y-axis scale. Labels give the recorded target or mechanism. P values and FDRs are from the original lineage-adjusted regression with donor-clustered standard errors; the box plots display unadjusted response distributions. None of the displayed compounds passed FDR <0.05 across the 1,265 REP.1M tests. Cell counts and group means were verified against the saved statistical results. Exploratory export: `figures/chapter3/ch03_supp_prism_rep1m_boxplots.pdf`; not assigned a coded figure number or inserted into Results.


## main - Fig 3.5. Drug screens recover expected controls but identify no broadly shared LOY sensitivity.

<a id="fig-3-5"></a>

**Figure:** `ch03_main_depmap_drug_response` | **Working ID:** CH3-S3.3.3-F02 | **Results:** 3.3.1.3 | **Status:** Selected nine-panel figure; incorporates S06 as e and g–i, with SR-1664 moved to f, 8 October 2026. Historical selected folder retains 3.4.3 numbering. | **Chapter number:** Fig 3.5

[Selected PDF](figures/chapter3/selected/CH3_S3_4_3_F02_depmap_drug_response/CH3_S3_4_3_F02_depmap_drug_response.pdf) | [Selected PNG](figures/chapter3/selected/CH3_S3_4_3_F02_depmap_drug_response/CH3_S3_4_3_F02_depmap_drug_response.png) | [Selected SVG](figures/chapter3/selected/CH3_S3_4_3_F02_depmap_drug_response/CH3_S3_4_3_F02_depmap_drug_response.svg)

**Fig 3.5. Drug screens recover expected controls but identify no broadly shared LOY sensitivity.**

**a,** Observed versus expected −log₁₀ *P* for 6,413 compound–dataset tests across seven screens. Expected *P* values are determined from test ranks within each screen, assuming a uniform distribution under the null hypothesis of no LOY-associated response difference. **b–d,** Pan-cancer associations in PRISM REP.1M, PRISM secondary and GDSC2, in response SD units. Dashed lines mark the BH FDR 0.05 reference; none met within-dataset FDR < 0.05. **e,** BRAF/MDM2 controls: the matching inhibitor with the smallest nominal *P* and altered/unaltered counts. Colour indicates −log₁₀ *P*. **f,** SR-1664 responses in mature B-cell neoplasms (PRISM REP.1M; 16 Y-positive/eight LOY, blue/red). Higher values indicate lower sensitivity. Boxes show medians and interquartile ranges (IQR); whiskers reach observations within 1.5 IQR. Δ is the LOY-minus-Y-positive mean difference; rank-normalised OLS *P* and primary-disease-search FDR are shown. **g,** Mixed-model associations for 352 compounds in ≥ 2 screening families; PRISM counts as one family. **h,i,** Mechanism (218) and target (385) classes with ≥ 5 compounds. Effects are mean signed response *t* statistics; *P* uses 2,000 within-lineage permutations. In g–i, dashed lines mark nominal *P* = 0.05; none met within-family FDR < 0.05. Labels identify nominal leads and ATR candidates. Negative effects indicate greater LOY sensitivity. Screens share models and compounds.

[Source and reproduction notes](figures/chapter3/selected/CH3_S3_4_3_F02_depmap_drug_response/selection_notes.md)

## Retired main — Standalone core LOY replication figure

CH3-S3.3.1-F02 was merged into panel f of CH3-S3.3.1-F01 at the author’s request on 7 October 2026. The former exports and caption are preserved in `archive/depmap_f02_into_f01_20261007_153616/`. The active caption is the combined F01 entry above.

## supp — Coverage and detection limits of LOY dependency searches within driver-defined cancer groups

**Figure:** `ch03_supp_depmap_driver_context_detection` | **Working ID:** CH3-S3.3.3-S05 | **Editing copy:** Results 3.3 | **Status:** Removed from Results_edit.md for now at the author’s request on 8 October 2026; caption and exports retained for possible reinstatement.

[PDF](figures/chapter3/ch03_supp_depmap_driver_context_detection.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_driver_context_detection.svg) | [PNG](figures/chapter3/ch03_supp_depmap_driver_context_detection.png)

**Coverage and detection limits of LOY dependency searches within driver-defined cancer groups.** **a,** Numbers of driver-gene/cancer-type pairs meeting each minimum count of driver-altered models per Y group, at lineage, disease and subtype levels. The five-per-group gate yields 31 lineage pairs; only 13 also have five driver-unaltered models per Y group. **b,** Y-group sizes in the 68 distinct eligible comparisons across 22 cancer types and 14 driver genes; 60 have fewer than ten models in the smaller altered Y group. Colours distinguish comparisons that can also test whether the Y effect differs by driver status; open points flag reduced information after accounting for library. **c,** Simulated effect sizes reaching 80% detection, in SD of between-line variation, under the indicated candidate and discovery thresholds. Crosses mark comparisons that did not reach 80% over the evaluated shifts, reflecting limitations of the selected rank-based test. The raw-Chronos comparison is shown for context and has calibration concerns in small groups. **d,** Simulated increases in sample size needed for 80% detection of each effect size, retaining the observed design and current testing thresholds; these are conditional projections, not required sample sizes for a new experiment. **e,** Observed effects for the four established LOY dependencies across the same comparisons. The focused family contains four genes across 68 comparisons; 161/272 tests meet BH FDR <0.05. Recovery is 35/88 (40%) with five lines in the smaller altered group and 30/32 (94%) with at least ten. Comparisons share models and are not independent power measurements. **f,** Observed versus expected P values for 414,052 tests among altered lines and 146,136 tests of a difference by driver status. No test meets its prespecified discovery threshold. The dashed line indicates the first-discovery boundary for the within-altered family. This supplementary figure retains the HPC Figure 14 artwork and annotations. Driver alteration includes qualifying mutations, amplifications or deep deletions; an effect within an altered group does not establish driver specificity.

[Source and reproduction notes](figures/chapter3/source_data/depmap_driver_context_detection/README.md)


## Retired supp — Drug-response controls and pooled tests of LOY-associated sensitivity

**Figure:** `ch03_supp_depmap_drug_controls_pooling` | **Working ID:** CH3-S3.3.2-S06 | **Editing copy:** Results 3.3.2 | **Status:** Merged into CH3-S3.3.3-F02 panels e and g–i on 8 October 2026. The following caption and exports are historical; use the combined main-figure caption above.

[Selected PDF](figures/chapter3/selected/CH3_S3_3_2_S06_depmap_drug_controls_pooling/ch03_supp_depmap_drug_controls_pooling.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_drug_controls_pooling.svg) | [PNG](figures/chapter3/ch03_supp_depmap_drug_controls_pooling.png)

**Drug-response controls and pooled tests of LOY-associated sensitivity.** **a,** Positive-control associations for BRAF alteration with BRAF-inhibitor sensitivity and TP53 alteration with MDM2-inhibitor resistance. Each cell shows the eligible matching inhibitor with the smallest nominal P value in that dataset; n gives genotype-altered/unaltered model counts. Colour indicates −log₁₀ P from lineage-adjusted regression with donor-clustered errors. All 12 selected controls had the expected direction and P <0.01; no eligible control was tested in PRISM REP.300. **b,** LOY-associated response differences for 352 compounds tested in at least two screening families, with PRISM screens treated as one family for eligibility. Responses were standardised within each source and analysed jointly using mixed models adjusted for lineage and source, with a cell-line random intercept. P values are from likelihood-ratio tests. **c,d,** Associations for 218 mechanism-of-action classes (c) and 385 target classes (d), each containing at least five compounds. The x axis shows the mean signed, lineage-adjusted drug-response t statistic; P values use 2,000 within-lineage permutations. Negative values in b–d indicate greater sensitivity in Y-negative models. Dashed lines mark nominal P = 0.05. Labels in b–d identify the four smallest nominal P values, including ties; coincident classes share a label. Five compounds annotated to ATR are additionally labelled in b, alongside the ATR target class in d. ATR mechanism groups contained 1–3 compounds each and were not tested in c. Benjamini–Hochberg correction was applied separately to the pooled-compound, mechanism and target families; none met FDR <0.05. Screening datasets overlap in models, and compounds can belong to multiple classes.

[Source and reproduction notes](figures/chapter3/source_data/depmap_drug_controls_pooling/README.md)


## supp - Supp Fig 3.1. Y-status classification, dependency robustness and genetic context analyses.

<a id="supp-fig-3-1"></a>

**Figure:** `ch03_supp_depmap_combined_checks` | **Working ID:** CH3-S3.4.1-S01 | **Results:** 3.3.1 | **Status:** Selected and inserted into Results_edit.md, 8 October 2026. | **Chapter number:** Supp Fig 3.1

**Combines:** CH3-S3.4.1-S01; CH3-S3.4.2-S01 panel a (panel b omitted); CH3-S3.3.3-S04 panels a and c; CH3-S3.3.2-S02 panels b and c. Former combined panels d (confounder check) and f (approximate interaction detection limits) were removed on 8 October 2026. The combined figure retains provisional ID CH3-S3.4.1-S01; the other source identifiers are retired from Results_edit.md.

[PDF](figures/chapter3/ch03_supp_depmap_combined_checks.pdf) | [Editable SVG](figures/chapter3/ch03_supp_depmap_combined_checks.svg) | [PNG](figures/chapter3/ch03_supp_depmap_combined_checks.png)

[Selected PDF](figures/chapter3/selected/CH3_S3_4_1_S01_depmap_combined_checks/ch03_supp_depmap_combined_checks.pdf)

**Supp Fig 3.1. Y-status classification, dependency robustness and genetic context analyses.**

**a,** Y-status classification of 1,143 male DepMap models using expression of seven Y-linked genes and DNA copy-number evidence. Definite RNA calls take priority; DNA supplies calls when RNA is unavailable. Models with intermediate RNA profiles, partial DNA loss without RNA, neither data type, or discordant RNA–DNA calls are set aside (grey). **b,** Chronos effects (LOY minus Y-positive) across models using alternative LOY definitions, lineage or disease adjustments, additional covariates and donor selection. Error bars show donor-clustered robust 95% confidence intervals. **c,** Random-forest prediction of the four established LOY dependencies (DDX3X, EIF1AX, RPS4X and ZFX) and other tested genes (grey). Predictor groups are added cumulatively: lineage and screening library; Y-status category; seven Y-linked gene expression measurements; and expression of 1,000 other genes selected by strongest absolute correlation with the dependency score within training folds. COMPOSITE averages the four standardised dependency scores. R² is from ten-fold donor-grouped cross-validation repeated three times. **d,** Numbers of dependency or drug-response features stably selected by Boruta with observed LOY status labels (black), versus mean selections after shuffling LOY status labels within cancer lineage (grey). **e,** ARHGEF2 dependency in 473 models, rank-normalised and adjusted for lineage and screening library, grouped into low, middle and high chromosomal instability (CIN) tertiles. The DepMap CIN score measures the fraction of the assessed genome differing from its predominant allele-specific copy-number state. Tertiles are for display; FDR refers to the interaction between LOY status and the continuous CIN score. **f,** Raw CIB1 Chronos scores in bladder/urinary tract models (ten Y-positive/eight LOY); FDR covers the rank-normalised screen of 6,089 genes in this subset. Lower scores indicate stronger dependency. Boxes show medians and interquartile ranges (IQR); whiskers reach observations within 1.5 IQR. Blue/red in a,e,f denote Y-positive/LOY.

[Source data and reproduction notes](figures/chapter3/source_data/depmap_combined_checks/README.md)

## supp — Exploratory RCOR3 dependency comparisons

**Figure:** `ch03_supp_rcor3_exploratory` | **Status:** Exploratory follow-up, 8 October 2026; not assigned a working figure number or inserted into Results.

[PDF](figures/chapter3/ch03_supp_rcor3_exploratory.pdf) | [Editable SVG](figures/chapter3/ch03_supp_rcor3_exploratory.svg) | [PNG](figures/chapter3/ch03_supp_rcor3_exploratory.png)

**Exploratory RCOR3 dependency comparisons.** **a,** Unadjusted RCOR3 Chronos scores in 338 Y-positive and 230 LOY models. The displayed effect (LOY minus Y-positive), P value and FDR come from the saved lineage- and library-adjusted pan-cancer limma analysis with donor blocking and correction across 17,205 genes. MAVER1 is the only model with DepMap RCOR3 dependency probability ≥0.5; it is Y-positive. **b,** Unadjusted RNAi dependency scores in 155 Y-positive and 151 LOY models. Displayed statistics come from the saved lineage-adjusted, donor-blocked analysis, corrected across 16,817 genes. CRISPR and RNAi scores are on different assay scales. **c,** Descriptive RCOR3 Chronos distributions in the 14 lineages containing at least five measured models per Y group. Parentheses give Y-positive/LOY counts; no new within-lineage tests were performed. Blue and red denote Y-positive and LOY models. Points represent individual models; boxes show medians and interquartile ranges, and whiskers extend to the furthest values within 1.5 interquartile ranges. Lower scores indicate stronger dependency; dashed lines mark zero. This follow-up was prompted by RCOR3 having the highest average shadow-win fraction among 200 random non-hit genes (81.7%), while failing the requirement of ≥80% under each of two seeds. RCOR3 did not meet the five-dependent-model eligibility criterion for the selective-gene context screen.

[Source data and reproduction notes](figures/chapter3/source_data/rcor3_exploratory/README.md)


## supp — Nominal dependency associations in driver-altered cancer subsets

**Figure:** `ch03_supp_depmap_driver_nominal_candidates` | **Working ID:** CH3-S3.3.1-S07 (provisional) | **Results:** 3.3.1 | **Status:** Merged into CH3-S3.3.2-F03 panels d and e on 8 October 2026; this caption is retained for historical exports.

[Selected PDF](figures/chapter3/selected/CH3_S3_3_1_S07_depmap_driver_nominal_candidates/ch03_supp_depmap_driver_nominal_candidates.pdf) | [Editable SVG](figures/chapter3/selected/CH3_S3_3_1_S07_depmap_driver_nominal_candidates/ch03_supp_depmap_driver_nominal_candidates.svg) | [PNG](figures/chapter3/selected/CH3_S3_3_1_S07_depmap_driver_nominal_candidates/ch03_supp_depmap_driver_nominal_candidates.png)

**INSM1 and ATRX associations do not pass correction across the full search.** **a,** INSM1 dependency in lung models grouped by MYC status: altered, n = 7/13; unaltered, n = 21/31 (Y-positive/LOY). **b,** ATRX dependency in melanoma grouped by CDKN2A status: altered, n = 20/5; unaltered, n = 9/0. Points show raw Chronos scores for male models (blue, Y-positive; red, LOY); lower values indicate stronger dependency. Boxes show medians and interquartile ranges; whiskers extend to observations within 1.5 interquartile ranges. Brackets compare Y groups; Δ is the library-adjusted Chronos difference (LOY minus Y-positive). P values use library-adjusted limma tests of rank-normalised scores. FDR correction covers 414,052 within-altered tests; the MYC-unaltered P value is unadjusted. The CDKN2A-unaltered comparison was not estimable because it contained no LOY models. Neither post hoc candidate (P <10⁻⁵) passed FDR <0.05.

[Source data and reproduction notes](figures/chapter3/source_data/depmap_driver_nominal_candidates/README.md)
