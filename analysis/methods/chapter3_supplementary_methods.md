# Chapter 3 — Supplementary Materials and Methods

## S3.1 Molecular profiling sample inventory and genomic assay validation

### S3.1.1 Populations and samples used for molecular profiling

| Cell line and populations | RNA-seq libraries and collection | Other profiling and sampling relationships |
|---|---|---|
| HT1080 sgCTRL and Y-depleted | Six per population: three in batch 1 (August 2025) and three in batch 2 (May 2026, nine months later) | Repeated sampling of one established population per condition. Baseline lpWGS used one sample per population; ATAC-seq used passages 11, 12 and 21 from each population. |
| HT1080 parental | Three in batch 1 | Contextual RNA-seq comparison; excluded from primary sgCTRL-versus-Y-depleted testing. |
| A549 sgCTRL and Y-depleted | Four per population (one batch) | Baseline lpWGS used one sample per population. |
| A549 parental | Three in the same A549 batch | Contextual RNA-seq comparison; excluded from primary sgCTRL-versus-Y-depleted testing. |
| DU145 sgCTRL and twice-nucleofected CP + Q53 + P47 | Three per population at P12, P13 and P14 (one batch) | Repeated sampling of one established population per condition; baseline lpWGS used one sample per population. |
| HT1080 UTY and KDM5D rescue comparison: sgCTRL, Y-depleted, empty vector, UTY, KDM5D and combined UTY–KDM5D | Three per population on matched collection occasions in batch 2 | The six batch 2 sgCTRL and Y-depleted libraries were also used in the baseline comparison. The four transduced populations contributed 12 additional libraries. The combined population descended from the sampled UTY population. |

The 32 baseline and parental libraries and 18 libraries in the rescue comparison therefore represent 44 unique RNA-seq libraries. Library counts do not denote independent population derivations. Original and corrected library identifiers, sampling information and analysis membership are provided in the [RNA-seq sample table](../../tables/chapter3/molecular_samples/rna_library_map.csv).

### S3.1.2 ddPCR validation and descriptive population comparisons

#### Cell-mixture assay check

KP4 and H2452 cells were mixed before DNA extraction at H2452 fractions of 0, 0.25, 0.5, 0.75 and 1. AMELY/AMELX ddPCR used 50 or 100 ng DNA per reaction with technical triplicates. Agreement between measured ratios and expected H2452 fractions was summarized separately by DNA input using Pearson correlation, linear regression with an intercept and mean absolute differences in percentage points. All available measurements were included; technical wells were not independent mixtures or extractions.

#### HT1080 and A549 guide and successive-nucleofection comparisons

HT1080 and A549 copy estimates were grouped by guide combination and nucleofection round and displayed as medians, interquartile ranges and individual measurements. HT1080 combined transfection experiments 1 and 2; the twice-nucleofected sample had one AMELY result and no DDX3Y result. Repeat assays were not independent population derivations, and no inferential tests were performed.

#### DU145 comparisons

DU145 one- and two-nucleofection populations were compared using separate target/reference ratios; KDM5D replaced DDX3Y in the twice-nucleofected assessment. RPS4Y1 was selected to assess distal Yp retention after lpWGS. Changes were described relative to the relevant control.

### S3.1.3 Extended low-pass WGS processing and DU145 regional comparisons

#### Bin-level filtering and numerical floor

In the QDNAseq 100-kb bin analysis of T2T-CHM13v2.0 alignments (Section 3.2.2.4), residual, blacklist, mappability and non-N-base-content filters were disabled, and sex chromosomes were retained. Counts were corrected for GC content and mappability, median-normalized and smoothed to reduce outlier-bin effects. Copy-number profiles were segmented using circular binary segmentation on square-root-transformed values, followed by normalization of the segmented profiles. Genome-wide profiles used segmented log2 relative copy number. Y-region summaries used median back-transformed, unsegmented corrected values from bins assigned by midpoint. A numerical floor of 0.001 was applied before log transformation and labelled in displays; it was not an experimental detection limit.

Baseline lpWGS included one sample per condition, and DNA and RNA were compared at the population level. Centromeric and Yq12 signals were limited by mappability; segment boundaries were not validated breakpoints, and bulk profiles did not resolve individual-cell karyotypes.

#### Y-PAR-masked regional atlases in Figs 3.4–3.6

The HT1080, A549 and DU145 regional atlases used lpWGS realigned with the Y copies of PAR masked. Male-specific Y bins required ≥60% non-N sequence, ≥30% mappability and ≥50 primary MAPQ ≥37 control reads; bins crossing regional boundaries, the centromere or Yq12 were excluded. Combined X/Y PAR coverage mapped to X was normalized to autosomal depth and then sgCTRL. PAR bins required ≥20 control reads and ≥50 reads per 100 kb; PAR1 and PAR2 retained 24 and four bins, respectively. Yq12 was compressed for display, and X-PAR bins were projected into the corresponding HG002 T2T-Y PAR spans. These atlas-specific alignment and display filters differ from the preceding genome-wide processing description. Full processing and dosage-inference details are preserved in the [regional-atlas source notes](../../figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md).

For the DU145 atlas, distal (~2.46–4.80 Mb) and proximal (~4.80–10.57 Mb) male-specific Yp were summarized separately; this division is not a precise structural breakpoint. The reference stock has one X and a rearranged Y ([ATCC](https://www.atcc.org/products/htb-81)), but catalogue cytogenetics do not establish the structure of these experimental populations. An intact XY-to-X model predicts 50% PAR retention. Local sgCTRL X:Yq dosage (~1.8:1) is closer to an XXY-equivalent model, predicting ~67% retention if X dosage is unchanged. Using Yq as the Y-dosage baseline and allowing for the higher X dosage in the depleted pool gives a conditional expectation of ~73–74%. These expectations assume complete loss of Y-derived PAR; bulk lpWGS does not establish the local structural karyotype.

#### DU145 treatment and recurrent regional-change comparisons

Drug-condition profiles were compared with sgCTRL within distal Yp (approximately 2.46–4.80 Mb), proximal Yp (6.00–10.00 Mb) and euchromatic Yq (10.883–27.450 Mb). Treatment-to-control ratios of regional copy-number summaries were multiplied by 100 to express retained signal. The broader male-specific Yp summary used 2.458–10.566 Mb. Region-specific ratios were not interpreted as absolute cell fractions.

Seven edited DU145 profiles were compared with the same sgCTRL using arm-median segment-log2 differences across 100-kb bins assigned by CHM13 midpoint. The fraction of bins with differences ≤−0.2 summarized lower dosage. The low-depth KU-only profile was displayed but excluded from the six-profile recurrence summary. Shared controls and unresolved ancestry limited interpretation as independent recurrent events.

## S3.2 Extended RNA-seq quality control, PCA, and sensitivity analyses

### S3.2.1 Provider RNA integrity and library quality control

Provider QC supplied RNA integrity measures: RINe for A549 and batch 1 HT1080, RQN for batch 2 HT1080 and RIN for DU145. RNA size summaries used DV200, except batch 2 HT1080, which used the fraction spanning 200–10,000 nucleotides.

For the baseline analyses, gene-level read counts were generated using featureCounts v2.1.1 with the GENCODE v49 primary-assembly annotation, assigning reads to exons grouped by gene identifier. Counting was unstranded. Paired-end input was specified, but reads were counted individually because paired-fragment counting was not enabled. Regional Y-expression analyses used GRCh38 gene counts, in which pseudoautosomal genes were assigned to chromosome X. These counts represent combined X- and Y-derived expression and cannot distinguish the chromosome of origin. HG002 T2T-Y coordinates were used for genomic illustrations only.

QC across all 44 libraries used MultiQC v1.34, samtools v1.23.1, Qualimap v2.3 and Picard CollectRnaSeqMetrics v3.1.0. Metrics covered read quality and depth, mapping, duplication, genomic read distribution, gene-body coverage, mitochondrial reads and rRNA. Qualimap used paired-end, non-strand-specific mode with GENCODE v49; Picard used STRAND=NONE.

Gene-body coverage across 100 relative positions was divided by each library’s mean coverage. Mitochondrial abundance was the percentage of mapped alignment records on chrM; rRNA abundance was the percentage of aligned bases in rRNA intervals. All 44 libraries were retained and displayed in **Supp Fig** 3.1 and **Supp Fig** 3.2.

[![Supp Fig 3.1. RNA-seq quality control across 44 in-house libraries.](../../figures/chapter3/ch03_methods_rna_library_qc.png)](../../figures/chapter3/ch03_methods_rna_library_qc.pdf)

**Supp Fig 3.1. RNA-seq quality control across 44 in-house libraries.**

Libraries comprise A549 (11), batch 1 HT1080 (9), batch 2 HT1080 including UTY/KDM5D rescue populations (18), and DU145 (6), ordered identically across panels. **a,** Trimmed read pairs entering STAR. **b,** STAR alignment classes as percentages of input pairs. **c,** Qualimap exonic, intronic and intergenic fractions of classified reads. **d,** fastp duplication estimates. **e,** Mitochondrial alignments as a percentage of mapped alignment records (blue), and rRNA-associated bases as a percentage of Picard aligned bases (orange); denominators differ. **f,** Qualimap 5′/3′ coverage bias; the dashed line denotes equal coverage. Condition colours apply in a,d,f. P, parental; C, sgCTRL; Yd, Y-depleted; EV, empty vector; U, UTY; K, KDM5D; UK, combined UTY–KDM5D rescue. Indices identify samples or collection occasions, not independent population derivations. All 44 libraries are shown.

<!-- Caption source: ../../Chapter3_figure_captions.md; figure: ch03_methods_rna_library_qc. -->

[![Supp Fig 3.2. Read quality, input RNA and gene-body coverage across 44 RNA-seq libraries.](../../figures/chapter3/ch03_methods_rna_input_quality_coverage.png)](../../figures/chapter3/ch03_methods_rna_input_quality_coverage.pdf)

**Supp Fig 3.2. Read quality, input RNA and gene-body coverage across 44 RNA-seq libraries.**

All 44 libraries, sample order and condition abbreviations follow **Supp Fig** 3.1. **a,** Bases with Phred scores ≥ 30 before (grey) and after (blue) fastp processing. **b,c,** Adapter-trimmed (b) and retained reads (c), each as a percentage of raw reads; trimmed reads can be retained. **d,** Input RNA integrity: RINe for A549/batch 1 HT1080, RQN for batch 2 HT1080, and provider-reported RIN for DU145. **e,** DV200 (percentage of RNA >200 nt), except for batch 2 HT1080, which shows the fraction within 200–10,000 nt. **f,** Qualimap coverage across 100 relative gene-body positions, divided by each library's mean. Curves represent libraries; dashed lines mark mean-normalized coverage of one.

<!-- Caption source: ../../Chapter3_figure_captions.md; figure: ch03_methods_rna_input_quality_coverage. -->

### S3.2.2 Extended baseline RNA-seq normalization, PCA, and network details

#### Normalization and additional feature sets

TMM normalization used expressed autosomal and non-PAR X-linked genes: 25,526 for HT1080, 22,792 for A549 and 18,767 for DU145. A549 filtering used edgeR::filterByExpr with minimum count 10 and minimum total count 15. Its primary test family additionally included 20 non-PAR Y-linked genes and 31 PAR features, giving 22,843 nuclear genes. CPM for additional genes used the primary effective library sizes. HT1080 Y-linked, PAR and mitochondrial genes were assessed separately using the same effective library sizes.

For HT1080, robust empirical Bayes moderation omitted an additional mean–variance trend. Three parental libraries per line provided context for HT1080 and A549 and were excluded from primary normalization and testing. HT1080 parental samples were collected in batch 1. DU145 retained legacy LNCaP identifiers in source records. A549 descriptive effects were differences between group-mean log2CPM values, distinct from the RINe-adjusted model coefficients. Complete primary results and model provenance are linked in the [gene-level results tables](../../tables/chapter3/a549_primary_nuclear_RINe/README.md).

#### Principal-component analysis and quality-control associations

PCA used TMM-normalized log2 CPM, centred by gene without variance scaling. HT1080 and A549 used prior count 2; their core PCAs retained 25,526 and 22,792 autosomal/non-PAR X genes, respectively. DU145 PCA used its primary autosomal/non-PAR X gene set. A549 PCA was not adjusted for RINe. Parental-inclusive PCAs used separate TMM normalizations across all 15 HT1080 or 11 A549 libraries while retaining the corresponding core gene sets.

The HT1080 all-library PCA used limma::removeBatchEffect to remove an additive batch effect while preserving population, assuming a common batch effect for parental cells sampled only in batch 1. This adjustment was for visualization. Within-batch PCAs used unadjusted expression with the all-15-library TMM factors and gene set. Core-library correlations used unadjusted expression, with average-linkage clustering on 1 minus Pearson r.

A549 PC1 associations with RINe, exonic fraction and mitochondrial fraction were summarized by Spearman correlation across eight primary libraries. Exonic fraction was the exonic count divided by the sum of exonic, intronic and intergenic Qualimap counts; mitochondrial fraction used mitochondrial gene counts divided by total assigned counts.

#### Enrichment and protein-association networks

A549 pathway analyses used unadjusted rankings and a condition-only model, without the RINe covariate used in its primary gene-level analysis. For the two-batch HT1080 comparison, consistent programs required overall FDR <0.05 and matching normalized enrichment-score directions across both batches.

HT1080 STRING v12.0 networks used scores ≥700 and genes with |mean log2CPM difference| ≥0.5, matching directions across batches and all leave-one-pair-out calculations, without a significance filter. Higher- and lower-expression networks were clustered separately using MCODE. Metascape analyzed the stated HT1080 and DU145 differential-expression lists against expressed-gene backgrounds, using hypergeometric tests and Benjamini–Hochberg correction.

## S3.3 Extended ATAC-seq quality control, normalization sensitivities, and external dataset reanalyses

### S3.3.1 HT1080 ATAC-seq quality control and sensitivity analyses

#### Quality control metrics

ATAC-seq QC assessed fragment counts and lengths, TSS enrichment and FRiP. Whole fragments were used for peak calling and counting; +4/−5 Tn5 offsets were applied only to insertion-based QC and footprinting. TSS enrichment compared central ±50-bp insertion coverage with the outer 100-bp flanks of a ±2-kb window. FRiP counted unique fragments overlapping sample or consensus peaks. Depth-matched QC used 13,595,988 fragments per library; differential analysis used full-depth counts from all six libraries.

All six libraries showed nucleosomal fragment-length patterns and enrichment at transcription start sites ([**Supp Fig** 3.3a,b,e](../../figures/chapter3/ch03_methods_ht1080_atac_qc.pdf)). Library depth, FRiP and depth-matched peak counts are shown in **Supp Fig** 3.3c; between-library accessibility correlations are shown in **Supp Fig** 3.3d.

[![Supp Fig 3.3. HT1080 ATAC-seq library quality.](../../figures/chapter3/ch03_methods_ht1080_atac_qc.png)](../../figures/chapter3/ch03_methods_ht1080_atac_qc.pdf)

**Supp Fig 3.3. HT1080 ATAC-seq library quality.**

Six libraries from one sgCTRL and one Y-depleted population, each sampled at passages 11, 12 and 21. **a,** Fragment-length distributions. **b,** Tn5 insertion profiles around transcription start sites (TSSs), normalized to flanks; enrichment compares the central 100 bp with outer flanks. **c,** Fragment counts, fraction overlapping each library's own peaks (FRiP), and peak counts before/after downsampling to 13.6 million fragments. **d,** Spearman correlations across 84,510 regions, excluding Y and PAR. **e,** Coverage around 34,174 TSSs per million non-Y fragments, with shared row order and colour scale capped at the 99th percentile. Blue/red denote sgCTRL (C)/Y-depleted (L); solid/dashed/dotted lines denote passages 11/12/21. Passage samples are repeated measurements of established populations.

<!-- Caption source: ../../Chapter3_figure_captions.md; figure: ch03_methods_ht1080_atac_qc. -->

#### Covariate, normalization, and genomic-context sensitivities

A sensitivity analysis additionally included passage as a categorical covariate (~ passage + condition), using the same 84,510 regions and TMM normalization; the condition-only model remained primary. PCA used log2 CPM with prior count 2, centred by region without variance scaling. Similarity between samples was assessed using Spearman correlations.

Normalization sensitivities used non-peak 10-kb bins with at least 30 counts across libraries, or csaw v1.44.0 loess offsets for the abundance–fold-change relationship. Differential testing and correction were repeated, with loess-adjusted results also used for pathway and promoter comparisons. Arm-level ATAC and RNA effects were compared with GRCh38-aligned HT1080 D0 lpWGS, excluding centromeres.

#### Motif, TSS-coverage, and regulatory-domain sensitivity details

In HOMER v5.1 motif testing (Section 3.2.3.3), backgrounds had FDR ≥0.20 and were matched by GC content and baseline logCPM, with up to ten controls per target without replacement. Benjamini–Hochberg correction covered motifs in both directions. Separate promoter and distal analyses used ATAC FDR <0.10, with correction within each region class. In chromVAR v1.32.0 and motifmatchr v1.32.0 analyses across 84,506 regions with usable sequence, the motif threshold was 5 × 10⁻⁵ with 50 background iterations, and shared motifs were interpreted at the transcription-factor-family level.

For RNA–ATAC promoter integration, ATAC FDR <0.05 and loess-adjusted effects were evaluated as sensitivities alongside the primary ATAC FDR <0.10 threshold. A copy-number sensitivity retained promoters with |Y-depleted-minus-sgCTRL log2 copy ratio| <0.2.

TSS coverage used strand-oriented 50-bp bins within ±3 kb, normalized per million non-Y fragments. Each gene’s TSS was selected by maximal pooled central ±100-bp coverage. RNA-defined groups were compared with three non-DE controls per gene without differential promoters, sampled without replacement within pooled TSS-coverage deciles; 200 additional draws assessed matching sensitivity. Per-gene log2 ratios used pseudocount 0.1. Two-sided Mann–Whitney tests compared target and control distributions, with Benjamini–Hochberg correction across four comparisons.

In the joint gene-program and Polycomb/rGREAT analyses, fgseaMultilevel assessed Hallmark, Reactome and GO Biological Process enrichment using ranked promoter effects, promoter test statistics and RNA effects in the same integrated gene universe. Benjamini–Hochberg correction was applied within each collection and ranking for fgsea, and separately to region- and gene-level tests within each direction for rGREAT v2.12.2. Sensitivity analyses used alternative matching, loess-adjusted accessibility and alternative region-to-gene assignments.

### S3.3.2 Extended regional DNA–RNA, DepMap co-expression, and cross-line sensitivity analyses

#### Regional expression and copy-number window definitions

Expressed genes were mapped to CHM13v2.0 by stable Ensembl identifier or unique same-chromosome symbol. DNA–RNA correlations described genomic units across the profiled populations.

For HT1080, mean log2CPM differences were averaged across 1-Mb windows and chromosome arms, weighting collection batches equally. Windows required at least five mapped genes and arms at least 20. DNA-window effects were differences between population medians of corrected unsegmented log2 values across usable 100-kb bins; arm effects were median paired-bin differences excluding centromeres. Pearson correlations covered 1,583 non-Y windows and 41 arms.

For A549 and DU145, gene-level DNA differences were calculated from median segmented log2 copy-number values in overlapping 100-kb bins and compared with Y-depleted-minus-sgCTRL mean log2CPM differences. Regional comparisons summarized mean RNA differences and median bin-level DNA differences within 10-Mb windows containing at least ten mapped genes and chromosome arms containing at least 20. Pearson correlations were calculated across non-Y genes, windows and arms. Additional regional comparisons included male-specific Yp and Yq, excluding pseudoautosomal regions, the centromere and Yq12.

#### DepMap GRHL2 co-expression overlap

The 200 genes most positively correlated with GRHL2 expression across 1,719 DepMap 26Q1 cell lines were selected using Pearson correlations of log2(TPM + 1), excluding GRHL2 itself. Over-representation among DU145 genes with lower expression at FDR <0.05 was assessed using a one-sided Fisher’s exact test against the 13,843 genes represented in both datasets. The reported P value was nominal.

#### Cross-line expression and abundance sensitivity analyses

In the cross-line comparison of HT1080, A549 and DU145 (Section 3.2.3.6), the saved comparison retained the earlier condition-only A549 results; these are distinct from the RINe-adjusted primary gene tests. The common effect ranking (mean log₂ CPM difference) also differs from the moderated-t ranking used for the primary HT1080 program analysis. Gene-set FDRs retained correction across all eligible sets within each background, including when displaying selected sets or summarizing Hallmark, Reactome and GO:BP alone. Shared pathway themes were assigned from leading-edge gene membership.

Post hoc sensitivity analyses required mean CPM ≥2, ≥5 or ≥10 in every background, using sgCTRL expression or the average across both conditions. Effect correlations, significant-gene overlaps and concordant large-effect counts were repeated, with nominal permutation tests. A separate sensitivity analysis expanded the interferon/inflammatory theme to include TNF-alpha/NF-κB and IL-6/JAK/STAT3 genes. [Complete sensitivity results and reproduction details](../../research/molecular_characterisation/crossline_abundance_sensitivity_2026-09-27/README.md) provide the full settings.

### S3.3.3 Reanalysis of external Y-loss and Y-paralog expression datasets

Eighteen external contrasts were reanalyzed from public reads or expression matrices using study-specific designs. Sources, samples and processing are listed in the [external dataset table](../../tables/chapter3/Table_CH3_external_expression_datasets.md). Several contrasts shared samples or came from the same study.

#### Bulk RNA-seq from constitutional karyotype and Y-linked knockout studies

Raw bulk RNA-seq reads were downloaded for [GSE126712](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126712), [GSE207112](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207112), [GSE264742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264742), [ERP156313/E-MTAB-13673](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13673) and [GSE301297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301297). Reads were aligned using two-pass STAR, and gene counts were generated using featureCounts v2.1.1 with GENCODE v49 and PAR-aware annotation. Sample metadata were joined to run accessions and study sample annotations. Quality assessment included read assignment, library size, sex-chromosome marker expression, PCA and sample clustering. Low-expression filtering used edgeR filterByExpr, followed by TMM normalization and voom-limma modelling. Each dataset was analyzed separately using a study-specific design.

Primary lymphocytes from [Zhang et al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32071206/) were modelled by karyotype, extracting the XO-minus-XY contrast (13 libraries in each group). The trophoblast comparison from [Ahern et al. (2022)](https://pubmed.ncbi.nlm.nih.gov/36161909/) used a karyotype-plus-differentiation-round design, excluding the recorded PCA outlier XO9_TBL_3. The neural-crest study ([Ahern et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11947172/)) compared XO and XY in the two male donor panels, adjusting for donor panel and differentiation round.

For the fetal-tissue study ([Suntharalingham et al., 2025](https://doi.org/10.1038/s42003-025-07699-4)), a tissue-adjusted XO-minus-XY contrast was fitted across placenta, kidney, liver, muscle, pancreas and skin using karyotype and tissue terms. The UTX/UTY knockout study ([Akiyama et al., 2026](https://pubmed.ncbi.nlm.nih.gov/41906541/)) contributed undifferentiated SEES3 UTY-knockout-versus-WT and double-knockout-versus-WT contrasts, each with two libraries per group. These comparisons assess loss of individual paralogs rather than chromosome loss.

#### A549 clone and MB49 comparisons

The A549 data from [Schlüter et al. (2026, preprint)](https://doi.org/10.64898/2026.06.02.729627) were reanalyzed from released count data. Bulk RSEM estimated counts were extracted from manuscript Table 3, rounded to integers and analyzed with DESeq2 using an unpaired LOY-versus-ROY design (four clones per group). Expression filtering used filterByExpr with min.total.count = 200; chrY was excluded from the downstream DE analysis. The multiome comparison used RNA UMI counts aggregated by clone from the released A549 single-cell 10X multiome object, followed by the corresponding DESeq2 comparison. These are two assays of the same clone panel. The source study lists [figshare](https://figshare.com/s/0bb9ea2ee2e0ea53f74c) for A549 multiomic objects and [GHGAS27424019412227](https://data.ghga.de/study/GHGAS27424019412227) for controlled-access raw data.

The MB49 study ([Abdel-Hafiz et al., 2023](https://doi.org/10.1038/s41586-023-06234-x)), deposited as [GSE229233](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229233), contributed separate cell-line and tumour comparisons, each with five Y-negative and five Y-positive samples. The external comparison imported local DESeq2 effect tables. Both locally generated STAR/featureCounts counts and reanalysis of author-supplied RSEM counts are present in the processing records. For cross-species gene-set matching, mouse symbols were upper-cased and matched to human symbols; this was a symbol-based approximation rather than a curated orthology map.

#### Released Page-lab expression matrices

Fibroblast and lymphoblastoid comparisons from [San Roman et al. (2024)](https://doi.org/10.1016/j.xgen.2023.100462) were reanalyzed using released gene-count matrices from [Zenodo](https://doi.org/10.5281/zenodo.10042162). Counts from the authors' kallisto/tximport workflow were rounded for DESeq2. The X0-versus-XY comparisons included 23 versus 14 fibroblast samples and 17 versus 17 lymphoblastoid samples. Both fitted designs included library-preparation batch and karyotype group.

#### Naturally Y-mosaic single-cell datasets

Single-cell RNA-seq (scRNA-seq) data from the HT1080 spheroid study ([Olson et al., 2026](https://pubmed.ncbi.nlm.nih.gov/41961858/)), [GSE313572](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313572), were reanalyzed from released count matrices. Only the control sample GSM9370527 contributed to the external comparison. The local configuration specified 3,000–12,000 detected features and a mitochondrial fraction no greater than 17.5%. Y status was assigned using RPS4Y1, EIF1AY, DDX3Y, KDM5D, UTY, USP9Y, ZFY and TMSB4Y. Counts from 1,185 Y-negative and 2,142 Y-positive cells were aggregated within four random partitions and fitted with TMM voom-limma using partition and Y-status terms. These partitions are pseudoreplicates of one sample; gene-level P values from that fit were not used in the downstream AUC test.

Tahoe-100M ([Zhang et al., 2025, preprint](https://doi.org/10.1101/2025.02.20.639398)) data were obtained from the public [Tahoe-100M release](https://huggingface.co/datasets/tahoebio/Tahoe-100M) and reanalyzed for naturally Y-mosaic cultures. The external comparison retained A549, J82, SK-MEL-2 and HCT15 DMSO contrasts under the recorded inclusion rule (MOSAIC call and PAR median log2 difference ≤ −0.5). The analysis used strict Y-status labels and required at least 3,000 UMIs per cell and at least 25 cells per Y-status group per well. Within-well pseudobulk comparisons included 22, 22, 19 and 6 wells, respectively.

### S3.3.4 Reanalysis of other chromosome-loss models

#### Datasets and experimental units

Four public datasets were analyzed separately. [GSE336056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336056) comprised RPTEC single-cell RNA-seq (scRNA-seq) and hashtag-oligonucleotide (HTO) data after sgNTC, sgYq or sg3p treatment at days 5 and 10 after editing. [GSE336576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336576) comprised bulk RNA-seq data from three RPTEC 3p-loss clones (C5, C7 and C10), each sampled in triplicate at day 18 and day 90 after editing, with three unedited control libraries. Both datasets came from [Dahiya et al. (2026, preprint)](https://doi.org/10.64898/2026.07.10.737056). [GSE222379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222379) contributed bulk RNA-seq data for nine contrasts after ReDACT-mediated removal of extra 1q, 7p or 8q copies in A2058, A2780, AGS, HCT116, MCF10A and MCF10A-HRAS backgrounds ([Girish et al., 2023](https://doi.org/10.1126/science.adg4521)). [GSE150686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150686) contributed bulk RNA-seq data for RPE1 monosomy 10+18 and monosomy 13 comparisons with p53-knockout controls, and a 19p-loss comparison with p53-knockdown controls ([Chunduri et al., 2021](https://doi.org/10.1038/s41467-021-25288-x)). Thus, the bulk RNA-seq collection contained 14 contrasts. [Dataset and analysis tables](../../tables/chapter3/arm_loss_comparison/README.md) specify sample membership, aggregation units and shared controls.

#### Single-cell loss classification and differential expression

HTO-assigned singlets with at least 10,000 UMIs were retained. Within the sgYq group, Yq-null cells had zero Yq-panel UMIs and at least five Yp-panel UMIs; Yq-intact cells had at least two Yq-panel UMIs and at least five Yp-panel UMIs. Other cells were excluded from this comparison. Within sg3p, loss calls required a 3p RNA log₂ ratio below log₂(0.8) and an absolute 3q ratio below 0.25, relative to same-day non-targeting controls; loss-called cells were compared with the remaining sg3p-treated cells. These were RNA-based classifications, not independent DNA measurements. [Classification and regional-expression source tables](../../figures/chapter3/source_data/arm_loss_comparison/README.md) retain the supporting checks. Differential expression used MAST v1.36.0, adjusting for detected-gene count and HTO sample, with likelihood-ratio tests and Benjamini–Hochberg correction. Pooled analyses were primary, supplemented by day-specific fits. Cells were subsamples of the same source cultures and did not represent independent experiments. Enrichment ranked genes by MAST effects, excluding Y, pseudoautosomal and unannotated genes, plus 3p genes from the 3p comparison. Full settings are provided in the [Yq analysis](../../../RNAseq/analysis/arm_loss_bundle/code/gse336056/de_yq_null_vs_intact.R), [3p classification](../../../RNAseq/analysis/arm_loss_bundle/code/gse336056/chr3p_scan.R) and [3p analysis](../../../RNAseq/analysis/arm_loss_bundle/code/gse336056/de_3p_loss_vs_intact.R) scripts.

#### Bulk differential expression and gene-set comparisons

Bulk analyses used deposited gene-count tables, with Ensembl identifiers standardized and duplicated PAR-Y entries removed. ReDACT libraries were summed within clone and RPTEC libraries within clone and time point; controls and RPE1 samples were retained separately. RPE1 samples represented repeated sampling of established lines. DESeq2 fitted an unpaired group-only model for each contrast, retaining genes with at least ten total counts. Unshrunken fold changes and Wald statistics were retained, with Benjamini–Hochberg correction. Early and late RPTEC contrasts shared controls. Bulk enrichment ranked genes by Wald statistics after excluding genes in the lost arm or chromosome; full settings are provided in the [bulk analysis script](../../../RNAseq/analysis/arm_loss_bundle/code/arm_loss_meta/02_de_gsea.R). Gene-set FDR correction was performed separately within each collection and contrast. The day-specific RPTEC display used the saved Hallmark results, corrected across all eligible sets at each time point (49 at days 5 and 10; 50 at days 18 and 90). These separate contrasts were not a longitudinal or condition-by-time interaction test.

## S3.4 Extended HT1080 UTY and KDM5D rescue analytical methods

### S3.4.1 Frozen 100-gene LOY signature and rank-scoring details

Alongside the 150-gene primary LOY signature (Section 3.2.4), the originally frozen 100-gene signature used the same batch 1 passage-pair consistency and |mean log2CPM difference| ≥0.5 rules without the FDR <0.20 requirement, taking the 50 largest effects in each direction. All 150 primary-signature genes, and 75 of the original 100 genes, passed the batch 2 primary expression filter; all 100 were retained for rank scoring.

Batch 2 rescue libraries were scored using within-sample percentile ranks across the 23,903-gene discovery background, following the single-sample rank-based approach of singscore (Foroutan et al., 2018), calculated directly on the frozen background. Scores were the mean percentile rank of higher-in-LOY genes minus that of lower-in-LOY genes. Score contrasts were fitted with collection-occasion adjustment. Percent recovery was scaled from the empty-vector mean (0%) to sgCTRL (100%). The 95% confidence intervals were estimated by sampling from the fitted model and recalculating percent recovery.

### S3.4.2 Matched control selection for promoter-class comparisons

Among 225 down-regulated genes represented in the rescue analysis (RNA FDR <0.05, |log₂ fold change| ≥1, and matching batch directions), 58 had at least one lost promoter peak at ATAC FDR <0.10 and 167 had no detected promoter loss. This expanded set was distinct from the 39 genes meeting FDR <0.05 in both assays. For the matched display, each of the 58 genes was paired with two controls using Euclidean distance in standardized baseline RNA log₂ fold change and rescue mean expression. Matching was greedy without replacement, with gene order randomized using seed 20260922, yielding 116 controls shared across contrasts. [Source tables and reproduction details](../../figures/chapter3/source_data/addback_atac_overlap/README.md) provide the full integration and matching settings.
