# Measured HT1080 chromatin sets: RNA enrichment

Completed 14 September 2026. Custom GSEA/CAMERA extension of the [measured reference atlas](CHROMATIN_REFERENCE_ATLAS.md), with complete primary tests, prespecified sensitivities, gene-level evidence and publication figures. Two large-set estimates remain explicitly unestimable in one passage sensitivity after the prescribed retry.

## Main interpretation

**The measured HT1080 reference does not support broad derepression of Polycomb-marked or co-marked promoters in LOY.** The primary GSEA instead identifies a lower-RNA tail among several reference chromatin categories. These are modest, mixed distributions, with important period differences. Neither fixed-correlation nor estimated-correlation CAMERA supports a combined discovery among the 20 custom sets.

This is an association between RNA differences in our established populations and baseline chromatin measured in another laboratory's HT1080 cultures. It does not measure a gain or loss of H3K27me3, PRC2 activity, HDAC activity, or a causal consequence of LOY.

## Primary results

All 20 frozen custom sets were tested in the original 25,461-gene RNA universe. Five have combined GSEA BH q < 0.05, all with negative NES:

| Reference category | Genes | GSEA NES | Custom BH q | Mean RNA difference, log2 CPM |
|---|---:|---:|---:|---:|
| H3K27me3 detected; H3K4me3 not detected | 527 | -1.530 | 0.000384 | -0.0619 |
| All H3K27me3-positive promoters | 845 | -1.368 | 0.00163 | -0.0312 |
| H4K20me1-positive promoters | 770 | -1.346 | 0.00445 | -0.0666 |
| H3K27me3/H2AK119ub1 jointly positive promoters | 318 | -1.445 | 0.00646 | -0.0635 |
| H2AK119ub1-positive promoters | 2,502 | -1.178 | 0.0149 | -0.0356 |

These overlapping categories are not five independent mechanistic findings. The 318 H3K4me3/H3K27me3 **co-marked** genes have NES +0.860, BH q = 1 and mean difference +0.0197. That category is distinct from the 318 H3K27me3/H2AK119ub1 jointly positive genes; equal sizes do not mean identical membership. Acetylated-promoter, H3K9me3/H4K20me3 and lamin categories show no combined GSEA BH discovery under the primary definition. Neither small high-coverage gene-body set passes the primary threshold.

The BY dependence sensitivity retains four of the five GSEA findings; H2AK119ub1 alone no longer passes 0.05. BH/BY values are conditional on valid gene-set P values and do not repair the experimental replication limitation.

**CAMERA gives a weaker conclusion.** No combined set passes custom BH q < 0.05 with either fixed rho = 0.01 or estimated nonnegative within-set correlation. Minimum combined adjusted values are approximately 0.784 and 0.935, respectively. Many adjusted q values are identical because of the BH step-up calculation; the underlying raw P values differ and are exported. For the broad H3K27me3 set, CAMERA's combined direction is Up while GSEA's is Down; neither CAMERA result is significant. GSEA ranks descriptive effect sizes, whereas CAMERA evaluates moderated expression-model statistics relative to the complement with correlation adjustment. Their disagreement is a reason to retain the full distributions and avoid describing uniform chromatin-state regulation.

## Period consistency and gene drivers

The broad H3K27me3 set changes enrichment direction: August NES +1.216, custom BH q = 0.0310; May NES -1.441, q = 0.000097. Its mean differences are +0.00827 and -0.0707 log2 CPM. The combined negative enrichment therefore does not represent a consistent broad shift across both collection periods.

The H3K27me3-detected/H3K4me3-undetected subset is more consistent in enrichment direction: August NES -1.307, q = 0.0310; May NES -1.606, q = 0.000011. Its August mean remains close to zero (-0.00331), compared with -0.121 in May. Negative-tail enrichment is not equivalent to every gene, or even the central portion of the distribution, declining. In fact, **50.1% of this subset and 51.5% of the full H3K27me3 set have positive combined differences**. Their median differences are approximately +0.00015 and +0.00950. The lower-RNA enrichment reflects the weighted tails rather than a majority-wide decrease.

Abundance differs substantially between categories. Median CPM is 0.509 in the H3K27me3-detected/H3K4me3-undetected subset, 1.10 in the full H3K27me3 set and 7.64 in the co-marked set, versus 2.85 across the entire ranking. In the first subset, 76.5% of genes have mean CPM below 2. This motivates the prior-count and abundance/length sensitivities; it does not by itself demonstrate that the signal is technical.

Illustrative leading-edge genes include PLA2G2D, CPNE4, CLSTN2, RYR2 and TMEM63C. Their combined differences range from roughly -2.5 to -3.1 log2 CPM, with lower RNA in both periods. CLSTN2 and RYR2 have mean CPM of 8.71 and 10.95, whereas PLA2G2D and CPNE4 are low abundance (0.46 and 0.27 CPM). The complete leading-edge export includes low-abundance genes, inconsistent genes and exact identifiers; these examples do not define the tested sets. In particular, the unnamed ENSG00000262681 representative is strongly May-weighted and should not be presented as a stable shared-period driver.

No mean, confidence interval or P value treats the genes or six passage comparisons as independently generated LOY populations. Sample/passage differences are descriptive measurements of these same populations.

Post hoc influence diagnostics remove the one or five genes with the largest absolute RNA differences from each set and ranking, then recalculate the unnormalized ES without assigning new P values. All five primary GSEA findings retain negative ES after removing five genes. Those five genes account for approximately 1.9% to 8.3% of total within-set hit weight. The largest chromosome-arm contribution to each leading edge is only 5.9% to 10.0%; this is not a single-arm cluster by simple gene counts. Nevertheless, chromosome-arm sensitivity cannot exclude finer-scale or sample-specific copy-number effects. Among these five sets, only H4K20me1 has a negative *set mean* in every one of the six matched passage comparisons.

## Relationship to the earlier epigenetic findings

The original Mikkelsen enrichments remain valid as results for their original source-defined gene sets and original global multiple-testing family. Their leading edges have limited support from this particular measured HT1080 promoter definition:

| Original source set | Leading-edge genes | HT1080 H3K27me3 promoter overlap | HT1080 H3K4me3/H3K27me3 co-marked overlap |
|---|---:|---:|---:|
| MIKKELSEN_MEF_HCP_WITH_H3K27ME3 | 86 | 10 | 1 |
| MIKKELSEN_MCV6_HCP_WITH_H3K27ME3 | 70 | 6 | 1 |
| PRC2_EED_DN.V1_UP | 46 | 5 | 1 |
| SENESE_HDAC2_TARGETS_DN | 47 | 1 | 0 |

These are descriptive overlaps, not new enrichment tests. **Low absolute coverage is not absence of relative enrichment:** only 3.32% of the full ranked background belongs to the primary H3K27me3 set. The Mikkelsen overlaps of 11.6% and 8.6% are approximately 3.5-fold and 2.6-fold above that background rate. This is partial reference support, without a new significance claim or adjustment for abundance/CpG composition. It does not establish that most leading-edge genes carry the mark in HT1080 or that the mark changes after LOY.

Low coverage may reflect cell context, promoter choice, peak sensitivity and culture differences; it is not proof that the remaining genes lack epigenetic regulation. Under the primary definition, CHRNB2, CHRNA7, MMP1, MMP3 and ZMYND8 have no qualifying promoter peak in this panel. The earlier Mikkelsen and HDAC-response findings therefore should not be relabelled as measured H3K27me3 changes or direct HDAC occupancy in our cells.

## Sensitivity analyses

The configuration was frozen before the new tests. It includes promoter windows of +/-1 and +/-2.5 kb, relaxed peak thresholds, any annotated TSS, exclusion of the 467 unevaluable primary promoters, prior-count and normalization alternatives, six leave-one-matched-passage-out rankings, and descriptive abundance/length/arm residualization. Primary gene-body definitions remain fixed under promoter alternatives. A set below 15 genes is retained as ineligible in the 20-set correction family; it is never truncated or filled with additional genes.

The +/-1 and +/-2.5 kb, any-TSS, evaluable-background, prior-count and RLE comparisons retain the main negative GSEA patterns and no co-marked enrichment. Relaxed peak thresholds retain four primary findings but weaken H2AK119ub1 alone (BH q approximately 0.224). After descriptive abundance/length residualization, the H3K27me3-detected/H3K4me3-undetected subset (q = 0.00648), broad H3K27me3 (q = 0.0269) and H4K20me1 (q = 0.0269) remain below 0.05. H3K27me3/H2AK119ub1 joint promoters weaken to q = 0.0672 and H2AK119ub1 alone to q = 0.234. This shows that some of the association depends on gene properties; the covariate ranking is not an independently validated biological effect.

Adding chromosome arm to the abundance/length residualization retains H3K27me3-detected/H3K4me3-undetected promoters (NES -1.453, q = 0.00262), broad H3K27me3 (NES -1.297, q = 0.0309) and H4K20me1 (NES -1.269, q = 0.0427). H3K27me3/H2AK119ub1 joint promoters have q = 0.0591 and H2AK119ub1 alone q = 0.301. Thus the three retained associations are not explained by these fitted covariates alone, while a global chromatin-activity interpretation remains unsupported.

All five primary findings keep a negative NES in all six leave-one-matched-passage-out rankings. Four remain below custom BH q < 0.05 throughout: the maximum q values are 0.000709 for H3K27me3-detected/H3K4me3-undetected promoters, 0.00460 for broad H3K27me3, 0.00958 for H4K20me1 and 0.0423 for H3K27me3/H2AK119ub1 joint promoters. H2AK119ub1 alone passes in four of six omissions, with maximum q = 0.202. These are robustness checks of the combined ranking, not independent replication; they do not negate the broad H3K27me3 August/May direction reversal. The co-marked set has no GSEA discovery in any scenario (minimum within-scenario q = 0.974).

Each sensitivity reports its own explicitly labelled 20-set BH/BY values. An additional pooled BH spans all 340 sensitivity records for GSEA and 140 per CAMERA correlation mode, including ineligible/unestimable records as unrejected. The combined primary family is kept separate. Under pooled sensitivity BH, the H3K27me3-detected/H3K4me3-undetected and H4K20me1 sets remain below 0.05 with negative NES in every sensitivity; broad H3K27me3 remains below 0.05 but retains its August sign reversal. Sensitivities are not alternative opportunities to select the smallest q value. None of the 318 estimable CAMERA tests passes within-scenario custom BH q < 0.05.

**Numerical limit in one passage sensitivity.** Both large H3K4me3 sets (10,032 and 9,714 genes) remain unestimable when May passage 1 is omitted, despite the prespecified 100,000-permutation retry. Their P values and NES remain missing; both hypotheses stay in the 20-set correction denominator as unrejected. The same two sets required a retry when August passage 2 was omitted, and both resolved to P = 1. All primary and period-specific GSEA estimates are finite. These are large-set directional-null normalization limits, not missing RNA samples or failed primary analyses. The [computational notes](provenance/chromatin_rna/COMPUTATION_NOTES.md) explain the faster per-set retry implementation and deterministic seeds.

Residualization uses four-degree natural splines of mean logCPM and log10 counted-feature length across all ranked genes, optionally adding chromosome arm. It produces a descriptive residual ranking. It can remove biology as well as technical structure and supplies no matched copy-number adjustment. CpG content was unavailable in the frozen covariates, and specimen-relevant DNA copy-number correction remains unresolved. CpG adjustment and matched copy-number correction were not performed.

## Methods and deliverables

Primary sets come directly from the frozen GSE261834 NEED-seq atlas: GENCODE v49, hg38, a selected TSS with a 4,001-bp promoter window, narrow-peak minus-log10(q) >=5 and broad-peak >=2, with ENCODE blacklist bases subtracted. More than 20% promoter masking gives an unevaluable annotation. Co-marking requires both marks within the same TSS window, without claiming same-cell or same-nucleosome bivalency. Lamin categories use deposited peaks, not validated LAD states.

The full original ranking is retained for primary comparability, including genes whose promoters are unevaluable; the 24,994-gene evaluable-promoter sensitivity changes that background explicitly. Broad custom sets are retained in full: 14 exceed the original MSigDB maximum size of 500. Here the maximum is universe size minus one, requiring a nonempty complement; the minimum is 15. Set membership was defined before inspecting the new RNA associations.

GSEA uses fgsea 1.36.2, weighted exponent 1, standard two-sided enrichment, nPermSimple = 10,000, eps = 1e-10 and deterministic seeds derived from 42. A prespecified numerical-failure retry uses nPermSimple = 100,000. Equal RNA metrics are ordered by gene symbol; ties in August/May are retained without artificial jitter. CAMERA uses limma 3.66.0 and the saved voom matrix, precision weights, normalization exposures and sample order. The design has six passage/day block intercepts and separate August/May LOY effects; the combined contrast averages the two effects equally and has four residual degrees of freedom. [fgsea documentation](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html), [limma CAMERA documentation](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf).

The four publication figures provide all-set GSEA/CAMERA results, six prespecified GSEA curves, complete effect distributions and the sensitivity matrix. Each has vector PDF/SVG and 600-dpi PNG exports, with captions and source data:

| Figure | Vector PDF |
|---|---|
| All 20 custom sets, GSEA and both CAMERA variants | [All-set results](outputs/figures/chromatin_rna/HT1080_chromatin_RNA_all_sets.pdf) |
| Six prespecified GSEA curves and leading edges | [GSEA curves](outputs/figures/chromatin_rna/HT1080_chromatin_RNA_GSEA_curves.pdf) |
| Complete RNA effect distributions | [Distributions](outputs/figures/chromatin_rna/HT1080_chromatin_RNA_effect_distributions.pdf) |
| Reference-definition and RNA-covariate sensitivities | [Sensitivity](outputs/figures/chromatin_rna/HT1080_chromatin_RNA_sensitivity.pdf) |

Open the [standalone illustrated report](outputs/report/chromatin_rna_enrichment.html) for all figures and download formats. [Captions](outputs/captions/chromatin_rna.md), the [primary results with effect summaries](outputs/tables/chromatin_rna/primary_results_with_effects.csv), [all GSEA results](outputs/tables/chromatin_rna/fgsea_results.csv), [all CAMERA results](outputs/tables/chromatin_rna/camera_results.csv), [leading-edge gene evidence](outputs/tables/chromatin_rna/leading_edge_gene_evidence.csv), [passage effects](outputs/tables/chromatin_rna/set_passage_effects.csv), [set/leading-edge overlap](outputs/tables/chromatin_rna/set_and_leading_edge_overlap.csv) and [original-source overlap](outputs/tables/chromatin_rna/original_epigenetic_leading_edge_reference_overlap.csv) accompany the plots. The [figure manifest](outputs/figures/chromatin_rna/manifest.csv) records paths, dimensions and hashes.

The run contains 360 GSEA records: 357 estimable, one below the 15-gene minimum (H3K4me1 at +/-1 kb), and the two documented numerical failures. CAMERA contains 320 records: 318 estimable and the same undersized window/set excluded under each correlation mode. Independent checks reproduced all BH/BY and pooled-sensitivity adjustments, all 20 primary ES/leading edges, and 120 primary/period CAMERA results from reconstructed sample blocks, gene-ID indices and explicit voom weights. All 14 frozen input hashes remained unchanged. See [scientific validation](provenance/chromatin_rna/validation.json) and [CAMERA validation](provenance/chromatin_rna/camera_validation.json).

Run `python run_chromatin_rna.py` from this analysis directory. Four isolated R partitions distribute the frozen tests; deterministic test seeds do not depend on worker scheduling. Completed calculations are reused only with unchanged frozen inputs. The serial finalizer merges the partitions, applies pooled sensitivity correction and assembles evidence. Counts, source peaks/bigWigs, original RNA models and original global GSEA/CAMERA results are not modified.

The next biological step is independent HT1080 H3K27me3/active-mark validation at these categories and selected loci, followed by group-level integration of our ATAC effects. Exact RNA-ATAC sample pairing remains unknown, so no sample-wise multimodal correlation is justified.
