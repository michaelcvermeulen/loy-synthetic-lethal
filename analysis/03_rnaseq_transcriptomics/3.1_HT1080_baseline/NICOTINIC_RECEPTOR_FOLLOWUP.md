# Nicotinic receptor follow-up in HT1080 LOY

The current data support selective remodeling of nicotinic-subunit RNA abundance, with clearer increases in CHRNA7 and CHRNB2, low-abundance increases in CHRNA6 and CHRNB4, and decreases in CHRNA1, CHRNA9 and CHRNA10. They do not establish increased receptor activity or one five-subunit receptor complex.

| Gene | Mean CPM | August Δlog2 CPM | May Δlog2 CPM | Combined Δlog2 CPM | Same-direction matches / 6 |
|---|---|---|---|---|---|
| CHRNA6 | 0.41 | 0.28 | 2.07 | 1.18 | 5 |
| CHRNA7 | 3.30 | 0.38 | 1.47 | 0.93 | 5 |
| CHRNB2 | 5.89 | 0.60 | 1.67 | 1.13 | 5 |
| CHRNB4 | 0.28 | 1.31 | 2.94 | 2.12 | 6 |
| CHRNE | 0.93 | 0.40 | 0.63 | 0.52 | 4 |

The five-member UP_01 module has MCODE score 5, ten association edges and ten edges above the STRING physical-score threshold of 700. Its top annotation, Reactome acetylcholine binding and downstream events, has raw P = 8.47e−17 and global annotation BH q = 5.51e−12. This quantifies functional coherence of the selected genes. It is not a test of RNA change or evidence that all five products assemble together in HT1080. STRING physical edges can encode co-complex evidence and transferred evidence; pair-specific experimental provenance is not resolved by the cached combined-score files. [STRING scoring documentation](https://version-12-0.string-db.org/help/scores/).

## Abundance and temporal consistency

CHRNA7 and CHRNB2 have appreciably more assigned reads than CHRNA6 and CHRNB4. In May, CHRNA7 averages 2.12 CPM in sgCTRL and 5.99 in LOY; CHRNB2 averages 3.39 and 10.80 CPM. Both are higher in five of six passage/day comparisons and remain positive after leaving out any one match. Their effect estimates change by about 0.01 log2 units with the smaller prior count or alternative RLE normalization.

CHRNA6 averages 0.41 CPM overall and includes zero-read controls. CHRNB4 averages 0.28 CPM and is higher in all six matches. Its May controls have 0–4 reads, while LOY samples have 12–22. All five genes remain eligible at the saved 50- and 100-total-read sensitivity thresholds. However, reducing the logCPM prior count from 2 to 0.25 changes the combined CHRNA6 estimate from +1.18 to +1.95, and CHRNB4 from +2.12 to +3.29. Direction is more stable than the exact magnitude at these low counts. Prior-count effects are not confidence intervals.

The mean module contrast is +0.60 in August and +1.75 in May. This is descriptive temporal variation. The larger CHRNB4 May contrast partly reflects lower sgCTRL abundance: its mean LOY abundance is 0.52 CPM in August and 0.40 in May, while controls fall from 0.17 to 0.03 CPM. A larger contrast therefore does not imply increasing LOY abundance over time. The same established populations were sampled in both periods; the six matches are not six independently generated populations.

## The rest of the receptor family changes in different directions

CHRNA1 is lower by −1.33 log2 CPM and CHRNA9 by −1.21, both in all six matches. CHRNA10 is lower by −0.62, with a stronger August effect and near-zero May mean. CHRNA3 and CHRNA5 show small positive combined means but opposite period signs. CHRNA4 is below the expression filter, and CHRNB3 has no assigned reads. The full family, including excluded low-expression genes, is retained in the figure and source tables. All-zero genes are displayed as having no assigned reads; tiny floating-point residuals in the original effect table are not interpreted biologically.

## What could the subunits assemble into?

| Candidate assembly | Relevant observations | Interpretation |
|---|---|---|
| α7-containing receptor | CHRNA7 +0.93; CHRFAM7A +0.51; RIC3 and TMEM35A below the expression filter | A testable receptor hypothesis, with assembly and gene-specific measurement unresolved. |
| α6β4-containing receptor | CHRNA6 +1.18 and CHRNB4 +2.12, both low abundance | A plausible pairing; not evidence that a functional channel exists in these cells. |
| β2-containing receptor | CHRNB2 +1.13; CHRNA6 increases, while CHRNA4 is very low and CHRNB3 has no reads | The RNA data do not identify the exact β2-containing receptor subtype. |
| Muscle-type receptor | CHRNE +0.52, CHRNA1 −1.33, CHRNB1 +0.21 and CHRND +0.43 | The muscle-type subunits do not show coordinated induction. |

α7 can form a homopentamer ([structural study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8135066/)); human α6β4 has been structurally characterized as two α6 plus three β4 subunits ([Su et al., 2025](https://www.nature.com/articles/s41467-025-58333-0)). CHRNE is the epsilon component of adult muscle-type receptors containing α1, β1 and δ ([NCBI RefSeq summary](https://www.ncbi.nlm.nih.gov/gene/1145/)). These are distinct known assemblies, so the five-node MCODE clique should not be drawn as one demonstrated heteropentamer. Gene-level CPM is not a measurement of protein abundance or subunit stoichiometry.

CHRFAM7A is a partial duplication of CHRNA7. Coexpression suppressed α7-mediated currents in heterologous experiments ([Araud et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3162115/)). Here CHRFAM7A also increases modestly (+0.51; 0.42 CPM). We cannot infer its net functional effect in HT1080, or resolve shared versus gene-specific reads from this count matrix. Gene-specific exon/junction inspection and orthogonal assays are needed before treating the CHRNA7 signal as a validated α7 change. Low RIC3/TMEM35A RNA is further context for assembly uncertainty, not proof that assembly is impossible; the NACHO/TMEM35A requirement was demonstrated in a different experimental setting ([Gu et al., 2016](https://pubmed.ncbi.nlm.nih.gov/26875622/)).

## Reconciliation with GSEA and the chromatin signatures

The original broader cholinergic sets are not globally significant: acetylcholine-gated channel complex has NES +1.20, global BH q = 0.620; cholinergic synaptic transmission has NES +1.11, q = 0.693; response to nicotine has q = 0.804; response to acetylcholine has q = 0.908. Fixed-correlation CAMERA also fails global q < 0.05 for these sets. The module's highly significant annotation does not override these whole-ranking tests.

The specific Reactome annotation used for UP_01 has only 11 represented ranking genes and was below the original GSEA minimum of 15. It was therefore not tested in that GSEA analysis; it is not a failed enrichment test. The network annotation used a separate, documented size range and annotation-testing family. No thresholds have been relaxed in this follow-up.

CHRNB2 occurs in the leading edges of both enriched Mikkelsen H3K27me3 source signatures and both CNS axonogenesis GO sets. It is the only UP_01 gene in the leading edges of any combined global-BH < 0.05 set. This links one receptor-subunit RNA change to the previously observed developmental/chromatin-associated pattern, but does not establish H3K27me3 loss at CHRNB2 in HT1080 or chromatin regulation of the whole module.

CHAT is very low and below the expression filter, SLC5A7 has no assigned reads, and SLC18A3 is below the filter with opposite period effects. This analysis supplies no positive evidence for an induced endogenous acetylcholine-production/release program. It does not measure ligand availability in the culture medium or exclude alternative sources of acetylcholine.

## Practical follow-up

1. Validate CHRNA7 and CHRNB2 RNA first, with CHRNA9 as a decreasing family comparator. Include CHRFAM7A-specific measurements and gene-specific CHRNA7 coverage. Validate CHRNA6/CHRNB4 with assays that work at their observed low abundance.
2. Establish surface receptor expression and an acetylcholine-evoked response before assigning channel activity. Combine functional readouts with genetic perturbation of CHRNA7 or CHRNB2 and suitable receptor-selectivity controls; a single drug response would not identify the receptor subtype.
3. If the RNA pattern and receptor function are confirmed, test whether the relevant subtype contributes to the LOY-associated phenotype, then examine the CHRNB2 chromatin hypothesis with locus-specific chromatin measurements. Independent LOY derivations or addback are required for stronger causal attribution.

## Deliverables and reproducibility

[Publication PDF](outputs/figures/nicotinic_receptors/HT1080_nicotinic_receptor_followup.pdf) · [600 dpi PNG](outputs/figures/nicotinic_receptors/HT1080_nicotinic_receptor_followup.png) · [Editable SVG](outputs/figures/nicotinic_receptors/HT1080_nicotinic_receptor_followup.svg) · [Figure caption](outputs/captions/nicotinic_receptor_followup.md).

Source data: [gene summary](outputs/tables/nicotinic_receptors/gene_summary.csv), [all sample counts and expression](outputs/tables/nicotinic_receptors/sample_expression.csv), [matched differences](outputs/tables/nicotinic_receptors/matched_differences.csv), [original enrichment results](outputs/tables/nicotinic_receptors/original_cholinergic_enrichments.csv), [gene-set coverage](outputs/tables/nicotinic_receptors/cholinergic_set_coverage.csv), and [links to enriched signatures](outputs/tables/nicotinic_receptors/UP_module_links_to_enriched_sets.csv).

Run `extract_nicotinic_followup.R` and then `nicotinic_followup.py` from any working directory. The frozen 01_ingest/02_expression models, original enrichment results and STRING tables are reused. Counts are assigned reads; no count files, models, gene-level tests, normalization or original enrichment statistics are changed. This follow-up is descriptive and post hoc.
