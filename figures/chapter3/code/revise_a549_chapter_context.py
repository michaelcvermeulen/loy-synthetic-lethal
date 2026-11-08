"""Apply the author-approved chapter-context revision, preserving adjacent sections."""
from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[3]
p=ROOT/'Chapter3_results.md'
text=p.read_text(encoding='utf-8')
start=text.index('#### 3.4.6.2 ');end=text.index('#### 3.4.6.3 ',start)
section='''#### 3.4.6.2 Transcriptional characterisation of A549

##### RNA-seq sampling design and global expression structure

Having established Y depletion in A549 cells, we next examined its transcriptional consequences. RNA-seq included four Y-depleted and four sgCTRL samples collected in a single batch, with three parental samples included for context. The primary comparison used the eight Y-depleted and sgCTRL samples and an unpaired analysis.

PCA showed substantial variation within both conditions, without clear separation between Y-depleted and sgCTRL samples ([global expression figure, panel a](figures/chapter3/ch03_main_a549_rna_global_structure.pdf)). PC1 and PC2 explained 50.4% and 13.5% of the variance. PC1 was strongly associated with the proportion of assigned reads mapping to mitochondrial genes (Spearman ρ = −0.93), indicating that variation in library composition accompanied the dominant expression pattern ([QC supplementary figure, panel d](figures/chapter3/ch03_supp_a549_rna_qc.pdf)).

Among 22,792 expressed autosomal and non-pseudoautosomal X-linked genes, none met FDR <0.05 in the unpaired voom–limma comparison ([global expression figure, panel b](figures/chapter3/ch03_main_a549_rna_global_structure.pdf)). Including mitochondrial read fraction as a covariate likewise yielded no significant genes. We therefore examined coordinated expression patterns alongside the individual-gene estimates.

[![Global transcriptional variation in A549 cells](figures/chapter3/ch03_main_a549_rna_global_structure.png)](figures/chapter3/ch03_main_a549_rna_global_structure.pdf)

**Figure (main; draft).** Global transcriptional variation in A549 cells. [Full caption](Chapter3_figure_captions.md#main---global-transcriptional-variation-in-a549-cells). [Supplementary figure](figures/chapter3/ch03_supp_a549_rna_qc.pdf).

##### Pseudoautosomal expression and retained X-linked paralogs

Having established depletion of Y-linked transcripts (Section 3.4.5.2), we next examined pseudoautosomal genes and the retained X-linked paralogs. Mean PAR1 expression was lower for PLCXD1 (41%), GTPBP6 (36%), CSF2RA (51%) and SLC25A6 (35%), with all four reductions retaining their direction when individual samples were omitted. PLCXD1, CSF2RA and GTPBP6 met correction across the 31 tested PAR features (adjusted P = 0.014, 0.014 and 0.018), whereas SLC25A6 did not (adjusted P = 0.062); none met correction across all expressed nuclear genes ([PAR and paralog figure, panel a](figures/chapter3/ch03_main_a549_rna_sex_chromosomes.pdf)). The PAR2 genes SPRY3 and VAMP7 showed smaller, less consistent differences ([PAR supplementary figure](figures/chapter3/ch03_supp_a549_rna_par_expression.pdf)).

The eight X-linked counterparts of the Y genes assessed in Section 3.4.5.2 showed small mean log2CPM differences, ranging from −0.29 for EIF1AX to +0.09 for RPS4X, with no coordinated transcriptional increase ([PAR and paralog figure, panel b](figures/chapter3/ch03_main_a549_rna_sex_chromosomes.pdf)). The broader paralog set used for cross-line comparison is considered in Section 3.4.6.4.

[![Pseudoautosomal and X-paralog expression in A549 cells](figures/chapter3/ch03_main_a549_rna_sex_chromosomes.png)](figures/chapter3/ch03_main_a549_rna_sex_chromosomes.pdf)

**Figure (main; draft).** Pseudoautosomal and X-paralog expression in A549 cells. [Full caption](Chapter3_figure_captions.md#main---y-linked-and-pseudoautosomal-expression-in-y-depleted-a549-cells). [Supplementary figure](figures/chapter3/ch03_supp_a549_rna_par_expression.pdf).

##### Broader transcriptional programmes

We ranked genes by their mean log2CPM difference to assess coordinated expression patterns. Across 23,457 MSigDB sets, 1,366 met a global adjusted P <0.05, predominantly with negative enrichment (1,300 sets). Eight Hallmark sets were negatively enriched: both MYC-target sets, E2F targets, the G2M checkpoint, mTORC1 signalling, oxidative phosphorylation, the p53 pathway and TNFA signalling via NFKB ([pathway figure, panel a](figures/chapter3/ch03_main_a549_rna_pathways.pdf)).

Reactome and Gene Ontology results extended these themes to DNA replication and repair, ribosome biogenesis, tRNA metabolism and mitochondrial gene expression. Positively enriched sets included fibrin-clot formation, protein activation cascades and GPCR ligand binding. Several coagulation-related sets shared FGA, FGB and FGG among their leading-edge genes, indicating overlapping contributions to these annotations (panel a).

Support for the enrichment patterns depended on the treatment of correlation between genes. CAMERA with a fixed inter-gene correlation of 0.01 supported six of the eight Hallmark sets at global adjusted P <0.05, whereas the p53 and TNFA/NFKB sets did not meet this threshold. When nonnegative correlations were estimated from the A549 data, no sets passed global correction (panel b). These results identify exploratory expression patterns, particularly lower cell-cycle-associated expression, rather than demonstrated changes in proliferation or pathway activity. The full Hallmark comparison is shown in the [supplementary figure](figures/chapter3/ch03_supp_a549_rna_hallmark_complete.pdf); comparisons with other engineered and public LOY models follow in Section 3.4.6.4.

[![Pathway enrichment in A549 is sensitive to inter-gene correlation](figures/chapter3/ch03_main_a549_rna_pathways.png)](figures/chapter3/ch03_main_a549_rna_pathways.pdf)

**Figure (main; draft).** Pathway enrichment in A549 is sensitive to inter-gene correlation. [Full caption](Chapter3_figure_captions.md#main---pathway-enrichment-in-a549-is-sensitive-to-inter-gene-correlation). [Supplementary figure](figures/chapter3/ch03_supp_a549_rna_hallmark_complete.pdf).

'''
p.write_text(text[:start]+section+text[end:],encoding='utf-8')

p=ROOT/'Chapter3_methods.md';text=p.read_text(encoding='utf-8')
marker='### 3.3.15 HT1080 UTY and KDM5D re-expression RNA analysis'
methods='''#### A549 baseline RNA-seq analysis

The A549 comparison comprised four Y-depleted and four sgCTRL libraries collected in one batch. Three parental libraries were included in a separate contextual analysis. Passage pairing and independent population derivation were not established in the retained sample records; the primary analysis was unpaired and did not treat the libraries as confirmed independent editing replicates.

Counts used the GENCODE v49 annotation. Expression filtering used edgeR::filterByExpr with minimum count 10 and minimum total count 15, followed by restriction to autosomal and non-PAR X-linked features. TMM normalisation was fitted to the 22,792 retained genes in the eight primary libraries. CPM and log2CPM (prior count 2) for other genes used the same effective library sizes. PCA was centred by gene without variance scaling. The parental-inclusive PCA used a separate normalisation across all 11 libraries and the same primary gene set.

Differential expression used limma-voom with a condition-only design and robust empirical Bayes moderation, with Y-depleted minus sgCTRL as the contrast. Benjamini–Hochberg correction spanned the primary test family. A sensitivity model included the fraction of assigned reads mapping to mitochondrial genes as a covariate. Descriptive effects were differences between group-mean log2CPM values, distinct from model coefficients. Leave-one-sample-out checks repeated normalisation and effect estimation after each library was omitted.

PAR features were tested in a secondary model fitted across expressed nuclear genes, retaining the primary effective library sizes. P values were corrected both across all expressed nuclear genes and separately across 31 expressed PAR-overlapping features. Percentage expression differences were calculated from arithmetic mean CPM. PAR-assigned counts represent total pseudoautosomal expression and do not distinguish X- and Y-derived transcripts. The eight-pair A549 display used the counterparts of RPS4Y1, ZFY, USP9Y, DDX3Y, UTY, NLGN4Y, KDM5D and EIF1AY; the cross-line analysis considered an expanded candidate set including PRKX and TMSB4X.

Gene-set analysis used mean log2CPM differences, one representative per gene symbol and frozen MSigDB 2026.1.Hs sets containing 15–500 represented genes. fgseaMultilevel tested 23,457 sets with seed 42, 10,000 simple permutations and a nominal P-value floor of 10⁻¹⁰. CAMERA used the same gene-symbol universe, normalisation and condition-only design, with either fixed inter-gene correlation 0.01 or estimated nonnegative correlations. Each analysis applied Benjamini–Hochberg correction across all tested collections. Representative Reactome and GO Biological Process sets were selected after inspection to illustrate themes alongside all eight significant Hallmark sets; they were not a prespecified test family. Complete tables, including gene membership and leading edges, were retained. These analyses assess ranked expression and model-dependent enrichment rather than directly measured pathway activity.

*Source: `RNAseq/analysis/A549/SEPT_18_2026_REANALYZE`, relative to the projects root. RNA preparation, sequencing details and unresolved collection/derivation links remain to be completed from laboratory records.*

'''
assert '#### A549 baseline RNA-seq analysis' not in text
text=text.replace(marker,methods+marker,1);p.write_text(text,encoding='utf-8')

p=ROOT/'Chapter3_figure_captions.md';text=p.read_text(encoding='utf-8')
def replace_caption(stem,caption):
 global text
 match=re.search(r'^## [^\n]+\n(?:(?!^## ).)*?\*\*Figure:\*\* `'+re.escape(stem)+r'`(?:(?!^## ).)*',text,re.M|re.S)
 if match is None: raise RuntimeError(stem)
 block=match.group(0)
 # Keep the heading stable so existing caption anchors continue to resolve.
 begin=block.index('\n\n**',block.index('[PDF]'))+2
 finish=block.index('\n\n[Source and reproduction notes]',begin)
 new=block[:begin]+caption+block[finish:]
 text=text[:match.start()]+new+text[match.end():]
replace_caption('ch03_main_a549_rna_sex_chromosomes','**Pseudoautosomal and X-paralog expression in A549 cells.** **a,** Selected PAR1 genes relative to mean sgCTRL CPM. Annotations give percentage changes and adjusted P values from unpaired voom–limma tests across 31 PAR features; none met genome-wide correction. Points show four samples per condition and horizontal marks show means. Blue circles denote sgCTRL and vermilion triangles Y-depleted samples; the dashed line marks 100% of mean sgCTRL. **b,** Mean log2CPM differences for eight X-linked counterparts. NLGN4X abundance was below 0.04 CPM in both groups.')
replace_caption('ch03_main_a549_rna_pathways','**Pathway enrichment in A549 is sensitive to inter-gene correlation.** **a,** The eight significant Hallmark sets (H), followed by selected Reactome (R) and GO Biological Process (GO) sets illustrating negative and positive enrichment. Rankings use mean log2CPM differences; blue indicates negative and orange positive enrichment. **b,** Benjamini–Hochberg-adjusted P values across 23,457 MSigDB sets for fgsea and CAMERA using fixed (0.01) or estimated nonnegative inter-gene correlations. Bold values meet adjusted P <0.05; shading represents −log10 adjusted P on a fixed 0–8 scale. The additional sets were selected to illustrate themes, not as an independent test family. Several nominal fgsea P values reached the analysis floor of 10⁻¹⁰.')
p.write_text(text,encoding='utf-8')

p=ROOT/'figures/chapter3/a549_rna_thesis_notes.md';text=p.read_text(encoding='utf-8')
text=text.replace('and the Y expression panel shows log2(CPM+1); these measures are not interchangeable.','whereas X-counterpart effects use differences in mean log2CPM; these measures are not interchangeable.')
text+='''
## Chapter-context revision, 24 September 2026

Removed the duplicate eight-gene Y-expression panel, already presented in Results 3.4.5.2. The sex-chromosome draft now contains PAR1 (a) and X counterparts (b). Earlier exports and captions are archived. The main pathway display now includes all eight significant Hallmarks plus nine illustrative Reactome/GO sets; exact labels, membership, leading edges and statistics are in `source_data/a549_rna_thesis/main_pathway_display.csv`. Selection is illustrative and post hoc, with no new tests or multiplicity family.

The full A549 fgsea table contains 1,366 significant sets (1,300 negative, 66 positive). The cross-line section reports 1,360 within its common-set subset; these counts have different denominators. Leading edges of fibrin-clot formation and protein activation cascade include FGA, FGB and FGG. The eight-pair A549 display and expanded cross-line paralog set are explicitly distinguished in the text and Methods.

Public A549 comparisons remain in 3.4.6.4. The outdated placeholder for a second public comparison in 3.4.6.2 was removed. A549-specific normalisation, statistical models, PAR testing and sample-design limitations are now documented in Chapter3_methods.md. No new biological analyses or LaTeX builds were performed.
'''
p.write_text(text,encoding='utf-8')
print('Updated A549 Results, Methods, captions and figure notes.')
