# Measured HT1080 chromatin-reference atlas

Completed 14 September 2026. Step 1 only: reference annotation and gene sets. No new GSEA, CAMERA, RNA differential analysis or ATAC integration was run.

## Deliverables

- [Publication figure (PDF)](outputs/figures/chromatin_atlas/HT1080_measured_chromatin_reference.pdf), [PNG](outputs/figures/chromatin_atlas/HT1080_measured_chromatin_reference.png), [SVG](outputs/figures/chromatin_atlas/HT1080_measured_chromatin_reference.svg).
- [Standalone illustrated report](outputs/report/chromatin_reference_atlas.html).
- [Compact promoter matrix: all 25,461 ranked genes](outputs/tables/chromatin_atlas/HT1080_promoter_marks_RNA_universe.csv).
- [Full ranked-gene feature matrix](outputs/tables/chromatin_atlas/HT1080_chromatin_features_RNA_universe.csv), including gene-body and signal summaries.
- [All 78,899 annotated genes](outputs/tables/chromatin_atlas/HT1080_chromatin_features_all_genes.csv.gz), including excluded/unknown features and Y genes as context.
- [Twenty gene-set definitions and sizes](outputs/tables/chromatin_atlas/gene_set_catalogue.csv), [symbol GMT](outputs/tables/chromatin_atlas/HT1080_GSE261834_measured_symbols.gmt), [gene-ID GMT](outputs/tables/chromatin_atlas/HT1080_GSE261834_measured_gene_ids.gmt), [membership](outputs/tables/chromatin_atlas/gene_set_membership_RNA_universe.csv).
- [Window, threshold and alternative-TSS sensitivity](outputs/tables/chromatin_atlas/gene_set_size_sensitivity.csv), [uncertain ranked promoters](outputs/tables/chromatin_atlas/RNA_genes_with_uncertain_promoter_annotation.csv).

## What the atlas contains

Eight standard NEED-seq histone tracks and two lamin peak tracks from GSE261834 were annotated against the frozen GENCODE v49 gene table. All 25,461 original representative gene IDs and symbols were preserved. There are 24,994 evaluable ranked promoters and 467 with more than 20% of the promoter covered by the hg38 blacklist. Unknown calls are retained, never converted to absence.

| mark       |   Any peak overlap |   >=10% coverage |   >=50% coverage |
|:-----------|-------------------:|-----------------:|-----------------:|
| H3K4me3    |              10032 |             7442 |                0 |
| H3K27ac    |               1673 |              743 |                0 |
| H3K4me1    |                 16 |                0 |                0 |
| H3K27me3   |                845 |              379 |                6 |
| H2AK119ub1 |               2502 |             1125 |                6 |
| H3K9me3    |               1247 |              530 |                3 |
| H4K20me3   |               1645 |              802 |                2 |
| H4K20me1   |                770 |              331 |                1 |
| LaminB1    |               3144 |             2221 |                0 |
| LaminB2    |               1769 |              957 |                0 |

Selected co-marked promoter sets:

|                       |   n_RNA_universe |
|:----------------------|-----------------:|
| H3K4me3 + H3K27me3    |              318 |
| H3K27me3 + H2AK119ub1 |              318 |
| H3K4me3 + H3K27ac     |             1524 |
| H3K9me3 + H4K20me3    |              390 |

The two 318-gene sets are separate definitions and must not be assumed to contain the same genes. Broad gene-body coverage sets and individual mark sets are also provided. Fourteen of the 20 sets exceed the original fgsea maximum of 500 genes; six fit the existing 15–500 limits. Sets were not truncated to force eligibility. Before step 2, prespecify size handling for the custom collection and assess the 467 unevaluable promoters as a background sensitivity.

## Frozen definitions and provenance

Inputs, peak thresholds and TSS rules were hashed before annotation. Only RNA gene IDs/names were read to define the universe; effect values did not determine genes, windows or thresholds. The input manifest includes the exact peak/bigWig files, blacklist, GTF, gene annotation and source metadata. Positive labels refer to deposited reference peaks under the chosen filter, not replicate-supported differential chromatin.

GENCODE contains 509,650 transcripts at 404,269 distinct gene-associated TSSs. One transcript per gene was selected by MANE_Select, then Ensembl_canonical, with deterministic fallback rules frozen in configuration (no fallback was needed). Transcript-level coordinates and all distinct-TSS peak calls are retained. The primary window is [TSS0−2000, TSS0+2001) in BED coordinates, after chromosome clipping. Minus-strand TSS0 is GTF transcript end minus one. The sensitivity windows use half-widths 1,000 and 2,500 bp. This symmetric convention includes the TSS base; the old `promoters(upstream=2500,downstream=2500)` definition contained 5,000 rather than 5,001 bp.

Narrow-peak minus-log10(q) must be >=5; broad-peak values must be >=2. The prespecified relaxed comparison uses >=3 and >=1.3, respectively. Filtering uses peak q-values within each source analysis; these are not RNA p-values. Blacklist bases are subtracted from both peaks and queried regions. Any nonzero overlap defines the main positive call; exact covered bases, covered fraction and stricter 10%/50% coverage summaries remain available. Gene-body overlap is reported separately from promoter overlap.

Co-marked sets require two marks in the same selected promoter window. The alternative-TSS analysis also requires both marks in the same TSS window, preventing false co-marking by combining marks at distant alternative promoters. This remains bulk co-marking, not proof of simultaneous marks in one cell or nucleosome.

Signal summaries integrate original bigWig intervals exactly over unmasked regions. Both represented-base mean and mean with unrepresented bases treated as zero are exported, together with represented fraction. Chromosome absence remains missing. All ten files have the expected GRCh38 chr1/chrX lengths and agree on shared sequence lengths. Eight histone tracks are described as RPKM in GEO; absolute values are not comparable between antibodies. **The lamin bigWigs are signed, approximately 20-kb signals, inconsistent with treating them as promoter-resolution RPKM.** Their normalization is unresolved at file level; they are retained as separately labelled coarse context. Lamin gene sets use the deposited peak files, not those coarse signal values, and are not claimed to reconstruct the paper's LAD states.

The publication describes replicate merging and GEO associates the processed histone files with the first sample accession. Exact per-file pooling and original assay-control relationships remain unresolved. This atlas therefore makes no independent-replicate support claim. No K562, protocol-comparison/MNase, imputed ENCODE or other-cell-line tracks enter these sets. Primary source: [GSE261834](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261834), [NEED-seq methods paper](https://academic.oup.com/nar/article/53/2/gkae1317/7958221).

## Interpretation and limitations

The H3K4me1 source file contains 1,905 deposited peaks, of which 293 pass the primary threshold; only 16 ranked promoters overlap a retained peak. The relaxed threshold identifies 92 such promoters. This is a limited reference feature, not evidence that other genes lack enhancer regulation. Distal overlapping H3K4me1/H3K27ac intervals outside all annotated ±2-kb TSS windows are exported only as candidate regulatory regions (minimum residual length 100 bp); no unsupported enhancer-to-gene assignments were made.

The deposited broadPeak files mostly contain short intervals (filtered medians approximately 244–335 bp). The broadPeak format alone does not establish broad heterochromatin domains. Coverage and continuous histone signal are retained to expose this limitation.

Candidate co-marked ranked genes number 260, 318 and 336 at ±1, ±2 and ±2.5 kb; the relaxed peak thresholds give 1,002 and the any-annotated-TSS definition gives 356. This demonstrates substantial threshold dependence. The old 432-gene list covered a different annotation/background; its 209 expressed genes are not directly comparable to 318 in the current ranked universe.

CHRNB2, CHRNA7, MMP1, MMP3 and ZMYND8 have no qualifying peak from the ten primary reference peak tracks at the selected promoter. This does not establish absence of chromatin regulation, and reference cultures can differ from ours. It does mean that Mikkelsen leading-edge membership is not independently corroborated as H3K27me3 promoter occupancy for CHRNB2 by this particular primary peak definition. Inspect quantitative signal, threshold/TSS sensitivity and independent ChIP before drawing a stronger conclusion.

Our ATAC samples are successive passages with uncertain RNA partners. No sample matching or population-level validation was inferred here. These are baseline reference features; measuring LOY-associated histone changes requires chromatin assays in the CTRL and LOY populations.

## Validation and rerun

Validation compares frozen hashes and identities, strand-aware TSS coordinates, independently reconstructed peak/blacklist overlaps and exact basewise bigWig signal summaries. See [validation](provenance/chromatin_atlas/validation.json) and [signal checks](provenance/chromatin_atlas/signal_validation.json). Plotting inputs and complete sensitivity membership are retained.

Run `python run_chromatin_atlas.py` from this analysis directory. It uses Windows Python/R for annotation and the existing WSL ATAC Python environment for bigWigs, then validates and regenerates the report. No package installation or original input modification is required. Original RNA/GSEA and source ATAC outputs remain unchanged.

## Figure caption

**Measured HT1080 chromatin reference from GSE261834.** (a) Numbers of genes with at least one unmasked base of a filtered reference peak in the selected promoter. Eight histone marks and two lamin peak tracks are shown. (b) Selected intersections require both marks in the same promoter window. (c) All 25,461 RNA-ranked genes are ordered by their eight-mark binary pattern and gene ID, without using RNA effect direction. Colour indicates the fraction with a detected peak within consecutive bins of up to 50 genes; grey denotes bins with no evaluable promoters. (d) Co-marked gene counts for promoter half-widths of 1, 2 and 2.5 kb; the dotted line marks the primary definition. Primary promoters include the TSS base plus 2,000 bp on either side (4,001 bp before clipping/masking), selected using MANE_Select then Ensembl_canonical tags. Peak filters are minus-log10(q) >=5 for the deposited narrow peaks and >=2 for broad peaks; ENCODE hg38 blacklist bases are subtracted. Promoters with more than 20% masked bases are unevaluable (467 ranked genes); 24,994 are evaluable. Counts are descriptive and gene sets overlap. No hypothesis tests, error bars or independent replicate counts are implied. Source tracks may pool replicates. Co-marking does not establish bivalency within the same cell or nucleosome. The data describe reference HT1080 cultures, not chromatin changes following LOY.

