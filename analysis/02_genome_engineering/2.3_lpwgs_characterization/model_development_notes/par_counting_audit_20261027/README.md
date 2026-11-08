# PAR counting audit — 27 October 2026

## Finding

The DNA PAR percentages in the current regional atlases should not be interpreted as reliable Y-specific retention. A confirmed upstream bin-filter mismatch, combined with PAR mapping ambiguity, requires correction before quantitative DNA summaries are finalised. The ratio calculations and transfer of source values into the new figures are arithmetically correct; the saved eligibility flags they inherited are not the original quality mask.

No figure exports, captions, Results prose or input analysis objects were changed during this audit. The current source-notes README is flagged with this finding. Files here are investigative outputs only.

## Confirmed DNA evidence

1. `WGS/t2t_bins_100kb.rds` has 22,682 of 31,185 bins marked usable, including only 58 of 625 chromosome-Y bins. Every fully contained PAR1/PAR2 bin is marked `use=FALSE`.
2. The feature metadata inside **all three final QDNAseq objects** (`HT1080_QDNAseq_cn.rds`, `A549_QDNAseq_cn.rds`, `LNCaP_QDNAseq_cn.rds`) instead marks **all 31,185 bins usable**, including all 625 Y bins. Their mappability, non-N-base and blacklist annotations exactly match the original annotation; the rejected bins have been re-enabled. This affects 8,503 bins genome-wide, so it is not solely a PAR display issue.
3. PAR1 bins currently drawn in the atlases have mappability values of approximately **0.7–35.4%**. PAR2 bins have approximately **3.4–4.4%**. They are below the intended mappability threshold (55% in the retained annotation script).
4. The three source notebooks call `applyFilters(..., residual=FALSE, blacklist=FALSE, mappability=FALSE, bases=FALSE, chromosomes=NA)`, with comments assuming that this preserves the custom annotation's `use` flag. The [QDNAseq function implementation](https://rdrr.io/bioc/QDNAseq/src/R/applyFilters.R) constructs a new mask and assigns it to `binsToUse(object)`; the all-disabled call is consistent with the all-TRUE mask observed in the final objects. The [package documentation](https://bioconductor.org/packages/release/bioc/manuals/QDNAseq/man/QDNAseq.pdf) also specifies that FALSE disables mappability/base filtering. The exact historical installed function has not been extracted from the HPC environment, but the flag mismatch in the saved objects is directly verified.
5. BAM headers checked for HT1080-sgCTRL-D0, A549-CTRL2, A549-LOY and LNCaP-sgCTRL record bwa-mem2 alignment to `/global/home/hpc5173/data/chm13v2.0.fa`. The local download script obtains the unmasked-PAR analysis-set file of that name. [T2T documentation](https://github.com/marbl/CHM13/blob/master/README.md#analysis-set) distinguishes this reference from the version with Y PARs hard-masked. X–Y sequence similarity means assignment to a Y PAR coordinate alone does not establish Y origin. A mask applied during plotting cannot resolve read-origin ambiguity.

The original annotation predates the final QDNAseq objects (11 May versus 13–21 May 2026), and the annotation fields match, supporting comparison to the original mask rather than a subsequently changed set of mappability values.

The current figures additionally mask centromere/Yq12 and reject invalid control denominators. Those checks do not rescue the original lost bin-quality mask. The mappability annotation was available in the original RDS files but absent from the slim figure-source CSVs; the preceding atlas work should have cross-checked the original mask before describing the PAR bins as usable.

## Exclusion-only sensitivity check

Applying the original annotation mask to the **already computed** bin-level ratios produces the following descriptive medians. This is not a rerun of GC/mappability correction, normalisation or segmentation and is not a replacement analysis.

| Model and region | Current bins → original-mask bins | Current median → mask-only median (% of sgCTRL) |
|---|---:|---:|
| HT1080 male-specific Yp | 72 → 7 | 43.19 → 49.87 |
| HT1080 euchromatic Yq | 138 → 51 | 0.228 → 0.229 |
| A549 male-specific Yp | 73 → 7 | 2.01 → 3.71 |
| A549 euchromatic Yq | 148 → 51 | 1.34 → 2.58 |
| DU145 distal male-specific Yp | 23 → 2 | 51.08 → 50.57 |
| DU145 proximal male-specific Yp | 49 → 5 | 5.29 → 2.08 |
| DU145 euchromatic Yq | 125 → 51 | 1.77 → 1.96 |

All PAR summary groups become empty under the original mask. The broad non-PAR patterns persist in this limited check, but several estimates change and Yp support becomes sparse, particularly the two remaining DU145 distal-Yp bins. This does not establish that a full corrected analysis would leave all conclusions unchanged.

PAR2 is especially unstable in the current figures: HT1080 has only three eligible bins with ratios of about 82%, 0.1% and 130%; A549's three depleted PAR2 values are all at the numerical floor; DU145 has only one valid-control bin. These are not persuasive quantitative estimates of retained PAR2 DNA.

## RNA evidence and remaining uncertainty

The displayed PAR genes (`SLC25A6`, `CD99`, `VAMP7`) use chrX annotations in all three models. They are shared PAR genes, not Y-specific expression readouts. Retained expression is therefore not, by itself, evidence that Y-linked PAR DNA remains. WGS uses T2T coordinates; RNA counts use GRCh38 annotations, with the gene map illustrating the corresponding Y loci.

`check_rna_par_values.py` compares the 78 PAR observations from 26 RNA samples with their recorded featureCounts source tables, library totals and percentage calculations. Results are saved in `rna_PAR_checks.json` and `rna_PAR_verified_observations.csv`. This establishes extraction/arithmetic consistency only. It does not validate the historical RNA alignment or prove that all X/Y PAR reads were assigned without bias.

An older proposed PAR re-quantification workflow exists at `RNAseq/code/par_requant/`. Its README reports a historical near-zero PAR-count problem and proposes primary-alignment counting against X-only PAR annotations. That proposal alone does not prove which procedure produced the current `STAR_v2_PARfix` data. In particular, a folder name and retained featureCounts header do not establish that a Y-PAR-masked reference was used; the merge script preserves the original header. The current plots contain substantial PAR counts, so the old zero-count diagnosis must not be applied to them without confirming provenance.

## Recommended follow-up

1. Correct and verify preservation of the intended QDNAseq quality mask before GC/mappability correction, normalisation and segmentation. Assess the consequences genome-wide as well as across Y; do not simply drop PAR points and assume the upstream issue is resolved.
2. Until valid Y-specific PAR quantification is established, retain PAR locations on the ideogram but mark PAR DNA as not reliably quantifiable in the main Y-retention plots, without retention percentages. Keep PAR RNA explicitly labelled as shared X/Y expression.
3. If total PAR dosage is a scientific objective, use a separately validated PAR-aware approach that combines X/Y homologous signal or remaps to an appropriately masked reference. Such a measurement is total PAR dosage, not intrinsically Y-specific retention.
4. Confirm RNA reference masking and the actual re-quantification procedure from the original STAR/BAM/run records before interpreting precise PAR expression ratios mechanistically.

## Files and reproduction

- `extract_bin_masks.R`: reads the original annotation and final QDNAseq feature metadata without modifying them.
- `original_bin_annotation.csv`, `*_final_feature_metadata.csv`: exported metadata used for the comparison.
- `*_atlas_bins_with_original_mask.csv`: joins the adopted figure-bin audit with the original quality mask.
- `regional_mask_sensitivity.csv`: the regional comparison above, including mappability ranges.
- `check_rna_par_values.py`, `rna_PAR_checks.json`, `rna_PAR_verified_observations.csv`: raw-count, CPM and relative-expression checks.

Use `C:/Program Files/R/R-4.5.2/bin/Rscript.exe --vanilla extract_bin_masks.R PROJECTS_ROOT AUDIT_DIRECTORY` to regenerate the metadata exports. Use the bundled Python with the thesis `.tmp/qi_figure_python` path for the RNA check. No analysis rerun or manuscript integration is performed by these scripts.
