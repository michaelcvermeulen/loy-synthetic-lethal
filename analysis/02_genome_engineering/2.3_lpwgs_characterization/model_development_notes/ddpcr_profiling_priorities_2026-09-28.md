# Additional ddPCR profiling priorities for Chapter 3

28 September 2026. Review of the model-generation text and ddPCR figures, together with the downstream molecular-characterisation and add-back sections. These are gaps in the documented evidence; some measurements may already exist outside the chapter. No manuscript or analysis changes were made during this review.

## Main recommendation

Prioritise RPS4Y1 profiling in DU145 and HT1080, followed by DNA profiling of the actual populations used for RNA-seq, ATAC-seq and add-back. These measurements would resolve more consequential gaps than adding further measurements to already well-characterised guide conditions.

## Ranked sample priorities

| Priority | Samples to profile | Additional ddPCR | Gap addressed |
|---|---|---|---|
| **1** | **DU145 sgCTRL, 1x and 2x pooled-guide populations**, ideally including DNA matched to the RNA-seq collections | **RPS4Y1**, alongside AMELY and a Yq target | The chapter explicitly leaves this unfinished. AMELY is depleted but lies outside the retained distal-Yp interval. RPS4Y1 would directly test that retention and connect it to persistent RPS4Y1/ZFY expression. |
| **2** | **HT1080 sgCTRL and 2x CP + Q53 population**, plus available August/May molecular-sampling material | **RPS4Y1 and DDX3Y**; retain AMELY for comparison | The guide figure has no DDX3Y measurement for 2x and only one AMELY measurement. RPS4Y1 would independently validate residual Yp seen by lpWGS and RNA-seq, including whether it persisted across sampling periods. |
| **3** | **All six HT1080 add-back populations**, preferably at the RNA collection points | **RPS4Y1, AMELY and DDX3Y** | The untransduced Y-depleted and empty-vector populations differ in residual Y-linked expression as well as thousands of other transcripts. DNA measurements would establish whether the add-back populations also differ in retained Y material. |
| **4** | **DU145 pharmacological series:** sgCTRL, 1x DMSO, KU alone, both panobinostat doses and both combinations; include 2x without drug | **RPS4Y1 across the series**, and **DDX3Y in 2x** | Tests whether panobinostat improves depletion of the retained distal segment, rather than predominantly reducing already-depleted regions. Also fills the explicit "not measured" cell in the treatment figure. |
| **5** | **A549 sgCTRL and 2x population actually used downstream** | **RPS4Y1**, with AMELY and DDX3Y on the same DNA | A549 already has substantial DNA evidence. Its bigger gap is linking ddPCR, lpWGS and RNA-seq samples. A matched three-locus profile would strengthen that connection and the contrast with regional retention in HT1080/DU145. |
| **6** | **HepG2 sgCTRL and edited population** | **RPS4Y1**, optionally KDM5D, plus confirmation of existing targets | Current evidence is one measurement each for AMELY and DDX3Y. Distal-Yp profiling would test whether HepG2 also retains a segment that AMELY misses. Lower priority because this model was not carried forward. |

For H2452, recover the existing five-guide results first. The text says RPS4Y1, AMELY, KDM5D and DDX3Y were assessed, so the current gap appears to be transferring results rather than necessarily running more assays.

## Additional DNA checks for the HT1080 ATAC interpretation

Profile genomic DNA from the HT1080 ATAC sampling populations for chromosome 18 dosage. The text infers a possible chromosome-wide gain from ATAC/RNA, whereas earlier lpWGS does not show it. Several separated chromosome-18 loci would test that inference.

A secondary check could target the chromosome-19 region around 52-53 Mb, with flanking loci, where local chromatin closing and a later deletion remain alternative explanations. Confirm the reference assembly and interval against the plotted source data before selecting assays.

## Sample selection and interpretation

- Prioritise archived DNA or pellets closest to the actual sequencing collections, with their corresponding sgCTRL samples. Freshly expanded stocks establish present-day status but do not settle historical sample differences.
- For add-back samples, do not interpret UTY or KDM5D DNA assays as endogenous Y dosage unless their amplicons exclude the introduced constructs. RPS4Y1, AMELY and DDX3Y provide a cleaner starting panel.
- Repeat wells improve measurement precision; additional independently generated populations address reproducibility of the editing outcome.
- Near-zero measurements, especially in A549, would benefit from documented assay background, detection limits and uncertainty, consistent with [dMIQE guidance](https://doi.org/10.1093/clinchem/hvaa125).
- Even a complete four-locus profile measures regional dosage; it does not by itself determine the fraction of cells lacking the entire chromosome.

## Minimal first pass

If resources are limited, start with:

1. DU145 distal Yp: RPS4Y1 in sgCTRL and the established edited populations.
2. HT1080 2x and downstream samples: RPS4Y1 and the missing DDX3Y measurement, with matched controls.
3. The six HT1080 add-back populations: a common RPS4Y1, AMELY and DDX3Y panel.

## Sources reviewed

- [Chapter 3 Results](../../Chapter3_results.md): Sections 3.4.6, relevant DNA/RNA/ATAC comparisons in 3.4.7, and add-back characterisation in 3.4.8.
- [HT1080 guide-screening and ddPCR figure](../../figures/chapter3/selected/CH3_S3_4_6_1_F01_ht1080_guide_screening_ddpcr/CH3_S3_4_6_1_F01_ht1080_guide_screening_ddpcr.pdf).
- [DU145 ddPCR figure](../../figures/chapter3/selected/CH3_S3_4_6_3_F01_du145_ddpcr/CH3_S3_4_6_3_F01_du145_ddpcr.pdf).
- [DU145 treatment-comparison figure](../../figures/chapter3/selected/CH3_S3_4_6_5_F01_du145_all_treatments/CH3_S3_4_6_5_F01_du145_all_treatments.pdf).
- [A549 drafting notes](a549_results_drafting_notes.md).
- [DU145 treatment drafting notes](du145_drug_optimisation_notes.md), including the later author clarification that treatment conditions represented independent nucleofections.
- [Chapter 3 figure captions](../../Chapter3_figure_captions.md).
- Huggett et al. (2020), [The Digital MIQE Guidelines Update](https://doi.org/10.1093/clinchem/hvaa125).
