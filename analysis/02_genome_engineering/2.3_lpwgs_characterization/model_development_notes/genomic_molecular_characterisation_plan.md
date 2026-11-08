# Adopted cell-line-first molecular Results structure

Adopted by the author, 17 September 2026. The author requested completing HT1080 RNA-seq, ATAC-seq and integration before moving to A549, DU145 and a cross-line comparison. This document defines the adopted structure; figure numbering remains provisional. It does not report new findings or replace the source checks in the existing Results plan.

## Adopted title

### 3.4.7 Genomic and molecular characterisation of Y-depleted models

The current 3.4.6 structure places lpWGS alongside model development and ddPCR in each background. Retain that genomic foundation there, then briefly identify the populations carried into molecular profiling here. This lets the new section begin with HT1080 RNA-seq while avoiding repeated genome-validation narratives. Refer back to regional Y retention and relevant autosomal changes when interpreting expression or accessibility.

A single introductory sample-to-assay map should show the retained populations, controls, passages/sampling windows, RNA-seq coverage and HT1080-only ATAC-seq coverage. Distinguish in-house populations from public models and editing events from repeated sampling.

### 3.4.7.1 Transcriptional and chromatin characterisation of HT1080

**RNA-seq.** Establish the design and sample relationships, then present global expression structure, Y-linked/PAR expression, retained X paralogs, autosomal differential expression and supported pathways. Develop persistence or variation across sampling windows within the original pools, rather than implying independent editing replicates. Include informative uncertain or discordant effects. Use the final source-matched baseline branch.

**ATAC-seq.** Introduce the available sample design and essential quality evidence, then present differential accessibility, genomic distribution of affected regions and representative loci. Include motif or regulator analyses only where retained analyses support them. Distinguish accessibility from histone modification or direct factor binding.

**Integration.** Relate RNA and accessibility effects using explicit peak-to-gene links and eligible backgrounds. Show concordant, opposing and mixed-peak examples. Incorporate lpWGS copy-number context at the same loci and separate group-effect comparisons from truly specimen-matched measurements. Use relevant public annotations or datasets to test specific interpretations, identifying which evidence is new analysis versus existing annotation. Public data should answer a stated question, not become a catch-all subsection.

**Main displays:** sample structure/PCA; gene effects and Y/X-paralog expression; pathway summaries; accessibility effects and regional annotations; RNA-ATAC effect comparisons; paired RNA, accessibility and copy-number locus views. Several figures can develop this sequence.

**Supplement:** full library/sample QC, all windows and contrasts, passage/parental sensitivities, full RNA and peak results, complete enrichments, peak-to-gene alternatives, background/normalisation checks and additional concordant and discordant loci.

### 3.4.7.2 Transcriptional characterisation of A549 and comparison with public LOY models

Repeat the interpretable RNA sequence: design and genomic context, sample structure, Y/PAR/X-paralog responses, autosomal gene effects and pathways. Report effect sizes and uncertainty even where few genes pass significance thresholds. Include the saved paired/unpaired sensitivity.

Then introduce the public A549 spontaneous-LOY clone datasets. Establish the different derivation and sampling designs, compare continuous expression effects and pathway-level evidence, and show both agreement and discordance. Keep clone composition, reused clones across modalities, estimator sensitivity and copy-number context explicit. Do not presume public models replicate the engineered-pool response.

**Main displays:** A549 expression overview; selected dosage/pathway effects; public-versus-in-house design and effect comparisons.

**Supplement:** complete RNA diagnostics and contrasts; public modalities, clone omissions, gene-set methods, genomic exclusions and the full comparison tables.

### 3.4.7.3 Transcriptional characterisation of DU145

Use the same RNA framework, linking the profiled population to its documented editing history and genomic state. Present Y-linked/PAR and X-paralog expression, autosomal effects and pathways. Integrate relevant copy-number differences without assuming all expression changes are caused by chromosome loss. Do not imply ATAC-seq was performed in this background.

**Main displays:** expression overview; dosage-related effects; supported pathways and informative loci with copy-number context.

**Supplement:** complete QC, model contrasts, gene/pathway tables and copy-number-sensitive comparisons. Retain historical LNCaP labels only in provenance.

### 3.4.7.4 Shared and background-specific molecular responses

Ask which responses recur across HT1080, A549 and DU145 and which are restricted or uncertain. Begin with Y-linked depletion and retained X-paralog/PAR responses; then compare autosomal expression effects and pathway responses on a defined shared gene universe. Compare continuous effect sizes and directions before thresholded overlap. A gene significant in one background but not another does not by itself establish a background-specific effect.

Use the HT1080 chromatin data to inform interpretation of its transcriptional changes, while clearly identifying that chromatin accessibility was not measured across all three lines. Compare public evidence with the appropriate in-house background before extending a conclusion across models.

**Main displays:** X-paralog/PAR effect matrix or forest; pairwise RNA effect scatterplots; pathway-by-model matrix; selected shared and divergent gene programmes with uncertainty.

**Supplement:** UpSet plots with explicit thresholds and universes; complete pairwise effects; alternative gene/pathway backgrounds; effect-direction, expression-level and copy-number sensitivities. Avoid using DEG counts alone to rank response magnitude.

## Follow-on section

Section 3.4.8: Transcriptional responses to UTY and KDM5D re-expression. Keeping add-back separate lets the baseline characterisation and cross-line synthesis first establish the states whose reversibility is tested. Present construct and vector controls, reversal/persistence/overshoot, and direct construct comparisons. There is no matched add-back ATAC dataset in the current source plan. Molecular reversal is distinct from functional dependency rescue.

## How this maps to the existing plan

- Former 3.4.7 per-line RNA, 3.4.9 HT1080 ATAC and 3.4.10 integration become the line-specific accounts above.
- Former 3.4.11 public A549 comparison moves into the A549 account.
- Former 3.4.8 cross-model RNA becomes the closing synthesis.
- Former 3.4.12 add-back remains a separate follow-on section, now numbered 3.4.8.

Retain the source register and readiness checks for S03/S18-S28 in `Chapter3_results_plan.md`. No figure-count ceiling applies. Keep complete diagnostics and alternate analyses in supplementary series, with results essential to interpreting the main claim in the main account. Apply the thesis figure style and main/supp naming when producing figures.

## Immediate drafting unit

Start with the HT1080 RNA-seq design and global sample structure: identify the exact input/output branch, specimen-to-pool map, sampling windows, contrasts and available QC. Build the first RNA figure and its claim/evidence table before writing numerical findings. Complete specimen-linked lpWGS interpretation in 3.4.6 as the model-development foundation.
