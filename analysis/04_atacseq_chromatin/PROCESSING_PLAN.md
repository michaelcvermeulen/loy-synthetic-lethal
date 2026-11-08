# ATACseq reprocessing and integration plan

Date: 2026-09-08

Status: Approved by the user on 2026-09-08; see RUN_STATUS.md for execution milestones and results. The user specified that the six samples should be treated as independent replicates (3 CTRL and 3 LOY). RNA selection is confirmed and group-effect integration is complete; direct sample matching and lpWGS selection remain pending. The initial audit below records what was known during planning; the completed six-BAM audit is in qc/source_bam_audit.csv.

RNA selection update (user instruction): use the HT1080 combined-batch RNAseq result included in Ben's thesis, `~ batch + treatment`. The matching primary result is `A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080_BATCHES/DE/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv`, superseding the preliminary voom v113 candidate. The source notebook verifies DESeq2 and LOY/CTRL direction; its annotation is EnsDb v86. Six CTRL and six LOY samples contribute to the contrast, with three parental samples also used in original model estimation. The thesis manuscript itself was not located. Snapshots, exact source hashes, mapping accounting, and the completed integration are in config/RNAseq/ and integration/RNAseq/.

## Objective

Generate fresh peaks from the six sample-level ATAC BAMs, quantify LOY versus CTRL accessibility, and integrate the results with RNA expression and low-pass WGS copy number. Keep all new scripts, configuration, logs, tables, and figures within this directory. Reference source BAMs by path instead of duplicating them.

## What the initial audit established

- The previous sample sheet identifies HT1080 CTRL_11, CTRL_12, CTRL_21, LOY_11, LOY_12, and LOY_21; its documentation says the suffix denotes passage. Per the user's 2026-09-08 instruction, treat each sample as an independent replicate of its condition. Retain passage as metadata.
- All 20 ATAC BAM files could be opened for reading in the initial inventory; each has a BAI. This is not yet a full BAM integrity check.
- Candidate primary inputs are `../filtered/{sample}.bam`. The inspected CTRL_11 header records GRCh38/no-alt alignment with Bowtie2 2.5.4, proper-pair/MAPQ >=30 filtering, mitochondrial exclusion, and Picard duplicate removal. Check all six headers and surviving alignments before adopting this set.
- The previous plan calls `rep_bam/` a byte-identical duplicate of `filtered/`. The inspected CTRL_11 headers agree, but full-file identity has not been reverified in this audit.
- The inspected nf-core CTRL_REP1 header has different alignment settings and MAPQ >=1 filtering. Its REP-to-passage mapping is unconfirmed. Do not combine these alternate processing outputs with the six candidate inputs as additional replicates.
- The August 27 analysis reused old peaks. Its results and QC are historical comparisons, not fresh measurements for this run.
- A candidate RNA result exists at `A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080_BATCHES/voom_combined_v113/DE/ht1080_batches_LOY_v_CTRL_voom.csv`. Previous documentation describes two batches, 6 CTRL versus 6 LOY, and Ensembl v113 annotation. Confirm that this is the intended final result and obtain the full tested-gene universe.
- WGS contains `hg38_bam/HT1080-LOY-D0.bam` and `hg38_bam/HT1080-sgCTRL-D0.bam`, plus T2T counterparts. Existing HT1080 QDNAseq figures describe 100 kb T2T bins. Confirm assay matching and use verified hg38 copy-number results, or derive them from the hg38 BAMs when the integration stage is authorized.
- Windows R 4.5.1 and 4.5.2 directories exist. `samtools` and `macs3` were not found on the current PATH; WSL is present as an executable, but its distributions/tools have not been checked. This does not establish that the tools are absent elsewhere.

## Processing stages

### 1. Freeze inputs, references, and design

Create a manifest containing sample, condition, passage, clone/line identity, culture/preparation batch, sequencing batch, biological unit, BAM/index paths, and RNA/WGS links. Record file sizes and checksums during execution; validate BAM completeness, sorting, indexes, read groups, sequence dictionaries, pairing, and prior filtering/shift history.

Use hg38/GRCh38 throughout the main analysis. Pin the exact reference sequence dictionary, gene annotation/TSS list, chromosome sizes, blacklist, software versions, and configuration. Harmonize chromosome names and explicitly convert BED's zero-based half-open coordinates when needed. Verify assembly identity before any overlap.

Confirmed analysis decision (user instruction, 2026-09-08): treat CTRL_11, CTRL_12, and CTRL_21 as three independent CTRL replicates, and LOY_11, LOY_12, and LOY_21 as three independent LOY replicates. Use primary design `~ condition`, with CTRL as the reference and LOY versus CTRL as the contrast. Retain passage for QC and a secondary `~ factor(passage) + condition` sensitivity analysis; this adjustment does not change the primary independent-replicate design. Record clone/line and batch information when available without inferring independent clone derivation from the requested analysis assumption.

Deliverable: input audit, reference manifest, sample sheet, and explicit design decision.

### 2. BAM preparation and initial QC

Validate and consistently apply primary, properly paired, high-quality alignment filters, excluding secondary/supplementary alignments, failed-QC reads, duplicate fragments, and mitochondrial reads. Ensure both mates satisfy the chosen fragment eligibility rule; read-level MAPQ filtering can leave orphaned mates. Avoid repeating coordinate-based duplicate removal on already deduplicated input. Remove blacklist-overlapping fragments using a documented pair-aware rule and retain usable canonical nuclear chromosomes.

Keep fragment coordinates unshifted for the proposed paired-fragment peak calling and counting. Generate a separate Tn5 insertion representation for TSS profiles or later footprinting, applying the correction once after checking input provenance.

Measure usable fragments, pairing/orphans, MAPQ distributions, insert-length distributions and nucleosomal structure, TSS enrichment, chromosome coverage, and filtering losses. Report reads and fragments separately. Original mitochondrial contamination, alignment rate, and pre-deduplication complexity must come from upstream BAMs or valid original metrics; filtered BAMs cannot reconstruct discarded molecules. Mark unavailable metrics as unavailable.

Compare QC across samples and passages. Use the named implementation and matching reference for TSS metrics, with contextual quality flags rather than automatic sample exclusion from one threshold. ENCODE lists FRiP and reference-dependent TSS enrichment among its ATAC criteria: https://www.encodeproject.org/atac-seq/.

Deliverable: per-sample QC table, fragment/TSS plots, and documented sample inclusion decisions.

### 3. Fresh peak calling

Proposed primary caller: version-pinned MACS3 narrow peaks in paired-fragment mode (`BAMPE`), with initial q-value cutoff 0.01 and summit output. Use the verified effective genome size and `--keep-dup all` only after duplicate handling has been validated. BAMPE uses actual fragments and does not take a nonzero shift; insertion-site processing is a separate representation. See the official documentation: https://macs3-project.github.io/MACS/docs/callpeak.html.

Call every sample separately with identical parameters. CTRL ATAC is a biological comparison group, not MACS genomic input. Condition-pooled calls can assist reproducibility assessment and browser display, but pooled files never replace sample-level counts.

Assess within-condition peak overlap and rank reproducibility. Produce separate permissive ranked calls if needed for pairwise IDR, documenting ranking, threshold, and pair-combination rules. If pseudoreplicates are used, partition whole read pairs into disjoint sets with fixed seeds; they do not add biological replication.

Calculate FRiP against each sample's own peaks and the shared peak universe using explicit fragment denominators. Compare peak discovery at matched fragment depth as a diagnostic while retaining full-depth data for primary differential analysis. Peak counts or FRiP alone cannot establish global chromatin closing.

Deliverable: fresh sample peaks/summits, peak statistics, FRiP, and reproducibility assessment.

### 4. Shared peak universe, counting, and coverage

Proposed consensus rule: support in at least 2 of the 3 samples within either condition, followed by union across conditions. This retains condition-specific accessibility. Adjust the rule explicitly if QC changes sample inclusion.

Construct deterministic, nonredundant 501 bp summit-centered regions, record sample support, resolve overlapping candidates with an explicit ranking/tie rule, and recheck blacklist and chromosome boundaries after resizing. Evaluate peak ranks in a way that limits domination by the deepest library. Preserve original variable-width calls as reference output.

Count eligible paired fragments once per overlapping feature against exactly the same coordinates in every sample. Implementation decision: allow an integer assignment to each overlapping peak when a fragment spans multiple regions; report extra assignments explicitly and count each fragment only once across the peak union for FRiP. Preserve integer raw counts and stable peak IDs. Generate fresh per-sample coverage tracks with documented scaling; normalized tracks are for display, not differential-model inputs.

Keep chrY as a separate LOY-associated signal check; exclude it from main normalization and primary differential testing. Retain PAR peaks in the complete consensus but exclude them from the primary model, using the sourced GRCh38 BED in config/. Require high mapping quality and blacklist exclusion for other X/Y-homologous regions; these filters do not guarantee resolution of all homologous loci. ATAC chrY depletion is not an estimate of the fraction of Y-negative cells.

Deliverable: consensus BED, membership table, raw count matrix, coverage tracks, PCA/MDS, correlations, and separate chrY summaries.

### 5. Differential accessibility

Proposed primary model: edgeR quasi-likelihood on the common peak matrix, with design-aware abundance filtering, explicit normalization, and LOY minus CTRL as the positive effect direction. Use the confirmed primary design `~ condition` with three independent replicates per group; report passage adjustment as a secondary sensitivity analysis. Report all tested peaks, effect sizes, uncertainty, and BH FDR; use FDR <0.05 for the main significance table. Report an explicitly defined minimum-effect test or effect-size subset separately.

Predeclare peak-based TMM and a background-bin normalization sensitivity analysis. Examine GC effects, composition, copy-number instability, and behavior of putatively stable bins before interpreting either normalization. Where possible, use copy-number-stable autosomal bins for the background arm. Neither library-size scaling nor comparable TSS/FRiP proves an absolute genome-wide accessibility shift without an external calibration or defensible invariant reference. Do not choose the headline method by how many desired hits it produces. ATAC normalization sensitivity is documented in this primary benchmark: https://pmc.ncbi.nlm.nih.gov/articles/PMC9701614/.

Use csaw sliding windows as corroboration, with region-level FDR after window aggregation. Distinguish the full primary peak-based results from a subset supported by reproducible peaks and concordant window-based evidence. Method agreement is robustness evidence, not independent experimental validation. csaw documentation: https://bioconductor.org/packages/release/bioc/html/csaw.html.

Deliverable: full and significant differential tables, MA/volcano plots, model diagnostics, normalization sensitivity, and robustness comparisons.

### 6. Annotation and regulatory interpretation

Annotate promoters (proposed TSS +/-2 kb), gene bodies, and distal peaks using a pinned annotation compatible with RNA analysis. Preserve all supported peak-to-gene links and their evidence type; nearest-gene assignments remain hypotheses. Define a gene-level promoter summary for integration so genes with multiple peaks are not counted as independent observations.

Analyze pathways for gained and lost accessibility using the tested, eligible gene universe. Where appropriate use region-based enrichment or matching to address gene length and number-of-peaks bias. Test motif enrichment with accessible backgrounds matched for GC, width, and baseline accessibility. Interpret motifs as sequence enrichment; footprinting/TF activity is an optional later extension contingent on depth and bias correction.

Deliverable: annotated peaks, promoter gene summaries, enrichment tables, motif results, and locus figures.

### 7. RNAseq and lpWGS integration

Start the assembly and sample-matching audit early; complete integrated interpretation after the ATAC results are stable.

RNA: use full differential tables and tested-gene membership; harmonize stable gene identifiers, annotation, effect direction, and coordinates. Compare promoter accessibility and RNA changes, retain discordant results, and assess overlap against genes evaluable in both assays. If samples are not directly matched, perform group-effect integration rather than presenting sample-wise correlations.

lpWGS: attach per-sample or per-line segment values, uncertainty, assembly, and segment overlap to peaks and genes. Prefer matching hg38 segments. If only T2T segments are available, recompute from available hg38 WGS BAMs or perform validated liftover with explicit accounting for split/unmapped intervals. Renaming chromosomes cannot convert genome assemblies.

Report copy-number-stable and copy-number-altered regions separately. Relative QDNAseq depth is not automatically absolute copy number. D0 line-level WGS may flag plausible dosage confounding but cannot automatically correct every later-passage ATAC sample. Consider copy-number offsets only with suitable matching and a justified model; do not subtract log2 copy number mechanically, particularly in deleted or low-mappability regions.

Produce a joint table containing peak/gene IDs, ATAC effect and FDR, RNA effect and FDR, copy-number context, assay matching level, and linkage confidence. Prioritize concordant regulatory changes in copy-number-stable regions while preserving dosage-associated and uncertain findings.

## Execution and deliverables

Use numbered, restartable stages with a single configuration, pinned software, captured commands/session information, and explicit checkpoints. Execution uses local Ubuntu-24.04 WSL for BAM/peak processing and Windows R 4.5.2 for statistical work; environment details are in README.md. Temporary-heavy work uses Linux temporary storage and final deliverables are kept here.

Planned subdirectories: `config/`, `code/`, `logs/`, `qc/`, `bam/`, `peaks/`, `counts/`, `tracks/`, `differential/`, `annotation/`, `motifs/`, `integration/`, and `report/`. Create them as stages are implemented. The final report should include methods, sample accounting, QC limitations, statistical scope, and regenerable figures/tables.

## Decisions still needed

1. Direct RNA-to-ATAC sample relationships, if sample-wise integration is desired. Dataset selection is resolved; the current comparison uses group effects.
2. Intended lpWGS dataset, whether D0 line-level data are appropriate, and whether passage-matched WGS exists.
3. Any external calibration/spike-in information for later interpretation. The current analysis uses relative normalization and does not establish absolute genome-wide accessibility changes.
