# A549-CTRL2 Y-PAR-masked alignment pilot

One Y-retaining A549 control sample was recovered from its existing BAM and realigned with BWA-MEM2 2.2.1. The comparison below tests PAR read assignment; it does **not** estimate Y-specific retention or compare a depleted population with its control. No other sample or draft-2 figure was updated.

Recovered and validated **6,581,127 read pairs** (13,162,254 primary read ends). Read-name, mate, forward-sequence and quality fingerprints agree between the original BAM, paired FASTQs and new BAM. No singleton reads were lost. The input is the previously processed/trimmed sequence stored in the BAM, not newly recovered raw sequencer output.

## PAR mapping and coverage

| Region | Original high-MAPQ / all PAR read ends | Masked high-MAPQ / all PAR read ends | Original → masked PAR/autosomal depth |
|---|---:|---:|---:|
| PAR1 | 1,712 / 10,346 (16.5%) | 9,122 / 10,100 (90.3%) | 0.175 → 0.932 |
| PAR2 | 113 / 1,507 (7.5%) | 1,328 / 1,474 (90.1%) | 0.079 → 0.981 |

High MAPQ means ≥37, matching the original QDNAseq read-counting threshold. Counts combine X- and Y-PAR primary alignments; after masking, homologous reads are expected to align to X. The full tables also provide MAPQ 0 and 30 sensitivity results.

Masking substantially restores high-MAPQ PAR signal in this control, supporting ambiguous X/Y-PAR mapping as the main reason for its low filtered PAR coverage on the original reference. This pilot establishes a workable BAM-to-FASTQ realignment route. A matched Y-depleted sample is still needed to assess changes in combined PAR dosage; the masked alignment cannot identify the X or Y origin of homologous PAR reads.

![PAR masking pilot](A:/Dropbox/LOY_paralogs/projects/WGS/reanalysis/Y_PAR_masked_20261027/comparison/A549_CTRL2_PAR_masking_pilot.png)

The combined coverage diagnostic sums aligned bases across X and Y PAR intervals, divides by the X-PAR span, then normalizes to the median depth of 21,446 usable autosomal 100-kb bins. X and Y PAR spans differ in length; this is a descriptive combined-signal measure, not an absolute copy-number call. Overlapping paired reads can both contribute. No GC correction or masked-reference mappability model has been fitted for this pilot, and the old unmasked-reference PAR exclusion mask is not reused to assess newly mapped PAR reads.

## Tracking the original Y-PAR reads

| Original location | Read ends followed | Remapped to corresponding X-PAR | Previously MAPQ <37, now X-PAR MAPQ ≥37 | Unmapped after realignment |
|---|---:|---:|---:|---:|
| chrY_PAR1 | 5,416 | 5,117 | 3,850 | 63 |
| chrY_PAR2 | 727 | 677 | 632 | 0 |

All original PAR primary read ends were followed through the complete new alignment. Detailed destinations and before/after MAPQs are in `comparison/PAR_read_reassignment.csv`.

## Whole-library checks

| Metric | Original | Y-PAR masked |
|---|---:|---:|
| Primary read ends | 13,162,254 | 13,162,254 |
| Unmapped primary read ends | 423,099 | 423,144 |
| Properly paired primary read ends | 12,204,542 | 12,204,514 |
| Median autosomal depth, MAPQ ≥30 | 0.53737 | 0.53736 |
| Median autosomal depth, MAPQ ≥37 | 0.53699 | 0.53697 |

Across the common autosomal bins, old/new MAPQ ≥37 depth has Pearson r = 0.999999; 97.53% of bins have exactly unchanged aligned-base totals. The median absolute depth change is 0.000000, its 95th percentile is 0.000000, and its maximum is 0.027280.

## Reference and alignment provenance

- Original BAM: `../../t2t_bam/A549-CTRL2.bam`. It includes paired primary reads with sequences and qualities; no primary hard clipping, QC-failure flags or duplicate flags were found.
- The original BAM SHA256 is unchanged after the pilot; the check is recorded in `qc/source_bam_unchanged.json`.
- Local source reference: `H:/reference/chm13v2.0.fa.gz`. Only chrY PAR1 [0, 2,458,320) and PAR2 [62,122,809, 62,460,029) were replaced with Ns; 2,795,540 reference positions were masked.
- The prepared reference matches every contig of the official `chm13v2.0_maskedY.fa.gz` analysis set by sequence MD5. The downloaded official compressed-file MD5 is `bd90ddb80c86af7fcfeefe7a0909b175`. This uses the original mitochondrial sequence, not the separate rCRS replacement reference.
- All non-Y bases remain unchanged relative to the local source FASTA. Reference chromosome names and lengths match the original BAM; original BAM sequence dictionary tags do not independently provide a checksum of the historical HPC FASTA.
- BWA-MEM2 2.2.1 and samtools 1.23.1. Mapping uses eight workers and fixed `-K 160000000`, the original 16-worker default batch size. Name collation changes the input order; minor placement differences outside PAR can therefore occur. Autosomal stability is checked above.
- Reads are not trimmed again. All primary reads, including unmapped and MAPQ-0 reads, are recovered; secondary/supplementary records are omitted to avoid duplicating sequence records.
- The index is stored in the dedicated WSL scratch folder `/home/micha/codex-work/wgs-ypar-pilot-20261027/reference/`. An extra 16-GiB task-specific swap file was used for index construction and is removed by the index script afterward.

## Files

- `fastq/A549-CTRL2.R1.fastq.gz` and `fastq/A549-CTRL2.R2.fastq.gz`: recovered pairs.
- `alignment/A549-CTRL2.YPARmasked.bam` and `.bam.bai`: new alignment and index.
- `qc/`: original/new read QC, fingerprints, paired-FASTQ verification, headers and flagstats.
- `reference/`: official PAR BED, prepared reference, checksum evidence and masking provenance.
- `comparison/`: quantitative tables, read reassignment, coverage diagnostics and review plots.
- `runtime/environment_explicit.txt`: exact software package build list.
- `logs/`: extraction, reference preparation, indexing and alignment logs.

Reproduction: `setup_runtime.sh`, `prepare_reference.py`, `verify_official_reference.py`, `build_index.sh`, `recover_pilot.sh`, `align_pilot.sh`, then `compare_pilot.py` and `render_pilot_report.py`. The pilot extraction used `recover_from_collated.sh` to reuse a completed name-collated BAM; the recovered-read fingerprint validates its equivalence. Scripts use separate output files and refuse to overwrite completed FASTQs/BAMs.

Primary documentation: [T2T analysis-set definitions](https://github.com/marbl/CHM13#analysis-set), [samtools FASTQ recovery](https://www.htslib.org/doc/samtools-fasta.html), [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2).

## Matched A549 LOY extension (27 October 2026)

At the author's request, `recover_align_A549_LOY.sh` recovered 6,557,196 paired reads from the complete original A549-LOY BAM and realigned them with the same reference, index and settings as the control pilot. The original BAM checksum is unchanged, and the primary-read fingerprints match through BAM-to-FASTQ-to-BAM conversion. The sample-specific files in `fastq/`, `alignment/`, `qc/` and `logs/` retain this provenance. The original control-pilot results above remain unchanged.

The matched-pair PAR diagnostic is produced by `compare_a549_pair.py` in `a549_pair_comparison/`. The quality-filtered QDNAseq rerun and figure comparison are staged in `../../../SEPT_11_2026_THESIS/research/model_development/a549_ypar_realign_20261027/`. Preview exports for A549 Fig 3.12 are in `../../../SEPT_11_2026_THESIS/figures/chapter3/drafts/a549_ypar_masked_20261027/`; these do not replace the selected figure or draft 2 caption. PAR coverage is combined X/Y signal and is displayed separately from male-specific Y retention.

## Four-sample HT1080/DU145 extension (27 October 2026)

`run_atlas_four.sh` completed the four authorised atlas samples: HT1080-sgCTRL-D0, HT1080-LOY-D0, LNCaP-sgCTRL and LNCaP-LOY-2X. The legacy LNCaP filenames identify the DU145 samples used by the existing atlas. The six additional DU145 treatment samples and A549-BC14 were not realigned. All 26,677,931 primary read pairs were preserved through paired-FASTQ recovery and realignment; the four original BAM checksums are unchanged. New BAMs and indexes are in `alignment/`, recovered FASTQs in `fastq/`, and sample-specific validation records in `qc/`.

This extension reused the verified reference and index, with the same BWA-MEM2 version and fixed 160-million-base batch size. Five mapping workers per sample allowed two concurrent alignments within the available resources. Matched-pair QDNAseq analyses were rerun for both the original and new alignments, and coverage diagnostics are in `atlas_pair_coverage/`.

The [completed analysis and comparison](../../../SEPT_11_2026_THESIS/research/model_development/ht1080_du145_ypar_20261027/README.md) links the HT1080 Fig 3.10 and DU145 Fig 3.14 PDF/PNG/SVG previews, captions and source snapshots. They use the broader 100-kb display established after the A549 bin-filter comparison, retaining each line's own ideogram and RNA data. X-mapped PAR bins are shown directly in panel a as combined X/Y signal. Selected figures and draft 2 remain unchanged pending review.
