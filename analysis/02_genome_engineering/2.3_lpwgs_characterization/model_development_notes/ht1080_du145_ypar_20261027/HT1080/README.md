# HT1080: Y-PAR-masked realignment and atlas preview

Both authorised samples are realigned and validated. Preview exports use the line-specific Fig 3.10 template, with the broader 100-kb Y-bin selection and X-PAR tracks requested after the A549 comparison. Selected exports and draft 2 prose/captions are unchanged.

## Regional comparison

| Region | Bins | Original alignment, same bins | New alignment |
|---|---:|---:|---:|
| MSY Yp | 67 | 42.62% | 42.65% |
| euchromatic Yq | 101 | 0.00% | 0.00% |

These columns use the same bins and the same broader QDNAseq analysis on the original and masked alignments. They isolate the alignment update from the earlier bin-filter change. The DU145 original baseline is refitted from its two atlas samples, rather than silently compared with a different treatment cohort. 

## Validation and interpretation

- QDNAseq GC/mappability correction, normalisation, smoothing, segmentation and segmented normalisation were rerun on the matched pair. The broad reference mask uses non-N sequence >=60% and mappability >=30%, excluding PAR from the correction fit. Displayed MSY bins must lie wholly within the specified regions, have finite estimates in both samples, control corrected CN >0.001, and >=50 primary MAPQ >=37 control reads. Centromere and Yq12 are excluded. No minimum LOY count or effect-direction filter is applied, and no positive numerator floor is imposed. Zero/near-zero bins do not establish absolute DNA absence.
- PAR primary-read depth excludes secondary, supplementary, QC-fail and marked-duplicate records, uses MAPQ >=37, and is normalised by the median depth of the same 21,446 strict-mask autosomal bins in each sample. PAR bins start at the X-PAR boundaries; shorter terminal bins are retained with >=20 control reads and a density equivalent to >=50 reads/100 kb. These tracks do not reuse the old unmasked-reference PAR mappability annotation or apply QDNAseq PAR correction.
- Original BAM checksums are unchanged. Every primary read name, mate, sequence and quality is preserved through recovery and realignment. Recovered FASTQs contain the original BAM's processed reads, not raw sequencer FASTQs.
- Same verified reference/index and BWA-MEM2 2.2.1 algorithm, with fixed `-K 160000000`; five mapping workers per sample allow two concurrent jobs within the available CPU/memory budget (A549 used eight workers). Sorting/indexing uses samtools 1.23.1. No index rebuild was needed.
- Primary-read fingerprints, paired-FASTQ checks, BAM integrity and QDNAseq masks were verified. All 132 RNA observations and the ideogram rendering code are unchanged. The line's own source tables are used throughout.
- Bin eligibility depends on reference properties and control support, not LOY effect direction. Difficulty-region overlaps remain in annotations; the broader display does not imply every complex sequence is fully resolved.
- PAR tracks report X-mapped combined X/Y coverage, not Y-PAR retention. `combined_PAR_comparison.csv` gives whole-region old/new diagnostics; the plotted lines are medians of bin ratios, a different summary.
- `autosomal_count_comparison.csv`, `read_validation.csv`, the figure source-data validation, and `mask_checkpoints.csv` contain the checks. A reference-threshold or bin-size choice is not a cell-fraction estimator.

## Outputs

PDF/PNG/SVG, caption, code and source snapshots: `A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/drafts/ht1080_du145_ypar_20261027/ht1080`.

The four-sample alignment driver is `WGS/reanalysis/Y_PAR_masked_20261027/run_atlas_four.sh`. Fit the original/masked pair with `fit_pairs.R HT1080 original|masked`; run `atlas_pair_coverage.py HT1080 original|masked` under WSL; then `prepare_inputs.py HT1080`, `build_renderer.py`, `make_atlas_preview.py HT1080 --output-dir <preview folder>`, and `compare_and_finalize.py HT1080` from the projects directory.
