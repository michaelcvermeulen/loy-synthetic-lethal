# A549 Fig 3.12 Y-PAR-masked alignment preview

Status: preview only. Selected figure exports and draft 2 prose/captions were not replaced.

Both A549-CTRL2 and A549-LOY were recovered from their complete original BAMs and realigned using the same verified Y-PAR-masked reference, BWA-MEM2 2.2.1 settings and samtools 1.23.1 sorting. Every primary read name, mate, sequence and quality was preserved. These are recovered fastp-processed sequences, not the original sequencer FASTQs. The reusable native-WSL index was not rebuilt.

## Regional DNA comparison

| Region | Current figure | Corrected filter, old alignment | Corrected filter, masked alignment |
|---|---:|---:|---:|
| MSY Yp | 2.01% (73 bins) | 3.77% (7 bins) | 3.77% (7 bins) |
| euchromatic Yq | 1.34% (148 bins) | 2.62% (51 bins) | 2.73% (51 bins) |

The old-filtered and new-filtered columns use the identical 58 bins and the linear ratio of corrected copy-number estimates. The current figure used its historical set of bins and a 0.001 numerator floor. The new linear plot does not introduce that numerator floor; a control estimate >0.001 is still required. Subnormal negative roundoff values (approximately -1.65e-309) are displayed as zero, with raw corrected values retained in the source tables. Thus the first-to-last comparison includes filtering/normalisation changes as well as realignment; the middle-to-last comparison isolates the alignment update with the corrected pipeline held constant.

## PAR interpretation

After masking, the combined PAR1 signal is 67.41% of control and PAR2 is 63.26%. These values combine homologous X/Y DNA, now mapped to X; they are not Y-PAR retention percentages. The old figure's Y-only PAR1/PAR2 percentages therefore are not a like-for-like comparison with this inset. All masked Y-PAR intervals have zero direct mapped coverage by design and are excluded rather than plotted as biological loss.

## Analysis and validation

- QDNAseq 1.46.0 counts use the historical settings (paired reads, MAPQ >=37; defaults otherwise, including supplementary alignments); identical settings were retained for comparison. The separate PAR depth diagnostic uses primary mapped reads, excludes secondary/supplementary/QC-fail/marked-duplicate records and excludes deletions/skips from autosomal bed coverage.
- The original annotation remains valid for retained bins because only Y-PAR bases changed, and every bin overlapping those intervals already fails the mask. The original low mappability of duplicated PAR sequence is not repurposed as a quality score for the masked-reference PAR inset.
- GC/mappability correction, normalisation, outlier smoothing, square-root CBS segmentation and segmented normalisation were rerun for the matched pair with seed 20261027. No copy-number calling was required for this descriptive figure.
- Read fingerprints match for both samples; quality masks are preserved at each analysis checkpoint; 58 Y bins and all 88 RNA measurements are verified; selected exports are unchanged.
- See `regional_comparison.csv`, `raw_count_alignment_comparison.csv`, `mask_checkpoints.csv`, `masked_bins.csv`, and the preview's source-data `validation.json`. PAR source tables are in `A:/Dropbox/LOY_paralogs/projects/WGS/reanalysis/Y_PAR_masked_20261027/a549_pair_comparison`.

## Outputs

Preview PDF/PNG/SVG and a standalone draft caption: `A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/drafts/a549_ypar_masked_20261027`.

Reproduction order, from the projects directory: `recover_align_A549_LOY.sh` (WSL, one-time; refuses overwrite), `count_new_bams.R` for each new BAM (one-time), `fit_new_pair.R`, `compare_a549_pair.py` (WSL), `prepare_preview.py`, generated `make_a549_ypar_preview.py A549 --output-dir <preview folder>`, then `finish_preview_report.py`. The scripts reuse the completed control alignment and never overwrite original BAMs or selected figures.
