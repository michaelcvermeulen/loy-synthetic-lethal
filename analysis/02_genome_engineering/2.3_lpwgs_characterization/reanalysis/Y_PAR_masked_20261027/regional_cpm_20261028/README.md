# Regional read totals and CPM — six realigned atlas samples

Values below are **primary read-end counts (CPM)** from the completed Y-PAR-masked BAMs. Both numerator and denominator require MAPQ >=37 and exclude unmapped, secondary, supplementary, QC-failed and marked-duplicate records (`-F 0xF04`). R1 and R2 are counted separately. The denominator is all eligible mapped primary read ends across all reference contigs, including chrM; it is not the number of pairs or QDNAseq's used-read total.

CPM = regional eligible read ends / total eligible mapped read ends x 1,000,000.

## Whole-region counts and CPM

| Sample | Male-specific Yp | PAR1 (X-mapped) | Euchromatic Yq | PAR2 (X-mapped) |
|---|---:|---:|---:|---:|
| A549 sgCTRL | 8,611 (776.05) | 9,122 (822.10) | 14,455 (1,302.72) | 1,328 (119.68) |
| A549 Y-depleted | 259 (23.19) | 6,199 (555.07) | 451 (40.38) | 852 (76.29) |
| HT1080 sgCTRL | 14,104 (1,073.93) | 9,387 (714.76) | 22,740 (1,731.51) | 1,528 (116.35) |
| HT1080 Y-depleted | 5,879 (474.78) | 6,087 (491.58) | 208 (16.80) | 775 (62.59) |
| DU145 sgCTRL | 12,339 (984.27) | 13,976 (1,114.85) | 14,807 (1,181.14) | 1,352 (107.85) |
| DU145 Y-depleted | 3,352 (301.01) | 9,330 (837.85) | 457 (41.04) | 944 (84.77) |

These are totals over the full stated intervals, without the figure's bin mappability or control-count exclusions. Male-specific Yp is chrY:2,458,321–10,565,750. Euchromatic Yq is chrY:10,883,086–27,449,937, matching the figure's Yq region. PAR1 and PAR2 are counted on chrX:1–2,394,410 and chrX:153,925,835–154,259,566, respectively. Coordinates here are one-based inclusive. PAR values are combined X/Y signal after Y-PAR masking, not Y-specific PAR counts.

The complete CSV also contains Yq including Yq12 (excluding PAR2), Yq12 alone, and DU145's distal/proximal Yp subdivisions. Reads overlapping a subdivision boundary may count in each subdivision; the combined region counts each eligible alignment once.

## CPM in Y-depleted cells as a percentage of sgCTRL

| Cell line | Male-specific Yp | PAR1 | Euchromatic Yq | PAR2 |
|---|---:|---:|---:|---:|
| A549 | 2.99% | 67.52% | 3.10% | 63.74% |
| HT1080 | 44.21% | 68.78% | 0.97% | 53.79% |
| DU145 | 30.58% | 75.15% | 3.47% | 78.60% |

## Comparison with the figures

| Cell / region | Whole-region CPM ratio | Plotted-bin CPM ratio | Median raw bin CPM ratio | Figure median |
|---|---:|---:|---:|---:|
| A549 / Male-specific Yp | 2.99% | 2.94% | 2.14% | 2.10% |
| A549 / PAR1 (X-mapped) | 67.52% | 67.52% | 67.43% | 66.83% |
| A549 / Euchromatic Yq | 3.10% | 3.15% | 2.68% | 2.69% |
| A549 / PAR2 (X-mapped) | 63.74% | 63.74% | 63.93% | 63.65% |
| HT1080 / Male-specific Yp | 44.21% | 44.06% | 43.89% | 42.65% |
| HT1080 / PAR1 (X-mapped) | 68.78% | 68.78% | 67.85% | 67.23% |
| HT1080 / Euchromatic Yq | 0.97% | 1.00% | 0.00% | 0.00% |
| HT1080 / PAR2 (X-mapped) | 53.79% | 53.79% | 53.60% | 53.20% |
| DU145 / Male-specific Yp | 30.58% | 30.66% | 26.23% | split in figure |
| DU145 / PAR1 (X-mapped) | 75.15% | 75.15% | 72.40% | 72.17% |
| DU145 / Euchromatic Yq | 3.47% | 3.53% | 2.14% | 2.20% |
| DU145 / PAR2 (X-mapped) | 78.60% | 78.60% | 72.49% | 71.95% |
| DU145 / Distal MSY Yp | 51.23% | 50.95% | 51.69% | 54.38% |
| DU145 / Proximal MSY Yp | 16.01% | 16.13% | 6.01% | 6.11% |

Whole-region totals, totals restricted to plotted bins and medians across bins answer different questions. The Y figure lines are medians of QDNAseq-corrected bin ratios; PAR lines are medians of autosome-normalised bin-depth ratios. CPM uses a library-wide eligible-read denominator, without GC/mappability correction. Therefore the percentages need not equal the plotted labels.

DU145 proximal Yp illustrates the distinction: its whole-region CPM ratio is 16.01%, its plotted-bin union ratio is 16.13%, and its median raw bin ratio is 6.01%, compared with the figure median of 6.11%. Of the 44 plotted proximal bins, 13 have raw retention above 20%; these account for 85.1% of the summed per-bin Y-depleted counts. Their higher counts raise the pooled total relative to the median. Per-bin counts are supplied as a diagnostic and can double-count reads crossing bin boundaries; regional and plotted-union totals above use the multi-region iterator to count each alignment once per union.

HT1080 has nonzero Yq read totals even though its corrected median bin value is zero (displayed as <1%); the plot does not claim complete absence of Yq reads.

## Denominators

| Sample ID | Eligible mapped primary read ends |
|---|---:|
| A549-CTRL2 | 11,096,006 |
| A549-LOY | 11,167,991 |
| HT1080-sgCTRL-D0 | 13,133,058 |
| HT1080-LOY-D0 | 12,382,503 |
| LNCaP-sgCTRL | 12,536,170 |
| LNCaP-LOY-2X | 11,135,664 |

## Files and validation

- `regional_reads_CPM_wide.csv`: compact six-sample table with both counts and CPM.
- `whole_region_reads_and_cpm.csv`: the four requested regions in long format.
- `regional_reads_and_cpm_all_definitions.csv`: whole regions, additional Yq/Yp definitions and exact plotted-bin unions.
- `CPM_ratios_vs_figure_medians.csv`, `figure_comparison_summary.csv` and `per_bin_raw_CPM_ratios.csv`: comparison diagnostics.
- `provenance.json` and `beds/`: filters, reference intervals, figure-input checksums and counted interval unions.

All six BAMs passed samtools quickcheck. X-PAR counts exactly reproduce the earlier coverage calculations; high-quality Y-mapped PAR counts are zero as expected from the reference mask. Plotted-union counts do not exceed their containing whole regions. BAM sizes and modification times are unchanged. These calculations do not modify the figures or draft 2.

Reproduction: run `count_regional_cpm.py` under the existing WSL runtime, then `render_cpm_summary.py` using Python. Legacy LNCaP sample IDs correspond to DU145.
