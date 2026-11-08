# Supp Fig 3.11: acrocentric short-arm support audit

The five masked arms have **zero bins passing the existing full quality mask**. The original heatmap had nevertheless included every nominal bin on these arms. This supports retaining the display mask, but the reason is unreliable arm-level measurement, not absence of sequencing reads.

| Arm | Bins included originally | Mappability >=55% | Full quality mask | Median annotated mappability |
|---|---:|---:|---:|---:|
| 13p | 165 | 0 | 0 | 1.13% |
| 14p | 114 | 2 | 0 | 8.59% |
| 15p | 172 | 4 | 0 | 0.84% |
| 21p | 111 | 1 | 0 | 4.22% |
| 22p | 142 | 7 | 0 | 7.52% |
| 5p | 483 | 456 | 334 | 91.66% |
| 8p | 453 | 418 | 379 | 91.21% |
| 8q | 1010 | 990 | 828 | 91.34% |
| 20p | 280 | 262 | 243 | 91.41% |

The saved mask requires non-N sequence >=60%, annotated mappability >=55%, and difficulty-region overlap <20%. Every bin on the five acrocentric short arms fails the difficulty-overlap criterion; only 0-7 bins per arm meet the mappability criterion alone. These thresholds predate the heatmap masking and were not chosen using the observed arm differences.

## Read counts are present

| Arm | Control bins with counted reads | Control total cached read counts | Median counts per 100-kb bin | Median / usable-autosomal median |
|---|---:|---:|---:|---:|
| 13p | 68/165 | 3,722 | 0 | 0.000 |
| 14p | 91/114 | 13,783 | 46.5 | 0.103 |
| 15p | 66/172 | 11,197 | 0 | 0.000 |
| 21p | 76/111 | 6,569 | 8 | 0.018 |
| 22p | 104/142 | 10,260 | 21 | 0.046 |
| 5p | 472/483 | 307,879 | 645 | 1.427 |
| 8p | 442/453 | 155,327 | 344 | 0.761 |
| 8q | 998/1010 | 601,636 | 601 | 1.330 |
| 20p | 267/280 | 87,371 | 321 | 0.710 |

These are cached QDNAseq alignment counts before GC/mappability correction and bin-quality exclusion, using the original read-counting settings. They are not sequencing depth, absolute copy number, or independent observations of an arm-level change. Ratios here only compare median bin counts with the same sample's usable autosomal-bin median. All eight libraries are tabulated in `arm_read_support_all_samples.csv`.

## Why the old heatmap still contained strong colours

The historical `wgs_bins.csv` has `use=TRUE` for all bins, including these arms. Its arm summarizer therefore included 165, 114, 172, 111 and 142 bins for 13p, 14p, 15p, 21p and 22p respectively. This is consistent with the previously diagnosed QDNAseq filter-reset problem. The corrected staged DU145 output marks all these bins excluded and their CN/segment values missing in every library.

The existing heatmap uses medians of saved segment differences. Having many nominal bins (often sharing segment values) does not establish many independent reliable measurements. This audit does not prove that every signal on these arms is artefactual, or that biological changes are absent; it shows that these data do not support interpreting them as well-resolved arm-level CN changes.

The grey RNA cells use a separate availability rule (at least 20 genes per arm) and should not be used as evidence of WGS coverage or mappability. No inference of zero genes is made from missing RNA summaries.

## Scope and reproduction

- No figure, caption, manuscript, source arm table or historical analysis was changed in this audit.
- Reproduced original arm bin counts for every arm and all seven heatmap DNA rows; checked the corrected exclusion in all eight DU145 libraries.
- Arm assignment follows the existing CHM13v2.0 cytoband/midpoint rule used by the figure source.
- Run `export_counts.R` with the existing R 4.5.2/QDNAseq environment, then `audit_arms.py` with the existing pandas environment.
- See `arm_bin_support.csv`, `arm_read_support_all_samples.csv`, `target_and_comparator_bin_details.csv`, and `input_manifest.json` for the numerical evidence and input checksums.
- The staged corrected WGS results remain separate from adopted draft-2 numerical figures; this check does not adopt that reanalysis.
