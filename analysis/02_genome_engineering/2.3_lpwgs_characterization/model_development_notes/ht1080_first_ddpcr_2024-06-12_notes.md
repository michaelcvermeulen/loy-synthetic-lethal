# First HT1080 ddPCR batch: 12 June 2024

Reviewed 17 September 2026. Source run: `A:/Dropbox/LOY_paralogs/data/ddPCR/JUNE12_24_AMELX/`.

## Identity and scope

The author confirmed on 17 September 2026 that **HT1080 CN means the CN guide**, superseding the earlier provisional CP assumption. Use CN throughout this batch. Raw labels remain in source snapshots and the well-level table. WT is the source comparator; it is not relabelled as a non-targeting-guide control. Neither CP nor CP + Q53 occurs in this batch. Assay date does not establish nucleofection date, number of transfections or linkage to the earlier RT-qPCR specimens.

## Measurements and checks

Results use instrument AMELY/AMELX **concentration ratios** (`Conc_ratio`), verified by dividing the two target concentrations in `data_table.csv` and by comparison with its AMELY `Ratio` field. The historical `ratio` field instead divides Y-positive by X-positive droplet counts, including double-positive droplets in both counts. These quantities are retained separately. No new droplet gating or model fitting was performed.

| Condition | Source wells | Evaluable wells | Mean concentration ratio | Range |
|---|---|---:|---:|---:|
| WT | A07, B07, H06 | 2 of 3 | 0.640785 | 0.629451-0.652120 |
| CN | C07, D07, E07 | 3 of 3 | 0.939060 | 0.890976-0.992253 |
| Q53 | F07, G07, H07 | 3 of 3 | 0.305769 | 0.302464-0.309252 |

B07 was an instrument No Call with 1,926 accepted droplets and zero positives; its missing ratio is not zero and was excluded from quantitative summaries. All eight evaluable HT1080 wells are marked manually analysed in the export and contain 14,710-19,034 accepted droplets. Well counts do not establish independent biological replication; no inferential statistics were performed. The rendered interpretation is lower AMELY relative to AMELX with Q53, with no corresponding reduction versus WT in the CN-treated samples in this run. The WT ratio is below one; its cause is unresolved. Ratios are not converted to percentages of Y-negative cells, and a single locus does not demonstrate whole-chromosome loss.

[Well-level data](../../figures/chapter3/source_data/ht1080_first_ddpcr/well_level_measurements.csv), [condition summary](../../figures/chapter3/source_data/ht1080_first_ddpcr/condition_summary.csv) and [source hashes](../../figures/chapter3/source_data/ht1080_first_ddpcr/source_manifest.json) preserve the measurements and source identity. Source snapshots are saved alongside these files.
