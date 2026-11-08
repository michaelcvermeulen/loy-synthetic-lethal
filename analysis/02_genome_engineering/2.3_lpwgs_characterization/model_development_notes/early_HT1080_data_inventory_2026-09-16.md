# Early HT1080 experiments — data inventory

## Author-confirmed chronology

On 16 September 2026, the author clarified that guide testing began in HT1080 because of its ease of use. Guides were first tested individually, and pooling followed. A subsequent clarification established that HT1080 settings were already published in the Lonza Knowledge Hub, allowing direct initiation of editing without a separate initial settings optimisation/QC phase. This supersedes the earlier wording that optimisation was performed in-house, and any implication that the work began with pooled guides or simultaneous testing in A549 and HT1080. The exact published protocol and programme used remain to be attached to the Methods record.

## Located data

The requested directory, `A:/Dropbox/LOY_paralogs/projects/ddPCR/data`, contains six curated files. These include early HT1080 ddPCR outcomes, but no explicit Lonza programme, pulse, kit or delivery-efficiency fields were located in the CSVs or the `LOY_clones.xlsx` workbook. The linked June 2024 raw-run text/CSV files also did not supply nucleofector settings. This finding does not establish that the settings are absent elsewhere in the project or bench records.

| ddPCR run date | HT1080 records | Source |
|---|---|---|
| 12 June 2024 | WT: A07, B07 and H06 (B07 has no usable ratio); CN: C07–E07; Q53: F07–H07. Q53 alone is explicitly labelled. | `projects/ddPCR/data/summary_all_INFO_TABLES.csv`, lines 3–24; original `data/ddPCR/JUNE12_24_AMELX/` |
| 3 December 2024 | CP LOY pool and Q53 LOY pool, plus C4, C7 and C24 clone labels. The master table assigns DNA date 26 November 2024 to the single-guide pools. | Both summary CSVs |
| 26 January 2025 | Historical C5 record, excluded from current analysis at the author's request. | Both summary CSVs and `26_JAN_2025.csv` |
| 24 and 27 February 2025 | CP, Q53 and CP+Q53 populations, Version 2, with controls/parental records; named combined-guide populations V2A/V2B in the all-INFO summary. | Both summary CSVs |

These are assay dates and table labels, not verified nucleofection dates or counts of independent editing experiments. Individual well entries must not be relabelled as independent biological replicates. The older rounds-analysis script interpreted `Version` as editing round; the author has corrected this: Version 2 is a repeated transfection experiment, and only MIX/MIX represents CP + Q53 delivered twice. C5 is excluded from current analysis. Population lineage still requires specimen records.

## CN label resolved on 17 September 2026

The author confirmed that the June 2024 label **HT1080 CN is the CN guide**. This supersedes the temporary assumption that it meant CP and resolves the identity question recorded below. Preserve CN in all derived tables and prose; do not classify it as CTRL using the summary notebook's parsing rule. See [the first-batch source record](ht1080_first_ddpcr_2024-06-12_notes.md).

## Historical label issue (resolved above)

The June 2024 raw exports label samples `HT1080 CN`. The summary notebook (`projects/ddPCR/analysis/ddPCR_summary_all.Rmd`, line 381) automatically assigns any sample matching `CN` to `CTRL`. This is a parsing rule, not independent evidence that CN was a non-targeting control. Given the separately documented centromeric guide named CN, its meaning in this experiment requires author confirmation before assigning a control role or normalising to it. No data or classification was changed.

## Raw source and measurement definition

`A:/Dropbox/LOY_paralogs/data/ddPCR/JUNE12_24_AMELX/` contains the June run's plate map, `data_table.csv`, `INFO_TABLE.csv`, droplet-amplitude exports and QuantaSoft files. The analysis route is `projects/ddPCR/analysis/ddPCR_12JUNE2024.R`.

The June summary's `ratio` is calculated from positive-droplet categories as `(AMELY + double-positive)/(AMELX + double-positive)`. Its separate `Conc_ratio` is calculated from the instrument concentration estimates, AMELY/AMELX. These should not be interchanged or interpreted directly as a percentage of Y-negative cells. No depletion effect size or inferential analysis was computed for this inventory.

## Next record needed

Locate the published Lonza HT1080 protocol and experimental record linking specimens to the programme used, kit, RNP amounts, cell inputs, recovery measurements and editing dates. A separate optimisation dataset is no longer expected for this initial stage, following the author's clarification. Confirm which original populations gave rise to the later single-guide and pooled-guide specimens; the June CN identity is now resolved as noted above.
