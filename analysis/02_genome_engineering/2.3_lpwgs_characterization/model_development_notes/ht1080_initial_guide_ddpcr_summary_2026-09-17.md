# HT1080 ddPCR: guide comparison, repeat transfections and sequential nucleofection

## Author-confirmed interpretation

Version 2 denotes a repeated transfection experiment. It must not be converted to a count of sequential nucleofections. CP, Q53 and CP + Q53 in that comparison are single-nucleofection populations. Only MIX/MIX represents CP + Q53 delivered twice. Sample C5 is excluded from current analysis, figure summaries and Results. Earlier interpretations are archived.

## Data used

The [main figure table](../../figures/chapter3/source_data/ht1080_ddpcr_nucleofection_rounds/main_measurements.csv) contains 23 measurements, all checked against raw instrument target/reference concentrations. The [sample mapping](../../figures/chapter3/source_data/ht1080_ddpcr_nucleofection_rounds/sample_round_mapping.csv) separates specimen identity and sequential nucleofection count. The main figure uses medians and interquartile ranges to match the requested reference plot type, with all individual measurements visible. No inference treats measurement counts as independent transfections.

| Condition | Sequential nucleofections | AMELY/AMELX range (n) | DDX3Y/RNaseP range (n) |
|---|---:|---|---|
| CP | 1 | 0.299-0.302 (2) | 0.090-0.101 (2) |
| Q53 | 1 | 0.315-0.329 (3) | 0.152 (1) |
| CP + Q53 | 1 | 0.166-0.207 (7) | 0.040-0.050 (4) |
| CP + Q53, MIX/MIX | 2 | 0.112 (1) | Unavailable |
| sgCTRL | 1 | 0.935-0.957 (2) | 0.481 (1) |

The February measurements use I6, I4, I7/I8 and H6 on 24 and 27 February 2025. The March measurements use J8, K1 and L2 for single-nucleofection CP + Q53 and K5 for MIX/MIX, all on 9 March 2025. The K5 ratio is 30.9185371398926 / 276.748962402344 = approximately 0.111721, from well C08. No DDX3Y measurement is present for K5 in that run. An older plotting script aliases November D0_LOY to MIX/MIX; that later specimen is not included without a verified lineage and assay-context link. Assay dates are not nucleofection dates.

## Repeatability across transfection experiments

The [repeatability table](../../figures/chapter3/source_data/ht1080_ddpcr_nucleofection_rounds/repeatability_measurements.csv) contains eight AMELY/AMELX measurements. Initial-experiment CP B5 had ratios 0.394-0.448 (two measurements), and Q53 B6 had 0.543 (one measurement), assayed on 3 December 2024. Repeat-experiment CP I6 had 0.299-0.302 (two measurements), and Q53 I4 had 0.315-0.329 (three measurements), assayed on 24/27 February 2025. Each experiment has one identified population per guide. These show the direction and variability of the response across separate experiments; no formal repeatability coefficient, control-normalised effect or significance test is calculated. Matched initial controls and exact recovery intervals remain unavailable. The older June CN/Q53 test is retained as its own supplementary batch and is not pooled into either group.

## Raw-data checks and exclusions

Use concentration ratios, not the historical positive-droplet-count ratios. Existing raw exports, gating and status flags remain unchanged. Manually analysed wells below 10,000 accepted droplets were retained; no new exclusion threshold was imposed. June WT B07 remains No Call. H6 is the HT1080 sgCTRL; the master table assigns I5 to A549 despite an erroneous HT1080 assignment in some all-INFO rows. Original source snapshots and inventories are retained as provenance, including records excluded from current analysis. C5 is listed in an explicit exclusion table and is absent from plotted and summary measurements. Ratios are not percentages of Y-negative cells, and locus-specific signals do not establish complete chromosome loss.

## Clone follow-up retained separately

The early all-INFO inventory additionally contains 17 clone assay records: 10 on 29 October 2024 (C2, C4, C7, C8 and C24), three on 3 December (C4, C7 and C24), two on 7 January 2025 (C24), and two on 27 January (C7). These are available in [the early inventory](../../figures/chapter3/source_data/ht1080_initial_guide_ddpcr/early_ht1080_inventory_including_clones.csv); they are not pooled with bulk-guide measurements or attributed to a guide without a derivation record.

Instrument checks show October C24 AMELY/AMELX concentration ratios of approximately 0-0.000362 and other clones approximately 0.917-1.049; December C24 remained zero while C4 and C7 were 0.945 and 0.978. January C24 RNaseP reference entries carry CHECK flags, so the all-INFO zero DDX3Y ratios require review rather than being interpreted as validated measurements of complete Y loss. January C7 has instrument concentration ratios of 0.502 for DDX3Y/RNaseP and 1.095 for AMELY/AMELX. The all-INFO inventory mislabels its C01 target as AMELY; the raw instrument export and master summary identify DDX3Y. These checks do not assign clone guide origin or establish whole-chromosome status.
