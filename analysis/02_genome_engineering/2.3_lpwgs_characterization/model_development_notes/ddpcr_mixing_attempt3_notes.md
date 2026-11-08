# ddPCR cell-mixture check — attempt3

## Author-confirmed procedure

On 16 September 2026, the author explained that this experiment preceded HT1080 guide testing. KP4 and H2452 cells were mixed in known proportions of 0, 0.25, 0.5, 0.75 and 1 **before DNA extraction**. DNA was then extracted and analysed in triplicate with the AMELX/Y TaqMan ddPCR assay. This supersedes older general notes describing a mixture of extracted DNA.

Original data: `A:/Dropbox/LOY_paralogs/data/ddPCR/control_test/attempt3/ratios.csv`. An unchanged copy is stored in `../../figures/chapter3/source_data/ddpcr_mixing_attempt3/ratios.csv`.

## Recorded measurements

The file contains two series labelled `50uL` and `100uL` in `Conc`. The author subsequently specified that these denote DNA concentrations of 50 and 100 ng/µL in the ddPCR reaction. The `Expected` fraction denominator remains to be explicitly confirmed. Rep1–Rep3 are technical-replicate slots, as labelled in the associated analysis script. There are 24 numeric measurements and six NA entries across 30 slots; missing entries are not zeros. The reason for each missing entry has not been established.

**Author clarification, 16 September 2026:** the earlier provisional recollection was superseded by the instruction to treat `50uL` as 50 ng/µL and `100uL` as 100 ng/µL for the ddPCR reaction. This is the author-supplied concentration mapping used in the figure, caption and Methods; it was not inferred from the legacy CSV units. The meaning of `Expected` was not explicitly confirmed in that reply.

| Expected | 50 ng/µL mean | Available wells | 100 ng/µL mean | Available wells |
|---|---:|---:|---:|---:|
| 0 | 0 | 3 | 0 | 2 |
| 0.25 | 0.2063 | 3 | 0.1997 | 3 |
| 0.5 | 0.3580 | 1 | 0.4190 | 3 |
| 0.75 | 0.6880 | 2 | 0.7867 | 3 |
| 1 | 0.9905 | 2 | 0.9615 | 2 |

Descriptive calculations from all nonmissing replicate values, computed separately by the recorded `Conc` label:

| Series | n | Pearson r | OLS slope | OLS intercept |
|---|---:|---:|---:|---:|
| 50 ng/µL | 11 | 0.9929547694 | 0.9761052632 | −0.0275000000 |
| 100 ng/µL | 13 | 0.9882356056 | 1.0194545455 | −0.0374965035 |

Calculations use `Expected` as x and each available Rep1/Rep2/Rep3 entry as y. Pearson r is `sum((x−mean(x))*(y−mean(y))) / sqrt(sum((x−mean(x))^2)*sum((y−mean(y))^2))`; the ordinary least-squares slope is the cross-product sum divided by the x sum of squares, and the intercept is `mean(y) − slope*mean(x)`. No hypothesis test was performed. Technical wells do not establish independent mixture/extraction replication. High correlation describes tracking across the series and does not establish exact agreement or an assay detection limit. The intermediate mixtures show deviations from the identity line.

## Source checks and open points

- `projects/ddPCR/analysis/ddPCR_attempt3_control.R` reads this ratios file for its regression/correlation plot and labels replicates as technical.
- The script also processes droplet files and calculates `(AMELY + double-positive)/AMELX-only`; it then loads the supplied ratios file without writing that calculated table to it. Therefore this script does not establish the exact producer or ratio definition of the supplied values. Do not silently treat the droplet formula as a verified AMELY/AMELX concentration ratio.
- The linked plate file is named `AMELX 07May2024.csv`; it has AMELX/VIC and AMELY/FAM targets and coded sample names (A1, A2, etc.). The filename suggests a May 2024 run but does not by itself identify mixture/extraction dates.
- Five wells are excluded in the script (B07, A07, B08, E09, H07). A well-to-mixture/replicate mapping is needed before connecting exclusions to the six NA entries.
- Confirm the mixture-fraction direction, exact assay identifier, DNA extraction method, reaction volume, ratio calculation and missing-well reasons before finalising quantitative Methods and a calibration caption.
- This experiment tests known mixtures of these cell lines. It does not independently establish a universal conversion of bulk Y dosage to Y-negative cell fractions in all engineered backgrounds.

## Placement

Use this as the opening assay check in 3.4.6, followed by HT1080 individual-guide testing using published Lonza settings, pooling and expansion to other backgrounds. The author clarified that no separate initial settings optimisation/QC phase was needed. Put the pre-extraction mixing procedure in Methods 3.3.11. A future figure should show individual wells and the identity line, keeping the two DNA concentrations in separate panels. Captions belong in `../../Chapter3_figure_captions.md`.

The Results opening and Methods 3.3.11 have now been drafted in that order. Numerical correlations are shown descriptively in the figure and retained in these notes; the exact ratio producer remains to be verified; the Results states only the observed increase across the known-mixture series.


## Agreement assessment, 29 September 2026

The author clarified on 23 September that inputs were 50 and 100 ng total DNA per reaction, Expected is the H2452 Y-positive cell fraction, and Observed is the measured Y/X ratio. These clarifications supersede the earlier concentration labels and denominator uncertainty above.

Using all 24 nonmissing technical measurements, mean absolute error relative to the nominal fraction was 0.042375 (4.2375 percentage points); mean signed error was -0.032375. At 50 ng (11 measurements), MAE was 0.0396364; at 100 ng (13 measurements), MAE was 0.0446923. Maximum absolute error was 0.142, in the sole available 50-ng measurement of the 50% mixture. These statistics weight each available well equally, retain zero-mixture wells, exclude missing values and apply no calibration. They describe the complete cell-mixing, extraction and measurement workflow; technical wells are not independent mixture preparations. Reproduction: figures/chapter3/source_data/ddpcr_mixing_attempt3/calculate_agreement.py, with per-measurement and summary CSVs alongside. Results now reports agreement and supports quantitative use under the tested mixture conditions.
