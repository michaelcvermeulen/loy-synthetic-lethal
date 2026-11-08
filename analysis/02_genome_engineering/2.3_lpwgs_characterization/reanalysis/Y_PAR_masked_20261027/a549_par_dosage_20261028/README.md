# A549 PAR retention: checking the starting X:Y dosage

The observed approximately 63–67% combined PAR retention is compatible with the local control having approximately twice as much X-derived as Y-derived DNA. A 50% expectation assumes equal starting X and Y contributions to PAR and unchanged X dosage. Neither male origin nor the cell-line name establishes that assumption in these cancer cultures.

## Local checks

The same fixed, strict reference mask was applied to both completed A549 alignments, independently of their observed counts: 1,178 non-PAR X bins, seven male-specific Yp bins and 51 euchromatic Yq bins. PAR and region-boundary overlaps were excluded. Primary-read depth uses MAPQ >=37 and excludes unmapped, secondary, supplementary, QC-failed and marked-duplicate records. It is divided by the same median depth over 21,446 autosomal bins used for the PAR plots. The existing broad-filter QDNAseq fit was also summarised over these identical fixed bins; it was not refitted.

| Region | Control raw relative depth | Y-depleted raw relative depth | Control corrected relative DNA | Y-depleted corrected relative DNA |
|---|---:|---:|---:|---:|
| Non-PAR X | 0.6860 | 0.6845 | 0.6508 | 0.6581 |
| Male-specific Yp | 0.3866 | 0.0066 | 0.3496 | 0.0133 |
| Euchromatic Yq | 0.3596 | 0.0028 | 0.3402 | 0.0092 |

These are relative signals, with the autosomal baseline around one, not absolute chromosome-copy calls. X dosage is effectively unchanged. Control X:Y signal is approximately 1.8–1.9:1 across these independent raw/corrected summaries.

Assuming representative non-PAR X and Y dosage reflects their respective PAR contributions, and removal of the Y contribution with unchanged X, predicted remaining PAR is X/(X+Y). The raw-depth estimates predict 64.0–65.6%; corrected estimates predict 65.1–65.7%. An illustrative two-X-plus-one-Y dosage model predicts 2/3 = 66.7%. The exact absolute karyotype and individual-cell distribution are not identified by this bulk-depth check.

The observed PAR bin-median retention is 66.83% for PAR1 and 63.65% for PAR2. Whole-region autosome-normalised depth gives 67.41% and 63.26%, respectively. The result therefore does not require approximately 30% residual Y material to explain departure from a 50% expectation. Low residual male-specific Y signal is reported separately and does not prove absence of every Y-derived sequence or structural PAR retention.

## External context

[ATCC's A549 record](https://www.atcc.org/products/ccl-185) describes an aneuploid, hypotriploid stock with modal chromosome number 66. Its historical cytogenetic description reports two X and two Y chromosomes in most cells, with one or both Y chromosomes lost in 40% of the cells examined; the record notes that the cytogenetics concern initial seed stock. This supports checking the current samples rather than assuming a fixed catalogue karyotype. It does **not** establish a two-X/one-Y karyotype for these particular cultures; that is an illustrative interpretation of the local relative-depth evidence.

## Reproduction

Run `../a549_par_dosage_check.py` under the existing WSL runtime. `nonPAR_XY_dosage.csv`, `per_bin_dosage.csv`, `PAR_dosage.csv`, `dosage_model_predictions.csv`, and `provenance.json` preserve the numerical results and source definitions. No alignments, fitted models, figures or chapter text were changed.
