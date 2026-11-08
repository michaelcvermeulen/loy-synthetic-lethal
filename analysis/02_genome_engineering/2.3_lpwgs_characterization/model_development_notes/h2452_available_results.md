# H2452 available results — 17 September 2026

## RT-qPCR excluded from the active thesis

Author decision, 17 September 2026: do not use the H2452 RT-qPCR findings. The interpretation was removed from Results, and both figure exports and their captions were retired to `archive/h2452_rtqpcr_retired_20260917_191609`. Raw measurements and processing records are preserved. This is an author-directed exclusion; no specific technical failure has been established.

## Early ddPCR

The master tables contain H2452 individual CP, CN and Q53 measurements and controls in December 2024, as well as parental/control observations in October 2024 and January 2025. These must not be confused with edited-population validation from the later five-guide retry.

The 3 December 2024 instrument export (`A:/Dropbox/LOY_paralogs/data/ddPCR/AMELX_03DEC2024_CRAIG/data_table.csv`) records AMELY/AMELX concentration ratios of 1.1303 and 1.0794 for control specimens A1/A5, 1.3390 and 1.1142 for CP specimens A3/A7, 1.0624 and 0.9679 for CN specimens A4/B1, and 1.1616 for Q53 specimen A6. The specimen assignments are recorded in `summary_table_ddPCR_edit.csv`. These individual-guide measurements show no substantial target depletion relative to the available controls. They are not independent biological replicate counts, and no inferential test was performed.

**Orientation issue to resolve before combining later runs:** in the 6 December export, the instrument ratio is attached to the AMELX row and is AMELX/AMELY according to the literal target labels. A375 has zero/near-zero signal in that row and strong signal in the row labelled AMELY, raising a target-label/channel assignment question. The 7 January parental H2452 export likewise reports RNaseP/DDX3Y on the RNaseP row. These fields cannot be pooled as direct Y/reference ratios without auditing assay orientation. Legacy master-table droplet ratios are also distinct from instrument concentration ratios. No new combined ddPCR percentage is reported here.

## Later retry and validation scope

The author previously confirmed a subsequent H2452 experiment using a five-guide system. Its exact composition and outcome have not yet been linked to an assay record. No edited H2452 five-guide group appears in the two reviewed master ddPCR tables. No H2452 lpWGS comparison was located in the WGS file inventory. Neither point establishes that the later experiment failed or that no other records exist.

The verified early ddPCR batch shows no substantial individual-guide AMELY depletion and remains separate from the later protocol development. A successful five-guide depletion result requires its own specimen-linked measurements.

## Author-directed pending Results draft

The author stated on 17 September 2026 that the five-guide retry data will be available soon and requested a draft comparing sgCTRL and LOY ratios for RPS4Y1, AMELY, KDM5D and DDX3Y. Results 3.4.6.4 now contains a manuscript-style comparison with explicit placeholders. These measurements are pending, not completed findings. Do not infer the target/reference denominator, number of replicates, guide identities, editing rounds or direction of effects from the requested group labels.

| Target | Reference gene | sgCTRL ratio | Edited ratio | Measurements per group | Change versus sgCTRL |
|---|---|---|---|---|---|
| RPS4Y1 | Pending | Pending | Pending | Pending | Pending |
| AMELY | Pending | Pending | Pending | Pending | Pending |
| KDM5D | Pending | Pending | Pending | Pending | Pending |
| DDX3Y | Pending | Pending | Pending | Pending | Pending |

Once verified, percent reduction can be calculated as 100 × (1 − edited/control) for each target using the appropriate group summaries. Retain increases if observed, define the sampling unit, and interpret changes as locus-dosage changes rather than a measured Y-negative cell fraction.


## Retrieved ddPCR tables

The current local master table contains 33 H2452 entries across October 2024, December 2024 and January 2025. [Extracted master rows](../../figures/chapter3/source_data/h2452_ddpcr_inventory/h2452_master_rows_unharmonised.csv) preserve the original ratios and labels; these are not a harmonised analysis table. [Verified 3 December concentration ratios](../../figures/chapter3/source_data/h2452_ddpcr_inventory/h2452_03Dec2024_verified_concentration_ratios.csv) were recalculated from AMELY and AMELX concentrations and agree with instrument Ratio values for all seven measurements. The specimen/guide assignments come from the edited master table. Concentration ratios remain unscaled because stock-specific reference copy number has not been confirmed here.

A fresh search of local ddPCR CSV records did not locate the later five-guide measurements. Their four-target Results placeholders remain pending. Repeat-run orientation questions described above remain unresolved; do not combine the historical master ratio fields with the verified concentration ratios.


## Five-guide data transfer status

Author clarification, 17 September 2026: the later five-guide H2452 ddPCR assays are complete, but the results are not yet on this computer. The author plans to add them tomorrow. Treat the four-target values as awaiting local transfer, not as experiments still to perform. Do not substitute the early individual-guide measurements for this completed five-guide experiment.


## Author-confirmed selection and delivery optimisation

Author account, 17 September 2026:

- H2452 was deliberately selected as a mesothelioma model because of the adverse association between tumour LOY and patient outcome, retained Y-linked expression, and resistance to selected X/Y paralog perturbations.
- A GFP plasmid was delivered using 16 distinct Lonza Nucleofector electroporation programs. Flow cytometry measured GFP-positive fractions, and the program with the highest fraction was selected for guide delivery.
- The initial guide experiment was unsuccessful; none of the assessed guides produced effective Y depletion.
- Y targeting was repeated with a pooled mixture of five guides. This is five-guide pooling, not five successive nucleofections or a fivefold dose. The follow-up ddPCR assays are complete and await local transfer.

The program ID, GFP-positive fractions, gating strategy, assay timing, viability measurements and replicate structure have not been supplied. Do not infer these details, or equate reporter positivity with Cas9 editing efficiency. The five-guide composition and exact retry delivery history remain to be linked to records.

Selection evidence: Qi et al. 2023, Fig. 4B-C/Table S4 (doi:10.1016/j.cell.2023.06.006) reports the mesothelioma progression-free-survival association. The source article was checked at https://ecor.mgh.harvard.edu/content/files/Rheinbay%20paper.pdf. Public baseline expression is recorded in `figures/chapter3/source_data/cell_line_profiles/selected_models_y_expression.csv` (nine detected Y-linked genes). This remains distinct from the excluded in-house H2452 RT-qPCR.

The saved 26Q1 dependency table (`figures/chapter3/source_data/model_dependency_comparisons/rna/dependency_comparison.csv`, ACH-000092) has Chronos effects DDX3X +0.792, EIF1AX -0.103, ZFX -0.186 and RPS4X -2.974. Accordingly, Results specify low sensitivity to the first three genes rather than resistance to all X/Y paralogs. These scores describe the current analysis; the author account establishes the historical model-selection rationale.
