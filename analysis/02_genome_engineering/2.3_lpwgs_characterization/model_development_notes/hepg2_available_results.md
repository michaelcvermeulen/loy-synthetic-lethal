# HepG2 available results — 17 September 2026

The edited summary table records a pooled HepG2 CP + Q53 + P47 condition with one nucleofection, alongside CTRL1 and parental cells. This is a recorded assignment, not a new direct author confirmation. No second-nucleofection HepG2 group was found in the reviewed summary tables.

## Instrument concentration ratios

Ratios below were recomputed as target concentration/reference concentration from matched wells in the original data-table exports. Each assay has one well per condition; no between-replicate uncertainty or inferential test is available.

| Assay | Recorded run | Parental | sgCTRL | CP + Q53 + P47 | Reduction versus sgCTRL |
|---|---|---:|---:|---:|---:|
| AMELY/AMELX | 7 November 2025 | 0.9873268043 | 1.0051006256 | 0.1296609868 | 87.0997% |
| DDX3Y/RNaseP | 13 November 2025 | 0.3183035029 | 0.3186712693 | 0.0101347990 | 96.8197% |

Reduction = 100 × (1 − edited ratio/control ratio). These are reductions in target/reference dosage ratios, not measured percentages of Y-negative cells. Reference copy numbers in the experimental HepG2 stock are not independently established here; ratios were not converted to absolute copies per cell.

Sources:

- `A:/Dropbox/LOY_paralogs/data/ddPCR/AMELX Craig_NOV7_2025_2025-11-07-15-53/data_table.csv`: B01 sgCTRL, C01 MIX and D01 parental; both channels marked OK.
- `A:/Dropbox/LOY_paralogs/data/ddPCR/AMELX Craig_NOV_13_2025-11-12-22-50/data_table.csv`: A03 parental, B03 MIX and C03 CTRL; edited and control channels marked Manual, parental OK. No regating was performed. The folder timestamp refers to 12 November; the master table records the run as 13 November.
- `A:/Dropbox/LOY_paralogs/summary_table_ddPCR_edit.csv`: guide/round/population annotation.

The master summary's legacy `ratio` field uses positive-droplet-count ratios and is numerically different: those values were not substituted for concentration ratios. Existing `LOY_effectiveness_HEPG2_LNCAP` and `LOY_reduction_HEPG2_LNCAP` plots in the shared code directory derive from that older measure. Their historical labels and scale need updating before thesis reuse.

## Available evidence and gaps

The author subsequently confirmed on 17 September 2026: “we just did the ddPCR for HEPG2”. HepG2 characterisation therefore comprised ddPCR only; RT-qPCR and lpWGS are outside the experimental scope, not missing validation datasets to recover. These measurements support depletion at two Y-linked loci. The assays were run on different dates and should not be treated as replicate editing experiments. DNA collection dates, recovery time and the exact specimen relationship between the two runs remain unresolved. Complete or uniform chromosome loss is not established by these two measurements.

Results subsection 3.4.6.5 now reports the verified two-locus depletion and the ddPCR-only scope. No new HepG2 figure has been generated.
