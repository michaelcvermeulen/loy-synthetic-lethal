# Guides undetected in both screens

Audit date: 20 October 2026. The analysis uses six HT1080 libraries and seven retained A549 libraries; all A549 P3 libraries remain excluded.

Strict non-detection means zero raw counts in every retained D0 and endpoint library. Both original matrices contain the same 77,441 unique guide IDs and matching gene annotations.

| Scope | Guides with zero reads throughout |
|---|---:|
| HT1080 | 236 |
| A549 | 307 |
| Both screens | 223 (0.288% of the library) |
| HT1080 only | 13 |
| A549 only | 84 |

DDX3X_4640_TCTCGGTTCCTTAAATGAGG is one of the shared undetected guides. There are 222 other such guides: the complete shared set comprises 222 gene-targeting guides across 221 genes and 1 non-targeting guide. This is a guide-level limitation and does not mean that all guides targeting the affected genes were absent.

[Full shared-zero list and raw counts](guides_undetected_in_both_screens.csv) · [Either-screen list](guides_undetected_in_either_screen.csv) · [Affected-gene coverage](affected_gene_coverage.csv) · [Selected genes](selected_gene_coverage.csv) · [Source hashes and definitions](summary.json).

## Separate endpoint-only losses

75 HT1080 guides and 123 A549 guides had at least one D0 read but zero reads across all retained endpoints within that screen. These are separate from the guides never detected. Neither category alone establishes a biological cause. [Counts](detected_D0_zero_endpoint_guides.csv).

## Exploratory sequence observation

Of the 223 shared zero-read guides, 196 begin with TCTC or AGACG; 480 of all 77,441 reference guides have either prefix. DDX3X guide 4640 begins with TCTC. This is a post hoc descriptive check, not a mechanistic test. Shared-zero sequences are unique within the reference, so duplicate sequence rows do not explain this set. These count tables cannot distinguish physical library underrepresentation from sequence/reference or read-processing problems. No synthesis, cloning, mapping or trimming cause is assigned, and no count matrix, figure or hit analysis is altered. [Prefix counts](exploratory_sequence_prefixes.csv).
