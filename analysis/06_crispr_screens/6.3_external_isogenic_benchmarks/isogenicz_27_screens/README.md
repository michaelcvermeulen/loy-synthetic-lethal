# Schukken et al. 2026 aneuploidy screens — drugZ and IsogenicZ rescoring (22 October 2026)

This is an exploratory analysis requested by the author. Nothing here is cited in the thesis. It rescored the paired aneuploidy screens of Schukken et al. (bioRxiv 2026, doi 10.64898/2026.04.26.720636) with the HT1080/external-screen procedure. The counts are in `crispr_screen/Schukken_2026_aneuploid_screens/` (GitHub commit `4db9338`; see that folder's README).

## Design

- **26 comparisons, each aneuploid against the control screened in the same batch.** Pairings follow the authors' `MegaA_Analysis_6_KMS.R` (druggable-genome library, batches 1–7) and `WGL_all_Analysis_3_KMS.R` (Brunello).
- **Brunello (7):**
  - **Restored-disomy designs (5), reoriented in the script on 22 October 2026 (17:30):** the parental line carries the 1q or 7q gain and is placed in the aneuploid position, with the disomic derivative as reference. These are named `<line>_<arm>_gain_vs_Di` (A2058 1q and 7q, A2780 1q, AGS 1q, MCF10A 1q).
  - Vaco trisomy 13 vs control
  - HCT116 trisomy 13 vs control
- **Druggable genome (19):** see `contrast_summary.csv`.
- **Skipped:** the genome-wide DLD1 trisomy 13, because its control file (`count-DLD1_Cont-i.f.csv`) is not in the repository.
- **One population per arm, one D0 and one D10 library, no replicate endpoints.** Every comparison is therefore 1 vs 1:
  - "matched endpoint drugZ" compares the D10 libraries;
  - IsogenicZ normalises each D10 to its own arm's D0.
  - Both use the pinned IsogenicZ upstream (commit `e9327ab`, SHA-256 checked), pseudocount 5 and empirical-Bayes window 500.
  - Negative normZ means depletion in the aneuploid line.
  - With n = 1 per arm, FDRs rest entirely on the empirical-Bayes noise model and are less reliable than for the two-replicate Desjardins/Feng screens.

## Input repairs (all recorded in `contrast_summary.csv` and the script)

- **Excel date-converted symbols** (e.g. DEC1 → 45261, MARCH1 → 44986) were recovered from the gene prefix of the guide ID (`GENE-SEQ` in Brunello, `GENE_domain` in the focused library). When the two files disagreed, the non-numeric label was used. No numeric gene labels remain.
- **Control guides were removed under every naming used:**
  - Non-Targeting Control
  - `neg_N` (A2058 7qDi) and `neg1`–`neg100` (focused library)
  - blank-ID controls in the HCT116 WGL control file.
- **CDK12 was removed from every comparison**, because the authors note CDK12 contamination (over 40 million reads) in one batch. CDK12-removed controls were used where the authors provided them.
- **Duplicate coordinate IDs.** Coordinate-style IDs shared by readthrough genes (e.g. APITD1/APITD1-CORT) were numbered in file order, not summed.
- **Column order.** `count-MCF10A_WT.csv` lists D10 before D0; columns are matched by name, not position.

## Results

- **IsogenicZ against drugZ:** most comparisons agree (ρ 0.85–1.0). The exceptions are DLD1 trisomies 4 and 2/18 (batch 2; ρ 0.55; 0 drugZ calls against 28–43 IsogenicZ calls; D0 depth about twofold lower than the merged control), A2058 7q (0.68) and A2780 trisomy 8 (0.66).
- **Possibly unreliable comparisons:**
  - The A2780 druggable-genome batches have low D10 agreement between arms (ρ 0.34–0.65).
  - Vaco trisomy 13 (846 depleted, 0 enriched) and DLD1 12/17 (1 depleted, 151 enriched) are strongly one-sided, which points to a global growth difference between the arms.
- **UBE2H** (the paper's lead gene):
  - depleted at FDR <0.25 in SNU-1 clone 114 (trisomy 9/19; normZ −3.7, FDR 0.008 by both methods);
  - also in HCT116 trisomy 5 by drugZ only (FDR 0.14);
  - its normZ is negative in 17 of 26 comparisons (median −0.36).
- **HT1080 candidates (current orientation):**
  - EIF1AX, MED12 and PELO are depleted only in the one-sided Vaco trisomy 13 comparison.
  - CREBBP is depleted in HCT116 trisomies 3 and 8, SREBF1 in the DLD1 2/18 repeat, and UBE2G2 in two SNU-1 clones.
  - ATXN7L3 and TAF5L are depleted in A2780 1q gain, and TAF6L in AGS 1q gain (IsogenicZ only).
  - DDX3X and SREBF2 are never depleted.
  - Before the reorientation, SUPT20H scored in MCF10A as the disomic arm; that call no longer applies.
- **Recurrent hits and themes:** see the 24 October summary in chat. The recurrent hits include KMT2D, TP53, FBXW7, APC/C subunits, PRPF19/CDC5L and PKMYT1. A2058 7q gain is led by mitochondrial translation and OXPHOS genes.

## Outputs

These are in `figures/chapter4/source_data/schukken_2026_rescored/`:
- per-comparison `<id>_IsogenicZ.csv` and `<id>_matched_endpoint_drugZ.csv`
- `schukken_scores_long.csv.gz`
- `focus_gene_scores.csv`
- `contrast_summary.csv` (files, columns, guide counts, label repairs, D0/D10 agreement, depth, call counts)
- `manifest.json` (input hashes)

Reproduce with Anaconda Python: `python -I run_schukken_screens.py`.
