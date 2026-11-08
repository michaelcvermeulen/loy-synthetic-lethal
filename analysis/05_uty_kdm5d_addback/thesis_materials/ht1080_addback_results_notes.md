# HT1080 add-back Results source notes

## Final source

The frozen thesis source is `RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026`. It supersedes the prior designation of `reannotate_v105` for Chapter 3 numerical add-back results. The older notebooks and `SEPT_2026_REANALYSIS` remain provenance and sensitivity sources.

The analysis reads the verified full count matrix and GENCODE v49 annotation produced from the original 18 featureCounts tables, the corrected sample manifest and relabelling record, and the August-only frozen signature from `HT1080/SEPT_11_2026_REANALYZE`. SHA256 hashes are recorded in `provenance/input_manifest.csv`; the earlier audit table includes hashes of the original count files.

## Author-confirmed experimental structure

- Empty-vector, UTY, KDM5D and combined UTY–KDM5D populations came from separate transductions of the same starting Y-depleted pool.
- Sample indices 1–3 identify matched collection occasions across all six populations.
- Each group represents one established population, so the three RNA samples do not provide independent transduction replication.
- Western blotting was completed for UTY and KDM5D in each sample and showed the expected abundance. Original images, lane identities, loading controls and exposure information remain to be linked before final figure insertion.

## Retained analysis

- Primary background: 23,992 expressed autosomal/non-PAR X genes.
- Primary model: limma-voom, `~ 0 + group + collection_occasion`, robust empirical Bayes.
- Complete output: 15 pairwise contrasts, TREAT at log2 fold change 0.5, unblocked sensitivity, plasmid factorial terms and leave-one-occasion-out sensitivity.
- Baseline reversal: 100-gene August-only frozen signature; 50 genes per direction. All 100 were present in May counts and 75 passed the May primary filter.
- Gene-set inference: directional `fry`, corrected across 12 prespecified signature-arm tests.
- Gene-level distance-to-control summaries are descriptive and retain control crossings; no binary rescue class is inferred.

## Central results used in prose

- Primary selected-gene counts versus vector: UTY 1,155; KDM5D 0; combined 2,153. Corresponding TREAT counts: 240, 0 and 197.
- Direct combined versus UTY: 0 selected genes. Combined versus KDM5D: 2,868.
- Frozen signature score differences versus vector: UTY −0.108 (95% CI −0.149 to −0.067; P = 1.6e−4); KDM5D 0.031 (−0.011 to 0.072; P = 0.13); combined −0.085 (−0.127 to −0.044; P = 9.7e−4).
- Combined versus UTY score difference: 0.023 (−0.018 to 0.064; P = 0.25).
- Both frozen signature arms changed in the corrective direction for UTY and combined versus vector; neither passed corrected significance for KDM5D versus vector.
- Of the 75 May-filtered signature genes, 49 were descriptively closer to sgCTRL after UTY, 59 after combined re-expression and 26 after KDM5D.

- Naive Y-depleted versus vector moved the frozen score away from sgCTRL: +0.114 (0.073 to 0.156; P = 1.1e−4); the August-up arm increased (fry FDR = 9.8e−6).

## Pathway analysis (added 21 September 2026)

- CAMERA (`code/run_camera_report.R`) and correlation/specificity sensitivity (`code/run_pathway_sensitivity.R`; outputs in `outputs/camera_sensitivity/`).
- UTY versus vector at inter-gene correlation 0.01: 19 Hallmark, 356 Reactome and 275 GO BP sets at FDR <0.05. Naive Y-depleted versus vector: 18, 546 and 538.
- Signed-z Spearman correlation with naive Y-depleted versus vector: UTY 0.93/0.88/0.79; combined 0.96/0.95/0.88; KDM5D 0.24/0.47/0.15. UTY sets shared in the same direction: 14/19, 311/356, 212/275; none in the opposite direction.
- The median estimated inter-gene correlation was 0.034. At 0.03, UTY versus vector retained 6/228/2 sets; with per-set estimated correlation, no set passed in any contrast.
- There were 17 small UTY-only sets that remained significant at both fixed settings, mainly IL-6 family, collagen-receptor, angiogenesis and coagulation sets (see `uty_sets_not_shared_with_naive_LOY.csv`). They are hypothesis-generating only.
- Interpretation: pathway shifts mark the empty-vector/KDM5D pair as distinct from the other Y-depleted populations. Pathway shifts are not UTY-specific, whereas the frozen-signature shift is.

## Interpretation limits

August-only selection removes May-dependent gene selection but does not create an independently derived validation model. The vector and naive Y-depleted populations differ in transduction/selection history and regional Y expression. The four plasmid groups are one population each. Shared empty-vector estimates induce dependence among trajectory summaries. RNA and protein abundance do not establish catalytic activity. There is no add-back ATAC dataset, and molecular reversal does not establish functional dependency rescue.

