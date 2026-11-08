# Cell-line profile and dependency comparison figures

## Scope

The author requested an HT1080 profile covering origin, cancer type, doubling time, patient details, mutations, Y expression and specific dependencies, followed by a comparison of HT1080, A549, DU145 and H2452 against XY and LOY lines for DDX3X, EIF1AX, RPS4X and ZFX. These figures support model selection in Results 3.4.6. The author confirmed ZFX as the fourth comparison gene on 17 September 2026; the earlier UTX/KDM6A versions were archived before the figures, data and captions were corrected.

Use baseline public-data measurements for model selection. Published cell-line properties and public profiles do not establish the corresponding measurements in the author's individual stocks or edited populations. Keep the experimental LOY screen for Chapter 4 unless the author requests its inclusion here.

## Background sources

Curated fields and claim-specific links are saved in [ht1080_background.json](../../figures/chapter3/source_data/cell_line_profiles/ht1080_background.json). DSMZ reports a doubling time of approximately 30 h under its stated conditions. Cellosaurus also cites an earlier 26-h estimate. Do not combine these estimates into an apparent experimental interval or label them as the author's own measurements.

## Data release reconciliation

The existing project scripts point to `A:/DEPMAP/DepMap`. Its `README.txt` identifies **DepMap Public 25Q3**; its gene-effect matrix and model metadata were saved in January 2026. The saved `depmap_tcga_loy/data/depmap_loy_classifications.RDS` contains 896 male-derived models (558 XY and 338 LOY), created in March 2026 by the expression/GMM caller in `01_Boruta_DepMap.Rmd`. This is not a verified 26Q1 input bundle.

The author supplied the HPC source `/global/home/hpc5173/projects/PRISM/data/CCLE/26Q1/` and expected transfer completion within approximately 30 minutes on 17 September 2026. A local directory `A:/DEPMAP/DepMap/26Q1` is present but was empty when checked after this message; its relationship to the HPC source is not yet verified. Wait for complete files and inspect the release README and supplied updated Y calls before generating the quantitative panels. Do not relabel 25Q3 files or merge their statistics into 26Q1 plots.

**Local path supplied by the author:** `H:/PRISM/data/CCLE/26Q1`. The release is organised by assay under `raw/depmap_public_26q1/`, with manifest and analysis documentation under `docs/`. The Model.csv and protein-coding expression files matched the HPC SHA-256 manifest on initial inspection. CRISPR and derived analysis files were still in transit. The supplied documentation identifies the reviewed calls at `derived/sex_chromosome_review_v1/sex_chromosome_calls_reviewed.csv` and the later vulnerability outputs under `derived/y_dosage_vulnerabilities_v1/final/`. Inspect these before choosing the plotted classification.

The local 26Q1 registry resolves HT1080 to ACH-000054, A549 to ACH-000681, DU145 to ACH-000979 and NCI-H2452 to ACH-000092. The HT1080 `PrimaryOrMetastasis` field says `Metastatic`, whereas repository descriptions identify the tumour as arising adjacent to the acetabulum. The figure therefore states the pelvic tumour origin without assigning primary/metastatic status. This conflict is retained for later source reconciliation. PRKY is absent from the protein-coding expression header; do not plot it as zero.

The portable extraction script is [prepare_cell_line_profiles.py](../../figures/chapter3/code/prepare_cell_line_profiles.py). It requires an explicit data directory, release label, Y-call CSV, call-column name and description of the caller. It checks model identities/default expression profiles, saves the inclusion table and source hashes, and detects changes to files while they are read. It does not infer new LOY calls. Use a fresh staging output directory, inspect the extracted tables, and then copy approved source data to the canonical figure source folder. The source folder's `ht1080_background.json` is independent of the data release.

The working HT1080 dependency display selects up to 16 strong dependencies (gene effect < -0.5), restricting to genes below -0.5 in at most 20% of other screened models. It ranks their differences from the other-model median. The author requested expansion from eight to 16 on 17 September 2026; 103 genes qualify under the unchanged criteria. No common-essential control list was available in the supplied bundle, so no such exclusion is claimed. This is a descriptive display choice for the new profile, not a historical selection method or evidence of a newly discovered dependency. Preserve the full ranking to make this choice reviewable. No inferential statistics are produced by the extraction or plotting scripts.

## Completed first drafts — 17 September 2026

Four figures were prepared from verified 26Q1 data: the main HT1080 profile, a supplementary four-model Y-expression heatmap, the main RNA-defined dependency comparison and a supplementary DNA-supported comparison. Active captions are in `Chapter3_figure_captions.md`; figure-specific notes contain reproduction instructions and source links. Earlier captions and the manifest were archived before each addition.

The main dependency comparison uses the executed RNA caller's 352 Y-expressing and 219 Y-depleted donor representatives with integrated CRISPR profiles. All four selected lines belong to the Y-expressing group. It is labelled by RNA state rather than asserting physical XY/LOY karyotypes. It pools partitions descriptively and reports observed mean differences. The DNA supplement uses the exact executed primary cohort (81 Y-retaining and 70 LOY), with source-matched adjusted effects and FDR-adjusted P values. A549, DU145 and H2452 are shown outside that stricter reference cohort, preserving their reasons for exclusion.

The `cohort_audit.csv` that arrived with the dependency results contains the reviewed calls and inclusion fields used by that analysis; this provided the necessary source before the standalone `sex_chromosome_calls_reviewed.csv` completed its transfer. There is no need to substitute the older 896-model classification file. The four selected lines' raw dependency scores and Y-status evidence are frozen in `figures/chapter3/source_data/model_dependency_comparisons/`.

HT1080's 26Q1 Chronos values are DDX3X +0.27, EIF1AX −1.15, RPS4X −2.53 and ZFX +0.10 (rounded). These values qualify any broad claim of resistance to all X-paralog knockouts. No Results prose or Word document was modified in this figure task.

### Original layout outline

The revised profile adds a DSMZ reference micrograph and expands the baseline dependency panel to 16 genes. Y-linked expression is ordered from Yp to Yq using frozen GRCh38 gene coordinates and linked to a full-length chromosome diagram. The four-model heatmap shares that order and diagram. The original reference image, attribution, coordinates and retrieval provenance are saved with the figure sources; superseded exports and captions are archived.

- **main: HT1080 profile.** Compact origin/patient/culture and reported-genotype panels, followed by Y-gene expression and selected baseline genetic dependencies with their reference distributions. Identify the source of each quantitative panel.
- **main: selected models in the LOY dependency landscape.** Four gene panels, XY before LOY, NPG blue/vermilion, the four selected models labelled directly and distinguished by shape. Display cohort sizes and source-matched effects/P values when available. Gene effects are on the Chronos scale, where more negative values indicate stronger dependency.
- Preserve the executed caller's actual meaning (RNA-based Y status, DNA-based LOY or integrated classification). Do not imply a normal diploid XY karyotype from an expression-only call. Account for missing/uncertain calls explicitly.
- Save plotted values, full selection criteria, source hashes, code and central captions. Use Open Sans and editable PDF/SVG exports with a PNG preview.
