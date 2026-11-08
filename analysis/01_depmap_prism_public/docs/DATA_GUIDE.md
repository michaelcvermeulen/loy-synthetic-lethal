# Data guide: supplied DepMap and drug-response collection

Audited 2026-09-14. This guide is based on the local files, their actual headers
and row identifiers, and their supplied READMEs. No biological analysis or LOY
classification was run.

## Location and organization

The collection remains at `/global/home/hpc5173/projects/PRISM/data/CCLE/26Q1`.
The redundant nested `26Q1/26Q1` directory has been removed. All 197 original
files, totaling 21,468,298,379 bytes (19.99 GiB), were retained.

```text
26Q1/
├── raw/
│   ├── depmap_public_26q1/
│   │   ├── metadata/        # Model, ModelCondition, Gene, profiles and screen maps
│   │   ├── crispr/          # Chronos effects, dependency probabilities, screen effects
│   │   ├── expression/      # Counts, gene TPM and transcript TPM
│   │   ├── copy_number/     # WGS, WES, segments and portal-transformed matrix
│   │   ├── mutations/       # Variant tables, MAF, damaging/hotspot matrices
│   │   ├── fusions/
│   │   ├── signatures/     # Genomic signatures, subtypes and microsatellites
│   │   └── documentation/  # Original 26Q1 README
│   ├── drug_response/
│   │   ├── prism/24Q2/primary/
│   │   ├── prism/20Q2/{primary,secondary}/
│   │   ├── ctrp/           # Published releases plus original ZIP archives
│   │   ├── gdsc/
│   │   ├── sanger_combinations/
│   │   └── ccle_np24/
│   ├── proteomics/harmonized/
│   ├── rnai/combined_demeter2/
│   ├── gene_sets/ssgsea_24Q2/
│   ├── epigenomics/global_chromatin_profiling/
│   └── studies/{depmap_nextgen,sanger_organoid}/
├── derived/{ctrp,sanger_organoid}/
├── results/sanger_organoid/figures/
├── references/sanger_organoid/
├── scripts/imported/ctrp/
├── docs/
└── archive/duplicate_downloads/
```

Here, `raw/` means the supplied source tables, many of which were already
processed by the data providers. It does not imply raw sequencing reads.
Existing CTRP mappings/AUC tables and organoid LOY results are in `derived/`;
the organoid figures are in `results/`. These inherited analyses were not rerun.

Sixteen byte-identical extra copies (1.82 GiB) were moved to the duplicate archive.
They were not deleted. `data_file_manifest.csv` records the old path, new path,
size, SHA-256, and canonical counterpart of each duplicate. Old paths are relative
to the former nested `26Q1/` directory; new paths are relative to the collection root.
Those mappings also specify how to reverse the relocation.

## Main DepMap tables

The following dimensions are observed locally. Feature counts exclude identifiers
and flags; total column counts are recorded in `table_summary.csv`.

| Table | Observed coverage | Intended role |
|---|---|---|
| `Model.csv` | 2,154 unique models, 49 columns | Model annotations and main join registry |
| `Gene.csv` | 44,083 rows, 56 columns | HGNC symbols, genomic locations and identifier mappings |
| `ModelCondition.csv` | 2,826 conditions covering 2,044 models | Model-condition mapping |
| `OmicsProfiles.csv` | 4,819 profiles covering 2,044 models | Profile metadata across assay types |
| `CRISPRGeneEffect.csv` | 1,208 models × 18,531 genes | Integrated Chronos effects |
| `CRISPRGeneDependency.csv` | 1,208 models × 18,531 genes | Gene dependency probabilities; a different measure from effects |
| `ScreenGeneEffect.csv` | 1,458 screens × 18,531 genes | Screen-level effects; map screens before model-level comparisons |
| `OmicsExpressionRawReadCountHumanAllGenes.csv` | 1,775 profiles × 60,649 genes; 1,719 models | Raw gene counts |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 1,775 profiles × 19,215 genes; 1,719 models | Log-transformed protein-coding gene expression |
| `OmicsExpressionTPMLogp1HumanAllGenes.csv` | 1,775 profiles × 60,649 genes | All-gene expression |
| `OmicsExpressionTranscriptTPMLogp1HumanAllGenes.csv` | 1,775 profiles; 237,018 total columns | Transcript expression; this is a much larger table |
| `OmicsCNGeneWGS.csv` | 1,132 profiles × 19,955 genes; 1,118 models | WGS gene-level copy number in linear scale |
| `OmicsCNGeneMC_WES.csv` | 1,904 profiles × 37,605 genes; 1,698 mapped models | Older WES copy number, indexed by model condition |
| `OmicsSomaticMutations.csv` | 1,172,688 variant records; 1,968 models | Long mutation table |
| Damaging / hotspot mutation matrices | 3,044 profile rows each; 1,968 models | Mutation feature matrices |
| `OmicsGlobalSignatures.csv` | 3,044 profile rows; 1,968 models | MSI, LOH, WGD, CIN, ploidy and aneuploidy fields |

The README describes `CRISPRGeneEffect.csv` as Chronos-integrated, copy-number
corrected, scaled, and screen-quality corrected. It describes WGS copy number as
linear-scale. The portal-transformed copy-number matrix has a different gene set;
do not assume the WGS, WES and portal matrices have identical columns or units.

## Joins and profile selection

- `ModelID` uses `ACH-` identifiers and joins to `Model.csv`. `PatientID` identifies
  the donor and can be shared by multiple models; it is not a replacement key.
- CRISPR matrices carry the `ACH-` identifier in their unnamed first column.
- The downloaded RNA expression and WGS gene matrices start with an unnamed
  numeric row index, then `SequencingID`, `ModelConditionID`, `ModelID`,
  `IsDefaultEntryForMC`, and `IsDefaultEntryForModel`. These are six metadata
  columns, not gene measurements. Their actual headers should take precedence
  over generic profile-column descriptions in the README.
- For one RNA or WGS profile per model, filter `IsDefaultEntryForModel == "Yes"`.
  RNA has 1,719 default and 56 nondefault rows; WGS has 1,118 default and 14
  nondefault rows. Neither default subset repeats a ModelID.
- WES has `ModelConditionID` and `IsDefaultEntryForMC`, without a direct ModelID
  column. Join through `ModelCondition.csv`; retain condition provenance and
  decide explicitly how to handle multiple conditions for one model.
- `ScreenID` uses `SC-` identifiers; use `CRISPRScreenMap.csv` and
  `ScreenSequenceMap.csv` to map screens. `pDNA` in the sequence map is a control
  label and is excluded from model counts.
- Many gene columns use `SYMBOL (EntrezID)`. Preserve the complete original name
  and a parsed identifier in analysis outputs, since symbol parsing can collapse
  distinct columns.
- Fusions, mutation records and copy-number segments are long tables. Repeated
  ModelIDs in these tables are expected and must not be removed indiscriminately.

## LOY-relevant coverage

The model registry contains 1,143 Male, 869 Female, 119 Unknown, and 23 blank Sex
annotations. These are source metadata, not chromosome-Y measurements.

| Available assays | Models | Marked Male | Prostate lineage |
|---|---:|---:|---:|
| CRISPR effects | 1,208 | 658 | 10 |
| Gene expression | 1,719 | 905 | 13 |
| WGS copy number | 1,118 | 608 | 6 |
| CRISPR + default gene expression | 1,140 | 619 | 10 |
| CRISPR + default gene expression + default WGS | 852 | 472 | 6 |

There are 15 prostate-lineage models in the registry. Seven are annotated
`Non-Cancerous`; lineage membership alone does not define a prostate cancer cohort.
See `prostate_model_availability.csv` for individual models and assay availability.

All 12 Y genes listed in the existing project's X/Y paralog table occur in the
all-gene raw-count matrix. Ten occur in the CRISPR effect and WGS gene matrices;
`PRKY` and `TXLNGY` are absent from those two feature sets. This establishes feature
availability only, not expression quality or LOY status. No LOY calls were made.

## Other collections and cross-release issues

| Collection | Observed main table | Join considerations |
|---|---|---|
| PRISM 24Q2 extended primary | 6,790 treatment rows × 919 model columns | Orientation is treatments × models; 4 ACH IDs are absent from the 26Q1 registry |
| PRISM 20Q2 secondary curves | 753,778 rows; 737 model IDs | Retain concentration/replicate and compound identifiers; 2 model IDs are absent from 26Q1 |
| CTRP v2.0 mapped AUC | 387,130 rows; 841 model IDs | Existing derived mappings are supplied; inspect match-tier fields and imported QC |
| GDSC dose response | 387,626 rows; 973 ACH IDs | `ARXSPAN_ID` carries ACH IDs; 10 are absent from the 26Q1 registry |
| Harmonized CCLE mass spectrometry | 375 models; 12,559 total columns | Protein identifiers require the supplied mapping table |
| Harmonized CCLE RPPA | 899 model IDs; 145 total columns | One model ID is absent from the 26Q1 registry |
| DepMap NextGen study | 1,592 screen metadata rows; 1,357 models | A separate publication dataset; 134 IDs are absent from the 26Q1 registry |
| Sanger organoid study | Publications, supplementary tables and inherited LOY outputs | Keep publication IDs and prior analysis provenance separate from DepMap model IDs |

These are separate dataset versions even though they were downloaded together.
The supplied PRISM 24Q2 README describes the extended matrix as incorporating the
previous primary screen. The supplied CTRP documentation describes v2.1 and v2.2
as subsets of v2.0. Do not count overlapping releases as independent replication.

Unmatched IDs are recorded in `data_profile.json`; they were not silently dropped
or reassigned. Availability tables are anchored to the current 2,154-model registry,
so models outside that registry appear in the unmatched lists rather than as extra
rows. One legacy NP24 CSV requires Windows-1252 for its header/sample preview;
the data were not transcoded.

## Completeness and validation

The supplied DepMap README has 85 exact file section headings. This collection
contains 32 of them; 53 are absent. The absent files include supplemental CRISPR
read counts, guide maps, QC reports, inferred parameters, and alternate expression
quantifications. `release_completeness.csv` lists every expected filename.

Every supplied file was hashed before and after relocation; all 197 SHA-256
digests matched. No partial-transfer extensions or zero-byte files were present.
Headers and first records were inspected for CSV/MAF inputs. The selected tables
in `table_summary.csv` were scanned through their complete contents for row counts
and identifier coverage; the remaining files have inventory/header inspection only.
Numeric matrix distributions and biological plausibility were not tested.

These checks establish unchanged local content and the stated parsing checks.
Official upstream checksums were unavailable, so they do not certify that every
browser download matches its publisher's original file.

## Existing code and next analysis prerequisites

The user subsequently supplied six Qi et al. supplementary workbooks. They are
organized in `raw/studies/qi_2023/supplementary_tables/` under the collection root.
See [the Qi table catalog](qi_2023/README.md): S1A contains 5,041 male TCGA cases
with Y annotations and S1B contains 5,406 female TCGA cases with X annotations.
S6 contains gene-level dependency statistics. No individual cell-line X/Y status
table was found. Three TCGA annotation extracts are in the working project's
`data/processed/qi_2023/`; the catalog/extraction ran successfully on SLURM.

`00_utils.R` now sources `config/data_paths.R` for the organized DepMap, PRISM,
GDSC and proteomics paths. The project's `data/depmap/26Q1/` is a flat symlink view,
so existing calls such as `file.path(DEPMAP_DIR, "Gene.csv")` resolve correctly.
Run the notebook from the project root.

After the original inventory, the user supplied
`/global/home/hpc5173/projects/PRISM/data/CCLE/26Q1/depmap_loy_classifications.RDS`.
It is linked at `data/depmap_loy_classifications.RDS`, the path the existing notebook
expects. The data frame contains 896 unique ModelIDs, `Y_score1`, `LOY_binary`,
and three Oncotree fields. There are 338 LOY and 558 XY labels; all 896 IDs match
Male entries in the 26Q1 Model registry, and no Y scores are missing. The file has
no independent X-status field or extra provenance attributes establishing a Qi
et al. source. Its classification method and source release remain unverified;
the XY label alone does not establish measured X copy number or inactivation.
See `loy_classification_file_audit.json`. This subsequently supplied file is outside
the original 197-file relocation manifest.

The notebook's count loader still needs review: its metadata-exclusion list does not remove the unnamed
numeric CSV row-index column observed in 26Q1. No differential analysis was run.

`DEPMAP_ALT`, `TCGA_DIR` and `LOY_DIR` in `00_utils.R` retain their old Windows
locations because their intended external datasets have not been identified.
Use `RNAI_FILE` and `SSGSEA_FILE` in the new configuration for the corresponding
downloaded legacy tables. The imported CTRP script is preserved verbatim and
also retains Windows paths; port it before rerunning its analysis.

No `Rscript` executable was on PATH during the initial audit. A later lightweight
inspection of the supplied RDS used
`/global/home/hpc5173/miniconda3/envs/seurat/bin/Rscript` successfully. The full R
analysis and path configuration have not been executed. Python profiling completed
with pandas; its source and run log are included for reproducibility.

## Additions 2026-10-04: public inputs for the LOY extension stages

Twelve files were downloaded on 2026-10-04 by
`scripts/loy_reanalysis_v1/fetch_extension_inputs.sh` for stages 25 to 30 of
`loy_reanalysis_v1`. They are outside the original 197-file relocation manifest and
`docs/SHA256SUMS`; their SHA-256 sums are in `docs/SHA256SUMS_extension_inputs`
(check with `sha256sum --check docs/SHA256SUMS_extension_inputs`). Each file was
verified against the publisher's checksum at download (figshare MD5; the SHA-256 in
the Git LFS pointer for BeatAML2) and is stored read-only. Source URLs and checksums
are also recorded in `scripts/loy_reanalysis_v1/extension_parameters.json`.

| Location | Files | Source and licence |
|---|---|---|
| `raw/studies/sanger_organoid/figshare/` | `CNV-segments_organoids.tsv` (PURPLE copy-number segments, chr1 to chrY, 256 organoids plus longitudinal samples), `data_CNV-genes_ALL-organoids.csv` (MSI, WGD, focal and arm-level calls), `data_covariates_ALL-organoids.csv`, `control-genes_organoids.csv`, `position_genes.csv` | Herranz-Ors et al., Nature 2026, figshare doi 10.6084/m9.figshare.28339340.v1, CC BY 4.0 |
| `raw/drug_response/beataml2/` | `beataml_wv1to4_clinical.xlsx` (942 specimens, 805 patients), `beataml_probit_curve_fits_v4_dbgap.txt` (drug response curve fits, 166 inhibitors), `beataml_waves1to4_norm_exp_dbgap.txt` (normalised RNA-seq, 707 samples), `beataml_wes_wv1to4_mutations_dbgap.txt`, `beataml_drug_families.xlsx`, `beataml_waves1to4_sample_mapping.xlsx`, `LICENSE.txt` | Bottomly et al., Cancer Cell 2022 (BeatAML2), `github.com/biodev/beataml2.0_data` at commit 6e403e1, CC BY 4.0 |

These are local integrity checks against the publishers' checksums; they do not
verify the biological content. The organoid supplementary tables already in
`raw/studies/sanger_organoid/supplementary_tables/` are unchanged.

## Additions 2026-10-05: DepMap NextGen inputs for stages 31 to 33

Six files were copied on 2026-10-05 by
`scripts/loy_reanalysis_v1/stage_nextgen_inputs.sh` from a folder of browser
downloads (`/global/home/hpc5173/projects/PRISM/data/nextgen/`, dated 2026-08-10,
left as it was). The DepMap portal answers scripts with a verification page, so
these files cannot be fetched by script. They are outside the original 197-file
relocation manifest and `docs/SHA256SUMS`; their SHA-256 sums are in
`docs/SHA256SUMS_nextgen_inputs` (check with
`sha256sum --check docs/SHA256SUMS_nextgen_inputs`), and their sources and
checksums are recorded under `nextgen_inputs` in
`scripts/loy_reanalysis_v1/extension_parameters.json`. All are stored read-only.

| Location | Files | Source, licence and check |
|---|---|---|
| `raw/studies/depmap_nextgen/` | `next_gen_expression.csv` (log2 TPM+1, protein-coding genes, 309 NextGen models), `next_gen_copy_number.csv` (gene-level copy number relative to ploidy, 291 models) | DepMap portal file set "NextGen Model Manuscript 2026" (Neiswender et al., Nature 2026). The portal publishes no checksum; the SHA-256 is that of the copy supplied |
| `raw/studies/depmap_nextgen/figshare/` | `model_metadata.csv` (1,999 models, 314 of them NextGen) | figshare 10.6084/m9.figshare.29472362.v1, CC BY 4.0; matches the MD5 figshare publishes |
| `raw/studies/depmap_nextgen/supplementary_tables/` | `41586_2026_10843_MOESM3_ESM.xlsx` (Supplementary Table 1: NextGen model metadata with recorded sex, a data dictionary, media formulations) | Neiswender et al., Nature 2026, 10.1038/s41586-026-10843-7; SHA-256 of the copy supplied |
| `references/depmap_nextgen/` | `s41586-026-10843-7.pdf` (the paper, open access), `29472362.zip` (the figshare record as downloaded, 36 files) | As above; SHA-256 of the copies supplied |

The other files of the same portal file set were already in
`raw/studies/depmap_nextgen/` and are unchanged; six files in the download folder
are byte-identical to them. `screen_metadata.csv` exists in two versions: the
portal copy already in the collection is the one the analysis uses, and the copy
inside the figshare archive differs from it only in three descriptive columns
(`OncotreeCode`, `PlateCoating`, `HasCNData`).

These are local integrity checks. Only `model_metadata.csv` is checked against a
publisher's checksum, and none of this verifies the biological content.
