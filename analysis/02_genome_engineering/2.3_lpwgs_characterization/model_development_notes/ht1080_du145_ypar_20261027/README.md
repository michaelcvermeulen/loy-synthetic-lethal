> **Adoption update — 28 October 2026:** The three realignment atlases are now adopted into draft 2 as Figs 3.4–3.6. See the [current source and reproduction record](../../../figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md). The text and checks below describe the earlier analysis/preview state. The separate-track supplement has not been realigned.

# HT1080 and DU145: completed realignments and atlas previews

The four samples authorised on 27 October 2026 are complete. The six additional DU145 treatment samples and A549-BC14 were not processed. The legacy `LNCaP` filenames identify the DU145 samples used by the existing atlas.

| Cell line | Original BAM sample ID | Read pairs |
|---|---|---:|
| HT1080 | HT1080-sgCTRL-D0 | 7,125,059 |
| HT1080 | HT1080-LOY-D0 | 6,664,446 |
| DU145 | LNCaP-sgCTRL | 6,826,656 |
| DU145 | LNCaP-LOY-2X | 6,061,770 |

Every primary read name, mate, sequence and quality was preserved through BAM-to-FASTQ recovery and realignment. Original BAM checksums are unchanged. BAM integrity, reference chromosome lengths, paired FASTQs, QDNAseq mask preservation, source data and figure exports passed their checks.

## Atlas previews

- [HT1080 Fig 3.10 PDF](A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/drafts/ht1080_du145_ypar_20261027/ht1080/ch03_main_ht1080_lpwgs_chrY_atlas.pdf) — [caption](A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/drafts/ht1080_du145_ypar_20261027/ht1080/preview_caption.md) and [analysis](A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/research/model_development/ht1080_du145_ypar_20261027/HT1080/README.md).
- [DU145 Fig 3.14 PDF](A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/drafts/ht1080_du145_ypar_20261027/du145/ch03_main_du145_lpwgs_chrY_atlas.pdf) — [caption](A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/drafts/ht1080_du145_ypar_20261027/du145/preview_caption.md) and [analysis](A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/research/model_development/ht1080_du145_ypar_20261027/DU145/README.md).

PNG and SVG exports, numerical tables, code and validation accompany each PDF. Each figure retains its own ideogram and RNA data, including the pre-existing DU145 4.8-Mb Yp split. Both use the broader 100-kb Y-bin display selected after the A549 comparison, with X-mapped PAR measurements directly in panel a. The selected figures, central captions and draft 2 have not been replaced.

## Effect of realignment on male-specific Y

| Cell line | Region | 100-kb bins | Original alignment | New alignment |
|---|---|---:|---:|---:|
| HT1080 | MSY Yp | 67 | 42.62% | 42.65% |
| HT1080 | euchromatic Yq | 101 | 0.00% | 0.00% |
| DU145 | distal MSY Yp | 23 | 54.37% | 54.38% |
| DU145 | proximal MSY Yp | 44 | 6.11% | 6.11% |
| DU145 | euchromatic Yq | 101 | 2.20% | 2.20% |

The old/new columns compare the same displayed bins under the same broad-filter, matched-pair analysis. They separate the alignment update from the earlier expansion of the bin display. Zero medians are descriptive at this sequencing depth and do not establish complete DNA absence; the plots label values below 1% as `<1%`.

## PAR interpretation and methods

Only the Y copies of PAR are masked in the alignment reference; the PAR measurements are shown. Teal tracks use X-mapped combined X/Y signal, normalised to the median primary-read depth of the same 21,446 autosomal bins per sample. They cannot measure Y-specific PAR retention. X-PAR positions are projected into labelled PAR spans, not treated as an exact coordinate liftover.

Male-specific Y bins use the rerun QDNAseq correction with reference non-N sequence >=60%, mappability >=30%, valid estimates in both samples and >=50 primary MAPQ >=37 control reads per 100-kb bin. Generic difficulty overlap remains annotated but does not automatically exclude the bin. Centromere, Yq12 and boundary-crossing bins are excluded. No minimum LOY count or effect-direction filter is used. PAR uses direct primary-read depth because the old PAR mappability annotation does not apply after Y-PAR masking.

The four-sample driver is `../../../../WGS/reanalysis/Y_PAR_masked_20261027/run_atlas_four.sh`; it reuses the verified A549 reference and index. Five BWA-MEM2 workers per sample and a fixed 160-million-base batch permit two concurrent alignments. The `fit_pairs.R`, `prepare_inputs.py`, `build_renderer.py`, `compare_and_finalize.py` and `finish_summary.py` scripts record the analysis and figure steps. Per-cell reports contain reproduction commands and detailed support/coverage diagnostics.
