# A549 X-PAR tracks and relaxed bin filter

Revised preview only. Selected figure exports and draft 2 text/captions are unchanged; the prior preview is preserved separately.

## What changed

X-PAR1 and X-PAR2 coverage now appears directly in panel a. The hatched PAR exclusions and separate inset were removed. The completed Y-PAR-masked BAMs were reused: the change is to the displayed PAR measurement, not to the reference or alignments. Centromere/Yq12 exclusions remain. All 88 RNA observations and the ideogram are unchanged.

## Filter comparison and choice

The original difficulty-overlap rule, rather than the 55% mappability threshold alone, caused most exclusions. Lowering mappability to 30% while retaining difficulty overlap <20% still leaves 58 male-specific Y bins. The chosen moderate relaxation keeps mappability >=30%, difficulty overlap <50% and non-N bases >=60%, with at least 50 primary MAPQ >=37 control reads per full 100-kb Y bin. All 81 chosen Y bins actually have at least 56 such control reads. Eligibility is not selected on LOY counts or effect direction.

| Filter | Yp bins | Yq bins | Median Yp (% control) | Median Yq (% control) |
|---|---:|---:|---:|---:|
| Mappability >=55%, difficulty <20% | 7 | 51 | 3.77 | 2.73 |
| Mappability >=30%, difficulty <50% — chosen | 10 | 71 | 2.57 | 2.76 |
| Mappability >=30%, difficulty <80% | 21 | 76 | 2.06 | 2.76 |
| Mappability >=30%, no difficulty-overlap exclusion | 67 | 100 | 2.10 | 2.69 |

All rows additionally apply the same anatomical exclusions and control-read support requirement. The chosen filter adds 23 Y bins (58 to 81; ten Yp and 71 Yq) while retaining a difficulty-overlap safeguard. More permissive alternatives were examined as sensitivity analyses, not chosen to achieve a particular retention value. The main strong-depletion conclusion persists; the Yp median is sensitive to which small set of bins is retained. This comparison supports a practical preview threshold, not a formally optimised or independently validated filter.

QDNAseq was refitted separately for each alternative on the same matched raw counts, preserving its specified mask throughout and using seed 20261027. The unchanged BAM read-count settings include the historical QDNAseq defaults; the additional displayed-bin support check and PAR diagnostic use primary reads. PAR bins were excluded from the QDNAseq correction fits because their original mappability scores describe the duplicated reference. The autosomal normaliser for the separate PAR depth ratios is held fixed across these displays.

## PAR values and coordinates

All 28 X-PAR bins pass the observed-control-coverage rule. The bin-median ratios are 66.83% for PAR1 and 63.65% for PAR2. The earlier inset used whole-region aggregate coverage ratios (67.41% and 63.26%); the small change reflects using medians of bin ratios for the new track. These are combined X/Y signals. They are displayed at fractional positions within labelled PAR spans and must not be read as exact X-to-Y coordinate liftovers.

## Files and reproduction

- `refit_filters.R`: staged QDNAseq fits; original RDS objects and global annotation were not changed.
- `compare_filters.py`: compares retained bins and results; `filter_result_comparison.csv` and `filter_bin_audit.csv` preserve all alternatives.
- `prepare_preview.py`: derives an isolated renderer, preserving b/c code exactly.
- X-PAR and primary Y control coverage: `A:/Dropbox/LOY_paralogs/projects/WGS/reanalysis/Y_PAR_masked_20261027/a549_panel_tracks`; produced by `a549_panel_tracks.py` without realignment.
- New PDF/PNG/SVG, caption and source snapshots: `A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/drafts/a549_xpar_relaxed_20261027`.
- Reproduce from the projects directory: `refit_filters.R`, `compare_filters.py`, `prepare_preview.py`, generated `make_a549_xpar_preview.py A549 --output-dir <new preview folder>`, then `finalize.py`.
