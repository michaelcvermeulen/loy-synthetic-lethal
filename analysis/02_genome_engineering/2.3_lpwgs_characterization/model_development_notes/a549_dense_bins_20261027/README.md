> **Adoption update — 28 October 2026:** The three realignment atlases are now adopted into draft 2 as Figs 3.4–3.6. See the [current source and reproduction record](../../../figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md). The text and checks below describe the earlier analysis/preview state. The separate-track supplement has not been realigned.

# A549: broader 100-kb track versus 50-kb bins

Recommendation: use the expanded 100-kb preview. It contains 67 male-specific Yp bins instead of ten, with 100 euchromatic Yq bins, and preserves stronger per-bin support than the 50-kb version. Both complete three-panel previews are available; neither replaces selected exports or draft 2 text.

## What was changed

The generic difficulty-overlap annotation is retained as metadata rather than used as a blanket exclusion. Reference mappability >=30%, non-N sequence >=60%, anatomical exclusion of centromere/Yq12, a valid control denominator and observed high-quality control-read support remain required. No selection uses the LOY count or effect direction. The 100-kb display requires at least 50 primary MAPQ >=37 control reads per Y bin; actual minimum is 54 in Yp and 50 in Yq. This is a deliberately broader descriptive display, not a claim that every complex region is fully resolved.

## Resolution comparison

| Bin size / region | Bins | Span (Mb) | Median control reads | Bins with no LOY reads | Median corrected ratio |
|---|---:|---:|---:|---:|---:|
| 100 kb / MSY Yp | 67 | 6.70 | 123 | 20.9% | 2.10% |
| 100 kb / euchromatic Yq | 100 | 10.00 | 140 | 20.0% | 2.69% |
| 50 kb / MSY Yp | 132 | 6.60 | 64 | 50.0% | 0.76% |
| 50 kb / euchromatic Yq | 194 | 9.70 | 71 | 43.3% | 2.24% |

The broader 100-kb track spans 6.70 Mb of the 8.11-Mb male-specific Yp interval (about 83%). The 50-kb version spans 6.60 Mb, so it adds spatial subdivisions rather than more genomic coverage or sequencing information. The minimum control-read threshold is scaled from 50 to 25 to compare similar read densities. At 50 kb, the median control count falls from 123 to 64 on Yp and half the retained Yp bins contain no primary LOY reads. Keeping a 50-read minimum at 50 kb retains 91 Yp bins but covers only 4.55 Mb; that sensitivity is included in `resolution_comparison.csv`.

The 50-kb bin-median Yp estimate falls to 0.76%, compared with 2.10% at 100 kb. This is a resolution/selection/fitting sensitivity in sparse data, not evidence of additional biological Y loss. The expanded 100-kb version provides a denser track without that increase in zero-read bins.

## Validation

- Exact 50-kb GC/non-N, mappability and difficulty annotations were recomputed from the existing aligned reference and original interval sources, not copied from parent bins.
- Recombining 50-kb annotations reproduces the frozen 100-kb annotation in all 31,156 bins outside masked Y-PAR intervals (tolerance 1e-5 percentage points; see `resolution_validation.json`).
- Summing 50-kb QDNAseq read-start counts reproduces the 100-kb counts exactly in all 31,185 genomic bins and both samples. Primary aligned-base coverage on Y is also exactly conserved when reaggregated.
- The same completed BAMs, MAPQ threshold and QDNAseq read-count settings were used. The 50-kb analysis was fitted independently with seed 20261027. PAR measurements retain their fixed autosomal depth normalisers.
- The ideogram and all 88 RNA observations are unchanged. Both figures were visually inspected; selected exports remain unchanged by checksum.

## Reproduction and files

- `WGS/reanalysis/Y_PAR_masked_20261027/build_50kb_comparison.py` builds 50-kb annotation and primary-read coverage under its `a549_50kb_comparison/` folder, using the existing WSL runtime/reference. No realignment is required.
- `fit_50kb.R` counts and fits the matched pair at 50 kb. The broader 100-kb fit is the existing `a549_par_track_relaxed_20261027/map30_no_difficulty` result.
- Run `prepare_previews.py 100` and `prepare_previews.py 50`, then the generated `make_a549_dense_preview.py A549 --output-dir <variant folder>` for each, followed by `compare_resolution.py` and `finalize.py`.
- Complete PDF/PNG/SVG previews, captions and source snapshots: `A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/figures/chapter3/drafts/a549_dense_bins_20261027/100kb` and `/50kb`.
