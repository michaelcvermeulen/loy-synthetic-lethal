# A549 Results 3.4.6.2 — source checks and drafting status

17 September 2026. The author requested development of the A549 subsection after adopting the integrated model-development structure. Existing A549 ddPCR summaries and lpWGS outputs were reviewed; no alignment, copy-number model or statistical test was rerun.

## Author-confirmed derivation

The author confirmed: “lpWGS was 2x”. The sequenced edited A549 population therefore received two successive CP + Q53 nucleofections. This clarification does not establish the recovery interval, identify it as a particular ddPCR specimen, or establish the editing history of every downstream RNA/screen sample. Those links remain open. The retained WGS script maps the comparison to A549-CTRL2 and A549-LOY; A549-BC14 is not included in this account.

## ddPCR values used

Source: [main group summary](../../figures/chapter3/source_data/a549_ddpcr_nucleofection_rounds/main_group_summary.csv), checked against [individual measurements](../../figures/chapter3/source_data/a549_ddpcr_nucleofection_rounds/main_measurements.csv). All ten group medians, quartiles, observation counts and reference-copy multiplications were independently recomputed and agreed. The [figure notes](../../figures/chapter3/ch03_main_a549_ddpcr_nucleofection_rounds_notes.md) document the raw exports, specimen mapping and exclusions.

| Condition | AMELY median copies/cell | Measurements | DDX3Y median copies/cell | Measurements |
|---|---:|---:|---:|---:|
| CP, 1x | 0.901086 | 2 | 0.875893 | 2 |
| Q53, 1x | 0.976561 | 2 | 0.968781 | 2 |
| CP + Q53, 1x | 0.058926 | 10 | 0.087788 | 9 |
| CP + Q53, 2x | 0.014522 | 4 | 0.007863 | 3 |
| sgCTRL | 1.032946 | 10 | 1.019731 | 10 |

The control condition is CTRL2 in source tables. Two I1 AMELY values remain excluded according to the retained analysis. Measurements span January–March 2025; their number does not establish independent derivations. AMELY 2x comprises L3 and L6; DDX3Y 2x comprises L6. No paired longitudinal effect, synergy, statistical superiority or Y-negative cell fraction is inferred. The numerical contrast supports a descriptive account of greater target depletion with pooled guides and lower values in the two-nucleofection samples.

## lpWGS values used

Canonical source reviewed: `../../../WGS/figures/A549_YposYneg/`, produced by `../../../WGS/code/figures_YposYneg.R` from the A549 QDNAseq table. Frozen copies of its four numerical summary tables are in [a549_lpwgs_source_snapshot](a549_lpwgs_source_snapshot/). The regional atlas PNG was visually inspected. The working Results link points to the existing atlas PDF; this is not yet a newly styled thesis export or an adopted final caption.

Raw chromosome-Y/autosomal ratios were independently recomputed from `bases / length` in each sample's `WGS/coverage/<sample>.genome.mosdepth.summary.txt`, using the sum of covered bases divided by the sum of chromosome lengths for chr1–chr22 as the autosomal denominator. They match `figure_raw_depth_ratios.tsv`:

| Sample | Mean autosomal depth | Mean Y depth | Y/autosomal ratio |
|---|---:|---:|---:|
| A549-CTRL2 | 0.460088561 | 0.150754653 | 0.327664424 |
| A549-LOY | 0.460878362 | 0.010147610 | 0.022017978 |

The control/edited ratio is approximately 14.9. The older `coverage_summary.tsv` uses rounded Y-depth entries and yields slightly different ratios; the draft uses the exact bases/length calculation. These raw mean-depth ratios are distinct from corrected regional QDNAseq copy-number ratios and from ddPCR estimated copies per cell.

The regional table records median corrected per-bin values of 0.348347 versus 0.007560 for male-specific Yp, and 0.329206 versus the exported 0.001 floor for euchromatic Yq. The Results describe this regional depletion qualitatively. The “mappable regions” aggregate in the plotting script includes PARs and should not be labelled male-specific. The separate MSY summary uses segmented values and a span that includes centromeric bins; it is not interchangeable with the per-region unsegmented summaries.

The producer uses 100-kb bins, QDNAseq correction/normalisation, and a `use` mask at plotting. The saved table producer applies `pmax(value, 1e-3)` before log transformation: **0.001 is a floor for values at or below that level, not proof of an exactly zero measurement**. The plotting threshold at log2 = -3 is an analysis/display convention, not an experimentally validated detection limit. Legacy captions claiming the male-specific region is entirely absent or attributing all non-Y differences to clonal drift were not adopted. No uniform complete-loss, breakpoint or causality claim is made.

## Context and outstanding work

A549's lung carcinoma origin is supported by [ATCC CCL-185](https://www.atcc.org/products/ccl-185), checked 17 September 2026. The culture/delivery selection rationale comes from the author's experimental account. Public Y-status calls must not replace measurements of the experimental stock.

Still required: recovery intervals; exact population-to-ddPCR-to-lpWGS-to-RNA/screen links; evaluation of regional non-Y differences; confirmation of reference FASTA/version and bin metadata for final Methods; and final thesis-style lpWGS panels and central captions. Existing April 2025 RT-qPCR plots have not been used to introduce new expression findings while their raw input remains unavailable.
