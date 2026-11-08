# A8 Pipeline Implementation Log

## Task 1 — Manifest and FCS loading

- RED: `.\codex\flow_env\python.exe -m unittest spectral.a8_pipeline.tests.test_manifest_fcs -v` failed with the expected missing-module error.
- GREEN: the same command completed with exit code 0; 3 tests passed in 0.178 seconds.
- Confirmed F11 contains 15,000 events, 253 parameters, hard-coded `eGFP-A`/`mCherry-A`, and matrix name `mcherry_GFP_matrix`.

## Task 2 — Acquisition segmentation and time QC

- RED 1: the time-QC test module failed with the expected missing-module error.
- GREEN 1: segmentation, Time reset, zero-MAD, and real-segment tests passed.
- RED 2: the new `time_qc` behavior test failed with the expected missing-function import error.
- GREEN 2: `.\codex\flow_env\python.exe -m unittest discover -s spectral/a8_pipeline/tests -v` completed with exit code 0; 8 tests passed.
- F11: one 15,000-event segment, median binned rate 192.69 events/s, no bins flagged.
- G11: two 10,000-event segments. Median binned rates were 397.23 and 383.43 events/s. The predefined rule flagged three 200-event Run 2 bins (first two bins after restart and bin 38); these remain subject to diagnostic-plot review and time-QC-off sensitivity analysis.

## Task 3 — Saturation decoding

- RED: the saturation test module failed with the expected missing-module error.
- GREEN: `.\codex\flow_env\python.exe -m unittest discover -s spectral/a8_pipeline/tests -v` completed with exit code 0; 12 tests passed.
- The decoded mask OR matched `Saturated` event-for-event in all five FCS files.
- Any/spectral/ImgB3 saturation percentages were F11 56.120/0.400/56.087 and G11 47.025/0.930/46.795.
- Control spectral saturation was 0.430% mCherry, 1.710% GFP, and 0.000% unstained.

## Task 4 — Scatter and pulse-singlet gates

- RED: the gating test module failed with the expected missing-module error.
- Initial GREEN attempt exposed a robustness failure: one high-leverage pulse outlier tilted the soft-loss slope to 1.243 instead of 1.0.
- Root cause: the ordinary initial fit inflated the robust-loss scale. A deterministic Theil–Sen fit on at most 2,000 ordered observations recovered the intended relationship; the targeted regression test then passed.
- The reviewed broad common scatter polygon retains 98.253% of F11, 96.890% of G11 Run 1, and 96.430% of G11 Run 2 before other exclusions. It excludes the sparse low-FSC/low-SSC tail without being called a viability gate.
- After time QC, spectral-saturation exclusion, and the cell polygon, the 3.5-MAD pulse gate retained 13,998/14,678 F11 events (95.37%) and 18,101/18,575 G11 events (97.45%).
- `.\codex\flow_env\python.exe -m unittest discover -s spectral/a8_pipeline/tests -v` completed with exit code 0; 16 tests passed.

## Task 5 — Reporter model, quadrant classification, and competition statistics

- RED: the reporter/statistics test module failed with the expected missing-module error.
- GREEN: `.\codex\flow_env\python.exe -m unittest discover -s spectral/a8_pipeline/tests -v` completed with exit code 0; 20 tests passed in 1.599 seconds.

### Primary-gated input to the pooled fit

Time QC, spectral-saturation exclusion, the common cell polygon, and the 3.5-MAD pulse
gate were applied in that order, reproducing the Task 4 counts exactly.

| Sample | Raw | After time/saturation/polygon | Primary singlets |
|---|---:|---:|---:|
| F11 | 15,000 | 14,678 | 13,998 |
| G11 | 20,000 | 18,575 | 18,101 |

Pooled primary events used for threshold estimation: 32,099.

### Fitted pooled two-Gaussian models

Transform `arcsinh(raw / 10000)`, seed 20260828, 8 initializations.

| Reporter | Mean 1 | Mean 2 | Sigma 1 | Sigma 2 | Weight 1 | Weight 2 | Separation | Log-likelihood |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `eGFP-A` | 1.415652 | 6.527939 | 0.917146 | 0.750221 | 0.473988 | 0.526012 | 6.101665 | -61377.7436 |
| `mCherry-A` | -2.907617 | 4.872710 | 0.935163 | 1.754918 | 0.276588 | 0.723412 | 5.533243 | -76573.7365 |

### Thresholds and validation outcome

| Reporter | Threshold (transformed) | Threshold (raw units) | Validation |
|---|---:|---:|---|
| `eGFP-A` | 4.186710 | 328,953.53 | Accepted |
| `mCherry-A` | -0.273371 | -2,767.88 | Accepted, flagged for review |

Both fits satisfy the three predefined criteria: every component weight is at or above
0.01, separation is at or above two pooled standard deviations, and a unique
equal-posterior root lies between the component means.

The negative raw `mCherry-A` threshold is expected rather than erroneous. These are
unmixed outputs, 28.08% of pooled `mCherry-A` values are at or below zero, and the
mCherry-negative mode sits near `arcsinh` -3.2, which is about -122,000 raw units. Both
thresholds also fall in the observed density minimum of their pooled histograms, so
neither is an artifact of the model disagreeing with the visible distribution.

### Review flag: intermediate mCherry shoulder and the double-positive population

The pooled `mCherry-A` distribution is not cleanly two-component. A broad shoulder spans
roughly `arcsinh` 1.6 to 4.9 between the negative mode and the main positive mode. The
two-component model absorbs that shoulder into a wide positive component (sigma 1.754918),
so the accepted threshold classifies the shoulder as mCherry-positive.

Classification with the accepted thresholds:

| Sample | GFP-only | mCherry-only | Double-positive | Double-negative |
|---|---:|---:|---:|---:|
| F11 | 3,183 | 6,836 | 3,685 | 294 |
| G11 Run 1 | 2,510 | 3,894 | 2,668 | 272 |
| G11 Run 2 | 2,368 | 3,642 | 2,471 | 276 |
| G11 combined | 4,878 | 7,536 | 5,139 | 548 |

The double-positive population has the same GFP intensity as the GFP-only population but
a markedly lower mCherry intensity than the mCherry-only population, in both dye
orientations:

| Sample | Population | Median `arcsinh` GFP | Median `arcsinh` mCherry |
|---|---|---:|---:|
| F11 | GFP-only | 6.51 | -3.06 |
| F11 | Double-positive | 6.52 | 3.08 |
| F11 | mCherry-only | 1.30 | 5.73 |
| G11 | GFP-only | 6.59 | -3.08 |
| G11 | Double-positive | 6.52 | 3.55 |
| G11 | mCherry-only | 1.26 | 6.05 |

92.4% of F11 and 87.1% of G11 double-positive events lie below `arcsinh` 4.86 on
mCherry, inside the shoulder. Because the effect tracks GFP intensity and appears at the
same GFP level in both orientations, it is more consistent with GFP spillover or an
unmixing residual on the mCherry axis than with genuine co-expression.

This is recorded as an unresolved ambiguity, not corrected here. The predefined model was
not modified, and no threshold was chosen by reference to the resulting LOY:CTRL values.
Task 6 image review and the Task 7 threshold sensitivity family are the appropriate
places to resolve whether the shoulder belongs with the negative or the positive
population. Until then the double-positive counts above should be treated as a
classification limitation, which matters because double-positive events are excluded from
the assignable LOY/CTRL denominator and are 26.3% of F11 and 28.4% of G11 primary events.

### Cross-check limitation

The three single-stain control FCS files do not contain `eGFP-A`, `mCherry-A`, or
`Autofluorescence_001-A`; only the two mixed recordings carry the hard-coded
`mcherry_GFP_matrix` outputs. The controls therefore cannot anchor the reporter
thresholds in unmixed space and can only be used on raw detector parameters.

## Task 6 — TIFF mapping and review montages

- RED: `.\codex\flow_env\python.exe -m unittest spectral.a8_pipeline.tests.test_images -v` failed with the expected missing-module error.
- GREEN: `.\codex\flow_env\python.exe -m unittest discover -s spectral/a8_pipeline/tests -v` completed with exit code 0; 28 tests passed in 2.177 seconds.
- Verified image coverage: F11 8,776/15,000 and G11 20,000/20,000, matching the design document.

### Deviation 1: the filename suffix is zero-based, not one-based

The spec and plan both state that TIFF filenames encode one-based FCS event
numbers. They do not. G11 contains exactly 20,000 files numbered 00000000
through 00019999 with no gaps, so a one-based reading is arithmetically
impossible, and the file `..._00000000.tiff` would map to event 0. F11 ranges
from 00000007 to 00014999, which is also consistent with a zero-based index.

Had the plan been followed literally, `build_image_manifest` would have rejected
G11 outright as containing an out-of-range event, and every F11 montage would
have been silently misattributed by one event.

`parse_event_number` returns the literal suffix, which keeps the behaviour the
plan's test specified. `build_image_manifest` converts it explicitly, so
`event_number` is one-based everywhere else in the pipeline and `image_index` is
the zero-based suffix.

### Deviation 2: `ImageSequence.Iterator` returns stale frames

The plan specifies `PIL.Image.open` plus `ImageSequence.Iterator`. For these
multi-frame float32 TIFFs the iterator yields identical pixel data for all six
frames: every frame reported min -0.008268, max 0.009471, mean 0.000000. Reading
with an explicit `seek` followed by `load` returns six genuinely different
frames. `read_tiff_channels` does the latter and asserts a common frame shape.
This is very likely the defect behind the existing `f11_cluster_image_preview`
and `f11_cluster_image_preview_corrected` pair.

### Frame-to-channel mapping

The TIFFs contain no channel metadata, so the mapping was derived from the data.
Each frame's peak was correlated against the per-channel `Max Intensity`
parameters in the FCS, which are populated for 100% of events.

| Frame | Imaging channel | Spearman rho |
|---:|---|---:|
| 1 | `FSC` | 0.939 |
| 2 | `SSC (Imaging)` | 0.967 |
| 3 | `GFP_reporter eGFP` | 0.944 |
| 5 | `mCherry_reporter mCherry` | 0.950 |

Frame 0 is the light-loss image, where cells appear dark against a bright field;
its peak statistic is background-dominated so the correlation is weak and it is
not asserted. Frame 4 responds to both reporters and is not confidently
assigned. `DEFAULT_CHANNEL_MAP` therefore uses frame 1 as grayscale structure,
frame 3 as green, and frame 5 as red.

Frame 5 correlates with `eGFP-A` at -0.058, so the mCherry imaging channel is
effectively independent of GFP and can serve as a second opinion on the spectral
mCherry axis.

### Deviation 3: montage display normalization

The plan specifies fixed 1st/99th percentile normalization. More than half of
every frame is background, so a 1st-percentile black point renders the noise
floor as a visible haze and makes the montages unreadable. The default is now
the 50th/99.5th percentile, computed once across all tiles in a montage and
recorded in the montage header. This is a display choice only and changes no
quantitative FCS measurement.

### Review of the Task 5 double-positive flag

Task 5 recorded that the double-positive population has the same GFP intensity
as GFP-only cells with intermediate mCherry, and flagged spillover as the
leading explanation. Image review revises that.

The mCherry imaging channel is independent of GFP, so it is an orthogonal test.
Comparing double-positives against GFP-only cells matched on `FSC-A` within
deciles, double-positives carry consistently more imaging mCherry:

| Sample | Size-matched imaging-mCherry ratio, double-positive / GFP-only | Range across deciles |
|---|---:|---|
| F11 | 1.38x | 1.35 - 1.45 |
| G11 | 1.88x | 1.64 - 2.13 |

The effect is present in all ten deciles of both samples, so it is not a size
artifact and not consistent with pure spillover, which would leave an
independent channel at background.

The montages show what the numbers cannot. Several double-positive events are
plainly two cells in one image: G11 events 263 and 1662 each show a separate
green cell and a separate red cell within a single event, and event 816 shows
the same. The pulse-singlet gate does remove many doublets, and the pulse
non-singlet row is full of obvious ones such as events 948 and 1688, but it does
not remove all of them.

`Delta CoM (GFP_reporter eGFP/mCherry_reporter mCherry)` is not elevated for
double-positives, at a median of 0.19 against 0.22 and 0.27 for the single
classes in G11. That parameter is only computed for 10.8% of F11 and 16.8% of
G11 primary events, so it is underpowered and does not overturn the direct
visual evidence.

The double-positive population is therefore a mixture of unresolved GFP/mCherry
doublets and single cells with genuine low mCherry, not a single artifact with
one cause. It remains excluded from the assignable LOY/CTRL denominator and
remains a limitation, at 26.3% of F11 and 28.4% of G11 primary events. Resolving
the split is a Task 7 sensitivity question, specifically the singlet-limit family
at 3.0, 3.5, and 4.0 scaled MAD.

### F11 image coverage is not missing at random

F11 coverage is 58.51% overall, but the missing images are structured:

- Coverage among scatter-kept events is 59.55%; among scatter-excluded events it is 0.00%, for all 262 such events.
- Events with an image have median `FSC-A` 75,924,880; events without have 105,022,420.
- Coverage by acquisition decile ranges from 45% to 64%.

Missing F11 images are still never used as an exclusion, as required. But F11
montages are drawn from a biased subset that under-represents large events and
contains no scatter-excluded events at all, so no F11 image montage supports a
quantitative claim about population frequencies. G11 has complete coverage and
carries any image-based argument.

### Artifacts

Written under `spectral/a8_pipeline/montages/`, repository-local because the
external output root still requires the Task 8 filesystem approval:

- `montage_F11_full.png` and `montage_G11_full.png`, up to 50 images per population.
- `montage_F11_review.png` and `montage_G11_review.png`, 8 per population for direct inspection.
- `image_montage_manifest.csv`, 955 rows recording sample, population, event number, image filename, raw `eGFP-A` and `mCherry-A`, both transformed threshold distances, and the any/spectral/ImgB3 saturation flags.
- `frame_identity_sheet.png`, the six-frame contact sheet supporting the channel mapping.

Populations covered: scatter-kept, scatter-excluded, pulse singlets, pulse
non-singlets, the four reporter quadrants, the 50 events nearest each reporter
threshold, and spectral-saturation events. Selection uses
`np.random.default_rng(20260828)`. Populations with fewer than 50 available
images are reported at their true count rather than padded; F11 scatter-excluded
yields 0 images and F11 GFP boundary yields 11.

## Task 7 — Sensitivity engine, SVG diagnostics, and reports

- RED: `.\codex\flow_env\python.exe -m unittest spectral.a8_pipeline.tests.test_sensitivity_report -v` failed with the expected missing-module error.
- GREEN: `.\codex\flow_env\python.exe -m unittest discover -s spectral/a8_pipeline/tests` completed with exit code 0; 38 tests passed in 2.753 seconds.
- Created `sensitivity.py`, `report.py`, `templates/report.html.j2`; extended `svg.py` with five figure writers.
- A preview run over the real export wrote 16 artifacts to `spectral/a8_pipeline/preview_report/`: Markdown, HTML, three CSVs, and 11 SVG figures. None are empty and every SVG parses as XML.

### Sensitivity matrix

22 scenarios across the seven required families, with exactly one baseline. The
threshold family carries all nine combinations of GFP and mCherry shifts at plus
or minus 10% of the transformed distance between component means, including the
unshifted centre as an internal check. The G11 Run 2 scenario is skipped for F11
because F11 has a single segment.

### Baseline competition result

| Sample | LOY reporter | GFP-only | mCherry-only | Double-pos | Double-neg | Assignable | LOY fraction | LOY:CTRL | log2 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| F11 | GFP | 3,183 | 6,836 | 3,685 | 294 | 10,019 | 0.3177 | 0.4656 | -1.1028 |
| G11 | mCherry | 4,878 | 7,536 | 5,139 | 548 | 12,414 | 0.6071 | 1.5449 | +0.6275 |

G11 segments agree closely with each other: Run 1 gives a LOY fraction of 0.6081
and Run 2 gives 0.6060, so the acquisition restart is not a material problem.

### The reciprocal dye swap does not reproduce

This is the most important result so far and it is a finding, not a defect in the
pipeline.

F11 and G11 disagree in direction: log2 ratio -1.1028 against +0.6275, a gap of
1.7303. The two recordings are the same competition read with the reporters
swapped, so they should agree.

They disagree because the measurement tracks the reporter rather than the
genotype. In both orientations the mCherry-labelled arm is the larger one:
mCherry-only over GFP-only is 2.148 in F11 and 1.545 in G11, regardless of which
genotype carries mCherry. Any per-sample LOY:CTRL from a single orientation is
therefore reporting reporter bias as much as biology.

### Where the disagreement comes from

The double-positive population identified in Tasks 5 and 6 is the dominant
contributor. Those events have GFP intensity identical to GFP-only cells, so
classifying them as double-positive removes GFP-side cells from the assignable
denominator in both orientations and depresses the GFP arm twice.

Reclassifying them, as a diagnostic only:

| Treatment of double-positives | F11 log2 | G11 log2 | Disagreement | Dye-swap mean log2 |
|---|---:|---:|---:|---:|
| Excluded (primary) | -1.103 | +0.628 | 1.730 | -0.238 |
| Counted as GFP | +0.007 | -0.411 | **0.417** | -0.202 |
| Counted as mCherry | -1.725 | +1.378 | 3.102 | -0.174 |
| Split evenly | -0.788 | +0.440 | 1.228 | -0.174 |

Counting double-positives as GFP reduces the disagreement roughly fourfold. That
is consistent with the Task 6 evidence that these are GFP cells carrying either
an unresolved mCherry partner or genuine low mCherry, but it is recorded here as
a diagnostic and is **not adopted**. Choosing it because it improves concordance
would violate the constraint against selecting gates by their effect on the
result. The primary analysis still excludes double-positives.

The dye-swap mean is far more stable than either arm: -0.238, -0.202, -0.174 and
-0.174 across all four treatments, a spread of 0.064 against 2.7 for the
individual arms. The exploratory dye-swap summary is therefore about log2 -0.2,
or a LOY:CTRL near 0.87, and that summary is the only competition number here
that survives the double-positive ambiguity. It still rests on one acquisition
per orientation and remains technical rather than biological replication.

### Sensitivity families

Maximum absolute change against the baseline:

| Family | F11 LOY fraction | F11 log2 | G11 LOY fraction | G11 log2 |
|---|---:|---:|---:|---:|
| time | 0.0000 | 0.0000 | 0.0011 | 0.0069 |
| segment | 0.0000 | 0.0000 | 0.0010 | 0.0061 |
| cell | 0.0037 | 0.0248 | 0.0009 | 0.0056 |
| singlet | 0.0014 | 0.0095 | 0.0040 | 0.0243 |
| **saturation** | **0.1767** | **1.5039** | **0.0982** | **0.5766** |
| threshold | 0.0156 | 0.1054 | 0.0155 | 0.0944 |
| cofactor | 0.0010 | 0.0069 | 0.0008 | 0.0051 |

Every gating family is small. Time QC, segment choice, polygon scale, singlet
limit and transform cofactor all move the LOY fraction by less than 0.005, and
the threshold family by about 0.016.

The saturation family dominates by an order of magnitude, exactly as the Task 3
analysis predicted. The `any` rule requires `Saturated == 0` and so discards
roughly half of each recording on the strength of `ImgB3 (790)` imaging
clipping. It moves the F11 log2 ratio by 1.50, which is comparable to the entire
dye-swap disagreement. This confirms that a blanket saturation exclusion is not a
defensible primary rule for a spectral fluorescence ratio, and the primary
spectral-only rule is retained.

### Note on the sensitivity summary

`summarise_families` takes the baseline LOY fraction and log2 ratio as explicit
arguments rather than reading a baseline row, so the same helper works for
per-sample subsets where the baseline row may have been filtered out. Deltas are
computed over every row in a family, including the baseline row itself, so a
family that contains the baseline reports a floor of zero rather than omitting it.

## Task 8 — CLI orchestration, integration dry run, full analysis, and verification

- RED: `.\codex\flow_env\python.exe -m unittest spectral.a8_pipeline.tests.test_integration -v` failed with the expected missing-module error.
- GREEN: `.\codex\flow_env\python.exe -m unittest discover -s spectral/a8_pipeline/tests -v` completed with exit code 0; 49 tests passed in 211.030 seconds.
- Created `run_analysis.py`, `tests/test_integration.py`, and `README.md`.

### Integration test structure

The integration tests share one dry run through `setUpClass`. A single run hashes
9.55 GB across 28,784 source files twice, at about 48 seconds per pass, so
running the pipeline once per test method would have dominated the suite. The
11 integration tests complete in 214 seconds.

### Step 5: repository-local CLI dry run

```
.\codex\flow_env\python.exe -m spectral.a8_pipeline.run_analysis
    --input 'A:\flow_jo_AUG28'
    --output '.\spectral\a8_pipeline\dry_run_output'
    --config '.\spectral\a8_pipeline\default_config.json' --dry-run
```

Exit code 0. Validation, segmentation, saturation agreement, gating and model
fitting all succeeded; staging was discarded and no output directory was created.

### Defect found and fixed during Step 5

The first dry run completed the entire analysis and hash re-verification, then
failed on the final `shutil.rmtree` of the staging directory with
`PermissionError: [WinError 32]`. The repository lives inside a Dropbox folder,
and the sync client had opened the just-written directory. Windows raises rather
than blocking.

`_retry_filesystem` now wraps both the staging removal and the staging-to-output
rename with eight attempts and a linear backoff, and raises a message naming a
sync client or virus scanner as the likely cause if it still fails. The rename
needed the same protection as the removal: an external run whose rename failed
would leave a complete analysis stranded in a hidden staging directory.

### Step 6: diagnostics review and configuration lock

Reviewed the scatter density and polygon, both pulse relationships, the pooled
reporter densities with their thresholds, the binned event-rate plots, the
saturation composition, and both image montages.

No value in `default_config.json` was changed. The polygon, singlet limit, time-QC
parameters, cofactor, seed, and validity criteria all survived review. The
reporter thresholds are fitted at runtime and frozen into the output's
`analysis_config.json` rather than hard-coded into the input configuration, so
the run stays reproducible without pre-committing a fitted value.

The one deliberate non-configuration change was the montage display
normalization, moved from the 1st/99th percentile to the 50th/99.5th in Task 6
because more than half of every frame is background. That is a display choice and
changes no measurement.

The mCherry threshold was **not** adjusted despite the double-positive ambiguity.
Adjusting it would have meant choosing a gate by its effect on LOY:CTRL. The
predefined threshold sensitivity family already brackets it at plus and minus 10%
of the component separation.

### Step 7: external analysis

```
.\codex\flow_env\python.exe -m spectral.a8_pipeline.run_analysis
    --input 'A:\flow_jo_AUG28'
    --output 'A:\flow_jo_AUG28\programmatic_analysis'
    --config '.\spectral\a8_pipeline\default_config.json'
```

Exit code 0. Artifacts were staged in `A:\flow_jo_AUG28\.programmatic_analysis_staging`,
source hashes were re-verified as unchanged, and the staging directory was then
renamed to `A:\flow_jo_AUG28\programmatic_analysis`. 30 files, 26.7 MB.

Both runs reproduced the same numbers: F11 13,998 and G11 18,101 primary
singlets, `eGFP-A` threshold 328,953.53 raw, `mCherry-A` threshold -2,767.88 raw,
dye-swap mean log2 ratio -0.2376.

### Step 8: verification

27 of 27 independent checks passed, recomputed from the sources rather than read
back from the pipeline's own summaries:

- All 13 required files present and nonempty; 11 SVG figures and 4 montage PNGs, none empty.
- Event labels unique by sample and event number; counts 15,000 and 20,000; one-based and contiguous.
- Gate counts nested from raw through time QC, saturation, polygon, and singlets, reconciling to source event counts.
- F11 and G11 both expose `eGFP-A`, `mCherry-A` and `Autofluorescence_001-A`, both with matrix `mcherry_GFP_matrix`.
- G11 segments 10,000/10,000 before QC; F11 a single 15,000-event segment.
- Decoded saturation OR matches the `Saturated` flag for every event in both recordings.
- Image coverage 8,776/15,000 and 20,000/20,000.
- Input manifest covers all 28,784 sources; no source missing or resized; all 8 FCS, CVW and FlowJo files re-hashed and matching.
- Frozen thresholds and fitted mixture parameters recorded in `analysis_config.json`.
- Report states counting precision, the viability limitation, and the dye-swap failure.

### Unresolved limitations carried into the report

1. The reciprocal dye swap does not reproduce: F11 log2 -1.1028 against G11 +0.6275. In both orientations the mCherry-labelled arm is larger regardless of genotype. Only the dye-swap mean, about -0.2, is stable across double-positive treatments.
2. Double-positive events are 26.3% of F11 and 28.4% of G11 primary events and are excluded from the assignable denominator. Image review shows a mixture of unresolved GFP/mCherry doublets and single cells with genuine low mCherry.
3. The saturation family dominates every other sensitivity family. A blanket `Saturated == 0` rule discards roughly half of each recording on `ImgB3 (790)` imaging clipping.
4. `ImgB3` saturation is genotype-associated, not reporter-associated: LOY 78.6% against CTRL 41.4% in F11, and LOY 59.9% against CTRL 39.4% in G11. Blanket saturation exclusion therefore preferentially removes LOY cells, which is the specific hazard the saturation analysis warned about and the reason that family moves the result so much.
5. F11 image coverage is 8,776/15,000 and not missing at random: no scatter-excluded event has an image, and events without images are larger.
6. The single-stain controls lack the unmixed reporter parameters and cannot anchor thresholds.
7. One acquisition per dye orientation and no culture replicates. This is a technical and QC analysis, not evidence of biological reproducibility.

### Final state

All eight tasks are complete. 49 tests pass. The analysis is at
`A:\flow_jo_AUG28\programmatic_analysis` and no source FCS, CVW, TIFF or FlowJo
file was modified at any point.
