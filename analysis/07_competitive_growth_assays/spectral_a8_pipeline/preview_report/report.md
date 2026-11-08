# A8 spectral competition analysis (preview)

Preview run of the Task 7 sensitivity engine on the August 28 export. Artifacts are repository-local; the external output root awaits Task 8 approval.

## Interpretation of intervals

Wilson intervals describe counting precision for the events acquired. They are not biological replicate uncertainty: there is one acquisition per dye orientation.

## Tables

### Competition results

| sample | scope | loy_reporter | loy | ctrl | assignable | loy_fraction | ratio | log2_ratio | wilson_low | wilson_high |
|---|---|---|---|---|---|---|---|---|---|---|
| F11 | all | GFP | 3183 | 6836 | 10019 | 0.317696 | 0.465623 | -1.10277 | 0.308651 | 0.326881 |
| F11 | Run 1 | GFP | 3183 | 6836 | 10019 | 0.317696 | 0.465623 | -1.10277 | 0.308651 | 0.326881 |
| G11 | all | mCherry | 7536 | 4878 | 12414 | 0.607057 | 1.5449 | 0.627509 | 0.598433 | 0.615614 |
| G11 | Run 1 | mCherry | 3894 | 2510 | 6404 | 0.608057 | 1.55139 | 0.633566 | 0.596039 | 0.619946 |
| G11 | Run 2 | mCherry | 3642 | 2368 | 6010 | 0.60599 | 1.53801 | 0.621062 | 0.593572 | 0.618272 |

### Sensitivity family summary

| sample | family | scenarios | max_abs_delta_loy_fraction | max_abs_delta_log2_ratio |
|---|---|---|---|---|
| F11 | time | 2 | 0 | 0 |
| F11 | segment | 1 | 0 | 0 |
| F11 | cell | 2 | 0.00374306 | 0.0248346 |
| F11 | singlet | 2 | 0.00142639 | 0.00948211 |
| F11 | saturation | 3 | 0.176673 | 1.50392 |
| F11 | threshold | 9 | 0.0156246 | 0.105422 |
| F11 | cofactor | 2 | 0.0010411 | 0.00692304 |
| G11 | time | 2 | 0.00114522 | 0.0069228 |
| G11 | segment | 2 | 0.00100098 | 0.00605129 |
| G11 | cell | 2 | 0.000927879 | 0.00560953 |
| G11 | singlet | 2 | 0.0040245 | 0.0242971 |
| G11 | saturation | 3 | 0.0982379 | 0.576614 |
| G11 | threshold | 9 | 0.0154972 | 0.0944179 |
| G11 | cofactor | 2 | 0.000844367 | 0.00510484 |

### Saturation composition (%)

| sample | any_saturated | spectral_saturated | imgb3_saturated |
|---|---|---|---|
| F11 | 56.12 | 0.4 | 56.0867 |
| G11 | 47.025 | 0.93 | 46.795 |

## Figures

![pulse fsc F11](figures/pulse_fsc_F11.svg)

*pulse fsc F11* — `figures/pulse_fsc_F11.svg`

![pulse fsc G11](figures/pulse_fsc_G11.svg)

*pulse fsc G11* — `figures/pulse_fsc_G11.svg`

![reporter eGFP A](figures/reporter_eGFP_A.svg)

*reporter eGFP A* — `figures/reporter_eGFP_A.svg`

![reporter mCherry A](figures/reporter_mCherry_A.svg)

*reporter mCherry A* — `figures/reporter_mCherry_A.svg`

![saturation composition](figures/saturation_composition.svg)

*saturation composition* — `figures/saturation_composition.svg`

![scatter F11](figures/scatter_F11.svg)

*scatter F11* — `figures/scatter_F11.svg`

![scatter G11](figures/scatter_G11.svg)

*scatter G11* — `figures/scatter_G11.svg`

![sensitivity delta F11](figures/sensitivity_delta_F11.svg)

*sensitivity delta F11* — `figures/sensitivity_delta_F11.svg`

![sensitivity delta G11](figures/sensitivity_delta_G11.svg)

*sensitivity delta G11* — `figures/sensitivity_delta_G11.svg`

![time qc F11](figures/time_qc_F11.svg)

*time qc F11* — `figures/time_qc_F11.svg`

![time qc G11](figures/time_qc_G11.svg)

*time qc G11* — `figures/time_qc_G11.svg`

## Limitations

- There was no viability dye, so scatter-excluded events are not dead cells.
- Double-positive events are 26-28% of primary events and are excluded from the assignable denominator.
- F11 image coverage is 58.51% and is not missing at random.
- Single-stain controls lack the unmixed reporter parameters and cannot anchor thresholds.
