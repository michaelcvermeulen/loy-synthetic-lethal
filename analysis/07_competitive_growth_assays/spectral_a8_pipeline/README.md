# A8 spectral competition analysis pipeline

Read-only analysis of the BD FACSDiscover A8 reciprocal dye-swap pair in
`A:\flow_jo_AUG28`. Nothing here modifies a source FCS, CVW, TIFF or FlowJo file.

## Running it

```powershell
.\codex\flow_env\python.exe -m spectral.a8_pipeline.run_analysis `
    --input 'A:\flow_jo_AUG28' `
    --output 'A:\flow_jo_AUG28\programmatic_analysis' `
    --config '.\spectral\a8_pipeline\default_config.json'
```

Add `--dry-run` to validate and compute everything without creating the output
directory. Run from the repository root so that `spectral` is importable.

Artifacts are generated into a sibling `.programmatic_analysis_staging`
directory and renamed into place only after source hashes are re-verified. An
existing output directory is never deleted or overwritten; the run exits with an
error instead.

## Tests

```powershell
.\codex\flow_env\python.exe -m unittest discover -s spectral/a8_pipeline/tests -v
```

The integration tests run the pipeline once for the whole class. A full run
hashes 9.55 GB across 28,784 source files twice, so per-method runs would
dominate the suite.

## Modules

| Module | Responsibility |
|---|---|
| `fcs_io.py` | Read FCS events and metadata without modifying the source. |
| `manifest.py` | Inventory sources, hash them, and verify they are unchanged. |
| `time_qc.py` | Detect acquisition segments and flag anomalous time bins. |
| `saturation.py` | Reinterpret the saturation bitmasks into named detector flags. |
| `gating.py` | Scatter polygon and robust pulse-singlet gates. |
| `reporters.py` | Two-component mixture fitting and quadrant classification. |
| `statistics.py` | Competition metrics and Wilson intervals. |
| `images.py` | Map TIFFs to events and render review montages. |
| `sensitivity.py` | The fixed sensitivity matrix and its evaluation. |
| `svg.py` | Dependency-free SVG figures. |
| `report.py` | Markdown and HTML reporting. |
| `run_analysis.py` | CLI orchestration, staging, and verification. |

## Things that are easy to get wrong

**Image filename suffixes are zero-based.** `..._00000007.tiff` is event index 7,
which is one-based event number 8. G11 has exactly 20,000 files numbered
`00000000`-`00019999`, so a one-based reading is impossible. Every table the
pipeline writes uses one-based `event_number`.

**Multi-frame TIFFs need an explicit seek.** `PIL.ImageSequence.Iterator` returns
stale pixel data for these float32 stacks; all six frames come back identical.
`read_tiff_channels` seeks and loads each frame.

**TIFF frames carry no channel metadata.** The mapping was established by
correlating each frame against the fully populated per-channel `Max Intensity`
parameters: frame 1 is FSC, frame 2 is SSC, frame 3 is GFP, frame 5 is mCherry.
Frames 0 and 4 are not confidently assigned.

**The reporter parameters are already unmixed.** `eGFP-A`, `mCherry-A` and
`Autofluorescence_001-A` come from the FACSChorus matrix `mcherry_GFP_matrix`.
Never apply conventional compensation to them. The single-stain controls do not
contain these parameters and cannot anchor the thresholds.

**Saturation is dominated by an imaging detector.** Roughly half of each mixed
recording has `Saturated == 1`, almost entirely from `ImgB3 (790)`. Spectral
detector saturation is under 1%. Excluding every saturated event is available as
a sensitivity scenario but is not a defensible primary rule.

**Scatter-excluded is not dead.** No viability dye was acquired.

## Known limitations of the current result

- The reciprocal dye swap does not reproduce. F11 and G11 disagree in direction,
  and in both orientations the mCherry-labelled arm is the larger one regardless
  of genotype. Only the dye-swap mean is stable.
- Double-positive events are 26-28% of primary events and are excluded from the
  assignable denominator. Image review shows they are a mixture of unresolved
  GFP/mCherry doublets and single cells with genuine low mCherry.
- F11 image coverage is 8,776/15,000 and is not missing at random.
- One acquisition per dye orientation, no culture replicates. This is a
  technical and QC analysis.
