# `07_competitive_growth_assays/`

Reciprocal two-colour (`eGFP` / `mCherry`) competitive growth assay pipeline quantified on the `BD FACSDiscover A8` spectral flow cytometer (`Chapter 4, Section 4.3.2–4.3.4`; `Fig 4.5–4.6`, `Supp Fig 4.12`).

- [`spectral_a8_pipeline/`](spectral_a8_pipeline/) — Modular Python package for spectral unmixing QC, gating, 9-point mixing calibration, and log-odds fitness scoring.
- [`comp_assay.py`](comp_assay.py), [`competition_pipeline.py`](competition_pipeline.py), [`flowkit_gating.py`](flowkit_gating.py), [`qc2_20260821_pipeline.py`](qc2_20260821_pipeline.py) — Top-level gating and competition analysis scripts.
- [`outputs/`](outputs/) — Calibration, gating QC, and competition ratio summary tables.
