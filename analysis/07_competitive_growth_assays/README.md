# `07_competitive_growth_assays/` — Reciprocal Two-Colour Spectral Flow-Cytometry Competition Validation

**Thesis Reference:** Chapter 4, §4.2.5, §4.3.2–4.3.4, and Supplementary Methods S4.5 (`Fig 4.5–4.6`, `Supp Fig 4.12`).

This module contains the analytical pipeline (`spectral_a8_pipeline/` and top-level calibration scripts) for the reciprocal two-colour (`eGFP` / `mCherry`) competitive growth assay quantified by spectral flow cytometry (`BD FACSDiscover A8` and `CytoFLEX S`), designed to validate candidate LOY synthetic lethalities in direct co-culture without single-cell cloning artefacts or metabolic viability dye confounders.

---

## Assay Architecture & Rationale

1. **Reciprocal Dye-Swap Co-Culture Design (`Fig 4.5a`):**
   - Stable lentiviral `eGFP` and `mCherry` reporter derivatives were established for both `HT1080 sgCTRL` and `HT1080 LOY` (`CP + Q53 2x`) populations.
   - Every Cas9 RNP perturbation (`NTC`, `OR2J2` cutting control, `DDX3X`, `TADA2B`, `TAF6L`, `SREBF1`) is nucleofected into individual reporter-labelled populations, which are then mixed 1:1 in **both reciprocal colour orientations**:
     - **Mix A (Orientation 1):** `LOY-eGFP + CTRL-mCherry`
     - **Mix B (Orientation 2):** `LOY-mCherry + CTRL-eGFP`
     - Plus **same-genotype controls** (`CTRL-eGFP + CTRL-mCherry` and `LOY-eGFP + LOY-mCherry`) to isolate fluorophore-intrinsic fitness effects.
2. **Spectral Unmixing & Autofluorescence Exclusion (`Fig 4.5b–d`, `SATURATION_ANALYSIS.md`):**
   - Unmixing on the `BD FACSDiscover A8` explicitly models cellular autofluorescence (`AF`) as a separate spectral endmember and flags detector-saturated events, eliminating false double-negative/double-positive spillover across passages.
3. **Log-Odds Fitness Scoring (`spectral_a8_pipeline/statistics.py`):**
   - For each well and timepoint $t$, the `LOY`-to-`CTRL` ratio is transformed to a baseline-referenced log$_2$ odds shift:
     $$\Delta\log_2\!\left(\frac{\mathrm{LOY}}{\mathrm{CTRL}}\right)_t = \log_2\!\left(\frac{p_{\mathrm{LOY},t}}{p_{\mathrm{CTRL},t}}\right) - \log_2\!\left(\frac{p_{\mathrm{LOY},\mathrm{D0}}}{p_{\mathrm{CTRL},\mathrm{D0}}}\right)$$
   - Subtracting the matched unperturbed (`NTC` / `OR2J2`) co-culture trajectory yields the net gene-specific `LOY` fitness penalty ($\Delta\Delta\log_2\text{ odds}$).

---

## Directory Structure

```text
07_competitive_growth_assays/
├── README.md
├── comp_assay_plan.md                       # Experimental design, guide sequences & plate layouts
├── SATURATION_ANALYSIS.md                   # BD FACSDiscover A8 detector saturation & unmixing audit
├── comp_assay.py                            # Core gating, mixing calibration & log-odds pipeline
├── competition_pipeline.py                  # Automated multi-run competition processor
├── autofluorescence_analysis.py             # Spectral autofluorescence detector selection & spillover check
├── flowkit_gating.py                        # Programmatic FlowKit/GatingML polygon & logicle gating
├── qc2_20260821_pipeline.py                 # 9-point reciprocal mixing calibration & stain-index pipeline
├── mix_5050_ratios.csv / qc_master_table.csv # Calibration & baseline QC summary tables
└── spectral_a8_pipeline/                    # Modular BD FACSDiscover A8 package (D0, D4, longitudinal & QC reports)
```

---

## Headline Findings

1. **9-Point Reciprocal Mixing Calibration & Baseline Characterization (`Fig 4.5e–h`, `Supp Fig 4.12`):**
   - Across 9-point volumetric titration series (10:90 to 90:10) in both `LOY-eGFP : CTRL-mCherry` and `LOY-mCherry : CTRL-eGFP` orientations, spectral gating recovers expected mixture fractions with **$R^2 > 0.996$** and near-unity slopes (`Fig 4.5e,f`).
   - Multiplex `ddPCR` confirms preserved Y depletion in the reporter-transduced `LOY` lines (`Fig 4.5g`), and genomic PCR confirms on-target editing across the RNP-nucleofected competition populations (`Fig 4.5h`).
2. **Baseline Proliferative Advantage of Unperturbed LOY Cells & Selective Fitness Penalties After *DDX3X* and *TADA2B* Knockout (`Fig 4.6`, `Supp Fig 4.12`):**
   - In unperturbed (`NTC` and `OR2J2`) co-cultures from `D0` to `D4`, Y-depleted HT1080 cells increase in relative abundance against `sgCTRL` in both colour orientations (`Fig 4.6a`), revealing a baseline co-culture fitness advantage of `LOY` cells.
   - Relative to the unperturbed baseline trajectory, ***DDX3X*** knockout and ***TADA2B*** (SAGA subunit) knockout induce consistent **negative $\Delta\Delta\log_2(\mathrm{LOY}/\mathrm{CTRL})$ fitness shifts across both reciprocal colour orientations (`Mix A` and `Mix B`)** from `D0` to `D4` (`Fig 4.6b,c`), whereas *SREBF1* shows a smaller, orientation-dependent shift over this early window.
