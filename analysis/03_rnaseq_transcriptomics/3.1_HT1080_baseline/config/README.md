# Core experimental design record

Confirmed by the user on 11 September 2026:

- Include HT1080 from both August 2025 and May 2026.
- Sample suffixes 1, 2 and 3 are successive passage numbers. LOY and sgCTRL with the same passage number were matched by passage/day within each period.
- Passages come from the same starting treatment populations.
- August and May sampled the same established populations again; they trace to the same original nucleofected pool, not independently generated LOY experiments.

The core has 12 RNA libraries: one established CTRL population and one established LOY population, each measured at three passages in each of two sampling periods. There is one population-generation history. Passage measurements are repeated biological observations, not independent population-generation experiments or repeated sequencing lanes from a single library.

`sample_manifest_draft.csv` preserves the raw matrix's 15-column order, including three August parental samples for context. The core rows have confirmed population history. Parental relationships, actual harvest dates and optional addback intervention units remain unconfirmed; their missing metadata do not block the core descriptive analysis.

Identifiers are analysis labels, not claimed laboratory barcodes:

| Field | Meaning |
|---|---|
| `sample_key` | Unique RNA-library key matching the count matrix |
| `batch` | Original batch label, also used as the sampling-period label |
| `passage_number` | Recorded passage order within period; not a known elapsed time |
| `passage_day_block` | CTRL/LOY timing match within period; not an independent replicate |
| `within_batch_population_id` | Population trajectory within a period, retained for plotting |
| `starting_population_id` | Established treatment population, shared across periods |
| `generation_experiment_id` | Common nucleofection origin across the core samples |
| `independent_unit_id` | Unit for an independently generated LOY comparison; the core has one shared origin ID |
| `culture_id` | Specific culture/vessel lineage, left unassigned without a laboratory identifier |
| `measurement_type` | Serial passage measurement of an established population |
| `between_batch_relation` | Whether the same established population was sampled again |

The analysis will estimate matched expression differences, display all passage measurements and assess repeatability. It will not use the six timing matches as six independent LOY-generation replicates. The two sampling periods are not an independent training/validation pair. Any optional intervention analysis needs its own experimental-unit definition.
