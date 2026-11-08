# HT1080 molecular integration — Supp Fig 4.7

Updated 28 October 2026 in response to the author's request to clarify the different meanings of colour and the biological row groups.

Current exports: [PDF](../../figures/chapter4/ch04_supp_ht1080_molecular_integration.pdf), [editable SVG](../../figures/chapter4/ch04_supp_ht1080_molecular_integration.svg), [600-dpi PNG](../../figures/chapter4/ch04_supp_ht1080_molecular_integration.png). Canvas: 9.6 × 10.5 inches; embedded Open Sans throughout.

## Display changes

- Screen columns use teal for negative values (resistors; relative enrichment in LOY) and purple for positive values (sensitizers; relative depletion in LOY). Baseline and add-back columns use blue for lower and red for higher expression/accessibility. Each panel has two separately labelled horizontal colourbars.
- Panel a displays seven columns, matching panel b. Only the untransduced-LOY-versus-EV column was omitted from the display. Its data remain in the unchanged original tables.
- Panel a contains three labelled groups, separated by subtle horizontal rules: screen sensitizer complexes; downstream regulated target programmes; screen resistor pathways. All 15 pathway rows and 17 gene rows are retained in their original order.
- Numerical cell values and FDR-based emphasis are preserved. Small FDRs use compact e notation in the new renderer. Light text on dark cells improves contrast; it carries no additional significance meaning.
- Display limits are unchanged: ±8 for signed −log10(FDR), ±8 for screen −normZ, and ±1.5 for gene-level log2 fold change. Printed values remain unclipped. Different tests and correction families remain distinct.

This is a display-only refinement of the saved 21 October integration tables. No input table, analysis, FDR, or underlying comparison was changed. The baseline and add-back results retain the interpretation limits documented in the [original source notes](../ch4_molecular_integration_20261021/README.md).

## Reproduction and validation

Run `figures/chapter4/code/make_ht1080_molecular_integration.py`, then `figures/chapter4/code/validate_ht1080_molecular_integration.py` using the local figure Python dependencies. The renderer reads the original `figure_panel_a_values.csv` and `figure_panel_b_values.csv` from `figures/chapter4/source_data/ch4_molecular_integration_20261021/`. It writes display tables and checksums to [the new provenance directory](../../figures/chapter4/source_data/ht1080_molecular_integration_supp_20261028/).

[Validation](../../figures/chapter4/source_data/ht1080_molecular_integration_supp_20261028/validation.json) confirms seven source files are byte-identical, all 105 retained pathway cells and 119 gene cells match the saved tables, four colourbars and three row categories are present, PDF fonts are embedded Open Sans, SVG text remains editable, PNG is 600 dpi, and no text overlaps or lies outside the page. The complete preview was visually inspected.

## Captions and numbering

The [central caption](../../Chapter4_figure_captions.md#supp-fig-4-7) and draft 2 inline caption are synchronised. Results had already been renumbered by the author; only this figure block was replaced, and all surrounding Results prose was verified unchanged. Central captions, the figure register and remaining chapter references now follow that sequence. Former main 4.5–4.8 are 4.4–4.7; the Discussion figure moves from 4.9 to 4.8. Former supplementary 4.7–4.11 are 4.8–4.12. Provisional public-analysis Fig 4.0/Supp Fig 4.0a–d retain their existing labels.

Previous Markdown, figure exports, source tables and relevant scripts are archived at `archive/ht1080_molecular_integration_supp_20261028_153527/`. The precise path is also stored in [archive_path.txt](archive_path.txt). The one-time integration script checks the archived hashes before writing; do not rerun it over subsequent author edits. Full-chapter DOCX/PDF exports were not regenerated during this figure revision.
