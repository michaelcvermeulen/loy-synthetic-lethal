# Repository Guidelines

## Project Structure & Data Organization

This repository stores flow-cytometry assay data rather than application source code. `comp_assay_qc_run_AUG_13_2026/` is the experiment bundle:

- `*.fcs` files are raw acquisitions for control (`CTRL`) and loss-of-Y (`LOY`) samples, with GFP, mCherry, and mixed conditions.
- `comp_assay_qc_run.xit` is the associated instrument/analysis workspace.
- `ExpSummaryForAPI.xml` records the plate, tube names, and signal-statistics settings.
- `Backup/` is reserved for instrument-generated backups.

Keep files from a run together. Place future runs in separate, date-stamped directories; do not mix acquisitions between runs.

## Validation and Local Commands

There is no build system or automated test suite. Use PowerShell for basic checks:

```powershell
Get-ChildItem .\comp_assay_qc_run_AUG_13_2026\*.fcs
[xml](Get-Content .\comp_assay_qc_run_AUG_13_2026\ExpSummaryForAPI.xml) | Out-Null
Get-FileHash .\comp_assay_qc_run_AUG_13_2026\*.fcs -Algorithm SHA256
```

These commands inventory acquisitions, confirm that the summary is well-formed XML, and produce hashes for transfer-integrity checks. Open `.fcs` and `.xit` files only with compatible cytometry software.

## Naming and Data-Handling Conventions

Preserve instrument-generated filenames and binary contents. Existing names follow `<plate>-<condition>-<well>.fcs`, for example `01-LOY_GFP-A4.fcs`. Use uppercase `CTRL` and `LOY`, retain channel labels (`GFP`, `mCherry`), and use well identifiers such as `A1` or `B8`. Avoid ambiguous suffixes such as `(1)` for new files; document why existing duplicates are retained. Format hand-edited XML with two-space indentation and UTF-8 encoding.

## Testing Changes

After any metadata change, parse the XML and verify every `<Tube Name>` maps to the intended `.fcs` filename. For added acquisitions, record file count, sizes, and SHA-256 hashes. Never rewrite raw `.fcs` files as part of validation.

## Commits and Pull Requests

No Git history is available here to infer an established convention. Use short, imperative commits such as `Add August QC run metadata`. Pull requests should identify the run and instrument, summarize files added or changed, explain duplicates or renames, and include validation results. Do not commit patient identifiers, credentials, temporary exports, or unrelated analysis output.
