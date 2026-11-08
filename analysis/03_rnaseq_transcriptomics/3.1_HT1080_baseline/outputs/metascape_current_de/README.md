# HT1080 current differential-expression lists for Metascape

These lists come from the current HT1080 limma-voom batch + condition model.
The contrast is Y-depleted minus sgCTRL. Genes were selected at overall BH FDR
< 0.05 and absolute log2 fold change >= 1, with matching effect directions in
August 2025 and May 2026.

Counts: 229 UP, 433 DOWN, 662 total. The custom
background contains 25461 unique tested gene symbols. Upload the UP
and DOWN text files as separate lists and use the tested-background text file as
the custom background. Each text file contains one HGNC-style symbol per line
with no header.

The fold-change threshold is a post-hoc filter. The reported FDR applies to the
non-zero-effect test, not a formal test that absolute log2FC exceeds 1. The
analysis is restricted to autosomal and non-PAR X genes; Y-linked and PAR genes
were excluded upstream. These samples are serial passages of one established
population per condition, so the list does not demonstrate reproducibility
across independently generated Y-depleted populations.

`HT1080_current_selected_DE_genes_with_statistics.csv` preserves the model
statistics and period effects. `manifest.json` records input/output hashes,
counts and the exact selection rule.

Reproduce from this analysis directory:

```powershell
python -X utf8 export_current_metascape_lists.py
```
