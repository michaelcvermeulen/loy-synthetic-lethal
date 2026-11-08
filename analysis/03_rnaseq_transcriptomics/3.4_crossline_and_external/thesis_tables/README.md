# Cross-line gene-set results

These tables preserve the saved enrichment results behind Results 3.3.4.4 after the 69-row supplementary heatmap was removed from the narrative on 8 October 2026.

- `primary_collection_results.tsv`: all 4,634 commonly tested Hallmark, Reactome and GO:BP sets, with collection identifiers, NES, original within-background BH values and shared/divergent classification.
- `shared_and_divergent_gene_sets.tsv`: the 46 shared and 23 divergent sets from the expanded 22,961-set comparison, including the original statistics, leading-edge genes and sensitivity annotations.

Shared means FDR <0.05 in at least two backgrounds and a matching enrichment direction in all three. Divergent means significance in at least two backgrounds in opposing directions. FDRs retain the original correction across all eligible sets within each background; they were not recalculated for these tables or the seven selected figure rows.

Source: `RNAseq/analysis/crossline_SEPT_2026/outputs/tables/gene_sets_all_common.csv` and `gene_sets_shared_divergent.csv`, relative to the projects directory. Regenerate these copies with `figures/chapter3/code/make_crossline_concise.py` from the thesis root. The input hashes are in `figures/chapter3/source_data/crossline_concise/input_manifest.json`.

The comparison uses mean-expression effects and effect-ranked enrichment. Its earlier condition-only A549 analysis differs from the current RINe-adjusted gene tests; HT1080 periods are equally weighted, and this ranking differs from the primary moderated-t enrichment analysis. These are descriptive comparisons between established populations and do not independently establish pathway activity or causal effects of Y loss.
