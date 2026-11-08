# DU145 section supporting tables

These are unchanged saved analysis results used to check the revision of Results 3.4.6.3 on 2 October 2026. No statistical model was refitted.

- [Primary gene-level results](DE_paired_voom_primary.csv): passage-blocked limma-voom; BH across 18,767 expressed autosomal/non-PAR X-linked genes. There are 147 FDR <0.05 genes. The current 114-gene reporting subset additionally uses |log2FC| ≥0.5 and agreement across three passage pairs; legacy thesis_DE flags in the original table use an older cutoff, so use figures/chapter3/code/du145_reporting.py for the current subset.
- [Y/PAR results](Y_PAR_gene_table.csv): percentages use arithmetic mean CPM; BH_all_expressed_nuclear corrects across 18,808 nuclear genes, while BH_Y_PAR_family corrects across 41 tested Y/PAR features (16 non-PAR Y and 25 PAR). Untested features are retained with missing statistics. PAR counts do not distinguish X- and Y-derived transcripts.
- [Highlighted default-CAMERA results](highlighted_CAMERA_results.csv) and [full default-CAMERA table](camera_all_collections.csv): fixed inter-gene correlation 0.01, passage-blocked model, BH_across_collections reported. The four highlighted signatures are significantly lower; 763 sets pass global FDR <0.05 overall.
- [Complete GSEA results](fgsea_all_collections.csv): effect-ranked fgsea, 23,252 eligible sets, full-collection BH and copy-number-exclusion sensitivity. There are 154 significant sets, of which 139 retain significance and direction after excluding genes in altered intervals.
- [Copy-number-restricted gene-level refit](../ch03_supp_du145_copy_number_stable_DE.csv): 16,719 genes, 91 at FDR <0.05 and 76 meeting the current effect-size and passage-consistency criteria. This is sensitivity analysis, not an additional main discovery set.

Primary source: RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE/outputs/tables, relative to the projects root. The regional plots retain their saved source values. Current captions: [Chapter3_figure_captions.md](../../../Chapter3_figure_captions.md).
