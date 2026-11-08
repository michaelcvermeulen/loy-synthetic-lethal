# Expanded screen GSEA — 21 October 2026

Author request: run GSEA on more gene sets than Hallmark. Exploratory; not in the Results prose.

**Scope (my call).** I used Hallmark, Reactome, GO Biological Process and GO Cellular Component from the same frozen MSigDB 2026.1.Hs resource as the Chapter 3 RNA analyses. These are the three Chapter 3 primary programme collections plus GO:CC, which holds protein complexes such as SAGA and SWI/SNF. Sets of 15–500 represented genes were eligible: 5,880 sets for the full rankings. Rankings and symbol mapping are reused unchanged from the [Hallmark analysis](../ht1080_hallmark_screen_gsea_20261021/README.md): −normZ, positive NES = sensitizer. The five rankings are HT1080 drugZ, HT1080 IsogenicZ, both HT1080 rankings without 1,526 common essentials, and A549 primary P4 drugZ. The primary FDR is BH across all eligible sets of the four collections, per ranking; within-collection BH is saved as a secondary column. fgsea 1.36.2 used the same settings and seed as before.

**Counts at global FDR ≤0.05.** HT1080 drugZ: 158 sets (58 sensitizer, 100 resistor). HT1080 IsogenicZ: 162. A549 drugZ: none. Sixty sets are significant in all four HT1080 rankings with the same sign. See [robust_in_all_four_ht1080_rankings.csv](../../figures/chapter4/source_data/screen_expanded_gsea_20261021/robust_in_all_four_ht1080_rankings.csv).

**Robust resistor sets (53).** These are dominated by mitochondrial translation and ribosome sets (NES −3.4; Reactome mitochondrial translation FDR 3×10⁻³¹), complex I biogenesis and the respiratory chain (NES −2.7 to −2.9), and ATP synthesis. They also include the cytosolic ribosome and translation (GOCC ribosome NES −2.6; Reactome translation NES −2.5), tRNA aminoacylation, nucleotide biosynthesis and interstrand cross-link repair. The Hallmark OXPHOS result is the same biology seen more broadly. **Caveat:** cytosolic translation and other core essential machinery are also resistor. LOY had lower essential-gene ROC AUC (0.907 versus 0.947), so part of the resistor signal may be weaker overall dropout in LOY rather than a mitochondrion-specific change. Removing common essentials does not remove the signal, and mitochondrial translation is stronger than cytosolic translation (NES −3.4 versus −2.5), but the analysis cannot fully separate the two.

**Robust sensitizer sets (7).**
- GOCC SAGA complex: NES 2.2, FDR 0.0037 (IsogenicZ 0.0018). The leading edge is exactly the six SAGA hits (TAF6L, TADA2B, SUPT20H, ATXN7L3, TADA1, TAF5L), so the set-level test restates the gene-level calls rather than adding independent evidence.
- SWI/SNF (GOCC SWI/SNF complex NES 2.1, FDR 0.0083), Reactome ATP-dependent chromatin remodelers and esBAF formation. The leading edge goes beyond the hit list: SMARCA4, SMARCE1, ARID1A, DPF2, ARID2, SMARCB1, BCL7A, SMARCD1, SMARCC1/2. This is new set-level support for a BAF chromatin-remodelling sensitizer theme alongside SAGA.
- Regulation of stem-cell population maintenance (TAF6L, SIRT6, TAF5L, SMARCA4, SMARCE1, ARID1A), and pattern specification and segmentation (led by MED12 and MTF2).

**Themes discussed earlier.**
- **Cholesterol/SREBP:** sets are positive (sensitizer direction) but not significant, for example Reactome SREBF/SREBP gene-expression activation NES 1.5, FDR 0.31 (IsogenicZ 0.22). The two SREBF hits do not extend to a pathway-level signal.
- **MYC:** Hallmark MYC targets V1 drops to FDR 0.12 under the larger correction (0.050 in the Hallmark-only family), and flips sign with IsogenicZ.
- **Mediator and HAT complexes:** not significant.

The Hallmark-only supplement and its FDR family are unchanged. The Results statement that six Hallmark sets met FDR ≤0.05 refers to the 50-set family.

**Outputs.** [Draft figure](../../figures/chapter4/review/ch04_supp_screen_expanded_gsea_draft.pdf): top 12 sets per direction among HT1080 drugZ global FDR ≤0.05, drugZ columns only. Source data in [source_data/screen_expanded_gsea_20261021](../../figures/chapter4/source_data/screen_expanded_gsea_20261021/): all fits with leading edges, the HT1080 significant-set wide table with sensitivity columns, the robust-set table, the selected-theme table, counts and a manifest. Reproduce with `run_expanded_screen_gsea.R`, then `summarise.R`.

**Author decision, 21 October 2026:** keep the draft figure as is: top 12 sets per direction among HT1080 drugZ global FDR ≤0.05, HT1080 and A549 drugZ columns. It is not restricted to the 60 robust sets, and redundant terms are not collapsed. The figure remains a review draft (`review/ch04_supp_screen_expanded_gsea_draft.*`), with no caption entry or Results placement yet.

**Author decision, 21 October 2026:** the draft now shows HT1080 only, with NES and FDR printed. This supersedes the two-column layout above.
