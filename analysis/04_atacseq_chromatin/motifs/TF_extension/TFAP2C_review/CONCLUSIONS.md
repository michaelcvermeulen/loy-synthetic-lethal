TFAP2C-focused review

Conclusion: LOY is associated with selective loss of accessibility at AP-2 motif-containing regulatory regions, accompanied by lower TFAP2C RNA and a modest shift in pooled predicted occupancy. This supports TFAP2C as a candidate contributor, but sequence motifs do not distinguish it from related AP-2 proteins.

All three TFAP2C models agree in direction. HOMER lost-peak enrichment is 1.89-2.12-fold, with FDR 2.0e-15 to 2.8e-7; no primary gained enrichment passes 0.05. The earlier JASPAR2022 CMH screen supports the same direction, with different motif versions and target definitions.

All LOY chromVAR scores are below all CTRL scores for every TFAP2C model. Primary limma FDR is 0.058-0.067; passage sensitivity FDR is 0.000125-0.000299, while Welch FDR is 0.255-0.276. Omitting any one sample preserves negative mean effects. Direction is robust; the significance threshold depends on statistical model. Technical batch is unconfirmed.

Across the union of 20,306 motif-positive peaks, 298 lose and 14 gain accessibility at FDR <0.05 and |log2FC| >1. Median accessibility log2FC is -0.0048; 50.6% of peaks have negative effects. These results do not imply uniform closure of all AP-2 regions.

TOBIAS predicts 2.1-6.1% fewer bound sites, depending on motif model, with both losses and gains at individual sites. Aggregate corrected cleavage profiles largely overlap and do not show a clear condition-specific loss of central protection. Its pooled-site significance is not a biological-replicate test and should not be treated as an independent validation experiment. Footprint-derived scores depend on local accessibility as well as cleavage patterns.

99.14% of TFAP2C-positive peaks also contain another AP-2 model somewhere in the same peak. This is peak-level overlap, not proof of identical binding positions. TFAP2A RNA is unchanged statistically; TFAP2B is absent from the selected RNA table. TFAP2C RNA log2FC is -0.747 (FDR 0.0117), approximately 40.4% lower. Both annotated TFAP2C promoter peaks are nonsignificant in primary ATAC.

At the strict DA cutoff, lost motif-containing promoter peaks link to 39 RNA-mapped genes: 17 have significantly lower RNA at FDR <0.05, none higher. With both assays at FDR <0.10 and no ATAC effect-size filter, 32 are RNA-down and 2 RNA-up (BMP7 and SLC9A2). These are positional candidate links, not experimentally established TFAP2C targets or pathway enrichment.

No copy-number adjustment has been applied. The assays are compared at group level; the three ATAC-based methods reuse the same data and are not independent biological experiments.
