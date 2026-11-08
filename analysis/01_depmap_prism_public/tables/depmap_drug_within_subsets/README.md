# Drug-response associations within subsets

Author-requested reporting revision removes drug interactions while retaining the existing direct LOY-versus-Y-positive fits within binary cancer and genomic subsets. This table contains 495,079 finite within-subset tests across 121 subsets. The original coefficients, nominal P values and within-class FDRs are unchanged. The only FDR <0.05 association is higher response / lower sensitivity to SR-1664 in LOY mature B-cell neoplasm models in PRISM REP.1M.

The original comparison set required five measured responses per Y group both inside and outside the subset. No models were refitted and no additional comparisons were introduced. Four continuous genomic interaction contexts have no corresponding binary subset comparison. Original HPC outputs remain the historical analysis record.

Reproduce with research/depmap_loy_reanalysis_v1/summarise_drug_within_subsets.py from the thesis root. Methods and Results describe the retained within-subset analysis.
