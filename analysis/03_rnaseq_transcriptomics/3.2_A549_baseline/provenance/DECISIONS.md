# A549 SEPT_18_2026_REANALYZE: decisions

## 27 September 2026: TF regulators and shared-target TF modules (tf_modules_v1), rules fixed before scoring

The author asked for the HT1080 tf_modules_v1 analysis (HT1080 `PLAN_20260924_tf_modules.md` and `provenance/DECISIONS.md`, 24 Sep 2026) to be applied to A549. HT1080 results were known when these rules were written; no A549 tf_modules_v1 score had been computed. The rules below must not be changed after A549 results are seen. Any later deviation is added as a dated note with its reason.

**Unchanged from HT1080**
- Priors: the HT1080 `resources/tf_modules_v1` files (CollecTRI current = frozen v1, DoRothEA A-C, TFActProfiler all and fibroblast), copied into `resources/tf_modules_v1` and verified against the HT1080 SHA-256 manifest. The A549 `resources/collectri_network.csv` is byte-identical to HT1080's.
- Methods: decoupleR ULM, MLM, WSUM (norm_wsum, 1,000 permutations, seed 42) and VIPER (analytic, pleiotropy=FALSE), minsize 5. Significance family: ULM(t) with BH across TFs per contrast. Consensus score for ranking and display only. limma::fry per regulon as a correlation-aware flag, not part of the rule.
- R1, R3, R4 and R6 as written for HT1080, including the AP1/NFKB complex-subunit self-edge definition for R6.
- Tiers: robust = R1-R6 all true; supported = R1 and R2 plus at least 3 of R3-R6.
- Module stage: nodes are robust TFs if at least 8, otherwise robust plus supported; weighted leading-contributor overlap, average-linkage clustering by best silhouette, bootstrap and edge-drop stability, degree-preserving rewired and size-matched nulls, expression-matched target-union permutation plus fry, confounder programmes (the same eight Hallmark sets and the 14 van den Brink 2017 Fig. 1 genes).

**A549-specific (author choice, 27 Sep 2026: line-specific substitutes with comparable tiers)**
- Design: four sgCTRL and four LOY libraries, one sampling batch, pairing and derivation uncertain (18 Sep 2026 decision). Overall input = moderated t for conditionLOY from the saved unpaired `~ condition` voom fit (`02_expression.rds` `de`, residual df 6).
- **R2 substitute:** there is no second sampling period. R2 = ULM on the descriptive group-mean difference with each of the eight libraries omitted (the saved `omit_*` effects, the basis of the existing stable-omission rule) has the same sign as the overall ULM(t) score in all 8 omissions.
- **R5 substitute:** 8 leave-one-sample-out refits of the `~ condition` voom model (7 libraries each), ULM(t) same sign in all 8. The full refit must reproduce the saved moderated t (max difference < 1e-8).
- R2 and R5 both draw on sample omission, so in A549 they are correlated. Neither is independent replication, and neither is equivalent to the HT1080 cross-period agreement. A549 tiers are therefore labelled as stability tiers, not replication tiers.
- Drift flag: CTRL_vs_parental moderated t from a new voom fit on all 11 libraries (`~ 0 + group`, the saved context TMM exposures, the primary A549 gene universe). It is a flag only.
- Copy number: arm-level rule as in HT1080 (median QDNAseq log2_seg for LOY minus sgCTRL at least 0.15 in absolute value; `WGS/A549_QDNAseq_cn_df.rds`; hg38 cytoband centromeres). The WGS specimens are not matched to the RNA libraries, and the A549 README defers coordinate integration; the flag is a screen, not a dosage account.
- Line-restricted similarity sensitivity (HT1080 used its 662 directionally consistent genes): the saved A549 descriptive candidates (0.5 logCPM plus stable omission direction).
- Annotation references: the A549 descriptive candidates and the HT1080 module 1 driver targets. HT1080 STRING modules and programme communities do not exist for A549 and are not used.
- Co-activity residuals remove `~ treatment` only.

**Shortfall rule (author choice, 27 Sep 2026)**
The node rule is not relaxed. With fewer than 3 nodes, clustering is skipped, per-regulator tables are still written, and a shortfall table states the reason. If clustering yields only singletons, no module is reported.
