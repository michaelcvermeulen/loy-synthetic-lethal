# DU145 (filed as LNCAP) SEPT_23_2026_REANALYZE: decisions

## 27 September 2026: TF regulators and shared-target TF modules (tf_modules_v1), rules fixed before scoring

The author asked for the HT1080 tf_modules_v1 analysis (HT1080 `PLAN_20260924_tf_modules.md` and `provenance/DECISIONS.md`, 24 Sep 2026) to be applied to DU145. HT1080 results, and the earlier DU145 CollecTRI ULM(logFC) results including the GRHL2/ZEB1 self-edge circularity, were known when these rules were written; no DU145 tf_modules_v1 score had been computed. The rules below must not be changed after DU145 results are seen. Any later deviation is added as a dated note with its reason.

**Unchanged from HT1080**
- Priors: the HT1080 `resources/tf_modules_v1` files (CollecTRI current = frozen v1, DoRothEA A-C, TFActProfiler all and fibroblast), copied into `resources/tf_modules_v1` and verified against the HT1080 SHA-256 manifest. The DU145 `resources/collectri_network.csv` is byte-identical to HT1080's.
- Methods: decoupleR ULM, MLM, WSUM (norm_wsum, 1,000 permutations, seed 42) and VIPER (analytic, pleiotropy=FALSE), minsize 5. Significance family: ULM(t) with BH across TFs per contrast. Consensus score for ranking and display only. limma::fry per regulon as a correlation-aware flag, not part of the rule.
- R1, R3, R4 and R6 as written for HT1080, including the AP1/NFKB complex-subunit self-edge definition for R6. R6 removes self-edges, so the known GRHL2 self-edge circularity is handled by the rule rather than by a special case.
- Tiers: robust = R1-R6 all true; supported = R1 and R2 plus at least 3 of R3-R6.
- Module stage: as in HT1080 (see the A549 and HT1080 DECISIONS entries for the full list).

**DU145-specific (author choice, 27 Sep 2026: line-specific substitutes with comparable tiers)**
- Design: three sgCTRL and three LOY libraries at matched passages P12, P13 and P14, one May 2026 window, no parental libraries. Overall input = moderated t for conditionLOY from the saved paired `~ passage + condition` voom fit (`02_expression.rds` `de`, residual df 2).
- **R2 substitute:** R2 = ULM on each saved passage-pair logCPM difference (`pair_diff`, P12, P13, P14) has the same sign as the overall ULM(t) score in all 3 pairs. This is the regulon-level analogue of the DU145 thesis DE rule (same direction in all three passage pairs).
- **R5 substitute:** 3 leave-one-pair-out refits of the paired voom model (4 libraries each, residual df 1), ULM(t) same sign in all 3. The full refit must reproduce the saved moderated t (max difference < 1e-8).
- Passages are repeated sampling of one population per arm. R2 and R5 measure stability across passages, not independent editing replication, and are not equivalent to the HT1080 cross-period agreement.
- Drift flag: not applicable (no parental libraries).
- Copy number: the established DU145 gene-level rule (`05_copy_number.rds` `CN_changed`: |LOY-2X minus sgCTRL| >= 0.2 median QDNAseq log2_seg on CHM13v2.0), not the HT1080 arm rule, because this line already has a documented gene-level screen and a focal 20p deletion that an arm median would dilute. A TF is flagged if half or more of its leading contributors sit on one arm or if its ULM(t) loses sign or p < 0.05 after removing CN-changed targets. Specimens were not matched to the RNA libraries.
- Line-restricted similarity sensitivity (HT1080 used its 662 directionally consistent genes): genes with FDR < 0.05 and the same direction in all three pairs, without the fold-change threshold.
- Annotation references: that consistent set, the 23 thesis DE genes, and the HT1080 module 1 driver targets.
- Co-activity residuals remove `~ passage + treatment`.

**Shortfall rule (author choice, 27 Sep 2026)**
The node rule is not relaxed. With fewer than 3 nodes, clustering is skipped, per-regulator tables are still written, and a shortfall table states the reason. If clustering yields only singletons, no module is reported.
