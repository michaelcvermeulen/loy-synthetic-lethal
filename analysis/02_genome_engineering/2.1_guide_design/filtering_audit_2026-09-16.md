# Reconstruction of the filtering behind top_hits

Read-only audit of the four original CSVs in this folder, the related source scripts and the saved site-level off-target results. Original tables and thesis figures were not changed.

## Main finding

The evidence supports this reconstruction:

**14,723 candidate 20-nt guides → 6,901 with zero exact non-Y matches and at least 10 exact Y matches → a 215-guide shortlist with zero reported one- or two-mismatch non-Y sites → ordered by decreasing exact Y matches.**

The exact-match filter is explicitly present in the code and reproduces 6,901 rows. The one-/two-mismatch exclusion is a strong inference from every retained guide and its saved site records, but the original command that created the shortlist was not found. Therefore, it is not yet proven that applying those two mismatch exclusions alone to all 6,901 guides yields exactly 215.

## What is confirmed

| Stage | Setting or observation | Evidence |
|---|---|---|
| Starting guide table | 14,723 distinct 20-nt guides | `candidate_guides.csv` |
| Exact specificity | `sumOther == 0` | `build_off_target_tables.R:4`, `:59`; `complile_final_table.R:4` |
| Minimum Y abundance | `sumY > 9`, equivalent to **at least 10**, summing forward and reverse counts | Same explicit code; 6,901 rows retained |
| Shortlist input | Read `data/CRISPRchromosome/guides_of_interest.csv` and iterate over `g$guides_of_interest` | `complile_final_table.R:5`, `:114`, `:116` |
| Saved shortlist | 215 distinct guides, all contained in the 6,901-guide set | Exact sequence join; every exact-match count agrees |
| One-/two-mismatch counts | `mismatch_1 == 0` and `mismatch_2 == 0` for all 215 | Both CSV and corresponding site-level RDS |
| Three-/four-mismatch sites | Retained; reported counts extend to 3,941 and 13,840, respectively | `top_hits.csv` |
| Ordering | All rows are in decreasing `chrY_exact_hits` order; ties occur | Checked directly; ordering command not recovered |
| Condensed export | Same 215 guides in the same order, with identical shared columns; drops the two coordinate-list columns | `complile_final_table.R:168` and direct file comparison |

In this particular candidate table, all 6,901 guides with `sumOther == 0` already have `sumY >= 10`; the abundance condition removes no additional rows after the exact-specificity filter. Eleven final guides have exactly 10 Y matches, and 48 have 20 or fewer. A `>20` guide-level threshold cannot reproduce this shortlist.

The previously audited **7,470** entries are 23-nt protospacer–PAM sequences from a separate saved-object/filter branch. They must not replace the 6,901-guide denominator here.

## Most likely missing filtering command

Conceptually, after joining the full mismatch summaries to the eligible guide table:

```r
# Reconstruction, not a recovered original command.
eligible <- candidate_guides |>
  dplyr::filter(sumOther == 0, sumY > 9)

shortlist <- eligible_with_mismatch_counts |>
  dplyr::filter(total_1 == 0, total_2 == 0) |>
  dplyr::arrange(dplyr::desc(sumY))
```

`eligible_with_mismatch_counts` denotes the unavailable full joined table. The source calls the aggregate mismatch columns `total_1`–`total_4`; the exported table calls them `mismatch_1`–`mismatch_4`. Additional criteria cannot be ruled out without the original shortlist-generation command and the mismatch results for excluded guides. No explicit DeepSpCas9 or CFD cutoff for producing the 215 rows was recovered. The saved table includes a maximum site CFD of 1 and a cumulative CFD as high as 949.04.

The `CFD >= 0.8` condition at `complile_final_table.R:135` defines the **chrY_CFD_80 annotation**. It is not a guide-retention filter: the enclosing loop still returns one row for each input guide.

## Resolved: the general CFD columns already exclude chrY

`A:/Dropbox/LOY_paralogs/data/guides_of_interest_off_targets_list.RDS` contains 215 guide-specific tables and **92,870 site rows**. Their chromosome fields cover chr1–22, chrX and chrM, with **zero chrY rows**. The saved site CFD values reproduce all 215 `mismatch_mean_CFD`, `mismatch_max_CFD` and `mismatch_sum_CFD` values to floating-point precision (largest absolute difference <5 × 10^-13). The one-, two- and three-mismatch counts also agree exactly.

Thus, for this saved dataset, `mismatch_sum_CFD` is already a **non-Y** cumulative score. The later `plot_top_hits_locations.R`/`.Rmd` calculation `pmax(mismatch_sum_CFD - chrY_off_sum_CFD, 0)` subtracts an unrelated chrY total from the non-Y total. That subtraction should not be used to reconstruct the original shortlist or interpret non-Y burden. This resolves the earlier uncertainty about the scope of the general CFD input. It does not validate PAM/orientation handling or the biological calibration of the original scores.

There is one count discrepancy: guide `ACACGCATGCTCAAATTGAA` has 785 four-mismatch sites in `top_hits`, but 786 in the RDS. The extra site is mitochondrial (CFD 0.0285205). This is consistent with the nuclear-only aggregation implied by the 115-column slice in `build_off_target_tables.R:99`; the exact historical column order was not recovered. Its contribution is included in the saved CFD summaries. All other mismatch-count comparisons agree.

## Distinguish the later 20-row plotting subset

`top_hits_armcounts.csv` contains five anchors plus 15 other guides. Its guide order is reproduced exactly by taking the five saved anchors and the top 15 non-anchor guides under the later score, **without a p-arm restriction**. All 15 non-anchor rows have zero p-arm hits. The current plotting script instead has `parm_only <- TRUE` and six anchors, so its current settings do not reproduce this older 20-row export. This plotting subset is downstream of the 215-row `top_hits.csv`; it did not generate that shortlist.

## Files needed for an exact historical reproduction

The source refers to `data/CRISPRchromosome/guides_of_interest.csv` and the full per-chromosome mismatch inputs under `data/CRISPRchromosome/`. These were not found in the inspected local tree. The recovered site-level RDS contains only the 215 retained guides, so it cannot establish how many excluded guides also pass a proposed mismatch rule.

## Reproduce this audit

```bash
Rscript --vanilla code/audit_top_hits.R /path/to/sgRNA_guide_development /path/to/guides_of_interest_off_targets_list.RDS
```

The script writes small comparison tables under [filter_audit/](filter_audit/), without altering source data. It uses base R and reruns no genome search or CFD prediction.

Primary source code: [build_off_target_tables.R](A:/Dropbox/LOY_paralogs/code/build_off_target_tables.R:59), [complile_final_table.R](A:/Dropbox/LOY_paralogs/code/find_unique_seq/complile_final_table.R:4), and [later plotting script](A:/Dropbox/LOY_paralogs/code/plot_top_hits_locations.R:31).
