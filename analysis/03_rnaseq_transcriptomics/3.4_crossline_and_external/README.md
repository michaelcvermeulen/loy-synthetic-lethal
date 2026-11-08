# Cross-line RNA-seq comparison (HT1080, A549, DU145), 23 September 2026

Numerical source for Chapter 3 Results 3.4.7.4 ("Shared and background-specific molecular responses"). Author decisions are recorded in `SEPT_11_2026_THESIS/AGENTS.md`.

This analysis reads the saved outputs of the three thesis reanalyses and compares them. It never recomputes a per-line statistic: effects, FDRs, NES and BH are carried over unchanged, and `validate.R` checks this.
- HT1080 `SEPT_11_2026_REANALYZE`, using its `combined` equal-weight two-period contrast.
- A549 `SEPT_18_2026_REANALYZE`.
- DU145 `LNCAP/SEPT_23_2026_REANALYZE`.

All three share GENCODE v49 gene IDs, the `mean_difference` effect measure (group-mean log2 CPM difference) and effect-ranked fgsea on the frozen MSigDB 2026.1.Hs resources.

## Shared rule

A gene, gene set or regulator is **shared** when it has BH < 0.05 in at least two lines and the same sign in all three. It is **divergent** when it is significant in at least two lines with opposite signs. Continuous genome-wide concordance is reported alongside, so weak global agreement is not hidden.

## Stages

1. **`01_ingest`.** Loads inputs, hashes them (`outputs/tables/input_hashes.csv`) and defines the common universes: 16,165 primary genes by exact `gene_id`, and 22,961 gene sets tested in all three lines.
2. **`02_sex_chromosome`.**
   - Y inventory, PAR genes and region sums, and X paralogs.
   - Percent of control is the arithmetic ratio of mean CPM on each line's primary TMM scaling. It is computed per sampling period and then averaged, so HT1080 weights its two periods equally.
   - This replaces the plan's 100·2^mean_difference, because per-sample CPM is available for all three lines; it matches the DU145 prose values.
   - HT1080 CPM comes from `outputs/models/02_expression.rds`.
3. **`03_gene_level`.** Pairwise correlations, then two separate questions:
   - **Magnitude co-occurrence** among genes large (|effect| ≥ 0.5) in all three lines, tested against 1,000 independent within-line shuffles.
   - **Sign agreement** among those genes, tested with an exact binomial against 1/4.
   - Also: significant-gene overlap, and shared and divergent genes.
4. **`04_gene_sets`.** Classifies every commonly tested set.
   - **Themes** come from the union of leading-edge genes in the lines where the set is significant, not from set names. The reference groups are fixed in advance: Hallmark interferon/inflammatory; Hallmark G2/M, mitotic spindle and E2F; the DepMap GRHL2 neighbourhood; PRC2/H3K27me3 targets (`BENPORATH_PRC2_TARGETS`, `BENPORATH_ES_WITH_H3K27ME3`); ribosomal proteins; replication-dependent histones.
   - **Positional screen** in every line: a set is flagged when at least 50% of a line's leading edge lies on one arm. The DU145 copy-number rerun is carried over. HT1080 and A549 have no GENCODE-keyed copy-number table and are not copy-number-screened.
   - **HT1080 period sensitivity:** the AUG and MAY contrasts are checked separately.
   - **Post hoc sensitivity (disclosed):** a widened inflammatory reference that adds Hallmark TNFA/NFKB and IL6/JAK/STAT3. It was added after inspecting results, because many unassigned shared sets had cytokine leading edges just below the 20% threshold.
5. **`05_regulators`.** CollecTRI target-pattern scores under the same rule.
6. **`06_external`.** Replication across the in-house lines and 13 external datasets listed in `config/external_datasets.csv`, grouped as A (somatic LOY), B (45,X vs 46,XY) and C (not LOY: UTY-KO, UTX/UTY-DKO and sign-flipped XYY).
   - **Metric:** AUC that a member gene's Y-loss-polarity effect exceeds a non-member's. AUC < 0.5 means members are lower with Y absent.
   - **Significance:** an empirical two-sided P from 500 size-matched random sets.
   - **Polarity** is verified with a male-specific Y gene wherever one is present (7 of 16 datasets). Most external DE tables excluded Y genes before testing, so their direction rests on the contrast label; the status is recorded per dataset.
   - The mouse MB49 join uses a naive upper-case symbol match.

## Key results (validated run)

- **Gene level.** Correlations are 0.11, −0.02 and −0.06.
  - Large effects co-occur more than chance: 90 genes against a null of 13 (P = 0.001).
  - Their signs agree less than chance: 14 of 90, 16%, against 25% (binomial P = 0.038).
  - The shared genes are PXDN, THBD, BST2, CALD1, GLB1L2 and APLP2, all lower. THBD's DU145 decrease is copy-number-concordant (20p).
- **Gene sets.** NES correlations are about zero.
  - 46 sets are shared, all lower.
    - By primary theme: 28 interferon/inflammatory, 15 other, 1 each GRHL2 epithelial, G2/M and ribosomal.
    - Under the widened reference, 35 are interferon/inflammatory.
    - None is positionally concentrated, and 45 of 46 have the same sign in both HT1080 periods.
  - 23 sets are divergent: 10 G2/M (A549 down, DU145 up), 2 PRC2/H3K27me3 target sets and others.
  - The BMI1_DN / MEL18_DN sets are named for polycomb proteins, but their shared leading edges are secreted inflammatory genes, not PRC2 targets.
- **Regulators shared:** JUN, NFKB and STAT1 target patterns are all lower.

## Reproduce

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' --vanilla run.R
```

`provenance/validation.json` records the checks; `provenance/code_hashes.csv` and `sessionInfo.txt` record the environment. Thesis figures are built separately in `SEPT_11_2026_THESIS/figures/chapter3/code/`.
