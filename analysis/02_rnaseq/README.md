# 02_rnaseq/

Bulk RNA-seq across the three engineered isogenic systems (A549, HT1080,
LNCaP), plus a cross-cell-line integrated analysis. Each cell line gets
three sub-numbered notebooks: **base** (canonical DE + GSEA + TF activity),
**extended** (robustness — filter sweep, voomQW, treat, edgeR-QL, camera,
paralog deep-dive, GSVA/PROGENy pathway scoring, lpWGS CN integration,
cross-framework robust DEGs), and **VIPER** (master-regulator inference via
CollecTRI + TCGA-SARC ARACNe regulons).

Pipeline: STAR two-pass + featureCounts (v2.1.1, GENCODE v49) →
`limma`–voom for DE → `fgsea` + `decoupleR` for enrichment / TF activity →
`VIPER` / `msVIPER` for master regulators.

**Counts pipeline (current).** All six base + extended notebooks load from
`projects/RNAseq/counts/STAR_v2_PARfix/<batch>/`. STAR v2 was re-aligned
against a PAR-Y-masked chrY reference + featureCounts re-quantified to
recover chrX-PAR reads (SHOX, CSF2RA, SLC25A6, AKAP17A, etc.) that the
first STAR run silently dropped at multi-mapper filtering. MSY counts are
unchanged across the two pipelines. See `projects/RNAseq/code/par_requant/`
for the upstream re-quantification recipe.

**Sex-chromosome filter policy (current).** Each notebook defines
`y_panel` + `excluded_gene_ids` (chrY ∪ PAR) + `par_gene_names` in §1.4.5
(extended) / §2.1.5 (base). The two filter sets are applied at different
boundaries:

- **DE flow** (limma-voom topTable, treat, edgeR-QL, eBayes): only chrY
  (MSY) is stripped from the DGEList. **PAR genes stay** so they get
  logFC + FDR in §4 volcanos like every other gene. The new §3.5 (base) /
  §4.0 (extended) Y-region volcanos colour-code these by region (PAR1 /
  PAR2 / MSY-Yp / MSY-Yq).
- **Pathway flow** (fgsea, camera, decoupleR ULM, DoRothEA, GSVA / ssGSEA,
  PROGENy, msVIPER): both chrY and PAR are stripped at the input boundary.
  `prep_ranks_t()` / `prep_ranks()` filter cascades into rank-based calls;
  `lcpm_lc` / camera voom / decoupleR `deg_mat` are filtered inline at
  build time.

## Planned analyses

| #     | File                                       | Description |
|-------|--------------------------------------------|-------------|
| 2.0   | `2.0_pipeline_comparison.Rmd`              | QC: compares OLD (HISAT2 / GENCODE v47 / featureCounts v2.0.1) vs NEW (STAR / GENCODE v49 / fc 2.1.1) count tables for HT1080 + A549. Library sizes, per-sample logCPM correlation, per-gene scatter, sample-correlation heatmap, DE concordance with logFC Pearson + Jaccard, automated verdict. **migrated** |
| 2.1.1 | `2.1.1_rnaseq_A549_base.Rmd`               | A549 base — paired limma–voom DE + fgsea + DoRothEA/decoupleR TFs + pooled and sgCTRL-vs-PARENTAL contrasts + msVIPER + cross-source TF consistency. Loads STAR_v2_PARfix counts; §3.5 Y-region volcano (PAR1/PAR2/MSY-Yp/MSY-Yq). **migrated** |
| 2.1.2 | `2.1.2_rnaseq_A549_extended.Rmd`           | A549 extended — robustness layer (filter sweep, voomQW, treat, edgeR-QL, camera × paired/unpaired/pooled designs; X/Y paralog deep-dive ±PARENTAL; GSVA/PROGENy/decoupleR sample scoring; lpWGS CN integration + cutoff sweep + per-design CN-flagged DEG sensitivity; cross-design × cross-framework robust-DEG panel; §11 LOY transcriptomics summary with PAR/MSY/Yp/Yq signatures). Loads STAR_v2_PARfix counts; §4.0 Y-region volcano. **migrated** |
| 2.1.3 | `2.1.3_rnaseq_A549_viper.Rmd`              | A549 VIPER — master-regulator inference. TODO. |
| 2.2.1 | `2.2.1_rnaseq_HT1080_base.Rmd`             | HT1080 base — paired limma–voom DE + fgsea + DoRothEA/decoupleR TFs + pooled and sgCTRL-vs-PARENTAL contrasts + msVIPER + cross-source TF consistency. Loads STAR_v2_PARfix counts; §3.5 Y-region volcano. **migrated** |
| 2.2.2 | `2.2.2_rnaseq_HT1080_extended.Rmd`         | HT1080 extended — robustness layer (filter sweep, voomQW, treat, edgeR-QL, camera × paired/unpaired/pooled designs; X/Y paralog deep-dive ±PARENTAL; GSVA/PROGENy/decoupleR sample scoring; lpWGS CN integration + cutoff sweep + per-design CN-flagged DEG sensitivity; cross-design × cross-framework robust-DEG panel; §11 LOY transcriptomics summary). Loads STAR_v2_PARfix counts; §4.0 Y-region volcano. **migrated** |
| 2.2.3 | `2.2.3_rnaseq_HT1080_viper.Rmd`            | HT1080 VIPER — master-regulator inference. TODO. |
| 2.3.1 | `2.3.1_rnaseq_LNCaP_base.Rmd`              | LNCaP base — paired limma–voom DE + fgsea + DoRothEA/decoupleR TFs + msVIPER (CollecTRI + TCGA-PRAD ARACNe) + cross-method TF consistency. No PARENTAL arm collected for LNCaP. Loads STAR_v2_PARfix counts; §3.5 Y-region volcano. **migrated** |
| 2.3.2 | `2.3.2_rnaseq_LNCaP_extended.Rmd`          | LNCaP extended — robustness layer. Same scaffold as HT1080/A549 but only **paired + unpaired** designs (no PARENTAL arm → no pooled-Y+ design). §10.4 robust-DEG panel uses 6 cells (voom/treat/edgeRQL × paired/unpaired). lpWGS CN integration uses the `sgCTRL` + `LOY-2X` arms of the LNCaP QDNAseq panel. Loads STAR_v2_PARfix counts; §4.0 Y-region volcano. **migrated** |
| 2.3.3 | `2.3.3_rnaseq_LNCaP_viper.Rmd`             | LNCaP VIPER. TODO. |
| 2.4   | `2.4_rnaseq_combined.Rmd`                  | Cross-cell-line integration — combined PCA (raw + batch-corrected), DEG overlap (strict Venn/UpSet + ranked Spearman/scatter), Hallmark / DoRothEA / decoupleR pathway-TF NES overlap (heatmap + Venn + pairwise scatter + sign-agrees-all-3 concordant calls). §10 msVIPER cross-line is a stub until per-cell VIPER notebooks land. **migrated** |

## Output layout

Each cell-line notebook writes outputs under a same-named subdirectory next
to the Rmd, namespaced so renders never collide:

```
analysis/02_rnaseq/
├── 2.2.1_rnaseq_HT1080_base.Rmd
├── HT1080/                          ← outputs from 2.2.1 + 2.2.2 + 2.2.3
│   ├── DE/                          ← from base
│   ├── GSEA/                        ← from base
│   ├── TF_analysis/                 ← from base
│   ├── figures_LOY_sgCTRL/          ← base figures (PNG + PDF)
│   ├── DE_extended/                 ← from extended (later)
│   ├── PATHWAY_SCORES/              ← from extended (later)
│   ├── CN_DEG/                      ← from extended (later)
│   ├── figures_LOY_sgCTRL_extended/ ← extended figures (later)
│   └── VIPER/                       ← from VIPER notebook (later)
└── A549/  LNCaP/  combined/         ← parallel structure once those are migrated
```

CSV/RDS outputs are gitignored repo-wide; PNG/PDF figures are tracked.

## Data path

Every notebook resolves data via a single `DATA_ROOT` constant near the top
of the setup chunk (default: `A:/Dropbox/LOY_paralogs/projects`). Update
that single line if data moves to a different mount or, post-publication,
to a public mirror of the GEO deposit.

## Existing upstream work (sources for the migration)

- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/A549/` — 3 `.Rmd` (base, extended, viper).
- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/HT1080/` — 3 `.Rmd`.
- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/LNCAP/` — 3 `.Rmd`.
- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/combined/` — `analysis_combined.Rmd`.
- `a:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/_shared/` — `viper_helpers.R`, `dcc_barplot.R`, `dcc_expression.R` (sourced by VIPER notebooks).
- `a:/Dropbox/LOY_paralogs/projects/RNAseq_LOY/` — older/exploratory outputs and DE tables; audit for any analyses to bring forward vs. supersede.
