## Expanded all-collection GSEA

[All-collection results and figures](all_collections/report.html): 23,696 eligible sets; 3,281 overall global-FDR discoveries. The earlier three-collection analysis is retained.

# HT1080 RNA programme analysis — current implementation

Open [the searchable report](report.html). Decisions are in [the agreed plan](../../PLAN_20260918.md).

## Completed

- Core batch + condition DE; coherent period-specific/interaction fit; separate August parental comparison.
- Gene consistency, Hallmark/Reactome/GO-BP GSEA, programme consistency and leading-edge groups.
- STRING/MCODE modules, CollecTRI activity, promoter and motif evidence from existing ATAC results.
- All-theme overview, five provisional detailed-review dossiers, thesis-style overview figure, independent numerical validation.

## Main results

- 662 directionally consistent genes (229 higher, 433 lower).
- 220 directionally consistent enriched pathways; 55 significant in both periods; 88 pathway groups including isolated terms.
- 60 directionally consistent TF patterns; five STRING association modules.
- 12 sgCTRL-versus-parental genes at both cutoffs; five overlap the primary gene set.
- Two interaction-significant genes in the factorial model versus 451 in the earlier passage-aware fit: major model sensitivity.
- CAMERA sensitivity: 17 discoveries at fixed correlation 0.01, zero with estimated nonnegative correlation. GSEA interpretation remains exploratory.

## Remaining interpretive work

Review theme labels and final 3–5 selections; inspect locus-level evidence and supported enhancer links; select broader signatures for follow-up; prepare detailed thesis figures and prose. RNA–ATAC specimen matching remains unconfirmed; current integration is at group level. The Y/PAR branch remains separate and unchanged.

## Reproduce

From this analysis directory, run with Rscript --vanilla and Python:

1. run_programme_analysis_v1.R
2. build_programme_networks_v1.py
3. check_programme_camera_v1.R
4. integrate_programmes_ATAC_v1.py
5. plot_programme_analysis_v1.py
6. report_programme_analysis_v1.py

Original inputs remain unchanged; input/code hashes and validation are in provenance/. All CSV tables retain full numerical precision. Parentals are independently filtered and normalized, and no PCA batch-adjusted values enter DE.
