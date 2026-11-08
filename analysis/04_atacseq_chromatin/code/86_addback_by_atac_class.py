"""Does promoter closing mark the genes that UTY add-back restores? (HT1080, 22 Sep 2026)

Question: among genes with lower expression in Y-depleted cells (thesis RNA DE, 'down'), are those that also lost
promoter accessibility (>= 1 lost promoter peak at ATAC FDR < 0.10; the 58 concordant genes) restored more by add-back
than down-regulated genes whose promoters did not close, once baseline effect size and abundance are accounted for?

Inputs (no refit):
  thesis integration universe (13,300 genes; SEPT_11_2026_THESIS figure source data, integration_genes.csv)
  add-back limma-voom contrasts (RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026/outputs/de/), each vs the empty
    vector LXHA: UTY, KDM5D, DUAL (UTY + KDM5D) and untransduced Y-depleted (LOY; construct-free reference)
  polycomb set membership (msigdbr 26.1.1 C2:CGP; exported by Rscript --vanilla to polycomb_set_membership.csv)
Restoration for a down-regulated gene = add-back log2FC vs vector (positive = towards sgCTRL).
Tests, per contrast:
  1. OLS (HC3 SEs): restoration ~ closed + baseline RNA log2FC + add-back AveExpr, DE-down genes tested in add-back.
  2. Matched: each closed gene vs its K = 2 nearest non-closed DE-down genes (standardised baseline log2FC and
     AveExpr; greedy, without replacement, closed genes processed in random order, seed 20260922); Wilcoxon signed-rank
     on closed - mean(matched) differences.
  3. Genome-wide: add-back log2FC ~ baseline RNA log2FC + median promoter ATAC log2FC + AveExpr (HC3), all universe
     genes tested in add-back.
  4. Polycomb: restoration ~ closed + ES_H3K27me3 + closed:ES_H3K27me3 + covariates, DE-down genes.
Contrasts 5-8 compare populations directly: three differ in UTY, and KDM5D vs untransduced (no UTY difference)
tests whether closed genes simply differ more between any two populations. The absolute-shift comparison
(|log2FC| closed vs matched) is reported for every contrast. BH across all contrasts within each test family. Add-back and baseline RNA are separate experiments and
the three add-back samples per group are repeated collections of one population, not independent transductions.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
PROJ = ROOT.parents[1]
OUT = ROOT / "integration/thesis_RNA_SEPT11/addback_by_atac_class"
OUT.mkdir(parents=True, exist_ok=True)
GENES = PROJ / "SEPT_11_2026_THESIS/figures/chapter3/source_data/ht1080_rna_atac_integration/integration_genes.csv"
DE = PROJ / "RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026/outputs/de"
CONTRASTS = {"UTY": "UTY_v_LXHA", "KDM5D": "KDM5D_v_LXHA", "UTY+KDM5D": "DUAL_v_LXHA", "Untransduced Y-depleted": "LOY_v_LXHA",
             # direct comparisons between populations with and without UTY, and a UTY-free control pair
             "UTY vs KDM5D": "UTY_v_KDM5D", "UTY vs untransduced": "UTY_v_LOY", "UTY+KDM5D vs untransduced": "DUAL_v_LOY",
             "KDM5D vs untransduced (no UTY difference)": "KDM5D_v_LOY"}
K, SEED = 2, 20260922

g = pd.read_csv(GENES)
pc = pd.read_csv(OUT / "polycomb_set_membership.csv")
for s, col in [("BENPORATH_ES_WITH_H3K27ME3", "ES_H3K27me3"), ("BENPORATH_PRC2_TARGETS", "PRC2")]:
    g[col] = g.entrez_id.isin(pc[pc.gs_name == s].ncbi_gene).astype(int)
g["closed"] = (g.lost_010 > 0).astype(int)
trd = pd.read_csv(ROOT / "differential/trended_normalisation/trended_all.csv", usecols=["peak_id", "logFC"]).set_index("peak_id").logFC
g["atac_trended"] = [np.median(trd.reindex(str(s).split(";")).dropna()) for s in g.peak_ids]
ab = {}
for lab, c in CONTRASTS.items():
    d = pd.read_csv(DE / f"de_{c}.csv")[["gene_id", "logFC", "FDR", "AveExpr"]]
    ab[lab] = d.set_index("gene_id")
ave = ab["UTY"].AveExpr
g = g[g.ensembl.isin(ave.index)].copy()
g["AveExpr"] = ave.reindex(g.ensembl).to_numpy()
for lab in CONTRASTS:
    g[f"ab_{lab}"] = ab[lab].logFC.reindex(g.ensembl).to_numpy()
    g[f"abFDR_{lab}"] = ab[lab].FDR.reindex(g.ensembl).to_numpy()
down = g[g.rna_class == "down"].copy()

# matching (once; the same matched sets are used for every contrast)
rng = np.random.default_rng(SEED)
z = down[["rna_logFC", "AveExpr"]].apply(lambda s: (s - s.mean()) / s.std())
closed_idx = down.index[down.closed == 1].to_numpy(); pool = list(down.index[down.closed == 0])
pairs = []
for i in rng.permutation(closed_idx):
    dist = np.sqrt(((z.loc[pool] - z.loc[i]) ** 2).sum(axis=1))
    pick = list(dist.nsmallest(K).index)
    pairs += [(i, j) for j in pick]
    pool = [p for p in pool if p not in pick]
pairs = pd.DataFrame(pairs, columns=["closed_gene", "matched_gene"])
bal = pd.DataFrame({
    "set": ["closed", "matched controls", "all non-closed down"],
    "n": [len(closed_idx), pairs.matched_gene.nunique(), int((down.closed == 0).sum())],
    "median_baseline_logFC": [down.loc[closed_idx, "rna_logFC"].median(), down.loc[pairs.matched_gene, "rna_logFC"].median(),
                              down[down.closed == 0].rna_logFC.median()],
    "median_AveExpr": [down.loc[closed_idx, "AveExpr"].median(), down.loc[pairs.matched_gene, "AveExpr"].median(),
                       down[down.closed == 0].AveExpr.median()]})
bal.to_csv(OUT / "matching_balance.csv", index=False)

rows, gw, pcr, per_gene = [], [], [], []
for lab in CONTRASTS:
    y = f"ab_{lab}"
    d = down.rename(columns={y: "restore"})
    fit = smf.ols("restore ~ closed + rna_logFC + AveExpr", data=d).fit(cov_type="HC3")
    diffs = [d.loc[i, "restore"] - d.loc[pairs[pairs.closed_gene == i].matched_gene, "restore"].mean() for i in closed_idx]
    w = stats.wilcoxon(diffs)
    cl, ot = d[d.closed == 1], d[d.closed == 0]
    mt = d.loc[pairs.matched_gene]
    rows.append(dict(contrast=lab, n_closed=len(cl), n_not_closed=len(ot), n_matched=len(mt),
                     median_restore_closed=cl.restore.median(), median_restore_not_closed=ot.restore.median(),
                     median_restore_matched=mt.restore.median(),
                     pct_towards_closed=100 * (cl.restore > 0).mean(), pct_towards_matched=100 * (mt.restore > 0).mean(),
                     pct_towards_FDR05_closed=100 * ((cl.restore > 0) & (cl[f"abFDR_{lab}"] < .05)).mean(),
                     pct_towards_FDR05_matched=100 * ((mt.restore > 0) & (mt[f"abFDR_{lab}"] < .05)).mean(),
                     ols_closed_coef=fit.params["closed"], ols_closed_CI_low=fit.conf_int().loc["closed", 0],
                     ols_closed_CI_high=fit.conf_int().loc["closed", 1], ols_closed_P=fit.pvalues["closed"],
                     ols_baseline_coef=fit.params["rna_logFC"],
                     matched_median_diff=float(np.median(diffs)), matched_wilcoxon_P=w.pvalue,
                     median_abs_closed=cl.restore.abs().median(), median_abs_matched=mt.restore.abs().median(),
                     abs_matched_wilcoxon_P=stats.wilcoxon([abs(d.loc[i, "restore"]) - d.loc[pairs[pairs.closed_gene == i].matched_gene, "restore"].abs().mean()
                                                           for i in closed_idx]).pvalue))
    f2 = smf.ols(f"Q('{y}') ~ rna_logFC + atac_median_logFC + AveExpr", data=g).fit(cov_type="HC3")
    # sensitivity: promoter ATAC under trended (loess) normalisation, which removes the abundance-dependent trend
    f2t = smf.ols(f"Q('{y}') ~ rna_logFC + atac_trended + AveExpr", data=g).fit(cov_type="HC3")
    gw.append(dict(contrast=lab, n=int(f2.nobs), atac_coef=f2.params["atac_median_logFC"],
                   atac_CI_low=f2.conf_int().loc["atac_median_logFC", 0], atac_CI_high=f2.conf_int().loc["atac_median_logFC", 1],
                   atac_P=f2.pvalues["atac_median_logFC"], rna_coef=f2.params["rna_logFC"], rna_P=f2.pvalues["rna_logFC"],
                   atac_trended_coef=f2t.params["atac_trended"], atac_trended_CI_low=f2t.conf_int().loc["atac_trended", 0],
                   atac_trended_CI_high=f2t.conf_int().loc["atac_trended", 1], atac_trended_P=f2t.pvalues["atac_trended"]))
    f3 = smf.ols("restore ~ closed * ES_H3K27me3 + rna_logFC + AveExpr", data=d).fit(cov_type="HC3")
    for term in ["closed", "ES_H3K27me3", "closed:ES_H3K27me3"]:
        pcr.append(dict(contrast=lab, term=term, coef=f3.params[term], CI_low=f3.conf_int().loc[term, 0],
                        CI_high=f3.conf_int().loc[term, 1], P=f3.pvalues[term]))
    for grp_, sub in d.groupby(["closed", "ES_H3K27me3"]):
        pcr.append(dict(contrast=lab, term=f"median_restore closed={grp_[0]} ES_H3K27me3={grp_[1]}",
                        coef=sub.restore.median(), n=len(sub)))
res = pd.DataFrame(rows)
res["ols_closed_BH"] = stats.false_discovery_control(res.ols_closed_P)
res["matched_BH"] = stats.false_discovery_control(res.matched_wilcoxon_P)
res["abs_matched_BH"] = stats.false_discovery_control(res.abs_matched_wilcoxon_P)
res.to_csv(OUT / "restoration_closed_vs_not.csv", index=False)
gw = pd.DataFrame(gw); gw["atac_BH"] = stats.false_discovery_control(gw.atac_P)
gw["atac_trended_BH"] = stats.false_discovery_control(gw.atac_trended_P); gw.to_csv(OUT / "genomewide_atac_term.csv", index=False)
pcr = pd.DataFrame(pcr)
m = pcr.P.notna()
pcr.loc[m, "BH"] = stats.false_discovery_control(pcr.loc[m, "P"]); pcr.to_csv(OUT / "polycomb_interaction.csv", index=False)
keep = ["entrez_id", "ensembl", "gene_name", "chrom", "rna_logFC", "AveExpr", "closed", "ES_H3K27me3", "PRC2", "lost_010",
        "atac_median_logFC"] + [c for c in down.columns if c.startswith("ab")]
down[keep].to_csv(OUT / "down_genes_addback.csv", index=False)
pairs.assign(closed_name=down.loc[pairs.closed_gene, "gene_name"].to_numpy(),
             matched_name=down.loc[pairs.matched_gene, "gene_name"].to_numpy()).to_csv(OUT / "matched_pairs.csv", index=False)
json.dump(dict(DE_down_in_addback=len(down), closed=int(down.closed.sum()), K=K, seed=SEED), open(OUT / "summary.json", "w"), indent=2)
pd.set_option("display.width", 250)
print(bal.round(3).to_string(index=False)); print(res.round(4).T.to_string()); print(gw.round(4).to_string(index=False))
print(pcr.round(4).to_string(index=False))
