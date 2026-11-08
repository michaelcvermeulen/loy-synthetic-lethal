"""Task 3: DepMap 26Q1 association of DU145 focus genes with a Y-expression score, the GRHL2
co-expression neighbourhood test, and GRHL2-ZEB1 context.

Inputs (read-only, hash-verified release):
  H:/PRISM/data/CCLE/26Q1/raw/depmap_public_26q1/expression/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
  H:/PRISM/data/CCLE/26Q1/raw/depmap_public_26q1/metadata/Model.csv
Run:  python 02_depmap_Y_association.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

OUT = os.path.dirname(os.path.abspath(__file__))
DM = "H:/PRISM/data/CCLE/26Q1/raw/depmap_public_26q1"
DU = ("A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE/"
      "outputs/tables")
RNG = np.random.default_rng(20260923)

FOCUS_DOWN = ("GRHL2 TACSTD2 AP1M2 CLDN7 EPS8L1 SCEL MAL2 TMC4 ST14 OVOL1 S100A14 LCN2 CXCL1 "
              "BST2 LTB SFTPA1 NAPSA HOXB9 PKIA CCR9 PXDN").split()
FOCUS_UP = "ZEB1 TAF9 F13A1 COL4A2 CREB3L1 RSPO3 ADAM12 CNTN1 TMEM45A".split()
EMT = "CDH1 CDH2 VIM ZEB1 ZEB2 SNAI1 SNAI2 TWIST1 ESRP1 ESRP2 OVOL2".split()
IFN = "STAT1 IRF1 IFI6 IFI44 ISG15 CXCL10".split()
GENES = list(dict.fromkeys(FOCUS_DOWN + FOCUS_UP + EMT + IFN))
PANEL = {g: ("focus_down" if g in FOCUS_DOWN else "focus_up" if g in FOCUS_UP
             else "EMT" if g in EMT else "IFN") for g in GENES}
YG = ["DDX3Y", "KDM5D", "UTY", "USP9Y", "EIF1AY", "ZFY", "RPS4Y1"]

# ------------------------------------------------------------------ load
expr = pd.read_csv(f"{DM}/expression/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
expr = expr[expr.IsDefaultEntryForModel.astype(str) == "Yes"].copy()
assert not expr.ModelID.duplicated().any()
gene_cols = [c for c in expr.columns if " (" in c]
X = expr[gene_cols].astype(np.float32)
X.columns = [c.split(" (")[0] for c in gene_cols]
X = X.loc[:, ~X.columns.duplicated()]
X.index = expr.ModelID.values
meta = pd.read_csv(f"{DM}/metadata/Model.csv").set_index("ModelID")
meta = meta.reindex(X.index)
print("default profiles:", X.shape)

# ------------------------------------------------------------------ Y score (male lines)
male = meta.Sex.eq("Male").values
Xm = X[male]
mm = meta[male].copy()
mm["Y_score"] = Xm[YG].mean(axis=1)
mm["lineage"] = mm.OncotreeLineage.fillna("Unknown")
lin_n = mm.lineage.value_counts()
print("male default lines:", male.sum())
# distribution / DU145 context
ctx = mm[["StrippedCellLineName", "lineage", "Y_score"]].copy()
for gname in YG + ["GRHL2", "ZEB1", "CDH1", "VIM"]:
    ctx[gname] = Xm[gname]
ctx.sort_values("Y_score").to_csv(f"{OUT}/task3_depmap_male_Yscore_by_line.csv")
q = mm.Y_score.quantile([0.05, 0.1, 0.25, 0.5, 0.75, 0.9]).round(3)
female_Y = X[meta.Sex.eq("Female").values][YG].mean(axis=1)
Y_LOW = float(female_Y.quantile(0.99))  # Y-low = at/below the female 99th percentile
mm["Y_low"] = (mm.Y_score <= Y_LOW).astype(int)
dist = dict(n_male=int(male.sum()), female_Y_p99=Y_LOW,
            male_Y_low_n=int(mm.Y_low.sum()), male_Y_low_frac=float(mm.Y_low.mean()),
            **{f"male_Y_q{int(k*100)}": v for k, v in q.items()},
            female_Y_median=float(female_Y.median()))
prost = ctx[ctx.lineage == "Prostate"].sort_values("Y_score")
print(prost.round(2).to_string())
prost.to_csv(f"{OUT}/task3_depmap_prostate_lines.csv")

# ------------------------------------------------------------------ (a) per-gene associations
lineages_30 = [l for l, n in lin_n.items() if n >= 30]
rows = []
for gname in GENES:
    if gname not in Xm.columns:
        rows.append(dict(gene=gname, panel=PANEL[gname], in_depmap=False))
        continue
    y = Xm[gname].values.astype(float)
    r = dict(gene=gname, panel=PANEL[gname], in_depmap=True,
             male_median_log2TPM1=float(np.median(y)),
             male_frac_expr_TPMge1=float((y >= 1).mean()))
    rho, p = stats.spearmanr(mm.Y_score, y)
    r.update(spearman_all_male=rho, spearman_all_male_P=p, n_all_male=len(y))
    d = pd.DataFrame(dict(g=y, Y=mm.Y_score.values, Ylow=mm.Y_low.values,
                          lin=mm.lineage.values))
    fit = smf.ols("g ~ Y + C(lin)", data=d).fit()
    r.update(OLS_Y_coef=fit.params["Y"], OLS_Y_SE=fit.bse["Y"], OLS_Y_t=fit.tvalues["Y"],
             OLS_Y_P=fit.pvalues["Y"])
    fit2 = smf.ols("g ~ Ylow + C(lin)", data=d).fit()
    r.update(OLS_Ylow_coef=fit2.params["Ylow"], OLS_Ylow_P=fit2.pvalues["Ylow"])
    for lin in lineages_30 + ["Prostate"]:
        s = d[d.lin == lin]
        if len(s) >= 5 and s.g.std() > 0:
            rr, pp = stats.spearmanr(s.Y, s.g)
        else:
            rr, pp = np.nan, np.nan
        r[f"rho_{lin}"] = rr
        r[f"P_{lin}"] = pp
        r[f"n_{lin}"] = len(s)
    rows.append(r)
A = pd.DataFrame(rows)
# direction expected under "DU145 effect is a generic Y-loss effect":
#   DU145-down gene -> positive association with Y score; DU145-up -> negative
du = pd.read_csv(f"{DU}/DU145_LOY_vs_sgCTRL_primary_effects.csv",
                 usecols=["gene_id", "gene_name", "mean_difference"])
de = pd.read_csv(f"{DU}/DE_paired_voom_primary.csv", usecols=["gene_id", "logFC", "FDR"])
du = du.merge(de, on="gene_id", how="left")
dmd = du.drop_duplicates("gene_name").set_index("gene_name").mean_difference
A["DU145_mean_difference"] = A.gene.map(dmd)
A["expected_sign_Y_coef"] = -np.sign(A.DU145_mean_difference)
A["OLS_Y_coef_matches_expected"] = np.sign(A.OLS_Y_coef) == A.expected_sign_Y_coef
A["OLS_Y_FDR_panel"] = np.nan
ok = A.OLS_Y_P.notna()
from statsmodels.stats.multitest import multipletests
A.loc[ok, "OLS_Y_FDR_panel"] = multipletests(A.loc[ok, "OLS_Y_P"], method="fdr_bh")[1]
A.to_csv(f"{OUT}/task3a_depmap_Yscore_association_focus_genes.csv", index=False)

# genome-wide partial association (Frisch-Waugh: residualise on lineage) for set tests
L = pd.get_dummies(mm.lineage, drop_first=False).values.astype(float)
Q, _ = np.linalg.qr(L)


def resid(M):
    return M - Q @ (Q.T @ M)


Yr = resid(mm.Y_score.values[:, None].astype(float))[:, 0]
Gm = Xm.values.astype(float)
Gr = resid(Gm)
df_res = len(Yr) - Q.shape[1] - 1
beta = (Gr.T @ Yr) / (Yr @ Yr)
res = Gr - np.outer(Yr, beta)
se = np.sqrt((res ** 2).sum(0) / df_res / (Yr @ Yr))
tval = np.where(se > 0, beta / se, np.nan)
GW = pd.DataFrame(dict(gene=Xm.columns, Y_coef=beta, Y_t=tval,
                       male_mean=Gm.mean(0)))
GW = GW[GW.male_mean > 0.5]  # drop essentially unexpressed genes
GW.to_csv(f"{OUT}/task3a_depmap_genomewide_Y_partial.csv", index=False)
# sanity: the check on the Y genes themselves (should be huge t)
print(GW[GW.gene.isin(YG)].to_string())
GWnoY = GW[~GW.gene.isin(YG + ["NLGN4Y", "TXLNGY", "PRKY", "TMSB4Y", "PCDH11Y", "TBL1Y",
                               "AMELY", "SRY", "TSPY1", "RBMY1A1", "VCY", "BPY2", "CDY1",
                               "DAZ1", "HSFY1", "XKRY", "PRY", "NLGN4Y"])]
uni = GWnoY.merge(du.drop_duplicates("gene_name"), left_on="gene", right_on="gene_name")
settests = []
for label, mask in [("DU145_FDR05_down", (uni.FDR < 0.05) & (uni.logFC < 0)),
                    ("DU145_FDR05_up", (uni.FDR < 0.05) & (uni.logFC > 0))]:
    a, b = uni.loc[mask, "Y_t"].dropna(), uni.loc[~mask, "Y_t"].dropna()
    U, p = stats.mannwhitneyu(a, b)
    settests.append(dict(set=label, n=len(a), median_Y_t_set=a.median(),
                         median_Y_t_rest=b.median(), AUC=U / (len(a) * len(b)), P=p))
rho_gw, p_gw = stats.spearmanr(uni.mean_difference, uni.Y_t)
settests.append(dict(set="all genes: spearman(DU145 mean_difference, DepMap Y_t)",
                     n=len(uni), AUC=np.nan, P=p_gw, median_Y_t_set=rho_gw))
ST = pd.DataFrame(settests)
ST.to_csv(f"{OUT}/task3a_depmap_Y_partial_settests.csv", index=False)
print(ST.to_string())

# ------------------------------------------------------------------ (b) GRHL2 neighbourhood
Xa = X.values.astype(np.float64)
keep = Xa.std(0) > 0
Xa = Xa[:, keep]
genes_all = X.columns[keep]
Z = (Xa - Xa.mean(0)) / Xa.std(0)
n = Z.shape[0]


def cor_with(gene):
    v = Z[:, genes_all.get_loc(gene)]
    return pd.Series(Z.T @ v / n, index=genes_all)


rg = cor_with("GRHL2").drop("GRHL2")
rz = cor_with("ZEB1").drop("ZEB1")
top_g = rg.sort_values(ascending=False).head(200)
top_z = rz.sort_values(ascending=False).head(200)
pd.DataFrame(dict(gene=top_g.index, r_with_GRHL2=top_g.values)).to_csv(
    f"{OUT}/task3b_GRHL2_top200_neighbourhood.csv", index=False)

duu = du.drop_duplicates("gene_name").set_index("gene_name")
univ = duu.index.intersection(genes_all)
duu = duu.loc[univ]
sig_down = set(duu.index[(duu.FDR < 0.05) & (duu.logFC < 0)])
sig_up = set(duu.index[(duu.FDR < 0.05) & (duu.logFC > 0)])


def nb_test(nb, label, direction_set, dir_label):
    inset = duu.index.isin(nb)
    a, b = duu.mean_difference[inset], duu.mean_difference[~inset]
    U, p = stats.mannwhitneyu(a, b)
    k = len(set(nb) & direction_set & set(univ))
    tab = [[k, inset.sum() - k], [len(direction_set) - k,
                                  len(univ) - inset.sum() - len(direction_set) + k]]
    OR, pf = stats.fisher_exact(tab, alternative="greater")
    return dict(neighbourhood=label, n_in_DU145_universe=int(inset.sum()),
                universe=len(univ), median_md_nb=a.median(), median_md_rest=b.median(),
                AUC_nb_higher=U / (len(a) * len(b)), wilcoxon_P=p,
                overlap_set=dir_label, overlap_k=k, set_size=len(direction_set),
                expected_k=inset.sum() * len(direction_set) / len(univ),
                fisher_OR=OR, fisher_P_greater=pf,
                overlap_genes=";".join(sorted(set(nb) & direction_set)))


NB = [nb_test(top_g.index, "GRHL2_top200", sig_down, "DU145_FDR05_down"),
      nb_test(top_z.index, "ZEB1_top200", sig_up, "DU145_FDR05_up"),
      nb_test(top_z.index, "ZEB1_top200_vs_down", sig_down, "DU145_FDR05_down")]
# calibration: top-200 neighbourhoods of 500 random expressed genes
cand = [g for g in genes_all if g in univ and X[g].mean() > 1 and g != "GRHL2"]
rand = RNG.choice(cand, 500, replace=False)
Rm = Z[:, [genes_all.get_loc(g) for g in rand]]
Cmat = (Z.T @ Rm) / n  # genes x 500
null_auc, null_k = [], []
md = duu.mean_difference
for j, g in enumerate(rand):
    s = pd.Series(Cmat[:, j], index=genes_all).drop(g).sort_values(ascending=False).head(200)
    inset = md.index.isin(s.index)
    a, b = md[inset], md[~inset]
    U, _ = stats.mannwhitneyu(a, b)
    null_auc.append(U / (len(a) * len(b)))
    null_k.append(len(set(s.index) & sig_down))
null_auc, null_k = np.array(null_auc), np.array(null_k)
NB[0]["null_frac_random_nb_AUC_le_obs"] = float((null_auc <= NB[0]["AUC_nb_higher"]).mean())
NB[0]["null_frac_random_nb_overlap_ge_obs"] = float((null_k >= NB[0]["overlap_k"]).mean())
NB[0]["null_AUC_median"] = float(np.median(null_auc))
NB[0]["null_overlap_median"] = float(np.median(null_k))
NB = pd.DataFrame(NB)
NB.to_csv(f"{OUT}/task3b_neighbourhood_enrichment.csv", index=False)
print(NB.T.to_string())
# which of the DU145 down genes are the strongest GRHL2 correlates
fd = duu.loc[duu.index.isin(sig_down)].copy()
fd["r_with_GRHL2"] = rg.reindex(fd.index)
fd["rank_in_GRHL2_corr"] = rg.rank(ascending=False).reindex(fd.index)
fd.sort_values("r_with_GRHL2", ascending=False).to_csv(
    f"{OUT}/task3b_DU145_down_genes_GRHL2_correlation.csv")

# ------------------------------------------------------------------ (c) GRHL2 vs ZEB1
cz = []
for label, mask in [("all_default_lines", np.ones(len(X), bool)),
                    ("male_lines", male),
                    ("female_lines", meta.Sex.eq("Female").values),
                    ("prostate_lines", meta.OncotreeLineage.eq("Prostate").values)]:
    a, b = X.loc[mask, "GRHL2"], X.loc[mask, "ZEB1"]
    rs, ps = stats.spearmanr(a, b)
    rp, pp = stats.pearsonr(a, b)
    cz.append(dict(subset=label, n=int(mask.sum()), spearman=rs, spearman_P=ps,
                   pearson=rp, pearson_P=pp))
# lineage-adjusted (all lines)
dd = pd.DataFrame(dict(g=X.GRHL2.values, z=X.ZEB1.values,
                       lin=meta.OncotreeLineage.fillna("Unknown").values))
f = smf.ols("g ~ z + C(lin)", data=dd).fit()
cz.append(dict(subset="all lines, OLS GRHL2 ~ ZEB1 + lineage", n=len(dd),
               pearson=f.params["z"], pearson_P=f.pvalues["z"]))
CZ = pd.DataFrame(cz)
CZ.to_csv(f"{OUT}/task3c_GRHL2_ZEB1_correlation.csv", index=False)
print(CZ.to_string())
pd.Series(dist).to_csv(f"{OUT}/task3_Yscore_distribution.csv", header=["value"])
print(pd.Series(dist))
print("lineages >=30 male lines:", {l: int(lin_n[l]) for l in lineages_30},
      "Prostate:", int(lin_n.get("Prostate", 0)))
