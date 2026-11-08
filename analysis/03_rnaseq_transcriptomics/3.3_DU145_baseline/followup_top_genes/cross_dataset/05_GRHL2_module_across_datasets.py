"""Module-level test: is the DepMap GRHL2 top-200 co-expression neighbourhood (task3b) shifted
down in each other dataset? Wilcoxon rank-sum of the dataset's Y-loss-polarity logFC for module
genes vs all other tested genes. Also the DU145 FDR<0.05-down genes as a set.
Symbol-level join (mouse: upper-cased symbols). Run after 01 and 02.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

OUT = os.path.dirname(os.path.abspath(__file__))
AN = "A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis"
DU = f"{AN}/LNCAP/SEPT_23_2026_REANALYZE/outputs/tables"
mod = set(pd.read_csv(f"{OUT}/task3b_GRHL2_top200_neighbourhood.csv").gene)
mod.add("GRHL2")
du = pd.read_csv(f"{DU}/DE_paired_voom_primary.csv", usecols=["gene_name", "logFC", "FDR"])
du_down = set(du.gene_name[(du.FDR < 0.05) & (du.logFC < 0)])

D = [  # key, path, symbol col, lfc col, polarity, upper-case symbols?
    ("DU145_LOY_primary", f"{DU}/DU145_LOY_vs_sgCTRL_primary_effects.csv", "gene_name",
     "mean_difference", 1, False),
    ("HT1080_LOY_primary", f"{AN}/HT1080/SEPT_11_2026_REANALYZE/outputs/tables/HT1080_LOY_vs_sgCTRL_primary_effects.csv",
     "gene_name", "mean_difference", 1, False),
    ("A549_LOY_primary", f"{AN}/A549/SEPT_18_2026_REANALYZE/outputs/tables/A549_LOY_vs_sgCTRL_primary_effects.csv",
     "gene_name", "mean_difference", 1, False),
    ("A549_clones_LOYvROY", f"{AN}/A549_clonal_LOY/DE/clonal_LOY_v_ROY_DESeq2.csv", "gene_name", "logFC", 1, False),
    ("A549_multiome_pb_LOYvROY", f"{AN}/A549_multiome_pseudobulk/DE/multiome_pb_LOY_v_ROY_DESeq2.csv", "gene_name", "logFC", 1, False),
    ("MB49_cellline_YnegvYpos", f"{AN}/GSE229233/DE/gse229233_cellline_Yneg_v_Ypos_DESeq2.csv", "gene_name", "log2FoldChange", 1, True),
    ("MB49_tumor_YnegvYpos", f"{AN}/GSE229233/DE/gse229233_tumor_Yneg_v_Ypos_DESeq2.csv", "gene_name", "log2FoldChange", 1, True),
    ("GSE126712_lymph_XOvXY", f"{AN}/GSE126712/DE/gse126712_XO_v_XY_DE.csv", "gene_name", "logFC", 1, False),
    ("GSE207112_troph_XOvXY", f"{AN}/GSE207112/DE/gse207112_XO_v_XY_DE_blocked.csv", "gene_name", "logFC", 1, False),
    ("GSE264742_NCC_XOvXY", f"{AN}/GSE264742/DE/gse264742_XO_v_XY_DE_blocked.csv", "gene_name", "logFC", 1, False),
    ("ERP156313_fetal_XOvXY", f"{AN}/ERP156313/DE/erp156313_consensus_XO_v_XY_DE.csv", "gene_name", "logFC", 1, False),
    ("PageLab_Fib_X0vXY", f"{AN}/PageLab_SCA/DE/trackA_X0_vs_XY_Fib.csv", "gene", "log2FoldChange", 1, False),
    ("PageLab_LCL_X0vXY", f"{AN}/PageLab_SCA/DE/trackA_X0_vs_XY_LCL.csv", "gene", "log2FoldChange", 1, False),
    ("GSE301297_undiff_UTYKOvWT", f"{AN}/GSE301297/DE/gse301297_SEES3undiff_UTYKO_v_WT_DE.csv", "gene_name", "logFC", 1, False),
    ("XYY_integrated_flipped", f"{AN}/integrated_XYY_XY/DE/integrated_consensus_XYY_v_XY_DE.csv", "gene_name", "logFC", -1, False),
]
# ---- null: top-200 DepMap co-expression neighbourhoods of 300 random expressed genes
DM = "H:/PRISM/data/CCLE/26Q1/raw/depmap_public_26q1"
ex = pd.read_csv(f"{DM}/expression/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
ex = ex[ex.IsDefaultEntryForModel.astype(str) == "Yes"]
gc = [c for c in ex.columns if " (" in c]
X = ex[gc].astype(float)
X.columns = [c.split(" (")[0] for c in gc]
X = X.loc[:, ~X.columns.duplicated()]
X = X.loc[:, X.std() > 0]
Z = ((X - X.mean()) / X.std(ddof=0)).values
genes = X.columns
rng = np.random.default_rng(7)
cand = [g for g in genes[(X.mean() > 1).values] if g != "GRHL2"]
seeds = rng.choice(cand, 300, replace=False)
Cm = Z.T @ Z[:, [genes.get_loc(g) for g in seeds]] / Z.shape[0]
NULL = []
for j, g in enumerate(seeds):
    s_ = pd.Series(Cm[:, j], index=genes).drop(g).sort_values(ascending=False).head(200)
    NULL.append(set(s_.index) | {g})
del X, Z, Cm, ex

rows = []
for key, path, sc, lc, pol, up in D:
    d = pd.read_csv(path)
    d["sym"] = d[sc].astype(str).str.replace(r"\.\d+$", "", regex=True)
    if up:
        d["sym"] = d.sym.str.upper()
    d = d.dropna(subset=[lc]).drop_duplicates("sym")
    x = pol * d[lc]
    for setname, S in [("GRHL2_top200_module", mod), ("DU145_FDR05_down", du_down)]:
        if key.startswith("DU145") and setname == "DU145_FDR05_down":
            continue
        m = d.sym.isin(S)
        a, b = x[m], x[~m]
        U, p = stats.mannwhitneyu(a, b)
        auc = U / (len(a) * len(b))
        r = dict(dataset=key, gene_set=setname, n_set_tested=int(m.sum()),
                 n_rest=int((~m).sum()), median_set=a.median(), median_rest=b.median(),
                 frac_set_below_zero=float((a < 0).mean()), AUC_set_higher=auc, P=p)
        if setname == "GRHL2_top200_module":
            na = []
            for S0 in NULL:
                m0 = d.sym.isin(S0)
                if m0.sum() < 20:
                    continue
                U0, _ = stats.mannwhitneyu(x[m0], x[~m0])
                na.append(U0 / (m0.sum() * (~m0).sum()))
            na = np.array(na)
            r.update(null_n=len(na), null_frac_AUC_le_obs=float((na <= auc).mean()),
                     null_frac_absdev_ge_obs=float((abs(na - 0.5) >= abs(auc - 0.5)).mean()),
                     null_AUC_sd=float(na.std()))
        rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/task2b_GRHL2_module_and_DU145down_across_datasets.csv", index=False)
pd.set_option("display.width", 250)
print(R.round(4).to_string())
