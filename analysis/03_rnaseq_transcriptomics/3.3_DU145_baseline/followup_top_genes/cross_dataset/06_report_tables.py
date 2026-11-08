"""Builds the markdown tables embedded in REPORT.md from the task CSVs (writes report_tables.md)."""
import os
import numpy as np
import pandas as pd

OUT = os.path.dirname(os.path.abspath(__file__))


def f(x, d=2, sign=True):
    if pd.isna(x):
        return "–"
    return f"{x:+.{d}f}" if sign else f"{x:.{d}f}"


def fdr(x):
    if pd.isna(x):
        return "–"
    return f"{x:.1e}" if x < 0.001 else f"{x:.3f}"


def cpm(x):
    if pd.isna(x):
        return "–"
    return f"{x:.1f}" if x >= 1 else f"{x:.2f}"


def md(df):
    cols = list(df.columns)
    s = "| " + " | ".join(cols) + " |\n|" + "|".join(["---"] * len(cols)) + "|\n"
    for _, r in df.iterrows():
        s += "| " + " | ".join(str(r[c]) for c in cols) + " |\n"
    return s


parts = []
# ---- Task 1
t = pd.read_csv(f"{OUT}/task1_inhouse_lines_focus_genes.csv")
order = {"focus_down": 0, "focus_up": 1, "EMT": 2, "IFN": 3}
t["o"] = t.panel.map(order)
t = t.sort_values(["o", "DU145_mean_difference"])
T1 = pd.DataFrame({
    "gene": t.gene_name, "panel": t.panel,
    "DU145 CTRL CPM": t.DU145_CTRL_mean_CPM.map(cpm),
    "DU145 md (FDR)": [f"{f(a)} ({fdr(b)})" for a, b in zip(t.DU145_mean_difference, t.DU145_FDR)],
    "HT1080 CPM": t.HT1080_allgene_mean_CPM.map(cpm),
    "HT1080 md": t.HT1080_allgene_mean_difference.map(f),
    "HT1080 batchDE logFC (FDR)": [f"{f(a)} ({fdr(b)})" if pd.notna(a) else "not tested"
                                   for a, b in zip(t.HT1080_batchDE_logFC, t.HT1080_batchDE_FDR)],
    "A549 CPM": t.A549_allgene_mean_CPM.map(cpm),
    "A549 md": t.A549_allgene_mean_difference.map(f),
    "A549 unpaired logFC (FDR)": [f"{f(a)} ({fdr(b)})" if pd.notna(a) else "not tested"
                                  for a, b in zip(t.A549_unpairedDE_logFC, t.A549_unpairedDE_FDR)],
})
parts.append("T1\n" + md(T1))
g = pd.read_csv(f"{OUT}/task1_genomewide_concordance.csv")
g["sign_agreement"] = g.sign_agreement.map(lambda x: f"{x:.2f}")
g["spearman"] = g.spearman.map(lambda x: f"{x:+.3f}")
parts.append("T1b\n" + md(g))
# ---- Task 2
w = pd.read_csv(f"{OUT}/task2_public_datasets_wide.csv")
parts.append("T2\n" + md(w.fillna("nt")))
s = pd.read_csv(f"{OUT}/task2_public_datasets_summary.csv").drop(columns=["description"])
parts.append("T2s\n" + md(s))
m = pd.read_csv(f"{OUT}/task2b_GRHL2_module_and_DU145down_across_datasets.csv")
m = m[m.gene_set == "GRHL2_top200_module"]
M = pd.DataFrame({"dataset": m.dataset, "n module genes": m.n_set_tested,
                  "median logFC module / rest": [f"{a:+.3f} / {b:+.3f}" for a, b in zip(m.median_set, m.median_rest)],
                  "AUC (module higher)": m.AUC_set_higher.map(lambda x: f"{x:.3f}"),
                  "Wilcoxon P": m.P.map(lambda x: f"{x:.1e}"),
                  "random-module null: frac AUC <= obs": m.null_frac_AUC_le_obs.map(lambda x: f"{x:.3f}")})
parts.append("T2m\n" + md(M))
m2 = pd.read_csv(f"{OUT}/task2b_GRHL2_module_and_DU145down_across_datasets.csv")
m2 = m2[m2.gene_set == "DU145_FDR05_down"]
M2 = pd.DataFrame({"dataset": m2.dataset, "n": m2.n_set_tested,
                   "median set / rest": [f"{a:+.3f} / {b:+.3f}" for a, b in zip(m2.median_set, m2.median_rest)],
                   "AUC": m2.AUC_set_higher.map(lambda x: f"{x:.3f}"),
                   "P": m2.P.map(lambda x: f"{x:.2g}")})
parts.append("T2d\n" + md(M2))
# ---- Task 3
a = pd.read_csv(f"{OUT}/task3a_depmap_Yscore_association_focus_genes.csv")
a["o"] = a.panel.map(order)
a = a.sort_values(["o"])
lin_cols = [c for c in a.columns if c.startswith("rho_") and c != "rho_Prostate"]
A = pd.DataFrame({
    "gene": a.gene, "panel": a.panel, "DU145 md": a.DU145_mean_difference.map(f),
    "male median log2(TPM+1)": a.male_median_log2TPM1.map(lambda x: f"{x:.2f}"),
    "Spearman vs Y (all male, n=905)": [f"{r:+.2f} ({p:.1g})" for r, p in zip(a.spearman_all_male, a.spearman_all_male_P)],
    "OLS Y coef ± SE (P) | lineage": [f"{c:+.3f} ± {s_:.3f} ({p:.2g})" for c, s_, p in zip(a.OLS_Y_coef, a.OLS_Y_SE, a.OLS_Y_P)],
    "OLS Y-low coef (P)": [f"{c:+.2f} ({p:.2g})" for c, p in zip(a.OLS_Ylow_coef, a.OLS_Ylow_P)],
    "sign as expected": a.OLS_Y_coef_matches_expected.map({True: "yes", False: "no"}),
    "within-lineage rho range (12 lineages ≥30)": [f"{np.nanmin(r):+.2f} to {np.nanmax(r):+.2f}" for r in a[lin_cols].values],
    "Prostate rho (n=13, P)": [f"{r:+.2f} ({p:.2g})" for r, p in zip(a.rho_Prostate, a.P_Prostate)],
})
parts.append("T3a\n" + md(A))
open(f"{OUT}/report_tables.md", "w", encoding="utf-8").write("\n\n".join(parts))
print("wrote report_tables.md")
