"""Task 4: CollecTRI targets of GRHL2 and ZEB1 vs DU145 DE results (plus context numbers used in
REPORT.md for the GRHL2 co-expression module).
Run:  python 03_collectri_targets.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

OUT = os.path.dirname(os.path.abspath(__file__))
BASE = "A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE"
net = pd.read_csv(f"{BASE}/resources/collectri_network.csv")
eff = pd.read_csv(f"{BASE}/outputs/tables/DU145_LOY_vs_sgCTRL_primary_effects.csv",
                  usecols=["gene_id", "gene_name", "CTRL_mean_CPM", "LOY_mean_CPM",
                           "mean_difference"])
de = pd.read_csv(f"{BASE}/outputs/tables/DE_paired_voom_primary.csv",
                 usecols=["gene_id", "logFC", "FDR", "thesis_DE"])
du = eff.merge(de, on="gene_id").drop_duplicates("gene_name").set_index("gene_name")

rows, summ = [], []
for tf in ["GRHL2", "ZEB1", "NKX2-1"]:
    t = net[net.source == tf].drop_duplicates("target").copy()
    t = t.join(du, on="target")
    t["tested_in_DU145"] = t.logFC.notna()
    t["DU145_FDR05"] = t.FDR < 0.05
    # For a TF that goes DOWN (GRHL2, DU145 logFC -3.95) an activating target (mor +1) is
    # expected DOWN and a repressed target (mor -1) expected UP; for ZEB1 (weakly UP,
    # logFC +0.33, FDR 0.48) the expectation is reversed.
    tf_dir = np.sign(du.loc[tf, "logFC"])
    t["expected_sign"] = tf_dir * np.sign(t.mor)
    t["agrees_with_mor"] = np.where(t.tested_in_DU145,
                                    np.sign(t.logFC) == t.expected_sign, np.nan)
    t.insert(0, "TF", tf)
    rows.append(t)
    t["self_edge"] = t.target == tf
    tt = t[t.tested_in_DU145 & ~t.self_edge]  # summary excludes the TF's own self-edge
    sig = tt[tt.DU145_FDR05]
    # signed-target score: mean of mor*logFC over tested targets, with a one-sample Wilcoxon
    sc = (np.sign(tt.mor) * tt.logFC)
    w = stats.wilcoxon(sc) if len(sc) >= 5 else None
    summ.append(dict(TF=tf, note="self-edge excluded from counts/score", TF_logFC=du.loc[tf, "logFC"], TF_FDR=du.loc[tf, "FDR"],
                     n_targets=len(t), n_activating=int((t.mor > 0).sum()),
                     n_repressive=int((t.mor < 0).sum()), n_tested=len(tt),
                     n_FDR05=len(sig), n_FDR05_agree=int(sig.agrees_with_mor.sum()),
                     n_tested_agree=int(tt.agrees_with_mor.sum()),
                     mean_mor_x_logFC=sc.mean(), wilcoxon_P=w.pvalue if w else np.nan,
                     FDR05_targets=";".join(f"{r.target}({r.logFC:+.2f},mor{int(r.mor):+d})"
                                            for r in sig.itertuples())))
T = pd.concat(rows)
T.to_csv(f"{OUT}/task4_collectri_GRHL2_ZEB1_targets.csv", index=False)
S = pd.DataFrame(summ)
S.to_csv(f"{OUT}/task4_collectri_summary.csv", index=False)
pd.set_option("display.width", 250)
print(S.T.to_string())
print(T[["TF", "target", "mor", "CTRL_mean_CPM", "logFC", "FDR", "agrees_with_mor"]]
      .round(3).to_string())

# context: DU145 effects on the canonical GRHL2-module / epithelial core genes
core = ("GRHL2 RAB25 ESRP1 ESRP2 CDH1 EPCAM CLDN4 CLDN7 OVOL1 OVOL2 SPINT1 PRSS8 KDF1 IRF6 CRB3 "
        "MARVELD3 ST14 TACSTD2 KRT8 KRT18 KRT19 ZEB1 ZEB2 VIM CDH2 FN1").split()
C = du.reindex(core)[["CTRL_mean_CPM", "mean_difference", "logFC", "FDR"]]
C.to_csv(f"{OUT}/task4_context_DU145_epithelial_core.csv")
print(C.round(3).to_string())
