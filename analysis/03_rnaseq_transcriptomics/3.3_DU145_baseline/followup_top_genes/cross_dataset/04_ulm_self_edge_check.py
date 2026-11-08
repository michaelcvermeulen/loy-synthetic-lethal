"""Sensitivity: does the DU145 CollecTRI ULM score for GRHL2 / ZEB1 / NKX2-1 survive removal of
self-edges and of the TF-TF edge GRHL2<->ZEB1? ULM = t of slope in lm(stat ~ mor-weight) over
the gene universe (decoupleR definition), stat = DU145 mean_difference (one row per symbol).
Run: python 04_ulm_self_edge_check.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

OUT = os.path.dirname(os.path.abspath(__file__))
BASE = "A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE"
net = pd.read_csv(f"{BASE}/resources/collectri_network.csv")
du = pd.read_csv(f"{BASE}/outputs/tables/DU145_LOY_vs_sgCTRL_primary_effects.csv",
                 usecols=["gene_name", "mean_difference", "CTRL_mean_CPM"])
du = du.sort_values("CTRL_mean_CPM", ascending=False).drop_duplicates("gene_name")
stat = du.set_index("gene_name").mean_difference


def ulm(tf, drop=()):
    t = net[(net.source == tf) & (~net.target.isin(drop))].drop_duplicates("target")
    w = pd.Series(0.0, index=stat.index)
    tt = t[t.target.isin(stat.index)]
    w[tt.target.values] = tt.mor.values
    r = stats.linregress(w.values, stat.values)
    tval = r.slope / r.stderr
    p = 2 * stats.t.sf(abs(tval), len(stat) - 2)
    return dict(TF=tf, dropped=";".join(drop) or "none", n_targets_in_universe=len(tt),
                ULM_t=tval, P=p)


rows = []
for tf, drops in [("GRHL2", [(), ("GRHL2",), ("GRHL2", "ZEB1")]),
                  ("ZEB1", [(), ("ZEB1",), ("GRHL2",)]),
                  ("NKX2-1", [(), ("NKX2-1",), ("NKX2-1", "SFTPA1")])]:
    for d in drops:
        rows.append(ulm(tf, d))
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/task4_ulm_self_edge_sensitivity.csv", index=False)
print(R.to_string())
