"""Tasks 1 + 2: look up DU145 focus genes in the HT1080 / A549 thesis reanalyses and in the
locally available public / clonal LOY and sex-chromosome datasets.

Read-only on every input. Writes CSVs next to this script.
Run:  python 01_inhouse_and_public_lookup.py
"""
import os
import re
import numpy as np
import pandas as pd

OUT = os.path.dirname(os.path.abspath(__file__))
AN = "A:/Dropbox/LOY_paralogs/projects/RNAseq/analysis"
DU = f"{AN}/LNCAP/SEPT_23_2026_REANALYZE/outputs/tables"
HT = f"{AN}/HT1080/SEPT_11_2026_REANALYZE/outputs/tables"
A5 = f"{AN}/A549/SEPT_18_2026_REANALYZE/outputs/tables"

FOCUS_DOWN = ("GRHL2 TACSTD2 AP1M2 CLDN7 EPS8L1 SCEL MAL2 TMC4 ST14 OVOL1 S100A14 LCN2 CXCL1 "
              "BST2 LTB SFTPA1 NAPSA HOXB9 PKIA CCR9 PXDN").split()
FOCUS_UP = "ZEB1 TAF9 F13A1 COL4A2 CREB3L1 RSPO3 ADAM12 CNTN1 TMEM45A".split()
EMT = "CDH1 CDH2 VIM ZEB1 ZEB2 SNAI1 SNAI2 TWIST1 ESRP1 ESRP2 OVOL2".split()
IFN = "STAT1 IRF1 IFI6 IFI44 ISG15 CXCL10".split()
GENES = list(dict.fromkeys(FOCUS_DOWN + FOCUS_UP + EMT + IFN))
PANEL = {g: ("focus_down" if g in FOCUS_DOWN else "focus_up" if g in FOCUS_UP
             else "EMT" if g in EMT else "IFN") for g in GENES}


def stable(x):
    return x.astype(str).str.replace(r"\.\d+$", "", regex=True)


# ---------------------------------------------------------------- DU145 reference
de = pd.read_csv(f"{DU}/DE_paired_voom_primary.csv")
eff = pd.read_csv(f"{DU}/DU145_LOY_vs_sgCTRL_primary_effects.csv")
du = eff[["gene_id", "gene_name", "CTRL_mean_CPM", "LOY_mean_CPM", "mean_difference",
          "expressed"]].merge(de[["gene_id", "logFC", "FDR", "thesis_DE"]], on="gene_id",
                              how="left")
ref = du[du.gene_name.isin(GENES)].copy()
dup = ref.gene_name.duplicated(keep=False)
if dup.any():  # keep the representative with highest expression
    ref = ref.sort_values("CTRL_mean_CPM", ascending=False).drop_duplicates("gene_name")
missing = sorted(set(GENES) - set(ref.gene_name))
print("DU145 genes missing from primary effects table:", missing)
# pull missing genes' IDs from the all-gene effects file if present
if missing:
    allg = pd.read_csv(f"{DU}/DU145_LOY_vs_sgCTRL_all_gene_effects.csv")
    extra = allg[allg.gene_name.isin(missing)]
    keep = [c for c in ["gene_id", "gene_name", "CTRL_mean_CPM", "LOY_mean_CPM",
                        "mean_difference", "expressed"] if c in extra.columns]
    ref = pd.concat([ref, extra[keep].assign(logFC=np.nan, FDR=np.nan, thesis_DE=False)])
ref["panel"] = ref.gene_name.map(PANEL)
ref["stable_id"] = stable(ref.gene_id)
ref = ref.rename(columns={c: f"DU145_{c}" for c in ref.columns
                          if c not in ("gene_id", "gene_name", "panel", "stable_id")})

# ---------------------------------------------------------------- Task 1: in-house lines
ht = pd.read_csv(f"{HT}/HT1080_LOY_vs_sgCTRL_primary_effects.csv",
                 usecols=["gene_id", "mean_CPM", "expressed", "mean_difference",
                          "AUG_2025", "MAY_2026", "matches_agreeing_with_mean",
                          "LOO_all_same_direction"])
ht = ht.rename(columns={c: f"HT1080_{c}" for c in ht.columns if c != "gene_id"})
hb = pd.read_csv(f"{HT}/batchDE_batch_condition__combined.csv",
                 usecols=["gene_id", "logFC", "PValue", "FDR"])
hb = hb.rename(columns={"logFC": "HT1080_batchDE_logFC", "PValue": "HT1080_batchDE_P",
                        "FDR": "HT1080_batchDE_FDR"})
hp = pd.read_csv(f"{HT}/batchDE_passage_period__combined.csv",
                 usecols=["gene_id", "logFC", "FDR"])
hp = hp.rename(columns={"logFC": "HT1080_passageDE_logFC", "FDR": "HT1080_passageDE_FDR"})
a5 = pd.read_csv(f"{A5}/A549_LOY_vs_sgCTRL_primary_effects.csv",
                 usecols=["gene_id", "mean_CPM", "CTRL_mean_CPM", "LOY_mean_CPM", "expressed",
                          "mean_difference", "LOO_all_same_direction"])
a5 = a5.rename(columns={c: f"A549_{c}" for c in a5.columns if c != "gene_id"})
au = pd.read_csv(f"{A5}/DE_unpaired_voom_exploratory.csv",
                 usecols=["gene_id", "logFC", "PValue", "FDR"])
au = au.rename(columns={"logFC": "A549_unpairedDE_logFC", "PValue": "A549_unpairedDE_P",
                        "FDR": "A549_unpairedDE_FDR"})
t1 = ref.copy()
for tab in (ht, hb, hp, a5, au):
    t1 = t1.merge(tab, on="gene_id", how="left")
t1["HT1080_in_primary_table"] = t1["HT1080_mean_difference"].notna()
t1["A549_in_primary_table"] = t1["A549_mean_difference"].notna()
# genes filtered out of a primary table: take CPM / descriptive difference from all-gene tables
for line, path in [("HT1080", f"{HT}/HT1080_LOY_vs_sgCTRL_all_gene_effects.csv"),
                   ("A549", f"{A5}/A549_LOY_vs_sgCTRL_all_gene_effects.csv")]:
    ag = pd.read_csv(path, usecols=["gene_id", "mean_CPM", "mean_difference",
                                    "exclusion_reason"]).set_index("gene_id")
    t1[f"{line}_exclusion_reason"] = t1.gene_id.map(ag.exclusion_reason)
    t1[f"{line}_allgene_mean_CPM"] = t1.gene_id.map(ag.mean_CPM)
    t1[f"{line}_allgene_mean_difference"] = t1.gene_id.map(ag.mean_difference)
    # practical expression call used in the report: mean CPM >= 1 across all libraries
    t1[f"{line}_CPM_ge1"] = t1[f"{line}_allgene_mean_CPM"] >= 1
t1["DU145_CPM_ge1"] = (t1.DU145_CTRL_mean_CPM + t1.DU145_LOY_mean_CPM) / 2 >= 1


def sign_agree(a, b):
    return np.where(a.isna() | b.isna(), np.nan, np.sign(a) == np.sign(b))


t1["HT1080_same_sign_as_DU145"] = sign_agree(t1.DU145_mean_difference, t1.HT1080_mean_difference)
t1["A549_same_sign_as_DU145"] = sign_agree(t1.DU145_mean_difference, t1.A549_mean_difference)
t1.to_csv(f"{OUT}/task1_inhouse_lines_focus_genes.csv", index=False)

# genome-wide context: concordance of the DU145 FDR<0.05 set with HT1080 / A549
g = du.merge(ht, on="gene_id").merge(a5, on="gene_id").merge(hb, on="gene_id", how="left")
sig = g[g.FDR < 0.05]
rows = []
for line, col in [("HT1080", "HT1080_mean_difference"), ("A549", "A549_mean_difference")]:
    for label, s in [("DU145_FDR05", sig),
                     ("DU145_FDR05_and_CPMge1_in_other", sig[sig[f"{line}_mean_CPM"] >= 1]),
                     ("DU145_thesis_DE", g[g.thesis_DE == True])]:
        s = s.dropna(subset=[col])
        agree = (np.sign(s.mean_difference) == np.sign(s[col])).mean()
        rho = s[["mean_difference", col]].corr(method="spearman").iloc[0, 1]
        rows.append(dict(line=line, gene_set=label, n=len(s), sign_agreement=agree,
                         spearman=rho))
    rows.append(dict(line=line, gene_set="all_shared_primary_genes", n=len(g),
                     sign_agreement=(np.sign(g.mean_difference) == np.sign(g[col])).mean(),
                     spearman=g[["mean_difference", col]].corr(method="spearman").iloc[0, 1]))
pd.DataFrame(rows).to_csv(f"{OUT}/task1_genomewide_concordance.csv", index=False)
print(pd.DataFrame(rows))

# ---------------------------------------------------------------- Task 2: public datasets
# Every entry: key -> (file, id-column, join-type, logFC col, FDR col, polarity, description)
# polarity = +1 if positive logFC means "higher when Y is lost/absent"; -1 if reversed.
P = []


def add(key, group, path, join, lfc, fdr, pol, desc, gene_col=None, sep=","):
    P.append(dict(key=key, group=group, path=path, join=join, lfc=lfc, fdr=fdr, pol=pol,
                  desc=desc, gene_col=gene_col))


# A. bona fide Y-loss vs Y-present (somatic LOY)
add("A549_clones_LOYvROY", "A_LOY", f"{AN}/A549_clonal_LOY/DE/clonal_LOY_v_ROY_DESeq2.csv",
    "id", "logFC", "FDR", 1, "A549 4 LOY vs 4 ROY clones (DESeq2, unpaired; LOY-ROY)")
add("A549_multiome_pb_LOYvROY", "A_LOY",
    f"{AN}/A549_multiome_pseudobulk/DE/multiome_pb_LOY_v_ROY_DESeq2.csv", "symbol", "logFC",
    "FDR", 1, "A549 multiome pseudobulk LOY vs ROY cells (DESeq2; LOY-ROY)", "gene_name")
add("MB49_cellline_YnegvYpos", "A_LOY",
    f"{AN}/GSE229233/DE/gse229233_cellline_Yneg_v_Ypos_DESeq2.csv", "mouse", "log2FoldChange",
    "padj", 1, "GSE229233 mouse MB49 bladder cancer cell line Y- vs Y+ (5v5, DESeq2)", "gene_name")
add("MB49_tumor_YnegvYpos", "A_LOY",
    f"{AN}/GSE229233/DE/gse229233_tumor_Yneg_v_Ypos_DESeq2.csv", "mouse", "log2FoldChange",
    "padj", 1, "GSE229233 MB49 subcutaneous tumour Y- vs Y+ (5v5, DESeq2; includes host stroma/immune)",
    "gene_name")
# B. constitutional 45,X vs 46,XY (X-dose matched, Y absent; also PAR dose 1 vs 2)
add("GSE126712_lymph_XOvXY", "B_XOvXY", f"{AN}/GSE126712/DE/gse126712_XO_v_XY_DE.csv", "id",
    "logFC", "FDR", 1, "GSE126712 primary lymphocytes 45,X vs 46,XY (limma-voom)")
add("GSE207112_troph_XOvXY", "B_XOvXY", f"{AN}/GSE207112/DE/gse207112_XO_v_XY_DE_blocked.csv",
    "id", "logFC", "FDR", 1, "GSE207112 isogenic iPSC-trophoblast XO vs XY (blocked by round)")
add("GSE264742_NCC_XOvXY", "B_XOvXY", f"{AN}/GSE264742/DE/gse264742_XO_v_XY_DE_blocked.csv",
    "id", "logFC", "FDR", 1, "GSE264742 isogenic iPSC-neural crest XO vs XY (panel+round blocked)")
add("ERP156313_fetal_XOvXY", "B_XOvXY",
    f"{AN}/ERP156313/DE/erp156313_consensus_XO_v_XY_DE.csv", "id", "logFC", "FDR", 1,
    "ERP156313 fetal tissues 45,X vs 46,XY, tissue-adjusted consensus (limma-voom)")
add("PageLab_Fib_X0vXY", "B_XOvXY", f"{AN}/PageLab_SCA/DE/trackA_X0_vs_XY_Fib.csv", "symbol",
    "log2FoldChange", "padj", 1, "Page lab SCA fibroblasts 45,X vs 46,XY (DESeq2)", "gene")
add("PageLab_LCL_X0vXY", "B_XOvXY", f"{AN}/PageLab_SCA/DE/trackA_X0_vs_XY_LCL.csv", "symbol",
    "log2FoldChange", "padj", 1, "Page lab SCA LCLs 45,X vs 46,XY (DESeq2)", "gene")
# C. not LOY: single Y-gene KO, Y gain (sign flipped so + = 'less Y'), shown for context only
add("GSE301297_undiff_UTYKOvWT", "C_other",
    f"{AN}/GSE301297/DE/gse301297_SEES3undiff_UTYKO_v_WT_DE.csv", "id", "logFC", "FDR", 1,
    "GSE301297 male hESC (SEES3) UTY-KO vs WT, undifferentiated (single Y gene KO, not LOY)")
add("GSE301297_undiff_DKOvWT", "C_other",
    f"{AN}/GSE301297/DE/gse301297_SEES3undiff_DKO_v_WT_DE.csv", "id", "logFC", "FDR", 1,
    "GSE301297 male hESC UTX/UTY double KO vs WT, undifferentiated (not LOY)")
add("XYY_integrated_flipped", "C_other",
    f"{AN}/integrated_XYY_XY/DE/integrated_consensus_XYY_v_XY_DE.csv", "id", "logFC", None, -1,
    "Integrated iPSC+NSC XYY vs XY consensus, SIGN FLIPPED (+ = lower with fewer Y); Y gain, not LOY")

mouse_note = []
long_rows = []
for p in P:
    d = pd.read_csv(p["path"])
    if p["join"] == "id":
        d["stable_id"] = stable(d["gene_id"])
        m = ref[["gene_name", "stable_id"]].merge(d, on="stable_id", how="left",
                                                   suffixes=("", "_src"))
    elif p["join"] == "symbol":
        col = p["gene_col"]
        d["_sym"] = d[col].astype(str).str.replace(r"\.\d+$", "", regex=True)
        d = d.sort_values(p["lfc"], key=lambda s: -s.abs()).drop_duplicates("_sym")
        m = ref[["gene_name"]].merge(d, left_on="gene_name", right_on="_sym", how="left", suffixes=("", "_src"))
    else:  # mouse: naive symbol case-folding ortholog match
        d["_sym"] = d["gene_name"].astype(str).str.upper()
        d = d.drop_duplicates("_sym")
        m = ref[["gene_name"]].merge(d, left_on="gene_name", right_on="_sym", how="left",
                                     suffixes=("", "_mouse"))
    for _, r in m.iterrows():
        lfc = r.get(p["lfc"], np.nan)
        fdr = r.get(p["fdr"], np.nan) if p["fdr"] else np.nan
        if p["key"].startswith("XYY"):
            fdr = np.nanmin([r.get("fdr_iPSC", np.nan), r.get("fdr_NSC", np.nan)]) \
                if pd.notna(lfc) else np.nan
        long_rows.append(dict(dataset=p["key"], group=p["group"], description=p["desc"],
                              gene_name=r["gene_name"],
                              logFC_Yloss_polarity=(p["pol"] * lfc) if pd.notna(lfc) else np.nan,
                              FDR=fdr, tested=pd.notna(lfc)))
L = pd.DataFrame(long_rows)
L = L.merge(ref[["gene_name", "panel", "DU145_mean_difference", "DU145_FDR"]], on="gene_name")
L["same_sign_as_DU145"] = np.where(L.tested, np.sign(L.logFC_Yloss_polarity) ==
                                   np.sign(L.DU145_mean_difference), np.nan)
L["sig_FDR05"] = L.FDR < 0.05
L.to_csv(f"{OUT}/task2_public_datasets_long.csv", index=False)


def cell(r):
    if not r.tested:
        return "nt"
    s = f"{r.logFC_Yloss_polarity:+.2f}"
    if r.FDR < 0.05:
        s += "*"
    return s


L["cell"] = L.apply(cell, axis=1)
W = L.pivot(index="gene_name", columns="dataset", values="cell")
W = W.reindex(GENES)[[p["key"] for p in P]]
W.to_csv(f"{OUT}/task2_public_datasets_wide.csv")

# per-dataset summary of the 30 DU145 focus genes
S = []
foc = FOCUS_DOWN + FOCUS_UP
for p in P:
    x = L[(L.dataset == p["key"]) & L.gene_name.isin(foc) & L.tested]
    xs = x[x.sig_FDR05]
    S.append(dict(dataset=p["key"], group=p["group"], description=p["desc"],
                  focus_tested=len(x), focus_same_sign=int(x.same_sign_as_DU145.sum()),
                  focus_sig=len(xs), focus_sig_same_sign=int(xs.same_sign_as_DU145.sum()),
                  focus_sig_opposite=int((~xs.same_sign_as_DU145.astype(bool)).sum())))
S = pd.DataFrame(S)
S.to_csv(f"{OUT}/task2_public_datasets_summary.csv", index=False)
print(S.to_string())

# membership in crossdataset strict DEG lists (context only)
sd = pd.read_csv(f"{AN}/crossdataset_LOY/DEG_overlap/strict_DEG_lists_long.csv")
sd[sd.gene.isin(GENES)].to_csv(f"{OUT}/task2_crossdataset_strictDEG_membership.csv", index=False)
print(sd[sd.gene.isin(GENES)].to_string())
