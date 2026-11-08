"""Supplementary figure: abundance trend and gene-set context of HT1080 accessibility changes.

Sources (ATACseq/SEPT_08_2026_REPROCESS; nothing refitted here):
  differential/trended_normalisation/ (88_trended_normalisation_sensitivity.R)
  integration/thesis_RNA_SEPT11/promoter_accessibility_GSEA/ (87_..., part A)
  integration/thesis_RNA_SEPT11/distal_nearest_gene_sets/ (87_..., part B; rGREAT 2.12.2)
a = median log2FC by baseline logCPM decile, primary TMM vs trended (loess) normalisation, with call counts;
b = Hallmark promoter GSEA NES, primary vs trended ranking, for sets named in either;
c = Hallmark NES, trended promoter ATAC vs RNA (same 13,300 genes);
d = distal lost and gained regions: top 20 and all 10 sets passing both GREAT tests, ranked by gene-based FDR;
e = distal gained regions: sets passing both tests, with and without chr18 regions.
"""
from pathlib import Path
import sys, shutil, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_atac_gene_set_context"
RUN = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM = "ch03_supp_ht1080_atac_gene_set_context"
LOWER, HIGHER, GREY = "#3C5488", "#E64B35", "#555555"

DATA.mkdir(parents=True, exist_ok=True)
srcs = {"trended_summary.csv": RUN / "differential/trended_normalisation/summary.csv",
        "trended_logFC_by_logCPM_decile.csv": RUN / "differential/trended_normalisation/trended_logFC_by_logCPM_decile.csv",
        "library_signal_to_noise.csv": RUN / "differential/trended_normalisation/library_signal_to_noise.csv",
        "fgsea_wide.csv": RUN / "integration/thesis_RNA_SEPT11/promoter_accessibility_GSEA/fgsea_wide.csv",
        "rGREAT_enrichment.csv": RUN / "integration/thesis_RNA_SEPT11/distal_nearest_gene_sets/rGREAT_enrichment.csv"}
for k, v in srcs.items():
    shutil.copyfile(v, DATA / k)
pd.DataFrame([dict(file=k, source=str(v), sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k, v in srcs.items()]
             ).to_csv(DATA / "source_manifest.csv", index=False)
summ = pd.read_csv(DATA / "trended_summary.csv")
dec = pd.read_csv(DATA / "trended_logFC_by_logCPM_decile.csv")
snr = pd.read_csv(DATA / "library_signal_to_noise.csv")
w = pd.read_csv(DATA / "fgsea_wide.csv")
gr = pd.read_csv(DATA / "rGREAT_enrichment.csv")

fig = plt.figure(figsize=(183 / 25.4, 175 / 25.4))


def letter(x, y, ch):
    fig.text(x, y, ch, fontsize=12, weight="semibold", va="bottom")


def nice(s):
    s = s.replace("HALLMARK_", "").replace("GOBP_", "").replace("REACTOME_", "").replace("_", " ").lower()
    for a, b in [("tnfa", "TNF-α"), ("nfkb", "NF-κB"), ("il6", "IL-6"), ("jak stat3", "JAK/STAT3"), ("kras", "KRAS"),
                 ("e2f", "E2F"), ("myc", "MYC"), ("g2m", "G2M"), ("il2 stat5", "IL-2/STAT5"), ("tgf beta", "TGF-β"),
                 ("smad2 smad3 smad4", "SMAD2/3/4"), ("smad2 3 smad4", "SMAD2/3/SMAD4"), (" dn", " down"), ("interferon", "IFN")]:
        s = s.replace(a, b)
    s = s[0].upper() + s[1:]
    return SHORT.get(s, s)


SHORT = {"Response to type ii IFN": "Response to type II interferon",
         "Regulation of cellular response to vascular endothelial growth factor stimulus": "Response to VEGF (regulation)",
         "Regulation of cellular response to growth factor stimulus": "Response to growth factor (regulation)",
         "Downregulation of SMAD2/3/SMAD4 transcriptional activity": "SMAD2/3/4 activity, downregulation",
         "Transcriptional activity of SMAD2/3/4 heterotrimer": "SMAD2/3/4 heterotrimer activity",
         "SMAD2/3/4 heterotrimer regulates transcription": "SMAD2/3/4 heterotrimer transcription",
         "Negative regulation of wnt signaling pathway": "Negative regulation of WNT signalling",
         "Signaling by tgfb family members": "TGF-β family signalling",
         "Signaling by TGF-β receptor complex": "TGF-β receptor complex signalling",
         "Runx2 regulates bone development": "RUNX2 regulates bone development",
         "IFN gamma response": "IFN-γ response", "IFN alpha response": "IFN-α response"}


# ---- a: abundance trend ----
a = fig.add_axes([.075, .675, .24, .264])
x = dec.dec.to_numpy()
a.plot(x, dec.median_logFC_primary, "-o", color=GREY, ms=3, lw=1, label="Primary (TMM)")
a.plot(x, dec.median_logFC_trended, "-o", color="#00A087", ms=3, lw=1, label="Trended (loess)")
a.axhline(0, color="#BBBBBB", lw=.6)
a.set_xticks([1, 5, 10]); a.set_xlabel("Baseline accessibility (logCPM decile)", fontsize=6.4)
a.set_ylabel("Median log2FC (Y-depleted / sgCTRL)", fontsize=6.4); a.tick_params(labelsize=6)
a.legend(loc="upper left", fontsize=5.8, frameon=False, handlelength=1.4, borderpad=0)
p = summ[(summ.arm == "trended") & (summ.FDR == .05)].iloc[0]
tx = snr.groupby("condition").tss_enrichment.mean()
a.text(.98, .04, f"FDR < 0.05, lost / gained\nPrimary {p.primary_lost} / {p.primary_gained}\nTrended {p.lost} / {p.gained}\n"
       f"Mean TSS enrichment\nsgCTRL {tx['CTRL']:.1f}, Y-depleted {tx['LOY']:.1f}",
       transform=a.transAxes, ha="right", va="bottom", fontsize=5.5, color=GREY, linespacing=1.15)
letter(.01, .95, "a")

# ---- b: Hallmark NES primary vs trended ----
b = fig.add_axes([.46, .675, .19, .264])
h = w[w.collection == "Hallmark"].copy()
pick = h[(h.padj_ATAC_median < .05) | (h.padj_ATAC_trended < .05)].copy()
pick = pick.sort_values("NES_ATAC_trended")
for i, r in enumerate(pick.itertuples()):
    b.plot([r.NES_ATAC_median, r.NES_ATAC_trended], [i, i], color="#CCCCCC", lw=1, zorder=1)
    b.scatter(r.NES_ATAC_median, i, s=14, color="white", edgecolor=GREY, lw=.8, zorder=2)
    b.scatter(r.NES_ATAC_trended, i, s=16, color="#00A087" if r.padj_ATAC_trended < .05 else "white", edgecolor="#00A087", lw=.8, zorder=3)
b.axvline(0, color="#999999", lw=.6, ls=(0, (4, 3)))
b.set_yticks(range(len(pick)), [nice(s) for s in pick.pathway], fontsize=5.9)
b.set_ylim(len(pick) - .5, -.8); b.set_xlim(-2.3, 2.3)
b.tick_params(axis="y", length=0); b.tick_params(axis="x", labelsize=6); b.spines["left"].set_visible(False)
b.set_xlabel("Promoter-accessibility NES", fontsize=6.4)
b.legend(handles=[Line2D([], [], ls="", marker="o", ms=3.8, mfc="white", mec=GREY, label="Primary"),
                  Line2D([], [], ls="", marker="o", ms=3.8, color="#00A087", label="Trended (filled FDR < 0.05)")],
         loc="lower center", bbox_to_anchor=(.4, 1.0), ncol=1, fontsize=5.6, frameon=False, handletextpad=.2, borderpad=0)
pick[["pathway", "NES_ATAC_median", "padj_ATAC_median", "NES_ATAC_trended", "padj_ATAC_trended"]].to_csv(DATA / "panel_b.csv", index=False)
letter(.30, .95, "b")

# ---- c: Hallmark NES trended ATAC vs RNA ----
c = fig.add_axes([.76, .675, .21, .264])
c.axhline(0, color="#DDDDDD", lw=.6); c.axvline(0, color="#DDDDDD", lw=.6)
sig = h.padj_ATAC_trended < .05
c.scatter(h.NES_RNA[~sig], h.NES_ATAC_trended[~sig], s=10, color="#BBBBBB", lw=0)
c.scatter(h.NES_RNA[sig], h.NES_ATAC_trended[sig], s=16, color="#00A087", lw=0)
LAB = {"HALLMARK_INTERFERON_GAMMA_RESPONSE": (-1.0, -.45), "HALLMARK_KRAS_SIGNALING_UP": (-.25, -2.2),
       "HALLMARK_IL6_JAK_STAT3_SIGNALING": (-2.3, -2.3), "HALLMARK_E2F_TARGETS": (.35, 1.55)}
for r in h[h.pathway.isin(LAB)].itertuples():
    c.annotate(nice(r.pathway), (r.NES_RNA, r.NES_ATAC_trended), xytext=LAB[r.pathway], textcoords="data", fontsize=5.3,
               color=GREY, arrowprops=dict(arrowstyle="-", color="#999999", lw=.5, shrinkA=0, shrinkB=2))
from scipy.stats import spearmanr
rho = spearmanr(h.NES_RNA, h.NES_ATAC_trended)
c.text(.03, .97, f"Hallmark sets (n = {len(h)})\nSpearman ρ = {rho.statistic:.2f}", transform=c.transAxes, va="top", fontsize=5.6, color=GREY)
c.set_xlabel("RNA NES (same 13,300 genes)", fontsize=6.4); c.set_ylabel("Promoter-accessibility NES (trended)", fontsize=6.4)
c.tick_params(labelsize=6); c.set_xlim(-2.4, 2.4); c.set_ylim(-2.4, 2.4)
letter(.68, .95, "c")

# ---- d: both directions, ranked by gene-based FDR across collections ----
import textwrap
selected=[]
for arm, left, width, bottom, height, colour, heading in [
    ('lost',.325,.095,.065,.49,LOWER,'Closing: top 20 of 31 sets'),
    ('gained',.815,.085,.34,.215,HIGHER,'Opening: all 10 sets')]:
    tab=gr[(gr.arm==arm)&(gr.BH<.05)&(gr.BH_hyper<.05)&(gr.fold_enrichment>1)].sort_values(['BH_hyper','BH','id']).head(20)
    selected.append(tab)
    da=fig.add_axes([left,bottom,width,height])
    da.set_xscale('log');da.set_xlim(1,16);da.set_xticks([1,2,4,8,16],[1,2,4,8,16]);da.minorticks_off()
    da.set_ylim(len(tab)-.5,-.5);da.set_yticks([]);da.spines['left'].set_visible(False);da.tick_params(labelsize=5)
    da.axvline(1,color='#BBBBBB',lw=.5)
    for j,r in enumerate(tab.itertuples()):
        da.plot([1,r.fold_enrichment],[j,j],color=colour,lw=1)
        da.scatter(r.fold_enrichment,j,color=colour,s=12,zorder=3)
        da.text(-.08,j,'\n'.join(textwrap.wrap(nice(r.id),35)),transform=da.get_yaxis_transform(),ha='right',va='center',fontsize=5.2,linespacing=1.05)
        for xpos,val in [(1.20,r.BH),(1.83,r.BH_hyper)]:
            da.text(xpos,j,f'{val:.1e}' if val<.001 else f'{val:.3f}',transform=da.get_yaxis_transform(),fontsize=4.8,va='center',ha='center')
    da.set_xlabel('Fold enrichment',fontsize=5.8)
    fig.text(left-.29,.595,heading,fontsize=6.4,weight='semibold',color=colour)
    for xpos,txt in [(1.20,'Region'),(1.83,'Gene')]:
        da.text(xpos,1.025,txt+'\nFDR',transform=da.transAxes,ha='center',fontsize=4.8)
lost=selected[0]
pd.concat(selected).to_csv(DATA/'panel_d.csv',index=False)
letter(.01,.615,'d')

# ---- e: distal gained, with vs without chr18 ----
e = fig.add_axes([.79, .065, .18, .185])
both = gr[(gr.arm == "gained") & (gr.BH < .05) & (gr.BH_hyper < .05) & (gr.fold_enrichment > 1)].sort_values("p_value")
ex = gr[gr.arm == "gained_excl_chr18"].set_index(["collection", "id"])
rows = []
for r in both.itertuples():
    k = (r.collection, r.id)
    rows.append(dict(id=r.id, FE_all=r.fold_enrichment, BH_all=r.BH, hyper_all=r.BH_hyper,
                     FE_excl=ex.loc[k, "fold_enrichment"] if k in ex.index else np.nan,
                     BH_excl=ex.loc[k, "BH"] if k in ex.index else np.nan,
                     hyper_excl=ex.loc[k, "BH_hyper"] if k in ex.index else np.nan))
pe = pd.DataFrame(rows).iloc[::-1].reset_index(drop=True)
pe.to_csv(DATA / "panel_e.csv", index=False)
for i, r in pe.iterrows():
    e.plot([r.FE_all, r.FE_excl if np.isfinite(r.FE_excl) else 1], [i, i], color="#DDDDDD", lw=1, zorder=1)
    e.scatter(r.FE_all, i, s=16, color=HIGHER, zorder=3)
    if np.isfinite(r.FE_excl):
        both_ok = (r.BH_excl < .05) and (r.hyper_excl < .05)
        e.scatter(r.FE_excl, i, s=16, color=HIGHER if both_ok else "white", edgecolor=HIGHER, marker="s", lw=.8, zorder=3)
    else:
        e.text(.83, i, "< 5 genes", va="center", fontsize=5.0, color=GREY)
e.set_xscale("log"); e.set_xlim(.8, 20); e.set_xticks([1, 2, 5, 10], ["1", "2", "5", "10"]); e.minorticks_off()
e.axvline(1, color="#999999", lw=.6, ls=(0, (4, 3)))
e.set_yticks(range(len(pe)), [nice(s) for s in pe.id], fontsize=5.6); e.tick_params(axis="y", length=0); e.tick_params(axis="x", labelsize=6)
e.set_ylim(len(pe) - .5, -.6)
e.set_xlabel("GREAT fold enrichment", fontsize=6.4)
e.legend(handles=[Line2D([], [], ls="", marker="o", ms=4, color=HIGHER, label="All (331 regions)"),
                  Line2D([], [], ls="", marker="s", ms=4, mfc="white", mec=HIGHER, label="Excluding chr18 (open: not both significant)")],
         loc="lower right", bbox_to_anchor=(1, 1.0), ncol=2, fontsize=5.4, frameon=False, handletextpad=.2, columnspacing=.6, borderpad=0)
letter(.53, .275, "e")

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=.035); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=200, bbox_inches="tight", pad_inches=.035)
plt.close(fig)
print("written", STEM, len(pick), len(lost), len(pe))
