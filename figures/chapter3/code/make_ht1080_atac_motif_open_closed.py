"""Supplementary figure: motif evidence at opening vs closing chromatin in Y-depleted HT1080 cells.

Sources (ATACseq/SEPT_08_2026_REPROCESS/motifs/TF_extension; TF extension 2026-09-08, no refit):
  concordance/TF_ATAC_footprint_RNA.csv  - per JASPAR2024 motif: HOMER lost/gained, chromVAR, TOBIAS BINDetect
  chromVAR/motifs_zscores.csv            - per-library chromVAR deviation z-scores
TF RNA is re-joined to the thesis RNA result (SEPT_11_2026_REANALYZE gene_consistency.csv); the RNA columns
in the concordance table come from the superseded DESeq2 branch and are not used.
Motif rows: motifs with FDR < 0.05 in HOMER (either direction, BH across 879 motifs x 2 directions) or chromVAR
(limma, independent-replicate design), one representative per Vierstra archetype (lowest FDR in any test).
Inputs are copied into source_data/ht1080_atac_motif_open_closed with SHA256 hashes.
"""
from pathlib import Path
import sys, shutil, hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_atac_motif_open_closed"
TFX = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS/motifs/TF_extension"
RNA = ROOT.parent / "RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1/tables/gene_consistency.csv"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from thesis_style import apply_style
from atac_volcano import spread
apply_style(8)
STEM = "ch03_supp_ht1080_atac_motif_open_closed"
LOWER, HIGHER, GREY = "#3C5488", "#E64B35", "#555555"
DIVERGE = LinearSegmentedColormap.from_list("div", [LOWER, "#FFFFFF", HIGHER])
HALO = [pe.withStroke(linewidth=1.5, foreground="white")]
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]

# ---- freeze inputs ----
DATA.mkdir(parents=True, exist_ok=True)
files = {"TF_ATAC_footprint_RNA.csv": TFX / "concordance/TF_ATAC_footprint_RNA.csv",
         "chromVAR_motifs_zscores.csv": TFX / "chromVAR/motifs_zscores.csv",
         "thesis_RNA_gene_consistency.csv": RNA}
pd.DataFrame([dict(file=k, source=str(v), sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k, v in files.items()]
             ).to_csv(DATA / "source_manifest.csv", index=False)
for k, v in files.items():
    shutil.copyfile(v, DATA / k)

t = pd.read_csv(DATA / "TF_ATAC_footprint_RNA.csv").drop_duplicates("motif_id")
assert len(t) == 879
z = pd.read_csv(DATA / "chromVAR_motifs_zscores.csv").set_index("id")
rna = pd.read_csv(DATA / "thesis_RNA_gene_consistency.csv")
rna_by = {s.upper(): r for s, r in zip(rna.gene_name, rna.itertuples())}

HL, HG, CV = "HOMER_lost_BH_FDR_both_directions", "HOMER_gained_BH_FDR_both_directions", "chromVAR_FDR"
sig = t[(t[HL] < .05) | (t[HG] < .05) | (t[CV] < .05)].copy()
sig["best"] = sig[[HL, HG, CV]].min(axis=1)
rows = sig.sort_values(["best", "motif_id"]).drop_duplicates("archetype").copy()
rows["log2FE_lost"] = np.log2(rows.HOMER_lost_fold_enrichment.clip(lower=.25))
rows["log2FE_gained"] = np.log2(rows.HOMER_gained_fold_enrichment.clip(lower=.25))
rows = rows.sort_values(["log2FE_gained", "log2FE_lost"], ascending=[False, True]).reset_index(drop=True)
# TF RNA (thesis branch): the motif's own TF; dimers use the first member; case-insensitive symbol match
def tf_rna(name):
    r = rna_by.get(name.split("::")[0].upper())
    return (np.nan, np.nan, False, False) if r is None else (r.logFC, r.FDR, bool(r.directionally_consistent), True)
rows[["tf_rna_log2FC", "tf_rna_FDR", "tf_rna_DE", "tf_expressed"]] = [tf_rna(n) for n in rows.motif_name]
for period, prefix in [('August','AUG'),('May','MAY')]:
    values=[]
    for name in rows.motif_name:
        r=rna_by.get(name.split('::')[0].upper())
        fc=getattr(r,prefix+'_logFC',np.nan); fdr=getattr(r,prefix+'_FDR',np.nan)
        values.append((fc,fdr,bool(np.isfinite(fdr) and fdr<.05 and abs(fc)>=1)))
    rows[['tf_rna_'+period+'_log2FC','tf_rna_'+period+'_FDR','tf_rna_'+period+'_DE']]=values
rows["label"] = [f"{n.upper()} · {a}" if not a.startswith("unmapped") else n.upper() for n, a in zip(rows.motif_name, rows.archetype)]  # display only; RNA lookup keeps JASPAR names
rows.to_csv(DATA / "panel_a_motif_rows.csv", index=False)

fig = plt.figure(figsize=(183 / 25.4, 200 / 25.4))


def letter(x, y, ch):
    fig.text(x, y, ch, fontsize=12, weight="semibold", va="bottom")


# a  evidence grid
COLS = [("log2FE_lost", HL, "HOMER\nclosing", "log₂ FE", 1.5),
        ("log2FE_gained", HG, "HOMER\nopening", "log₂ FE", 1.5),
        ("chromVAR_delta_z", CV, "chromVAR\nΔz", "Δz", 18),
        ("LOY_CTRL_change", None, "TOBIAS\nΔscore", "score Δ", .25),
        ("tf_rna_log2FC", "tf_rna_FDR", "RNA\ncombined", "log₂FC", 3)]
# Retain HOMER closing/opening alongside chromVAR, footprint and RNA evidence.
COLS += [('tf_rna_August_log2FC','tf_rna_August_FDR','RNA\nAugust','log₂FC',3),
         ('tf_rna_May_log2FC','tf_rna_May_FDR','RNA\nMay','log₂FC',3)]
a = fig.add_axes([.20, .47, .49, .43])
n = len(rows)
for j, (val, fdr, head, unit, lim) in enumerate(COLS):
    homer = val.startswith("log2FE_")
    norm = Normalize(0, lim) if homer else Normalize(-lim, lim)
    cmap = LinearSegmentedColormap.from_list(val, ["white", LOWER if val == "log2FE_lost" else HIGHER]) if homer else DIVERGE
    for i, r in rows.iterrows():
        v = r[val]
        missing=not np.isfinite(v)
        col='#D9D9D9' if missing else cmap(norm(np.clip(v,0 if homer else -lim,lim)))
        a.add_patch(plt.Rectangle((j-.49,i-.49),.98,.98,facecolor=col,edgecolor='white',lw=.3))
        if not missing and fdr is not None:
            f=r[fdr]
            if np.isfinite(f) and f<.05:
                a.scatter(j,i,s=9,color='#111111',lw=0,zorder=4)
                if val.startswith('tf_rna_') and abs(v)>=1:
                    a.scatter(j,i,s=40,facecolor='none',edgecolor='#111111',lw=.7,zorder=4)
    a.text(j, -1.0, head, ha="center", va="bottom", fontsize=6.2, weight="semibold", linespacing=1.1)
    a.text(j, n - .2, (f"0–{lim:g}\n{unit}" if homer else f"±{lim:g}\n{unit}"), ha="center", va="top", fontsize=5.3, color=GREY, linespacing=1.1)
a.set(xlim=(-.6, len(COLS) - .4), ylim=(n - .5, -.6))
a.set_xticks([]); a.set_yticks(range(n), rows.label, fontsize=6)
a.tick_params(axis="y", length=0, pad=3)
for sp in a.spines.values(): sp.set_visible(False)
for j in range(len(COLS) - 1):
    a.axvline(j + .5, color="#EEEEEE", lw=.6, zorder=0)
a.axvline(1.5, color="#CCCCCC", lw=.8, zorder=0)
a.axvline(3.5, color="#CCCCCC", lw=.8, zorder=0)
letter(.02, .925, "a")
cbx = fig.add_axes([.215, .415, .12, .008])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-1, 1), cmap=DIVERGE), cax=cbx, orientation="horizontal")
cb.set_ticks([-1, 0, 1], labels=["Lower", "0", "Higher"]); cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=5.2, length=2)
leg = [Line2D([], [], ls='', marker='o', markersize=3, color='#111111', label='FDR < 0.05'),
       Line2D([], [], ls='', marker=r'$\odot$', markersize=7, color='#111111', label='RNA: FDR < 0.05 and |log₂FC| ≥ 1'),
       plt.Rectangle((0,0),1,1,facecolor='#D9D9D9',label='Unavailable RNA')]
fig.legend(handles=leg,loc='upper left',bbox_to_anchor=(.36,.43),fontsize=5.4,handletextpad=.3,borderpad=0,labelspacing=.4,frameon=False)
fig.text(.43,.372,'TOBIAS: colour only; pooled-site tests',fontsize=5.3,color=GREY)
for xpos,col,label in [(.12,LOWER,'HOMER closing'),(.27,HIGHER,'HOMER opening')]:
    cx=fig.add_axes([xpos,.367,.10,.006])
    cm=LinearSegmentedColormap.from_list(label,['white',col])
    cbh=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1.5),cmap=cm),cax=cx,orientation='horizontal')
    cbh.set_ticks([0,1.5],labels=['≤0','1.5']);cbh.outline.set_visible(False);cbh.ax.tick_params(labelsize=5,length=2)
    cx.set_title(label+' (log₂ FE)',fontsize=5,pad=3)

# b  chromVAR per-library z-scores, same rows
b = fig.add_axes([.735, .47, .19, .43])
zraw = z.loc[rows.motif_id, SAMPLES].to_numpy()
zm = zraw - zraw.mean(axis=1, keepdims=True)   # row-centred: each motif relative to its own mean across libraries
vmax = float(np.nanpercentile(np.abs(zm), 98))
im = b.imshow(zm, cmap=DIVERGE, vmin=-vmax, vmax=vmax, aspect="auto")
b.set_xticks(range(6), [f"{'C' if s.startswith('CTRL') else 'L'} p{s.split('_')[1]}" for s in SAMPLES], rotation=90,
             fontsize=6)
for t_, s in zip(b.get_xticklabels(), SAMPLES):
    t_.set_color(LOWER if s.startswith("CTRL") else HIGHER)
b.set_yticks([]); b.tick_params(length=0); b.set_ylim(n - .5, -.6)
for sp in b.spines.values(): sp.set_visible(False)
b.axvline(2.5, color="white", lw=1.2)
b.set_title("chromVAR deviation z\n(row-centred)", fontsize=6.6, weight="semibold", pad=6, linespacing=1.1)
cb2 = fig.colorbar(im, cax=fig.add_axes([.94, .60, .009, .15])); cb2.outline.set_visible(False)
cb2.ax.tick_params(labelsize=5.4, length=2); cb2.set_label("z minus motif mean", fontsize=5.6)
letter(.705, .925, "b")

# d  TOBIAS differential footprint volcano (pooled conditions)
d = fig.add_axes([.15, .06, .70, .28])
lp = -np.log10(t.LOY_CTRL_pvalue.clip(lower=1e-300))
d.scatter(t.LOY_CTRL_change, lp, s=4, c="#CFCFCF", lw=0, zorder=1)
sel = t.motif_id.isin(rows.motif_id)
d.scatter(t.LOY_CTRL_change[sel], lp[sel], s=14, c=np.where(t.LOY_CTRL_change[sel] < 0, LOWER, HIGHER),
          edgecolors="#222222", linewidths=.3, zorder=3)
# Labels: the six panel-a motifs with the strongest footprint loss, in an ordered column right of the cluster;
# plus every panel-a motif with a positive change.
neg = t[sel & (t.LOY_CTRL_change < -.12)].assign(lp=lp).sort_values("lp", ascending=False).drop_duplicates("motif_name").head(6)
ys = spread(neg.lp.to_numpy(), 7, lo=60, hi=146)
for r, ty in zip(neg.itertuples(), ys):
    d.annotate(r.motif_name.upper(), (r.LOY_CTRL_change, r.lp), xytext=(-.085, ty), fontsize=5.6, color=LOWER,
               ha="left", va="center", path_effects=HALO, arrowprops=dict(arrowstyle="-", color="#999999", lw=.3, shrinkB=1))
for r in t[sel & (t.LOY_CTRL_change > 0)].itertuples():
    d.annotate(r.motif_name.upper(), (r.LOY_CTRL_change, -np.log10(max(r.LOY_CTRL_pvalue, 1e-300))), xytext=(6, 4),
               textcoords="offset points", fontsize=5.6, color=HIGHER, ha="left", path_effects=HALO)
d.axvline(0, color="#999999", lw=.5, zorder=0)
d.set_ylim(0, 150)
assert lp.max() < 150
d.set(xlabel="Footprint score change (Y-depleted − sgCTRL)", ylabel="−log₁₀ P (pooled sites)")
d.text(.03, .03, "Filled: motifs in panel a\nPooled-condition BINDetect;\nnot a replicate-aware test",
       transform=d.transAxes, va="bottom", ha="left", fontsize=5.6, color=GREY, linespacing=1.3)
letter(.02, .335, "c")

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=180)
plt.close(fig)
json.dump({"motif_rows": rows.label.tolist(), 
           "chromVAR_sig": int((t[CV] < .05).sum()), "zscore_vmax_p98": vmax}, open(DATA / "figure_checks.json", "w"), indent=1)
print("written", STEM, n, "rows")
