"""Main HT1080 ATAC-seq figure from frozen source data; run from any directory.

Inputs are produced by prepare_ht1080_atac_figures.py (Windows) and
extract_ht1080_atac_locus_coverage.py (WSL, pyBigWig). No model is refitted.
"""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle, FancyArrow
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_atac_accessibility"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
sys.path.insert(0, str(Path(__file__).resolve().parent))
from atac_volcano import select_labels, draw_volcano
apply_style(8)
STEM = "ch03_main_ht1080_atac_accessibility"
COL = {"CTRL": "#3C5488", "LOY": "#E64B35"}
LABEL = {"CTRL": "sgCTRL", "LOY": "Y-depleted"}
DIR = {"gained": "#E64B35", "lost": "#3C5488", "ns": "#C6C6C6"}
FEAT = {"promoter": ("Promoter (TSS ±2 kb)", "#00A087"), "gene_body": ("Gene body", "#91D1C2"),
        "intergenic": ("Intergenic", "#B09C85")}
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
GREY = "#555555"

qc = pd.read_csv(DATA / "library_qc.csv").set_index("sample").loc[SAMPLES]
pc = pd.read_csv(DATA / "pca_scores.csv").set_index("sample").loc[SAMPLES]
var = pd.read_csv(DATA / "pca_variance.csv").set_index("component").variance_fraction
da = pd.read_csv(DATA / "differential_peaks.csv")
feat = pd.read_csv(DATA / "feature_enrichment_fisher.csv")
loci = pd.read_csv(DATA / "locus_definitions.csv")
cov = pd.read_csv(DATA / "locus_coverage_50bp.csv")
genes = pd.read_csv(DATA / "locus_genes_unlabelled.csv")
n_gain, n_lost = (da.direction == "gained").sum(), (da.direction == "lost").sum()
assert (len(da), n_gain, n_lost) == (84510, 95, 663)

fig = plt.figure(figsize=(183 / 25.4, 215 / 25.4))
a = fig.add_axes([.10, .695, .18, .21])
b = fig.add_axes([.37, .695, .175, .21])
c = fig.add_axes([.62, .695, .36, .21])
d = fig.add_axes([.20, .445, .31, .14])
e = fig.add_axes([.65, .445, .33, .14])

def letter(ax, ch, dx=-.075, dy=.035):
    pos = ax.get_position()
    fig.text(pos.x0 + dx, pos.y1 + dy, ch, fontsize=12, weight="semibold", va="bottom")

# a  chrY signal within the ATAC libraries
y = qc.chrY_per_million_nonY
x = np.array([0, 1, 2, 3.5, 4.5, 5.5])
a.bar(x, y, color=[COL[s[:4].rstrip("_")] for s in SAMPLES], width=.78)
for xi, yi in zip(x, y):
    a.text(xi, yi + 25, f"{yi:,.0f}", ha="center", va="bottom", fontsize=5.6)
a.set_xticks(x, [f"p{s.split('_')[1]}" for s in SAMPLES], fontsize=6, rotation=90)
a.set(ylim=(0, 1650), xlim=(-.6, 6.1), ylabel="chrY fragments per 10⁶\nnon-Y fragments")
a.set_yticks([0, 500, 1000, 1500])
for cx, cond in [(1, "CTRL"), (4.5, "LOY")]:
    a.text(cx, -.2, LABEL[cond], transform=a.get_xaxis_transform(), ha="center", va="top", fontsize=6.8, color=COL[cond])
ratio = y[3:].mean() / y[:3].mean()
a.text(.76, .24, f"Y-depleted mean\n{100 * ratio:.1f}% of sgCTRL", transform=a.transAxes, ha="center", fontsize=6.3, color=GREY)
letter(a, "a", dx=-.09)

# b  PCA with passage-matched pairs joined
for p in ["11", "12", "21"]:
    s0, s1 = pc.loc[f"CTRL_{p}"], pc.loc[f"LOY_{p}"]
    b.annotate("", xy=(s1.PC1, s1.PC2), xytext=(s0.PC1, s0.PC2),
               arrowprops=dict(arrowstyle="-|>", color="#AAAAAA", lw=.7, shrinkA=4, shrinkB=4, mutation_scale=6))
for s, r in pc.iterrows():
    b.scatter(r.PC1, r.PC2, s=30, c=COL[r.condition], edgecolors="white", linewidths=.4, zorder=3)
    b.annotate(f"p{r.passage}", (r.PC1, r.PC2), xytext=(4, 4), textcoords="offset points", fontsize=6.5, color=COL[r.condition])
b.set(xlabel=f"PC1 ({100 * var['PC1']:.1f}%)", ylabel=f"PC2 ({100 * var['PC2']:.1f}%)", xlim=(-140, 130), ylim=(-110, 110))
b.text(.02, .03, "84,510 tested regions", transform=b.transAxes, ha="left", va="bottom", fontsize=6.3, color=GREY)
letter(b, "b", dx=-.07)

# c  volcano; named-gene labels only (LOC-labelled version is supplementary)
da["y"] = -np.log10(da.PValue)
lab = select_labels(da, n_p=10, n_ext=3, forced=dict(zip(loci.peak_id, loci.label)), include_loc=False)
lab[["peak_id", "text", "direction", "logFC", "PValue", "FDR", "annotation", "nearest_tss_distance_signed"]].to_csv(
    DATA / "volcano_labels_main.csv", index=False)
n_treat, cut = draw_volcano(c, da, lab, label_gap=.37, headroom=1.22)
c.set_xticks([-4, -2, 0, 2])
letter(c, "c", dx=-.065)

# d  genomic-feature composition vs tested non-significant background
bars = [("ns", "Tested, not significant"), ("gained", "Gained"), ("lost", "Lost")]
for i, (key, name) in enumerate(bars):
    sub = feat[feat.direction == ("lost" if key == "ns" else key)]
    total = int(sub.background_total.iloc[0] if key == "ns" else sub.foreground_total.iloc[0])
    left = 0
    for f in ["promoter", "gene_body", "intergenic"]:
        r = sub[sub.feature == f].iloc[0]
        pct = r.background_pct if key == "ns" else r.foreground_pct
        d.barh(i, pct, left=left, color=FEAT[f][1], height=.66, edgecolor="white", linewidth=.5)
        d.text(left + pct / 2, i, f"{pct:.0f}%", ha="center", va="center", fontsize=6.3,
               color="white" if f == "promoter" else "#222222")
        left += pct
    d.text(-2, i, f"{name}\n(n = {total:,})", ha="right", va="center", fontsize=6.8,
           color=DIR[key] if key != "ns" else "#222222", linespacing=1.15)
d.set(xlim=(0, 100), ylim=(2.5, -.6), yticks=[])
d.set_xticks([0, 25, 50, 75, 100]); d.set_xlabel("Regions (%)")
d.spines["left"].set_visible(False)
minbh = feat.fisher_BH.min()
d.text(.5, -.50, f"Fisher's exact test vs tested background: all BH-adjusted P ≥ {minbh:.2f}",
       transform=d.transAxes, ha="center", fontsize=6.3, color=GREY)
d.legend(handles=[Patch(facecolor=v[1], label=v[0]) for v in FEAT.values()], loc="lower center",
         bbox_to_anchor=(.35, 1.03), ncol=3, fontsize=6.5, handlelength=1, columnspacing=1, handletextpad=.4)
letter(d, "d", dx=-.185)

# e  effect-size distribution of the FDR < 0.05 regions
bins = np.arange(-4.0, 2.6, .1)
hist_stats = {}
for key in ["lost", "gained"]:
    x = da.loc[da.direction == key, "logFC"]
    e.hist(x, bins=bins, color=DIR[key], edgecolor="white", linewidth=.25)
    hist_stats[key] = dict(n=len(x), median=float(x.median()), abs_gt1=int((x.abs() > 1).sum()),
                           min=float(x.min()), max=float(x.max()))
for v in [-1, 1]:
    e.axvline(v, color="#888888", ls=(0, (4, 3)), lw=.6, zorder=0)
e.set(xlabel="log₂ fold change (Y-depleted / sgCTRL)", ylabel="Regions (FDR < 0.05)", xlim=(-4.2, 2.7))
top = e.get_ylim()[1]
for key, xpos, ha in [("lost", .03, "left"), ("gained", .97, "right")]:
    s_ = hist_stats[key]
    e.text(xpos, .97, f"{key} (n = {s_['n']})\nmedian {s_['median']:+.2f}\n|log₂FC| > 1: {s_['abs_gt1']}",
           transform=e.transAxes, ha=ha, va="top", fontsize=6.2, color=DIR[key], linespacing=1.3)
json.dump(hist_stats, open(DATA / "logFC_histogram_stats.json", "w"), indent=1)
pd.DataFrame({"bin_left": bins[:-1],
              **{k: np.histogram(da.loc[da.direction == k, "logFC"], bins=bins)[0] for k in ["lost", "gained"]}}
             ).to_csv(DATA / "logFC_histogram_bins.csv", index=False)
letter(e, "e", dx=-.075)

# f  coverage tracks at representative promoter peaks
ncol, x0, width, gap = len(loci), .085, .895, .032
w = (width - gap * (ncol - 1)) / ncol
fig.text(x0 - .075, .335, "f", fontsize=12, weight="semibold", va="bottom")
for j, r in enumerate(loci.itertuples()):
    left = x0 + j * (w + gap)
    sub = cov[cov.peak_id == r.peak_id]
    ymax = sub.mean_fragment_cpm.max() * 1.08
    for k, cond in enumerate(["CTRL", "LOY"]):
        ax = fig.add_axes([left, .235 - k * .072, w, .062])
        ax.axvspan(r.peak_start, r.peak_end, color="#EEEEEE", zorder=0, lw=0)
        mat = []
        for s in [s for s in SAMPLES if s.startswith(cond)]:
            t = sub[sub["sample"] == s]
            xs = (t.bin_start + t.bin_end) / 2
            ax.plot(xs, t.mean_fragment_cpm, color=COL[cond], lw=.45, alpha=.7)
            mat.append(t.mean_fragment_cpm.to_numpy())
        ax.fill_between(xs, np.mean(mat, axis=0), color=COL[cond], alpha=.35, lw=0)
        ax.set(xlim=(r.window_start, r.window_end), ylim=(0, ymax), xticks=[])
        ax.set_yticks([0, round(ymax * .8, 1)])
        ax.tick_params(axis="y", labelsize=5.8, length=2, pad=1)
        ax.spines["bottom"].set_visible(False)
        if j == 0:
            ax.set_ylabel(LABEL[cond], fontsize=6.6, color=COL[cond], labelpad=2)
        if k == 0:
            ax.set_title(r.label, fontsize=7.5, pad=12, weight="semibold", style="italic")
            ax.text(.5, 1.03, f"{r.direction}  log₂FC {r.logFC:+.2f}  FDR {r.FDR:.3f}",
                    transform=ax.transAxes, ha="center", va="bottom", fontsize=5.9, color=DIR[r.direction])
    g = fig.add_axes([left, .03, w, .092])
    g.set(xlim=(r.window_start, r.window_end), ylim=(-3.0, .8)); g.axis("off")
    gsub = genes[genes.peak_id == r.peak_id].reset_index(drop=True)
    for i, gr in gsub.iterrows():
        yy = -i * 1.05
        s0, s1 = max(gr.start, r.window_start), min(gr.end, r.window_end)
        g.plot([s0, s1], [yy, yy], color="#444444", lw=1.1, solid_capstyle="butt")
        if r.window_start <= gr.tss <= r.window_end:
            sign = 1 if gr.strand == "+" else -1
            g.plot([gr.tss, gr.tss], [yy, yy + .45], color="#444444", lw=.7)
            g.annotate("", xy=(gr.tss + sign * 700, yy + .45), xytext=(gr.tss, yy + .45),
                       arrowprops=dict(arrowstyle="-|>", color="#444444", lw=.7, mutation_scale=5, shrinkA=0, shrinkB=0))
        lx = min(max((s0 + s1) / 2, r.window_start + 1800), r.window_end - 1800)
        g.text(lx, yy - .18, gr.symbol, ha="center", va="top", fontsize=5.8, style="italic", color="#444444")
    # scale bar and coordinates
    g.plot([r.window_end - 2000, r.window_end], [-2.75, -2.75], color="#222222", lw=.8)
    g.text(r.window_end - 1000, -2.7, "2 kb", ha="center", va="bottom", fontsize=5.6)
    g.text(r.window_start, -2.7, f"{r.chrom}:{r.window_start / 1e6:.2f} Mb", ha="left", va="bottom", fontsize=5.6, color=GREY)
fig.text(.025, .20, "Fragments per million non-Y fragments", rotation=90, ha="center", va="center", fontsize=6.6)

handles = [Patch(facecolor=COL[k], label=LABEL[k]) for k in ["CTRL", "LOY"]]
handles += [Line2D([], [], ls="", marker="o", markersize=3.5, color=DIR["gained"], label="Gained (FDR < 0.05)"),
            Line2D([], [], ls="", marker="o", markersize=3.5, color=DIR["lost"], label="Lost (FDR < 0.05)"),
            Line2D([], [], ls="", marker="o", markersize=4.5, markerfacecolor="none", markeredgecolor="#222222", label="TREAT |log₂FC| > 0.5")]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .995), ncol=5, fontsize=7, handlelength=1.1, columnspacing=1.6)

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=180)
plt.close(fig)
json.dump({"regions_tested": len(da), "gained": int(n_gain), "lost": int(n_lost), "treat_fdr005": n_treat,
           "fdr_raw_P_cutoff": float(10 ** -cut), "chrY_ratio_LOY_over_CTRL": float(ratio),
           "pca_variance_PC1_PC2": [float(var["PC1"]), float(var["PC2"])], "volcano_labels": int(len(lab)), "logFC_hist": hist_stats,
           "font": "Open Sans", "size_mm": [183, 215], "png_dpi": 600, "models_refitted": False},
          open(DATA / "main_figure_checks.json", "w"), indent=1)
print("written", STEM)
