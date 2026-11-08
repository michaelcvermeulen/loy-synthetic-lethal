"""Draft: HT1080 RNA global structure (F01) and expression changes (F02) as one figure.

Panel order: a sampling design, b period-adjusted PCA, c assigned reads,
d correlation heatmap, e volcano, f August/May agreement. Reads only the frozen
source_data tables; no analysis is recomputed and no existing export is changed.
Run from any working directory.
"""
from pathlib import Path
import sys, json, hashlib
from itertools import groupby
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from adjustText import adjust_text
from matplotlib import patheffects as pe

FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
GS = FIG / "source_data/ht1080_rna_global_structure"
EC = FIG / "source_data/ht1080_rna_expression_changes"
OUT = FIG / "drafts"
sys.path.insert(0, str(THESIS / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(7.5); np.random.seed(42)
STEM = "ch03_main_ht1080_rna_structure_and_expression_draft"
COL = {"CTRL": "#3C5488", "LOY": "#E64B35", "PARENTAL": "#888888"}
GREY = "#C8C8C8"; MUTED = "#555555"; REF = "#888888"
SHAPE = {"AUG_2025": "o", "MAY_2026": "^"}
TICK, ANNOT, GENE = 7, 6.5, 6          # pt; axis labels use the 7.5-pt base
HALO = [pe.withStroke(linewidth=1.6, foreground="white")]   # keeps gene labels legible over points

# ---- frozen inputs -------------------------------------------------------
pc, qc, cr, design = [pd.read_csv(GS / f"Fig01_design_QC_{x}.csv") for x in "abcd"]
plot_pc = pd.read_csv(GS / "panel_a_adjusted_all15.csv")
var = pd.read_csv(GS / "panel_a_variance.csv")
var = var[var.view == "adjusted"].set_index("component").variance_fraction
d = pd.read_csv(EC / "all_genes.csv")
assert len(plot_pc) == 15 and len(qc) == 15 and len(cr) == 144 and len(design) == 12 and len(d) == 25526

qc_order = [f"{p}{i}" for p in ["AC", "AL", "AP", "ML", "MC"] for i in [1, 2, 3]]
qc = qc.set_index("short").loc[qc_order].reset_index()
order = [f"{p}{i}" for p in ["AC", "AL", "MC", "ML"] for i in [1, 2, 3]]
matrix = cr.pivot(index="sample_2", columns="sample_1", values="r").loc[order, order]
assert np.allclose(matrix, matrix.T) and np.allclose(np.diag(matrix), 1)
dist = 1 - matrix.to_numpy(); dist = (dist + dist.T) / 2; np.fill_diagonal(dist, 0)
tree = linkage(squareform(dist), method="average", optimal_ordering=True)
order = [order[i] for i in leaves_list(tree)[::-1]]
saved_order = pd.read_csv(GS / "correlation_cluster_order.csv")["sample"].tolist()
assert order == saved_order, "Clustering must reproduce the saved heatmap order"
matrix = matrix.loc[order, order]

s = d.directionally_consistent.astype(bool)
assert (s & (d.logFC < 0)).sum() == 433 and (s & (d.logFC > 0)).sum() == 229
cols = np.where(s, np.where(d.logFC > 0, COL["LOY"], COL["CTRL"]), GREY)
volcano_lab = pd.concat([d[s & (np.sign(d.logFC) == g)].sort_values("FDR").head(12) for g in [-1, 1]])
agree = d[s & ~d.gene_name.str.startswith("ENSG")].copy()
agree["gap"] = (agree.AUG_logFC - agree.MAY_logFC).abs()
agree["minabs"] = agree[["AUG_logFC", "MAY_logFC"]].abs().min(axis=1)
agree = agree[agree.gap <= .5]
agree = pd.concat([agree[np.sign(agree.logFC) == g].sort_values(["minabs", "FDR", "gene_id"],
                   ascending=[False, True, True]).head(6) for g in [-1, 1]])
# Label selections must match the approved F02 display.
assert set(volcano_lab.gene_id) == set(pd.read_csv(EC / "labelled_genes.csv").gene_id)
assert set(agree.gene_id) == set(pd.read_csv(EC / "panel_b_labelled_genes.csv").gene_id)

# ---- layout (mm, origin top-left) ---------------------------------------
W, H = 183, 133
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
def ax_mm(x, y, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)
def letter(x, y, t):
    fig.text(x / W, 1 - y / H, t, fontsize=10, weight="semibold", va="top")

R1, R2, AH, AH2 = 11, 72, 44, 50   # row tops and axis heights
a = ax_mm(21, R1, 42, AH)        # design
b = ax_mm(77, R1, 43, AH)        # PCA
c = ax_mm(134, R1, 47, AH)       # assigned reads
hm_x, cell = 16.5, 3.3           # heatmap
dh = ax_mm(hm_x, R2, 12 * cell, 12 * cell)
e = ax_mm(77, R2, 43, AH2)       # volcano
f = ax_mm(134, R2, 47, 47)       # agreement
for ax, t, x in [(a, "a", 1), (b, "b", 64), (c, "c", 121)]: letter(x, R1 - 3, t)
for ax, t, x in [(dh, "d", 1), (e, "e", 64), (f, "f", 121)]: letter(x, R2 - 3, t)

# shared legend
handles = [Patch(facecolor=COL[k], label=l) for k, l in
           [("CTRL", "sgCTRL"), ("LOY", "Y-depleted"), ("PARENTAL", "Parental")]]
handles += [Line2D([], [], color="#333333", marker=m, lw=0, markersize=4.5, label=l)
            for m, l in [("o", "Aug 2025"), ("^", "May 2026")]]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, 1 - .6 / H), ncol=5,
           handlelength=1, handletextpad=.5, columnspacing=1.8, fontsize=7.5, borderaxespad=0)

# a: sampling design (spacing denotes order, not elapsed time)
for tr, y in [("CTRL", 1), ("LOY", 0)]:
    sub = design[design.treatment == tr].sort_values("occasion")
    a.plot(sub.occasion, [y] * len(sub), c=COL[tr], lw=.9, zorder=1)
    for r in sub.itertuples():
        a.scatter(r.occasion, y, c=COL[tr], marker=SHAPE[r.batch], s=22, lw=0, zorder=3)
        a.text(r.occasion, y + .17, str(r.passage_number), ha="center", fontsize=ANNOT, color=COL[tr])
a.set(xlim=(.5, 7.5), ylim=(-.45, 1.6))
a.set_xticks([2, 6], ["Aug 2025", "May 2026"], fontsize=TICK)
a.set_yticks([1, 0], ["sgCTRL", "Y-depleted"], fontsize=TICK)
for y, lab in [(1, "2×"), (0, "2× CP + Q53")]:
    a.annotate(lab, xy=(0, y), xycoords=a.get_yaxis_transform(), xytext=(-3, -9),
               textcoords="offset points", ha="right", va="center", fontsize=ANNOT - .5, color=REF)
a.tick_params(length=0)
for sp in a.spines.values(): sp.set_visible(False)
a.set_xlabel("Successive passages within each period", fontsize=ANNOT, labelpad=3)
a.text(.5, 1.0, "Same established populations", transform=a.transAxes, ha="center",
       va="top", fontsize=ANNOT, color=MUTED)
a.text(.5, .5, "12 libraries; six passage matches", transform=a.transAxes, ha="center",
       va="center", fontsize=ANNOT, color=MUTED)

# b: period-adjusted PCA, all 15 libraries
for r in plot_pc.itertuples():
    b.scatter(r.PC1, r.PC2, c=COL[r.treatment], marker=SHAPE[r.batch], s=20,
              edgecolors="white", linewidths=.35, zorder=3)
    dx, dy, ha = 4, 3, "left"
    if r.treatment == "PARENTAL" and r.passage_number == 2: dx, dy, ha = -5, 6, "right"
    if r.batch == "AUG_2025" and r.treatment == "CTRL" and r.passage_number in [2, 3]: dx, dy, ha = -4, -6, "right"
    if r.batch == "MAY_2026" and r.treatment == "LOY" and r.passage_number == 3: dx, dy, ha = -4, -5, "right"
    b.annotate(str(r.passage_number), (r.PC1, r.PC2), xytext=(dx, dy), textcoords="offset points",
               ha=ha, fontsize=ANNOT - .5, color=COL[r.treatment])
b.set(xlabel=f"PC1 ({100 * var['PC1']:.1f}%)", ylabel=f"PC2 ({100 * var['PC2']:.1f}%)",
      xlim=(-110, 94), ylim=(-95, 92))
b.set_xticks([-100, -50, 0, 50]); b.set_yticks([-50, 0, 50])
b.tick_params(labelsize=TICK)
b.text(.03, .03, "15 libraries; 25,526 genes", transform=b.transAxes, fontsize=ANNOT - .5, color=MUTED)

# c: assigned reads, source panel order
c.bar(np.arange(15), qc.assigned_fraction * 100, color=[COL[t] for t in qc.treatment], width=.74)
c.set_xticks(np.arange(15), qc.short, rotation=90, fontsize=ANNOT - .3)
c.tick_params(axis="x", length=2, pad=2); c.tick_params(axis="y", labelsize=TICK)
c.set(ylabel="Assigned reads (%)", ylim=(0, 70), xlim=(-.7, 14.3))
c.set_yticks([0, 20, 40, 60])

# d: unadjusted core correlation heatmap with dendrogram and annotation strips
im = dh.imshow(matrix.to_numpy(), cmap="magma", vmin=float(matrix.min().min()), vmax=1,
               interpolation="nearest", aspect="equal")
dh.set_yticks(range(12), order, fontsize=ANNOT - .3); dh.yaxis.tick_right()
dh.set_xticks([]); dh.tick_params(length=0, pad=2)
for sp in dh.spines.values(): sp.set_visible(False)
top, side = R2, 12 * cell
den = ax_mm(1, top, 8.5, side)
dendrogram(tree, orientation="left", no_labels=True, ax=den, color_threshold=0,
           above_threshold_color="#666666", link_color_func=lambda _: "#666666")
for coll in den.collections: coll.set_linewidth(.5)
den.set_ylim(0, 120); den.axis("off")
meta = pc.set_index("short").loc[order]
per = ax_mm(10, top, 3.2, side); con = ax_mm(13.7, top, 1.8, side)
start = 0
for batch, grp in groupby(meta.batch):
    n = len(list(grp))
    lab, shade = {"AUG_2025": ("Aug 2025", "#EEEEEE"), "MAY_2026": ("May 2026", "#DADADA")}[batch]
    per.add_patch(Rectangle((0, start), 1, n, facecolor=shade, edgecolor="white", lw=.4))
    per.text(.5, start + n / 2, lab, ha="center", va="center", rotation=90, fontsize=ANNOT - .5)
    start += n
for i, tr in enumerate(meta.treatment):
    con.add_patch(Rectangle((0, i), 1, 1, facecolor=COL[tr], edgecolor="white", lw=.3))
for axx in (per, con): axx.set(xlim=(0, 1), ylim=(12, 0)); axx.axis("off")
cb = fig.colorbar(im, cax=ax_mm(hm_x + 4, top + side + 3.2, side - 8, 1.8), orientation="horizontal")
cb.set_ticks([.94, .96, .98, 1]); cb.ax.tick_params(labelsize=ANNOT, length=1.5, pad=1.5)
cb.outline.set_visible(False); cb.set_label("Pearson r", fontsize=ANNOT, labelpad=2)

# e: volcano
e.scatter(d.logFC, -np.log10(d.FDR), c=cols, s=np.where(s, 4, 2), alpha=.55, lw=0, rasterized=True)
for x in [-1, 1]: e.axvline(x, c=REF, ls="--", lw=.5)
e.axhline(-np.log10(.05), c=REF, ls="--", lw=.5)
e.set(xlabel="log$_2$FC (Y-depleted / sgCTRL)", ylabel="−log$_{10}$(BH-adjusted P)",
      xlim=(-7.5, 6.5), ylim=(-.05, 7.4))
e.set_xticks([-6, -3, 0, 3, 6]); e.tick_params(labelsize=TICK)
e.text(.02, .02, "433 lower", transform=e.transAxes, fontsize=ANNOT, color=COL["CTRL"], weight="semibold")
e.text(.98, .02, "229 higher", transform=e.transAxes, fontsize=ANNOT, color=COL["LOY"],
       weight="semibold", ha="right")
texts = [e.text(r.logFC, -np.log10(r.FDR) + .08, r.gene_name, fontsize=GENE, path_effects=HALO,
                color=COL["CTRL"] if r.logFC < 0 else COL["LOY"]) for r in volcano_lab.itertuples()]
hit = s & (d.FDR < .01)       # avoid the dense upper points, not the whole cloud
adjust_text(texts, x=d.loc[hit, "logFC"].to_numpy(), y=-np.log10(d.loc[hit, "FDR"].to_numpy()), ax=e,
            iter_lim=500, expand=(1.1, 1.3), force_text=(.3, .5), force_static=(.2, .3),
            arrowprops=dict(arrowstyle="-", color=REF, lw=.35))

# f: August versus May effects, square with equal scaling
f.scatter(d.AUG_logFC, d.MAY_logFC, c=GREY, s=2, alpha=.35, lw=0, rasterized=True)
f.scatter(d.loc[s, "AUG_logFC"], d.loc[s, "MAY_logFC"], c=cols[s], s=5, alpha=.75, lw=0, rasterized=True)
lim = 10.5   # common limits; leaves label columns clear of the data
f.plot([-lim, lim], [-lim, lim], ls="--", c=REF, lw=.5)
f.axvline(0, c="#BBBBBB", lw=.4); f.axhline(0, c="#BBBBBB", lw=.4)
f.set(xlabel="August log$_2$FC", ylabel="May log$_2$FC", xlim=(-lim, lim), ylim=(-lim, lim))
f.set_xticks([-8, -4, 0, 4, 8]); f.set_yticks([-8, -4, 0, 4, 8]); f.tick_params(labelsize=TICK)
f.set_aspect("equal", adjustable="box")
r_all = d.AUG_logFC.corr(d.MAY_logFC)
f.text(.03, .97, f"r = {r_all:.2f} (all genes)", transform=f.transAxes, va="top", fontsize=ANNOT)
for g in [-1, 1]:
    lab = agree[np.sign(agree.logFC) == g].sort_values("MAY_logFC")
    # Label columns sit off the diagonal: blue above it on the left, vermilion below it on the right.
    ys = np.linspace(-6.4, -1.2, len(lab)) if g < 0 else np.linspace(1.0, 6.2, len(lab))
    for r, y in zip(lab.itertuples(), ys):
        f.annotate(r.gene_name, xy=(r.AUG_logFC, r.MAY_logFC), xytext=(-6.9 if g < 0 else 6.6, y),
                   ha="right" if g < 0 else "left", va="center", fontsize=GENE, path_effects=HALO,
                   color=COL["CTRL"] if g < 0 else COL["LOY"],
                   arrowprops=dict(arrowstyle="-", color=REF, lw=.35), zorder=6)

# ---- export ---------------------------------------------------------------
OUT.mkdir(exist_ok=True)
records = []
for ext in ["pdf", "svg", "png"]:
    p = OUT / f"{STEM}.{ext}"
    fig.savefig(p, dpi=600, facecolor="white")
    records.append({"file": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
fig.savefig(OUT / f"{STEM}_preview.png", dpi=200, facecolor="white")
plt.close(fig)
checks = {"status": "draft; not selected for the thesis", "size_mm": [W, H], "font": "Open Sans",
          "panel_map": {"a": "F01 d", "b": "F01 a", "c": "F01 b", "d": "F01 c", "e": "F02 a", "f": "F02 b"},
          "sources": ["source_data/ht1080_rna_global_structure", "source_data/ht1080_rna_expression_changes"],
          "heatmap_order_matches_saved": True, "label_sets_match_saved": True,
          "pearson_r_all_genes": round(float(r_all), 4), "analysis_recomputed": False, "exports": records}
(OUT / f"{STEM}_render_checks.json").write_text(json.dumps(checks, indent=2))
print(json.dumps({k: checks[k] for k in ["size_mm", "pearson_r_all_genes"]}))
