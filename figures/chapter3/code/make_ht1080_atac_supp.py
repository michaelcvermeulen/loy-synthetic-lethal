"""Supplementary HT1080 ATAC-seq QC and sensitivity figures from frozen source data.

Inputs come from prepare_ht1080_atac_figures.py. No model is refitted.
TSS profiles use the pipeline's own normalisation (mean of the outer 100 bp on
each flank), so the plotted centre matches the reported TSS enrichment score.
"""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_atac_accessibility"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
COL = {"CTRL": "#3C5488", "LOY": "#E64B35"}
LABEL = {"CTRL": "sgCTRL", "LOY": "Y-depleted"}
DIR = {"gained": "#E64B35", "lost": "#3C5488", "ns": "#C6C6C6"}
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
STYLE = {"11": "-", "12": (0, (4, 2)), "21": (0, (1, 1.3))}
GREY = "#555555"
cond = lambda s: s.split("_")[0]
short = lambda s: f"{'C' if cond(s) == 'CTRL' else 'L'} p{s.split('_')[1]}"


def letter(fig, ax, ch, dx=-.07, dy=.03):
    pos = ax.get_position()
    fig.text(pos.x0 + dx, pos.y1 + dy, ch, fontsize=12, weight="semibold", va="bottom")


def save(fig, stem):
    (FIG / "previews").mkdir(exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        target, stage = FIG / f"{stem}.{ext}", FIG / f"{stem}.staging.{ext}"
        fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
    fig.savefig(FIG / "previews" / f"{stem}_preview.png", dpi=180)
    plt.close(fig)


def sample_legend(fig, y):
    h = [Line2D([], [], color=COL[c], lw=1.2, label=LABEL[c]) for c in ["CTRL", "LOY"]]
    h += [Line2D([], [], color="#444444", lw=1, ls=STYLE[p], label=f"passage {p}") for p in STYLE]
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(.5, y), ncol=5, fontsize=7, handlelength=2.2, columnspacing=1.6)


# ---------------------------------------------------------------- QC figure
qc = pd.read_csv(DATA / "library_qc.csv").set_index("sample").loc[SAMPLES]
frag = pd.read_csv(DATA / "fragment_lengths.csv")
tss = pd.read_csv(DATA / "tss_profiles.csv")
corr = pd.read_csv(DATA / "spearman_correlations.csv", index_col=0).loc[SAMPLES, SAMPLES]
md = pd.read_csv(DATA / "matched_depth_peak_counts.csv").set_index("sample").loc[SAMPLES]

fig = plt.figure(figsize=(183 / 25.4, 225 / 25.4))
a = fig.add_axes([.09, .715, .38, .19]); b = fig.add_axes([.60, .715, .38, .19])
sample_legend(fig, .995)
checks = {}
for s in SAMPLES:
    f = frag[frag["sample"] == s]
    f = f[f.length <= 800]
    a.plot(f.length, 100 * f.fragments / frag.loc[frag["sample"] == s, "fragments"].sum(), color=COL[cond(s)],
           ls=STYLE[s.split("_")[1]], lw=.8)
    t = tss[tss["sample"] == s].insertions.to_numpy().astype(float)
    flank = np.r_[t[:100], t[-100:]].mean()
    score = t[1950:2050].mean() / flank
    assert abs(score - qc.loc[s, "tss_enrichment"]) < 1e-6, s
    offs = tss.loc[tss["sample"] == s, "offset"].to_numpy()
    smooth = np.convolve(t / flank, np.ones(25) / 25, mode="same")
    b.plot(offs[12:-12], smooth[12:-12], color=COL[cond(s)], ls=STYLE[s.split("_")[1]], lw=.8)
checks["tss_scores_reproduced"] = True
a.set(xlabel="Fragment length (bp)", ylabel="Fragments (%)", xlim=(0, 800))
for x, lab in [(100, "NFR"), (200, "mono"), (400, "di")]:
    a.text(x, a.get_ylim()[1] * .97, lab, ha="center", va="top", fontsize=6.3, color=GREY)
b.set(xlabel="Distance from TSS (bp)", ylabel="Tn5 insertions\n(relative to flanks)", xlim=(-2000, 2000))
b.set_xticks([-2000, -1000, 0, 1000, 2000])
b.text(.98, .96, "TSS enrichment\n" + "\n".join(f"{short(s)}: {qc.loc[s, 'tss_enrichment']:.1f}" for s in SAMPLES),
       transform=b.transAxes, ha="right", va="top", fontsize=6.2, color=GREY, linespacing=1.3)
letter(fig, a, "a"); letter(fig, b, "b")

# c  per-library metrics
metrics = [("accepted_fragments", "Retained fragments (M)", 1e-6, "{:.1f}"),
           ("frip_own", "FRiP, own peaks", 1, "{:.2f}"),
           ("n_peaks", "Peaks called (thousands)", 1e-3, "{:.1f}"),
           ("matched_depth_peaks", "Peaks at 13.6 M fragments\n(thousands)", 1e-3, "{:.1f}")]
qc["matched_depth_peaks"] = md.matched_depth_peaks
xs = np.array([0, 1, 2, 3.5, 4.5, 5.5])
w = .185
for i, (col, name, scale, fmt) in enumerate(metrics):
    ax = fig.add_axes([.07 + i * (w + .06), .465, w, .145])
    vals = qc[col] * scale
    ax.bar(xs, vals, color=[COL[cond(s)] for s in SAMPLES], width=.8)
    for x, v in zip(xs, vals):
        ax.text(x, v * 1.01, fmt.format(v), ha="center", va="bottom", fontsize=5.3, rotation=90)
    ax.set_xticks(xs, [f"p{s.split('_')[1]}" for s in SAMPLES], fontsize=6, rotation=90)
    ax.set_ylim(0, vals.max() * 1.3); ax.set_title(name, fontsize=7, weight="semibold", pad=4)
    ax.tick_params(axis="y", labelsize=6.5)
    if i == 0:
        letter(fig, ax, "c", dx=-.06, dy=.035)
    checks[col] = [float(v) for v in vals]

# d  Spearman correlation heatmap with clustering
order = list(corr.index)
dist = 1 - corr.to_numpy(); np.fill_diagonal(dist, 0); dist = (dist + dist.T) / 2
tree = linkage(squareform(dist), method="average", optimal_ordering=True)
order = [order[i] for i in leaves_list(tree)[::-1]]
m = corr.loc[order, order]
d = fig.add_axes([.165, .105, .165, .20])
im = d.imshow(m.to_numpy(), cmap="magma", vmin=float(m.min().min()), vmax=1, aspect="auto")
d.set_xticks(range(6), [short(s) for s in order], rotation=90, fontsize=6.3)
d.set_yticks(range(6), [short(s) for s in order], fontsize=6.3); d.yaxis.tick_right()
d.tick_params(length=0, pad=2)
for sp in d.spines.values(): sp.set_visible(False)
for i, s in enumerate(order):
    fig.add_artist(Rectangle((d.get_position().x0 - .014, d.get_position().y1 - (i + 1) * d.get_position().height / 6),
                             .01, d.get_position().height / 6, transform=fig.transFigure, color=COL[cond(s)], ec="white", lw=.4))
dend = fig.add_axes([.055, d.get_position().y0, .09, d.get_position().height])
dendrogram(tree, orientation="left", no_labels=True, ax=dend, color_threshold=0, above_threshold_color="#666666",
           link_color_func=lambda _: "#666666")
dend.set_ylim(0, 60); dend.axis("off")
cb = fig.colorbar(im, cax=fig.add_axes([.395, .135, .01, .14])); cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=6.3, length=2); cb.set_label("Spearman ρ (log-CPM)", fontsize=6.8)
letter(fig, d, "d", dx=-.13, dy=.058)
checks["correlation_order"] = order
checks["correlation_offdiag_range"] = [float(m.min().min()), float(m.where(~np.eye(6, dtype=bool)).max().max())]

# e  TSS-centred coverage heatmaps (rows sorted by mean signal across all six libraries)
tz = np.load(DATA / "tss_heatmap_matrix.npz")
mats = {s_: tz[s_] for s_ in SAMPLES}
row_order = np.argsort(-np.mean([mats[s_].mean(axis=1) for s_ in SAMPLES], axis=0))
vmax = float(np.percentile(np.concatenate([mats[s_].ravel() for s_ in SAMPLES]), 99))
n_bins = mats[SAMPLES[0]].shape[1]; xs_kb = (np.arange(n_bins) + .5) * 6 / n_bins - 3
x0, cw, gap = .51, .058, .008
prof_max = max(mats[s_].mean(axis=0).max() for s_ in SAMPLES) * 1.08
for k, s_ in enumerate(SAMPLES):
    left = x0 + k * (cw + gap)
    pa = fig.add_axes([left, .305, cw, .05])
    pa.plot(xs_kb, mats[s_].mean(axis=0), color=COL[cond(s_)], lw=.8)
    pa.set(xlim=(-3, 3), ylim=(0, prof_max), xticks=[])
    pa.set_title(short(s_), fontsize=6.6, color=COL[cond(s_)], pad=2)
    pa.tick_params(axis="y", labelsize=5.6, length=2, pad=1)
    if k:
        pa.set_yticklabels([])
    else:
        pa.set_ylabel("Mean", fontsize=6.3, labelpad=2)
    ha_ = fig.add_axes([left, .095, cw, .20])
    hm = ha_.imshow(mats[s_][row_order], aspect="auto", cmap="magma", vmin=0, vmax=vmax,
                    extent=(-3, 3, len(row_order), 0), interpolation="antialiased")
    ha_.set_xticks([-3, 0, 3], ["−3", "0", "+3"], fontsize=5.8)
    for t_, al in zip(ha_.get_xticklabels(), ["left", "center", "right"]):
        t_.set_ha(al)
    ha_.tick_params(axis="x", length=2, pad=1)
    ha_.set_yticks([])
    for sp in ha_.spines.values(): sp.set_visible(False)
    if k == 0:
        ha_.set_ylabel(f"{len(row_order):,} TSSs\n(sorted by mean signal)", fontsize=6.3, labelpad=3)
        letter(fig, pa, "e", dx=-.06, dy=.013)
fig.text(x0 + 3 * (cw + gap) - gap / 2, .05, "Distance from TSS (kb)", ha="center", fontsize=7)
cb2 = fig.colorbar(hm, cax=fig.add_axes([x0 + 6 * (cw + gap) + .004, .125, .009, .14]), extend="max")
cb2.outline.set_visible(False); cb2.ax.tick_params(labelsize=5.8, length=2)
cb2.set_label("Fragments per 10⁶\nnon-Y fragments", fontsize=6.2, labelpad=2)
checks["tss_heatmap"] = {"sites": int(len(row_order)), "vmax_p99": vmax,
                         "centre_to_flank": {s_: float(mats[s_][:, n_bins // 2 - 1:n_bins // 2 + 1].mean()
                                                   / np.r_[mats[s_][:, :2], mats[s_][:, -2:]].mean()) for s_ in SAMPLES}}
save(fig, "ch03_supp_ht1080_atac_qc")

# ------------------------------------------------------- sensitivity figure
da = pd.read_csv(DATA / "differential_peaks.csv")
summ = pd.read_csv(DATA / "analysis_summary.csv").set_index("arm")
csaw = pd.read_csv(DATA / "csaw_region_counts.csv").set_index("direction").regions
fig = plt.figure(figsize=(183 / 25.4, 160 / 25.4))
axes = [fig.add_axes([.10, .55, .33, .33]), fig.add_axes([.60, .55, .33, .33]),
        fig.add_axes([.24, .09, .25, .33]), fig.add_axes([.73, .09, .25, .33])]
sens = {}
for ax, key, name in [(axes[0], "background", "Background-bin TMM"), (axes[1], "passage", "Passage-adjusted")]:
    for dkey in ["ns", "lost", "gained"]:
        q = da[da.direction == dkey]
        ax.scatter(q.logFC, q[f"{key}_logFC"], s=1.5 if dkey == "ns" else 4, c=DIR[dkey], lw=0,
                   alpha=.4 if dkey == "ns" else .85, rasterized=True, zorder=1 if dkey == "ns" else 2)
    ax.plot([-5, 5], [-5, 5], color="#888888", lw=.5, ls=(0, (4, 3)), zorder=0)
    rho = spearmanr(da.logFC, da[f"{key}_logFC"]).statistic
    row = summ.loc[f"{key}_TMM"]
    assert abs(rho - row.spearman_logFC_vs_primary) < 1e-9
    ax.set(xlim=(-4.5, 4.5), ylim=(-4.5, 4.5), xlabel="Primary log₂FC", ylabel=f"{name} log₂FC")
    ax.set_aspect("equal")
    ax.text(.03, .97, f"Spearman ρ = {rho:.4f}\n{int(row.gained):,} gained, {int(row.lost):,} lost\n"
            f"{int(row.shared_significant_same_direction)} of 758 primary calls retained",
            transform=ax.transAxes, va="top", fontsize=6.2, color=GREY, linespacing=1.35)
    ax.set_title(name, fontsize=7.5, weight="semibold", pad=4)
    sens[key] = dict(rho=float(rho), gained=int(row.gained), lost=int(row.lost))
letter(fig, axes[0], "a", dx=-.085, dy=.04); letter(fig, axes[1], "b", dx=-.085, dy=.04)

c = axes[2]
arms = [("Primary\n(edgeR, peaks)", int(summ.loc["primary_TMM", "lost"]), int(summ.loc["primary_TMM", "gained"])),
        ("Background-bin\nTMM", int(summ.loc["background_TMM", "lost"]), int(summ.loc["background_TMM", "gained"])),
        ("Passage-\nadjusted", int(summ.loc["passage_TMM", "lost"]), int(summ.loc["passage_TMM", "gained"])),
        ("csaw windows\n(merged regions)", int(csaw["lost"]), int(csaw["gained"]))]
for i, (name, l, g) in enumerate(arms):
    c.barh(i, -np.log10(l), color=DIR["lost"], height=.6); c.barh(i, np.log10(g), color=DIR["gained"], height=.6)
    c.text(-np.log10(l) - .08, i, f"{l:,}", ha="right", va="center", fontsize=6.2, color=DIR["lost"])
    c.text(np.log10(g) + .08, i, f"{g:,}", ha="left", va="center", fontsize=6.2, color=DIR["gained"])
c.axvline(0, color="#222222", lw=.6)
c.set_yticks(range(len(arms)), [a_[0] for a_ in arms], fontsize=6.4)
c.set(xlim=(-5.6, 5.6), ylim=(len(arms) - .4, -.6), xlabel="Regions at FDR < 0.05 (log₁₀ scale)")
c.set_xticks([-4, -2, 0, 2, 4], ["10⁴", "10²", "1", "10²", "10⁴"])
c.spines["left"].set_visible(False); c.tick_params(axis="y", length=0)
c.text(.25, 1.02, "lost", transform=c.transAxes, ha="center", fontsize=6.8, color=DIR["lost"])
c.text(.75, 1.02, "gained", transform=c.transAxes, ha="center", fontsize=6.8, color=DIR["gained"])
letter(fig, c, "c", dx=-.20, dy=.04)
h = [Line2D([], [], ls="", marker="o", markersize=3.5, color=DIR[k], label=l) for k, l in
     [("gained", "Primary gained"), ("lost", "Primary lost"), ("ns", "Primary not significant")]]
fig.legend(handles=h, loc="upper center", bbox_to_anchor=(.5, .995), ncol=3, fontsize=7, columnspacing=1.6)

# d  how many primary FDR < 0.05 calls each check retains (moved from the main figure)
e = axes[3]
tiers = [("Primary FDR < 0.05", da.direction),
         ("TREAT |log₂FC| > 0.5", da.direction.where(da["TREAT_FDR_lfc0.5"] < .05, "ns")),
         ("csaw region agrees", da.direction.where(da.csaw_significant_same_direction, "ns")),
         ("Same call, background-bin TMM", da.direction.where(da.background_direction == da.direction, "ns")),
         ("Same call, passage-adjusted", da.direction.where(da.passage_direction == da.direction, "ns"))]
tier_rows = []
for i, (name, series) in enumerate(tiers):
    series = series.where(da.direction != "ns", "ns")
    g, l = int((series == "gained").sum()), int((series == "lost").sum())
    tier_rows.append(dict(tier=name, gained=g, lost=l))
    e.barh(i, -l, color=DIR["lost"], height=.62); e.barh(i, g, color=DIR["gained"], height=.62)
    e.text(-l - 12, i, f"{l}", ha="right", va="center", fontsize=6.2, color=DIR["lost"])
    e.text(g + 12, i, f"{g}", ha="left", va="center", fontsize=6.2, color=DIR["gained"])
e.axvline(0, color="#222222", lw=.6)
e.set_yticks(range(len(tiers)), [t[0] for t in tiers], fontsize=6.4)
e.set(xlim=(-830, 260), ylim=(len(tiers) - .4, -.6), xlabel="Primary calls retained (lost | gained)")
e.set_xticks([-600, -300, 0, 200], ["600", "300", "0", "200"])
e.spines["left"].set_visible(False); e.tick_params(axis="y", length=0)
pd.DataFrame(tier_rows).to_csv(DATA / "evidence_tier_counts.csv", index=False)
letter(fig, e, "d", dx=-.225, dy=.04)
checks["tiers"] = tier_rows
save(fig, "ch03_supp_ht1080_atac_sensitivity")
checks["sensitivity"] = sens; checks["arms"] = arms
json.dump(checks, open(DATA / "supp_figure_checks.json", "w"), indent=1, default=str)
print("written supp figures")
