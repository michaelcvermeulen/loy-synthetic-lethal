"""Supplementary TSS-centred accessibility at promoters of RNA-defined gene groups (HT1080).

Companion to ch03_supp_ht1080_atac_differential_heatmaps (region-centred). Inputs
(source_data/ht1080_atac_accessibility), from prepare_ht1080_atac_tss_groups.py:
  tss_groups_genes.csv, tss_groups_tests.csv, tss_groups_control_draws.csv,
  tss_groups_ratio_by_passage.csv; coverage from tss_transcript_matrix.npz.
Profiles and heatmaps show the first (seeded) control draw; P values and median
differences on the plot summarise 200 independent control draws.
"""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_atac_accessibility"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM = "ch03_supp_ht1080_atac_tss_groups"
COL = {"CTRL": "#3C5488", "LOY": "#E64B35"}
LABEL = {"CTRL": "sgCTRL", "LOY": "Y-depleted"}
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
GROUPS = [("down", "Down-regulated"), ("down_no_lost_peak", "Down-regulated,\nno called lost peak"),
          ("up", "Up-regulated"), ("promoter_lost", "Lost promoter peak\n(ATAC-selected)")]
TARGET, CONTROL = "#222222", "#9A9A9A"
PSEUDO, RATIO_LIM = .1, 2.0
DIVERGE = LinearSegmentedColormap.from_list("lost_gained", ["#3C5488", "#FFFFFF", "#E64B35"])

genes = pd.read_csv(DATA / "tss_groups_genes.csv", dtype={"entrez_id": str})
tests = pd.read_csv(DATA / "tss_groups_tests.csv").set_index("group")
draws = pd.read_csv(DATA / "tss_groups_control_draws.csv").set_index("group")
passage = pd.read_csv(DATA / "tss_groups_ratio_by_passage.csv")
z = np.load(DATA / "tss_transcript_matrix.npz")
mats = {s: z[s] for s in SAMPLES}
nb = mats[SAMPLES[0]].shape[1]
half_kb = 3
xs = (np.arange(nb) + .5) * 2 * half_kb / nb - half_kb
C = slice(nb // 2 - 2, nb // 2 + 2)
ctrl = np.mean([mats[s] for s in SAMPLES[:3]], axis=0)
loy = np.mean([mats[s] for s in SAMPLES[3:]], axis=0)
ratio = np.log2((loy + PSEUDO) / (ctrl + PSEUDO))


def rows(group):
    d = genes[genes.group == group]
    return d.row.to_numpy()[np.argsort(-d.pooled_centre.to_numpy(), kind="stable")]


def fmt_p(p):
    if p >= 1e-3:
        return f"P = {p:.2g}"
    m, e = f"{p:.1e}".split("e")
    return f"P = {m}×10$^{{{int(e)}}}$"


fig = plt.figure(figsize=(183 / 25.4, 200 / 25.4))

# Average profiles are retained in main Figure F09d; keep central values for checks.
pw, px0, pstep = .19, .075, .232
stats = {}
for group, _ in GROUPS:
    for grp in [group, group+'_ref']:
        idx=rows(grp)
        stats[grp]={cond:float(np.mean([mats[key][idx].mean(0) for key in keys],axis=0)[C].mean())
                    for cond,keys in [('CTRL',SAMPLES[:3]),('LOY',SAMPLES[3:])]}

# b  heatmaps: pooled sgCTRL, pooled Y-depleted and log2 ratio; genes above their controls
vmax = float(np.percentile(np.concatenate([ctrl[genes.row].ravel(), loy[genes.row].ravel()]), 99))
cw, cgap, bgap = .052, .005, .232
top, hb, hgap = .625, .13, .012
fig.text(.015, top + .075, "a", fontsize=12, weight="semibold")
for gi, (g, name) in enumerate(GROUPS):
    x0 = px0 + gi * bgap
    fig.text(x0+.085, top+.042, name+f"\nn = {int(tests.loc[g, 'targets'])} genes",ha="center",va="bottom",fontsize=6,weight="semibold")
    for bi, grp in enumerate([g, f"{g}_ref"]):
        y = top - (bi + 1) * hb - bi * hgap
        idx = rows(grp)
        for k, (lab, m) in enumerate([("sgCTRL", ctrl), ("Y-depl.", loy), ("log₂ L/C", ratio)]):
            ax = fig.add_axes([x0 + k * (cw + cgap) + (.006 if k == 2 else 0), y, cw, hb])
            if k == 2:
                im_r = ax.imshow(m[idx], aspect="auto", cmap=DIVERGE, vmin=-RATIO_LIM, vmax=RATIO_LIM,
                                 extent=(-3, 3, len(idx), 0), interpolation="antialiased")
            else:
                im = ax.imshow(m[idx], aspect="auto", cmap="magma", vmin=0, vmax=vmax,
                               extent=(-3, 3, len(idx), 0), interpolation="antialiased")
            ax.set_yticks([]); ax.set_xticks([])
            for sp in ax.spines.values(): sp.set_visible(False)
            if bi == 0:
                ax.set_title(lab, fontsize=6.2, pad=2,
                             color={0: COL["CTRL"], 1: COL["LOY"], 2: "#222222"}[k])
            else:
                ax.set_xticks([-3, 0, 3], ["−3", "0", "+3"], fontsize=5.6)
                for t_, al in zip(ax.get_xticklabels(), ["left", "center", "right"]):
                    t_.set_ha(al)
                ax.tick_params(axis="x", length=2, pad=1)
            if k == 0:
                ax.set_ylabel("Genes" if bi == 0 else "Controls", fontsize=6.2, labelpad=2,
                              color=TARGET if bi == 0 else "#666666")
fig.text(px0 + 1.5 * bgap + pw / 2, top - 2 * hb - hgap - .036, "Distance from TSS (kb)", ha="center", fontsize=6.8)
cb = fig.colorbar(im, cax=fig.add_axes([.958, top - hb + .005, .008, hb - .01]), extend="max")
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.6, length=2)
cb.set_label("CPM", fontsize=5.9, labelpad=1)
cb2 = fig.colorbar(im_r, cax=fig.add_axes([.958, top - 2 * hb - hgap + .005, .008, hb - .01]), extend="both")
cb2.outline.set_visible(False); cb2.ax.tick_params(labelsize=5.6, length=2)
cb2.set_label("log₂ L/C", fontsize=5.9, labelpad=1)

# c  per-gene central log2 ratio, genes vs matched controls
cx = fig.add_axes([.085, .065, .50, .205])
rng = np.random.default_rng(1)
for gi, (g, name) in enumerate(GROUPS):
    for off, grp, colour in [(-.19, g, TARGET), (.19, f"{g}_ref", CONTROL)]:
        v = genes.loc[genes.group == grp, "log2_ratio"].to_numpy()
        pos = gi + off
        cx.scatter(pos + rng.uniform(-.1, .1, len(v)), np.clip(v, -3.4, 2.4), s=1.2, color=colour, alpha=.35,
                   lw=0, rasterized=True)
        cx.boxplot(v, positions=[pos], widths=.28, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="none", edgecolor=colour, lw=.9), medianprops=dict(color=colour, lw=1.4),
                   whiskerprops=dict(color=colour, lw=.8), capprops=dict(color=colour, lw=.8))
    d = draws.loc[g]
    lab = (f"Δ {d.median_difference_median:+.2f}\n({d.median_difference_q025:+.2f} to "
           f"{d.median_difference_q975:+.2f})").replace("-", "−")
    cx.text(gi, 2.55, f"{lab}\n{fmt_p(d.MW_P_median)}",
            ha="center", va="bottom", fontsize=5.7, linespacing=1.2)
cx.axhline(0, color="#888888", lw=.6, ls=(0, (4, 3)), zorder=0)
cx.set_xticks(range(4), [n for _, n in GROUPS], fontsize=6.1)
cx.set_xlim(-.55, 3.55)
cx.set_ylim(-3.5, 3.9)
cx.set_yticks([-3, -2, -1, 0, 1, 2])
cx.set_ylabel("Central 200 bp log₂ (Y-depleted / sgCTRL)", fontsize=6.6)
cx.legend(handles=[Patch(facecolor="none", edgecolor=TARGET, label="Genes"),
                   Patch(facecolor="none", edgecolor=CONTROL, label="Matched not-DE controls")],
          loc="lower left", fontsize=6, frameon=False, ncol=2, handlelength=1.2, borderaxespad=.2)
fig.text(.015, .29, "b", fontsize=12, weight="semibold")

# d  group ratio of central coverage in each passage pair
dx = fig.add_axes([.69, .105, .28, .165])
pas = ["p11", "p12", "p21"]
for gi, (g, name) in enumerate(GROUPS):
    for off, grp, colour, mk in [(-.14, g, TARGET, "o"), (.14, f"{g}_ref", CONTROL, "s")]:
        vals = [passage[(passage.group == grp) & (passage.passage == p)].ratio.iloc[0] for p in pas]
        dx.plot([gi + off] * 3, vals, color=colour, lw=.6, zorder=1)
        dx.scatter([gi + off] * 3, vals, s=11, marker=mk, color=colour, zorder=2, lw=0)
dx.axhline(1, color="#888888", lw=.6, ls=(0, (4, 3)), zorder=0)
dx.set_xticks(range(4), ["Down", "Down, no\ncalled loss", "Up", "Lost\npromoter"], fontsize=5.9)
dx.set_xlim(-.5, 3.5)
dx.set_ylim(.3, 1.2)
dx.set_ylabel("Central coverage ratio,\nY-depleted / sgCTRL", fontsize=6.4)
dx.tick_params(axis="y", labelsize=6)
dx.legend(handles=[Line2D([], [], ls="", marker="o", markersize=3.5, color=TARGET, label="Genes"),
                   Line2D([], [], ls="", marker="s", markersize=3.5, color=CONTROL, label="Controls")],
          loc="lower left", fontsize=5.9, frameon=False, handletextpad=.2, borderaxespad=.2)
dx.set_title("One point per passage pair (p11, p12, p21)", fontsize=6.1, pad=3)
fig.text(.625, .29, "c", fontsize=12, weight="semibold")

# Remove the vacated top row while preserving remaining panel sizes.
fig.set_size_inches(183/25.4,145/25.4)
scale=200/145
for ax in fig.axes:
    box=ax.get_position()
    ax.set_position([box.x0,box.y0*scale,box.width,box.height*scale])
for text in fig.texts:
    x,y=text.get_position();text.set_position((x,y*scale))

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target_, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=.04); stage.replace(target_)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=180, bbox_inches="tight", pad_inches=.04)
plt.close(fig)
json.dump({"centre_signal_plotted_draw": stats, "vmax_p99": vmax, "ratio_pseudocount": PSEUDO,
           "ratio_limits": RATIO_LIM, "per_gene_points_clipped": int(((genes.log2_ratio < -3.4) | (genes.log2_ratio > 2.4)).sum())},
          open(DATA / "tss_groups_checks.json", "w"), indent=1)
print("written", STEM)
