"""Supplementary heatmaps of coverage and motif content at HT1080 differential regions.

Inputs (source_data/ht1080_atac_accessibility):
  region_heatmap_regions.csv, region_heatmap_motif_CMH.csv  <- prepare_ht1080_atac_region_sets.R
  region_heatmap_matrix.npz                                  <- extract_ht1080_atac_region_matrix.py (WSL)
  region_heatmap_homer_stats.csv                             <- prepare_ht1080_atac_region_heatmaps.py
Each differential group is shown with its own GC + accessibility-matched controls
(3 per target). Panel b shows the stratified CMH tests of motif presence; the
seven motifs were chosen from the HOMER screen.
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
STEM = "ch03_supp_ht1080_atac_differential_heatmaps"
COL = {"CTRL": "#3C5488", "LOY": "#E64B35"}
LABEL = {"CTRL": "sgCTRL", "LOY": "Y-depleted"}
DIRC = {"lost": "#3C5488", "gained": "#E64B35"}
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
GROUPS = [("lost", "Lost", "#3C5488"), ("lost_ref", "Matched\ncontrols", "#555555"),
          ("gained", "Gained", "#E64B35"), ("gained_ref", "Matched\ncontrols", "#555555")]
FAMILIES = ["TFAP2", "EBF", "ZIC", "TEAD", "RUNX", "GATA", "SOX6"]
PSEUDO, RATIO_LIM = .1, 2.0
short = lambda s: f"{'C' if s.startswith('CTRL') else 'L'} p{s.split('_')[1]}"
DIVERGE = LinearSegmentedColormap.from_list("lost_gained", ["#3C5488", "#FFFFFF", "#E64B35"])
PRESENCE = LinearSegmentedColormap.from_list("presence", ["#F2F2F2", "#222222"])

regions = pd.read_csv(DATA / "region_heatmap_regions.csv")
z = np.load(DATA / "region_heatmap_matrix.npz")
assert list(z["peak_id"]) == list(regions.peak_id)
cmh = pd.read_csv(DATA / "region_heatmap_motif_CMH.csv")
homer = pd.read_csv(DATA / "region_heatmap_homer_stats.csv")
REP_NAME = homer.drop_duplicates("family").set_index("family").motif_name.str.upper()
mats = {s: z[s] for s in SAMPLES}
nb = mats[SAMPLES[0]].shape[1]
xs = (np.arange(nb) + .5) * 4 / nb - 2
ctrl = np.mean([mats[s] for s in SAMPLES[:3]], axis=0)
loy = np.mean([mats[s] for s in SAMPLES[3:]], axis=0)
ratio = np.log2((loy + PSEUDO) / (ctrl + PSEUDO))
vmax = float(np.percentile(np.concatenate([mats[s].ravel() for s in SAMPLES]), 99))

# Row order: lost most-negative first, gained most-positive first, controls by mean signal.
order = {}
for g, _, _ in GROUPS:
    idx = np.flatnonzero(regions.group == g)
    lfc = regions.logFC.to_numpy()[idx]
    if g == "lost":
        idx = idx[np.argsort(lfc, kind="stable")]
    elif g == "gained":
        idx = idx[np.argsort(-lfc, kind="stable")]
    else:
        idx = idx[np.argsort(-(ctrl[idx].mean(1) + loy[idx].mean(1)), kind="stable")]
    order[g] = idx

fig = plt.figure(figsize=(183 / 25.4, 205 / 25.4))

# a  mean profiles per group
prof_axes = [fig.add_axes([.07 + i * .153, .775, .128, .14]) for i in range(4)]
ymax, stats = 0, {}
for ax, (g, name, colour) in zip(prof_axes, GROUPS):
    idx = order[g]
    for cond, keys in [("CTRL", SAMPLES[:3]), ("LOY", SAMPLES[3:])]:
        per_rep = np.array([mats[s][idx].mean(0) for s in keys])
        ax.fill_between(xs, per_rep.min(0), per_rep.max(0), color=COL[cond], alpha=.2, lw=0)
        ax.plot(xs, per_rep.mean(0), color=COL[cond], lw=1)
        ymax = max(ymax, per_rep.max())
        stats.setdefault(g, {})[cond] = float(per_rep.mean(0)[nb // 2 - 2:nb // 2 + 2].mean())
    title = {"lost": "Lost", "lost_ref": "Controls for lost", "gained": "Gained",
             "gained_ref": "Controls for gained"}[g]
    ax.set_title(f"{title}\nn = {len(idx):,}", fontsize=6.8, weight="semibold", color=colour, pad=3, linespacing=1.15)
    ax.set_xticks([-2, 0, 2], ["−2", "0", "+2"], fontsize=6.3)
    ax.tick_params(axis="y", labelsize=6.3)
    ax.text(.97, .95, f"{stats[g]['LOY'] / stats[g]['CTRL']:.2f}×", transform=ax.transAxes, ha="right",
            va="top", fontsize=6.3, color="#555555")
for ax in prof_axes:
    ax.set_ylim(0, ymax * 1.08)
for ax in prof_axes[1:]:
    ax.set_yticklabels([])
prof_axes[0].set_ylabel("Mean fragments per 10⁶\nnon-Y fragments", fontsize=6.8)
fig.text(.07 + 1.5 * .153 + .064, .735, "Distance from region centre (kb)", ha="center", fontsize=6.8)
fig.text(.02, .935, "a", fontsize=12, weight="semibold")

# b  CMH odds ratios (targets vs their own matched controls), log scale
bx = fig.add_axes([.765, .755, .20, .16])
for i, fam in enumerate(FAMILIES):
    for dy, direction in [(-.17, "lost"), (.17, "gained")]:
        r = cmh[(cmh.family == fam) & (cmh.direction == direction)].iloc[0]
        sig = r.BH_FDR < .05
        bx.plot([r.CI_low, r.CI_high], [i + dy, i + dy], color=DIRC[direction], lw=.9, solid_capstyle="butt")
        bx.scatter(r.MH_odds_ratio, i + dy, s=16, marker="o", zorder=3, linewidths=.8,
                   facecolors=DIRC[direction] if sig else "white", edgecolors=DIRC[direction])
bx.axvline(1, color="#888888", lw=.6, ls=(0, (4, 3)), zorder=0)
bx.set_xscale("log")
bx.set(xlim=(.2, 8), ylim=(len(FAMILIES) - .4, -.6))
bx.set_xticks([.25, .5, 1, 2, 4], ["0.25", "0.5", "1", "2", "4"], fontsize=6)
bx.minorticks_off()
bx.set_yticks(range(len(FAMILIES)), [REP_NAME[f] for f in FAMILIES], fontsize=6.4, style="normal")
bx.tick_params(axis="y", length=0)
bx.spines["left"].set_visible(False)
bx.set_xlabel("Odds ratio vs matched controls\n(CMH, 95% CI)", fontsize=6.5, linespacing=1.15)
fig.text(.695, .935, "b", fontsize=12, weight="semibold")
bh = [Line2D([], [], ls="", marker="o", markersize=4, color=DIRC["lost"], label="Lost"),
      Line2D([], [], ls="", marker="o", markersize=4, color=DIRC["gained"], label="Gained"),
      Line2D([], [], ls="", marker="o", markersize=4, markerfacecolor="white", markeredgecolor="#555555",
             label="BH FDR ≥ 0.05")]
bx.legend(handles=bh, loc="lower center", bbox_to_anchor=(.42, 1.0), ncol=3, fontsize=5.9, handletextpad=.2,
          columnspacing=.8, borderpad=0)

# c  heatmaps + ratio + motif strips; each target block sits directly above its controls
heights = {"lost": .145, "lost_ref": .145, "gained": .07, "gained_ref": .07}
gap_after = {"lost": 0, "lost_ref": .045, "gained": 0, "gained_ref": 0}
top = .595
x0, cw, gap = .125, .072, .007
sw, sgap = .016, .003
strip_x0 = x0 + 7 * (cw + gap) + .012 + .016
fig.text(.02, top + .045, "c", fontsize=12, weight="semibold")
fig.text(strip_x0 + 3.5 * (sw + sgap) - sgap / 2, top + .045, "Motif present\n(% of regions)",
         ha="center", va="bottom", fontsize=6.2, linespacing=1.15)
y = top
for g, name, colour in GROUPS:
    h = heights[g]
    y -= h
    idx = order[g]
    for k, s in enumerate(SAMPLES + ["ratio"]):
        ax = fig.add_axes([x0 + k * (cw + gap) + (.012 if s == "ratio" else 0), y, cw, h])
        if s == "ratio":
            im_r = ax.imshow(ratio[idx], aspect="auto", cmap=DIVERGE, vmin=-RATIO_LIM, vmax=RATIO_LIM,
                             extent=(-2, 2, len(idx), 0), interpolation="antialiased")
        else:
            im = ax.imshow(mats[s][idx], aspect="auto", cmap="magma", vmin=0, vmax=vmax,
                           extent=(-2, 2, len(idx), 0), interpolation="antialiased")
        ax.set_yticks([]); ax.set_xticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        if g == "lost":
            ax.set_title("log₂ L/C" if s == "ratio" else short(s), fontsize=6.8, pad=3,
                         color="#222222" if s == "ratio" else COL["CTRL" if s.startswith("CTRL") else "LOY"])
        if g == "gained_ref":
            ax.set_xticks([-2, 0, 2], ["−2", "0", "+2"], fontsize=5.8)
            for t_, al in zip(ax.get_xticklabels(), ["left", "center", "right"]):
                t_.set_ha(al)
            ax.tick_params(axis="x", length=2, pad=1)
    direction = g.replace("_ref", "")
    for m_i, fam in enumerate(FAMILIES):
        sx = fig.add_axes([strip_x0 + m_i * (sw + sgap), y, sw, h])
        present = regions[fam].to_numpy()[idx].astype(float)[:, None]
        sx.imshow(present, aspect="auto", cmap=PRESENCE, vmin=0, vmax=1, interpolation="nearest",
                  extent=(0, 1, len(idx), 0))
        sx.set_xticks([]); sx.set_yticks([])
        for sp in sx.spines.values(): sp.set_visible(False)
        if g == "lost":
            sx.set_title(REP_NAME[fam], fontsize=5.9, rotation=90, style="normal", pad=3)
        # Target percentages are bold where the CMH test vs their own controls has BH FDR < 0.05.
        bold = g == direction and cmh[(cmh.family == fam) & (cmh.direction == direction)].BH_FDR.iloc[0] < .05
        sx.text(.5, -.02, f"{100 * present.mean():.0f}", transform=sx.transAxes, ha="center", va="top",
                fontsize=5.3, color="#222222" if bold else "#777777", weight="semibold" if bold else "normal")
    fig.text(x0 - .012, y + h / 2, f"{name}\nn = {len(idx):,}", ha="right", va="center", fontsize=6.6,
             color=colour, weight="semibold", linespacing=1.2)
    y -= gap_after[g] + (.012 if g in ("lost", "gained") else 0)
fig.text(x0 + 3 * (cw + gap), y - .025, "Distance from region centre (kb)", ha="center", fontsize=7)
cb = fig.colorbar(im, cax=fig.add_axes([.875, .40, .01, .15]), extend="max")
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.8, length=2)
cb.set_label("Fragments per 10⁶\nnon-Y fragments", fontsize=6.3)
cb2 = fig.colorbar(im_r, cax=fig.add_axes([.875, .14, .01, .15]), extend="both")
cb2.outline.set_visible(False); cb2.ax.tick_params(labelsize=5.8, length=2)
cb2.set_label("log₂ (Y-depleted / sgCTRL)", fontsize=6.3)

handles = [Line2D([], [], color=COL[c], lw=1.2, label=f"{LABEL[c]} mean") for c in ["CTRL", "LOY"]]
handles.append(Patch(facecolor="#999999", alpha=.3, label="Range across three libraries"))
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.38, .995), ncol=3, fontsize=7, columnspacing=1.8)

# Remove 8 mm from the gap above panel c without shrinking any data panel.
old_height, new_height = 205., 197.
for axis in fig.axes:
    pos = axis.get_position()
    shift = 8. if pos.y0 >= .70 else 0.
    axis.set_position([pos.x0, (pos.y0*old_height-shift)/new_height,
                       pos.width, pos.height*old_height/new_height])
for text in fig.texts:
    tx, ty = text.get_position()
    text.set_position((tx, (ty*old_height-(8. if ty >= .70 else 0.))/new_height))
fig.set_size_inches(183/25.4, new_height/25.4)
(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target_, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target_)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=180)
plt.close(fig)
json.dump({"groups": {g: int(len(order[g])) for g in order}, "centre_signal": stats, "vmax_p99": vmax,
           "ratio_pseudocount": PSEUDO, "ratio_limits": RATIO_LIM,
           "ratio_clipped_fraction": float((np.abs(ratio) > RATIO_LIM).mean())},
          open(DATA / "region_heatmap_checks.json", "w"), indent=1)
print("written", STEM)
