"""Main figure: expression responses by baseline promoter-accessibility class (HT1080).

Source: ATACseq/SEPT_08_2026_REPROCESS/code/86_addback_by_atac_class.py (outputs in
integration/thesis_RNA_SEPT11/addback_by_atac_class/). Nothing is refitted; tables are frozen into
source_data/ht1080_addback_by_atac_class/.
a = per-gene expression differences, closed (58) vs matched non-closed DE-down genes (116), for the direct population contrasts;
b = OLS 'closed' coefficient (baseline log2FC + AveExpr adjusted, HC3) across all eight contrasts.
The former panel c remains in historical exports; no models are refitted.
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
DATA = FIG / "source_data/ht1080_addback_by_atac_class"
SRC = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS/integration/thesis_RNA_SEPT11/addback_by_atac_class"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM = "ch03_main_ht1080_addback_by_atac_class"
LOWER, GREY, MATCH = "#3C5488", "#555555", "#A8A8A8"

# Read frozen tables only; preserve the original analysis and source manifest.
res = pd.read_csv(DATA / "restoration_closed_vs_not.csv").set_index("contrast")
dn = pd.read_csv(DATA / "down_genes_addback.csv")
pairs = pd.read_csv(DATA / "matched_pairs.csv")
matched_names = set(pairs.matched_name)

fig = plt.figure(figsize=(10.0, 4.5))


def letter(x, y, ch):
    fig.text(x, y, ch, fontsize=12, weight="semibold", va="bottom")


def fmt_p(p):
    return f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}".replace("e-0", "e-")


# ---- a ----
PANEL_A = [("UTY vs untransduced", "UTY vs\nuntransduced"), ("UTY+KDM5D vs untransduced", "UTY+KDM5D vs\nuntransduced"),
           ("UTY vs KDM5D", "UTY vs\nKDM5D"), ("KDM5D vs untransduced (no UTY difference)", "KDM5D vs\nuntransduced")]
a = fig.add_axes([.07, .15, .43, .72])
rng = np.random.default_rng(1)
out = []
for i, (key, lab) in enumerate(PANEL_A):
    col = f"ab_{key}"
    cl = dn[dn.closed == 1][col]; mt = dn[dn.gene_name.isin(matched_names) & (dn.closed == 0)][col]
    for off, vals, c, fc in [(-.19, cl, LOWER, LOWER), (.19, mt, MATCH, MATCH)]:
        a.scatter(i + off + rng.uniform(-.07, .07, len(vals)), vals, s=6, color=fc, lw=0, alpha=.75, zorder=2, rasterized=True)
        q1, q2, q3 = np.percentile(vals, [25, 50, 75])
        a.plot([i + off - .13, i + off + .13], [q2, q2], color="black", lw=1.2, zorder=4)
        a.plot([i + off, i + off], [q1, q3], color="black", lw=.7, zorder=4)
    r = res.loc[key]
    a.text(i, 8.3, f"Δ = {r.ols_closed_coef:+.2f}\nFDR {fmt_p(r.ols_closed_BH)}", ha="center", va="top", fontsize=8,
           color="black" if r.ols_closed_BH < .05 else GREY, linespacing=1.1)
    out.append(dict(contrast=key, closed_median=cl.median(), matched_median=mt.median(), ols_closed=r.ols_closed_coef, BH=r.ols_closed_BH))

vals_all = pd.concat([dn[f"ab_{k}"] for k, _ in PANEL_A])
assert vals_all.min() > -4.6 and vals_all.max() < 7.3, "panel a would clip"
a.axhline(0, color="#999999", lw=.6, ls=(0, (4, 3)))
a.axvspan(2.5, 3.5, color="#F4F4F4", zorder=0)
a.set_xticks(range(len(PANEL_A)), [l for _, l in PANEL_A], fontsize=8)
a.set_xlim(-.55, len(PANEL_A) - .45); a.set_ylim(-4.6, 8.5)
a.set_ylabel("Expression difference (log2)", fontsize=9)
a.tick_params(axis="y", labelsize=8); a.tick_params(axis="x", length=0)
a.text(3, 7.3, "no UTY\ndifference", ha="center", va="top", fontsize=8, color=GREY, linespacing=1.05)
a.legend(handles=[Line2D([], [], ls="", marker="o", ms=4, color=LOWER, label="Promoter loss (58)"),
                  Line2D([], [], ls="", marker="o", ms=4, color=MATCH, label="Matched controls (116)")],
         loc="upper left", bbox_to_anchor=(0, 1.1), ncol=2, fontsize=8, frameon=False, handletextpad=.2, borderpad=0)
letter(.01, .92, "a")

# ---- b: adjusted estimates, vector comparisons first ----
ORDER = [("UTY", "UTY vs vector"), ("UTY+KDM5D", "UTY+KDM5D vs vector"),
         ("KDM5D", "KDM5D vs vector"), ("Untransduced Y-depleted", "Untransduced vs vector"),
         ("UTY vs untransduced", "UTY vs untransduced"), ("UTY+KDM5D vs untransduced", "UTY+KDM5D vs untransduced"),
         ("UTY vs KDM5D", "UTY vs KDM5D"),
         ("KDM5D vs untransduced (no UTY difference)", "KDM5D vs untransduced")]


def forest(ax, est, lo, hi, bh, xlabel, xlim):
    for i, (key, _) in enumerate(ORDER):
        sig = bh[key] < .05
        ax.plot([lo[key], hi[key]], [i, i], color=LOWER, lw=.9)
        ax.scatter(est[key], i, s=16, color=LOWER if sig else "white", edgecolor=LOWER, lw=.8, zorder=3)
        ax.text(xlim[1], i, fmt_p(bh[key]), ha="right", va="center", fontsize=8, color="black" if sig else GREY)
    ax.axvline(0, color="#888888", lw=.6, ls=(0, (4, 3)))
    ax.axhspan(-.6, 3.5, color="#F7F7F7", zorder=0)
    ax.set_ylim(len(ORDER) - .4, -.6); ax.set_xlim(*xlim)
    ax.set_yticks(range(len(ORDER)), [l for _, l in ORDER], fontsize=8)
    ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", labelsize=8); ax.spines["left"].set_visible(False)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.text(xlim[1], -.95, "FDR", ha="right", va="bottom", fontsize=8, color=GREY)


b = fig.add_axes([.755, .15, .22, .72])
forest(b, res.ols_closed_coef, res.ols_closed_CI_low, res.ols_closed_CI_high, res.ols_closed_BH,
       "Promoter loss − no loss, adjusted (log2)", (-1.0, 1.6))
letter(.545, .92, "b")
res.loc[[k for k, _ in ORDER], ["ols_closed_coef", "ols_closed_CI_low", "ols_closed_CI_high", "ols_closed_BH"]].to_csv(DATA / "panel_b_final_display.csv")

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=200)
plt.close(fig)
print("written", STEM)
