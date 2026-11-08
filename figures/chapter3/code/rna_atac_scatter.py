"""Shared drawing for peak-level RNA-ATAC link scatters (main integration figure and supplement).

Categories come from prepare_ht1080_rna_atac_peak_scatter.R: neither, peak_only, RNA_only,
concordant, opposite. Keeping the styling here makes the promoter panel identical in both figures.
"""
import numpy as np
from matplotlib.lines import Line2D

LOWER, HIGHER, PEAK_ONLY = "#3C5488", "#E64B35", "#9A9A9A"


def draw_link_points(ax, d):
    for cat, kw in [("neither", dict(s=1.2, c="#D2D2D2", zorder=1)),
                    ("peak_only", dict(s=3.5, c=PEAK_ONLY, alpha=.35, zorder=1.5))]:
        q = d[d.category == cat]
        ax.scatter(q.rna_logFC, q.logFC, lw=0, rasterized=True, **kw)
    q = d[d.category == "RNA_only"]
    ax.scatter(q.rna_logFC, q.logFC, s=3.5, c=np.where(q.rna_class == "up", HIGHER, LOWER), alpha=.35, lw=0,
               rasterized=True, zorder=2)
    q = d[d.category == "concordant"]
    ax.scatter(q.rna_logFC, q.logFC, s=7, c=np.where(q.rna_class == "up", HIGHER, LOWER), edgecolors="#111111",
               linewidths=.3, zorder=4)
    q = d[d.category == "opposite"]
    ax.scatter(q.rna_logFC, q.logFC, s=8, marker="D", c="#222222", edgecolors="white", linewidths=.3, zorder=5)
    ax.axhline(0, color="#999999", lw=.5, zorder=0); ax.axvline(0, color="#999999", lw=.5, zorder=0)


def link_legend_handles():
    return [Line2D([], [], ls="", marker="o", markersize=3, color="#D2D2D2", label="Neither"),
            Line2D([], [], ls="", marker="o", markersize=3.5, color=PEAK_ONLY, alpha=.6, label="Peak only"),
            Line2D([], [], ls="", marker="o", markersize=3.5, color=LOWER, alpha=.5, label="RNA DE only, lower"),
            Line2D([], [], ls="", marker="o", markersize=3.5, color=HIGHER, alpha=.5, label="RNA DE only, higher"),
            Line2D([], [], ls="", marker="o", markersize=3.4, markerfacecolor=LOWER, markeredgecolor="#111111",
                   label="Concordant, lower"),
            Line2D([], [], ls="", marker="o", markersize=3.4, markerfacecolor=HIGHER, markeredgecolor="#111111",
                   label="Concordant, higher"),
            Line2D([], [], ls="", marker="D", markersize=3, markerfacecolor="#222222", markeredgecolor="white",
                   label="Opposite")]
