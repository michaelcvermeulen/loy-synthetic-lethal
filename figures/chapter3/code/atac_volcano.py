"""Shared HT1080 ATAC volcano drawing for the main and supplementary figures.

Labels follow a fixed rule applied to FDR < 0.05 regions: the n_p lowest P values,
the n_ext most extreme log2FC in each direction, plus any forced peak IDs.
Label text is the nearest knownGene TSS gene (positional only for non-promoter
regions); bold marks promoter-overlapping regions.
"""
import numpy as np
import matplotlib.patheffects as pe

DIR = {"gained": "#E64B35", "lost": "#3C5488", "ns": "#C6C6C6"}


def select_labels(da, n_p, n_ext, forced=None, include_loc=True):
    sig = da[da.direction != "ns"]
    pool = sig if include_loc else sig[~sig.nearest_symbol.fillna("LOC").str.startswith("LOC")]
    chosen = (set(pool.nsmallest(n_p, "PValue").peak_id) | set(pool.nsmallest(n_ext, "logFC").peak_id)
              | set(pool.nlargest(n_ext, "logFC").peak_id))
    forced = forced or {}
    chosen |= set(forced)
    lab = sig[sig.peak_id.isin(chosen)].copy()
    lab["text"] = lab.peak_id.map(forced).fillna(lab.nearest_symbol)
    lab["promoter"] = lab.annotation == "promoter"
    return lab


def spread(targets, gap, lo, hi, iterations=500):
    """1-D label spreading: keep order (descending), push neighbours apart
    symmetrically until every gap >= `gap`, clamped to [lo, hi]."""
    ys = np.array(targets, dtype=float)
    for _ in range(iterations):
        moved = False
        for k in range(1, len(ys)):
            deficit = gap - (ys[k - 1] - ys[k])
            if deficit > 1e-9:
                ys[k - 1] += deficit / 2; ys[k] -= deficit / 2; moved = True
        ys = np.clip(ys, lo, hi)
        if not moved:
            break
    return ys


def draw_volcano(ax, da, lab, fontsize=5.9, point_scale=1.0, xlim=(-5.4, 3.6), label_gap=.4, headroom=1.1, lost_split=None, right_lo=5.0):
    da = da.assign(y=-np.log10(da.PValue))
    for key in ["ns", "lost", "gained"]:
        q = da[da.direction == key]
        ax.scatter(q.logFC, q.y, s=(2.2 if key == "ns" else 5) * point_scale, c=DIR[key], lw=0,
                   alpha=.5 if key == "ns" else .85, rasterized=True, zorder=1 if key == "ns" else 2)
    treat = da[da["TREAT_FDR_lfc0.5"] < .05]
    ax.scatter(treat.logFC, treat.y, s=11 * point_scale, facecolors="none", edgecolors="#222222",
               linewidths=.45, zorder=3)
    cut = -np.log10(da.loc[da.FDR < .05, "PValue"].max())
    ax.axhline(cut, color="#888888", ls=(0, (4, 3)), lw=.6, zorder=0)
    ax.set(xlim=xlim, ylim=(0, da.y.max() * headroom))
    ax.set_xlabel("log₂ fold change (Y-depleted / sgCTRL)")
    ax.set_ylabel("−log₁₀ P value")
    halo = [pe.withStroke(linewidth=1.6, foreground="white")]
    lab = lab.assign(y=-np.log10(lab.PValue))
    # Deterministic placement: order labels by point height, spread them to a
    # minimum vertical gap, then park each just outside the cloud's edge at that
    # height so leader lines stay short and do not cross.
    ymax = ax.get_ylim()[1]
    renderer = ax.figure.canvas.get_renderer()
    sig_pts = da[da.direction != "ns"]
    sig_px = ax.transData.transform(np.c_[sig_pts.logFC, sig_pts.y])
    # Optional: lost labels whose point lies right of `lost_split` go into the empty
    # gap between the two clouds (above `right_lo`), halving the left-hand stack.
    groups = [("lost", -1, lab[lab.direction == "lost"]), ("gained", 1, lab[lab.direction == "gained"])]
    if lost_split is not None:
        lost = lab[lab.direction == "lost"]
        groups = [("lost", -1, lost[lost.logFC <= lost_split]), ("lost", 1, lost[lost.logFC > lost_split]),
                  groups[1]]
    for key, sign, sub in groups:
        sub = sub.sort_values("y", ascending=False)
        lo = right_lo if (key == "lost" and sign > 0) else label_gap
        ys = spread(sub.y.to_numpy(), label_gap, lo=lo, hi=ymax * .97)
        for r, yt in zip(sub.itertuples(), ys):
            # Clear only nearby points; distant outliers should not push the label away.
            band = da[((da.y - yt).abs() < .3) & ((da.logFC - r.logFC) * sign < .8)]
            band = band[band.direction == key] if key == "lost" and sign > 0 else band
            side = band.logFC.min() if sign < 0 else band.logFC.max()
            side = side if np.isfinite(side) else r.logFC
            xt = min(side, r.logFC) - .18 if sign < 0 else max(side, r.logFC) + .18
            ann = ax.annotate(r.text, (r.logFC, r.y), xytext=(xt, yt), fontsize=fontsize, color=DIR[key],
                              style="italic", weight="semibold" if r.promoter else "normal",
                              ha="right" if sign < 0 else "left", va="center", path_effects=halo, zorder=5,
                              arrowprops=dict(arrowstyle="-", color="#888888", lw=.35, shrinkA=1, shrinkB=1.5))
            # If the label box covers any significant point, step outward past it.
            for _ in range(8):
                box = ann.get_window_extent(renderer).expanded(1.0, 1.0).padded(2.5)
                hit = sig_px[(sig_px[:, 0] > box.x0) & (sig_px[:, 0] < box.x1) &
                             (sig_px[:, 1] > box.y0) & (sig_px[:, 1] < box.y1)]
                if not len(hit):
                    break
                edge_px = hit[:, 0].min() - 4 if sign < 0 else hit[:, 0].max() + 4
                xt = ax.transData.inverted().transform((edge_px, 0))[0]
                ann.set_position((xt, yt))
    n_lost, n_gain = (da.direction == "lost").sum(), (da.direction == "gained").sum()
    ax.text(0, 1.02, f"{n_lost} lost", transform=ax.transAxes, va="bottom", fontsize=7, color=DIR["lost"],
            weight="semibold")
    ax.text(1, 1.02, f"{n_gain} gained", transform=ax.transAxes, va="bottom", ha="right", fontsize=7,
            color=DIR["gained"], weight="semibold")
    ax.text(.99, cut / ax.get_ylim()[1] - .015, "FDR 0.05", transform=ax.transAxes, ha="right", va="top", fontsize=5.8,
            color="#888888")
    return len(treat), cut
