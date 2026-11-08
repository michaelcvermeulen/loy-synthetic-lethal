"""Main HT1080 RNA-ATAC promoter integration figure (thesis RNA branch).

Inputs (source_data/ht1080_rna_atac_integration) come from
prepare_ht1080_rna_atac_integration.R and extract_ht1080_atac_locus_coverage.py
(WSL). The integration itself is ATACseq/.../code/80_integrate_thesis_RNA.R.
Panel c reuses the supplementary promoter peak links (source_data/ht1080_rna_atac_peak_scatter,
RNA FDR 0.05 / ATAC FDR 0.10; prepare_ht1080_rna_atac_peak_scatter.R).
ATAC significance is peak FDR < 0.10 in integration panels; panel d loss groups use FDR < 0.05. RNA uses the thesis criteria.
Counts and correlations are descriptive group-level comparisons; nothing is refitted.
"""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_rna_atac_integration"
PEAKS = FIG / "source_data/ht1080_rna_atac_peak_scatter"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from thesis_style import apply_style
from atac_volcano import spread
from rna_atac_scatter import draw_link_points, link_legend_handles
apply_style(8)
STEM = "ch03_main_ht1080_rna_atac_integration"
COL = {"CTRL": "#3C5488", "LOY": "#E64B35"}
LABEL = {"CTRL": "sgCTRL", "LOY": "Y-depleted"}
LOWER, HIGHER = "#3C5488", "#E64B35"
GREY = "#555555"
LABEL_FS = 5.0  # scatter gene labels (reduced at author request, 22 Sep 2026)
HALO = [pe.withStroke(linewidth=1.5, foreground="white")]
SAMPLES = ["CTRL_11", "CTRL_12", "CTRL_21", "LOY_11", "LOY_12", "LOY_21"]
CAT = {"tested_promoter": ("Tested promoter peak", "#444444"), "no_promoter_peak": ("No promoter peak", "#9E9E9E"),
       "no_knownGene_TSS": ("No knownGene TSS", "#C9C9C9"), "no_ID_or_name_match": ("No ID or name match", "#E6E6E6")}
S = 230 / 250  # vertical scale for the lower section after the figure grew from 230 to 250 mm
UNNAMED = r"^(LOC|ENSG)"
ATAC_FDR = 0.10
HIGHLIGHT = "ONECUT2"  # author request, 22 Sep 2026: grey label; RNA up, promoter peak opening but FDR 0.145
N_EXTREME = 9  # extra grey labels for extreme non-joint links in panel c
GREY_OFFSET = {"MMP1": (0, 30), "ACSS1": (-4, -10)}  # per-gene label offsets in points, (dx, dy)
LEFT_SLOTS = 20  # labels in the lower-left ordered column

g = pd.read_csv(DATA / "integration_genes.csv")
att = pd.read_csv(DATA / "DE_gene_attrition.csv")
cors = pd.read_csv(DATA / "continuous_agreement.csv").set_index("subset")
loci = pd.read_csv(DATA / "locus_definitions.csv")
lpeaks = pd.read_csv(DATA / "locus_peaks.csv")
lgenes = pd.read_csv(DATA / "locus_genes.csv")
cov = pd.read_csv(DATA / "locus_coverage_50bp.csv")
expr = pd.read_csv(DATA / "locus_RNA_logCPM.csv")
prom = pd.read_csv(PEAKS / "promoter_peak_links.csv")
atac_loci = pd.read_csv(FIG / "source_data/ht1080_atac_accessibility/locus_definitions.csv")[["peak_id", "label"]]
assert len(g) == 13300 and (g.joint_010 == "concordant").sum() == 61 and len(att) == 662
assert set(pd.read_csv(PEAKS / "source_manifest.csv")[["RNA_FDR", "ATAC_FDR"]].drop_duplicates().itertuples(index=False)) == {(0.05, 0.10)}

fig = plt.figure(figsize=(183 / 25.4, 250 / 25.4))


def letter(x, y, ch):
    fig.text(x, y, ch, fontsize=12, weight="semibold", va="bottom")


def label_column(ax, rows, x_text, ha, colour_of, gap, lo, hi, weight_of=None):
    """Ordered label column beside a cluster: labels keep point order, spread to a minimum gap."""
    rows = rows.sort_values("y", ascending=False)
    ys = spread(rows.y.to_numpy(), gap, lo=lo, hi=hi)
    for r, ty in zip(rows.itertuples(), ys):
        ax.annotate(r.text, (r.x, r.y), xytext=(x_text, ty), fontsize=LABEL_FS, style="italic", ha=ha, va="center",
                    color=colour_of(r), weight=weight_of(r) if weight_of else "normal", path_effects=HALO, zorder=8,
                    arrowprops=dict(arrowstyle="-", color="#999999", lw=.3, shrinkA=1, shrinkB=1))


# a  where the 662 DE genes go
a = fig.add_axes([.105, .895, .185, .06])
for i, (dd, name, colour) in enumerate([("down", "Lower", LOWER), ("up", "Higher", HIGHER)]):
    sub = att[att.direction == dd]
    left = 0
    for key, (_, c_) in CAT.items():
        n = int((sub.category == key).sum())
        a.barh(i, n, left=left, color=c_, height=.62, edgecolor="white", linewidth=.5)
        if n >= 25:
            a.text(left + n / 2, i, f"{n}", ha="center", va="center", fontsize=5.8,
                   color="white" if key in ("tested_promoter", "no_promoter_peak") else "#222222")
        left += n
    a.text(-8, i, f"{name}\nn = {len(sub)}", ha="right", va="center", fontsize=6.3, color=colour, linespacing=1.1)
a.set(ylim=(1.6, -.6), yticks=[], xlim=(0, 440))
a.set_xticks([0, 200, 400]); a.tick_params(axis="x", labelsize=6.3)
a.set_xlabel("DE genes (thesis criteria)", fontsize=6.6, labelpad=1)
a.spines["left"].set_visible(False)
a.legend(handles=[Patch(facecolor=c_, label=n) for n, c_ in CAT.values()], loc="upper center", bbox_to_anchor=(.35, -.95),
         ncol=2, fontsize=5.6, handlelength=.9, columnspacing=.7, handletextpad=.3, borderpad=0)
a.text(0, 1.2, f"{int((att.category == 'tested_promoter').sum())}/662 DE genes have a tested promoter",
       transform=a.transAxes, fontsize=5.8, color=GREY)
letter(.02, .965, "a")

# b  promoter change rate by RNA class (gene level; RNA thesis criteria, ATAC peak FDR < 0.10)
b = fig.add_axes([.105, .715, .185, .085])
groups = [("down", "Lower"), ("not_DE", "Not DE"), ("up", "Higher")]
rates = []
for i, (k, name) in enumerate(groups):
    sub = g[g.rna_class == k]
    for j, (atac_k, colour) in enumerate([("lost", LOWER), ("gained", HIGHER)]):
        n = int((sub.atac_class_010 == atac_k).sum()); pct = 100 * n / len(sub)
        rates.append(dict(rna_class=k, atac_class_010=atac_k, n=n, total=len(sub), pct=pct))
        x = i + (j - .5) * .36
        b.bar(x, pct, width=.34, color=colour)
        b.text(x, pct + .5, f"{n}", ha="center", va="bottom", fontsize=5.3)
b.set_xticks(range(3), [f"{n}\n({(g.rna_class == k).sum():,})" for k, n in groups], fontsize=5.9)
for t_, col_ in zip(b.get_xticklabels(), [LOWER, "#222222", HIGHER]):
    t_.set_color(col_)
b.set_ylabel("Genes with a significant\npromoter peak (%)", fontsize=6.3)
b.set(ylim=(0, 30)); b.tick_params(axis="y", labelsize=6.3)
b.legend(handles=[Patch(facecolor=LOWER, label="Lost peak"), Patch(facecolor=HIGHER, label="Gained peak")],
         loc="upper right", fontsize=5.8, handlelength=1, borderpad=0)
pd.DataFrame(rates).to_csv(DATA / "promoter_rate_by_RNA_class.csv", index=False)
letter(.02, .802, "b")

# c  peak level, promoter links (RNA thesis criteria, ATAC FDR < 0.10)
c = fig.add_axes([.40, .725, .42, .24])
draw_link_points(c, prom)
c.set(xlim=(-5.8, 5.8), ylim=(-3.2, 3.4), xlabel="RNA log₂FC (Y-depleted / sgCTRL)",
      ylabel="Promoter peak ATAC log₂FC")
assert prom.rna_logFC.abs().max() < 5.8 and prom.logFC.between(-3.2, 2.9).all()  # no clipping
x = prom[prom.category.isin(["concordant", "opposite"])].assign(mag=lambda t: (t.rna_logFC * t.logFC).abs())
x = x.sort_values("mag", ascending=False).drop_duplicates("entrez")
loc = prom.merge(atac_loci, on="peak_id"); loc = loc[loc.gene_name == loc.label]
f_loci = prom[prom.gene_name.isin(loci.label)].sort_values("FDR").drop_duplicates("gene_name")  # panel-e loci
up = pd.concat([x[x.rna_class == "up"], f_loci[(f_loci.rna_logFC > 0) & (f_loci.gene_name != HIGHLIGHT)]]).drop_duplicates("gene_name")
label_column(c, up.assign(x=up.rna_logFC, y=up.logFC, text=up.gene_name), 3.4, "left",
             lambda r: "#222222" if r.category == "opposite" else HIGHER, .26, -3.0, 2.8)
down_fixed = pd.concat([loc[loc.logFC < 0], f_loci[f_loci.rna_logFC < 0]]).drop_duplicates("gene_name")
# Grey labels (author request, 22 Sep 2026): ONECUT2 plus the most extreme named links that are not jointly
# significant but significant in one assay (peak only or RNA only), ranked by distance from the origin on axes
# scaled to the data range; one per gene.
labelled0 = set(up.gene_name) | set(down_fixed.gene_name) | set(loc.label)
ext = prom[prom.category.isin(["peak_only", "RNA_only"]) & ~prom.gene_name.str.match(UNNAMED)
           & ~prom.gene_name.isin(labelled0)].copy()
ext["r"] = np.hypot(ext.rna_logFC / prom.rna_logFC.abs().max(), ext.logFC / prom.logFC.abs().max())
ext = ext.sort_values("r", ascending=False).drop_duplicates("entrez").head(N_EXTREME)
grey = pd.concat([prom[prom.gene_name == HIGHLIGHT].head(1), ext]).drop_duplicates("gene_name")
# Grey links inside the lower-left cluster join the ordered left column; the rest sit beside their points.
in_cluster = (grey.rna_logFC < 0) & (grey.logFC < -.8)
grey_col, grey_pt = grey[in_cluster], grey[~in_cluster]
down_top = x[(x.rna_class == "down") & ~x.gene_name.str.match(UNNAMED) & ~x.gene_name.isin(down_fixed.gene_name)]
left_rows = pd.concat([down_fixed, grey_col.assign(grey=True),
                       down_top.head(LEFT_SLOTS - len(down_fixed) - len(grey_col))])
left_rows["grey"] = left_rows["grey"].eq(True)
label_column(c, left_rows.assign(x=left_rows.rna_logFC, y=left_rows.logFC, text=left_rows.gene_name), -4.0, "right",
             lambda r: "#666666" if r.grey else ("#222222" if r.category != "concordant" else LOWER), .23, -3.05, 2.0)
for r in loc[loc.logFC > 0].itertuples():  # DCC: peak only, just outside the RNA criteria
    c.annotate(r.label, (r.rna_logFC, r.logFC), xytext=(-.8, 1.75), fontsize=LABEL_FS, style="italic", ha="right",
               va="center", color="#222222", path_effects=HALO, zorder=8,
               arrowprops=dict(arrowstyle="-", color="#999999", lw=.3, shrinkA=1, shrinkB=1))
for r in grey_pt.itertuples():
    dx, dy = GREY_OFFSET.get(r.gene_name, (8 if r.rna_logFC > 0 else -8, 6 if r.logFC > 0 else -6))
    c.annotate(r.gene_name, (r.rna_logFC, r.logFC), xytext=(dx, dy), textcoords="offset points", fontsize=LABEL_FS,
               style="italic", ha="center" if dx == 0 else ("left" if dx > 0 else "right"), va="center", color="#666666", path_effects=HALO,
               zorder=8, arrowprops=dict(arrowstyle="-", color="#AAAAAA", lw=.3, shrinkA=1, shrinkB=1))
grey[["gene_name", "rna_logFC", "logFC", "FDR", "rna_FDR", "category"]].to_csv(DATA / "panel_c_grey_labels.csv", index=False)
both = prom[prom.category.isin(["concordant", "opposite"])]
cnt = lambda rc, pc: int(((both.rna_class == rc) & (both.peak_class == pc)).sum())
c.text(.02, .98, "\n".join(["Promoter peak–gene links (RNA FDR < 0.05, ATAC FDR < 0.10)",
                            f"concordant  down/lost {cnt('down', 'lost')}   up/gained {cnt('up', 'gained')}",
                            f"opposite  up/lost {cnt('up', 'lost')}   down/gained {cnt('down', 'gained')}",
                            f"Spearman ρ  all {spearmanr(prom.rna_logFC, prom.logFC).statistic:.2f} "
                            f"({len(prom):,})   both {spearmanr(both.rna_logFC, both.logFC).statistic:.2f} ({len(both)})"]),
       transform=c.transAxes, ha="left", va="top", fontsize=5.4, color=GREY, linespacing=1.3)
letter(.34, .965, "c")
c.legend(handles=link_legend_handles(), loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=5.8,
         handletextpad=.3, labelspacing=.7, borderpad=0)

import runpy
runpy.run_path(str(Path(__file__).with_name('make_ht1080_joint_promoter_genes.py')))
jg = g[g.joint_010.isin(['concordant', 'opposite', 'mixed_promoter'])].sort_values('rna_logFC').reset_index(drop=True)

# d  TSS-centred coverage in RNA-defined groups (source: S16 panel a).
TDATA = FIG / "source_data/ht1080_atac_accessibility"
tgenes = pd.read_csv(TDATA / "tss_groups_genes.csv")
ttests = pd.read_csv(TDATA / "tss_groups_tests.csv").set_index("group")
tz = np.load(TDATA / "tss_transcript_matrix.npz")
tnb = tz[SAMPLES[0]].shape[1]
txs = (np.arange(tnb) + .5) * 6 / tnb - 3
tcentre = slice(tnb // 2 - 2, tnb // 2 + 2)
tgroups = [("down", "Lower expression"), ("down_no_lost_peak", "Lower expression, no called loss"),
           ("up", "Higher expression"), ("promoter_lost", "Lost promoter peak (ATAC-selected)")]
taxes = []; tstats = []; tymax = 0
for j, (group, title) in enumerate(tgroups):
    ax = fig.add_axes([.10, .52-j*.15, .20, .095]); taxes.append(ax)
    ratios = {}
    for grp, ls in [(group, "-"), (group+"_ref", (0,(3,2)))]:
        ids = tgenes.loc[tgenes.group==grp, "row"].to_numpy()
        means = {}
        for cond, keys in [("CTRL", SAMPLES[:3]), ("LOY", SAMPLES[3:])]:
            mat = np.array([tz[key][ids].mean(0) for key in keys])
            if grp==group: ax.fill_between(txs,mat.min(0),mat.max(0),color=COL[cond],alpha=.2,lw=0)
            ax.plot(txs,mat.mean(0),color=COL[cond],lw=1 if grp==group else .8,ls=ls,alpha=1 if grp==group else .4)
            means[cond] = float(mat.mean(0)[tcentre].mean()); tymax=max(tymax,float(mat.max()))
        ratios[grp] = means['LOY']/means['CTRL']
    ax.set_title(title+f"\nn = {int(ttests.loc[group,'targets'])} genes, {int(ttests.loc[group,'controls'])} controls",fontsize=5.8,weight='semibold',pad=3)
    ax.text(.97,.96,f"Genes {ratios[group]:.2f}x\nControls {ratios[group+'_ref']:.2f}x",transform=ax.transAxes,ha='right',va='top',fontsize=5.4,color=GREY)
    ax.set_xticks([-3,0,3],['-3','TSS','+3']); ax.tick_params(labelsize=5.6)
    tstats.append(dict(group=group,genes_ratio=ratios[group],controls_ratio=ratios[group+'_ref']))
for ax in taxes: ax.set_ylim(0,tymax*1.05)
fig.text(.027,.335,'Mean fragments per million non-Y fragments',rotation=90,ha='center',va='center',fontsize=6.1)
fig.text(.20,.043,'Distance from TSS (kb)',ha='center',fontsize=6.1)
fig.legend(handles=[Line2D([],[],color=COL['CTRL'],label='sgCTRL'),Line2D([],[],color=COL['LOY'],label='Y-depleted'),
                    Line2D([],[],color=GREY,label='Genes'),Line2D([],[],color=GREY,ls='--',alpha=.4,label='Matched controls')],
           loc='lower left',bbox_to_anchor=(.07,.653),ncol=2,fontsize=5.2,frameon=False,handlelength=1.3,columnspacing=.7)
letter(.02,.685,'d')
pd.DataFrame(tstats).to_csv(DATA/'panel_d_TSS_profiles.csv',index=False)

NL_ = chr(10)
# e  representative loci (3 rows x 2 columns): promoter coverage + per-sample RNA
E_TOP, E_ROW, E_X = .64, .205, (.40, .71)   # compact layout (author request: less white space)
tw, rw = .19, .05
PEAK_BAND = {"lost": "#C9D0E3", "gained": "#F3CFC8", "ns": "#E6E6E6"}   # at ATAC peak FDR < 0.10
for k_, r in enumerate(loci.itertuples()):
    x0 = E_X[k_ % 2]
    top = E_TOP - (k_ // 2) * E_ROW               # top of this cell (header baseline area)
    sub = cov[cov.peak_id == r.peak_id]
    ymax = sub.mean_fragment_cpm.max() * 1.1
    colour = HIGHER if r.rna_logFC > 0 else LOWER
    fig.text(x0, top, r.label, fontsize=7, weight="semibold", style="italic", color=colour, va="top")
    fig.text(x0 + .012 + .0075 * len(r.label), top - .001, r._3, fontsize=5.6, color=GREY, va="top")
    fig.text(x0, top - .013, f"RNA log₂FC {r.rna_logFC:+.2f} (FDR {r.rna_FDR:.1e})"
             + NL_ + f"promoter ATAC log₂FC {r.atac_median_logFC:+.2f} (peak FDR {r.atac_min_peak_FDR:.2g})",
             fontsize=5.3, color=GREY, va="top", linespacing=1.2)
    t_h, t_gap, hd = .040, .004, .040
    for k, cond in enumerate(["CTRL", "LOY"]):
        ax = fig.add_axes([x0, top - hd - (k + 1) * t_h - k * t_gap, tw, t_h])
        for p in lpeaks[lpeaks.locus == r.peak_id].itertuples():
            cls = "ns" if p.FDR >= ATAC_FDR else ("gained" if p.logFC > 0 else "lost")
            ax.axvspan(p.start, p.end, color=PEAK_BAND[cls], lw=0, zorder=0)
        mat = []
        for s in [s for s in SAMPLES if s.startswith(cond)]:
            t = sub[sub["sample"] == s]; xs = (t.bin_start + t.bin_end) / 2
            ax.plot(xs, t.mean_fragment_cpm, color=COL[cond], lw=.45, alpha=.7)
            mat.append(t.mean_fragment_cpm.to_numpy())
        ax.fill_between(xs, np.mean(mat, axis=0), color=COL[cond], alpha=.35, lw=0)
        ax.set(xlim=(r.window_start, r.window_end), ylim=(0, ymax), xticks=[])
        ax.set_yticks([0, round(ymax * .8, 1)]); ax.tick_params(axis="y", labelsize=5.3, length=2, pad=1)
        ax.spines["bottom"].set_visible(False)
        ax.text(.01, .95, LABEL[cond], transform=ax.transAxes, fontsize=5.4, color=COL[cond], va="top")
    gy = top - hd - 2 * t_h - t_gap - .046
    gx = fig.add_axes([x0, gy, tw, .044]); gx.axis("off")
    gx.set(xlim=(r.window_start, r.window_end), ylim=(-2.9, .9))
    gsub = lgenes[lgenes.locus == r.peak_id].reset_index(drop=True)
    gsub = pd.concat([gsub[gsub.symbol == r.label], gsub[gsub.symbol != r.label]]).reset_index(drop=True).head(2)
    for i, gr in gsub.iterrows():
        yy = -i * 1.05
        s0, s1 = max(gr.start, r.window_start), min(gr.end, r.window_end)
        gx.plot([s0, s1], [yy, yy], color="#444444", lw=1.0, solid_capstyle="butt")
        if r.window_start <= gr.tss <= r.window_end:
            sign = 1 if gr.strand == "+" else -1
            gx.plot([gr.tss, gr.tss], [yy, yy + .45], color="#444444", lw=.6)
            gx.annotate("", xy=(gr.tss + sign * 700, yy + .45), xytext=(gr.tss, yy + .45),
                        arrowprops=dict(arrowstyle="-|>", color="#444444", lw=.6, mutation_scale=4.5, shrinkA=0, shrinkB=0))
        lx = min(max((s0 + s1) / 2, r.window_start + 1500), r.window_end - 1500)
        gx.text(lx, yy - .18, str(gr.symbol), ha="center", va="top", fontsize=5.2, style="italic", color="#444444")
    gx.plot([r.window_end - 2000, r.window_end], [-2.65, -2.65], color="#222222", lw=.7)
    gx.text(r.window_end - 1000, -2.6, "2 kb", ha="center", va="bottom", fontsize=5.1)
    gx.text(r.window_start, -2.6, f"{r.chrom}:{r.tss / 1e6:.2f} Mb", ha="left", va="bottom", fontsize=5.1, color=GREY)
    rx = fig.add_axes([x0 + tw + .03, top - hd - 2 * t_h - t_gap, rw, 2 * t_h + t_gap])
    ex = expr[expr.gene_name == r.label]
    for i, cond in enumerate(["CTRL", "LOY"]):
        for period, marker, dx in [("AUG_2025", "o", -.12), ("MAY_2026", "^", .12)]:
            v = ex[(ex.condition == cond) & (ex.period == period)].logCPM
            rx.scatter(np.full(len(v), i + dx), v, s=8, marker=marker, color=COL[cond], lw=0, zorder=3)
        rx.plot([i - .3, i + .3], [ex[ex.condition == cond].logCPM.mean()] * 2, color="#222222", lw=.8, zorder=4)
    rx.set_xticks([0, 1], ["C", "L"], fontsize=5.6)
    for t_, cond in zip(rx.get_xticklabels(), ["CTRL", "LOY"]):
        t_.set_color(COL[cond])
    rx.set_xlim(-.6, 1.6); rx.tick_params(axis="y", labelsize=5.3, length=2, pad=1)
    rx.set_ylabel("RNA log₂CPM", fontsize=5.4, labelpad=1)
letter(.345, E_TOP + .012, "e")
# Legend for panel e only, drawn with the marks used in the panel (tracks and per-library RNA points).
e_handles = [Patch(facecolor=COL["CTRL"], alpha=.35, edgecolor=COL["CTRL"], label="sgCTRL"),
             Patch(facecolor=COL["LOY"], alpha=.35, edgecolor=COL["LOY"], label="Y-depleted"),
             Line2D([], [], ls="", marker="o", markersize=3.5, color="#666666", label="Aug 2025"),
             Line2D([], [], ls="", marker="^", markersize=3.5, color="#666666", label="May 2026"),
             Patch(facecolor=PEAK_BAND["lost"], label="Lost peak"), Patch(facecolor=PEAK_BAND["gained"], label="Gained peak"),
             Patch(facecolor=PEAK_BAND["ns"], label="Tested peak, n.s.")]
fig.legend(handles=e_handles, loc="lower left", bbox_to_anchor=(.40, E_TOP + .012), ncol=7, fontsize=5.6,
           handlelength=1.1, handletextpad=.3, columnspacing=.9, borderpad=0)
(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=180)
plt.close(fig)
json.dump({"universe": len(g), "joint": int(len(jg)), "joint_chr19": int((jg.chrom == "chr19").sum()),
                      "panel_c_labels": up.gene_name.tolist() + left_rows.gene_name.tolist() + loc[loc.logFC > 0].label.tolist(),
           "panel_c_grey_labels": grey.gene_name.tolist(),
           "joint_010": jg.joint_010.value_counts().to_dict(), "atac_fdr": 0.10,
           "DE_attrition": att.groupby(["direction", "category"]).size().rename("n").reset_index().to_dict("records"),
           "loci": loci[["label", "class"]].to_dict("records"), "font": "Open Sans", "size_mm": [183, 250]},
          open(DATA / "figure_checks.json", "w"), indent=1)
print("written", STEM)
