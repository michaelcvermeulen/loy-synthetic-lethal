"""Supplementary peak-level RNA-ATAC scatter by peak type, RNA FDR 0.05 / ATAC FDR 0.10 (HT1080).

Inputs from prepare_ht1080_rna_atac_peak_scatter.R. One point per peak-gene link:
a = promoter peaks linked to every overlapping TSS (+/- 2 kb); b = gene-body peaks and
c = intergenic peaks, each linked to the nearest expressed gene's TSS (exploratory).
Descriptive; nothing is refitted. The earlier two-panel version (promoter + all peaks)
is archived in archive/peak_scatter_two_panel_2026-09-22/.
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
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
# Optional args: RNA_FDR ATAC_FDR (default 0.05 0.10; author decision 22 Sep 2026). Other thresholds read the matching suffixed
# source_data folder and write an exploratory figure (not the supplementary figure).
RNA_FDR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
ATAC_FDR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.10
SUFFIX = "" if (RNA_FDR, ATAC_FDR) == (0.05, 0.10) else f"_rna{round(RNA_FDR * 100):03d}_atac{round(ATAC_FDR * 100):03d}"
DATA = FIG / f"source_data/ht1080_rna_atac_peak_scatter{SUFFIX}"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from thesis_style import apply_style
from atac_volcano import spread
from rna_atac_scatter import draw_link_points, link_legend_handles
apply_style(8)
HALO = [pe.withStroke(linewidth=1.6, foreground="white")]
STEM = "ch03_supp_ht1080_rna_atac_peak_scatter" if not SUFFIX else f"ch03_exploratory_ht1080_rna_atac_peak_scatter{SUFFIX}"
LOWER, HIGHER, GREY, PEAK_ONLY = "#3C5488", "#E64B35", "#555555", "#9A9A9A"
TOP_DOWN = 4  # down/lost genes labelled per distal panel

prom = pd.read_csv(DATA / "promoter_peak_links.csv")
near = pd.read_csv(DATA / "nearest_gene_peak_links.csv")
panels = [("promoter", "Promoter (TSS ± 2 kb)", prom),
          ("gene_body", "Gene body (nearest gene)", near[near.annotation == "gene_body"]),
          ("intergenic", "Intergenic (nearest gene)", near[near.annotation == "intergenic"])]
ATAC_LOCI = pd.read_csv(FIG / "source_data/ht1080_atac_accessibility/locus_definitions.csv")[["peak_id", "label"]]
LOCUS_TEXT = {"DCC": (-.4, 2.6), "CNTN4": (-2.1, -2.85), "INSR": (.35, -2.45), "SMPD3": (.35, -2.05)}
xlim = max(abs(d.rna_logFC).max() for _, _, d in panels) + .4
ylim = max(abs(d.logFC).max() for _, _, d in panels) + .3

fig = plt.figure(figsize=(183 / 25.4, 92 / 25.4))
axes = [fig.add_axes([.065 + i * .318, .15, .265, .64]) for i in range(3)]
checks, cor_rows, label_rows = {}, [], []
for ax, (key, title, d), ch in zip(axes, panels, "abc"):
    draw_link_points(ax, d)
    ax.set(xlim=(-xlim, xlim), ylim=(-ylim, ylim + 1.2), xlabel="RNA log₂FC (Y-depleted / sgCTRL)")
    ax.set_xticks([-6, -3, 0, 3, 6])
    if ch == "a":
        ax.set_ylabel("Peak ATAC log₂FC")
    ax.set_title(title, fontsize=7, weight="semibold", pad=4)

    # Labels (fixed rule, one per gene at its strongest link): every up-regulated gene with a significant
    # peak, every down-regulated gene with a gained peak, and - in the distal panels - the TOP_DOWN
    # down/lost genes by |RNA log2FC x peak log2FC|.
    x = d[d.category.isin(["concordant", "opposite"])].copy()
    x["mag"] = (x.rna_logFC * x.logFC).abs()
    x = x.sort_values("mag", ascending=False).drop_duplicates("entrez")
    groups = [(x[(x.rna_class == "up") & (x.logFC > 0)], 3.0, "left", HIGHER),
              (x[(x.rna_class == "up") & (x.logFC < 0)], 3.0, "left", "#222222"),
              (x[(x.rna_class == "down") & (x.logFC > 0)], -3.0, "right", "#222222")]
    if key != "promoter":
        named = ~x.gene_name.str.match(r"^(LOC|ENSG)")  # unnamed loci are not picked for the top-down labels
        groups.append((x[(x.rna_class == "down") & (x.logFC < 0) & named].head(TOP_DOWN), -3.6, "right", LOWER))
    for sub, tx, ha, colour in groups:
        sub = sub.sort_values("logFC", ascending=False)
        ys = spread(sub.logFC.to_numpy(), .3, lo=-ylim + .15, hi=ylim + .6)
        for r, ty in zip(sub.itertuples(), ys):
            ax.annotate(r.gene_name, (r.rna_logFC, r.logFC), xytext=(tx, ty), fontsize=5.5, style="italic",
                        ha=ha, va="center",
                        color=colour, path_effects=HALO, zorder=6,
                        arrowprops=dict(arrowstyle="-", color="#999999", lw=.3, shrinkA=1, shrinkB=1.2))
            label_rows.append(dict(panel=key, gene_name=r.gene_name, category=r.category, rna_class=r.rna_class,
                                   rna_logFC=r.rna_logFC, peak_logFC=r.logFC, peak_id=r.peak_id))
    if key == "promoter":  # ATAC main-figure loci, labelled wherever they fall (DDIT4-AS1 is not expressed)
        ref = d.merge(ATAC_LOCI, on="peak_id")
        ref = ref[ref.gene_name == ref.label]
        for r in ref.itertuples():
            tx, ty = LOCUS_TEXT[r.label]
            ax.annotate(r.label, (r.rna_logFC, r.logFC), xytext=(tx, ty), fontsize=5.5, style="italic", color="#222222",
                        ha="left" if tx > r.rna_logFC else "right", va="center", path_effects=HALO, zorder=8,
                        arrowprops=dict(arrowstyle="-", color="#555555", lw=.35, shrinkA=1, shrinkB=1))
        checks["atac_figure_loci"] = ref[["label", "rna_logFC", "logFC", "rna_FDR", "category"]].to_dict("records")

    both = d[d.category.isin(["concordant", "opposite"])]
    rho_all = spearmanr(d.rna_logFC, d.logFC).statistic
    rho_both = spearmanr(both.rna_logFC, both.logFC).statistic if len(both) > 2 else np.nan
    cor_rows += [dict(panel=key, subset="all", n=len(d), spearman=rho_all),
                 dict(panel=key, subset="both significant", n=len(both), spearman=rho_both)]
    n = lambda rc, pc: (int(((both.rna_class == rc) & (both.peak_class == pc)).sum()),
                        int(both[(both.rna_class == rc) & (both.peak_class == pc)].entrez.nunique()))
    dl, ug, ul, dg = n("down", "lost"), n("up", "gained"), n("up", "lost"), n("down", "gained")
    txt = "\n".join([f"{len(d):,} links · both significant (links/genes):",
                     f"concordant  down/lost {dl[0]}/{dl[1]}   up/gained {ug[0]}/{ug[1]}",
                     f"opposite  up/lost {ul[0]}/{ul[1]}   down/gained {dg[0]}/{dg[1]}",
                     f"Spearman ρ  all {rho_all:.2f}   both {rho_both:.2f} (n = {len(both)})"])
    ax.text(.03, .97, txt, transform=ax.transAxes, ha="left", va="top", fontsize=5.3, color=GREY, linespacing=1.35)
    fig.text(ax.get_position().x0 - .05, .855, ch, fontsize=12, weight="semibold")
    checks[key] = {"links": len(d), "peaks": int(d.peak_id.nunique()), "genes": int(d.entrez.nunique()),
                   "categories": d.category.value_counts().to_dict()}
    assert d.rna_logFC.abs().max() < xlim and d.logFC.abs().max() < ylim  # no clipping
for ax in axes[1:]:
    ax.tick_params(axis="y", labelleft=True)
pd.DataFrame(cor_rows).to_csv(DATA / "spearman_by_peak_type.csv", index=False)
pd.DataFrame(label_rows).to_csv(DATA / "labels_by_peak_type.csv", index=False)

h = link_legend_handles()
fig.legend(handles=h, loc="upper center", bbox_to_anchor=(.5, .995), ncol=7, fontsize=6.2, handletextpad=.25, columnspacing=.9)
fig.text(.5, .92, f"RNA: FDR < {RNA_FDR:.2f}, |log₂FC| ≥ 1, same direction in both periods   ·   ATAC peak: FDR < {ATAC_FDR:.2f}",
         ha="center", fontsize=6.6, color=GREY)

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=180)
plt.close(fig)
json.dump(checks, open(DATA / "figure_checks.json", "w"), indent=1, default=float)
print("written", STEM)
