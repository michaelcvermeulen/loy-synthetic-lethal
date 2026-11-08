"""Supplementary figure: over-representation in RNA-ATAC agreeing quadrants and polycomb-target promoters (HT1080).

Sources (ATACseq/SEPT_08_2026_REPROCESS/integration/thesis_RNA_SEPT11; scripts 80-83 in that run's code/):
  RNA_ATAC_promoter_genes.csv, quadrant_ORA_rna025_atac025/ and quadrant_ORA_rna050_atac050/,
  polycomb_matched/polycomb_matched_tests.csv and _K1.csv.
Inputs are copied into source_data/ht1080_quadrant_ORA with SHA256 hashes; nothing is refitted here.
"""
from pathlib import Path
import sys, json, shutil, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_quadrant_ORA"
SRC = ROOT.parent / "ATACseq/SEPT_08_2026_REPROCESS/integration/thesis_RNA_SEPT11"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM = "ch03_supp_ht1080_quadrant_ORA"
LOWER, HIGHER, GREY = "#3C5488", "#E64B35", "#555555"
QCOL = {"UP_OPEN": HIGHER, "DOWN_CLOSED": LOWER, "UP_CLOSED": "#8C6BB1", "DOWN_OPEN": "#E3A33B", "none": "#D2D2D2"}
QLAB = {"UP_OPEN": "UP-OPEN", "DOWN_CLOSED": "DOWN-CLOSED", "UP_CLOSED": "UP-CLOSED", "DOWN_OPEN": "DOWN-OPEN"}
PC_SETS = {"BENPORATH_ES_WITH_H3K27ME3": "ES H3K27me3", "BENPORATH_PRC2_TARGETS": "PRC2 targets"}

# ---- freeze inputs ----
DATA.mkdir(parents=True, exist_ok=True)
files = {"RNA_ATAC_promoter_genes.csv": SRC / "RNA_ATAC_promoter_genes.csv",
         "ORA_all_rna025_atac025.csv": SRC / "quadrant_ORA_rna025_atac025/ORA_all.csv",
         "quadrant_counts_rna050_atac050.json": SRC / "quadrant_ORA_rna050_atac050/summary.json",
         "PRC2_membership_by_quadrant.csv": SRC / "quadrant_ORA_rna025_atac025/PRC2_membership_by_quadrant.csv",
         "polycomb_matched_tests_K3.csv": SRC / "polycomb_matched/polycomb_matched_tests.csv",
         "polycomb_matched_tests_K1.csv": SRC / "polycomb_matched/polycomb_matched_tests_K1.csv"}
manifest = []
for name, src in files.items():
    shutil.copyfile(src, DATA / name)
    manifest.append(dict(file=name, source=str(src), sha256=hashlib.sha256(src.read_bytes()).hexdigest()))
pd.DataFrame(manifest).to_csv(DATA / "source_manifest.csv", index=False)

g = pd.read_csv(DATA / "RNA_ATAC_promoter_genes.csv")
ora = pd.read_csv(DATA / "ORA_all_rna025_atac025.csv")
q050 = json.load(open(DATA / "quadrant_counts_rna050_atac050.json"))["counts"]
memb = pd.read_csv(DATA / "PRC2_membership_by_quadrant.csv")
k3 = pd.read_csv(DATA / "polycomb_matched_tests_K3.csv"); k1 = pd.read_csv(DATA / "polycomb_matched_tests_K1.csv")
T1, T2 = .25, .5
g["quadrant"] = np.select([(g.rna_logFC > T1) & (g.atac_median_logFC > T1), (g.rna_logFC < -T1) & (g.atac_median_logFC < -T1),
                           (g.rna_logFC > T1) & (g.atac_median_logFC < -T1), (g.rna_logFC < -T1) & (g.atac_median_logFC > T1)],
                          ["UP_OPEN", "DOWN_CLOSED", "UP_CLOSED", "DOWN_OPEN"], "none")
n025 = g.quadrant.value_counts()
assert (n025["UP_OPEN"], n025["DOWN_CLOSED"]) == (391, 573)

fig = plt.figure(figsize=(183 / 25.4, 215 / 25.4))


def letter(x, y, ch):
    fig.text(x, y, ch, fontsize=12, weight="semibold", va="bottom")


# a  quadrants
a = fig.add_axes([.075, .665, .30, .285])
for q in ["none", "UP_CLOSED", "DOWN_OPEN", "UP_OPEN", "DOWN_CLOSED"]:
    s = g[g.quadrant == q]
    a.scatter(s.rna_logFC, s.atac_median_logFC, s=1.2 if q == "none" else 3, c=QCOL[q], lw=0,
              alpha=1 if q == "none" else .7, rasterized=True, zorder=1 if q == "none" else 2)
for v, ls in [(T1, (0, (3, 2))), (T2, (0, (1, 1.5)))]:
    for sgn in (-1, 1):
        a.axvline(sgn * v, color="#888888", lw=.5, ls=ls, zorder=0); a.axhline(sgn * v, color="#888888", lw=.5, ls=ls, zorder=0)
a.set(xlim=(-5.8, 5.8), ylim=(-3.2, 3.2), xlabel="RNA log₂FC (Y-depleted / sgCTRL)",
      ylabel="Promoter ATAC log₂FC (median of tested peaks)")
assert g.rna_logFC.abs().max() < 5.8 and g.atac_median_logFC.abs().max() < 3.2
for q, (x, y, ha) in {"UP_OPEN": (5.5, 2.9, "right"), "DOWN_CLOSED": (-5.5, -2.9, "left"),
                      "UP_CLOSED": (5.5, -2.9, "right"), "DOWN_OPEN": (-5.5, 2.9, "left")}.items():
    a.text(x, y, f"{QLAB[q]}\n{n025[q]:,} ({q050[q]:,})", ha=ha, va="top" if y > 0 else "bottom", fontsize=6.4,
           color=QCOL[q], weight="semibold", linespacing=1.2)
a.text(.5, 1.02, "Genes per quadrant at |log₂FC| > 0.25 (> 0.5); 13,300 genes", transform=a.transAxes,
       ha="center", fontsize=6.2, color=GREY)
letter(.01, .962, "a")

# b  ORA dot plots (universe background, thresholds 0.25/0.25)
cond = ora[ora.background != "universe"].set_index(["quadrant", "collection", "gs_name"]).BH_FDR
uni = ora[ora.background == "universe"]
NORM = Normalize(0, 12); CMAP = plt.get_cmap("viridis")
def pretty(name):
    for p in ["HALLMARK_", "GOBP_", "REACTOME_"]:
        name = name.replace(p, "")
    s = name.replace("_", " ").lower()
    return "\n".join(__import__("textwrap").wrap(s, 43))
rows_out = []
for k, (q, box) in enumerate([("DOWN_CLOSED", [.60, .70, .12, .25]), ("UP_OPEN", [.60, .415, .12, .185])]):
    s = uni[(uni.quadrant == q) & (uni.BH_FDR < .05)].sort_values("BH_FDR")
    s = s.groupby("collection", sort=False).head(4)
    order = {"Hallmark": 0, "GO:BP": 1, "CGP": 2, "Reactome": 3}
    s = s.assign(o=s.collection.map(order)).sort_values(["o", "BH_FDR"]).reset_index(drop=True)
    ax = fig.add_axes(box)
    for i, r in s.iterrows():
        cf = cond.get((q, r.collection, r.gs_name), np.nan)
        ax.scatter(r.fold_enrichment, i, s=6 + 1.6 * r.Count, color=CMAP(NORM(min(-np.log10(r.BH_FDR), 12))),
                   edgecolors="#111111" if cf < .05 else "none", linewidths=.9, zorder=3)
        rows_out.append(dict(quadrant=q, collection=r.collection, gs_name=r.gs_name, Count=r.Count, query_n=r.query_n,
                             fold_enrichment=r.fold_enrichment, BH_FDR=r.BH_FDR, conditional_BH_FDR=cf))
    ax.set_yticks(range(len(s)), [f"{pretty(n)}" for n in s.gs_name], fontsize=5.2)
    for t_, col in zip(ax.get_yticklabels(), s.collection):
        t_.set_color({"Hallmark": "#222222", "GO:BP": "#444444", "CGP": "#6A3D9A", "Reactome": "#444444"}[col])
    for col in s.collection.unique():
        idx = s.index[s.collection == col]
        ax.text(1.04, (idx.min() + idx.max()) / 2, col, transform=ax.get_yaxis_transform(), fontsize=5.6, color=GREY,
                va="center", rotation=90)
        if idx.min() > 0:
            ax.axhline(idx.min() - .5, color="#DDDDDD", lw=.5, zorder=0)
    ax.set(ylim=(len(s) - .5, -.5), xscale="log", xlim=(1.8, 20))
    ax.set_xticks([2, 5, 10, 20], ["2", "5", "10", "20"], fontsize=6); ax.minorticks_off()
    ax.set_xlabel("Fold enrichment vs all genes", fontsize=6.6)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.set_title(f"{QLAB[q]} ({int(s.query_n.iloc[0])} genes; top 4 sets per collection)", fontsize=7, weight="semibold",
                 color=QCOL[q], pad=4)
pd.DataFrame(rows_out).to_csv(DATA / "panel_b_sets.csv", index=False)
cax = fig.add_axes([.905, .83, .009, .10])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, extend="max")
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.6, length=2); cb.set_label("−log₁₀ BH FDR", fontsize=5.8)
bh = [Line2D([], [], ls="", marker="o", markersize=np.sqrt(6 + 1.6 * n), color="#999999", label=str(n)) for n in (10, 30, 60)]
bh.append(Line2D([], [], ls="", marker="o", markersize=5, markerfacecolor="white", markeredgecolor="#111111",
                 markeredgewidth=.9, label="FDR < 0.05 vs\nRNA-direction\nbackground"))
fig.legend(handles=bh, loc="upper left", bbox_to_anchor=(.885, .80), fontsize=5.6, title="Genes", title_fontsize=5.8,
           labelspacing=.9, handletextpad=.4, borderpad=0, frameon=False)
letter(.39, .962, "b")

# c  polycomb-target membership by RNA x ATAC direction
c = fig.add_axes([.10, .42, .17, .145])
m = memb[memb.gene_set == "BENPORATH_ES_WITH_H3K27ME3"].set_index(["rna", "atac"])
rows_, cols_ = ["up", "flat", "down"], ["closed", "flat", "open"]
mat = np.array([[m.loc[(r, cl), "pct"] for cl in cols_] for r in rows_])
im = c.imshow(mat, cmap="magma_r", vmin=0, vmax=16, aspect="auto")
for i, r in enumerate(rows_):
    for j, cl in enumerate(cols_):
        c.text(j, i, f"{mat[i, j]:.1f}%\n{int(m.loc[(r, cl), 'members'])}/{int(m.loc[(r, cl), 'n']):,}", ha="center",
               va="center", fontsize=5.4, color="white" if mat[i, j] > 8 else "#222222", linespacing=1.15)
c.set_xticks(range(3), ["closed", "flat", "open"], fontsize=6.2); c.set_yticks(range(3), ["up", "flat", "down"], fontsize=6.2)
c.set_xlabel("Promoter ATAC (|log₂FC| > 0.25)", fontsize=6.4); c.set_ylabel("RNA (|log₂FC| > 0.25)", fontsize=6.4)
c.tick_params(length=0)
for sp in c.spines.values(): sp.set_visible(False)
c.set_title("ES H3K27me3 targets (% of genes)", fontsize=6.8, weight="semibold", pad=4)
cb2 = fig.colorbar(im, cax=fig.add_axes([.285, .43, .008, .125])); cb2.outline.set_visible(False)
cb2.ax.tick_params(labelsize=5.4, length=2)
letter(.01, .60, "c")

# d  matched-background forest plot
ROWS = [("gene | ATAC closed (|log2FC| > 0.25) vs ALL flat, UNMATCHED | all genes", "Closed genes, unmatched"),
        ("gene | ATAC closed (|log2FC| > 0.25) vs matched flat | all genes", "Closed genes (> 0.25), matched"),
        ("gene | ATAC closed (|log2FC| > 0.25) vs matched flat | RNA-flat genes", "  RNA-flat genes only"),
        ("gene | ATAC closed (|log2FC| > 0.50) vs matched flat | all genes", "Closed genes (> 0.5), matched"),
        ("gene | ATAC open (|log2FC| > 0.25) vs matched flat | all genes", "Open genes (> 0.25), matched"),
        ("peak | lost promoter peaks FDR < 0.05 vs matched non-DA promoter peaks", "Lost peaks, FDR < 0.05"),
        ("peak | lost promoter peaks FDR < 0.10 vs matched non-DA promoter peaks", "Lost peaks, FDR < 0.10"),
        ("peak | gained promoter peaks FDR < 0.10 vs matched non-DA promoter peaks", "Gained peaks, FDR < 0.10")]
d = fig.add_axes([.30, .06, .55, .25])
fr = []
for i, (cmp, name) in enumerate(ROWS):
    for dy, (set_, col) in zip((-.16, .16), [("BENPORATH_ES_WITH_H3K27ME3", "#6A3D9A"), ("BENPORATH_PRC2_TARGETS", "#1B9E77")]):
        r = k3[(k3.comparison == cmp) & (k3.gene_set == set_)].iloc[0]
        r1 = k1[(k1.comparison == cmp) & (k1.gene_set == set_)]
        lo, hi = max(r.CI_low, .1), min(r.CI_high, 30)
        d.plot([lo, hi], [i + dy, i + dy], color=col, lw=.9, solid_capstyle="butt", zorder=2)
        d.scatter(r.MH_OR, i + dy, s=18, color=col if r.BH_FDR < .05 else "white", edgecolors=col, linewidths=.8, zorder=3)
        if len(r1) and "UNMATCHED" not in cmp and np.isfinite(r1.MH_OR.iloc[0]):
            d.scatter(min(r1.MH_OR.iloc[0], 30), i + dy, s=9, marker="|", color=col, zorder=3)
        fr.append(dict(row=name, gene_set=set_, n_targets=r.targets, OR_K3=r.MH_OR, CI_low=r.CI_low, CI_high=r.CI_high,
                       BH_FDR_K3=r.BH_FDR, OR_K1=r1.MH_OR.iloc[0] if len(r1) else np.nan))
    d.text(32, i, f"n = {int(r.targets):,}", fontsize=5.4, color=GREY, va="center")
d.axvline(1, color="#888888", lw=.6, ls=(0, (4, 3)), zorder=0)
d.axhline(.5, color="#DDDDDD", lw=.6); d.axhline(4.5, color="#DDDDDD", lw=.6)
d.set_xscale("log"); d.set(xlim=(.1, 30), ylim=(len(ROWS) - .4, -.6))
d.set_xticks([.1, .3, 1, 3, 10, 30], ["0.1", "0.3", "1", "3", "10", "30"], fontsize=6); d.minorticks_off()
d.set_yticks(range(len(ROWS)), [n for _, n in ROWS], fontsize=6)
d.tick_params(axis="y", length=0); d.spines["left"].set_visible(False)
d.set_xlabel("Odds ratio (Mantel–Haenszel, 95% CI)", fontsize=6.6)
pd.DataFrame(fr).to_csv(DATA / "panel_d_forest.csv", index=False)
dh = [Line2D([], [], ls="", marker="o", markersize=4, color="#6A3D9A", label="ES H3K27me3 targets"),
      Line2D([], [], ls="", marker="o", markersize=4, color="#1B9E77", label="PRC2 targets"),
      Line2D([], [], ls="", marker="o", markersize=4, markerfacecolor="white", markeredgecolor="#555555", label="BH FDR ≥ 0.05"),
      Line2D([], [], ls="", marker="|", markersize=5, color="#555555", label="K = 1 matching")]
d.legend(handles=dh, loc="lower center", bbox_to_anchor=(.45, 1.01), ncol=4, fontsize=5.6, handletextpad=.3,
         columnspacing=.9, borderpad=0)
d.text(0, -.2, "Matched on baseline logCPM × GC × CpG o/e (K = 3 controls; bars show K = 1)",
       transform=d.transAxes, fontsize=5.6, color=GREY)
letter(.01, .345, "d")

(FIG / "previews").mkdir(exist_ok=True)
for ext in ["pdf", "svg", "png"]:
    target, stage = FIG / f"{STEM}.{ext}", FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600, facecolor="white"); stage.replace(target)
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=180)
plt.close(fig)
print("written", STEM)
