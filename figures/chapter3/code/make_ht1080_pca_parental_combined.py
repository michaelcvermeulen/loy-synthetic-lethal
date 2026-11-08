"""Merge current HT1080 S01/S02 display panels without recomputing analyses."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

FIG = Path(__file__).resolve().parents[1]
ROOT = FIG.parents[1]
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)

STEM = "ch03_supp_ht1080_expression_structure_parental"
PCA = FIG / "source_data/ht1080_rna_parental_batch_pca"
PAR = FIG / "source_data/ht1080_rna_programme_overview"
OUT = FIG / "source_data/ht1080_pca_parental_combined"
OUT.mkdir(parents=True, exist_ok=True)
INPUTS = ([PCA / f"PCA_{v}.csv" for v in ("unadjusted", "adjusted", "batch1", "batch2")]
          + [PCA / "PCA_variance.csv", PCA / "PCA_within_batch_variance.csv"]
          + [PAR / f for f in ("DE_CTRL_vs_parental.csv", "heatmap_b_parental.csv",
                              "heatmap_layout.json", "parental_GSEA_provenance.json")])
hashes = {}
for source in INPUTS:
    target = OUT / source.name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if target.exists():
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, source
    else:
        shutil.copy2(source, target)
    hashes[str(source.relative_to(ROOT))] = digest

COL = {"CTRL": "#3C5488", "LOY": "#E64B35", "PARENTAL": "#888888"}
SHAPE = {"AUG_2025": "o", "MAY_2026": "^"}
variance = pd.concat([pd.read_csv(OUT / "PCA_variance.csv"),
                      pd.read_csv(OUT / "PCA_within_batch_variance.csv")], ignore_index=True)
fig = plt.figure(figsize=(11.8, 8.0))
positions = [[.07, .665, .175, .26], [.31, .665, .175, .26],
             [.55, .665, .175, .26], [.79, .665, .175, .26]]

def letter_at(letter, x, y):
    fig.text(x, y, letter, weight="semibold", fontsize=12)

for pos, view, letter, title in zip(
        positions, ["unadjusted", "adjusted", "batch1", "batch2"], "abcd",
        ["All libraries: unadjusted", "All libraries: period-adjusted",
         "August 2025", "May 2026"]):
    ax = fig.add_axes(pos)
    data = pd.read_csv(OUT / f"PCA_{view}.csv")
    assert len(data) == {"unadjusted": 15, "adjusted": 15, "batch1": 9, "batch2": 6}[view]
    v = variance[variance.view == view].set_index("component").variance_fraction
    for row in data.itertuples():
        ax.scatter(row.PC1, row.PC2, c=COL[row.treatment], marker=SHAPE[row.batch],
                   s=28, edgecolors="white", linewidths=.4, zorder=3)
        dx, dy, ha = 4, 4, "left"
        if view == "unadjusted" and row.batch == "MAY_2026" and row.treatment == "CTRL" and row.passage_number == 2:
            dx, dy, ha = -5, 5, "right"
        if view == "unadjusted" and row.batch == "AUG_2025" and row.treatment == "CTRL" and row.passage_number == 2:
            dx, dy, ha = -6, 3, "right"
        if view == "unadjusted" and row.batch == "AUG_2025" and row.treatment == "LOY" and row.passage_number == 1:
            dx, dy = 4, -9
        if view == "adjusted" and row.batch == "AUG_2025" and row.treatment == "CTRL" and row.passage_number in [2, 3]:
            dx, dy, ha = -6, -6, "right"
        if view == "adjusted" and row.batch == "MAY_2026" and row.treatment == "LOY" and row.passage_number == 3:
            dx, dy, ha = -5, -6, "right"
        if view == "adjusted" and row.treatment == "PARENTAL" and row.passage_number == 2:
            dx, dy, ha = -7, 8, "right"
        if view == "batch1" and row.treatment == "PARENTAL" and row.passage_number == 2:
            dx, dy = 5, -9
        if view == "batch2" and row.treatment == "LOY" and row.passage_number == 2:
            dx, dy, ha = -5, -6, "right"
        ax.annotate(str(row.passage_number), (row.PC1, row.PC2), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, color=COL[row.treatment], fontsize=7)
    ax.set_xlabel(f"PC1 ({100*v['PC1']:.1f}%)")
    ax.set_ylabel(f"PC2 ({100*v['PC2']:.1f}%)")
    if view == "unadjusted":
        ax.set(xlim=(-95, 132), ylim=(-108, 90))
        ax.set_xticks([-75, -25, 25, 75, 125])
        ax.set_yticks([-100, -50, 0, 50])
    elif view == "adjusted":
        ax.set(xlim=(-110, 92), ylim=(-88, 90))
        ax.set_xticks([-100, -50, 0, 50])
        ax.set_yticks([-50, 0, 50])
    else:
        ax.margins(x=.18, y=.20)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    letter_at(letter, pos[0] - .043, pos[1] + pos[3] + .014)
    fig.text(pos[0], pos[1] + pos[3] + .014, title, fontsize=8.5)
    ax.text(.02, .025, f"{len(data)} libraries", transform=ax.transAxes,
            fontsize=7, color="#555555")

handles = [Patch(facecolor=COL[t], label=l) for t, l in
           [("CTRL", "sgCTRL"), ("LOY", "Y-depleted"), ("PARENTAL", "Parental")]]
handles += [Line2D([], [], color="#333333", marker=m, lw=0, markersize=5, label=l)
            for m, l in [("o", "Aug 2025"), ("^", "May 2026")]]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.52, .987),
           ncol=5, fontsize=8, handlelength=1, columnspacing=1.7)

# Former S02 panel a: saved sgCTRL-versus-parental gene-level results.
ax = fig.add_axes([.065, .165, .265, .36])
parent = pd.read_csv(OUT / "DE_CTRL_vs_parental.csv")
sig = (parent.FDR < .05) & (parent.logFC.abs() >= 1)
lab = parent[sig].sort_values("logFC")
assert sig.sum() == 12
assert set(lab.gene_name) == {"HHIP", "PXDN", "CSMD2", "DOCK3", "LAT2", "FBLN1",
                              "IL2RB", "CECR2", "MARCKS", "ABLIM1", "MLPH", "ADARB2"}
ax.scatter(parent.logFC, -np.log10(parent.FDR), s=4, c="#C8C8C8", alpha=.4,
           lw=0, rasterized=True)
ax.scatter(lab.logFC, -np.log10(lab.FDR), s=20,
           c=np.where(lab.logFC > 0, COL["LOY"], COL["CTRL"]), lw=0)
for x in [-1, 1]:
    ax.axvline(x, c="#888888", ls="--", lw=.6)
ax.axhline(-np.log10(.05), c="#888888", ls="--", lw=.6)
ax.set(xlim=(-5.2, 6), ylim=(-.05, 3),
       xlabel="log₂ fold change (sgCTRL / parental)", ylabel="−log₁₀(BH-adjusted P)")
# Fixed label positions keep the twelve selected names legible in all exports.
label_positions = {
    "HHIP": (-1.8, 1.67), "PXDN": (.2, 1.62), "DOCK3": (.6, 1.92),
    "CSMD2": (1.1, 2.18), "LAT2": (1.5, 1.18), "FBLN1": (2.55, 2.44),
    "IL2RB": (2.75, 1.70), "CECR2": (3.5, 1.18), "MARCKS": (3.55, 2.20),
    "ABLIM1": (4.7, 2.49), "ADARB2": (4.85, 1.56), "MLPH": (5.25, 2.03),
}
for row in lab.itertuples():
    ax.annotate(row.gene_name, (row.logFC, -np.log10(row.FDR)),
                xytext=label_positions[row.gene_name], ha="center", va="center", fontsize=7,
                arrowprops=dict(arrowstyle="-", color="#888888", lw=.5,
                                shrinkA=2, shrinkB=3))
letter_at("e", .018, .55)
fig.text(.065, .55, "sgCTRL versus parental: expression", fontsize=8.5)

# Former S02 panel b: retain all 17 selected sets, order, NES and saved FDR.
b = pd.read_csv(OUT / "heatmap_b_parental.csv")
assert len(b) == 17 and (b.FDR < .05).all()
vmin, vmax = json.loads((OUT / "heatmap_layout.json").read_text())["shared_NES_limits"]
cmap = LinearSegmentedColormap.from_list("NES", [COL["CTRL"], "#F7F7F7", COL["LOY"]])
norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
ax = fig.add_axes([.587, .085, .047, .44])
im = ax.imshow(b[["NES"]].to_numpy(), cmap=cmap, norm=norm, aspect="auto")

def pathway_label(row):
    x = row.pathway
    for pre in ["HALLMARK_", "GOBP_", "REACTOME_"]:
        x = x.removeprefix(pre)
    return "\n".join(textwrap.wrap(f"[{row.group}] " + x.replace("_", " ").lower(), 32))

ax.set_yticks(range(len(b)), [pathway_label(q) for q in b.itertuples()], fontsize=7)
ax.tick_params(length=0, pad=4)
ax.spines[:].set_visible(False)
for i, row in enumerate(b.itertuples()):
    if i and b.group.iloc[i] != b.group.iloc[i-1]:
        ax.axhline(i-.5, color="white", lw=2)
    ax.text(0, i, f"{row.NES:.2f}", ha="center", va="center", fontsize=7,
            color="white" if abs(row.NES) > vmax*.58 else "#222222")
    ax.text(1.16, i, f"{row.FDR:.2g}", transform=ax.get_yaxis_transform(), va="center", fontsize=7)
ax.set_xticks([])
ax.text(.5, 1.015, "NES", transform=ax.transAxes, fontsize=7, ha="center")
ax.text(1.16, 1.015, "FDR", transform=ax.transAxes, fontsize=7)
letter_at("f", .345, .55)
fig.text(.38, .55, "sgCTRL versus parental: GSEA", fontsize=8.5)
cb = fig.colorbar(im, cax=fig.add_axes([.425, .042, .209, .009]), orientation="horizontal")
cb.set_label("Normalised enrichment score (NES)", fontsize=7, labelpad=2)
cb.ax.tick_params(labelsize=7, length=2)
cb.outline.set_visible(False)

# Former F01b: retain the exact saved matrix and dendrogram order.
from ht1080_expression_heatmaps import draw_saved_correlation
draw_saved_correlation(fig, FIG)

for ext in ("pdf", "svg", "png"):
    stage = FIG / f"{STEM}.staging.{ext}"
    fig.savefig(stage, dpi=600)
    stage.replace(FIG / f"{STEM}.{ext}")
fig.savefig(FIG / "previews" / f"{STEM}_preview.png", dpi=170)
plt.close(fig)

for source in INPUTS:
    assert hashlib.sha256(source.read_bytes()).hexdigest() == hashes[str(source.relative_to(ROOT))]
provenance = {
    "date": "2026-10-09", "working_id": "CH3-S3.3.6.1-S01", "analysis_recomputed": False,
    "panels": {"a": "all-15 unadjusted PCA", "b": "all-15 period-adjusted PCA",
               "c": "August PCA (9 libraries)", "d": "May PCA (6 libraries)",
               "e": "former S02a: sgCTRL versus parental volcano (12 qualifying genes)",
               "f": "former S02b: 17 selected GSEA sets",
               "g": "former F01b: saved expression correlations and clustering"},
    "size_inches": [11.8, 8.0], "NES_limits": [vmin, vmax], "input_sha256": hashes,
    "PCA_coordinates_and_variance_unchanged": True,
    "DE_and_GSEA_values_and_FDR_families_unchanged": True,
    "qualifying_gene_labels": lab.gene_name.tolist(),
    "GSEA_pathway_order": b.pathway.tolist(),
}
(OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"figure": STEM, "PCA_libraries": [15, 15, 9, 6],
                  "gene_labels": len(lab), "GSEA_rows": len(b), "analysis_recomputed": False}))
