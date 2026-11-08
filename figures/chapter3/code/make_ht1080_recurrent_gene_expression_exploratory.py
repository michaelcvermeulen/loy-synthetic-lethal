"""Exploratory sample-level expression figure for recurrent UP Metascape genes."""
from pathlib import Path
import hashlib
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
SOURCE = FIG / "source_data" / "ht1080_recurrent_gene_expression"
sys.path.insert(0, str(THESIS / "skills" / "thesis-figure-style" / "scripts"))
from thesis_style import apply_style


apply_style(8)
BLUE = "#3C5488"
RED = "#E64B35"
PARENTAL = "#777777"
GREY = "#A8A8A8"
COLOURS = {"Parental": PARENTAL, "sgCTRL": BLUE, "Y-depleted": RED}
MARKERS = {"August": "o", "May": "^"}

expression = pd.read_csv(SOURCE / "sample_expression.csv")
statistics = pd.read_csv(SOURCE / "batch_aware_statistics.csv").set_index("gene_name")

columns = [
    ("Neurite morphogenesis", ["WNT7A", "PAX2", "CHRNB2"]),
    ("Epithelial/cardiovascular development", ["GATA4", "BMP7", "IHH"]),
    ("Synapse organisation/signalling", ["NLGN4X", "MAPK8IP2", "NRXN2"]),
]
conditions = ["Parental", "sgCTRL", "Y-depleted"]
xmap = {condition: index for index, condition in enumerate(conditions)}
passage_offsets = {1: -0.09, 2: 0.0, 3: 0.09}
batch_offsets = {"August": -0.025, "May": 0.025}

fig, axes = plt.subplots(3, 3, figsize=(183 / 25.4, 165 / 25.4), sharex=True)
fig.subplots_adjust(left=0.10, right=0.98, bottom=0.15, top=0.76, wspace=0.26, hspace=0.43)

displayed = []
for column_index, (module, genes) in enumerate(columns):
    for row_index, gene in enumerate(genes):
        ax = axes[row_index, column_index]
        data = expression[expression["gene"] == gene].copy()

        # Connect the six confirmed within-period, passage-matched sgCTRL/LOY pairs.
        for (batch, passage), pair in data[data["condition"].isin(["sgCTRL", "Y-depleted"])].groupby(
            ["batch_label", "passage_number"], sort=True
        ):
            if set(pair["condition"]) == {"sgCTRL", "Y-depleted"}:
                pair = pair.set_index("condition").loc[["sgCTRL", "Y-depleted"]]
                xs = [
                    xmap[condition] + passage_offsets[int(passage)] + batch_offsets[batch]
                    for condition in ["sgCTRL", "Y-depleted"]
                ]
                ax.plot(xs, pair["log2_CPM"], color="#C8C8C8", lw=0.55, zorder=1)

        for record in data.itertuples():
            x = (
                xmap[record.condition]
                + passage_offsets[int(record.passage_number)]
                + batch_offsets[record.batch_label]
            )
            ax.scatter(
                x, record.log2_CPM, s=25, marker=MARKERS[record.batch_label],
                facecolor=COLOURS[record.condition], edgecolor="white", linewidth=0.55,
                zorder=3
            )
            displayed.append(
                {
                    "gene": gene, "sample_key": record.sample_key,
                    "condition": record.condition, "batch": record.batch_label,
                    "passage_number": record.passage_number,
                    "plot_x": x, "log2_CPM": record.log2_CPM,
                }
            )

        stat = statistics.loc[gene]
        ax.text(0.02, 1.02, gene, transform=ax.transAxes, va="bottom", ha="left",
                fontsize=8, weight="semibold", clip_on=False)
        ax.text(
            0.98, 1.02,
            f"log₂FC {stat.LOY_vs_sgCTRL_log2FC:+.2f}\nFDR {stat.LOY_vs_sgCTRL_FDR:.2g}",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=5.6,
            color="#444444", linespacing=0.92, clip_on=False
        )
        ax.set_xlim(-0.35, 2.35)
        ax.grid(axis="y", color="#E7E7E7", lw=0.5)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", length=0)
        if column_index == 0:
            ax.set_ylabel("log₂ CPM")
        if row_index == 2:
            ax.set_xticks(range(3), conditions, rotation=22, ha="right")
        else:
            ax.tick_params(labelbottom=False)

for column_index, (module, _) in enumerate(columns):
    position = axes[0, column_index].get_position()
    fig.text((position.x0 + position.x1) / 2, 0.805, module,
             ha="center", va="bottom", fontsize=7.2, weight="semibold")

fig.text(0.07, 0.965, "Expression of recurrent genes across HT1080 populations",
         fontsize=11, weight="semibold", va="top")
fig.text(
    0.07, 0.925,
    "All-15-library contextual TMM normalisation; lines connect passage-matched sgCTRL and Y-depleted samples",
    fontsize=7.1, color="#555555", va="top"
)

condition_handles = [
    Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=COLOURS[name],
           markeredgecolor="white", markersize=6, label=name)
    for name in conditions
]
batch_handles = [
    Line2D([0], [0], marker=marker, linestyle="none", markerfacecolor="#555555",
           markeredgecolor="white", markersize=6, label=batch)
    for batch, marker in MARKERS.items()
]
fig.legend(condition_handles + batch_handles, conditions + list(MARKERS),
           loc="upper center", bbox_to_anchor=(0.64, 0.885), frameon=False,
           ncol=5, fontsize=6.5, handletextpad=0.35, columnspacing=0.9)
fig.text(
    0.07, 0.040,
    "Points are individual RNA-seq libraries; shapes identify sampling period and small horizontal offsets identify passage.",
    fontsize=6.0, color="#555555", va="bottom"
)
fig.text(
    0.07, 0.020,
    "Parental samples were collected only in August and are contextual because matching to the edited populations is not established.",
    fontsize=6.0, color="#555555", va="bottom"
)

stem = "ch03_exploratory_ht1080_recurrent_gene_expression"
exports = []
for extension in ("pdf", "svg", "png"):
    target = FIG / f"{stem}.{extension}"
    staging = FIG / f"{stem}.staging.{extension}"
    fig.savefig(staging, dpi=600)
    staging.replace(target)
    exports.append({"file": target.name, "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
preview = FIG / "previews" / f"{stem}_preview.png"
preview.parent.mkdir(exist_ok=True)
fig.savefig(preview, dpi=170)
plt.close(fig)

pd.DataFrame(displayed).to_csv(SOURCE / "displayed_points.csv", index=False)
(SOURCE / "figure_provenance.json").write_text(json.dumps({
    "status": "exploratory; not selected for the thesis figure set",
    "expression": "TMM log2 CPM with prior count 2 from an all-15-library contextual normalization",
    "statistics": "Batch + condition limma-voom model on the 12 core sgCTRL and Y-depleted libraries",
    "parental": "August only; contextual and not passage-matched to edited populations",
    "exports": exports,
}, indent=2) + "\n", encoding="utf-8")
print(json.dumps(exports, indent=2))
