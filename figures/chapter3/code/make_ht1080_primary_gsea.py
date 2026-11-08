"""Create the selected primary HT1080 GSEA figure for Chapter 3."""
from pathlib import Path
import hashlib
import json
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd


FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
TABLE = (
    THESIS.parent
    / "RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE"
    / "outputs/programme_analysis_v1/tables/programme_consistency.csv"
)
sys.path.insert(0, str(THESIS / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style


apply_style(8)
BLUE = "#3C5488"
RED = "#E64B35"
LIGHT_BLUE = "#AEBBD3"
LIGHT_RED = "#F2B2A7"
TEXT = "#242424"
CMAP = LinearSegmentedColormap.from_list("signed_nes", [BLUE, "#F7F7F7", RED])


def clean_name(value):
    for prefix in ("HALLMARK_", "REACTOME_", "GOBP_"):
        value = value.removeprefix(prefix)
    value = value.replace("_", " ").lower()
    return "\n".join(textwrap.wrap(value, width=43))


def export(fig, selected, counts):
    stem = "ch03_main_ht1080_primary_gsea"
    source = FIG / "source_data" / "ht1080_primary_gsea"
    source.mkdir(parents=True, exist_ok=True)
    selected.to_csv(source / "displayed_pathways.csv", index=False)
    counts.to_csv(source / "discovery_counts.csv", index=False)

    records = []
    for ext in ("pdf", "svg", "png"):
        output = FIG / f"{stem}.{ext}"
        staging = FIG / f"{stem}.staging.{ext}"
        fig.savefig(staging, dpi=600)
        staging.replace(output)
        records.append(
            {"file": output.name, "sha256": hashlib.sha256(output.read_bytes()).hexdigest()}
        )
    preview = FIG / "previews" / f"{stem}_preview.png"
    preview.parent.mkdir(exist_ok=True)
    fig.savefig(preview, dpi=170)
    plt.close(fig)

    provenance = {
        "source_table": str(TABLE),
        "font": "Open Sans",
        "testing_family": "Hallmark, Reactome and GO Biological Process",
        "eligibility": "overall FDR < 0.05 and matching NES directions in August and May",
        "display_selection": (
            "Four smallest-overall-FDR eligible sets per collection, plus the three "
            "smallest-overall-FDR positively enriched sets; duplicates removed."
        ),
        "exports": records,
    }
    (source / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(stem)


data = pd.read_csv(TABLE)
eligible = data[data["directionally_consistent"].astype(bool)].copy()
eligible["direction"] = np.where(eligible["NES"] < 0, "Lower", "Higher")

assert len(data) == 4868
assert len(eligible) == 220
assert (eligible["NES"] < 0).sum() == 209
assert (eligible["NES"] > 0).sum() == 11
assert eligible["both_periods_FDR05"].astype(bool).sum() == 55

selected = pd.concat(
    [
        group.sort_values(["FDR", "pathway"]).head(4)
        for _, group in eligible.groupby("analysis_collection", sort=True)
    ]
    + [eligible[eligible["NES"] > 0].sort_values(["FDR", "pathway"]).head(3)],
    ignore_index=True,
).drop_duplicates("pathway")
selected = selected.sort_values(["NES", "FDR", "pathway"]).reset_index(drop=True)

collection_order = ["Hallmark", "Reactome", "GO Biological Process"]
counts = []
for collection in collection_order:
    subset = eligible[eligible["analysis_collection"] == collection]
    for direction in ("Lower", "Higher"):
        group = subset[subset["direction"] == direction]
        counts.append(
            {
                "collection": collection,
                "direction": direction,
                "overall_and_direction_matched": len(group),
                "significant_in_both_periods": int(group["both_periods_FDR05"].astype(bool).sum()),
            }
        )
counts = pd.DataFrame(counts)

fig = plt.figure(figsize=(183 / 25.4, 205 / 25.4))
# Stack the count summary over the heatmap. This leaves enough width for full
# pathway names at the final 183-mm page width.
ax_count = fig.add_axes([0.18, 0.715, 0.68, 0.145])
ax_heat = fig.add_axes([0.49, 0.115, 0.31, 0.50])

# Panel a: discovery counts. Pale bars are the directionally persistent total;
# saturated overlays are the subset significant in both periods.
y = np.arange(len(collection_order))
for i, collection in enumerate(collection_order):
    row = counts[counts["collection"] == collection]
    lower = int(row.loc[row.direction == "Lower", "overall_and_direction_matched"].iloc[0])
    higher = int(row.loc[row.direction == "Higher", "overall_and_direction_matched"].iloc[0])
    lower_both = int(row.loc[row.direction == "Lower", "significant_in_both_periods"].iloc[0])
    higher_both = int(row.loc[row.direction == "Higher", "significant_in_both_periods"].iloc[0])
    ax_count.barh(i, -lower, height=0.54, color=LIGHT_BLUE, edgecolor="none")
    ax_count.barh(i, -lower_both, height=0.54, color=BLUE, edgecolor="none")
    ax_count.barh(i, higher, height=0.54, color=LIGHT_RED, edgecolor="none")
    ax_count.barh(i, higher_both, height=0.54, color=RED, edgecolor="none")
    if lower:
        ax_count.text(-lower - 4, i, str(lower), ha="right", va="center", fontsize=7, color=BLUE)
    if higher:
        ax_count.text(higher + 4, i, str(higher), ha="left", va="center", fontsize=7, color=RED)

ax_count.axvline(0, color="#777777", lw=0.65)
ax_count.set_yticks(y, collection_order)
ax_count.invert_yaxis()
ax_count.set_xlim(-175, 25)
ax_count.set_xticks([-150, -100, -50, 0], ["150", "100", "50", "0"])
ax_count.set_xlabel("Number of enriched gene sets")
ax_count.tick_params(axis="y", length=0)
ax_count.spines[["top", "right", "left"]].set_visible(False)
ax_count.grid(axis="x", color="#E6E6E6", lw=0.5, zorder=0)
ax_count.set_axisbelow(True)
ax_count.text(0.02, -0.28, "Lower in Y-depleted", transform=ax_count.transAxes,
              color=BLUE, ha="left", va="top", fontsize=7, weight="semibold")
ax_count.text(0.98, -0.28, "Higher", transform=ax_count.transAxes,
              color=RED, ha="right", va="top", fontsize=7, weight="semibold")

# Panel b: representative pathways across the overall and period-specific models.
matrix = selected[["NES", "AUG_2025_NES", "MAY_2026_NES"]].to_numpy()
pvalues = selected[["FDR", "AUG_2025_FDR", "MAY_2026_FDR"]].to_numpy()
limit = np.ceil(np.abs(matrix).max() * 10) / 10
image = ax_heat.imshow(
    matrix, aspect="auto", cmap=CMAP,
    norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
)
labels = []
codes = {"Hallmark": "H", "Reactome": "R", "GO Biological Process": "GO"}
for row in selected.itertuples():
    labels.append(f"[{codes[row.analysis_collection]}]  {clean_name(row.pathway)}")
ax_heat.set_yticks(np.arange(len(selected)), labels, fontsize=6.2)
ax_heat.set_xticks([0, 1, 2], ["Overall", "August", "May"])
ax_heat.tick_params(axis="both", length=0)
for i, j in zip(*np.where(pvalues < 0.05)):
    colour = "white" if abs(matrix[i, j]) > limit * 0.55 else TEXT
    ax_heat.text(j, i, "*", ha="center", va="center", fontsize=9, color=colour, weight="semibold")
for i, value in enumerate(selected["FDR"]):
    ax_heat.text(3.0, i, f"{value:.1e}", va="center", ha="left", fontsize=5.9, clip_on=False)
ax_heat.text(3.0, -1.05, "Overall\nFDR", va="bottom", ha="left", fontsize=6.1, clip_on=False)
ax_heat.spines[:].set_visible(False)

# Separate the lower- and higher-enrichment representatives.
n_lower = int((selected["NES"] < 0).sum())
ax_heat.axhline(n_lower - 0.5, color="#777777", lw=0.7)
ax_heat.text(-0.02, 1.025, "Representative pathways", transform=ax_heat.transAxes,
             fontsize=7.4, weight="semibold", ha="left")

cax = fig.add_axes([0.925, 0.235, 0.014, 0.27])
colourbar = fig.colorbar(image, cax=cax)
colourbar.ax.set_title("NES", fontsize=6.5, pad=6)
colourbar.ax.tick_params(labelsize=6, length=2)

# Figure hierarchy and concise on-canvas statistical summary.
fig.text(0.09, 0.955, "Persistent transcriptional programmes after Y depletion",
         fontsize=11, weight="semibold", ha="left", va="top")
fig.text(
    0.09, 0.902,
    "Hallmark, Reactome and GO Biological Process  |  4,868 tested  |  "
    "220 persistent: 209 lower, 11 higher  |  55 significant in both periods",
    fontsize=7.2, color="#4D4D4D", ha="left", va="top"
)
fig.text(0.04, 0.865, "a", fontsize=12, weight="semibold")
fig.text(0.04, 0.62, "b", fontsize=12, weight="semibold")
fig.text(
    0.09, 0.047,
    "Saturated bars: FDR < 0.05 in August and May; pale extent: overall FDR < 0.05 with matching directions.",
    fontsize=6.2, color="#555555", ha="left", va="bottom"
)
fig.text(
    0.09, 0.025,
    "* FDR < 0.05 in the indicated contrast.  H, Hallmark; R, Reactome; GO, Gene Ontology Biological Process.",
    fontsize=6.2, color="#555555", ha="left", va="bottom"
)

export(fig, selected, counts)
