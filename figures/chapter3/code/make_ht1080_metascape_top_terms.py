"""Plot the top Metascape terms for higher- and lower-expression HT1080 genes."""
from pathlib import Path
import hashlib
import json
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
INPUT = FIG / "source_data" / "ht1080_metascape_current"
OUTPUT = FIG / "source_data" / "ht1080_metascape_top_terms"
OUTPUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(THESIS / "skills" / "thesis-figure-style" / "scripts"))
from thesis_style import apply_style


apply_style(8)
BLUE = "#3C5488"
VERMILION = "#E64B35"
GREY = "#555555"


def wrap_label(value, width=48):
    value = str(value).replace(" HCP WITH H3K27ME3", " HCP H3K27me3")
    return "\n".join(textwrap.wrap(value, width=width))


def q_text(log10q):
    qvalue = 10 ** float(log10q)
    return f"{qvalue:.1e}"


def load(direction):
    data = pd.read_csv(INPUT / f"{direction}_nodes.csv")
    if "n_hits" not in data:
        data["n_hits"] = data["hits"].fillna("").str.split("|").map(
            lambda values: sum(bool(value) for value in values)
        )
    data = data.sort_values(
        ["mlog10p", "log10q", "description"], ascending=[False, True, True]
    ).head(20).copy()
    data["direction"] = direction
    data["rank"] = np.arange(1, len(data) + 1)
    data["q_value"] = 10 ** data["log10q"]
    return data


up = load("up")
down = load("down")
up.to_csv(OUTPUT / "higher_top20.csv", index=False)
down.to_csv(OUTPUT / "lower_top20.csv", index=False)


fig = plt.figure(figsize=(183 / 25.4, 235 / 25.4))
ax_up = fig.add_axes([0.40, 0.555, 0.50, 0.33])
ax_down = fig.add_axes([0.40, 0.115, 0.50, 0.33])


def draw(axis, frame, colour, title, total_terms):
    shown = frame.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(shown))
    axis.barh(y, shown["mlog10p"], color=colour, alpha=0.88, height=0.68)
    axis.set_yticks(y, [wrap_label(value) for value in shown["description"]], fontsize=6.1)
    axis.tick_params(axis="y", length=0, pad=4)
    axis.set_xlabel(r"Metascape enrichment $-\log_{10}P$")
    axis.grid(axis="x", color="#E5E5E5", lw=0.5)
    axis.set_axisbelow(True)
    axis.spines[["top", "right", "left"]].set_visible(False)
    maximum = max(frame["mlog10p"].max(), 1)
    axis.set_xlim(0, maximum + 5.3)
    for index, row in shown.iterrows():
        axis.text(
            row["mlog10p"] + 0.20,
            index,
            f"{int(row['n_hits'])} genes; q={q_text(row['log10q'])}",
            fontsize=5.7,
            va="center",
            ha="left",
            color=GREY,
        )
    axis.text(
        0,
        1.045,
        title,
        transform=axis.transAxes,
        color=colour,
        fontsize=8.5,
        weight="semibold",
        ha="left",
    )
    axis.text(
        1,
        1.045,
        f"Top 20 of {total_terms} exported terms",
        transform=axis.transAxes,
        color=GREY,
        fontsize=6.3,
        ha="right",
    )


draw(ax_up, up, VERMILION, "Higher expression in Y-depleted HT1080", 124)
draw(ax_down, down, BLUE, "Lower expression in Y-depleted HT1080", 76)

fig.text(
    0.08,
    0.965,
    "Top Metascape terms associated with Y depletion",
    fontsize=11,
    weight="semibold",
    va="top",
)
fig.text(
    0.08,
    0.925,
    "Over-representation analysis of 229 higher- and 433 lower-expression genes",
    fontsize=7.2,
    color=GREY,
    va="top",
)
fig.text(0.04, 0.89, "a", fontsize=12, weight="semibold")
fig.text(0.04, 0.45, "b", fontsize=12, weight="semibold")
fig.text(
    0.08,
    0.047,
    "Terms are ranked independently in each direction by nominal enrichment P. "
    "Numbers beside bars give input-list hits and Metascape-reported q values.",
    fontsize=6.2,
    color=GREY,
    va="bottom",
)
fig.text(
    0.08,
    0.025,
    "Related terms may share genes and are not statistically independent.",
    fontsize=6.2,
    color=GREY,
    va="bottom",
)

stem = "ch03_supp_ht1080_metascape_top_terms"
exports = []
for extension in ("pdf", "svg", "png"):
    target = FIG / f"{stem}.{extension}"
    staging = FIG / f"{stem}.staging.{extension}"
    fig.savefig(staging, dpi=600)
    staging.replace(target)
    exports.append(
        {"file": target.name, "sha256": hashlib.sha256(target.read_bytes()).hexdigest()}
    )
preview = FIG / "previews" / f"{stem}_preview.png"
preview.parent.mkdir(exist_ok=True)
fig.savefig(preview, dpi=160)
plt.close(fig)

provenance = {
    "source": str(INPUT),
    "selection": "Top 20 exported terms in each direction by descending -log10 nominal P; deterministic ties by log10 q and description.",
    "gene_lists": {"higher": 229, "lower": 433, "background": 25461},
    "exported_terms": {"higher": 124, "lower": 76, "lower_connected_in_network": 70},
    "font": "Open Sans",
    "exports": exports,
}
(OUTPUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
print(json.dumps(exports, indent=2))
