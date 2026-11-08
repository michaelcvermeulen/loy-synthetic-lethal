"""Selected supplementary composite of Metascape terms and recurrent UP genes."""
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
TERM_DATA = FIG / "source_data" / "ht1080_metascape_top_terms"
GENE_DATA = FIG / "source_data" / "ht1080_metascape_up_gene_modules_exploratory"
OUTPUT = FIG / "source_data" / "ht1080_metascape_terms_and_gene_modules"
OUTPUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(THESIS / "skills" / "thesis-figure-style" / "scripts"))
from thesis_style import apply_style


apply_style(7.5)
BLUE = "#3C5488"
VERMILION = "#E64B35"
GREY = "#555555"
PALETTE = ["#3C5488", "#E64B35", "#00A087", "#F39B7F", "#4DBBD5", "#8491B4", "#7E6148"]


up = pd.read_csv(TERM_DATA / "higher_top20.csv")
down = pd.read_csv(TERM_DATA / "lower_top20.csv")
genes = pd.read_csv(GENE_DATA / "displayed_top40_gene_module_counts.csv")
modules = pd.read_csv(GENE_DATA / "module_key.csv").sort_values("community")
module_ids = modules["community"].astype(int).tolist()
module_names = dict(zip(module_ids, modules["theme"]))
module_colours = dict(zip(module_ids, PALETTE))


def wrap(value, width=39):
    value = str(value).replace(" HCP WITH H3K27ME3", " HCP H3K27me3")
    return "\n".join(textwrap.wrap(value, width=width))


def q_text(log10q):
    return f"{10 ** float(log10q):.1e}"


fig = plt.figure(figsize=(290 / 25.4, 210 / 25.4))
ax_up = fig.add_axes([0.175, 0.555, 0.285, 0.315])
ax_down = fig.add_axes([0.175, 0.125, 0.285, 0.315])
ax_gene = fig.add_axes([0.585, 0.15, 0.385, 0.72])


def draw_terms(axis, frame, colour, title):
    shown = frame.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(shown))
    axis.barh(y, shown["mlog10p"], color=colour, alpha=0.88, height=0.68)
    axis.set_yticks(y, [wrap(value) for value in shown["description"]], fontsize=5.5)
    axis.tick_params(axis="y", length=0, pad=3)
    axis.set_xlabel(r"Metascape enrichment $-\log_{10}P$", fontsize=6.5)
    axis.grid(axis="x", color="#E5E5E5", lw=0.45)
    axis.set_axisbelow(True)
    axis.spines[["top", "right", "left"]].set_visible(False)
    maximum = max(frame["mlog10p"].max(), 1)
    axis.set_xlim(0, maximum + 5.6)
    for index, row in shown.iterrows():
        axis.text(
            row["mlog10p"] + 0.16, index,
            f"{int(row['n_hits'])} genes; q={q_text(row['log10q'])}",
            fontsize=4.8, va="center", ha="left", color=GREY
        )
    axis.text(
        0, 1.035, title, transform=axis.transAxes,
        color=colour, fontsize=7.2, weight="semibold", ha="left"
    )


draw_terms(ax_up, up, VERMILION, "Higher expression in Y-depleted HT1080")
draw_terms(ax_down, down, BLUE, "Lower expression in Y-depleted HT1080")

# Recreate grouped row positions for the recurrent-gene view.
y_positions = []
group_bounds = []
cursor = 0.0
for module, group in genes.groupby("dominant_module", sort=True):
    start = cursor
    for _ in range(len(group)):
        y_positions.append(cursor)
        cursor += 1.0
    group_bounds.append((int(module), start, cursor - 1.0))
    cursor += 0.85
genes = genes.copy()
genes["y"] = y_positions

left = pd.Series(0.0, index=genes.index)
for module in module_ids:
    values = genes[f"module_{module}_term_hits"]
    ax_gene.barh(
        genes["y"], values, left=left, height=0.70,
        color=module_colours[module], edgecolor="white", linewidth=0.22,
        label=module_names[module]
    )
    left = left + values

ax_gene.set_yticks(genes["y"], genes["gene"], fontsize=5.7)
ax_gene.tick_params(axis="y", length=0, pad=3)
ax_gene.invert_yaxis()
ax_gene.set_xlabel("Number of higher-expression Metascape terms containing the gene", fontsize=6.5)
ax_gene.grid(axis="x", color="#E5E5E5", lw=0.45)
ax_gene.set_axisbelow(True)
ax_gene.spines[["top", "right", "left"]].set_visible(False)
ax_gene.set_xlim(0, genes["total_term_hits"].max() + 6.2)
for row in genes.itertuples():
    tie = "†" if row.dominant_tie else ""
    ax_gene.text(row.total_term_hits + 0.3, row.y, f"{row.total_term_hits}{tie}",
                 va="center", fontsize=5.1, color=GREY)
for module, start, end in group_bounds:
    ax_gene.text(
        0.15, start - 0.62, module_names[module],
        ha="left", va="center", fontsize=5.2,
        color=module_colours[module], weight="semibold"
    )
    if end < genes["y"].max():
        ax_gene.axhline(end + 0.92, color="#D8D8D8", lw=0.5)
ax_gene.text(
    0, 1.035, "Genes recurring across higher-expression pathways",
    transform=ax_gene.transAxes, fontsize=7.2, weight="semibold", ha="left"
)

fig.text(0.06, 0.955, "Metascape pathways and their recurrent contributing genes",
         fontsize=11, weight="semibold", va="top")
fig.text(
    0.06, 0.918,
    "Over-representation analysis of 229 higher- and 433 lower-expression genes; 25,461-gene tested background",
    fontsize=7, color=GREY, va="top"
)
fig.text(0.025, 0.875, "a", fontsize=11, weight="semibold")
fig.text(0.025, 0.445, "b", fontsize=11, weight="semibold")
fig.text(0.505, 0.875, "c", fontsize=11, weight="semibold")

handles, labels = ax_gene.get_legend_handles_labels()
fig.legend(
    handles, labels, loc="lower left", bbox_to_anchor=(0.575, 0.035),
    ncol=2, frameon=False, fontsize=5.2,
    title="UP pathway-community contribution", title_fontsize=5.8,
    handlelength=1.4, columnspacing=1.1
)
fig.text(
    0.06, 0.025,
    "a,b, Top 20 terms per direction ranked by nominal enrichment P; labels give input-list hits and Metascape q.  "
    "c, Top 40 genes ranked by recurrence across UP terms and grouped by dominant pathway community.",
    fontsize=5.6, color=GREY, va="bottom"
)
fig.text(
    0.06, 0.009,
    "Pathway recurrence reflects overlapping annotations rather than expression magnitude. † Equal maxima in multiple modules.",
    fontsize=5.6, color=GREY, va="bottom"
)

stem = "ch03_supp_ht1080_metascape_terms_and_gene_modules"
exports = []
for extension in ("pdf", "svg", "png"):
    target = FIG / f"{stem}.{extension}"
    staging = FIG / f"{stem}.staging.{extension}"
    fig.savefig(staging, dpi=600)
    staging.replace(target)
    exports.append({"file": target.name, "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
preview = FIG / "previews" / f"{stem}_preview.png"
fig.savefig(preview, dpi=150)
plt.close(fig)

(OUTPUT / "provenance.json").write_text(json.dumps({
    "status": "selected supplementary figure CH3-S3.4.6.1-S04",
    "term_panels": str(TERM_DATA),
    "gene_panel": str(GENE_DATA),
    "layout": "290 x 210 mm landscape composite",
    "exports": exports,
}, indent=2) + "\n", encoding="utf-8")
print(json.dumps(exports, indent=2))
