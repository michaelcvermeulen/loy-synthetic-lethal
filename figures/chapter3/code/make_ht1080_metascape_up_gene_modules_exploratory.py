"""Exploratory view of recurrent genes across UP Metascape pathway communities."""
from collections import Counter, defaultdict
from pathlib import Path
import hashlib
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
SOURCE = FIG / "source_data" / "ht1080_metascape_current" / "up_nodes.csv"
OUTDATA = FIG / "source_data" / "ht1080_metascape_up_gene_modules_exploratory"
OUTDATA.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(THESIS / "skills" / "thesis-figure-style" / "scripts"))
from thesis_style import apply_style


apply_style(8)
PALETTE = ["#3C5488", "#E64B35", "#00A087", "#F39B7F", "#4DBBD5", "#8491B4", "#7E6148"]
GREY = "#555555"


nodes = pd.read_csv(SOURCE)
module_table = (
    nodes[["community", "theme"]]
    .drop_duplicates()
    .sort_values("community")
    .reset_index(drop=True)
)
assert len(module_table) == 7
module_ids = module_table["community"].astype(int).tolist()
module_names = dict(zip(module_table["community"].astype(int), module_table["theme"]))
module_colours = dict(zip(module_ids, PALETTE))

gene_module_counts = defaultdict(Counter)
for row in nodes.itertuples():
    for gene in str(row.hits).split("|"):
        if gene and gene != "nan":
            gene_module_counts[gene][int(row.community)] += 1

records = []
for gene, counts in gene_module_counts.items():
    total = sum(counts.values())
    maximum = max(counts.values())
    # Community number is used only as a deterministic tie-break for display.
    dominant = min(module for module, value in counts.items() if value == maximum)
    record = {
        "gene": gene,
        "total_term_hits": total,
        "dominant_module": dominant,
        "dominant_theme": module_names[dominant],
        "dominant_term_hits": maximum,
        "dominant_fraction": maximum / total,
        "dominant_tie": sum(value == maximum for value in counts.values()) > 1,
    }
    for module in module_ids:
        record[f"module_{module}_term_hits"] = counts.get(module, 0)
    records.append(record)

all_genes = pd.DataFrame(records).sort_values(
    ["total_term_hits", "gene"], ascending=[False, True]
).reset_index(drop=True)
top = all_genes.head(40).copy()
top = top.sort_values(
    ["dominant_module", "total_term_hits", "gene"], ascending=[True, False, True]
).reset_index(drop=True)

all_genes.to_csv(OUTDATA / "all_gene_module_counts.csv", index=False)
top.to_csv(OUTDATA / "displayed_top40_gene_module_counts.csv", index=False)
module_table.assign(colour=module_table.community.map(module_colours)).to_csv(
    OUTDATA / "module_key.csv", index=False
)

# Insert visible gaps between dominant-module groups.
y_positions = []
group_bounds = []
cursor = 0.0
for module, group in top.groupby("dominant_module", sort=True):
    start = cursor
    for _ in range(len(group)):
        y_positions.append(cursor)
        cursor += 1.0
    end = cursor - 1.0
    group_bounds.append((module, start, end))
    cursor += 0.85
top["y"] = y_positions

fig = plt.figure(figsize=(183 / 25.4, 230 / 25.4))
ax = fig.add_axes([0.28, 0.225, 0.66, 0.635])
left = pd.Series(0.0, index=top.index)
for module in module_ids:
    values = top[f"module_{module}_term_hits"]
    ax.barh(
        top["y"], values, left=left, height=0.70,
        color=module_colours[module], edgecolor="white", linewidth=0.25,
        label=module_names[module]
    )
    left = left + values

ax.set_yticks(top["y"], top["gene"], fontsize=6.5)
ax.tick_params(axis="y", length=0, pad=4)
ax.invert_yaxis()
ax.set_xlabel("Number of higher-expression Metascape terms containing the gene")
ax.grid(axis="x", color="#E5E5E5", lw=0.5)
ax.set_axisbelow(True)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.set_xlim(0, top["total_term_hits"].max() + 6.5)

for row in top.itertuples():
    tie = "†" if row.dominant_tie else ""
    ax.text(
        row.total_term_hits + 0.35, row.y, f"{row.total_term_hits}{tie}",
        va="center", ha="left", fontsize=5.9, color=GREY
    )

for module, start, end in group_bounds:
    ax.text(
        0.15, start - 0.62, module_names[module],
        ha="left", va="center", fontsize=5.9,
        color=module_colours[module], weight="semibold"
    )
    if end < top["y"].max():
        ax.axhline(end + 0.92, color="#D8D8D8", lw=0.55)

fig.text(
    0.08, 0.965,
    "Genes recurring across higher-expression Metascape pathways",
    fontsize=11, weight="semibold", va="top"
)
fig.text(
    0.08, 0.925,
    "Top 40 of 105 represented genes across 124 enriched terms; rows grouped by dominant pathway community",
    fontsize=7.1, color=GREY, va="top"
)
handles, labels = ax.get_legend_handles_labels()
fig.legend(
    handles, labels, loc="lower left", bbox_to_anchor=(0.08, 0.055),
    ncol=3, frameon=False, fontsize=5.8, title="Pathway-community contribution",
    title_fontsize=6.3, handlelength=1.4, columnspacing=1.2
)
fig.text(
    0.08, 0.027,
    "Bar length counts enriched terms containing each gene; it does not measure expression magnitude.",
    fontsize=5.8, color=GREY, va="bottom"
)
fig.text(
    0.08, 0.011,
    "Group headings identify the dominant module. † Equal maxima in multiple modules; the lower-numbered module determines row placement.",
    fontsize=5.8, color=GREY, va="bottom"
)

stem = "ch03_exploratory_ht1080_metascape_up_gene_modules"
exports = []
for extension in ("pdf", "svg", "png"):
    target = FIG / f"{stem}.{extension}"
    staging = FIG / f"{stem}.staging.{extension}"
    fig.savefig(staging, dpi=600)
    staging.replace(target)
    exports.append({"file": target.name, "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
preview = FIG / "previews" / f"{stem}_preview.png"
preview.parent.mkdir(exist_ok=True)
fig.savefig(preview, dpi=160)
plt.close(fig)

provenance = {
    "status": "exploratory; not selected for the thesis figure set",
    "source": str(SOURCE),
    "terms": int(len(nodes)),
    "represented_genes": int(len(all_genes)),
    "displayed_genes": int(len(top)),
    "module_definition": "Louvain pathway communities from the current UP Metascape term-similarity network",
    "gene_assignment": "Module with the largest number of containing terms; ties resolved by lower community number for row placement only",
    "ranking": "Descending total number of enriched terms containing the gene; gene symbol resolves ties",
    "exports": exports,
}
(OUTDATA / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
print(json.dumps(exports, indent=2))
