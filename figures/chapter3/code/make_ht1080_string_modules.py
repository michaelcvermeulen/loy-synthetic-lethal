"""Plot current HT1080 STRING/MCODE modules in the thesis figure style."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
ANALYSIS = THESIS.parent / "RNAseq" / "analysis" / "HT1080" / "SEPT_11_2026_REANALYZE"
STRING = ANALYSIS / "outputs" / "programme_analysis_v1" / "tables" / "STRING"
SOURCE = FIG / "source_data" / "ht1080_string_modules"
SOURCE.mkdir(parents=True, exist_ok=True)

import sys
sys.path.insert(0, str(THESIS / "skills" / "thesis-figure-style" / "scripts"))
from thesis_style import apply_style

apply_style(8)

PALETTE = {
    "DOWN_01": "#3C5488",
    "DOWN_02": "#00A087",
    "DOWN_03": "#8491B4",
    "UP_01": "#E64B35",
    "UP_02": "#F39B7F",
}
TITLES = {
    "DOWN_01": "Interferon-response module",
    "DOWN_02": "IL12-family receptor module",
    "DOWN_03": "IGFBP3–MMP module",
    "UP_01": "BMP-response module",
    "UP_02": "Neurexin–neuroligin module",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fmt_effect(value: float) -> str:
    return f"{value:+.2f}".replace("-", "−")


def load_data():
    modules = pd.read_csv(STRING / "modules.csv")
    effects = pd.read_csv(STRING / "all_symbol_mapping_and_effects.csv")
    effects = effects.loc[effects["directionally_consistent"], ["gene_name", "logFC", "FDR", "direction"]]
    association = pd.read_csv(STRING / "all_candidate_association_edges.csv")
    physical = pd.read_csv(STRING / "all_candidate_physical_edges.csv")
    physical_pairs = {tuple(sorted((row.gene_a, row.gene_b))) for row in physical.itertuples()}
    return modules, effects, association, physical_pairs


def circular_positions(members: list[str]) -> dict[str, np.ndarray]:
    n = len(members)
    angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, n, endpoint=False)
    if n == 4:
        angles += np.pi / 4
    return {gene: np.array([np.cos(angle), np.sin(angle)]) for gene, angle in zip(members, angles)}


def draw_module(ax, row, effects, association, physical_pairs):
    members = str(row.members).split(";")
    effect = effects.set_index("gene_name").loc[members]
    graph = nx.Graph()
    graph.add_nodes_from(members)
    subset = association[association["gene_a"].isin(members) & association["gene_b"].isin(members)]
    graph.add_weighted_edges_from(subset[["gene_a", "gene_b", "score"]].itertuples(index=False, name=None))
    pos = circular_positions(members)

    for u, v, data in graph.edges(data=True):
        pair = tuple(sorted((u, v)))
        physical = pair in physical_pairs
        width = 0.65 + 1.35 * (float(data["weight"]) - 700) / 300
        ax.plot(
            [pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
            color="#59666D" if physical else "#AEB8BD",
            lw=width,
            ls="-" if physical else (0, (2.2, 2.2)),
            alpha=0.86 if physical else 0.72,
            zorder=1,
        )

    sizes = [78 + 38 * abs(float(effect.loc[gene, "logFC"])) for gene in members]
    ax.scatter(
        [pos[g][0] for g in members], [pos[g][1] for g in members],
        s=sizes, c=PALETTE[row.module_id], edgecolor="#333333", linewidth=0.65, zorder=3,
    )
    for gene in members:
        point = pos[gene]
        label_point = point * 1.48
        ax.text(label_point[0], label_point[1] + 0.055, gene, fontsize=7.4,
                fontstyle="italic", ha="center", va="center", color="#202020")
        ax.text(label_point[0], label_point[1] - 0.12, fmt_effect(effect.loc[gene, "logFC"]),
                fontsize=6.2, ha="center", va="center", color="#555555")

    code = row.module_id.replace("DOWN_", "D").replace("UP_", "U")
    ax.text(0.0, 1.22, f"{code}  {TITLES[row.module_id]}", transform=ax.transAxes,
            fontsize=8.2, weight="semibold", color=PALETTE[row.module_id], va="top")
    ax.text(0.0, 1.08, f"{len(members)} genes  |  MCODE score {row.mcode_score:.1f}",
            transform=ax.transAxes, fontsize=6.4, color="#555555", va="top")
    ax.set_xlim(-1.70, 1.70)
    ax.set_ylim(-1.65, 1.65)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    modules, effects, association, physical_pairs = load_data()
    assert len(modules) == 5
    assert set(modules["module_id"]) == set(PALETTE)
    assert effects["gene_name"].is_unique

    figure = plt.figure(figsize=(183 / 25.4, 158 / 25.4), facecolor="white")
    grid = figure.add_gridspec(6, 2, left=0.055, right=0.97, bottom=0.18, top=0.735,
                              width_ratios=[1.08, 0.92], hspace=1.00, wspace=0.28)
    down_axes = [figure.add_subplot(grid[0:2, 0]), figure.add_subplot(grid[2:4, 0]),
                 figure.add_subplot(grid[4:6, 0])]
    up_axes = [figure.add_subplot(grid[0:3, 1]), figure.add_subplot(grid[3:6, 1])]

    for ax, (_, row) in zip(down_axes, modules.loc[modules.direction.eq("DOWN")].iterrows()):
        draw_module(ax, row, effects, association, physical_pairs)
    for ax, (_, row) in zip(up_axes, modules.loc[modules.direction.eq("UP")].iterrows()):
        draw_module(ax, row, effects, association, physical_pairs)

    figure.text(0.055, 0.955, "Reference protein associations among persistent expression changes",
                fontsize=11.0, weight="semibold", va="top")
    figure.text(0.055, 0.910,
                "500 of 662 selected genes mapped to STRING v12.0  |  association score ≥0.7  |  five MCODE modules",
                fontsize=7.2, color="#555555", va="top")
    figure.text(0.040, 0.820, "a", fontsize=12, weight="semibold")
    figure.text(0.080, 0.820, "Lower in Y-depleted HT1080", fontsize=9.2, weight="semibold")
    figure.text(0.565, 0.820, "b", fontsize=12, weight="semibold")
    figure.text(0.605, 0.820, "Higher in Y-depleted HT1080", fontsize=9.2, weight="semibold")

    legend = figure.add_axes([0.08, 0.055, 0.84, 0.075])
    legend.axis("off")
    legend.plot([0.02, 0.09], [0.66, 0.66], color="#59666D", lw=1.6)
    legend.text(0.105, 0.66, "STRING physical support ≥0.7", va="center", fontsize=6.7)
    legend.plot([0.36, 0.43], [0.66, 0.66], color="#AEB8BD", lw=1.2, ls=(0, (2.2, 2.2)))
    legend.text(0.445, 0.66, "Other STRING association ≥0.7", va="center", fontsize=6.7)
    for x, value in zip([0.72, 0.80, 0.90], [1, 2, 4]):
        legend.scatter(x, 0.66, s=78 + 38 * value, c="#D7DDE0", edgecolor="#555555", linewidth=0.5)
        legend.text(x, 0.18, str(value), ha="center", fontsize=6.3)
    legend.text(0.81, 1.05, "|log₂FC|", ha="center", fontsize=6.7)
    legend.set_xlim(0, 1)
    legend.set_ylim(0, 1.2)
    figure.text(0.055, 0.018,
                "Numbers beneath genes are model log₂ fold changes (Y-depleted versus sgCTRL). Layout positions aid display only.",
                fontsize=6.0, color="#555555")

    stem = "ch03_main_ht1080_string_modules"
    exports = []
    for extension in ["pdf", "svg", "png"]:
        target = FIG / f"{stem}.{extension}"
        figure.savefig(target, dpi=600)
        exports.append({"file": target.name, "sha256": sha256(target)})
    preview = FIG / "previews" / f"{stem}_preview.png"
    preview.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(preview, dpi=180)
    plt.close(figure)

    module_members = (modules.assign(gene_name=modules["members"].str.split(";"))
                      .explode("gene_name").drop(columns="members"))
    source_nodes = module_members.merge(effects.drop(columns="direction"), on="gene_name", how="left")
    source_nodes.to_csv(SOURCE / "module_nodes.csv", index=False)
    selected = set(source_nodes["gene_name"])
    source_edges = association[association["gene_a"].isin(selected) & association["gene_b"].isin(selected)].copy()
    source_edges["physical_support"] = [tuple(sorted((a, b))) in physical_pairs
                                        for a, b in source_edges[["gene_a", "gene_b"]].itertuples(index=False)]
    membership = source_nodes.set_index("gene_name")["module_id"].to_dict()
    source_edges = source_edges[source_edges.apply(lambda r: membership[r.gene_a] == membership[r.gene_b], axis=1)]
    source_edges["module_id"] = source_edges["gene_a"].map(membership)
    source_edges.to_csv(SOURCE / "module_edges.csv", index=False)
    modules.to_csv(SOURCE / "modules.csv", index=False)

    provenance = {
        "analysis": "Current batch + condition directionally consistent set: overall BH FDR <0.05, post-hoc |log2FC| >=1, matching August and May directions.",
        "mapping": "500 of 662 selected symbols mapped exactly to STRING preferred symbols (331 DOWN; 169 UP).",
        "network": "STRING v12.0 association score >=700; MCODE settings inherited from the frozen analysis. All five detected modules and all within-module edges are shown. Solid edges also have physical-network score >=700.",
        "layout": "Deterministic circular layouts within modules. Positions and between-module distances have no biological meaning.",
        "source_files": {str(path.relative_to(ANALYSIS)): sha256(path) for path in [
            STRING / "modules.csv", STRING / "all_symbol_mapping_and_effects.csv",
            STRING / "all_candidate_association_edges.csv", STRING / "all_candidate_physical_edges.csv"]},
        "exports": exports,
        "script_sha256": sha256(Path(__file__)),
    }
    (SOURCE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"modules": len(modules), "nodes": len(source_nodes),
                      "edges": len(source_edges), "exports": exports}, indent=2))


if __name__ == "__main__":
    main()
