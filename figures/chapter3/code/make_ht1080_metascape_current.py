"""Render the 19 September 2026 HT1080 Metascape rerun in thesis style."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import textwrap
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from adjustText import adjust_text
import networkx as nx
import numpy as np
import pandas as pd


FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
RESEARCH = THESIS / "research"
SOURCE = FIG / "source_data" / "ht1080_metascape_current"
SOURCE.mkdir(parents=True, exist_ok=True)

import sys
sys.path.insert(0, str(THESIS / "skills" / "thesis-figure-style" / "scripts"))
from thesis_style import apply_style

apply_style(8)

NS = {"x": "http://www.cs.rpi.edu/XGMML"}
NPG = ["#3C5488", "#E64B35", "#00A087", "#F39B7F", "#4DBBD5", "#8491B4", "#7E6148"]
GREY = "#B8B8B8"

THEME_NAMES = {
    "up": {
        "GPCR ligand binding": "GPCR and neuroactive-ligand signalling",
        "regulation of anatomical structure morphogenesis": "Morphogenesis and chemotaxis",
        "epithelium development": "Epithelial and cardiovascular development",
        "neuron projection morphogenesis": "Neurite morphogenesis",
        "synapse organization": "Synapse organisation and neuronal signalling",
        "MIKKELSEN MCV6 HCP WITH H3K27ME3": "Polycomb-associated promoter signatures",
        "LIM MAMMARY STEM CELL UP": "Epithelial and stem-cell signatures",
    },
    "down": {
        "regulation of system process": "Adhesion and vascular-system processes",
        "extracellular matrix": "Stromal ECM and collagen turnover",
        "external encapsulating structure": "Stromal ECM and collagen turnover",
        "hormone binding": "Hormone and peptide-receptor binding",
        "cell morphogenesis": "Cell and neuron morphogenesis",
        "SWEET LUNG CANCER KRAS DN": "KRAS and CDH1-associated signatures",
        "WALLACE PROSTATE CANCER RACE UP": "Cancer and ESR1-reference signatures",
        "positive regulation of MAPK cascade": "RTK, PI3K and MAPK signalling",
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wrap(value: str, width: int) -> str:
    return "\n".join(textwrap.wrap(str(value), width))


def shorten(value: str, limit: int = 33) -> str:
    value = value.replace("MIKKELSEN ", "").replace(" HCP WITH H3K27ME3", " HCP H3K27me3")
    return value if len(value) <= limit else value[: limit - 3] + "..."


def parse_export(direction: str) -> tuple[pd.DataFrame, nx.Graph, str, Path]:
    archive = RESEARCH / f"{direction.upper()}_METASCAPE.zip"
    with ZipFile(archive) as bundle:
        data = bundle.read("Enrichment_GO/GONetwork.xgmml")
    root = ET.fromstring(data)
    date_node = root.find(".//{http://purl.org/dc/elements/1.1/}date")
    export_date = date_node.text if date_node is not None else "unknown"

    rows = []
    lookup = {}
    for node in root.findall("x:node", NS):
        attrs = {a.attrib["name"]: a.attrib.get("value") for a in node.findall("x:att", NS)}
        lookup[node.attrib["id"]] = attrs["canonicalName"]
        rows.append(
            {
                "name": attrs["canonicalName"],
                "description": attrs["Description"],
                "mlog10p": -float(attrs["LogP"]),
                "log10q": float(attrs.get("Log(q-value)", "nan")),
                "hits": attrs.get("Hits", ""),
                "metascape_group": attrs.get("GROUP_ID"),
            }
        )
    nodes = pd.DataFrame(rows)
    graph = nx.Graph()
    graph.add_nodes_from(nodes["name"])
    for edge in root.findall("x:edge", NS):
        source = lookup[edge.attrib["source"]]
        target = lookup[edge.attrib["target"]]
        if source != target:
            graph.add_edge(source, target)

    assert graph.number_of_nodes() == len(nodes)
    assert nodes["name"].is_unique
    return nodes, graph, export_date, archive


def annotate_communities(nodes: pd.DataFrame, graph: nx.Graph, direction: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    communities = nx.community.louvain_communities(graph, seed=42)
    membership = {node: index + 1 for index, community in enumerate(communities) for node in community}
    nodes = nodes.copy()
    nodes["community"] = nodes["name"].map(membership)

    summaries = []
    for community, group in nodes.groupby("community", sort=True):
        ordered = group.sort_values(["mlog10p", "description"], ascending=[False, True])
        top_term = ordered.iloc[0]["description"]
        label = THEME_NAMES[direction].get(top_term, "Other reference signatures")
        hit_counts = Counter(ordered["hits"].str.split("|").explode())
        hit_counts.pop("", None)
        recurring = [gene for gene, _ in sorted(hit_counts.items(), key=lambda pair: (-pair[1], pair[0]))[:4]]
        summaries.append(
            {
                "community": community,
                "theme": label,
                "terms": len(group),
                "best_mlog10p": ordered.iloc[0]["mlog10p"],
                "representative_pathways": "; ".join(shorten(x) for x in ordered.head(2)["description"]),
                "recurring_hit_genes": ", ".join(recurring),
                "top_term": top_term,
            }
        )
    summary = pd.DataFrame(summaries)
    named = summary["top_term"].isin(THEME_NAMES[direction])
    named_ids = summary.loc[named, "community"].tolist()
    color_map = {community: NPG[index] for index, community in enumerate(named_ids)}
    nodes["theme"] = nodes["community"].map(summary.set_index("community")["theme"])
    nodes["color"] = nodes["community"].map(color_map).fillna(GREY)
    summary["color"] = summary["community"].map(color_map).fillna(GREY)
    summary["displayed_theme"] = np.where(named, summary["theme"], "Other isolated reference signatures")
    return nodes, summary


def component_layout(graph: nx.Graph) -> dict[str, np.ndarray]:
    components = sorted(nx.connected_components(graph), key=lambda part: (-len(part), sorted(part)[0]))
    positions: dict[str, np.ndarray] = {}
    for index, component in enumerate(components):
        subgraph = graph.subgraph(component)
        if len(component) > 1:
            local = nx.spring_layout(subgraph, seed=42, k=1.1 / np.sqrt(len(component)), iterations=500)
        else:
            local = {next(iter(component)): np.zeros(2)}
        if index == 0:
            center = np.array([0.0, 0.0])
            scale = 1.0
        else:
            column = (index - 1) % 2
            row = (index - 1) // 2
            center = np.array([1.55 + 0.62 * column, 0.85 - 0.56 * row])
            scale = 0.18 if len(component) < 5 else 0.28
        for node, point in local.items():
            positions[node] = point * scale + center
    return positions


def community_layout(nodes: pd.DataFrame, graph: nx.Graph, direction: str) -> dict[str, np.ndarray]:
    # Display anchors only: communities and all exported connections are unchanged.
    centers = {
        'up': [(-0.80, 0.00), (-0.72, -0.68), (0.00, 0.12),
               (0.78, -0.25), (0.12, -0.72), (0.65, 0.72), (-0.64, 0.72)],
        'down': [(-0.65, 0.64), (-0.68, -0.34), (0.68, 0.69),
                 (0.12, -0.68), (-0.79, -0.85), (0.80, -0.62), (0.10, 0.08)],
    }
    positions = {}
    visible = nodes.loc[nodes['color'].ne(GREY)]
    for index, (_, group) in enumerate(visible.groupby('community', sort=True)):
        sub = graph.subgraph(group['name'])
        local = nx.spring_layout(sub, seed=42, iterations=700, k=0.65)
        radius = 0.15 + 0.23 * np.sqrt(len(group) / 33)
        center = np.asarray(centers[direction][index])
        for name, point in local.items():
            positions[name] = center + point * radius
    return positions


def collapse_summary(summary: pd.DataFrame) -> pd.DataFrame:
    named = summary[summary["color"].ne(GREY)].copy()
    grey = summary[summary["color"].eq(GREY)]
    if len(grey):
        aggregate = {
            "displayed_theme": "Other isolated reference signatures",
            "terms": int(grey["terms"].sum()),
            "best_mlog10p": float(grey["best_mlog10p"].max()),
            "representative_pathways": "; ".join(grey.nlargest(2, "best_mlog10p")["top_term"].map(shorten)),
            "recurring_hit_genes": "See source table",
            "color": GREY,
        }
        named = pd.concat([named, pd.DataFrame([aggregate])], ignore_index=True)
    return named


def draw_table(fig: plt.Figure, summary: pd.DataFrame) -> None:
    table = collapse_summary(summary)
    axis = fig.add_axes([0.045, 0.045, 0.91, 0.25])
    axis.set_axis_off()
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    x_positions = [0.025, 0.315, 0.395, 0.47, 0.79]
    headers = ["Theme", "Terms", "Best\n−log₁₀P", "Representative pathways", "Recurring hit genes"]
    for x, label in zip(x_positions, headers):
        axis.text(x, 0.99, label, fontsize=6.3, weight="semibold", va="top")
    axis.plot([0, 1], [0.85, 0.85], color="#333333", lw=0.55)
    row_height = 0.76 / len(table)
    for index, row in enumerate(table.itertuples()):
        y = 0.81 - index * row_height
        axis.scatter(0.009, y - 0.018, s=20, color=row.color)
        values = [
            wrap(row.displayed_theme, 28),
            str(row.terms),
            f"{row.best_mlog10p:.1f}",
            str(row.representative_pathways).replace("; ", "\n"),
            wrap(row.recurring_hit_genes, 24),
        ]
        for x, value in zip(x_positions, values):
            axis.text(x, y, value, fontsize=5.9, va="top")


def render(direction: str) -> tuple[plt.Figure, dict]:
    nodes, graph, export_date, archive = parse_export(direction)
    nodes, summary = annotate_communities(nodes, graph, direction)
    position = community_layout(nodes, graph, direction)
    nodes['displayed'] = nodes['name'].isin(position)
    nodes["x"] = nodes["name"].map(lambda name: position.get(name, [np.nan, np.nan])[0])
    nodes["y"] = nodes["name"].map(lambda name: position.get(name, [np.nan, np.nan])[1])
    shown = nodes.loc[nodes['displayed']].copy()
    display_graph = graph.subgraph(shown['name'])

    figure = plt.figure(figsize=(183 / 25.4, 230 / 25.4))
    axis = figure.add_axes([0.05, 0.34, 0.90, 0.57])
    axis.set_axis_off()
    nx.draw_networkx_edges(display_graph, position, ax=axis, width=0.40, edge_color="#B9B9B9", alpha=0.38)
    sizes = 12 + shown["mlog10p"].to_numpy() * 4
    nx.draw_networkx_nodes(
        display_graph,
        position,
        ax=axis,
        nodelist=shown["name"].tolist(),
        node_size=sizes,
        node_color=shown["color"].tolist(),
        edgecolors="#555555",
        linewidths=0.4,
    )
    theme_labels = []
    for row in summary.loc[summary["color"].ne(GREY)].itertuples():
        group = nodes[nodes["community"].eq(row.community)]
        theme_labels.append(axis.text(
            group["x"].mean(),
            group["y"].max() + 0.10,
            wrap(row.theme, 24),
            fontsize=8,
            weight="semibold",
            ha="center",
            va="center",
            bbox={"facecolor": "white", "alpha": 0.70, "edgecolor": "none", "pad": 1},
        ))
    np.random.seed(42)
    adjust_text(
        theme_labels,
        ax=axis,
        x=shown["x"].to_numpy(),
        y=shown["y"].to_numpy(),
        arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.45},
        iter_lim=300,
    )
    axis.margins(0.10)

    condition = "Higher" if direction == "up" else "Lower"
    figure.text(
        0.055,
        0.958,
        (
            f"{condition} in Y-depleted HT1080\n"
            f"{len(shown)} displayed terms | {display_graph.number_of_edges()} connections"
        ),
        fontsize=10,
        weight="semibold",
        linespacing=1.8,
    )
    handles = [axis.scatter([], [], s=12 + value * 4, c=GREY, edgecolors="#555555", linewidths=0.4) for value in [5, 10, 15]]
    figure.legend(
        handles=handles,
        labels=["5", "10", "15"],
        title=r"$-\log_{10}P$",
        loc="lower center",
        bbox_to_anchor=(0.5, 0.30),
        ncol=3,
        frameon=False,
        fontsize=7,
        title_fontsize=7,
    )
    draw_table(figure, summary)
    figure.text(
        0.055,
        0.019,
        ("Six isolated terms omitted from network; retained in source tables. " if direction == 'down' else '')
        + "Positions aid display; connections represent exported term similarity.",
        fontsize=5.8,
        color="#555555",
    )

    stem = f"ch03_main_ht1080_metascape_current_{direction}_network_v2"
    exports = []
    for extension in ["pdf", "svg", "png"]:
        target = FIG / f"{stem}.{extension}"
        staging = FIG / f"{stem}.staging.{extension}"
        figure.savefig(staging, dpi=600)
        staging.replace(target)
        exports.append({"file": target.name, "sha256": sha256(target)})
    figure.savefig(FIG / "previews" / f"{stem}_preview.png", dpi=150)

    nodes.to_csv(SOURCE / f"{direction}_nodes.csv", index=False)
    pd.DataFrame(graph.edges(), columns=["source", "target"]).to_csv(SOURCE / f"{direction}_edges.csv", index=False)
    summary.to_csv(SOURCE / f"{direction}_communities.csv", index=False)
    record = {
        "direction": direction,
        "nodes": len(nodes),
        "displayed_nodes": len(shown),
        "omitted_isolated_nodes": int((~nodes['displayed']).sum()),
        "edges": graph.number_of_edges(),
        "components": nx.number_connected_components(graph),
        "communities": len(summary),
        "export_date": export_date,
        "source_archive": str(archive),
        "source_sha256": sha256(archive),
        "exports": exports,
    }
    return figure, record


def main() -> None:
    figures = []
    records = []
    for direction in ["up", "down"]:
        figure, record = render(direction)
        figures.append(figure)
        records.append(record)

    combined = FIG / "ch03_supp_ht1080_metascape_current_networks_v2.pdf"
    with PdfPages(combined) as pdf:
        for figure in figures:
            pdf.savefig(figure)
            plt.close(figure)
    manifest = {
        "analysis": (
            "Metascape rerun using the current limma-voom batch + condition lists: "
            "overall BH FDR <0.05, post-hoc absolute log2FC >=1, matching period directions; "
            "25,461-symbol custom tested background."
        ),
        "network": (
            "All nodes and unique edges are retained in source tables. Six grey isolated DOWN terms "
            "are omitted from display. Communities are NetworkX Louvain communities with seed 42. "
            "Layout uses curated community anchors and seeded within-community spring layouts; "
            "positions have no biological meaning."
        ),
        "runs": records,
        "combined_pdf": {"file": combined.name, "sha256": sha256(combined)},
        "script_sha256": sha256(Path(__file__)),
    }
    (SOURCE / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
