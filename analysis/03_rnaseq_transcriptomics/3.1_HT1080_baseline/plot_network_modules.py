"""Publication figures and a browsable evidence report for fresh HT1080 modules."""
from __future__ import annotations

import html
import json
import math
import re
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd
from adjustText import adjust_text
from PIL import Image

from network_modules_analysis import ROOT, TABLES, PROV, digest

FIGURES = ROOT / "outputs/figures/network_modules"
CAPTIONS = ROOT / "outputs/captions/network_modules_captions.md"
COLORS = ["#0072B2", "#009E73", "#CC79A7", "#D55E00", "#6A51A3", "#8C6D31", "#56B4E9", "#B88A00"]


def short_term(name):
    # Display changes only; exact source name and membership remain in tables.
    return re.sub(r"^(GOBP_|REACTOME_)", "", name).replace("_", " ").lower()


def fmt_q(value):
    if not np.isfinite(value):
        return "—"
    return f"{value:.1e}" if value < .001 else f"{value:.3f}"


def symbol_size(effect):
    return 12 + 8 * abs(effect)


def separate_remaining_labels(ax, labels):
    """Resolve residual bounding-box intersections deterministically in pixels."""
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    bounds = ax.get_window_extent(renderer)
    for _ in range(150):
        changed = False
        for i, a in enumerate(labels):
            for b in labels[i+1:]:
                ra = a.get_window_extent(renderer).expanded(1.05, 1.18)
                rb = b.get_window_extent(renderer).expanded(1.05, 1.18)
                if not ra.overlaps(rb):
                    continue
                ca, cb = (ra.y0+ra.y1)/2, (rb.y0+rb.y1)/2
                distance = (ra.height+rb.height)/2 + 1.5 - abs(ca-cb)
                sign = 1 if ca >= cb else -1
                for text, box, step in [(a, ra, sign*distance/2), (b, rb, -sign*distance/2)]:
                    pixel = ax.transData.transform(text.get_position())
                    step = np.clip(step, bounds.y0+1-box.y0, bounds.y1-1-box.y1)
                    pixel[1] += step
                    text.set_position(ax.transData.inverted().transform(pixel))
                changed = True
        if not changed:
            return
    raise RuntimeError("Label placement did not converge without overlap")


def module_positions(graph, members):
    # Use explicit order, fixed initialization and unweighted geometry. STRING
    # score controls edge width, not placement or MCODE topology after filtering.
    sub = nx.Graph()
    sub.add_nodes_from(sorted(members))
    sub.add_edges_from(sorted(tuple(sorted((a, b))) for a, b in graph.subgraph(members).edges()))
    positions = nx.spring_layout(sub, seed=42, iterations=500, weight=None, scale=1)
    # Each module has a fixed display cell. Geometry is not a biological axis.
    xy = np.array([positions[g] for g in members])
    xy[:, 0] = (xy[:, 0] - xy[:, 0].mean()) * .80
    xy[:, 1] = (xy[:, 1] - xy[:, 1].mean()) * .70
    positions = dict(zip(members, xy))
    return positions, xy


def plot_module(ax, row, graph, genes, label, color, position_rows, page):
    members = row.members.split(";")
    positions, xy = module_positions(graph, members)
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.16, 1.18)
    ax.axis("off")
    ax.text(0, 1.08, f"{row.module_id.replace('DOWN_', 'D').replace('UP_', 'U')}  {label}",
            transform=ax.transAxes, fontsize=6.5, fontweight="bold", color=color,
            ha="left", va="bottom", clip_on=False)
    for a, b in sorted(graph.subgraph(members).edges()):
        attrs = graph[a][b]
        physical = attrs["physical_supported"]
        width = .3 + .35 * (attrs["association_score"] - 700) / 300
        ax.plot([positions[a][0], positions[b][0]], [positions[a][1], positions[b][1]],
                color="#596A72" if physical else "#B8BEC2", linewidth=width,
                linestyle="-" if physical else (0, (2.5, 2.2)), alpha=.8,
                zorder=1)
    effects = [genes.loc[g, "mean_difference"] for g in members]
    points = ax.scatter(xy[:, 0], xy[:, 1], s=[symbol_size(v) for v in effects],
                        c=color, edgecolors="#343D42", linewidths=.35, zorder=3)
    labels = [ax.text(positions[g][0], positions[g][1], g, fontsize=5.8,
                      fontstyle="italic", ha="center", va="center", zorder=5)
              for g in members]
    np.random.seed(42)
    adjust_text(labels, x=xy[:, 0], y=xy[:, 1], objects=points, ax=ax,
                force_text=(.35, .55), force_static=(.3, .5), force_pull=(.015, .015),
                expand=(1.18, 1.45), max_move=(8, 12), iter_lim=600,
                ensure_inside_axes=True, prevent_crossings=True)
    separate_remaining_labels(ax, labels)
    for gene, text in zip(members, labels):
        delta = np.linalg.norm(ax.transData.transform(text.get_position()) - ax.transData.transform(positions[gene]))
        if delta > 5:
            ax.annotate("", xy=positions[gene], xytext=text.get_position(),
                        arrowprops=dict(arrowstyle="-", color="#8D969B", lw=.25,
                                        shrinkA=2, shrinkB=3), zorder=2)
        position_rows.append({"page": page, "module_id": row.module_id, "gene_name": gene,
                              "node_x": positions[gene][0], "node_y": positions[gene][1],
                              "label_x": text.get_position()[0], "label_y": text.get_position()[1]})
    return labels


def draw_page(rows, direction, page_number, total_modules, graphs, genes, labels, config, position_rows):
    full_height = config["figure_height_mm"]
    unused_network_mm = (4-math.ceil(len(rows)/2))*.155*full_height
    unused_table_mm = (8-len(rows))*.103*.193*full_height
    height = full_height-unused_network_mm-unused_table_mm
    table_height_mm = .193*full_height-unused_table_mm
    top_y = lambda y, removed=0: 1-((1-y)*full_height-removed)/height
    bottom_y = lambda y: y*full_height/height
    fig = plt.figure(figsize=(config["figure_width_mm"]/25.4, height/25.4), facecolor="white")
    direction_title = "Lower RNA in LOY" if direction == "DOWN" else "Higher RNA in LOY"
    shown = "Highest-scoring eight modules" if page_number == 1 else f"Module atlas · page {page_number}"
    fig.text(.035, top_y(.978), f"{direction_title}: protein-association modules", fontsize=9, fontweight="bold", va="top")
    fig.text(.035, top_y(.953), f"{shown} · {len(rows)} of {total_modules} modules · STRING v12.0 · score ≥ 0.7", fontsize=6.5, color="#505B60", va="top")
    fig.text(.011, top_y(.926), "a", fontsize=8, fontweight="bold")
    plotted_labels = []
    for i, row in enumerate(rows.itertuples()):
        col, grid_row = i % 2, i // 2
        ax = fig.add_axes([.055 + col*.485, top_y(.785-grid_row*.155), .41, .121*full_height/height])
        plotted_labels.extend(plot_module(ax, row, graphs[direction], genes,
                                         labels[row.module_id], COLORS[i], position_rows,
                                         f"{direction}_{page_number}"))
    # The same size mapping is used for both directions and every atlas page.
    size_handles = [plt.scatter([], [], s=symbol_size(v), facecolors="#D4DADF", edgecolors="#343D42", linewidths=.35) for v in [1, 2, 4]]
    fig.legend(size_handles, ["1", "2", "4"], title="|mean log$_2$ RNA difference|", loc="center left",
               bbox_to_anchor=(.035, top_y(.282, unused_network_mm)), ncol=3, frameon=False, fontsize=6, title_fontsize=6,
               columnspacing=1.1, handletextpad=.45)
    line_handles = [Line2D([0], [0], color="#596A72", lw=.65),
                    Line2D([0], [0], color="#B8BEC2", lw=.65, linestyle=(0, (2.5, 2.2)))]
    fig.legend(line_handles, ["Physical support ≥ 0.7", "Other STRING association"],
               loc="center left", bbox_to_anchor=(.52, top_y(.282, unused_network_mm)), frameon=False, fontsize=6)
    fig.text(.011, top_y(.242, unused_network_mm), "b", fontsize=8, fontweight="bold")
    table_ax = fig.add_axes([.035, bottom_y(.055), .93, table_height_mm/height])
    table_ax.axis("off")
    table_ax.set_xlim(0, 1)
    table_ax.set_ylim(0, table_height_mm)
    table_y = lambda y: table_height_mm-(1-y)*.193*full_height
    cols = [.00, .075, .13, .62, .72]
    for x, header in zip(cols, ["Module", "n", "Top source annotation", "BH q", "Largest RNA differences"]):
        table_ax.text(x, table_y(.96), header, fontsize=6, fontweight="bold", va="top")
    table_ax.axhline(table_y(.99), color="#35434B", lw=.65)
    table_ax.axhline(table_y(.88), color="#35434B", lw=.4)
    for i, row in enumerate(rows.itertuples()):
        y = table_y(.82-i*.103)
        module_id = row.module_id.replace("DOWN_", "D").replace("UP_", "U")
        table_ax.text(cols[0], y, module_id, fontsize=6, fontweight="bold", color=COLORS[i], va="center")
        table_ax.text(cols[1], y, str(row.nodes), fontsize=6, va="center")
        source = short_term(row.pathway) if isinstance(row.pathway, str) else "No eligible annotation"
        wrapped = textwrap.wrap(source, width=51)
        if len(wrapped)>2:
            wrapped = [wrapped[0], textwrap.shorten(" ".join(wrapped[1:]), width=50, placeholder="…")]
        table_ax.text(cols[2], y, "\n".join(wrapped), fontsize=5.5, va="center", linespacing=1.05)
        table_ax.text(cols[3], y, fmt_q(row.BH_all_modules_and_terms), fontsize=5.7, va="center")
        names = sorted(row.members.split(";"), key=lambda g: (-abs(genes.loc[g,"mean_difference"]),g))
        member_text = ", ".join(names[:3]) + (f" +{len(names)-3}" if len(names)>3 else "")
        table_ax.text(cols[4], y, "\n".join(textwrap.wrap(member_text, width=28)), fontsize=5.5,
                      fontstyle="italic", va="center", linespacing=1.05)
    table_ax.axhline(.025*.193*full_height, color="#35434B", lw=.65)
    fig.text(.035, bottom_y(.039), "MCODE: local implementation; modules ranked by topology. Only within-module edges shown.", fontsize=5.7, color="#505B60")
    fig.text(.035, bottom_y(.025), "Annotation q spans all 42 modules × 6,203 terms; exploratory. RNA changes and reference edges do not measure protein activity.", fontsize=5.5, color="#505B60")
    return fig, plotted_labels


def caption(direction, rows, coverage):
    word = "decreased" if direction == "DOWN" else "increased"
    return (f"**Protein-association modules among genes with {word} RNA in established HT1080 LOY populations.** "
        f"The eight highest-scoring of {coverage['modules']} modules are displayed ({int(rows.nodes.sum())} genes). "
        "Candidates satisfy the existing absolute combined log2-expression difference ≥ 0.5, period-direction agreement and six leave-one-match-out direction checks. "
        "Combined effects equally weight August 2025 and May 2026. "
        "a, Human STRING v12.0 associations with scores ≥ 0.7 were clustered using a local Python implementation of MCODE "
        "(degree cutoff 2; node-score cutoff 0.2; k-core 2; haircut enabled; fluff disabled; maximum depth 100). "
        "Colors distinguish modules; node sizes scale with absolute combined RNA differences. Solid edges also have STRING physical-network scores ≥ 0.7; "
        "dashed edges denote other STRING associations. Only within-module edges are shown; layout distances carry no biological meaning. "
        "b, Top GO biological-process/Reactome annotations and representative genes. Annotation BH adjustment covers all 260,526 module–term tests "
        "across both directions, using 14,564 expressed STRING-mapped genes as background. These exploratory annotations share knowledge with the network resource. "
        "The complete atlas and candidate-background sensitivity are supplied separately. Six timing matches repeatedly sample the same established populations; "
        "RNA changes and reference associations do not establish altered protein abundance, binding or pathway activity.")


def make_report(modules, coverage, labels, figure_rows, config):
    body = []
    body.append("<h1>HT1080 LOY: fresh protein-association modules</h1><p>14 September 2026 · SEPT_11_2026_REANALYZE</p>")
    body.append("<p><strong>The current candidate genes form 24 lower-RNA modules (150 genes) and 18 higher-RNA modules (86 genes).</strong> "
                "The leading lower-RNA modules include interferon-response genes, an integrin–collagen–matrix group and an ankyrin–spectrin scaffold. "
                "The higher-RNA network includes receptor, Hedgehog-associated, Eph–ephrin and other programs. These describe selected groups of RNA changes, "
                "not measured protein complexes or biochemical pathway activation.</p>")
    body.append("<p>The updated matrix module includes MMP1, SPP1, collagens and integrins. A different higher-RNA matrix module contains MMP14/MMP15, "
                "COL2A1/COL6A2 and other genes (UP_17, atlas page 3). The presence of both directions reinforces selective remodeling rather than a uniform decrease "
                "of every matrix gene. Similarly, opposing receptor modules contain different subunits. Source labels such as heart development, photoreceptor "
                "maintenance or immune pathways do not demonstrate acquisition of those cell identities in HT1080.</p>")
    body.append("<h2>Coverage and selection</h2><p>The frozen analysis has 4,892 candidate gene loci. Applying its existing representative-symbol mapping gives "
                "4,879 unique-symbol candidates (2,379 lower and 2,500 higher); 2,288 map exactly to STRING preferred names. No fuzzy aliases, external connector "
                "proteins or legacy gene-level significance thresholds were introduced. Unmapped and unclustered genes remain in the evidence tables. "
                "Module ranks reflect density × module size, not annotation P values. The main figures show the first eight ranks in each direction; the atlas retains all ranks.</p>")
    body.append(coverage.to_html(index=False, border=0))
    legacy = pd.read_csv(TABLES / "legacy_DOWN_members_current_evidence.csv")
    body.append(f"<p>For comparison with the supplied historical figure, all {len(legacy)} pictured lower-RNA genes map to the current expression table, "
                f"{int(legacy.mean_difference.lt(0).sum())} remain lower and {int(legacy.descriptive_candidate.sum())} meet the current descriptive-candidate criteria. "
                "This supports continuity of the gene-level pattern. Module membership can differ because this extension uses a new candidate list and STRING "
                "network rather than the historical Metascape network and annotations.</p>")
    body.append("<h2>Main figures</h2>")
    for row in figure_rows:
        if row["page"] != 1:
            continue
        stem = row["stem"]
        body.append(f'<h3>{html.escape(row["title"])}</h3><p><a href="../figures/network_modules/{stem}.pdf">PDF</a> · '
                    f'<a href="../figures/network_modules/{stem}.svg">SVG</a> · <a href="../figures/network_modules/{stem}.png">600 dpi PNG</a></p>')
        body.append("<svg" + (FIGURES / f"{stem}.svg").read_text(encoding="utf-8").split("<svg",1)[1])
        direction = row["direction"]
        cap = caption(direction, modules[(modules.direction==direction) & modules.shown_main], coverage.set_index("direction").loc[direction].to_dict())
        body.append(f'<p class="caption">{html.escape(cap.replace("**", ""))}</p>')
    body.append("<h2>Complete module atlas</h2><p><a href='../figures/network_modules/HT1080_LOY_network_module_atlas.pdf'>Download all six atlas pages as one PDF</a>. "
                "Every module has its full membership below; all page exports are in the figure manifest.</p>")
    for direction in ("DOWN", "UP"):
        sub = modules[modules.direction==direction].copy()
        sub["label"] = sub.module_id.map(labels)
        sub["physical_support"] = sub.physical_edges.astype(str)+"/"+sub.edges.astype(str)
        sub["annotation_q"] = sub.BH_all_modules_and_terms.map(fmt_q)
        sub["candidate_background_q"] = sub.candidate_background_BH.map(fmt_q)
        sub["source_annotation"] = sub.pathway.map(short_term)
        body.append(f"<h3>{direction}</h3>" + sub[["module_id","label","nodes","mcode_score","physical_support", "source_annotation",
                    "annotation_q","candidate_background_q","members"]].round(3).to_html(index=False,border=0))
    body.append("<h2>What the network supports</h2><p>The lower-RNA interferon module is strongly connected in the association network but has only 1 of 45 "
                "edges supported at the same threshold in the physical network. By contrast, the ankyrin–spectrin module has physical support for 13 of 14 edges. "
                "This is why functional association and physical interaction are drawn differently. Even physical support may represent membership in the same complex "
                "rather than direct pairwise binding. Neither edge type was measured in these HT1080 samples.</p>")
    body.append("<p>The strict physical-only analysis identifies 12 lower-RNA and 11 higher-RNA modules; memberships are exported separately. "
                "The broad association modules must not be called independently validated protein complexes. Functional annotations share information sources "
                "with STRING, so annotation enrichment is descriptive support, not independent confirmation of module biology.</p>")
    body.append("<p>These are networks of genes meeting the RNA candidate criteria, not a map of all chromatin-regulating proteins. "
                "HDAC1, HDAC2, ZMYND8, EED, EZH2 and SUZ12 have absolute combined RNA differences below the 0.5 selection threshold; "
                "KDM5D and UTY are outside the autosomal/non-PAR-X primary universe. Their absence from these modules does not test their regulatory involvement.</p>")
    body.append("<h2>Methods and limits</h2><p>Both networks use the existing descriptive candidate flag and the same signed mean log2 RNA differences as the primary analysis. "
                "The plotted network is the STRING v12.0 full association network, threshold 700/1000; the separately cached physical network supplies physical support "
                "at threshold 700. Network files were cached on 9 September 2026 and verified against their recorded MD5 hashes; current inputs are SHA-256 frozen. "
                "Edges are undirected, duplicate pairs collapsed, self-loops excluded. Clustering uses binary topology after filtering; edge widths vary with association scores. "
                "Node area in points squared is 12 + 8 × absolute RNA difference. Layout is seeded Fruchterman–Reingold (seed 42), calculated independently for each module, "
                "with fixed geometry across copies of the same page. Local alphabetical tie/neighbor order can differ from Cytoscape node order. "
                "This was not a Metascape submission or an execution of the Cytoscape MCODE app.</p>")
    body.append("<p>Functional annotation uses frozen MSigDB 2026.1.Hs GO:BP and Reactome membership, restricted to 10–1,000 genes in the 14,564-gene expressed "
                "STRING-mapped background. All 6,203 eligible terms are tested against all 42 modules, including zero-overlap tests. One-sided hypergeometric P values "
                "receive BH correction over 260,526 tests. A second, separate family uses the 1,306 DOWN or 982 UP mapped candidates as the background, "
                "retaining the same terms and modules. This addresses how much annotation exceeds the directional input-list composition. "
                "The displayed top term has at least two overlapping genes; all tests are exported. These annotation q values are separate from the original GSEA/CAMERA "
                "families and are not gene-level or population-level significance. MCODE topology scores are not P values.</p>")
    body.append("<p>Candidate selection already requires period and omission direction agreement, so recurrence within these selected modules is partly built into selection. "
                "Module means and normalization/omission diagnostics are descriptive. There is one established population per condition, sampled again in May; "
                "the networks cannot isolate a general causal LOY effect. No RNA, chromatin or protein experiment was added.</p>")
    body.append("<h2>Source data and reproducibility</h2><ul>")
    for filename, title in [("modules.csv","All module memberships and statistics"), ("candidate_node_evidence.csv","Every mapped candidate, RNA evidence and clustering status"),
                            ("all_symbol_mapping_and_effects.csv","Full expression-universe mapping, including unmapped symbols"),
                            ("module_annotations_all_tests.csv.gz","All annotation tests and both backgrounds"),
                            ("strict_physical_modules.csv","Physical-only MCODE modules"),
                            ("module_expression_sensitivity.csv","Module expression sensitivity"),
                            ("legacy_DOWN_members_current_evidence.csv","Original example's members joined to current expression evidence"),
                            ("node_and_label_positions.csv","Plot positions")]:
        body.append(f'<li><a href="../tables/network_modules/{filename}">{title}</a></li>')
    body.append("</ul><p><a href='../figures/network_modules/manifest.csv'>Figure manifest</a> · "
                "<a href='../captions/network_modules_captions.md'>Captions</a> · "
                "<a href='../../provenance/network_modules/analysis_summary.json'>Analysis summary</a></p>")
    body.append("<p>From the RNAseq repository root:</p><pre>python analysis/HT1080/SEPT_11_2026_REANALYZE/network_modules_analysis.py\n"
                "python analysis/HT1080/SEPT_11_2026_REANALYZE/plot_network_modules.py</pre>")
    body.append("<p>Methods: <a href='https://pmc.ncbi.nlm.nih.gov/articles/PMC149346/'>Bader and Hogue, 2003</a>; "
                "<a href='https://github.com/BaderLab/MCODE/tree/a5657b9281b709ecd09f9747ef110f7f156c6cdb'>pinned MCODE reference implementation</a>; "
                "<a href='https://mcode.readthedocs.io/en/latest/Tuning.html'>MCODE parameters</a>; "
                "<a href='https://version-12-0.string-db.org/help/scores/'>STRING association and physical evidence</a>.</p>")
    document = "<!doctype html><html><head><meta charset='utf-8'><title>HT1080 LOY network modules</title><style>" \
        "body{font-family:Arial,sans-serif;max-width:1250px;margin:35px auto;padding:0 25px;line-height:1.55;color:#24313b}" \
        "h1,h2,h3{color:#17334a}table{border-collapse:collapse;width:100%;font-size:12px;display:block;overflow-x:auto}" \
        "th,td{padding:7px;border-bottom:1px solid #d6dce0;text-align:left}th{background:#eef2f4}svg{width:100%;height:auto;max-width:950px}" \
        ".caption{font-size:13px}pre{background:#eef2f4;padding:15px}a{color:#076593}</style></head><body>" + "\n".join(body) + "</body></html>"
    (ROOT / "outputs/report/network_modules.html").write_text(document, encoding="utf-8")


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    config = json.loads((ROOT / "config/network_modules.json").read_text())
    labels = json.loads((ROOT / "config/network_module_labels.json").read_text(encoding="utf-8"))
    font = font_manager.findfont("Arial", fallback_to_default=False)
    plt.rcParams.update({"font.family":"Arial", "font.size":6, "pdf.fonttype":42, "ps.fonttype":42,
                         "mathtext.fontset":"custom", "mathtext.rm":"Arial", "mathtext.it":"Arial:italic", "mathtext.bf":"Arial:bold",
                         "svg.fonttype":"none", "svg.hashsalt":"HT1080_network_modules_42", "axes.unicode_minus":True})
    modules = pd.read_csv(TABLES / "modules.csv")
    coverage = pd.read_csv(TABLES / "coverage.csv")
    nodes = pd.read_csv(TABLES / "candidate_node_evidence.csv").set_index("gene_name")
    graphs = {d:nx.read_graphml(ROOT / f"outputs/models/network_modules/STRING_association_{d}.graphml") for d in ("DOWN","UP")}
    assert set(modules.module_id) == set(labels)
    positions, manifest, label_checks, caption_blocks = [], [], [], []
    with PdfPages(FIGURES / "HT1080_LOY_network_module_atlas.pdf") as atlas:
        for direction in ("DOWN","UP"):
            sub = modules[modules.direction==direction].sort_values("rank")
            for first in range(0,len(sub),config["atlas_modules_per_page"]):
                page = first//config["atlas_modules_per_page"]+1
                selected = sub.iloc[first:first+config["atlas_modules_per_page"]]
                stem = (f"HT1080_LOY_network_modules_{direction}" if page==1 else f"HT1080_LOY_network_atlas_{direction}_page{page:02d}")
                fig, text_labels = draw_page(selected,direction,page,len(sub),graphs,nodes,labels,config,positions)
                fig.canvas.draw()
                renderer = fig.canvas.get_renderer()
                overlap_count = 0
                for i,a in enumerate(text_labels):
                    for b in text_labels[i+1:]:
                        if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)):
                            overlap_count += 1
                label_checks.append({"page":stem,"gene_labels":len(text_labels),"overlapping_label_pairs":overlap_count})
                assert overlap_count == 0, f"Overlapping gene labels on {stem}"
                for ext in ("pdf","svg","png"):
                    fig.savefig(FIGURES / f"{stem}.{ext}", dpi=config["png_dpi"], facecolor="white")
                atlas.savefig(fig)
                plt.close(fig)
                manifest.append({"stem":stem,"direction":direction,"page":page,"title":f"{direction}: module ranks {first+1}–{first+len(selected)}",
                                 "modules":";".join(selected.module_id),"genes":int(selected.nodes.sum()),
                                 "width_mm":config["figure_width_mm"],"height_mm":float(fig.get_size_inches()[1]*25.4),"png_dpi":config["png_dpi"],
                                 "review_status":"pending_visual_review"})
                if page==1:
                    caption_blocks.append(f"## {stem}\n\n" + caption(direction,selected,coverage.set_index("direction").loc[direction].to_dict()))
                print(f"Rendered {stem}: {len(text_labels)} gene labels; overlap pairs {overlap_count}",flush=True)
    pd.DataFrame(positions).to_csv(TABLES / "node_and_label_positions.csv",index=False)
    pd.DataFrame(manifest).to_csv(FIGURES / "manifest.csv",index=False)
    CAPTIONS.write_text("# HT1080 LOY protein-association module captions\n\n"+"\n\n".join(caption_blocks)+
        "\n\n## Complete atlas\n\nThe six-page atlas displays all 24 lower-RNA and 18 higher-RNA modules in descending MCODE-score order. "
        "Page 1 of each direction matches its main figure. Module identifiers remain fixed throughout; colors distinguish modules within a page and repeat on later pages. "
        "Methods, visual encodings, annotation correction and experimental limitations are the same as in the main captions.\n",encoding="utf-8")
    make_report(modules,coverage,labels,manifest,config)
    checks = {"font_file":font,"font_size_gene_labels_pt":5.8,"main_canvas_mm":[config["figure_width_mm"],config["figure_height_mm"]],
              "atlas_page_heights_mm":[row["height_mm"] for row in manifest],
              "label_checks":label_checks,"source_modules":len(modules),"source_genes":int(modules.nodes.sum()),
              "figure_pages":len(manifest),"matplotlib_version":matplotlib.__version__,"networkx_version":nx.__version__}
    for row in manifest:
        with Image.open(FIGURES / f"{row['stem']}.png") as im:
            assert abs(im.info['dpi'][0]-600)<1
            assert im.width==round(config["figure_width_mm"]/25.4*600) or abs(im.width-config["figure_width_mm"]/25.4*600)<1
    (PROV / "figure_export_checks.json").write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2),flush=True)


if __name__ == "__main__":
    main()
