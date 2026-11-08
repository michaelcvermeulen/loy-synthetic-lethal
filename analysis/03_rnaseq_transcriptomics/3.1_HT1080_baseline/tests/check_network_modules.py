"""Validate source preservation, graph memberships, statistics and figure exports."""
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote, urlparse

import networkx as nx
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from network_modules_analysis import digest, check_algorithm
from plot_network_modules import module_positions

REPO = ROOT.parents[2]
T = ROOT / "outputs/tables/network_modules"
F = ROOT / "outputs/figures/network_modules"
P = ROOT / "provenance/network_modules"
config = json.loads((ROOT / "config/network_modules.json").read_text())
summary = json.loads((P / "analysis_summary.json").read_text())
check_algorithm(config["mcode"])
for row in json.loads((P / "frozen_scientific_input_hashes.json").read_text()):
    assert digest(REPO / row["file"]) == row["sha256"], row["file"]
for row in json.loads((ROOT / "resources/network_modules/reference_sources.json").read_text()):
    assert digest(ROOT / "resources/network_modules" / row["file"]) == row["sha256"]

modules = pd.read_csv(T / "modules.csv")
nodes = pd.read_csv(T / "candidate_node_evidence.csv").fillna({"module_id": ""})
coverage = pd.read_csv(T / "coverage.csv").set_index("direction")
expression = pd.read_csv(ROOT / "outputs/tables/HT1080_LOY_vs_sgCTRL_primary_effects.csv").set_index("gene_id")
assert len(modules) == 42 and len(nodes) == 2288
assert nodes.gene_name.is_unique and nodes.gene_id.is_unique
assert nodes.descriptive_candidate.all() and nodes.mean_difference.abs().ge(.5).all()
assert nodes.period_direction_agrees.all() and nodes.LOO_all_same_direction.all()
for field in ["mean_difference", "AUG_2025", "MAY_2026", "mean_CPM", "prior025_mean_difference", "RLE_mean_difference"]:
    assert np.allclose(nodes[field], expression.loc[nodes.gene_id, field]), field
module_members = set()
for direction in ("DOWN", "UP"):
    graph = nx.read_graphml(ROOT / f"outputs/models/network_modules/STRING_association_{direction}.graphml")
    edge_table = pd.read_csv(T / f"{direction}_association_edges.csv")
    assert graph.number_of_edges() == len(edge_table) == coverage.loc[direction, "association_edges"]
    assert nx.number_of_selfloops(graph) == 0
    assert set(graph) == set(nodes.loc[nodes.direction.eq(direction), "gene_name"])
    assert edge_table.association_score.ge(700).all()
    for row in edge_table.itertuples():
        assert graph[row.gene_a][row.gene_b]["physical_supported"] == row.physical_supported
        assert row.physical_supported == pd.notna(row.physical_score)
        if row.physical_supported:
            assert row.physical_score >= 700
    for row in modules[modules.direction.eq(direction)].itertuples():
        members = set(row.members.split(";"))
        assert not members & module_members
        module_members |= members
        sub = graph.subgraph(members)
        assert len(sub) == row.nodes and sub.number_of_edges() == row.edges
        assert nx.is_connected(sub) and len(nx.k_core(sub, k=2)) == len(sub)
        assert np.isclose(nx.density(sub)*len(sub), row.mcode_score)
        assert sum(d["physical_supported"] for _, _, d in sub.edges(data=True)) == row.physical_edges
        assert members == set(nodes.loc[nodes.module_id.eq(row.module_id), "gene_name"])
    # Layout must survive different node/edge insertion orders, not just seed reuse.
    members = modules[modules.direction.eq(direction)].iloc[0].members.split(";")
    altered = nx.Graph()
    altered.add_nodes_from(reversed(list(graph)))
    altered.add_edges_from((b,a) for a,b in reversed(list(graph.edges())))
    a, _ = module_positions(graph, members)
    b, _ = module_positions(altered, members)
    assert all(np.array_equal(a[g], b[g]) for g in members), "Layout order dependence"
assert len(module_members) == 236

annotations = pd.read_csv(T / "module_annotations_all_tests.csv.gz")
assert len(annotations) == 2 * 42 * 6203
assert annotations.groupby(["background", "module_id"]).size().eq(6203).all()
assert annotations.P_value.between(0,1).all()
assert annotations.BH_all_modules_and_terms.between(0,1).all()
assert annotations.loc[annotations.overlap_genes.eq(0), "P_value"].eq(1).all()
assert annotations.overlap_genes.le(annotations.module_genes).all()
assert annotations.overlap_genes.le(annotations.term_genes).all()
for background, frame in annotations.groupby("background"):
    assert len(frame) == summary["annotation_tests_per_background"]
    ordered = frame.sort_values("P_value")
    assert np.diff(ordered.BH_all_modules_and_terms).min() >= -1e-12
    assert (frame.BH_all_modules_and_terms + 1e-12 >= frame.P_value).all()
    if background.startswith("expressed"):
        assert frame.background_genes.eq(14564).all()
    else:
        assert frame.loc[frame.direction.eq("DOWN"), "background_genes"].eq(1306).all()
        assert frame.loc[frame.direction.eq("UP"), "background_genes"].eq(982).all()

manifest = pd.read_csv(F / "manifest.csv")
assert len(manifest) == 6
assert len(PdfReader(F / "HT1080_LOY_network_module_atlas.pdf").pages) == 6
positions = pd.read_csv(T / "node_and_label_positions.csv")
assert len(positions) == 236 and set(positions.gene_name) == module_members
fonts, export_hashes = set(), []
for row in manifest.itertuples():
    pdf = PdfReader(F / f"{row.stem}.pdf")
    assert len(pdf.pages) == 1
    page = pdf.pages[0]
    assert abs(float(page.mediabox.width)*25.4/72-183) < .01
    assert abs(float(page.mediabox.height)*25.4/72-row.height_mm) < .01
    assert row.height_mm <= 247
    text_chunks = []
    page.extract_text(visitor_text=lambda text, cm, tm, font, size: text_chunks.append(text))
    text = "\n".join(text_chunks)
    expected = set()
    for mid in row.modules.split(";"):
        expected |= set(modules.set_index("module_id").loc[mid, "members"].split(";"))
    for gene in expected:
        assert re.search(r"(?<![A-Za-z0-9])"+re.escape(gene)+r"(?![A-Za-z0-9])", text), (row.stem,gene)
    assert len(expected) == row.genes
    for ref in page["/Resources"]["/Font"].values():
        font = ref.get_object()
        for item in font.get("/DescendantFonts", [font]):
            f = item.get_object()
            descriptor = f["/FontDescriptor"].get_object()
            assert any(k in descriptor for k in ["/FontFile", "/FontFile2", "/FontFile3"])
            assert "Arial" in str(f.get("/BaseFont")), str(f.get("/BaseFont"))
            fonts.add(str(f.get("/BaseFont")))
    svg = (F / f"{row.stem}.svg").read_text(encoding="utf-8")
    assert "<text" in svg and "Arial" in svg
    with Image.open(F / f"{row.stem}.png") as image:
        assert abs(image.width-183/25.4*600) < 1.1
        assert abs(image.height-row.height_mm/25.4*600) < 1.1
        assert all(abs(d-600)<1 for d in image.info["dpi"])
    for ext in ("pdf", "svg", "png"):
        path = F / f"{row.stem}.{ext}"
        export_hashes.append({"path":path.relative_to(ROOT).as_posix(), "sha256":digest(path)})
plot_checks = json.loads((P / "figure_export_checks.json").read_text())
assert all(r["overlapping_label_pairs"]==0 for r in plot_checks["label_checks"])
report = ROOT / "outputs/report/network_modules.html"
soup = BeautifulSoup(report.read_text(encoding="utf-8"), "html.parser")
assert len(soup.find_all("svg")) == 2
local_links = 0
for a in soup.find_all("a", href=True):
    target = a["href"]
    if target.startswith(("#", "https:", "http:", "data:")):
        continue
    path = unquote(urlparse(target).path)
    assert (report.parent / path).exists(), target
    local_links += 1
for path in [ROOT / "outputs/captions/network_modules_captions.md", report, F / "HT1080_LOY_network_module_atlas.pdf"]:
    export_hashes.append({"path":path.relative_to(ROOT).as_posix(),"sha256":digest(path)})

# Independent preservation snapshots from the completed GSEA-curve task.
preserved = 0
for row in pd.read_csv(ROOT / "provenance/GSEA_CURVES_frozen_inputs.csv").itertuples():
    assert digest(ROOT / row.path) == row.sha256, row.path
    preserved += 1
for row in pd.read_csv(ROOT / "provenance/GSEA_CURVES_prior_figure_hashes.csv").to_dict("records"):
    for ext in ("pdf", "svg", "png"):
        assert digest(ROOT / row[ext]) == row[ext+"_sha256"], row[ext]

result = {"validated_utc":datetime.now(timezone.utc).isoformat(), "modules":42,"module_genes":236,
          "mapped_candidates":2288,"annotation_tests_per_background":260526,
          "figure_pages":6,"individual_page_exports":18,"atlas_pages":6,"broken_local_links":0,
          "local_report_links_checked":local_links,"embedded_fonts":sorted(fonts),
          "prior_scientific_inputs_unchanged":preserved,"prior_figure_exports_unchanged":51,
          "known_graph_algorithm_checks":"passed","layout_insertion_order_check":"passed",
          "gene_label_overlap_pairs":0,"all_plotted_genes_found_in_PDF_text":True}
(P / "validation.json").write_text(json.dumps(result,indent=2))
pd.DataFrame(export_hashes).to_csv(P / "figure_caption_report_hashes.csv",index=False)
print(json.dumps(result,indent=2))
