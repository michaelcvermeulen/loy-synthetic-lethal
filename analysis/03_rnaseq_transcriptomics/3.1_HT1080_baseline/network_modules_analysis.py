"""Fresh STRING/MCODE module analysis of the frozen HT1080 expression results.

Local Python implementation of the published MCODE algorithm for simple,
undirected graphs, with lexical tie/neighbor ordering for reproducibility.
Reference: Bader & Hogue (2003), doi:10.1186/1471-2105-4-2; algorithm behavior
checked against BaderLab/MCODE commit a5657b9281b709ecd09f9747ef110f7f156c6cdb.
This is not a Metascape submission or a Cytoscape-app execution.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
STRING = REPO / "analysis/_shared/STRING"
TABLES = ROOT / "outputs/tables/network_modules"
PROV = ROOT / "provenance/network_modules"
CONFIG = ROOT / "config/network_modules.json"


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def write_csv(df, name):
    df.to_csv(TABLES / name, index=False)


def mcode(graph, params):
    """Density of highest neighborhood k-core times k, then seed expansion.

    Closed neighborhoods include the seed. Nodes with degree < 2 score zero.
    Seen status persists after unsuccessful clusters and haircut, as in the
    reference implementation. Fluff and self-loops are deliberately unsupported.
    """
    assert not params["fluff"]
    assert not graph.is_directed() and nx.number_of_selfloops(graph) == 0
    scores = {}
    for node in sorted(graph):
        if graph.degree(node) < params["degree_cutoff"]:
            scores[node] = 0.0
            continue
        neighborhood = graph.subgraph([node, *graph.neighbors(node)])
        cores = nx.core_number(neighborhood)
        k = max(cores.values())
        highest = neighborhood.subgraph([n for n in cores if cores[n] >= k])
        scores[node] = k * nx.density(highest)

    seen, clusters = set(), []
    for seed in sorted(graph, key=lambda n: (-scores[n], n)):
        if seed in seen:
            continue
        members = set()
        threshold = scores[seed] * (1 - params["node_score_cutoff"])

        def expand(node, depth):
            if node in seen:
                return
            seen.add(node)
            if depth > params["max_depth"]:
                return
            # Degree-zero/one nodes have no expansion neighborhood in MCODE.
            if graph.degree(node) < 2:
                return
            for neighbor in sorted(graph.neighbors(node)):
                if neighbor not in seen and scores[neighbor] >= threshold:
                    members.add(neighbor)
                    expand(neighbor, depth + 1)

        expand(seed, 1)
        if not members:
            continue
        members.add(seed)
        sub = graph.subgraph(members)
        if not nx.k_core(sub, k=params["k_core"]):
            continue
        if params["haircut"]:
            members = set(nx.k_core(sub, k=2))
        sub = graph.subgraph(members)
        clusters.append({"seed": seed, "members": sorted(members),
                         "nodes": len(sub), "edges": sub.number_of_edges(),
                         "density": nx.density(sub),
                         "mcode_score": len(sub) * nx.density(sub)})
    clusters.sort(key=lambda c: (-c["mcode_score"], -c["nodes"], c["seed"]))
    return clusters, scores


def check_algorithm(params):
    # Exact expected outcomes distinguish core density from full-neighborhood
    # density, preserve dense clusters joined by a bridge, and prune a tail.
    clique = nx.complete_graph(["A", "B", "C", "D"])
    clique.add_edges_from([("D", "E"), ("E", "F")])
    clusters, scores = mcode(clique, params)
    assert [c["members"] for c in clusters] == [["A", "B", "C", "D"]]
    assert math.isclose(scores["D"], 3.0)
    assert math.isclose(clusters[0]["mcode_score"], 4.0)
    graph = nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(3))
    graph.add_edge(3, 4)
    clusters, _ = mcode(graph, params)
    assert {frozenset(c["members"]) for c in clusters} == {
        frozenset(range(4)), frozenset(range(4, 7))}
    assert not mcode(nx.path_graph(8), params)[0]
    assert not mcode(nx.star_graph(8), params)[0]
    assert len(mcode(nx.cycle_graph(5), params)[0]) == 1
    reversed_graph = nx.Graph()
    reversed_graph.add_edges_from(reversed(list(clique.edges())))
    assert mcode(clique, params)[0] == mcode(reversed_graph, params)[0]


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(adjusted, 1)
    return out


def load_inputs():
    config = json.loads(CONFIG.read_text())
    files = [CONFIG, Path(__file__),
             ROOT / "outputs/tables/HT1080_LOY_vs_sgCTRL_primary_effects.csv",
             ROOT / "outputs/tables/symbol_mapping_representatives.csv",
             ROOT / "resources/msigdb_all_membership.csv.gz",
             ROOT / "resources/msigdb_all_set_metadata.csv",
             STRING / "checksums.csv"]
    checksums = pd.read_csv(STRING / "checksums.csv")
    for name in ["9606.protein.info.v12.0.txt.gz", "9606.protein.links.v12.0.txt.gz",
                 "9606.protein.physical.links.v12.0.txt.gz"]:
        path = STRING / name
        assert digest(path, "md5") == checksums.set_index("file").loc[name, "md5"]
        files.append(path)
    hashes = [{"file": p.relative_to(REPO).as_posix(), "sha256": digest(p)} for p in files]
    scientific_hashes = [x for x in hashes if not x["file"].endswith((".py", ".json"))]
    frozen = PROV / "frozen_scientific_input_hashes.json"
    if frozen.exists():
        assert json.loads(frozen.read_text()) == scientific_hashes, "Frozen scientific input changed"
    else:
        frozen.write_text(json.dumps(scientific_hashes, indent=2))
    (PROV / "run_input_hashes.json").write_text(json.dumps(hashes, indent=2))
    p = pd.read_csv(files[2])
    reps = pd.read_csv(files[3])
    data = p[p.gene_id.isin(reps.gene_id)].copy().sort_values("gene_name")
    assert len(data) == 25461 and data.gene_name.is_unique
    required = (data.mean_difference.abs().ge(.5) & data.period_direction_agrees & data.LOO_all_same_direction)
    assert np.array_equal(required, data.descriptive_candidate)
    info = pd.read_csv(STRING / "9606.protein.info.v12.0.txt.gz", sep="\t")
    assert info.preferred_name.is_unique and info["#string_protein_id"].is_unique
    data["string_id"] = data.gene_name.map(info.set_index("preferred_name")["#string_protein_id"])
    data["direction"] = np.where(data.mean_difference < 0, "DOWN", "UP")
    data["mapping_status"] = np.where(data.string_id.notna(), "exact_preferred_symbol", "unmapped")
    write_csv(data, "all_symbol_mapping_and_effects.csv")
    return config, data, info, len(p), int(p.descriptive_candidate.sum())


def load_edges(info, candidate_symbols, filename, score_cutoff):
    lookup = info.set_index("#string_protein_id").preferred_name.to_dict()
    pieces = []
    for frame in pd.read_csv(STRING / filename, sep=" ", chunksize=500000):
        frame = frame[frame.combined_score >= score_cutoff].copy()
        frame["a"] = frame.protein1.map(lookup)
        frame["b"] = frame.protein2.map(lookup)
        frame = frame[frame.a.isin(candidate_symbols) & frame.b.isin(candidate_symbols) & frame.a.ne(frame.b)]
        if len(frame):
            pairs = [tuple(sorted(pair)) for pair in zip(frame.a, frame.b)]
            pieces.append(pd.DataFrame({"gene_a": [x[0] for x in pairs],
                                        "gene_b": [x[1] for x in pairs],
                                        "score": frame.combined_score.to_numpy()}))
    edges = pd.concat(pieces, ignore_index=True)
    return edges.groupby(["gene_a", "gene_b"], as_index=False).score.max()


def annotate_modules(modules, data, config):
    meta = pd.read_csv(ROOT / "resources/msigdb_all_set_metadata.csv").fillna("")
    meta = meta[meta.gs_subcollection.isin(config["annotation_collections"])]
    wanted = set(meta.gs_name)
    universe = set(data.loc[data.string_id.notna(), "gene_name"])
    terms = defaultdict(set)
    with gzip.open(ROOT / "resources/msigdb_all_membership.csv.gz", "rt", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["pathway"] in wanted and row["symbol"] in universe:
                terms[row["pathway"]].add(row["symbol"])
    terms = {k: v for k, v in sorted(terms.items())
             if config["annotation_min_size"] <= len(v) <= config["annotation_max_size"]}
    term_names = list(terms)
    term_arrays = list(terms.values())
    sizes = np.array([len(v) for v in term_arrays])
    direction_sets = {d: set(data.loc[data.descriptive_candidate & data.string_id.notna() & data.direction.eq(d), "gene_name"])
                      for d in ("DOWN", "UP")}
    families = []
    for background in ("expressed_STRING_mapped", "direction_candidate_STRING_mapped"):
        frames = []
        for module in modules:
            genes = set(module["members"])
            bg = universe if background.startswith("expressed") else direction_sets[module["direction"]]
            bg_sizes = sizes if background.startswith("expressed") else np.array([len(t & bg) for t in term_arrays])
            overlap = np.array([len(genes & t) for t in term_arrays])
            pval = hypergeom.sf(overlap-1, len(bg), bg_sizes, len(genes))
            assert np.isfinite(pval).all()
            frames.append(pd.DataFrame({"module_id": module["module_id"], "direction": module["direction"],
                "background": background, "background_genes": len(bg), "module_genes": len(genes),
                "pathway": term_names, "term_genes": bg_sizes, "overlap_genes": overlap,
                "P_value": pval, "overlap_symbols": [";".join(sorted(genes & t)) for t in term_arrays]}))
        family = pd.concat(frames, ignore_index=True)
        family["BH_all_modules_and_terms"] = bh(family.P_value)
        families.append(family)
    annotation = pd.concat(families, ignore_index=True)
    annotation.to_csv(TABLES / "module_annotations_all_tests.csv.gz", index=False, compression="gzip")
    display = annotation[(annotation.background == "expressed_STRING_mapped") & (annotation.overlap_genes >= 2)]
    top = display.sort_values(["module_id", "BH_all_modules_and_terms", "P_value", "pathway"]).groupby("module_id").head(1)
    sensitivity = annotation[annotation.background == "direction_candidate_STRING_mapped"].set_index(["module_id", "pathway"])
    top = top.copy()
    top["candidate_background_BH"] = [sensitivity.loc[(r.module_id, r.pathway), "BH_all_modules_and_terms"] for r in top.itertuples()]
    write_csv(top, "module_top_annotations.csv")
    write_csv(meta[meta.gs_name.isin(term_names)], "annotation_term_metadata.csv")
    return annotation, top, len(terms)


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    PROV.mkdir(parents=True, exist_ok=True)
    config, data, info, gene_count, gene_candidates = load_inputs()
    check_algorithm(config["mcode"])
    candidates = data[data.descriptive_candidate & data.string_id.notna()].copy()
    print(f"Selected {data.descriptive_candidate.sum()} unique-symbol candidates; {len(candidates)} exact STRING mappings", flush=True)
    edges = load_edges(info, set(candidates.gene_name), "9606.protein.links.v12.0.txt.gz", config["association_score"])
    physical = load_edges(info, set(candidates.gene_name), "9606.protein.physical.links.v12.0.txt.gz", config["physical_score"])
    edges = edges.rename(columns={"score": "association_score"}).merge(
        physical.rename(columns={"score": "physical_score"}), on=["gene_a", "gene_b"], how="left")
    edges["physical_supported"] = edges.physical_score.notna()
    direction_map = candidates.set_index("gene_name").direction.to_dict()
    edges["direction_a"] = edges.gene_a.map(direction_map)
    edges["direction_b"] = edges.gene_b.map(direction_map)
    write_csv(edges, "all_candidate_edges_including_cross_direction.csv")
    modules, node_rows, coverage, physical_modules, omission_rows = [], [], [], [], []
    graph_dir = ROOT / "outputs/models/network_modules"
    graph_dir.mkdir(parents=True, exist_ok=True)
    for direction in ("DOWN", "UP"):
        selected = candidates[candidates.direction == direction]
        symbols = set(selected.gene_name)
        ed = edges[edges.direction_a.eq(direction) & edges.direction_b.eq(direction)].copy()
        graph = nx.Graph()
        graph.add_nodes_from(sorted(symbols))
        for row in ed.itertuples():
            graph.add_edge(row.gene_a, row.gene_b, association_score=int(row.association_score),
                           physical_supported=bool(row.physical_supported),
                           physical_score=int(row.physical_score) if pd.notna(row.physical_score) else 0)
        cluster_list, scores = mcode(graph, config["mcode"])
        assigned = {}
        for i, cluster in enumerate(cluster_list, 1):
            mid = f"{direction}_{i:02d}"
            cluster.update(module_id=mid, direction=direction, rank=i,
                           shown_main=i <= config["main_modules_per_direction"])
            sub = graph.subgraph(cluster["members"])
            cluster["physical_edges"] = sum(attr["physical_supported"] for _, _, attr in sub.edges(data=True))
            cluster["physical_edge_fraction"] = cluster["physical_edges"] / cluster["edges"]
            assert len(nx.k_core(sub, k=2)) == len(sub)
            for gene in cluster["members"]:
                assert gene not in assigned
                assigned[gene] = mid
            modules.append(cluster)
            effect = selected.set_index("gene_name").loc[cluster["members"]]
            for key in ["mean_difference", "AUG_2025", "MAY_2026", "mean_CPM", "matches_agreeing_with_mean"]:
                cluster[f"mean_{key}"] = float(effect[key].mean())
            for scenario in ["prior025_mean_difference", "RLE_mean_difference", *[x for x in effect if x.startswith("omit_")]]:
                omission_rows.append({"module_id": mid, "scenario": scenario, "genes": len(effect),
                                      "mean_difference": effect[scenario].mean(),
                                      "fraction_original_direction": (effect[scenario] * (-1 if direction == "DOWN" else 1) > 0).mean()})
        for gene in sorted(graph):
            row = {"gene_name": gene, "direction": direction, "module_id": assigned.get(gene, ""),
                   "node_mcode_score": scores[gene], "degree": graph.degree(gene),
                   "status": "module_member" if gene in assigned else ("connected_unclustered" if graph.degree(gene) else "no_within_direction_edge")}
            node_rows.append(row)
            graph.nodes[gene].update(row)
        nx.write_graphml(graph, graph_dir / f"STRING_association_{direction}.graphml")
        strict = nx.Graph()
        strict.add_nodes_from(sorted(symbols))
        strict_ed = physical[physical.gene_a.isin(symbols) & physical.gene_b.isin(symbols)]
        strict.add_edges_from(zip(strict_ed.gene_a, strict_ed.gene_b))
        strict_clusters, _ = mcode(strict, config["mcode"])
        for i, cluster in enumerate(strict_clusters, 1):
            physical_modules.append({**cluster, "direction": direction, "physical_module_id": f"{direction}_PHYSICAL_{i:02d}"})
        nx.write_graphml(strict, graph_dir / f"STRING_physical_{direction}.graphml")
        coverage.append({"direction": direction, "candidate_unique_symbols": int((data.descriptive_candidate & data.direction.eq(direction)).sum()),
                         "STRING_mapped_candidates": len(symbols), "connected_genes": sum(graph.degree(n)>0 for n in graph),
                         "association_edges": len(ed), "modules": len(cluster_list), "module_genes": len(assigned),
                         "strict_physical_connected_genes": sum(strict.degree(n)>0 for n in strict),
                         "strict_physical_edges": strict.number_of_edges(), "strict_physical_modules": len(strict_clusters)})
        write_csv(ed, f"{direction}_association_edges.csv")
        print(coverage[-1], flush=True)
    node_df = pd.DataFrame(node_rows).merge(data, on=["gene_name", "direction"], validate="one_to_one")
    write_csv(node_df, "candidate_node_evidence.csv")
    write_csv(pd.DataFrame(coverage), "coverage.csv")
    write_csv(pd.DataFrame(omission_rows), "module_expression_sensitivity.csv")
    physical_df = pd.DataFrame(physical_modules)
    physical_df["members"] = physical_df.members.map(lambda x: ";".join(x))
    write_csv(physical_df, "strict_physical_modules.csv")
    annotation, top, term_count = annotate_modules(modules, data, config)
    module_df = pd.DataFrame(modules)
    module_df["members"] = module_df.members.map(lambda x: ";".join(x))
    module_df = module_df.merge(top[["module_id", "pathway", "overlap_genes", "BH_all_modules_and_terms", "candidate_background_BH"]], on="module_id", how="left")
    write_csv(module_df, "modules.csv")
    # Quantify recurrence of the pictured legacy members; this is an overlay,
    # not a reuse of the legacy DEG selection or of its module annotation P.
    legacy = REPO / "analysis/HT1080/reannotate_v105/metascape/ppi_network_combinedrun_DOWN_nodes.csv"
    if legacy.exists():
        old = pd.read_csv(legacy)
        merged = old.merge(data, left_on="gene", right_on="gene_name", how="left", suffixes=("_legacy", "_current"))
        write_csv(merged, "legacy_DOWN_members_current_evidence.csv")
    summary = {"completed_utc": datetime.now(timezone.utc).isoformat(), "eligible_gene_loci": gene_count,
               "candidate_gene_loci": gene_candidates, "unique_symbols": len(data),
               "candidate_unique_symbols": int(data.descriptive_candidate.sum()),
               "STRING_mapped_background": int(data.string_id.notna().sum()), "annotation_terms": term_count,
               "annotation_tests_per_background": int(len(annotation)/2), "coverage": coverage,
               "algorithm": "Local Python MCODE implementation; lexical ordering; no Cytoscape or Metascape execution",
               "checks": "Known-graph core/haircut tests; insertion-order invariance; frozen input hashes; source MD5; exact mappings; selection flags; module 2-cores; unique membership; finite annotation tests"}
    (PROV / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print(module_df[["module_id", "nodes", "mcode_score", "pathway"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
