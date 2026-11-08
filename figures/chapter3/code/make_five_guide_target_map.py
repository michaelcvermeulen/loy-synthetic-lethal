"""Render the five-guide target map with verified repeat intervals.

Uses the portable thesis Open Sans theme. Coordinates and repeat annotations
are frozen beside the figure; no network or genome search is required.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import numpy as np

THESIS = Path(__file__).resolve().parents[3]
CHAPTER = THESIS / "figures/chapter3"
SOURCE = CHAPTER / "source_data/guide_target_map_five"
sys.path.insert(0, str(THESIS / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style

STEM = "ch03_main_guide_target_map_stacked"
WIDTH, HEIGHT = 183, 214
LEFT, PLOT_WIDTH = 37, 138
LENGTH_BP = 62460029
ORDER = ["Q53", "CN", "CP", "P47", "P51"]
FAMILY = {"Q53": "DYZ1_HSat3", "CN": "DYZ3", "CP": "DYZ3", "P47": "TSPY", "P51": "TSPY"}
REPEAT_LABEL = {"DYZ1_HSat3": "DYZ1 / HSat3", "DYZ3": "DYZ3", "TSPY": "TSPY composite"}
GUIDE_LABEL = {"DYZ1_HSat3": "DYZ1 / HSat3 (Yq12)", "DYZ3": "DYZ3 alpha-satellite", "TSPY": "TSPY composite (Yp)"}
EXPECTED = {"Q53": 532, "CN": 53, "CP": 54, "P47": 47, "P51": 51}
INK, MUTED, LINE = "#222222", "#555C60", "#B8C3C6"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(name):
    with (SOURCE / name).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=CHAPTER)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    source_manifest = json.loads((SOURCE / "source_provenance.json").read_text(encoding="utf-8"))
    for filename, expected_hash in source_manifest["files"].items():
        assert digest(SOURCE / filename) == expected_hash, filename
    summary = {r["guide_id"]: r for r in read_csv("guide_summary.csv")}
    sites = read_csv("exact_Y_repeat_sites.csv")
    regions = read_csv("repeat_regions.csv")
    cytobands = read_csv("T2T_chrY_cytobands.csv")
    pars = read_csv("PAR_regions.csv")
    colours = json.loads((SOURCE / "repeat_colours.json").read_text(encoding="utf-8"))["palette"]
    intervals = {name: [] for name in REPEAT_LABEL}
    for row in regions:
        start, end = int(row["start_0based"]), int(row["end_0based"])
        assert 0 <= start < end <= LENGTH_BP
        intervals[row["repeat_class"]].append((start, end))
    assert {key: len(value) for key, value in intervals.items()} == {"DYZ1_HSat3": 43, "DYZ3": 1, "TSPY": 16}
    cen = np.array(intervals["DYZ3"][0]) / 1e6
    q12 = min(int(r["start"]) for r in cytobands if r["band"] == "q12") / 1e6
    par2_start = min(int(r["start_0based"]) for r in pars if r["region"] == "PAR2") / 1e6
    positions, jitter, overlap_counts = {}, {}, {}
    for index, guide in enumerate(ORDER):
        subset = sorted((r for r in sites if r["guide_id"] == guide), key=lambda r: int(r["start_1based"]))
        assert len(subset) == EXPECTED[guide] == int(summary[guide]["exact_Y_matches"])
        starts = [int(r["start_1based"]) for r in subset]
        assert len(starts) == len(set(starts))
        assert min(starts) == int(summary[guide]["min_position_bp"])
        assert max(starts) == int(summary[guide]["max_position_bp"])
        passed = 0
        for row in subset:
            assert row["guide"] == summary[guide]["sequence"]
            start, end = int(row["start_1based"]) - 1, int(row["end_1based"])
            assert end - start == 20
            assert row["canonical_NGG"] == "True" and row["pam"].endswith("GG")
            passed += any(start < b and end > a for a, b in intervals[FAMILY[guide]])
        assert passed == EXPECTED[guide]
        overlap_counts[guide] = passed
        # One-based inclusive intervals become zero-based half-open intervals;
        # the midpoint is only a plotting coordinate, not an inferred cut site.
        positions[guide] = np.array([(int(r["start_1based"]) - 1 + int(r["end_1based"])) / 2e6 for r in subset])
        jitter[guide] = np.random.default_rng(20260916 + index).uniform(.16, .84, len(subset))
    assert sum(EXPECTED.values()) == len(sites) == 737

    family = apply_style(base_size=8.5)
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4))
    guide_texts = []

    def axis(x, y, width, height):
        return fig.add_axes([x / WIDTH, 1 - (y + height) / HEIGHT, width / WIDTH, height / HEIGHT])

    def text(x, y, label, size=8.5, weight="normal", color=INK, ha="left"):
        return fig.text(x / WIDTH, 1 - y / HEIGHT, label, ha=ha, va="center", fontsize=size,
                        fontweight=weight, color=color)

    def panel(letter, x, y, label):
        text(x, y, letter, 11, "semibold")
        text(x + 6, y, label, 9.3, "semibold")

    def repeat_bars(ax, key, low=0, height=1, alpha=1):
        for start, end in intervals[key]:
            ax.add_patch(Rectangle((start / 1e6, low), (end - start) / 1e6, height,
                                   facecolor=colours[key], alpha=alpha, edgecolor="none"))

    def decorate_track(ax, guide, xlim, ticks, show_ticklabels=False, genome_context=False):
        ax.set(xlim=xlim, ylim=(0, 1), yticks=[])
        if genome_context:
            ax.axvspan(q12, par2_start, color="#000000", alpha=.028, lw=0, zorder=0)
            for r in pars:
                ax.axvspan(int(r["start_0based"]) / 1e6, int(r["end_0based"]) / 1e6,
                           color="#DCE9EC", alpha=.6, lw=0, zorder=0)
        repeat_bars(ax, FAMILY[guide], alpha=.13)
        if genome_context:
            for edge in cen:
                ax.axvline(edge, color="#7E888D", lw=.4, alpha=.7, zorder=1)
        ax.scatter(positions[guide], jitter[guide], s=4.3, c=colours[FAMILY[guide]],
                   alpha=.88, linewidths=0, zorder=3)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set(color=LINE, linewidth=.6)
        ax.set_xticks(ticks)
        ax.tick_params(axis="x", labelbottom=show_ticklabels, length=2.2, pad=2, labelsize=7.5)
        ax.tick_params(axis="y", length=0)

    text(8, 6, "a", 11, "semibold")
    ax = axis(LEFT, 14, PLOT_WIDTH, 10)
    stain_colours = {"gneg": "white", "gpos50": "#979EA4", "gvar": "#D8DDE1"}
    for row in cytobands:
        start, end = int(row["start"]) / 1e6, int(row["end"]) / 1e6
        if row["stain"] == "acen":
            vertices = [(start, .25), (start, .72), (end, .485)] if row["band"].startswith("p") else [(start, .485), (end, .25), (end, .72)]
            ax.add_patch(Polygon(vertices, facecolor=colours["DYZ3"], edgecolor=MUTED, linewidth=.4))
        else:
            ax.add_patch(Rectangle((start, .25), end - start, .47, facecolor=stain_colours[row["stain"]], edgecolor=MUTED, linewidth=.4))
    for row in pars:
        start, end = int(row["start_0based"]) / 1e6, int(row["end_0based"])/1e6
        ax.add_patch(Rectangle((start, .25), end - start, .47, facecolor="#DCE9EC", edgecolor=MUTED, linewidth=.4))
    ax.text(0, 1.0, "PAR1", fontsize=7, color=MUTED, ha="left")
    ax.text(6, 1.0, "Yp", fontsize=8, ha="center")
    ax.text(LENGTH_BP / 1e6, 1.0, "PAR2", fontsize=7, color=MUTED, ha="right")
    ax.text(37, 1.0, "Yq", fontsize=8, ha="center")
    ax.annotate("Centromere", xy=(np.mean(cen), .78), xytext=(15, 1.0), fontsize=7, color=MUTED,
                arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=.5), annotation_clip=False)
    ax.text((q12 + par2_start) / 2, .485, "Yq12 heterochromatin", fontsize=7.5, color=MUTED, ha="center", va="center")
    ax.set(xlim=(0, LENGTH_BP / 1e6), ylim=(0, 1.2))
    ax.set_axis_off()

    for row_index, key in enumerate(REPEAT_LABEL):
        y = 27 + row_index * 5
        text(LEFT - 2, y + 1.5, REPEAT_LABEL[key], 7.7, color=colours[key], ha="right")
        ax = axis(LEFT, y, PLOT_WIDTH, 3)
        repeat_bars(ax, key)
        ax.set(xlim=(0, LENGTH_BP / 1e6), ylim=(0, 1))
        ax.set_axis_off()

    for row_index, guide in enumerate(ORDER):
        y = 47 + row_index * 18
        color = colours[FAMILY[guide]]
        text(LEFT - 2, y, guide, 10, "semibold", color=color, ha="right")
        guide_texts.append(text(LEFT, y, summary[guide]["sequence"], 8.2))
        label = GUIDE_LABEL[FAMILY[guide]].replace("alpha", "\u03b1")
        text(84, y, label, 7.8, color=color)
        text(175, y, f"{EXPECTED[guide]} exact matches", 8, color=MUTED, ha="right")
        ax = axis(LEFT, y + 3, PLOT_WIDTH, 11)
        decorate_track(ax, guide, (0, LENGTH_BP / 1e6), np.arange(0, 61, 10), row_index == len(ORDER) - 1, True)
    text(LEFT + PLOT_WIDTH / 2, 144, "Position on HG002 T2T-Y (Mb)", 8.5, ha="center")

    zooms = [dict(letter="b", label_x=8, title="Centromere: DYZ3", x=21, width=65, key="DYZ3",
                  xlim=(10.54, 10.91), ticks=[10.6, 10.7, 10.8, 10.9], guides=["CN", "CP"]),
             dict(letter="c", label_x=97, title="Yp: TSPY repeats", x=111, width=64, key="TSPY",
                  xlim=(5.8, 10.4), ticks=[6, 7, 8, 9, 10], guides=["P47", "P51"])]
    for zoom in zooms:
        panel(zoom["letter"], zoom["label_x"], 154, zoom["title"])
        ax = axis(zoom["x"], 162, zoom["width"], 3)
        repeat_bars(ax, zoom["key"])
        ax.set(xlim=zoom["xlim"], ylim=(0, 1))
        ax.set_axis_off()
        for row_index, guide in enumerate(zoom["guides"]):
            assert ((positions[guide] >= zoom["xlim"][0]) & (positions[guide] <= zoom["xlim"][1])).all()
            y = 169 + 13 * row_index
            text(zoom["x"] - 2, y + 4, guide, 8.5, "semibold", color=colours[zoom["key"]], ha="right")
            ax = axis(zoom["x"], y, zoom["width"], 8)
            decorate_track(ax, guide, zoom["xlim"], zoom["ticks"], row_index == 1)
        text(zoom["x"] + zoom["width"] / 2, 199, "Position (Mb)", 8, ha="center")

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for item in guide_texts:
        assert "\n" not in item.get_text()
        right_mm = item.get_window_extent(renderer).x1 / fig.bbox.width * WIDTH
        assert right_mm < 82, (item.get_text(), right_mm)
    for item in fig.texts:
        bounds = item.get_window_extent(renderer)
        assert fig.bbox.contains(bounds.x0, bounds.y0) and fig.bbox.contains(bounds.x1, bounds.y1), item.get_text()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    title = "main - Five Y-targeting guides and their repeat regions"
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(args.output_dir / f"{STEM}.{suffix}", dpi=args.dpi, metadata={"Title": title})
    preview = args.output_dir / "previews"
    preview.mkdir(exist_ok=True)
    fig.savefig(preview / f"{STEM}.png", dpi=180)
    plt.close(fig)
    render = {"figure": STEM, "size_mm": [WIDTH, HEIGHT], "dpi": args.dpi, "font": family,
              "matplotlib": matplotlib.__version__, "guide_order": ORDER, "exact_site_counts": EXPECTED,
              "repeat_overlap_counts": overlap_counts, "repeat_colours": colours,
              "source_manifest_sha256": digest(SOURCE / "source_provenance.json"), "script_sha256": digest(__file__),
              "coordinate_display": "Midpoints of exact 20-nt match intervals; 1 Mb = 1,000,000 bp. No inferred cleavage positions.",
              "repeat_display": "All chrY annotations for the three repeat classes, including intervals with no target for these guides.",
              "jitter": "Deterministic vertical display offsets only, seed 20260916 plus guide index; identical between full map and zooms.",
              "exports": {f"{STEM}.{s}": digest(args.output_dir / f"{STEM}.{s}") for s in ("pdf", "svg", "png")}}
    (SOURCE / "render_provenance.json").write_text(json.dumps(render, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output_dir / f"{STEM}.pdf"), "counts": EXPECTED,
                      "verified_repeat_overlaps": overlap_counts, "size_mm": [WIDTH, HEIGHT]}, indent=2))


if __name__ == "__main__":
    main()
