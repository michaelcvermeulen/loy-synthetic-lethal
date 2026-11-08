"""Render the audited guide-discovery workflow, with editable vector text.

Run from any working directory; all defaults resolve from this file.
This visualises saved results, not a newly executed guide-search algorithm.
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
from matplotlib.patches import Rectangle, FancyArrowPatch


THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style

CHAPTER = THESIS / "figures/chapter3"
STEM = "ch03_main_guide_design"
WIDTH, HEIGHT = 183, 236  # mm; explicit final physical size
INK = "#222222"
MUTED = "#555C60"
LINE = "#B8C3C6"
TEAL = "#008B78"  # structural emphasis, not a condition/guide colour mapping
PALE = "#E9F3F0"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=CHAPTER)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    source_dir = CHAPTER / "source_data/guide_design"
    count_path = source_dir / "sequence_filter_counts.csv"
    guide_path = CHAPTER / "source_data/guide_selection_initial_five/guide_summary.csv"
    with count_path.open(encoding="utf-8-sig", newline="") as handle:
        stages = {r["stage"]: r for r in csv.DictReader(handle)}
    with guide_path.open(encoding="utf-8-sig", newline="") as handle:
        guides = {r["guide_id"]: r for r in csv.DictReader(handle)}
    order = ["Q53", "CN", "CP", "P47", "P51"]
    for name in order:
        g = guides[name]
        assert int(g["exact_Y_matches"]) == int(g["forward_Y_matches"]) + int(g["reverse_Y_matches"])
        assert int(g["exact_non_Y_matches"]) == 0
        assert len(g["sequence"]) == 20

    family = apply_style(base_size=9)
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, WIDTH), ylim=(HEIGHT, 0))
    ax.set_axis_off()
    text_artists = []

    def text(x, y, label, size=9, weight="normal", color=INK, ha="left", va="center", **kwargs):
        item = ax.text(x, y, label, fontsize=size, fontweight=weight, color=color,
                       ha=ha, va=va, linespacing=1.35, **kwargs)
        text_artists.append(item)
        return item

    def line(x1, y1, x2, y2, color=LINE, lw=.65, **kwargs):
        ax.plot([x1, x2], [y1, y2], color=color, lw=lw, **kwargs)

    def arrow(x1, y1, x2, y2, dashed=False, color=MUTED):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=7.5, linewidth=.8, color=color,
                     linestyle=(0, (2.5, 2)) if dashed else "solid",
                     shrinkA=0, shrinkB=0))

    def panel(letter, y, title):
        text(8, y, letter, 11, "semibold")
        text(14, y, title, 9.5, "semibold")

    text(8, 8, "main \u2014 Discovery of Y-chromosome deletion guides", 11, "semibold")
    text(8, 16, "Repetitive Y-specific targets for chromosome shredding and centromere targeting", 8.5, color=MUTED)

    # a. Search unit and orientation. Diagram lengths are schematic.
    panel("a", 28, "Enumerate repeated, PAM-compatible chrY sequences")
    text(10, 44, "HG002 T2T-Y", 10, "semibold")
    text(10, 51, "62,460,029 bp", 8.5, color=MUTED)
    text(10, 58, "CP086569.2", 8.5, color=MUTED)
    text(10, 66, "Count identical 23-mers", 8, color=MUTED)
    arrow(52, 51, 66, 51)
    text(115, 37, "23-nt protospacer\u2013PAM sequence", 9, "semibold", ha="center")
    for y, reverse in [(43, False), (55, True)]:
        spacer_x, pam_x = (88, 74) if reverse else (74, 140)
        ax.add_patch(Rectangle((spacer_x, y), 66, 8, facecolor=PALE, edgecolor="none"))
        ax.add_patch(Rectangle((pam_x, y), 14, 8, facecolor=TEAL, edgecolor="none"))
        text(spacer_x + 33, y + 4, "20-nt protospacer", 8.5, ha="center")
        text(pam_x + 7, y + 4, "CCN" if reverse else "NGG", 8.5, "semibold", color="white", ha="center")
        text(70.5, y + 4, "5\u2032", 8.5, ha="right")
        text(157, y + 4, "3\u2032", 8.5)
    text(115, 69, "Recognise either PAM orientation", 8.5, color=MUTED, ha="center")

    # b. Counts refer exclusively to saved, distinct 23-mers.
    panel("b", 80, "Filter recurrent Y target sequences")
    boxes = [
        (8, "recurrent", "Recurrent,\nPAM-compatible", "\u226518 stored Y occurrences"),
        (66, "Y_restricted", "No exact matches\noutside chrY", "chr1\u201322, chrX and chrM"),
        (124, "high_copy", "High repeat\nabundance", ">20 occurrences on chrY"),
    ]
    for x, key, heading, criterion in boxes:
        ax.add_patch(Rectangle((x, 89), 51, 46, facecolor="white", edgecolor=LINE, linewidth=.65))
        line(x, 89, x + 51, 89, TEAL, 1.6)
        text(x + 25.5, 97, heading, 8.8, "semibold", ha="center")
        text(x + 25.5, 111, f"{int(stages[key]['count']):,}", 20, "semibold", ha="center")
        text(x + 25.5, 121, "distinct 23-mers", 8.5, ha="center")
        text(x + 25.5, 130, criterion, 7.9, color=MUTED, ha="center")
    arrow(60, 111, 65, 111)
    arrow(118, 111, 123, 111)

    # The dashed link deliberately denotes a conceptual handoff, not a
    # recovered row-wise mapping from 7,470 target sequences to 20-nt guides.
    arrow(149.5, 136, 149.5, 142, dashed=True)
    panel("c", 146, "Assess candidate guides (20-nt spacers)")
    criteria = [
        (8, "Targeting pattern", "Exact-match abundance\nand regional distribution\nacross both orientations"),
        (66, "Near-match profile", "1\u20134 mismatches and CFD\nSummarise chrY\nnear-matches separately"),
        (124, "Predicted activity", "DeepSpCas9 predictions\nwhere available"),
    ]
    for x, heading, body in criteria:
        line(x, 153, x + 51, 153, TEAL, .9)
        text(x, 158, heading, 8.8, "semibold")
        text(x, 164, body, 8.2, va="top")
    arrow(91.5, 179, 91.5, 184, dashed=True)

    # d. Initial selection is user-confirmed; these are not algorithmic top ranks.
    panel("d", 189, "Initial selection: five complementary candidates")
    centres = [24.7, 58.1, 91.5, 124.9, 158.3]
    regions = [(8, 41.4, "Yq repeats"), (44.4, 105.2, "Centromere"), (110.2, 175, "Yp")]
    for left, right, label in regions:
        text((left + right) / 2, 198, label, 8.8, color=MUTED, ha="center")
        line(left, 202, right, 202)
    for name, x in zip(order, centres):
        text(x, 208, name, 10, "semibold", ha="center")
        text(x, 217, f"{int(guides[name]['exact_Y_matches']):,}", 15, "semibold", ha="center")
    text(91.5, 226, "Exact chrY matches; all five have zero exact non-Y matches", 8.5, ha="center")
    text(91.5, 232, "Counts sum forward- and reverse-orientation matches", 8, color=MUTED, ha="center")

    # At physical export size, catch accidental canvas overflow before writing.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    for item in text_artists:
        box = item.get_window_extent(renderer)
        if box.x0 < canvas.x0 or box.x1 > canvas.x1 or box.y0 < canvas.y0 or box.y1 > canvas.y1:
            raise RuntimeError(f"Text outside figure canvas: {item.get_text()}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    title = "main - Discovery of Y-chromosome deletion guides"
    for suffix in ("pdf", "svg", "png"):
        metadata = {"Title": title}
        fig.savefig(args.output_dir / f"{STEM}.{suffix}", dpi=args.dpi, metadata=metadata)
    # Small reviewer copy, separate from thesis exports.
    preview = args.output_dir / "previews"
    preview.mkdir(exist_ok=True)
    fig.savefig(preview / f"{STEM}.png", dpi=180)
    plt.close(fig)
    provenance = {
        "figure": STEM,
        "status": "draft for author review",
        "scope": "Results 3.4.4; sequence discovery through initial candidate selection",
        "size_mm": [WIDTH, HEIGHT], "raster_dpi": args.dpi,
        "font": family, "matplotlib": matplotlib.__version__,
        "colour_semantics": "Teal is structural emphasis only; no condition or guide colour mapping",
        "inputs": [{"path_relative_to_thesis": str(p.relative_to(THESIS)).replace('\\', '/'), "sha256": sha256(p)} for p in (count_path, guide_path)],
        "script_sha256": sha256(Path(__file__)),
        "source_audit": "research/sgRNA_design/sgRNA_design_source_audit_2026-09-16.md",
        "inference_limits": [
            "23-mer counts and 20-nt guide counts are distinct units; dashed links are conceptual.",
            "Five initial candidates are user-confirmed, not reconstructed automated top ranks.",
            "Scoring methods are shown qualitatively pending input/context reconciliation.",
            "Exact reference matches do not establish cleavage, off-target absence or chromosome loss."
        ],
        "exports": {f"{STEM}.{s}": sha256(args.output_dir / f"{STEM}.{s}") for s in ("pdf", "svg", "png")},
    }
    (source_dir / "render_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"outputs": str(args.output_dir), "size_mm": [WIDTH, HEIGHT], "font": family, "counts": {k: int(v["count"]) for k, v in stages.items()}}, indent=2))


if __name__ == "__main__":
    main()
