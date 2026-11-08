"""Render the revised guide-design figure from audited source tables.

Preserves the compact, untitled layout of the author's Illustrator edit.
Three panels show the sequence scan, guide filtering and five named candidates.
The dashed handoff denotes the shortlist with verified mismatch annotations;
the saved records do not cover all 6,901 exact-match-eligible guides.
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
STEM = "ch03_main_guide_design_edit"
WIDTH, HEIGHT = 183, 157  # mm
INK, MUTED, LINE, TEAL, PALE = "#222222", "#555C60", "#B8C3C6", "#008B78", "#E9F3F0"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=CHAPTER)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    source = CHAPTER / "source_data/guide_design_edit"
    count_path = source / "filtering_counts.csv"
    guide_path = source / "named_guide_decisions.csv"
    manifest_path = source / "run_manifest.json"
    counts = {r["stage"]: int(r["guides"]) for r in read_csv(count_path)}
    guides = {r["guide_id"]: r for r in read_csv(guide_path)}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for path in (count_path, guide_path, source / "retained_guides.csv"):
        assert sha256(path) == manifest["output_sha256"][path.name], path
    assert counts["two_mismatch_filter"] == manifest["counts"]["retained"]
    assert counts["required_mismatch_counts_available"] == counts["two_mismatch_filter"]
    assert manifest["filters"] == {"minimum_exact_y_matches": 10, "maximum_exact_non_y_matches": 0,
                                  "maximum_non_y_1_mismatch_sites": 0, "maximum_non_y_2_mismatch_sites": 0}
    order = ["Q53", "CN", "CP", "P47", "P51"]
    for name in order:
        g = guides[name]
        assert g["decision"] == "retained"
        assert len(g["guide"]) == 20 and set(g["guide"]) <= set("ACGT")
        assert int(g["exact_y_matches"]) == int(g["exact_y_forward"]) + int(g["exact_y_reverse"])
        assert all(int(g[field]) == 0 for field in ("exact_non_y_matches", "non_y_mm1", "non_y_mm2"))

    family = apply_style(base_size=9)
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, WIDTH), ylim=(HEIGHT, 0))
    ax.set_axis_off()
    texts = []

    def text(x, y, label, size=9, weight="normal", color=INK, ha="left", va="center"):
        item = ax.text(x, y, label, fontsize=size, fontweight=weight, color=color,
                       ha=ha, va=va, linespacing=1.3)
        texts.append(item)
        return item

    def line(x1, y1, x2, y2, color=LINE, lw=.65):
        ax.plot([x1, x2], [y1, y2], color=color, lw=lw)

    def arrow(x1, y1, x2, y2, dashed=False):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=7.5,
                     linewidth=.8, color=MUTED, linestyle=(0, (2.5, 2)) if dashed else "solid",
                     shrinkA=0, shrinkB=0))

    def panel(letter, y, label):
        text(8, y, letter, 11, "semibold")
        text(14, y, label, 9.3, "semibold")

    def box(x, y, height=30):
        ax.add_patch(Rectangle((x, y), 51, height, facecolor="white", edgecolor=LINE, linewidth=.65))
        line(x, y, x + 51, y, TEAL, 1.4)

    # Author's compact layout: no on-canvas overall title.
    panel("a", 6, "Scan T2T-Y for PAM-compatible guide sequences")
    text(8, 13, "HG002 T2T-Y | 62,460,029 bp", 8.5, "semibold")
    text(48, 18, "23-nt window", 8, "semibold", color=TEAL, ha="center")
    # Each small tile denotes one nucleotide. The repeated tracks are the
    # same reference segment at successive scan steps, not different loci.
    track_x, tile_width, tile_count = 17.5, 2, 35
    for step, y in enumerate((24, 30, 36)):
        window_x = track_x + (3 + step) * tile_width
        text(11, y, str(step + 1), 7.5, color=MUTED, ha="center")
        for index in range(tile_count):
            within_window = 3 + step <= index < 26 + step
            ax.add_patch(Rectangle((track_x + index * tile_width, y - 1.65), tile_width, 3.3,
                                  facecolor=PALE if within_window else "#F4F6F6",
                                  edgecolor=LINE, linewidth=.3))
        ax.add_patch(Rectangle((window_x, y - 2.3), 23 * tile_width, 4.6,
                              facecolor="none", edgecolor=TEAL, linewidth=1))
        text(15, y, "...", 7, color=MUTED, ha="right")
        text(89, y, "...", 7, color=MUTED)
        if step < 2:
            line(window_x, y + 2.3, window_x + tile_width, y + 3.7, TEAL, .7)
            line(window_x + 46, y + 2.3, window_x + 48, y + 3.7, TEAL, .7)
    arrow(25, 42, 80, 42)
    text(52.5, 46.5, "Slide window along the Y sequence", 8, color=MUTED, ha="center")

    arrow(94, 30, 102, 30)
    text(140, 17, "Retain PAM-compatible windows", 8.5, "semibold", ha="center")
    for y, reverse in [(23, False), (34, True)]:
        spacer_x, pam_x = (121, 109) if reverse else (109, 161)
        ax.add_patch(Rectangle((spacer_x, y), 52, 6.5, facecolor=PALE, edgecolor="none"))
        ax.add_patch(Rectangle((pam_x, y), 12, 6.5, facecolor=TEAL, edgecolor="none"))
        text(spacer_x + 26, y + 3.25, "20-nt protospacer", 8, ha="center")
        text(pam_x + 6, y + 3.25, "CCN" if reverse else "NGG", 8, "semibold", color="white", ha="center")
        text(107, y + 3.25, "5\u2032", 8, ha="right")
        text(175, y + 3.25, "3\u2032", 8)
    text(140, 46.5, "Both strands; reference orientation", 7.6, color=MUTED, ha="center")

    panel("b", 55, "Candidate screening and the guide shortlist (20-nt)")
    text(14, 61, "Assess Y-target abundance and genome-wide specificity", 8.2, color=MUTED)
    stages = [(8, counts["input"], "Candidate table", "distinct guides", ["Both orientations"]),
              (66, counts["exact_y_abundance_filter"], "Exact-match filters", "eligible guides", ["\u226510 exact Y matches", "0 exact non-Y matches"]),
              (124, counts["two_mismatch_filter"], "Guide shortlist", "shortlisted guides", ["No recorded non-Y sites", "with 0, 1 or 2 mismatches"])]
    for x, count, heading, unit, detail in stages:
        box(x, 66, 30)
        text(x + 25.5, 70.5, heading, 8.6, "semibold", ha="center")
        text(x + 25.5, 78, f"{count:,}", 19, "semibold", ha="center")
        text(x + 25.5, 84, unit, 8.2, ha="center")
        text(x + 25.5, 91, "\n".join(detail), 8, color=MUTED, ha="center")
    arrow(60, 80, 65, 80)
    arrow(118, 80, 123, 80, dashed=True)
    text(8, 101, "Shortlisted guides ranked by decreasing exact Y-match count", 8, color=MUTED)

    panel("c", 110, "Five guides selected for initial testing")
    centres = [23.5, 57.5, 91.5, 125.5, 159.5]
    ax.add_patch(Rectangle((6, 117), 171, 29, facecolor="white", edgecolor=LINE, linewidth=.75))
    ax.add_patch(Rectangle((6, 117), 171, 7, facecolor="#F1F5F4", edgecolor="none"))
    line(6, 117, 177, 117, LINE, .75)
    line(6, 124, 177, 124)
    for left, right, label in [(6, 41, "Yq repeats"), (41, 109, "Centromere"), (109, 177, "Yp")]:
        text((left + right) / 2, 120.5, label, 8.8, "semibold", color=MUTED, ha="center")
    for divider in (41, 109):
        line(divider, 117, divider, 146, LINE, .85)
    for name, x in zip(order, centres):
        text(x, 128, name, 10, "semibold", ha="center")
        sequence = guides[name]["guide"]
        text(x, 134, sequence, 6.5, ha="center")
        text(x, 140.5, f"{int(guides[name]['exact_y_matches']):,}", 15, "semibold", ha="center")
    text(91.5, 151, "Exact Y matches (both orientations); guide sequences shown 5\u2032 to 3\u2032, without PAM", 8, ha="center")

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for item in texts:
        bounds = item.get_window_extent(renderer)
        if not (fig.bbox.x0 <= bounds.x0 <= bounds.x1 <= fig.bbox.x1 and
                fig.bbox.y0 <= bounds.y0 <= bounds.y1 <= fig.bbox.y1):
            raise RuntimeError(f"Text outside canvas: {item.get_text()}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(args.output_dir / f"{STEM}.{suffix}", dpi=args.dpi,
                    metadata={"Title": "main - Y-chromosome guide filtering and the top-hits table"})
    preview = args.output_dir / "previews"
    preview.mkdir(exist_ok=True)
    fig.savefig(preview / f"{STEM}.png", dpi=180)
    plt.close(fig)
    provenance = {
        "figure": STEM, "size_mm": [WIDTH, HEIGHT], "dpi": args.dpi,
        "font": family, "matplotlib": matplotlib.__version__,
        "inputs": [{"path_relative_to_thesis": p.relative_to(THESIS).as_posix(), "sha256": sha256(p)}
                   for p in (count_path, guide_path, manifest_path)],
        "script_sha256": sha256(__file__),
        "scope": "Illustrated sliding-window discovery, guide filtering and five initial candidates grouped by Y region",
        "interpretation": ["Panel a tiles are schematic nucleotide positions, not reference sequence letters or a chromosome-scale map.",
                           "The three scan rows depict the same reference segment at successive window positions.",
                           "The 14,723-guide candidate table is not asserted to enumerate every PAM-compatible window in T2T-Y.",
                           "Mismatch annotations exist for the retained 215-guide subset; the dashed link is limited to available annotations.",
                           "No assertion of a new genome-wide mismatch search over all eligible guides.",
                           "Named historical guides pass the same filters; they are not the five highest ranked candidates."],
        "exports": {f"{STEM}.{s}": sha256(args.output_dir / f"{STEM}.{s}") for s in ("pdf", "svg", "png")}}
    (source / "render_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figure": str(args.output_dir / f"{STEM}.pdf"), "counts": manifest["counts"],
                      "size_mm": [WIDTH, HEIGHT]}, indent=2))


if __name__ == "__main__":
    main()
