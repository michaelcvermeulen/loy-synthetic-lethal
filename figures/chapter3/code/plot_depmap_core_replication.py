"""Compact display of saved 7 October HPC estimates; no models are refitted."""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
EXPORT = ROOT / "research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07"
FIGURES = ROOT / "figures/chapter3"
DATA = FIGURES / "source_data/depmap_core_replication"
STEM = "ch03_main_depmap_core_replication"
STYLE = ROOT / "skills/thesis-figure-style/scripts/thesis_style.py"


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    spec = importlib.util.spec_from_file_location("thesis_style", STYLE)
    style = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(style)
    font = style.apply_style(9)
    DATA.mkdir(parents=True, exist_ok=True)
    inputs = {
        "fig13d_paired_genes_across_sources.csv": EXPORT / "results/loy_reanalysis_v1/figure_source_data/fig13d_paired_genes_across_sources.csv",
        "gene_pan_results.csv": EXPORT / "derived/loy_reanalysis_v1/final/gene_pan_results.csv",
        "organoid_pan_results.csv": EXPORT / "derived/loy_reanalysis_v1/final/organoid_pan_results.csv",
        "nextgen_paired_genes.csv": EXPORT / "derived/loy_reanalysis_v1/final/nextgen_paired_genes.csv",
    }
    records = read_csv(inputs["fig13d_paired_genes_across_sources.csv"])
    genes = ["DDX3X", "EIF1AX", "RPS4X", "ZFX"]
    sources = ["Cell lines, CRISPR", "Sanger organoids, CRISPR", "NextGen, CRISPR"]
    tables = {
        sources[0]: {r["symbol"]: r for r in read_csv(inputs["gene_pan_results.csv"]) if r["symbol"] in genes},
        sources[1]: {r["symbol"]: r for r in read_csv(inputs["organoid_pan_results.csv"]) if r["symbol"] in genes},
        sources[2]: {r["symbol"]: r for r in read_csv(inputs["nextgen_paired_genes.csv"])},
    }
    plotted = []
    for gene in genes:
        for source in sources:
            matches = [r for r in records if r["symbol"] == gene and r["source"] == source]
            if not matches:
                assert gene == "EIF1AX" and source == sources[1]
                continue
            assert len(matches) == 1
            row = dict(matches[0])
            row["display_fdr"] = tables[source][gene]["permutation_fdr" if source == sources[2] else "fdr"]
            row["fdr_family"] = {sources[0]: "17,205 genes", sources[1]: "16,547 genes", sources[2]: "4 prespecified genes; permutation P"}[source]
            plotted.append(row)
    assert len(plotted) == 11
    lookup = {(r["symbol"], r["source"]): r for r in plotted}
    colours = {sources[0]: "#444444", sources[1]: "#009E87", sources[2]: "#8570B2"}
    labels = ["DepMap cell lines", "Sanger organoids", "DepMap NextGen"]
    fig, axes = plt.subplots(2, 2, figsize=(7.45, 4.15), sharex=True)
    for index, (ax, gene) in enumerate(zip(axes.flat, genes)):
        ax.axvline(0, color="#929292", lw=.7, ls=(0, (3, 3)), zorder=0)
        ax.text(-.02, 1.1, "abcd"[index], transform=ax.transAxes, weight="semibold", fontsize=11)
        ax.text(.09, 1.1, gene, transform=ax.transAxes, fontstyle="italic", fontsize=10)
        ax.text(.92, 1.1, "FDR", transform=ax.transAxes, ha="right", fontsize=8)
        for y, source in zip([2, 1, 0], sources):
            row = lookup.get((gene, source))
            if row is None:
                ax.text(-1.3, y, "Not in library", va="center", fontsize=8, color="#777777")
                ax.text(.78, y, "—", va="center", ha="center", fontsize=8)
                continue
            effect, low, high, fdr = [float(row[k]) for k in ("smd", "smd_low", "smd_high", "display_fdr")]
            assert low <= effect <= high
            colour = colours[source]
            ax.errorbar(effect, y, xerr=[[effect-low], [high-effect]], fmt="o", markersize=5,
                        color=colour, markerfacecolor=colour if fdr < .05 else "white",
                        markeredgewidth=1, elinewidth=1, capsize=2.2, zorder=3)
            qtext = f"{fdr:.1e}" if fdr < .001 else f"{fdr:.3f}"
            ax.text(.78, y, qtext, va="center", ha="center", fontsize=8)
        ax.set_yticks([2, 1, 0], labels)
        ax.tick_params(axis="y", length=0, pad=5)
        ax.set_ylim(-.55, 2.65)
        ax.set_xlim(-2.3, 1.08)
        ax.set_xticks([-2, -1, 0])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_bounds(-2.3, .18)
    fig.supxlabel("Standardised dependency difference (Y-negative − Y-positive)", y=.16, fontsize=9)
    fig.text(.5, .087, "Negative values indicate stronger dependency in Y-negative models.", ha="center", fontsize=8)
    fig.text(.5, .025, "FDR: genome-wide for cell lines and Sanger; four-gene permutation family for NextGen.", ha="center", fontsize=7.7)
    fig.subplots_adjust(left=.16, right=.985, top=.90, bottom=.265, wspace=.76, hspace=.74)
    for extension in ("pdf", "svg", "png"):
        fig.savefig(FIGURES / f"{STEM}.{extension}", dpi=300)
    plt.close(fig)
    with (DATA / "plotted_estimates.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(plotted[0]))
        writer.writeheader(); writer.writerows(plotted)
    shutil.copy2(inputs["fig13d_paired_genes_across_sources.csv"], DATA / "fig13d_paired_genes_across_sources.csv")
    manifest = {name: {"source": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for name, path in inputs.items()}
    (DATA / "source_manifest.json").write_text(json.dumps({"font": font, "inputs": manifest}, indent=2), encoding="utf-8")
    print(json.dumps({"figure": STEM, "estimates": len(plotted), "font": font}))


if __name__ == "__main__":
    main()
