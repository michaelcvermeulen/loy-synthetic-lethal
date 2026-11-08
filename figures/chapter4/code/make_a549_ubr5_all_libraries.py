"""Display UBR5 in all 11 A549 libraries, including the excluded culture."""
from pathlib import Path
import sys
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent / ".figure_python_deps"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import StrMethodFormatter
import pymupdf

sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)

FIGS = ROOT / "figures/chapter4"
DATA = FIGS / "source_data/a549_ubr5_all_libraries_20261020"
DATA.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT.parent / "crispr_screen/a549_loy_screen1/CRISPR_screen/Screen1/Mageck/A549_full.count.txt"
MASK = FIGS / "source_data/a549_p3_rescue_20261020/top250_union_excluded_guides.csv"
EXCLUDED = "CTRL_P3_R1"
SAMPLES = ["CTRL_D0_1", "CTRL_D0_2", "LOY_D0_1", "LOY_D0_2",
           "CTRL_P3_R2", "LOY_P3_R1", "LOY_P3_R2", EXCLUDED,
           "CTRL_P4_R2", "LOY_P4_R1", "LOY_P4_R2"]
LABELS = ["CTRL\nT1", "CTRL\nT2", "LOY\nT1", "LOY\nT2",
          "CTRL\nR2", "LOY\nR1", "LOY\nR2", "LOY*",
          "CTRL\nR2", "LOY\nR1", "LOY\nR2"]
X = np.array([0, 1, 2, 3, 4.6, 5.6, 6.6, 7.6, 9.2, 10.2, 11.2])
BLUE, RED = "#3C5488", "#E64B35"
COLORS = [RED if s.startswith("LOY") or s == EXCLUDED else BLUE for s in SAMPLES]


def main():
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE, MASK]}
    d = pd.read_csv(SOURCE, sep="\t", index_col="sgRNA")
    e = d.loc[d.Gene.eq("UBR5")].copy()
    e["guide_id"] = e.index.str.split("_").str[1]
    e = e.sort_values("guide_id")
    assert len(e) == 4 and e.guide_id.tolist() == ["38607", "38608", "38609", "38610"]
    assert set(SAMPLES) == set(d.columns) - {"Gene"}
    raw = e[SAMPLES].T.copy()
    raw.columns = e.guide_id
    raw.index.name = "recorded_sample"
    raw.to_csv(DATA / "UBR5_raw_counts_all_libraries.csv")
    table = raw.copy()
    table.insert(0, "condition", ["LOY (author-confirmed)" if s == EXCLUDED else s.split("_")[0] for s in table.index])
    table.insert(1, "status", ["excluded; shown for review" if s == EXCLUDED else "retained for sensitivity analysis" for s in table.index])
    table.to_csv(DATA / "UBR5_raw_counts_annotated.csv")

    # Match the exploratory NTC scaling, calculated for the excluded sample only
    # for visualization. It is never added back to any inference or gene ranking.
    excluded_guides = set(pd.read_csv(MASK).sgRNA)
    pooled = {g: d[[f"{g}_D0_1", f"{g}_D0_2"]].sum(axis=1) for g in ["CTRL", "LOY"]}
    ref = sum(v / v.sum() * 1e6 for v in pooled.values()) / 2
    stable = (d.Gene.eq("Non-Targeting-Control") &
              d[SAMPLES[:4]].ge(50).all(axis=1) &
              ~d.index.isin(excluded_guides) & ref.gt(0))
    assert stable.sum() == 925 and not e.index.isin(excluded_guides).any()
    factors = d.loc[stable, SAMPLES].div(ref.loc[stable], axis=0).median()
    assert factors.gt(0).all()
    norm = e[SAMPLES].div(factors, axis=1).T
    norm.columns = e.guide_id
    norm.index.name = "recorded_sample"
    norm.to_csv(DATA / "UBR5_NTC_scaled_abundance_all_libraries.csv")
    factors.rename("raw_reads_per_normalized_unit").to_csv(DATA / "normalization_factors.csv")
    saved = pd.read_csv(FIGS / "source_data/a549_p3_sensitivity_20261020/normalization_factors.csv")
    saved = saved.query("excluded_guides == 342 and normalization == 'ntc'").set_index("sample")
    retained = [s for s in SAMPLES if s != EXCLUDED]
    np.testing.assert_allclose(factors.loc[retained], saved.loc[retained, "raw_counts_per_unit"])

    outputs = {}
    for mode, values, ylabel, ylim, ticks in [
        ("raw_counts", raw, "Raw mapped reads", 4000, [0, 1000, 2000, 3000, 4000]),
        ("ntc_scaled", norm, "NTC-scaled guide abundance", 18, [0, 3, 6, 9, 12, 15, 18]),
    ]:
        assert values.min().min() >= 0 and values.max().max() < ylim
        fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.15), sharey=True)
        fig.subplots_adjust(left=.085, right=.99, top=.90, bottom=.18,
                            wspace=.15, hspace=.62)
        for ax, guide, letter in zip(axes.flat, values.columns, "abcd"):
            ax.axvspan(3.8, 8.4, color="#F4F4F4", zorder=0)
            ax.axvline(3.8, color=".75", lw=.5)
            ax.axvline(8.4, color=".75", lw=.5)
            for i, s in enumerate(SAMPLES):
                val = float(values.loc[s, guide])
                marker = "X" if s == EXCLUDED else ("D" if "_D0_" in s else "o")
                ax.scatter(X[i], val, marker=marker, color=COLORS[i], s=29 if s == EXCLUDED else 20,
                           linewidth=.4, edgecolor="white", zorder=3, clip_on=False)
                text = f"{int(val):,}" if mode == "raw_counts" else f"{val:.2f}"
                ax.annotate(text, (X[i], val), xytext=(0, 5), textcoords="offset points",
                            ha="center", va="bottom", fontsize=6.5, color=COLORS[i])
            ax.set(xticks=X, xticklabels=LABELS, yticks=ticks, ylim=(0, ylim),
                   xlim=(-.5, 11.7), title=f"UBR5 guide {guide}")
            ax.tick_params(axis="x", labelsize=6.8, pad=5, length=0)
            ax.tick_params(axis="y", labelsize=7)
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
            for label, col in zip(ax.get_xticklabels(), COLORS):
                label.set_color(col)
            for x, label in [(1.5, "D0"), (6.1, "P3"), (10.2, "P4")]:
                ax.text(x, -.26, label, transform=ax.get_xaxis_transform(), ha="center",
                        fontsize=8, fontweight="semibold")
            ax.text(-.08, 1.08, letter, transform=ax.transAxes, fontsize=12, fontweight="semibold")
            ax.grid(axis="y", color="#E9EBEE", lw=.4, zorder=0)
            ax.set_axisbelow(True)
        fig.supylabel(ylabel, x=.012, fontsize=10)
        heading = "A549 UBR5: all libraries"
        qualifier = "Raw counts; linear scale" if mode == "raw_counts" else "NTC median scaling; linear scale"
        fig.suptitle(f"{heading}  |  {qualifier}", x=.52, y=.98, fontsize=11)
        fig.legend(handles=[
            Line2D([], [], marker="D", ls="", color=".4", label="D0 technical library"),
            Line2D([], [], marker="o", ls="", color=".4", label="P3/P4 culture sample"),
            Line2D([], [], marker="X", ls="", color=RED, label="Excluded LOY culture*")],
            loc="lower center", bbox_to_anchor=(.53, .06), ncol=3, frameon=False,
            fontsize=7.5, handletextpad=.4, columnspacing=1.3)
        fig.text(.085, .04, "* Recorded as CTRL_P3_R1; author-confirmed LOY. Shown for review only; remains excluded.",
                 fontsize=7.5, color=".3")
        tail = ("Library depths differ; no normalization or pseudocount applied." if mode == "raw_counts" else
                "Scaling matches the 342-guide sensitivity branch; no pseudocount. Units are not total-count CPM.")
        fig.text(.085, .016, tail, fontsize=7.5, color=".3")
        stem = f"ch04_supp_a549_ubr5_all_libraries_{mode}_review"
        for ext in ["pdf", "svg", "png"]:
            fig.savefig(FIGS / f"{stem}.{ext}", dpi=600)
        preview = FIGS / "review" / f"{stem}_preview.png"
        fig.savefig(preview, dpi=170)
        digest = hashlib.sha256(preview.read_bytes()).hexdigest()[:12]
        final_preview = preview.with_name(f"{stem}_{digest}.png")
        preview.replace(final_preview)
        plt.close(fig)
        with pymupdf.open(FIGS / f"{stem}.pdf") as pdf:
            page = pdf[0]
            assert all("OpenSans" in f[3].replace(" ", "") for f in page.get_fonts())
            outside = [w[4] for w in page.get_text("words")
                       if w[0] < 0 or w[1] < 0 or w[2] > page.rect.width or w[3] > page.rect.height]
            assert not outside, outside
        outputs[mode] = {"preview": str(final_preview), "pdf": str(FIGS / f"{stem}.pdf"),
                         "all_pdf_text_within_page": True, "font_open_sans": True,
                         "linear_y_limits": [0, ylim], "plotted_values": 44}
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == v for p, v in before.items())
    manifest = {"sources": before, "outputs": outputs, "guide_order": e.guide_id.tolist(),
                "samples": SAMPLES, "raw_values_preserved": True, "sources_unchanged": True,
                "status": "Review only. Excluded culture displayed but not reintroduced into analysis.",
                "NTC_normalization_guides": 925, "retained_scaling_matches_sensitivity_analysis": True,
                "normalization": "Median raw count / common pooled-D0 CPM reference over eligible NTC guides; "
                                 "exclude frozen 342-guide mask and require >=50 counts in all four D0 libraries.",
                "D0_samples": "Technical libraries; not independent cultures.",
                "P3_P4_samples": "Repeated samples of the same continuing cultures."}
    (DATA / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(table.to_string())
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
