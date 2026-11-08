"""Publication candidate panels from frozen exploratory DepMap analyses.

No gene re-selection, fitting, new significance tests or manuscript integration.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / "figures/chapter4"
SOURCE = OUT / "source_data"
DATA = SOURCE / "ddx3x_publication_panels"
sys.path.insert(0, str(THESIS / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style

GROUPS = ["XY", "LOY", "Female"]
COLOURS = {"XY": "#3C5488", "LOY": "#E64B35", "Female": "#008577",
           "XY_low_DDX3Y": "#71A2C1", "XY_high_DDX3Y": "#3C5488"}
LABELS = {"XY": "Y-retaining", "LOY": "Y-depleted", "Female": "Female"}
LINESTYLES = {"XY": "-", "LOY": "--", "Female": ":"}
CMAP = LinearSegmentedColormap.from_list(
    "signed_correlation", ["#3C5488", "#F7F7F7", "#C96C2D"])
INPUTS = []
EXPORTS = []


def read(rel):
    path = SOURCE / rel
    INPUTS.append(path)
    return pd.read_csv(path)


def fmt(value):
    text = f"{value:.2f}"
    return "0.00" if text == "-0.00" else text.replace("-", "−")


def panel(fig, letter, x, y):
    fig.text(x, y, letter, fontweight="semibold", fontsize=12)


def export(fig, stem, info):
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT / f"{stem}.{ext}", dpi=600 if ext == "png" else 300)
    fig.savefig(OUT / "review" / f"{stem}_preview.png", dpi=155)
    dimensions = fig.get_size_inches().tolist()
    plt.close(fig)
    with fitz.open(OUT / f"{stem}.pdf") as doc:
        page = doc[0]
        assert all("OpenSans" in font[3] for font in page.get_fonts())
        spans = [s for b in page.get_text("dict")["blocks"] if "lines" in b
                 for line in b["lines"] for s in line["spans"] if s["text"].strip()]
        outside = [s["text"] for s in spans
                   if not page.rect.contains(fitz.Rect(s["bbox"]))]
        overlaps = []
        for i, a in enumerate(spans):
            for b in spans[i+1:]:
                ra, rb = fitz.Rect(a["bbox"]), fitz.Rect(b["bbox"])
                area = (ra & rb).get_area()
                if area > 1 and area / min(ra.get_area(), rb.get_area()) > .12:
                    overlaps.append([a["text"], b["text"]])
        assert not outside, outside
        assert not overlaps, overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), colorspace=fitz.csGRAY).save(
            OUT / "review" / f"{stem}_greyscale.png")
    record = dict(stem=stem, dimensions_inches=dimensions, font="Open Sans",
                  text_outside_page=outside, text_overlaps=overlaps, **info)
    (DATA / f"{stem}_render_manifest.json").write_text(
        json.dumps(record, indent=2), encoding="utf-8")
    EXPORTS.append(record)
    print(stem, "passed font, overlap and page-boundary checks", flush=True)


def prepare():
    selected = read("ddx3x_genomewide_y_pattern/selected_partners.csv")
    summary = read("xy_paralog_genomewide_y_pattern/anchor_summary.csv")
    sensitivity = read("xy_paralog_genomewide_y_pattern/threshold_sensitivity_counts.csv")
    same = []
    for anchor in summary.anchor:
        source = read(f"xy_paralog_genomewide_y_pattern/{anchor}_all_gene_correlations.csv.gz")
        rows = source.set_index("gene").loc[selected.gene].reset_index()
        assert len(rows) == 31 and rows.eligible.all()
        same.append(rows)
    same = pd.concat(same, ignore_index=True)
    ddx = same.loc[same.anchor.eq("DDX3X")].set_index("gene").loc[selected.gene]
    for col in ["xy_r", "loy_r", "female_r"]:
        np.testing.assert_allclose(ddx[col], selected[col], atol=1e-10, rtol=0)
    assert len(same) == 496
    assert summary.selected_partners.sum() == 32
    assert summary.loc[summary.anchor.eq("DDX3X"), "selected_partners"].item() == 31
    assert same.loc[~same.anchor.eq("DDX3X"), "selected"].sum() == 0
    same.to_csv(DATA / "same_31_partners_all_16_anchors.csv", index=False)
    selected.to_csv(DATA / "displayed_ddx3x_partners.csv", index=False)
    summary.to_csv(DATA / "genomewide_pattern_counts.csv", index=False)
    sensitivity.to_csv(DATA / "pattern_threshold_sensitivity.csv", index=False)
    residuals = read("ddx3x_partner_adjustment/residuals.csv")
    residuals = residuals.loc[residuals.analysis.eq("lineage_efficacy")]
    stats = read("ddx3x_partner_adjustment/correlations.csv")
    stats = stats.loc[stats.analysis.eq("lineage_efficacy")]
    for row in stats.itertuples():
        values = residuals.loc[residuals.group.eq(row.group)]
        assert len(values) == row.n_models
        np.testing.assert_allclose(np.corrcoef(values.DDX3X, values[row.partner])[0, 1],
                                   row.pearson_r, atol=1e-10, rtol=0)
    residuals.to_csv(DATA / "representative_scatter_residuals.csv", index=False)
    stats.to_csv(DATA / "representative_correlations_and_intervals.csv", index=False)
    return selected, summary, sensitivity, same, residuals, stats


def clean_heat(ax, rows, cols):
    ax.set_xticks(np.arange(-.5, cols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, rows, 1), minor=True)
    ax.grid(which="minor", color="white", lw=.6)
    ax.tick_params(axis="both", which="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def main_panel(selected, summary, sensitivity, residuals, stats):
    apply_style(8)
    fig = plt.figure(figsize=(7.8, 8.0))
    panel(fig, "a", .024, .964)
    handles = [Line2D([], [], marker="o", linestyle=LINESTYLES[g], lw=1, markersize=4,
                      color=COLOURS[g], label=f"{LABELS[g]} (n = {n})")
               for g, n in zip(GROUPS, [324, 229, 472])]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.53, .991),
               ncol=3, frameon=False, fontsize=8, handletextpad=.4, columnspacing=1.4)
    for j, partner in enumerate(["IGBP1", "GBF1", "IK"]):
        ax = fig.add_axes([.09+j*.31, .711, .25, .216])
        for g in GROUPS:
            frame = residuals.loc[residuals.group.eq(g)]
            row = stats.loc[stats.group.eq(g) & stats.partner.eq(partner)].iloc[0]
            ax.scatter(frame.DDX3X, frame[partner], s=5.5, alpha=.24,
                       c=COLOURS[g], edgecolors="none", rasterized=False)
            grid = np.linspace(frame.DDX3X.min(), frame.DDX3X.max(), 80)
            ax.plot(grid, row.residual_intercept+row.residual_slope*grid,
                    c=COLOURS[g], ls=LINESTYLES[g], lw=1.1, zorder=4)
        ax.set_xlim(-1.6, 1.25)
        low, high = residuals[partner].min(), residuals[partner].max()
        ax.set_ylim(low-.08*(high-low), high+.48*(high-low))
        ax.set_xticks([-1, 0, 1])
        ax.tick_params(labelsize=7.4)
        ax.set_title(partner, fontsize=9, fontweight="semibold", pad=6)
        for k, g in enumerate(GROUPS):
            row = stats.loc[stats.group.eq(g) & stats.partner.eq(partner)].iloc[0]
            text = f"r = {fmt(row.pearson_r)} [{fmt(row.bootstrap_ci_low)}, {fmt(row.bootstrap_ci_high)}]"
            ax.text(.025, .967-k*.103, text, transform=ax.transAxes,
                    color=COLOURS[g], fontsize=7.0, va="top",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.84, pad=.6))
        if j == 0:
            ax.set_ylabel("Partner gene-effect residual", fontsize=8)
    fig.text(.52, .661, "DDX3X gene-effect residual", ha="center", fontsize=8.5)

    panel(fig, "b", .024, .615)
    ax = fig.add_axes([.145, .082, .285, .491])
    mat = selected[["xy_r", "loy_r", "female_r"]].to_numpy()
    assert np.abs(mat).max() <= .5
    im = ax.imshow(mat, cmap=CMAP, vmin=-.5, vmax=.5, aspect="auto")
    ax.set_yticks(range(31), selected.gene)
    ax.set_xticks(range(3), ["Y-retaining", "Y-depleted", "Female"])
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", labelsize=7.8, pad=6)
    ax.tick_params(axis="y", labelsize=7.0, pad=5)
    clean_heat(ax, 31, 3)
    for i, gene in enumerate(selected.gene):
        if gene in ["IGBP1", "GBF1", "IK", "PSMD8"]:
            ax.get_yticklabels()[i].set_fontweight("semibold")
        for j in range(3):
            val = mat[i, j]
            ax.text(j, i, fmt(val), ha="center", va="center", fontsize=7.0,
                    color="white" if abs(val) > .38 else "#242424")
    cbax = fig.add_axes([.18, .051, .21, .009])
    cb = fig.colorbar(im, cax=cbax, orientation="horizontal", ticks=[-.5, 0, .5])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7.3, length=2, pad=2)
    cb.set_label("Adjusted Pearson r with DDX3X", fontsize=7.8, labelpad=2)

    panel(fig, "c", .50, .615)
    ax = fig.add_axes([.657, .158, .28, .395])
    sens = sensitivity.loc[np.isclose(sensitivity.max_abs_xy_r, .15) &
                           np.isclose(sensitivity.min_delta, .20)]
    counts = sens.pivot(index="anchor", columns="min_positive_r",
                        values="selected_partners").loc[summary.anchor, [.20, .25, .30]]
    values = counts.to_numpy()
    ax.imshow(values, cmap="magma_r", vmin=0, vmax=int(values.max()), aspect="auto")
    ax.set_yticks(range(16), counts.index)
    ax.set_xticks(range(3), ["0.20", "0.25", "0.30"])
    ax.xaxis.tick_top()
    ax.tick_params(axis="both", labelsize=8.1, pad=6)
    clean_heat(ax, 16, 3)
    ax.add_patch(Rectangle((.5, -.5), 1, 16, fill=False, edgecolor="#222222",
                           lw=1.1, clip_on=False))
    for i in range(16):
        if counts.index[i] == "DDX3X":
            ax.get_yticklabels()[i].set_fontweight("semibold")
        for j in range(3):
            ax.text(j, i, str(values[i, j]), ha="center", va="center", fontsize=8.4,
                    color="white" if values[i, j] > 35 else "#222222")
    ax.text(.5, 1.11, "Minimum r in LOY and female", ha="center",
            transform=ax.transAxes, fontsize=8.2)
    fig.text(.785, .119, "Number of matching partners", ha="center", fontsize=8.2)
    fig.text(.785, .091, "17,786 eligible partners per anchor", ha="center",
             fontsize=7.5, color="#555555")
    fig.text(.785, .063, "Fixed: |r in XY| ≤ 0.15; each Δr ≥ 0.20", ha="center",
             fontsize=7.4, color="#555555")
    fig.text(.785, .036, "Outlined column: original descriptive rule", ha="center",
             fontsize=7.4, color="#555555")
    export(fig, "ch04_main_ddx3x_context_publication_draft",
           dict(partners=31, anchors=16, no_new_tests=True,
                scatter_intervals="Saved 95% donor-bootstrap correlation intervals.",
                counts="Descriptive threshold counts; not significant-difference counts.",
                color_limits=[-.5, .5], no_color_clipping=True))


def companion():
    apply_style(8)
    expr = read("ddx3y_expression_dependency/expression_dependency_residuals.csv")
    expr = expr.loc[expr.analysis.eq("primary")]
    estat = read("ddx3y_expression_dependency/ddx3x_expression_associations.csv")
    row = estat.loc[estat.analysis.eq("primary")].iloc[0]
    ps = read("psmd8_ddx3y_followup/plotted_residuals.csv")
    pst = read("psmd8_ddx3y_followup/correlations_and_sensitivities.csv")
    pst = pst.loc[pst.analysis.eq("primary")].set_index("group")
    interactions = read("proteasome_ddx3y_context/interaction_results.csv")
    # Preserve both original correction families in the caption and provenance.
    old = read("ddx3y_expression_dependency/expression_partner_interactions.csv")
    interaction = interactions.loc[interactions.analysis.eq("primary") &
                                   interactions.gene.eq("PSMD8")].iloc[0]
    original = old.loc[old.analysis.eq("primary") & old.gene.eq("PSMD8")].iloc[0]
    np.testing.assert_allclose(interaction.p_value, original.p_value, atol=1e-10)
    expr.to_csv(DATA / "expression_scatter_residuals.csv", index=False)
    ps.to_csv(DATA / "psmd8_scatter_residuals.csv", index=False)
    pst.to_csv(DATA / "psmd8_correlations_with_intervals.csv")
    fig = plt.figure(figsize=(7.8, 5.35))
    panel(fig, "a", .025, .958)
    ax = fig.add_axes([.103, .61, .335, .30])
    x, y = expr.expression_residual, expr.dependency_residual
    assert len(x) == row.n_models == 305
    np.testing.assert_allclose(np.corrcoef(x, y)[0, 1], row.partial_r, atol=1e-10)
    ax.scatter(x, y, c=COLOURS["XY"], alpha=.4, s=8, edgecolors="none")
    grid = np.linspace(x.min(), x.max(), 100)
    ax.plot(grid, row.effect*grid, c="#222222", lw=1.1)
    ax.margins(x=.06, y=.18)
    ax.set_xlabel("DDX3Y expression residual", fontsize=8)
    ax.set_ylabel("DDX3X gene-effect residual", fontsize=8)
    p_text = f"{row.p_value:.2e}".replace("e-0", "e−")
    ax.text(.03, .97, f"r = {fmt(row.partial_r)} [{fmt(row.r_ci_low)}, {fmt(row.r_ci_high)}]\n"
            f"P = {p_text}; n = {int(row.n_models)}", transform=ax.transAxes, va="top",
            fontsize=7.6, bbox=dict(facecolor="white", edgecolor="none", alpha=.88, pad=1))
    panel(fig, "b", .51, .958)
    ax = fig.add_axes([.727, .61, .23, .30])
    groups = ["XY_low_DDX3Y", "XY_high_DDX3Y", "LOY", "Female"]
    labels = ["XY, lower DDX3Y", "XY, higher DDX3Y", "Y-depleted", "Female"]
    ax.axvline(0, color="#999999", lw=.65, ls="--")
    for i, g in enumerate(groups):
        r = pst.loc[g]
        ax.plot([r.r_ci_low, r.r_ci_high], [i, i], color=COLOURS[g], lw=1.4)
        ax.plot(r.pearson_r, i, "o", color=COLOURS[g], ms=4.5)
        ax.text(.99, i-.26, f"r = {fmt(r.pearson_r)}; n = {int(r.n_models)}",
                transform=ax.get_yaxis_transform(), ha="right", fontsize=7.2, color="#444444")
    ax.set(yticks=range(4), yticklabels=labels, ylim=(3.65, -.65),
           xlim=(-.5, .5), xticks=[-.5, 0, .5])
    ax.tick_params(axis="y", length=0, pad=6, labelsize=8)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("DDX3X–PSMD8 adjusted r\n95% donor-bootstrap interval", fontsize=8)
    panel(fig, "c", .025, .453)
    paired = ps.loc[ps.group.isin(groups[:2])]
    xlo, xhi = paired.ddx3x_residual.min(), paired.ddx3x_residual.max()
    ylo, yhi = paired.psmd8_residual.min(), paired.psmd8_residual.max()
    xlim = (xlo-.06*(xhi-xlo), xhi+.06*(xhi-xlo))
    ylim = (ylo-.06*(yhi-ylo), yhi+.28*(yhi-ylo))
    for j, g in enumerate(groups[:2]):
        ax = fig.add_axes([.103+j*.506, .135, .335, .265])
        frame = ps.loc[ps.group.eq(g)]
        r = pst.loc[g]
        x, y = frame.ddx3x_residual, frame.psmd8_residual
        assert len(frame) == r.n_models
        np.testing.assert_allclose(np.corrcoef(x, y)[0, 1], r.pearson_r, atol=1e-10)
        ax.scatter(x, y, s=9, alpha=.47, color=COLOURS[g], edgecolors="none")
        slope, intercept = np.polyfit(x, y, 1)
        grid = np.linspace(x.min(), x.max(), 100)
        ax.plot(grid, intercept+slope*grid, color=COLOURS[g], lw=1.25)
        ax.set(xlim=xlim, ylim=ylim)
        assert x.min() > xlim[0] and x.max() < xlim[1]
        assert y.min() > ylim[0] and y.max() < ylim[1]
        ax.set_xlabel("DDX3X gene-effect residual", fontsize=8)
        ax.set_ylabel("PSMD8 gene-effect residual", fontsize=8)
        ax.set_title(labels[j], fontsize=8.6, pad=5)
        ax.text(.025, .965, f"r = {fmt(r.pearson_r)}; n = {int(r.n_models)}",
                transform=ax.transAxes, va="top", fontsize=7.6,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=1))
    fig.text(.52, .027, f"Continuous DDX3Y × PSMD8: P = {interaction.p_value:.4f}; "
             f"proteasome-set FDR = {interaction.fdr_bh_set:.3f}",
             ha="center", fontsize=8, color="#444444")
    export(fig, "ch04_main_ddx3y_psmd8_publication_draft",
           dict(expression_models=305, psmd8_strata=[153, 152, 216, 448],
                expression_test="Saved donor-clustered regression test.",
                interaction_p=float(interaction.p_value),
                interaction_fdr_proteasome_65=float(interaction.fdr_bh_set),
                interaction_fdr_original_31=float(original.fdr_bh_31_within_analysis),
                intervals="Saved nominal 95% donor-bootstrap intervals.",
                no_new_tests=True))


def comparison(selected, summary, same):
    apply_style(8)
    anchors = summary.anchor.tolist()
    genes = selected.gene.tolist()
    matrices = [same.loc[same.anchor.eq(a)].set_index("gene").loc[
        genes, ["xy_r", "loy_r", "female_r"]].to_numpy() for a in anchors]
    matrix = np.column_stack(matrices)
    assert matrix.shape == (31, 48) and np.isfinite(matrix).all()
    assert np.abs(matrix).max() <= .5
    fig = plt.figure(figsize=(13.2, 7.5))
    ax = fig.add_axes([.075, .13, .91, .73])
    im = ax.imshow(matrix, cmap=CMAP, vmin=-.5, vmax=.5, aspect="auto")
    ax.set_yticks(range(31), genes)
    ax.set_xticks(range(48), ["XY", "LOY", "F"]*16)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", labelsize=6.2, pad=4)
    ax.tick_params(axis="y", labelsize=8.4, pad=6)
    clean_heat(ax, 31, 48)
    for k, anchor in enumerate(anchors):
        for i in range(31):
            for j in range(3):
                val = matrices[k][i, j]
                ax.text(3*k+j, i, fmt(val), ha="center", va="center", fontsize=6.2,
                        color="white" if abs(val) > .38 else "#303030")
        ax.text(k*3+1, 1.055, anchor, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=8.2,
                fontweight="semibold" if anchor == "DDX3X" else "normal")
        if k > 0:
            ax.axvline(k*3-.5, color="#777777", lw=.75)
    ax.add_patch(Rectangle((-.5, -.5), 3, 31, fill=False,
                           edgecolor="#222222", lw=1.2, clip_on=False))
    cbax = fig.add_axes([.37, .065, .29, .012])
    cb = fig.colorbar(im, cax=cbax, orientation="horizontal", ticks=[-.5, 0, .5])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=2, labelsize=8)
    cb.set_label("Adjusted Pearson r: anchor versus partner", fontsize=8.5, labelpad=3)
    fig.text(.075, .954, "XY: Y-retaining male   •   LOY: Y-depleted male   •   F: recorded female",
             fontsize=8.5)
    export(fig, "ch04_supp_ddx3x_gametolog_comparison_publication_draft",
           dict(partners=31, anchors=16, cells=1488, color_limits=[-.5, .5],
                no_color_clipping=True, no_new_tests=True,
                selection="31 partners selected on DDX3X correlations; conditional comparison."))


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    (OUT / "review").mkdir(exist_ok=True)
    selected, summary, sensitivity, same, residuals, stats = prepare()
    main_panel(selected, summary, sensitivity, residuals, stats)
    companion()
    comparison(selected, summary, same)
    manifest = dict(
        purpose="Draft publication candidate panels; no final figure assignment or Results integration.",
        dependency_release="DepMap Public 25Q3",
        annotations="Previously joined Chapter 3 26Q1 annotations.",
        selection="Frozen descriptive DDX3X rule: LOY and female r>=.25, abs(XY r)<=.15, both delta>=.20; >=100 complete models/group.",
        specificity="31 matches for DDX3X, one TGIF2LX–OR51A2 match, zero for 14 other anchors. Not genome-wide uniqueness or a formal specificity test.",
        same_partner_comparison="The fixed 31 DDX3X-selected partners have zero pattern matches for any other tested anchor; this is conditional on DDX3X selection.",
        cohorts="Main: 324 XY, 229 LOY, 472 female for DDX3X. Companion RNA subset: 305 XY, 216 LOY, 448 female. Female does not imply verified XX karyotype.",
        validation="Every main scatter and DDX3X matrix reproduced from saved values to 1e-10; PSMD8 and expression scatter correlations reproduced; no hidden clipping.",
        sources={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in sorted(set(INPUTS + [Path(__file__)]))},
        exports=EXPORTS)
    (DATA / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("All publication candidate figures and source tables saved.", flush=True)


if __name__ == "__main__":
    main()
