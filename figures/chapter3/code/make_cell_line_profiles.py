"""Thesis cell-line profiles from frozen source tables; no analysis in plotting.

Run with --background-only while quantitative sources are being reconciled.
The complete figures require figure_config.json and the CSVs listed below.
"""
from pathlib import Path
import argparse
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
import numpy as np
import pandas as pd

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style

BLUE, RED = "#3C5488", "#E64B35"
COLOURS = {"XY": BLUE, "LOY": RED}
LINES = ["HT1080", "A549", "DU145", "H2452"]
MARKERS = dict(zip(LINES, ["o", "s", "^", "D"]))
GENES = ["DDX3X", "EIF1AX", "RPS4X", "ZFX"]


def heading(ax, letter, title):
    ax.text(-.14, 1.05, letter, transform=ax.transAxes,
            weight="semibold", fontsize=11, va="bottom")
    ax.set_title(title, loc="left", pad=13, fontsize=9, weight="semibold")


def background(fig, rect, meta):
    """Reusable background strip with publication/repository properties."""
    ax = fig.add_axes(rect)
    ax.set_axis_off()
    ax.text(0, .98, "a", weight="semibold", fontsize=11, va="top")
    ax.text(.035, .98, f"{meta['cell_line']}  |  {meta['reported_disease']}",
            fontsize=11, weight="semibold", va="top")
    ax.text(1, .98, f"RRID:{meta['rrid']}", fontsize=7.5,
            color="#666666", ha="right", va="top")
    ax.plot([0, 1], [.85, .85], lw=.6, color="#B9B9B9")
    cols = [0, .365, .70]
    titles = ["Origin and donor", "Growth and karyotype", "Reported alterations"]
    for x, title in zip(cols, titles):
        ax.text(x, .78, title, fontsize=8.7, weight="semibold", va="top")
    origin = [
        "Pelvic fibrosarcoma",
        "Adjacent to the acetabulum",
        f"{meta['donor_age_years']}-year-old male",
        "No prior chemotherapy/radiotherapy",
        f"Established {meta['established_year']}; reported {meta['first_report_year']}",
    ]
    culture = [
        "Adherent monolayer",
        "Doubling time: ~30 h (DSMZ)",
        "Pseudodiploid; modal count 46",
        "ATCC CCL-121 / DSMZ ACC 315",
    ]
    for x, lines in [(cols[0], origin), (cols[1], culture)]:
        for i, s in enumerate(lines):
            ax.text(x, .63-i*.115, s, fontsize=7.8, va="top")
    for i, variant in enumerate(meta["reported_alterations"]):
        y = .63-i*.175
        text = variant["gene"] + " " + variant["alteration"]
        ax.text(cols[2], y, text, fontsize=8.2, weight="semibold", va="top")
        ax.text(cols[2], y-.078, variant["zygosity"].lower(),
                fontsize=7.4, color="#666666", va="top")
    ax.text(0, -.02, "Published/repository properties; doubling time depends on culture conditions.",
            fontsize=7.1, color="#666666", va="top")
    ax.set_xlim(0, 1)
    ax.set_ylim(-.1, 1)
    return ax


def export(fig, out, stem):
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=600)
    fig.savefig(out / f"{stem}_preview.png", dpi=170)
    plt.close(fig)


def ordered_y_genes(source):
    """Use frozen GRCh38 gene coordinates for order and chromosome positions."""
    coords = pd.read_csv(source / 'reference_material/y_gene_coordinates_grch38.csv')
    assert coords.assembly.eq('GRCh38').all() and coords.chromosome.eq('Y').all()
    return coords.sort_values('start').reset_index(drop=True)


def chromosome_diagram(fig, rect, source, coords):
    """Full GRCh38 chromosome to scale, linked to equally spaced expression rows."""
    bands = pd.read_csv(source / 'reference_material/chrY_cytobands_grch38.csv')
    length = bands.end.max()
    ax = fig.add_axes(rect)
    n = len(coords)
    ax.set_xlim(0, 1)
    ax.set_ylim(n-.5, -.5)
    ax.set_axis_off()
    y = lambda bp: -.5 + bp / length * (n-.0)
    left, right = .24, .36
    shades = {'gneg':'#FFFFFF','gpos25':'#DEDEDE','gpos50':'#B9B9B9',
              'gpos75':'#909090','gpos100':'#686868','gvar':'#EAEAEA'}
    for band in bands.itertuples():
        y0, y1 = y(band.start), y(band.end)
        if band.stain == 'acen':
            if str(band.band).startswith('p'):
                xy = [(left,y0),(right,y0),((left+right)/2,y1)]
            else:
                xy = [((left+right)/2,y0),(right,y1),(left,y1)]
            ax.add_patch(Polygon(xy, facecolor='#707070', edgecolor='#777777', lw=.45))
        else:
            ax.add_patch(Rectangle((left,y0),right-left,y1-y0,
                                   facecolor=shades[band.stain],edgecolor='#777777',lw=.45))
    centromere = bands.loc[bands.stain.eq('acen')]
    for label, bp in [('Yp',centromere.start.min()/2),
                      ('Yq',(centromere.end.max()+length)/2)]:
        ax.text(.13,y(bp),label,rotation=90,ha='center',va='center',fontsize=7,color='#555555')
    ax.text(.30,-.75,'chrY',ha='center',va='bottom',fontsize=7.3,clip_on=False)
    ax.text(.30,n-.35,'GRCh38\n57.2 Mb',ha='center',va='top',fontsize=6.1,color='#555555',clip_on=False)
    for i, row in enumerate(coords.itertuples()):
        midpoint = (row.start - 1 + row.end) / 2
        gy = y(midpoint)
        ax.plot([right,.44,.94],[gy,gy,i],lw=.55,color='#A8A8A8',zorder=1)
        ax.plot([left,right],[gy,gy],lw=.8,color=BLUE,zorder=3)
    return ax


def background_with_photo(fig, source, meta):
    fig.text(.045,.972,f"{meta['cell_line']}  |  {meta['reported_disease']}",
             fontsize=11,weight='semibold',va='top')
    fig.text(.955,.969,f"RRID:{meta['rrid']}",fontsize=7.5,color='#666666',ha='right',va='top')
    fig.text(.045,.925,'a',fontsize=11,weight='semibold',va='top')
    fig.text(.079,.923,'Origin and background',fontsize=9,weight='semibold',va='top')
    groups = [
        (.881,'Origin and donor',[
            'Pelvic tumour adjacent to the acetabulum',
            f"{meta['donor_age_years']}-year-old man; no prior chemo-/radiotherapy",
            f"Established {meta['established_year']}; first reported {meta['first_report_year']}"]),
        (.787,'Growth and karyotype',[
            'Adherent monolayer; doubling time ~30 h (DSMZ)',
            'Pseudodiploid; modal chromosome count 46']),
        (.711,'Reported genotype',[
            f"{v['gene']} {v['alteration']} ({v['zygosity'].lower()})"
            for v in meta['reported_alterations']]
            + ([meta['reported_tp53']['display']] if 'reported_tp53' in meta else [])),
    ]
    for pos, title, lines in groups:
        fig.text(.045,pos,title,fontsize=8.1,weight='semibold',va='top')
        for i, text in enumerate(lines):
            fig.text(.045,pos-.025-i*.018,text,fontsize=7.5,va='top')
    fig.text(.555,.925,'b',fontsize=11,weight='semibold',va='top')
    fig.text(.589,.923,'Reference culture morphology',fontsize=9,weight='semibold',va='top')
    image_file = source / 'reference_material/HT-1080_Bild_2.jpg'
    ax=fig.add_axes([.555,.677,.40,.206])
    ax.imshow(plt.imread(image_file),cmap='gray',vmin=0,vmax=255,
              interpolation='none',aspect='equal')
    ax.set_axis_off()
    fig.text(.555,.662,'© DSMZ · HT-1080, ACC 315',fontsize=7.1,color='#555555',va='top')
    fig.text(.555,.642,'Reference image; original scale bar: 50 µm',fontsize=7.1,color='#555555',va='top')


def profile(source, out, config, meta):
    expr = pd.read_csv(source / "ht1080_y_expression.csv")
    deps = pd.read_csv(source / "ht1080_dependency_context.csv")
    coords = ordered_y_genes(source)
    expr = expr.set_index('gene').loc[coords.gene].reset_index()
    fig = plt.figure(figsize=(183/25.4, 200/25.4))
    background_with_photo(fig, source, meta)
    fig.text(.045,.585,'c',fontsize=11,weight='semibold',va='bottom')
    fig.text(.079,.585,'Y-linked expression in chromosome order',fontsize=9,weight='semibold',va='bottom')
    chromosome_diagram(fig,[.035,.16,.205,.385],source,coords)
    ax = fig.add_axes([.30,.16,.195,.385])
    yy = np.arange(len(expr))
    ax.barh(yy, expr["expression"], color=BLUE, height=.60)
    ax.set_yticks(yy, expr["gene"], fontsize=8)
    ax.set_ylim(len(expr)-.5,-.5)
    ax.set_xlabel("Expression [log₂(TPM + 1)]")
    ax.set_xlim(0, max(expr.expression)*1.18)
    for y, val in zip(yy, expr.expression):
        ax.text(val+.08, y, f"{val:.2f}", va="center", fontsize=7)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax = fig.add_axes([.69, .16, .27, .385])
    fig.text(.595,.585,'d',fontsize=11,weight='semibold',va='bottom')
    fig.text(.629,.585,'Selected genetic dependencies',fontsize=9,weight='semibold',va='bottom')
    nref = sorted(deps.n_reference.unique())
    nref_text = f'{nref[0]:,}' if len(nref) == 1 else f'{min(nref):,}–{max(nref):,}'
    ax.text(0, 1.008, f'Other models: n = {nref_text}', transform=ax.transAxes,
            fontsize=7, color='#666666', va='bottom')
    yy = np.arange(len(deps))
    for y, row in enumerate(deps.itertuples()):
        ax.plot([row.reference_q25, row.reference_q75], [y, y],
                color="#BDBDBD", lw=3, solid_capstyle='butt', zorder=1)
    ax.scatter(deps.reference_median, yy, marker="|", s=100,
               color="#777777", label="Reference median", zorder=3)
    ax.scatter(deps.gene_effect, yy, s=27, color=BLUE,
               label="HT1080", zorder=4)
    ax.set_yticks(yy, deps.gene, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Chronos gene effect")
    ax.axvline(0, lw=.5, color="#BDBDBD", zorder=0)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper left", bbox_to_anchor=(-.05, -.15),
              ncol=2, fontsize=6.7, handletextpad=.3, columnspacing=.8)
    fig.text(.045, .065, f"{config['release_label']}  |  Baseline public cell-line profiles", fontsize=7.3, color="#666666")
    fig.text(.045, .045, 'Background and image: repository properties; doubling time depends on culture conditions.',
             fontsize=7.1,color='#666666')
    export(fig, out, "ch03_main_ht1080_profile")


def expression_heatmap(source, out, config):
    data = pd.read_csv(source / 'selected_models_y_expression.csv')
    coords = ordered_y_genes(source)
    genes = coords.gene.tolist()
    values = data.pivot(index='gene', columns='cell_line', values='expression').reindex(index=genes, columns=LINES)
    fig, ax = plt.subplots(figsize=(183/25.4, 112/25.4))
    fig.subplots_adjust(left=.30, right=.83, top=.86, bottom=.15)
    chromosome_diagram(fig,[.04,.15,.19,.71],source,coords)
    im = ax.imshow(values, cmap='magma', vmin=0, vmax=10, aspect='auto', interpolation='nearest')
    ax.set_xticks(range(4), LINES, fontsize=9)
    ax.set_yticks(range(len(genes)), genes, fontsize=8.5)
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', length=0, pad=6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(len(genes)):
        for j in range(4):
            val = values.iloc[i,j]
            ax.text(j, i, f'{val:.2f}' if pd.notna(val) else 'NA', ha='center', va='center',
                    color='white' if val < 5.5 else '#222222', fontsize=8)
    ax.set_xticks(np.arange(-.5,4,1), minor=True)
    ax.set_yticks(np.arange(-.5,len(genes),1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.65)
    ax.tick_params(which='minor', length=0)
    cbax = fig.add_axes([.865,.31,.018,.43])
    cb = fig.colorbar(im, cax=cbax, ticks=[0,2,4,6,8,10])
    cb.set_label('Expression [log₂(TPM + 1)]', fontsize=8)
    cb.outline.set_visible(False)
    fig.text(.30, .065, f"{config['release_label']}  |  Default RNA profile for each model", fontsize=7.5, color='#555555')
    export(fig, out, 'ch03_supp_selected_models_y_expression')


def comparison(source, out, config):
    data = pd.read_csv(source / "dependency_comparison.csv")
    stats_path = source / "dependency_statistics.csv"
    stats = pd.read_csv(stats_path).set_index("gene") if stats_path.exists() else None
    if data.duplicated(["ModelID", "gene"]).any():
        raise ValueError("Duplicate model/gene entries")
    fig, axs = plt.subplots(2, 2, figsize=(183/25.4, 169/25.4))
    fig.subplots_adjust(left=.20, right=.97, top=.89, bottom=.15, hspace=.75, wspace=.68)
    for ax, gene, letter in zip(axs.flat, GENES, "abcd"):
        sub = data.loc[data.gene.eq(gene)]
        title = gene
        title_y = 1.14 if config.get('statistics_kind') == 'dna' else 1.085
        ax.text(-.14,title_y,letter,transform=ax.transAxes,fontsize=11,weight='semibold',va='bottom')
        ax.text(0,title_y,title,transform=ax.transAxes,fontsize=9,weight='semibold',va='bottom')
        for y, group in enumerate(["XY", "LOY"]):
            vals = sub.loc[sub.LOY_status.eq(group), "gene_effect"].dropna().to_numpy()
            parts = ax.violinplot(vals, positions=[y], vert=False,
                                  widths=.72, showextrema=False)
            for body in parts["bodies"]:
                body.set_facecolor(COLOURS[group])
                body.set_edgecolor(COLOURS[group])
                body.set_alpha(.28)
            ax.boxplot(vals, positions=[y], orientation="horizontal", widths=.15,
                       showfliers=False, patch_artist=True,
                       boxprops={"facecolor": COLOURS[group], "edgecolor": COLOURS[group]},
                       whiskerprops={"color": COLOURS[group], "linewidth": .7},
                       capprops={"color": COLOURS[group], "linewidth": .7},
                       medianprops={"color": "white", "linewidth": 1})
        labels = [f"{config.get('group_labels', {}).get(group, group)}\n(n = {sub.loc[sub.LOY_status.eq(group), 'gene_effect'].notna().sum()})"
                  for group in ["XY", "LOY"]]
        for y, line in enumerate(LINES, start=2):
            rows = sub.loc[sub.cell_line.eq(line)]
            if len(rows) != 1:
                raise ValueError(f"Expected one {line}/{gene} row, got {len(rows)}")
            row = rows.iloc[0]
            ax.scatter(row.gene_effect, y, marker=MARKERS[line], s=32,
                       color=COLOURS.get(row.LOY_status, "#888888"),
                       edgecolor="#222222", linewidth=.55, zorder=5)
            ax.annotate(f'{row.gene_effect:.2f}', (row.gene_effect,y), xytext=(6,0),
                        textcoords='offset points',va='center',fontsize=7)
            labels.append(line)
        ax.axhline(1.55, color="#D6D6D6", lw=.6)
        ax.axvline(0, color="#BDBDBD", lw=.5, zorder=0)
        ax.set_yticks(range(6), labels, fontsize=7.6)
        ax.set_ylim(5.65, -.7)
        ax.set_xlabel("Chronos gene effect", fontsize=8)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.margins(x=.22)
        if stats is not None and gene in stats.index:
            ax.text(0, 1.015, str(stats.loc[gene, "plot_annotation"]),
                    transform=ax.transAxes, fontsize=7, color="#555555")
    fig.text(.045, .070, config['release_label'], fontsize=7.2, color="#555555")
    fig.text(.045, .047, config['comparison_label'], fontsize=7.2, color="#555555")
    fig.text(.045, .024, "More negative gene effects indicate stronger dependency. Axis ranges vary by gene.", fontsize=7.2, color="#555555")
    export(fig, out, config.get('output_stem', 'ch03_main_selected_models_paralog_dependencies'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=THESIS / "figures/chapter3/source_data/cell_line_profiles")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--background-only", action="store_true")
    parser.add_argument("--profile-only", action="store_true")
    parser.add_argument("--comparison-only", action="store_true")
    args = parser.parse_args()
    apply_style(base_size=8)
    if args.comparison_only:
        config = json.loads((args.source_dir / 'figure_config.json').read_text(encoding='utf-8'))
        comparison(args.source_dir, args.output_dir, config)
        return
    meta = json.loads((args.source_dir / "ht1080_background.json").read_text(encoding="utf-8"))
    if args.background_only:
        fig = plt.figure(figsize=(183/25.4, 64/25.4))
        background(fig, [.04, .12, .92, .82], meta)
        export(fig, args.output_dir, "ch03_main_ht1080_background_draft")
    else:
        config = json.loads((args.source_dir / "figure_config.json").read_text(encoding="utf-8"))
        profile(args.source_dir, args.output_dir, config, meta)
        expression_heatmap(args.source_dir, args.output_dir, config)
        if not args.profile_only:
            comparison(args.source_dir, args.output_dir, config)


if __name__ == "__main__":
    main()
