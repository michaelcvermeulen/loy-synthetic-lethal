"""Render the main HT1080 overview and companion QC supplement.

Run prepare_ht1080_overview_qc.py first. Scientific sources are immutable.
Outputs: PDF, editable-text SVG, 600-dpi PNG and 160-dpi review PNG.
"""
from pathlib import Path
import json
import os
import sys
import argparse

THESIS = Path(__file__).resolve().parents[3]
DEPS = THESIS.parent / '.figure_python_deps'
if DEPS.is_dir():
    sys.path.insert(0, str(DEPS))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_figure1_qc'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch, Polygon, ConnectionPatch
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm, Normalize
from matplotlib import font_manager
from matplotlib.cm import ScalarMappable
from make_ht1080_screen_schematic import draw_schematic
from make_ht1080_ddpcr_passage import draw_panel as draw_ddpcr_panel

STYLE = THESIS / 'skills/thesis-figure-style'
sys.path.insert(0, str(STYLE / 'scripts'))
from thesis_style import apply_style, condition_colours
FAMILY = apply_style(base_size=9)
COLOURS = condition_colours()
BLUE, RED = COLOURS['Y-retaining'], COLOURS['Y-depleted']
GREY, PALE, DARK = '#727984', '#F3F5F7', '#222222'
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2', 'POST_LOY_REP1', 'POST_LOY_REP2']
LABELS = ['D0 CTRL', 'D0 LOY', 'CTRL R1', 'CTRL R2', 'LOY R1', 'LOY R2']
SAMPLE_COLOURS = [BLUE, RED, BLUE, BLUE, RED, RED]
SUMMARY = json.loads((DATA / 'summary.json').read_text())
REVIEW = OUT / 'review'
REVIEW.mkdir(exist_ok=True)


def letter(fig, x, y, label):
    fig.text(x, y, label, fontsize=13, fontweight='semibold', va='top')


def box(ax, x, y, w, h, text='', colour=GREY, fill='white', size=9):
    p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.006,rounding_size=0.012',
                      edgecolor=colour, facecolor=fill, linewidth=.9, clip_on=False)
    ax.add_patch(p)
    if text:
        ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=size, linespacing=1.35, color=DARK)
    return p


def arrow(ax, start, end, colour=GREY, dashed=False, style='-|>', scale=9):
    p = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=scale,
                        lw=.85, color=colour, linestyle='--' if dashed else '-',
                        shrinkA=1, shrinkB=1)
    ax.add_patch(p)


def y_ideogram(ax):
    """Horizontal reference chromosome above the heatmaps; genomic lengths preserved."""
    source = THESIS / 'figures/chapter3/source_data/ht1080_y_expression_confirmation/T2T_chrY_cytobands.csv'
    bands = pd.read_csv(source)
    bands = bands.loc[bands.chrom.eq('chrY')]
    colours = {'gneg': 'white', 'gpos50': '#999999', 'gvar': '#AFC0CD'}
    for band in bands.itertuples():
        start, end = band.start / 1e6, band.end / 1e6
        if band.stain == 'acen':
            points = ([(start, .1), (start, .9), (end, .5)] if band.band.startswith('p')
                      else [(start, .5), (end, .1), (end, .9)])
            ax.add_patch(Polygon(points, facecolor='#777777', edgecolor='#555555', lw=.4))
        else:
            ax.add_patch(Rectangle((start, .1), end-start, .8,
                                   facecolor=colours[band.stain], edgecolor='#555555', lw=.4))
    length = bands.end.max()/1e6
    centromere = bands.loc[bands.stain.eq('acen')]
    ax.set(xlim=(0, length), ylim=(0, 1))
    ax.axis('off')
    ax.text(-.025, .5, 'Y', transform=ax.transAxes, ha='right', va='center', fontsize=8.5)
    ax.text(centromere.start.min()/2e6, 1.3, 'Yp', ha='center', va='bottom', fontsize=7.5)
    ax.text((centromere.end.max()/1e6+length)/2, 1.3, 'Yq', ha='center', va='bottom', fontsize=7.5)


def workflow(fig):
    ax = fig.add_axes([.045, .495, .91, .465])
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    letter(fig, .016, .975, 'a')
    # Stylised vector 4D instrument and two established populations.
    ax.text(.08, .975, 'HT1080', ha='center', va='top', fontweight='semibold')
    ax.add_patch(Rectangle((.037, .83), .084, .065, facecolor='#D9E2EB', edgecolor=BLUE, lw=.8))
    ax.add_patch(Rectangle((.050, .839), .058, .045, facecolor='#F4F8FB', edgecolor=BLUE, lw=.6))
    for y in [.785, .74]:
        box(ax, .021, y, .12, .040, colour=GREY, fill='#DCE2E8')
        ax.add_patch(Rectangle((.028, y+.005), .084, .028, facecolor=BLUE, edgecolor='none'))
    ax.text(.08, .725, 'Lonza 4D-Nucleofector\nY-targeting / control RNP', ha='center', va='top', fontsize=8.5)
    arrow(ax, (.15, .81), (.21, .81))
    ax.text(.30, .975, 'Expansion', ha='center', va='top', fontweight='semibold')
    box(ax, .22, .83, .16, .075, 'CTRL', BLUE, '#F0F3F8')
    box(ax, .22, .70, .16, .075, 'Y-depleted', RED, '#FCF1EC')
    arrow(ax, (.39, .86), (.45, .82), BLUE)
    arrow(ax, (.39, .74), (.45, .78), RED)
    ax.text(.54, .975, 'Library delivery', ha='center', va='top', fontweight='semibold')
    box(ax, .45, .69, .19, .22, 'Brunello transduction\n77,441 sgRNAs\nPuromycin selection', size=8.5)
    ax.text(.54, .655, 'Each genotype separately', ha='center', va='top', fontsize=8, color=GREY)
    # Separate D0 sampling and post-selection outgrowth replication.
    arrow(ax, (.64, .805), (.71, .805))
    ax.text(.84, .975, 'Replicate outgrowth', ha='center', va='top', fontweight='semibold')
    for y, colour, label in [(.855, BLUE, 'CTRL'), (.725, RED, 'LOY')]:
        ax.text(.70, y, label, ha='right', va='center', color=colour, fontweight='semibold', fontsize=8.5)
        for x, rep in [(.77, 'R1'), (.91, 'R2')]:
            box(ax, x-.046, y-.037, .092, .074, rep, colour, '#F8F9FB', 8.5)
    ax.text(.84, .65, '2 outgrowth replicates / genotype', ha='center', va='top', fontsize=8)
    arrow(ax, (.54, .69), (.54, .515), dashed=True)
    box(ax, .49, .435, .10, .075, 'D0', GREY, PALE, 8.5)
    ax.text(.54, .41, '1 sample\nper genotype', ha='center', va='top', fontsize=8)
    # Shared timeline denotes parallel passaging, not pooling.
    arrow(ax, (.95, .60), (.95, .515), colour=GREY)
    ax.text(.81, .56, 'Four passages', ha='center', va='bottom', fontsize=9)
    for x, p in zip([.66, .75, .84, .93], ['P1', 'P2', 'P3', 'P4']):
        ax.plot(x, .48, 'o', color=GREY, markersize=3.5)
        ax.text(x, .445, p, ha='center', va='top', fontsize=8)
    arrow(ax, (.66, .48), (.95, .48))
    ax.text(.81, .30, '32 million cells / culture / passage\n9.3–9.8 cumulative doublings\n4 endpoint samples', ha='center', va='center', fontsize=8.5, linespacing=1.4)
    # Reference Y chromosome replaces the former molecular-confirmation heading.
    ideogram_ax = ax.inset_axes([.055, .565, .395, .028])
    y_ideogram(ideogram_ax)
    # Inset heatmaps remain separate measurements and scales.
    dna = pd.read_csv(DATA / 'dna_inset.csv').pivot(index='condition', columns='region', values='median_relative_cn')
    d = dna.loc[['sgCTRL', 'LOY'], ['MSY Yp', 'MSY Yq (euchromatic)']].to_numpy()
    r = pd.read_csv(DATA / 'rna_inset.csv').set_index('gene_name')
    dna_ax = ax.inset_axes([.055, .325, .105, .16])
    di = dna_ax.imshow(d, cmap='viridis', vmin=0, vmax=.55, aspect='auto')
    dna_ax.set_xticks([0, 1], ['Yp', 'Yq'], fontsize=8)
    dna_ax.set_yticks([0, 1], ['CTRL', 'LOY'], fontsize=8)
    for lab, colour in zip(dna_ax.get_yticklabels(), [BLUE, RED]): lab.set_color(colour)
    dna_ax.tick_params(length=0, pad=3)
    dna_ax.set_title('DNA', fontsize=8.5, pad=4)
    for i in range(2):
        for j in range(2):
            value = f'{d[i,j]:.3f}' if d[i,j] < .00101 else f'{d[i,j]:.2f}'
            dna_ax.text(j, i, value, ha='center', va='center', fontsize=8, color='white' if d[i,j] < .3 else DARK)
    rna_ax = ax.inset_axes([.205, .325, .245, .16])
    rna_values = r[['sgCTRL', 'LOY']].to_numpy().T
    ri = rna_ax.imshow(rna_values, cmap='magma', vmin=0, vmax=100, aspect='auto')
    for i, row in enumerate(rna_values):
        for j, value in enumerate(row):
            rna_ax.add_patch(Rectangle((j-.5, i-.5), 1, 1, fill=False,
                                      edgecolor='white', linewidth=.8))
            label = f'{value:.0f}' if i == 0 else f'{value:.1f}'
            rna_ax.text(j, i, label, ha='center', va='center', fontsize=8,
                        color=DARK if value >= 65 else 'white')
    rna_ax.set_xticks(range(len(r)), r.index, rotation=55, ha='right', fontsize=8, fontstyle='italic')
    rna_ax.set_yticks([]); rna_ax.tick_params(length=0, pad=2)
    rna_ax.set_title('RNA', fontsize=8.5, pad=4,
                     bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1})
    # Map each equally spaced RNA column to its actual T2T-Y gene midpoint.
    gene_source = THESIS / 'figures/chapter3/source_data/ht1080_y_expression_confirmation/gene_mapping_plotted.csv'
    genes = pd.read_csv(gene_source).set_index('gene_name').loc[r.index]
    for j, gene in enumerate(genes.itertuples()):
        midpoint_mb = (gene.t2t_start + gene.t2t_end) / 2e6
        ax.add_artist(ConnectionPatch(
            xyA=(midpoint_mb, .1), coordsA=ideogram_ax.transData,
            xyB=(j, -.5), coordsB=rna_ax.transData,
            color='#808080', alpha=.2, linewidth=.7,
            zorder=1, clip_on=False))
    cb1 = fig.colorbar(di, cax=ax.inset_axes([.055, .165, .105, .018]), orientation='horizontal', ticks=[0,.5])
    cb1.ax.tick_params(labelsize=7, length=2, pad=1)
    cb1.set_label('Relative CN', fontsize=7.5, labelpad=1)
    cb2 = fig.colorbar(ri, cax=ax.inset_axes([.27, .165, .18, .018]), orientation='horizontal', ticks=[0,50,100])
    cb2.ax.tick_params(labelsize=7, length=2, pad=1)
    cb2.set_label('% of CTRL expression', fontsize=7.5, labelpad=1)
    for a in [dna_ax, rna_ax]:
        for spine in a.spines.values(): spine.set_visible(False)
    # Sample processing sequence, followed by the endpoint-only contrast.
    ax.text(.006, .036, 'Each D0 /\nendpoint sample:', fontsize=8.5, va='center')
    seq = [(.23, .09, 'gDNA'), (.40, .20, 'sgRNA PCR + indexing'), (.60, .13, 'Sequencing'), (.77, .13, 'Guide counts'), (.93, .12, 'drugZ')]
    for x, w, text in seq:
        box(ax, x-w/2, .003, w, .065, text, GREY, PALE, 8.5)
    for a, b in zip(seq[:-1], seq[1:]): arrow(ax, (a[0]+a[1]/2+.001, .0355), (b[0]-b[1]/2-.001, .0355))
    ax.text(.93, -.03, 'Endpoints: LOY vs CTRL', ha='center', va='top', fontsize=8)
    return ax


def molecular_confirmation(fig, rect):
    """Preserve the Chapter 3 DNA/RNA inset alongside the new screen schematic."""
    ax=fig.add_axes(rect, xlim=(0,1), ylim=(0,1))
    ax.axis('off')
    # Reference Y chromosome replaces the former molecular-confirmation heading.
    ideogram_ax = ax.inset_axes([.02, .92, .96, .045])
    y_ideogram(ideogram_ax)
    # Inset heatmaps remain separate measurements and scales.
    dna = pd.read_csv(DATA / 'dna_inset.csv').pivot(index='condition', columns='region', values='median_relative_cn')
    d = dna.loc[['sgCTRL', 'LOY'], ['MSY Yp', 'MSY Yq (euchromatic)']].to_numpy()
    r = pd.read_csv(DATA / 'rna_inset.csv').set_index('gene_name')
    dna_ax = ax.inset_axes([.02, .57, .23, .20])
    di = dna_ax.imshow(d, cmap='viridis', vmin=0, vmax=.55, aspect='auto')
    dna_ax.set_xticks([0, 1], ['Yp', 'Yq'], fontsize=8)
    dna_ax.set_yticks([0, 1], ['CTRL', 'LOY'], fontsize=8)
    for lab, colour in zip(dna_ax.get_yticklabels(), [BLUE, RED]): lab.set_color(colour)
    dna_ax.tick_params(length=0, pad=3)
    dna_ax.set_title('DNA', fontsize=8.5, pad=4)
    for i in range(2):
        for j in range(2):
            value = f'{d[i,j]:.3f}' if d[i,j] < .00101 else f'{d[i,j]:.2f}'
            dna_ax.text(j, i, value, ha='center', va='center', fontsize=8, color='white' if d[i,j] < .3 else DARK)
    rna_ax = ax.inset_axes([.35, .57, .63, .20])
    rna_values = r[['sgCTRL', 'LOY']].to_numpy().T
    ri = rna_ax.imshow(rna_values, cmap='magma', vmin=0, vmax=100, aspect='auto')
    for i, row in enumerate(rna_values):
        for j, value in enumerate(row):
            rna_ax.add_patch(Rectangle((j-.5, i-.5), 1, 1, fill=False,
                                      edgecolor='white', linewidth=.8))
            label = f'{value:.0f}' if i == 0 else f'{value:.1f}'
            rna_ax.text(j, i, label, ha='center', va='center', fontsize=8,
                        color=DARK if value >= 65 else 'white')
    rna_ax.set_xticks(range(len(r)), r.index, rotation=55, ha='right', fontsize=8, fontstyle='italic')
    rna_ax.set_yticks([]); rna_ax.tick_params(length=0, pad=2)
    rna_ax.set_title('RNA', fontsize=8.5, pad=4,
                     bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1})
    # Map each equally spaced RNA column to its actual T2T-Y gene midpoint.
    gene_source = THESIS / 'figures/chapter3/source_data/ht1080_y_expression_confirmation/gene_mapping_plotted.csv'
    genes = pd.read_csv(gene_source).set_index('gene_name').loc[r.index]
    for j, gene in enumerate(genes.itertuples()):
        midpoint_mb = (gene.t2t_start + gene.t2t_end) / 2e6
        ax.add_artist(ConnectionPatch(
            xyA=(midpoint_mb, .1), coordsA=ideogram_ax.transData,
            xyB=(j, -.5), coordsB=rna_ax.transData,
            color='#808080', alpha=.2, linewidth=.7,
            zorder=1, clip_on=False))
    cb1 = fig.colorbar(di, cax=ax.inset_axes([.02, .19, .23, .025]), orientation='horizontal', ticks=[0,.5])
    cb1.ax.tick_params(labelsize=7, length=2, pad=1)
    cb1.set_label('Relative CN', fontsize=7.5, labelpad=1)
    cb2 = fig.colorbar(ri, cax=ax.inset_axes([.52, .19, .46, .025]), orientation='horizontal', ticks=[0,50,100])
    cb2.ax.tick_params(labelsize=7, length=2, pad=1)
    cb2.set_label('% of CTRL expression', fontsize=7.5, labelpad=1)
    for a in [dna_ax, rna_ax]:
        for spine in a.spines.values(): spine.set_visible(False)
    return ax


def plot_doublings(ax):
    d = pd.read_csv(DATA / 'population_doublings.csv')
    for group, colour in [('CTRL', BLUE), ('LOY', RED)]:
        for rep, marker, ls in [('Rep1', 'o', '-'), ('Rep2', '^', '--')]:
            v = d.loc[d.Group.eq(group) & d.Replicate.eq(rep)].sort_values('Passage')
            ax.plot(range(1,5), v['Cumulative Doublings'], color=colour, marker=marker,
                    ls=ls, ms=4.0, lw=1.4, markerfacecolor='white', label=f'{group} R{rep[-1]}')
    ax.set(xlim=(.85,4.15), ylim=(0,10.4), xlabel='Screening passage', ylabel='Cumulative population doublings')
    ax.set_xticks([1,2,3,4], ['P1','P2','P3','P4'])
    ax.set_yticks([0,2,4,6,8,10])
    ax.legend(loc='upper left', fontsize=8, handlelength=1.8, borderaxespad=.1)


def plot_roc(ax):
    for group, colour in [('CTRL', BLUE), ('LOY', RED)]:
        d = pd.read_csv(DATA / f'roc_curve_{group}.csv')
        ax.plot(d.false_positive_rate, d.true_positive_rate, color=colour, lw=1.6,
                label=f'{group}  AUC {SUMMARY["roc"][group]["auc"]:.3f}')
    ax.plot([0,1], [0,1], color='#B4B8BF', lw=.8, ls='--', zorder=0)
    ax.set(xlim=(0,1), ylim=(0,1.02), xlabel='False-positive rate', ylabel='True-positive rate')
    ax.set_xticks([0,.5,1]); ax.set_yticks([0,.5,1])
    ax.legend(loc='lower right', fontsize=8.5, handlelength=1.5)
    ax.text(.98, .37, 'Essential / non-essential genes\nCTRL: 989 / 740\nLOY: 994 / 736',
            ha='right', va='bottom', transform=ax.transAxes, fontsize=8, linespacing=1.4)


def plot_matrix(ax, fig):
    d = pd.read_csv(DATA / 'spearman_guide_counts.csv', index_col=0).loc[SAMPLES, SAMPLES]
    norm = Normalize(vmin=.83, vmax=1)
    im = ax.imshow(d, cmap='magma', norm=norm, aspect='equal')
    for i in range(6):
        for j in range(6):
            ax.text(j,i,f'{d.iloc[i,j]:.3f}', ha='center', va='center', fontsize=7.7,
                    color=DARK if d.iloc[i,j] > .94 else 'white')
    ax.set_xticks(range(6), LABELS, rotation=55, ha='right', fontsize=8)
    ax.set_yticks(range(6), LABELS, fontsize=8)
    for labs in [ax.get_xticklabels(),ax.get_yticklabels()]:
        for lab,c in zip(labs,SAMPLE_COLOURS): lab.set_color(c)
    ax.tick_params(length=0, pad=3)
    ax.axhline(1.5, color='white', lw=1); ax.axvline(1.5, color='white', lw=1)
    for sp in ax.spines.values(): sp.set_visible(False)
    cb = fig.colorbar(im, cax=ax.inset_axes([0,1.075,1,.045]), orientation='horizontal', ticks=[.85,.9,.95,1])
    cb.ax.tick_params(labelsize=8, length=2, pad=2)
    cb.ax.xaxis.set_ticks_position('top'); cb.ax.xaxis.set_label_position('top')
    cb.set_label('Spearman correlation', fontsize=9, labelpad=4)


def main_figure():
    fig = plt.figure(figsize=(9.6,234/25.4))
    width_mm, height_mm = 9.6*25.4, 234
    draw_schematic(fig, rect=[10/width_mm,151/height_mm,224/width_mm,(224*64/183)/height_mm],
                   include_note=False, compact=True)
    letter(fig,.016,.985,'a')
    molecular_confirmation(fig,[14/width_mm,148/height_mm,80/width_mm,35/height_mm])
    # Retain the compact a-d arrangement and add two readable QC panels below.
    row_y=(.09*203.2+60)/height_mm
    row_h=.285*203.2/height_mm
    a1 = fig.add_axes([.075,row_y,.235,row_h]); plot_doublings(a1)
    a2 = fig.add_axes([.395,row_y,.235,row_h]); plot_depletion(a2)
    a3 = fig.add_axes([.744,row_y,.225,row_h]); plot_matrix(a3,fig)
    for x,l in [(.02,'b'),(.34,'c'),(.68,'d')]: letter(fig,x,(.416*203.2+60)/height_mm,l)
    plot_pca(fig.add_axes([.075,14/height_mm,.235,44/height_mm]))
    plot_reference_densities(fig.add_axes([.395,14/height_mm,.574,44/height_mm]))
    for x,l in [(.02,'e'),(.34,'f')]: letter(fig,x,62/height_mm,l)
    return fig


def plot_mapping(fig):
    d = pd.read_csv(DATA / 'sequencing_qc.csv', index_col=0).loc[SAMPLES]
    positions = np.arange(6)
    a1 = fig.add_axes([.085,.82,.225,.11])
    a2 = fig.add_axes([.085,.72,.225,.075])
    for ax,col in [(a1,'percent_mapped'),(a2,'percent_detected')]:
        if col == 'percent_mapped':
            ax.bar(positions,d[col],color=SAMPLE_COLOURS,width=.65,alpha=.9)
        else:
            ax.scatter(positions,d[col],c=SAMPLE_COLOURS,s=26,zorder=3)
        ax.axvline(1.5,color='#BCC2C9',lw=.8,ls='--')
        ax.set_xlim(-.6,5.6)
    a1.set_ylim(0,102); a1.set_yticks([0,50,100]); a1.set_ylabel('Mapped reads (%)')
    a1.set_xticks([])
    for x,v in enumerate(d.percent_mapped): a1.text(x,v+2,f'{v:.1f}',ha='center',va='bottom',fontsize=8)
    a2.set_ylim(98.8,100); a2.set_yticks([98.8,99.4,100]); a2.set_ylabel('Guides\ndetected (%)',fontsize=8.5,labelpad=0)
    a2.set_xticks(positions,LABELS,rotation=45,ha='right',fontsize=8)
    for x,v in enumerate(d.percent_detected): a2.text(x,v+.035,f'{v:.2f}',ha='center',va='bottom',fontsize=7.5)
    # Dots make the explicitly narrowed guide-detection scale unambiguous.
    return a1,a2


def plot_counts(ax, log):
    values = [log[s].to_numpy() for s in SAMPLES]
    viol = ax.violinplot(values, positions=range(6), widths=.78, showmedians=True, showextrema=False, points=120)
    for body,c in zip(viol['bodies'],SAMPLE_COLOURS): body.set(facecolor=c, edgecolor=c, alpha=.58, linewidth=.7)
    viol['cmedians'].set(color=DARK,linewidth=1)
    ax.set_xticks(range(6),LABELS,rotation=45,ha='right',fontsize=8)
    ax.set(ylabel='Guide abundance, log₂(CPM + 1)',ylim=(-.25,16),xlim=(-.6,5.6))
    ax.set_yticks([0,4,8,12,16])
    ax.axvline(1.5,color='#BCC2C9',lw=.8,ls='--')
    zeros = log[SAMPLES].eq(0).sum()
    for x,s in enumerate(SAMPLES): ax.text(x,1.015,str(zeros[s]),ha='center',fontsize=8,transform=ax.get_xaxis_transform())
    ax.text(.02,1.08,'Zero-count guides',transform=ax.transAxes,va='bottom',fontsize=8)
    ax.text(.98,.93,'77,441 guides / sample',transform=ax.transAxes,ha='right',fontsize=8)


def plot_scatter(ax, log, group, norm):
    s=SUMMARY['regression'][group]
    x,y=[log[s[col]].to_numpy() for col in ['x_sample','y_sample']]
    is_d0=group=='D0'
    limit=15 if is_d0 else 8.25
    assert max(x.max(),y.max()) < limit
    hb=ax.hexbin(x,y,gridsize=95,extent=(0,limit,0,limit),mincnt=1,cmap='magma',norm=norm,linewidths=0,rasterized=True)
    ax.plot([0,limit],[0,limit],color='#9AA1A8',ls='--',lw=.9,label='Identity')
    colour=DARK if is_d0 else BLUE if group=='CTRL' else RED
    ax.plot([0,limit],[s['intercept'],s['intercept']+limit*s['slope']],color=colour,lw=1.35,label='Linear fit')
    x_label='D0 CTRL, log₂(CPM + 1)' if is_d0 else 'R1, log₂(CPM + 1)'
    y_label='D0 LOY, log₂(CPM + 1)' if is_d0 else 'R2, log₂(CPM + 1)'
    ax.set(xlim=(0,limit),ylim=(0,limit),xlabel=x_label,ylabel=y_label)
    ticks=[0,5,10,15] if is_d0 else [0,2,4,6,8]
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    ax.text(.025,.975,'D0: CTRL vs LOY' if is_d0 else group,transform=ax.transAxes,va='top',fontweight='semibold',color=colour)
    ax.text(.025,.87,f'Pearson r = {s["pearson_r"]:.3f}\nn = 77,441 guides',transform=ax.transAxes,va='top',fontsize=8.5)
    ax.legend(loc='lower right',fontsize=8,handlelength=2)
    return hb


def plot_gini(ax):
    d=pd.read_csv(DATA/'gini_per_sample.csv',index_col=0).loc[SAMPLES]
    ax.bar(range(6),d.GiniIndex,color=SAMPLE_COLOURS,width=.65,alpha=.9)
    ax.set(xlim=(-.6,5.6),ylim=(0,.14),ylabel='MAGeCK Gini index (log counts)')
    ax.set_xticks(range(6),LABELS,rotation=45,ha='right',fontsize=8)
    ax.set_yticks([0,.04,.08,.12])
    ax.axvline(1.5,color='#BCC2C9',lw=.8,ls='--')
    for x,value in enumerate(d.GiniIndex):
        ax.text(x,value+.003,f'{value:.5f}',ha='center',va='bottom',fontsize=8,rotation=45)


def plot_depletion(ax):
    d=pd.read_csv(DATA/'reference_gene_depletion.csv')
    combos=[('CTRL','Essential'),('CTRL','Non-Essential'),('LOY','Essential'),('LOY','Non-Essential')]
    values=[d.loc[d.category.eq(cat),f'LFC_{g}'].dropna().to_numpy() for g,cat in combos]
    viol=ax.violinplot(values,positions=[0,1,2.4,3.4],widths=.82,showmedians=True,showextrema=False,points=120)
    for body,c in zip(viol['bodies'],[BLUE,BLUE,RED,RED]): body.set(facecolor=c,edgecolor=c,alpha=.58,lw=.7)
    viol['cmedians'].set(color=DARK,linewidth=1.1)
    ax.axhline(0,color='#9AA1A8',ls='--',lw=.8)
    ax.set(xlim=(-.65,4.05),ylabel='Gene log₂(endpoint / D0)')
    lo=min(map(np.min,values)); hi=max(map(np.max,values))
    ax.set_ylim(np.floor(lo)-.3,np.ceil(hi)+1.4)
    ax.set_xticks([0,1,2.4,3.4],['Ess.','Non-ess.','Ess.','Non-ess.'],fontsize=7.5)
    for x,v in zip([0,1,2.4,3.4],values): ax.text(x,np.ceil(hi)+.2,f'{len(v)}',ha='center',fontsize=7.5)
    ax.text(.02,.98,'n genes',transform=ax.transAxes,va='top',fontsize=7)
    ax.text(.5,-.2,'CTRL',transform=ax.get_xaxis_transform(),ha='center',color=BLUE,fontweight='semibold')
    ax.text(2.9,-.2,'LOY',transform=ax.get_xaxis_transform(),ha='center',color=RED,fontweight='semibold')


def plot_ntc(ax):
    d=pd.read_csv(DATA/'ntc_guide_endpoint_lfc.csv').LFC
    s=SUMMARY['ntc']
    edges=np.arange(np.floor(d.min()*5)/5-.1,np.ceil(d.max()*5)/5+.21,.10)
    ax.axvspan(s['q025'],s['q975'],color='#E8EDF3',zorder=0)
    ax.hist(d,bins=edges,color='#858E9D',edgecolor='white',linewidth=.35,zorder=2)
    ax.axvline(0,color=DARK,ls='--',lw=.9,zorder=3)
    ax.set(xlabel='NTC guide log₂ fold change\n(LOY / CTRL)',ylabel='Number of guides')
    ax.set_xlim(min(edges),max(edges))
    ax.text(.98,.96,f'n = {s["n"]:,}\nMedian = {s["median"]:.3f}\nSD = {s["sd"]:.3f}\nCentral 95%: {s["q025"]:.2f} to {s["q975"]:.2f}',
            ha='right',va='top',transform=ax.transAxes,fontsize=7.2,linespacing=1.3)


def plot_reference_densities(ax):
    table=pd.read_csv(DATA/'endpoint_reference_density_guides.csv')
    curves=pd.read_csv(DATA/'endpoint_reference_density_curves.csv')
    meta=json.loads((DATA/'endpoint_reference_density_manifest.json').read_text())
    ntc=table.loc[table.included & table.category.eq('Non-targeting'),'LFC']
    # Low-opacity reference distributions sit behind the NTC histogram and curve.
    for category,colour,ls in [('Essential','#8473AC','-'),('Non-essential','#B79549','--')]:
        ax.fill_between(curves.LFC,curves[category],color=colour,alpha=.14,zorder=1,linewidth=0)
        ax.plot(curves.LFC,curves[category],color=colour,alpha=.55,lw=1.1,ls=ls,zorder=2,
                label=f'{category} guides (n = {meta["groups"][category]["included_guides"]:,})')
    edges=np.arange(np.floor(ntc.min()*10)/10,np.ceil(ntc.max()*10)/10+.1001,.1)
    ax.hist(ntc,bins=edges,density=True,color='#858E9D',alpha=.22,edgecolor='white',linewidth=.25,zorder=3)
    ax.plot(curves.LFC,curves['Non-targeting'],color=DARK,lw=1.5,zorder=4,
            label=f'Non-targeting guides (n = {len(ntc):,})')
    ax.axvline(0,color='#9AA1A8',ls=':',lw=.7,zorder=0)
    ax.set(xlabel='Guide log₂ fold change (LOY / CTRL endpoints)',ylabel='Density',
           xlim=(curves.LFC.min(),curves.LFC.max()),ylim=(0,curves['Non-targeting'].max()*1.17))
    ax.set_xticks([-8,-4,0,4,8])
    handles,labels=ax.get_legend_handles_labels()
    ax.legend([handles[i] for i in [2,0,1]],[labels[i] for i in [2,0,1]],loc='upper left',
              fontsize=8,handlelength=1.8,labelspacing=.3)
    s=meta['groups']['Non-targeting']
    ax.text(.985,.96,f'NTC median = {s["median"]:.3f}\nSD = {s["sd"]:.3f}\nCentral 95%: {s["q025"]:.2f} to {s["q975"]:.2f}',
            ha='right',va='top',transform=ax.transAxes,fontsize=7.5,linespacing=1.15)


def plot_pca(ax):
    d=pd.read_csv(DATA/'pca_sample_scores.csv').set_index('sample').loc[SAMPLES]
    pct=pd.read_csv(DATA/'pca_variance.csv').variance_percent.to_numpy()
    # Replicate scores are close; separate labels vertically with short leaders.
    offsets={'D0_CTRL':(-3,10),'D0_LOY':(-3,-14),
             'POST_CTRL_REP1':(5,13),'POST_CTRL_REP2':(5,-13),
             'POST_LOY_REP1':(5,13),'POST_LOY_REP2':(5,-13)}
    for sample,label,colour in zip(SAMPLES,LABELS,SAMPLE_COLOURS):
        row=d.loc[sample]
        ax.scatter(row.PC1,row.PC2,s=30,marker='D' if row.stage=='D0' else 'o',
                   color=colour,edgecolor='white',linewidth=.5,zorder=3)
        ax.annotate(label,(row.PC1,row.PC2),xytext=offsets[sample],textcoords='offset points',
                    fontsize=7,color=colour,ha='right' if sample.startswith('D0') else 'left',
                    va='center',arrowprops={'arrowstyle':'-','color':colour,'lw':.5,'alpha':.55})
    ax.set(xlabel=f'PC1 ({pct[0]:.1f}%)',ylabel=f'PC2 ({pct[1]:.1f}%)',xlim=(-245,425),ylim=(-130,130))
    ax.set_xticks([-200,0,200,400]); ax.set_yticks([-100,0,100])
    ax.tick_params(labelsize=8)


def supplementary_figure():
    fig=plt.figure(figsize=(9.6,9.8))
    log=pd.read_csv(DATA/'guide_log2_cpm.csv')
    plot_mapping(fig)
    plot_counts(fig.add_axes([.405,.72,.235,.21]),log)
    plot_gini(fig.add_axes([.735,.72,.235,.21]))
    d=fig.add_axes([.075,.355,.235,.22])
    d0_norm=LogNorm(vmin=1,vmax=10000)
    d0=plot_scatter(d,log,'D0',d0_norm)
    # Retain the complete original scatter and add the post-exclusion comparison.
    d.get_legend().remove()
    d.legend(loc='upper left',fontsize=7.5,handlelength=1.5,
             bbox_to_anchor=(.02,.69),labelspacing=.3,borderpad=.1)
    filtered_data=OUT/'source_data/ht1080_d0_jackpot_removed'
    filtered=pd.read_csv(filtered_data/'retained_guide_coordinates.csv.gz')
    stats=pd.read_csv(filtered_data/'correlation_summary.csv').set_index('analysis').loc['exclude_QC173']
    x=filtered.D0_CTRL_log2_CPM_plus1.to_numpy()
    y=filtered.D0_LOY_log2_CPM_plus1.to_numpy()
    assert len(x)==77268 and np.isclose(np.corrcoef(x,y)[0,1],stats.pearson_recalculated_CPM)
    inset=d.inset_axes([.58,.17,.40,.34],zorder=5)
    inset.set_facecolor('white')
    inset.hexbin(x,y,gridsize=95,extent=(0,15,0,15),mincnt=1,cmap='magma',
                 norm=d0_norm,linewidths=0,rasterized=True)
    assert max(x.max(),y.max())<15
    inset.plot([0,15],[0,15],color='#9AA1A8',ls='--',lw=.6)
    inset.plot([0,15],[stats.intercept,stats.intercept+15*stats.slope],color=DARK,lw=.8)
    inset.set(xlim=(0,8.25),ylim=(0,8.25),xticks=[0,2,4,6,8],yticks=[0,2,4,6,8])
    inset.tick_params(labelsize=6,length=2,pad=1)
    inset.set_xlabel('CTRL',fontsize=6.5,labelpad=0)
    inset.set_ylabel('LOY',fontsize=6.5,labelpad=0)
    inset.text(0,1.05,f'173 removed; r = {stats.pearson_recalculated_CPM:.3f}',
               transform=inset.transAxes,va='bottom',fontsize=6,clip_on=False)
    assert d0.get_array().max()<=d0_norm.vmax
    d0cb=fig.colorbar(d0,cax=fig.add_axes([.075,.615,.235,.009]),orientation='horizontal',ticks=[1,100,10000])
    d0cb.set_label('D0 guides / hexagon',fontsize=8,labelpad=4)
    d0cb.ax.xaxis.set_label_position('top')
    d0cb.ax.tick_params(labelsize=7.5,length=2,pad=2)
    norm=LogNorm(vmin=1,vmax=3000)
    e=fig.add_axes([.405,.355,.235,.22]); f=fig.add_axes([.735,.355,.235,.22])
    h1=plot_scatter(e,log,'CTRL',norm); h2=plot_scatter(f,log,'LOY',norm)
    assert max(h1.get_array().max(),h2.get_array().max()) <= norm.vmax
    cb=fig.colorbar(h2,cax=fig.add_axes([.405,.615,.565,.009]),orientation='horizontal',ticks=[1,10,100,1000])
    cb.set_label('Endpoint guides / hexagon',fontsize=8,labelpad=4)
    cb.ax.xaxis.set_label_position('top')
    cb.ax.tick_params(labelsize=7.5,length=2,pad=2)
    plot_roc(fig.add_axes([.075,.075,.235,.18]))
    # Reuse the frozen passage estimates; both condition facets form panel h.
    ddpcr=OUT/'source_data/ht1080_ddpcr_passage'
    wells=pd.read_csv(ddpcr/'plotted_wells.csv')
    summary=pd.read_csv(ddpcr/'plotted_summary.csv')
    left=fig.add_axes([.425,.075,.2625,.18])
    right=fig.add_axes([.7075,.075,.2625,.18],sharey=left)
    draw_ddpcr_panel(left,wells,summary,'CTRL')
    draw_ddpcr_panel(right,wells,summary,'LOY',show_y=False)
    fig.text(.6975,.022,'Screening passage',ha='center',va='bottom',fontsize=9)
    assert left.get_ylim()==right.get_ylim()
    for x,y,l in [(.02,.965,'a'),(.35,.965,'b'),(.68,.965,'c'),(.02,.60,'d'),(.35,.60,'e'),(.68,.60,'f'),(.02,.278,'g'),(.35,.278,'h')]: letter(fig,x,y,l)
    return fig


def save(fig,stem):
    for fmt in ['pdf','svg','png']:
        fig.savefig(OUT/f'{stem}.{fmt}',dpi=600,metadata={'Creator':'Thesis figure draft; source data and code supplied'} if fmt=='pdf' else None)
    fig.savefig(REVIEW/f'{stem}_preview.png',dpi=160)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    scope=parser.add_mutually_exclusive_group()
    scope.add_argument('--supp-only',action='store_true')
    scope.add_argument('--main-only',action='store_true')
    args=parser.parse_args()
    if not args.supp_only:
        save(main_figure(),'ch04_main_ht1080_crispr_screen_workflow')
    if not args.main_only:
        save(supplementary_figure(),'ch04_supp_ht1080_crispr_screen_qc')
    version={'font':FAMILY,'regular_font':font_manager.findfont(font_manager.FontProperties(family=FAMILY)),
             'matplotlib':matplotlib.__version__,'numpy':np.__version__,'pandas':pd.__version__,
             'main_size_inches':[9.6,234/25.4],'supp_size_inches':[9.6,9.8],'supp_panel_rows':[3,3,2],'png_dpi':600,
             'scope':'Main c shows essential/non-essential depletion; e shows saved PCA; f shows matched endpoint guide-LFC densities. Supplementary d includes a D0 inset after removal of the 173 QC-excluded guides (all 164 jackpot guides plus nine further D0 outliers) and CPM recalculation; g retains ROC; h contains the two frozen ddPCR passage facets. Hit calling unchanged.'}
    (DATA/'render_manifest.json').write_text(json.dumps(version,indent=2),encoding='utf-8')
    print(json.dumps(version,indent=2))


if __name__=='__main__': main()
