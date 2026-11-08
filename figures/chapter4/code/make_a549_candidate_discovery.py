"""Render an A549 review draft in the HT1080 Figure 2 visual format."""
from pathlib import Path
import sys
import os
import json
import hashlib
import itertools
import math
THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/a549_candidate_discovery'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

STEM = 'ch04_main_a549_candidate_discovery_draft'
W, H = 9.6, 8.35
STRONG, EXPLORATORY, NS = '#7A3E85', '#BC8A26', '#B8B8B8'
POSITIVE, NEGATIVE = '#629EBB', '#C97769'
AREA = 15.


def fdr_text(x):
    return f'{x:.3g}' if x < .01 else f'{x:.3f}'


def axes(fig, x, y, w, h):
    return fig.add_axes([x/W, y/H, w/W, h/H])


def text(fig, x, y, value, **kw):
    return fig.text(x/W, y/H, value, **kw)


def letter(fig, x, y, value):
    return text(fig, x, y, value, fontsize=13, fontweight='semibold', va='top')


def guide_panel(ax, genes, guides, density, marks, panel, slots=None):
    slots = slots or len(genes)
    grid, curve = density.guide_log2FC, density.display_density_peak_scaled
    xlim = (grid.min(), grid.max())
    for j, gene in enumerate(genes):
        base = slots-1-j
        sub = guides.loc[guides.Gene.eq(gene)]
        ax.fill_between(grid, base+.06, base+.06+curve*.76, color='#EAEAEA', lw=0)
        ax.plot(grid, base+.06+curve*.76, color='#A7A7A7', lw=.6)
        ax.hlines(base+.06, *xlim, color='#DEDEDE', lw=.4)
        ax.vlines(0, base+.06, base+.86, color='#C3C3C3', lw=.45, ls=(0, (2, 2)))
        for r in sub.loc[sub.lfc_displayed].itertuples():
            ax.vlines(r.mageck_reported_log2FC, base+.06, base+.86, color='#333333', lw=.85)
            marks.append(dict(panel=panel, gene=gene, sgRNA=r.sgRNA,
                              x=r.mageck_reported_log2FC, y=base+.06, height=.8))
        label = gene + ('*' if not sub.lfc_displayed.all() else '')
        ax.text(-.045, base+.47, label, transform=ax.get_yaxis_transform(),
                va='center', ha='right', fontsize=7.5)
        ax.text(1.045, base+.47, fdr_text(sub.drugz_directional_fdr.iloc[0]),
                transform=ax.get_yaxis_transform(), va='center', fontsize=7.3)
    ax.text(1.045, 1.018, 'drugZ FDR', transform=ax.transAxes, fontsize=6.8, va='bottom')
    ax.set(xlim=xlim, ylim=(-.08, slots+.02), yticks=[])
    ax.spines['left'].set_visible(False)
    ax.set_xticks([-8, -4, 0, 4])
    ax.tick_params(axis='x', labelsize=7, length=2, pad=2)
    ax.set_xlabel('Guide log₂FC (LOY / CTRL)', fontsize=7.2, labelpad=3)


def layout(nodes, pairs):
    core = nodes.index[nodes.depmap_degree.add(nodes.string_degree).gt(0)].tolist()
    strip = nodes.index[nodes.depmap_degree.add(nodes.string_degree).eq(0)].tolist()
    graph = nx.Graph()
    graph.add_nodes_from(core)
    for e in pairs.loc[pairs.depmap_edge | pairs.string_edge].itertuples():
        wt = (4*e.pearson_r if e.depmap_edge and e.pearson_r > 0 else .15 if e.depmap_edge else 0)
        wt += .35*e.string_score if e.string_edge else 0
        graph.add_edge(e.gene_a, e.gene_b, weight=wt)
    # Pack disconnected components into compact rows. Shared coordinates show
    # which source supplies each edge; distances have no quantitative meaning.
    pos = dict(zip(['MRPL44','TRA2A','FGF12','MYH9','CTSH','PTPRN'],
                   [(x,.92) for x in [.02,.21,.41,.61,.80,.99]]))
    pos.update({'UBR5':(.21,.70),'CHGA':(.99,.70)})
    for centre, middle, left, right in [
        (.15,'EIF1AX','AK9','DDX42'), (.49,'PNPLA4','PLA2G15','WWC3'),
        (.85,'SLC22A10','TMEM136','ZNF714')]:
        pos[middle]=(centre,.47)
        pos[left]=(centre-.105,.27)
        pos[right]=(centre+.105,.27)
    for centre, left, right in [(.12,'GSR','AIFM2'),(.37,'ANXA3','MCEMP1'),
                                (.62,'BMPR1A','GREM1'),(.87,'NUP155','THOC7')]:
        pos[left]=(centre-.065,.04)
        pos[right]=(centre+.065,.04)
    assert set(pos) == set(core)
    return core, strip, pos


def node(ax, gene, xy, row, side, strip=False):
    absent = side == 'depmap' and not row.depmap_available
    colour, face = ('#888888', '#D6D6D6') if absent else (STRONG, '#EDE1EF')
    ax.scatter(*xy, s=AREA*abs(row.normZ), facecolor=face, edgecolor=colour, lw=1, zorder=4)
    suffix = ('*' if not row.depmap_available else '') if (strip or side == 'depmap') else ''
    if (strip or side == 'string') and not row.string_available:
        suffix += '†'
    radius = np.sqrt(AREA*abs(row.normZ))/2
    return ax.annotate(gene+suffix, xy, xytext=(radius+2.8, 0) if strip else (0, radius+1.5),
        textcoords='offset points', fontsize=6.9 if strip else 7.2, ha='left' if strip else 'center',
        va='center' if strip else 'bottom', zorder=5,
        bbox=dict(facecolor='white', edgecolor='none', alpha=.35, pad=.1))


def main():
    apply_style(8)
    dz = pd.read_csv(DATA / 'all_gene_scores.csv')
    named = pd.read_csv(DATA / 'main_named_hits.csv')
    guides = pd.read_csv(DATA / 'selected_guides.csv')
    density = pd.read_csv(DATA / 'ntc_density.csv')
    nodes = pd.read_csv(DATA / 'candidates.csv').set_index('GENE')
    pairs = pd.read_csv(DATA / 'pairwise_evidence.csv')
    manifest = json.loads((DATA / 'analysis_manifest.json').read_text())
    fig = plt.figure(figsize=(W, H))

    letter(fig, .14, H-.08, 'a')
    ax = axes(fig, .70, H-1.98, 2.18, 1.61)
    for cl, colour, size in [('FDR ≥ 0.25', NS, 3), ('0.05–0.25', EXPLORATORY, 10),
                             ('FDR < 0.05', STRONG, 14)]:
        sub = dz.loc[dz.display_class.eq(cl)]
        ax.scatter(sub.rank_synth, sub.normZ, s=size, color=colour, lw=0,
                   rasterized=cl == 'FDR ≥ 0.25', zorder=2 if cl == 'FDR ≥ 0.25' else 3)
    ax.axhline(0, color='#888888', lw=.6, zorder=0)
    ax.set(xscale='log', xlim=(.8, 25000), ylim=(-5.9, 4.5), yticks=[-4, -2, 0, 2, 4])
    ax.set_xticks([1, 100, 10000], ['1', '100', '10,000'])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(labelsize=7.5, length=2.5, pad=2)
    ax.set_xlabel('Depletion rank (log scale)', fontsize=8, labelpad=3)
    ax.set_ylabel('drugZ normZ', fontsize=8, labelpad=3)
    ax.text(.04, .95, '19,112 genes', transform=ax.transAxes, fontsize=7.5, va='top')
    indexed = dz.set_index('GENE')
    labels = {'L3MBTL2': (1.2, -3.6), 'UBR5': (1.2, -2.5), 'TM7SF2': (1.2, -1.4),
              'CLYBL': (20, -5.45), 'EIF1AX': (65, -4.25)}
    label_artists = []
    for gene, xy in labels.items():
        p = indexed.loc[gene]
        label_artists.append(ax.annotate(gene, (p.rank_synth, p.normZ), xytext=xy,
            fontsize=6.6, color='#555555', ha='left', va='center',
            arrowprops=dict(arrowstyle='-', color='#A0A0A0', lw=.45, shrinkA=1, shrinkB=2)))
    text(fig, .70, H-2.51, 'FDR <0.25: 73 depleted; 0 enriched', fontsize=7.4)
    letter(fig, .14, H-2.56, 'b')
    marks = []
    bx = axes(fig, 1.10, H-3.67, 1.78, .91)
    guide_panel(bx, ['DDX3X', 'EIF1AX', 'ZFX', 'RPS4X'], guides, density, marks, 'b')
    letter(fig, 3.55, H-.08, 'c')
    for col, start in enumerate([0, 10]):
        cx = axes(fig, 4.37+col*2.91, H-3.67, 1.56, 3.21)
        guide_panel(cx, named.GENE.iloc[start:start+10].tolist(),
                    guides, density, marks, f'c{col+1}', slots=10)
    # Shared reference and statistics key, without repeating long explanations.
    fdr_handles = [Line2D([], [], color=c, marker='o', ms=4, lw=0, label=s) for c, s in [
        (STRONG, 'FDR <0.05'), (EXPLORATORY, '0.05–0.25'), (NS, 'FDR ≥0.25')]]
    fig.legend(handles=fdr_handles, loc='upper left', bbox_to_anchor=(.065, (H-4.04)/H),
        ncol=3, fontsize=7.2, columnspacing=.8, handletextpad=.3, borderaxespad=0)
    fig.legend(handles=[Line2D([], [], color='#333333', marker='|', ms=6, lw=0, label='Guide'),
        Line2D([], [], color='#AAAAAA', lw=2, label='Non-targeting reference (n = 997)')],
        loc='upper left', bbox_to_anchor=(.46, (H-4.04)/H), ncol=2, fontsize=7.2,
        handletextpad=.35, columnspacing=.8, borderaxespad=0)
    text(fig, 9.35, H-4.32, '* One guide omitted per marked gene (no P4 reads)', ha='right', fontsize=6.8, color='#666666')

    core, strip, pos = layout(nodes, pairs)
    letter(fig, .14, H-4.46, 'd')
    n_string = int(pairs.string_edge.sum())
    n_depmap = int(pairs.depmap_edge.sum())
    text(fig, .69, H-4.61, f'STRING: score ≥0.4 · {n_string} associations', fontsize=8)
    text(fig, 5.23, H-4.61, f'DepMap: |r| >0.2 · {n_depmap} connection' + ('s' if n_depmap != 1 else ''), fontsize=8)
    network_labels, net_axes = [], []
    for side, x in [('string', .18), ('depmap', 5.02)]:
        na = axes(fig, x, 1.65, 4.38, 2.00)
        na.set(xlim=(-.08, 1.08), ylim=(-.08, 1.14))
        na.axis('off')
        edges = pairs.loc[pairs.string_edge if side == 'string' else pairs.depmap_edge]
        for e in edges.itertuples():
            xx, yy = zip(pos[e.gene_a], pos[e.gene_b])
            if side == 'string':
                colour, width, ls = '#A3A7AA', .25+2.1*e.string_score, '-'
            else:
                colour = POSITIVE if e.pearson_r > 0 else NEGATIVE
                width = .35+2.3*abs(e.pearson_r)
                ls = '-' if e.pearson_r > 0 else (0, (3, 2))
            na.plot(xx, yy, color=colour, lw=width, ls=ls, alpha=.62, zorder=1)
        network_labels.extend([node(na, g, pos[g], nodes.loc[g], side) for g in core])
        net_axes.append(na)
    text(fig, .36, 1.49, f'Candidates without a displayed edge in either view ({len(strip)})',
         fontsize=7.2, color='#555555')
    strip_ax = axes(fig, .18, .36, 9.22, 1.02)
    strip_ax.set(xlim=(0, 1), ylim=(0, 1))
    strip_ax.axis('off')
    strip_ax.axhline(1, color='#DDDDDD', lw=.65)
    rows = math.ceil(len(strip)/8)
    for i, gene in enumerate(strip):
        xy = (.015+(i%8)/8, 1-(i//8+.55)/rows)
        network_labels.append(node(strip_ax, gene, xy, nodes.loc[gene], 'string', strip=True))
    text(fig, .36, .19, '* DepMap unavailable    † No verified STRING mapping', fontsize=6.7, color='#555555')
    text(fig, 9.35, .19, 'Area: |normZ|    Blue: positive r    Dashed red: negative r',
         ha='right', fontsize=6.7, color='#555555')
    text(fig, .36, .015, 'P4: one CTRL, two LOY. Endpoint contrasts; not adjusted for D0. Original saved analysis.',
         fontsize=7.1, color='#555555', va='bottom')

    table = pd.DataFrame(marks)
    assert len(table) == 90 and table.sgRNA.nunique() == 86
    table.to_csv(DATA / 'displayed_tick_positions.csv', index=False)
    # Data and geometry checks: identical candidates, positions and node areas.
    assert len(net_axes[0].lines) == n_string and len(net_axes[1].lines) == n_depmap
    for i, g in enumerate(core):
        a, b = [na.collections[i] for na in net_axes]
        assert np.array_equal(a.get_offsets(), b.get_offsets())
        assert np.array_equal(a.get_sizes(), b.get_sizes())
        assert np.isclose(a.get_sizes()[0], AREA*abs(nodes.loc[g, 'normZ']))
    assert set(core) | set(strip) == set(nodes.index) and not set(core) & set(strip)
    pd.DataFrame([{'GENE':g, 'region':'network', 'x':pos[g][0], 'y':pos[g][1]} for g in core] +
        [{'GENE':g, 'region':'no_edge_strip', 'x':.015+(i%8)/8, 'y':1-(i//8+.55)/rows}
         for i,g in enumerate(strip)]).to_csv(DATA / 'node_layout.csv', index=False)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    overlaps = [[a.get_text(), b.get_text()] for i, a in enumerate(network_labels)
        for b in network_labels[i+1:]
        if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer))]
    print('Network label overlaps:', overlaps, flush=True)
    assert not overlaps, overlaps
    review = OUT / 'review'
    review.mkdir(exist_ok=True)
    for fmt in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{fmt}', dpi=600 if fmt == 'png' else 300)
    preview = review / f'{STEM}_preview.png'
    fig.savefig(preview, dpi=160)
    digest = hashlib.sha256(preview.read_bytes()).hexdigest()[:12]
    hashed = review / f'{STEM}_{digest}.png'
    hashed.write_bytes(preview.read_bytes())
    record = dict(stem=STEM, preview=str(hashed), inches=[W, H], font='Open Sans',
        panels=['a','b','c','d'], strong_depleted=len(named), network_candidates=len(nodes),
        network_core=len(core), no_edge_strip=len(strip), string_edges=n_string,
        depmap_edges=n_depmap, identical_network_coordinates=True,
        network_label_overlaps=overlaps, guide_marks=len(table), guide_unique_marks=86,
        guide_xlim=manifest['guide_xlim'], node_area_factor=AREA,
        layout='Compact rows of union-graph connected components; distances are illustrative',
        modules='No community shading; insufficient positive DepMap connectivity',
        source_manifest_sha256=hashlib.sha256((DATA/'analysis_manifest.json').read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (DATA/'render_manifest.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    plt.close(fig)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
