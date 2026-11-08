"""Display the three saved DU145 MCODE modules; no reanalysis.

The former selected-pair panels d-f were retired on 9 October 2026.
Their exports, display tables and renderer are preserved in the thesis archive.
"""
from pathlib import Path
import hashlib
import json
import sys
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
INPUT = FIG / 'source_data/du145_metascape_current'
SOURCE = FIG / 'source_data/du145_metascape_ppi_modules'
STEM = 'ch03_supp_du145_metascape_ppi_modules'
NS = {'x': 'http://www.cs.rpi.edu/XGMML'}
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8.5)

PANELS = [
    dict(direction='up', cluster=1, letter='a', colour='#00A087',
         direction_label='Higher in Y-depleted', theme='Collagen-associated',
         positions={'P3H2': (0, .78), 'COL4A2': (.78, .22),
                    'COL12A1': (.52, -.70), 'TAF9': (-.52, -.70),
                    'COL4A1': (-.78, .22)}),
    dict(direction='down', cluster=1, letter='b', colour='#8491B4',
         direction_label='Lower in Y-depleted', theme='Immune-associated',
         positions={'CSF1': (.62, .62), 'SYK': (.62, -.62),
                    'B2M': (-.62, -.62), 'SIRPA': (-.62, .62)}),
    dict(direction='down', cluster=2, letter='c', colour='#B89139',
         direction_label='Lower in Y-depleted', theme='Chemokine-associated',
         positions={'CCR9': (-.62, .62), 'CCL22': (.62, .62),
                    'CXCL1': (0, -.70)}),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_modules():
    node_rows, edge_rows, sources = [], [], []
    for direction in ['up', 'down']:
        folder = INPUT / direction / 'Enrichment_PPI'
        path = folder / 'xgmml/MyList_MCODE_ALL.xgmml'
        sources.extend([path, folder / 'MCODE.csv', folder / 'GO_MCODE_Top3.csv'])
        root = ET.parse(path).getroot()
        symbols, clusters = {}, {}
        for n in root.findall('x:node', NS):
            attrs = {a.get('name'): a.get('value') for a in n.findall('x:att', NS)}
            symbol, cluster = attrs['Symbol'], int(attrs['MCODE_CLUSTER_ID'])
            symbols[n.get('id')] = symbol
            clusters[symbol] = cluster
            node_rows.append(dict(direction=direction, cluster=cluster,
                                  symbol=symbol, entrez_id=attrs['Gene'],
                                  mcode_score=float(attrs['MCODE_SCORE']),
                                  mcode_type=attrs['MCODE_TYPE']))
        seen = set()
        for e in root.findall('x:edge', NS):
            a, b = sorted([symbols[e.get('source')], symbols[e.get('target')]])
            assert clusters[a] == clusters[b], 'Unexpected cross-module edge'
            assert (a, b) not in seen, 'Unexpected duplicate exported edge'
            seen.add((a, b))
            attrs = {x.get('name'): x.get('value') for x in e.findall('x:att', NS)}
            edge_rows.append(dict(direction=direction, cluster=clusters[a],
                                  source=a, target=b, exported_attributes=json.dumps(attrs)))
        table = pd.read_csv(folder / 'MCODE.csv')
        exported = {(r['symbol'], r['cluster']) for r in node_rows if r['direction'] == direction}
        assert exported == set(zip(table.Symbol, table.Cluster))
    return pd.DataFrame(node_rows), pd.DataFrame(edge_rows), sources


def main():
    SOURCE.mkdir(exist_ok=True)
    nodes, edges, inputs = load_modules()
    assert len(nodes) == 12 and len(edges) == 17
    nodes['display_type'] = 'MCODE_module'
    edges['display_type'] = 'MCODE_module'
    for spec in PANELS:
        for table in [nodes, edges]:
            mask = (table.direction == spec['direction']) & (table.cluster == spec['cluster'])
            table.loc[mask, 'panel'] = spec['letter']
    height = 72
    fig = plt.figure(figsize=(183 / 25.4, height / 25.4), facecolor='white')
    summaries = []
    for i, spec in enumerate(PANELS):
        direction, cluster = spec['direction'], spec['cluster']
        n = nodes[nodes.panel == spec['letter']]
        e = edges[edges.panel == spec['letter']]
        pos = spec['positions']
        assert set(n.symbol) == set(pos)
        left = .014 + (i % 3) * .328
        top, bottom = 10, 59
        header_y = 4
        ax = fig.add_axes([left, 1 - bottom / height, .318, (bottom - top) / height])
        ax.set(xlim=(-1.39, 1.39), ylim=(-1.12, 1.12), aspect='equal')
        ax.axis('off')
        for row in e.itertuples():
            p, q = pos[row.source], pos[row.target]
            ax.plot([p[0], q[0]], [p[1], q[1]], color='#ADB2B9', lw=1.25, zorder=1)
        for symbol, (x, y) in pos.items():
            ax.scatter([x], [y], s=145, c=spec['colour'], edgecolors='white', linewidths=.8, zorder=3)
            dx = 9 if x > .1 else (-9 if x < -.1 else 0)
            dy = 0 if abs(x) > .1 else (10 if y > 0 else -10)
            ax.annotate(symbol, (x, y), xytext=(dx, dy), textcoords='offset points',
                        ha='left' if dx > 0 else ('right' if dx < 0 else 'center'),
                        va='bottom' if dy > 0 else ('top' if dy < 0 else 'center'),
                        fontsize=8.5, color='#242424', zorder=4)
        fig.text(left, 1 - header_y / height, spec['letter'], ha='left', va='top', fontsize=11, weight='semibold')
        fig.text(left + .031, 1 - (header_y + .3) / height, spec['direction_label'], ha='left', va='top', fontsize=8.4)
        fig.text(left + .159, 1 - 60 / height, spec['theme'], ha='center', va='top', fontsize=8.4, weight='semibold')
        score = float(n.mcode_score.iloc[0])
        detail = f'{len(n)} proteins; {len(e)} links  |  MCODE {cluster}'
        fig.text(left + .159, 1 - 65 / height, detail,
                 ha='center', va='top', fontsize=7.5, color='#555555')
        summaries.append(dict(direction=direction, cluster=cluster, panel=spec['letter'],
                              display_type='MCODE_module', theme=spec['theme'],
                              proteins=len(n), links=len(e), mcode_score=score,
                              colour=spec['colour']))
        for symbol, (x, y) in pos.items():
            mask = (nodes.direction == direction) & (nodes.symbol == symbol)
            nodes.loc[mask, ['x', 'y']] = [x, y]
    outputs = []
    for extension in ['pdf', 'svg', 'png']:
        path = FIG / f'{STEM}.{extension}'
        fig.savefig(path, dpi=350, facecolor='white')
        outputs.append(path)
    preview = FIG / 'previews' / f'{STEM}_preview.png'
    fig.savefig(preview, dpi=180, facecolor='white')
    plt.close(fig)
    nodes.to_csv(SOURCE / 'nodes.csv', index=False)
    edges.to_csv(SOURCE / 'edges.csv', index=False)
    pd.DataFrame(summaries).to_csv(SOURCE / 'modules.csv', index=False)
    provenance = dict(
        figure=STEM, date='2026-10-09', status='Selected as CH3-S3.4.6.3-S09; three-panel layout adopted on 2026-10-09',
        method='Display all 12 nodes and 17 edges in the saved UP and DOWN MCODE_ALL networks; no PPI or enrichment analysis rerun.',
        inputs={str(p.relative_to(THESIS)): sha(p) for p in inputs},
        outputs={str(p.relative_to(THESIS)): sha(p) for p in outputs},
        layout='Fixed compact module positions; positions carry no biological distance meaning.',
        node_area='Uniform; not expression magnitude or connectivity.',
        edges='Uniform undirected grey lines; database associations of mixed evidence; not measured in DU145.',
        colour='Categorical MCODE module identity; expression direction is labelled explicitly.',
        retired_panels='Selected pairs formerly shown as d-f; archived in archive/DU145_S09_remove_pairs_20261009.',
        modules=summaries, font='Bundled Open Sans', size_mm=[183, height])
    (SOURCE / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'figure': STEM, 'modules': summaries, 'outputs': [str(p) for p in outputs]}, indent=2))


if __name__ == '__main__':
    main()
