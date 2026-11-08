"""Panels h/i: unaltered gel pixels and frozen, verified clone ddPCR values."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from PIL import Image
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parents[1]
SOURCE = OUT / 'source_data/ht1080_ko_validation_20261022'
DATA = OUT / 'source_data/competition_workflow_qc'
GEL_SOURCE = OUT / 'source_data/ht1080_author_ko_gel_20261022/author_supplied_ko_gel.png'
GEL_BOX = (0, 420, 709, 575)
CLONES = ['LOY1', 'LOY3', 'LOY4', 'LOY5', 'LOY7', 'LOY8',
          'CTRL1', 'CTRL3', 'CTRL4', 'CTRL6', 'CTRL7', 'CTRL8']
ASSAYS = ['RPS4Y1 (Yp)', 'AMELY/AMELX (Yp)', 'KDM5D (Yq)']
COLOURS = ['#2a78d6', '#eb6834', '#1baf7a']
MARKERS = ['o', 's', '^']


def draw_validation(fig, width, height):
    def label(x, y, value, **kw):
        return fig.text(x/width, y/height, value, **kw)

    full = Image.open(GEL_SOURCE)
    gel = full.crop(GEL_BOX)
    gel.save(DATA / 'panel_h_gel_native.png')
    x, w = 63, 58
    h = w * gel.height / gel.width
    # Keep the cropped band region directly beneath its existing lane labels.
    y = 48.5 - h
    ax = fig.add_axes([x/width, y/height, w/width, h/height])
    ax.imshow(gel, interpolation='none')
    ax.axis('off')
    label(60, 64, 'h', fontsize=9, weight='semibold', va='center')
    label(91, 64, 'Genomic PCR', fontsize=7.7, weight='semibold', ha='center', va='center')
    # Same labelled lanes as the grant panel, in the author-supplied image.
    lane_x = [[110, 159, 209, 257], [306, 356, 404, 452], [500, 549, 596, 646]]
    def px(v):
        return x + (v-GEL_BOX[0])/gel.width*w
    for gene, lanes in zip(['DDX3X', 'TADA2B', 'TAF6L'], lane_x):
        label(px(np.mean(lanes)), 59.7, gene, fontsize=6.2, ha='center', va='center')
        fig.add_artist(Line2D([px(lanes[0]-18)/width, px(lanes[-1]+18)/width],
                            [58.1/height]*2, transform=fig.transFigure, color='#444444', lw=.5))
        label(px(np.mean(lanes[2:])), 55.9, 'KO', fontsize=5.8, ha='center', va='center')
        fig.add_artist(Line2D([px(lanes[2]-18)/width, px(lanes[-1]+18)/width],
                            [54.4/height]*2, transform=fig.transFigure, color='#444444', lw=.5))
        for lane, name in zip(lanes, ['CTL', 'LOY', 'CTL', 'LOY']):
            label(px(lane), 51.7, name, fontsize=6, rotation=90, ha='center', va='center')

    d = pd.read_csv(SOURCE / 'display_values.csv')
    reference = pd.read_csv(SOURCE / 'DDX3X_clones_grant_source.csv')
    check = d.merge(reference, on=['clone', 'assay'], validate='one_to_one', suffixes=('', '_source'))
    assert len(d) == len(check) == 36 and set(d.label) == set(CLONES)
    assert np.array_equal(check.dosage, check.dosage_source)
    assert np.allclose(d.dosage, d.ratio/d.exp_ratio, atol=1e-12, rtol=0)
    ax = fig.add_axes([136/width, 13/height, 44/width, 47/height])
    positions = list(range(13, 7, -1)) + list(range(6, 0, -1))
    yy = dict(zip(CLONES, positions))
    ax.axvline(1, color='#777777', lw=.6, zorder=0)
    ax.axvline(0, color='#BBBBBB', lw=.5, zorder=0)
    for value in positions:
        ax.axhline(value, color='#EEEEEE', lw=.4, zorder=0)
    for j, assay in enumerate(ASSAYS):
        s = d.loc[d.assay.eq(assay)]
        ax.scatter(s.dosage, s.label.map(yy)+[.22, 0, -.22][j], s=8,
                   marker=MARKERS[j], facecolors='white', edgecolors=COLOURS[j],
                   linewidths=.6, zorder=3)
    ax.set_yticks(positions, CLONES)
    ax.set(xlim=(-.045, 1.56), ylim=(.35, 14.4), xticks=[0, .5, 1, 1.5], xlabel='Relative Y dosage')
    ax.tick_params(axis='y', labelsize=5.9, length=0, pad=2)
    ax.tick_params(axis='x', labelsize=6.2, length=1.6, pad=1.2)
    ax.xaxis.label.set_size(6.5)
    ax.xaxis.labelpad=2
    ax.spines['left'].set_visible(False)
    ax.text(0, 13.95, 'LOY (n = 6)', fontsize=6.1, weight='semibold')
    ax.text(0, 7.1, 'Y-retaining (n = 6)', fontsize=6.1, weight='semibold')
    label(125, 64, 'i', fontsize=9, weight='semibold', va='center')
    label(158, 64, 'Clone Y status', fontsize=7.7, weight='semibold', ha='center', va='center')
    handles=[Line2D([], [], marker=m, color=c, markerfacecolor='white', markersize=3,
                    linestyle='none', label=assay) for m, c, assay in zip(MARKERS, COLOURS, ASSAYS)]
    fig.legend(handles=handles, loc='lower right', bbox_to_anchor=(181/width, .1/height),
               ncol=3, fontsize=6.1, columnspacing=.7, handletextpad=.3, borderaxespad=0)
    record={'gel_source_crop_ltrb': list(GEL_BOX), 'gel_size_px': list(gel.size),
            'gel_pixels_unchanged': True, 'gel_source': str(GEL_SOURCE),
            'gel_labels': 'Author-supplied gene/KO/lane labels retypeset in Open Sans; lane order unchanged.',
            'clone_count': 12, 'measurement_count': 36, 'values_match_frozen_grant_csv': True,
            'sources_sha256': {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [GEL_SOURCE, SOURCE/'display_values.csv', SOURCE/'DDX3X_clones_grant_source.csv']}}
    (DATA/'validation_panels_manifest.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
