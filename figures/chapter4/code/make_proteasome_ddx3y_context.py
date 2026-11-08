"""Full annotated proteasome set and a focused guide-representation supplement."""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/proteasome_ddx3y_context'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

GROUPS = ['XY_low_DDX3Y', 'XY_high_DDX3Y', 'LOY', 'Female']
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2',
           'POST_LOY_REP1', 'POST_LOY_REP2']
FOLLOWUP = ['PSMD8', 'PSMD11', 'PSMD2', 'PSMC1', 'PSMC2', 'PSMC4']


def audit_export(fig, stem, dimensions, info):
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{stem}.{ext}', dpi=600 if ext == 'png' else 300)
    fig.savefig(OUT / 'review' / f'{stem}_preview.png', dpi=145)
    plt.close(fig)
    with fitz.open(OUT / f'{stem}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside = [s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps = []
        for i, a in enumerate(spans):
            for b in spans[i + 1:]:
                ra, rb = fitz.Rect(a['bbox']), fitz.Rect(b['bbox'])
                area = (ra & rb).get_area()
                if area > 1 and area / min(ra.get_area(), rb.get_area()) > .12:
                    overlaps.append([a['text'], b['text']])
        assert not outside, outside
        assert not overlaps, overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.1, 1.1), colorspace=fitz.csGRAY).save(
            OUT / 'review' / f'{stem}_greyscale.png')
    result = dict(stem=stem, dimensions_inches=dimensions, font='Open Sans',
                  text_outside_page=outside, text_overlaps=overlaps, **info)
    (DATA / f'{stem}_render_manifest.json').write_text(
        json.dumps(result, indent=2), encoding='utf-8')
    print(stem, 'passed text and font checks', flush=True)


def main_figure():
    apply_style(8)
    data = pd.read_csv(DATA / 'combined_evidence.csv').sort_values('display_order')
    n = len(data)
    assert n == 66
    fig = plt.figure(figsize=(10.6, 12.4))
    bottom, height = .082, .817
    heat = fig.add_axes([.123, bottom, .257, height])
    forest = fig.add_axes([.465, bottom, .20, height], sharey=heat)
    drug = fig.add_axes([.792, bottom, .145, height], sharey=heat)
    corr = data[GROUPS].to_numpy(float)
    cmap = LinearSegmentedColormap.from_list('correlation', ['#3C5488', '#F7F7F7', '#C96C2D'])
    cmap.set_bad('#D5D5D5')
    im = heat.imshow(np.ma.masked_invalid(corr), cmap=cmap, vmin=-.5, vmax=.5, aspect='auto')
    heat.set(yticks=np.arange(n), yticklabels=data.gene.tolist(), ylim=(n-.5, -.5))
    heat.set_xticks(range(4), ['XY lower\nDDX3Y\nn = 153', 'XY higher\nDDX3Y\nn = 152',
                              'LOY\n\nn = 216', 'Female\n\nn = 448'])
    heat.xaxis.tick_top()
    heat.tick_params(axis='both', length=0, labelsize=7.7)
    heat.tick_params(axis='x', pad=6)
    heat.tick_params(axis='y', pad=6)
    heat.set_xticks(np.arange(-.5, 4, 1), minor=True)
    heat.set_yticks(np.arange(-.5, n, 1), minor=True)
    heat.grid(which='minor', color='white', lw=.65)
    heat.tick_params(which='minor', bottom=False, top=False, left=False)
    for spine in heat.spines.values():
        spine.set_visible(False)
    forest.axvline(0, color='#888888', ls='--', lw=.65)
    drug.axvline(0, color='#888888', ls='--', lw=.65)
    for ax in [forest, drug]:
        ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        ax.spines['left'].set_visible(False)
    assert data.ci_low.min() > -.15 and data.ci_high.max() < .15
    forest.set(xlim=(-.15, .15), xticks=[-.1, 0, .1])
    drug.set(xlim=(-3, 3), xticks=[-2, 0, 2])
    forest.set_xlabel('Interaction coefficient', fontsize=8.4)
    drug.set_xlabel('HT1080 drugZ normZ', fontsize=8.4)
    for i, row in enumerate(data.itertuples()):
        focus = row.gene == 'PSMD8'
        if focus:
            heat.get_yticklabels()[i].set_fontweight('semibold')
            for ax in [forest, drug]:
                ax.axhspan(i-.5, i+.5, color='#F0E8F5', zorder=0)
        for j in range(4):
            r = corr[i, j]
            label = 'NA' if not np.isfinite(r) else f'{r:.2f}'.replace('-0.00', '0.00')
            heat.text(j, i, label, ha='center', va='center', fontsize=7.1,
                      color='white' if np.isfinite(r) and abs(r) > .38 else '#222222')
        if np.isfinite(row.effect):
            color = '#7A3E85' if focus else '#666666'
            forest.plot([row.ci_low, row.ci_high], [i, i], color=color, lw=.8)
            forest.plot(row.effect, i, 'o', color=color, ms=3.2)
            forest.text(1.075, i, f'{row.fdr_bh_set:.3f}',
                        transform=forest.get_yaxis_transform(), ha='left', va='center',
                        fontsize=7.2, color='#333333')
        else:
            forest.text(0, i, 'NA', ha='center', va='center', fontsize=7.2, color='#888888')
            forest.text(1.075, i, 'NA', transform=forest.get_yaxis_transform(),
                        ha='left', va='center', fontsize=7.2, color='#888888')
        if np.isfinite(row.normZ):
            drug.plot(row.normZ, i, marker='D' if focus else 'o', ms=3.5 if focus else 3.0,
                      color='#7A3E85' if focus else '#666666', linestyle='none')
        else:
            drug.text(0, i, 'NA', ha='center', va='center', fontsize=7.2, color='#888888')
    for annotation, rows in data.groupby('annotation_group', sort=False):
        lo, hi = rows.display_order.min(), rows.display_order.max()
        short = annotation.replace(' annotation', '').replace('Other associated', 'Other associated proteins')
        heat.text(-.35, (lo+hi)/2, short, transform=heat.get_yaxis_transform(),
                  rotation=90, ha='center', va='center', fontsize=8.0, color='#666666')
        if lo > 0:
            for ax in [heat, forest, drug]:
                ax.axhline(lo-.5, color='#888888', lw=.65)
    # Short headings separate meaning without placing a manuscript caption on the canvas.
    for x, letter, label in [(.065, 'a', 'DDX3X co-dependency'),
                              (.425, 'b', 'DDX3Y interaction'),
                              (.75, 'c', 'HT1080 screen')]:
        fig.text(x, .965, letter, fontsize=12, fontweight='semibold')
        fig.text(x+.04, .965, label, fontsize=9.2, fontweight='semibold')
    forest.text(1.075, 1.012, 'FDR', transform=forest.transAxes, ha='left',
                va='bottom', fontsize=8)
    drug.text(.5, 1.012, 'All FDR ≥ 0.25', transform=drug.transAxes,
              ha='center', va='bottom', fontsize=7.6, color='#555555')
    cbax = fig.add_axes([.17, .043, .16, .009])
    cb = fig.colorbar(im, cax=cbax, orientation='horizontal', ticks=[-.5, 0, .5])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7.4, length=2)
    cb.set_label('Adjusted Pearson r', fontsize=8, labelpad=2)
    fig.text(.566, .018, 'Negative: slope decreases\nwith higher DDX3Y RNA',
             ha='center', fontsize=7.7, color='#555555')
    fig.text(.865, .018, 'Negative: LOY-depleted\nPositive: LOY-enriched',
             ha='center', fontsize=7.7, color='#555555')
    audit_export(fig, 'ch04_supp_proteasome_ddx3y_context_draft', [10.6, 12.4],
                 dict(rows=n, correlation_cells=int(corr.size),
                      primary_interaction_tests=int(data.effect.notna().sum()),
                      correlation_limits=[-.5, .5], forest_limits=[-.15, .15],
                      drugz_limits=[-3, 3], no_clipping=True,
                      effect_units='DDX3X effect / partner XY SD / DDX3Y log2(TPM+1) unit.',
                      highlighted_gene='PSMD8; emphasis denotes follow-up focus, not significance.',
                      group_labels='GO-derived annotation groups; not a new curated structural assignment.'))


def guide_figure():
    apply_style(8.5)
    guides = pd.read_csv(DATA / 'ht1080_guide_evidence.csv')
    guides = guides.loc[guides.gene.isin(FOLLOWUP)].copy()
    guides['gene_order'] = guides.gene.map({g: i for i, g in enumerate(FOLLOWUP)})
    guides = guides.sort_values(['gene_order', 'guide_number'])
    assert len(guides) == 24 and not guides.all_samples_zero.any()
    cpm = guides[[s+'_cpm' for s in SAMPLES]].to_numpy(float)
    values = np.log2(cpm + 1)
    vmax = float(np.ceil(values.max()))
    fig = plt.figure(figsize=(9.3, 7.0))
    heat = fig.add_axes([.15, .15, .52, .70])
    effect = fig.add_axes([.755, .15, .19, .70], sharey=heat)
    im = heat.imshow(values, cmap='magma', vmin=0, vmax=vmax, aspect='auto')
    heat.set_yticks(range(24), [f'{r.gene} g{r.guide_number}' for r in guides.itertuples()])
    heat.set_xticks(range(6), ['D0\nCTRL', 'D0\nLOY', 'Endpoint\nCTRL 1',
                              'Endpoint\nCTRL 2', 'Endpoint\nLOY 1', 'Endpoint\nLOY 2'])
    heat.xaxis.tick_top()
    heat.tick_params(axis='both', length=0, pad=6, labelsize=8)
    heat.set_xticks(np.arange(-.5, 6, 1), minor=True)
    heat.set_yticks(np.arange(-.5, 24, 1), minor=True)
    heat.grid(which='minor', color='white', lw=.6)
    heat.tick_params(which='minor', bottom=False, top=False, left=False)
    for spine in heat.spines.values():
        spine.set_visible(False)
    effect.spines['left'].set_visible(False)
    effect.tick_params(axis='y', which='both', left=False, labelleft=False)
    effect.axvline(0, color='#777777', ls='--', lw=.7)
    effect.set_xlabel('Endpoint contrast\nLOY − CTRL', fontsize=8.5)
    contrast_limit = float(max(2, np.ceil(guides.endpoint_loy_minus_ctrl.abs().max())))
    effect.set_xlim(-contrast_limit, contrast_limit)
    effect.set_xticks([-contrast_limit, 0, contrast_limit])
    for i, row in enumerate(guides.itertuples()):
        for j in range(6):
            heat.text(j, i, f'{cpm[i,j]:.2g}', ha='center', va='center', fontsize=7.8,
                      color='white' if values[i,j] < vmax*.62 else '#222222')
        effect.plot(row.endpoint_loy_minus_ctrl, i, 'o', ms=4,
                    color='#7A3E85' if row.gene == 'PSMD8' else '#666666')
        if i % 4 == 0 and i > 0:
            heat.axhline(i-.5, color='#AAAAAA', lw=.9)
            effect.axhline(i-.5, color='#AAAAAA', lw=.9)
    fig.text(.04, .94, 'a', fontsize=12, fontweight='semibold')
    fig.text(.15, .94, 'Guide abundance (cell labels: CPM)', fontsize=9.5, fontweight='semibold')
    fig.text(.72, .94, 'b', fontsize=12, fontweight='semibold')
    fig.text(.78, .94, 'Guide contrast', fontsize=9.5, fontweight='semibold')
    cbax = fig.add_axes([.30, .065, .24, .016])
    cb = fig.colorbar(im, cax=cbax, orientation='horizontal')
    cb.outline.set_visible(False)
    cb.set_label('Colour: log₂(CPM + 1)', fontsize=8.5, labelpad=3)
    cb.ax.tick_params(labelsize=8, length=2)
    audit_export(fig, 'ch04_supp_proteasome_ht1080_guides_draft', [9.3, 7.0],
                 dict(genes=FOLLOWUP, guides=24, sample_measurements=144,
                      heatmap_labels='CPM; colour log2(CPM+1).',
                      colour_limits=[0, vmax], guide_contrast_limits=[-contrast_limit, contrast_limit],
                      contrast='Mean endpoint log2(CPM+1), LOY minus CTRL.',
                      no_guides_excluded=True, no_clipping=True))


if __name__ == '__main__':
    main_figure()
    guide_figure()
