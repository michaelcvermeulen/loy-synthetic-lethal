"""Stacked guide-level endpoint contrasts against an empirical NTC background.

This is a descriptive display of frozen MAGeCK values, not a residual model or
a replacement for the primary drugZ analysis. Run from any working directory.
"""
from pathlib import Path
import hashlib
import json
import os
import sys

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_guide_rugs'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

SCREEN = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
INPUTS = {
    'mageck_guides': SCREEN / 'test/POST_LOY_vs_CTRL.sgrna_summary.txt',
    'raw_counts': SCREEN / 'count/ht1080_screen2_R1.count.txt',
    'drugz_genes': SCREEN / 'results/drugZ/LOY_vs_CTRL_unfiltered.txt',
}
STEM = 'ch04_supp_ht1080_guide_consistency_draft'
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2',
           'POST_LOY_REP1', 'POST_LOY_REP2']
PANELS = [
    ['DDX3X', 'EIF1AX', 'RPS4X', 'ZFX', 'SREBF1', 'SREBF2', 'TIMM8A', 'MED12'],
    ['TAF6L', 'TADA2B', 'SUPT20H', 'ATXN7L3', 'TADA1', 'TAF5L', 'ATP6V0E1', 'USP7'],
]
ENRICHED = {'ATP6V0E1', 'USP7'}
XLIM = (-7.2, 5.8)


def fdr_label(value):
    return f'{value:.3g}' if value < .01 else f'{value:.3f}'


def prepare():
    counts = pd.read_csv(INPUTS['raw_counts'], sep='\t')
    guides = pd.read_csv(INPUTS['mageck_guides'], sep='\t')
    scores = pd.read_csv(INPUTS['drugz_genes'], sep='\t').set_index('GENE')
    joined = counts.merge(guides[['sgrna', 'Gene', 'LFC', 'FDR']],
                          left_on=['sgRNA', 'Gene'], right_on=['sgrna', 'Gene'],
                          how='left', validate='one_to_one')
    assert len(joined) == 77441 and joined.sgRNA.is_unique
    assert joined.LFC.notna().all()
    joined['endpoint_read_total'] = joined[SAMPLES[2:]].sum(axis=1)
    joined['endpoint_detected_samples'] = joined[SAMPLES[2:]].gt(0).sum(axis=1)
    joined['all_samples_zero'] = joined[SAMPLES].sum(axis=1).eq(0)
    joined['lfc_displayed'] = joined.endpoint_read_total.gt(0) & np.isfinite(joined.LFC)

    selected = joined.loc[joined.Gene.isin(sum(PANELS, []))].copy()
    selected = selected.sort_values(['Gene', 'sgRNA'])
    selected['guide_number'] = selected.groupby('Gene').cumcount() + 1
    selected['sequence_5prime_to_3prime'] = selected.sgRNA.str.rsplit('_', n=1).str[-1]
    selected['display_panel'] = selected.Gene.map({g: chr(97+i) for i,p in enumerate(PANELS) for g in p})
    selected['display_row'] = selected.Gene.map({g: j+1 for p in PANELS for j,g in enumerate(p)})
    selected['drugz_normZ'] = selected.Gene.map(scores.normZ)
    selected['drugz_fdr_direction'] = np.where(selected.Gene.isin(ENRICHED), 'enrichment', 'depletion')
    selected['drugz_directional_fdr'] = [scores.loc[r.Gene, 'fdr_supp' if r.Gene in ENRICHED else 'fdr_synth']
                                        for r in selected.itertuples()]
    selected = selected.rename(columns={'LFC': 'mageck_reported_log2FC', 'FDR': 'mageck_guide_FDR'})
    selected = selected.drop(columns='sgrna').sort_values(['display_panel', 'display_row', 'guide_number'])
    assert len(selected) == 64 and selected.groupby('Gene').size().eq(4).all()
    assert selected.lfc_displayed.sum() == 62
    selected.to_csv(DATA / 'selected_guides.csv', index=False)

    ntc = joined.loc[joined.Gene.eq('Non-Targeting-Control')].copy()
    assert len(ntc) == 1000 and ntc.lfc_displayed.sum() == 999
    ntc.to_csv(DATA / 'non_targeting_guides.csv', index=False)
    values = ntc.loc[ntc.lfc_displayed, 'LFC'].to_numpy()
    # Explicit normal-reference bandwidth; no fitted null or normality assumption.
    scale = min(values.std(ddof=1), np.subtract(*np.percentile(values, [75, 25])) / 1.34)
    bandwidth = .9 * scale * len(values) ** (-.2)
    grid = np.linspace(*XLIM, 2601)
    density = np.exp(-.5 * ((grid[:, None] - values[None, :]) / bandwidth) ** 2).mean(axis=1)
    density /= bandwidth * np.sqrt(2 * np.pi)
    pd.DataFrame({'guide_log2FC': grid, 'kde_density': density,
                  'display_density_peak_scaled': density / density.max()}).to_csv(DATA / 'ntc_density.csv', index=False)
    row_summary = selected.groupby(['display_panel', 'display_row', 'Gene'], sort=False).apply(
        lambda d: pd.Series({
            'guides_detected': int(d.lfc_displayed.sum()),
            'negative_LFC_guides': int((d.lfc_displayed & d.mageck_reported_log2FC.lt(0)).sum()),
            'positive_LFC_guides': int((d.lfc_displayed & d.mageck_reported_log2FC.gt(0)).sum()),
            'drugz_normZ': d.drugz_normZ.iloc[0],
            'drugz_fdr_direction': d.drugz_fdr_direction.iloc[0],
            'drugz_directional_fdr': d.drugz_directional_fdr.iloc[0],
        }), include_groups=False).reset_index()
    row_summary.to_csv(DATA / 'gene_summary.csv', index=False)
    manifest = {
        'inputs': {k: {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                   for k, p in INPUTS.items()},
        'description': 'Saved MAGeCK endpoint LOY-versus-CTRL guide log2 fold changes; descriptive display only.',
        'contrast': 'Two endpoint outgrowth cultures per genotype; D0 excluded from both endpoint contrasts.',
        'ntc_total': 1000, 'ntc_with_endpoint_reads': 999, 'ntc_excluded_zero_endpoint_reads': 1,
        'density': {'method': 'Gaussian KDE', 'bandwidth': float(bandwidth),
                    'bandwidth_rule': '0.9 * min(sample SD, IQR/1.34) * n**(-1/5)',
                    'same_curve_each_row': True, 'height': 'scaled to common peak; no density y axis',
                    'ntc_observed_range': [float(values.min()), float(values.max())]},
        'selected_genes': PANELS, 'selected_guides': 64, 'displayed_guides': 62,
        'excluded_selected': selected.loc[~selected.lfc_displayed, ['sgRNA', 'mageck_reported_log2FC']].to_dict('records'),
        'x_limits': XLIM, 'primary_gene_analysis': 'Original drugZ directional FDRs; no recalculation',
        'selection': 'Four established gametologs; SREBF1/2; TIMM8A/MED12; six SAGA candidates; two enriched hits at FDR<0.05.',
        'normalisation_caveat': 'MAGeCK used non-targeting controls; their centring is not independent QC.',
    }
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    return selected, grid, density / density.max(), manifest


def render(selected, grid, curve, manifest):
    apply_style(9)
    fig = plt.figure(figsize=(9.6, 6.15))
    marks = []
    for panel, genes in enumerate(PANELS):
        x = .113 + panel * .5
        ax = fig.add_axes([x, .195, .306, .715])
        fig.text(.014 + panel*.5, .973, chr(97+panel), fontsize=13, fontweight='semibold', va='top')
        fig.text(x+.317, .929, 'drugZ FDR', fontsize=7.5, ha='left', va='bottom')
        for j, gene in enumerate(genes):
            base = 7-j
            data = selected.loc[selected.Gene.eq(gene)]
            ax.fill_between(grid, base+.045, base+.045+curve*.79, color='#EBEBEB', linewidth=0, zorder=0)
            ax.plot(grid, base+.045+curve*.79, color='#AAAAAA', lw=.8, zorder=1)
            ax.hlines(base+.045, *XLIM, color='#D8D8D8', lw=.5, zorder=1)
            ax.vlines(0, base+.045, base+.92, color='#C4C4C4', lw=.55, ls=(0,(2,2)), zorder=1)
            for r in data.itertuples():
                if not r.lfc_displayed:
                    continue
                # Height offsets expose near-overlapping LFCs without horizontal jitter.
                top = base + .38 + .17*(r.guide_number-1)
                ax.plot([r.mageck_reported_log2FC]*2, [base+.045,top], color='#303030', lw=1.05, zorder=3)
                ax.plot(r.mageck_reported_log2FC, top, 'o', color='#303030', ms=2.4, mew=0, zorder=4)
                marks.append({'sgRNA':r.sgRNA, 'x':r.mageck_reported_log2FC, 'top':top, 'panel':chr(97+panel)})
            missing = not data.lfc_displayed.all()
            ax.text(-.035, base+.48, gene+('†' if missing else ''), transform=ax.get_yaxis_transform(),
                    ha='right', va='center', fontsize=9,
                    fontweight='semibold' if panel == 1 and j < 6 else 'normal')
            fdr = data.drugz_directional_fdr.iloc[0]
            ax.text(1.035, base+.48, fdr_label(fdr), transform=ax.get_yaxis_transform(),
                    ha='left', va='center', fontsize=8)
        # Subtle spaces in reading groups without changing row heights.
        for y in ([3.99, 1.99] if panel == 0 else [1.99]):
            ax.axhline(y, color='white', lw=2, zorder=2)
        ax.set_xlim(*XLIM)
        ax.set_ylim(-.12, 8.02)
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.set_xticks([-6,-4,-2,0,2,4])
        ax.tick_params(axis='x',labelsize=8, length=3, pad=3)
        ax.set_xlabel('Guide log₂ fold change (LOY / CTRL)', fontsize=8.5, labelpad=4)
        ax.text(0,-.135,'Depleted in LOY', transform=ax.transAxes, fontsize=7.5, va='top')
        ax.text(1,-.135,'Enriched in LOY', transform=ax.transAxes, fontsize=7.5, va='top',ha='right')
    handles = [Line2D([],[], color='#303030', marker='o', lw=1, ms=2.8,
                      label='Detected guide (saved MAGeCK estimate)'),
               Line2D([],[], color='#AAAAAA', lw=2,
                      label='Non-targeting guide distribution (n = 999)')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(.018,.035),
               ncol=2, fontsize=8, handlelength=1.7, columnspacing=2.3, borderaxespad=0)
    fig.text(.025,.010,'† One guide has no reads; its fold change is omitted. Tick heights are staggered to reveal overlapping guides.', fontsize=7.5)
    fig.canvas.draw()
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            bb=t.get_window_extent(fig.canvas.get_renderer())
            assert bb.x0 >= 0 and bb.y0 >= 0 and bb.x1 <= fig.bbox.width+.1 and bb.y1 <= fig.bbox.height+.1, t.get_text()
    for ext, dpi in [('pdf',300), ('svg',300), ('png',600)]:
        fig.savefig(OUT/f'{STEM}.{ext}', dpi=dpi)
    review=OUT/'review';review.mkdir(exist_ok=True)
    fig.savefig(review/f'{STEM}_preview.png', dpi=160)
    plt.close(fig)
    pd.DataFrame(marks).to_csv(DATA/'displayed_tick_positions.csv',index=False)
    manifest.update({'figure_inches':[9.6,6.15],'font':'Open Sans','png_dpi':600,
                     'render_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    (DATA/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({'stem':STEM,'genes':16,'guides_plotted':len(marks),'background_guides':999,
                      'KDE_bandwidth':manifest['density']['bandwidth']}))


if __name__ == '__main__':
    render(*prepare())
