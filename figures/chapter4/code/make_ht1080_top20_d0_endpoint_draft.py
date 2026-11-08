"""Top-20 D0/endpoint draft: omit DDX3X guide 2 from display, retain analyses."""
from pathlib import Path
import hashlib
import json
import os
import sys

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_top20_d0_endpoint_draft'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pymupdf as fitz

sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(9)

STEM = 'ch04_supp_ht1080_top20_d0_endpoint_draft'
SCREEN = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
SOURCES = {
    'counts': SCREEN / 'count/ht1080_screen2_R1.count.txt',
    'library': SCREEN / 'brunello/broadgpp-brunello-library-corrected.txt',
    'drugz_scores': SCREEN / 'results/drugZ/LOY_vs_CTRL_unfiltered.txt',
    'existing_profiles': OUT / 'source_data/ht1080_candidate_discovery/guide_profiles.csv',
}
GENES = []  # Filled from the saved directional drugZ depletion ranking.
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2',
           'POST_LOY_REP1', 'POST_LOY_REP2']
COLORS = {'CTRL': '#3C5488', 'LOY': '#E64B35'}
EXCLUDED_GUIDE = 'DDX3X_4638_GCACCACCATAAACCACGCA'


def main():
    counts = pd.read_csv(SOURCES['counts'], sep='\t')
    assert len(counts) == 77441 and counts.sgRNA.is_unique
    cpm = counts[SAMPLES].div(counts[SAMPLES].sum(), axis=1) * 1e6
    scores = pd.read_csv(SOURCES['drugz_scores'], sep='\t')
    scores = scores.loc[scores.GENE.ne('Non-Targeting-Control')].set_index('GENE')
    top = scores.sort_values('rank_synth').head(20)
    GENES[:] = top.index.tolist()
    assert top.rank_synth.tolist() == list(range(1, 21)) and top.fdr_synth.lt(.05).all()
    selected = counts.loc[counts.Gene.isin(GENES)].sort_values(['Gene', 'sgRNA']).copy()
    selected['guide_number'] = selected.groupby('Gene').cumcount() + 1
    library = pd.read_csv(SOURCES['library'], sep='\t')
    pairs = set(zip(library.gene, library.Seq))
    assert selected.groupby('Gene').size().eq(4).all()
    assert all((r.Gene, r.sgRNA.rsplit('_', 1)[-1]) in pairs for r in selected.itertuples())
    records = []
    for idx, row in selected.iterrows():
        for sample in SAMPLES:
            records.append(dict(gene=row.Gene, sgrna=row.sgRNA, guide_number=row.guide_number,
                                sample=sample, genotype='CTRL' if 'CTRL' in sample else 'LOY',
                                stage='D0' if sample.startswith('D0') else 'Endpoint',
                                replicate=0 if sample.startswith('D0') else int(sample[-1]),
                                raw_count=int(row[sample]), cpm=float(cpm.loc[idx, sample]),
                                plotted_cpm_plus_1=float(cpm.loc[idx, sample] + 1)))
    all_values = pd.DataFrame(records)
    all_values['display_included'] = all_values.sgrna.ne(EXCLUDED_GUIDE)
    assert len(all_values) == 480
    excluded = all_values.loc[~all_values.display_included].copy()
    assert len(excluded) == 6 and excluded.guide_number.eq(2).all()
    assert excluded.gene.eq('DDX3X').all()
    excluded['exclusion_reason'] = 'Author-requested display exclusion of unusually abundant DDX3X D0 guide; same guide omitted at both stages in both genotypes.'
    excluded.to_csv(DATA / 'excluded_guide_observations.csv', index=False)
    all_values.to_csv(DATA / 'all_guide_abundances.csv', index=False)
    values = all_values.loc[all_values.display_included].copy()
    assert len(values) == 474 and values.sgrna.nunique() == 79
    assert values.plotted_cpm_plus_1.between(0, 80).all()
    existing = pd.read_csv(SOURCES['existing_profiles']).set_index(['sgrna', 'sample'])
    common = values.set_index(['sgrna', 'sample']).index.intersection(existing.index)
    check = values.set_index(['sgrna', 'sample']).loc[common]
    assert np.array_equal(check.raw_count, existing.loc[common, 'raw_count'])
    assert np.allclose(check.cpm, existing.loc[common, 'cpm'], rtol=1e-13, atol=1e-13)
    values.to_csv(DATA / 'plotted_guide_abundances.csv', index=False)
    key = selected[['Gene', 'sgRNA', 'guide_number']].copy()
    key['display_included'] = key.sgRNA.ne(EXCLUDED_GUIDE)
    key.to_csv(DATA / 'guide_key.csv', index=False)
    scores.loc[GENES].to_csv(DATA / 'displayed_gene_scores.csv')

    fig = plt.figure(figsize=(12.4, 12.8))
    handles = [Line2D([], [], marker='s', color='#666666', markerfacecolor='white',
                      lw=0, ms=4.5, label='D0'),
               Line2D([], [], marker='o', color='#666666', lw=0, ms=4.5, label='Endpoint R1'),
               Line2D([], [], marker='^', color='#666666', lw=0, ms=4.5, label='Endpoint R2'),
               Line2D([], [], color='#888888', lw=.9, label='Same guide: D0 to endpoint mean')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.51, .998),
               ncol=4, fontsize=8, handletextpad=.5, columnspacing=1.8)
    endpoints = []
    for i, gene in enumerate(GENES):
        row, col = divmod(i, 4)
        left = .055 + col * .239
        bottom = (10.48 - row * 2.38) / 12.8
        ax = fig.add_axes([left, bottom, .21, 1.40 / 12.8])
        fdr = float(scores.loc[gene, 'fdr_synth'])
        label = f'{i + 1}. {gene}' + ('*' if gene == 'DDX3X' else '')
        ax.text(.5, 1.30, label, ha='center', transform=ax.transAxes,
                fontsize=9.5, fontweight='semibold')
        ax.text(.5, 1.15, f'Endpoint drugZ FDR = {fdr:.3g}', ha='center',
                transform=ax.transAxes, fontsize=8, color='#555555')
        for genotype, start in [('CTRL', 0), ('LOY', 2.4)]:
            color = COLORS[genotype]
            sub = values.loc[values.gene.eq(gene) & values.genotype.eq(genotype)]
            for number, guide in sub.groupby('guide_number'):
                d0 = float(guide.loc[guide.stage.eq('D0'), 'plotted_cpm_plus_1'].iloc[0])
                ep = guide.loc[guide.stage.eq('Endpoint')].sort_values('replicate')
                mean = float(ep.cpm.mean() + 1)
                offset = (int(number) - 2.5) * .055
                ax.plot([start + offset, start + 1 + offset], [d0, mean],
                        color=color, alpha=.36, lw=.9, zorder=1)
                ax.scatter(start + offset, d0, s=21, marker='s', facecolors='white',
                           edgecolors=color, linewidths=.8, zorder=3)
                for r in ep.itertuples():
                    x = start + 1 + offset + (-.035 if r.replicate == 1 else .035)
                    ax.scatter(x, r.plotted_cpm_plus_1, s=22,
                               marker='o' if r.replicate == 1 else '^',
                               color=color, edgecolors='white', linewidths=.3, zorder=4)
                endpoints.append(dict(gene=gene, genotype=genotype, guide_number=int(number),
                                      D0_cpm_plus_1=d0, endpoint_mean_cpm_plus_1=mean))
            ax.text(start + .5, 76, genotype, ha='center', va='center',
                    fontsize=8.5, color=color, clip_on=False)
        ax.axhline(1, color='#DADADA', lw=.55, zorder=0)
        ax.axvline(1.7, color='#EEEEEE', lw=.7, zorder=0)
        ax.set_xlim(-.3, 3.7)
        ax.set_ylim(0, 80)
        ax.set_yticks([0, 20, 40, 60, 80])
        ax.set_xticks([0, 1, 2.4, 3.4], ['D0', 'Endpoint', 'D0', 'Endpoint'])
        ax.tick_params(labelsize=8, length=2.5, pad=3)
        if col == 0:
            ax.set_ylabel('Guide abundance (CPM + 1)', fontsize=8)
        else:
            ax.tick_params(axis='y', labelleft=False)
    fig.text(.5, .033, 'One D0 sample per genotype; two endpoint cultures. Lines connect each guide to its endpoint mean.',
             ha='center', fontsize=8, color='#555555')
    fig.text(.5, .016, '* DDX3X guide 2 omitted from display at both stages. CPM denominators and original drugZ FDRs are unchanged.',
             ha='center', fontsize=8, color='#555555')
    pd.DataFrame(endpoints).to_csv(DATA / 'guide_connection_values.csv', index=False)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600 if ext == 'png' else 300)
    fig.savefig(OUT / 'review' / f'{STEM}_preview.png', dpi=160)
    plt.close(fig)

    with fitz.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        text = page.get_text()
        assert all(g in text for g in GENES)
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        overlaps = []
        for i, a in enumerate(spans):
            assert page.rect.contains(fitz.Rect(a['bbox'])), a['text']
            for b in spans[i + 1:]:
                intersection = fitz.Rect(a['bbox']) & fitz.Rect(b['bbox'])
                if not intersection.is_empty and intersection.width > .5 and intersection.height > .5:
                    overlaps.append([a['text'], b['text']])
        assert not overlaps, overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.6, 1.6), colorspace=fitz.csGRAY).save(
            OUT / 'review' / f'{STEM}_greyscale.png')
    manifest = {name: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                for name, path in SOURCES.items()}
    manifest['script'] = {'path': str(Path(__file__)), 'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    summary = dict(genes=GENES, n_guides=79, D0_observations=158, endpoint_observations=316,
                   connections=158, existing_profile_values_verified=len(common), scale='linear CPM + 1',
                   shared_ylim=[0, 80], all_guides_retained=False, excluded_guide=EXCLUDED_GUIDE,
                   excluded_observations=6, source_guide_count=80,
                   filtering='Display-only exclusion of DDX3X guide 2 in both stages and genotypes; CPM still uses all 77,441 guides in each sample.',
                   baseline_samples_per_genotype=1, endpoint_cultures_per_genotype=2,
                   annotations='Saved endpoint drugZ directional depletion FDR, not D0-to-endpoint tests.',
                   interpretation='Descriptive relative abundance; D0 library composition is distorted by honeypotting. No absolute growth interpretation.',
                   text_overlaps=0, font='Open Sans', integrated_into_manuscript=False)
    (DATA / 'validation.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
