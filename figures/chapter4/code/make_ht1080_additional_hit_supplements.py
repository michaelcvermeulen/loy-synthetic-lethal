"""Three standalone Section 4.3.1.2 supplements from frozen local outputs.

No screen refits, new enrichment tests, changes to hit calls or guide exclusions.
"""
from pathlib import Path
import hashlib
import json
import os
import sys
import textwrap

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_additional_hit_supplements'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import pymupdf as fitz
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8.5)

SCREEN = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
DEP = OUT / 'source_data/ht1080_depmap_codependency'
CTRL, LOY = '#3C5488', '#E64B35'
STRONG, EXPLORATORY, NS = '#7A3E85', '#BC8A26', '#A7A7A7'
SOURCES = {
    'drugZ': SCREEN / 'results/drugZ/LOY_vs_CTRL_unfiltered.txt',
    'mageck_genes': SCREEN / 'test/POST_LOY_vs_CTRL.gene_summary.txt',
    'mageck_guides': SCREEN / 'test/POST_LOY_vs_CTRL.sgrna_summary.txt',
    'counts': SCREEN / 'count/ht1080_screen2_R1.count.txt',
    'partners': DEP / 'top_partners.csv',
    'pairwise': DEP / 'seed_pairwise_pearson.csv',
    'pairwise_n': DEP / 'seed_pairwise_n.csv',
    'depmap_manifest': DEP / 'analysis_manifest.json',
}
DZ = pd.read_csv(SOURCES['drugZ'], sep='\t').set_index('GENE')
MAG = pd.read_csv(SOURCES['mageck_genes'], sep='\t')
DEPLETED = set(DZ.index[DZ.fdr_synth.lt(.25)])
ENRICHED = DZ.loc[DZ.fdr_supp.lt(.25)].sort_values('normZ', ascending=False)
SCREEN_ALIASES = {'SPRING1': 'C12orf49', 'HAPSTR1': 'C16orf72', 'VIRMA': 'KIAA1429'}
assert len(DEPLETED) == 48 and len(ENRICHED) == 9
RECORD = {}


def colour(q):
    return STRONG if q < .05 else EXPLORATORY if q < .25 else NS


def fdr(gene):
    r = DZ.loc[SCREEN_ALIASES.get(gene, gene)]
    return float(r.fdr_synth if r.normZ < 0 else r.fdr_supp)


def label(fig, letter, x, y):
    fig.text(x, y, letter, fontsize=12, fontweight='semibold', va='top')


def save(fig, stem, required):
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{stem}.{ext}', dpi=600 if ext == 'png' else 300)
    fig.savefig(OUT / 'review' / f'{stem}_preview.png', dpi=160)
    plt.close(fig)
    with fitz.open(OUT / f'{stem}.pdf') as doc:
        page = doc[0]
        assert all(s in page.get_text() for s in required)
        assert all('OpenSans' in f[3] for f in page.get_fonts()), page.get_fonts()
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside, overlaps = [], []
        for i, a in enumerate(spans):
            if not page.rect.contains(fitz.Rect(a['bbox'])):
                outside.append(a['text'])
            for b in spans[i + 1:]:
                rect = fitz.Rect(a['bbox']) & fitz.Rect(b['bbox'])
                if not rect.is_empty and rect.width > .5 and rect.height > .5:
                    overlaps.append([a['text'], b['text']])
        page.get_pixmap(matrix=fitz.Matrix(1.6, 1.6), colorspace=fitz.csGRAY).save(
            OUT / 'review' / f'{stem}_greyscale.png')
        RECORD[stem] = dict(page_inches=[page.rect.width / 72, page.rect.height / 72],
                           outside_text=outside, text_overlaps=overlaps, font='Open Sans')


def metascape():
    runs = {}
    for name, rel, cutoff in [('strict', 'depleted', .25),
                             ('broad', 'depleted/FDR_less_than_0.4/all.t1kkj1m1w', .4)]:
        root = SCREEN / 'results/metascape' / rel
        for key, path in [('terms', root / 'Enrichment_GO/_FINAL_GO.csv'),
                          ('annotation', root / 'metascape_result.xlsx'),
                          ('report', root / 'AnalysisReport.html')]:
            SOURCES[f'metascape_{name}_{key}'] = path
        a = pd.read_excel(root / 'metascape_result.xlsx', sheet_name='Annotation',
                          usecols=['MyList', 'Gene ID', 'Gene Symbol'])
        genes = set(a.MyList)
        assert genes == set(MAG.loc[MAG['neg|fdr'].lt(cutoff), 'id'])
        a['drugZ_depleted_fdr_025'] = a.MyList.isin(DEPLETED)
        a.to_csv(DATA / f'metascape_{name}_input.csv', index=False)
        t = pd.read_csv(root / 'Enrichment_GO/_FINAL_GO.csv')
        t['q_value'] = 10 ** t['Log(q-value)']
        t['p_value'] = 10 ** t.LogP
        converted_to_input = dict(zip(a['Gene Symbol'], a.MyList))
        t['input_gene_symbols'] = t.Hits.apply(lambda x: '|'.join(converted_to_input[g] for g in x.split('|')))
        t['drugZ_depleted_overlap'] = t.input_gene_symbols.apply(lambda x: len(set(x.split('|')) & DEPLETED))
        assert t['#GeneInHitList'].eq(len(genes)).all()
        assert t['#TotalGeneInLibrary'].eq(30238).all()
        assert all(set(s.split('|')) <= genes for s in t.input_gene_symbols)
        t.to_csv(DATA / f'metascape_{name}_all_terms.csv', index=False)
        reps = t.loc[t.FirstInGroupByLogP.eq(1)].sort_values('LogP').copy()
        reps.to_csv(DATA / f'metascape_{name}_all_cluster_representatives.csv', index=False)
        runs[name] = reps
        RECORD[f'metascape_{name}'] = dict(input_n=len(genes), mageck_depletion_fdr=cutoff,
            primary_drugZ_overlap=len(genes & DEPLETED), term_count=len(t),
            cluster_count=len(reps), q_lt_005_representatives=int(reps.q_value.lt(.05).sum()),
            background_n=30238, database_updated='2025-07-01')

    fig = plt.figure(figsize=(10.2, 10.7))
    for letter, name, bottom, height, number in [('a', 'strict', .745, .195, 5),
                                               ('b', 'broad', .362, .305, 10)]:
        sub = runs[name].head(number).copy()
        sub.to_csv(DATA / f'metascape_panel_{letter}_terms.csv', index=False)
        title = ('MAGeCK depletion FDR < 0.25: 23 genes (21 also drugZ hits)' if name == 'strict'
                 else 'MAGeCK depletion FDR < 0.40: 86 genes (40 also drugZ hits)')
        label(fig, letter, .025, bottom + height + .044)
        fig.text(.068, bottom + height + .04, title, fontsize=9.5, fontweight='semibold', va='top')
        ax = fig.add_axes([.47, bottom, .33, height])
        yy = np.arange(len(sub))
        vals = -sub['Log(q-value)'].to_numpy()
        ax.axvline(-np.log10(.05), lw=.8, ls='--', color='#888888', zorder=0)
        ax.hlines(yy, 0, vals, color='#D5D5D5', lw=1.7)
        for y, value, (_, r) in zip(yy, vals, sub.iterrows()):
            ax.scatter(value, y, s=30, facecolor=STRONG if r.q_value < .05 else 'white',
                       edgecolor=STRONG if r.q_value < .05 else '#999999', lw=.9, zorder=3)
            ax.text(-.035, y, textwrap.fill(r.Description.strip(), 48), ha='right', va='center',
                    transform=ax.get_yaxis_transform(), fontsize=8.5,
                    fontweight='semibold' if r.q_value < .05 else 'normal')
            for pos, value in [(1.10, int(r['#GeneInGOAndHitList'])),
                               (1.29, int(r.drugZ_depleted_overlap))]:
                ax.text(pos, y, str(value), transform=ax.get_yaxis_transform(),
                        ha='center', va='center', fontsize=8.5)
        ax.text(1.10, 1.05, 'Input\nhits', transform=ax.transAxes, ha='center', fontsize=8)
        ax.text(1.29, 1.05, 'drugZ\nhits', transform=ax.transAxes, ha='center', fontsize=8)
        ax.set(xlim=(-.13, 8), ylim=(len(sub) - .5, -.5), yticks=[], xticks=[0, 2, 4, 6, 8],
               xlabel='−log₁₀(Metascape q value)')
        ax.spines['left'].set_visible(False)

    reps = runs['broad'].loc[runs['broad'].q_value.lt(.05)]
    genes = list(dict.fromkeys(g for hits in reps.Hits for g in hits.split('|')))
    membership = pd.DataFrame([[int(g in r.Hits.split('|')) for g in genes]
                              for r in reps.itertuples()], index=reps.GO, columns=genes)
    membership.to_csv(DATA / 'metascape_panel_c_membership.csv')
    DZ.loc[genes].to_csv(DATA / 'metascape_panel_c_gene_scores.csv')
    label(fig, 'c', .025, .291)
    fig.text(.068, .287, 'Genes contributing to broad-list clusters with q < 0.05',
             fontsize=9.5, fontweight='semibold', va='top')
    ax = fig.add_axes([.35, .163, .61, .100])
    for y, (_, row) in enumerate(membership.iterrows()):
        ax.axhspan(y - .45, y + .45, color='#F6F6F6', zorder=0)
        for x, g in enumerate(genes):
            if row[g]:
                ax.scatter(x, y, s=35, color=colour(float(DZ.loc[g, 'fdr_synth'])),
                           edgecolor='white', linewidth=.3)
    names = ['STAGA complex', 'Stem-cell population maintenance',
             'SREBF / miR33 lipid homeostasis', 'Epigenetic gene regulation']
    ax.set(xticks=np.arange(len(genes)), xticklabels=genes, yticks=np.arange(len(reps)),
           yticklabels=names, xlim=(-.5, len(genes) - .5), ylim=(len(reps) - .5, -.5))
    ax.tick_params(axis='both', length=0, labelsize=8)
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center')
    for spine in ax.spines.values(): spine.set_visible(False)
    for tick, g in zip(ax.get_xticklabels(), genes):
        if DZ.loc[g, 'fdr_synth'] < .05: tick.set_fontweight('semibold')
    handles = [Line2D([], [], marker='o', color=c, lw=0, markersize=5, label=txt)
               for c, txt in [(STRONG, 'drugZ depletion FDR < 0.05'),
                              (EXPLORATORY, '0.05 ≤ FDR < 0.25'), (NS, 'FDR ≥ 0.25')]]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .038),
               ncol=3, fontsize=8, columnspacing=1.5, handletextpad=.4)
    fig.text(.5, .019, 'Saved Metascape runs; genome background = 30,238 genes. Dashed line: q = 0.05.',
             ha='center', fontsize=8, color='#555555')
    save(fig, 'ch04_supp_ht1080_metascape_overrepresentation_draft', ['STAGA', 'MAGeCK', 'SREBF1'])


def neighbourhoods():
    partners = pd.read_csv(SOURCES['partners'])
    pair = pd.read_csv(SOURCES['pairwise'], index_col=0)
    n = pd.read_csv(SOURCES['pairwise_n'], index_col=0)
    seeds = ['SREBF1', 'SREBF2', 'WTAP', 'UBE2G2', 'ITGB5', 'USP7']
    comparisons = {'SREBF1': 'SREBF2', 'SREBF2': 'SREBF1', 'USP7': 'TP53'}
    records = []
    for seed in seeds:
        rows = [dict(gene=seed, pearson_r=np.nan, n_models=np.nan, selection='Seed')]
        top = partners.loc[partners.seed.eq(seed) & partners.selection.eq('Top positive')].sort_values('rank_positive').head(5)
        assert len(top) == 5
        rows += [dict(gene=r.partner, pearson_r=r.pearson_r, n_models=int(r.n_models),
                      selection='Top five positive') for r in top.itertuples()]
        if seed in comparisons:
            g = comparisons[seed]
            rows.append(dict(gene=g, pearson_r=float(pair.loc[seed, g]), n_models=int(n.loc[seed, g]),
                             selection='Selected comparison'))
        for row in rows:
            screen_gene = SCREEN_ALIASES.get(row['gene'], row['gene'])
            row.update(seed=seed, screen_gene=screen_gene, normZ=float(DZ.loc[screen_gene, 'normZ']),
                       directional_fdr=fdr(row['gene']))
        records.extend(rows)
    table = pd.DataFrame(records)
    table.to_csv(DATA / 'hit_neighbourhoods_plotted_values.csv', index=False)
    fig = plt.figure(figsize=(10.2, 8.7))
    for i, seed in enumerate(seeds):
        rr, cc = divmod(i, 2)
        left = .10 + cc * .49
        bottom = .715 - rr * .295
        label(fig, chr(97 + i), left - .075, bottom + .244)
        fig.text(left - .025, bottom + .24, seed, fontsize=10, fontweight='semibold', va='top')
        sub = table.loc[table.seed.eq(seed)].reset_index(drop=True)
        rax = fig.add_axes([left, bottom, .138, .199])
        zax = fig.add_axes([left + .215, bottom, .112, .199])
        for j, r in sub.iterrows():
            if r.selection == 'Selected comparison':
                for ax in [rax, zax]: ax.axhspan(j - .43, j + .43, color='#F2F2F2', zorder=0)
            if np.isfinite(r.pearson_r):
                rax.hlines(j, 0, r.pearson_r, color='#C9C9C9', lw=1)
                rax.scatter(r.pearson_r, j, s=22, color='#525252',
                            marker='D' if r.selection == 'Selected comparison' else 'o', zorder=3)
                rax.text(1.06, j, f'{r.pearson_r:.2f}', transform=rax.get_yaxis_transform(),
                         va='center', fontsize=8)
            zax.hlines(j, 0, r.normZ, color='#D5D5D5', lw=1)
            zax.scatter(r.normZ, j, s=27, color=colour(r.directional_fdr), zorder=3)
            zax.text(1.12, j, f'{r.directional_fdr:.3g}', transform=zax.get_yaxis_transform(),
                     va='center', fontsize=8)
        rax.axvline(0, color='#BBBBBB', lw=.6, zorder=0)
        zax.axvline(0, color='#BBBBBB', lw=.6, zorder=0)
        rax.set(xlim=(-.55, .74), xticks=[-.5, 0, .5], ylim=(6.6, -.6),
                 yticks=np.arange(len(sub)), yticklabels=[g + ('*' if g == 'ETNK1' else '') for g in sub.gene],
                 xlabel='DepMap Pearson r')
        rax.tick_params(axis='y', length=0, pad=5, labelsize=8)
        for tick, q in zip(rax.get_yticklabels(), sub.directional_fdr):
            if q < .05: tick.set_fontweight('semibold')
        zax.set(xlim=(-6.4, 6.4), xticks=[-6, 0, 6], ylim=(6.6, -.6), yticks=[], xlabel='drugZ normZ')
        zax.text(1.12, 1.035, 'FDR', transform=zax.transAxes, fontsize=8)
        for ax in [rax, zax]:
            ax.spines['left'].set_visible(False)
            ax.tick_params(axis='x', labelsize=8, pad=3)
    handles = [Line2D([], [], color=c, marker='o', lw=0, markersize=5, label=t)
               for c, t in [(STRONG, 'drugZ FDR < 0.05'), (EXPLORATORY, '0.05 ≤ FDR < 0.25'),
                             (NS, 'FDR ≥ 0.25')]]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.50, .029), ncol=3,
               fontsize=8, columnspacing=1.8, handletextpad=.4)
    fig.text(.5, .014, 'Top five positive partners; shaded diamonds are selected comparisons. n = 1,186 models (* ETNK1: n = 355).',
             ha='center', fontsize=7.8, color='#555555')
    RECORD['neighbourhoods'] = dict(seeds=seeds, rows=len(table), top_partners=30,
        selected_comparisons=3, partners_meeting_drugz_fdr_025=int(table.loc[table.selection.eq('Top five positive'), 'directional_fdr'].lt(.25).sum()),
        srebs_pearson_r=float(pair.loc['SREBF1', 'SREBF2']),
        usp7_tp53_pearson_r=float(pair.loc['USP7', 'TP53']),
        method='Saved Pearson correlations, unadjusted, DepMap Public 25Q3; not a LOY-specificity test')
    save(fig, 'ch04_supp_ht1080_additional_hit_neighbourhoods_draft', seeds + ['TP53', 'ETNK1*'])


def enriched():
    counts = pd.read_csv(SOURCES['counts'], sep='\t')
    samples = list(counts.columns[2:])
    assert len(counts) == 77441 and counts.sgRNA.is_unique
    cpm = counts[samples].div(counts[samples].sum(), axis=1) * 1e6
    guides = pd.read_csv(SOURCES['mageck_guides'], sep='\t')
    genes = list(ENRICHED.index)
    selected = counts.loc[counts.Gene.isin(genes)].sort_values(['Gene', 'sgRNA']).copy()
    selected['guide_number'] = selected.groupby('Gene').cumcount() + 1
    assert len(selected) == 36 and selected.groupby('Gene').size().eq(4).all()
    rows = []
    for idx, row in selected.iterrows():
        for sample in samples:
            rows.append(dict(gene=row.Gene, sgrna=row.sgRNA, guide_number=int(row.guide_number),
                sample=sample, genotype='CTRL' if 'CTRL' in sample else 'LOY',
                stage='D0' if sample.startswith('D0') else 'Endpoint',
                replicate=0 if sample.startswith('D0') else int(sample[-1]),
                raw_count=int(row[sample]), cpm=float(cpm.loc[idx, sample]),
                cpm_plus_1=float(cpm.loc[idx, sample] + 1)))
    values = pd.DataFrame(rows)
    assert len(values) == 216 and values.cpm_plus_1.between(0, 70).all()
    values.to_csv(DATA / 'enriched_hits_guide_abundances.csv', index=False)
    selected[['Gene', 'sgRNA', 'guide_number']].to_csv(DATA / 'enriched_hits_guide_key.csv', index=False)
    ENRICHED.to_csv(DATA / 'enriched_hit_scores.csv')
    effects = guides.loc[guides.Gene.isin(genes)].copy()
    assert len(effects) == 36 and effects.sgrna.is_unique
    assert set(effects.sgrna) == set(selected.sgRNA)
    effects.to_csv(DATA / 'enriched_hits_saved_guide_lfc.csv', index=False)
    ntc = guides.loc[guides.Gene.eq('Non-Targeting-Control'), ['sgrna', 'LFC']]
    assert len(ntc) == 1000
    ntc.to_csv(DATA / 'enriched_panel_b_ntc_lfc.csv', index=False)
    # Gaussian KDE, Scott bandwidth; only for descriptive display of the NTC distribution.
    x = np.linspace(-4, 4, 801)
    null = ntc.LFC.to_numpy()
    bw = float(np.std(null, ddof=1) * len(null) ** (-.2))
    density = np.exp(-.5 * ((x[:, None] - null[None, :]) / bw) ** 2).mean(axis=1) / (bw * np.sqrt(2 * np.pi))
    pd.DataFrame({'lfc': x, 'density': density}).to_csv(DATA / 'enriched_panel_b_null_density.csv', index=False)
    fig = plt.figure(figsize=(10.2, 11.5))
    label(fig, 'a', .025, .98)
    label(fig, 'b', .46, .98)
    ax = fig.add_axes([.13, .715, .21, .237])
    y = np.arange(9)
    ax.hlines(y, 0, ENRICHED.normZ, color='#D7D7D7', lw=1)
    ax.scatter(ENRICHED.normZ, y, c=[colour(q) for q in ENRICHED.fdr_supp], s=30, zorder=3)
    ax.set(yticks=y, yticklabels=genes, ylim=(8.5, -.5), xlim=(0, 5.1), xticks=[0, 2, 4], xlabel='drugZ normZ (enrichment)')
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    for j, q in enumerate(ENRICHED.fdr_supp):
        ax.text(1.05, j, f'{q:.3g}', transform=ax.get_yaxis_transform(), va='center', fontsize=8)
    ax.text(1.05, 1.03, 'FDR', transform=ax.transAxes, fontsize=8)
    for tick, q in zip(ax.get_yticklabels(), ENRICHED.fdr_supp):
        if q < .05: tick.set_fontweight('semibold')
    ax = fig.add_axes([.58, .715, .36, .237])
    for i, gene in enumerate(genes):
        curve = .64 * density / density.max()
        ax.fill_between(x, i + .33, i + .33 - curve, color='#E5E5E5')
        ax.plot(x, i + .33 - curve, color='#9E9E9E', lw=.75)
        vals = effects.loc[effects.Gene.eq(gene), 'LFC']
        ax.vlines(vals, i - .3, i + .32, color=colour(float(ENRICHED.loc[gene, 'fdr_supp'])), lw=1)
    ax.axvline(0, color='#AAAAAA', ls=':', lw=.7)
    ax.set(xlim=(-4, 4), xticks=[-4, -2, 0, 2, 4], yticks=y, yticklabels=genes,
           ylim=(8.55, -.55), xlabel='Endpoint guide log₂ fold change (LOY / sgCTRL)')
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    fig.text(.76, .966, 'Guide ticks over 1,000 non-targeting controls', ha='center', fontsize=8)
    label(fig, 'c', .025, .659)
    handles = [Line2D([], [], color=CTRL, lw=1, label='sgCTRL'),
               Line2D([], [], color=LOY, lw=1, label='Y-depleted'),
               Line2D([], [], marker='s', color='#666666', markerfacecolor='white', lw=0, ms=4, label='D0'),
               Line2D([], [], marker='o', color='#666666', lw=0, ms=4, label='Endpoint R1'),
               Line2D([], [], marker='^', color='#666666', lw=0, ms=4, label='Endpoint R2')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.52, .665), ncol=5,
               fontsize=8, columnspacing=1.9, handletextpad=.4)
    connections = []
    for i, gene in enumerate(genes):
        rr, cc = divmod(i, 3)
        left, bottom = .075 + cc * .315, .467 - rr * .197
        ax = fig.add_axes([left, bottom, .267, .132])
        ax.text(.5, 1.19, gene, ha='center', transform=ax.transAxes, fontsize=9,
                fontweight='semibold' if ENRICHED.loc[gene, 'fdr_supp'] < .05 else 'normal')
        for xpos, name, col in [(.20, 'sgCTRL', CTRL), (.80, 'Y-depleted', LOY)]:
            ax.text(xpos, 1.02, name, transform=ax.transAxes, ha='center', va='bottom',
                    fontsize=7.5, color=col)
        for gt, start in [('CTRL', 0), ('LOY', 2.4)]:
            c = CTRL if gt == 'CTRL' else LOY
            sub = values.loc[values.gene.eq(gene) & values.genotype.eq(gt)]
            for gn, guide in sub.groupby('guide_number'):
                d0 = float(guide.loc[guide.stage.eq('D0'), 'cpm_plus_1'].iloc[0])
                ep = guide.loc[guide.stage.eq('Endpoint')].sort_values('replicate')
                mean = float(ep.cpm_plus_1.mean())
                offset = (int(gn) - 2.5) * .055
                ax.plot([start + offset, start + 1 + offset], [d0, mean], color=c, lw=.8, alpha=.40)
                ax.scatter(start + offset, d0, s=18, marker='s', facecolors='white', edgecolors=c, linewidths=.7, zorder=3)
                for r in ep.itertuples():
                    ax.scatter(start + 1 + offset + (-.035 if r.replicate == 1 else .035),
                               r.cpm_plus_1, marker='o' if r.replicate == 1 else '^',
                               color=c, s=19, edgecolor='white', lw=.25, zorder=4)
                connections.append(dict(gene=gene, genotype=gt, guide_number=int(gn),
                                        d0_cpm_plus_1=d0, endpoint_mean_cpm_plus_1=mean))
        ax.axhline(1, lw=.5, color='#DADADA', zorder=0)
        ax.axvline(1.7, lw=.5, color='#EEEEEE', zorder=0)
        ax.set(xlim=(-.3, 3.7), ylim=(0, 70), yticks=[0, 20, 40, 60],
               xticks=[0, 1, 2.4, 3.4], xticklabels=['D0', 'Endpoint', 'D0', 'Endpoint'])
        ax.tick_params(axis='both', labelsize=8, length=2.5, pad=3)
        if cc == 0: ax.set_ylabel('Guide abundance (CPM + 1)', fontsize=8)
        else: ax.tick_params(axis='y', labelleft=False)
    pd.DataFrame(connections).to_csv(DATA / 'enriched_hits_guide_connections.csv', index=False)
    fig.text(.5, .023, 'One D0 sample and two endpoint cultures per genotype. Lines connect each guide to its endpoint mean.',
             ha='center', fontsize=8, color='#555555')
    RECORD['enriched_hits'] = dict(genes=genes, guides=36, abundance_observations=216,
        baseline_observations=72, endpoint_observations=144, connections=72,
        ntc_guides=1000, ntc_kde_bandwidth=bw, kde_display_range=[-4, 4],
        ntc_outside_display=int(np.sum((null < -4) | (null > 4))),
        guide_lfc_range=[float(effects.LFC.min()), float(effects.LFC.max())],
        max_cpm_plus_1=float(values.cpm_plus_1.max()), raw_zero_count=int(values.raw_count.eq(0).sum()),
        denominator='All 77,441 guides per sample; no exclusions or renormalisation',
        guide_lfc='Saved MAGeCK endpoint LFC; not drugZ normZ or a new D0-corrected score')
    save(fig, 'ch04_supp_ht1080_enriched_hit_guides_draft', genes)


def main():
    metascape()
    neighbourhoods()
    enriched()
    aliases = [dict(display_gene=g, screen_gene=s, source=url) for g, s, url in [
        ('SPRING1', 'C12orf49', 'https://www.ncbi.nlm.nih.gov/gene/79794'),
        ('HAPSTR1', 'C16orf72', 'https://www.ncbi.nlm.nih.gov/gene/29035'),
        ('VIRMA', 'KIAA1429', 'https://www.ncbi.nlm.nih.gov/gene/25962')]]
    pd.DataFrame(aliases).to_csv(DATA / 'verified_screen_symbol_aliases.csv', index=False)
    # Ready-to-use inputs for a future Metascape rerun, not yet submitted or analysed.
    for name, genes in [('drugz_depleted_48', DZ.loc[DZ.fdr_synth.lt(.25)].sort_values('normZ').index),
                        ('drugz_enriched_9', ENRICHED.index),
                        ('drugz_tested_gene_background', DZ.index[DZ.index != 'Non-Targeting-Control'])]:
        (DATA / f'future_metascape_{name}.txt').write_text('\n'.join(genes) + '\n', encoding='utf-8')
    SOURCES['script'] = Path(__file__)
    manifest = {key: dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                for key, p in SOURCES.items()}
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    (DATA / 'validation.json').write_text(json.dumps(RECORD, indent=2), encoding='utf-8')
    for name, record in RECORD.items():
        if name.startswith('ch04_'):
            assert not record['outside_text'] and not record['text_overlaps'], (name, record)
    print(json.dumps(RECORD, indent=2))


if __name__ == '__main__':
    main()
