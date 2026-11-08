"""Restyle the original full-48 Metascape results, excluding GPEC rescoring.

All numerical enrichment results are from GO_AllLists.csv. Final-table cluster
assignments are used only to organise the display. No enrichment tests are run.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import sys
import textwrap

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_drugz_metascape_48'
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
apply_style(8.5)

ROOT = THESIS / 'research/chapter4_metascape_drugz_48_2026-10-18'
REPORT = ROOT / 'report'
STEM = 'ch04_supp_ht1080_drugz_metascape_draft'
SOURCES = {
    'zip': ROOT / 'all.t4l6xpy52.zip',
    'annotation': REPORT / 'metascape_result.xlsx',
    'original_full_list_enrichment': REPORT / 'Enrichment_GO/GO_AllLists.csv',
    'final_table_for_cluster_ids_only': REPORT / 'Enrichment_GO/_FINAL_GO.csv',
    'evidence_counts': REPORT / 'Evidence/Evidence.csv',
    'report': REPORT / 'AnalysisReport.html',
    'drugz': THESIS.parent / 'crispr_screen/ht1080_loy_screen_2/results/drugZ/LOY_vs_CTRL_unfiltered.txt',
    'prepared_input': OUT / 'source_data/ht1080_additional_hit_supplements/future_metascape_drugz_depleted_48.txt',
    'prepared_background': OUT / 'source_data/ht1080_additional_hit_supplements/future_metascape_drugz_tested_gene_background.txt',
    'script': Path(__file__),
}
STRONG, EXPLORATORY = '#7A3E85', '#BC8A26'


def main():
    scores = pd.read_csv(SOURCES['drugz'], sep='\t').set_index('GENE')
    expected = set(scores.index[scores.fdr_synth.lt(.25)])
    annotation = pd.read_excel(SOURCES['annotation'], sheet_name='Annotation',
                               usecols=['MyList', 'Gene ID', 'Type', 'Tax ID'])
    assert len(annotation) == 48 and annotation.MyList.is_unique and annotation['Gene ID'].is_unique
    assert set(annotation.MyList) == expected == set(SOURCES['prepared_input'].read_text().splitlines())
    annotation.rename(columns={'MyList': 'input_symbol', 'Gene ID': 'entrez_id'}, inplace=True)
    annotation.merge(scores.reset_index(), left_on='input_symbol', right_on='GENE',
                     validate='one_to_one').to_csv(DATA / 'input_48_verified.csv', index=False)
    gene_map = dict(zip(annotation.entrez_id.astype(str), annotation.input_symbol))
    original = pd.read_csv(SOURCES['original_full_list_enrichment'])
    final = pd.read_csv(SOURCES['final_table_for_cluster_ids_only'])
    evidence = pd.read_csv(SOURCES['evidence_counts'])
    assert len(original) == 133 and original.GO.is_unique and final.GO.is_unique
    assert original['#GeneInHitList'].eq(48).all()
    assert original['#TotalGeneInLibrary'].eq(18975).all()
    assert set(original.GeneList) == {'MyList'} and set(original.GO) == set(final.GO)
    original['source_row'] = np.arange(len(original))
    original['input_gene_symbols'] = original.GeneID.apply(lambda x: '|'.join(gene_map[g] for g in x.split('|')))
    assert original.apply(lambda r: len(r.input_gene_symbols.split('|')) == r['#GeneInGOAndHitList'], axis=1).all()
    original['approx_q_value'] = 10 ** original['Log(q-value)']
    # q values are approximate because this original table rounds log10(q).
    # Values printed as -1.3 cannot resolve which side of 0.05 the exact q lies.
    original['q_boundary_rounding_ambiguous'] = original['Log(q-value)'].eq(-1.3)
    original['reported_q_clearly_lt_005'] = original['Log(q-value)'].lt(-1.35)
    original = original.merge(final[['GO', 'GROUP_ID']], on='GO', validate='one_to_one')
    original.to_csv(DATA / 'full_48_original_enrichment.csv', index=False)
    # Preserve source ranking when rounded P values are tied.
    reps = original.sort_values(['LogP', 'source_row'], kind='stable').drop_duplicates('GROUP_ID').copy()
    assert len(reps) == 13
    assert not reps.q_boundary_rounding_ambiguous.any()
    reps.to_csv(DATA / 'panel_a_cluster_representatives.csv', index=False)
    selected = reps.loc[reps.reported_q_clearly_lt_005].copy()
    assert len(selected) == 6
    genes = list(dict.fromkeys(g for hits in selected.input_gene_symbols for g in hits.split('|')))
    assert len(genes) == 21
    membership = pd.DataFrame([[int(g in row.input_gene_symbols.split('|')) for g in genes]
                              for row in selected.itertuples()], index=selected.GO, columns=genes)
    membership.to_csv(DATA / 'panel_b_gene_membership.csv')
    scores.loc[genes].to_csv(DATA / 'panel_b_original_drugz_scores.csv')
    comparison = original[['GO', 'Description', 'LogP', 'Log(q-value)', '#GeneInHitList']].merge(
        final[['GO', 'LogP', 'Log(q-value)', '#GeneInHitList', 'EvidenceCutoff']], on='GO',
        suffixes=('_original_full_48', '_final_prioritised'), validate='one_to_one')
    comparison.to_csv(DATA / 'original_vs_prioritised_audit.csv', index=False)
    evidence.to_csv(DATA / 'saved_evidence_counts.csv', index=False)

    fig = plt.figure(figsize=(10.2, 9.4))
    fig.text(.025, .98, 'a', fontsize=12, fontweight='semibold', va='top')
    fig.text(.07, .976, 'Metascape enrichment of all 48 drugZ depleted hits',
             fontsize=10, fontweight='semibold', va='top')
    ax = fig.add_axes([.39, .505, .20, .425])
    ax.axvline(-np.log10(.05), color='#888888', lw=.8, ls='--', zorder=0)
    values = -reps['Log(q-value)'].to_numpy()
    ax.hlines(np.arange(13), 0, values, color='#D0D0D0', lw=1.5, zorder=1)
    for y, (_, r) in enumerate(reps.iterrows()):
        ax.scatter(-r['Log(q-value)'], y, s=32,
                   facecolor=STRONG if r.reported_q_clearly_lt_005 else 'white',
                   edgecolor=STRONG if r.reported_q_clearly_lt_005 else '#999999', lw=.8, zorder=3)
        ax.text(-.045, y, textwrap.fill(r.Description.strip(), 44), ha='right', va='center',
                transform=ax.get_yaxis_transform(), fontsize=8.4,
                fontweight='semibold' if r.reported_q_clearly_lt_005 else 'normal')
        ax.text(1.13, y, str(int(r['#GeneInGOAndHitList'])), ha='center', va='center',
                transform=ax.get_yaxis_transform(), fontsize=8)
        q = float(r.approx_q_value)
        qtext = f'{q:.2g}' if q < .01 else f'{q:.3f}'.rstrip('0').rstrip('.')
        ax.text(1.47, y, qtext, ha='center', va='center', transform=ax.get_yaxis_transform(), fontsize=8)
        names = ', '.join(r.input_gene_symbols.split('|'))
        ax.text(1.87, y, textwrap.fill(names, 32, break_long_words=False, break_on_hyphens=False),
                ha='left', va='center', transform=ax.get_yaxis_transform(), fontsize=7.7,
                linespacing=1.0, color='#444444')
    ax.text(1.13, 1.025, 'n', transform=ax.transAxes, ha='center', fontsize=8)
    ax.text(1.47, 1.025, 'q ≈', transform=ax.transAxes, ha='center', fontsize=8)
    ax.text(1.87, 1.025, 'Contributing genes', transform=ax.transAxes, ha='left', fontsize=8)
    ax.set(xlim=(-.1, 7.1), ylim=(12.5, -.5), yticks=[], xticks=[0, 2, 4, 6],
           xlabel='−log₁₀(q), original full-list output')
    ax.spines['left'].set_visible(False)

    short = {
        'CORUM:6664': 'STAGA complex',
        'WP2011': 'SREBF / miR33 lipid homeostasis',
        'M145': 'p53 downstream pathway',
        'M11896': 'Lindgren bladder cancer cluster 2A DN',
        'GO:2000036': 'Stem-cell population maintenance',
        'GO:0120545': 'Nucleic-acid conformation isomerase activity',
    }
    fig.text(.025, .419, 'b', fontsize=12, fontweight='semibold', va='top')
    fig.text(.07, .415, 'Contributing genes in the six representative terms with q < 0.05',
             fontsize=9.5, fontweight='semibold', va='top')
    ax = fig.add_axes([.365, .18, .6, .188])
    for y, (go, row) in enumerate(membership.iterrows()):
        ax.axhspan(y - .45, y + .45, color='#F5F5F5', zorder=0)
        for x, gene in enumerate(genes):
            if row[gene]:
                strong = scores.loc[gene, 'fdr_synth'] < .05
                ax.scatter(x, y, s=35, color=STRONG if strong else EXPLORATORY,
                           marker='o' if strong else 'D', edgecolor='white', lw=.3, zorder=3)
    ax.set(yticks=np.arange(len(selected)), yticklabels=[short[g] for g in selected.GO],
           xticks=np.arange(len(genes)), xticklabels=genes, xlim=(-.5, len(genes)-.5),
           ylim=(len(selected)-.5, -.5))
    ax.tick_params(axis='both', length=0, labelsize=8)
    plt.setp(ax.get_xticklabels(), rotation=90, ha='center')
    for tick, gene in zip(ax.get_xticklabels(), genes):
        if scores.loc[gene, 'fdr_synth'] < .05: tick.set_fontweight('semibold')
    for spine in ax.spines.values(): spine.set_visible(False)
    handles = [Line2D([], [], color=STRONG, marker='o', lw=0, ms=5, label='drugZ depletion FDR < 0.05'),
               Line2D([], [], color=EXPLORATORY, marker='D', lw=0, ms=5, label='0.05 ≤ drugZ depletion FDR < 0.25')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.53, .052), ncol=2,
               fontsize=8, handletextpad=.5, columnspacing=2)
    fig.text(.5, .031, 'Custom background: 18,975 genes in the enrichment tables. Dashed line: q = 0.05.',
             ha='center', fontsize=8, color='#555555')
    fig.text(.5, .014, 'Values retain the original export precision; automatic prioritised-subset statistics are not plotted.',
             ha='center', fontsize=8, color='#555555')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600 if ext == 'png' else 300)
    fig.savefig(OUT / 'review' / f'{STEM}_preview.png', dpi=160)
    plt.close(fig)

    with fitz.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all(g in page.get_text() for g in genes)
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside, overlaps = [], []
        for i, a in enumerate(spans):
            if not page.rect.contains(fitz.Rect(a['bbox'])): outside.append(a['text'])
            for b in spans[i+1:]:
                rect = fitz.Rect(a['bbox']) & fitz.Rect(b['bbox'])
                if not rect.is_empty and rect.width > .5 and rect.height > .5:
                    overlaps.append([a['text'], b['text']])
        page.get_pixmap(matrix=fitz.Matrix(1.6, 1.6), colorspace=fitz.csGRAY).save(
            OUT / 'review' / f'{STEM}_greyscale.png')

    record = dict(input_genes=48, exact_drugz_input_match=True, unique_entrez=48,
        report_custom_background_genes=19048, effective_enrichment_background_genes=18975,
        prepared_background_symbols=19112,
        full_background_membership_verifiable_from_export=False,
        all_original_terms=133, displayed_clusters=13, displayed_q_lt_005_clusters=6,
        original_terms_clearly_below_q_005=int(original.reported_q_clearly_lt_005.sum()),
        terms_with_rounded_q_at_boundary=int(original.q_boundary_rounding_ambiguous.sum()),
        contributing_genes=21, panel_a_gene_column=True,
        panel_a_gene_column_entries=int(reps['#GeneInGOAndHitList'].sum()),
        final_subset_term_counts={str(int(k)):int(v) for k,v in final['#GeneInHitList'].value_counts().items()},
        evidence_subset_sizes={str(c):int(evidence.TotalEvidenceCounts.ge(c).sum()) for c in [1,2,3]},
        statistics_source='GO_AllLists.csv only, all terms use 48 input genes; retain rounded LogP and Log(q-value)',
        clustering='Saved final GROUP_ID only; representative reselected by original full-list LogP and source-row tie-break',
        no_new_tests_or_q_values=True, no_gene_filtering=True,
        database_updated='2026-07-01', font='Open Sans', page_inches=[10.2,9.4],
        text_overlaps=overlaps, text_outside_page=outside)
    (DATA / 'validation.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    manifest={key:dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for key,p in SOURCES.items()}
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    assert not overlaps and not outside, (overlaps,outside)
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
