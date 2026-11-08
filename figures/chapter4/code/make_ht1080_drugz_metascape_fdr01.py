"""Import and plot the user's 34-gene Metascape export (19 October 2026).

Retain original full-list statistics; final-table GROUP_ID is display metadata.
No enrichment tests, drugZ refits or multiple-testing corrections are performed.
"""
from pathlib import Path, PurePosixPath
import csv
import hashlib
import io
import json
import os
import re
import shutil
import stat
import sys
import textwrap
import zipfile

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_drugz_metascape_34_fdr01'
ROOT = THESIS / 'research/chapter4_metascape_drugz_34_fdr01_2026-10-19'
REPORT = ROOT / 'report'
DOWNLOAD = Path(r'C:\Users\micha\Downloads\all.tjmy5fqrd.zip')
STEM = 'ch04_supp_ht1080_drugz_metascape_fdr01_draft'
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

STRONG, EXPLORATORY = '#7A3E85', '#BC8A26'
SOURCES = {
    'zip': ROOT / DOWNLOAD.name,
    'annotation': REPORT / 'metascape_result.xlsx',
    'original_full_list_enrichment': REPORT / 'Enrichment_GO/GO_AllLists.csv',
    'final_table_for_cluster_ids_only': REPORT / 'Enrichment_GO/_FINAL_GO.csv',
    'evidence_counts': REPORT / 'Evidence/Evidence.csv',
    'report': REPORT / 'AnalysisReport.html',
    'drugz': THESIS.parent / 'crispr_screen/ht1080_loy_screen_2/results/drugZ/LOY_vs_CTRL_unfiltered.txt',
    'prepared_input': OUT / 'source_data/ht1080_additional_hit_supplements/drugz_depleted_fdr_0p1_34.txt',
    'prepared_background': OUT / 'source_data/ht1080_additional_hit_supplements/future_metascape_drugz_tested_gene_background.txt',
    'previous_full_48_results': OUT / 'source_data/ht1080_drugz_metascape_48/full_48_original_enrichment.csv',
    'script': Path(__file__),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def import_report():
    """Validate every member before extracting; never execute report content."""
    ROOT.mkdir(parents=True, exist_ok=True)
    archived = SOURCES['zip']
    if not archived.exists():
        shutil.copy2(DOWNLOAD, archived)
    elif DOWNLOAD.exists():
        assert sha(archived) == sha(DOWNLOAD), 'Existing import differs from download'
    entries = []
    seen = set()
    with zipfile.ZipFile(archived) as z:
        assert sum(x.file_size for x in z.infolist()) < 100_000_000
        for item in z.infolist():
            path = PurePosixPath(item.filename)
            assert not path.is_absolute() and '..' not in path.parts
            assert '\\' not in item.filename and ':' not in item.filename
            assert not stat.S_ISLNK(item.external_attr >> 16)
            key = item.filename.casefold()
            assert key not in seen, 'Duplicate member path'
            seen.add(key)
            target = REPORT.joinpath(*path.parts).resolve()
            assert target.is_relative_to(REPORT.resolve())
            entries.append((item, target))
        for item, target in entries:
            if item.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                payload = z.read(item)
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists():
                    assert target.read_bytes() == payload, f'Import changed: {target}'
                else:
                    target.write_bytes(payload)
    manifest = dict(import_date='2026-10-19', original_download=str(DOWNLOAD),
        archived_zip=str(archived), sha256=sha(archived),
        entries=len(entries), uncompressed_bytes=sum(i.file_size for i, _ in entries),
        paths_and_symlinks_validated=True,
        files=[dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=sha(p))
               for i,p in entries if not i.is_dir()])
    (ROOT / 'import_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')


def prepare_tables():
    scores = pd.read_csv(SOURCES['drugz'], sep='\t').set_index('GENE')
    expected = set(scores.index[scores.fdr_synth.lt(.1) & scores.normZ.lt(0)])
    annotation = pd.read_excel(SOURCES['annotation'], sheet_name='Annotation',
                              usecols=['MyList', 'Gene ID', 'Type', 'Tax ID'])
    assert len(annotation) == 34 and annotation.MyList.is_unique and annotation['Gene ID'].is_unique
    assert set(annotation.MyList) == expected == set(SOURCES['prepared_input'].read_text().splitlines())
    annotation.rename(columns={'MyList': 'input_symbol', 'Gene ID': 'entrez_id'}, inplace=True)
    annotation.merge(scores.reset_index(), left_on='input_symbol', right_on='GENE',
                     validate='one_to_one').to_csv(DATA / 'input_34_verified.csv', index=False)
    gene_map = dict(zip(annotation.entrez_id.astype(str), annotation.input_symbol))
    original = pd.read_csv(SOURCES['original_full_list_enrichment'])
    final = pd.read_csv(SOURCES['final_table_for_cluster_ids_only'])
    evidence = pd.read_csv(SOURCES['evidence_counts'])
    assert len(original) == 80 and original.GO.is_unique and final.GO.is_unique
    assert original['#GeneInHitList'].eq(34).all()
    assert original['#TotalGeneInLibrary'].eq(18975).all()
    assert set(original.GeneList) == {'MyList'} and set(original.GO) == set(final.GO)
    original['source_row'] = np.arange(len(original))
    original['input_gene_symbols'] = original.GeneID.apply(lambda x: '|'.join(gene_map[g] for g in x.split('|')))
    assert original.apply(lambda r: len(r.input_gene_symbols.split('|')) == r['#GeneInGOAndHitList'], axis=1).all()
    assert original.apply(lambda r: set(r.input_gene_symbols.split('|')) == set(r.Hits.split('|')), axis=1).all()
    original['approx_q_value'] = 10 ** original['Log(q-value)']
    original['q_boundary_rounding_ambiguous'] = original['Log(q-value)'].eq(-1.3)
    assert not original.q_boundary_rounding_ambiguous.any()
    original['reported_q_clearly_lt_005'] = original['Log(q-value)'].lt(-1.35)
    original = original.merge(final[['GO', 'GROUP_ID']], on='GO', validate='one_to_one')
    original.to_csv(DATA / 'full_34_original_enrichment.csv', index=False)
    reps = original.sort_values(['LogP', 'source_row'], kind='stable').drop_duplicates('GROUP_ID').copy()
    assert len(reps) == 11
    reps.to_csv(DATA / 'panel_a_cluster_representatives.csv', index=False)
    selected = reps.loc[reps.reported_q_clearly_lt_005].copy()
    assert len(selected) == 2 and set(selected.GO) == {'CORUM:6664', 'WP2011'}
    genes = list(dict.fromkeys(g for hits in selected.input_gene_symbols for g in hits.split('|')))
    assert len(genes) == 8
    membership = pd.DataFrame([[int(g in row.input_gene_symbols.split('|')) for g in genes]
                              for row in selected.itertuples()], index=selected.GO, columns=genes)
    membership.to_csv(DATA / 'panel_b_gene_membership.csv')
    scores.loc[genes].to_csv(DATA / 'panel_b_original_drugz_scores.csv')
    audit = original[['GO', 'Description', 'LogP', 'Log(q-value)', '#GeneInHitList']].merge(
        final[['GO', 'LogP', 'Log(q-value)', '#GeneInHitList', 'EvidenceCutoff']], on='GO',
        suffixes=('_original_full_34', '_final_prioritised'), validate='one_to_one')
    audit.to_csv(DATA / 'original_vs_prioritised_audit.csv', index=False)
    evidence.to_csv(DATA / 'saved_evidence_counts.csv', index=False)
    old = pd.read_csv(SOURCES['previous_full_48_results'])
    columns = ['GO', 'Description', 'LogP', 'Log(q-value)', 'approx_q_value',
               '#GeneInHitList', '#GeneInGOAndHitList', 'input_gene_symbols']
    comparison = old[columns].merge(original[columns], on='GO', how='outer',
        suffixes=('_fdr025_48', '_fdr01_34'), indicator='reported_in', validate='one_to_one')
    comparison.to_csv(DATA / 'full_list_fdr025_vs_fdr01_comparison.csv', index=False)
    html = SOURCES['report'].read_text(encoding='utf-8')
    plain = re.sub(r'\s+', ' ', re.sub('<[^>]+>', ' ', html))
    assert '19048 genes' in plain and '2026-07-01' in plain
    return scores, original, final, evidence, reps, selected, genes, membership


def draw(scores, reps, selected, genes, membership):
    # More compact than the 48-gene figure: two supported representative groups.
    fig = plt.figure(figsize=(10.2, 7.4))
    fig.text(.025, .975, 'a', fontsize=12, fontweight='semibold', va='top')
    ax = fig.add_axes([.39, .325, .20, .61])
    ax.axvline(-np.log10(.05), color='#888888', lw=.8, ls='--', zorder=0)
    ax.hlines(np.arange(len(reps)), 0, -reps['Log(q-value)'], color='#D0D0D0', lw=1.5)
    for y, (_, r) in enumerate(reps.iterrows()):
        strong = r.reported_q_clearly_lt_005
        ax.scatter(-r['Log(q-value)'], y, s=32, facecolor=STRONG if strong else 'white',
                   edgecolor=STRONG if strong else '#999999', lw=.8, zorder=3)
        ax.text(-.045, y, textwrap.fill(r.Description.strip(), 43), ha='right', va='center',
                transform=ax.get_yaxis_transform(), fontsize=8.2,
                fontweight='semibold' if strong else 'normal', linespacing=1.05)
        ax.text(1.13, y, str(int(r['#GeneInGOAndHitList'])), ha='center', va='center',
                transform=ax.get_yaxis_transform(), fontsize=8)
        q = float(r.approx_q_value)
        qtext = f'{q:.2g}' if q < .01 else f'{q:.3f}'.rstrip('0').rstrip('.')
        ax.text(1.47, y, qtext, ha='center', va='center', transform=ax.get_yaxis_transform(), fontsize=8)
        names = ', '.join(r.input_gene_symbols.split('|'))
        ax.text(1.87, y, textwrap.fill(names, 32, break_long_words=False, break_on_hyphens=False),
                ha='left', va='center', transform=ax.get_yaxis_transform(), fontsize=7.7,
                linespacing=1.0, color='#444444')
    for x, name, ha in [(1.13, 'n', 'center'), (1.47, 'q ≈', 'center'), (1.87, 'Contributing genes', 'left')]:
        ax.text(x, 1.025, name, transform=ax.transAxes, ha=ha, fontsize=8)
    ax.set(xlim=(-.1, 5.65), ylim=(len(reps)-.5, -.5), yticks=[], xticks=[0, 1, 2, 3, 4, 5],
           xlabel='−log₁₀(Metascape q)')
    ax.spines['left'].set_visible(False)

    fig.text(.025, .24, 'b', fontsize=12, fontweight='semibold', va='top')
    ax = fig.add_axes([.39, .167, .54, .075])
    short = {'CORUM:6664': 'STAGA complex', 'WP2011': 'SREBF / miR33 lipid homeostasis'}
    for y, (_, row) in enumerate(membership.iterrows()):
        ax.axhspan(y-.45, y+.45, color='#F5F5F5', zorder=0)
        for x, gene in enumerate(genes):
            if row[gene]:
                strong = scores.loc[gene, 'fdr_synth'] < .05
                ax.scatter(x, y, s=40, color=STRONG if strong else EXPLORATORY,
                           marker='o' if strong else 'D', edgecolor='white', lw=.3, zorder=3)
    ax.set(yticks=np.arange(len(selected)), yticklabels=[short[g] for g in selected.GO],
           xticks=np.arange(len(genes)), xticklabels=genes, xlim=(-.5,len(genes)-.5),
           ylim=(len(selected)-.5, -.5))
    ax.tick_params(axis='both', length=0, labelsize=8)
    plt.setp(ax.get_xticklabels(), rotation=0, ha='center')
    for tick, gene in zip(ax.get_xticklabels(), genes):
        if scores.loc[gene, 'fdr_synth'] < .05:
            tick.set_fontweight('semibold')
    for spine in ax.spines.values():
        spine.set_visible(False)
    handles = [Line2D([], [], color=STRONG, marker='o', lw=0, ms=5, label='drugZ depletion FDR < 0.05'),
               Line2D([], [], color=EXPLORATORY, marker='D', lw=0, ms=5, label='0.05 ≤ drugZ depletion FDR < 0.1')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.53, .065), ncol=2,
               fontsize=8, handletextpad=.5, columnspacing=2)
    (OUT / 'review').mkdir(exist_ok=True)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600 if ext == 'png' else 300)
    fig.savefig(OUT / 'review' / f'{STEM}_preview.png', dpi=160)
    plt.close(fig)


def validate_pdf(genes, reps):
    with fitz.open(OUT / f'{STEM}.pdf') as doc:
        page = doc[0]
        text = page.get_text()
        assert all(g in text for g in genes)
        assert all(g in text for row in reps.input_gene_symbols for g in row.split('|'))
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside, overlaps = [], []
        for i, a in enumerate(spans):
            if not page.rect.contains(fitz.Rect(a['bbox'])):
                outside.append(a['text'])
            for b in spans[i+1:]:
                rect = fitz.Rect(a['bbox']) & fitz.Rect(b['bbox'])
                if not rect.is_empty and rect.width > .5 and rect.height > .5:
                    overlaps.append([a['text'], b['text']])
        page.get_pixmap(matrix=fitz.Matrix(1.6,1.6), colorspace=fitz.csGRAY).save(
            OUT / 'review' / f'{STEM}_greyscale.png')
    return overlaps, outside


def main():
    import_report()
    scores, original, final, evidence, reps, selected, genes, membership = prepare_tables()
    draw(scores, reps, selected, genes, membership)
    overlaps, outside = validate_pdf(genes, reps)
    record = dict(input_genes=34, drugz_direction='depleted', drugz_fdr_cutoff=.1,
        exact_drugz_input_match=True, unique_entrez=34, report_custom_background_genes=19048,
        effective_enrichment_background_genes=18975, prepared_background_symbols=19112,
        full_background_membership_verifiable_from_export=False, all_original_terms=len(original),
        displayed_clusters=len(reps), displayed_q_lt_005_clusters=len(selected),
        original_terms_clearly_below_q_005=int(original.reported_q_clearly_lt_005.sum()),
        terms_with_rounded_q_at_boundary=int(original.q_boundary_rounding_ambiguous.sum()),
        contributing_genes=len(genes), panel_a_gene_column=True,
        panel_a_gene_column_entries=int(reps['#GeneInGOAndHitList'].sum()),
        final_subset_term_counts={str(int(k)):int(v) for k,v in final['#GeneInHitList'].value_counts().items()},
        evidence_subset_sizes={str(c):int(evidence.TotalEvidenceCounts.ge(c).sum()) for c in [1,2,3]},
        statistics_source='GO_AllLists.csv only, all terms use 34 input genes; retain rounded LogP and Log(q-value)',
        clustering='Saved final GROUP_ID only; representative reselected by original full-list LogP and source-row tie-break',
        no_new_tests_or_q_values=True, no_gene_filtering=True,
        database_updated='2026-07-01', font='Open Sans', page_inches=[10.2,7.4],
        text_overlaps=overlaps, text_outside_page=outside)
    (DATA / 'validation.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    manifest = {key:dict(path=str(p), sha256=sha(p)) for key,p in SOURCES.items()}
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    exports = {ext:dict(path=str(OUT / f'{STEM}.{ext}'), sha256=sha(OUT / f'{STEM}.{ext}'))
               for ext in ['pdf','svg','png']}
    (DATA / 'export_manifest.json').write_text(json.dumps(exports, indent=2), encoding='utf-8')
    assert not overlaps and not outside, (overlaps, outside)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
