"""Display every Y-annotated Brunello guide using frozen HT1080 screen results."""
from pathlib import Path
import sys, os, json, hashlib, math

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_y_guides'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.backends.backend_pdf import PdfPages
import pymupdf as fitz
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(9)
SCREEN = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
FILES = {
    'original_library': SCREEN / 'brunello/broadgpp-brunello-library-contents.txt',
    'corrected_library': SCREEN / 'brunello/broadgpp-brunello-library-corrected.txt',
    'counts': SCREEN / 'count/ht1080_screen2_R1.count.txt',
    'guide_results': SCREEN / 'test/POST_LOY_vs_CTRL.sgrna_summary.txt',
    'drugz': SCREEN / 'results/drugZ/LOY_vs_CTRL_unfiltered.txt',
}
CTRL, LOY = '#3C5488', '#E64B35'
REVIEW = OUT / 'review'
REVIEW.mkdir(exist_ok=True)
STEM = 'ch04_supp_ht1080_y_guides'
KEY_GENES = ['DDX3Y', 'EIF1AY', 'RPS4Y1', 'RPS4Y2', 'ZFY', 'KDM5D', 'UTY', 'USP9Y']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare():
    original = pd.read_csv(FILES['original_library'], sep='\t')
    annotation = original.loc[original['Genomic Sequence'].astype(str).str.startswith('NC_000024.')].copy()
    assert len(annotation) == annotation['sgRNA Target Sequence'].nunique() == 188
    annotation = annotation.rename(columns={'Target Gene Symbol': 'annotation_gene', 'sgRNA Target Sequence': 'sequence'})
    lib = pd.read_csv(FILES['corrected_library'], sep='\t').rename(columns={'sgRNAID': 'guide_id', 'Seq': 'sequence'})
    selected = lib.loc[lib.sequence.isin(annotation.sequence)].copy()
    selected = selected.merge(annotation, on='sequence', validate='one_to_one')
    assert selected.gene.eq(selected.annotation_gene).all()
    counts = pd.read_csv(FILES['counts'], sep='\t')
    assert counts.sgRNA.is_unique and len(counts) == 77441
    samples = counts.columns[2:].tolist()
    endpoints = samples[2:]
    totals = counts[samples].sum()
    counts = counts.rename(columns={'sgRNA': 'guide_id', 'Gene': 'count_gene'})
    selected = selected.merge(counts, on='guide_id', validate='one_to_one')
    guide = pd.read_csv(FILES['guide_results'], sep='\t').rename(columns={
        'sgrna': 'guide_id', 'Gene': 'mageck_gene',
        **{c: 'mageck_' + c for c in ['control_count', 'treatment_count', 'control_mean', 'treat_mean',
            'LFC', 'control_var', 'adj_var', 'score', 'p.low', 'p.high', 'p.twosided', 'FDR', 'high_in_treatment']}})
    selected = selected.merge(guide, on='guide_id', validate='one_to_one')
    assert selected.gene.eq(selected.count_gene).all() and selected.gene.eq(selected.mageck_gene).all()
    selected = selected.sort_values(['gene', 'guide_id']).reset_index(drop=True)
    selected['guide_number'] = selected.groupby('gene').cumcount() + 1
    assert len(selected) == 188 and selected.groupby('gene').size().eq(4).all()
    selected['all_samples_zero'] = selected[samples].sum(axis=1).eq(0)
    selected['endpoint_zero'] = selected[endpoints].sum(axis=1).eq(0)
    selected['endpoint_samples_detected'] = selected[endpoints].gt(0).sum(axis=1)
    selected['endpoint_detected_any'] = ~selected.endpoint_zero
    selected['endpoint_detected_all4'] = selected.endpoint_samples_detected.eq(4)
    selected['D0_detected_both'] = selected[samples[:2]].gt(0).all(axis=1)
    selected['D0_zero_both'] = selected[samples[:2]].sum(axis=1).eq(0)
    for genotype, cols in [('sgCTRL', endpoints[:2]), ('LOY', endpoints[2:])]:
        selected[f'endpoint_{genotype}_raw_total'] = selected[cols].sum(axis=1)
        selected[f'endpoint_{genotype}_cultures_detected'] = selected[cols].gt(0).sum(axis=1)
        selected[f'endpoint_detected_{genotype}'] = selected[f'endpoint_{genotype}_raw_total'].gt(0)
    for sample in samples:
        selected[sample + '_CPM'] = selected[sample] / totals[sample] * 1e6
    # Descriptive baseline contrast: one sample per genotype, +1 CPM pseudocount.
    # A guide absent from both samples has no observed baseline effect estimate.
    selected['D0_log2FC_CPM_plus_1'] = (
        np.log2(selected.D0_LOY_CPM + 1) - np.log2(selected.D0_CTRL_CPM + 1)
    ).mask(selected.D0_zero_both)
    selected['endpoint_mean_sgCTRL_CPM'] = selected[[s + '_CPM' for s in endpoints[:2]]].mean(axis=1)
    selected['endpoint_mean_LOY_CPM'] = selected[[s + '_CPM' for s in endpoints[2:]]].mean(axis=1)
    drugz_all = pd.read_csv(FILES['drugz'], sep='\t')
    drugz = drugz_all.loc[drugz_all.GENE.isin(selected.gene)].copy().rename(columns={'GENE': 'gene'})
    assert len(drugz) == 47 and drugz.gene.is_unique
    selected = selected.merge(drugz.rename(columns={c: 'drugz_' + c for c in drugz if c != 'gene'}), on='gene', validate='many_to_one')
    genes = KEY_GENES + sorted(set(selected.gene) - set(KEY_GENES))
    gene_summary = drugz.set_index('gene').loc[genes].reset_index()
    gene_summary['guides_total'] = 4
    by_gene = selected.groupby('gene')
    gene_summary['guides_with_endpoint_reads'] = gene_summary.gene.map(by_gene.endpoint_zero.apply(lambda x: int((~x).sum())))
    gene_summary['median_observed_guide_LFC'] = gene_summary.gene.map(selected.loc[~selected.endpoint_zero].groupby('gene').mageck_LFC.median())
    for col in ['endpoint_mean_sgCTRL_CPM', 'endpoint_mean_LOY_CPM']:
        gene_summary[col] = gene_summary.gene.map(by_gene[col].mean())
    front = ['gene', 'guide_number', 'guide_id', 'sequence', 'endpoint_zero', 'all_samples_zero',
             'endpoint_mean_sgCTRL_CPM', 'endpoint_mean_LOY_CPM', 'endpoint_detected_any',
             'endpoint_detected_sgCTRL', 'endpoint_detected_LOY', 'endpoint_detected_all4', 'D0_detected_both',
             'D0_log2FC_CPM_plus_1', 'mageck_LFC', 'mageck_FDR',
             'drugz_normZ', 'drugz_fdr_synth', 'drugz_fdr_supp']
    selected = selected[front + [c for c in selected if c not in front]]
    selected.to_csv(DATA / 'all_y_guides.csv', index=False)
    detection_table(selected).to_csv(DATA / 'guide_detection.csv', index=False)
    gene_summary.to_csv(DATA / 'y_gene_summary.csv', index=False)
    rows = []
    for _, row in selected.iterrows():
        for sample in samples:
            rows.append({'gene': row.gene, 'guide_number': row.guide_number, 'guide_id': row.guide_id, 'sample': sample,
                'timepoint': 'D0' if sample.startswith('D0') else 'Endpoint',
                'genotype': 'sgCTRL' if 'CTRL' in sample else 'LOY',
                'culture': '' if sample.startswith('D0') else sample[-1],
                'raw_count': row[sample], 'CPM': row[sample + '_CPM'], 'CPM_plus_1': row[sample + '_CPM'] + 1,
                'log2_CPM_plus_1': np.log2(row[sample + '_CPM'] + 1)})
    long = pd.DataFrame(rows)
    long.to_csv(DATA / 'all_y_guides_sample_abundance.csv', index=False)
    observed = selected.loc[~selected.endpoint_zero]
    summary = {
        'y_annotated_genes': len(gene_summary), 'y_annotated_guides': len(selected),
        'endpoint_observed_guides': len(observed), 'endpoint_positive_guide_LFC': int(observed.mageck_LFC.gt(0).sum()),
        'endpoint_negative_guide_LFC': int(observed.mageck_LFC.lt(0).sum()),
        'endpoint_detected_all4': int(selected.endpoint_detected_all4.sum()),
        'endpoint_absent_sgCTRL': int((~selected.endpoint_detected_sgCTRL).sum()),
        'endpoint_absent_LOY': int((~selected.endpoint_detected_LOY).sum()),
        'D0_LFC_available_guides': int(selected.D0_log2FC_CPM_plus_1.notna().sum()),
        'D0_LFC_min': float(selected.D0_log2FC_CPM_plus_1.min()),
        'D0_LFC_max': float(selected.D0_log2FC_CPM_plus_1.max()),
        'median_observed_guide_LFC': float(observed.mageck_LFC.median()),
        'minimum_drugz_depletion_fdr': float(gene_summary.fdr_synth.min()),
        'minimum_drugz_enrichment_fdr': float(gene_summary.fdr_supp.min()),
        'depleted_genes_fdr025': int(gene_summary.fdr_synth.lt(.25).sum()),
        'enriched_genes_fdr025': int(gene_summary.fdr_supp.lt(.25).sum()),
        'zero_endpoint_guides': selected.loc[selected.endpoint_zero, ['gene', 'guide_number', 'guide_id', 'all_samples_zero']].to_dict('records'),
        'gene_order': genes, 'sample_library_totals': {k: int(v) for k, v in totals.items()},
        'maximum_endpoint_CPM_plus_1': float(long.loc[long.timepoint.eq('Endpoint'), 'CPM_plus_1'].max()),
    }
    (DATA / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return selected, gene_summary, long, summary


def detection_table(guides):
    columns = ['gene', 'guide_number', 'endpoint_mean_sgCTRL_CPM', 'endpoint_mean_LOY_CPM',
        'endpoint_zero', 'endpoint_detected_any', 'endpoint_detected_sgCTRL', 'endpoint_detected_LOY',
        'endpoint_detected_all4', 'endpoint_samples_detected',
        'endpoint_sgCTRL_cultures_detected', 'endpoint_LOY_cultures_detected',
        'endpoint_sgCTRL_raw_total', 'endpoint_LOY_raw_total',
        'D0_detected_both', 'D0_zero_both', 'D0_CTRL', 'D0_LOY', 'D0_CTRL_CPM', 'D0_LOY_CPM', 'D0_log2FC_CPM_plus_1',
        'all_samples_zero', 'mageck_LFC', 'guide_id', 'sequence']
    return guides[columns]


def overview(guides, genes):
    fig = plt.figure(figsize=(16.8, 11.8))
    order = genes.gene.tolist()
    indexed = guides.set_index(['gene', 'guide_number'])
    values = np.array([[indexed.loc[(g, n), 'mageck_LFC'] for n in range(1, 5)] for g in order])
    missing = np.array([[indexed.loc[(g, n), 'endpoint_zero'] for n in range(1, 5)] for g in order])
    values[missing] = np.nan
    ax = fig.add_axes([.065, .135, .355, .78])
    cmap = LinearSegmentedColormap.from_list('guide_lfc', ['#72518C', '#FAFAFA', '#319581'])
    cmap.set_bad('#D9D9D9')
    cpm_cmap = plt.get_cmap('magma_r')
    lfc_norm, cpm_norm = Normalize(-3, 3), Normalize(0, 80)
    rgba = np.zeros((len(order), 12, 4))
    abundance_values = np.array([[[indexed.loc[(g, n), f'endpoint_mean_{condition}_CPM']
        for condition in ['sgCTRL', 'LOY']] for n in range(1, 5)] for g in order])
    for i in range(len(order)):
        for j in range(4):
            rgba[i, 3*j] = cmap(lfc_norm(values[i, j])) if not missing[i, j] else (.85, .85, .85, 1)
            for k in range(2):
                val = abundance_values[i, j, k]
                rgba[i, 3*j+k+1] = cpm_cmap(cpm_norm(val)) if val > 0 else (.85, .85, .85, 1)
    ax.imshow(rgba, aspect='auto', interpolation='nearest')
    ax.set_xticks(range(12), ['LFC', 'sgCTRL', 'LOY'] * 4, fontsize=8)
    ax.xaxis.tick_top()
    for j, label in enumerate(ax.get_xticklabels()):
        label.set_color(['#222222', CTRL, LOY][j % 3])
    for j in range(4):
        ax.text(3*j+1, 1.036, f'Guide {j+1}', transform=ax.get_xaxis_transform(),
            ha='center', va='bottom', fontsize=9, weight='semibold')
    ax.set_yticks(range(len(order)), order, fontsize=8.6)
    ax.tick_params(length=0)
    for label, gene in zip(ax.get_yticklabels(), order):
        if gene in KEY_GENES: label.set_weight('semibold')
    ax.set_xticks(np.arange(-.5, 12, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(order), 1), minor=True)
    ax.grid(which='minor', color='white', lw=1.2)
    ax.tick_params(which='minor', length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    for x in [2.5, 5.5, 8.5]: ax.axvline(x, color='white', lw=3.2)
    for i in range(len(order)):
        for j in range(4):
            val = values[i, j]
            ax.text(3*j, i, '*' if np.isnan(val) else f'{val:.2f}', ha='center', va='center', fontsize=8,
                color='white' if np.isfinite(val) and abs(val) > 2.2 else '#222222')
            for k in range(2):
                val = abundance_values[i, j, k]
                # Preserve exact zeros; never round a detected guide to displayed zero.
                label = '0*' if val == 0 else f'{val:.1f}' if val >= .1 else f'{val:.2g}'
                ax.text(3*j+k+1, i, label, ha='center', va='center', fontsize=8,
                    color='white' if val > 40 else '#222222')
    cax = fig.add_axes([.065, .08, .17, .013])
    cb = fig.colorbar(ScalarMappable(norm=lfc_norm, cmap=cmap), cax=cax, orientation='horizontal', ticks=[-3, -2, -1, 0, 1, 2, 3])
    cb.set_label('Guide log₂ fold change (LOY / sgCTRL)', labelpad=4)
    cb.outline.set_visible(False)
    cpm_ax = fig.add_axes([.25, .08, .17, .013])
    cb_cpm = fig.colorbar(ScalarMappable(norm=cpm_norm, cmap=cpm_cmap), cax=cpm_ax,
        orientation='horizontal', ticks=[0, 20, 40, 60, 80])
    cb_cpm.set_label('Mean endpoint abundance (CPM)', labelpad=4)
    cb_cpm.outline.set_visible(False)
    fig.text(.065, .974, 'Endpoint MAGeCK LFC and abundance', fontsize=10, weight='semibold')
    fig.text(.015, .98, 'a', fontsize=14, weight='semibold')

    # Baseline counts are individual samples, not replicate means. Use a separate
    # labelled abundance scale to expose the lower D0 guide representation.
    dx = fig.add_axes([.445, .135, .35, .78], sharey=ax)
    d0_values = np.array([[[indexed.loc[(g, n), sample + '_CPM']
        for sample in ['D0_CTRL', 'D0_LOY']] for n in range(1, 5)] for g in order])
    d0_lfc = np.array([[indexed.loc[(g, n), 'D0_log2FC_CPM_plus_1'] for n in range(1, 5)] for g in order])
    d0_norm = Normalize(0, 10)
    d0_rgba = np.zeros((len(order), 12, 4))
    for i in range(len(order)):
        for j in range(4):
            d0_rgba[i, 3*j] = cmap(lfc_norm(d0_lfc[i, j])) if np.isfinite(d0_lfc[i, j]) else (.85, .85, .85, 1)
            for k in range(2):
                val = d0_values[i, j, k]
                d0_rgba[i, 3*j+k+1] = cpm_cmap(d0_norm(val)) if val > 0 else (.85, .85, .85, 1)
    dx.imshow(d0_rgba, aspect='auto', interpolation='nearest')
    dx.set_xticks(range(12), ['LFC', 'sgCTRL', 'LOY'] * 4, fontsize=8)
    dx.xaxis.tick_top()
    for j, label in enumerate(dx.get_xticklabels()): label.set_color(['#222222', CTRL, LOY][j % 3])
    for j in range(4):
        dx.text(3*j+1, 1.036, f'Guide {j+1}', transform=dx.get_xaxis_transform(),
            ha='center', va='bottom', fontsize=9, weight='semibold')
    dx.set_xticks(np.arange(-.5, 12, 1), minor=True)
    dx.grid(which='minor', color='white', lw=1.2)
    dx.tick_params(axis='y', which='both', length=0, labelleft=False)
    dx.tick_params(axis='x', which='both', length=0)
    for spine in dx.spines.values(): spine.set_visible(False)
    for x in [2.5, 5.5, 8.5]: dx.axvline(x, color='white', lw=3.2)
    for i in range(len(order)):
        for j in range(4):
            val = d0_lfc[i, j]
            dx.text(3*j, i, '*' if np.isnan(val) else f'{val:.2f}'.replace('-0.00', '0.00'), ha='center', va='center', fontsize=8,
                color='white' if np.isfinite(val) and abs(val) > 2.2 else '#222222')
            for k in range(2):
                val = d0_values[i, j, k]
                label = '0*' if val == 0 else f'{val:.2f}' if val >= .1 else f'{val:.2g}'
                dx.text(3*j+k+1, i, label, ha='center', va='center', fontsize=8,
                    color='white' if val > 5 else '#222222')
    d0_cax = fig.add_axes([.445, .08, .35, .013])
    d0_cb = fig.colorbar(ScalarMappable(norm=d0_norm, cmap=cpm_cmap), cax=d0_cax,
        orientation='horizontal', ticks=[0, 2, 4, 6, 8, 10])
    d0_cb.set_label('D0 abundance (CPM; separate scale)', labelpad=4)
    d0_cb.outline.set_visible(False)
    fig.text(.445, .974, 'D0 LFC and abundance · one sample per genotype', fontsize=10, weight='semibold')
    fig.text(.43, .98, 'b', fontsize=14, weight='semibold')
    fig.text(.62, .956, 'D0 LFC = log₂[(LOY CPM + 1)/(sgCTRL CPM + 1)]', fontsize=8, ha='center')

    bx = fig.add_axes([.82, .135, .074, .78], sharey=ax)
    yy = np.arange(len(order))
    bx.axvline(0, color='#888888', lw=.7)
    bx.hlines(yy, 0, genes.normZ, color='#CBCBCB', lw=1)
    bx.scatter(genes.normZ, yy, color='#666666', s=17, zorder=3)
    bx.set_xlim(-2, 4)
    bx.set_xticks([-2, 0, 2, 4])
    bx.set_xlabel('drugZ normZ', labelpad=4)
    bx.tick_params(axis='y', which='both', length=0, labelleft=False)
    bx.spines['left'].set_visible(False)
    fig.text(.82, .974, 'Gene-level drugZ', fontsize=10, weight='semibold')
    fig.text(.802, .98, 'c', fontsize=14, weight='semibold')
    for i, row in genes.iterrows():
        bx.text(1.45, i, f'{row.fdr_synth:.3f}', transform=bx.get_yaxis_transform(), fontsize=8, ha='center', va='center')
        bx.text(2.08, i, f'{row.fdr_supp:.3f}', transform=bx.get_yaxis_transform(), fontsize=8, ha='center', va='center')
    bx.text(1.45, 1.015, 'Depletion\nFDR', transform=bx.transAxes, fontsize=8, ha='center', va='bottom')
    bx.text(2.08, 1.015, 'Enrichment\nFDR', transform=bx.transAxes, fontsize=8, ha='center', va='bottom')
    for axis in [ax, dx, bx]: axis.axhline(7.5, color='#999999', lw=.8)
    fig.text(.82, .09, 'Negative: depleted in LOY\nPositive: enriched in LOY', fontsize=8.5, va='top')
    fig.text(.065, .019, '0*: zero reads in the D0 sample or both endpoint cultures of that genotype. LFC *: absent from both genotypes at that timepoint.', fontsize=8)
    fig.text(.065, .005, 'Endpoint LFC: saved MAGeCK; D0 LFC: descriptive CPM + 1 contrast. Abundance colour scales differ. Original gene-level FDRs retained; no D0 adjustment.', fontsize=8)
    return fig


def abundance(guides, subset, page_number):
    fig, axes = plt.subplots(6, 4, figsize=(9.6, 11.8))
    fig.subplots_adjust(left=.075, right=.985, top=.933, bottom=.085, wspace=.19, hspace=.70)
    offsets = {'POST_CTRL_REP1': -.24, 'POST_CTRL_REP2': -.08, 'POST_LOY_REP1': .08, 'POST_LOY_REP2': .24}
    for i, ax in enumerate(axes.flat):
        if i >= len(subset):
            ax.set_visible(False)
            continue
        gene = subset[i]
        rows = guides.loc[guides.gene.eq(gene)].sort_values('guide_number')
        for sample, offset in offsets.items():
            ax.scatter(rows.guide_number + offset, rows[sample + '_CPM'] + 1,
                color=CTRL if 'CTRL' in sample else LOY, marker='o' if sample.endswith('1') else '^',
                s=20, edgecolor='white', lw=.35, zorder=3)
        ax.axhline(1, color='#D0D0D0', lw=.6)
        zero = set(rows.loc[rows.endpoint_zero, 'guide_number'])
        ax.set_xticks([1, 2, 3, 4], [f'{n}*' if n in zero else str(n) for n in range(1, 5)])
        ax.set_xlim(.5, 4.5)
        ax.set_ylim(0, 80)
        ax.set_yticks([0, 20, 40, 60, 80])
        ax.tick_params(labelsize=8, length=2.5, pad=2)
        if i % 4: ax.tick_params(labelleft=False)
        ax.set_title(gene, fontsize=9, pad=4, weight='semibold' if gene in KEY_GENES else 'normal')
    handles = [Line2D([], [], color=c, marker='o', lw=0, ms=5, label=g) for g, c in [('sgCTRL', CTRL), ('LOY', LOY)]]
    handles += [Line2D([], [], color='#666666', marker=m, lw=0, ms=5, label=g) for g, m in [('R1', 'o'), ('R2', '^')]]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.075, .99), ncol=4, columnspacing=1.5)
    fig.text(.98, .973, f'Endpoint abundance · {page_number}/2', ha='right', fontsize=9)
    fig.supylabel('Guide abundance (CPM + 1)', x=.013, fontsize=10)
    fig.supxlabel('Guide number', y=.041, fontsize=10)
    fig.text(.075, .015, '* No endpoint reads. R1/R2: outgrowth cultures from one established population per genotype.', fontsize=8)
    return fig


def export(fig, suffix, report):
    stem = STEM + suffix
    report.savefig(fig)
    fig.savefig(OUT / f'{stem}.pdf', dpi=300)
    fig.savefig(OUT / f'{stem}.svg', dpi=300)
    fig.savefig(OUT / f'{stem}.png', dpi=600)
    fig.savefig(REVIEW / f'{stem}_preview.png', dpi=160)
    plt.close(fig)
    return stem


def workbook(guides, genes, long):
    absent = guides.loc[guides.endpoint_zero]
    absent_labels = '; '.join(f'{r.gene} guide {r.guide_number}' for r in absent.itertuples())
    all_zero_labels = '; '.join(f'{r.gene} guide {r.guide_number}' for r in guides.loc[guides.all_samples_zero].itertuples())
    notes = [
        ('Scope', 'All 188 unique guides annotated to NC_000024.* (chromosome Y), 47 genes, four guides each.'),
        ('Selection', 'Supplied original Brunello genomic annotation; exclusive Y specificity has not been established. PAR/shared targets annotated on X are not included by this rule.'),
        ('Primary results', 'Unfiltered saved drugZ gene-level endpoint results; original full-screen depletion and enrichment FDRs, no new correction or subset testing.'),
        ('Guide effects', 'Saved MAGeCK endpoint LOY versus CTRL log2FC; original non-targeting-guide normalisation. Distinct from drugZ normZ.'),
        ('Counts', 'Raw integrated sgRNA read counts. These are not target-gene expression or Y-chromosome DNA copy number.'),
        ('CPM', 'Raw counts / total counts across all 77,441 guides in the sample x 1,000,000. Individual-culture abundance pages use linear CPM + 1; overview endpoint abundance columns use mean CPM and D0 abundance columns use single-sample CPM, both without a pseudocount.'),
        ('Overview abundance', 'Per-guide sgCTRL and LOY columns show arithmetic mean CPM across the two endpoint cultures in each genotype, without a pseudocount. Zero means no raw reads in either culture.'),
        ('Overview D0', 'Middle panel b shows individual D0 sgCTRL and LOY CPM for all four guides per gene, one sample per genotype; no pseudocount. Its 0-10 colour scale is separate from the endpoint 0-80 scale. Gene-level drugZ is now panel c.'),
        ('D0 LFC', 'D0_log2FC_CPM_plus_1 = log2[(D0 LOY CPM + 1)/(D0 sgCTRL CPM + 1)]. This descriptive contrast uses a +1 CPM pseudocount and one baseline sample per genotype. Blank when both D0 samples have zero reads; no P value or FDR. It is distinct from the saved MAGeCK endpoint estimator and must not be subtracted from that LFC as a matched-normalisation baseline correction.'),
        ('Detection filters', 'Guide detection sheet: endpoint_zero TRUE selects guides with zero reads in all four endpoints; endpoint_detected_sgCTRL/LOY records any reads within each genotype; endpoint_detected_all4 records at least one read in every endpoint. No low-count threshold or guide removal is applied.'),
        ('Baseline detection', 'D0_detected_both requires at least one raw read in each D0 sample; D0_zero_both means no reads in either baseline. Baseline raw counts/CPM are provided for inspection.'),
        ('Design', 'Two endpoint outgrowth cultures from one established population per genotype, one D0 baseline sample per genotype.'),
        ('Baseline', 'D0 is displayed in overview panel b and provided in the table for context; it does not enter the saved endpoint comparisons.'),
        ('Guide numbering', 'Guide IDs sorted alphabetically within each gene; sequence and full ID retained.'),
        ('No endpoint reads', f'{absent_labels}; software LFC retained in table but omitted from heatmap. All six samples zero: {all_zero_labels}.'),
        ('Readme', 'See the source manifest and adjacent source notes for reproduction and exact input hashes.'),
    ]
    with pd.ExcelWriter(DATA / 'ht1080_all_y_guides.xlsx', engine='openpyxl') as writer:
        for name, frame in [('Gene summary', genes), ('Guide detection', detection_table(guides)), ('All 188 guides', guides), ('Sample abundance', long), ('Readme', pd.DataFrame(notes, columns=['Item', 'Description']))]:
            frame.to_excel(writer, sheet_name=name, index=False)
            sheet = writer.sheets[name]
            sheet.freeze_panes = 'C2' if name != 'Readme' else 'A2'
            sheet.auto_filter.ref = sheet.dimensions
            from openpyxl.styles import Font, PatternFill
            from openpyxl.utils import get_column_letter
            for cell in sheet[1]:
                cell.font = Font(name='Open Sans', bold=True, color='FFFFFF')
                cell.fill = PatternFill('solid', fgColor='3C5488')
            for colnum, col in enumerate(frame, 1):
                sheet.column_dimensions[get_column_letter(colnum)].width = min(48, max(14, len(col) + 2))
            if name == 'Readme': sheet.column_dimensions['B'].width = 130


def main():
    guides, genes, long, summary = prepare()
    workbook(guides, genes, long)
    stems = []
    with PdfPages(OUT / f'{STEM}_report.pdf') as report:
        stems.append(export(overview(guides, genes), '_overview', report))
        for part in range(2):
            subset = genes.gene.tolist()[part * 24:(part + 1) * 24]
            stems.append(export(abundance(guides, subset, part + 1), f'_abundance_{part + 1}', report))
    doc = fitz.open(OUT / f'{STEM}_report.pdf')
    assert len(doc) == 3
    all_text = '\n'.join(page.get_text() for page in doc)
    assert all(g in all_text for g in genes.gene)
    fonts = sorted({f[3] for page in doc for f in page.get_fonts()})
    assert all('OpenSans' in f for f in fonts)
    gene_pages = '\n'.join(doc[i].get_text() for i in [1, 2])
    assert all(g in gene_pages for g in genes.gene)
    assert len(long) == 1128 and long.loc[long.timepoint.eq('Endpoint')].shape[0] == 752
    assert summary['maximum_endpoint_CPM_plus_1'] < 80
    assert guides[['D0_CTRL_CPM', 'D0_LOY_CPM']].max().max() < 10
    assert guides.D0_log2FC_CPM_plus_1.abs().max() < 3
    assert guides.D0_log2FC_CPM_plus_1.isna().equals(guides.D0_zero_both)
    # Independent join to the pre-existing full-screen CPM export, if available.
    old = OUT / 'source_data/ht1080_figure1_qc/guide_log2_cpm.csv'
    independent = 'not_checked'
    if old.exists():
        prior = pd.read_csv(old)
        if 'sgRNA' in prior and 'D0_CTRL' in prior:
            joined = guides.merge(prior, left_on='guide_id', right_on='sgRNA', validate='one_to_one', suffixes=('', '_prior'))
            for sample in summary['sample_library_totals']:
                assert np.allclose(np.log2(joined[sample + '_CPM'] + 1), joined[sample + '_prior'])
            independent = 'all_1128_values_agree'
    manifest = {
        'inputs': {key: {'path': str(path), 'sha256': sha(path)} for key, path in FILES.items()},
        'selection': 'Original Brunello Genomic Sequence starts NC_000024.; exact sequence join to corrected library, then unique guide ID joins.',
        'CPM_denominator': 'All 77,441 guide rows in each sample',
        'guide_effect': 'Saved MAGeCK guide LFC; NTC normalised; endpoint LOY versus CTRL; no D0 adjustment',
        'gene_statistics': 'Saved unfiltered drugZ results and original full-screen FDRs, including all guides',
        'plots': stems, 'font': fonts, 'png_dpi': 600,
        'figure_inches': {'overview': [16.8, 11.8], 'abundance_pages': [9.6, 11.8]},
        'overview_abundance': 'Arithmetic mean endpoint CPM within sgCTRL and LOY separately; two cultures each; no pseudocount; zero shown as 0*',
        'overview_D0': 'Panel b: individual sgCTRL and LOY D0 CPM, one sample each, no pseudocount; 0-10 colour scale separate from endpoint 0-80 scale; drugZ moved to panel c',
        'D0_LFC': 'log2[(D0 LOY CPM + 1)/(D0 sgCTRL CPM + 1)]; descriptive one-sample-per-genotype contrast; both-zero guides omitted/NaN; shared -3 to +3 LFC colour scale; distinct from saved MAGeCK endpoint normalisation',
        'detection_filters': 'Count > 0; exact zero and per-genotype/per-culture flags only, no arbitrary low-count cutoff or guide removal',
        'validation': {'guides': 188, 'genes': 47, 'unique_sequences': 188, 'endpoint_marks': 752,
            'overview_finite_LFC': 186, 'overview_no_endpoint_reads': 2, 'report_pages': 3,
            'overview_endpoint_CPM_cells': 376, 'overview_D0_CPM_cells': 376,
            'D0_finite_LFC': int(guides.D0_log2FC_CPM_plus_1.notna().sum()),
            'D0_missing_LFC': int(guides.D0_log2FC_CPM_plus_1.isna().sum()),
            'D0_zero_cells': int(guides[['D0_CTRL', 'D0_LOY']].eq(0).sum().sum()),
            'D0_colour_scale_clipping': False,
            'all_genes_present_in_detail_pages': True, 'endpoint_plot_clipping': False,
            'independent_existing_logCPM_check': independent},
    }
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'summary': summary, 'validation': manifest['validation'], 'report': str(OUT / f'{STEM}_report.pdf')}, indent=2))


if __name__ == '__main__':
    main()
