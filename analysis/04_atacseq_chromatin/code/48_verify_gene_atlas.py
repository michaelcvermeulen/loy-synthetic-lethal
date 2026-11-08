"""Independently check gene selection, peak measurements and final atlas files."""
from pathlib import Path
import json
import math
import re
import zipfile
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'integration/gene_atlas'
selection = json.loads((OUT / 'selection_summary.json').read_text())
genes = pd.read_csv(OUT / 'data/selected_genes.csv', dtype={'entrez_id': str})
reasons = pd.read_csv(OUT / 'data/selection_reasons.csv', dtype={'entrez_id': str})
peaks = pd.read_csv(ROOT / 'annotation/annotated_differential_peaks.csv').set_index('peak_id')
links = pd.read_csv(ROOT / 'integration/feature_overlaps/peak_gene_link_audit.csv', dtype={'entrez_id': str})
strong = set(peaks[(peaks.FDR < selection['strong_FDR']) & (peaks.logFC.abs() > selection['strong_abs_logFC'])].index)
assert set(reasons.loc[reasons.reason == 'Strong DAR', 'qualifying_peak_id']) == strong
expected = set(links.loc[links.peak_id.isin(strong), 'entrez_id'])
combined = pd.read_csv(ROOT / 'integration/feature_scatters/feature_gene_effects.csv', dtype={'entrez_id': str}).assign(RNA_result='Combined')
batch = pd.read_csv(ROOT / 'integration/batch_RNA_overlaps/all_gene_results.csv', dtype={'entrez_id': str}).rename(columns={'batch': 'RNA_result'})
effects = pd.concat([combined, batch], ignore_index=True)
scope = selection['RNA_scope']
source_set = ['Combined', 'AUG_2025', 'MAY_2026'] if scope == 'all' else (['Combined'] if scope == 'combined' else ['AUG_2025', 'MAY_2026'])
joint = effects[(effects.RNA_result.isin(source_set)) & (effects.minimum_peak_FDR < .1) & (effects.rna_FDR < .05)]
expected.update(joint.entrez_id)
assert genes.entrez_id.is_unique and set(genes.entrez_id) == expected
for row in joint.itertuples():
    linked = links[(links.feature == row.feature) & (links.entrez_id == row.entrez_id) & (links.atac_FDR < .1)]
    actual = reasons[(reasons.entrez_id == row.entrez_id) & (reasons.feature == row.feature) & (reasons.RNA_result == row.RNA_result)]
    assert set(actual.qualifying_peak_id) == set(linked.peak_id)
raw = pd.read_csv(ROOT / 'counts/peak_counts.tsv', sep='\t').set_index('peak_id')
norm = pd.read_csv(ROOT / 'differential/normalization_factors.csv')
norm = norm[norm.arm == 'primary_TMM'].set_index('sample')
windows = pd.read_csv(OUT / 'data/plot_windows.csv', dtype={'entrez_id': str})
measurements = 0
for gene in genes.itertuples():
    saved = pd.read_csv(OUT / f'data/{gene.entrez_id}_peak_measurements.csv').set_index('peak_id')
    wanted = set(links.loc[links.entrez_id == gene.entrez_id, 'peak_id'])
    assert saved.index.is_unique and set(saved.index) == wanted
    assert np.allclose(saved.logFC, peaks.loc[saved.index, 'logFC'], rtol=1e-12, atol=1e-14)
    assert np.allclose(saved.FDR, peaks.loc[saved.index, 'FDR'], rtol=1e-12, atol=1e-14)
    for sample in raw.columns:
        assert np.array_equal(saved[sample + '_raw_count'], raw.loc[saved.index, sample])
        calculated = raw.loc[saved.index, sample] * 1e6 / norm.loc[sample, 'effective_size']
        assert np.allclose(saved[sample + '_TMM_CPM'], calculated, rtol=1e-12, atol=1e-14)
    measurements += len(saved)
    overview = windows[(windows.entrez_id == gene.entrez_id) & (windows.view == 'Gene overview')]
    assert len(overview) == 1 and overview.iloc[0].start <= gene.start - 1 and overview.iloc[0].end >= gene.end
    closeups = windows[(windows.entrez_id == gene.entrez_id) & (windows.view != 'Gene overview')]
    for pid in reasons.loc[reasons.entrez_id == gene.entrez_id, 'qualifying_peak_id'].unique():
        peak = peaks.loc[pid]
        assert ((closeups.chrom == peak.chrom) & (closeups.start <= peak.start) & (closeups.end >= peak.end)).any()
render = pd.read_csv(OUT / 'data/render_audit.csv')
assert len(render) == len(windows) and set(render.window_id) == set(windows.window_id)
for row in render.itertuples():
    data = np.load(OUT / f'coverage/{row.window_id}.npz')
    assert data['values'].max() < row.shared_coverage_ymax
    assert len(data['samples']) == 6
final = json.loads((OUT / 'ANALYSIS_COMPLETE.json').read_text())
assert final['genes'] == len(genes) and final['views'] == len(windows)
result = {'status': 'passed', 'genes': len(genes), 'views': len(windows), 'distinct_qualifying_peaks': reasons.qualifying_peak_id.nunique(),
          'gene_peak_measurement_rows_checked': measurements,
          'checks': ['Exact strong-DAR and RNA-overlap union reconstructed from original inputs',
                     'All overlap inclusion reasons reproduce the correct feature and RNA result',
                     'Every linked peak raw count, TMM CPM, log2FC and FDR matches its original source',
                     'Every selected gene is fully included in its overview; every qualifying peak has a close-up',
                     'All six coverage tracks share the displayed limit and no plotted coverage value is clipped']}
(OUT / 'verification.json').write_text(json.dumps(result, indent=2))
archive = OUT / 'gene_diagrams_PDFs.zip'
with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as bundle:
    for gene in genes.sort_values(['symbol', 'entrez_id']).itertuples():
        symbol = re.sub(r'[^A-Za-z0-9_.-]', '_', str(gene.symbol))
        bundle.write(OUT / f'figures/{gene.entrez_id}.pdf', f'{symbol}_{gene.entrez_id}.pdf')
    for name in ['METHODS.txt', 'selection_summary.json', 'verification.json', 'data/selected_genes.csv', 'data/selection_reasons.csv']:
        bundle.write(OUT / name, name)
with zipfile.ZipFile(archive) as bundle:
    assert bundle.testzip() is None
    assert sum(name.endswith('.pdf') for name in bundle.namelist()) == len(genes)
report = ROOT / 'report/gene_atlas.html'
text = report.read_text(encoding='utf-8')
snippet = '<p id="bulk-pdf-download"><a href="../integration/gene_atlas/gene_diagrams_PDFs.zip">Download all 251 gene PDF books as one ZIP</a> | <a href="../integration/gene_atlas/verification.json">Data verification</a></p>'
snippet = snippet.replace('251', str(len(genes)))
if 'id="bulk-pdf-download"' not in text:
    text = text.replace('<h2>Data and methods</h2>', snippet + '<h2>Data and methods</h2>')
report.write_text(text, encoding='utf-8')
final['independent_data_validation'] = 'passed'
final['PDF_archive'] = 'gene_diagrams_PDFs.zip'
(OUT / 'ANALYSIS_COMPLETE.json').write_text(json.dumps(final, indent=2))
print(json.dumps(result, indent=2))
