"""Publish the separate-batch overlap outputs without changing earlier analyses."""
from pathlib import Path
import csv
import html
import json
import re
import math

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'integration/batch_RNA_overlaps'

def read(name):
    with (OUT / name).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def table(rows, columns):
    return '<table><thead><tr>' + ''.join('<th>' + html.escape(label) + '</th>' for _, label in columns) + '</tr></thead><tbody>' + ''.join(
        '<tr>' + ''.join('<td>' + html.escape(str(row[key])) + '</td>' for key, _ in columns) + '</tr>' for row in rows) + '</tbody></table>'

def figure(stem):
    path = '../integration/batch_RNA_overlaps/' + stem
    return f'<figure><a href="{path}.png"><img src="{path}.png" alt="{html.escape(stem)}"></a><figcaption><a href="{path}.pdf">Vector PDF</a> | <a href="{path}.png">300 dpi PNG</a></figcaption></figure>'

# Independent source-value and fixed-ATAC-coordinate checks across every exported gene.
with (ROOT / 'integration/feature_scatters/feature_gene_effects.csv').open(newline='', encoding='utf-8-sig') as handle:
    previous = {(r['feature'], r['entrez_id']): r for r in csv.DictReader(handle)}
source_rows = {}
for batch in ('AUG_2025', 'MAY_2026'):
    source_rows[batch] = {r['gene_id']: r for r in read(f'sources/ht1080_{batch}_LOY_v_CTRL_DESeq2.csv')}
for row in read('all_gene_results.csv'):
    original = source_rows[row['batch']][row['rna_ensembl_id']]
    assert row['entrez_id'] == original['entrezid']
    for target, source in [('rna_logFC', 'logFC'), ('rna_FDR', 'FDR')]:
        assert math.isclose(float(row[target]), float(original[source]), rel_tol=1e-12, abs_tol=1e-300)
    prior = previous.get((row['feature'], row['entrez_id']))
    if prior:
        for field in ('atac_median_logFC', 'minimum_peak_FDR'):
            assert math.isclose(float(row[field]), float(prior[field]), rel_tol=1e-12, abs_tol=1e-14)
with (OUT / 'verification.txt').open('a', encoding='utf-8') as handle:
    handle.write('\nPASS: Python independently verifies all exported RNA effects/FDR against source and ATAC coordinates against prior integration.\n')
stats = read('overlap_summary.csv')
native = [r for r in stats if r['universe_type'] == 'batch-specific']
common = [r for r in stats if r['universe_type'] == 'common-to-both']
for row in stats:
    row['spearman_all_genes'] = f"{float(row['spearman_all_genes']):.3f}"
columns = [('batch', 'RNA batch'), ('feature', 'ATAC feature'), ('eligible_genes', 'Eligible genes'),
           ('ATAC_genes', 'ATAC significant'), ('RNA_genes', 'RNA significant'), ('Both', 'Overlap'),
           ('Up_up', 'Up/up'), ('Down_down', 'Down/down'), ('Opposing', 'Opposing'), ('spearman_all_genes', 'All-gene Spearman rho')]
rep = read('replication_summary.csv')
support = read('combined_down70_support_summary.csv')
prom = {r['batch']: r for r in native if r['feature'] == 'promoter'}
pr = next(r for r in rep if r['feature'] == 'promoter')
intro = (f"Promoter overlap: August 2025 has {prom['AUG_2025']['Both']} genes "
         f"({prom['AUG_2025']['Down_down']} down/down and {prom['AUG_2025']['Up_up']} up/up); "
         f"May 2026 has {prom['MAY_2026']['Both']} genes "
         f"({prom['MAY_2026']['Down_down']} down/down and {prom['MAY_2026']['Up_up']} up/up). "
         f"On the common eligible background, {pr['overlap_in_both_batches']} genes overlap ATAC in both RNA batches, "
         f"including {pr['down_down_in_both_batches']} down/down in both.")
parts = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>ATAC overlap with separate RNA-seq batches</title>',
         '<style>body{font-family:Arial,sans-serif;max-width:1400px;margin:35px auto;padding:0 25px;color:#24313b;line-height:1.5}h1,h2{color:#182c3d}table{border-collapse:collapse;width:100%;font-size:14px;margin:20px 0}th,td{border:1px solid #d8dfe4;padding:8px;text-align:left}th{background:#edf3f5}img{max-width:100%;height:auto}figure{margin:25px 0}a{color:#006da0}.note{background:#eef4f8;padding:16px}code{overflow-wrap:anywhere}</style><body>',
         '<h1>ATAC overlap with August 2025 and May 2026 RNA-seq separately</h1>',
         '<p><strong>ATAC FDR &lt;0.10; RNA FDR &lt;0.05. LOY versus CTRL; no fold-change filter.</strong></p>',
         '<p>' + intro + '</p>',
         '<p>Both RNA results come from the batch-stratified section of the HT1080_BATCHES notebook. Each DESeq2 model uses <code>~ treatment</code> and only three LOY plus three CTRL samples from that batch. Parental and add-back samples do not enter these separate fits. Existing effect estimates and adjusted P-values are reused unchanged.</p>',
         '<p class="note">These are separately fitted results within the combined notebook, using its initial gene filter (at least 10 counts in at least 3 combined-input samples). They are not asserted to reproduce every historical standalone RNA analysis. Each batch has its own DESeq2 fit and eligible gene background. The main ATAC model and peak-to-gene links are unchanged.</p>',
         '<h2>Overlap and concordant directions</h2>', table(native, columns),
         '<p>ATAC membership means any linked peak passes FDR. Up/up and down/down use the median effect across all tested linked ATAC peaks and the RNA effect. This median has no gene-level FDR. Counts of significant gained and lost linked peaks are included in the gene tables.</p>',
         figure('feature_overlap_grid'), '<h2>Promoter accessibility versus RNA</h2>', figure('promoter_scatters'),
         '<p>Up to five jointly significant genes per concordant direction per batch are labelled, selected by the smaller absolute effect across the two assays. Grey points are genes significant in neither assay. The two panels share axes and retain every eligible gene.</p>',
         '<h2>All genomic features</h2>', figure('feature_scatter_grid'),
         '<p>Promoter links use TSS +/-2 kb, gene-body links use overlaps, and intergenic links use the nearest TSS as an exploratory assignment. A gene may occur in multiple feature classes. Intergenic assignments do not establish enhancer targets.</p>',
         '<h2>Comparison using genes eligible in both batches</h2>', table(common, columns),
         table(rep, [('feature', 'ATAC feature'), ('common_eligible', 'Common eligible'), ('overlap_in_both_batches', 'Overlap in both batches'), ('down_down_in_both_batches', 'Down/down in both'), ('up_up_in_both_batches', 'Up/up in both')]),
         '<h2>What happens to the combined-model 70 down/down genes?</h2>',
         table(support, [('batch', 'RNA batch'), ('original_genes', 'Original genes'), ('eligible', 'Eligible in batch'), ('RNA_down_significant', 'Still significant down/down'), ('RNA_down_any_FDR', 'Negative RNA effect, any FDR'), ('ZNF_ZFP_down_significant', 'Significant down/down among 18 ZNF/ZFP genes')]),
         '<p>The combined-model 70 genes were selected using these same RNA batches, so this is a descriptive consistency check, not independent validation. Lower significance in one batch can reflect power as well as effect differences; it does not by itself demonstrate a batch-by-treatment interaction.</p>',
         '<h2>Gene lists and reproducibility</h2><ul>']
for batch in ('AUG_2025', 'MAY_2026'):
    for feature in ('promoter', 'gene_body', 'intergenic'):
        name = f'{batch}_{feature}_overlap_genes.csv'
        parts.append(f'<li><a href="../integration/batch_RNA_overlaps/{name}">{batch}: {feature} overlap genes</a></li>')
for name, label in [('all_gene_results.csv', 'All eligible genes and their effects'), ('common_gene_batch_comparison.csv', 'Both batches side by side'),
                    ('combined_down70_batch_support.csv', 'Batch support for each original down/down gene'), ('source_manifest.csv', 'Source files and SHA256 checksums'),
                    ('gene_mapping_audit.csv', 'All gene mapping and eligibility decisions'), ('METHODS.txt', 'Complete methods'), ('verification.txt', 'Verification')]:
    parts.append(f'<li><a href="../integration/batch_RNA_overlaps/{name}">{label}</a></li>')
parts += ['</ul><p>Effects are compared at group level; direct ATAC/RNA sample matching is not established. No passage term is included in these RNA fits, and no copy-number adjustment has been applied. Shared peaks, clustered genes and variable peak counts limit gene independence; no overlap-enrichment P-value is assigned.</p>',
          '<p><a href="ATACseq_report.html">Main ATAC report</a></p></body></html>']
(ROOT / 'report/batch_RNA_overlaps.html').write_text('\n'.join(parts), encoding='utf-8')
snippet = '<section id="batch-rna-overlaps"><h2>ATAC overlap with each RNA-seq batch separately</h2><p>ATAC FDR &lt;0.10 and RNA FDR &lt;0.05. ' + intro + ' <a href="batch_RNA_overlaps.html">Separate-batch overlap plots, scatter plots and gene lists</a>.</p></section>'
main = ROOT / 'report/ATACseq_report.html'
content = main.read_text(encoding='utf-8')
if '<section id="batch-rna-overlaps">' in content:
    content = re.sub(r'<section id="batch-rna-overlaps">.*?</section>', lambda _: snippet, content, flags=re.S)
else:
    content = content.replace('</body>', snippet + '\n</body>')
main.write_text(content, encoding='utf-8')
# Verify every local link and image in the new report resolves.
for link in re.findall(r'(?:href|src)="([^"]+)"', '\n'.join(parts)):
    assert (ROOT / 'report' / link).resolve().exists(), link
(OUT / 'ANALYSIS_COMPLETE.json').write_text(json.dumps({'status': 'complete', 'report': 'report/batch_RNA_overlaps.html',
    'ATAC_FDR': 0.1, 'RNA_FDR': 0.05, 'RNA_model': '~ treatment separately in each batch', 'source_fits_reused': True,
    'validation': 'source hashes, gene mapping, set partitions, directions, exact70 accounting, report links passed'}, indent=2), encoding='utf-8')
print(intro)
print('Saved report/batch_RNA_overlaps.html')
