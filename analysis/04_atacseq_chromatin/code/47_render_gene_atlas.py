"""Render reproducible gene diagrams, per-gene PDFs, and a searchable HTML atlas."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import csv
import hashlib
import html
import json
import math
import re
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'integration/gene_atlas'
CTRL = '#0072B2'
LOY = '#D55E00'
PURPLE = '#7B3294'
SAMPLES = ['CTRL_11', 'CTRL_12', 'CTRL_21', 'LOY_11', 'LOY_12', 'LOY_21']
RNA_NAMES = {'Combined': 'Combined', 'AUG_2025': 'Aug 2025', 'MAY_2026': 'May 2026'}
DATA = {}

def esc(value):
    return html.escape(str(value))

def fmt(value, digits=3):
    if pd.isna(value):
        return 'NA'
    return f'{float(value):.{digits}g}'

def init_worker():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False})
    DATA['genes'] = pd.read_csv(OUT / 'data/selected_genes.csv', dtype={'entrez_id': str}).set_index('entrez_id')
    DATA['models'] = pd.read_csv(OUT / 'data/gene_models.csv', dtype={'gene_id': str})
    DATA['exons'] = {k: v for k, v in pd.read_csv(OUT / 'data/merged_gene_exons.csv', dtype={'gene_id': str}).groupby('gene_id')}
    DATA['windows'] = pd.read_csv(OUT / 'data/plot_windows.csv', dtype={'entrez_id': str})
    DATA['reasons'] = pd.read_csv(OUT / 'data/selection_reasons.csv', dtype={'entrez_id': str}).fillna('')
    DATA['links'] = pd.read_csv(OUT / 'data/all_selected_gene_peak_links.csv', dtype={'entrez_id': str})
    DATA['peaks'] = pd.read_csv(ROOT / 'annotation/annotated_differential_peaks.csv').set_index('peak_id')
    DATA['rna'] = pd.read_csv(OUT / 'data/RNA_source_rows.csv', dtype={'entrezid': str})
    DATA['effects'] = pd.read_csv(OUT / 'data/selected_gene_RNA_ATAC_effects.csv', dtype={'entrez_id': str})
    DATA['counts'] = pd.read_csv(OUT / 'data/peak_counts.csv').set_index('peak_id')
    norm = pd.read_csv(ROOT / 'differential/normalization_factors.csv')
    eff = norm[norm.arm == 'primary_TMM'].set_index('sample').loc[SAMPLES, 'effective_size']
    DATA['cpm'] = DATA['counts'][SAMPLES].div(eff, axis=1) * 1e6
    DATA['code_hash'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

def model_layout(chrom, left, right, target):
    models = DATA['models']
    models = models[(models.chrom == chrom) & (models.start - 1 < right) & (models.end > left)].copy()
    models['plot_start'] = np.maximum(left, models.start - 1)
    models = models.sort_values(['plot_start', 'end', 'gene_id'])
    lane_ends = []
    layout = []
    span = right - left
    for row in models.itertuples():
        lo = max(left, row.start - 1)
        hi = min(right, row.end)
        label_hi = lo + max(6, len(str(row.symbol))) * span * .006
        lane = next((i for i, end in enumerate(lane_ends) if end < lo - .012 * span), len(lane_ends))
        end = max(hi, label_hi)
        if lane == len(lane_ends):
            lane_ends.append(end)
        else:
            lane_ends[lane] = end
        layout.append((row, lane))
    return layout, max(1, len(lane_ends))

def draw_gene_models(ax, layout, lanes, left, right, target):
    span = right - left
    for gene, lane in layout:
        y = lanes - lane - .65
        lo, hi = max(left, gene.start - 1), min(right, gene.end)
        color = PURPLE if gene.gene_id == target else '#758592'
        ax.plot([lo, hi], [y, y], color=color, lw=.7)
        exons = DATA['exons'].get(gene.gene_id)
        if exons is not None:
            for exon in exons[(exons.start0 < right) & (exons.end > left)].itertuples():
                a, b = max(left, exon.start0), min(right, exon.end)
                ax.add_patch(Rectangle((a, y - .12), b - a, .24, color=color, linewidth=0))
        points = np.arange(lo + .016 * span, hi, .04 * span)
        if len(points):
            ax.scatter(points, np.full(len(points), y), marker='>' if gene.strand == '+' else '<', s=8, color=color, linewidths=0)
        label = str(gene.symbol) + (' *' if gene.gene_id == target else '')
        near_right = lo > right - max(8, len(label)) * span * .006
        ax.text(right - .003 * span if near_right else lo + .002 * span, y + .19, label, fontsize=6.4, color=color,
                ha='right' if near_right else 'left',
                weight='bold' if gene.gene_id == target else 'normal', clip_on=True)
    ax.set_ylim(0, lanes + .3)
    ax.set_yticks([])
    ax.set_ylabel('Genes', rotation=0, ha='right', va='center', labelpad=9)
    ax.set_xlabel('hg38 / GRCh38 genomic position (bp; 0-based axis)')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x):,}'))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.spines['left'].set_visible(False)

def rna_rows(gene_id):
    r = DATA['rna']
    output = []
    for name in RNA_NAMES:
        rows = r[(r.entrezid == gene_id) & (r.RNA_result == name)]
        if len(rows) == 1 and rows.iloc[0].entrez_multiplicity == 1:
            row = rows.iloc[0]
            output.append([RNA_NAMES[name], fmt(row.logFC), fmt(row.FDR)])
        else:
            output.append([RNA_NAMES[name], 'NA', 'ambiguous' if len(rows) else 'not tested'])
    return output

def draw_window(gene_id, gene, window):
    left, right = int(window.start), int(window.end)
    coverage = np.load(OUT / f'coverage/{window.window_id}.npz')
    assert list(coverage['samples']) == SAMPLES
    vals = coverage['values']
    starts, ends = coverage['starts'], coverage['ends']
    assert starts[0] == left and ends[-1] == right and vals.shape == (6, len(starts))
    xedges = np.r_[starts, ends[-1]]
    ymax = max(.05, float(vals.max()) * 1.06)
    region_peaks = DATA['peaks'][(DATA['peaks'].chrom == window.chrom) & (DATA['peaks'].start < right) & (DATA['peaks'].end > left)]
    gene_links = DATA['links'][DATA['links'].entrez_id == gene_id]
    linked_ids = set(gene_links.peak_id)
    selection = DATA['reasons'][DATA['reasons'].entrez_id == gene_id]
    qualifying_ids = set(selection.qualifying_peak_id)
    qualifying = region_peaks[region_peaks.index.isin(qualifying_ids)].sort_values(['FDR', 'start'])
    layout, lanes = model_layout(window.chrom, left, right, gene_id)
    model_height = max(1.0, .23 * lanes)
    heights = [1.25] + [.58] * 6 + [.9, model_height]
    fig = plt.figure(figsize=(13.2, 4.7 + sum(heights)))
    grid = fig.add_gridspec(9, 1, left=.085, right=.735, top=.855, bottom=.135, height_ratios=heights, hspace=.24)
    axes = []
    for i in range(9):
        axes.append(fig.add_subplot(grid[i], sharex=axes[0] if axes else None))
    for ax in axes:
        ax.set_xlim(left, right)
    for ax in axes[:-1]:
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
    mean_ax = axes[0]
    mean_ax.stairs(vals[:3].mean(axis=0), xedges, color=CTRL, lw=1.0, label='CTRL mean')
    mean_ax.stairs(vals[3:].mean(axis=0), xedges, color=LOY, lw=1.0, label='LOY mean')
    mean_ax.set_ylabel('Group\nmeans', rotation=0, ha='right', va='center', labelpad=10)
    mean_ax.legend(loc='upper right', ncol=2, fontsize=7, frameon=False)
    mean_ax.set_ylim(0, ymax)
    for i, sample in enumerate(SAMPLES):
        ax = axes[i + 1]
        color = CTRL if i < 3 else LOY
        ax.stairs(vals[i], xedges, color=color, fill=True, linewidth=.5, alpha=.8)
        ax.set_ylim(0, ymax)
        ax.set_yticks([0, ymax / 2])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.2g}'))
        ax.tick_params(axis='y', labelsize=6)
        ax.set_ylabel(sample, rotation=0, ha='right', va='center', color=color, fontsize=8, labelpad=8)
    for ax in axes[:7]:
        for peak in qualifying.itertuples():
            ax.axvspan(peak.start, peak.end, color=PURPLE, alpha=.065, linewidth=0)
    # Local DAR effects include all tested peaks; gene-linked peaks receive black outlines.
    effect_ax = axes[7]
    lim = max(.5, float(region_peaks.logFC.abs().max()) * 1.15) if len(region_peaks) else 1
    for pid, peak in region_peaks.iterrows():
        color = (LOY if peak.logFC > 0 else CTRL) if peak.FDR < .1 else '#BEC6CC'
        effect_ax.add_patch(Rectangle((peak.start, min(0, peak.logFC)), peak.end - peak.start,
                                     abs(peak.logFC), facecolor=color,
                                     edgecolor='#222222' if pid in linked_ids else color,
                                     linewidth=.65 if pid in linked_ids else .25))
    effect_ax.axhline(0, color='#555555', lw=.5)
    effect_ax.set_ylim(-lim, lim)
    effect_ax.set_ylabel('Peak\nlog2FC', rotation=0, ha='right', va='center', labelpad=10)
    effect_ax.yaxis.set_major_locator(MaxNLocator(3))
    for number, (pid, peak) in enumerate(qualifying.head(3).iterrows(), 1):
        mean_ax.text((peak.start + peak.end) / 2, ymax * .78, f'P{number}', ha='center', fontsize=7,
                     bbox={'facecolor': 'white', 'alpha': .8, 'edgecolor': 'none', 'pad': 1})
    draw_gene_models(axes[8], layout, lanes, left, right, gene_id)
    symbol = str(gene.symbol)
    fig.text(.085, .966, f'{symbol}  |  {window.view}', weight='bold', fontsize=17, va='top')
    fig.text(.085, .935, f'Entrez {gene_id}  |  {window.chrom}:{left:,}-{right:,}  |  gene strand {gene.strand}', fontsize=10)
    reason_text = []
    if gene.strong_DAR:
        reason_text.append('Strong DAR')
    for name in RNA_NAMES:
        if bool(gene['overlap_' + name]):
            reason_text.append('Overlap: ' + RNA_NAMES[name])
    fig.text(.085, .911, '  |  '.join(reason_text), fontsize=9, color=PURPLE)
    fig.text(.085, .884, f'ATAC fragment coverage per million non-Y fragments; shared range 0-{ymax:.3g}; plot bins {int(coverage["bin_width"])} bp', fontsize=8)
    # Right-hand panel: RNA estimates and original peak-model measurements.
    info = fig.add_axes([.78, .68, .20, .16]); info.axis('off')
    info.text(0, 1.05, 'RNA: LOY / CTRL', weight='bold', fontsize=10)
    t = info.table(cellText=rna_rows(gene_id), colLabels=['Result', 'log2FC', 'FDR'], loc='upper left',
                   cellLoc='left', colLoc='left', colWidths=[.43, .25, .32], bbox=[0, .43, 1, .56])
    t.auto_set_font_size(False); t.set_fontsize(8)
    for cell in t.get_celld().values():
        cell.set_edgecolor('#DCE2E6'); cell.set_linewidth(.5)
    info.text(0, .27, 'Gene links: ' + str(gene.features).replace(';', ', '), fontsize=8, wrap=True)
    info.text(0, .09, f'{len(qualifying)} qualifying peak(s) in this view', fontsize=8)
    for number, (pid, peak) in enumerate(qualifying.head(3).iterrows(), 1):
        y = .535 - (number - 1) * .145
        ax = fig.add_axes([.805, y, .165, .095])
        values = DATA['cpm'].loc[pid, SAMPLES].to_numpy(float)
        jitter = np.array([-.12, 0, .12])
        for group, offset, color in [(values[:3], 0, CTRL), (values[3:], 1, LOY)]:
            ax.scatter(offset + jitter, group, s=20, color=color, edgecolor='white', linewidth=.3, zorder=3)
            ax.plot([offset - .23, offset + .23], [group.mean()] * 2, color=color, lw=1.1)
        ax.set_xlim(-.5, 1.5); ax.set_ylim(bottom=0)
        ax.set_xticks([0, 1], ['CTRL', 'LOY'], fontsize=7)
        ax.tick_params(axis='y', labelsize=7)
        ax.set_ylabel('TMM peak CPM', fontsize=7)
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.set_title(f'P{number}: log2FC {peak.logFC:+.2f} | FDR {peak.FDR:.3g}\n{pid}', fontsize=7.5, loc='left')
    fig.text(.78, .14, 'Dots: individual replicates\nBars: group means\n\nAll qualifying peaks and raw counts\nare listed on the gene page.', fontsize=8, va='bottom')
    footer = ('Peak track: orange = gained, blue = lost at FDR <0.10; grey = other tested peaks. Black outlines = links to this gene.\n'
              'Purple shading = selection-supporting peaks. Gene models merge exons across isoforms; * marks the selected gene.\n'
              'Coverage uses non-Y fragment CPM; peak measurements use primary-model TMM normalization. No copy-number adjustment.\n'
              'Intergenic nearest-TSS links are exploratory. Primary ATAC model: ~ condition; positive log2FC means higher in LOY.')
    fig.text(.085, .025, footer, fontsize=7.5, va='bottom', linespacing=1.4)
    audit = {'window_id': window.window_id, 'entrez_id': gene_id, 'chrom': window.chrom, 'start': left, 'end': right,
             'display_bin_bp': int(coverage['bin_width']), 'shared_coverage_ymax': ymax, 'samples': 6,
             'local_tested_peaks': len(region_peaks), 'local_qualifying_peaks': len(qualifying), 'gene_models_shown': len(layout)}
    return fig, audit

STYLE = '''body{font:15px/1.5 Arial,sans-serif;color:#24313b;max-width:1450px;margin:30px auto;padding:0 25px}h1,h2{color:#172f41}a{color:#006ea1}img{width:100%;height:auto}table{border-collapse:collapse;width:100%;font-size:13px;margin:15px 0}th,td{border-bottom:1px solid #d9e1e6;padding:7px;text-align:left}th{background:#edf3f5;position:sticky;top:0}input,select{font:inherit;padding:8px;margin:5px}figure{margin:20px 0}.note{background:#edf4f8;padding:14px}.scroll{overflow-x:auto}code{overflow-wrap:anywhere}'''

def df_table(df):
    return '<div class="scroll">' + df.to_html(index=False, escape=True, border=0, na_rep='NA', float_format=lambda v: f'{v:.4g}') + '</div>'

def render_gene(gene_id):
    gene = DATA['genes'].loc[gene_id]
    windows = DATA['windows'][DATA['windows'].entrez_id == gene_id]
    reasons = DATA['reasons'][DATA['reasons'].entrez_id == gene_id]
    links = DATA['links'][DATA['links'].entrez_id == gene_id].copy()
    peaks = DATA['peaks'].loc[links.peak_id.unique()].copy()
    peaks['qualifying_for_selection'] = peaks.index.isin(reasons.qualifying_peak_id)
    evidence = links.groupby('peak_id').feature.agg(lambda v: ';'.join(sorted(set(v))))
    peaks['gene_link_feature'] = peaks.index.map(evidence)
    peaks = peaks.sort_values(['FDR', 'start'])
    table = peaks[['chrom', 'start', 'end', 'gene_link_feature', 'logFC', 'FDR', 'qualifying_for_selection']].join(
        DATA['counts'][SAMPLES].rename(columns={s: s + '_raw_count' for s in SAMPLES})).join(
        DATA['cpm'].rename(columns={s: s + '_TMM_CPM' for s in SAMPLES}))
    table.to_csv(OUT / f'data/{gene_id}_peak_measurements.csv')
    audits = []
    with PdfPages(OUT / f'figures/{gene_id}.pdf', metadata={'Title': f'{gene.symbol}: ATAC locus diagrams', 'Author': 'ATACseq analysis'}) as book:
        for window in windows.itertuples(index=False):
            fig, audit = draw_window(gene_id, gene, window)
            fig.savefig(OUT / f'figures/{window.window_id}.png', dpi=300, facecolor='white')
            book.savefig(fig, facecolor='white')
            plt.close(fig)
            audits.append(audit)
    effects = DATA['effects'][DATA['effects'].entrez_id == gene_id]
    effect_cols = ['RNA_result', 'feature', 'rna_logFC', 'rna_FDR', 'atac_median_logFC', 'minimum_peak_FDR', 'tested_peaks']
    view = [f'<!doctype html><html lang="en"><meta charset="utf-8"><title>{esc(gene.symbol)} ATAC locus</title><style>{STYLE}</style><body>',
            f'<p><a href="../../../report/gene_atlas.html">Back to searchable atlas</a></p><h1>{esc(gene.symbol)} — ATAC locus diagrams</h1>',
            f'<p>Entrez {gene_id}; {gene.chrom}:{gene.start:,}-{gene.end:,} ({gene.strand} strand; gene annotation coordinates are 1-based closed).</p>',
            f'<p><a href="../figures/{gene_id}.pdf">Download all views as vector PDF</a> | <a href="../data/{gene_id}_peak_measurements.csv">All linked peaks: raw counts, TMM CPM and DAR statistics</a></p>',
            '<h2>Why this gene is included</h2>', df_table(reasons[['feature', 'reason', 'RNA_result', 'qualifying_peak_id']]),
            '<h2>RNA and ATAC effect summaries</h2>', df_table(effects[effect_cols]),
            '<p>Empty effect summaries mean no eligible RNA integration row. Source RNA estimates, when uniquely mapped, are shown in the diagrams. RNA absence is not evidence of no expression.</p>',
            '<p class="note">All six ATAC replicates share one coverage scale within each view. The scale can differ between views and genes. Coverage is normalized per million non-Y fragments; peak dots use the primary edgeR TMM effective library sizes. These are distinct normalizations. Original FDR values are reused, not recalculated within each locus.</p>']
    for window in windows.itertuples(index=False):
        view += [f'<h2>{esc(window.view)}</h2><figure><a href="../figures/{window.window_id}.png"><img loading="lazy" src="../figures/{window.window_id}.png" alt="{esc(gene.symbol)} {esc(window.view)}"></a><figcaption><a href="../figures/{window.window_id}.png">300 dpi PNG</a></figcaption></figure>']
    display = table.reset_index()[['peak_id', 'gene_link_feature', 'logFC', 'FDR', 'qualifying_for_selection']]
    view += ['<h2>All linked tested peaks</h2>', df_table(display),
             '<p>Every selection-supporting peak is covered by a close-up. Gene overviews include the complete selected gene and all selection-supporting peaks. Nearby genes are shown for context; overlap or proximity does not establish regulation. The gene track shows merged exons across isoforms, not a single transcript.</p></body></html>']
    (OUT / f'genes/{gene_id}.html').write_text('\n'.join(view), encoding='utf-8')
    (OUT / f'data/{gene_id}_render_audit.json').write_text(json.dumps({'code_sha256': DATA['code_hash'], 'views': audits}, indent=2))
    return gene_id, len(audits)

def publish_atlas():
    genes = DATA['genes'].reset_index()
    summary = json.loads((OUT / 'selection_summary.json').read_text())
    rows = []
    for row in genes.sort_values('symbol').to_dict('records'):
        gene_id = row['entrez_id']
        assert (OUT / f'genes/{gene_id}.html').exists()
        flags = []
        if row['strong_DAR']:
            flags.append('Strong DAR')
        flags.extend(RNA_NAMES[name] for name in RNA_NAMES if row['overlap_' + name])
        rows.append(f'<tr data-strong="{int(row["strong_DAR"])}" data-feature="{row["features"]}" data-rna="' +
                    ';'.join(name for name in RNA_NAMES if row['overlap_' + name]) + '">' +
                    f'<td><a href="../integration/gene_atlas/genes/{gene_id}.html">{esc(row["symbol"])}</a></td><td>{gene_id}</td><td>{row["chrom"]}</td>' +
                    f'<td>{esc(row["features"])}</td><td>{esc("; ".join(flags))}</td><td>{row["strongest_linked_FDR"]:.3g}</td><td>{row["qualifying_peaks"]}</td>' +
                    f'<td><a href="../integration/gene_atlas/figures/{gene_id}.pdf">PDF</a></td></tr>')
    intro = f'{summary["selected_genes"]} genes; {summary["windows"]} locus views. {summary["overlap_genes"]} genes meet ATAC/RNA overlap criteria; {summary["strong_DAR_genes"]} genes link to {summary["strong_DAR_peaks"]} strong DARs. Categories overlap.'
    text = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>ATAC gene-locus atlas</title>', f'<style>{STYLE}</style><body>',
            '<h1>ATAC gene-locus atlas</h1><p>' + intro + '</p>',
            f'<p><strong>Strong DAR:</strong> peak FDR &lt;{summary["strong_FDR"]} and |log2FC| &gt;{summary["strong_abs_logFC"]}. <strong>ATAC/RNA overlap:</strong> any linked ATAC peak FDR &lt;0.10 and RNA gene FDR &lt;0.05, without an effect-size or direction filter. RNA scope: {esc(summary["RNA_scope"])}.</p>',
            '<p>Each gene has a full-locus overview and close-ups covering every qualifying peak. Figures show six ATAC replicates, group means, peak log2 fold changes, merged-exon gene models, RNA results from the combined and two separate fits, and per-replicate TMM peak measurements. Download vector PDFs and 300 dpi PNGs from the gene pages.</p>',
            '<p class="note">Promoter, gene-body and intergenic links are included and marked separately. Intergenic nearest-TSS links remain exploratory. Coverage tracks use non-Y fragment CPM; DAR testing and peak dots use TMM normalization. All samples share the coverage y-axis within each view. The primary ATAC model treats samples as independent replicates; no copy-number adjustment is applied.</p>',
            '<p><input id="search" type="search" placeholder="Search gene, Entrez ID or chromosome" aria-label="Search genes">',
            '<select id="strong" aria-label="Gene selection"><option value="">All selected genes</option><option value="1">Strong DAR genes</option></select>',
            '<select id="feature" aria-label="Link type"><option value="">All link types</option><option>promoter</option><option>gene_body</option><option>intergenic</option></select>',
            '<select id="rna" aria-label="RNA overlap support"><option value="">Any RNA support</option><option value="Combined">Combined overlap</option><option value="AUG_2025">August 2025 overlap</option><option value="MAY_2026">May 2026 overlap</option></select></p>',
            '<p id="count"></p><table id="atlas"><thead><tr><th>Gene</th><th>Entrez</th><th>Chr</th><th>Selection link types</th><th>Selection support</th><th>Minimum linked peak FDR</th><th>Qualifying peaks</th><th>Download</th></tr></thead><tbody>',
            *rows, '</tbody></table>',
            '<p>Feature and RNA filters mark gene-level inclusion and need not refer to the same peak; exact feature-by-result evidence is listed on each gene page. Minimum linked FDR is a peak statistic, not a calibrated gene-level FDR.</p>',
            '<h2>Data and methods</h2><ul>']
    for name, label in [('selected_genes.csv', 'Selected genes'), ('selection_reasons.csv', 'Every inclusion reason and supporting peak'),
                        ('RNA_source_rows.csv', 'RNA source rows'), ('all_selected_gene_peak_links.csv', 'All selected gene–peak links'),
                        ('plot_windows.bed', 'Plot regions, BED'), ('source_manifest.csv', 'Input checksums'), ('coverage_source_manifest.csv', 'BigWig checksums'),
                        ('render_audit.csv', 'Rendering audit')]:
        text.append(f'<li><a href="../integration/gene_atlas/data/{name}">{label}</a></li>')
    text += ['</ul><p><a href="../integration/gene_atlas/METHODS.txt">Complete methods</a> | <a href="ATACseq_report.html">Main ATAC report</a></p>',
             '''<script>const rows=[...document.querySelectorAll('#atlas tbody tr')];function filter(){let n=0;const q=document.getElementById('search').value.toLowerCase(),s=document.getElementById('strong').value,f=document.getElementById('feature').value,r=document.getElementById('rna').value;for(const row of rows){const show=row.textContent.toLowerCase().includes(q)&&(!s||row.dataset.strong===s)&&(!f||row.dataset.feature.split(';').includes(f))&&(!r||row.dataset.rna.split(';').includes(r));row.hidden=!show;if(show)n++;}document.getElementById('count').textContent=n+' of '+rows.length+' genes shown';}for(const id of ['search','strong','feature','rna'])document.getElementById(id).addEventListener('input',filter);filter();</script></body></html>''']
    (ROOT / 'report/gene_atlas.html').write_text('\n'.join(text), encoding='utf-8')
    snippet = '<section id="gene-locus-atlas"><h2>ATAC gene-locus diagrams</h2><p>' + intro + ' <a href="gene_atlas.html">Open the searchable gene atlas, coverage diagrams and peak measurements</a>.</p></section>'
    main = ROOT / 'report/ATACseq_report.html'
    content = main.read_text(encoding='utf-8')
    if '<section id="gene-locus-atlas">' in content:
        content = re.sub(r'<section id="gene-locus-atlas">.*?</section>', lambda _: snippet, content, flags=re.S)
    else:
        content = content.replace('</body>', snippet + '\n</body>')
    main.write_text(content, encoding='utf-8')
    print(intro, flush=True)

def validate_and_publish():
    genes = DATA['genes']
    windows = DATA['windows']
    audits = []
    for gene_id in genes.index:
        result = json.loads((OUT / f'data/{gene_id}_render_audit.json').read_text())
        assert result['code_sha256'] == DATA['code_hash']
        audits.extend(result['views'])
        assert (OUT / f'figures/{gene_id}.pdf').stat().st_size > 1000
    assert len(audits) == len(windows)
    assert set(a['window_id'] for a in audits) == set(windows.window_id)
    for window in windows.itertuples():
        assert (OUT / f'figures/{window.window_id}.png').stat().st_size > 1000
        for pid in window.qualifying_peak_ids.split(';'):
            p = DATA['peaks'].loc[pid]
            assert window.chrom == p.chrom and window.start <= p.start and window.end >= p.end
    pd.DataFrame(audits).to_csv(OUT / 'data/render_audit.csv', index=False)
    methods = ['Selection thresholds and RNA scope are recorded in selection_summary.json. No selected gene is omitted.',
               'ATAC/RNA overlap is gene-level membership based on any linked significant peak and the original RNA FDR. There is no joint FDR or effect-size filter for overlap.',
               'Strong DAR gene assignment reuses all original promoter, gene-body and intergenic nearest-TSS links, including genes unavailable in RNA.',
               'Primary ATAC edgeR QL ~ condition statistics are reused. Primary TMM effective library sizes give peak CPM = count / effective_size * 1e6; all six replicates are shown, with arithmetic group means.',
               'Coverage comes from the six original 25 bp fragment-span BigWigs normalized by accepted non-Y fragment counts. No Tn5 shift, footprint correction or new coverage normalization is applied.',
               'Coverage is displayed as exact averages over at least 25 bp bins, approximately 1800 bins maximum per view. Missing BigWig bases are zero. Extraction verifies the signal integral is preserved.',
               'All six replicates and the group-mean track use the same zero-based y-axis limit within a view, without clipping. Y-axis scales may differ between views and genes.',
               'Gene overview includes the complete gene plus all selection-supporting peaks and 10 kb flanks. Every qualifying peak receives a 5 kb-flanked close-up; overlapping close-ups are merged.',
               'Plot coordinates are hg38, 0-based half-open. Gene model input is 1-based closed and is converted explicitly. Exons are unioned across isoforms, not a single canonical transcript.',
               'All same-annotation neighboring genes intersecting a window are shown. The selected gene is purple with an asterisk. A close-up may be distal to the selected gene and therefore not show its body.',
               'The DAR-effect track includes all tested peaks in the window, with black outlines for links to the selected gene. Qualifying peaks are shaded purple. Top three qualifying peaks by FDR have numbered peak-CPM panels; complete measurements for every linked peak are in CSV.',
               'RNA results are existing combined DESeq2 ~ batch + treatment and separate AUG_2025/MAY_2026 ~ treatment fits. Ambiguous Entrez mappings are not summarized as a single RNA effect.',
               'Intergenic nearest-TSS assignments do not establish enhancer regulation. Gene-body overlap also does not establish a regulatory relationship. A peak can link to multiple genes.',
               'No copy-number correction or causal attribution is made. RNA and ATAC are integrated at group level; direct sample matching has not been established.',
               'Figures: matplotlib, individual 300 dpi PNGs and a multipage vector PDF per gene. Gene atlas HTML is searchable locally and has no external dependency.']
    (OUT / 'METHODS.txt').write_text('\n'.join(methods), encoding='utf-8')
    publish_atlas()
    paths = [ROOT / 'report/gene_atlas.html'] + list((OUT / 'genes').glob('*.html'))
    link_count = 0
    for path in paths:
        for link in re.findall(r'(?:href|src)="([^"]+)"', path.read_text(encoding='utf-8')):
            assert (path.parent / link).resolve().exists(), (path, link)
            link_count += 1
    complete = {'status': 'complete', 'genes': len(genes), 'views': len(windows), 'local_links_verified': link_count,
                'validation': 'All selected genes, PDFs, PNGs, plot windows, qualifying-peak coverage, source-value extraction and local links passed',
                'render_code_sha256': DATA['code_hash'], 'matplotlib': matplotlib.__version__, 'numpy': np.__version__, 'pandas': pd.__version__}
    (OUT / 'ANALYSIS_COMPLETE.json').write_text(json.dumps(complete, indent=2))
    print(json.dumps(complete, indent=2), flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--genes', nargs='+')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--publish-only', action='store_true')
    args = parser.parse_args()
    init_worker()
    jobs = list(DATA['genes'].index)
    if args.genes:
        jobs = [i for i in jobs if i in args.genes or str(DATA['genes'].loc[i, 'symbol']) in args.genes]
        assert jobs, args.genes
    if args.resume:
        remaining = []
        for gene_id in jobs:
            path = OUT / f'data/{gene_id}_render_audit.json'
            if not path.exists() or json.loads(path.read_text()).get('code_sha256') != DATA['code_hash']:
                remaining.append(gene_id)
        jobs = remaining
    if not args.publish_only:
        with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker) as executor:
            pending = {executor.submit(render_gene, gene_id): gene_id for gene_id in jobs}
            for i, future in enumerate(as_completed(pending), 1):
                gene_id, n = future.result()
                print(f'Rendered {i}/{len(jobs)}: {DATA["genes"].loc[gene_id, "symbol"]} ({n} views)', flush=True)
    if not args.genes:
        validate_and_publish()

if __name__ == '__main__':
    main()
