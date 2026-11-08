"""Select gene loci using existing DAR and RNA results; preserve every inclusion reason."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'integration/gene_atlas'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--strong-fdr', type=float, default=.01)
    parser.add_argument('--strong-lfc', type=float, default=1.)
    parser.add_argument('--rna-scope', choices=['all', 'combined', 'separate'], default='all')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for folder in ['data', 'figures', 'genes', 'coverage']:
        (OUT / folder).mkdir(exist_ok=True)
    peak_path = ROOT / 'annotation/annotated_differential_peaks.csv'
    link_path = ROOT / 'integration/feature_overlaps/peak_gene_link_audit.csv'
    peaks = pd.read_csv(peak_path, dtype={'nearest_entrez_id': str})
    assert peaks.peak_id.is_unique
    peaks['strong_DAR'] = (peaks.FDR < args.strong_fdr) & (peaks.logFC.abs() > args.strong_lfc)
    links = pd.read_csv(link_path, dtype={'entrez_id': str})[
        ['feature', 'peak_id', 'chrom', 'entrez_id', 'atac_logFC', 'atac_FDR', 'link_evidence']].drop_duplicates()
    assert not links.duplicated(['feature', 'peak_id', 'entrez_id']).any()
    reasons = []
    strong = links[links.peak_id.isin(peaks.loc[peaks.strong_DAR, 'peak_id'])].copy()
    for row in strong.itertuples():
        reasons.append({'entrez_id': row.entrez_id, 'feature': row.feature, 'reason': 'Strong DAR', 'RNA_result': '', 'qualifying_peak_id': row.peak_id})
    combined_path = ROOT / 'integration/feature_scatters/feature_gene_effects.csv'
    batch_path = ROOT / 'integration/batch_RNA_overlaps/all_gene_results.csv'
    combined = pd.read_csv(combined_path, dtype={'entrez_id': str}).assign(RNA_result='Combined')
    batch = pd.read_csv(batch_path, dtype={'entrez_id': str}).rename(columns={'batch': 'RNA_result'})
    effects = pd.concat([combined, batch], ignore_index=True)
    sources = ['Combined', 'AUG_2025', 'MAY_2026'] if args.rna_scope == 'all' else (['Combined'] if args.rna_scope == 'combined' else ['AUG_2025', 'MAY_2026'])
    overlap = effects[(effects.RNA_result.isin(sources)) & (effects.minimum_peak_FDR < .1) & (effects.rna_FDR < .05)].copy()
    for row in overlap.itertuples():
        qualifying = links[(links.entrez_id == row.entrez_id) & (links.feature == row.feature) & (links.atac_FDR < .1)]
        assert len(qualifying)
        for p in qualifying.itertuples():
            reasons.append({'entrez_id': row.entrez_id, 'feature': row.feature, 'reason': 'ATAC/RNA overlap', 'RNA_result': row.RNA_result, 'qualifying_peak_id': p.peak_id})
    reasons = pd.DataFrame(reasons).drop_duplicates()
    selected = set(reasons.entrez_id)
    genes_path = ROOT / 'config/genes_knownGene.tsv'
    genes = pd.read_csv(genes_path, sep='\t', dtype={'gene_id': str}).rename(columns={'gene_id': 'entrez_id'})
    assert genes.entrez_id.is_unique and selected.issubset(set(genes.entrez_id))
    symbols = dict(zip(effects.entrez_id, effects.rna_symbol))
    for row in peaks.itertuples():
        if row.nearest_entrez_id not in symbols and pd.notna(row.nearest_symbol):
            symbols[row.nearest_entrez_id] = row.nearest_symbol
    model_path = OUT / 'data/gene_models.csv'
    if model_path.exists():
        for row in pd.read_csv(model_path, dtype={'gene_id': str}).itertuples():
            if row.gene_id not in symbols or pd.isna(symbols[row.gene_id]):
                symbols[row.gene_id] = row.symbol
    genes['symbol'] = genes.entrez_id.map(symbols).fillna(genes.entrez_id)
    atlas = genes[genes.entrez_id.isin(selected)].copy()
    atlas['strong_DAR'] = atlas.entrez_id.isin(strong.entrez_id)
    for label in ['Combined', 'AUG_2025', 'MAY_2026']:
        atlas['overlap_' + label] = atlas.entrez_id.isin(overlap.loc[overlap.RNA_result == label, 'entrez_id'])
        atlas['promoter_overlap_' + label] = atlas.entrez_id.isin(overlap.loc[(overlap.RNA_result == label) & (overlap.feature == 'promoter'), 'entrez_id'])
    atlas['features'] = atlas.entrez_id.map(reasons.groupby('entrez_id').feature.agg(lambda v: ';'.join(sorted(set(v)))))
    atlas['strongest_linked_FDR'] = atlas.entrez_id.map(links.groupby('entrez_id').atac_FDR.min())
    atlas['qualifying_peaks'] = atlas.entrez_id.map(reasons.groupby('entrez_id').qualifying_peak_id.nunique())
    atlas = atlas.sort_values(['strongest_linked_FDR', 'symbol', 'entrez_id'])
    windows = []
    peak_index = peaks.set_index('peak_id')
    sizes = dict(line.split() for line in (ROOT / 'config/bam.chrom.sizes').read_text().splitlines())
    sizes = {k: int(v) for k, v in sizes.items()}
    for row in atlas.itertuples():
        chosen = peak_index.loc[reasons.loc[reasons.entrez_id == row.entrez_id, 'qualifying_peak_id'].unique()].sort_values('start')
        assert set(chosen.chrom) == {row.chrom}
        # Gene annotation is 1-based closed; plotting and BED exports use 0-based half-open intervals.
        left = max(0, min(row.start - 1, int(chosen.start.min())) - 10000)
        right = min(sizes[row.chrom], max(row.end, int(chosen.end.max())) + 10000)
        windows.append({'window_id': f'{row.entrez_id}_overview', 'entrez_id': row.entrez_id, 'symbol': row.symbol,
                        'chrom': row.chrom, 'start': left, 'end': right, 'view': 'Gene overview',
                        'qualifying_peak_ids': ';'.join(chosen.index)})
        # Every qualifying peak receives a close-up; overlapping 10 kb windows are merged.
        clusters = []
        for pid, peak in chosen.iterrows():
            lo = max(0, int(peak.start) - 5000); hi = min(sizes[row.chrom], int(peak.end) + 5000)
            if clusters and lo <= clusters[-1][1]:
                clusters[-1][1] = max(hi, clusters[-1][1]); clusters[-1][2].append(pid)
            else:
                clusters.append([lo, hi, [pid]])
        for index, (lo, hi, ids) in enumerate(clusters, 1):
            windows.append({'window_id': f'{row.entrez_id}_peak{index}', 'entrez_id': row.entrez_id, 'symbol': row.symbol,
                            'chrom': row.chrom, 'start': lo, 'end': hi, 'view': f'Peak close-up {index}', 'qualifying_peak_ids': ';'.join(ids)})
    windows = pd.DataFrame(windows)
    atlas.to_csv(OUT / 'data/selected_genes.csv', index=False)
    reasons.to_csv(OUT / 'data/selection_reasons.csv', index=False)
    windows.to_csv(OUT / 'data/plot_windows.csv', index=False)
    windows[['chrom', 'start', 'end', 'window_id']].to_csv(OUT / 'data/plot_windows.bed', sep='\t', index=False, header=False)
    peaks[peaks.strong_DAR].to_csv(OUT / 'data/strong_DAR_peaks.csv', index=False)
    links[links.entrez_id.isin(selected)].to_csv(OUT / 'data/all_selected_gene_peak_links.csv', index=False)
    effects[effects.entrez_id.isin(selected)].to_csv(OUT / 'data/selected_gene_RNA_ATAC_effects.csv', index=False)
    # Store RNA statistics even when a selected strong-DAR gene is not assay-eligible.
    rnafiles = {'Combined': ROOT / 'config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv',
                **{b: ROOT / f'integration/batch_RNA_overlaps/sources/ht1080_{b}_LOY_v_CTRL_DESeq2.csv' for b in ['AUG_2025', 'MAY_2026']}}
    rnalist = []
    for name, path in rnafiles.items():
        r = pd.read_csv(path, dtype={'entrezid': str})
        r['entrez_multiplicity'] = r.groupby('entrezid').entrezid.transform('size')
        rnalist.append(r[r.entrezid.isin(selected)].assign(RNA_result=name))
    pd.concat(rnalist, ignore_index=True).to_csv(OUT / 'data/RNA_source_rows.csv', index=False)
    for fname in ['peak_counts.tsv', 'logCPM_primary.tsv']:
        c = pd.read_csv(ROOT / 'counts' / fname, sep='\t')
        c[c.peak_id.isin(links.loc[links.entrez_id.isin(selected), 'peak_id'])].to_csv(OUT / 'data' / fname.replace('.tsv', '.csv'), index=False)
    allpaths = [peak_path, link_path, combined_path, batch_path, genes_path, *rnafiles.values(), ROOT / 'differential/normalization_factors.csv']
    manifest = []
    for path in allpaths:
        manifest.append({'path': str(path.relative_to(ROOT)), 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()})
    pd.DataFrame(manifest).to_csv(OUT / 'data/source_manifest.csv', index=False)
    summary = {'strong_FDR': args.strong_fdr, 'strong_abs_logFC': args.strong_lfc, 'overlap_ATAC_FDR': .1,
               'overlap_RNA_FDR': .05, 'RNA_scope': args.rna_scope, 'selected_genes': len(atlas),
               'strong_DAR_peaks': int(peaks.strong_DAR.sum()), 'strong_DAR_genes': int(atlas.strong_DAR.sum()),
               'overlap_genes': overlap.entrez_id.nunique(), 'windows': len(windows),
               'overlap_genes_by_RNA_result': overlap.groupby('RNA_result').entrez_id.nunique().to_dict(),
               'overlap_genes_by_feature': overlap.groupby('feature').entrez_id.nunique().to_dict()}
    (OUT / 'selection_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
