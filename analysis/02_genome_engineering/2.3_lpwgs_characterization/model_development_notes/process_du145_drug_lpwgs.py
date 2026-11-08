"""Summarise existing DU145 treatment lpWGS outputs; no alignment or refitting."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
WGS = ROOT.parent / 'WGS'
OUT = ROOT / 'figures/chapter3/source_data/du145_drug_optimisation/lpwgs'
LABELS = {
    'LNCaP-sgCTRL': 'sgCTRL', 'LNCaP-LOY-2X': 'LOY 2x',
    'LNCaP-DMSO-LOY-1X': 'DMSO vehicle', 'LNCaP-30': 'Panobinostat 30 nM',
    'LNCaP-70': 'Panobinostat 75 nM', 'LNCaP-KU': 'KU-55933',
    'LNCaP-KU-30': 'KU-55933 + panobinostat 30 nM',
    'LNCaP-KU-70': 'KU-55933 + panobinostat 75 nM',
}

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = [WGS / 'LNCaP_QDNAseq_cn_df.rds',
              WGS / 'code/coverage_report_LNCaP_QDNAseq.Rmd',
              WGS / 'code/coverage_report_LNCaP_QDNAseq.html']
    r = f"d <- readRDS('{inputs[0].as_posix()}')\nwrite.csv(d[d$chromosome == 'chrY',], '{(OUT / 'chrY_qdnaseq_bins.csv').as_posix()}', row.names=FALSE)\n"
    subprocess.run(['Rscript', '--vanilla', '-'], input=r, text=True, check=True)
    d = pd.read_csv(OUT / 'chrY_qdnaseq_bins.csv')
    assert set(d['sample']) == set(LABELS)
    rows, regional = [], []
    for sample, label in LABELS.items():
        p = WGS / 'coverage' / f'{sample}.genome.mosdepth.summary.txt'
        inputs.append(p)
        shutil.copy2(p, OUT / p.name)
        t = pd.read_csv(p, sep='\t').set_index('chrom')
        a = t.loc[[f'chr{i}' for i in range(1, 23)]]
        auto = a.bases.sum() / a.length.sum()
        y = t.loc['chrY', 'bases'] / t.loc['chrY', 'length']
        rows.append(dict(sample=sample, condition=label, autosome_depth=auto,
                         chrY_depth=y, chrY_autosome_ratio=y / auto))
        g = d[d['sample'].eq(sample) & d['use'].eq(True)]
        for name, lo, hi in [('distal_MSY_Yp_summary', 2.46, 4.8),
                             ('proximal_Yp_summary', 6, 10),
                             ('euchromatic_Yq_summary', 13, 27.449937)]:
            x = g[g.midpoint.ge(lo) & g.midpoint.lt(hi)]
            regional.append(dict(sample=sample, condition=label, region=name,
                start_Mb=lo, end_Mb=hi, n_bins=len(x),
                median_corrected_copy_ratio=x.copynumber.median(),
                median_segmented_copy_ratio=(2 ** x.log2_seg).median()))
    table = pd.DataFrame(rows)
    vehicle = table.loc[table.condition.eq('DMSO vehicle'), 'chrY_autosome_ratio'].iloc[0]
    table['percent_lower_than_vehicle'] = 100 * (1 - table.chrY_autosome_ratio / vehicle)
    table.to_csv(OUT / 'raw_depth_ratios.csv', index=False)
    pd.DataFrame(regional).to_csv(OUT / 'regional_summary.csv', index=False)
    pd.read_html(inputs[2])[0].to_csv(OUT / 'saved_report_read_counts.csv', index=False)
    manifest = {
        'inputs': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs],
        'notes': ['Historical LNCaP labels refer to DU145; 70 is displayed as author-confirmed 75 nM.',
                  'One saved WGS sample per condition; ddPCR same-specimen matching is unresolved.',
                  'Raw depth ratios include repeats and PARs; not absolute copies or Y-negative cell fractions.',
                  'Regional windows are descriptive summaries, not newly inferred deletion boundaries.',
                  'Unsegmented copynumber is retained without flooring; segmented ratios inherit the original 0.001 floor.',
                  'All eight samples exist in saved RDS and HTML; stale low-depth KU comment is unsupported by current summaries.']}
    (OUT / 'provenance.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(table.to_string(index=False))

if __name__ == '__main__':
    main()
