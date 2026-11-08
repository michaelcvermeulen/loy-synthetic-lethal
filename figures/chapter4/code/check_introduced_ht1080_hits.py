"""Tabulate saved results for candidates introduced in Chapter 4; no refitting."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
FIG = ROOT / 'SEPT_11_2026_THESIS/figures/chapter4'
OUT = FIG / 'source_data/ht1080_introduced_hits_20261022'
SOURCES = {
    'primary_drugZ': ROOT / 'crispr_screen/ht1080_loy_screen_2/results/drugZ/LOY_vs_CTRL_unfiltered.txt',
    'MAGeCK': ROOT / 'crispr_screen/ht1080_loy_screen_2/test/POST_LOY_vs_CTRL.gene_summary.txt',
    'IsogenicZ': FIG / 'source_data/ht1080_isogenicz_20261020/mask173_detected_IsogenicZ_genes.tsv',
    'saved_comparison': FIG / 'source_data/ht1080_candidate_discovery/method_comparison.csv',
}
GENES = ['DDX3X', 'TIMM8A', 'MED12', 'SCLY', 'UBE2G2', 'ITGB5', 'SREBF2', 'SREBF1', 'LLGL1',
         'TAF6L', 'TADA2B', 'SUPT20H', 'ACSL4', 'ATXN7L3', 'TADA1',
         'KMT2D', 'TAF5L', 'KAT2A']

def main():
    d = pd.read_csv(SOURCES['primary_drugZ'], sep='\t', index_col=0)
    m = pd.read_csv(SOURCES['MAGeCK'], sep='\t', index_col=0)
    i = pd.read_csv(SOURCES['IsogenicZ'], sep='\t', index_col=0)
    saved = pd.read_csv(SOURCES['saved_comparison'], index_col=0)
    for df in [d, m, i, saved]:
        assert df.index.is_unique and set(GENES).issubset(df.index)
    rows = []
    for gene in GENES:
        x, y, z = d.loc[gene], m.loc[gene], i.loc[gene]
        assert np.isclose(x.fdr_synth, saved.loc[gene, 'fdr_synth'])
        assert np.isclose(y['neg|fdr'], saved.loc[gene, 'neg|fdr'])
        rows.append(dict(
            gene=gene, role='functional context only' if gene == 'KAT2A' else 'candidate',
            drugZ_rank=int(x.rank_synth), drugZ_normZ=x.normZ, drugZ_FDR=x.fdr_synth,
            MAGeCK_rank=int(y['neg|rank']), MAGeCK_FDR=y['neg|fdr'], MAGeCK_lfc=y['neg|lfc'],
            IsogenicZ_rank=int(z.rank_synth), IsogenicZ_normZ=z.normZ, IsogenicZ_FDR=z.fdr_synth,
        ))
    table = pd.DataFrame(rows)
    for method in ['drugZ', 'MAGeCK', 'IsogenicZ']:
        table[method + '_FDR_lt_005'] = table[method + '_FDR'] < .05
        table[method + '_FDR_lt_025'] = table[method + '_FDR'] < .25
    assert (table.loc[table.role.eq('candidate'), 'drugZ_FDR'] < .25).all()
    assert (table.drugZ_normZ < 0).all() and (table.IsogenicZ_normZ < 0).all()
    OUT.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT / 'candidate_method_support.csv', index=False)
    manifest = dict(
        sources={name: dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                 for name, path in SOURCES.items()},
        interpretation='Saved depletion-tail FDRs with original testing families; no refitting or subset FDR adjustment. Agreement across methods is not independent replication.',
        IsogenicZ_branch='mask173_detected: 173 D0 outliers and all-library-undetected guides excluded; no 50-read D0 filter',
        KAT2A='Checked as functional context for SAGA, not a nominated hit.',
    )
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(table[['gene', 'drugZ_FDR', 'MAGeCK_FDR', 'IsogenicZ_FDR']].to_string(index=False))

if __name__ == '__main__':
    main()
