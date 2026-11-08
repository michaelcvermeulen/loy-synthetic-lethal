"""Audit strict zero-read guides in the retained HT1080 and A549 screens."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT.parent / '.figure_python_deps'))
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
SOURCES = {
    'HT1080': ROOT.parent / 'crispr_screen/ht1080_loy_screen_2/count/ht1080_screen2_R1.count.txt',
    'A549': ROOT.parent / 'crispr_screen/a549_loy_screen1/CRISPR_screen/Screen1/Mageck/A549_full.count.txt',
}
SAMPLES = {
    'HT1080': ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2', 'POST_LOY_REP1', 'POST_LOY_REP2'],
    'A549': ['CTRL_D0_1', 'CTRL_D0_2', 'LOY_D0_1', 'LOY_D0_2', 'CTRL_P4_R2', 'LOY_P4_R1', 'LOY_P4_R2'],
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    hashes = {line: digest(path) for line, path in SOURCES.items()}
    tables = {line: pd.read_csv(path, sep='\t').set_index('sgRNA').sort_index()
              for line, path in SOURCES.items()}
    ht, a549 = tables['HT1080'], tables['A549']
    pd.testing.assert_index_equal(ht.index, a549.index)
    pd.testing.assert_series_equal(ht.Gene, a549.Gene)
    assert len(ht) == 77441 and ht.index.is_unique
    audit = ht[['Gene']].copy()
    audit['sequence'] = audit.index.str.rsplit('_', n=1).str[-1]
    assert audit.sequence.str.fullmatch('[ACGT]{20}').all()
    metrics = []
    endpoint_loss = []
    for line, raw in tables.items():
        samples = SAMPLES[line]
        assert not any('_P3_' in sample for sample in samples)
        x = raw[samples]
        assert np.isfinite(x).all().all() and (x >= 0).all().all()
        assert np.equal(x, np.floor(x)).all().all()
        d0 = [sample for sample in samples if 'D0' in sample]
        end = [sample for sample in samples if 'D0' not in sample]
        zero = x.eq(0).all(axis=1)
        zero_d0 = x[d0].eq(0).all(axis=1)
        zero_end = x[end].eq(0).all(axis=1)
        assert np.array_equal(zero, zero_d0 & zero_end)
        assert np.array_equal(zero, x.sum(axis=1).eq(0))
        audit[f'{line}_zero_all'] = zero
        audit[f'{line}_zero_D0'] = zero_d0
        audit[f'{line}_zero_endpoint'] = zero_end
        audit[f'{line}_D0_total_reads'] = x[d0].sum(axis=1)
        audit[f'{line}_endpoint_total_reads'] = x[end].sum(axis=1)
        for sample in samples:
            audit[f'{line}__{sample}'] = x[sample]
        metrics.append({
            'screen': line, 'n_retained_libraries': len(samples), 'n_guides': len(x),
            'zero_all_libraries': int(zero.sum()), 'zero_all_percent': float(zero.mean() * 100),
            'zero_all_D0_libraries': int(zero_d0.sum()),
            'zero_all_endpoint_libraries': int(zero_end.sum()),
            'detected_at_D0_but_zero_all_endpoints': int((~zero_d0 & zero_end).sum()),
            'zero_at_D0_but_detected_at_endpoint': int((zero_d0 & ~zero_end).sum()),
        })
        losses = audit.loc[~zero_d0 & zero_end, ['Gene', 'sequence']].copy()
        losses.insert(0, 'screen', line)
        losses['D0_total_reads'] = x.loc[losses.index, d0].sum(axis=1)
        losses['endpoint_total_reads'] = 0
        endpoint_loss.append(losses.reset_index())
    both = audit.HT1080_zero_all & audit.A549_zero_all
    either = audit.HT1080_zero_all | audit.A549_zero_all
    common = audit.loc[both].copy()
    common.to_csv(OUT / 'guides_undetected_in_both_screens.csv', index_label='sgRNA')
    audit.loc[either].to_csv(OUT / 'guides_undetected_in_either_screen.csv', index_label='sgRNA')
    pd.DataFrame(metrics).to_csv(OUT / 'screen_zero_read_summary.csv', index=False)
    pd.concat(endpoint_loss, ignore_index=True).to_csv(OUT / 'detected_D0_zero_endpoint_guides.csv', index=False)
    gene_coverage = audit.groupby('Gene').agg(
        total_guides=('sequence', 'size'),
        HT1080_zero_guides=('HT1080_zero_all', 'sum'),
        A549_zero_guides=('A549_zero_all', 'sum'))
    gene_coverage['shared_zero_guides'] = both.groupby(audit.Gene).sum()
    gene_coverage.loc[gene_coverage.HT1080_zero_guides.gt(0) | gene_coverage.A549_zero_guides.gt(0)].to_csv(
        OUT / 'affected_gene_coverage.csv')
    ntc = common.Gene.eq('Non-Targeting-Control')
    # Post hoc descriptive sequence check; no mechanistic or inferential claim.
    prefix = audit.sequence.str.startswith('TCTC') | audit.sequence.str.startswith('AGACG')
    prefix_records = []
    for label, mask in [('TCTC', audit.sequence.str.startswith('TCTC')),
                        ('AGACG', audit.sequence.str.startswith('AGACG')),
                        ('either_prefix', prefix)]:
        prefix_records.append({'prefix': label, 'n_in_full_library': int(mask.sum()),
                               'n_shared_zero': int((mask & both).sum()),
                               'n_other_guides': int((mask & ~both).sum())})
    pd.DataFrame(prefix_records).to_csv(OUT / 'exploratory_sequence_prefixes.csv', index=False)
    focus = ['DDX3X', 'EIF1AX', 'ZFX', 'RPS4X', 'MED12', 'CDK9', 'SREBF1', 'SREBF2', 'TADA2B', 'TAF6L']
    gene_coverage.loc[focus].to_csv(OUT / 'selected_gene_coverage.csv')
    assert common.loc['DDX3X_4640_TCTCGGTTCCTTAAATGAGG', ['HT1080_zero_all', 'A549_zero_all']].all()
    assert hashes == {line: digest(path) for line, path in SOURCES.items()}
    summary = {
        'date': '2026-10-20', 'definition': 'Exactly zero raw counts in every retained library, including D0 and endpoint.',
        'sources': {line: {'path': str(path), 'sha256': hashes[line]} for line, path in SOURCES.items()},
        'samples': SAMPLES, 'per_screen': metrics,
        'shared_zero_guides': int(both.sum()), 'shared_zero_percent': float(both.mean() * 100),
        'other_shared_zero_guides_besides_DDX3X_4640': int(both.sum() - 1),
        'HT1080_only_zero_guides': int((audit.HT1080_zero_all & ~audit.A549_zero_all).sum()),
        'A549_only_zero_guides': int((~audit.HT1080_zero_all & audit.A549_zero_all).sum()),
        'either_screen_zero_guides': int(either.sum()),
        'shared_zero_gene_targeting_guides': int((~ntc).sum()),
        'shared_zero_target_genes': int(common.loc[~ntc, 'Gene'].nunique()),
        'shared_zero_non_targeting_guides': int(ntc.sum()),
        'max_shared_zero_guides_per_target_gene': int(common.loc[~ntc].groupby('Gene').size().max()),
        'shared_zero_with_duplicate_sequence': int(audit.sequence.duplicated(keep=False).loc[both].sum()),
        'exploratory_prefix_counts': prefix_records,
        'interpretation': 'Non-detection in existing count matrices; physical absence, mapping/reference failure and other technical causes are not distinguished. Absence at D0 is not measured loss during outgrowth.',
        'excluded_samples': 'All A549 P3 samples remain excluded. No absent CTRL P4 replicate is imputed.',
        'validation': {'guide_ids_and_gene_annotations_match': True, 'integer_nonnegative_counts': True,
                       'all_zero_equals_zero_sum': True, 'source_hashes_unchanged': True},
        'software': {'pandas': pd.__version__, 'numpy': np.__version__},
    }
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    text = f"""# Guides undetected in both screens

Audit date: 20 October 2026. The analysis uses six HT1080 libraries and seven retained A549 libraries; all A549 P3 libraries remain excluded.

Strict non-detection means zero raw counts in every retained D0 and endpoint library. Both original matrices contain the same 77,441 unique guide IDs and matching gene annotations.

| Scope | Guides with zero reads throughout |
|---|---:|
| HT1080 | {metrics[0]['zero_all_libraries']} |
| A549 | {metrics[1]['zero_all_libraries']} |
| Both screens | {int(both.sum())} ({both.mean()*100:.3f}% of the library) |
| HT1080 only | {summary['HT1080_only_zero_guides']} |
| A549 only | {summary['A549_only_zero_guides']} |

DDX3X_4640_TCTCGGTTCCTTAAATGAGG is one of the shared undetected guides. There are {int(both.sum()-1)} other such guides: the complete shared set comprises {int((~ntc).sum())} gene-targeting guides across {common.loc[~ntc,'Gene'].nunique()} genes and {int(ntc.sum())} non-targeting guide. This is a guide-level limitation and does not mean that all guides targeting the affected genes were absent.

[Full shared-zero list and raw counts](guides_undetected_in_both_screens.csv) · [Either-screen list](guides_undetected_in_either_screen.csv) · [Affected-gene coverage](affected_gene_coverage.csv) · [Selected genes](selected_gene_coverage.csv) · [Source hashes and definitions](summary.json).

## Separate endpoint-only losses

{metrics[0]['detected_at_D0_but_zero_all_endpoints']} HT1080 guides and {metrics[1]['detected_at_D0_but_zero_all_endpoints']} A549 guides had at least one D0 read but zero reads across all retained endpoints within that screen. These are separate from the guides never detected. Neither category alone establishes a biological cause. [Counts](detected_D0_zero_endpoint_guides.csv).

## Exploratory sequence observation

Of the {int(both.sum())} shared zero-read guides, {int((prefix & both).sum())} begin with TCTC or AGACG; {int(prefix.sum())} of all 77,441 reference guides have either prefix. DDX3X guide 4640 begins with TCTC. This is a post hoc descriptive check, not a mechanistic test. Shared-zero sequences are unique within the reference, so duplicate sequence rows do not explain this set. These count tables cannot distinguish physical library underrepresentation from sequence/reference or read-processing problems. No synthesis, cloning, mapping or trimming cause is assigned, and no count matrix, figure or hit analysis is altered. [Prefix counts](exploratory_sequence_prefixes.csv).
"""
    (OUT / 'README.md').write_text(text, encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k not in ['sources', 'samples']}, indent=2))
    print('\nSelected-gene coverage:\n' + gene_coverage.loc[focus].to_string())


if __name__ == '__main__':
    main()
