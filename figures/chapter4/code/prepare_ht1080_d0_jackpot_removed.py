"""Prepare the descriptive D0 correlation after the frozen jackpot exclusions."""
from pathlib import Path
import hashlib
import json
import sys

THESIS = Path(__file__).resolve().parents[3]
DEPS = THESIS.parent / '.figure_python_deps'
if DEPS.is_dir():
    sys.path.insert(0, str(DEPS))
import numpy as np
import pandas as pd

OUT = THESIS / 'figures/chapter4/source_data/ht1080_d0_jackpot_removed'
COUNTS = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2/count/ht1080_screen2_R1.count.txt'
FROZEN = THESIS / 'figures/chapter4/source_data/ht1080_figure1_qc/d0_outlier_guides.csv'
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2',
           'POST_LOY_REP1', 'POST_LOY_REP2']


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(COUNTS, sep='\t')
    assert len(counts) == 77441 and counts.sgRNA.is_unique
    raw = counts[SAMPLES].astype(float)
    cpm = raw.div(raw.sum(), axis=1) * 1e6
    log = np.log2(cpm[SAMPLES[:2]] + 1)
    d0_mean = cpm[SAMPLES[:2]].mean(axis=1)
    end_mean = cpm[SAMPLES[2:]].mean(axis=1)
    jackpot = (d0_mean > 200) & (d0_mean / end_mean.replace(0, np.nan) > 50)
    z = (log - log.mean()) / log.std(ddof=1)
    outliers = z.abs().gt(3).any(axis=1)
    assert jackpot.sum() == 164 and outliers.sum() == 173
    assert not (jackpot & ~outliers).any()
    assert set(counts.loc[outliers, 'sgRNA']) == set(pd.read_csv(FROZEN).sgRNA)
    rows = []
    for name, mask in [('all_guides', pd.Series(False, index=counts.index)),
                       ('exclude_jackpot164', jackpot), ('exclude_QC173', outliers)]:
        retained = raw.loc[~mask, SAMPLES[:2]]
        recpm = retained.div(retained.sum(), axis=1) * 1e6
        relog = np.log2(recpm + 1)
        x, y = relog.to_numpy().T
        slope, intercept = np.polyfit(x, y, 1)
        r = float(np.corrcoef(x, y)[0, 1])
        # Independently check the Pearson calculation using centred vectors.
        xc, yc = x - x.mean(), y - y.mean()
        assert np.isclose(r, np.dot(xc, yc) / np.sqrt(np.dot(xc, xc) * np.dot(yc, yc)), atol=1e-12)
        rows.append(dict(analysis=name, excluded=int(mask.sum()), retained=len(retained),
                         pearson_recalculated_CPM=r,
                         pearson_original_CPM=float(log.loc[~mask].corr().iloc[0, 1]),
                         spearman=float(retained.corr(method='spearman').iloc[0, 1]),
                         slope=float(slope), intercept=float(intercept)))
        # The inset uses the 173-guide QC set, which adds five near-threshold
        # jackpot guides (D0/endpoint ratio 37-42) missed by the 164-guide rule.
        if name == 'exclude_QC173':
            assert len(retained) == 77268
            data = counts.loc[~mask, ['sgRNA', 'Gene']].copy()
            data[['D0_CTRL_log2_CPM_plus1', 'D0_LOY_log2_CPM_plus1']] = relog.to_numpy()
            data.to_csv(OUT / 'retained_guide_coordinates.csv.gz', index=False)
            excluded = counts.loc[mask, ['sgRNA', 'Gene']].copy()
            excluded['in_jackpot164'] = jackpot[mask].to_numpy()
            excluded['D0_mean_CPM'] = d0_mean[mask].to_numpy()
            excluded['D0_over_endpoint_CPM'] = (d0_mean / end_mean.replace(0, np.nan))[mask].to_numpy()
            excluded.to_csv(OUT / 'excluded_qc173_guides.csv', index=False)
    assert round(rows[0]['pearson_recalculated_CPM'], 3) == .957
    assert all(round(row['pearson_recalculated_CPM'], 3) == .927 for row in rows[1:])
    pd.DataFrame(rows).to_csv(OUT / 'correlation_summary.csv', index=False)
    manifest = {
        'method': 'Pearson correlation of log2(CPM + 1); CPM denominators recalculated after removal, separately for each D0 library. No minimum-read filter; retained zero-count guides included.',
        'jackpot_rule': 'Mean D0 CPM >200 and mean D0 / mean endpoint CPM >50, reproducing the existing 164-guide diagnostic set.',
        'qc173_rule': 'Guides with |z| >3 for D0 log2(CPM + 1) in either D0 library; contains all 164 jackpot guides and matches ht1080_figure1_qc/d0_outlier_guides.csv.',
        'inset_exclusion': 'exclude_QC173 (author decision 22 October 2026): also removes five jackpot guides (KEAP1, RET, CEBPA, KMT2B, SMARCB1; mean D0 CPM 2,180-12,665, D0/endpoint ratio 37-42) that the ratio >50 cutoff missed.',
        'scope': 'Descriptive D0 QC inset only; no screen candidate scores or existing exclusions changed.',
        'results': rows,
        'sources': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [COUNTS, FROZEN, Path(__file__)]},
    }
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
