"""Freeze DepMap positions of the parental HT1080 and A549 models for gametologs and screen hits.

Descriptive only: DepMap Public 25Q3 Chronos values for the parental lines are placed in the
distribution of all models and of the Chapter 3 Y-retaining / Y-depleted male cohorts. The
HT1080 per-arm guide changes reuse the Fig. 1c essential-gene QC convention.
"""
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
from scipy.stats import spearmanr

FIGS = THESIS / 'figures/chapter4'
SD = FIGS / 'source_data'
OUT = SD / 'depmap_model_positions'
DEPMAP = Path('A:/DEPMAP/DepMap')
EFFECT = DEPMAP / 'CRISPRGeneEffect.csv'
COHORT = SD / 'xy_paralog_ranked_partners_y_status/model_cohort_and_exclusions.csv'
HT_SCORES = SD / 'ht1080_candidate_discovery/all_gene_scores.csv'
A5_SCORES = SD / 'a549_candidate_discovery/all_gene_scores.csv'
SAGA = SD / 'ht1080_candidate_discovery/saga_reference_scores.csv'
COUNTS = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2/count/ht1080_screen2_R1.count.txt'
QC173 = SD / 'ht1080_figure1_qc/d0_outlier_guides.csv'
HT, A5 = 'ACH-000054', 'ACH-000681'
GAMETOLOGS = ['DDX3X', 'EIF1AX', 'ZFX', 'RPS4X']
ALIASES = {'ICK': 'CILK1'}  # screen library symbol -> DepMap 25Q3 symbol
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2', 'POST_LOY_REP1', 'POST_LOY_REP2']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sources = [EFFECT, DEPMAP / 'README.txt', COHORT, HT_SCORES, A5_SCORES, SAGA, COUNTS, QC173]
    before = {str(p): sha(p) for p in sources}

    ge = pd.read_csv(EFFECT, index_col=0)
    ge.columns = ge.columns.str.split(' ').str[0]
    assert {HT, A5} <= set(ge.index) and len(ge) == 1186
    coh = pd.read_csv(COHORT)
    coh = coh[coh.included.eq(True)]
    group = coh.set_index('ModelID').group
    assert (group == 'XY').sum() == 324 and (group == 'LOY').sum() == 229
    assert HT in group.index and group[HT] == 'XY' and A5 not in group.index

    ht = pd.read_csv(HT_SCORES).set_index('GENE')
    a5 = pd.read_csv(A5_SCORES).set_index('GENE')
    ht_hits = ht[(ht.directional_fdr < .05) & (ht.normZ < 0)].sort_values('normZ').index.tolist()
    a5_hits = a5[(a5.screen_fdr < .05) & (a5.direction == 'depleted')].sort_values('normZ').index.tolist()
    assert len(ht_hits) == 28 and len(a5_hits) == 19
    saga = set(pd.read_csv(SAGA).gene)

    panels = ([('a', g) for g in GAMETOLOGS] +
              [('b', g) for g in ht_hits if g not in GAMETOLOGS] +
              [('c', g) for g in a5_hits if g not in GAMETOLOGS])
    rows, missing = [], []
    for panel, g in panels:
        sym = ALIASES.get(g, g)
        if sym not in ge.columns:
            missing.append(dict(panel=panel, gene=g, reason='absent from DepMap 25Q3 gene-effect matrix'))
            continue
        v = ge[sym].dropna()
        xy, loy = v[v.index.isin(group[group == 'XY'].index)], v[v.index.isin(group[group == 'LOY'].index)]
        rows.append(dict(
            panel=panel, gene=g, depmap_symbol=sym, saga_reference=g in saga,
            HT1080_hit=g in ht_hits, A549_hit=g in a5_hits,
            n_models=len(v), all_median=v.median(), n_XY=len(xy), XY_q1=xy.quantile(.25), XY_median=xy.median(),
            XY_q3=xy.quantile(.75), n_LOY=len(loy), LOY_q1=loy.quantile(.25), LOY_median=loy.median(),
            LOY_q3=loy.quantile(.75),
            HT1080_chronos=v[HT], A549_chronos=v[A5],
            # Percent of models with a more negative (more dependent) effect than the parental line.
            HT1080_pct_all=(v < v[HT]).mean() * 100, A549_pct_all=(v < v[A5]).mean() * 100,
            HT1080_pct_XY=(xy < v[HT]).mean() * 100, A549_pct_XY=(xy < v[A5]).mean() * 100,
            HT1080_normZ=ht.normZ.get(g, np.nan), HT1080_fdr=ht.directional_fdr.get(g, np.nan),
            A549_normZ=a5.normZ.get(g, np.nan), A549_fdr=a5.screen_fdr.get(g, np.nan)))
    pos = pd.DataFrame(rows)
    assert [m['gene'] for m in missing] == ['NARR']

    # All-model values for plotted genes, with cohort labels.
    long = (ge[pos.depmap_symbol].rename(columns=dict(zip(pos.depmap_symbol, pos.gene)))
            .rename_axis('ModelID').reset_index().melt('ModelID', var_name='gene', value_name='chronos').dropna())
    long['cohort'] = long.ModelID.map(group).fillna('other_or_excluded')

    # HT1080 screen, per arm: 173-guide exclusion, CPM recalculated, >=50 raw D0 reads within the arm,
    # guide log2[(mean endpoint CPM + 1)/(D0 CPM + 1)] averaged per gene (Fig. 1c QC convention).
    c = pd.read_csv(COUNTS, sep='\t')
    raw = c[SAMPLES].astype(float)
    cpm = raw.div(raw.sum()) * 1e6
    lg = np.log2(cpm[SAMPLES[:2]] + 1)
    outlier = ((lg - lg.mean()) / lg.std(ddof=1)).abs().gt(3).any(axis=1)
    assert set(c.loc[outlier, 'sgRNA']) == set(pd.read_csv(QC173).sgRNA)
    keep = ~outlier
    r = raw[keep].div(raw[keep].sum()) * 1e6
    gene = c.loc[keep, 'Gene']
    arm_rows, ntc_rows = {}, []
    for arm in ['CTRL', 'LOY']:
        ok = raw.loc[keep, f'D0_{arm}'] >= 50
        lfc = np.log2((r[[f'POST_{arm}_REP1', f'POST_{arm}_REP2']].mean(axis=1) + 1) / (r[f'D0_{arm}'] + 1))
        arm_rows[f'{arm}_lfc'] = lfc[ok].groupby(gene[ok]).mean()
        arm_rows[f'{arm}_n_guides'] = lfc[ok].groupby(gene[ok]).size()
        ntc = lfc[ok & gene.eq('Non-Targeting-Control')]
        ntc_rows.append(dict(arm=arm, n_guides=len(ntc), q05=ntc.quantile(.05), median=ntc.median(),
                             q95=ntc.quantile(.95)))
    arm = pd.DataFrame(arm_rows)
    d = pos[pos.HT1080_hit][['gene', 'saga_reference', 'HT1080_chronos', 'HT1080_normZ']].copy()
    d = d.join(arm, on='gene')
    d['plotted'] = d[['CTRL_lfc', 'LOY_lfc']].notna().all(axis=1)
    dp = d[d.plotted]
    stats = {}
    for a in ['CTRL', 'LOY']:
        rho, p = spearmanr(dp.HT1080_chronos, dp[f'{a}_lfc'])
        stats[a] = dict(spearman_rho=float(rho), p_two_sided=float(p), n_genes=int(len(dp)))

    pos.to_csv(OUT / 'gene_positions.csv', index=False)
    pd.DataFrame(missing).to_csv(OUT / 'genes_not_in_depmap.csv', index=False)
    long.to_csv(OUT / 'plotted_model_values.csv.gz', index=False)
    d.to_csv(OUT / 'ht1080_hit_arm_changes.csv', index=False)
    pd.DataFrame(ntc_rows).to_csv(OUT / 'ht1080_ntc_arm_reference.csv', index=False)
    assert all(sha(p) == h for p, h in before.items())
    manifest = dict(
        purpose='Descriptive placement of parental HT1080 and A549 in DepMap Public 25Q3 Chronos distributions.',
        cohorts='Chapter 3 primary Y-status calls (DepMap 26Q1 eligibility) matched to 25Q3 by ModelID: 324 XY, 229 LOY '
                'male models. A549 has an unresolved primary call and is not in either cohort. Grey distribution: all 1,186 models.',
        percentile='Percent of models with a more negative Chronos value than the parental line.',
        hits='Depleted at drugZ FDR <0.05 (HT1080 directional FDR; A549 screen FDR). DDX3X (HT1080 hit) and EIF1AX '
             '(A549 hit) are shown once, in panel a. ICK is matched to its current symbol CILK1; NARR is absent from 25Q3.',
        arm_changes='173-guide exclusion, CPM recalculated, >=50 raw D0 reads within the arm, guide '
                    'log2[(mean endpoint CPM + 1)/(D0 CPM + 1)] averaged per gene; genes lacking an eligible guide in '
                    'either arm are not plotted.',
        spearman=stats, missing=missing, n_genes={k: int((pos.panel == k).sum()) for k in 'abc'},
        sources=before)
    (OUT / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({k: manifest[k] for k in ['n_genes', 'spearman', 'missing']}, indent=2))
    print('panel d not plotted:', d.loc[~d.plotted, 'gene'].tolist())


if __name__ == '__main__':
    main()
