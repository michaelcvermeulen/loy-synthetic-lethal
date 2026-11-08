"""Score the Schukken et al. (2026) paired aneuploidy screens with IsogenicZ and matched endpoint drugZ.

Exploratory analysis requested 22 October 2026. Each screen has one population per arm, with one D0 and one
D10 library, so every contrast is one aneuploid versus one control endpoint (IsogenicZ: each endpoint normalised
to its own arm's D0). Pairings follow the authors' own batches in Code/MegaA_Analysis_6_KMS.R and
Code/WGL_all_Analysis_3_KMS.R (commit 4db9338). The procedure and pinned IsogenicZ copy are those of
research/external_isogenicz_20261022. Negative normZ = depleted in the aneuploid line.
Run with Anaconda Python: python -I run_schukken_screens.py
"""
from pathlib import Path
import hashlib, importlib.util, json, logging, re, shutil, warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
TH = HERE.parents[1]
DATA = TH.parent / 'crispr_screen/Schukken_2026_aneuploid_screens/repo/Data'
OUT = TH / 'figures/chapter4/source_data/schukken_2026_rescored'
PINNED = TH / 'research/external_isogenicz_20261022/upstream'
UPSTREAM_SHA = '939667ec51e7c602e4302678a4870f85bb847c756367622645105486a657beb9'
FOCUS = ['UBE2H', 'SREBF1', 'SREBF2', 'TAF6L', 'TADA2B', 'SUPT20H', 'ATXN7L3', 'TADA1', 'TAF5L', 'MED12', 'CREBBP',
         'DDX3X', 'EIF1AX', 'PELO', 'UBE2G2', 'ACSL4', 'TIMM8A']
B, F = 'Brunello genome-wide', 'druggable-genome focused'
# (contrast, library, batch, aneuploid file, control file)
DESIGN = [
    # Restored-disomy designs: the parental line carries the gain and is the MORE aneuploid arm, so it is
    # placed in the "aneuploid" position and the disomic derivative is the reference.
    ('A2058_1q_gain_vs_Di', B, 'A2058', 'count-A2058_WT.csv', 'count-A2058_1qDi.csv'),
    ('A2058_7q_gain_vs_Di', B, 'A2058', 'count-A2058_WT.csv', 'count-A2058_7qDi.csv'),
    ('A2780_1q_gain_vs_Di', B, 'A2780 WGL', 'count-A2780-Cont1-d0.d10.csv', 'count-A2780-1qDi-d0.d10.csv'),
    ('AGS_1q_gain_vs_Di', B, 'AGS', 'count-AGS_Control.csv', 'count-AGS_1qDi.csv'),
    ('MCF10A_1q_gain_vs_Di', B, 'MCF10A', 'count-MCF10A_WT.csv', 'count-MCF10A_1qDi.csv'),
    ('Vaco_Ts13', B, 'Vaco432', 'count-Vaco-Ts13-i.f.csv', 'count-Vaco-Cont-i.f.csv'),
    ('HCT116_Ts13_WGL', B, 'HCT116 WGL', 'Ryan-Count_HCT116_Ts13_i.f.txt', 'Ryan-Count_HCT116_Cont1_i.f.head.txt'),
    ('HCT116_Ts8', F, 'MegaA 1', 'count-HCT116-Ts8-i.f.csv', 'count-HCT116-Cont1-i.f.csv'),
    ('DLD1_Ts13', F, 'MegaA 1', 'count-DLD1-Ts13-i.f.csv', 'count-DLD1-Cont1-i.f.csv'),
    ('DLD1_Ts8', F, 'MegaA 1', 'count-DLD1-Ts8-i.f.csv', 'count-DLD1-Cont1-i.f.csv'),
    ('HCT116_Ts3', F, 'MegaA 2', 'count-HCT116-Ts3-i.f.csv', 'count-HCT116-Cont2-i.f.csv'),
    ('HCT116_Ts5', F, 'MegaA 2', 'count-HCT116-Ts5-i.f.csv', 'count-HCT116-Cont2-i.f.csv'),
    ('DLD1_Ts4', F, 'MegaA 2', 'count-DLD1-Ts4-i.f.csv', 'count-DLD1-Cont2_P10.13-i.f-wo-cdk12.csv'),
    ('DLD1_Ts2.18', F, 'MegaA 2', 'count-DLD1-Ts2.18-i.f.csv', 'count-DLD1-Cont2_P10.13-i.f-wo-cdk12.csv'),
    ('DLD1_Ts10.21', F, 'MegaA 3', 'count-DLD1-Ts10.21-i.f.csv', 'count-DLD1-Cont3-i.f.csv'),
    ('DLD1_Ts5.15', F, 'MegaA 3', 'count-DLD1-Ts5.15-i.f.csv', 'count-DLD1-Cont3-i.f.csv'),
    ('A2780_Ts2', F, 'MegaA 4', 'count-A2780-Ts2-i.f.csv', 'count-A2780-Cont1-i.f.csv'),
    ('A2780_Ts10', F, 'MegaA 4', 'count-A2780-Ts10-i.f.csv', 'count-A2780-Cont1-i.f.csv'),
    ('A2780_Ts8', F, 'MegaA 5', 'count-A2780-Ts8.csv', 'count-A2780-Cont2.csv'),
    ('A2780_Ts18', F, 'MegaA 5', 'count-A2780-Ts18.csv', 'count-A2780-Cont2.csv'),
    ('SNU1_C12', F, 'MegaA 6', 'count-Snu-1_clone_12.csv', 'count-Snu-1_WT_C13.csv'),
    ('SNU1_C24', F, 'MegaA 6', 'count-Snu1C24D0_v_Snu1C23D10.csv', 'count-Snu-1_WT_C13.csv'),
    ('SNU1_C111', F, 'MegaA 6', 'count-Snu-1_clone_111.csv', 'count-Snu-1_WT_C13.csv'),
    ('SNU1_C114_Ts9.19', F, 'MegaA 6', 'count-Snu-1_clone_114.csv', 'count-Snu-1_WT_C13.csv'),
    ('DLD1_Ts2.18_rep', F, 'MegaA 7', 'count-DLD1__Ts2.18.csv', 'count-DLD1_Ctrl_3.csv'),
    ('DLD1_Ts12.17', F, 'MegaA 7', 'count-DLD1_Ts12.17.csv', 'count-DLD1_Ctrl_3.csv'),
]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_program():
    (HERE / 'upstream').mkdir(exist_ok=True)
    for name in ['isogenicz.py', 'README.md', 'LICENSE', 'commit.json']:
        dst = HERE / 'upstream' / name
        if not dst.exists():
            shutil.copyfile(PINNED / name, dst)
    assert sha(HERE / 'upstream/isogenicz.py') == UPSTREAM_SHA
    spec = importlib.util.spec_from_file_location('pinned_isogenicz', HERE / 'upstream/isogenicz.py')
    program = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(program)
    pd.options.mode.copy_on_write = False
    program.log_.setLevel(logging.WARNING)
    return program


def read_counts(name):
    """Return sgRNA-indexed frame with Gene, D0, D10 (column order and naming vary between files)."""
    p = DATA / name
    if name.startswith('Ryan-Count') and not name.endswith('.head.txt'):
        d = pd.read_csv(p, sep='\t', header=None, names=['sgRNA', 'Gene', 'initial', 'final'], dtype={'sgRNA': str})
    else:
        d = pd.read_csv(p, sep='\t' if p.suffix == '.txt' else ',', encoding='utf-8-sig', dtype={'sgRNA': str})
    d = d.loc[:, ~d.columns.str.startswith('Unnamed')]
    cols = list(d.columns)
    d0 = [c for c in cols[2:] if re.search(r'initial|D=?0$|D0$', c)]
    d10 = [c for c in cols[2:] if re.search(r'final$|D=?10$', c)]
    assert len(d0) == 1 and len(d10) == 1, (name, cols)
    out = d[['sgRNA', 'Gene', d0[0], d10[0]]].rename(columns={d0[0]: 'D0', d10[0]: 'D10'})
    # The HCT116 WGL control file leaves non-targeting guide IDs blank; they are removed later anyway.
    blank = out.sgRNA.isna()
    assert is_control(out.loc[blank, 'Gene'].astype(str)).all(), name
    out = out[~blank].copy()
    # Some files carry Excel date-converted symbols (e.g. DEC1 -> 45261); recover them from GENE-SEQ guide IDs.
    out['Gene'] = out.Gene.astype(str)
    num = out.Gene.str.fullmatch(r'\d+(\.0)?') & out.sgRNA.notna()
    dash = num & out.sgRNA.str.contains('-', na=False)        # Brunello GENE-SEQ IDs
    under = num & ~dash & out.sgRNA.str.contains('_', na=False)  # focused-library GENE_domain_exon IDs
    out.loc[dash, 'Gene'] = out.loc[dash, 'sgRNA'].str.rsplit('-', n=1).str[0]
    out.loc[under, 'Gene'] = out.loc[under, 'sgRNA'].str.split('_', n=1).str[0]
    # Coordinate-style IDs repeat for readthrough genes (e.g. APITD1 / APITD1-CORT); number repeats in file order.
    rep = out.groupby('sgRNA').cumcount()
    out.loc[rep > 0, 'sgRNA'] = out.loc[rep > 0, 'sgRNA'] + '#' + rep[rep > 0].astype(str)
    assert out.sgRNA.is_unique, name
    return out.set_index('sgRNA'), (d0[0], d10[0])


def is_control(g):
    # Exact control labels only, so genes such as NEGR1 are not caught.
    return g.str.contains(r'(?i)non.?targeting', regex=True) | g.str.fullmatch(r'(?i)neg_?\d*|control|ctrl|egfp|lacz|luciferase')


def run_one(program, row):
    cid, lib, batch, ane_f, ctl_f = row
    a, acols = read_counts(ane_f)
    c, ccols = read_counts(ctl_f)
    sh = a.index.intersection(c.index)
    ga, gc = a.loc[sh, 'Gene'], c.loc[sh, 'Gene']
    mismatch = ga != gc
    # Where labels still disagree, prefer the non-numeric (non-date-converted) label.
    gene = ga.where(~(mismatch & ga.str.fullmatch(r'\d+(\.0)?')), gc)
    assert (gene == gc).mean() > .99 or (gene == ga).mean() > .99, (cid, int(mismatch.sum()))
    reads = pd.DataFrame({'GENE': gene, 'T0_REF': c.loc[sh, 'D0'], 'T0_LES': a.loc[sh, 'D0'],
                          'REF1': c.loc[sh, 'D10'], 'LES1': a.loc[sh, 'D10']})
    S = ['T0_REF', 'T0_LES', 'REF1', 'LES1']
    reads[S] = reads[S].fillna(0).round().astype(int)
    # A guide is a control if either file labels it as one (A2058 7qDi uses neg_N where WT uses Non-Targeting Control).
    ctrl = is_control(reads.GENE) | is_control(c.loc[sh, 'Gene'].astype(str)) | is_control(a.loc[sh, 'Gene'].astype(str))
    keep = ~ctrl & reads.GENE.ne('CDK12') & reads[S].gt(0).any(axis=1)
    reads = reads[keep].copy()
    iso_norm = program.normalize_readcounts(reads, ['T0_LES'], ['LES1'], ['T0_REF'], ['REF1'])[['REF1', 'LES1']]
    plain = reads[['REF1', 'LES1']].div(reads[['REF1', 'LES1']].sum()) * 1e7
    out = {}
    for method, norm in [('IsogenicZ', iso_norm), ('matched_endpoint_drugZ', plain)]:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fc = program.calculate_unpaired_foldchange(reads, norm, ['REF1'], ['LES1'], 5)
            fc = program.empirical_bayes(fc, 500, len(fc), 'mean_fc', 'eb_std', 'zscore')[0]
            agg = fc.groupby('GENE').zscore.agg(sumZ='sum', numObs=lambda x: np.count_nonzero(x))
            res = program.calculate_drugz_score_unpaired(agg.copy(), 1)
        assert np.isfinite(res.normZ).all(), (cid, method)
        res.index.name = 'Gene'
        out[method] = res
    gini = lambda x: float((2 * np.arange(1, len(x) + 1) - len(x) - 1) @ np.sort(x) / (len(x) * np.sum(x)))
    qc = dict(id=cid, library=lib, batch=batch, aneuploid_file=ane_f, control_file=ctl_f,
              aneuploid_cols='; '.join(acols), control_cols='; '.join(ccols), shared_guides=len(sh),
              gene_label_mismatches_after_repair=int((ga != gc).sum()),
              numeric_gene_labels_remaining=int(reads.GENE.str.fullmatch(r'\d+(\.0)?').sum()),
              scored_guides=int(len(reads)), genes=int(reads.GENE.nunique()),
              d0_spearman=float(spearmanr(reads.T0_REF, reads.T0_LES).statistic),
              d10_spearman=float(spearmanr(reads.REF1, reads.LES1).statistic),
              gini_d0_ctl=gini(reads.T0_REF), gini_d0_ane=gini(reads.T0_LES),
              median_reads_d0_ctl=float(reads.T0_REF.median()), median_reads_d0_ane=float(reads.T0_LES.median()))
    return out, qc


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    program = load_program()
    longs, summ = [], []
    for row in DESIGN:
        out, qc = run_one(program, row)
        iso, mat = out['IsogenicZ'], out['matched_endpoint_drugZ']
        sh = iso.index.intersection(mat.index)
        qc.update(rho_iso_vs_drugz=float(spearmanr(iso.loc[sh, 'normZ'], mat.loc[sh, 'normZ']).statistic),
                  iso_dep_fdr025=int((iso.fdr_synth < .25).sum()), drugz_dep_fdr025=int((mat.fdr_synth < .25).sum()),
                  both_dep_fdr025=len(set(iso.index[iso.fdr_synth < .25]) & set(mat.index[mat.fdr_synth < .25])),
                  drugz_enr_fdr025=int((mat.fdr_supp < .25).sum()))
        summ.append(qc)
        for method, res in out.items():
            res.to_csv(OUT / f'{row[0]}_{method}.csv')
            longs.append(res[['normZ', 'fdr_synth', 'fdr_supp']].assign(id=row[0], library=row[1], method=method).reset_index())
        print(row[0], {k: qc[k] for k in ['scored_guides', 'd0_spearman', 'rho_iso_vs_drugz', 'drugz_dep_fdr025', 'iso_dep_fdr025']}, flush=True)
    pd.DataFrame(summ).to_csv(OUT / 'contrast_summary.csv', index=False)
    L = pd.concat(longs, ignore_index=True)
    L.to_csv(OUT / 'schukken_scores_long.csv.gz', index=False)
    L[L.Gene.isin(FOCUS)].to_csv(OUT / 'focus_gene_scores.csv', index=False)
    json.dump(dict(repo_commit='4db9338d90c90862fff569c6d3f5c15c5c13000e', upstream_sha256=UPSTREAM_SHA, pseudocount=5,
                   eb_window=500, design='one aneuploid vs one control endpoint per contrast; IsogenicZ uses each arm\'s D0',
                   removed='non-targeting controls, CDK12 (contaminated MegaA screen) and guides with no reads in any library',
                   pandas=pd.__version__, inputs={r[0]: [sha(DATA / r[3]), sha(DATA / r[4])] for r in DESIGN},
                   script_sha256=sha(__file__)), open(OUT / 'manifest.json', 'w'), indent=2)


if __name__ == '__main__':
    main()
