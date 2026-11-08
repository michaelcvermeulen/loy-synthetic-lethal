"""Prepare descriptive DepMap partner profiles for the archived gametolog set."""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
from prepare_ht1080_depmap_codependency import cross_cor,sha,MIN_N

THESIS=Path(__file__).resolve().parents[3]
DATA=THESIS/'figures/chapter4/source_data/xy_paralog_codependency'
REF=THESIS.parent/'pubs/summaries/_assets/2025_Decasien_XY_divergence/gametologs_code/gametologs-main/data/gametologs_in_genome.csv'
DEPMAP=Path('A:/DEPMAP/DepMap')
PREVIOUS=THESIS/'figures/chapter4/source_data/ht1080_depmap_codependency'


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    reference=pd.read_csv(REF,sep='\t')
    reference.to_csv(DATA/'reference_gametolog_pairs.csv',index=False)
    print('Reading the same DepMap 25Q3 gene-effect matrix...',flush=True)
    ge=pd.read_csv(DEPMAP/'CRISPRGeneEffect.csv',index_col=0)
    full_columns=ge.columns.tolist()
    ge.columns=ge.columns.str.replace(r' \(\d+\)$','',regex=True)
    assert ge.columns.is_unique and ge.index.is_unique
    cohort=pd.read_csv(PREVIOUS/'model_cohort.csv').set_index('ModelID')
    assert ge.index.tolist()==cohort.index.tolist()
    models=pd.read_csv(DEPMAP/'Model.csv').set_index('ModelID').reindex(ge.index)
    models[['CellLineName','OncotreeLineage','OncotreePrimaryDisease','Sex']].to_csv(DATA/'model_cohort.csv',index_label='ModelID')
    known={'DDX3X','EIF1AX','RPS4X','ZFX'}
    rows=[]
    for g,r in reference.drop_duplicates('common_name').set_index('common_name').iterrows():
        pairs=reference.loc[reference.common_name.eq(g),'pair'].tolist()
        counterparts=reference.loc[reference.pair.isin(pairs)&reference.common_name.ne(g),'common_name'].drop_duplicates().tolist()
        n=int(ge[g].notna().sum()) if g in ge.columns else 0
        variable=bool(ge[g].std()>1e-9) if g in ge.columns else False
        eligible=n>=MIN_N and variable
        rows.append(dict(gene=g,chromosome=r.Gametolog,counterparts=';'.join(counterparts),
                         reference_group=r.Group,core_four=g in known,depmap_present=g in ge.columns,
                         n_models=n,eligible_seed=eligible,
                         exclusion='' if eligible else ('No DepMap profile' if g not in ge.columns else f'Only {n} models; minimum {MIN_N}'),
                         depmap_column=full_columns[ge.columns.get_loc(g)] if g in ge.columns else ''))
    annotations=pd.DataFrame(rows).set_index('gene')
    annotations.to_csv(DATA/'gene_annotations.csv')
    seeds=annotations.index[annotations.eligible_seed].tolist()
    selected=ge[seeds]
    selected.to_csv(DATA/'seed_gene_effect.csv',index_label='ModelID')
    print(f'Computing {len(seeds)} seed profiles against {ge.shape[1]} genes...',flush=True)
    r,n=cross_cor(selected.to_numpy(),ge.to_numpy())
    np.savez_compressed(DATA/'genomewide_correlations.npz',seeds=np.array(seeds),genes=ge.columns.to_numpy(dtype=str),r=r,n=n)
    top=[]
    for i,g in enumerate(seeds):
        tab=pd.DataFrame(dict(seed=g,partner=ge.columns,pearson_r=r[i],n_models=n[i]))
        tab=tab.loc[tab.partner.ne(g)&tab.pearson_r.gt(0)].sort_values(['pearson_r','partner'],ascending=[False,True]).head(10).copy()
        tab['rank_positive']=np.arange(1,len(tab)+1)
        assert len(tab)==10
        top.append(tab)
    pd.concat(top,ignore_index=True).to_csv(DATA/'top_10_partners.csv',index=False)
    # Reproduce the existing DDX3X partner profile before making the new display.
    previous=np.load(PREVIOUS/'genomewide_correlations.npz')
    assert previous['genes'].tolist()==ge.columns.tolist()
    previous_seed=previous['seeds'].tolist().index('DDX3X')
    error=float(np.nanmax(np.abs(previous['r'][previous_seed]-r[seeds.index('DDX3X')])))
    assert error<1e-9,error
    checks=[]
    rng=np.random.default_rng(20261019)
    for i,j in zip(rng.integers(len(seeds),size=40),rng.integers(ge.shape[1],size=40)):
        pair=ge[[seeds[i],ge.columns[j]]].dropna()
        assert n[i,j]==len(pair)
        if len(pair)>=MIN_N and pair.std().min()>1e-9:
            checks.append(abs(float(np.corrcoef(pair.to_numpy().T)[0,1])-r[i,j]))
    assert max(checks)<1e-9
    source_hash=sha(DEPMAP/'CRISPRGeneEffect.csv')
    old=json.loads((PREVIOUS/'analysis_manifest.json').read_text())
    assert source_hash==old['sources'][str(DEPMAP/'CRISPRGeneEffect.csv')]['sha256']
    record=dict(release='DepMap Public 25Q3',models=len(ge),genomewide_genes=ge.shape[1],
                reference_pairs=int(reference.pair.nunique()),reference_unique_genes=len(annotations),
                reference_X=int(annotations.chromosome.eq('X').sum()),reference_Y=int(annotations.chromosome.eq('Y').sum()),
                eligible_seeds=len(seeds),eligible_X=int((annotations.chromosome.eq('X')&annotations.eligible_seed).sum()),
                eligible_Y=int((annotations.chromosome.eq('Y')&annotations.eligible_seed).sum()),
                exclusions=annotations.loc[~annotations.eligible_seed,['chromosome','n_models','exclusion']].reset_index().to_dict('records'),
                min_models=MIN_N,correlation='Unadjusted pairwise-complete Pearson of corrected Chronos gene effects',
                cohort='All models in local gene-effect matrix, all recorded sexes and lineages; no LOY restriction',
                partner_selection='Top ten positive correlations per seed, excluding self-correlations; no r or significance cutoff',
                scope='All unique genes in the archived gametolog reference table, not an exhaustive homology search',
                validation=dict(existing_ddx3x_profile_max_error=error,numpy_pair_checks=len(checks),numpy_max_error=max(checks)),
                sources={str(p):dict(sha256=sha(p),bytes=p.stat().st_size) for p in [REF,DEPMAP/'CRISPRGeneEffect.csv',DEPMAP/'Model.csv']})
    (DATA/'analysis_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k!='sources'},indent=2))


if __name__=='__main__':main()
