"""Extend the saved partner display to author-selected KIF1A, preserving old results."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from prepare_ht1080_depmap_codependency import cross_cor,sha,MIN_N,DEPMAP,DATA as PREVIOUS

DATA=PREVIOUS.parent/'ht1080_combined_dependency_context'

def main():
    source=DEPMAP/'CRISPRGeneEffect.csv'
    old=json.loads((PREVIOUS/'analysis_manifest.json').read_text())
    source_hash=sha(source)
    assert source_hash==old['sources'][str(source)]['sha256']
    ge=pd.read_csv(source,index_col=0)
    ge.columns=ge.columns.str.replace(r' \(\d+\)$','',regex=True)
    assert ge.columns.is_unique and ge.index.is_unique
    assert ge.index.tolist()==pd.read_csv(PREVIOUS/'model_cohort.csv').ModelID.tolist()
    seeds=['KIF1A','ATXN7L3']
    r,n=cross_cor(ge[seeds].to_numpy(),ge.to_numpy())
    cached=np.load(PREVIOUS/'genomewide_correlations.npz')
    assert cached['genes'].tolist()==ge.columns.tolist()
    old_r=cached['r'][cached['seeds'].tolist().index('ATXN7L3')]
    error=float(np.nanmax(np.abs(r[1]-old_r)))
    assert error<1e-9
    annot=pd.read_csv(PREVIOUS/'seed_annotations.csv').set_index('gene')
    tab=pd.DataFrame(dict(seed='KIF1A',partner=ge.columns,pearson_r=r[0],n_models=n[0]))
    tab=tab.loc[tab.partner.ne('KIF1A')&tab.pearson_r.notna()].sort_values(['pearson_r','partner'],ascending=[False,True]).reset_index(drop=True)
    tab['rank_positive']=np.arange(1,len(tab)+1)
    tab['partner_screen_direction']=tab.partner.map(annot.screen_direction).fillna('Not nominated')
    tab['partner_saga_module']=tab.partner.map(annot.saga_module).fillna('')
    top=tab.loc[tab.pearson_r.gt(0)].head(20).assign(selection='Top positive')
    errors=[]
    for q in top.itertuples():
        pair=ge[['KIF1A',q.partner]].dropna()
        assert len(pair)==q.n_models
        errors.append(abs(float(pair.corr().iloc[0,1])-q.pearson_r))
    assert max(errors)<1e-9
    top.to_csv(DATA/'additional_top_partners.csv',index=False)
    tab.to_csv(DATA/'kif1a_genomewide_partner_rankings.csv.gz',index=False)
    ge[['KIF1A']].to_csv(DATA/'kif1a_gene_effect.csv',index_label='ModelID')
    record=dict(seed='KIF1A',seed_n_models=int(ge.KIF1A.notna().sum()),release='DepMap Public 25Q3',models=len(ge),genomewide_genes=ge.shape[1],
        correlation='Unadjusted pairwise-complete Pearson',min_models=MIN_N,selection='Top ten positive partners displayed; self excluded; no r or significance cutoff',
        source=str(source),source_sha256=source_hash,reference_atxn7l3_max_error=error,top20_independent_check_max_error=max(errors),
        original_analysis_modified=False,primary_candidate_list_modified=False)
    (DATA/'additional_partner_analysis.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(top.head(10).to_string(index=False))
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
