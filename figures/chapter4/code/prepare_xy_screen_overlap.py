"""Combine frozen X-gametolog co-dependency profiles with HT1080 drugZ scores."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd

THESIS=Path(__file__).resolve().parents[3]
SRC=THESIS/'figures/chapter4/source_data'
DATA=SRC/'xy_screen_overlap'
LIB=THESIS.parent/'crispr_screen/ht1080_loy_screen_2/brunello/broadgpp-brunello-library-contents.txt'

def main():
    DATA.mkdir(exist_ok=True)
    frozen=np.load(SRC/'xy_paralog_codependency/genomewide_correlations.npz')
    annotation=pd.read_csv(SRC/'xy_paralog_codependency/gene_annotations.csv').set_index('gene')
    score=pd.read_csv(SRC/'ht1080_candidate_discovery/all_gene_scores.csv').set_index('GENE')
    hits=pd.read_csv(SRC/'ht1080_candidate_discovery/depleted_candidates.csv').sort_values(['rank_synth','GENE']).set_index('GENE')
    saga=pd.read_csv(SRC/'ht1080_candidate_discovery/saga_reference_scores.csv').set_index('gene')
    saga_genes=set(saga.index[saga.reference_20])
    seeds=annotation.index[annotation.chromosome.eq('X')&annotation.eligible_seed].tolist()
    assert len(seeds)==16 and len(hits)==48
    identifiers=pd.read_csv(SRC/'ht1080_depmap_codependency/depmap_gene_identifiers.csv').set_index('gene')
    library=pd.read_csv(LIB,sep='\t',usecols=['Target Gene ID','Target Gene Symbol']).drop_duplicates()
    library=library.loc[library['Target Gene Symbol'].isin(score.index)]
    lookup=library.groupby('Target Gene ID')['Target Gene Symbol'].agg(lambda s:sorted(set(s))).to_dict()
    mapping=[]
    for gene in frozen['genes']:
        entrez=int(identifiers.loc[gene,'entrez_id'])
        choices=lookup.get(entrez,[])
        symbol=gene if gene in score.index else choices[0] if len(choices)==1 else None
        method='Exact symbol' if gene in score.index else 'Shared Entrez ID' if symbol else 'No unique screen match'
        mapping.append(dict(depmap_gene=gene,entrez_id=entrez,screen_gene=symbol,mapping=method))
    mapping=pd.DataFrame(mapping).set_index('depmap_gene')
    assert mapping.loc['ARHGAP45','screen_gene']=='HMHA1'
    assert mapping.loc['ARHGAP45','entrez_id']==23526
    mapping.to_csv(DATA/'gene_identity_mapping.csv')
    hitset=set(hits.index)
    rankings={};overlap=[];members=[]
    for seed in seeds:
        i=frozen['seeds'].tolist().index(seed)
        tab=pd.DataFrame(dict(seed=seed,partner=frozen['genes'],pearson_r=frozen['r'][i],n_models=frozen['n'][i]))
        tab=tab.loc[tab.partner.ne(seed)&tab.pearson_r.gt(0)].sort_values(['pearson_r','partner'],ascending=[False,True]).copy()
        tab['partner_rank']=np.arange(1,len(tab)+1)
        tab['screen_gene']=tab.partner.map(mapping.screen_gene)
        tab['screen_available']=tab.screen_gene.notna()
        tab['depleted_hit']=tab.screen_gene.isin(hitset)
        tab['saga_reference']=tab.partner.isin(saga_genes)
        tab['normZ']=tab.screen_gene.map(score.normZ)
        tab['fdr_synth']=tab.screen_gene.map(score.fdr_synth)
        rankings[seed]=tab
        members.append(tab.head(100))
        for k in [10,25,50,100]:
            selected=tab.head(k)
            overlap.append(dict(seed=seed,top_k=k,depleted_hits=int(selected.depleted_hit.sum()),
                                screen_measured=int(selected.screen_available.sum()),
                                hit_genes=';'.join(selected.loc[selected.depleted_hit,'screen_gene'])))
    counts=pd.DataFrame(overlap).pivot(index='seed',columns='top_k',values='depleted_hits')
    order=sorted(seeds,key=lambda g:(-counts.loc[g,25],-counts.loc[g,10],-counts.loc[g,50],-counts.loc[g,100],g))
    counts=counts.loc[order]
    assert counts.loc['KDM5C'].tolist()==[3,5,5,5]
    counts.to_csv(DATA/'panel_a_overlap_counts.csv')
    pd.DataFrame(overlap).to_csv(DATA/'panel_a_overlap_details.csv',index=False)
    pd.concat(members,ignore_index=True).to_csv(DATA/'top100_partner_memberships.csv',index=False)

    r=pd.DataFrame(frozen['r'],index=frozen['seeds'],columns=frozen['genes']).reindex(index=order,columns=hits.index).T
    n=pd.DataFrame(frozen['n'],index=frozen['seeds'],columns=frozen['genes']).reindex(index=order,columns=hits.index).T
    selected=pd.DataFrame(False,index=hits.index,columns=order)
    for seed in order:
        if seed in r.index:r.loc[seed,seed]=np.nan
        for g in rankings[seed].head(25).screen_gene:
            if g in selected.index:selected.loc[g,seed]=True
    assert selected.to_numpy().sum()==counts[25].sum()
    assert r.loc['ZCCHC16'].isna().all()
    assert int(r.isna().to_numpy().sum())==17 # 16 absent + DDX3X self.
    assert n.where(r.notna()).min().min()==1046
    r.to_csv(DATA/'panel_b_correlations.csv',index_label='screen_gene')
    n.to_csv(DATA/'panel_b_model_counts.csv',index_label='screen_gene')
    selected.to_csv(DATA/'panel_b_top25_membership.csv',index_label='screen_gene')
    hits=hits.assign(saga_reference=hits.index.isin(saga_genes))
    hits.to_csv(DATA/'panel_b_screen_scores.csv')
    detail=rankings['KDM5C'].head(25).copy()
    assert detail.screen_available.all()
    assert detail.loc[detail.depleted_hit,'partner_rank'].tolist()==[4,8,9,13,23]
    detail.to_csv(DATA/'panel_c_kdm5c_top25.csv',index=False)
    sources=[SRC/'xy_paralog_codependency/genomewide_correlations.npz',SRC/'xy_paralog_codependency/gene_annotations.csv',
             SRC/'ht1080_candidate_discovery/all_gene_scores.csv',SRC/'ht1080_candidate_discovery/depleted_candidates.csv',
             SRC/'ht1080_candidate_discovery/saga_reference_scores.csv',SRC/'ht1080_depmap_codependency/depmap_gene_identifiers.csv',LIB]
    manifest=dict(release='DepMap Public 25Q3',primary_depleted_candidates=48,primary_fdr_threshold=.25,
                  seed_order=order,hit_order=hits.index.tolist(),analysis_refit=False,
                  correlation='Unadjusted pairwise-complete Pearson; minimum 200 shared models; all sexes and lineages.',
                  ranking='Positive partners, self excluded, descending r then symbol; rankings are not restricted to screen-measured genes.',
                  identity_mapping='Exact symbol, otherwise unique shared Entrez ID from the archived Brunello library; original screen symbols and scores retained.',
                  overlap='Descriptive counts, not enrichment tests; no new P values or FDRs. Denominator is top k partners; measured-screen coverage saved separately.',
                  missing='ZCCHC16 has no DepMap profile; DDX3X self-comparison masked. Neither is zero.',
                  panel_b_valid_pair_models=[1046,1186],panel_c_pair_models=sorted(detail.n_models.unique().tolist()),
                  kdm5c_screen=score.loc['KDM5C',['normZ','fdr_synth']].to_dict(),
                  kdm5c_hit_ranks=detail.loc[detail.depleted_hit,['partner','partner_rank','pearson_r','normZ','fdr_synth']].to_dict('records'),
                  atxn7l3_kdm5c_rank=int(rankings['KDM5C'].set_index('partner').loc['ATXN7L3','partner_rank']),
                  source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (DATA/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(counts.to_string())
    print('KDM5C detail: all 25 screen scores matched; ARHGAP45 maps to HMHA1 through Entrez 23526.')

if __name__=='__main__':main()
