"""DepMap 25Q3 co-dependencies for current drugZ hits and frozen SAGA members.

Primary: Pearson correlation, pairwise complete, at least 200 shared models.
Genome-wide partner rankings are descriptive; no co-dependency hit/FDR calls.
Sensitivity: correlations of per-gene OLS residuals after lineage + recorded sex.
"""
from pathlib import Path
import hashlib
import json
import sys
import re

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
import numpy as np
import pandas as pd

DATA = THESIS / 'figures/chapter4/source_data/ht1080_depmap_codependency'
CAND = DATA.parent / 'ht1080_candidate_discovery'
DEPMAP = Path('A:/DEPMAP/DepMap')
MIN_N = 200


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def cross_cor(x, y):
    """Pairwise-complete Pearson correlations and denominators (features in columns)."""
    mx=np.isfinite(x).astype(float); my=np.isfinite(y).astype(float)
    x=np.nan_to_num(x); y=np.nan_to_num(y)
    n=mx.T@my
    sx=x.T@my; sy=mx.T@y
    sxx=(x*x).T@my; syy=mx.T@(y*y)
    with np.errstate(invalid='ignore',divide='ignore'):
        num=x.T@y-sx*sy/n
        den=np.sqrt(np.maximum(sxx-sx*sx/n,0)*np.maximum(syy-sy*sy/n,0))
        r=num/den
    r[(n<MIN_N)|(den<=1e-12)]=np.nan
    return np.clip(r,-1,1),n.astype(int)


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    hit=pd.concat([pd.read_csv(CAND/'depleted_candidates.csv').assign(direction='Depleted'),
                   pd.read_csv(CAND/'enriched_candidates.csv').assign(direction='Enriched')],ignore_index=True)
    assert len(hit)==57 and hit.GENE.is_unique
    saga=pd.read_csv(CAND/'saga_reference_scores.csv')
    assert len(saga)==24 and saga.reference_20.sum()==20
    h=hit.set_index('GENE'); s=saga.set_index('gene')
    seeds=list(dict.fromkeys(hit.GENE.tolist()+saga.gene.tolist()))
    print('Reading local DepMap 25Q3 gene-effect matrix...',flush=True)
    ge=pd.read_csv(DEPMAP/'CRISPRGeneEffect.csv',index_col=0)
    gene_map=pd.DataFrame({'depmap_column':ge.columns})
    gene_map['gene']=gene_map.depmap_column.str.replace(r' \(\d+\)$','',regex=True)
    gene_map['entrez_id']=gene_map.depmap_column.str.extract(r'\((\d+)\)$').iloc[:,0]
    assert gene_map.gene.is_unique and ge.index.is_unique
    ge.columns=gene_map.gene
    assert all(re.fullmatch(r'ACH-\d+',v) for v in ge.index)
    model=pd.read_csv(DEPMAP/'Model.csv').set_index('ModelID').reindex(ge.index)
    assert model.OncotreeLineage.notna().all()
    model[['CellLineName','OncotreeLineage','OncotreePrimaryDisease','Sex']].to_csv(DATA/'model_cohort.csv',index_label='ModelID')
    rows=[]
    for g in seeds:
        rows.append({'gene':g,'screen_hit':g in h.index,
            'screen_direction':h.loc[g,'direction'] if g in h.index else 'Not nominated',
            'normZ':h.loc[g,'normZ'] if g in h.index else s.loc[g,'normZ'],
            'screen_fdr':h.loc[g,'directional_fdr'] if g in h.index else min(s.loc[g,'fdr_synth'],s.loc[g,'fdr_supp']),
            'saga_module':s.loc[g,'module'] if g in s.index else '',
            'saga_reference_20':bool(s.loc[g,'reference_20']) if g in s.index else False,
            'depmap_available':g in ge.columns,
            'depmap_n':int(ge[g].notna().sum()) if g in ge.columns else 0,
            'depmap_column':gene_map.set_index('gene').loc[g,'depmap_column'] if g in ge.columns else '',
            'mapping':'exact symbol' if g in ge.columns else 'absent; no substitution'})
    annot=pd.DataFrame(rows).set_index('gene')
    annot.to_csv(DATA/'seed_annotations.csv')
    gene_map.to_csv(DATA/'depmap_gene_identifiers.csv',index=False)
    present=[g for g in seeds if g in ge.columns]
    selected=ge[present]
    selected.to_csv(DATA/'seed_gene_effect.csv',index_label='ModelID')
    # Pairwise complete, not listwise deletion across the full genome.
    print(f'Correlating {len(present)} seeds against {ge.shape[1]} genes in {len(ge)} models...',flush=True)
    r,n=cross_cor(selected.to_numpy(),ge.to_numpy())
    np.savez_compressed(DATA/'genomewide_correlations.npz',seeds=np.array(present),genes=ge.columns.to_numpy(dtype=str),r=r,n=n)
    ix=ge.columns.get_indexer(present)
    raw=pd.DataFrame(r[:,ix],index=present,columns=present).reindex(index=seeds,columns=seeds)
    ns=pd.DataFrame(n[:,ix],index=present,columns=present).reindex(index=seeds,columns=seeds)
    raw.to_csv(DATA/'seed_pairwise_pearson.csv',index_label='gene')
    ns.to_csv(DATA/'seed_pairwise_n.csv',index_label='gene')
    # Descriptive nuisance-adjusted sensitivity; no claim of LOY specificity.
    cov=pd.get_dummies(model[['OncotreeLineage','Sex']].fillna('Unknown'),drop_first=True,dtype=float)
    design=np.column_stack([np.ones(len(model)),cov.to_numpy()])
    residual=selected.copy()
    ranks={}
    for g in present:
        ok=selected[g].notna().to_numpy()
        beta,_,rank,_=np.linalg.lstsq(design[ok],selected.loc[ok,g].to_numpy(),rcond=None)
        residual.loc[ok,g]=selected.loc[ok,g].to_numpy()-design[ok]@beta
        ranks[g]=int(rank)
    rr,rn=cross_cor(residual.to_numpy(),residual.to_numpy())
    adjusted=pd.DataFrame(rr,index=present,columns=present).reindex(index=seeds,columns=seeds)
    adjusted.to_csv(DATA/'seed_pairwise_lineage_sex_residual_pearson.csv',index_label='gene')
    residual.to_csv(DATA/'seed_lineage_sex_residuals.csv',index_label='ModelID')
    # Genome-wide full ranking and compact positive/negative selections.
    top=[]; full=[]; summary=[]
    for i,g in enumerate(present):
        tab=pd.DataFrame({'seed':g,'partner':ge.columns,'pearson_r':r[i],'n_models':n[i]})
        tab=tab.loc[tab.partner.ne(g)&tab.pearson_r.notna()].sort_values(['pearson_r','partner'],ascending=[False,True]).reset_index(drop=True)
        tab['rank_positive']=np.arange(1,len(tab)+1)
        tab['partner_screen_direction']=tab.partner.map(annot.screen_direction).fillna('Not nominated')
        tab['partner_saga_module']=tab.partner.map(annot.saga_module).fillna('')
        full.append(tab)
        pos=tab.loc[tab.pearson_r.gt(0)].head(20).assign(selection='Top positive')
        neg=tab.loc[tab.pearson_r.lt(0)].sort_values(['pearson_r','partner']).head(10).assign(selection='Top negative')
        top.extend([pos,neg])
        summary.append({'gene':g,'screen_direction':annot.loc[g,'screen_direction'],
            'saga_module':annot.loc[g,'saga_module'],'eligible_partners':len(tab),
            'positive_partners_r_ge_03':int(tab.pearson_r.ge(.3).sum()),
            'strongest_partner':tab.iloc[0].partner,'strongest_r':float(tab.iloc[0].pearson_r),
            'top_5':'; '.join(f'{a.partner} ({a.pearson_r:.3f})' for a in tab.head(5).itertuples())})
    pd.concat(full,ignore_index=True).to_csv(DATA/'genomewide_partner_rankings.csv.gz',index=False,float_format='%.10g')
    top=pd.concat(top,ignore_index=True)
    top.to_csv(DATA/'top_partners.csv',index=False)
    pd.DataFrame(summary).to_csv(DATA/'seed_partner_summary.csv',index=False)
    # Validate against numpy.corrcoef, with and without missing observations.
    rng=np.random.default_rng(20261015)
    errors=[]
    pairs=[('DDX3X','EIF1AX'),('TADA2B','TAF6L'),('SREBF1','SREBF2')]
    pairs += [(present[int(i)],ge.columns[int(j)]) for i,j in zip(rng.integers(len(present),size=120),rng.integers(ge.shape[1],size=120))]
    for g,p in pairs:
        complete=ge[[g,p]].dropna()
        i=present.index(g); j=ge.columns.get_loc(p)
        assert n[i,j]==len(complete)
        if len(complete)>=MIN_N and (complete.std()>1e-9).all():
            errors.append(abs(float(np.corrcoef(complete.to_numpy().T)[0,1])-r[i,j]))
    assert max(errors)<1e-9
    # Check surviving historical top-partner measurements independently of old hit calls.
    oldpath=THESIS.parent/'thesis/LOY_paralog/tables/codependency_partners_long.csv'
    old=pd.read_csv(oldpath)
    oldcheck=[]
    for row in old.itertuples():
        if row.hit in present and row.partner in ge.columns:
            oldcheck.append(abs(row.r-r[present.index(row.hit),ge.columns.get_loc(row.partner)]))
    five=['TAF6L','TADA2B','SUPT20H','TADA1','TAF5L']
    six=five+['ATXN7L3']
    sagapairs=[]
    for i,g in enumerate(six):
        for p in six[i+1:]:
            sagapairs.append({'gene1':g,'gene2':p,'pearson_r':raw.loc[g,p],'residual_r':adjusted.loc[g,p],'n_models':ns.loc[g,p]})
    pd.DataFrame(sagapairs).to_csv(DATA/'saga_hit_pair_summary.csv',index=False)
    manifest={'release':'DepMap Public 25Q3','release_evidence':'Local README.txt first line',
        'data_type':'CRISPRGeneEffect.csv: integrated, corrected Chronos gene-effect scores',
        'models':len(ge),'genes':ge.shape[1],'sex_counts':model.Sex.value_counts(dropna=False).to_dict(),
        'lineages':int(model.OncotreeLineage.nunique()),
        'cohort':'All models in the local CRISPRGeneEffect matrix; no expression or LOY-classification restriction.',
        'screen_hits_total':len(hit),'depleted_hits':int(hit.direction.eq('Depleted').sum()),'enriched_hits':int(hit.direction.eq('Enriched').sum()),
        'mapped_screen_hits':int(annot.loc[hit.GENE,'depmap_available'].sum()),
        'missing_screen_hits':annot.loc[hit.GENE].query('not depmap_available').index.tolist(),
        'saga_reference_members':20,'saga_alternative_members':4,'all_seed_genes':len(seeds),'available_seed_genes':len(present),
        'correlation':'Pearson, pairwise-complete observations; minimum 200 shared models; self-correlations excluded from partner rankings.',
        'multiple_testing':'No P-value or FDR threshold is used; correlation rankings and r>=0.3 counts are descriptive only. Screen FDR labels retain the original drugZ family.',
        'sensitivity':'Per-gene OLS residualisation on intercept + OncotreeLineage + recorded Sex (Male/Female/Unknown); fit using available models per gene, then pairwise-complete Pearson correlation. Not a LOY-specificity test.',
        'residual_design_rank_range':[min(ranks.values()),max(ranks.values())],
        'limitations':'Shared dependency profiles can reflect biology, lineage, sex, other context or residual screen artefacts; correlation is not direct binding or evidence of LOY-specific synthetic lethality. No LOY subgroup contrast performed.',
        'validation':{'numpy_corrcoef_pairs_checked':len(errors),'max_corrcoef_error':max(errors),
                      'historical_partner_pairs_checked':len(oldcheck),'max_historical_r_difference':max(oldcheck)},
        'sources':{str(p):{'sha256':sha(p),'bytes':p.stat().st_size} for p in [DEPMAP/'CRISPRGeneEffect.csv',DEPMAP/'Model.csv',DEPMAP/'README.txt',CAND/'depleted_candidates.csv',CAND/'enriched_candidates.csv',CAND/'saga_reference_scores.csv',oldpath]}}
    (DATA/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in manifest.items() if k not in ['sources']},indent=2))
    print(pd.DataFrame(summary).query("gene in ['DDX3X','TADA2B','TAF6L','ATXN7L3','SREBF1','SREBF2','USP7']").to_string(index=False))


if __name__=='__main__':
    main()
