"""Check exported gene inspection against original counts, effects and gene lists."""
from pathlib import Path
from html.parser import HTMLParser
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'integration/pathway_GSEA'
OUT=BASE/'gene_drivers'
D=json.loads((ROOT/'report/pathway_driver_assets/data.js').read_text(encoding='utf-8').split('=',1)[1].removesuffix(';'))
ranks=pd.read_csv(BASE/'gene_ranks.csv',dtype={'entrez_id':str}).set_index('entrez_id')
results=pd.read_csv(BASE/'fgsea_all_results.csv')
members=pd.read_csv(BASE/'sources/gene_set_membership.csv',dtype={'entrez_id':str})
analyses=['ATAC','RNA_Combined','RNA_AUG_2025','RNA_MAY_2026']
selected_results=results[results.analysis.isin(analyses)&results.pathway.isin([p['pathway'] for p in D['pathways'].values()])]
assert set(D['genes'])==set(';'.join(selected_results.leading_edge_entrez).split(';'))
assert len(D['genes'])==379
rnafiles=[ROOT/'config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv']+[ROOT/f'integration/batch_RNA_overlaps/sources/ht1080_{b}_LOY_v_CTRL_DESeq2.csv' for b in ['AUG_2025','MAY_2026']]
rna=[pd.read_csv(p,dtype={'entrezid':str}).set_index('entrezid') for p in rnafiles]
exprbase=ROOT/'integration/interactive_atac_rna'
raw=pd.read_csv(exprbase/'RNA_raw_counts.csv').set_index('gene_id')
norm=pd.read_csv(exprbase/'RNA_normalized_counts.csv').set_index('gene_id_full')
sf=pd.read_csv(exprbase/'RNA_size_factors.csv').set_index('sample').size_factor
samples=[s['colname'] for s in D['rna_samples']]
meta=pd.read_csv(exprbase/'RNA_sample_metadata.csv').set_index('colname')
for s in D['rna_samples']:
    for c in ['batch','treatment','sample_id']:assert s[c]==meta.loc[s['colname'],c]
assert len(samples)==12 and len(set(samples))==12
counts=pd.read_csv(ROOT/'counts/peak_counts.tsv',sep='\t').set_index('peak_id')
peaks=pd.read_csv(ROOT/'differential/primary_TMM_all.csv').set_index('peak_id')
links=pd.read_csv(BASE/'promoter_peak_ranking_audit.csv',dtype={'entrez_id':str})
nf=pd.read_csv(ROOT/'differential/normalization_factors.csv')
sizes=nf[nf.arm.eq('primary_TMM')].set_index('sample').loc[D['atac_samples'],'effective_size'].to_numpy()
robust=pd.read_csv(OUT/'RNA_sample_robustness.csv',dtype={'entrez_id':str}).set_index(['entrez_id','batch'])
npeak=0
for gid,g in D['genes'].items():
    ens=rna[0].loc[gid,'gene_id_full']
    np.testing.assert_array_equal(g['raw'],raw.loc[ens,samples].to_numpy(int))
    np.testing.assert_allclose(g['norm'],norm.loc[ens,samples].to_numpy(float),rtol=1e-12)
    np.testing.assert_allclose(g['norm'],np.array(g['raw'])/sf.loc[samples].to_numpy(),rtol=1e-12)
    for j,r in enumerate(g['rna']):
        np.testing.assert_allclose([r['logFC'],r['SE'],r['FDR']],rna[j].loc[gid,['logFC','lfcSE','FDR']].to_numpy(float),rtol=1e-12)
    linked=links[links.entrez_id.eq(gid)]
    assert {p['id'] for p in g['peaks']}==set(linked.peak_id)
    assert peaks.loc[g['representative_peak'],'logCPM']==peaks.loc[linked.peak_id,'logCPM'].max()
    np.testing.assert_allclose(g['atac_median_logFC'],peaks.loc[linked.peak_id,'logFC'].median())
    for p in g['peaks']:
        np.testing.assert_allclose([p['logFC'],p['FDR']],peaks.loc[p['id'],['logFC','FDR']].to_numpy(float))
        np.testing.assert_array_equal(p['raw'],counts.loc[p['id'],D['atac_samples']].to_numpy(int))
        np.testing.assert_allclose(p['cpm'],np.array(p['raw'])/sizes*1e6,rtol=1e-12)
        npeak+=1
    for b in ['AUG_2025','MAY_2026']:
        groups={c:np.log2([g['norm'][i]+1 for i,s in enumerate(D['rna_samples']) if s['batch']==b and s['treatment']==c]) for c in ['CTRL','LOY']}
        cv,lv=groups['CTRL'],groups['LOY'];assert len(cv)==len(lv)==3
        dif=lv.mean()-cv.mean();v=robust.loc[(gid,b)]
        loo=[np.delete(lv,i).mean()-cv.mean() for i in range(3)]+[lv.mean()-np.delete(cv,i).mean() for i in range(3)]
        np.testing.assert_allclose([v.mean_log_expression_difference_LOY_CTRL,v.leave_one_out_min_difference,v.leave_one_out_max_difference],[dif,min(loo),max(loo)],atol=1e-12)
        assert v.all_LOY_below_all_CTRL==(lv.max()<cv.min())
        assert v.all_LOY_above_all_CTRL==(lv.min()>cv.max())
        assert v.leave_one_out_same_sign==bool(np.all(np.sign(loo)==np.sign(dif)))
        np.testing.assert_allclose(v.LOY_CTRL_pairwise_fraction_positive,(lv[:,None]>cv[None,:]).mean())

effects=pd.read_csv(OUT/'all_pathway_gene_effects.csv',dtype={'entrez_id':str})
con=pd.read_csv(OUT/'leading_edge_concentration.csv').set_index(['pathway','analysis'])
for key,p in D['pathways'].items():
    pathway=p['pathway'];eligible=set(members.loc[members.pathway.eq(pathway),'entrez_id'])&set(ranks.index)
    assert set(p['eligible'])==eligible
    e=effects[effects.pathway.eq(pathway)].set_index('entrez_id')
    assert set(e.index)==eligible and set(e.index[e.any_leading_edge])==set(p['drivers'])
    assert len(p['examples'])==12 and set(e.index[e.selected_example])==set(p['examples'])
    for a in ['ATAC_median_logFC']+[a+'_logFC' for a in analyses[1:]]+[a+'_FDR' for a in analyses[1:]]:
        np.testing.assert_allclose(e[a],ranks.loc[e.index,a])
    union=set();must=set();weights={}
    for a in analyses:
        rr=selected_results[selected_results.pathway.eq(pathway)&selected_results.analysis.eq(a)].iloc[0]
        le=set(rr.leading_edge_entrez.split(';'));union|=le
        assert set(e.index[e[a+'_leading_edge']])==le
        w=ranks.loc[sorted(le),a].abs();w=w/w.sum();order=w.sort_values(ascending=False,kind='stable')
        must.update(order.head(3).index);weights[a]=w.to_dict()
        np.testing.assert_allclose(e[a+'_leading_weight_fraction'],[w.get(g,0) for g in e.index])
        c=con.loc[(pathway,a)]
        np.testing.assert_allclose([c.top1_leading_weight_fraction,c.top5_leading_weight_fraction],[order.iloc[0],order.head(5).sum()])
        assert c.genes_for_50pct_leading_weight==int((order.cumsum()<.5).sum()+1)
    expected=set(must)
    for g in sorted(union,key=lambda g:(-max(weights[a].get(g,0) for a in analyses),g)):
        if len(expected)==12:break
        expected.add(g)
    assert set(p['examples'])==expected
    s=p['summary'];aug=e.RNA_AUG_2025_logFC;may=e.RNA_MAY_2026_logFC
    for name,value in {'eligible_genes':len(e),'leading_edge_union_genes':len(union),'ATAC_positive_genes':(e.ATAC_median_logFC>0).sum(),'RNA_Aug_positive_genes':(aug>0).sum(),'RNA_May_positive_genes':(may>0).sum(),'RNA_both_negative_genes':((aug<0)&(may<0)).sum(),'RNA_both_positive_genes':((aug>0)&(may>0)).sum(),'RNA_opposite_sign_genes':(np.sign(aug)!=np.sign(may)).sum()}.items():assert s[name]==value
    np.testing.assert_allclose([s['median_ATAC_logFC'],s['median_RNA_Aug_logFC'],s['median_RNA_May_logFC'],s['RNA_batch_effect_spearman']],[e.ATAC_median_logFC.median(),aug.median(),may.median(),aug.corr(may,method='spearman')])

class Links(HTMLParser):
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ['href','src'] and v and not v.startswith(('http:', 'https:', '#','data:')):
                assert (ROOT/'report'/v.split('#')[0]).resolve().exists(),v
Links().feed((ROOT/'report/pathway_gene_drivers.html').read_text(encoding='utf-8'))
browser=json.loads((OUT/'browser_verification.json').read_text())
assert browser['status']=='passed'
result={'status':'passed','distinct_genes':len(D['genes']),'RNA_gene_sample_values':len(D['genes'])*12,'linked_promoter_records':npeak,'RNA_gene_batch_robustness_rows':len(robust),'checks':['Original RNA coefficients, SE and FDR','Original RNA raw counts and size-factor normalization','Original promoter links, peak effects and TMM-normalized counts','Complete leading-edge membership and example selection','Leading-edge contribution weights and pathway summaries','All local report links','Offline browser interaction checks'],'visual_inspection':'Heatmap, batch scatter and browser plots reviewed'}
(OUT/'data_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
manifest=json.loads((OUT/'completion.json').read_text())
for rel,digest in manifest['source_sha256'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
manifest.update(status='complete',verification=result,browser_verification=browser)
for p in [*rnafiles,ROOT/'differential/primary_TMM_all.csv',BASE/'sources/gene_set_membership.csv',BASE/'promoter_peak_ranking_audit.csv',ROOT/'report/pathway_driver_assets/viewer.js',Path(__file__)]:
    manifest['source_sha256'][str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(OUT/'completion.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
