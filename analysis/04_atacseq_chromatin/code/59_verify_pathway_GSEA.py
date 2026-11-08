"""Independent source-rank, enrichment-score, leading-edge, FDR and report-link checks."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'integration/pathway_GSEA'


def read(path,**kw):return pd.read_csv(ROOT/path,**kw)


def bh(p):
    p=np.array(p,float);ix=np.argsort(p);q=np.empty(len(p));q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q


class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):self.links.extend(v for k,v in attrs if k in ('href','src') and v)


def main():
    ranks=pd.read_csv(OUT/'gene_ranks.csv',dtype={'entrez_id':str}).set_index('entrez_id')
    a=read('integration/feature_scatters/feature_gene_effects.csv',dtype={'entrez_id':str})
    b=read('integration/batch_RNA_overlaps/all_gene_results.csv',dtype={'entrez_id':str})
    universe=set(a.loc[a.feature.eq('promoter'),'entrez_id'])
    for batch in ('AUG_2025','MAY_2026'):universe&=set(b.loc[b.feature.eq('promoter')&b.batch.eq(batch),'entrez_id'])
    assert universe==set(ranks.index) and len(universe)==11891
    links=read('integration/feature_overlaps/peak_gene_link_audit.csv',dtype={'entrez_id':str})
    links=links[links.feature.eq('promoter')&links.entrez_id.isin(universe)][['entrez_id','peak_id']].drop_duplicates()
    for model,column in [('primary','ATAC'),('passage','ATAC_passage')]:
        peaks=read(f'differential/{model}_TMM_all.csv')
        merged=links.merge(peaks,on='peak_id',validate='many_to_one')
        merged['score']=np.sign(merged.logFC)*np.sqrt(merged.F)
        score=merged.groupby('entrez_id').score.median().reindex(ranks.index)
        np.testing.assert_allclose(score,ranks[column],atol=1e-12,rtol=1e-12)
        if model=='primary':
            np.testing.assert_allclose(merged.groupby('entrez_id').logFC.median().reindex(ranks.index),ranks.ATAC_median_logFC,atol=1e-12)
            representative=merged.sort_values(['logCPM','peak_id'],ascending=[False,True]).drop_duplicates('entrez_id').set_index('entrez_id').reindex(ranks.index)
            assert representative.peak_id.tolist()==ranks.representative_peak.tolist()
            np.testing.assert_allclose(representative.score,ranks.ATAC_representative,atol=1e-12)
    for name,path in [('RNA_Combined','config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv'),
        ('RNA_AUG_2025','integration/batch_RNA_overlaps/sources/ht1080_AUG_2025_LOY_v_CTRL_DESeq2.csv'),
        ('RNA_MAY_2026','integration/batch_RNA_overlaps/sources/ht1080_MAY_2026_LOY_v_CTRL_DESeq2.csv')]:
        frame=read(path,dtype={'entrezid':str});frame=frame[frame.entrezid.isin(universe)].set_index('entrezid').reindex(ranks.index)
        assert not frame.index.duplicated().any()
        np.testing.assert_allclose(frame.t,ranks[name],atol=1e-12)
        np.testing.assert_allclose(frame.logFC,ranks[name+'_logFC'],atol=1e-12)
        assert np.array_equal(np.sign(frame.t),np.sign(frame.logFC))
    results=pd.read_csv(OUT/'fgsea_all_results.csv')
    for (analysis,collection),g in results.groupby(['analysis','collection']):
        assert np.isfinite(g[['pval','NES','ES','padj']].to_numpy()).all()
        np.testing.assert_allclose(bh(g.pval),g.padj,atol=1e-12)
    for analysis,g in results.groupby('analysis'):
        np.testing.assert_allclose(bh(g.pval),g.padj_global,atol=1e-12)
    sets=pd.read_csv(OUT/'sources/gene_set_membership.csv',dtype={'entrez_id':str}).groupby('pathway').entrez_id.agg(set).to_dict()
    max_es_error=0.;checked=0
    for analysis in ('ATAC','RNA_Combined','RNA_AUG_2025','RNA_MAY_2026'):
        ordered=ranks.sort_values(analysis,ascending=False,kind='stable')
        gene_ids=ordered.index.to_numpy();scores=ordered[analysis].to_numpy();n=len(scores)
        # gene_ranks.csv is already in the deterministic primary gene-ID tie order.
        for row in results[results.analysis.eq(analysis)].itertuples():
            hit=np.isin(gene_ids,list(sets[row.pathway]));positions=np.flatnonzero(hit)
            assert hit.sum()==row.size and 15<=row.size<=500
            weight=np.abs(scores[hit]);steps=np.full(n,-1/(n-hit.sum()),float);steps[hit]=weight/weight.sum()
            curve=np.cumsum(steps)
            imax=int(np.argmax(curve));imin=int(np.argmin(curve))
            es=curve[imax] if curve[imax]>-curve[imin] else curve[imin]
            max_es_error=max(max_es_error,float(abs(es-row.ES)))
            assert abs(es-row.ES)<1e-9,(analysis,row.pathway,es,row.ES)
            leading=set(gene_ids[positions[positions<=imax]]) if es>0 else set(gene_ids[positions[positions>imin]])
            stored=set(str(row.leading_edge_entrez).split(';'))
            assert leading==stored,(analysis,row.pathway,'leading-edge mismatch',leading^stored)
            checked+=1
    sizes=pd.read_csv(OUT/'pathway_sizes_ATAC.csv')
    filtered=sizes[~sizes.tested].merge(pd.read_csv(OUT/'sources/pathway_annotation.csv'),on='pathway')
    filtered.to_csv(OUT/'excluded_pathways.csv',index=False)
    pairs=pd.read_csv(OUT/'pathway_concordance.csv')
    edges=pd.read_csv(OUT/'shared_leading_edge_genes.csv',dtype={'entrez_id':str})
    result_index=results.set_index(['analysis','pathway'])
    for r in pairs.itertuples():
        aa=result_index.loc[('ATAC',r.pathway)];rr=result_index.loc[(r.RNA_analysis,r.pathway)]
        expected=set(str(aa.leading_edge_entrez).split(';'))&set(str(rr.leading_edge_entrez).split(';'))
        got=set(edges.loc[edges.RNA_analysis.eq(r.RNA_analysis)&edges.pathway.eq(r.pathway),'entrez_id'])
        assert got==expected and len(got)==r.shared_leading_edge_genes
    summary=pd.read_csv(OUT/'pathway_concordance_summary.csv')
    assert summary.loc[summary.RNA_analysis.isin(['RNA_Combined','RNA_MAY_2026']),'both_FDR005'].sum()==0
    assert summary.loc[summary.RNA_analysis.eq('RNA_AUG_2025'),'both_concordant_FDR005'].sum()==3
    assert summary.loc[summary.RNA_analysis.eq('RNA_AUG_2025'),'both_opposing_FDR005'].sum()==1
    page=ROOT/'report/pathway_GSEA.html';parser=Links();parser.feed(page.read_text(encoding='utf-8'))
    for link in parser.links:
        url=urlsplit(link)
        if not url.scheme and url.path:assert (page.parent/unquote(url.path)).is_file(),link
    result={'status':'passed','common_gene_universe':len(universe),'all_ATAC_RNA_ranks_match_sources':True,
            'independent_enrichment_scores_and_leading_edges_checked':checked,'max_ES_absolute_error':max_es_error,
            'BH_corrections_match':True,'leading_edge_intersections_checked':len(pairs),'all_local_report_links_exist':True,
            'excluded_Hallmark_sets':filtered.loc[filtered.collection.eq('Hallmark'),['pathway','eligible_genes']].to_dict('records'),
            'limitation':'Gene-set permutation tests do not explicitly correct shared-promoter or inter-gene correlation.'}
    (OUT/'verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    completion=json.loads((OUT/'completion.json').read_text());completion['status']='complete';completion['verification']=result
    (OUT/'completion.json').write_text(json.dumps(completion,indent=2),encoding='utf-8')
    (OUT/'ANALYSIS_COMPLETE.json').write_text(json.dumps(completion,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
