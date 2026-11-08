"""Independently verify the exclusion, ranks, GSEA statistics and report links."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'integration/pathway_GSEA';OUT=BASE/'ribosomal_sensitivity'


def bh(p):
    p=np.asarray(p,float);ix=np.argsort(p);out=np.empty(len(p));out[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return out


class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):self.links.extend(v for k,v in attrs if k in ('href','src') and v)


def main():
    base=pd.read_csv(BASE/'gene_ranks.csv',dtype={'entrez_id':str}).set_index('entrez_id')
    excluded=pd.read_csv(OUT/'excluded_genes.csv',dtype={'entrez_id':str})
    catalog=pd.read_csv(OUT/'sources/HGNC_ribosomal_protein_catalog.csv',dtype={'entrez_id':str})
    sets=pd.read_csv(BASE/'sources/gene_set_membership.csv',dtype={'entrez_id':str}).groupby('pathway').entrez_id.agg(set).to_dict()
    results=pd.read_csv(OUT/'fgsea_exclusion_results.csv')
    checked=0;max_error=0.;universes={}
    for mode in ['no_cytosolic','no_cytosolic_or_mitochondrial']:
        groups=catalog[catalog['class'].eq('cytosolic')] if mode=='no_cytosolic' else catalog
        mask=base.index.isin(groups.entrez_id)|base.symbol.isin(groups.symbol)
        removed=set(base.index[mask]);expected=set(excluded.loc[excluded['class'].eq('cytosolic'),'entrez_id']) if mode=='no_cytosolic' else set(excluded.entrez_id)
        assert removed==expected
        assert not base.loc[list(removed),'symbol'].str.match(r'^RPS6K').any()
        rank=pd.read_csv(OUT/f'gene_ranks_{mode}.csv',dtype={'entrez_id':str}).set_index('entrez_id')
        pd.testing.assert_frame_equal(rank,base.loc[~mask],check_exact=False,rtol=1e-12,atol=1e-12)
        universes[mode]=len(rank)
        for analysis in ['ATAC','RNA_Combined','RNA_AUG_2025','RNA_MAY_2026']:
            ordered=rank.sort_values(analysis,ascending=False,kind='stable')
            ids=ordered.index.to_numpy();scores=ordered[analysis].to_numpy();positions={g:i for i,g in enumerate(ids)};n=len(ids)
            frame=results[results.exclusion.eq(mode)&results.analysis.eq(analysis)]
            eligible={p for p,g in sets.items() if 15<=len(g&set(ids))<=500}
            assert eligible==set(frame.pathway)
            for r in frame.itertuples():
                pos=np.array(sorted(positions[g] for g in sets[r.pathway] if g in positions),dtype=int)
                assert len(pos)==r.size
                steps=np.full(n,-1/(n-len(pos)));weight=np.abs(scores[pos]);steps[pos]=weight/weight.sum();curve=np.cumsum(steps)
                imax=int(np.argmax(curve));imin=int(np.argmin(curve));es=curve[imax] if curve[imax]>-curve[imin] else curve[imin]
                error=float(abs(es-r.ES));max_error=max(max_error,error);assert error<1e-9,(mode,analysis,r.pathway,error)
                leading=set(ids[pos[pos<=imax]]) if es>0 else set(ids[pos[pos>imin]])
                assert leading==set(str(r.leading_edge_entrez).split(';'))
                assert not leading&removed
                checked+=1
            np.testing.assert_allclose(bh(frame.pval),frame.padj_global,atol=1e-12)
            for _,g in frame.groupby('collection'):np.testing.assert_allclose(bh(g.pval),g.padj,atol=1e-12)
    all_results=pd.read_csv(OUT/'all_scenarios_results.csv')
    for _,frame in all_results.groupby(['exclusion','analysis','collection']):
        frame=frame[frame.FDR_common_pathway_family.notna()]
        np.testing.assert_allclose(bh(frame.pval),frame.FDR_common_pathway_family,atol=1e-12)
    original=pd.read_csv(BASE/'fgsea_all_results.csv');original=original[original.analysis.isin(['ATAC','RNA_Combined','RNA_AUG_2025','RNA_MAY_2026'])]
    got=all_results[all_results.exclusion.eq('original')][original.columns]
    pd.testing.assert_frame_equal(got.reset_index(drop=True),original.reset_index(drop=True),check_exact=False,rtol=1e-12,atol=1e-12)
    completed=json.loads((BASE/'ANALYSIS_COMPLETE.json').read_text())
    assert hashlib.sha256((BASE/'fgsea_all_results.csv').read_bytes()).hexdigest()==completed['source_checksums']['integration\\pathway_GSEA\\fgsea_all_results.csv']
    focus=all_results[all_results.analysis.eq('ATAC')&all_results.pathway.isin(['HALLMARK_MYC_TARGETS_V1','HALLMARK_MYC_TARGETS_V2','HALLMARK_E2F_TARGETS'])]
    assert len(focus)==9 and focus.NES.gt(0).all() and focus.padj.lt(.05).all()
    page=ROOT/'report/ribosomal_GSEA_sensitivity.html';parser=Links();parser.feed(page.read_text(encoding='utf-8'))
    for link in parser.links:
        u=urlsplit(link)
        if not u.scheme and u.path:assert (page.parent/unquote(u.path)).is_file(),link
    result={'status':'passed','remaining_genes':universes,'HGNC_exclusion_verified':True,'retained_ranks_unchanged':True,'original_GSEA_results_unchanged':True,
            'independent_ES_and_leading_edge_checks':checked,'max_ES_absolute_error':max_error,'all_FDR_recalculations_match':True,'local_report_links_exist':True}
    (OUT/'verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    # Publish a link without replacing the original analysis or its figures.
    original_page=ROOT/'report/pathway_GSEA.html';html=original_page.read_text(encoding='utf-8')
    section='<section id="ribosomal-sensitivity"><h2>Ribosomal protein exclusion sensitivity</h2><p><a href="ribosomal_GSEA_sensitivity.html">Compare the original GSEA with cytosolic-only and cytosolic-plus-mitochondrial ribosomal protein exclusions</a>. Includes unchanged gene statistics, recomputed pathway FDR, explicit excluded-gene lists and size-filter audits.</p></section>'
    if 'id="ribosomal-sensitivity"' not in html:original_page.write_text(html.replace('</body>',section+'\n</body>'),encoding='utf-8')
    complete=json.loads((OUT/'completion.json').read_text());complete['status']='complete';complete['verification']=result
    complete['file_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'code/61_ribosomal_gene_sets.py',ROOT/'code/62_ribosomal_GSEA.R',ROOT/'code/63_report_ribosomal_GSEA.py',Path(__file__),OUT/'excluded_genes.csv']}
    (OUT/'ANALYSIS_COMPLETE.json').write_text(json.dumps(complete,indent=2),encoding='utf-8');(OUT/'completion.json').write_text(json.dumps(complete,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
