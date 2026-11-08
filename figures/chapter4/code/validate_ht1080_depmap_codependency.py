"""Validate co-dependency provenance, saved matrices, rankings and figure exports."""
from pathlib import Path
import sys,json,hashlib,xml.etree.ElementTree as ET
THESIS=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(THESIS.parent/'.figure_python_deps'))
import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/ht1080_depmap_codependency'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as fh:
        for b in iter(lambda:fh.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    manifest=json.loads((DATA/'analysis_manifest.json').read_text())
    for path,entry in manifest['sources'].items(): assert sha(Path(path))==entry['sha256'],path
    a=pd.read_csv(DATA/'seed_annotations.csv').set_index('gene')
    assert len(a)==75 and a.screen_hit.sum()==57 and a.depmap_available.sum()==74
    assert a.loc[~a.depmap_available].index.tolist()==['ZCCHC16']
    assert a.loc[a.depmap_available,'depmap_n'].eq(1186).all()
    x=pd.read_csv(DATA/'seed_gene_effect.csv',index_col=0)
    r=pd.read_csv(DATA/'seed_pairwise_pearson.csv',index_col=0)
    ref=x.corr()
    assert np.allclose(r.loc[x.columns,x.columns],ref,atol=1e-10)
    assert r.loc['ZCCHC16'].isna().all() and r['ZCCHC16'].isna().all()
    m=pd.read_csv(DATA/'model_cohort.csv',index_col=0).loc[x.index]
    cov=pd.get_dummies(m[['OncotreeLineage','Sex']].fillna('Unknown'),drop_first=True,dtype=float)
    design=np.column_stack([np.ones(len(m)),cov.to_numpy()])
    resid=pd.read_csv(DATA/'seed_lineage_sex_residuals.csv',index_col=0).loc[x.index,x.columns]
    assert np.max(np.abs(design.T@resid.to_numpy()))<1e-7
    adjusted=pd.read_csv(DATA/'seed_pairwise_lineage_sex_residual_pearson.csv',index_col=0)
    assert np.allclose(adjusted.loc[x.columns,x.columns],resid.corr(),atol=1e-10)
    top=pd.read_csv(DATA/'top_partners.csv')
    assert top.n_models.ge(200).all() and top.seed.ne(top.partner).all()
    assert top.loc[top.selection.eq('Top positive'),'pearson_r'].gt(0).all()
    assert top.loc[top.selection.eq('Top negative'),'pearson_r'].lt(0).all()
    result={'source_hashes_unchanged':True,'mapped_hits':56,'unavailable_hit':'ZCCHC16',
            'pairwise_matrix_matches_pandas_corr':True,'residuals_orthogonal_to_covariates':True,
            'exports':{}}
    cluster=json.loads((DATA/'hit_clustering_manifest.json').read_text())
    assert sha(DATA/'seed_pairwise_pearson.csv')==cluster['source_sha256']
    genes=cluster['input_order']
    tree=pd.read_csv(DATA/'hit_heatmap_linkage.csv')
    clusters={i:[g] for i,g in enumerate(genes)}
    assert len(genes)==56 and len(tree)==55
    # Independently check each saved merge distance directly from the original
    # inter-cluster correlations, and require the minimum eligible distance.
    for i,row in enumerate(tree.itertuples()):
        l,rr=int(row.left_cluster),int(row.right_cluster)
        direct=1-r.loc[clusters[l],clusters[rr]].to_numpy().mean()
        assert np.isclose(row.distance_1_minus_r,direct,atol=1e-12)
        candidates=[1-r.loc[clusters[u],clusters[v]].to_numpy().mean()
                    for u in clusters for v in clusters if u<v]
        assert np.isclose(direct,min(candidates),atol=1e-12)
        members=clusters.pop(l)+clusters.pop(rr)
        assert len(members)==row.size
        clusters[len(genes)+i]=members
    leaves=next(iter(clusters.values()))
    assert leaves+['ZCCHC16']==cluster['display_order']
    assert pd.read_csv(DATA/'hit_heatmap_order.csv').gene.tolist()==cluster['display_order']
    result['clustering']={'genes':56,'merges':55,'all_average_linkage_distances_verified':True,
                          'both_axes_use_same_tree':cluster['same_tree_for_rows_and_columns']}
    stems={'ch04_supp_ht1080_depmap_hit_codependencies_draft':1,
           'ch04_supp_ht1080_depmap_saga_codependencies_draft':1,
           'ch04_supp_ht1080_depmap_hit_partner_atlas':8,
           'ch04_supp_ht1080_depmap_saga_partner_atlas':3}
    for stem,pages in stems.items():
        pdf=fitz.open(OUT/f'{stem}.pdf'); assert len(pdf)==pages
        for i,page in enumerate(pdf):
            words=page.get_text('words')
            outside=[w[:5] for w in words if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
            assert not outside,(stem,i,outside)
            assert all('OpenSans' in f[3] for f in page.get_fonts()),(stem,page.get_fonts())
        alltext='\n'.join(p.get_text() for p in pdf)
        required=a.index[a.screen_hit].tolist() if 'hit_' in stem else a.index[a.saga_module.notna()].tolist()
        assert all(g in alltext for g in required),(stem,[g for g in required if g not in alltext])
        result['exports'][stem]={'pages':len(pdf),'out_of_bounds_words':0,'all_requested_genes_labelled':True,'all_fonts_open_sans':True}
        if pages==1:
            tree=ET.parse(OUT/f'{stem}.svg')
            assert len(tree.findall('.//{http://www.w3.org/2000/svg}text'))>50
            png=Image.open(OUT/f'{stem}.png'); assert min(png.info['dpi'])>590
            pdf[0].get_pixmap(matrix=fitz.Matrix(1.4,1.4),colorspace=fitz.csGRAY).save(OUT/'review'/f'{stem}_greyscale_check.png')
        pdf.close()
    (DATA/'export_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
