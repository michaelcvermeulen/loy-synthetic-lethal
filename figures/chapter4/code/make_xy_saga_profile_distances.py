"""Distances from X-linked gametolog profiles to prespecified SAGA hit centroids.

Uses the frozen 56-feature correlation profiles, not UMAP coordinates. Distances
have the same feature-count normalisation as the UMAP input graph.
"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
import pandas as pd

THESIS=Path(__file__).resolve().parents[3]
OUT=THESIS/'figures/chapter4'
SOURCE=OUT/'source_data/ht1080_hit_profile_umap'
DATA=OUT/'source_data/xy_saga_profile_distance'
DATA.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pymupdf as fitz
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
STEM='ch04_supp_xy_saga_profile_distance_draft'
FIVE=['TAF6L','TADA2B','SUPT20H','TADA1','TAF5L']
SIX=FIVE+['ATXN7L3']


def distances(profiles,centroid):
    delta=profiles.sub(centroid,axis=1)
    n=delta.notna().sum(axis=1)
    assert n.ge(55).all()
    return np.sqrt(56/n*delta.pow(2).sum(axis=1))


def main():
    z=np.load(SOURCE/'input_profiles.npz')
    p=pd.DataFrame(z['distance_r'],index=z['genes'],columns=z['features'])
    ann=pd.read_csv(SOURCE/'gene_annotations.csv',index_col=0)
    refpath=OUT/'source_data/xy_paralog_codependency/reference_gametolog_pairs.csv'
    ref=pd.read_csv(refpath)
    genes=ref.loc[ref.Gametolog.eq('X'),'common_name'].drop_duplicates().tolist()
    assert len(genes)==16 and set(genes).issubset(p.index)
    centres=pd.DataFrame({'five_hit_centroid':p.loc[FIVE].mean(axis=0),
                          'six_hit_centroid':p.loc[SIX].mean(axis=0)})
    support=pd.DataFrame({'five_hit_centroid':p.loc[FIVE].notna().sum(axis=0),
                         'six_hit_centroid':p.loc[SIX].notna().sum(axis=0)})
    assert centres.notna().all().all()
    assert support.five_hit_centroid.value_counts().to_dict()=={5:51,4:5}
    assert support.six_hit_centroid.value_counts().to_dict()=={6:50,5:6}
    centres.to_csv(DATA/'centroid_profiles.csv',index_label='feature')
    support.to_csv(DATA/'centroid_feature_support.csv',index_label='feature')
    p.loc[SIX].to_csv(DATA/'saga_member_profiles.csv',index_label='gene')
    p.loc[genes].to_csv(DATA/'x_gametolog_profiles.csv',index_label='gene')
    all_dist=pd.DataFrame({'five_hit_distance':distances(p,centres.five_hit_centroid),
                           'six_hit_distance':distances(p,centres.six_hit_centroid)},index=p.index)
    all_dist.index.name='gene'
    all_dist.to_csv(DATA/'all_gene_centroid_distances.csv')
    result=all_dist.loc[genes].copy()
    result['rank_X_five']=result.five_hit_distance.rank(method='min').astype(int)
    result['rank_X_six']=result.six_hit_distance.rank(method='min').astype(int)
    result['n_features']=p.loc[genes].notna().sum(axis=1)
    result['n_models']=ann.loc[genes,'n_models']
    # Additional audit: distances to individual members use the exact same
    # pairwise self-exclusion rule as the UMAP neighbour input.
    pair_rows=[]
    for g in genes:
        for s in SIX:
            valid=p.loc[g].notna() & p.loc[s].notna()
            d=np.sqrt(56/valid.sum()*((p.loc[g,valid]-p.loc[s,valid])**2).sum())
            pair_rows.append(dict(gene=g,saga_gene=s,profile_distance=float(d),n_features=int(valid.sum())))
    pair=pd.DataFrame(pair_rows)
    pair.to_csv(DATA/'individual_saga_member_distances.csv',index=False)
    mean_five=pair.loc[pair.saga_gene.isin(FIVE)].groupby('gene').profile_distance.mean()
    result['mean_distance_to_five_members']=mean_five
    mean_six=pair.groupby('gene').profile_distance.mean()
    result['mean_distance_to_six_members']=mean_six
    # Common-feature sensitivity: omit the DDX3X feature from every X gene so
    # all 16 comparisons use the same 55 observed coordinates.
    common=p.columns[p.columns!='DDX3X']
    assert p.loc[genes,common].notna().all().all()
    sensitivity=np.sqrt(56/55*p.loc[genes,common].sub(centres.loc[common,'five_hit_centroid'],axis=1).pow(2).sum(axis=1))
    result['five_hit_distance_common_55_features']=sensitivity
    result['rank_X_common_55_five']=sensitivity.rank(method='min').astype(int)
    sensitivity_six=np.sqrt(56/55*p.loc[genes,common].sub(centres.loc[common,'six_hit_centroid'],axis=1).pow(2).sum(axis=1))
    result['six_hit_distance_common_55_features']=sensitivity_six
    result['rank_X_common_55_six']=sensitivity_six.rank(method='min').astype(int)
    result=result.sort_values(['six_hit_distance','gene'])
    result.to_csv(DATA/'x_gametolog_saga_distances.csv')
    # Validate the centroid without matrix reductions and the distance without
    # pandas operations; no input correlations are refitted or modified.
    errors=[]
    for members,column in [(FIVE,'five_hit_centroid'),(SIX,'six_hit_centroid')]:
        for feature in p.columns:
            values=[float(p.loc[g,feature]) for g in members if np.isfinite(p.loc[g,feature])]
            errors.append(abs(sum(values)/len(values)-centres.loc[feature,column]))
        for gene in genes:
            values=[(float(p.loc[gene,f])-float(centres.loc[f,column]))**2 for f in p.columns if np.isfinite(p.loc[gene,f])]
            direct=(56/len(values)*sum(values))**.5
            errors.append(abs(direct-all_dist.loc[gene,column.replace('centroid','distance')]))
    assert max(errors)<1e-12
    assert result.index[0]=='KDM5C'
    manifest=dict(primary_cluster=SIX,sensitivity_cluster=FIVE,display='Six-hit reference only; five-hit sensitivity retained in source tables',features=56,genes=genes,
                  centroid='Feature-wise arithmetic mean of observed member correlations; self-comparisons excluded, not imputed',
                  primary_centroid_support='Six members at 50 features; five other members at each of the six self-feature coordinates',
                  metric='sqrt(56 / number of observed gene-centroid features * sum squared differences)',
                  query_features='55 for DDX3X after self exclusion; 56 for every other X gene',
                  n_models_range=[int(result.n_models.min()),int(result.n_models.max())],
                  coordinates='Original raw signed 56-feature profiles, not 2D UMAP distances',
                  refit=False,new_significance_tests=False,
                  independent_validation_max_error=max(errors),
                  same_X_ranking_with_ATXN7L3=bool(result.rank_X_five.equals(result.rank_X_six)),
                  closest_common_55=sensitivity_six.idxmin(),closest_mean_member=mean_six.idxmin(),
                  interpretation='Ranks describe similarity to prespecified SAGA profiles, not SAGA membership, enrichment or a test of LOY specificity',
                  circularity='The distance includes the SAGA-hit correlation features that define the reference group; this is not independent evidence of association',
                  source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [SOURCE/'input_profiles.npz',SOURCE/'gene_annotations.csv',refpath]})
    (DATA/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')

    apply_style(8.5)
    fig=plt.figure(figsize=(6.3,5.7))
    ax=fig.add_axes([.155,.16,.675,.75])
    y=np.arange(len(result))
    for i,(gene,row) in enumerate(result.iterrows()):
        if gene=='KDM5C':ax.axhspan(i-.44,i+.44,color='#795082',alpha=.08,zorder=0)
        ax.scatter(row.six_hit_distance,i,s=32,color='#795082',edgecolor='white',linewidth=.4,zorder=4)
    ax.set(yticks=y,yticklabels=result.index,ylim=(15.65,-.65),xlim=(1.10,2.03),xticks=[1.2,1.4,1.6,1.8,2.0])
    ax.set_xlabel('Profile distance to centroid (smaller = closer)',labelpad=7,fontsize=8.4)
    ax.tick_params(axis='y',length=0,pad=7)
    ax.get_yticklabels()[0].set_fontweight('semibold');ax.get_yticklabels()[0].set_color('#795082')
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x',color='#ECEEF1',lw=.65,zorder=0)
    bx=fig.add_axes([.86,.16,.125,.75],sharey=ax)
    bx.set_xlim(0,1);bx.axis('off')
    for i,(gene,row) in enumerate(result.iterrows()):
        bx.text(.5,i,f'{row.six_hit_distance:.3f}',ha='center',va='center',fontsize=8.4,fontweight='semibold' if gene=='KDM5C' else 'normal')
    bx.text(.5,1.035,'Distance',transform=bx.transAxes,ha='center',fontsize=8,fontweight='semibold')
    fig.text(.155,.955,'SAGA reference: all six screen hits',fontsize=9,fontweight='semibold')
    fig.text(.155,.035,'TAF6L · TADA2B · SUPT20H · TADA1 · TAF5L · ATXN7L3',fontsize=7.7,color='#555555')
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=180)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps=[]
        for i,a in enumerate(spans):
            for b in spans[i+1:]:
                ra,rb=fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
                area=(ra&rb).get_area()
                if area>1 and area/min(ra.get_area(),rb.get_area())>.12:overlaps.append([a['text'],b['text']])
        page.get_pixmap(matrix=fitz.Matrix(1.6,1.6),colorspace=fitz.csGRAY).save(OUT/'review'/f'{STEM}_greyscale.png')
    assert not outside,outside
    assert not overlaps,overlaps
    (DATA/'render_validation.json').write_text(json.dumps(dict(font='Open Sans',dimensions_inches=[6.3,5.7],
        png_dpi=600,text_outside_page=outside,text_overlaps=overlaps,rows=16,plotted_distance='six_hit_distance',saga_reference=SIX),indent=2),encoding='utf-8')
    print(result[['five_hit_distance','six_hit_distance','rank_X_five','rank_X_six','n_models']].round(6).to_string(),flush=True)
    print(json.dumps({k:manifest[k] for k in ['same_X_ranking_with_ATXN7L3','closest_common_55','closest_mean_member','independent_validation_max_error']},indent=2),flush=True)


if __name__=='__main__':main()
