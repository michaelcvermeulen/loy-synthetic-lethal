"""Thesis-style exploratory UMAP and offline gene-search viewer."""
from pathlib import Path
import base64,hashlib,json,sys,os
THESIS=Path(__file__).resolve().parents[3]
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/ht1080_hit_profile_umap'
STEM='ch04_supp_ht1080_hit_profile_umap_draft'
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
import pymupdf as fitz
from PIL import Image
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
COL={'Depleted':'#795082','Enriched':'#00877C','SAGA':'#A88F25','X':'#426EA0','Y':'#C27B38','Other':'#CBD0D5'}


def base(ax, e):
    ax.scatter(e.UMAP1,e.UMAP2,s=2,c=COL['Other'],alpha=.30,lw=0,rasterized=False,zorder=1)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.set(xticks=[],yticks=[],xlabel='UMAP 1',ylabel='UMAP 2')
    ax.set_aspect('equal',adjustable='box')
    lo=e.min();hi=e.max();margin=(hi-lo)*.09
    ax.set_xlim(lo.UMAP1-margin.UMAP1,hi.UMAP1+margin.UMAP1)
    ax.set_ylim(lo.UMAP2-margin.UMAP2,hi.UMAP2+margin.UMAP2)
    ax.xaxis.labelpad=1;ax.yaxis.labelpad=1


def labels(ax,e,genes,colours=None,size=7.2):
    genes=[g for g in genes if g in e.index]
    np.random.seed(42)
    texts=[ax.text(e.loc[g,'UMAP1'],e.loc[g,'UMAP2'],g,fontsize=size,
                   color=colours.get(g,'#333333') if colours else '#333333',
                   zorder=6,bbox=dict(facecolor='white',edgecolor='none',alpha=.65,pad=.25)) for g in genes]
    adjust_text(texts, x=e.loc[genes,'UMAP1'].to_numpy(),y=e.loc[genes,'UMAP2'].to_numpy(),
                ax=ax,expand=(1.12,1.28),force_text=(.35,.5),force_static=(.15,.3),
                force_pull=(.008,.008),iter_lim=700,ensure_inside_axes=True,
                arrowprops=dict(arrowstyle='-',lw=.4,color='#777777',alpha=.7))
    return texts


def legend(ax,items):
    handles=[Line2D([],[],marker=marker,ls='',mfc=COL[key],mec='white',mew=.4,ms=6,label=label)
             for key,marker,label in items]
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(0,1.065),ncol=len(items),
              fontsize=7.6,handletextpad=.3,columnspacing=1.0,borderaxespad=0)


def static():
    apply_style(8.3)
    ann=pd.read_csv(DATA/'gene_annotations.csv',index_col=0).fillna('')
    e=pd.read_csv(DATA/'embedding_primary.csv',index_col=0)
    assert ann.index.equals(e.index)
    fig=plt.figure(figsize=(12.4,6.6))
    axes=[fig.add_axes([.045,.105,.455,.79]),fig.add_axes([.54,.105,.455,.79])]
    label_texts=[]
    for ax in axes:base(ax,e)
    ax=axes[0]
    for direction in ['Depleted','Enriched']:
        sub=e.loc[ann.screen_direction.eq(direction)]
        ax.scatter(sub.UMAP1,sub.UMAP2,s=28,c=COL[direction],ec='white',lw=.45,zorder=4)
    hits=ann.loc[ann.screen_direction.ne('Not nominated')]
    top=hits.loc[hits.screen_direction.eq('Depleted')].sort_values('candidate_normZ').head(10).index.tolist()
    top+=hits.loc[hits.saga_reference_20.astype(bool)].index.tolist()
    top+=hits.loc[hits.screen_direction.eq('Enriched')].index.tolist()
    top=list(dict.fromkeys(top))
    label_texts+=labels(ax,e,top,{g:COL[hits.loc[g,'screen_direction']] for g in top})
    legend(ax,[('Depleted','o','LOY-depleted (47)'),('Enriched','o','LOY-enriched (9)')])
    ax=axes[1]
    for key,mask,marker in [('SAGA',ann.saga_reference_20.astype(bool),'o'),
                            ('X',ann.gametolog.eq('X'),'D'),('Y',ann.gametolog.eq('Y'),'^')]:
        sub=e.loc[mask]
        ax.scatter(sub.UMAP1,sub.UMAP2,s=32,c=COL[key],marker=marker,ec='white',lw=.45,zorder=4)
    focus=['KDM5C','KDM5D','KDM6A','UTY','DDX3X','DDX3Y','EIF1AX','EIF1AY','RPS4X','RPS4Y1','ZFX','ZFY',
           'TAF6L','TADA2B','SUPT20H','TADA1','TAF5L','ATXN7L3','SUPT7L','USP22']
    colours={g:COL[ann.loc[g,'gametolog']] if ann.loc[g,'gametolog'] in ['X','Y'] else '#806B16' for g in focus if g in ann.index}
    label_texts+=labels(ax,e,focus,colours)
    legend(ax,[('SAGA','o','SAGA'),('X','D','X-linked gametolog'),('Y','^','Y-linked gametolog')])
    fig.text(.012,.972,'a',fontsize=13,fontweight='semibold',va='top')
    fig.text(.515,.972,'b',fontsize=13,fontweight='semibold',va='top')
    fig.text(.05,.965,'HT1080 candidates',fontsize=9,fontweight='semibold',va='top')
    fig.text(.55,.965,'SAGA and X–Y gametologs',fontsize=9,fontweight='semibold',va='top')
    fig.text(.5,.022,'17,931 genes · 56 signed co-dependency correlations per gene · shared UMAP coordinates',
             fontsize=8.2,color='#555555',ha='center')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    overlaps=[]
    for i,a in enumerate(label_texts):
        for b in label_texts[i+1:]:
            if a.axes is b.axes and a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)):
                overlaps.append([a.get_text(),b.get_text()])
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{ext}',dpi=600)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=165)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        page.get_pixmap(matrix=fitz.Matrix(1.4,1.4),colorspace=fitz.csGRAY).save(OUT/'review'/f'{STEM}_greyscale.png')
    # Machado et al. severe deuteranopia approximation, applied in linear RGB.
    im=np.asarray(Image.open(OUT/'review'/f'{STEM}_preview.png').convert('RGB')).astype(float)/255
    linear=np.where(im<=.04045,im/12.92,((im+.055)/1.055)**2.4)
    mat=np.array([[.367322,.860646,-.227968],[.280085,.672501,.047413],[-.011820,.042940,.968881]])
    sim=np.clip(linear@mat.T,0,1)
    sim=np.where(sim<=.0031308,sim*12.92,1.055*sim**(1/2.4)-.055)
    Image.fromarray(np.uint8(np.clip(sim,0,1)*255)).save(OUT/'review'/f'{STEM}_deuteranopia.png')
    record=dict(stem=STEM,font='Open Sans',dimensions_inches=[12.4,6.6],png_dpi=600,
                input_genes=len(e),panels_same_coordinates=True,
                panel_a_labels=top,panel_b_labels=[g for g in focus if g in e.index],
                label_overlaps=overlaps,text_outside_page=outside,
                highlighted_saga=int(ann.saga_reference_20.astype(bool).sum()),
                highlighted_X=int(ann.gametolog.eq('X').sum()),highlighted_Y=int(ann.gametolog.eq('Y').sum()))
    (DATA/'render_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2),flush=True)
    assert not outside,outside
    assert not overlaps,overlaps


def sensitivity():
    apply_style(8)
    ann=pd.read_csv(DATA/'gene_annotations.csv',index_col=0).fillna('')
    from sklearn.neighbors import NearestNeighbors
    knn=np.load(DATA/'exact_knn.npz')['indices'][:,1:16]
    fig,axes=plt.subplots(2,2,figsize=(10.5,9.6))
    summaries=[]
    for ax,(name,title) in zip(axes.flat,[('primary','30 neighbours · seed 42'),('neighbours15','15 neighbours · seed 42'),
                                        ('neighbours50','50 neighbours · seed 42'),('seed43','30 neighbours · seed 43')]):
        e=pd.read_csv(DATA/f'embedding_{name}.csv',index_col=0)
        base(ax,e)
        for direction in ['Depleted','Enriched']:
            sub=e.loc[ann.screen_direction.eq(direction)]
            ax.scatter(sub.UMAP1,sub.UMAP2,s=19,c=COL[direction],ec='white',lw=.3,zorder=4)
        labels(ax,e,['KDM5C','KDM6A','DDX3X','TADA2B','ATXN7L3','USP7','SREBF1'],size=7)
        ax.set_title(title,loc='left',fontsize=9,pad=5)
        emb_knn=NearestNeighbors(n_neighbors=16).fit(e.to_numpy()).kneighbors(e.to_numpy(),return_distance=False)[:,1:]
        recall=np.array([len(set(a)&set(b))/15 for a,b in zip(knn,emb_knn)])
        summaries.append(dict(fit=name,median_15nn_recall=float(np.median(recall)),mean_15nn_recall=float(recall.mean()),
                              candidate_mean_recall=float(recall[ann.screen_direction.ne('Not nominated')].mean())))
    fig.subplots_adjust(left=.04,right=.99,top=.965,bottom=.055,hspace=.14,wspace=.06)
    fig.savefig(OUT/'review'/f'{STEM}_settings.pdf')
    fig.savefig(OUT/'review'/f'{STEM}_settings.png',dpi=150)
    plt.close(fig)
    pd.DataFrame(summaries).to_csv(DATA/'embedding_neighbour_preservation.csv',index=False)
    print('Neighbour-preservation diagnostics:',json.dumps(summaries),flush=True)


def interactive():
    z=np.load(DATA/'input_profiles.npz')
    e=pd.read_csv(DATA/'embedding_primary.csv',index_col=0)
    a=pd.read_csv(DATA/'gene_annotations.csv',index_col=0).fillna('')
    k=np.load(DATA/'exact_knn.npz')
    matrix=np.round(z['distance_r'],4)
    payload=dict(genes=e.index.tolist(),xy=np.round(e.to_numpy(),5).tolist(),features=z['features'].tolist(),
                 profiles=[[None if not np.isfinite(x) else float(x) for x in row] for row in matrix],
                 directions=a.screen_direction.tolist(),saga=a.saga_reference_20.astype(bool).astype(int).tolist(),
                 gametolog=a.gametolog.tolist(),n=a.n_models.tolist(),rms=np.round(a.profile_rms,4).tolist(),
                 neighbours=k['indices'][:,1:11].tolist(),distances=np.round(k['distances'][:,1:11],4).tolist())
    saga_distance=pd.read_csv(OUT/'source_data/xy_saga_profile_distance/x_gametolog_saga_distances.csv')
    payload['saga_distance']=saga_distance[['gene','five_hit_distance','six_hit_distance','rank_X_five','rank_X_six']].to_dict('records')
    template=(Path(__file__).with_name('ht1080_hit_profile_umap_viewer.html.in')).read_text(encoding='utf-8')
    font=(THESIS/'skills/thesis-figure-style/assets/fonts/OpenSans-Regular.ttf').read_bytes()
    template=template.replace('FONT_PAYLOAD',base64.b64encode(font).decode('ascii'))
    html=template.replace('/*DATA_PAYLOAD*/',json.dumps(payload,separators=(',',':'),allow_nan=False))
    (OUT/f'{STEM}.html').write_text(html,encoding='utf-8')
    # Retain the previously opened path as a fully populated offline viewer.
    # The non-browsable .html.in source can no longer be confused with this entry.
    Path(__file__).with_name('ht1080_hit_profile_umap_viewer.html').write_text(html,encoding='utf-8')


if __name__=='__main__':
    static()
    sensitivity()
    interactive()
