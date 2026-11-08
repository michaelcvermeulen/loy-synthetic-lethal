"""Cluster shared partner profiles for all reference X-linked gametolog genes."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pymupdf as fitz

THESIS=Path(__file__).resolve().parents[3]
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/xy_paralog_codependency'
REVIEW=OUT/'review'
STEM='ch04_supp_xy_paralog_codependency_partners_draft'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
CORE={'DDX3X','EIF1AX','RPS4X','ZFX'}


def average_linkage(distance):
    """Deterministic UPGMA on the supplied pairwise profile distances."""
    names=list(distance.index);n=len(names)
    assert np.isfinite(distance.to_numpy()).all()
    assert np.allclose(distance,distance.T)
    d={(i,j):float(distance.iloc[i,j]) for i in range(n) for j in range(i+1,n)}
    active={i:[i] for i in range(n)};sizes={i:1 for i in range(n)};merges=[]
    while len(active)>1:
        a,b=min(((i,j) for i in active for j in active if i<j),key=lambda p:(d[p],p))
        k=n+len(merges);height=d[a,b]
        for c in active:
            if c in (a,b):continue
            d[tuple(sorted((c,k)))]=(sizes[a]*d[tuple(sorted((a,c)))]+sizes[b]*d[tuple(sorted((b,c)))])/(sizes[a]+sizes[b])
        active[k]=active.pop(a)+active.pop(b);sizes[k]=sizes[a]+sizes[b]
        merges.append([a,b,height,sizes[k]])
    order=[names[i] for i in next(iter(active.values()))]
    assert len(set(order))==n and np.all(np.diff([m[2] for m in merges])>=-1e-10)
    return order,merges


def dendrogram(ax,names,order,merges,side):
    nodes={i:(order.index(g),0.) for i,g in enumerate(names)}
    for k,(a,b,h,size) in enumerate(merges):
        p1,h1=nodes[int(a)];p2,h2=nodes[int(b)]
        positions=[p1,p1,p2,p2];heights=[h1,h,h,h2]
        ax.plot(positions if side=='top' else heights,heights if side=='top' else positions,
                color='#737B84',lw=.6)
        nodes[len(names)+k]=((p1+p2)/2,h)
    limit=merges[-1][2]*1.05
    if side=='top':ax.set(xlim=(-.5,len(order)-.5),ylim=(0,limit))
    else:ax.set(xlim=(limit,0),ylim=(len(order)-.5,-.5))
    ax.axis('off')


def main():
    font=apply_style(8)
    annotation=pd.read_csv(DATA/'gene_annotations.csv').set_index('gene')
    seeds=annotation.index[annotation.chromosome.eq('X')&annotation.eligible_seed].tolist()
    assert len(seeds)==16 and annotation.loc[seeds,'chromosome'].eq('X').all()
    top=pd.read_csv(DATA/'top_10_partners.csv')
    top=top.loc[top.seed.isin(seeds)].copy()
    assert len(top)==160 and top.groupby('seed').size().eq(10).all()
    partners=sorted(top.partner.unique())
    source=np.load(DATA/'genomewide_correlations.npz')
    r=pd.DataFrame(source['r'],index=source['seeds'],columns=source['genes']).loc[seeds,partners].T
    n=pd.DataFrame(source['n'],index=source['seeds'],columns=source['genes']).loc[seeds,partners].T
    self_cells=[]
    for g in set(seeds)&set(partners):
        r.loc[g,g]=np.nan;self_cells.append(g)
    # Cluster similarities between correlation profiles, not the gene-effect scores.
    # Pairwise-complete profile distances omit masked cells without zero imputation.
    row_distance=(1-r.T.corr(min_periods=8)).clip(lower=0,upper=2)
    col_distance=(1-r.corr(min_periods=20)).clip(lower=0,upper=2)
    for distance in (row_distance,col_distance):
        for i in range(len(distance)):distance.iat[i,i]=0.
    row_order,row_tree=average_linkage(row_distance)
    col_order,col_tree=average_linkage(col_distance)
    shown=r.loc[row_order,col_order]
    counts=n.loc[row_order,col_order]
    selected=pd.DataFrame(False,index=row_order,columns=col_order)
    for q in top.itertuples():selected.loc[q.partner,q.seed]=True
    assert selected.to_numpy().sum()==160
    assert all(np.isclose(shown.loc[q.partner,q.seed],q.pearson_r) for q in top.itertuples())
    shown.to_csv(DATA/'x_member_shared_partner_correlations.csv',index_label='partner')
    counts.to_csv(DATA/'x_member_shared_partner_n_models.csv',index_label='partner')
    selected.to_csv(DATA/'x_member_top10_membership.csv',index_label='partner')
    top.to_csv(DATA/'x_member_top10_long.csv',index=False)
    for axis,order,tree,distance in [('row',row_order,row_tree,row_distance),('column',col_order,col_tree,col_distance)]:
        pd.DataFrame({'gene':order,'display_order':np.arange(1,len(order)+1)}).to_csv(DATA/f'x_member_{axis}_order.csv',index=False)
        pd.DataFrame(tree,columns=['left_cluster','right_cluster','profile_distance','size']).to_csv(DATA/f'x_member_{axis}_linkage.csv',index=False)
        distance.to_csv(DATA/f'x_member_{axis}_distances.csv')
    # Tall exploratory export retains every selected partner at readable label size.
    width=10.4;cell_height=.13;heat_height=len(partners)*cell_height;height=heat_height+2.40
    fig=plt.figure(figsize=(width,height))
    def rect(x,y,w,h):return [x/width,y/height,w/width,h/height]
    x,y,w=1.55,1.45,7.20
    ax=fig.add_axes(rect(x,y,w,heat_height))
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#C8CCD1')
    assert np.nanmax(np.abs(shown.to_numpy()))<.60
    im=ax.imshow(shown,cmap=cmap,vmin=-.60,vmax=.60,aspect='auto',interpolation='none')
    ax.set_xticks(range(len(col_order)),col_order,rotation=90,fontsize=8.5)
    ax.set_yticks(range(len(row_order)),row_order,fontsize=7.2)
    ax.yaxis.tick_right();ax.tick_params(length=0,pad=4,labelleft=False,labelright=True)
    for label in ax.get_xticklabels()+ax.get_yticklabels():
        if label.get_text() in CORE:label.set_fontweight('semibold')
    for spine in ax.spines.values():spine.set_visible(False)
    yy,xx=np.where(selected.to_numpy())
    colours=['white' if abs(shown.iloc[i,j])>.37 else '#252525' for i,j in zip(yy,xx)]
    ax.scatter(xx,yy,s=4,c=colours,marker='o',linewidths=0)
    topax=fig.add_axes(rect(x,y+heat_height+.06,w,.60))
    leftax=fig.add_axes(rect(.08,y,1.36,heat_height))
    dendrogram(topax,list(col_distance.index),col_order,col_tree,'top')
    dendrogram(leftax,list(row_distance.index),row_order,row_tree,'left')
    ax.set_xlabel('X-linked member of an X–Y pair',fontsize=9,labelpad=7)
    cb=fig.colorbar(im,cax=fig.add_axes(rect(1.55,.30,2.65,.10)),orientation='horizontal',ticks=[-.6,-.3,0,.3,.6])
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=8,length=2,pad=2)
    cb.ax.xaxis.set_label_position('top')
    cb.set_label('Pearson r',fontsize=8,labelpad=1)
    fig.legend(handles=[Line2D([],[],marker='o',color='#252525',ls='',ms=3,label='Top 10 for that column'),
                        Patch(facecolor='#C8CCD1',label='Self or <200 shared models')],
               loc='lower left',bbox_to_anchor=(4.60/width,.13/height),frameon=False,fontsize=8,ncol=1)
    fig.text(4.7/width,.07/height,'Bold labels: DDX3X, EIF1AX, RPS4X and ZFX.',fontsize=7.6,color='#555555')
    fig.canvas.draw()
    assert np.allclose(ax.images[0].get_array(),shown.to_numpy(),equal_nan=True)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(REVIEW/f'{STEM}_preview.png',dpi=120)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all(g in page.get_text() for g in seeds+partners)
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for s in line['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        assert not outside,outside
        page.get_pixmap(matrix=fitz.Matrix(.9,.9),colorspace=fitz.csGRAY).save(REVIEW/f'{STEM}_greyscale.png')
    record=dict(stem=STEM,font=font,width_inches=width,height_inches=height,seeds=seeds,clustered_column_order=col_order,
                selected_pairs=len(top),unique_partner_rows=len(partners),clustered_row_order=row_order,
                correlation_scale=[-.6,.6],no_colour_clipping=True,seed_model_range=[int(annotation.loc[seeds,'n_models'].min()),int(annotation.loc[seeds,'n_models'].max())],
                displayed_pair_model_range=[int(counts.where(shown.notna()).min().min()),int(counts.where(shown.notna()).max().max())],
                missing_cells=int(shown.isna().to_numpy().sum()),self_masked=sorted(self_cells),
                clustering='Both axes: average linkage on 1 minus Pearson correlation between their unscaled correlation profiles; pairwise-complete distances, no zero imputation; deterministic tie/leaf order; no discrete cluster cut',
                profile_distance_minimum_shared_features=dict(rows=8,columns=20),
                interpretation='Descriptive partner-selected profiles; neither a LOY-specificity test nor evidence of direct physical interaction. Strong partners selected in one column are shown in every column, including weak and negative correlations.',
                source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [DATA/'genomewide_correlations.npz',DATA/'top_10_partners.csv',DATA/'gene_annotations.csv']},
                text_outside_page=outside,analysis_values_preserved=True)
    (DATA/'render_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k not in ['clustered_row_order','source_sha256']},indent=2))


if __name__=='__main__':main()
