"""Three-panel draft linking X-gametolog co-dependencies to HT1080 drugZ hits."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle,Patch
from matplotlib.lines import Line2D
import pymupdf as fitz

THESIS=Path(__file__).resolve().parents[3]
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/xy_screen_overlap'
STEM='ch04_supp_xy_paralog_screen_overlap_draft'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
W,H=12.4,11.1
HIT='#7A3E85';GREY='#777F89';MISSING='#C8CCD1';GOLD='#BAA04B'

def rect(x,y,w,h):return [x/W,y/H,w/W,h/H]
def fdr_label(value):
    return f'{value:.1e}' if value<.001 else f'{value:.2g}' if value<.01 else f'{value:.3f}'.rstrip('0').rstrip('.')
def clean(ax):
    for side in ['top','right','left']:ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_linewidth(.6)
    ax.tick_params(axis='y',length=0)
    ax.tick_params(axis='x',labelsize=7.2,length=2,width=.6)
def marker(ax,x,y,hit,fdr,size=19):
    edge=HIT if hit else GREY
    ax.scatter(x,y,s=size,facecolor=edge if (not hit or fdr<.05) else 'white',edgecolor=edge,lw=.8,zorder=4)

def main():
    apply_style(8)
    counts=pd.read_csv(DATA/'panel_a_overlap_counts.csv',index_col=0)
    r=pd.read_csv(DATA/'panel_b_correlations.csv',index_col=0)
    selected=pd.read_csv(DATA/'panel_b_top25_membership.csv',index_col=0)
    scores=pd.read_csv(DATA/'panel_b_screen_scores.csv',index_col=0)
    detail=pd.read_csv(DATA/'panel_c_kdm5c_top25.csv')
    manifest=json.loads((DATA/'analysis_manifest.json').read_text())
    seeds=counts.index.tolist();genes=scores.index.tolist()
    assert r.columns.tolist()==seeds and r.index.tolist()==genes
    fig=plt.figure(figsize=(W,H))
    fig.text(.12/W,10.98/H,'a',fontsize=14,fontweight='semibold',va='top')
    fig.text(4.72/W,10.98/H,'b',fontsize=14,fontweight='semibold',va='top')
    fig.text(.12/W,7.05/H,'c',fontsize=14,fontweight='semibold',va='top')

    # a: exact hit counts, without an enrichment or significance interpretation.
    ax=fig.add_axes(rect(1.2,7.88,2.8,2.65))
    cmap_a=LinearSegmentedColormap.from_list('hit_count',['#FAF8FA','#B18EB7',HIT])
    ax.imshow(counts,cmap=cmap_a,vmin=0,vmax=5,aspect='auto',interpolation='none')
    ax.set_xticks(range(4),['10','25','50','100'],fontsize=8)
    ax.xaxis.tick_top();ax.xaxis.set_label_position('top')
    ax.set_xlabel('Top co-dependency partners',fontsize=8.5,labelpad=7)
    ax.set_yticks(range(16),seeds,fontsize=7.5)
    ax.tick_params(length=0,pad=4)
    ax.set_xticks(np.arange(-.5,4,1),minor=True)
    ax.set_yticks(np.arange(-.5,16,1),minor=True)
    ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',length=0)
    for spine in ax.spines.values():spine.set_visible(False)
    for i in range(16):
        for j in range(4):
            value=int(counts.iloc[i,j])
            ax.text(j,i,str(value),ha='center',va='center',fontsize=8,
                    fontweight='semibold' if value else 'normal',color='white' if value>=3 else '#777777' if value==0 else '#333333')
    ax.add_patch(Rectangle((-.49,-.49),3.98,.98,fill=False,ec=HIT,lw=1,clip_on=False))
    ax.get_yticklabels()[0].set_fontweight('semibold');ax.get_yticklabels()[0].set_color(HIT)
    fig.text(.55/W,7.65/H,'Cells: number of depleted HT1080 hits',fontsize=8,color='#555555')
    fig.text(.55/W,7.42/H,'Purple: depleted hit; bold partner: SAGA member',fontsize=7.7,color='#555555')
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=HIT,ms=4,label='FDR <0.05'),
                        Line2D([],[],marker='o',ls='',color=HIT,mfc='white',ms=4,label='FDR 0.05–<0.25')],
               loc='lower left',bbox_to_anchor=(.46/W,7.12/H),ncol=2,frameon=False,fontsize=7.5,columnspacing=1.2,handletextpad=.45)

    # b: all 48 depleted candidates, raw r and original screen scores.
    bx=fig.add_axes(rect(6.0,1.4,4.0,8.7))
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad(MISSING)
    assert np.nanmax(np.abs(r.to_numpy()))<.6
    im=bx.imshow(r,cmap=cmap,vmin=-.6,vmax=.6,aspect='auto',interpolation='none')
    bx.set_xticks(range(16),seeds,rotation=90,fontsize=7.3)
    bx.xaxis.tick_top();bx.tick_params(axis='x',length=0,pad=5)
    bx.set_yticks(range(48),genes,fontsize=7.1)
    bx.tick_params(axis='y',length=0,pad=11)
    for spine in bx.spines.values():spine.set_visible(False)
    bx.get_xticklabels()[0].set_fontweight('semibold');bx.get_xticklabels()[0].set_color(HIT)
    for i,label in enumerate(bx.get_yticklabels()):
        if scores.iloc[i].saga_reference:
            label.set_fontweight('semibold')
            bx.add_patch(Rectangle((-.87,i-.34),.22,.68,fc=GOLD,ec='none',clip_on=False))
        if genes[i]=='ZCCHC16':label.set_color('#777777')
    bx.add_patch(Rectangle((-.49,-.49),.98,47.98,fill=False,ec=HIT,lw=1.1,clip_on=False))
    yy,xx=np.where(selected.to_numpy())
    bx.scatter(xx,yy,s=6,c=['white' if abs(r.iloc[i,j])>.37 else '#252525' for i,j in zip(yy,xx)],lw=0,zorder=4)
    sx=fig.add_axes(rect(10.25,1.4,1.02,8.7))
    sx.set(xlim=(-12,1),ylim=(47.5,-.5),xticks=[-10,-5,0],yticks=[])
    sx.axvline(0,color='#A8ADB4',lw=.6)
    for i,q in enumerate(scores.itertuples()):
        sx.hlines(i,0,q.normZ,color=HIT,alpha=.35,lw=.65)
        marker(sx,q.normZ,i,True,q.fdr_synth,size=14)
    clean(sx);sx.set_xlabel('drugZ normZ',fontsize=8,labelpad=5)
    fig.text(10.76/W,10.28/H,'HT1080 screen',fontsize=8,ha='center')
    fx=fig.add_axes(rect(11.42,1.4,.72,8.7))
    fx.set(xlim=(0,1),ylim=(47.5,-.5));fx.axis('off')
    for i,q in enumerate(scores.itertuples()):fx.text(.5,i,fdr_label(q.fdr_synth),ha='center',va='center',fontsize=6.8,color=HIT if q.fdr_synth<.05 else '#555555')
    fig.text(11.78/W,10.28/H,'FDR',fontsize=8,ha='center')
    cb=fig.colorbar(im,cax=fig.add_axes(rect(6.0,.72,3.0,.12)),orientation='horizontal',ticks=[-.6,-.3,0,.3,.6])
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=7.4,length=2,pad=2)
    cb.set_label('DepMap Pearson r · unadjusted',fontsize=8,labelpad=3)
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color='#252525',ms=3,label="In that column's top 25"),
                        Patch(facecolor=GOLD,label='SAGA reference member'),Patch(facecolor=MISSING,label='No profile / self')],
               loc='lower left',bbox_to_anchor=(9.15/W,.31/H),frameon=False,fontsize=7.3,labelspacing=.45,handlelength=1)
    fig.text(10.15/W,.17/H,'Negative: depleted in LOY',fontsize=7.4,color='#555555')

    # c: every KDM5C top-25 partner, including screen non-hits.
    cx=fig.add_axes(rect(1.38,.85,.98,5.9))
    dx=fig.add_axes(rect(2.70,.85,1.12,5.9))
    ex=fig.add_axes(rect(4.02,.85,.52,5.9))
    labels=[f'{q.partner_rank}  {q.partner}' for q in detail.itertuples()]
    cx.set(xlim=(0,.47),ylim=(24.5,-.5),yticks=range(25),yticklabels=labels,xticks=[0,.2,.4])
    cx.tick_params(axis='y',labelsize=7.2,pad=3)
    dx.set(xlim=(-5.5,2.4),ylim=(24.5,-.5),yticks=[],xticks=[-5,0,2])
    dx.axvline(0,color='#A8ADB4',lw=.6)
    ex.set(xlim=(0,1),ylim=(24.5,-.5));ex.axis('off')
    for i,q in enumerate(detail.itertuples()):
        if q.depleted_hit:
            for target in (cx,dx,ex):target.axhspan(i-.5,i+.5,color='#F4EDF5',zorder=0)
        colour=HIT if q.depleted_hit else GREY
        cx.hlines(i,0,q.pearson_r,color=colour,alpha=.5,lw=.7)
        marker(cx,q.pearson_r,i,q.depleted_hit,q.fdr_synth,size=17)
        cx.text(q.pearson_r+.018,i,f'{q.pearson_r:.2f}',va='center',fontsize=6.8,color=colour)
        dx.hlines(i,0,q.normZ,color=colour,alpha=.5,lw=.7)
        marker(dx,q.normZ,i,q.depleted_hit,q.fdr_synth,size=17)
        ex.text(.5,i,fdr_label(q.fdr_synth),ha='center',va='center',fontsize=6.8,color=colour)
        label=cx.get_yticklabels()[i]
        label.set_color(colour if q.depleted_hit else '#333333')
        if q.saga_reference:label.set_fontweight('semibold')
    for target in (cx,dx):clean(target)
    cx.set_xlabel('Pearson r',fontsize=8,labelpad=4)
    dx.set_xlabel('drugZ normZ',fontsize=8,labelpad=4)
    fig.text(1.38/W,6.93/H,'KDM5C co-dependency',fontsize=8)
    fig.text(2.70/W,6.93/H,'HT1080 screen',fontsize=8)
    fig.text(4.28/W,6.93/H,'FDR',fontsize=8,ha='center')
    fig.text(.35/W,.20/H,'c: all 25 partners, ranked by KDM5C co-dependency',fontsize=7.7,color='#555555')

    fig.canvas.draw()
    assert np.allclose(bx.images[0].get_array(),r.to_numpy(),equal_nan=True)
    assert np.allclose(ax.images[0].get_array(),counts.to_numpy())
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=140)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all(g in page.get_text() for g in seeds+genes+detail.partner.tolist())
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps=[]
        for i,a in enumerate(spans):
            for b in spans[i+1:]:
                ra,rb=fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
                area=(ra&rb).get_area()
                if area>1 and area/min(ra.get_area(),rb.get_area())>.12:overlaps.append([a['text'],b['text'],area])
        page.get_pixmap(matrix=fitz.Matrix(1,1),colorspace=fitz.csGRAY).save(OUT/'review'/f'{STEM}_greyscale.png')
    files=list(DATA.glob('panel_*.csv'))
    record=dict(stem=STEM,dimensions_inches=[W,H],font='Open Sans',png_dpi=600,panels=['a: overlap counts','b: all 48 depleted hits across 16 X-linked genes','c: all KDM5C top-25 partners'],
                original_values_preserved=True,correlation_scale=[-.6,.6],colour_clipping=False,missing_correlations=int(r.isna().to_numpy().sum()),
                top25_membership_dots=int(selected.to_numpy().sum()),text_outside_page=outside,text_overlaps=overlaps,
                source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (DATA/'render_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k!='source_sha256'},indent=2))
    assert not outside,outside
    assert not overlaps,overlaps

if __name__=='__main__':main()
