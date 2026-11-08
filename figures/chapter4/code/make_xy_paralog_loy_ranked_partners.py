"""Show genome-wide top-ten partners selected and ranked within LOY males."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
import pymupdf as fitz

THESIS=Path(__file__).resolve().parents[3]
OUT=THESIS/'figures/chapter4'
BASE=OUT/'source_data/xy_paralog_codependency'
DATA=OUT/'source_data/xy_paralog_loy_ranked_partners'
STEM='ch04_supp_xy_paralog_loy_ranked_partners_draft'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

def main():
    apply_style(8)
    top=pd.read_csv(DATA/'loy_top10_partners.csv')
    assert len(top)==160 and top.notna().all().all()
    annotation=pd.read_csv(BASE/'gene_annotations.csv').set_index('gene')
    screen=pd.read_csv(OUT/'source_data/ht1080_depmap_codependency/seed_annotations.csv').fillna('').set_index('gene')
    core=['DDX3X','EIF1AX','RPS4X','ZFX']
    seeds=core+sorted(set(top.seed)-set(core))
    assert len(seeds)==16 and annotation.loc[seeds,'chromosome'].eq('X').all()
    tables=[]
    for i,g in enumerate(seeds):
        tab=top.loc[top.seed.eq(g)].sort_values('rank_positive').copy()
        assert len(tab)==10 and tab.rank_positive.tolist()==list(range(1,11))
        assert tab.loy_pearson_r.gt(0).all() and tab.loy_n_models.ge(100).all()
        assert tab.partner.ne(g).all()
        tab.insert(0,'seed_display_order',i+1)
        tab['partner_screen_direction']=tab.partner.map(screen.screen_direction).fillna('Not nominated')
        tab['partner_saga_reference']=tab.partner.map(screen.saga_reference_20).fillna(False).astype(bool)
        tables.append(tab)
    displayed=pd.concat(tables,ignore_index=True)
    displayed.to_csv(DATA/'ranked_partners_display.csv',index=False)

    width,height=8.2,8.3
    fig=plt.figure(figsize=(width,height))
    def rect(x,y,w,h):return [x/width,y/height,w/width,h/height]
    cmap=plt.get_cmap('RdBu_r')
    norm=Normalize(-1,1)
    axes=[]
    for block,y in enumerate([4.70,.95]):
        ax=fig.add_axes(rect(.57,y,7.48,3.25))
        selected=tables[block*8:(block+1)*8]
        ax.set_xlim(-.5,7.5);ax.set_ylim(9.5,-.5)
        ax.set_xticks(range(8),seeds[block*8:(block+1)*8],fontsize=9,fontweight='semibold')
        ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=7)
        ax.set_yticks(range(10),range(1,11),fontsize=8)
        ax.tick_params(axis='y',length=0,pad=4)
        ax.set_ylabel('LOY partner rank',fontsize=9,labelpad=5)
        for spine in ax.spines.values():spine.set_visible(False)
        for col,t in enumerate(selected):
            for row,p in enumerate(t.itertuples()):
                ax.add_patch(Rectangle((col-.5,row-.5),1,.5,facecolor='#f0f0f0',edgecolor='none'))
                name=p.partner+({'Depleted':'*','Enriched':'**'}.get(p.partner_screen_direction,''))
                ax.text(col,row-.25,name,ha='center',va='center',fontsize=8.5,
                        fontweight='semibold' if p.partner_saga_reference else 'normal',color='#222222')
                for group,offset in [('xy',-.5),('loy',0)]:
                    r=getattr(p,f'{group}_pearson_r')
                    rgb=np.array(cmap(norm(r))[:3])
                    linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
                    colour='#222222' if float(linear@np.array([.2126,.7152,.0722]))>.179 else 'white'
                    ax.add_patch(Rectangle((col+offset,row),.5,.5,facecolor=cmap(norm(r)),edgecolor='white',linewidth=.35))
                    value=0.0 if abs(r)<.005 else r
                    ax.text(col+offset+.25,row+.25,f'{value:.2f}',ha='center',va='center',fontsize=8,color=colour)
                ax.add_patch(Rectangle((col-.5,row-.5),1,1,facecolor='none',edgecolor='white',linewidth=.8))
        axes.append(ax)
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes(rect(2.65,.71,3.85,.10)),orientation='horizontal',ticks=[-1,-.5,0,.5,1])
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=8,length=2,pad=2)
    cb.set_label('Pearson r · within male Y-status group, unadjusted',fontsize=8,labelpad=2)
    fig.text(.57/width,.205/height,'Left: XY (Y-retaining); right: LOY. Partners selected and ranked within LOY males.',fontsize=8,color='#555555')
    fig.text(.57/width,.045/height,'* Depleted; ** enriched HT1080 hit; semibold partner: SAGA member.',fontsize=8,color='#555555')
    fig.canvas.draw()
    for block,ax in enumerate(axes):
        assert len(ax.patches)==8*10*4
        for col,t in enumerate(tables[8*block:8*(block+1)]):
            for row,p in enumerate(t.itertuples()):
                for band,group in [(1,'xy'),(2,'loy')]:
                    patch=ax.patches[(col*10+row)*4+band]
                    assert np.allclose(patch.get_facecolor(),cmap(norm(getattr(p,f'{group}_pearson_r'))))
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=160)
    plt.close(fig)

    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all(g in page.get_text() for g in set(seeds+displayed.partner.tolist()))
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps=[]
        for i,a in enumerate(spans):
            for b in spans[i+1:]:
                ra,rb=fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
                area=(ra&rb).get_area()
                if area>1 and area/min(ra.get_area(),rb.get_area())>.10:overlaps.append([a['text'],b['text']])
        assert not outside,outside
        assert not overlaps,overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.5,1.5),colorspace=fitz.csGRAY).save(OUT/'review'/f'{STEM}_greyscale.png')
    sources=[DATA/'loy_top10_partners.csv',DATA/'selected_partners_by_group_long.csv',DATA/'analysis_manifest.json',BASE/'gene_annotations.csv',OUT/'source_data/ht1080_depmap_codependency/seed_annotations.csv']
    record=dict(stem=STEM,seeds=seeds,seed_count=16,displayed_pairs=len(displayed),partners_per_gene=10,
                order='Four established gametolog dependencies first, followed by the other twelve genes alphabetically; no clustering.',
                layout='Two blocks of eight columns; rows are independent within-LOY partner ranks. Each cell: partner name above side-by-side r values, XY male left and LOY male right. Counts in source table and caption.',
                dimensions_inches=[width,height],font='Open Sans',pearson_scale=[-1,1],correlation='Unadjusted pairwise-complete Pearson, separately in XY and LOY males',
                model_count_range={s:[int(displayed[f'{s}_n_models'].min()),int(displayed[f'{s}_n_models'].max())] for s in ['xy','loy']},
                selection_group="LOY",min_paired_loy_models=100,xy_does_not_affect_selection=True,
                new_analysis=False,analysis='Genome-wide within-LOY positive Pearson ranking; XY values shown for the same partners. No significance tests.',text_outside_page=outside,text_overlaps=overlaps,
                source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (DATA/'ranked_partners_render_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':main()
