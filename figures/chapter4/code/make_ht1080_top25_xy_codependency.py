"""Top-25 HT1080 depleted hits versus six X-linked and two Y-linked profiles."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,Rectangle
import pymupdf as fitz

THESIS=Path(__file__).resolve().parents[3]
OUT=THESIS/'figures/chapter4'
SOURCE=OUT/'source_data/xy_screen_overlap'
DATA=OUT/'source_data/xy_screen_top25_focus'
STEM='ch04_supp_ht1080_top25_xy_codependency_draft'
X_SEEDS=['KDM5C','KDM6A','DDX3X','EIF1AX','RPS4X','ZFX']
SEEDS=['KDM5C','KDM5D','KDM6A','UTY','DDX3X','EIF1AX','RPS4X','ZFX']
XY_SOURCE=OUT/'source_data/xy_paralog_codependency'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

def main():
    DATA.mkdir(exist_ok=True)
    apply_style(8)
    scores=pd.read_csv(SOURCE/'panel_b_screen_scores.csv',index_col=0).sort_values('rank_synth').head(25)
    assert scores.rank_synth.tolist()==list(range(1,26)) and scores.fdr_synth.lt(.05).all()
    r=pd.read_csv(SOURCE/'panel_b_correlations.csv',index_col=0).loc[scores.index,X_SEEDS].copy()
    n=pd.read_csv(SOURCE/'panel_b_model_counts.csv',index_col=0).loc[scores.index,X_SEEDS].copy()
    frozen=np.load(XY_SOURCE/'genomewide_correlations.npz')
    for g in ['KDM5D','UTY']:
        i=frozen['seeds'].tolist().index(g)
        r[g]=pd.Series(frozen['r'][i],index=frozen['genes']).reindex(scores.index)
        n[g]=pd.Series(frozen['n'][i],index=frozen['genes']).reindex(scores.index)
    r=r[SEEDS];n=n[SEEDS]
    assert r.isna().to_numpy().sum()==1 and pd.isna(r.loc['DDX3X','DDX3X'])
    assert n[X_SEEDS].where(r[X_SEEDS].notna()).stack().dropna().eq(1186).all()
    assert n[['KDM5D','UTY']].eq(355).all().all()
    assert np.nanmax(np.abs(r.to_numpy()))<.6
    assert 'TAF5L' not in scores.index
    scores.to_csv(DATA/'screen_scores.csv')
    r.to_csv(DATA/'correlations.csv',index_label='screen_gene')
    n.to_csv(DATA/'model_counts.csv',index_label='screen_gene')
    cohort=pd.read_csv(XY_SOURCE/'model_cohort.csv').set_index('ModelID')
    effects=pd.read_csv(XY_SOURCE/'seed_gene_effect.csv',index_col=0)
    coverage=[]
    for g in SEEDS:
        available=effects[g].notna()
        sexes=cohort.loc[available,'Sex'].fillna('Unknown').value_counts()
        coverage.append(dict(gene=g,chromosome='Y' if g in ['KDM5D','UTY'] else 'X',
                             n_models=int(available.sum()),Male=int(sexes.get('Male',0)),
                             Female=int(sexes.get('Female',0)),Unknown=int(sexes.get('Unknown',0))))
    coverage=pd.DataFrame(coverage).set_index('gene')
    coverage.to_csv(DATA/'column_coverage.csv')
    assert effects.KDM5D.notna().equals(effects.UTY.notna())

    width,height=9.15,8.05
    fig=plt.figure(figsize=(width,height))
    def rect(x,y,w,h):return [x/width,y/height,w/width,h/height]
    ax=fig.add_axes(rect(1.45,1.10,5.60,6.25))
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#C8CCD1')
    im=ax.imshow(r,cmap=cmap,vmin=-.6,vmax=.6,aspect='auto',interpolation='none')
    ax.set_xticks(range(8),SEEDS,fontsize=8.8,fontweight='semibold')
    ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=19)
    for j,g in enumerate(SEEDS):
        ax.text(j,1.015,f'n={int(coverage.loc[g,"n_models"]):,}',transform=ax.get_xaxis_transform(),
                ha='center',va='bottom',fontsize=7,color='#666666')
    ax.get_xticklabels()[0].set_color('#7A3E85')
    labels=[f'{int(q.rank_synth)}  {q.Index}' for q in scores.itertuples()]
    ax.set_yticks(range(25),labels,fontsize=8.2)
    ax.tick_params(axis='y',length=0,pad=12)
    ax.set_xticks(np.arange(-.5,8,1),minor=True)
    ax.set_yticks(np.arange(-.5,25,1),minor=True)
    ax.grid(which='minor',color='white',linewidth=.6)
    ax.tick_params(which='minor',length=0)
    for spine in ax.spines.values():spine.set_visible(False)
    for i,q in enumerate(scores.itertuples()):
        if q.saga_reference:
            ax.get_yticklabels()[i].set_fontweight('semibold')
            ax.add_patch(Rectangle((-.68,i-.30),.12,.60,fc='#BAA04B',ec='none',clip_on=False))
        for j in range(8):
            value=r.iloc[i,j]
            text='—' if pd.isna(value) else '0.00' if abs(value)<.005 else f'{value:.2f}'
            ax.text(j,i,text,ha='center',va='center',fontsize=7.8,
                    color='white' if pd.notna(value) and abs(value)>.37 else '#333333')

    bx=fig.add_axes(rect(7.4,1.10,1.30,6.25))
    bx.set(xlim=(-12,.8),ylim=(24.5,-.5),yticks=[],xticks=[-10,-5,0])
    bx.axvline(0,color='#A8ADB4',lw=.7)
    for i,q in enumerate(scores.itertuples()):
        bx.hlines(i,q.normZ,0,color='#7A3E85',alpha=.40,lw=.7)
        bx.scatter(q.normZ,i,s=18,c='#7A3E85',lw=0,zorder=3)
    for side in ['top','right','left']:bx.spines[side].set_visible(False)
    bx.spines['bottom'].set_linewidth(.6)
    bx.tick_params(axis='x',labelsize=8,length=3,width=.6)
    bx.set_xlabel('drugZ normZ',fontsize=8.5,labelpad=5)
    fig.text(8.05/width,7.57/height,'HT1080 screen',fontsize=8.8,ha='center')
    cb=fig.colorbar(im,cax=fig.add_axes(rect(1.45,.54,3.65,.12)),orientation='horizontal',ticks=[-.6,-.3,0,.3,.6])
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=8,length=2,pad=2)
    cb.set_label('DepMap Pearson r · unadjusted',fontsize=8.5,labelpad=3)
    fig.legend(handles=[Patch(fc='#BAA04B',label='SAGA reference member'),Patch(fc='#C8CCD1',label='Self-comparison')],
               loc='lower left',bbox_to_anchor=(6.75/width,.12/height),frameon=False,fontsize=7.8,handlelength=1,labelspacing=.5)
    fig.canvas.draw()
    assert np.allclose(ax.images[0].get_array(),r.to_numpy(),equal_nan=True)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=170)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all(g in page.get_text() for g in SEEDS+scores.index.tolist())
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps=[]
        for i,a in enumerate(spans):
            for b in spans[i+1:]:
                ra,rb=fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
                area=(ra&rb).get_area()
                if area>1 and area/min(ra.get_area(),rb.get_area())>.12:overlaps.append([a['text'],b['text']])
        page.get_pixmap(matrix=fitz.Matrix(1.5,1.5),colorspace=fitz.csGRAY).save(OUT/'review'/f'{STEM}_greyscale.png')
    sources=[SOURCE/f for f in ['panel_b_correlations.csv','panel_b_model_counts.csv','panel_b_screen_scores.csv']]
    sources += [XY_SOURCE/f for f in ['genomewide_correlations.npz','model_cohort.csv','seed_gene_effect.csv']]
    record=dict(stem=STEM,dimensions_inches=[width,height],font='Open Sans',seeds=SEEDS,genes=scores.index.tolist(),
                selection='Exactly original HT1080 drugZ depletion ranks 1–25; six X-linked and two Y-linked profiles; X/Y chromatin-gene pairs adjacent; no clustering.',
                saga_rows=scores.index[scores.saga_reference].tolist(),omitted_saga_hit='TAF5L, depletion rank 37, outside the top-25 cutoff',
                r_scale=[-.6,.6],numeric_precision=2,finite_correlations=199,masked_cells=1,column_coverage=coverage.reset_index().to_dict('records'),
                cohort_limit='X columns use 1,186 models; Y columns use the same 355-model subset, including 192 Male, 149 Female and 14 Unknown. No Y-retaining restriction; comparisons are not matched-cohort estimates.',
                analysis_refit=False,original_values_preserved=True,colour_clipping=False,text_outside_page=outside,text_overlaps=overlaps,
                source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (DATA/'render_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k not in ['genes','source_sha256']},indent=2))
    assert not outside,outside
    assert not overlaps,overlaps

if __name__=='__main__':main()
