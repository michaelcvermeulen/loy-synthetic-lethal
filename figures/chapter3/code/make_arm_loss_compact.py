"""Condensed chromosome-loss comparison; display saved results without refitting."""
from pathlib import Path
import sys, json, hashlib, argparse, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
SRC=FIG/'source_data/arm_loss_comparison'
PREVIOUS=FIG/'source_data/arm_loss_consolidated'
OUT=FIG/'source_data/arm_loss_compact_20261030'
OUT.mkdir(exist_ok=True)
SNAP=OUT/'inputs'
SNAP.mkdir(exist_ok=True)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path,default=FIG/'selected/page_reduction_20261030')
args=parser.parse_args()
DEST=args.output_dir.resolve()
DEST.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
STEM='ch03_supp_arm_loss_comparison'
INPUTS=[]

def read(path,sep):
    INPUTS.append({'file':path.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    frozen=SNAP/path.name
    if frozen.exists():
        assert hashlib.sha256(frozen.read_bytes()).hexdigest()==INPUTS[-1]['sha256']
    else:
        shutil.copy2(path,frozen)
    return pd.read_csv(frozen,sep=sep)

hallmark=read(SRC/'main_a_hallmark.tsv','\t')
go=read(SRC/'main_b_go.tsv','\t')
programmes=read(PREVIOUS/'panel_d_programmes.csv',',')
times=read(PREVIOUS/'panel_f_four_times.csv',',')
assert len(hallmark)==49 and len(go)==3178
assert hallmark['class'].value_counts().get('shared')==16
assert hallmark['class'].value_counts().get('opposite',0)==0
ORDER=['RPTEC_early_sc_Yq','RPTEC_early_sc_3p','A2058_1q','A2780_1q_vsWT','A2780_1q_vsTK','AGS_1q_vsTK','MCF10A_1q','MCF10Ahr_1q','RPTEC_3p_early','RPTEC_3p_late','A2058_7p','A2058_8q','HCT116_8q','RPE1_mono1018','RPE1_mono13','RPE1_mono19p']
LABELS=['RPTEC: early Yq','RPTEC: early 3p','1q | A2058 [1]','1q | A2780 vs WT [2]','1q | A2780 vs TK [2]','1q | AGS vs TK','1q | MCF10A','1q | MCF10A-HRAS','3p | RPTEC day 18 [3]','3p | RPTEC day 90 [3]','7p | A2058 [1]','8q | A2058 [1]','8q | HCT116','10+18 | RPE1 p53 KO [4]','13 | RPE1 p53 KO [4]','19p | RPE1 p53 KD']
SETS={'HALLMARK_P53_PATHWAY':'p53','HALLMARK_E2F_TARGETS':'E2F','HALLMARK_G2M_CHECKPOINT':'G2/M','HALLMARK_MYC_TARGETS_V1':'MYC V1','HALLMARK_OXIDATIVE_PHOSPHORYLATION':'OXPHOS','RIBOSOMAL_PROTEIN':'Ribosomal','HALLMARK_INTERFERON_ALPHA_RESPONSE':'IFN-α','HALLMARK_INTERFERON_GAMMA_RESPONSE':'IFN-γ','HALLMARK_TNFA_SIGNALING_VIA_NFKB':'TNFα–NF-κB','HALLMARK_INFLAMMATORY_RESPONSE':'Inflammatory','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION':'EMT'}
TIME_SETS=['HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_MYC_TARGETS_V1','HALLMARK_P53_PATHWAY']
COLOURS={'neither':'#C8C8C8','3p_only':'#00A087','Yq_only':'#F39B7F','shared':'#3C5488','opposite':'#E64B35'}
MARKERS={'neither':'.','3p_only':'^','Yq_only':'s','shared':'o','opposite':'D'}
CMAP=LinearSegmentedColormap.from_list('signed',['#3C5488','#FFFFFF','#E64B35'])
CMAP.set_bad('#DDDDDD')

def scatter(ax,table,label,rho):
    ax.axhline(0,color='.84',lw=.5);ax.axvline(0,color='.84',lw=.5)
    ax.plot([-3.5,3.5],[-3.5,3.5],'--',color='.75',lw=.5)
    for cls in COLOURS:
        sub=table[table['class']==cls]
        ax.scatter(sub.NES_3p,sub.NES_Yq,c=COLOURS[cls],marker=MARKERS[cls],s=15 if len(table)<100 else (17 if cls=='opposite' else 6),linewidths=0,alpha=.65 if cls=='neither' else .95)
    ax.set(xlim=(-3.5,3.5),ylim=(-3.5,3.5),xticks=[-3,0,3],yticks=[-3,0,3],xlabel='NES, 3p loss vs remaining cells',ylabel='NES, Yq loss vs Yq-intact cells')
    ax.text(.03,.97,label,transform=ax.transAxes,ha='left',va='top',fontsize=8)
    ax.text(.98,.025,f'ρ = {rho:.2f}; {len(table):,} sets',transform=ax.transAxes,ha='right',va='bottom',fontsize=7.5)

def annotate(ax,table,key,text,xytext):
    row=table.set_index('pathway').loc[key]
    ax.annotate(text,(row.NES_3p,row.NES_Yq),xytext=xytext,fontsize=7.5,va='center',arrowprops={'arrowstyle':'-','lw':.5,'color':'.4'},bbox={'facecolor':'white','alpha':.95,'edgecolor':'none','pad':1})

fig=plt.figure(figsize=(8,6.7))
a=fig.add_axes([.085,.640,.390,.325])
b=fig.add_axes([.570,.640,.390,.325])
scatter(a,hallmark,'Hallmark',.71390);scatter(b,go,'GO:BP',.49870)
annotate(a,hallmark,'HALLMARK_P53_PATHWAY','p53',(-.1,2.9))
annotate(a,hallmark,'HALLMARK_E2F_TARGETS','E2F',(-3.25,-1.7))
annotate(b,go,'GOBP_RESPONSE_TO_TYPE_I_INTERFERON','Type I interferon',(-3.2,.3))
handles=[Line2D([],[],marker=MARKERS[k],color=COLOURS[k],lw=0,markersize=4,label=v) for k,v in [('shared','Both, same direction'),('opposite','Both, opposite'),('3p_only','3p only'),('Yq_only','Yq only'),('neither','Neither')]]
fig.text(.085,.560,'FDR <0.05:',fontsize=7.5,va='center')
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.58,.560),ncol=5,fontsize=7,frameon=False,columnspacing=.85,handletextpad=.3)

c=fig.add_axes([.255,.240,.695,.280])
nes=programmes.pivot(index='contrast',columns='pathway',values='NES').reindex(index=ORDER,columns=list(SETS))
q=programmes.pivot(index='contrast',columns='pathway',values='padj').reindex(index=ORDER,columns=list(SETS))
im=c.imshow(nes,cmap=CMAP,vmin=-4,vmax=4,aspect='auto',interpolation='none')
yy,xx=np.where(q.to_numpy()<.05);c.scatter(xx,yy,s=6,c='#222222')
c.set_xticks(range(len(SETS)),list(SETS.values()),rotation=48,ha='right',fontsize=7.5)
c.set_yticks(range(len(ORDER)),LABELS,fontsize=7.5)
c.axhline(1.5,color='white',lw=2);c.tick_params(length=0)
for spine in c.spines.values():spine.set_visible(False)
cax=fig.add_axes([.964,.300,.012,.160])
cb=fig.colorbar(im,cax=cax,ticks=[-4,0,4]);cb.ax.set_title('NES',fontsize=7.5,pad=6);cb.ax.tick_params(labelsize=7,pad=2)

d=fig.add_axes([.255,.040,.695,.090])
tnes=times.pivot(index='day',columns='pathway',values='NES').reindex(index=[5,10,18,90],columns=TIME_SETS)
tq=times.pivot(index='day',columns='pathway',values='padj').reindex(index=[5,10,18,90],columns=TIME_SETS)
assert tnes.shape==(4,6) and not tnes.isna().any().any()
d.imshow(tnes,cmap=CMAP,vmin=-4,vmax=4,aspect='auto',interpolation='none')
for yy in range(4):
    for xx in range(6):
        val=tnes.iloc[yy,xx]
        d.text(xx,yy,f'{val:+.2f}'+('*' if tq.iloc[yy,xx]<.05 else ''),ha='center',va='center',fontsize=7.5,color='white' if abs(val)>2.6 else '#222222')
d.axhline(1.5,color='white',lw=2)
d.set_xticks(range(6),[SETS[x] for x in TIME_SETS],fontsize=7.5)
d.set_yticks(range(4),['Day 5 (scRNA)','Day 10 (scRNA)','Day 18 (bulk)','Day 90 (bulk)'],fontsize=7.5)
d.tick_params(length=0)
for spine in d.spines.values():spine.set_visible(False)
d.text(1,1.10,'NES; * FDR <0.05',transform=d.transAxes,ha='right',fontsize=7)
for x,y,label in [(.025,.960,'a'),(.510,.960,'b'),(.025,.528,'c'),(.025,.140,'d')]:
    fig.text(x,y,label,fontsize=12,weight='semibold',va='bottom')

# Layout-only checks: preserve every point, matrix entry and annotation.
for ax,table in [(a,hallmark),(b,go)]:
    assert sum(len(artist.get_offsets()) for artist in ax.collections)==len(table)
    assert table.NES_3p.between(*ax.get_xlim()).all()
    assert table.NES_Yq.between(*ax.get_ylim()).all()
assert nes.shape==(16,11) and tnes.shape==(4,6)
assert np.allclose(c.images[0].get_array().filled(np.nan),nes.to_numpy(),equal_nan=True)
assert np.allclose(d.images[0].get_array().filled(np.nan),tnes.to_numpy(),equal_nan=True)
assert len(c.collections[0].get_offsets())==int((q.to_numpy()<.05).sum())
assert len(d.texts)==25  # 24 annotated cells plus the NES/FDR key.
fig.canvas.draw()
from matplotlib.text import Text
renderer=fig.canvas.get_renderer()
clipped=[]
for item in fig.findobj(Text):
    if not item.get_visible() or not item.get_text().strip():continue
    bounds=item.get_window_extent(renderer)
    if bounds.x0 < -1 or bounds.y0 < -1 or bounds.x1 > fig.bbox.width+1 or bounds.y1 > fig.bbox.height+1:
        clipped.append(item.get_text())
assert not clipped,clipped

for ext in ['pdf','svg','png']:
    stage=DEST/f'{STEM}.staging.{ext}'
    fig.savefig(stage,dpi=600,facecolor='white');stage.replace(DEST/f'{STEM}.{ext}')
fig.savefig(DEST/f'{STEM}_preview.png',dpi=200,facecolor='white')
preview=DEST/f'{STEM}_preview.png'
shutil.copy2(preview,DEST/f'{STEM}_preview_{hashlib.sha256(preview.read_bytes()).hexdigest()[:12]}.png')
plt.close(fig)
hallmark.to_csv(OUT/'panel_a_hallmark.tsv',sep='\t',index=False)
go.to_csv(OUT/'panel_b_go.tsv',sep='\t',index=False)
programmes[programmes.pathway.isin(SETS)].to_csv(OUT/'panel_c_programmes.csv',index=False)
times.to_csv(OUT/'panel_d_four_times.csv',index=False)
(OUT/'provenance.json').write_text(json.dumps({'date':'2026-10-30','script':str(Path(__file__)),'inputs':INPUTS,'statistical_refit':False,'font':'Open Sans','figure_inches':[8,6.7],'previous_figure_inches':[8,8.7],'height_reduction_percent':100*(1-6.7/8.7),'text_sizes_unchanged':True,'panel_changes':'Layout only: reduced top/bottom margins and row gutters; slightly shallower, wider scatter panels; wider heatmaps. All four panels, displayed values, labels, correlations and FDR annotations preserved.','FDR':'Original collection/contrast families, unchanged; display selection does not recalculate FDR.'},indent=2),encoding='utf-8')
print('Rendered four-panel chromosome-loss comparison from saved data.')
