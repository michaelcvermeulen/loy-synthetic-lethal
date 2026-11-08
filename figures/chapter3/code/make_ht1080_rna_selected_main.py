"""Main RNA results figures selected from the applied programme analysis."""
from pathlib import Path
import sys,hashlib,json,shutil,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from adjustText import adjust_text
FIG=Path(__file__).resolve().parents[1];THESIS=FIG.parents[1]
ROOT=THESIS.parent/'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE'
T=ROOT/'outputs/programme_analysis_v1/tables'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8);np.random.seed(42)
blue='#3C5488';red='#E64B35';grey='#C8C8C8'
cmap=LinearSegmentedColormap.from_list('signed',[blue,'#F7F7F7',red])
d=pd.read_csv(T/'gene_consistency.csv');g=pd.read_csv(T/'programme_consistency.csv');tf=pd.read_csv(T/'TF_consistency.csv')
def export(fig,stem,tables,selection):
 src=FIG/'source_data'/stem.removeprefix('ch03_main_');src.mkdir(parents=True,exist_ok=True)
 for name,table in tables.items():table.to_csv(src/f'{name}.csv',index=False)
 records=[]
 for ext in ['pdf','svg','png']:
  p=FIG/f'{stem}.{ext}';stage=FIG/f'{stem}.staging.{ext}'
  fig.savefig(stage,dpi=600);stage.replace(p)
  records.append(dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 fig.savefig(FIG/'previews'/f'{stem}_preview.png',dpi=160);plt.close(fig)
 (src/'provenance.json').write_text(json.dumps({'font':'Open Sans','source_folder':str(T),'selection':selection,'exports':records},indent=2))
 print(stem,flush=True)

fig,(a,b)=plt.subplots(1,2,figsize=(183/25.4,125/25.4))
fig.subplots_adjust(left=.10,right=.97,bottom=.15,top=.83,wspace=.40)
for ax,letter,title in [(a,'a','Overall expression differences'),(b,'b','Agreement across sampling periods')]:
 pos=ax.get_position();fig.text(pos.x0-.06,.905,letter,fontsize=12,weight='semibold');fig.text(pos.x0,.905,title,fontsize=7.5,weight='semibold')
s=d.directionally_consistent
cols=np.where(s,np.where(d.logFC>0,red,blue),grey)
a.scatter(d.logFC,-np.log10(d.FDR),c=cols,s=np.where(s,7,3),alpha=.55,lw=0,rasterized=True)
for x in [-1,1]:a.axvline(x,c='#888888',ls='--',lw=.6)
a.axhline(-np.log10(.05),c='#888888',ls='--',lw=.6)
a.set(xlabel='log₂FC (Y-depleted / sgCTRL)',ylabel='−log₁₀(BH-adjusted P)',xlim=(-7.5,6.5),ylim=(-.05,6.8))
selected=pd.concat([d[s&(np.sign(d.logFC)==sign)].sort_values('FDR').head(12) for sign in [-1,1]])
texts=[a.text(r.logFC,-np.log10(r.FDR)+.08,r.gene_name,fontsize=5.5,color=blue if r.logFC<0 else red) for r in selected.itertuples()]
adjust_text(texts,ax=a,iter_lim=170,expand=(1.2,1.3),force_text=(.4,.6),arrowprops=dict(arrowstyle='-',color='#888888',lw=.4))
b.scatter(d.AUG_logFC,d.MAY_logFC,c=grey,s=3,alpha=.35,lw=0,rasterized=True)
b.scatter(d.loc[s,'AUG_logFC'],d.loc[s,'MAY_logFC'],c=cols[s],s=8,alpha=.75,lw=0,rasterized=True)
lim=max(d.AUG_logFC.abs().max(),d.MAY_logFC.abs().max())+.3
b.plot([-lim,lim],[-lim,lim],ls='--',c='#888888',lw=.6)
b.axvline(0,c='#BBBBBB',lw=.5);b.axhline(0,c='#BBBBBB',lw=.5)
b.set(xlabel='August log₂FC',ylabel='May log₂FC',xlim=(-lim,lim),ylim=(-lim,lim))
b.text(.03,.97,f'All tested genes: r = {d.AUG_logFC.corr(d.MAY_logFC):.2f}',transform=b.transAxes,va='top',fontsize=7)
b.set_box_aspect(1)
b.set_aspect('equal',adjustable='box')
b.set_anchor('N')
# Label large effects shared across periods, using a deterministic display rule.
agree=d[s & ~d.gene_name.str.startswith('ENSG')].copy()
agree['period_effect_gap']=(agree.AUG_logFC-agree.MAY_logFC).abs()
agree['minimum_period_abs_logFC']=agree[['AUG_logFC','MAY_logFC']].abs().min(axis=1)
agree=agree[agree.period_effect_gap<=.5]
agree=pd.concat([agree[np.sign(agree.logFC)==sign].sort_values(
 ['minimum_period_abs_logFC','FDR','gene_id'],ascending=[False,True,True]).head(6) for sign in [-1,1]])
fig.canvas.draw()
for sign in [-1,1]:
 labelled=agree[np.sign(agree.logFC)==sign].sort_values('MAY_logFC')
 label_y=np.linspace(-6.5,-1.0,len(labelled)) if sign<0 else np.linspace(2.3,6.7,len(labelled))
 for r,y in zip(labelled.itertuples(),label_y):
  b.annotate(r.gene_name,xy=(r.AUG_logFC,r.MAY_logFC),
   xytext=(-5.7 if sign<0 else 5.3,y),ha='right' if sign<0 else 'left',va='center',
   fontsize=5.8,color=blue if sign<0 else red,
   arrowprops=dict(arrowstyle='-',color='#888888',lw=.45),zorder=6)
fig.canvas.draw()
box=b.get_window_extent()
assert abs(box.width-box.height)<.1, 'Panel b must have a square plotting area'
fig.text(.10,.97,'FDR < 0.05 and |log₂FC| ≥ 1: 433 lower / 229 higher',fontsize=8,va='top')
export(fig,'ch03_main_ht1080_rna_expression_changes',{'all_genes':d,'labelled_genes':selected,'panel_b_labelled_genes':agree},'Core autosomal/non-PAR X model. Panel a: 12 lowest-FDR qualifying genes per direction labelled. Panel b: square, equal unit scaling; six named qualifying genes per direction with largest minimum absolute period effect among those with August/May log2FC difference <=0.5 (ties: overall FDR then gene ID). This is a display rule, not a new inferential test. All genes shown in both panels; Pearson r uses all 25526 genes.')
if '--expression-only' in sys.argv:
 sys.exit(0)

sig=g[g.directionally_consistent]
selected=pd.concat([q.sort_values('FDR').head(4) for _,q in sig.groupby('analysis_collection')]+[sig[sig.NES>0].sort_values('FDR').head(3)]).drop_duplicates('pathway')
selected=selected.sort_values(['NES','FDR'])
ft=tf[tf.directionally_consistent].sort_values('FDR').head(12).sort_values('score')
fig=plt.figure(figsize=(183/25.4,205/25.4))
a=fig.add_axes([.55,.49,.28,.43]);b=fig.add_axes([.55,.09,.28,.29])
for ax,letter,title in [(a,'a','Pathway enrichment across periods'),(b,'b','Inferred transcription-factor activity')]:
 pos=ax.get_position();fig.text(.03,pos.y1+.027,letter,fontsize=12,weight='semibold');fig.text(.11,pos.y1+.027,title,fontsize=9,weight='semibold')
def nice(v):
 for pre in ['HALLMARK_','REACTOME_','GOBP_']:v=v.removeprefix(pre)
 return '\n'.join(textwrap.wrap(v.replace('_',' ').lower(),44))
def heat(ax,frame,values,pvals,names,bar_title):
 mat=frame[values].to_numpy();pv=frame[pvals].to_numpy();mx=max(abs(mat).max(),1)
 im=ax.imshow(mat,cmap=cmap,norm=TwoSlopeNorm(vmin=-mx,vcenter=0,vmax=mx),aspect='auto')
 ax.set_yticks(range(len(frame)),names,fontsize=6.6);ax.set_xticks(range(3),['Overall','August','May'],fontsize=7)
 ax.tick_params(axis='both',length=0)
 for i,j in zip(*np.where(pv<.05)):ax.text(j,i,'*',ha='center',va='center',fontsize=9,color='white' if abs(mat[i,j])>.6*mx else '#222222')
 for i,val in enumerate(pv[:,0]):ax.text(3.03,i,f'{val:.1e}',fontsize=6.2,va='center',clip_on=False)
 ax.text(3.03,-1,'Overall\nFDR',fontsize=6.2,clip_on=False)
 pos=ax.get_position();cb=fig.colorbar(im,cax=fig.add_axes([.95,pos.y0,.012,pos.height]));cb.ax.set_title(bar_title,fontsize=6.5,pad=6);cb.ax.tick_params(labelsize=6)
heat(a,selected,['NES','AUG_2025_NES','MAY_2026_NES'],['FDR','AUG_2025_FDR','MAY_2026_FDR'],[nice(s) for s in selected.pathway],'NES')
heat(b,ft,['score','AUG_2025score','MAY_2026score'],['FDR','AUG_2025FDR','MAY_2026FDR'],ft.source,'Score')
fig.text(.11,.025,'* FDR < 0.05; selected patterns shown, complete results in the supplement',fontsize=7,color='#555555')
export(fig,'ch03_main_ht1080_rna_programmes_regulators',{'pathways':selected,'TFs':ft},'Four lowest-overall-FDR consistent pathways per collection, plus three lowest-FDR higher-enrichment pathways; deduplicated. Twelve lowest-overall-FDR consistent TF patterns. Full three-collection and TF testing families retained.')
