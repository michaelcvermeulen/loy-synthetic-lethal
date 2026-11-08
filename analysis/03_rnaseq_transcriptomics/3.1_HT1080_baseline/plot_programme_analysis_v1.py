"""Thesis-style review figure for the implemented RNA programme plan."""
from pathlib import Path
import sys,json,hashlib,shutil,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from adjustText import adjust_text
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs/programme_analysis_v1';T=OUT/'tables'
THESIS=ROOT.parents[3]/'SEPT_11_2026_THESIS';FIG=THESIS/'figures/chapter3'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8);np.random.seed(42)
STEM='ch03_supp_ht1080_rna_programme_overview'
blue='#3C5488';red='#E64B35';grey='#C8C8C8'
cmap=LinearSegmentedColormap.from_list('signed',[blue,'#F7F7F7',red])
d=pd.read_csv(T/'gene_consistency.csv');parent=pd.read_csv(T/'DE_CTRL_vs_parental.csv')
g=pd.read_csv(T/'programme_consistency.csv');tf=pd.read_csv(T/'TF_consistency.csv')
fig=plt.figure(figsize=(183/25.4,235/25.4))
a=fig.add_axes([.105,.74,.355,.19]);b=fig.add_axes([.615,.74,.355,.19])
c=fig.add_axes([.57,.395,.355,.235]);f=fig.add_axes([.57,.055,.355,.235])
for ax,letter,title in [(a,'a','Period-specific expression effects'),(b,'b','sgCTRL versus parental (August)'),(c,'c','Leading pathway patterns'),(f,'d','Inferred transcription-factor activity')]:
 p=ax.get_position();fig.text(.025 if letter in 'acd' else .53,p.y1+.018,letter,fontsize=12,weight='semibold')
 fig.text(.105 if letter in 'acd' else .615,p.y1+.018,title,fontsize=8,weight='semibold')
s=d.directionally_consistent
a.scatter(d.AUG_logFC,d.MAY_logFC,c=grey,s=2,alpha=.35,lw=0,rasterized=True)
a.scatter(d.loc[s,'AUG_logFC'],d.loc[s,'MAY_logFC'],c=np.where(d.loc[s,'logFC']>0,red,blue),s=5,alpha=.7,lw=0,rasterized=True)
lim=max(d.AUG_logFC.abs().max(),d.MAY_logFC.abs().max())+.2
a.plot([-lim,lim],[-lim,lim],c='#777777',ls='--',lw=.6)
a.axhline(0,c='#BBBBBB',lw=.5);a.axvline(0,c='#BBBBBB',lw=.5)
a.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel='August log₂FC',ylabel='May log₂FC')
a.text(.04,.96,'662 directionally consistent\n229 higher / 433 lower',transform=a.transAxes,va='top',fontsize=6.5)
ps=(parent.FDR<.05)&(parent.logFC.abs()>=1)
b.scatter(parent.logFC,-np.log10(parent.FDR),s=3,c=grey,alpha=.35,lw=0,rasterized=True)
b.scatter(parent.loc[ps,'logFC'],-np.log10(parent.loc[ps,'FDR']),s=14,c=np.where(parent.loc[ps,'logFC']>0,red,blue),lw=0)
for q in [-1,1]:b.axvline(q,c='#888888',ls='--',lw=.6)
b.axhline(-np.log10(.05),c='#888888',ls='--',lw=.6)
b.set(xlim=(-5.2,6),ylim=(-.05,3.0),xlabel='log₂FC (sgCTRL / parental)',ylabel='−log₁₀(BH-adjusted P)')
lab=parent[ps].sort_values('logFC')
tx=[b.text(r.logFC+(.2 if i%2 else -.2),-np.log10(r.FDR)+.1+(i%4)*.13,r.gene_name,fontsize=5.8) for i,r in enumerate(lab.itertuples())]
adjust_text(tx,ax=b,target_x=lab.logFC.to_numpy(),target_y=-np.log10(lab.FDR.to_numpy()),expand=(1.35,1.5),force_text=(.6,.8),iter_lim=180,arrowprops=dict(arrowstyle='-',color='#888888',lw=.4))
def nice(s):
 for prefix in ['HALLMARK_','REACTOME_','GOBP_']:s=s.removeprefix(prefix)
 return '\n'.join(textwrap.wrap(s.replace('_',' ').lower(),48))
# Representation from each collection avoids allowing GO size alone to determine display.
selected=pd.concat([q[q.directionally_consistent].sort_values('FDR').head(4) for _,q in g.groupby('analysis_collection')])
selected=selected.sort_values(['NES','FDR'])
mat=selected[['NES','AUG_2025_NES','MAY_2026_NES']].to_numpy()
im=c.imshow(mat,cmap=cmap,norm=TwoSlopeNorm(vmin=-max(abs(mat).max(),1),vcenter=0,vmax=max(abs(mat).max(),1)),aspect='auto')
c.set_yticks(range(len(selected)),[nice(s) for s in selected.pathway],fontsize=6.8)
c.set_xticks(range(3),['Overall','August','May'],fontsize=7);c.tick_params(axis='both',length=0)
pv=selected[['FDR','AUG_2025_FDR','MAY_2026_FDR']].to_numpy()
for i,j in zip(*np.where(pv<.05)):c.text(j,i,'*',ha='center',va='center',fontsize=10,color='white' if abs(mat[i,j])>1.4 else '#222222')
cb=fig.colorbar(im,cax=fig.add_axes([.946,.395,.013,.235]));cb.ax.set_title('NES',fontsize=7,pad=6);cb.ax.tick_params(labelsize=6)
ft=tf[tf.directionally_consistent].sort_values('FDR').head(12).sort_values('score')
fm=ft[['score','AUG_2025score','MAY_2026score']].to_numpy();vmax=max(abs(fm).max(),1)
im=f.imshow(fm,cmap=cmap,norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax),aspect='auto')
f.set_yticks(range(len(ft)),ft.source,fontsize=7);f.set_xticks(range(3),['Overall','August','May'],fontsize=7);f.tick_params(axis='both',length=0)
fp=ft[['FDR','AUG_2025FDR','MAY_2026FDR']].to_numpy()
for i,j in zip(*np.where(fp<.05)):f.text(j,i,'*',ha='center',va='center',fontsize=10,color='white' if abs(fm[i,j])>vmax*.6 else '#222222')
cb=fig.colorbar(im,cax=fig.add_axes([.946,.055,.013,.235]));cb.ax.set_title('Score',fontsize=7,pad=6);cb.ax.tick_params(labelsize=6)
fig.text(.105,.18,'60 consistent TF patterns\n\n* FDR < 0.05\nFull statistics in source tables',fontsize=7,color='#555555',va='top')
fig.text(.105,.985,'HT1080 RNA programme analysis',fontsize=10,weight='semibold',va='top')
src=FIG/'source_data/ht1080_rna_programme_overview';src.mkdir(parents=True,exist_ok=True)
selected.to_csv(src/'panel_c_pathways.csv',index=False);ft.to_csv(src/'panel_d_TFs.csv',index=False)
for name in ['gene_consistency','DE_CTRL_vs_parental']:shutil.copy2(T/f'{name}.csv',src/f'{name}.csv')
records=[]
for ext in ['pdf','svg','png']:
 path=FIG/f'{STEM}.{ext}';fig.savefig(path,dpi=600);shutil.copy2(path,OUT/'figures'/path.name)
 records.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
fig.savefig(FIG/'previews'/f'{STEM}_preview.png',dpi=150);plt.close(fig)
(src/'provenance.json').write_text(json.dumps({'exports':records,'font':'Open Sans','dimensions_mm':[183,235],
 'pathway_display':'Four smallest overall FDR directionally consistent pathways per collection; not all significant pathways.',
 'TF_display':'Twelve smallest overall FDR directionally consistent TF target patterns.'},indent=2))
print(STEM)
