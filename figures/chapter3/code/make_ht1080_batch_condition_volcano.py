"""PAR-inclusive companion volcano; uses fitted results, never refits in Python."""
from pathlib import Path
import sys, json, hashlib, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

FIG=Path(__file__).resolve().parents[1]
THESIS=FIG.parents[1]
ROOT=THESIS.parent/'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8.5)
np.random.seed(42)
STEM='ch03_supp_ht1080_batch_condition_volcano'
SOURCE=ROOT/'outputs/tables/batchDE_batch_condition_with_PAR__combined.csv'
d=pd.read_csv(SOURCE)
assert len(d)==25553 and d.PAR.sum()==27
d['passes_cutoffs']=(d.FDR<.05)&(d.logFC.abs()>=1)
d['y']=-np.log10(d.FDR)
assert np.isfinite(d[['logFC','y']]).all().all()
purple='#8064A2';blue='#3C5488';red='#E64B35';grey='#C6C6C6'
fig,ax=plt.subplots(figsize=(183/25.4,185/25.4))
fig.subplots_adjust(left=.115,right=.975,bottom=.10,top=.84)
ax.scatter(d.logFC,d.y,s=5,c=grey,alpha=.42,lw=0,rasterized=True)
for sign,col in [(-1,blue),(1,red)]:
 q=d[d.passes_cutoffs & ~d.PAR & (np.sign(d.logFC)==sign)]
 ax.scatter(q.logFC,q.y,s=12,c=col,alpha=.8,lw=0,rasterized=True)
q=d[d.PAR]
ax.scatter(q.logFC,q.y,s=25,marker='D',c=purple,edgecolors='white',linewidths=.35,zorder=4)
for x in [-1,1]: ax.axvline(x,color='#888888',ls=(0,(4,3)),lw=.65,zorder=0)
ax.axhline(-np.log10(.05),color='#888888',ls=(0,(4,3)),lw=.65,zorder=0)
ax.set(xlim=(-7.5,6.0),ylim=(-.12,6.45),xlabel='log₂ fold change (Y-depleted / sgCTRL)',ylabel='−log₁₀(BH-adjusted P value)')
ax.set_xticks(np.arange(-6,7,2))
# Prespecified label rule: strongest adjusted P values and largest effects per direction,
# plus all expressed protein-coding PAR genes and all threshold-passing PAR features.
ids=[]
for sign in [-1,1]:
 q=d[d.passes_cutoffs & ~d.PAR & (np.sign(d.logFC)==sign)]
 ids+=q.sort_values(['FDR','gene_id']).head(18).gene_id.tolist()
 ids+=q.assign(magnitude=q.logFC.abs()).sort_values(['magnitude','gene_id'],ascending=[False,True]).head(4).gene_id.tolist()
ids+=d[d.PAR & ((d.gene_type=='protein_coding')|d.passes_cutoffs)].gene_id.tolist()
labels=d[d.gene_id.isin(ids)].copy()
texts=[]
for r in labels.itertuples():
 col=purple if r.PAR else (blue if r.logFC<0 else red)
 # Seed PAR labels away from the crowded central axis before collision adjustment.
 shift=-.7 if r.PAR and r.logFC<0 else (.7 if r.PAR else 0)
 texts.append(ax.text(r.logFC+shift,r.y+.07,r.gene_name,fontsize=6.5,color=col,
                      weight='semibold' if r.PAR else 'normal',zorder=6))
adjust_text(texts,x=d.loc[d.passes_cutoffs|d.PAR,'logFC'].to_numpy(),
            y=d.loc[d.passes_cutoffs|d.PAR,'y'].to_numpy(),
            target_x=labels.logFC.to_numpy(),target_y=labels.y.to_numpy(),
            ax=ax,expand=(1.12,1.2),force_text=(.35,.6),force_static=(.1,.2),
            force_pull=(.005,.01),iter_lim=180,
            arrowprops=dict(arrowstyle='-',color='#888888',lw=.4),min_arrow_len=2)
fig.text(.115,.955,'HT1080 | batch + condition',weight='semibold',fontsize=11)
fig.text(.115,.925,'FDR < 0.05 and |log₂FC| ≥ 1',fontsize=8.5)
fig.text(.975,.925,'437 lower  ·  232 higher',ha='right',fontsize=8.5)
handles=[Line2D([],[],ls='',marker='o',color=c,label=l,markersize=4) for c,l in
         [(blue,'Lower'),(red,'Higher'),(grey,'Below either cutoff')]]
handles.append(Line2D([],[],ls='',marker='D',color=purple,label='PAR genes',markersize=4))
fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.10,.905),ncol=4,fontsize=7.3,
           columnspacing=1.6,handletextpad=.5)
src=FIG/'source_data/ht1080_batch_condition_volcano';src.mkdir(parents=True,exist_ok=True)
shutil.copy2(SOURCE,src/SOURCE.name)
labels.to_csv(src/'labelled_genes.csv',index=False)
d[d.PAR].to_csv(src/'PAR_genes.csv',index=False)
for ext in ['pdf','svg','png']:
 target=FIG/f'{STEM}.{ext}'
 fig.savefig(target,dpi=600)
fig.savefig(FIG/'previews'/f'{STEM}_preview.png',dpi=160)
plt.close(fig)
meta={'source':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      'model':'limma-voom ~ batch + treatment, robust empirical Bayes',
      'tested_features':len(d),'PAR_features':int(d.PAR.sum()),'labels':len(labels),
      'FDR_cutoff':.05,'absolute_log2FC_cutoff':1,'font':'Open Sans',
      'note':'PAR-inclusive companion fit; BH across 25553 features. Original model unchanged.',
      'PAR_annotation':'Inherited GRCh38 PAR interval-overlap annotation. No allele-of-origin interpretation.',
      'exports':{ext:hashlib.sha256((FIG/f'{STEM}.{ext}').read_bytes()).hexdigest() for ext in ['pdf','svg','png']}}
(src/'provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print(json.dumps(meta,indent=2))
