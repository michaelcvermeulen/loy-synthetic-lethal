"""Exploratory combined HT1080 bulk-guide measurements; C5 and B6 excluded."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8.5)
STEM='ch03_supp_ht1080_ddpcr_combined_guides_preview'
OUT=FIG/'source_data/ht1080_ddpcr_combined_guides_preview'
OUT.mkdir(exist_ok=True)
source=FIG/'source_data/ht1080_initial_guide_ddpcr/bulk_population_measurements.csv'
d=pd.read_csv(source)
# One curated measurement table, avoiding duplicate rows from overlapping summaries.
d=d[d.condition.isin(['CP','Q53','CP + Q53','MIX/MIX','CTRL2']) & ~d.sample_id.isin(['C5','B6'])].copy()
assert d.concentration_ratio.notna().all()
assert not d.duplicated(['source_path','export_well','target']).any()
assert len(d)==28
assert not d.sample_id.isin(['C5','B6']).any()
# Recheck ratios from the recorded concentrations and raw instrument tables.
for path,g in d.groupby('source_path'):
 raw=pd.read_csv(path)
 col=next(c for c in raw if c.startswith('Conc(copies/'))
 for _,r in g.iterrows():
  y=raw[raw.Well.eq(r.export_well)&raw.Target.eq(r.target)]
  x=raw[raw.Well.eq(r.export_well)&raw.Target.str.upper().eq(r.reference.upper())]
  assert len(y)==len(x)==1
  assert np.isclose(float(y.iloc[0][col])/float(x.iloc[0][col]),r.concentration_ratio)
d['reference_copy_number']=d.assay.map({'AMELY/AMELX':1,'DDX3Y/RNaseP':2})
d['estimated_target_copies_per_cell']=d.concentration_ratio*d.reference_copy_number
d['plot_group']=d.condition.replace({'MIX/MIX':'CP + Q53 2x','CTRL2':'sgCTRL'})
d['point_class']=np.where(d.nucleofection_count.eq(2),'Two nucleofections',np.where(d.date.lt('2025-01-01'),'Earlier single nucleofection','Later single nucleofection'))
d.to_csv(OUT/'combined_measurements.csv',index=False)
pd.DataFrame([{'sample_id':x,'reason':'Excluded by author instruction; raw data retained'} for x in ['C5','B6']]).to_csv(OUT/'excluded_samples.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(183/25.4,108/25.4))
fig.subplots_adjust(left=.093,right=.987,bottom=.25,top=.85,wspace=.25)
groups=['CP','Q53','CP + Q53','CP + Q53 2x','sgCTRL'];summary=[]
markers={'Earlier single nucleofection':'s','Later single nucleofection':'o','Two nucleofections':'D'}
for ax,assay,letter in zip(axes,['AMELY/AMELX','DDX3Y/RNaseP'],['a','b']):
 for i,group in enumerate(groups):
  g=d[d.assay.eq(assay)&d.plot_group.eq(group)].sort_values(['date','export_well'])
  v=g.estimated_target_copies_per_cell.to_numpy()
  if not len(v):
   ax.text(i,.045,'NA',transform=ax.get_xaxis_transform(),ha='center',fontsize=8,color='#777777')
   summary.append(dict(assay=assay,group=group,n_measurements=0,median=None,q1=None,q3=None,n_double_nucleofection_measurements=0))
   continue
  med=float(np.median(v));q1,q3=np.quantile(v,[.25,.75])
  color='#3C5488' if group=='sgCTRL' else '#E64B35'
  ax.bar(i,med,width=.67,color=color,zorder=2)
  if group=='CP + Q53 2x':ax.bar(i,med,width=.67,facecolor='none',edgecolor='#87301E',linewidth=.7,hatch='///',zorder=3)
  ax.errorbar(i,med,yerr=[[med-q1],[q3-med]],fmt='none',ecolor='#333333',elinewidth=.8,capsize=2,zorder=4)
  offsets=np.linspace(-.18,.18,len(g)) if len(g)>1 else [0]
  for offset,(_,row) in zip(offsets,g.iterrows()):
   ax.scatter(i+offset,row.estimated_target_copies_per_cell,s=16,marker=markers[row.point_class],facecolor='white',edgecolor='#222222',linewidth=.65,zorder=5)
  y=max(v.max(),q3)+.058
  if .93<y<1.05:y=1.05
  ax.text(i,y,f'{med:.3f}',ha='center',va='bottom',fontsize=7.4,fontweight='semibold')
  ax.text(i,y+.06,f'n={len(v)}',ha='center',va='bottom',fontsize=6.8,color='#666666')
  summary.append(dict(assay=assay,group=group,n_measurements=len(v),median=med,q1=q1,q3=q3,n_double_nucleofection_measurements=int(g.nucleofection_count.eq(2).sum())))
 ax.axhline(1,color='#888888',lw=.8,ls=(0,(4,3)),zorder=1)
 ax.axvline(3.5,color='#DDDDDD',lw=.65)
 ax.set(xlim=(-.55,4.55),ylim=(0,1.35),xticks=range(5),xticklabels=['CP\n1\u00d7','Q53\n1\u00d7','CP + Q53\n1\u00d7','CP + Q53\n2\u00d7','sgCTRL'],yticks=np.arange(0,1.21,.2))
 ax.tick_params(axis='x',length=0,pad=8,labelsize=7.5)
 ax.set_title(assay.split('/')[0],pad=12,fontsize=10)
 ax.text(-.14,1.08,letter,transform=ax.transAxes,fontsize=12,fontweight='semibold')
 ax.grid(axis='y',color='#ECECEC',lw=.45);ax.set_axisbelow(True)
axes[0].set_ylabel('Estimated Y-linked target copies per cell',labelpad=9)
handles=[Line2D([],[],marker=markers[k],color='none',markerfacecolor='white',markeredgecolor='#222222',markersize=4,label=k) for k in markers]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.54,.068),ncol=3,frameon=False,fontsize=6.6,handletextpad=.4,columnspacing=1.1)
fig.text(.54,.027,'1\u00d7 = one nucleofection; 2\u00d7 = two successive nucleofections',ha='center',fontsize=7,color='#444444')
for ext in ['pdf','svg','png']:fig.savefig(FIG/(STEM+'.'+ext),dpi=600,facecolor='white')
fig.savefig(FIG/'previews'/(STEM+'.png'),dpi=170,facecolor='white');plt.close(fig)
pd.DataFrame(summary).to_csv(OUT/'combined_group_summary.csv',index=False)
meta={'status':'exploratory preview; does not replace approved main figure or Results','source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'raw_concentration_checks':'passed','n_measurements':len(d),'excluded_samples':['C5','B6'],'exclusion_basis':'author instruction; B6 excluded on 17 September 2026 without a specified technical QC reason','normalization':'AMELY/AMELX x 1; DDX3Y/RNaseP x 2; no control normalization','rounds':'CP + Q53 single and double nucleofections have separate bars. K5 is the sole double-nucleofection AMELY measurement; DDX3Y is unavailable for this group.','summary':'Measurement-weighted median/IQR; not independent biological replicates','point_shapes':markers}
(OUT/'validation.json').write_text(json.dumps(meta,indent=2)+'\n')
print(pd.DataFrame(summary).to_string(index=False))
