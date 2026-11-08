"""Plot H2452 initial-guide ddPCR exactly as recorded in the master table.
Author-selected combined 4 and 6 December display; exact recorded ratios preserved.
"""
from pathlib import Path
import argparse,sys,json,hashlib,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[3];FIG=ROOT/'figures/chapter3'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8.5)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--master',type=Path,default=Path('A:/Dropbox/LOY_paralogs/summary_table_ddPCR_edit.csv'))
args=parser.parse_args()
source=FIG/'source_data/h2452_ddpcr_initial_guides';source.mkdir(parents=True,exist_ok=True)
master=pd.read_csv(args.master);master['source_row_1based']=np.arange(len(master))+2
all_h=master[master.cell_line.eq('H2452')].copy();assert len(all_h)==33
all_h['condition']=all_h.Guides.fillna('Parental').replace({'CTRL2':'sgCTRL'})
runs=[('04_DEC_2024','4 December 2024'),('06_DEC_2024','6 December 2024')]
mask=all_h.ddPCR_date.isin([r[0] for r in runs]) & all_h.TargetName.eq('AMELY') & all_h.RefGene.eq('AMELX')
d=all_h[mask].copy();assert len(d)==22 and d.ratio.notna().all()
assert d[d.condition.ne('Parental')].Nucelofections.eq(1).all()
all_h['plot_status']=np.where(mask,'included in combined 4 and 6 December plot','other assay date; retained in source context table')
all_h.to_csv(source/'all_h2452_measurements.csv',index=False)
d.to_csv(source/'plotted_measurements.csv',index=False)
all_h[~mask].to_csv(source/'measurements_outside_selected_runs.csv',index=False)
shutil.copy2(args.master,source/'source_summary_table_ddPCR_edit.csv')
fig,axis=plt.subplots(figsize=(120/25.4,105/25.4))
axes=[axis]
fig.subplots_adjust(left=.16,right=.975,bottom=.23,top=.84)
groups=['CP','CN','Q53','sgCTRL','Parental'];colors={g:'#E64B35' for g in groups};colors.update(sgCTRL='#3C5488',Parental='#999999')
summary=[]
for ax in axes:
 run='04_DEC_2024 + 06_DEC_2024';title='H2452 | AMELY';letter=''
 for i,g in enumerate(groups):
  part=d[d.condition.eq(g)].sort_values(['ddPCR_date','source_row_1based']);values=part.ratio.to_numpy()
  if len(values):
   med=float(np.median(values));q1,q3=np.quantile(values,[.25,.75])
   ax.bar(i,med,width=.67,color=colors[g],edgecolor=colors[g],lw=.6,zorder=2)
   ax.errorbar(i,med,yerr=[[med-q1],[q3-med]],fmt='none',ecolor='#333333',elinewidth=.7,capsize=2,zorder=4)
   offsets=np.linspace(-.18,.18,len(values)) if len(values)>1 else [0]
   for date,marker in [('04_DEC_2024','o'),('06_DEC_2024','^')]:
    selected=part.ddPCR_date.eq(date).to_numpy()
    ax.scatter((i+np.array(offsets))[selected],values[selected],s=20,marker=marker,facecolor='white',edgecolor='#222222',lw=.7,zorder=5)
   y=max(values.max(),q3)+.05
   ax.text(i,y,f'{med:.3f}',ha='center',va='bottom',fontsize=8,fontweight='semibold')
   ax.text(i,y+.09,f'n={len(values)}',ha='center',va='bottom',fontsize=7.5,color='#666666')
  else:
   med=q1=q3=np.nan;ax.text(i,.055,'NA',ha='center',va='bottom',fontsize=7,color='#777777')
  summary.append(dict(run=run,condition=g,n_measurements=len(values),DNA_sample_ids=';'.join(dict.fromkeys(part.Sample.astype(str))),median=med,q1=q1,q3=q3,value='recorded summary-table ratio'))
 ax.axvline(2.5,color='#DDDDDD',lw=.65)
 ax.set(xlim=(-.6,4.6),ylim=(0,1.7),yticks=[0,.5,1,1.5])
 ax.set_xticks(range(5),['CP\n1\u00d7','CN\n1\u00d7','Q53\n1\u00d7','sgCTRL','Parental'])
 ax.tick_params(axis='x',length=0,pad=7,labelsize=8)
 ax.set_title(title,fontsize=10,pad=12)
 ax.legend(handles=[Line2D([],[],marker=marker,linestyle='none',markerfacecolor='white',markeredgecolor='#222222',markersize=4,label=label) for marker,label in [('o','Batch 1'),('^','Batch 2')]],loc='upper right',frameon=False,ncol=2,fontsize=7.3,handletextpad=.5,columnspacing=1)
 ax.grid(axis='y',color='#ECECEC',lw=.45);ax.set_axisbelow(True)
axes[0].set_ylabel('Recorded AMELY/AMELX ratio',labelpad=8)
fig.text(.54,.065,'1\u00d7 = one nucleofection',ha='center',fontsize=7,color='#444444')
stem='ch03_main_h2452_ddpcr_initial_guides'
for ext in ['pdf','svg','png']:fig.savefig(FIG/f'{stem}.{ext}',dpi=600,facecolor='white')
fig.savefig(FIG/'previews'/f'{stem}.png',dpi=170,facecolor='white');plt.close(fig)
pd.DataFrame(summary).to_csv(source/'group_summary.csv',index=False)
validation={'master_path':str(args.master),'master_sha256':hashlib.sha256(args.master.read_bytes()).hexdigest(),'n_total_H2452_rows':33,'n_plotted_rows':22,'rows_outside_selected_runs':11,'normalization':'None: exact summary-table ratio field; no copy-number scaling, control normalization or concentration-ratio substitution','aggregation':'median and IQR of all 4 and 6 December measurements within condition, as requested by author; measurement-weighted' ,'n_unit':'assay measurements, including repeated DNA specimens; not independent nucleofections','combined_runs':[x[0] for x in runs],'no_inferential_tests':True,'five_guide_data_included':False,'axis_limits':[0,1.7],'font':'Open Sans','dimensions_mm':[120,105],'historical_ratio_limitation':'Recorded ratio conventions/channel labels need run-specific audit; this display is descriptive and is not a harmonised concentration-ratio analysis','exports':{ext:hashlib.sha256((FIG/f'{stem}.{ext}').read_bytes()).hexdigest() for ext in ['pdf','svg','png']}}
(source/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
print('Exported',stem,'with',len(d),'observations; combined 4 and 6 December')
