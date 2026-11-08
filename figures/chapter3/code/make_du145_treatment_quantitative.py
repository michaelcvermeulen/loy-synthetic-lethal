"""Main and supplementary DU145 treatment figures from verified saved outputs."""
from pathlib import Path
import sys, json, hashlib, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
SRC = FIG / 'source_data/du145_drug_optimisation'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
SAMPLES = ['LNCaP-sgCTRL','LNCaP-DMSO-LOY-1X','LNCaP-30','LNCaP-70','LNCaP-LOY-2X','LNCaP-KU','LNCaP-KU-30','LNCaP-KU-70']
LABEL = dict(zip(SAMPLES,['sgCTRL','1× DMSO','1× pano 30 nM','1× pano 75 nM','2×, no drug','1× KU','1× KU + pano 30 nM','1× KU + pano 75 nM']))
COLOR = dict(zip(SAMPLES,['#3C5488','#777777','#E99B57','#E64B35','#705A8C','#4DBBD5','#00A087','#8491B4']))
MARKER = dict(zip(SAMPLES,['o','o','s','D','^','v','P','X']))
COND = dict(zip(SAMPLES,['sgCTRL','DMSO vehicle','Panobinostat 30 nM','Panobinostat 75 nM','LOY 2x','KU-55933','KU-55933 + panobinostat 30 nM','KU-55933 + panobinostat 75 nM']))

dd=pd.read_csv(SRC/'drug_vs_vehicle_and_sgCTRL.csv')
two_round=pd.read_csv(SRC/'two_round_ddpcr_comparator.csv').set_index('target').to_dict('index')
INPUTS=['drug_vs_vehicle_and_sgCTRL.csv','two_round_ddpcr_comparator.csv','lpwgs/Y_content_summary_depth.csv','lpwgs/Y_content_summary_regions.csv','lpwgs/plotted_chrY_bins.csv','lpwgs/plotted_distal_Yp_summary.csv','panobinostat_grouped_measurements.csv','panobinostat_grouped_summary.csv']
input_hashes={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in INPUTS}

def tag(fig,letter,y):
    fig.text(.025,y,letter,fontsize=12,fontweight='semibold',va='top')

def ddpanel(ax,target,samples,full=False):
    ref='AMELX' if target=='AMELY' else 'RNaseP'
    ys=[]
    for i,s in enumerate(samples):
        row=dd[dd.condition.eq(COND[s]) & dd.target.eq(target)]
        assert len(row)==1
        row=row.iloc[0];y=row.concentration_ratio;ys.append(y)
        ax.vlines(i,0,y,color=COLOR[s],lw=1.2,alpha=.8)
        ax.scatter(i,y,s=25,marker=MARKER[s],facecolor='white',edgecolor=COLOR[s],lw=1.2,zorder=3)
        if not full:
            ax.text(i,y+.004,f'{row.percent_reduction_vs_sgCTRL:.1f}%',ha='center',fontsize=7,color=COLOR[s])
    ax.set(title=f'{target} / {ref}',ylim=(0,.09),yticks=[0,.03,.06,.09],xlim=(-.5,len(samples)-.5))
    ax.set_xticks(range(len(samples)))
    if full:ax.set_xticklabels(['1× DMSO','Pano 30','Pano 75','KU','KU + 30','KU + 75'],rotation=55,ha='right',fontsize=7)
    else:ax.set_xticklabels(['1× DMSO','Pano\n30 nM','Pano\n75 nM'],fontsize=7)
    ax.axhline(ys[0],ls=':',lw=.7,color='#999999',zorder=0)
    ax.grid(axis='y',color='#EEEEEE',lw=.4);ax.set_axisbelow(True)

def save(fig,stem):
    destinations=[(FIG/f'{stem}.{ext}',600) for ext in ['pdf','svg','png']]+[(FIG/'previews'/f'{stem}.png',160)]
    for dest,dpi in destinations:
        tmp=dest.with_name(dest.stem+'.staging'+dest.suffix)
        fig.savefig(tmp,dpi=dpi,facecolor='white')
        for attempt in range(10):
            try:tmp.replace(dest);break
            except OSError:
                if attempt==9:raise
                time.sleep(.5)
    plt.close(fig)

fig=plt.figure(figsize=(183/25.4,265/25.4))
tag(fig,'a',.970)
fig.text(.10,.970,'DU145 ddPCR | all six single-round conditions',fontsize=9,fontweight='semibold',va='top')
ss=[SAMPLES[i] for i in [1,2,3,5,6,7]]
for i,t in enumerate(['AMELY','DDX3Y','KDM5D']):
    ax=fig.add_axes([.10+i*.30,.817,.25,.10]);ddpanel(ax,t,ss,full=True)
    if i==0:ax.set_ylabel('Target/reference ratio')
tag(fig,'b',.737)
for bottom,metric,title in [(.641,'percent_reduction_vs_sgCTRL','Overall target-dosage reduction (%) versus sgCTRL'),
                            (.545,'percent_reduction_vs_vehicle','Additional reduction of remaining signal (%) versus 1X DMSO')]:
    overall = metric == 'percent_reduction_vs_sgCTRL'
    cmap = plt.get_cmap('magma_r' if overall else 'RdBu_r')
    norm = matplotlib.colors.Normalize(vmin=0 if overall else -100, vmax=100)
    ax=fig.add_axes([.10,bottom,.77,.061]);ax.axis('off')
    fig.text(.10,bottom+.074,title,fontsize=7.5,fontweight='semibold')
    table=[]
    for t in ['AMELY','DDX3Y','KDM5D']:
        values=[f'{dd.loc[dd.condition.eq(COND[s]) & dd.target.eq(t),metric].iloc[0]:.1f}' for s in ss]
        if overall:
            values.append(f'{two_round[t]["reduction"]:.1f}' if t in two_round else 'Not measured')
        table.append([t]+values)
    headers=['Target','1X DMSO','Pano 30','Pano 75','KU','KU + 30','KU + 75']
    if overall: headers.append('2X*')
    tab=ax.table(cellText=table,colLabels=headers,cellLoc='center',bbox=[0,0,1,1])
    tab.auto_set_font_size(False);tab.set_fontsize(6.6)
    for (r,c),cell in tab.get_celld().items():
        cell.set_edgecolor('#DDDDDD');cell.set_linewidth(.4)
        if r==0:cell.set_facecolor('#F1F1F1')
        elif c>0:
            try:
                value=float(cell.get_text().get_text())
            except ValueError:
                cell.set_facecolor('#E7E7E7');cell.get_text().set_fontsize(5.8)
                continue
            rgba=cmap(norm(value));cell.set_facecolor(rgba)
            luminance=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
            cell.get_text().set_color('white' if luminance<.45 else '#202020')
    cax=fig.add_axes([.90,bottom,.012,.061])
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,
                    ticks=[0,50,100] if overall else [-100,0,100])
    cb.ax.tick_params(labelsize=6,pad=2,length=2);cb.outline.set_linewidth(.4)
    cb.ax.set_title('%',fontsize=6,pad=3)
fig.text(.10,.744,'All drug conditions: 1X pooled CP + Q53 + P47. *2X: no drug, separate February assays.',fontsize=6.3)
fig.text(.10,.013,'March: same-run sgCTRL / DMSO. 2X: assay-matched sgCTRL; no cross-run comparison with March DMSO.',fontsize=6.3,color='#666666')

ORDER=[SAMPLES[i] for i in [0,1,2,3,5,6,7,4]]
LABELS=['sgCTRL','1\u00d7 DMSO','Pano 30','Pano 75','KU','KU + 30','KU + 75','2\u00d7, no drug']
depth=pd.read_csv(SRC/'lpwgs/Y_content_summary_depth.csv').set_index('sample').loc[ORDER]
regions=pd.read_csv(SRC/'lpwgs/Y_content_summary_regions.csv')
bins=pd.read_csv(SRC/'lpwgs/plotted_chrY_bins.csv')
rows=[]
for sample in ORDER:
 for region,lo,hi in [('Distal Yp',2.46,4.8),('Proximal Yp',6.,10.),('Euchromatic Yq',10.883085,27.449937)]:
  vals=bins.loc[bins['sample'].eq(sample)&bins.midpoint.ge(lo)&bins.midpoint.lt(hi),'relative_cn']
  rows.append(dict(sample=sample,region=region,start_Mb=lo,end_Mb=hi,n_bins=len(vals),median_relative_cn=vals.median()))
regional=pd.DataFrame(rows)
previous=pd.read_csv(SRC/'lpwgs/plotted_distal_Yp_summary.csv')
for row in previous.itertuples():
 value=regional.loc[regional['sample'].eq(row.sample)&regional.region.eq('Distal Yp'),'median_relative_cn'].iloc[0]
 assert np.isclose(value,row.median_relative_cn,rtol=1e-12)
for row in regional[regional.region.eq('Euchromatic Yq')].itertuples():
 value=regions.loc[regions['sample'].eq(row.sample)&regions.region.eq('MSY Yq'),'median_relative_cn'].iloc[0]
 assert np.isclose(value,row.median_relative_cn,rtol=1e-12)
fig.text(.54,.503,'d',fontsize=12,fontweight='semibold',va='top')
fig.text(.60,.498,'Depletion by panobinostat exposure',fontsize=7.7,fontweight='semibold')
arm_rows=[]
for region,bottom,title in [('Whole Y',.365,'Whole-Y coverage'),('MSY Yp',.245,'Male-specific Yp'),('MSY Yq',.125,'Euchromatic Yq')]:
 arm=regions[regions.region.eq(region)].set_index('sample')
 control=float(depth.loc['LNCaP-sgCTRL','chrY_autosome_ratio']) if region=='Whole Y' else float(arm.loc['LNCaP-sgCTRL','median_relative_cn'])
 ax=fig.add_axes([.625,bottom,.295,.090])
 for group,samples in enumerate([[SAMPLES[1],SAMPLES[5]],[SAMPLES[2],SAMPLES[3],SAMPLES[6],SAMPLES[7]]]):
  vals=[]
  for sample,offset in zip(samples,np.linspace(-.12,.12,len(samples))):
   if region=='Whole Y':
    value=float(depth.loc[sample,'chrY_autosome_ratio']);lo=hi=np.nan
    prior=pd.read_csv(SRC/'panobinostat_grouped_measurements.csv')
    ref=prior.loc[prior.measure.eq('lpWGS')&prior.condition.eq(COND[sample])].iloc[0]
    assert np.isclose(value,ref.ratio,rtol=1e-12)
    assert np.isclose(100*value/control,ref.residual_percent,rtol=1e-12)
   else:
    row=arm.loc[sample];value=float(row.median_relative_cn);lo=row.start_Mb;hi=row.end_Mb
    check=bins.loc[bins['sample'].eq(sample)&bins.midpoint.ge(lo)&bins.midpoint.lt(hi),'relative_cn'].median()
    assert np.isclose(value,check,rtol=1e-12)
   depletion=100*(1-value/control);vals.append(depletion)
   arm_rows.append(dict(sample=sample,condition=COND[sample],region=region,panobinostat=bool(group),relative_cn=value,sgCTRL_relative_cn=control,depletion_percent=depletion,start_Mb=lo,end_Mb=hi))
   ax.scatter(group+offset,depletion,s=28,marker=MARKER[sample],color=COLOR[sample],edgecolors='white',linewidths=.4,zorder=3)
  med=float(np.median(vals));ax.plot([group-.20,group+.20],[med,med],color='#333333',lw=1.1)
  ax.text(group+.24,med,f'{med:.1f}%',fontsize=6.6,va='center')
 ax.set(xlim=(-.4,1.65),ylim=(75,100),yticks=[75,80,85,90,95,100],xticks=[0,1])
 ax.set_title(title,fontsize=8,fontweight='semibold',pad=5)
 ax.set_ylabel('Depletion (%)',fontsize=7)
 ax.set_xticklabels(['No panobinostat\nn = 2','Panobinostat\nn = 4'] if region=='MSY Yq' else [])
 ax.tick_params(labelsize=6.7);ax.grid(axis='y',color='#eeeeee',lw=.4);ax.set_axisbelow(True)
handles=[Line2D([],[],marker=MARKER[SAMPLES[i]],color=COLOR[SAMPLES[i]],ls='',markersize=4,label=label) for i,label in [(1,'DMSO'),(5,'KU'),(2,'Pano 30'),(3,'Pano 75'),(6,'KU + 30'),(7,'KU + 75')]]
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.77,.069),ncol=3,frameon=False,fontsize=6.1,columnspacing=.8,handletextpad=.3)
tag(fig,'c',.503)
fig.text(.24,.498,'Regional Y copy number',fontsize=9,fontweight='semibold')
heat=fig.add_axes([.24,.105,.23,.35])
values=regional.pivot(index='sample',columns='region',values='median_relative_cn').loc[ORDER,['Distal Yp','Proximal Yp','Euchromatic Yq']].values
norm=matplotlib.colors.LogNorm(vmin=.001,vmax=1)
im=heat.imshow(values,cmap='magma',norm=norm,aspect='auto',interpolation='nearest')
heat.set(ylim=(7.6,-.6),yticks=range(8),yticklabels=LABELS,xticks=range(3),xticklabels=['Distal\nYp','Proximal\nYp','Euchromatic\nYq'])
heat.tick_params(axis='x',length=0,labelsize=6.8,pad=5)
heat.tick_params(axis='y',length=0,pad=5,labelsize=7)
for i,label in enumerate(heat.get_yticklabels()):label.set_color(COLOR[ORDER[i]])
heat.spines[['top','right','left','bottom']].set_visible(False)
for i in range(8):
 for j in range(3):
  v=values[i,j];label='Floor' if v<=.0010000001 else f'{v:.3f}'
  heat.text(j,i,label,ha='center',va='center',fontsize=7.5,color='#222222' if norm(v)>.62 else 'white')
for y in np.arange(.5,7.5):heat.axhline(y,color='white',lw=.7)
for x in [.5,1.5]:heat.axvline(x,color='white',lw=.7)
cax=fig.add_axes([.49,.14,.012,.28]);cb=fig.colorbar(im,cax=cax,ticks=[.001,.01,.1,1])
cb.ax.set_yticklabels(['0.001','0.01','0.1','1']);cb.ax.tick_params(labelsize=6)
cb.set_label('Median relative copy number',fontsize=7,labelpad=3)
fig.text(.10,.047,'WGS: one sample per condition. Floor: numerical analysis floor (0.001).',fontsize=6.5,color='#555555')
stem='ch03_main_du145_drug_treatment_comparison'
save(fig,stem)
out=SRC/'quantitative_comparison';out.mkdir(exist_ok=True)
regional.to_csv(out/'regional_values.csv',index=False)
depth.to_csv(out/'whole_Y_values.csv')
pd.DataFrame(arm_rows).to_csv(out/'grouped_arm_depletion.csv',index=False)
pd.DataFrame(arm_rows).groupby(['region','panobinostat']).depletion_percent.agg(['count','median']).to_csv(out/'grouped_arm_summary.csv')
assert all(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h for f,h in input_hashes.items())
(out/'validation.json').write_text(json.dumps({'figure':stem,'font':'Open Sans','size_mm':[183,265],'source_hashes':input_hashes,'source_inputs_unchanged':True,'ddPCR':'Saved ratios and percent reductions, unchanged','distal_Yp_matches_previous':True,'euchromatic_Yq_matches_previous':True,'regional_intervals_Mb':{'distal_Yp':[2.46,4.8],'proximal_Yp':[6,10],'euchromatic_Yq':[10.883085,27.449937]},'biological_samples_per_WGS_condition':1,'grouped_arm_values_match_saved_bins':True,'panel_d':'Stacked whole-Y, male-specific Yp and euchromatic Yq depletion; each metric normalised to its sgCTRL value','whole_Y_matches_saved_grouped_values':True,'grouped_WGS_n':[2,4],'inferential_tests':False},indent=2),encoding='utf-8')
print('Rendered quantitative DU145 treatment comparison; saved inputs unchanged.')
