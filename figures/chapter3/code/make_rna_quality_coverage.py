"""Quality, trimming and coverage companion; frozen QC values, no reanalysis."""
from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
OUT = ROOT / 'figures/chapter3'
SRC = OUT / 'source_data/rna_quality_coverage'
d = pd.read_csv(SRC/'sample_metrics.tsv', sep='\t')
cov = pd.read_csv(SRC/'genebody_coverage.tsv', sep='\t')
provider=pd.read_csv(OUT/'source_data/provider_rna_qc/provider_qc.tsv',sep='\t')
assert set(provider.qc_id)==set(d.qc_id) and provider.qc_id.is_unique
d=d.merge(provider[['qc_id','integrity_type','integrity','size_metric','size_percent','mapping_status']],on='qc_id',validate='one_to_one',sort=False)
d.to_csv(SRC/'plotted_metrics_with_provider.tsv',sep='\t',index=False)
assert len(d) == 44 and d.qc_id.is_unique and d.plot_x.is_unique
assert set(cov.qc_id) == set(d.qc_id)
assert not cov.duplicated(['qc_id','pos']).any()
for _, t in cov.groupby('qc_id'):
    assert set(t.pos) == set(range(100))
    assert np.isfinite(t.cov_rel).all() and (t.cov_rel >= 0).all()
    assert np.isclose(t.cov_rel.mean(), 1, atol=1e-10)
for col in ['q30_raw','q30_trim','adapter_trimmed_pct','fastp_pass_pct']:
    assert np.isfinite(d[col]).all() and d[col].between(0,100).all()
assert np.allclose(d.fastp_pass_pct, 100*d.fastp_pass_reads/d.raw_reads)
cov = cov.merge(d[['qc_id','group','display_condition','display_label']],on='qc_id',validate='many_to_one')
cov.to_csv(SRC/'plotted_coverage.tsv',sep='\t',index=False)
colors = {'PARENTAL':'#8C8C8C','CTRL':'#3C5488','sgCTRL':'#3C5488','LOY':'#E64B35',
          'LOY_LXHA':'#C77C52','LOY_UTY_AB':'#00A087','LOY_KDM5D_AB':'#8491B4','LOY_KDM5D_UTY_AB':'#7E6148'}
groups = ['A549','HT1080','HT1080 ADDBACK','DU145']
group_labels = ['A549 · Aug 2025','HT1080 · Aug 2025','HT1080 · May 2026','DU145 · May 2026']
x = d.plot_x.to_numpy()
fig = plt.figure(figsize=(8.3,12.2))
gs = fig.add_gridspec(6,4,left=.105,right=.98,top=.925,bottom=.14,hspace=.95,wspace=.22,height_ratios=[1,1,1,1.2,1,1])
axs = [fig.add_subplot(gs[i,:]) for i in range(3)]
for i,ax in enumerate(axs):
    ax.set_xlim(-1,x[-1]+1)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True,color='#EAEAEA',lw=.4)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis='x',length=0)
    ax.set_xticks(x, d.display_label if i==2 else ['']*44,rotation=90,fontsize=7)
    ax.text(-.09,1.10,chr(97+i),transform=ax.transAxes,fontweight='semibold',fontsize=11)
    for g in ['HT1080','DU145']:
        xx=d.loc[d['group']==g,'plot_x'];ax.axvspan(xx.min()-.5,xx.max()+.5,color='#F3F4F5',zorder=0)
for g,label in zip(groups,group_labels):
    xx=d.loc[d['group']==g,'plot_x']
    axs[0].text((xx.min()+xx.max())/2,1.65,label,transform=axs[0].get_xaxis_transform(),ha='center',fontsize=8,fontweight='semibold')
axs[0].set_title('Base quality before and after filtering',loc='left',pad=8)
axs[0].vlines(x,d.q30_raw,d.q30_trim,color='#999999',lw=.7)
axs[0].scatter(x,d.q30_raw,marker='o',s=18,facecolors='white',edgecolors='#444444',lw=.7,label='Before',zorder=3)
axs[0].scatter(x,d.q30_trim,marker='^',s=20,color='#3C5488',label='After',zorder=3)
axs[0].set_ylabel('Q30 bases (%)');axs[0].set_ylim(92,100)
axs[0].legend(loc='lower right',bbox_to_anchor=(1,1.01),ncol=2,fontsize=7)
axs[1].set_title('Reads trimmed for adapters',loc='left',pad=8)
axs[1].bar(x,d.adapter_trimmed_pct,color=d.display_condition.map(colors),width=.78)
axs[1].set_ylabel('Raw reads (%)');axs[1].set_ylim(0,35)
axs[2].set_title('Read retention after filtering',loc='left',pad=8)
axs[2].scatter(x,d.fastp_pass_pct,c=d.display_condition.map(colors),s=19,zorder=3)
axs[2].set_ylabel('Raw reads retained (%)');axs[2].set_ylim(94,100)
coverage_axes=[]
ymax=np.ceil(cov.cov_rel.max()*10)/10
for j,(g,label) in enumerate(zip(groups,group_labels)):
    ax=fig.add_subplot(gs[3,j]);coverage_axes.append(ax)
    for r in d.loc[d['group']==g].itertuples():
        t=cov.loc[cov.qc_id==r.qc_id].sort_values('pos')
        ax.plot(t.pos,t.cov_rel,color=colors[r.display_condition],lw=.65,alpha=.75)
    ax.axhline(1,color='#888888',lw=.5,ls='--')
    ax.set_xlim(0,99);ax.set_ylim(0,ymax)
    ax.set_xticks([0,49.5,99],['5′','Mid','3′'])
    ax.set_title(label.replace(' · ','\n'),fontsize=8,pad=6)
    if j==0:ax.set_ylabel('Coverage / library mean')
    else:ax.tick_params(labelleft=False)
coverage_axes[0].text(-.42,1.43,'d',transform=coverage_axes[0].transAxes,fontsize=11,fontweight='semibold')
coverage_axes[0].text(0,1.43,'Gene-body coverage',transform=coverage_axes[0].transAxes,fontweight='semibold')
legend_labels={'PARENTAL':'Parental','CTRL':'sgCTRL','LOY':'Y-depleted','LOY_LXHA':'Empty vector','LOY_UTY_AB':'UTY','LOY_KDM5D_AB':'KDM5D','LOY_KDM5D_UTY_AB':'UTY + KDM5D'}
for row,key,title,ylabel,limits in [(4,'integrity','Input RNA integrity','Integrity score',(7.5,10.5)),(5,'size_percent','Input RNA size distribution','RNA fraction (%)',(75,100))]:
    ax=fig.add_subplot(gs[row,:]);ax.set_xlim(-1,x[-1]+1);ax.set_ylim(*limits)
    ax.set_title(title,loc='left',pad=24);ax.set_ylabel(ylabel)
    ax.text(-.09,1.10,chr(97+row),transform=ax.transAxes,fontweight='semibold',fontsize=11)
    ax.spines['bottom'].set_visible(False);ax.tick_params(axis='x',length=0)
    ax.set_xticks(x,d.display_label if row==5 else ['']*44,rotation=90,fontsize=7)
    ax.set_axisbelow(True);ax.yaxis.grid(True,color='#EAEAEA',lw=.4)
    for g in groups:
        t=d.loc[d['group']==g]
        ax.scatter(t.plot_x,t[key],s=20,c=t.display_condition.map(colors),marker='o',zorder=3)
        if g in ['HT1080','DU145']:ax.axvspan(t.plot_x.min()-.5,t.plot_x.max()+.5,color='#F3F4F5',zorder=0)
        label=('RQN' if g=='HT1080 ADDBACK' else 'RIN (reported)' if g=='DU145' else 'RINe') if row==4 else ('200–10,000 nt' if g=='HT1080 ADDBACK' else 'DV200')
        ax.text(t.plot_x.mean(),1.03,label,transform=ax.get_xaxis_transform(),ha='center',fontsize=7)
fig.legend([Line2D([0],[0],color=colors[k],lw=2) for k in legend_labels],list(legend_labels.values()),loc='lower center',bbox_to_anchor=(.54,.075),ncol=4,fontsize=7.5,columnspacing=1.3,handlelength=1.3)
fig.text(.105,.046,'P, parental; C, sgCTRL; Yd, Y-depleted; EV, empty vector; U, UTY; K, KDM5D; UK, UTY + KDM5D.',fontsize=7)
stem='ch03_supp_rna_quality_coverage_draft'
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{stem}.{ext}',dpi=300)
(SRC/'validation.json').write_text(json.dumps({'libraries':len(d),'coverage_rows':len(cov),'positions_per_library':100,'coverage_mean_per_library':1,'percentage_bounds_passed':True,'read_retention_recalculated':True,'sample_order':'identical to first QC figure','font':'Open Sans'},indent=2),encoding='utf-8')
print('Saved',stem,'with 44 libraries and 4400 coverage values.')
