"""Render the 44-library QC draft from frozen archive tables; no DE reanalysis."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
OUT = ROOT / 'figures/chapter3'
SRC = OUT / 'source_data/rna_all_samples_qc'
d = pd.read_csv(SRC / 'qc_metrics.tsv', sep='\t')
m = pd.read_csv(SRC / 'addback_sample_manifest.csv')
lookup = {(r.original_label, str(r.collection_occasion)): r.final_group for r in m.itertuples()}
names = {'CTRL':'C', 'sgCTRL':'C', 'PARENTAL':'P', 'LOY':'Yd',
         'LOY_LXHA':'EV', 'LOY_UTY_AB':'U', 'LOY_KDM5D_AB':'K', 'LOY_KDM5D_UTY_AB':'UK'}
groups = ['A549', 'HT1080', 'HT1080 ADDBACK', 'DU145']
order = ['PARENTAL','CTRL','sgCTRL','LOY','LOY_LXHA','LOY_UTY_AB','LOY_KDM5D_AB','LOY_KDM5D_UTY_AB']
d['display_condition'] = d.condition
for i,r in d.iterrows():
    if r['group'] == 'HT1080 ADDBACK':
        d.loc[i,'display_condition'] = lookup[(r.condition, r.qc_id.rsplit('-',1)[1])]
d['group_order'] = d['group'].map({g:i for i,g in enumerate(groups)})
d['condition_order'] = d.display_condition.map({g:i for i,g in enumerate(order)})
d = d.sort_values(['group_order','condition_order','qc_id']).reset_index(drop=True)
assert len(d)==44 and d.qc_id.is_unique and d.display_condition.notna().all()
colors = {'PARENTAL':'#8C8C8C','CTRL':'#3C5488','sgCTRL':'#3C5488','LOY':'#E64B35',
          'LOY_LXHA':'#C77C52','LOY_UTY_AB':'#00A087','LOY_KDM5D_AB':'#8491B4','LOY_KDM5D_UTY_AB':'#7E6148'}
x=[]; spans=[]; cursor=0
for g in groups:
    n=(d['group']==g).sum(); x.extend(range(cursor,cursor+n));spans.append((g,cursor,cursor+n-1));cursor+=n+1
x=np.array(x)
d['plot_x']=x
d['display_label']=[names[c]+q.rsplit('-',1)[1] for c,q in zip(d.display_condition,d.qc_id)]
d.to_csv(SRC/'plotted_metrics.tsv',sep='\t',index=False)
fig,axs=plt.subplots(6,1,figsize=(8.3,9.4),sharex=True)
fig.subplots_adjust(left=.10,right=.99,top=.90,bottom=.13,hspace=.55)
titles=['Sequencing depth','Alignment outcome','Genomic read distribution','Sequence duplication','Mitochondrial and rRNA abundance','Gene-body coverage bias']
for i,(ax,title) in enumerate(zip(axs,titles)):
    ax.set_title(title,loc='left',pad=7)
    ax.text(-.085,1.07,chr(97+i),transform=ax.transAxes,fontweight='semibold',fontsize=11)
    for g,l,r in spans:
        if g in ['HT1080','DU145']:ax.axvspan(l-.5,r+.5,color='#F3F4F5',zorder=0)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis='x',length=0)
    ax.set_xlim(-1,x[-1]+1)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True,color='#EAEAEA',lw=.4)
axs[0].bar(x,d.star_input_pairs/1e6,color=d.display_condition.map(colors),width=.78)
axs[0].set_ylabel('Read pairs (M)');axs[0].set_ylim(0,52)
for g,l,r in spans:
    label={'A549':'A549 · Aug 2025','HT1080':'HT1080 · Aug 2025','HT1080 ADDBACK':'HT1080 · May 2026','DU145':'DU145 · May 2026'}[g]
    axs[0].text((l+r)/2,1.55,label,ha='center',transform=axs[0].get_xaxis_transform(),fontsize=8,fontweight='semibold')
bottom=np.zeros(len(d))
for v,label,c in [(d.star_unique_pct,'Unique','#3C5488'),(d.star_multi_pct,'Multiple loci','#91B8B5'),(100-d.star_unique_pct-d.star_multi_pct,'Unmapped','#D5D5D5')]:
    assert np.isfinite(v).all() and (v>=-1e-5).all()
    axs[1].bar(x,v,bottom=bottom,width=.78,color=c,label=label);bottom+=v
axs[1].set_ylabel('Input pairs (%)');axs[1].set_ylim(0,100)
bottom=np.zeros(len(d))
for key,label,c in [('qm_exonic_pct','Exonic','#3C5488'),('qm_intronic_pct','Intronic','#E9A36F'),('qm_intergenic_pct','Intergenic','#B8B3C9')]:
    axs[2].bar(x,d[key],bottom=bottom,width=.78,color=c,label=label);bottom+=d[key]
assert np.allclose(bottom,100)
axs[2].set_ylabel('Classified reads (%)');axs[2].set_ylim(0,100)
axs[3].bar(x,d.dup_pct,color=d.display_condition.map(colors),width=.78)
axs[3].set_ylabel('Duplication (%)');axs[3].set_ylim(0,40)
axs[4].scatter(x-.12,d.chrM_pct,s=13,color='#3C5488',label='chrM alignments',zorder=3)
axs[4].scatter(x+.12,100*d.pic_ribosomal,s=15,marker='D',color='#E69F00',label='rRNA bases',zorder=3)
axs[4].set_ylabel('Abundance (%)');axs[4].set_ylim(0,20)
axs[5].scatter(x,d.qm_bias_53,s=17,c=d.display_condition.map(colors),zorder=3)
axs[5].axhline(1,color='#888888',lw=.6,ls='--')
axs[5].set_ylabel("5′ / 3′ coverage");axs[5].set_ylim(.95,1.4)
for i in [1,2,4]:axs[i].legend(loc='lower right',bbox_to_anchor=(1,1.01),ncol=3,fontsize=7,handlelength=1.1,columnspacing=1)
axs[-1].set_xticks(x,d.display_label,rotation=90,fontsize=7)
fig.text(.10,.055,'P, parental; C, sgCTRL; Yd, Y-depleted; EV, empty vector; U, UTY; K, KDM5D; UK, UTY + KDM5D.',fontsize=7.5)
fig.text(.10,.033,'One position per library across all panels; add-back labels follow the corrected sample manifest.',fontsize=7.5,color='#555555')
stem='ch03_supp_rna_all_samples_qc_draft'
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{stem}.{ext}',dpi=300)
print('Rendered',len(d),'libraries:',OUT/stem)
