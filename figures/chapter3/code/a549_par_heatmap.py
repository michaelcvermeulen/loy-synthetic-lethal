"""Compact distal-Yp schematic and sample-level PAR/MSY expression heatmap."""
from pathlib import Path
import sys,json,hashlib
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'figures/chapter3'
DATA=OUT/'source_data/a549_rna_thesis'
SOURCE=ROOT.parent/'RNAseq/analysis/A549/SEPT_18_2026_REANALYZE/outputs/tables'

def draw(par,xy,effects,par_de):
    genes=['PLCXD1','GTPBP6','PPP2R3B','CRLF2','CSF2RA','SLC25A6','ASMTL','AKAP17A','DHRSX','ZBED1','CD99','RPS4Y1','ZFY']
    data=pd.concat([par[par.PAR_region=='PAR1'],xy[xy.gene_name.isin(['RPS4Y1','ZFY'])]],ignore_index=True)
    cols=[f'A549-{g}-{i}' for g in ['CTRL','LOY'] for i in range(1,5)]
    cp=data.pivot(index='gene_name',columns='sample_key',values='CPM').loc[genes,cols]
    assert cp.shape==(13,8) and cp.notna().all().all()
    rel=cp.div(cp.iloc[:,:4].mean(axis=1),axis=0)*100
    means=rel.iloc[:,4:].mean(axis=1)
    top=max(150,int(np.ceil(rel.to_numpy().max()/50)*50))
    fig=plt.figure(figsize=(183/25.4,98/25.4))
    fig.text(.025,.95,'a',fontsize=12,weight='semibold')
    fig.text(.07,.95,'Distal Yp and PAR1 expression',fontsize=8,weight='semibold')
    fig.text(.775,.95,'b',fontsize=12,weight='semibold')
    fig.text(.82,.95,'X counterparts',fontsize=8,weight='semibold')
    y0=.24;h=.57
    ax=fig.add_axes([.28,y0,.275,h])
    im=ax.imshow(rel,cmap='magma',vmin=0,vmax=top,aspect='auto',interpolation='nearest')
    ax.set_yticks(range(13),genes,fontsize=6.6,fontstyle='italic')
    ax.set_xticks(range(8),[1,2,3,4]*2,fontsize=6.5)
    ax.tick_params(length=0,pad=3);ax.spines[:].set_visible(False)
    ax.axvline(3.5,color='white',lw=2);ax.axhline(10.5,color='white',lw=2)
    for x0,colour,label in [(-.5,'#3C5488','sgCTRL'),(3.5,'#E64B35','Y-depleted')]:
        ax.add_patch(Rectangle((x0,-1.02),4,.3,facecolor=colour,clip_on=False,lw=0))
        ax.text(x0+2,-1.45,label,ha='center',va='bottom',fontsize=7,color=colour)
    ax.set_xlabel('Sample',fontsize=7,labelpad=3)
    schem=fig.add_axes([.038,y0,.083,h]);schem.set(xlim=(0,1),ylim=(12.5,-.5));schem.axis('off')
    # Only the distal Yp segment is illustrated, not a full chromosome/copy-number track.
    shape=FancyBboxPatch((.47,-.45),.29,12.9,boxstyle='round,pad=0,rounding_size=.18',facecolor='#ECECEC',edgecolor='#777777',lw=.6)
    schem.add_patch(shape)
    parblock=Rectangle((.47,-.45),.29,10.95,facecolor='#B8DEDA',edgecolor='none');parblock.set_clip_path(shape);schem.add_patch(parblock)
    schem.plot([.47,.76],[10.5,10.5],color='#777777',lw=.6)
    schem.text(.3,4.8,'PAR1',rotation=90,ha='center',va='center',fontsize=7,color='#38766E')
    schem.text(.3,11.55,'MSY',rotation=90,ha='center',va='center',fontsize=6.4,color='#555555')
    for i in range(13):schem.plot([.77,1.0],[i,i],color='#999999',lw=.4)
    schem.text(.61,-1.05,'Yp tip',ha='center',va='bottom',fontsize=6.5)
    schem.annotate('',xy=(.61,13.5),xytext=(.61,12.5),arrowprops={'arrowstyle':'->','lw':.6,'color':'#777777'},annotation_clip=False)
    stats=fig.add_axes([.575,y0,.15,h]);stats.set(xlim=(-.5,1.5),ylim=(12.5,-.5));stats.axis('off')
    lookup=par_de.set_index('gene_name')
    for i,gene in enumerate(genes):
        stats.text(0,i,f'{means[gene]:.1f}%',ha='center',va='center',fontsize=6.5)
        val=lookup.loc[gene,'BH_PAR_family'] if gene in lookup.index else np.nan
        stats.text(1,i,f'{val:.3f}' if np.isfinite(val) else '—',ha='center',va='center',fontsize=6.3,weight='semibold' if val<.05 else 'normal')
    stats.text(0,-1.4,'Mean\nY-depleted',ha='center',va='bottom',fontsize=6.2)
    stats.text(1,-1.4,'PAR\nadj. P',ha='center',va='bottom',fontsize=6.2)
    c=fig.add_axes([.845,y0,.13,h])
    c.axvline(0,color='.65',lw=.7);c.scatter(effects.X_mean_difference,range(8),color='#333333',s=19)
    c.set_yticks(range(8),effects.X,fontstyle='italic',fontsize=6.5);c.invert_yaxis()
    c.set(xlim=(-.35,.15),xlabel='Mean log₂CPM\ndifference')
    c.xaxis.label.set_size(7);c.tick_params(axis='x',labelsize=6.5);c.set_xticks([-.3,-.1,.1])
    cb=fig.colorbar(im,cax=fig.add_axes([.28,.115,.275,.025]),orientation='horizontal')
    cb.ax.tick_params(labelsize=6.5);cb.set_label('Expression (% of mean sgCTRL)',fontsize=7,labelpad=2)
    fig.text(.04,.15,'Distal Yp\nNot to scale',fontsize=6.2,color='#555555')
    table=rel.copy();table['LOY_mean_percent']=means;table['PAR_BH']=par_de.set_index('gene_name').BH_PAR_family.reindex(genes)
    table.to_csv(DATA/'distal_Yp_PAR_heatmap.csv')
    return fig

if __name__=='__main__':
    sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
    from thesis_style import apply_style
    apply_style(8)
    args=[pd.read_csv(SOURCE/(name+'.csv')) for name in ['PAR_selected_gene_sample_abundance','XY_paralog_sample_abundance','XY_paralog_effects_and_length_proxy','DE_PAR_secondary_voom']]
    fig=draw(*args)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/('ch03_main_a549_rna_sex_chromosomes.'+ext),dpi=300)
    print('Saved compact PAR1/distal-Yp heatmap and X counterparts.')
