"""Six separate chromosome-Y lpWGS profiles on one shared linear scale."""
from pathlib import Path
import argparse, hashlib, json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path,default=FIG)
args=parser.parse_args()
OUT=args.output_dir.resolve();OUT.mkdir(parents=True,exist_ok=True)
DATA=OUT/'source_data/regional_y_atlases_normalised'
STEM='ch03_supp_regional_y_individual_tracks'
COL={'sgCTRL':'#3C5488','LOY':'#E64B35'}
LABEL={'sgCTRL':'sgCTRL','LOY':'Y-depleted'}
P1,C0,C1,Q0,P2,LENGTH=2.458320,10.565750,10.883085,27.449937,62.122809,62.460029
UPPER=1.3
fig=plt.figure(figsize=(183/25.4,224/25.4))
gs=fig.add_gridspec(6,1,left=.115,right=.955,top=.94,bottom=.075,hspace=.36)
audits=[];profiles=[];hashes={}
for i,cell in enumerate(['ht1080','a549','du145']):
    path=DATA/cell/'paired_Y_bins_and_exclusions.csv'
    hashes[str(path.relative_to(OUT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    pairs=pd.read_csv(path);d=pairs.loc[pairs.included].copy()
    regions=[('PAR1',0,P1),('MSY Yp',P1,C0),('euchromatic Yq',C1,Q0),('PAR2',P2,LENGTH)]
    if cell=='du145':regions[1:2]=[('distal MSY Yp',P1,4.8),('proximal MSY Yp',4.8,C0)]
    for j,condition in enumerate(['sgCTRL','LOY']):
        index=2*i+j;ax=fig.add_subplot(gs[index])
        values=np.exp2(d['log2_'+condition])
        for a,b in [(0,P1),(P2,LENGTH)]:ax.axvspan(a,b,color='#EAF6F7',zorder=0)
        for a,b in [(C0,C1),(Q0,P2)]:
            ax.axvspan(a,b,facecolor='#F0F0F0',edgecolor='#DDDDDD',hatch='///',lw=.4,zorder=0)
        # All included bins are represented; above-range values have explicit triangles.
        normal=values.le(UPPER)
        ax.scatter(d.loc[normal,'midpoint'],values[normal],s=5.5,color=COL[condition],alpha=.45,lw=0,zorder=2)
        ax.scatter(d.loc[~normal,'midpoint'],np.full((~normal).sum(),UPPER),s=12,
                   marker='^',color=COL[condition],lw=0,zorder=3)
        for region,a,b in regions:
            v=values[d.region.eq(region)]
            ax.plot([a,b],[v.median()]*2,color=COL[condition],lw=1.5,solid_capstyle='butt',zorder=4)
        ax.text(0,1.055,f'{cell.upper()}  |  {LABEL[condition]}',transform=ax.transAxes,
                ha='left',va='bottom',color=COL[condition],fontsize=8.5,fontweight='semibold')
        ax.text(-.10,1.14,chr(97+index),transform=ax.transAxes,ha='left',va='top',
                fontsize=12,fontweight='semibold')
        ax.text((Q0+P2)/2,.64,'Yq12 heterochromatin (unmappable, excluded)',
                ha='center',va='center',fontsize=7.1,color='#666666',
                bbox=dict(facecolor='#F0F0F0',edgecolor='none',alpha=.9,pad=2))
        if (~normal).any():
            ax.text((Q0+P2)/2,1.12,f'{(~normal).sum()} bins >1.3; maximum {values.max():.2f}',
                    fontsize=7,ha='center',va='center',color=COL[condition],
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.5))
        ax.set(xlim=(0,LENGTH),ylim=(-.025,1.37),yticks=[0,.5,1],ylabel='Relative CN')
        ax.set_xticks([0,10,20,30,40,50,60])
        ax.spines[['right','top']].set_visible(False)
        ax.tick_params(axis='x',labelbottom=index==5,length=2.5)
        if index==5:ax.set_xlabel('HG002 T2T-Y position (Mb; full chromosome scale)',fontsize=8)
        if index==0:
            ax.annotate('Centromere excluded',xy=((C0+C1)/2,.82),xytext=(12.7,1.10),
                        fontsize=6.6,color='#666666',arrowprops=dict(arrowstyle='-',lw=.5,color='#999999'))
        profiles.append(pd.DataFrame({'cell_line':cell.upper(),'condition':condition,
            'bin_id':d.bin_id,'midpoint_Mb':d.midpoint,'region':d.region,
            'corrected_relative_CN':values,'above_display_range':~normal}))
        audits.append({'panel':chr(97+index),'cell_line':cell.upper(),'condition':condition,
                       'n_bins':len(d),'n_above_1_3':int((~normal).sum()),'maximum':float(values.max())})
handles=[Line2D([],[],marker='o',ls='none',color='#888888',markersize=3,label='100-kb bin'),
         Line2D([],[],color='#666666',lw=1.5,label='Regional median'),
         Line2D([],[],marker='^',ls='none',color='#666666',markersize=4,label='Above displayed range')]
fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.955,.995),ncol=3,
           frameon=False,fontsize=7,handlelength=1.5,columnspacing=1.2)
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600,facecolor='white')
(OUT/'previews').mkdir(exist_ok=True)
fig.savefig(OUT/'previews'/f'{STEM}.png',dpi=160,facecolor='white');plt.close(fig)
pd.concat(profiles,ignore_index=True).to_csv(DATA/'individual_tracks_plotted_values.csv',index=False)
record={'panels':audits,'x_axis':'Uncompressed full HG002 T2T-Y chromosome',
        'y_axis':'Linear corrected relative copy number; common limits [-0.025, 1.37]',
        'above_range':'Every value above 1.3 represented by a triangle at 1.3; counts/maxima on plot; exact values in CSV',
        'selection':'Same paired eligible bins as the three main normalised atlases',
        'n_WGS_samples_per_population':1,'source_sha256':hashes,
        'source_files_unchanged':all(hashlib.sha256((OUT/p).read_bytes()).hexdigest()==v for p,v in hashes.items())}
(DATA/'individual_tracks_validation.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
assert record['source_files_unchanged']
print(json.dumps(audits,indent=2))
