"""Shared relative-CN scales for five public baseline WGS profiles."""
from pathlib import Path
import json
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import numpy as np
import pandas as pd
from make_cell_line_profiles import THESIS, BLUE, apply_style

OUT=THESIS/'figures/chapter3'
SOURCE=OUT/'source_data/selected_models_y_cn'
STEM='ch03_main_selected_models_Y_copy_number'
LINES=['HT1080','A549','DU145','H2452','HepG2']

def label(v):return '<0.001' if v<.001 else f'{v:.3f}'

def main():
    apply_style(base_size=8)
    seg=pd.read_csv(SOURCE/'chrY_wgs_segments.csv')
    cn=pd.read_csv(SOURCE/'y_gene_wgs_cn.csv')
    bands=pd.read_csv(SOURCE/'chrY_cytobands_grch38.csv')
    coords=pd.read_csv(SOURCE/'y_gene_coordinates_grch38.csv').sort_values('start')
    genes=coords.gene.tolist()
    mat=cn.pivot(index='cell_line',columns='gene',values='relative_copy_ratio').loc[LINES,genes]
    assert np.isfinite(mat.to_numpy()).all() and mat.max().max()<1
    assert seg.End.max()==26637000 and set(seg.cell_line)==set(LINES)
    fig=plt.figure(figsize=(183/25.4,190/25.4))
    fig.text(.05,.966,'a',fontsize=12,weight='semibold')
    fig.text(.09,.966,'Chromosome-Y copy-number segments',fontsize=10,weight='semibold')
    fig.text(.09,.939,'DepMap Public 26Q1 | Baseline WGS profiles | GRCh38',fontsize=7.5,color='#555555')
    ax=fig.add_axes([.22,.87,.70,.034]);ax.set_xlim(0,27);ax.set_ylim(0,1);ax.set_axis_off()
    shades={'gneg':'white','gpos50':'#B9B9B9','gvar':'#E8E8E8'}
    for b in bands.itertuples():
        left=b.start/1e6;right=min(b.end/1e6,27)
        if left>=27:continue
        if b.stain=='acen':
            xy=[(left,.15),(left,.85),(right,.5)] if b.band.startswith('p') else [(left,.5),(right,.85),(right,.15)]
            ax.add_patch(Polygon(xy,facecolor='#777777',edgecolor='#666666',lw=.6))
        else:ax.add_patch(Rectangle((left,.15),right-left,.7,facecolor=shades[b.stain],edgecolor='#777777',lw=.6))
    # GRCh38 PAR1: 1-based inclusive 10,001-2,781,479; plotting intervals are 0-based.
    par_start,par_end=.010,2.781479
    ax.add_patch(Rectangle((par_start,.15),par_end-par_start,.7,
                           facecolor='#4DBBD5',edgecolor='#007A87',lw=.65,alpha=.60))
    for x,text in [(1.4,'PAR1'),(6.4,'Yp'),(10.45,'Centromere'),(19,'Yq')]:
        ax.text(x,1.12,text,fontsize=7.2,ha='center',va='bottom',clip_on=False)
    for i,line in enumerate(LINES):
        ax=fig.add_axes([.22,.79-i*.070,.70,.052])
        ax.set_xlim(0,27);ax.set_ylim(-.06,1.02)
        ax.axvspan(26.637,27,facecolor='#EEEEEE',hatch='////',edgecolor='#CCCCCC',lw=0)
        ax.axvspan(10.3,10.6,color='#EEEEEE',lw=0,zorder=0)
        ax.axvspan(par_start,par_end,facecolor='#4DBBD5',alpha=.10,lw=0,zorder=0)
        ax.axvline(par_end,color='#007A87',lw=.55,ls=':',alpha=.7,zorder=0)
        for val in [0,.5,1]:ax.axhline(val,color='#E5E5E5',lw=.55,zorder=0)
        sub=seg.loc[seg.cell_line.eq(line)]
        for s in sub.itertuples():
            ax.hlines(s.RelCopyRatio,s.Start/1e6,s.End/1e6,color=BLUE,lw=2)
            if s.End-s.Start>200000:
                ax.text((s.Start+s.End)/2e6,s.RelCopyRatio+.10,label(s.RelCopyRatio),
                        ha='center',va='bottom',fontsize=7,color=BLUE)
        # Show segment transitions without creating interpolated measurements.
        for a,b in zip(sub.itertuples(),list(sub.itertuples())[1:]):
            if a.End==b.Start:ax.plot([a.End/1e6]*2,[a.RelCopyRatio,b.RelCopyRatio],color=BLUE,lw=.8)
        ax.text(-.115,.50,line,transform=ax.transAxes,ha='right',va='center',fontsize=9,weight='semibold')
        ax.set_yticks([0,.5,1],['0','0.5','1'],fontsize=6.5)
        ax.tick_params(axis='y',length=2,pad=3)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([0,5,10,15,20,25])
        if i<4:ax.tick_params(axis='x',bottom=False,labelbottom=False);ax.spines['bottom'].set_visible(False)
        else:ax.set_xlabel('chrY position (Mb; displayed interval 0\u201327 Mb)',fontsize=8)
        if line=='H2452':
            ax.annotate('2.93 Mb',(2.928,.607),xytext=(4.5,.88),fontsize=6.8,
                        arrowprops={'arrowstyle':'-','lw':.6,'color':'#777777'})
    fig.text(.037,.674,'Relative copy-number ratio',rotation=90,fontsize=8,ha='center',va='center')
    fig.text(.22,.441,'Segment calls end at 26.637 Mb; distal Yq is not reported.',fontsize=7.2,color='#555555')
    fig.text(.05,.391,'b',fontsize=12,weight='semibold')
    fig.text(.09,.391,'Gene-level copy number in chromosome order',fontsize=10,weight='semibold')
    ax=fig.add_axes([.22,.145,.70,.190])
    im=ax.imshow(mat,cmap='viridis',vmin=0,vmax=1,aspect='auto',interpolation='nearest')
    ax.set_yticks(range(5),LINES,fontsize=8.5)
    ax.set_xticks(range(9),genes,fontsize=7.2,rotation=45,ha='left')
    ax.xaxis.tick_top();ax.tick_params(axis='both',length=0,pad=5)
    for spine in ax.spines.values():spine.set_visible(False)
    for i in range(5):
        for j in range(9):
            value=mat.iloc[i,j]
            text='<0.001' if value<.001 else f'{value:.2f}'
            ax.text(j,i,text,ha='center',va='center',fontsize=7,
                    color='white' if value<.56 else '#222222')
    ax.set_xticks(np.arange(-.5,9,1),minor=True);ax.set_yticks(np.arange(-.5,5,1),minor=True)
    ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',length=0)
    cbax=fig.add_axes([.35,.099,.44,.013]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal')
    cb.set_label('Relative copy-number ratio',fontsize=7.3,labelpad=3)
    cb.ax.tick_params(labelsize=6.5,length=2);cb.outline.set_visible(False)
    fig.text(.09,.035,'A549: public Y-DNA and RNA profiles disagree. Values describe public models, not the study stocks.',fontsize=7,color='#555555')
    for ext,dpi,folder in [('pdf',600,OUT),('svg',600,OUT),('png',600,OUT),('png',180,OUT/'previews')]:
        target=folder/f'{STEM}.{ext}';stage=folder/f'{STEM}.staging.{ext}'
        fig.savefig(stage,dpi=dpi)
        for attempt in range(10):
            try:stage.replace(target);break
            except OSError:
                if attempt==9:raise
                time.sleep(.5)
    plt.close(fig)
    print(STEM)

if __name__=='__main__':main()
