"""Compact b-i panels from frozen flow, ratio, gel and clone ddPCR sources."""
from pathlib import Path
import json
import os
import sys

THESIS=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(THESIS.parent/'.figure_python_deps'))
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/competition_workflow_qc'
FLOW=OUT/'source_data/competition_workflow'
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from thesis_style import apply_style
from draw_competition_validation_panels import draw_validation

GREEN,CHERRY,GREY='#32965B','#D34F4F','#777E86'
WIDTH,HEIGHT=183,122


def render_lower():
    apply_style(7.5)
    fig=plt.figure(figsize=(WIDTH/25.4,HEIGHT/25.4))
    manifest=json.loads((FLOW/'flow_examples_manifest.json').read_text())
    controls=json.loads((FLOW/'flow_controls_manifest.json').read_text())
    by_label={r['label']:r for r in controls['controls']}
    rows=[by_label[k] for k in ('Unstained','GFP only','mCherry only')]+manifest['wells']
    weights=np.asarray(manifest['classifier']['weights'])
    intercept=manifest['classifier']['intercept']
    gfp_limit,cherry_limit=manifest['negative_rule']['limits']
    limits=(-8.5,12.5)
    def label(x,y,s,**kw):
        return fig.text(x/WIDTH,y/HEIGHT,s,**kw)
    for i,row in enumerate(rows):
        left=12+34.3*i
        bottom=83
        title_y=120
        ax=fig.add_axes([left/WIDTH,bottom/HEIGHT,31.8/WIDTH,31.8/HEIGHT])
        historical=row.get('historical_unstained',False)
        mixture=row.get('mixture')
        title=f'Mixture {mixture}' if mixture else row['label']
        subtitle='WT' if mixture else 'CTRL'
        label(left-2,title_y,chr(ord('b')+i),fontsize=9,fontweight='semibold',va='center')
        label(left+16.5,title_y,title,fontsize=7.1,fontweight='semibold',ha='center',va='center')
        label(left+15.9,title_y-3.7,subtitle,fontsize=6.4,color=GREY,ha='center',va='center')
        events=pd.read_csv(FLOW/row['events_csv'])
        values=events[['eGFP_arcsinh','mCherry_arcsinh']].to_numpy()
        assert len(values)==row['primary_singlets']
        assert ((values>=limits[0]) & (values<=limits[1])).all()
        ax.add_patch(Rectangle((limits[0],limits[0]),gfp_limit-limits[0],cherry_limit-limits[0],
                              facecolor='#F2F3F5',edgecolor='#A3A9B0',lw=.5,ls=(0,(3,2)),zorder=0))
        xx=np.linspace(*limits,500)
        ax.plot(xx,-(weights[0]*xx+intercept)/weights[1],color='#222222',lw=.55,ls=(0,(4,3)),zorder=3)
        if historical:
            ax.scatter(values[:,0],values[:,1],s=.75,c=GREY,alpha=.23,edgecolors='none',rasterized=True,zorder=2)
            label(left+15.9,76.8,'Historical reference',fontsize=6.3,color=GREY,ha='center',va='center')
            label(left+15.9,73.6,'Draft gate overlay',fontsize=6.3,color='#7954A1',ha='center',va='center')
        else:
            excluded=events.excluded.astype(bool).to_numpy()
            gfp=events.reporter_class.eq('GFP').to_numpy()
            for mask,colour in [(excluded,'#A3A9B0'),(~excluded & ~gfp,CHERRY),(~excluded & gfp,GREEN)]:
                ax.scatter(values[mask,0],values[mask,1],s=.85,c=colour,alpha=.32,edgecolors='none',rasterized=True,zorder=2)
            if mixture:
                gfp_genotype,cherry_genotype=('CTRL','LOY') if mixture==1 else ('LOY','CTRL')
                annotations=[(f'{gfp_genotype}–GFP {100*row["gfp_fraction"]:.1f}%',GREEN),
                             (f'{cherry_genotype}–mCherry {100*(1-row["gfp_fraction"]):.1f}%',CHERRY)]
            else:
                is_gfp=row['label']=='GFP only'
                percent=100*(row['gfp_fraction'] if is_gfp else 1-row['gfp_fraction'])
                annotations=[(f'{"GFP" if is_gfp else "mCherry"} {percent:.1f}%',GREEN if is_gfp else CHERRY),
                             (f'Excluded {100*row["negative_broad_fraction"]:.1f}%',GREY)]
            for y,(s,colour) in zip([76.8,73.6],annotations):
                label(left+15.9,y,s,fontsize=6.3,color=colour,ha='center',va='center')
        ax.set(xlim=limits,ylim=limits,xticks=[-5,0,5,10],yticks=[-5,0,5,10],
               xlabel='',ylabel='mCherry' if i==0 else '')
        ax.tick_params(labelsize=6.4,length=1.6,pad=1.3,labelleft=i==0)
        ax.xaxis.label.set_size(6.7);ax.xaxis.labelpad=1.2
        ax.yaxis.label.set_size(7);ax.yaxis.labelpad=2
    # Single recovery panel, with the reporter legend inside the empty upper left.
    ax=fig.add_axes([12/WIDTH,18/HEIGHT,43/WIDTH,35/HEIGHT])
    summary=pd.read_csv(OUT/'source_data/competition_known_ratios/summary.csv')
    assert len(summary)==10 and summary.technical_n.eq(2).all()
    ax.plot([0,100],[0,100],color='#A7ADB5',lw=.75,ls=(0,(3,3)),zorder=0)
    for orientation,colour,marker,style,name in [
            ('LOY_GFP',GREEN,'o','-','LOY–GFP'),('CTRL_GFP',CHERRY,'s','--','LOY–mCherry')]:
        d=summary.loc[summary.orientation.eq(orientation)].sort_values('nominal_percent')
        error=np.vstack([d.mean_percent-d.minimum_percent,d.maximum_percent-d.mean_percent])
        ax.errorbar(d.nominal_percent,d.mean_percent,yerr=error,color=colour,marker=marker,
                    ms=3,mfc='white',mew=.7,lw=.9,ls=style,elinewidth=.7,capsize=1.8,label=name)
    ax.set(xlim=(0,100),ylim=(0,100),xticks=[10,25,50,75,90],yticks=[0,25,50,75,100],
           xlabel='Nominal LOY (%)',ylabel='Measured LOY (%)')
    ax.tick_params(labelsize=6.4,length=1.6,pad=1.3)
    ax.xaxis.labelpad=2;ax.yaxis.labelpad=2
    ax.legend(loc='upper left',fontsize=6.5,handlelength=1.6,labelspacing=.2,
              handletextpad=.45,borderaxespad=.3)
    label(2,64,'g',fontsize=9,fontweight='semibold',va='center')
    label(31,64,'Known-ratio recovery',fontsize=7.7,fontweight='semibold',ha='center',va='center')
    label(96,69.5,'GFP · Flow axes: arcsinh(unmixed signal / 1000)',fontsize=6.5,color=GREY,ha='center',va='center')
    draw_validation(fig, WIDTH, HEIGHT)
    stem=DATA/'compact_b_to_i'
    for ext in ['pdf','svg']:
        fig.savefig(stem.with_suffix('.'+ext),dpi=600)
    plt.close(fig)
    return stem
