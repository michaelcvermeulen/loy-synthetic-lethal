"""All-15 before/after PCA and independently fitted within-period PCAs."""
from pathlib import Path
import sys,json
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[3];FIG=ROOT/"figures/chapter3";DATA=FIG/"source_data/ht1080_rna_parental_batch_pca"
sys.path.insert(0,str(ROOT/"skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM="ch03_supp_ht1080_parental_batch_adjusted_pca"
COL={"CTRL":"#3C5488","LOY":"#E64B35","PARENTAL":"#888888"};SHAPE={"AUG_2025":"o","MAY_2026":"^"}
variance=pd.concat([pd.read_csv(DATA/"PCA_variance.csv"),pd.read_csv(DATA/"PCA_within_batch_variance.csv")],ignore_index=True)
fig,axes=plt.subplots(2,2,figsize=(183/25.4,184/25.4))
fig.subplots_adjust(left=.095,right=.975,bottom=.09,top=.83,wspace=.38,hspace=.60)
for ax,view,letter,title in zip(axes.flat,["unadjusted","adjusted","batch1","batch2"],"abcd",["All libraries: TMM-normalised","All libraries: period-adjusted","Batch 1 (August 2025)","Batch 2 (May 2026)"]):
    data=pd.read_csv(DATA/f"PCA_{view}.csv"); assert len(data)=={"unadjusted":15,"adjusted":15,"batch1":9,"batch2":6}[view]
    v=variance[variance.view==view].set_index("component").variance_fraction
    for row in data.itertuples():
        ax.scatter(row.PC1,row.PC2,c=COL[row.treatment],marker=SHAPE[row.batch],s=28,edgecolors="white",linewidths=.4,zorder=3)
        dx,dy,ha=4,4,"left"
        if view=="unadjusted" and row.batch=="MAY_2026" and row.treatment=="CTRL" and row.passage_number==2: dx,dy,ha=-5,5,"right"
        if view=="unadjusted" and row.batch=="AUG_2025" and row.treatment=="CTRL" and row.passage_number==2: dx,dy,ha=-6,3,"right"
        if view=="unadjusted" and row.batch=="AUG_2025" and row.treatment=="LOY" and row.passage_number==1: dx,dy=4,-9
        if view=="adjusted" and row.batch=="AUG_2025" and row.treatment=="CTRL" and row.passage_number in [2,3]: dx,dy,ha=-6,-6,"right"
        if view=="adjusted" and row.batch=="MAY_2026" and row.treatment=="LOY" and row.passage_number==3: dx,dy,ha=-5,-6,"right"
        if view=="adjusted" and row.treatment=="PARENTAL" and row.passage_number==2: dx,dy,ha=-7,8,"right"
        if view=="batch1" and row.treatment=="PARENTAL" and row.passage_number==2: dx,dy=5,-9
        if view=="batch2" and row.treatment=="LOY" and row.passage_number==2: dx,dy,ha=-5,-6,"right"
        ax.annotate(str(row.passage_number),(row.PC1,row.PC2),xytext=(dx,dy),textcoords="offset points",ha=ha,color=COL[row.treatment],fontsize=6.5)
    ax.set_xlabel(f"PC1 ({100*v['PC1']:.1f}%)");ax.set_ylabel(f"PC2 ({100*v['PC2']:.1f}%)")
    if view=="unadjusted": ax.set(xlim=(-95,132),ylim=(-108,90));ax.set_xticks([-75,-25,25,75,125])
    elif view=="adjusted": ax.set(xlim=(-110,92),ylim=(-88,90));ax.set_xticks([-100,-50,0,50])
    else:
        ax.margins(x=.18,y=.20)
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4));ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    if view in ["unadjusted","adjusted"]: ax.set_yticks([-100,-50,0,50] if view=="unadjusted" else [-50,0,50])
    pos=ax.get_position();fig.text(pos.x0-.065,pos.y1+.027,letter,weight="semibold",fontsize=12);fig.text(pos.x0,pos.y1+.027,title,fontsize=8.5,weight="semibold")
    ax.text(.02,.02,f"{len(data)} libraries",transform=ax.transAxes,fontsize=7,color="#555555")
handles=[Patch(facecolor=COL[t],label=l) for t,l in [("CTRL","sgCTRL"),("LOY","Y-depleted"),("PARENTAL","Parental")]]
fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.54,.992),ncol=3,fontsize=8,handlelength=1,columnspacing=1.8)
handles=[Line2D([],[],color="#333333",marker=m,lw=0,markersize=5,label=l) for m,l in [("o","Aug 2025"),("^","May 2026")]]
fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.54,.948),ncol=2,fontsize=7,columnspacing=1.8)
for ext in ["pdf","svg","png"]:
    stage=FIG/f"{STEM}.staging.{ext}";fig.savefig(stage,dpi=600);stage.replace(FIG/f"{STEM}.{ext}")
fig.savefig(FIG/"previews"/f"{STEM}_preview.png",dpi=200)
plt.close(fig)
print(STEM)
