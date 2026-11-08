"""Render the verified 15-library contextual PCA in the thesis style."""
from pathlib import Path
import sys,json
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/"figures/chapter3"
DATA=FIG/"source_data/ht1080_rna_parental_pca"
sys.path.insert(0,str(ROOT/"skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM="ch03_supp_ht1080_parental_pca"
d=pd.read_csv(DATA/"PCA_parental_all15.csv")
v=pd.read_csv(DATA/"PCA_variance.csv").set_index("component").variance_fraction
assert len(d)==15 and d.treatment.value_counts().to_dict()=={"CTRL":6,"LOY":6,"PARENTAL":3}
COL={"CTRL":"#3C5488","LOY":"#E64B35","PARENTAL":"#888888"}
SHAPE={"AUG_2025":"o","MAY_2026":"^"}
fig,ax=plt.subplots(figsize=(125/25.4,113/25.4))
fig.subplots_adjust(left=.17,right=.96,bottom=.16,top=.77)
for row in d.itertuples():
    ax.scatter(row.PC1,row.PC2,s=35,marker=SHAPE[row.batch],c=COL[row.treatment],edgecolors="white",linewidths=.5,zorder=3)
    dx,dy,align=5,4,"left"
    if row.batch=="AUG_2025" and row.treatment=="CTRL" and row.passage_number==2: dx,dy,align=-6,3,"right"
    if row.batch=="AUG_2025" and row.treatment=="LOY" and row.passage_number==1: dx,dy=5,-9
    ax.annotate(str(row.passage_number),(row.PC1,row.PC2),xytext=(dx,dy),textcoords="offset points",ha=align,fontsize=7,color=COL[row.treatment])
ax.set(xlabel=f"PC1 ({100*v['PC1']:.1f}%)",ylabel=f"PC2 ({100*v['PC2']:.1f}%)",xlim=(-95,128),ylim=(-105,84))
ax.set_xticks([-75,-25,25,75,125]);ax.set_yticks([-100,-50,0,50])
fig.text(.565,.035,"15 libraries; 25,526 autosomal / non-PAR X genes",ha="center",fontsize=7,color="#555555")
handles=[Patch(facecolor=COL[t],label=l) for t,l in [("CTRL","sgCTRL"),("LOY","Y-depleted"),("PARENTAL","Parental")]]
fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.56,.99),ncol=3,fontsize=8,handlelength=1,columnspacing=1.2)
handles=[Line2D([],[],color="#333333",marker=m,linestyle="none",markersize=5,label=l) for m,l in [("o","Aug 2025"),("^","May 2026")]]
fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.56,.905),ncol=2,fontsize=7,columnspacing=1.8)
for ext in ["pdf","svg","png"]:
    stage=FIG/f"{STEM}.staging.{ext}";fig.savefig(stage,dpi=600);stage.replace(FIG/f"{STEM}.{ext}")
fig.savefig(FIG/"previews"/f"{STEM}_preview.png",dpi=220)
plt.close(fig)
(DATA/"render_checks.json").write_text(json.dumps({"libraries":15,"parental":3,"sgCTRL":6,"Y_depleted":6,"PC1_percent":100*v['PC1'],"PC2_percent":100*v['PC2'],"size_mm":[125,113],"font":"Open Sans","dpi":600},indent=2))
print(STEM)
