"""Restyle saved Fig01 panel data; run from any working directory."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from itertools import groupby
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "figures/chapter3"
DATA = FIG / "source_data/ht1080_rna_global_structure"
sys.path.insert(0, str(ROOT / "skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM = "ch03_main_ht1080_rna_global_structure"
COL = {"CTRL": "#3C5488", "LOY": "#E64B35", "PARENTAL": "#888888"}
SHAPE = {"AUG_2025": "o", "MAY_2026": "^"}
pc, qc, cr, design = [pd.read_csv(DATA / f"Fig01_design_QC_{x}.csv") for x in "abcd"]
assert len(pc)==12 and len(qc)==15 and len(cr)==144 and len(design)==12
plot_pc = pd.read_csv(DATA / "panel_a_adjusted_all15.csv")
plot_variance = pd.read_csv(DATA / "panel_a_variance.csv")
plot_variance = plot_variance[plot_variance.view == "adjusted"].set_index("component").variance_fraction
assert len(plot_pc)==15
order = [f"{prefix}{i}" for prefix in ["AC","AL","MC","ML"] for i in [1,2,3]]
qc_order = [f"{prefix}{i}" for prefix in ["AC","AL","AP","ML","MC"] for i in [1,2,3]]
qc = qc.set_index("short").loc[qc_order].reset_index()
matrix = cr.pivot(index="sample_2", columns="sample_1", values="r").loc[order,order]
assert np.allclose(matrix, matrix.T) and np.allclose(np.diag(matrix),1)
# Hierarchical clustering of existing correlations; no expression values recalculated.
distance = 1 - matrix.to_numpy()
distance = (distance + distance.T) / 2
np.fill_diagonal(distance, 0)
pd.DataFrame({"leaf_id":range(12),"sample":order}).to_csv(DATA/"correlation_linkage_leaf_ids.csv",index=False)
tree = linkage(squareform(distance), method="average", optimal_ordering=True)
order = [order[i] for i in leaves_list(tree)[::-1]]
matrix = matrix.loc[order, order]
pd.DataFrame({"heatmap_top_to_bottom": range(1,13), "sample": order}).to_csv(DATA/"correlation_cluster_order.csv",index=False)
pd.DataFrame(tree,columns=["child_1","child_2","distance","n_leaves"]).to_csv(DATA/"correlation_linkage.csv",index=False)

fig = plt.figure(figsize=(183/25.4,180/25.4))
a=fig.add_axes([.10,.61,.35,.26])
b=fig.add_axes([.60,.61,.37,.26])
c=fig.add_axes([.14,.145,.31,.30])
d=fig.add_axes([.64,.16,.33,.285])
def heading(ax, letter, title):
    pos=ax.get_position()
    title_y = .51 if letter in "cd" else pos.y1+.034
    fig.text(.035 if letter == "c" else pos.x0-.065,title_y,letter,fontsize=12,weight="semibold")
    fig.text(pos.x0,title_y,title,fontsize=9,weight="semibold")
heading(a,"a","")
heading(b,"b","")
heading(c,"c","")
heading(d,"d","Sampling design")
for row in plot_pc.itertuples():
    a.scatter(row.PC1,row.PC2,c=COL[row.treatment],marker=SHAPE[row.batch],s=30,edgecolors="white",linewidths=.4,zorder=3)
    dx,dy,align=5,4,"left"
    if row.treatment=="PARENTAL" and row.passage_number==2: dx,dy,align=-7,8,"right"
    if row.batch=="AUG_2025" and row.treatment=="CTRL" and row.passage_number in [2,3]: dx,dy,align=-6,-6,"right"
    if row.batch=="MAY_2026" and row.treatment=="LOY" and row.passage_number==3: dx,dy,align=-5,-6,"right"
    a.annotate(str(row.passage_number),(row.PC1,row.PC2),xytext=(dx,dy),textcoords="offset points",ha=align,fontsize=7,color=COL[row.treatment])
a.set(xlabel=f"PC1 ({100*plot_variance['PC1']:.1f}%)",ylabel=f"PC2 ({100*plot_variance['PC2']:.1f}%)",xlim=(-110,94),ylim=(-95,92))
a.set_xticks([-100,-50,0,50]); a.set_yticks([-75,-50,-25,0,25,50,75])
a.text(.02,.02,"15 libraries; 25,526 genes",transform=a.transAxes,fontsize=7,color="#555555")
# Keep the source assignment-panel order, including all parental libraries.
b.bar(np.arange(len(qc)),qc.assigned_fraction*100,color=[COL[t] for t in qc.treatment],width=.74)
b.set_xticks(np.arange(len(qc)),qc.short,rotation=90,fontsize=6.7)
b.set(ylabel="Assigned reads (%)",ylim=(0,70),xlim=(-.7,len(qc)-.3))
b.set_yticks([0,20,40,60]); b.set_xlabel("Period / population / passage",fontsize=7,labelpad=4)
# Use the original observed range, without clipping or changing correlation values.
im=c.imshow(matrix.to_numpy(),cmap="magma",vmin=float(matrix.min().min()),vmax=1,interpolation="nearest",aspect="equal")
c.set_xticks(range(12),order,rotation=90,fontsize=6.5)
c.set_yticks(range(12),order,fontsize=6.5)
c.tick_params(length=0,pad=3)
c.yaxis.tick_right()
for spine in c.spines.values(): spine.set_visible(False)
# A left-hand tree uses the same leaf order as the heatmap rows.
pos = c.get_position()
dendax = fig.add_axes([.035,pos.y0,pos.x0-.078,pos.height])
dendrogram(tree,orientation="left",no_labels=True,ax=dendax,
           color_threshold=0,above_threshold_color="#666666",
           link_color_func=lambda _: "#666666")
for collection in dendax.collections: collection.set_linewidth(.6)
dendax.set_ylim(0,120)
dendax.axis("off")
# Row annotations are keyed by sample ID and share the heatmap top-to-bottom order.
meta = pc.set_index("short").loc[order]
pos = c.get_position()
ann = fig.add_axes([pos.x0-.012,pos.y0,.008,pos.height])
for i, treatment in enumerate(meta.treatment):
    ann.add_patch(Rectangle((0,i),1,1,facecolor=COL[treatment],edgecolor="white",linewidth=.35))
ann.set(xlim=(0,1),ylim=(12,0)); ann.axis("off")
ann.text(.5,-.18,"Condition",ha="center",va="bottom",rotation=90,fontsize=6.3)
period = fig.add_axes([pos.x0-.036,pos.y0,.020,pos.height])
start = 0
for batch, samples in groupby(meta.batch):
    size = len(list(samples))
    label, shade = {"AUG_2025": ("Aug 2025", "#EEEEEE"), "MAY_2026": ("May 2026", "#DADADA")}[batch]
    period.add_patch(Rectangle((0,start),1,size,facecolor=shade,edgecolor="white",linewidth=.4))
    period.text(.5,start+size/2,label,ha="center",va="center",rotation=90,fontsize=6.5)
    start += size
period.set(xlim=(0,1),ylim=(12,0)); period.axis("off")
period.text(.5,-.18,"Period",ha="center",va="bottom",rotation=90,fontsize=6.3)
cbax=fig.add_axes([.14,.086,.27,.013]); cb=fig.colorbar(im,cax=cbax,orientation="horizontal")
cb.set_ticks([.94,.96,.98,1]); cb.ax.tick_params(labelsize=7,length=2); cb.outline.set_visible(False)
cb.set_label("Pearson r",fontsize=7,labelpad=3)
# Same established populations; spacing denotes order only.
for treatment,y in [("CTRL",1),("LOY",0)]:
    sub=design[design.treatment==treatment].sort_values("occasion")
    d.plot(sub.occasion,[y]*len(sub),c=COL[treatment],lw=1,zorder=1)
    for row in sub.itertuples():
        d.scatter(row.occasion,y,c=COL[treatment],marker=SHAPE[row.batch],s=32,zorder=3)
        d.text(row.occasion,y+.15,str(row.passage_number),ha="center",fontsize=7,color=COL[treatment])
d.set(xlim=(.5,7.5),ylim=(-.55,1.6))
d.set_xticks([2,6],["Aug 2025","May 2026"],fontsize=8)
d.set_yticks([1,0],["sgCTRL","Y-depleted"],fontsize=8)
# Secondary population labels; the meaning of 2x is defined only in the caption.
for y, label in [(1,"2\u00d7"),(0,"2\u00d7 CP + Q53")]:
    d.annotate(label,xy=(0,y),xycoords=d.get_yaxis_transform(),xytext=(-3,-12),textcoords="offset points",ha="right",va="center",fontsize=6.5,color="#888888")
d.tick_params(axis="both",length=0)
for spine in d.spines.values(): spine.set_visible(False)
d.set_xlabel("Successive passages within each period",fontsize=7,labelpad=8)
d.text(.5,1.02,"Same established populations",transform=d.transAxes,ha="center",fontsize=7,color="#555555")
d.text(.5,.03,"One population per condition\n12 libraries; six passage matches",transform=d.transAxes,ha="center",va="bottom",fontsize=7,color="#555555",linespacing=1.5)
handles=[Patch(facecolor=COL[t],label=l) for t,l in [("CTRL","sgCTRL"),("LOY","Y-depleted"),("PARENTAL","Parental")]]
fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.995),ncol=3,handlelength=1.1,columnspacing=2,fontsize=8)
handles=[Line2D([],[],color="#333333",marker=m,lw=0,markersize=5,label=l) for m,l in [("o","Aug 2025"),("^","May 2026")]]
a.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,1.17),ncol=2,fontsize=7,handletextpad=.4,columnspacing=1.5,borderpad=0)
(FIG/"previews").mkdir(exist_ok=True)
for ext in ["pdf","svg","png"]:
    target=FIG/f"{STEM}.{ext}"; stage=FIG/f"{STEM}.staging.{ext}"
    fig.savefig(stage,dpi=600,facecolor="white")
    stage.replace(target)
fig.savefig(FIG/"previews"/f"{STEM}_preview.png",dpi=180)
plt.close(fig)
checks={"panel_a_libraries":len(plot_pc),"panel_a_adjustment":"additive sampling period; population preserved; visualization only","panel_b_libraries":len(qc),"panel_c_correlations":len(cr),"panel_d_libraries":len(design),"correlation_range":[float(matrix.min().min()),float(matrix.max().max())],"assigned_percent_range":[float(qc.assigned_fraction.min()*100),float(qc.assigned_fraction.max()*100)],"font":"Open Sans","size_mm":[183,180],"png_dpi":600,"analysis_recomputed":False,"clustering":{"method":"average linkage","distance":"1 - Pearson r","optimal_leaf_ordering":True,"sample_order":order}}
(DATA/"render_checks.json").write_text(json.dumps(checks,indent=2))
print(json.dumps(checks))
