"""Plot the completed batch-aware DE extension; do not refit any model."""
from pathlib import Path
import sys,hashlib,json,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parent
THESIS=ROOT.parents[3]/"SEPT_11_2026_THESIS"
FIG=THESIS/"figures/chapter3"
DATA=ROOT/"outputs/tables"
sys.path.insert(0,str(THESIS/"skills/thesis-figure-style/scripts"))
from thesis_style import apply_style
apply_style(8)
STEM="ch03_supp_ht1080_batch_aware_de"
p=pd.read_csv(DATA/"batchDE_passage_period__combined.csv")
q=pd.read_csv(DATA/"batchDE_batch_condition__combined.csv").set_index("gene_id").loc[p.gene_id].reset_index()
aug=pd.read_csv(DATA/"batchDE_passage_period__AUG_2025.csv").set_index("gene_id").loc[p.gene_id]
may=pd.read_csv(DATA/"batchDE_passage_period__MAY_2026.csv").set_index("gene_id").loc[p.gene_id]
assert len(p)==25526 and (q.gene_id.to_numpy()==p.gene_id.to_numpy()).all()
sig=p.FDR<.05
colour=np.where(sig,np.where(p.logFC>0,"#E64B35","#3C5488"),"#BBBBBB")
fig,axes=plt.subplots(3,2,figsize=(183/25.4,229/25.4))
fig.subplots_adjust(left=.11,right=.975,top=.90,bottom=.065,hspace=.65,wspace=.42)
labels=["Combined population difference","Agreement with descriptive effects","Sensitivity to model choice","Period-specific effects","Voom mean-variance fit","Raw P-value distributions"]
for ax,letter,title in zip(axes.flat,"abcdef",labels):
 pos=ax.get_position();fig.text(pos.x0-.065,pos.y1+.018,letter,fontsize=12,weight="semibold");fig.text(pos.x0,pos.y1+.018,title,fontsize=8,weight="semibold")
ax=axes[0,0]
ax.scatter(p.logFC,-np.log10(p.FDR.clip(lower=np.finfo(float).tiny)),s=3,c=colour,alpha=.5,rasterized=True,linewidths=0)
ax.axhline(-np.log10(.05),c="#777777",ls="--",lw=.6)
ax.set(xlabel="Model log2 fold change",ylabel="-log10(BH-adjusted P)")
ax.text(.98,.96,f"FDR < 0.05\nHigher: {((p.logFC>0)&sig).sum():,}\nLower: {((p.logFC<0)&sig).sum():,}",transform=ax.transAxes,ha="right",va="top",fontsize=7)
def concord(ax,x,y,xlab,ylab):
 ax.scatter(x,y,s=2,c=colour,alpha=.3,linewidths=0,rasterized=True)
 lo=min(np.min(x),np.min(y));hi=max(np.max(x),np.max(y));ax.plot([lo,hi],[lo,hi],ls="--",lw=.6,c="#777777")
 ax.set(xlabel=xlab,ylabel=ylab)
 ax.text(.03,.96,f"Pearson r = {np.corrcoef(x,y)[0,1]:.3f}",transform=ax.transAxes,va="top",fontsize=7)
concord(axes[0,1],p.mean_difference,p.logFC,"Descriptive mean log2 difference","Passage-aware model log2 FC")
concord(axes[1,0],q.logFC,p.logFC,"Batch + condition model log2 FC","Passage-aware model log2 FC")
concord(axes[1,1],aug.logFC.to_numpy(),may.logFC.to_numpy(),"August model log2 FC","May model log2 FC")
ax=axes[2,0]
v=pd.read_csv(DATA/"batchDE_voom_points_passage_period.csv")
ax.scatter(v.mean_log2_count,v.sqrt_residual_SD,s=2,color="#999999",alpha=.2,linewidths=0,rasterized=True)
for name,col,label in [("passage_period","#3C5488","Passage-aware"),("batch_condition","#E64B35","Batch + condition")]:
 trend=pd.read_csv(DATA/f"batchDE_voom_trend_{name}.csv")
 ax.plot(trend.mean_log2_count,trend.sqrt_residual_SD,c=col,lw=1,label=label)
ax.set(xlabel="Mean log2 count",ylabel="Square root of residual SD")
ax.legend(fontsize=6.5,loc="upper right")
ax=axes[2,1]
for frame,col,label in [(p,"#3C5488","Passage-aware"),(q,"#E64B35","Batch + condition")]:
 ax.hist(frame.PValue,bins=np.linspace(0,1,31),histtype="step",lw=1,color=col,label=label)
ax.set(xlabel="Raw P value",ylabel="Genes (log scale)",yscale="log");ax.legend(fontsize=6.5,loc="upper right")
fig.legend(handles=[Line2D([],[],marker="o",ls="",color=c,label=l,markersize=4) for c,l in [("#3C5488","Lower in Y-depleted"),("#E64B35","Higher in Y-depleted"),("#BBBBBB","FDR >= 0.05")]],loc="upper center",bbox_to_anchor=(.54,.99),ncol=3,fontsize=7)
records=[]
for ext in ["pdf","svg","png"]:
 target=FIG/f"{STEM}.{ext}";stage=FIG/f"{STEM}.staging.{ext}"
 fig.savefig(stage,dpi=600);stage.replace(target)
 report=ROOT/"outputs/figures"/target.name;shutil.copy2(target,report)
 records.append({"file":target.name,"sha256":hashlib.sha256(target.read_bytes()).hexdigest(),"width_mm":183,"height_mm":229,"dpi":600,"font":"Open Sans"})
fig.savefig(FIG/"previews"/f"{STEM}_preview.png",dpi=150)
plt.close(fig)
pd.DataFrame(records).to_csv(ROOT/"outputs/figures/batchDE_figure_manifest.csv",index=False)
print(STEM)
