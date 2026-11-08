"""S18: joint RNA/promoter changes with chromosome 19 regional context. No model refitting."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
DATA=FIG/'source_data/ht1080_rna_atac_integration'
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
LOWER,HIGHER,GREY='#3C5488','#E64B35','#555555'
g=pd.read_csv(DATA/'integration_genes.csv')
# d  jointly significant genes (RNA thesis criteria, ATAC promoter peak FDR < 0.10)
supfig = plt.figure(figsize=(183 / 25.4, 205 / 25.4))
d = supfig.add_axes([.145, .065, .22, .82])
jg = g[g.joint_010.isin(["concordant", "opposite", "mixed_promoter"])].sort_values("rna_logFC").reset_index(drop=True)
row_col = np.where(jg.joint_010 != "concordant", "#222222", np.where(jg.rna_logFC > 0, HIGHER, LOWER))
for i, r in jg.iterrows():
    d.plot([r.rna_logFC, r.atac_median_logFC], [i, i], color="#BBBBBB", lw=.5, zorder=1)
d.scatter(jg.rna_logFC, jg.index, s=7, c=row_col, zorder=3)
d.scatter(jg.atac_median_logFC, jg.index, s=7, facecolors="white", edgecolors=row_col, linewidths=.6, zorder=3)
d.set_yticks(jg.index, jg.gene_name, fontsize=5.1, style="italic")
for t_, chrom, colour in zip(d.get_yticklabels(), jg.chrom, row_col):
    t_.set_color(colour)
    if chrom == "chr19":
        t_.set_weight("semibold")
d.set(ylim=(len(jg) - .4, -.6), xlim=(-4.2, 2.8), xlabel="log₂FC (Y-depleted / sgCTRL)")
d.set_xticks([-4, -2, 0, 2])
d.tick_params(axis="y", length=0, pad=2)
d.axvline(0, color="#999999", lw=.5)
d.legend(handles=[Line2D([], [], ls="", marker="o", markersize=3.5, color="#555555", label="RNA"),
                  Line2D([], [], ls="", marker="o", markersize=3.5, markerfacecolor="white", markeredgecolor="#555555",
                         label="Promoter ATAC (median)")],
         loc="lower center", bbox_to_anchor=(.4, 1.005), ncol=2, fontsize=5.1, handletextpad=.2, columnspacing=.8,
         borderpad=0)
n_conc = int((jg.joint_010 == "concordant").sum()); n_opp = int((jg.joint_010 == "opposite").sum())
d.set_title(f"Jointly significant genes (n = {len(jg)})", fontsize=6.4, weight="semibold", pad=30)
d.text(.4, 1.045, f"{n_conc} concordant, {n_opp} opposite (black) · bold = chr19 ({(jg.chrom == 'chr19').sum()})",
       transform=d.transAxes, ha="center", fontsize=5.0, color=GREY)



supfig.text(.015,.965,'a',fontsize=12,weight='semibold')
CND=FIG/'source_data/ht1080_atac_copy_number'
reg=pd.read_csv(CND/'region_CN_context.csv');reg=reg[reg.chrom=='chr19'].copy();reg['Mb']=(reg.start+reg.end)/2e6
bins=pd.read_csv(CND/'WGS_hg38_100kb_bins.csv');bins=bins[(bins.chrom=='chr19')&bins.use&bins.delta_segment_log2.notna()]
summary=reg.groupby(np.floor(reg.Mb).astype(int)).agg(n=('logFC','size'),median_log2FC=('logFC','median'),negative_fraction=('logFC',lambda v:(v<0).mean()))
summary.index.name='bin_start_Mb';summary.to_csv(DATA/'S18_chr19_1Mb_summary.csv')
axes=[]
for letter,box,limits,title in [('b',[.48,.57,.49,.30],(0,59.2),'Chromosome 19'),('c',[.48,.085,.49,.34],(50,59.2),'Distal chromosome 19')]:
 ax=supfig.add_axes(box);axes.append(ax)
 sub=reg[reg.Mb.between(*limits)]
 ax.scatter(sub.Mb,sub.logFC,s=3,color='#BEBEBE',lw=0,rasterized=True,label='Tested peaks')
 sig=sub[sub.FDR<.05];ax.scatter(sig.Mb,sig.logFC,s=9,color=LOWER,lw=0,zorder=3,label='ATAC FDR <0.05')
 ss=summary[(summary.index>=limits[0])&(summary.index<limits[1])]
 ax.plot(ss.index+.5,ss.median_log2FC,color='#222222',lw=1,marker='o',ms=2,label='1 Mb median')
 ax.scatter((bins.start0+bins.end0)/2e6,bins.delta_segment_log2,color='#8C6BB1',s=2,lw=0,label='Earlier D0 lpWGS')
 ax.axhline(0,color='#AAAAAA',lw=.5);ax.axvspan(52,53,color='#F5E9D6',alpha=.6,zorder=-1)
 ax.set(xlim=limits,ylim=(-3.2,3.0 if letter=='c' else 2.0))
 ax.set_xlabel('Position (Mb, hg38)',fontsize=6.2)
 ax.set_ylabel('ATAC / lpWGS difference (log2)',fontsize=6.2)
 ax.set_title(title,fontsize=7,weight='semibold',pad=8);ax.tick_params(labelsize=6)
 supfig.text(.415,box[1]+box[3]+.03,letter,fontsize=12,weight='semibold')
axes[0].legend(loc='upper right',ncol=2,fontsize=5.3,frameon=False,handlelength=1.1,columnspacing=.7)
axes[0].text(.03,.04,'52-53 Mb highlighted',transform=axes[0].transAxes,fontsize=5.8,color=GREY)
axes[1].set_xticks(range(50,60))
# Label unique linked peaks; shared promoter assignments are kept together.
zg=jg[(jg.chrom=='chr19')&jg.gene_name.str.match('ZNF|ZFP')].copy()
links=[]
for row in zg.itertuples():
 for peak in row.peak_ids.split(';'):
  if peak in set(reg.loc[reg.Mb.between(50,59.2),'peak_id']):links.append({'peak_id':peak,'gene':row.gene_name})
labs=pd.DataFrame(links).groupby('peak_id').gene.agg(lambda x:' / '.join(x)).reset_index().merge(reg[['peak_id','Mb','logFC']],on='peak_id').sort_values('Mb')
labs.to_csv(DATA/'S18_ZNF_peak_labels.csv',index=False)
for xpos,row in zip(np.linspace(50.7,57.4,len(labs)),labs.itertuples()):
 axes[1].annotate(row.gene,(row.Mb,row.logFC),xytext=(xpos,.8),rotation=65,ha='left',va='bottom',fontsize=5.0,color=LOWER,
                  arrowprops=dict(arrowstyle='-',color='#888888',lw=.45,shrinkA=1,shrinkB=2))
axes[1].text(.03,.035,'52-53 Mb: 22/26 peaks lower; median log2FC -0.53',transform=axes[1].transAxes,fontsize=5.7,color=GREY)
(FIG/'previews').mkdir(exist_ok=True)
for ext in ['pdf','svg','png']:
 target=FIG/f'ch03_supp_ht1080_joint_promoter_genes.{ext}';stage=FIG/f'ch03_supp_ht1080_joint_promoter_genes.staging.{ext}'
 supfig.savefig(stage,dpi=600,facecolor='white',bbox_inches='tight',pad_inches=.04);stage.replace(target)
supfig.savefig(FIG/'previews/ch03_supp_ht1080_joint_promoter_genes_preview.png',dpi=180,bbox_inches='tight',pad_inches=.04)
plt.close(supfig)
print('Written S18 with 63 genes, chromosome 19 overview and labelled zoom.')
