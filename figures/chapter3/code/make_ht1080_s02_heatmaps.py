"""Compact S02 pathway heatmaps using the previously selected saved results."""
from pathlib import Path
import sys,json,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
FIG=Path(__file__).resolve().parents[1];ROOT=FIG.parents[1]
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
SRC=FIG/'source_data/ht1080_rna_programme_overview'
a=pd.read_csv(SRC/'panel_c_pathways.csv')
b=pd.read_csv(SRC/'parental_GSEA_display.csv')
def group(p):return 'H' if p.startswith('HALLMARK_') else 'GO' if p.startswith('GOBP_') else 'R'
for table in [a,b]:
    table['group']=table.pathway.map(group)
    table['group_order']=table['group'].map({'H':0,'GO':1,'R':2})
a=a.sort_values(['group_order','NES','FDR']).reset_index(drop=True)
b=b.sort_values(['group_order','NES','FDR']).reset_index(drop=True)
vals=a[['NES','AUG_2025_NES','MAY_2026_NES']].to_numpy()
vmax=max(abs(vals).max(),b.NES.abs().max())
cmap=LinearSegmentedColormap.from_list('NES',['#3C5488','#F7F7F7','#E64B35'])
norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
fig=plt.figure(figsize=(10.5,6.3))
axes=[fig.add_axes([.31,.16,.18,.75]),fig.add_axes([.84,.16,.055,.75])]
def label(path):
    x=path
    for pre in ['HALLMARK_','GOBP_','REACTOME_']:x=x.removeprefix(pre)
    return '\n'.join(textwrap.wrap(f'[{group(path)}] '+x.replace('_',' ').lower(),38))
for ax,table,matrix in [(axes[0],a,vals),(axes[1],b,b[['NES']].to_numpy())]:
    im=ax.imshow(matrix,cmap=cmap,norm=norm,aspect='auto')
    ax.set_yticks(range(len(table)),[label(p) for p in table.pathway],fontsize=7)
    ax.tick_params(length=0,pad=4)
    ax.spines[:].set_visible(False)
    for i in range(1,len(table)):
        if table.group.iloc[i]!=table.group.iloc[i-1]:ax.axhline(i-.5,color='white',lw=2)
    for i in range(len(table)):
        for j in range(matrix.shape[1]):
            q=matrix[i,j];suffix=''
            if ax is axes[0] and table[['FDR','AUG_2025_FDR','MAY_2026_FDR']].iloc[i,j]<.05:suffix='*'
            ax.text(j,i,f'{q:.2f}{suffix}',ha='center',va='center',fontsize=6.5,color='white' if abs(q)>vmax*.58 else '#222222')
axes[0].set_xticks([0,1,2],['Overall','August','May'],fontsize=7)
axes[1].set_xticks([0],['August'],fontsize=7)
for i,q in enumerate(b.FDR):axes[1].text(1.24,i,f'{q:.2g}',transform=axes[1].get_yaxis_transform(),va='center',fontsize=7)
axes[1].text(1.24,1.015,'FDR',transform=axes[1].transAxes,fontsize=7)
fig.text(.018,.95,'a',fontsize=12,weight='semibold');fig.text(.055,.95,'Y-depleted versus sgCTRL',fontsize=9,weight='semibold')
fig.text(.53,.95,'b',fontsize=12,weight='semibold');fig.text(.56,.95,'sgCTRL versus parental',fontsize=9,weight='semibold')
cb=fig.colorbar(im,cax=fig.add_axes([.36,.052,.27,.02]),orientation='horizontal');cb.set_label('Normalised enrichment score (NES)',fontsize=7,labelpad=2);cb.ax.tick_params(labelsize=7,length=2);cb.outline.set_visible(False)
fig.text(.02,.095,'H: Hallmark   GO: GO Biological Process   R: Reactome',fontsize=7)
fig.text(.66,.065,'* FDR < 0.05 in a\nAll displayed sets in b: FDR < 0.05',fontsize=7)
stem='ch03_supp_ht1080_rna_programme_overview'
for ext in ['pdf','svg','png']:fig.savefig(FIG/f'{stem}.{ext}',dpi=300)
fig.savefig(FIG/'previews'/f'{stem}_preview.png',dpi=160)
a.to_csv(SRC/'heatmap_a_pathways.csv',index=False);b.to_csv(SRC/'heatmap_b_parental.csv',index=False)
(SRC/'heatmap_layout.json').write_text(json.dumps({'panels':{'a':'former c; same 12 pathways','b':'former e-g; same 17 pathways'},'size_inches':[10.5,6.3],'shared_NES_limits':[-float(vmax),float(vmax)],'analysis_recomputed':False},indent=2))
print('Saved compact two-panel S02.')
