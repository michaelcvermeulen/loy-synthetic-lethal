"""Standalone May-only RNA version of F09c; original F09 remains unchanged."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from adjustText import adjust_text

R=Path(__file__).resolve().parents[3];F=R/'figures/chapter3'
D=F/'source_data/ht1080_F09c_May';D.mkdir(exist_ok=True)
sys.path.insert(0,str(R/'skills/thesis-figure-style/scripts'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from thesis_style import apply_style
from rna_atac_scatter import draw_link_points,link_legend_handles
apply_style(8)
p=pd.read_csv(F/'source_data/ht1080_rna_atac_peak_scatter/promoter_peak_links.csv')
g=pd.read_csv(F/'source_data/ht1080_rna_atac_integration/integration_genes.csv')
r=pd.read_csv(R.parent/'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1/tables/gene_consistency.csv')
r['ensembl']=r.gene_id.str.split('.').str[0]
mapping=g[['entrez_id','ensembl']].merge(r[['ensembl','MAY_logFC','MAY_FDR']],on='ensembl',validate='one_to_one')
q=p.merge(mapping.rename(columns={'entrez_id':'entrez'}),on='entrez',how='left',validate='many_to_one')
assert len(q)==len(p) and q[['MAY_logFC','MAY_FDR']].notna().all().all()
q['combined_rna_logFC']=q.rna_logFC
q['combined_rna_FDR']=q.rna_FDR
q['combined_category']=q.category
q['rna_logFC']=q.MAY_logFC;q['rna_FDR']=q.MAY_FDR
de=(q.rna_FDR<.05)&(q.rna_logFC.abs()>=1)
q['rna_class']=np.where(de,np.where(q.rna_logFC>0,'up','down'),'not_DE')
q['category']=np.select([~de&q.peak_class.eq('ns'),~de,q.peak_class.eq('ns'),q.rna_class.eq('up')==q.peak_class.eq('gained')],['neither','peak_only','RNA_only','concordant'],default='opposite')
q.to_csv(D/'May_promoter_peak_links.csv',index=False)
both=q[q.category.isin(['concordant','opposite'])]
stats={'links':len(q),'genes':int(q.entrez.nunique()),'all_spearman':float(spearmanr(q.rna_logFC,q.logFC).statistic),'combined_all_spearman':float(spearmanr(q.combined_rna_logFC,q.logFC).statistic),'joint_links':len(both),'joint_spearman':float(spearmanr(both.rna_logFC,both.logFC).statistic),'categories':q.category.value_counts().to_dict(),'joint_directions':both.groupby(['rna_class','peak_class']).size().to_dict()}
stats['joint_directions']={'/'.join(k):int(v) for k,v in stats['joint_directions'].items()}
(D/'statistics.json').write_text(json.dumps(stats,indent=2))
fig,ax=plt.subplots(figsize=(7.2,5.1));fig.subplots_adjust(left=.11,right=.76,bottom=.13,top=.97)
draw_link_points(ax,q)
xl=max(5.8,float(q.rna_logFC.abs().max())+.5)
ax.set(xlim=(-xl,xl),ylim=(-3.2,3.7),xlabel='May RNA log₂FC (Y-depleted / sgCTRL)',ylabel='Promoter peak ATAC log₂FC')
ax.text(.02,.98,f"May RNA: FDR < 0.05 and |log₂FC| ≥ 1; ATAC: FDR < 0.10\nSpearman ρ = {stats['all_spearman']:.2f} ({len(q):,} links)\nBoth assays: ρ = {stats['joint_spearman']:.2f} ({len(both):,} links)",transform=ax.transAxes,va='top',fontsize=7)
named=both[~both.gene_name.str.match(r'^(LOC|ENSG)',na=False)].copy()
named['mag']=(named.rna_logFC*named.logFC).abs()
labels=pd.concat([named.sort_values('mag',ascending=False).drop_duplicates('entrez').head(20), named[named.rna_class.eq('up')].sort_values('mag',ascending=False).drop_duplicates('entrez')]).drop_duplicates('entrez')
labels[['gene_name','rna_logFC','rna_FDR','logFC','FDR','category']].to_csv(D/'labelled_genes.csv',index=False)
texts=[]
for row in labels.itertuples():
    color='#222222' if row.category=='opposite' else ('#3C5488' if row.rna_logFC<0 else '#E64B35')
    texts.append(ax.text(row.rna_logFC,row.logFC,row.gene_name,fontsize=6,color=color,fontstyle='italic'))
adjust_text(texts,ax=ax,expand=(1.4,1.6),force_text=(.5,.8),arrowprops=dict(arrowstyle='-',color='.55',lw=.4))
ax.legend(handles=link_legend_handles(),loc='upper left',bbox_to_anchor=(1.01,1),fontsize=6.5,frameon=False)
stem='ch03_supp_ht1080_F09c_May_RNA_draft'
for ext in ['pdf','svg','png']:fig.savefig(F/(stem+'.'+ext),dpi=250)
plt.close(fig)
(D/'README.md').write_text('May-only exploratory F09c. Same promoter peak–gene links and ATAC effects as original. May effects and FDR joined via integration Entrez-to-Ensembl mapping and current gene_consistency.csv. RNA DE uses May FDR <0.05 and absolute log2FC >=1, without requiring August direction agreement. ATAC FDR <0.10. Colours/classes recomputed for May. No model refit. Correlations descriptive: links share genes/peaks and are not independent biological replicates. Top 20 named joint links by absolute product of effects plus all jointly significant RNA-up genes labelled, one per gene. Combined comparison statistic retained. F09 not replaced. Reproduce with code/plot_ht1080_F09c_May.py.\n',encoding='utf-8')
print(json.dumps(stats,indent=2))
