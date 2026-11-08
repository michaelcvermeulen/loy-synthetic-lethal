"""Combine frozen parental comparisons; no models are refitted."""
from pathlib import Path
import hashlib,json,sys,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'figures/chapter3'
DATA=OUT/'source_data/a549_parental_comparisons';DATA.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
stem='ch03_supp_a549_parental_comparisons'
inputs={};tables={};thresholds=[]
for group in ['sgCTRL','LOY']:
    folder=OUT/f'drafts/A549_parental_vs_{group}_2026-09-25'
    for kind,file in [('DE','DE_condition_only.csv'),('GSEA','GSEA_Hallmark.csv'),('summary','summary.csv')]:
        source=folder/file;shutil.copy2(source,DATA/f'{group}_{file}')
        inputs[str(source.relative_to(ROOT))]=hashlib.sha256(source.read_bytes()).hexdigest()
        tables[group,kind]=pd.read_csv(source)
fig=plt.figure(figsize=(7.2,8.5))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.65],left=.11,right=.90,bottom=.08,top=.96,hspace=.38,wspace=.4)
for j,(group,label,color) in enumerate([('sgCTRL','sgCTRL','#3C5488'),('LOY','Y-depleted','#E64B35')]):
    ax=fig.add_subplot(gs[0,j]);d=tables[group,'DE'];sig=d['adj.P.Val']<.05
    ax.scatter(d.logFC,-np.log10(d['P.Value']),s=3,c='#cccccc',alpha=.7,lw=0,rasterized=True)
    for sign,c in [(-1,color),(1,'#666666')]:
        k=sig & (np.sign(d.logFC)==sign)
        ax.scatter(d.loc[k,'logFC'],-np.log10(d.loc[k,'P.Value']),s=7,c=c,lw=0,rasterized=True)
    tested=int(d['P.Value'].notna().sum())
    cutoff=.05*int(sig.sum())/tested
    assert np.array_equal((d['P.Value']<cutoff).to_numpy(),sig.to_numpy())
    height=-np.log10(cutoff)
    thresholds.append(dict(comparison=group,tested=tested,significant=int(sig.sum()),nominal_P_cutoff=cutoff,minus_log10_P=height))
    ax.axvline(0,c='.65',lw=.6);ax.axhline(height,c='.6',ls='--',lw=.6)
    ax.text(.98,height+.10,'FDR = 0.05',transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=6,color='.4')
    ax.set(xlabel=f'log2FC (parental / {label})',ylabel='-log10(nominal P)',xlim=(-12,12),ylim=(-.1,9))
    ax.set_title(f'Parental versus {label}',loc='left',fontsize=9)
    ax.text(-.2,1.08,'ab'[j],transform=ax.transAxes,fontsize=11,fontweight='semibold')
    ax.text(.04,.97,f'{sig.sum()} genes at FDR <0.05',transform=ax.transAxes,va='top',fontsize=7)
    from adjustText import adjust_text
    texts=[]
    for _,r in d[d.gene_name.isin(['IFIT2','OASL','ISG15'])].iterrows():texts.append(ax.text(r.logFC,-np.log10(r['P.Value']),r.gene_name,fontsize=6.5))
    adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',lw=.5,color='.5'))
g1=tables['sgCTRL','GSEA'].set_index('pathway');g2=tables['LOY','GSEA'].set_index('pathway')
keep=sorted(set(g1.index[g1.BH_across_collections<.05])|set(g2.index[g2.BH_across_collections<.05]))
nes=pd.DataFrame({'Parental / sgCTRL':g1.loc[keep,'NES'],'Parental / Y-depleted':g2.loc[keep,'NES']})
nes=nes.loc[nes.mean(axis=1).sort_values(ascending=False).index]
fdr=pd.DataFrame({'Parental / sgCTRL':g1.loc[nes.index,'BH_across_collections'],'Parental / Y-depleted':g2.loc[nes.index,'BH_across_collections']})
nes.to_csv(DATA/'Hallmark_NES_displayed.csv');fdr.to_csv(DATA/'Hallmark_global_FDR_displayed.csv')
# Reserve left-hand space for full pathway labels within the figure.
sub=gs[1,:].subgridspec(1,3,width_ratios=[1.7,2,.14],wspace=.15)
ax=fig.add_subplot(sub[0,1]);cax=fig.add_subplot(sub[0,2])
cmap=LinearSegmentedColormap.from_list('effect',['#3C5488','#FFFFFF','#E64B35'])
im=ax.imshow(nes,cmap=cmap,norm=TwoSlopeNorm(vmin=-3.2,vcenter=0,vmax=3.2),aspect='auto')
labels=[x.replace('HALLMARK_','').replace('_',' ') for x in nes.index]
ax.set_yticks(range(len(nes)),labels,fontsize=7)
ax.set_xticks([0,1],['Parental /\nsgCTRL','Parental /\nY-depleted'],fontsize=8)
for i in range(len(nes)):
    for j in range(2):ax.text(j,i,f'{nes.iloc[i,j]:.2f}'+('*' if fdr.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=7,color='white' if abs(nes.iloc[i,j])>2.1 else '#222222')
ax.set_title('Hallmark enrichment',loc='left',fontsize=9,pad=10)
ax.text(-1.05,1.035,'c',transform=ax.transAxes,fontsize=11,fontweight='semibold')
fig.colorbar(im,cax=cax,label='NES (positive: parental)')
fig.text(.11,.015,'* Global BH FDR <0.05. Heatmap colour represents enrichment direction, not sample identity.',fontsize=7)
pd.DataFrame(thresholds).to_csv(DATA/'volcano_FDR_thresholds.csv',index=False)
print(pd.DataFrame(thresholds).to_string(index=False))
exports={}
for ext in ['pdf','svg','png']:
    file=OUT/f'{stem}.{ext}';fig.savefig(file,dpi=300);exports[file.name]=hashlib.sha256(file.read_bytes()).hexdigest()
plt.close(fig)
(DATA/'provenance.json').write_text(json.dumps({'inputs_sha256':inputs,'exports_sha256':exports,'note':'Frozen exploratory comparisons, unpaired; positive effects are parental-high. No models refitted.'},indent=2),encoding='utf-8')
print('Wrote',stem,'with',len(nes),'Hallmark sets')
