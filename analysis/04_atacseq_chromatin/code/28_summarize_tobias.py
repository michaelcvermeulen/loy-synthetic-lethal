from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'motifs/TF_extension/TOBIAS'
assert (OUT/'ANALYSIS_COMPLETE.json').exists()

def label_points(ax,points):
    """Place labels without overlapping labels or covering another labeled point."""
    ax.figure.canvas.draw()
    renderer=ax.figure.canvas.get_renderer();bounds=ax.get_window_extent(renderer)
    anchors=[ax.transData.transform((x,y)) for _,x,y in points]
    occupied=[]
    if ax.get_legend() is not None:occupied.append(ax.get_legend().get_window_extent(renderer))
    offsets=sorted([(dx,dy) for dx in [6,-6,18,-18,30,-30,45,-45,60,-60]
                    for dy in [6,-8,18,-20,30,-32,42,-44,54,-56]],key=lambda p:p[0]**2+p[1]**2)
    for (label,x,y),anchor in zip(points,anchors):
        placed=False
        for dx,dy in offsets:
            position=ax.transData.inverted().transform(anchor+np.array([dx,dy])*ax.figure.dpi/72)
            text=ax.text(*position,label,fontsize=8,ha='left' if dx>0 else 'right',va='center',
                         zorder=4,bbox=dict(facecolor='white',edgecolor='none',alpha=.75,pad=.2))
            box=text.get_window_extent(renderer).padded(2)
            inside=box.x0>bounds.x0+2 and box.x1<bounds.x1-2 and box.y0>bounds.y0+2 and box.y1<bounds.y1-2
            if inside and not any(box.overlaps(b) for b in occupied) and not any(box.contains(*p) for p in anchors):
                occupied.append(box)
                ax.annotate('',xy=(x,y),xytext=position,arrowprops=dict(arrowstyle='-',color='0.4',lw=.5,shrinkA=4,shrinkB=3),zorder=2)
                placed=True;break
            text.remove()
        assert placed,f'No clear label position for {label}'
df=pd.read_csv(OUT/'BINDetect/bindetect_results.txt',sep='\t')
assert 'LOY_CTRL_change' in df and 'LOY_CTRL_pvalue' in df
assert df.motif_id.is_unique
df['LOY_CTRL_pvalue']=pd.to_numeric(df.LOY_CTRL_pvalue)
df['LOY_CTRL_change']=pd.to_numeric(df.LOY_CTRL_change)
valid=np.isfinite(df.LOY_CTRL_pvalue)&np.isfinite(df.LOY_CTRL_change)
assert df.loc[valid,'LOY_CTRL_pvalue'].between(0,1).all()
p=df.loc[valid,'LOY_CTRL_pvalue'].to_numpy();order=np.argsort(p)
q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
adjusted=np.empty(len(p));adjusted[order]=q
df['pooled_site_BH_FDR']=np.nan;df.loc[valid,'pooled_site_BH_FDR']=adjusted
cv=pd.read_csv(ROOT/'motifs/TF_extension/chromVAR/motifs_independent_limma.csv')
joined=df.merge(cv[['id','motif_name','archetype','mapping_status','LOY_minus_CTRL_zscore','adj.P.Val']],left_on='motif_id',right_on='id',how='left',validate='one_to_one')
assert joined['id'].notna().all()
joined.to_csv(OUT/'TF_scores_with_chromVAR.csv',index=False)
fig,ax=plt.subplots(figsize=(9,6))
plot=joined[valid]
ax.scatter(plot.LOY_CTRL_change,-np.log10(plot.LOY_CTRL_pvalue.clip(lower=1e-300)),
 c=np.where(plot.pooled_site_BH_FDR<.05,'#0072b2','#aaaaaa'),s=14,alpha=.6)
ax.axvline(0,color='grey',ls='--',lw=.8)
labels=plot.sort_values('pooled_site_BH_FDR').drop_duplicates('archetype').head(12)
ax.set(xlabel='TOBIAS differential score (LOY versus CTRL)',ylabel='-log10 pooled-site P-value',title='Predicted differential motif occupancy: TOBIAS')
ax.legend(handles=[Line2D([],[],marker='o',ls='',color='#0072b2',label='Pooled-site BH FDR <0.05'),
                   Line2D([],[],marker='o',ls='',color='#aaaaaa',label='Other tested motifs')],frameon=False,fontsize=8)
fig.text(.5,.015,'Condition pools have comparable depth. P-values describe motif-site score differences; they are not biological-replicate tests.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.04,1,1))
label_points(ax,[(r.motif_name,r.LOY_CTRL_change,-np.log10(max(r.LOY_CTRL_pvalue,1e-300))) for r in labels.itertuples()])
for ext in ['png','pdf']:fig.savefig(OUT/f'TOBIAS_differential_volcano.{ext}',dpi=220)
plt.close(fig)
# One representative per mapped group, selected by chromVAR FDR before examining agreement.
representatives=joined.sort_values('adj.P.Val').drop_duplicates('archetype')
fig,ax=plt.subplots(figsize=(9,6))
ax.scatter(representatives.LOY_minus_CTRL_zscore,representatives.LOY_CTRL_change,
 c=np.where(representatives['adj.P.Val']<.05,'#0072b2','#aaaaaa'),s=18,alpha=.65)
ax.axvline(0,color='grey',ls='--',lw=.8);ax.axhline(0,color='grey',ls='--',lw=.8)
ax.set(xlabel='chromVAR z-score difference (LOY - CTRL)',ylabel='TOBIAS differential score (LOY versus CTRL)',title='Motif accessibility versus predicted occupancy')
ax.legend(handles=[Line2D([],[],marker='o',ls='',color='#0072b2',label='chromVAR FDR <0.05'),
                   Line2D([],[],marker='o',ls='',color='#aaaaaa',label='Other representatives')],frameon=False,fontsize=8)
fig.text(.5,.015,'One representative motif per mapped group; unresolved motifs remain separate. Shared data mean agreement is not independent validation.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.04,1,1))
label_points(ax,[(r.motif_name,r.LOY_minus_CTRL_zscore,r.LOY_CTRL_change) for r in representatives.head(12).itertuples()])
for ext in ['png','pdf']:fig.savefig(ROOT/f'motifs/TF_extension/concordance/chromVAR_TOBIAS.{ext}',dpi=220)
rna=pd.read_csv(ROOT/'motifs/TF_extension/concordance/chromVAR_TF_RNA_all.csv')
combined=rna.merge(joined[['motif_id','LOY_CTRL_change','LOY_CTRL_pvalue','pooled_site_BH_FDR']],on='motif_id',how='left',validate='many_to_one')
homer=pd.read_csv(ROOT/'motifs/TF_extension/HOMER/known_motif_all_results.csv')
for direction in ['gained','lost']:
    h=homer[homer.arm.eq('primary_TMM')&homer.direction.eq(direction)][['motif_id','fold_enrichment','PValue','BH_FDR_both_directions']]
    h=h.rename(columns={c:f'HOMER_{direction}_{c}' for c in h if c!='motif_id'})
    combined=combined.merge(h,on='motif_id',how='left',validate='many_to_one')
combined.to_csv(ROOT/'motifs/TF_extension/concordance/TF_ATAC_footprint_RNA.csv',index=False)
joint_ids=combined.loc[(combined.chromVAR_FDR<.05)&((combined.HOMER_gained_BH_FDR_both_directions<.05)|(combined.HOMER_lost_BH_FDR_both_directions<.05)),'motif_id']
combined[combined.motif_id.isin(joint_ids)].to_csv(ROOT/'motifs/TF_extension/concordance/shared_HOMER_chromVAR_candidates_with_footprints_RNA.csv',index=False)
summary=dict(motifs=len(df),finite_comparisons=int(valid.sum()),pooled_site_FDR005=int((df.pooled_site_BH_FDR<.05).sum()),
 sample_level_inference=False,score_comparison='LOY_CTRL',predicted_occupancy=True)
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'METHODS.txt').write_text('TOBIAS 0.17.5 ATACorrect -> ScoreBigwig -> BINDetect. JASPAR2024 CORE vertebrates; motif site P=5e-5.\nUnshifted, filtered proper-pair BAMs were sampled consistently by read name to approximately equal depth per replicate, then pooled by condition. Actual retained counts are recorded. No second duplicate removal or manual Tn5 shift was applied.\nBias correction uses internal +4/-5 shifts and the verified LF-format hg38 genome; chrY/mitochondrial chromosomes are excluded. Footprint regions are the tested consensus excluding chrY/PAR. Default signal normalizations are retained.\nBINDetect compares LOY first with CTRL second. Its score and pooled-site P-values are retained; BH is added across finite motif comparisons. These P-values are not replicate-aware biological differential tests, and native highlighted flags are not FDR calls.\nFootprints predict occupancy; a motif-family signal does not identify a unique TF paralog or establish binding experimentally.\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
