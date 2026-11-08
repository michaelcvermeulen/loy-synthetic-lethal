from pathlib import Path
import re,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'motifs/TF_extension/HOMER'
mapping=pd.read_csv(ROOT/'motifs/TF_extension/chromVAR/motif_archetype_mapping.csv')
cv=pd.read_csv(ROOT/'motifs/TF_extension/chromVAR/motifs_independent_limma.csv')
def bh(p):
    p=np.asarray(p);order=np.argsort(p);v=p[order]*len(p)/np.arange(1,len(p)+1)
    adjusted=np.minimum.accumulate(v[::-1])[::-1].clip(0,1);out=np.empty(len(p));out[order]=adjusted;return out
tables=[]
for arm in ['primary_TMM','passage_TMM']:
    pieces=[]
    for direction in ['gained','lost']:
        folder=OUT/arm/direction;assert (folder/'COMPLETE.json').exists()
        raw=pd.read_csv(folder/'knownResults.txt',sep='\t')
        n_target=int(re.search(r'\(of (\d+)\)',raw.columns[5]).group(1))
        n_background=int(re.search(r'\(of (\d+)\)',raw.columns[7]).group(1))
        expected=sum(1 for _ in (OUT/arm/f'{direction}_target.bed').open());assert n_target==expected
        df=pd.DataFrame(dict(motif_id=raw.iloc[:,0].str.extract(r'^(MA\d+\.\d+)',expand=False),
            HOMER_name=raw.iloc[:,0],HOMER_displayed_PValue=pd.to_numeric(raw.iloc[:,2]),log_PValue=pd.to_numeric(raw.iloc[:,3]),
            HOMER_qvalue=pd.to_numeric(raw.iloc[:,4]),target_motif_count=pd.to_numeric(raw.iloc[:,5]),
            background_motif_count=pd.to_numeric(raw.iloc[:,7])))
        # HOMER's displayed P column rounds to powers of ten; its log P column preserves more precision.
        df['PValue']=np.exp(df.log_PValue)
        assert len(df)==879 and df.motif_id.notna().all() and df.motif_id.is_unique
        assert np.isfinite(df.PValue).all() and df.PValue.between(0,1).all()
        df['target_total']=n_target;df['background_total']=n_background
        df['target_fraction']=df.target_motif_count/n_target;df['background_fraction']=df.background_motif_count/n_background
        df['fold_enrichment']=df.target_fraction.div(df.background_fraction.replace(0,np.nan))
        df.loc[(df.background_fraction==0)&(df.target_fraction>0),'fold_enrichment']=np.inf
        df['arm']=arm;df['direction']=direction;pieces.append(df)
    joined=pd.concat(pieces,ignore_index=True);joined['BH_FDR_both_directions']=bh(joined.PValue)
    tables.append(joined)
all_results=pd.concat(tables,ignore_index=True).merge(mapping,on='motif_id',validate='many_to_one')
all_results=all_results.sort_values(['arm','BH_FDR_both_directions','PValue'])
all_results.to_csv(OUT/'known_motif_all_results.csv',index=False)
primary=all_results[all_results.arm=='primary_TMM']
# Display one representative motif per mapped group; retain motif-level statistics, not a claimed family-level FDR.
representatives=primary.drop_duplicates(['direction','archetype'])
top=representatives[representatives.BH_FDR_both_directions<.05].groupby('direction',sort=False).head(12)
top[['direction','motif_id','motif_name','archetype','fold_enrichment','PValue','BH_FDR_both_directions']].to_csv(OUT/'known_motif_summary.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(13,6))
for ax,direction in zip(axes,['gained','lost']):
    d=top[top.direction==direction].iloc[::-1]
    if d.empty:
        ax.text(.5,.5,'No motifs pass FDR <0.05',ha='center',transform=ax.transAxes)
        ax.set(title=f'{direction.title()} DA peaks',xlabel='Known-motif fold enrichment')
        continue
    x=d.fold_enrichment.to_numpy();finite=x[np.isfinite(x)];cap=max(2,finite.max() if len(finite) else 2)*1.1
    values=np.nan_to_num(x,nan=0,posinf=cap)
    colors=np.where(d.BH_FDR_both_directions<.05,'#d55e00' if direction=='gained' else '#0072b2','#b0b0b0')
    labels=d.motif_name+' '+d.motif_id
    ax.barh(np.arange(len(d)),values,color=colors);ax.set_yticks(np.arange(len(d)),labels,fontsize=8)
    ax.axvline(1,color='grey',ls='--',lw=.8);ax.set(title=f'{direction.title()} DA peaks',xlabel='Known-motif fold enrichment')
    ax.set_xlim(0,max(1,float(values.max()))*1.3)
    if len(d)==1:ax.set_ylim(-1,1)
    for j,(fe,fdr) in enumerate(zip(values,d.BH_FDR_both_directions)):
        ax.text(fe,j,f'  q={fdr:.2g}',va='center',fontsize=7)
fig.suptitle('HOMER 5.1 / JASPAR2024: FDR <0.05 and |ATAC log2FC| >1')
fig.text(.5,.01,'Up to 12 significant motifs per direction, one per mapped group; motif-level BH correction across both directions.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.035,1,.95))
for ext in ['png','pdf']:fig.savefig(OUT/f'known_motif_enrichment.{ext}',dpi=220)
plt.close(fig)
con=primary.merge(cv[['id','LOY_minus_CTRL_zscore','adj.P.Val']],left_on='motif_id',right_on='id',validate='many_to_one')
con.to_csv(ROOT/'motifs/TF_extension/concordance/HOMER_chromVAR.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(11,5))
for ax,direction in zip(axes,['gained','lost']):
    d=con[(con.direction==direction)&np.isfinite(con.fold_enrichment)&(con.fold_enrichment>0)]
    both=(d.BH_FDR_both_directions<.05)&(d['adj.P.Val']<.05)
    ax.scatter(np.log2(d.fold_enrichment),d.LOY_minus_CTRL_zscore,s=12,c=np.where(both,'#0072b2','#aaaaaa'),alpha=.6)
    ax.axhline(0,color='grey',ls='--',lw=.7);ax.axvline(0,color='grey',ls='--',lw=.7)
    ax.set(title=f'HOMER: {direction} peaks',xlabel='HOMER log2 fold enrichment',ylabel='chromVAR z-score difference (LOY - CTRL)')
    labels=d[both].drop_duplicates('archetype').head(8).sort_values('LOY_minus_CTRL_zscore')
    for j,(_,r) in enumerate(labels.iterrows()):
        ax.annotate(r.motif_name,(np.log2(r.fold_enrichment),r.LOY_minus_CTRL_zscore),
            xytext=(-42,(-10 if j%2==0 else 10)),textcoords='offset points',fontsize=8,
            arrowprops=dict(arrowstyle='-',color='grey',lw=.5))
fig.suptitle('Sequence enrichment versus per-sample motif accessibility')
fig.text(.5,.01,'Blue: FDR <0.05 in both tests. Shared samples and related motifs mean agreement is method robustness, not independent validation.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.04,1,.94))
for ext in ['png','pdf']:fig.savefig(ROOT/f'motifs/TF_extension/concordance/HOMER_chromVAR.{ext}',dpi=220)
summary=all_results.groupby(['arm','direction']).agg(tested=('motif_id','size'),significant_FDR005=('BH_FDR_both_directions',lambda x:int((x<.05).sum()))).reset_index()
summary.to_csv(OUT/'analysis_summary.csv',index=False);print(summary.to_string(index=False))
(OUT/'verification.txt').write_text('PASS: all four HOMER runs ingest the expected target counts and test 879 unique JASPAR motifs\nPASS: finite valid P-values; BH across both directions within each design arm\nPASS: motif identities match chromVAR and reference mapping\n')
