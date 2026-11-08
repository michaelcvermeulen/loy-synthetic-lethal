"""Export DA subsets and GC/accessibility matched consensus backgrounds."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'motifs/TF_extension'
OUT.mkdir(parents=True,exist_ok=True)
rng = np.random.default_rng(20260908)
cov = pd.read_csv(ROOT/'motifs/peak_covariates.csv')
accounts=[]
for arm in ['primary_TMM','passage_TMM']:
    df=pd.read_csv(ROOT/f'differential/{arm}_all.csv').merge(cov[['peak_id','GC','N_fraction']],on='peak_id',validate='one_to_one')
    df=df[df.N_fraction==0].copy()
    df['gc_bin']=np.minimum((df.GC/.05).astype(int),19)
    df['accessibility_bin']=pd.qcut(df.logCPM,10,labels=False,duplicates='drop')
    folder=OUT/'HOMER'/arm
    folder.mkdir(parents=True,exist_ok=True)
    selected=(df.FDR<.05)&(df.logFC.abs()>1)
    df.assign(requested_DA_subset=selected).to_csv(folder/'all_tested_peaks.csv',index=False)
    for direction in ['gained','lost']:
        q=df[selected & ((df.logFC>0) if direction=='gained' else (df.logFC<0))]
        b=df[df.FDR>=.2]
        qparts=[];bparts=[]
        ratios=[]
        for key,g in q.groupby(['gc_bin','accessibility_bin']):
            available=((b.gc_bin==key[0])&(b.accessibility_bin==key[1])).sum()
            if available: ratios.append(max(1,int(available)//len(g)))
        control_ratio=min([10]+ratios)
        for key,g in q.groupby(['gc_bin','accessibility_bin']):
            pool=b[(b.gc_bin==key[0])&(b.accessibility_bin==key[1])]
            # Up to ten controls per target without reuse within each direction.
            if pool.empty: continue
            keep=min(len(g),len(pool))
            qparts.append(g if keep==len(g) else g.iloc[rng.choice(len(g),keep,replace=False)])
            bparts.append(pool.iloc[rng.choice(len(pool),control_ratio*keep,replace=False)])
        qt=pd.concat(qparts) if qparts else q.iloc[:0]
        bt=pd.concat(bparts) if bparts else b.iloc[:0]
        assert set(qt.peak_id).isdisjoint(bt.peak_id)
        for data,label in [(qt,'target'),(bt,'background')]:
            data.to_csv(folder/f'{direction}_{label}.csv',index=False)
            data.assign(score=0,strand='+')[['chrom','start','end','peak_id','score','strand']].to_csv(folder/f'{direction}_{label}.bed',sep='\t',header=False,index=False)
        accounts.append(dict(arm=arm,direction=direction,eligible_target=len(q),matched_target=len(qt),unmatched_target=len(q)-len(qt),background=len(bt),controls_per_target=control_ratio))
    fig,axs=plt.subplots(1,2,figsize=(11,4.5))
    cols=np.where(selected,np.where(df.logFC>0,'#d55e00','#0072b2'),'#b5b5b5')
    axs[0].scatter(df.logFC,-np.log10(df.FDR.clip(lower=1e-300)),s=3,c=cols,alpha=.5,rasterized=True)
    axs[0].set(xlabel='ATAC log2 fold change (LOY / CTRL)',ylabel='-log10 FDR',title='DA volcano')
    axs[0].axhline(-np.log10(.05),color='grey',ls='--',lw=.7)
    for x in [-1,1]: axs[0].axvline(x,color='grey',ls='--',lw=.7)
    axs[1].scatter(df.logCPM,df.logFC,s=3,c=cols,alpha=.5,rasterized=True)
    axs[1].set(xlabel='Average logCPM',ylabel='ATAC log2 fold change (LOY / CTRL)',title='MA plot')
    for y in [-1,1]: axs[1].axhline(y,color='grey',ls='--',lw=.7)
    fig.suptitle(f'{arm}: colored peaks FDR <0.05 and |log2FC| >1')
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(folder/f'DA_volcano_MA.{ext}',dpi=200)
    plt.close(fig)
pd.DataFrame(accounts).to_csv(OUT/'HOMER/matching_accounting.csv',index=False)
peaks=pd.read_csv(ROOT/'differential/primary_TMM_all.csv')
peaks[['chrom','start','end','peak_id']].to_csv(OUT/'tested_consensus.bed',sep='\t',header=False,index=False)
files=list((ROOT/'config/TF_extension').glob('*.jaspar'))+list((ROOT/'config/TF_extension').glob('*.tsv'))
pd.DataFrame([dict(file=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]).to_csv(OUT/'reference_checksums.csv',index=False)
(OUT/'HOMER/METHODS.txt').write_text('Targets: edgeR FDR <0.05 and |log2FC| >1; no refitting. Passage adjustment remains sensitivity, not confirmed technical batch.\nControls: tested consensus FDR >=0.20, no N bases, same 5%-GC bins and logCPM deciles; up to 10 controls per target without replacement. Unsupported target strata are explicitly excluded.\nHOMER applies its additional sequence-composition normalization. Fold enrichment is target motif percentage / background motif percentage.\n',encoding='utf-8')
print(pd.DataFrame(accounts).to_string(index=False))
