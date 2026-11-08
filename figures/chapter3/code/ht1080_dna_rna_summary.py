"""CHM13 mapping and regional summaries of saved HT1080 expression effects."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

def summarize(thesis, out, paired, arms, cb, aggregate):
    root = thesis.parent
    rna = root/'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/tables/batchDE_batch_condition__combined.csv'
    annotation = root/'WGS/T2T_genes.gz'
    a = pd.read_csv(rna)
    assert len(a)==25526 and not a.gene_id.duplicated().any()
    assert not a.chrom.eq('chrY').any()
    ysource=rna.parent/'Y_PAR_all_gene_effects.csv'
    y=pd.read_csv(ysource)
    y=y[y.chrom.eq('chrY') & y.expressed & ~y.PAR].copy()
    a=pd.concat([a,y[a.columns.intersection(y.columns)]],ignore_index=True)
    assert not a.gene_id.duplicated().any()
    t = pd.read_csv(annotation,sep='\t',usecols=[0,1,2,12,20],low_memory=False)
    t.columns=['chrom','t2t_start','t2t_end','symbol','source_gene']
    t=t[t.source_gene.str.startswith('ENSG',na=False)].copy()
    t['stable']=t.source_gene.str.split('.').str[0]
    t=t.groupby(['chrom','stable'],sort=False).agg(t2t_start=('t2t_start','min'),t2t_end=('t2t_end','max'),symbol=('symbol','first')).reset_index()
    a['stable']=a.gene_id.str.split('.').str[0]
    a=a.merge(t,on=['chrom','stable'],how='left',validate='many_to_one')
    a['mapping']=np.where(a.t2t_start.notna(),'stable_ENSG_same_chrom','unmapped')
    unique=t[~t.duplicated(['chrom','symbol'],keep=False)].set_index(['chrom','symbol'])
    for i in a.index[a.t2t_start.isna()]:
        key=(a.at[i,'chrom'],a.at[i,'gene_name'])
        if key in unique.index:
            z=unique.loc[key]
            a.loc[i,['t2t_start','t2t_end','mapping']]=[z.t2t_start,z.t2t_end,'symbol_same_chrom_unique']
    a['mid']=(a.t2t_start+a.t2t_end)/2
    a['RNA']=a.mean_difference
    a['window_start']=(a.mid//1_000_000).astype('Int64')*1_000_000
    cen=cb[cb.stain.eq('acen')].groupby('chrom').agg(left=('start','min'),right=('end','max'))
    a['arm']=None
    p=a.mid.le(a.chrom.map(cen.left));q=a.mid.gt(a.chrom.map(cen.right))
    a.loc[p,'arm']=a.loc[p,'chrom'].str[3:]+'p'
    a.loc[q,'arm']=a.loc[q,'chrom'].str[3:]+'q'
    yq12=cb[cb.chrom.eq('chrY') & cb.band.eq('q12')].iloc[0]
    a['Yq12_excluded']=a.chrom.eq('chrY') & a.t2t_start.lt(yq12.end) & a.t2t_end.gt(yq12.start)
    a.to_csv(out/'RNA_gene_mapping.csv',index=False)
    mapped=a[a.t2t_start.notna() & ~a.Yq12_excluded].copy()
    assert np.isfinite(mapped.RNA).all()
    r=mapped.groupby(['chrom','window_start']).agg(RNA=('RNA','mean'),genes=('gene_id','size')).reset_index()
    b=aggregate(paired,1_000_000).rename(columns={'chromosome':'chrom','difference':'DNA'})
    windows=b.merge(r,on=['chrom','window_start'],how='left',validate='one_to_one')
    windows['included']=windows.genes.ge(5) | (windows.chrom.eq('chrY') & windows.genes.ge(1))
    windows.to_csv(out/'DNA_RNA_1Mb_all_windows.csv',index=False)
    windows=windows[windows.included].copy()
    ar=mapped.groupby('arm').agg(RNA=('RNA','mean'),genes=('gene_id','size')).reset_index()
    ya=[]
    for arm,left,right in [('Yp',0,cen.loc['chrY','left']),('Yq',cen.loc['chrY','right'],cb[cb.chrom.eq('chrY')].end.max())]:
        v=paired[paired.chromosome.eq('chrY') & paired.mid.gt(left) & paired.mid.le(right)]
        ya.append(dict(chromosome='chrY',arm=arm,n_bins=len(v),difference=v.difference.median(),low_mappability_arm=False))
    arms=pd.concat([arms,pd.DataFrame(ya)],ignore_index=True)
    ar=arms.rename(columns={'difference':'DNA'}).merge(ar,on='arm',how='left',validate='one_to_one')
    ar['included']=ar.genes.ge(20) | (ar.chromosome.eq('chrY') & ar.genes.ge(1))
    ar.to_csv(out/'DNA_RNA_all_arms.csv',index=False)
    ar=ar[ar.included].copy()
    stats=[]
    for name,z in [('window_1Mb',windows),('arm',ar)]:
        assert np.isfinite(z[['DNA','RNA']]).all().all()
        nony=z[~z['chrom' if name=='window_1Mb' else 'chromosome'].eq('chrY')]
        stats.append(dict(level=name,n=len(z),Pearson=z.DNA.corr(z.RNA),Spearman=z.DNA.corr(z.RNA,method='spearman'),nonY_n=len(nony),nonY_Pearson=nony.DNA.corr(nony.RNA)))
        z.to_csv(out/f'DNA_RNA_{name}_plotted.csv',index=False)
    pd.DataFrame(stats).to_csv(out/'DNA_RNA_correlations.csv',index=False)
    allw=pd.read_csv(out/'DNA_RNA_1Mb_all_windows.csv')
    sensitivity=[]
    for minimum in [1,5,10]:
        z=allw[allw.genes.ge(minimum) & ~allw.chrom.eq('chrY')]
        sensitivity.append(dict(minimum_genes=minimum,n=len(z),Pearson=z.DNA.corr(z.RNA)))
    pd.DataFrame(sensitivity).to_csv(out/'DNA_RNA_window_gene_count_sensitivity.csv',index=False)
    provenance=dict(inputs={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [rna,ysource,annotation]},
        RNA='Saved mean_difference: equal-weight two-period mean difference in TMM log2 CPM; expressed autosomal/non-PAR X genes plus expressed non-PAR Y genes; no DE refit',
        mapping=a.mapping.value_counts().to_dict(),Yq12_excluded_RNA_genes=a.loc[a.Yq12_excluded,'gene_name'].tolist(),
        mapping_rule='stable ENSG on same chromosome, then unique symbol on same chromosome; CHM13 transcript spans collapsed per stable gene',
        aggregation='RNA: mean mapped gene effect by midpoint; DNA windows: difference of per-population median corrected log2 values; DNA arms: median paired bin difference as in bar panel',
        minimum_genes={'window_1Mb_nonY':5,'arm_nonY':20,'Y_units':1},arm_centromeres='Excluded by midpoint, using same boundaries as DNA arm bars; Yq12 excluded from Y DNA',
        scope='Descriptive population comparison; RNA and WGS specimens not established as matched; genomic units not biological replicates',correlations=stats)
    return windows,ar,provenance
