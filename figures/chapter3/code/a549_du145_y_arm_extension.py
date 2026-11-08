"""Display-only 5-Mb tracks and Y arm extension for the eight-panel supplement."""
import numpy as np
import pandas as pd

def extend(fig,src,at,dt,cb,t,units,profiles,tables,inputs):
    yq=cb[cb.chrom.eq('chrY') & cb.band.eq('q12')].iloc[0]
    tracks={}; armplots={}; report={}; correlations=[]
    for line,ctrl,loy,path in [('A549','A549-CTRL2','A549-LOY',at/'Y_PAR_all_gene_effects.csv'),
                             ('DU145','LNCaP-sgCTRL','LNCaP-LOY-2X',dt/'DU145_LOY_vs_sgCTRL_all_gene_effects.csv')]:
        inputs.append(path)
        raw=pd.read_csv(fig/f'source_data/{line.lower()}_lpwgs_confirmation/qdnaseq_bins.csv')
        raw=raw[raw.use & raw['sample'].isin([ctrl,loy]) & np.isfinite(raw.log2)]
        paired=raw.pivot(index=['chromosome','start','end'],columns='sample',values='log2').dropna().reset_index()
        excluded=paired.chromosome.eq('chrY') & paired.start.le(yq.end) & paired.end.gt(yq.start)
        paired=paired[~excluded].copy()
        paired['window_start']=(((paired.start+paired.end)/2)//5_000_000).astype(int)*5_000_000
        b=paired.groupby(['chromosome','window_start']).agg(control=(ctrl,'median'),LOY=(loy,'median'),n_bins=('start','size'),source_start=('start','min'),source_end=('end','max')).reset_index()
        b['window_end']=np.minimum(b.window_start+5_000_000,b.source_end)
        b['mid']=(b.window_start+b.window_end)/2
        b['DNA']=b.LOY-b.control
        assert len(b)==630 and excluded.sum()==351
        b.to_csv(src/f'{line}_8panel_track_5Mb.csv',index=False);tracks[line]=b
        y=pd.read_csv(path)
        y=y[y.chrom.eq('chrY') & y.expressed & ~y.PAR].copy()
        y=y[['gene_id','gene_name','chrom','mean_difference']]
        y['stable']=y.gene_id.str.split('.').str[0]
        y=y.merge(t,on=['chrom','stable'],how='left',validate='many_to_one')
        y['mapping']=np.where(y.t2t_start.notna(),'stable_ENSG_same_chrom','unmapped')
        unique=t[~t.duplicated(['chrom','symbol'],keep=False)].set_index(['chrom','symbol'])
        for i in y.index[y.t2t_start.isna()]:
            key=(y.at[i,'chrom'],y.at[i,'gene_name'])
            if key in unique.index:
                z=unique.loc[key];y.loc[i,['t2t_start','t2t_end','mapping']]=[z.t2t_start,z.t2t_end,'symbol_same_chrom_unique']
        y['Yq12_excluded']=y.t2t_start.lt(yq.end) & y.t2t_end.gt(yq.start)
        y['mid']=(y.t2t_start+y.t2t_end)/2
        units(y,y.chrom,y.mid)
        y['included']=y.t2t_start.notna() & ~y.Yq12_excluded
        y.to_csv(src/f'{line}_8panel_Y_gene_mapping.csv',index=False)
        r=y[y.included].groupby('arm').agg(RNA=('mean_difference','mean'),genes=('gene_id','size'))
        dna=profiles[line]
        dna=dna[dna.chromosome.eq('chrY') & ~(dna.start.le(yq.end) & dna.end.gt(yq.start))]
        z=dna.groupby('arm').agg(DNA=('DNA','median'),bins=('DNA','size')).join(r).dropna()
        z=z[z.genes.ge(1)].reset_index().rename(columns={'arm':'unit'})
        assert set(z.unit)=={'chrYp','chrYq'}
        arms=pd.concat([tables[line]['arm'],z],ignore_index=True)
        assert len(arms)==43
        arms.to_csv(src/f'{line}_8panel_arm_Y_inclusive.csv',index=False);armplots[line]=arms
        correlations.append(dict(cell_line=line,n=43,Pearson=arms.DNA.corr(arms.RNA),nonY_n=41,nonY_Pearson=tables[line]['arm'].DNA.corr(tables[line]['arm'].RNA)))
        report[line]=dict(Y_mapping=y.mapping.value_counts().to_dict(),Yq12_excluded_genes=y.loc[y.Yq12_excluded,'gene_name'].tolist(),Y_arm_genes=z.set_index('unit').genes.to_dict(),excluded_track_bins=int(excluded.sum()))
    pd.DataFrame(correlations).to_csv(src/'8panel_Y_inclusive_correlations.csv',index=False)
    report['rules']={'track':'5-Mb windows; per-population medians of corrected unsegmented log2, then LOY minus control',
      'Yq12_BED':[int(yq.start),int(yq.end)],'Y_RNA':'expressed non-PAR Y genes, saved mean_difference on the same normalized scale; Yq12 genes excluded',
      'Y_DNA':'median paired saved segment log2 differences after Yq12 exclusion; same arm assignment as non-Y',
      'minimum_genes':'20 for non-Y arms, 1 for Y arms; non-Y values unchanged',
      'scope':'Y added only to arm panels d/h; gene and 10-Mb scatter panels retain their previous non-Y scope'}
    return tracks,armplots,report
