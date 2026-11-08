from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[3];F=R/'figures/chapter3';D=F/'source_data/ht1080_atac_copy_number';run=R.parent/'ATACseq/SEPT_08_2026_REPROCESS'
bands=pd.read_csv(F/'source_data/ht1080_profile_wgs_cn/reference_material/cytoBand_hg38.txt.gz',sep='\t',header=None,names=['chrom','start','end','band','stain'])
cen=bands[bands.stain.eq('acen')].groupby('chrom').agg(start=('start','min'),end=('end','max'))
def arms(df,start,end):
 mid=(df[start]+df[end])/2;res=pd.Series(index=df.index,dtype='object')
 for chrom in ['chr'+str(i) for i in range(1,23)]:
  mask=df.chrom.eq(chrom);a,b=cen.loc[chrom,['start','end']]
  res.loc[mask & (mid<a)]=chrom+'p';res.loc[mask & (mid>=b)]=chrom+'q'
 return res
lc=pd.read_csv(run/'counts/logCPM_primary.tsv',sep='\t');coord=lc.peak_id.str.extract(r'(chr[^:]+):(\d+)-(\d+)');lc['chrom']=coord[0];lc['start']=coord[1].astype(int);lc['end']=coord[2].astype(int);lc['arm']=arms(lc,'start','end')
out=pd.DataFrame()
for ps in ['11','12','21']:
 delta=lc['LOY_'+ps]-lc['CTRL_'+ps];out['ATAC_p'+ps]=delta.groupby(lc.arm).median()-delta.median()
out['ATAC_mean_shift']=out[['ATAC_p11','ATAC_p12','ATAC_p21']].mean(axis=1);out['ATAC_regions']=lc.groupby('arm').size()
g=pd.read_csv(F/'source_data/ht1080_rna_atac_integration/integration_genes.csv');ann=pd.read_csv(run/'config/genes_knownGene.tsv',sep='\t').rename(columns={'gene_id':'entrez_id'})
g=g.merge(ann[['entrez_id','start','end']],on='entrez_id',validate='one_to_one');g['arm']=arms(g,'start','end');out['RNA_median_logFC']=g.groupby('arm').rna_logFC.median()-g.rna_logFC.median();out['RNA_genes']=g.groupby('arm').size()
b=pd.read_csv(D/'WGS_hg38_100kb_bins.csv');b=b[b['use']&b.delta_segment_log2.notna()].copy();b['arm']=arms(b,'start0','end0');out['WGS_D0_median_segment_delta']=b.groupby('arm').delta_segment_log2.median();out['WGS_bins']=b.groupby('arm').size()
out.index.name='arm';out.to_csv(D/'arm_correlations_all.csv');out=out.dropna(subset=['ATAC_mean_shift','RNA_median_logFC','WGS_D0_median_segment_delta']);out.to_csv(D/'panel_b_arm_correlations.csv')
stats={}
for key,x,y in [('RNA_vs_WGS','WGS_D0_median_segment_delta','RNA_median_logFC'),('ATAC_vs_WGS','WGS_D0_median_segment_delta','ATAC_mean_shift'),('ATAC_vs_RNA','RNA_median_logFC','ATAC_mean_shift')]:
 v=spearmanr(out[x],out[y]);stats[key]=[float(v.statistic),float(v.pvalue)]
(D/'arm_correlation_statistics.json').write_text(json.dumps({'n':len(out),'tests':stats},indent=2));print(len(out),stats)
