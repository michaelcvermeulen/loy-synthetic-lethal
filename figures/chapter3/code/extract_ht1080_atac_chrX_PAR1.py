from pathlib import Path
import pyBigWig,numpy as np,csv
root=Path(__file__).resolve().parents[3]
run=root.parent/'ATACseq/SEPT_08_2026_REPROCESS'
out=root/'figures/chapter3/source_data/ht1080_atac_chrY'
with (out/'chrX_PAR1_10kb.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['sample','chrom','start','end','mean_fragment_CPM'])
 for sample in ['CTRL_11','CTRL_12','CTRL_21','LOY_11','LOY_12','LOY_21']:
  b=pyBigWig.open(str(run/'tracks'/f'{sample}.fragment_cpm.bw'))
  values=np.nan_to_num(b.values('chrX',10000,2781479,numpy=True))
  for offset in range(0,len(values),10000):
   w.writerow([sample,'chrX',10000+offset,min(20000+offset,2781479),float(values[offset:offset+10000].mean())])
  b.close()
print('Extracted chrX PAR1 for six libraries.')
