from pathlib import Path
import pyBigWig,numpy as np,csv
root=Path(__file__).resolve().parents[3]
run=root.parent/'ATACseq/SEPT_08_2026_REPROCESS'
out=root/'figures/chapter3/source_data/ht1080_atac_chrY'
out.mkdir(parents=True,exist_ok=True)
with (out/'coverage_10kb.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['sample','start','end','mean_fragment_CPM'])
 for sample in ['CTRL_11','CTRL_12','CTRL_21','LOY_11','LOY_12','LOY_21']:
  b=pyBigWig.open(str(run/'tracks'/f'{sample}.fragment_cpm.bw'))
  length=b.chroms()['chrY']
  # Chunked base-level reads: retain zero-coverage bases in bin means.
  for start in range(0,length,1000000):
   end=min(start+1000000,length)
   v=np.nan_to_num(b.values('chrY',start,end,numpy=True))
   for j in range(0,len(v),10000):
    w.writerow([sample,start+j,min(start+j+10000,end),float(v[j:j+10000].mean())])
  b.close();print(sample,flush=True)
