from pathlib import Path
import importlib.util,collections,json,csv,numpy as np,pysam
R=Path(__file__).resolve().parents[3];run=R.parent/'ATACseq/SEPT_08_2026_REPROCESS';out=R/'figures/chapter3/source_data/ht1080_atac_chrY'
spec=importlib.util.spec_from_file_location('prep',run/'code/02_prepare.py');prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)
start,end=2685900,2746400
class LocalBam:
 def __init__(self,b):self.b=b
 def fetch(self,until_eof=True):return self.b.fetch('chrX',start-100000,end+100000)
with (out/'CD99_unmasked_100bp.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['gene','chrom','gene_start','gene_end','sample','start','end','CPM'])
 for s in ['CTRL_11','CTRL_12','CTRL_21','LOY_11','LOY_12','LOY_21']:
  coverage=np.zeros(end-start+1,dtype=np.int64)
  with pysam.AlignmentFile(str(run.parent/'filtered'/f'{s}.bam'),'rb') as b:
   for _,__,chrom,a,z in prep.eligible_pairs(LocalBam(b),{'chrX'},{},collections.defaultdict(int)):
    a=max(a,start);z=min(z,end)
    if z>a:coverage[a-start]+=1;coverage[z-start]-=1
  values=np.cumsum(coverage[:-1])*1e6/json.loads((run/'qc'/f'{s}_preparation.json').read_text())['nonY_fragments']
  for i in range(0,len(values),100):w.writerow(['CD99','chrX',2690988,2741312,s,start+i,start+i+100,float(values[i:i+100].mean())])
  print(s,flush=True)
