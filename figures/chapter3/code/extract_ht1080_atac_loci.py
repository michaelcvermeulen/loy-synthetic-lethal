from pathlib import Path
import csv,numpy as np,pyBigWig
R=Path(__file__).resolve().parents[3];F=R/'figures/chapter3/source_data/ht1080_atac_chrY'
loci=[('P2RY8','chrX',1462572,1537187),('AKAP17A','chrX',1591604,1602532),('RPS4Y1','chrY',2841602,2932000),('DDX3Y','chrY',12904108,12920478)]
samples=['CTRL_11','CTRL_12','CTRL_21','LOY_11','LOY_12','LOY_21']
with (F/'locus_coverage_100bp.csv').open('w',newline='') as out:
 w=csv.writer(out);w.writerow(['gene','chrom','gene_start','gene_end','sample','start','end','CPM'])
 for sample in samples:
  b=pyBigWig.open(str(R.parent/'ATACseq/SEPT_08_2026_REPROCESS/tracks'/f'{sample}.fragment_cpm.bw'))
  for gene,chrom,gs,ge in loci:
   start=(gs-1)//100*100-5000;end=((ge+99)//100)*100+5000
   v=np.nan_to_num(b.values(chrom,start,end,numpy=True))
   for i in range(0,len(v),100):w.writerow([gene,chrom,gs,ge,sample,start+i,min(start+i+100,end),float(v[i:i+100].mean())])
  b.close()
print('Four loci extracted at 100 bp.')
