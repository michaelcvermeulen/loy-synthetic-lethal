"""Pair-consistent, approximately equal-depth sampling per replicate for TOBIAS."""
from pathlib import Path
import json
import pandas as pd
import pysam
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'motifs/TF_extension/TOBIAS'; OUT.mkdir(parents=True,exist_ok=True)
qc=pd.read_csv(ROOT/'qc/qc_metrics_master.csv')
target=int(qc.nonY_fragments.min())
chroms=[f'chr{i}' for i in range(1,23)]+['chrX']
records=[]
for i,row in qc.iterrows():
    source=ROOT/f'bam/{row["sample"]}.bam'
    dest=OUT/f'{row["sample"]}_matched.bam'
    done=dest.with_suffix('.complete.json')
    if done.exists(): records.append(json.loads(done.read_text()));continue
    if dest.exists(): raise RuntimeError(f'Incomplete output exists: {dest}')
    fraction=target/int(row.nonY_fragments)
    print('Sampling',row['sample'],'fraction',fraction,flush=True)
    args=['-b','-@','4','-o',str(dest),'--subsample-seed',str(20260908+i),'--subsample',str(fraction),str(source),*chroms]
    pysam.view(*args,catch_stdout=False)
    pysam.index('-@','4',str(dest))
    with pysam.AlignmentFile(dest,'rb') as bam:
        actual=sum(x.mapped for x in bam.get_index_statistics() if x.contig in chroms)
    assert actual%2==0 and abs(actual/2-target)/target<.002
    entry=dict(sample=row['sample'],condition=row.condition,source_nonY_fragments=int(row.nonY_fragments),
               target_fragments=target,sampling_probability=fraction,seed=20260908+i,actual_reads=actual,actual_fragments=actual//2)
    done.write_text(json.dumps(entry,indent=2)+'\n');records.append(entry)
    pd.DataFrame(records).to_csv(OUT/'replicate_pool_accounting.csv',index=False)
for condition in ['CTRL','LOY']:
    dest=OUT/f'{condition}_matched_pool.bam';done=dest.with_suffix('.complete.json')
    if done.exists():continue
    if dest.exists():raise RuntimeError(f'Incomplete pool exists: {dest}')
    sources=[str(OUT/f'{s}_matched.bam') for s in qc.loc[qc.condition==condition,'sample']]
    print('Merging',condition,flush=True)
    pysam.merge('-@','4','-r',str(dest),*sources,catch_stdout=False)
    pysam.index('-@','4',str(dest))
    with pysam.AlignmentFile(dest,'rb') as bam:
        actual=sum(x.mapped for x in bam.get_index_statistics() if x.contig in chroms)
    expected=sum(r['actual_reads'] for r in records if r['condition']==condition)
    assert actual==expected
    done.write_text(json.dumps(dict(condition=condition,reads=actual,fragments=actual//2,samples=sources),indent=2)+'\n')
fa=ROOT/'config/TF_extension/hg38_primary.fa'
assert Path(str(fa)+'.complete').exists()
pysam.faidx(str(fa))
with pysam.FastaFile(str(fa)) as genome,pysam.AlignmentFile(str(ROOT/'bam/CTRL_11.bam'),'rb') as bam:
    for chrom in chroms: assert genome.get_reference_length(chrom)==bam.get_reference_length(chrom)
(OUT/'POOL_COMPLETE.json').write_text(json.dumps(dict(target_fragments_per_replicate=target,sampling='samtools pair-consistent subsampling; actual depths recorded',reference_lengths_verified=True),indent=2)+'\n')
