from pathlib import Path
import sys,json,pysam
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'motifs/TF_extension/TOBIAS'
condition=sys.argv[1];assert condition in ('CTRL','LOY')
samples=[f'{condition}_{p}' for p in ['11','12','21']]
entries=[json.loads((out/f'{s}_matched.complete.json').read_text()) for s in samples]
dest=out/f'{condition}_matched_pool.bam';done=dest.with_suffix('.complete.json')
if done.exists():sys.exit(0)
assert not dest.exists(), 'Partial or running pool exists'
sources=[str(out/f'{s}_matched.bam') for s in samples]
pysam.merge('-@','4','-r',str(dest),*sources,catch_stdout=False)
pysam.index('-@','4',str(dest))
with pysam.AlignmentFile(dest,'rb') as bam:actual=sum(x.mapped for x in bam.get_index_statistics())
assert actual==sum(e['actual_reads'] for e in entries)
done.write_text(json.dumps(dict(condition=condition,reads=actual,fragments=actual//2,samples=sources),indent=2)+'\n')
fa=ROOT/'config/TF_extension/hg38_primary.fa'
if not Path(str(fa)+'.fai').exists():pysam.faidx(str(fa))
print(condition,'pool verified:',actual,'reads',flush=True)
