"""Recompute 50-kb annotation and primary-read track coverage from existing inputs."""
from pathlib import Path
import csv,gzip,json,math,subprocess,hashlib
import pysam

ROOT=Path(__file__).resolve().parent
WGS=ROOT.parents[1]
OUT=ROOT/'a549_50kb_comparison';OUT.mkdir(exist_ok=True)
REF=Path('/home/micha/codex-work/wgs-ypar-pilot-20261027/reference/chm13v2.0.YPARmasked.fa')
SIZE=50000
chroms=['chr'+str(i) for i in range(1,23)]+['chrX','chrY']
fasta=pysam.FastaFile(str(REF))
rows=[];bychrom={}
for chrom in chroms:
    sequence=fasta.fetch(chrom).upper()
    rr=[]
    for s in range(0,len(sequence),SIZE):
        e=min(s+SIZE,len(sequence));seq=sequence[s:e]
        n=sum(seq.count(base) for base in 'ACGT');gc=seq.count('G')+seq.count('C')
        row=dict(bin_id=f'{chrom[3:]}:{s+1}-{e}',chromosome=chrom[3:],start=s+1,end=e,
            bases=100*n/(e-s),gc=100*gc/max(n,1),mappability=0.,blacklist=0.)
        rr.append(row);rows.append(row)
    bychrom[chrom]=rr
    print('Sequence annotation:',chrom,len(rr),flush=True)
fasta.close()
for filename,field in [('mappability_T2T.gz','mappability'),('difficulty_regions_T2T.gz','blacklist')]:
    intervals={c:[] for c in chroms}
    with gzip.open(WGS/filename,'rt') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            c,s,e,*_=line.split();c=c if c.startswith('chr') else 'chr'+c
            if c in intervals:intervals[c].append((int(s),int(e)))
    for chrom in chroms:
        data=sorted(intervals.pop(chrom));merged=[]
        for s,e in data:
            if merged and s<=merged[-1][1]:merged[-1]=(merged[-1][0],max(e,merged[-1][1]))
            else:merged.append((s,e))
        rr=bychrom[chrom]
        for s,e in merged:
            for j in range(s//SIZE,min((e-1)//SIZE+1,len(rr))):
                row=rr[j]
                overlap=max(0,min(e,row['end'])-max(s,row['start']-1))
                row[field]+=100*overlap/(row['end']-row['start']+1)
        assert all(-1e-8<=r[field]<=100+1e-7 for r in rr)
    print('Completed interval annotation:',field,flush=True)
with (OUT/'annotation_50kb.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

# A meaningful coordinate/provenance check: recombining the exact 50-kb
# annotation must reproduce the frozen 100-kb annotation outside masked Y-PAR.
oldpath=WGS.parent/'SEPT_11_2026_THESIS/research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv'
errors=[];checked=0;max_errors={x:0. for x in ['bases','gc','mappability','blacklist']}
with oldpath.open(newline='') as f:
    for old in csv.DictReader(f):
        chrom='chr'+old['chromosome'];s=int(float(old['start']))-1;e=int(float(old['end']))
        if chrom=='chrY' and (s<2458320 or e>62122809):continue
        rr=bychrom[chrom][s//SIZE:(e-1)//SIZE+1]
        assert sum(r['end']-r['start']+1 for r in rr)==e-s
        calculated={}
        for field in ['bases','mappability','blacklist']:
            calculated[field]=sum(r[field]*(r['end']-r['start']+1) for r in rr)/(e-s)
        total_nonN=sum(r['bases']*(r['end']-r['start']+1) for r in rr)
        calculated['gc']=sum(r['gc']*r['bases']*(r['end']-r['start']+1) for r in rr)/total_nonN if total_nonN else 0
        for field,value in calculated.items():
            err=abs(value-float(old[field]));max_errors[field]=max(max_errors[field],err)
            if err>1e-5:errors.append((old['bin_id'],field,err))
        checked+=1
(OUT/'annotation_validation.json').write_text(json.dumps({'checked_100kb_bins':checked,'maximum_absolute_errors':max_errors,'mismatches':errors[:30],'n_mismatches':len(errors)},indent=2))
assert not errors,errors[:10]
print('50-kb annotation reproduces original 100-kb annotation:',checked,'bins',flush=True)

bed=OUT/'track_bins_50kb.bed'
regions=[('chrX',0,2394410,'PAR1_X'),('chrX',153925834,154259566,'PAR2_X'),('chrY',0,62460029,'Y')]
with bed.open('w') as f:
    for chrom,s,e,name in regions:
        for left in range(s,e,SIZE):
            right=min(left+SIZE,e)
            f.write(f'{chrom}\t{left}\t{right}\t{name}:{left}-{right}\n')
coverage=[]
for sample in ['A549-CTRL2','A549-LOY']:
    bam=ROOT/'alignment'/f'{sample}.YPARmasked.bam'
    target=OUT/f'{sample}.bedcov.tsv'
    with target.open('w') as f:
        subprocess.run([str(ROOT/'runtime/env/bin/samtools'),'bedcov','-j','-c','-Q','37','-G','0x800',str(bed),str(bam)],stdout=f,check=True)
    for line in target.read_text().splitlines():
        c,s,e,name,total,count=line.split('\t');s,e=int(s),int(e)
        coverage.append(dict(sample=sample,chromosome=c,start0=s,end=e,track=name.split(':')[0],
            aligned_bases=int(total),read_count=int(count),depth=int(total)/(e-s)))
with (OUT/'binned_coverage.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(coverage[0]));w.writeheader();w.writerows(coverage)
(OUT/'provenance.json').write_text(json.dumps({'bin_size_bp':SIZE,'reference':str(REF),'masking':'Y-PAR only, same aligned reference',
    'annotation_sources':{name:hashlib.sha256((WGS/name).read_bytes()).hexdigest() for name in ['mappability_T2T.gz','difficulty_regions_T2T.gz']},
    'coverage_MAPQ':37,'coverage_excluded_flags':'0xF04','coverage_control_filter':'Evaluated after counting; never filter on LOY abundance.'},indent=2))
print('50-kb annotation and coverage complete',flush=True)
