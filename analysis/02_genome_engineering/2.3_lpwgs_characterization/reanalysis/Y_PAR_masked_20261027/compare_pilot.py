"""One-sample descriptive PAR mapping/coverage comparison; no Y-retention inference."""
import bisect,collections,csv,hashlib,json,statistics,subprocess,sys
from pathlib import Path
import pysam

ROOT=Path(__file__).resolve().parent
SAMTOOLS=str(ROOT/'runtime/env/bin/samtools')
OLD=ROOT.parents[1]/'t2t_bam/A549-CTRL2.bam'
NEW=ROOT/'alignment/A549-CTRL2.YPARmasked.bam'
OUT=ROOT/'comparison';OUT.mkdir(exist_ok=True)
BED=ROOT/'reference/chm13v2.0_PAR.bed'
regions=[]
for line in BED.read_text().splitlines():
    c,s,e=line.split()[:3];s,e=int(s),int(e)
    regions.append((c,s,e,'PAR1' if s==0 else 'PAR2'))
EXCLUDE=0xF04 # unmapped, secondary, supplementary, QC failure, marked duplicates
def write_table(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def where(read):
    if read.is_unmapped:return 'unmapped'
    for c,s,e,p in regions:
        if read.reference_name==c and read.reference_start<e and read.reference_end>s:
            return f'{c}_{p}'
    return 'outside_PAR'
def key(read):return (read.query_name,1 if read.is_read1 else 2)

# Same autosomal eligibility mask for the within-sample before/after diagnostic.
annotation=ROOT.parents[2]/'SEPT_11_2026_THESIS/research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv'
autobed=OUT/'autosomal_usable_100kb.bed'
with annotation.open(newline='') as f,autobed.open('w') as o:
    count=0
    for row in csv.DictReader(f):
        if row['chromosome'] in map(str,range(1,23)) and row['use'].upper()=='TRUE':
            o.write(f"chr{row['chromosome']}\t{int(float(row['start']))-1}\t{int(float(row['end']))}\t{row['bin_id']}\n")
            count+=1
assert count==21446,count
coverage=[];autonorm=[];mapping=[];old_records={}
baseline_only='--baseline-only' in sys.argv[1:]
cache_path=OUT/'original_baseline_cache.json'
cache_signature={'bam_size':OLD.stat().st_size,'bam_mtime_ns':OLD.stat().st_mtime_ns,
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'PAR_bed_sha256':hashlib.sha256(BED.read_bytes()).hexdigest(),
    'bin_annotation_sha256':hashlib.sha256(annotation.read_bytes()).hexdigest()}
cached=False
if not baseline_only and cache_path.exists():
    cache=json.loads(cache_path.read_text())
    if cache['signature']==cache_signature:
        coverage=cache['coverage'];autonorm=cache['autonorm'];mapping=cache['mapping']
        old_records={(r['name'],r['mate']):r['data'] for r in cache['old_records']}
        cached=True
        print('Reusing validated original-alignment baseline.',flush=True)
alignments=[('original',OLD)] if baseline_only else [('original',OLD),('Y_PAR_masked',NEW)]
for label,path in alignments:
    if label=='original' and cached:continue
    print('Comparing:',label,flush=True)
    for q in (0,30,37):
        for chrom,start,end,par in regions:
            cmd=[SAMTOOLS,'coverage','-H','-q',str(q),'-Q','0','--ff',str(EXCLUDE),
                '-r',f'{chrom}:{start+1}-{end}',str(path)]
            result=subprocess.run(cmd,check=True,text=True,capture_output=True)
            fields=result.stdout.strip().split('\t')
            assert len(fields)==9,(cmd,result.stdout,result.stderr)
            row={'alignment':label,'region':par,'chromosome':chrom,'start_0based':start,'end':end,
                'minimum_MAPQ':q,'read_count':int(fields[3]),'covered_bases':int(fields[4]),
                'breadth_percent':float(fields[5]),'mean_depth':float(fields[6]),
                'mean_base_quality':float(fields[7]),'mean_mapping_quality':float(fields[8])}
            coverage.append(row)
    for q in (30,37):
        target=OUT/f'{label}.autosomal_bedcov.Q{q}.tsv'
        with target.open('w') as f:
            subprocess.run([SAMTOOLS,'bedcov','-j','-Q',str(q),'-G','0x800',str(autobed),str(path)],stdout=f,check=True)
        depths=[]
        for line in target.read_text().splitlines():
            c,s,e,name,total=line.split('\t')
            depths.append(int(total)/(int(e)-int(s)))
        assert len(depths)==count
        autonorm.append({'alignment':label,'minimum_MAPQ':q,'usable_autosomal_bins':count,
                        'median_autosomal_depth':statistics.median(depths)})
    with pysam.AlignmentFile(path,'rb',threads=4) as b:
        for chrom,start,end,par in regions:
            histogram=collections.Counter()
            for r in b.fetch(chrom,start,end):
                if r.flag&EXCLUDE:continue
                histogram[r.mapping_quality]+=1
                if label=='original':
                    assert key(r) not in old_records
                    old_records[key(r)]={'old_region':f'{chrom}_{par}','old_MAPQ':r.mapping_quality}
            for q,n in sorted(histogram.items()):
                mapping.append({'alignment':label,'region':f'{chrom}_{par}','MAPQ':q,'read_count':n})
if baseline_only:
    cache_path.write_text(json.dumps({'signature':cache_signature,'coverage':coverage,
        'autonorm':autonorm,'mapping':mapping,
        'old_records':[{'name':k[0],'mate':k[1],'data':v} for k,v in old_records.items()]}))
    print(f'Original baseline ready: {len(old_records):,} PAR read ends to follow.',flush=True)
    sys.exit(0)
write_table(OUT/'regional_coverage.csv',coverage)
write_table(OUT/'autosomal_normalization.csv',autonorm)
write_table(OUT/'PAR_MAPQ_distribution.csv',mapping)

normal={(r['alignment'],r['minimum_MAPQ']):r['median_autosomal_depth'] for r in autonorm}
combined=[]
for label in ('original','Y_PAR_masked'):
    for par in ('PAR1','PAR2'):
        xspan=next(e-s for c,s,e,p in regions if c=='chrX' and p==par)
        for q in (0,30,37):
            rr=[r for r in coverage if r['alignment']==label and r['region']==par and r['minimum_MAPQ']==q]
            bases=sum(r['mean_depth']*(r['end']-r['start_0based']) for r in rr)
            combined.append({'alignment':label,'region':par,'minimum_MAPQ':q,
                'X_plus_Y_read_count':sum(r['read_count'] for r in rr),
                'X_plus_Y_aligned_bases_per_X_PAR_base':bases/xspan,
                'relative_to_median_autosomal_depth':bases/xspan/normal[(label,q)] if q>0 else ''})
write_table(OUT/'combined_PAR_coverage.csv',combined)

# Follow each old PAR primary read end through the complete realignment.
transitions=collections.Counter();followed=0
with pysam.AlignmentFile(NEW,'rb',threads=4) as bam:
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary:continue
        old=old_records.get(key(r))
        if old is None:continue
        transitions[(old['old_region'],old['old_MAPQ'],where(r),r.mapping_quality)]+=1
        followed+=1
assert followed==len(old_records),(followed,len(old_records))
rows=[{'old_region':a,'old_MAPQ':b,'new_region':c,'new_MAPQ':d,'read_count':n}
      for (a,b,c,d),n in sorted(transitions.items())]
write_table(OUT/'PAR_read_reassignment.csv',rows)
summary={}
for c,s,e,p in regions:
    name=f'{c}_{p}'
    rr=[r for r in rows if r['old_region']==name]
    n=sum(r['read_count'] for r in rr)
    summary[name]={'original_PAR_read_ends':n,
        'new_X_PAR_read_ends':sum(r['read_count'] for r in rr if r['new_region']==f'chrX_{p}'),
        'old_MAPQ_ge_37':sum(r['read_count'] for r in rr if r['old_MAPQ']>=37),
        'new_MAPQ_ge_37_anywhere':sum(r['read_count'] for r in rr if r['new_MAPQ']>=37),
        'old_below37_new_X_PAR_ge37':sum(r['read_count'] for r in rr if r['old_MAPQ']<37 and r['new_MAPQ']>=37 and r['new_region']==f'chrX_{p}'),
        'new_unmapped':sum(r['read_count'] for r in rr if r['new_region']=='unmapped')}
(OUT/'read_reassignment_summary.json').write_text(json.dumps(summary,indent=2))
(OUT/'analysis_parameters.json').write_text(json.dumps({'sample':'A549-CTRL2','minimum_base_quality':0,
    'MAPQ_thresholds':[0,30,37],'excluded_flags':EXCLUDE,
    'coverage_definition':'Aligned-base depth; overlapping mates can both contribute; deletions/skips excluded from bedcov.',
    'normalizer':'Median depth of 21,446 original-mask usable autosomal 100-kb bins; same bins both alignments.',
    'interpretation':'Combined homologous PAR mapping/coverage diagnostic, not Y-specific retention, absolute copy number, or a LOY/control comparison.',
    'PAR_lengths_differ_between_X_and_Y':True,
    'old_PAR_read_ends_followed':followed},indent=2))
print(json.dumps(summary,indent=2),flush=True)
