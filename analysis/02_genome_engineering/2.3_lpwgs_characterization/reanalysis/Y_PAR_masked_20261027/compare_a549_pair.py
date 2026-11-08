"""Descriptive PAR coverage for the matched A549 pair, separate from Y-specific CN."""
import csv, hashlib, json, statistics, subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SAMTOOLS=str(ROOT/'runtime/env/bin/samtools')
OUT=ROOT/'a549_pair_comparison'; OUT.mkdir(exist_ok=True)
ANNOTATION=ROOT.parents[2]/'SEPT_11_2026_THESIS/research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv'
BED=ROOT/'reference/chm13v2.0_PAR.bed'
EXCLUDE=0xF04
regions=[]
for line in BED.read_text().splitlines():
    c,s,e=line.split()[:3]; s,e=int(s),int(e)
    regions.append((c,s,e,'PAR1' if s==0 else 'PAR2'))
def write_csv(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

assert (ROOT/'qc/REALIGNMENT_COMPLETE.txt').exists()
for sample in ['A549-CTRL2','A549-LOY']:
    bam=ROOT/'alignment'/f'{sample}.YPARmasked.bam'
    assert Path(str(bam)+'.bai').exists()
    subprocess.run([SAMTOOLS,'quickcheck',str(bam)],check=True)
# Coverage can be calculated while the complete read-content QC finishes.
# The figure/report requires both full read-content checks to pass before release.
autobed=OUT/'autosomal_usable_100kb.bed'
with ANNOTATION.open(newline='') as f,autobed.open('w') as o:
    n=0
    for r in csv.DictReader(f):
        if r['chromosome'] in map(str,range(1,23)) and r['use'].upper()=='TRUE':
            o.write(f"chr{r['chromosome']}\t{int(float(r['start']))-1}\t{int(float(r['end']))}\t{r['bin_id']}\n");n+=1
assert n==21446
coverage=[];norms=[];manifest=[]
for sample in ['A549-CTRL2','A549-LOY']:
    for alignment in ['original','Y_PAR_masked']:
        path=ROOT.parents[1]/'t2t_bam'/f'{sample}.bam' if alignment=='original' else ROOT/'alignment'/f'{sample}.YPARmasked.bam'
        print('Coverage:',sample,alignment,flush=True)
        manifest.append({'sample':sample,'alignment':alignment,'path':str(path),'size':path.stat().st_size})
        for chrom,start,end,par in regions:
            result=subprocess.run([SAMTOOLS,'coverage','-H','-q','37','-Q','0','--ff',str(EXCLUDE),
                '-r',f'{chrom}:{start+1}-{end}',str(path)],check=True,text=True,capture_output=True)
            f=result.stdout.strip().split('\t');assert len(f)==9
            coverage.append({'sample':sample,'alignment':alignment,'region':par,'chromosome':chrom,
                'start_0based':start,'end':end,'minimum_MAPQ':37,'read_count':int(f[3]),
                'covered_bases':int(f[4]),'breadth_percent':float(f[5]),'mean_depth':float(f[6])})
        target=OUT/f'{sample}.{alignment}.autosomal_bedcov.Q37.tsv'
        with target.open('w') as f:
            subprocess.run([SAMTOOLS,'bedcov','-j','-Q','37','-G','0x800',str(autobed),str(path)],stdout=f,check=True)
        depths=[int(t)/(int(e)-int(s)) for c,s,e,name,t in (line.split('\t') for line in target.read_text().splitlines())]
        assert len(depths)==n
        norms.append({'sample':sample,'alignment':alignment,'minimum_MAPQ':37,'n_autosomal_bins':n,
                      'median_autosomal_depth':statistics.median(depths)})
write_csv(OUT/'regional_coverage.csv',coverage)
write_csv(OUT/'autosomal_normalization.csv',norms)
normal={(r['sample'],r['alignment']):r['median_autosomal_depth'] for r in norms}
combined=[]
for sample in ['A549-CTRL2','A549-LOY']:
    for alignment in ['original','Y_PAR_masked']:
        for par in ['PAR1','PAR2']:
            rr=[r for r in coverage if r['sample']==sample and r['alignment']==alignment and r['region']==par]
            xlength=next(e-s for c,s,e,p in regions if c=='chrX' and p==par)
            depth=sum(r['mean_depth']*(r['end']-r['start_0based']) for r in rr)/xlength
            if alignment=='Y_PAR_masked':
                assert next(r['mean_depth'] for r in rr if r['chromosome']=='chrY')==0
            combined.append({'sample':sample,'alignment':alignment,'region':par,
                'X_plus_Y_read_count':sum(r['read_count'] for r in rr),'aligned_bases_per_X_PAR_base':depth,
                'relative_to_autosomal_median':depth/normal[(sample,alignment)]})
write_csv(OUT/'combined_PAR_coverage.csv',combined)
ratios=[]
for alignment in ['original','Y_PAR_masked']:
    for par in ['PAR1','PAR2']:
        ctrl=next(r for r in combined if r['sample']=='A549-CTRL2' and r['alignment']==alignment and r['region']==par)
        loy=next(r for r in combined if r['sample']=='A549-LOY' and r['alignment']==alignment and r['region']==par)
        ratios.append({'alignment':alignment,'region':par,'control_relative_depth':ctrl['relative_to_autosomal_median'],
            'LOY_relative_depth':loy['relative_to_autosomal_median'],
            'combined_PAR_percent_of_control':100*loy['relative_to_autosomal_median']/ctrl['relative_to_autosomal_median']})
write_csv(OUT/'PAR_ratios.csv',ratios)
(OUT/'provenance.json').write_text(json.dumps({'BAMs':manifest,'min_MAPQ':37,'excluded_flags':hex(EXCLUDE),
    'method':'Sum X and Y PAR aligned bases, divide by X PAR length and median usable-autosomal depth; compare LOY/control.',
    'interpretation':'Combined X/Y PAR signal; cannot assign to Y or estimate Y-PAR retention.',
    'reference_PAR_bed_sha256':hashlib.sha256(BED.read_bytes()).hexdigest(),
    'annotation_sha256':hashlib.sha256(ANNOTATION.read_bytes()).hexdigest()},indent=2))
print(ratios,flush=True)
