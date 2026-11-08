"""Read-only coverage diagnostics for the four authorised atlas samples."""
import argparse,csv,hashlib,json,shutil,statistics,subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('cell',choices=['HT1080','DU145']);ap.add_argument('mode',choices=['original','masked']);a=ap.parse_args()
SAMPLES={'HT1080':['HT1080-sgCTRL-D0','HT1080-LOY-D0'],'DU145':['LNCaP-sgCTRL','LNCaP-LOY-2X']}[a.cell]
OUT=ROOT/'atlas_pair_coverage'/a.cell/a.mode;OUT.mkdir(parents=True,exist_ok=True)
SAMTOOLS=str(ROOT/'runtime/env/bin/samtools')
BED=ROOT/'a549_panel_tracks/track_bins.bed'
AUTO=ROOT/'a549_pair_comparison/autosomal_usable_100kb.bed'
PAR=ROOT/'reference/chm13v2.0_PAR.bed'
coverage=[];normalizers=[];regional=[];manifest=[]
def save(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for sample in SAMPLES:
    bam=ROOT.parents[1]/'t2t_bam'/f'{sample}.bam' if a.mode=='original' else ROOT/'alignment'/f'{sample}.YPARmasked.bam'
    assert bam.exists() and Path(str(bam)+'.bai').exists()
    subprocess.run([SAMTOOLS,'quickcheck',str(bam)],check=True)
    manifest.append(dict(sample=sample,path=str(bam),size=bam.stat().st_size))
    print(a.cell,a.mode,sample,'regional coverage',flush=True)
    target=OUT/f'{sample}.tracks.bedcov.tsv'
    with target.open('w') as f:
        subprocess.run([SAMTOOLS,'bedcov','-j','-c','-Q','37','-G','0x800',str(BED),str(bam)],stdout=f,check=True)
    for line in target.read_text().splitlines():
        c,s,e,name,total,count=line.split('\t');s,e=int(s),int(e)
        coverage.append(dict(sample=sample,chromosome=c,start0=s,end=e,track=name.split(':')[0],
            aligned_bases=int(total),read_count=int(count),depth=int(total)/(e-s)))
    for line in PAR.read_text().splitlines():
        c,s,e=line.split();s,e=int(s),int(e)
        result=subprocess.run([SAMTOOLS,'coverage','-H','-q','37','-Q','0','--ff','0xF04','-r',f'{c}:{s+1}-{e}',str(bam)],check=True,text=True,capture_output=True)
        v=result.stdout.strip().split('\t');assert len(v)==9
        regional.append(dict(sample=sample,alignment=a.mode,region='PAR1' if s==0 else 'PAR2',chromosome=c,start0=s,end=e,
            read_count=int(v[3]),covered_bases=int(v[4]),breadth_percent=float(v[5]),mean_depth=float(v[6])))
    print(a.cell,a.mode,sample,'autosomal normalization',flush=True)
    target=OUT/f'{sample}.autosomal.bedcov.tsv'
    with target.open('w') as f:
        subprocess.run([SAMTOOLS,'bedcov','-j','-Q','37','-G','0x800',str(AUTO),str(bam)],stdout=f,check=True)
    depths=[int(total)/(int(e)-int(s)) for c,s,e,name,total in (line.split('\t') for line in target.read_text().splitlines())]
    assert len(depths)==21446
    normalizers.append(dict(sample=sample,alignment=a.mode,minimum_MAPQ=37,n_autosomal_bins=len(depths),median_autosomal_depth=statistics.median(depths)))
save(OUT/'binned_coverage.csv',coverage)
save(OUT/'autosomal_normalization.csv',normalizers)
save(OUT/'regional_PAR_coverage.csv',regional)
(OUT/'provenance.json').write_text(json.dumps({'cell':a.cell,'mode':a.mode,'BAMs':manifest,'MAPQ':37,
    'excluded_flags':'0xF04','autosomal_bin_mask':'Original strict mask, fixed across all samples and alignments',
    'roi_bed_sha256':hashlib.sha256(BED.read_bytes()).hexdigest(),'autosomal_bed_sha256':hashlib.sha256(AUTO.read_bytes()).hexdigest()},indent=2))
print(a.cell,a.mode,'coverage complete',flush=True)
