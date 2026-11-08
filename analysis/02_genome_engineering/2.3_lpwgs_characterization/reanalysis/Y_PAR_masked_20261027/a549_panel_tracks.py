"""Collect primary-read 100-kb Y and X-PAR coverage for the revised A549 preview."""
from pathlib import Path
import csv,json,subprocess

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'a549_panel_tracks';OUT.mkdir(exist_ok=True)
SAMTOOLS=ROOT/'runtime/env/bin/samtools'
regions=[('chrX',0,2394410,'PAR1_X'),('chrX',153925834,154259566,'PAR2_X'),('chrY',0,62460029,'Y')]
bed=OUT/'track_bins.bed'
with bed.open('w') as f:
    for c,s,e,name in regions:
        # PAR bins start at the PAR boundary, and the final shorter bin is
        # retained with its actual width recorded. No cross-boundary bins.
        for left in range(s,e,100000):
            right=min(left+100000,e)
            f.write(f'{c}\t{left}\t{right}\t{name}:{left}-{right}\n')
rows=[]
for sample in ['A549-CTRL2','A549-LOY']:
    bam=ROOT/'alignment'/f'{sample}.YPARmasked.bam'
    assert bam.exists() and Path(str(bam)+'.bai').exists()
    output=OUT/f'{sample}.bedcov.tsv'
    print('Binned primary-read coverage:',sample,flush=True)
    with output.open('w') as f:
        subprocess.run([str(SAMTOOLS),'bedcov','-j','-c','-Q','37','-G','0x800',str(bed),str(bam)],stdout=f,check=True)
    for line in output.read_text().splitlines():
        c,s,e,name,bases,count=line.split('\t')
        rows.append({'sample':sample,'chromosome':c,'start0':int(s),'end':int(e),
            'track':name.split(':')[0],'aligned_bases':int(bases),'read_count':int(count),
            'depth':int(bases)/(int(e)-int(s))})
with (OUT/'binned_coverage.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(OUT/'provenance.json').write_text(json.dumps({'MAPQ_min':37,'excluded_flags':'0xF04',
    'bin_size_bp':100000,'PAR_bins':'X coordinates, anchored at each PAR boundary; final bin may be shorter',
    'normalization':'Use the existing strict-mask autosomal median depth for each sample.',
    'interpretation':'X-mapped PAR coverage in these alignments is combined homologous X/Y signal.'},indent=2))
print('Coverage tables complete',flush=True)
