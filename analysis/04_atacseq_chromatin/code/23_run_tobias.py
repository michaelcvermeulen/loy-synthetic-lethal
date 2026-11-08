"""Run published TOBIAS stages on depth-matched condition pools."""
from pathlib import Path
import subprocess,json,time,sys,shutil,pysam,fcntl,os
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'motifs/TF_extension/TOBIAS'
BIN=Path('/tmp/atacseq-sept08-2026-tf-tools/tobias_venv/bin/TOBIAS')
reference=json.loads((ROOT/'config/TF_extension/reference_verified.json').read_text())
genome=Path(reference['linux_cache']);assert genome.exists()
peaks=ROOT/'motifs/TF_extension/tested_consensus.bed'
motifs=ROOT/'config/TF_extension/JASPAR2024_CORE_vertebrates.jaspar'
only_condition=sys.argv[1] if len(sys.argv)>1 else None
if only_condition:
    assert only_condition in ['CTRL','LOY']
    assert (OUT/f'{only_condition}_matched_pool.complete.json').exists()
else:
    assert (OUT/'POOL_COMPLETE.json').exists(), 'Footprint pools must be completed first'
assert BIN.exists()
cores=os.environ.get('TF_TOBIAS_CORES','8')
def run(stage,args,required):
    marker=OUT/f'{stage}.complete.json'
    if marker.exists():
        assert all(p.exists() and p.stat().st_size>0 for p in required)
        return
    cmd=[str(BIN),*map(str,args)]
    print('Running',stage,flush=True)
    (OUT/f'{stage}.command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    start=time.time()
    with (ROOT/f'logs/{stage}.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    assert all(p.exists() and p.stat().st_size>0 for p in required)
    marker.write_text(json.dumps(dict(elapsed_seconds=time.time()-start,command=cmd),indent=2)+'\n')
for condition in ([only_condition] if only_condition else ['CTRL','LOY']):
    # Separate condition jobs may overlap; the final comparison waits for each lock.
    condition_lock=(genome.parent/f'TOBIAS_{condition}.lock').open('w')
    fcntl.flock(condition_lock,fcntl.LOCK_EX)
    folder=OUT/condition
    corrected=folder/f'{condition}_corrected.bw'
    source_bam=OUT/f'{condition}_matched_pool.bam'
    cached_bam=genome.parent/source_bam.name
    if not (OUT/f'23_ATACorrect_{condition}.complete.json').exists():
        if not cached_bam.exists() or cached_bam.stat().st_size!=source_bam.stat().st_size:
            print('Caching verified pool',condition,flush=True);shutil.copyfile(source_bam,cached_bam)
        shutil.copyfile(Path(str(source_bam)+'.bai'),Path(str(cached_bam)+'.bai'))
        pysam.quickcheck(str(cached_bam))
        with pysam.AlignmentFile(str(cached_bam),'rb') as b:
            actual=sum(x.mapped for x in b.get_index_statistics())
        assert actual==json.loads(source_bam.with_suffix('.complete.json').read_text())['reads']
    run(f'23_ATACorrect_{condition}',['ATACorrect','--bam',cached_bam,
        '--genome',genome,'--peaks',peaks,'--blacklist',ROOT/'config/hg38-blacklist.v2.bed',
        '--drop-chroms','chrM','chrY','--outdir',folder,'--prefix',condition,'--cores',cores],[corrected,folder/f'{condition}_atacorrect.pdf'])
    run(f'23_ScoreBigwig_{condition}',['ScoreBigwig','--signal',corrected,'--regions',peaks,
        '--output',folder/f'{condition}_footprints.bw','--cores',cores],[folder/f'{condition}_footprints.bw'])
    fcntl.flock(condition_lock,fcntl.LOCK_UN)
    condition_lock.close()
if only_condition:sys.exit(0)
# BINDetect writes many per-motif files; compute locally in WSL, then preserve all outputs in the project.
marker=OUT/'23_BINDetect.complete.json'
required=['bindetect_results.txt','bindetect_figures.pdf']
if not marker.exists():
    start=time.time();cached_signals=[]
    for condition in ['LOY','CTRL']:
        src=OUT/f'{condition}/{condition}_footprints.bw';dst=genome.parent/src.name
        if not dst.exists() or dst.stat().st_size!=src.stat().st_size:shutil.copyfile(src,dst)
        cached_signals.append(dst)
    cached_output=genome.parent/'BINDetect'
    cmd=[str(BIN),*map(str,['BINDetect','--signals',*cached_signals,
        '--motifs',motifs,'--genome',genome,'--peaks',peaks,'--cond-names','LOY','CTRL',
        '--motif-pvalue','0.00005','--naming','id','--skip-excel','--outdir',cached_output,'--cores','8'])]
    (OUT/'23_BINDetect.command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    if not (cached_output/'NATIVE_COMPLETE.json').exists():
        print('Running 23_BINDetect on local cache',flush=True)
        with (ROOT/'logs/23_BINDetect.log').open('w') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        assert all((cached_output/f).stat().st_size>0 for f in required)
        (cached_output/'NATIVE_COMPLETE.json').write_text(json.dumps(dict(command=cmd),indent=2)+'\n')
    print('Preserving complete BINDetect outputs in project',flush=True)
    shutil.copytree(cached_output,OUT/'BINDetect',dirs_exist_ok=True)
    assert all((OUT/'BINDetect'/f).stat().st_size>0 for f in required)
    marker.write_text(json.dumps(dict(elapsed_seconds=time.time()-start,command=cmd,all_outputs_preserved=True),indent=2)+'\n')
else:
    assert all((OUT/'BINDetect'/f).stat().st_size>0 for f in required)
(OUT/'ANALYSIS_COMPLETE.json').write_text(json.dumps(dict(tool='TOBIAS',pooled_conditions=['LOY','CTRL'],
    sample_level_inference=False,predicted_occupancy=True),indent=2)+'\n')
print('TOBIAS stages complete',flush=True)
