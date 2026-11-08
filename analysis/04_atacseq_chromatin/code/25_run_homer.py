from pathlib import Path
import os,subprocess,json,time,re,shutil
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'motifs/TF_extension/HOMER'
BIN=Path('/tmp/atacseq-sept08-2026-tf-tools/homer_package/bin')
env=os.environ.copy();env['PATH']=str(BIN)+os.pathsep+env['PATH']
reference=json.loads((ROOT/'config/TF_extension/reference_verified.json').read_text())
genome=Path(reference['linux_cache']);assert genome.exists()
expected_sequences=dict(line.rstrip('\n').split('\t') for line in (ROOT/'config/TF_extension/tested_peak_sequences.tsv').open() if not line.startswith('peak_id\t'))
for arm in ['primary_TMM','passage_TMM']:
    for direction in ['gained','lost']:
        folder=OUT/arm/direction;folder.mkdir(exist_ok=True)
        marker=folder/'COMPLETE.json'
        if marker.exists():
            header=(folder/'knownResults.txt').read_text().splitlines()[0]
            expected=sum(1 for _ in (OUT/arm/f'{direction}_target.bed').open())
            observed=int(re.search(r'Target Sequences with Motif\(of (\d+)\)',header).group(1))
            if observed==expected:continue
            invalid=folder.with_name(folder.name+'_invalid_configuration')
            assert not invalid.exists()
            folder.rename(invalid);folder.mkdir()
            previous_log=ROOT/f'logs/25_HOMER_{arm}_{direction}.log'
            if previous_log.exists():shutil.copy2(previous_log,invalid/'invalid_run.log')
        args=[str(BIN/'findMotifsGenome.pl'),str(OUT/arm/f'{direction}_target.bed'),
            str(genome),str(folder),
            '-size','given','-bg',str(OUT/arm/f'{direction}_background.bed'),
            '-mknown',str(ROOT/'config/TF_extension/JASPAR2024_CORE_vertebrates.homer'),
            '-nomotif','-useOldBg','-nofacts','-keepFiles','-p','6',
            '-preparsedDir','/tmp/atacseq-sept08-2026-tf-tools/homer_preparsed']
        (folder/'command.json').write_text(json.dumps(args,indent=2)+'\n')
        print('HOMER',arm,direction,flush=True);start=time.time()
        with (ROOT/f'logs/25_HOMER_{arm}_{direction}.log').open('w') as log:
            subprocess.run(args,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        result=folder/'knownResults.txt'
        assert result.exists() and result.stat().st_size>100
        header=result.read_text().splitlines()[0]
        expected=sum(1 for _ in (OUT/arm/f'{direction}_target.bed').open())
        observed=int(re.search(r'Target Sequences with Motif\(of (\d+)\)',header).group(1))
        assert observed==expected, f'HOMER extracted {observed} of {expected} target sequences'
        checked=0
        for line in (folder/'seq.tsv').read_text().splitlines():
            sid,seq=line.split('\t')[:2];pid=sid[2:] if sid.startswith('BG') else sid
            assert seq.upper()==expected_sequences[pid], f'HOMER sequence mismatch for {pid}'
            checked+=1
        expected_bg=sum(1 for _ in (OUT/arm/f'{direction}_background.bed').open())
        assert checked==expected+expected_bg
        marker.write_text(json.dumps(dict(elapsed_seconds=time.time()-start,version='5.1',command=args),indent=2)+'\n')
print('HOMER stages complete',flush=True)
