from pathlib import Path
import shutil,json,hashlib,sys
import pandas as pd
import pysam
ROOT=Path(__file__).resolve().parents[1]
verification_file=ROOT/'config/TF_extension/reference_verified.json'
previous=json.loads(verification_file.read_text()) if verification_file.exists() else None
source=ROOT/'config/TF_extension/hg38_primary.fa'
assert Path(str(source)+'.LF_verified').exists()
cache=Path('/tmp/atacseq-sept08-2026-tf-tools/hg38_primary.fa')
print('Caching normalized FASTA in Linux storage',flush=True)
shutil.copyfile(source,cache)
pysam.faidx(str(cache))
expected=pd.read_csv(ROOT/'config/TF_extension/tested_peak_sequences.tsv',sep='\t')
with pysam.FastaFile(str(cache)) as fa:
    for row in expected.itertuples(index=False):
        chrom,coords=row.peak_id.split(':');start,end=map(int,coords.split('-'))
        actual=fa.fetch(chrom,start,end).upper()
        assert actual==row.sequence, f'Reference mismatch at {row.peak_id}'
    lengths=dict(zip(fa.references,fa.lengths))
with pysam.AlignmentFile(str(ROOT/'bam/CTRL_11.bam'),'rb') as bam:
    assert all(bam.get_reference_length(c)==n for c,n in lengths.items())
h=hashlib.sha256()
with cache.open('rb') as f:
    while block:=f.read(8*1024*1024):
        assert b'\r' not in block
        h.update(block)
info=dict(source=str(source),linux_cache=str(cache),sha256=h.hexdigest(),verified_peak_sequences=len(expected),
          BAM_chromosome_lengths_agree=True,line_endings='LF',sequence_reference='Independent BSgenome getSeq export')
if previous is not None:
    assert previous['sha256']==info['sha256'], 'Reference changed; inspect before reusing analysis outputs'
verification_file.write_text(json.dumps(info,indent=2)+'\n')
if previous is not None:
    print(json.dumps(info,indent=2),flush=True);sys.exit(0)
# Retain failed/interrupted outputs without presenting them as completed analyses.
for arm in ['primary_TMM','passage_TMM']:
    for direction in ['gained','lost']:
        folder=ROOT/f'motifs/TF_extension/HOMER/{arm}/{direction}'
        if folder.exists():
            dest=folder.with_name(direction+'_invalid_reference_or_configuration')
            assert not dest.exists();folder.rename(dest)
            log=ROOT/f'logs/25_HOMER_{arm}_{direction}.log'
            if log.exists():shutil.copy2(log,dest/'invalid_run.log')
folder=ROOT/'motifs/TF_extension/TOBIAS/CTRL'
if folder.exists():
    dest=folder.with_name('CTRL_interrupted_reference_repair');assert not dest.exists();folder.rename(dest)
    log=ROOT/'logs/23_ATACorrect_CTRL.log'
    if log.exists():shutil.copy2(log,dest/'interrupted_run.log')
print(json.dumps(info,indent=2),flush=True)
