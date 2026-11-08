"""Validate BAM->FASTQ recovery using all primary read sequences and qualities."""
import argparse, collections, gzip, hashlib, json
from pathlib import Path
import pysam

MOD=(1<<128)-1
def fingerprint(name,mate,seq,qual):
    data=(name+'\0'+str(mate)+'\0'+seq+'\0'+qual).encode('ascii')
    return int.from_bytes(hashlib.blake2b(data,digest_size=16).digest(),'little')
def update(d,name,mate,seq,qual):
    value=fingerprint(name,mate,seq,qual)
    d['records']+=1;d['bases']+=len(seq)
    d['hash_xor']^=value;d['hash_sum']=(d['hash_sum']+value)&MOD
def blank():return dict(records=0,bases=0,hash_xor=0,hash_sum=0)
def bam_qc(bam_path,out):
    counts=collections.Counter();fp=blank();lengths=collections.Counter()
    with pysam.AlignmentFile(bam_path,'rb',threads=4) as bam:
        header=bam.header.to_dict()
        for r in bam.fetch(until_eof=True):
            counts['alignment_records']+=1
            if r.is_secondary:counts['secondary']+=1;continue
            if r.is_supplementary:counts['supplementary']+=1;continue
            counts['primary']+=1
            counts['primary_unmapped']+=r.is_unmapped
            counts['primary_properly_paired']+=r.is_proper_pair
            counts['primary_duplicate_flag']+=r.is_duplicate
            counts['primary_qcfail_flag']+=r.is_qcfail
            counts['primary_hardclipped']+=any(op==5 for op,n in (r.cigartuples or []))
            mate=1 if r.is_read1 and not r.is_read2 else 2 if r.is_read2 and not r.is_read1 else 0
            counts[f'mate_{mate}']+=1
            if not r.is_paired:counts['primary_unpaired']+=1
            seq=r.get_forward_sequence();qual=r.qual
            if seq is None or qual is None:
                counts['missing_sequence_or_quality']+=1;continue
            if r.is_reverse:qual=qual[::-1]
            assert len(seq)==len(qual)
            lengths[len(seq)]+=1
            update(fp,r.query_name,mate,seq,qual)
            if counts['primary']%2000000==0:print('Checked primary reads:',counts['primary'],flush=True)
    result={'bam':str(bam_path),'header':header,'counts':dict(counts),'length_histogram':dict(lengths),'read_fingerprint':fp}
    Path(out).write_text(json.dumps(result,indent=2))
    for key in ('primary_hardclipped','missing_sequence_or_quality','primary_unpaired','mate_0','primary_qcfail_flag','primary_duplicate_flag'):
        assert counts[key]==0,(key,counts[key])
    assert counts['mate_1']==counts['mate_2']
    print('BAM recovery prerequisites pass:',counts['primary'],'primary records',flush=True)
def fastq_qc(r1,r2,original,out):
    fp=blank();pairs=0
    with pysam.FastxFile(r1) as f1,pysam.FastxFile(r2) as f2:
        it2=iter(f2)
        for a in f1:
            b=next(it2,None)
            assert b is not None and a.name==b.name,('unmatched FASTQ pair',a.name)
            update(fp,a.name,1,a.sequence,a.quality)
            update(fp,b.name,2,b.sequence,b.quality)
            pairs+=1
            if pairs%1000000==0:print('Checked FASTQ pairs:',pairs,flush=True)
        assert next(it2,None) is None
    expected=json.loads(Path(original).read_text())['read_fingerprint']
    result={'pairs':pairs,'read_fingerprint':fp,'matches_original_primary_reads':fp==expected}
    Path(out).write_text(json.dumps(result,indent=2))
    assert fp==expected,(fp,expected)
    print('FASTQ read names, sequences, qualities and pairing validated:',pairs,'pairs',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='mode',required=True)
    b=sub.add_parser('bam');b.add_argument('bam');b.add_argument('out')
    f=sub.add_parser('fastq');f.add_argument('r1');f.add_argument('r2');f.add_argument('original');f.add_argument('out')
    a=p.parse_args()
    if a.mode=='bam':bam_qc(a.bam,a.out)
    else:fastq_qc(a.r1,a.r2,a.original,a.out)
