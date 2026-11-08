"""Mask only official chrY PAR intervals; retain every other reference base."""
import argparse,gzip,hashlib,json,os
from pathlib import Path
import pysam

def prepare(source,bed,bam,dest):
    intervals={}
    for line in Path(bed).read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        fields=line.split();c=fields[0];s,e=map(int,fields[1:3])
        if c=='chrY':intervals.setdefault(c,[]).append((s,e))
    assert sorted(intervals.get('chrY',[]))==[(0,2458320),(62122809,62460029)],intervals
    with pysam.AlignmentFile(bam,'rb') as b:expected=dict(zip(b.references,b.lengths))
    source=Path(source);dest=Path(dest)
    assert not dest.exists(),'Refusing to overwrite completed reference'
    tmp=dest.with_suffix(dest.suffix+'.partial')
    digests={};name=None;pos=0;changed=0;yin=0
    md5src=hashlib.md5();sha=hashlib.sha256()
    with source.open('rb') as handle:
        while chunk:=handle.read(8*1024*1024):md5src.update(chunk)
    with gzip.open(source,'rt',encoding='ascii') as inp,tmp.open('wt',encoding='ascii',newline='\n') as output:
        for line in inp:
            if line.startswith('>'):
                if name is not None:
                    assert pos==expected[name],(name,pos,expected[name])
                    digests[name]={'length':pos,'original_sequence_md5':original.hexdigest(),'masked_sequence_md5':masked.hexdigest(),'changed_bases':changed}
                name=line[1:].split()[0];pos=0;changed=0
                assert name in expected,('unexpected reference contig',name)
                original=hashlib.md5();masked=hashlib.md5()
                output.write(line.rstrip()+'\n');continue
            seq=line.strip();original.update(seq.upper().encode('ascii'))
            new=seq
            for start,end in intervals.get(name,[]):
                a=max(0,start-pos);b=min(len(seq),end-pos)
                if b>a:
                    changed+=sum(x.upper()!='N' for x in seq[a:b]);yin+=b-a
                    new=new[:a]+'N'*(b-a)+new[b:]
            masked.update(new.upper().encode('ascii'))
            output.write(new+'\n');pos+=len(seq)
        assert name is not None and pos==expected[name]
        digests[name]={'length':pos,'original_sequence_md5':original.hexdigest(),'masked_sequence_md5':masked.hexdigest(),'changed_bases':changed}
    assert set(digests)==set(expected),('contig mismatch',set(expected)-set(digests))
    assert yin==2458320+(62460029-62122809),yin
    assert all(v['changed_bases']==0 for k,v in digests.items() if k!='chrY')
    assert all(v['original_sequence_md5']==v['masked_sequence_md5'] for k,v in digests.items() if k!='chrY')
    os.replace(tmp,dest)
    result={'source':str(source),'source_compressed_md5':md5src.hexdigest(),'masked_reference':str(dest),
        'mask_bed':str(bed),'mask_bed_sha256':hashlib.sha256(Path(bed).read_bytes()).hexdigest(),
        'masked_intervals_0based':intervals,'masked_span_bases':yin,'contigs':digests,
        'BAM_dictionary_names_lengths_match':True,'all_non_Y_sequences_unchanged':True}
    dest.with_suffix('.provenance.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='contigs'},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('bed');p.add_argument('bam');p.add_argument('dest')
    a=p.parse_args();prepare(a.source,a.bed,a.bam,a.dest)
