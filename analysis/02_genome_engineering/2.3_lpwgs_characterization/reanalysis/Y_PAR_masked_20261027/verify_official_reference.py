import gzip,hashlib,json,sys
from pathlib import Path
import pysam
official=Path(sys.argv[1]);provenance=Path(sys.argv[2]);out=Path(sys.argv[3])
known='bd90ddb80c86af7fcfeefe7a0909b175'
h=hashlib.md5()
with official.open('rb') as f:
    while chunk:=f.read(8*1024*1024):h.update(chunk)
assert h.hexdigest()==known,('Official reference compressed MD5 differs',h.hexdigest())
expected=json.loads(provenance.read_text())['contigs']
actual={}
with pysam.FastxFile(str(official)) as f:
    for record in f:
        seq=record.sequence
        actual[record.name]={'length':len(seq),'md5':hashlib.md5(seq.upper().encode('ascii')).hexdigest()}
assert set(actual)==set(expected)
for name,row in actual.items():
    assert row['length']==expected[name]['length'],name
    assert row['md5']==expected[name]['masked_sequence_md5'],('Reference sequence mismatch',name)
out.write_text(json.dumps({'official_gzip_md5':known,'all_contig_sequences_match_prepared_reference':True,'contigs':actual},indent=2))
print('Prepared reference matches every contig of the official masked-Y analysis set.',flush=True)
