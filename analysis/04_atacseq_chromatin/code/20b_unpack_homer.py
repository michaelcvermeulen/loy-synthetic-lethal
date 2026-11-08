from pathlib import Path
import tarfile, hashlib, shlex, re
ROOT=Path(__file__).resolve().parents[1]
archive=ROOT/'tools/TF_extension/homer-5.1.tar.bz2'
assert hashlib.sha256(archive.read_bytes()).hexdigest()=='58436010bfebe77abc867204b629e93622df09abe51acbe41fe92815468c0570'
dest=Path('/tmp/atacseq-sept08-2026-tf-tools/homer_package')
dest.mkdir(parents=True,exist_ok=True)
with tarfile.open(archive,'r:bz2') as t:
    t.extractall(dest,filter='data')
for line in (dest/'info/has_prefix').read_text().splitlines():
    old,mode,rel=shlex.split(line)
    path=dest/rel
    data=path.read_bytes()
    if mode=='text': data=data.replace(old.encode(),str(dest).encode())
    elif mode=='binary':
        new=str(dest).encode();assert len(new)<=len(old)
        # Relocate the whole NUL-terminated string; padding inside the path truncates its suffix.
        def relocate(match):
            chunk=match.group(0); replaced=chunk.replace(old.encode(),new)
            return replaced+b'\0'*(len(chunk)-len(replaced))
        data=re.sub(re.escape(old.encode())+b'[^\0]*',relocate,data)
    else:raise ValueError(mode)
    path.write_bytes(data)
print((dest/'info/index.json').read_text())
print('Relocated files:',len((dest/'info/has_prefix').read_text().splitlines()))
print('HOMER paths:',[str(p) for p in dest.rglob('findMotifsGenome.pl')])
