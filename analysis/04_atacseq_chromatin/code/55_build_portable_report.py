"""Export the complete viewer as local HTML/JS plus lossless 25 bp coverage tiles.

No BAMs, fragment-level records, Python runtime, or server are shipped.
Coverage tiles store unsigned 32-bit integer fragment-covered bases per 25 bp bin.
The browser sums these integers and applies the original non-Y CPM denominator.
"""
from pathlib import Path
import base64
import gzip
import hashlib
import importlib.util
import json
import math
import shutil
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'exports/HT1080_ATAC_RNA_portable'
AUDIT = ROOT / 'integration/portable_atac_rna'
TILE = 1_000_000


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False), encoding='utf-8')


def packet(path, key, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    compressed = gzip.compress(payload, compresslevel=6, mtime=0)
    text = 'window.receivePortableData(' + json.dumps(key) + ',' + json.dumps(base64.b64encode(compressed).decode()) + ');'
    path.write_text(text, encoding='ascii')


def main():
    started = time.time()
    AUDIT.mkdir(parents=True, exist_ok=True)
    assets = OUT / 'report/interactive_assets'
    assets.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'report/interactive_assets/data.js'
    data = json.loads(source.read_text(encoding='utf-8').split('=', 1)[1].rstrip(';\n'))
    spec = importlib.util.spec_from_file_location('viewer_server', ROOT / 'code/51_interactive_server.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    db = module.LocusData()
    eligible = {(row[0], name.split(':')[1]) for name, panel in data['panels'].items() for row in panel['rows']}
    required = {}
    for gid, feature in sorted(eligible):
        g = db.genes.loc[gid]
        p = db.peaks.loc[db.links[(gid, feature)]]
        assert (p.chrom == g.chrom).all()
        left = max(0, (min(int(g.start)-1, int(p.start.min()))-10000)//25*25)
        right = min(db.sizes[g.chrom], math.ceil((max(int(g.end), int(p.end.max()))+10000)/25)*25)
        required.setdefault(g.chrom, set()).update(range(left//TILE, (right-1)//TILE+1))
    print(f'Exporting {len(data["genes"]):,} genes, {len(eligible):,} gene/feature combinations, {sum(map(len, required.values())):,} coverage tiles.', flush=True)
    coverage_audit = []
    checks = []
    for chrom in sorted(required, key=lambda c: list(db.sizes).index(c)):
        size = db.sizes[chrom]
        n = math.ceil(size/25)
        matrix = np.empty((6,n), dtype='<u4')
        for s, sample in enumerate(module.SAMPLES):
            a = db.npz[sample][chrom+'_start'].astype(np.int64)
            b = db.npz[sample][chrom+'_end'].astype(np.int64)
            assert np.all(a >= 0) and np.all(b <= size) and np.all(b > a)
            sb, eb = a//25, b//25
            delta = np.bincount(sb, minlength=n+1)-np.bincount(eb, minlength=n+1)
            integral = np.cumsum(delta, dtype=np.int64)[:n]*25
            integral -= np.bincount(sb, weights=a%25, minlength=n+1)[:n].astype(np.int64)
            integral += np.bincount(eb, weights=b%25, minlength=n+1)[:n].astype(np.int64)
            total = int((b-a).sum())
            assert int(integral.sum()) == total
            assert integral.min() >= 0 and integral.max() <= np.iinfo(np.uint32).max
            matrix[s] = integral
            coverage_audit.append({'chrom':chrom,'sample':sample,'fragments':len(a),'fragment_base_integral':total,'tile_integral_matches':True,'max_bin_integral':int(integral.max())})
            del a,b,sb,eb,delta,integral
        for i in sorted(required[chrom]):
            payload = np.ascontiguousarray(matrix[:,i*TILE//25:min((i+1)*TILE//25,n)]).tobytes()
            path = OUT / f'report/portable_data/{chrom}/{i}.js'
            packet(path, f'{chrom}:{i}', payload)
            # Verify every stored tile survives serialization and decompression exactly.
            stored = json.loads(path.read_text().split(',',1)[1][:-2])
            assert gzip.decompress(base64.b64decode(stored)) == payload
        del matrix
        chromosome_genes = db.genes[db.genes.chrom == chrom]
        exon_groups = {gid:frame[['start0','end']].values.tolist() for gid,frame in db.exons[db.exons.chrom==chrom].groupby('gene_id', sort=False)}
        models = [{'id':gid,'symbol':r.symbol,'start':int(r.start)-1,'end':int(r.end),'strand':r.strand,'exons':exon_groups.get(gid,[])} for gid,r in chromosome_genes.iterrows()]
        peaks = [{'id':pid,'start':int(p.start),'end':int(p.end),'logFC':float(p.logFC),'FDR':float(p.FDR),
                  'raw':[int(x) for x in db.counts.loc[pid,module.SAMPLES]],'cpm':db.cpm.loc[pid,module.SAMPLES].round(7).tolist()} for pid,p in db.peaks[db.peaks.chrom==chrom].iterrows()]
        links = {gid+':'+f:db.links[(gid,f)] for gid,f in sorted(eligible) if data['genes'][gid]['chrom']==chrom}
        packet(OUT/f'report/portable_data/{chrom}.js', 'meta:'+chrom,
               json.dumps({'size':size,'genes':models,'peaks':peaks,'links':links},separators=(',',':'),allow_nan=False).encode())
        # Independent direct-fragment reference fixtures, one close-up and one large overview per chromosome.
        options = [(gid,f) for gid,f in sorted(eligible) if data['genes'][gid]['chrom']==chrom]
        chosen = max(options, key=lambda t: data['genes'][t[0]]['end']-data['genes'][t[0]]['start'])
        for view in ('best','overview'):
            checks.append({'gene':chosen[0],'feature':chosen[1],'view':view,'expected':db.locus(*chosen,view)})
        db.cache.clear()
        print(f'{chrom}: {len(required[chrom])} tiles, {len(links)} gene/feature links verified ({time.time()-started:.0f}s).', flush=True)
    for filename in ('plotly.min.js','data.js','viewer.js'):
        shutil.copy2(ROOT/'report/interactive_assets'/filename, assets/filename)
    shutil.copy2(ROOT/'code/portable_coverage.js', assets/'portable_coverage.js')
    settings = {'samples':module.SAMPLES,'denominators':[db.denominator[s] for s in module.SAMPLES],'tile_bp':TILE,'source_bin_bp':25}
    (assets/'portable_settings.js').write_text('window.ATAC_PORTABLE_SETTINGS='+json.dumps(settings)+';',encoding='utf-8')
    html = (ROOT/'report/interactive_ATAC_RNA.html').read_text(encoding='utf-8')
    html = html.replace('<script src="interactive_assets/viewer.js"></script>', '<script src="interactive_assets/portable_settings.js"></script><script src="interactive_assets/portable_coverage.js"></script><script src="interactive_assets/viewer.js"></script>')
    html = html.replace('The viewer reads the same accepted ATAC fragment arrays used to generate the existing 25 bp BigWigs.', 'This portable viewer uses losslessly precomputed 25 bp coverage from the same accepted ATAC fragments used for the original BigWigs. All track data load from the extracted package; no server is needed.')
    html = html.replace('<a href="batch_RNA_overlaps.html">Separate-batch report</a> · <a href="ATACseq_report.html">Main ATAC report</a> · ', '')
    html = html.replace('This viewer runs locally; no research data are sent to an external service.', 'Portable offline edition · Keep the extracted folders together. No installation, server or internet connection is needed. No research data are sent to an external service.')
    (OUT/'report/interactive_ATAC_RNA.html').write_text(html,encoding='utf-8')
    (OUT/'START_HERE.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=report/interactive_ATAC_RNA.html"><title>Open HT1080 ATAC/RNA explorer</title><body><h1>HT1080 ATAC/RNA explorer</h1><p><a href="report/interactive_ATAC_RNA.html">Open the interactive report</a></p><p>Extract the entire ZIP before opening. Keep the report and integration folders together.</p></body></html>',encoding='utf-8')
    # Include existing publication-quality locus diagrams, gene pages, and their linked measurement tables.
    shutil.copy2(ROOT/'report/gene_atlas.html',OUT/'report/gene_atlas.html')
    for subdir in ('figures','genes','data'):
        shutil.copytree(ROOT/f'integration/gene_atlas/{subdir}',OUT/f'integration/gene_atlas/{subdir}',dirs_exist_ok=True)
    for filename in ('RNA_normalization_verification.csv','RNA_size_factors.csv','RNA_sample_metadata.csv','RNA_expression_METHODS.txt'):
        dest = OUT/'integration/interactive_atac_rna'/filename
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/'integration/interactive_atac_rna'/filename,dest)
    # Remove main-report links from the copied static atlas, which is a standalone part of this package.
    atlas = OUT/'report/gene_atlas.html'
    atlas.write_text(atlas.read_text(encoding='utf-8').replace('href="ATACseq_report.html"','href="interactive_ATAC_RNA.html"').replace('Back to ATAC report','Back to interactive report'),encoding='utf-8')
    readme = '''HT1080 ATAC / RNA INTERACTIVE REPORT — PORTABLE EDITION

HOW TO OPEN
1. Download the ZIP to your computer.
2. Extract the entire ZIP (Windows: right-click > Extract All; macOS: double-click).
3. Open START_HERE.html in the extracted folder using Chrome, Edge, Firefox or Safari.
   Chrome was used for the offline browser verification. A current browser with
   DecompressionStream support is required; older browsers may not load tracks.
4. Keep report/ and integration/ together with START_HERE.html.

No Python, R, command window, server, BAM files or internet connection is required.
Do not open the HTML inside the ZIP or a cloud-storage preview. Download and extract first.

CONTENTS
All 14,358 selectable genes; combined, August 2025 and May 2026 RNA results;
promoter, gene-body and intergenic ATAC links; selectable FDR thresholds;
per-sample RNA counts and six ATAC coverage tracks; gene models; peak measurements;
full-gene views and individual-peak close-ups; 251 prepared gene diagram PDF books.
Scatter and track plots can be downloaded as PNG through the Plotly toolbar.

METHODS
Statistics and scatter coordinates are copied unchanged from the existing analysis.
Defaults: ATAC FDR < 0.10 and RNA FDR < 0.05, with no fold-change filter.
ATAC primary edgeR QL design: ~ condition. RNA: combined ~ batch + treatment;
separate batch results ~ treatment. RNA count plots share the combined-model
DESeq2 size factors and are not batch-corrected.
Promoter ATAC coordinates are medians over tested linked peaks, whereas ATAC
significance uses any linked peak; this is not a gene-level FDR.
Coverage stores exact integer fragment-covered bases in 25 bp bins, compressed
losslessly. Larger display bins sum these values and divide by their actual
width, then normalize per million original accepted non-Y fragments. All six
replicates and their group means share the same displayed coverage scale.
Gene models merge exons across isoforms. Peak-count plots use primary-model TMM
CPM, a different denominator from the coverage tracks. No CN adjustment or model
refit was performed. Intergenic nearest-TSS links are exploratory.

The package contains normalized coverage and research results; it does not contain
BAMs or individual fragment records. It does not transmit data to a remote service.
Build and verification details are in verification/.
'''
    (OUT/'README.txt').write_text(readme,encoding='utf-8')
    write_json(AUDIT/'coverage_integral_verification.json',coverage_audit)
    write_json(AUDIT/'reference_loci.json',checks)
    result = {'status':'built_pending_browser_verification','genes':len(data['genes']),'panels':len(data['panels']),
              'gene_feature_pairs':len(eligible),'coverage_tiles':sum(map(len,required.values())),'chromosomes':len(required),
              'coverage_encoding':'gzip/base64 little-endian uint32 fragment-covered bases per 25 bp bin, sample-major',
              'source_files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,ROOT/'code/51_interactive_server.py',Path(__file__),ROOT/'code/portable_coverage.js']},
              'seconds':round(time.time()-started,1)}
    write_json(AUDIT/'build_manifest.json',result)
    print(json.dumps(result,indent=2),flush=True)


if __name__ == '__main__':
    main()
