"""Assemble a student handoff: working reports, processed results, code and documentation.

Excludes alignments, per-fragment records, signal tracks, the reference genome, native
per-motif TOBIAS site files, software installations, logs and browser test profiles.

The interactive explorer and gene atlas are taken unchanged from the verified portable
package (exports/HT1080_ATAC_RNA_portable); the live versions need a local Python server
reading counts/*_fragments.npz and would not work for the recipient.

Usage: python code/build_student_handoff.py [--zip]
Refuses to overwrite an existing handoff folder.
"""
import argparse
import csv
import hashlib
import html
import json
import re
import shutil
import sys
import zipfile
from datetime import date
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
PORTABLE = ROOT / 'exports/HT1080_ATAC_RNA_portable'
NAME = 'HT1080_ATAC_student_handoff'
OUT = ROOT / 'exports' / NAME

# Never shipped, wherever they occur.
RAW_SUFFIXES = {'.bam', '.bai', '.cram', '.crai', '.bw', '.bigwig', '.bedgraph', '.bdg', '.fa', '.fasta',
                '.2bit', '.fastq', '.fq', '.npz', '.npy', '.pickle', '.pid', '.log'}
SKIP_DIR_NAMES = {'__pycache__', 'browser_test_profile', 'extracted_check', 'node_modules', '.git'}
SKIP_DIR_PATTERNS = (re.compile(r'_invalid_'), re.compile(r'_interrupted_'))
MAX_CONFIG_BYTES = 10_000_000

# Portable versions of these replace the live server-backed ones.
PORTABLE_WINS = {'report/interactive_ATAC_RNA.html', 'report/gene_atlas.html'}
PORTABLE_WINS_PREFIX = ('report/interactive_assets/', 'report/portable_data/')

REPORTS = [
    ('report/ATACseq_report.html', 'Main ATAC-seq report',
     'QC, peaks, differential accessibility, annotation, motifs/TF extension and RNA integration. Start here.'),
    ('report/interactive_ATAC_RNA.html', 'Interactive ATAC / RNA explorer',
     'Click any gene: RNA counts by sample, six ATAC coverage tracks, linked peaks. Portable edition; no server needed.'),
    ('report/gene_atlas.html', 'Gene-locus atlas',
     '251 prepared gene pages with locus diagrams and PDF books.'),
    ('report/batch_RNA_overlaps.html', 'Separate RNA-seq batch overlaps',
     'ATAC overlap with the August 2025 and May 2026 RNA batches analysed separately.'),
    ('report/down70_enrichment.html', 'Exact 70-gene enrichment',
     'GO and chromosome clustering for the 70 genes down in both assays (ATAC FDR 0.10, RNA FDR 0.05).'),
    ('report/pathway_GSEA.html', 'Pathway enrichment (GSEA)',
     'Hallmark and Reactome fgsea on ATAC and RNA ranks, with concordance.'),
    ('report/pathway_gene_drivers.html', 'E2F, MYC and interferon gene drivers',
     'Gene-level inspection of the leading pathways.'),
    ('report/ribosomal_GSEA_sensitivity.html', 'Ribosomal exclusion sensitivity',
     'GSEA repeated without cytosolic / mitochondrial ribosomal genes.'),
    ('report/RNA_interaction_CN.html', 'RNA interaction and copy-number context',
     'Batch-by-treatment interaction and lpWGS dosage context.'),
    ('report/TFAP2C_review.html', 'TFAP2C evidence review',
     'TFAP2C across HOMER, chromVAR, TOBIAS and RNA.'),
]

ATTR = re.compile(r"""(?:href|src|poster)\s*=\s*(?:"([^"]*)"|'([^']*)')""", re.I)
CSS_URL = re.compile(r"""url\(\s*['"]?([^'")]+)['"]?\s*\)""", re.I)


def rel(p, base=ROOT):
    return p.relative_to(base).as_posix()


def skipped_dir(parts):
    return any(x in SKIP_DIR_NAMES or any(r.search(x) for r in SKIP_DIR_PATTERNS) for x in parts)


def allowed(path, base=ROOT):
    parts = path.relative_to(base).parts
    if skipped_dir(parts[:-1]):
        return False
    name = path.name.lower()
    if any(name.endswith(s) for s in RAW_SUFFIXES) or name.endswith(('.fastq.gz', '.fq.gz', '.fa.gz')):
        return False
    return True


def local_refs(page):
    """Yield (raw reference, resolved path) for local links in an HTML/CSS file."""
    text = page.read_text(encoding='utf-8', errors='replace')
    found = {m.group(1) or m.group(2) or '' for m in ATTR.finditer(text)}
    if page.suffix.lower() in ('.html', '.htm', '.css'):
        found |= {m.group(1) for m in CSS_URL.finditer(text)}
    for ref in found:
        ref = ref.strip()
        parts = urlsplit(ref)
        if not ref or ref.startswith('#') or parts.scheme or parts.netloc:
            continue
        p = unquote(parts.path)
        # attribute-like fragments inside minified JS or template strings are not links
        if not p or p.startswith('/') or any(c in p for c in "'+(){}<>;,=\n") or ' ' in p.strip() and '.' not in p:
            continue
        yield ref, (page.parent / p).resolve()


def crawl(start_pages, base):
    """Follow local links from HTML pages. Returns (files, directory links, missing)."""
    queue, seen, dirs, missing = list(start_pages), set(), set(), []
    while queue:
        f = queue.pop()
        if f in seen:
            continue
        seen.add(f)
        if f.suffix.lower() not in ('.html', '.htm', '.css'):
            continue
        for ref, target in local_refs(f):
            try:
                target.relative_to(base)
            except ValueError:
                missing.append((rel(f, base), ref, 'outside package'))
                continue
            if target.is_file():
                queue.append(target)
            elif target.is_dir():
                dirs.add(target)
            else:
                missing.append((rel(f, base), ref, 'not found'))
    return seen, dirs, missing


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


class Builder:
    def __init__(self):
        self.copied = {}       # destination relpath -> source description
        self.conflicts = []
        self.excluded = []     # (relpath, reason) for files deliberately left out of included trees

    def add(self, src, dest_rel, origin):
        dest = OUT / dest_rel
        if dest_rel in self.copied:
            if dest.stat().st_size != src.stat().st_size or sha256(dest) != sha256(src):
                self.conflicts.append((dest_rel, self.copied[dest_rel], origin))
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        self.copied[dest_rel] = origin

    def add_tree(self, folder, origin, recursive=True, max_bytes=None):
        folder = ROOT / folder
        if not folder.exists():
            raise FileNotFoundError(folder)
        for src in sorted(folder.rglob('*') if recursive else folder.glob('*')):
            if not src.is_file():
                continue
            r = rel(src)
            if not allowed(src):
                if not skipped_dir(src.relative_to(ROOT).parts[:-1]):
                    self.excluded.append((r, 'raw/intermediate file type'))
                continue
            if max_bytes and src.stat().st_size > max_bytes:
                self.excluded.append((r, f'larger than {max_bytes // 1_000_000} MB'))
                continue
            if r in PORTABLE_WINS or r.startswith(PORTABLE_WINS_PREFIX):
                continue
            self.add(src, r, origin)


def build_index():
    rows = '\n'.join(
        f'<li><a href="{html.escape(path)}">{html.escape(title)}</a><span>{html.escape(note)}</span></li>'
        for path, title, note in REPORTS)
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HT1080 LOY ATAC-seq handoff</title>
<style>:root{{--blue:#0072b2;--ink:#203446;--muted:#637583;--line:#dbe4eb}}body{{margin:0;background:#f3f6f8;color:var(--ink);font:15px/1.5 system-ui,Arial,sans-serif}}main{{max-width:860px;margin:0 auto;padding:32px 16px 48px}}h1{{font-size:26px;margin:0 0 4px;letter-spacing:-.5px}}h2{{font-size:17px;margin:28px 0 8px}}p{{margin:6px 0}}.muted{{color:var(--muted)}}a{{color:var(--blue)}}ul.reports{{list-style:none;padding:0;margin:0;background:#fff;border:1px solid var(--line);border-radius:9px}}ul.reports li{{padding:12px 16px;border-bottom:1px solid var(--line)}}ul.reports li:last-child{{border-bottom:0}}ul.reports a{{font-weight:600;display:block}}ul.reports span{{color:var(--muted);font-size:13px}}code{{background:#e8eef2;padding:1px 5px;border-radius:4px;font-size:13px}}</style></head><body><main>
<h1>HT1080 LOY vs CTRL · ATAC-seq reprocessing</h1>
<p class="muted">Handoff package built {date.today().isoformat()}. Three independent CTRL and three independent LOY replicates; positive log fold change means gain in LOY.</p>
<h2>Reports</h2>
<ul class="reports">
{rows}
</ul>
<h2>Using this folder</h2>
<p>Extract the whole ZIP first, then open this page in Chrome, Edge, Firefox or Safari. Keep the folder structure intact: the reports load figures, tables and coverage data through relative paths. No internet, Python, R or local server is needed to read the reports.</p>
<p>Processed tables are beside the reports: <code>differential/</code>, <code>annotation/</code>, <code>counts/</code>, <code>peaks/</code>, <code>qc/</code>, <code>motifs/</code>, <code>integration/</code>. Scripts are in <code>code/</code>. See <a href="HANDOFF_README.txt">HANDOFF_README.txt</a> for what is and is not included, and <a href="README.md">README.md</a> for the full analysis description.</p>
</main></body></html>
'''


def build_readme(n_files, n_bytes, excluded_summary):
    return f'''HT1080 LOY vs CTRL ATAC-seq REPROCESSING — STUDENT HANDOFF
Built {date.today().isoformat()} by code/build_student_handoff.py. {n_files} files, {n_bytes / 1e9:.2f} GB.

HOW TO OPEN
1. Extract the entire ZIP to a local folder (not a cloud-storage preview).
2. Open START_HERE.html in Chrome, Edge, Firefox or Safari. It links all ten reports.
3. Keep every folder where it is. The reports find their figures, tables and coverage data
   through relative paths. No internet connection, Python, R or server is required.
4. On Windows, extract to a short location such as C:\\Data\\ or your Documents folder. Paths inside
   the package are up to about 110 characters; very deep destinations exceed the Windows path limit
   and extraction stops partway.

WHAT IS INCLUDED
report/         Ten HTML reports and their scripts/data. The interactive explorer and gene atlas are
                the verified portable edition (see PORTABLE_REPORT_README.txt and verification/).
integration/    RNA integration, overlaps, scatters, gene atlas, pathway GSEA, RNA x copy-number context.
differential/   edgeR / csaw results for every tested region, plus fitted model objects (models.rds).
annotation/     Peak annotation, promoter links, GO tables.
counts/         Peak count matrix, logCPM, background bins, csaw windows.
peaks/          Per-sample MACS3 narrowPeak/summit files, consensus membership, matched-depth recalls.
qc/             Per-sample preparation metrics, TSS profiles, QC figures.
motifs/         Motif enrichment, chromVAR, HOMER, concordance tables, TFAP2C review, and TOBIAS
                summary-level outputs (BINDetect results table, figures, HTML).
code/ config/   All analysis scripts, designs, parameters, annotations and environment locks.
README.md, PROCESSING_PLAN.md, RUN_STATUS.md, TF_EXTENSION_PLAN.md    Original project documentation.
MANIFEST.tsv    Every file with size and SHA-256.  EXCLUDED.tsv lists what was deliberately left out.

WHAT IS NOT INCLUDED
- Alignments: source BAMs (bam/) and the depth-matched / pooled BAMs under motifs/TF_extension/TOBIAS/.
- Per-fragment records (counts/*_fragments.npz) and bigWig signal tracks (tracks/, TOBIAS *.bw).
- The hg38 reference FASTA and tested peak sequences (config/TF_extension/).
- TOBIAS per-motif site folders ({excluded_summary['tobias_dirs']} folders, about 3 GB).
- HOMER/TOBIAS runs labelled invalid or interrupted in the project (excluded from results there too).
- Software installations (tools/), logs/, browser test profiles, and the earlier export build area.

CONSEQUENCES
- All reports, tables and figures work as delivered. The ten reports, the index and a gene page were
  opened from an extracted copy of the ZIP in headless Chrome with networking disabled: no script
  errors, failed loads or broken images, and the explorer loaded coverage tracks for a test gene.
- Two native third-party pages are included as produced and are not part of the ten reports:
  motifs/TF_extension/TOBIAS/BINDetect/bindetect_LOY_CTRL.html loads Highcharts from the internet and
  is blank offline; HOMER knownResults.html pages link to de novo / GO pages that these known-motif
  runs never wrote. Use bindetect_results.txt / .xlsx and the concordance tables instead.
- The live server version of the explorer (Start_interactive_report.cmd, code/51_interactive_server.py)
  is not runnable from this package; use report/interactive_ATAC_RNA.html instead.
- Scripts from BAM preparation through counting (code/02-04, 10, 22-23, 29) and any stage needing the
  reference genome cannot be rerun without the excluded inputs. Stages that start from counts/,
  differential/ and annotation/ (for example 05-07, 13-17, 57-69) have their inputs here, but paths
  to the external RNA-seq and lpWGS sources recorded in the scripts point to the original workstation.
- config/samples.tsv refers to BAM locations that are not part of this package.

Ask for the excluded data if you need to rerun upstream stages.
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--zip', action='store_true', help='also write exports/%s.zip and its SHA-256' % NAME)
    args = parser.parse_args()
    if OUT.exists():
        sys.exit(f'{OUT} already exists. Inspect and remove it before rebuilding.')
    if not PORTABLE.is_dir():
        sys.exit(f'Missing verified portable package: {PORTABLE}')
    b = Builder()

    # 1. Verified portable explorer + atlas, byte-for-byte.
    for src in sorted(PORTABLE.rglob('*')):
        if not src.is_file():
            continue
        r = rel(src, PORTABLE)
        if r == 'START_HERE.html':
            continue
        if r == 'README.txt':
            r = 'PORTABLE_REPORT_README.txt'
        b.add(src, r, 'portable package')

    # 2. Everything the live report pages link to.
    live_pages = [ROOT / p for p, _, _ in REPORTS]
    files, dir_links, missing = crawl(live_pages, ROOT)
    if missing:
        print('Unresolved links in live reports (informational):', *missing[:20], sep='\n  ')
    for src in sorted(files):
        r = rel(src)
        if r in PORTABLE_WINS or r.startswith(PORTABLE_WINS_PREFIX):
            continue
        if not allowed(src):
            sys.exit(f'A report links to an excluded file type: {r}')
        b.add(src, r, 'linked from report')
    for d in sorted(dir_links):
        b.add_tree(rel(d), 'directory linked from report')
    for extra in ('report/summary.json', 'report/verification.txt'):
        b.add(ROOT / extra, extra, 'report record')

    # 3. Processed results, code and documentation.
    for doc in ('README.md', 'PROCESSING_PLAN.md', 'RUN_STATUS.md', 'TF_EXTENSION_PLAN.md'):
        b.add(ROOT / doc, doc, 'documentation')
    for folder in ('code', 'qc', 'annotation', 'differential', 'counts', 'peaks'):
        b.add_tree(folder, 'processed results' if folder != 'code' else 'code')
    b.add_tree('config', 'configuration', max_bytes=MAX_CONFIG_BYTES)
    b.add_tree('motifs', 'processed results', recursive=False)
    b.add_tree('motifs/TF_extension', 'processed results', recursive=False)
    for folder in ('chromVAR', 'concordance', 'HOMER', 'TFAP2C_review'):
        b.add_tree(f'motifs/TF_extension/{folder}', 'processed results')
    b.add_tree('motifs/TF_extension/TOBIAS', 'TOBIAS summary', recursive=False)
    b.add_tree('motifs/TF_extension/TOBIAS/BINDetect', 'TOBIAS summary', recursive=False)
    for cond in ('CTRL', 'LOY'):
        b.add_tree(f'motifs/TF_extension/TOBIAS/{cond}', 'TOBIAS summary', recursive=False)
    tobias_dirs = sum(1 for d in (ROOT / 'motifs/TF_extension/TOBIAS/BINDetect').iterdir() if d.is_dir())
    b.excluded.append(('motifs/TF_extension/TOBIAS/BINDetect/*/', f'{tobias_dirs} per-motif site folders'))
    for folder in ('RNAseq', 'FDR_sensitivity', 'feature_overlaps', 'feature_scatters', 'promoter_FDR010_review',
                   'batch_RNA_overlaps', 'down70_ATAC010_RNA005', 'pathway_GSEA', 'RNA_interaction_CN'):
        b.add_tree(f'integration/{folder}', 'processed results')
    b.add_tree('integration/interactive_atac_rna', 'processed results', recursive=False)
    for whole in ('bam/', 'tracks/', 'tools/', 'logs/', 'integration/portable_atac_rna/', 'integration/gene_atlas/coverage/',
                  'exports/HT1080_ATAC_RNA_portable.zip', 'Start_interactive_report.cmd'):
        b.excluded.append((whole, 'not part of handoff'))

    # 4. Front page and handoff notes.
    (OUT / 'START_HERE.html').write_text(build_index(), encoding='utf-8')
    copied_files = [p for p in OUT.rglob('*') if p.is_file()]
    (OUT / 'HANDOFF_README.txt').write_text(
        build_readme(len(copied_files) + 3, sum(p.stat().st_size for p in copied_files), {'tobias_dirs': tobias_dirs}),
        encoding='utf-8')

    # 5. Verify the assembled folder on its own.
    problems, preexisting = [], []
    for path in OUT.rglob('*'):
        if path.is_file() and not allowed(path, OUT):
            problems.append(f'excluded file type present: {rel(path, OUT)}')
    pages = [OUT / 'START_HERE.html'] + [OUT / p for p, _, _ in REPORTS]
    for p in pages:
        if not p.is_file():
            problems.append(f'missing report: {rel(p, OUT)}')
    _, _, missing_out = crawl([p for p in pages if p.is_file()], OUT)
    # also check every HTML in the package, not only pages reachable from the index
    _, _, missing_all = crawl(sorted(OUT.rglob('*.html')), OUT)
    for page, ref, why in sorted(set(missing_out) | set(missing_all)):
        # A link that is equally unresolved in the source project is not a packaging defect
        # (HOMER's stock knownResults.html links to de novo/GO pages that known-motif runs never write).
        source_page = ROOT / page
        if source_page.is_file() and not (source_page.parent / unquote(urlsplit(ref).path)).exists():
            preexisting.append(f'{page} -> {ref}')
        else:
            problems.append(f'broken link: {page} -> {ref} ({why})')
    for dest_rel, first, second in b.conflicts:
        problems.append(f'content differs between {first} and {second}: {dest_rel}')
    portable_changed = [r for r in (rel(s, PORTABLE) for s in PORTABLE.rglob('*') if s.is_file())
                        if r not in ('START_HERE.html', 'README.txt') and sha256(PORTABLE / r) != sha256(OUT / r)]
    problems += [f'portable file altered: {r}' for r in portable_changed]

    with (OUT / 'EXCLUDED.tsv').open('w', newline='', encoding='utf-8') as handle:
        w = csv.writer(handle, delimiter='\t')
        w.writerow(['path', 'reason'])
        w.writerows(sorted(set(b.excluded)))
    with (OUT / 'MANIFEST.tsv').open('w', newline='', encoding='utf-8') as handle:
        w = csv.writer(handle, delimiter='\t')
        w.writerow(['path', 'bytes', 'sha256', 'origin'])
        for p in sorted(x for x in OUT.rglob('*') if x.is_file() and x.name != 'MANIFEST.tsv'):
            r = rel(p, OUT)
            w.writerow([r, p.stat().st_size, sha256(p), b.copied.get(r, 'handoff build')])

    by_top = {}
    for p in OUT.rglob('*'):
        if p.is_file():
            k = p.relative_to(OUT).parts[0]
            n, s = by_top.get(k, (0, 0))
            by_top[k] = (n + 1, s + p.stat().st_size)
    summary = {'built': date.today().isoformat(), 'files': sum(n for n, _ in by_top.values()),
               'bytes': sum(s for _, s in by_top.values()),
               'by_top_level': {k: {'files': n, 'MB': round(s / 1e6, 1)} for k, (n, s) in sorted(by_top.items())},
               'longest_relative_path': max(len(rel(p, OUT)) for p in OUT.rglob('*') if p.is_file()),
               'excluded_entries': len(set(b.excluded)),
               'links_unresolved_in_source_too': preexisting, 'problems': problems}
    print(json.dumps(summary, indent=1))
    if problems:
        sys.exit('Handoff folder assembled with problems; not packaging.')

    if args.zip:
        archive = OUT.parent / f'{NAME}.zip'
        if archive.exists():
            sys.exit(f'{archive} already exists.')
        stored = {'.png', '.jpg', '.jpeg', '.pdf', '.zip', '.gz', '.xlsx', '.tiff', '.rds'}
        with zipfile.ZipFile(archive, 'w', allowZip64=True) as z:
            for p in sorted(x for x in OUT.rglob('*') if x.is_file()):
                method = zipfile.ZIP_STORED if p.suffix.lower() in stored else zipfile.ZIP_DEFLATED
                z.write(p, f'{NAME}/{rel(p, OUT)}', compress_type=method, compresslevel=6)
        with zipfile.ZipFile(archive) as z:
            bad = z.testzip()
            names = len(z.namelist())
        if bad:
            sys.exit(f'Archive failed CRC check at {bad}')
        digest = sha256(archive)
        (archive.parent / f'{archive.name}.sha256').write_text(f'{digest}  {archive.name}\n')
        print(json.dumps({'archive': str(archive), 'GB': round(archive.stat().st_size / 1e9, 2),
                          'entries': names, 'sha256': digest}, indent=1))


if __name__ == '__main__':
    main()
