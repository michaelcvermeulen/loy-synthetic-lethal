"""One-time focused integration of the combined DepMap F01 figure."""
from pathlib import Path
from datetime import datetime
import hashlib
import json
import re
import shutil

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
STEM = 'ch03_main_depmap_y_status_dependencies'
CODE = 'CH3-S3.4.1-F01'
SELECTED = 'CH3_S3_4_1_F01_depmap_y_status_dependencies'
REL = f'figures/chapter3/selected/{SELECTED}/{SELECTED}'
TITLE = 'Y-status classification and associated gene dependencies in DepMap'
ANCHOR = 'main--y-status-classification-and-associated-gene-dependencies-in-depmap'
BOUNDARY = '### 3.4.4 A repeat-targeting CRISPR strategy for Y depletion'
old_stems = ['CH3_S3_4_1_F01_depmap_y_status_calls', 'CH3_S3_4_2_F01_depmap_pan_cancer_dependencies']
doc_names = ['Chapter3_results.md', 'Results_edit.md', 'Chapter3_figure_captions.md', 'Chapter3_source_notes.md']
raw = {name: (ROOT/name).read_bytes() for name in doc_names}
text = {name: value.decode('utf-8').replace('\r\n', '\n') for name, value in raw.items()}

caption_body = '''**a,** Distribution of the Y-expression score, the mean log₂(TPM + 1) of DDX3Y, EIF1AY, KDM5D, RPS4Y1, USP9Y, UTY and ZFY, in 711 female and 905 male default RNA profiles. Histograms use 0.2-unit bins; dashed lines mark score thresholds of 1.0 and 2.5. **b,** Y-expression score against the median DNA copy ratio of the same seven genes in 882 male lines with both measurements. The copy ratio is relative to the sample mean, so approximately 0.5 corresponds to one Y chromosome in a diploid line. Blue, red and grey denote integrated Y-positive, Y-negative and set-aside calls, respectively. **c,** Evidence supporting the calls in all 1,143 annotated male lines: 588 Y-positive, 383 Y-negative and 172 set aside. Counts are printed for segments containing at least 40 lines. Expression determines the integrated call, DNA supplies a call where RNA is unavailable, and definite disagreement sets a line aside. Copy ratios derive from whole-genome sequencing where available and otherwise from exome sequencing. Female lines define the expression reference and are not themselves classified. **d,** Genome-wide association of Y status with Chronos gene effect across 17,205 autosomal and X-linked genes, estimated with limma while adjusting for lineage and CRISPR library and blocking on donor. The x axis shows Y-negative minus Y-positive gene effect; negative values indicate stronger dependency in Y-negative lines. The y axis shows raw P values. Red marks genes passing Benjamini–Hochberg FDR <0.05 across all 17,205 tests, green marks other X-linked paralogs of Y-linked genes, and grey marks the remaining genes. **e,** Individual Chronos values for the four significant genes in 338 Y-positive (Y+) and 230 Y-negative (Y−) lines. Boxes show the median and interquartile range; whiskers extend to the most extreme observations within 1.5 interquartile ranges of the box. Each gene has its own y-axis scale. Printed Δ values and FDRs come from the adjusted model in d. These public Y-status calls do not establish the status of in-house stocks, and group-level gene-effect differences do not establish dependency in every individual line.'''
caption_block = f'''## main — {TITLE}

**Figure:** `{STEM}` | **Working ID:** {CODE} | **Results:** 3.4.1–3.4.2 | **Title:** {TITLE}.

[Selected PDF]({REL}.pdf) | [Selected SVG]({REL}.svg) | [Selected PNG]({REL}.png)

[PDF](figures/chapter3/{STEM}.pdf) | [Editable SVG](figures/chapter3/{STEM}.svg) | [PNG](figures/chapter3/{STEM}.png)

**{TITLE}.** {caption_body}

[Source and reproduction notes](figures/chapter3/selected/{SELECTED}/selection_notes.md)

'''
blocks = re.split(r'(?=^## )', text['Chapter3_figure_captions.md'], flags=re.M)
changed = []
for i, block in enumerate(blocks):
    if '**Working ID:** CH3-S3.4.1-F01' in block:
        blocks[i] = caption_block
        changed.append('first')
    elif '**Working ID:** CH3-S3.4.2-F01' in block:
        blocks[i] = ''
        changed.append('second')
assert sorted(changed) == ['first', 'second']
text['Chapter3_figure_captions.md'] = ''.join(blocks)

for name in ['Chapter3_results.md', 'Results_edit.md']:
    prefix, suffix = text[name].split(BOUNDARY, 1)
    if name == 'Results_edit.md':
        new_block = f'---\n\n**Figure {CODE} (main). {TITLE}.**\n\n[![{CODE}]({ROOT.as_posix()}/{REL}.png)]({ROOT.as_posix()}/{REL}.pdf)\n\n{caption_body}\n\n---'
        for identifier, replacement in [(CODE, new_block), ('CH3-S3.4.2-F01', '')]:
            pattern = rf'---\n\n\*\*Figure {re.escape(identifier)} \(main\).*?\n\n---'
            prefix, count = re.subn(pattern, lambda m: replacement, prefix, count=1, flags=re.S)
            assert count == 1, (name, identifier)
    else:
        new_block = f'[![Figure {CODE}. {TITLE}.]({REL}.png)]({REL}.pdf)\n\n**Figure {CODE} (main).** {TITLE}. [Full caption](Chapter3_figure_captions.md#{ANCHOR}).'
        for identifier, replacement in [(CODE, new_block), ('CH3-S3.4.2-F01', '')]:
            pattern = rf'^\[!\[Figure {re.escape(identifier)}\..*?\n\n\*\*Figure {re.escape(identifier)} \(main\).*?$'
            prefix, count = re.subn(pattern, lambda m: replacement, prefix, count=1, flags=re.S | re.M)
            assert count == 1, (name, identifier)
    prefix = prefix.replace('[Figure CH3-S3.4.2-F01, panels a and b]', '[Figure CH3-S3.4.1-F01, panels d and e]')
    for old_stem in old_stems:
        prefix = prefix.replace(f'figures/chapter3/selected/{old_stem}/{old_stem}', REL)
    prefix = re.sub(r'\n{4,}', '\n\n', prefix)
    assert 'CH3-S3.4.2-F01' not in prefix
    assert prefix.count('**Figure CH3-S3.4.1-F01 (main)') == 1
    assert len(re.findall(r'!\[', prefix)) == 11
    text[name] = prefix + BOUNDARY + suffix

note = f'''\n\n## DepMap main-figure combination — 5 October 2026

At the author's request, the former Y-status F01 (3.4.1) and dependency F01 (3.4.2) are combined as {CODE}. The first three panels retain Y-status classification; the former dependency panels a and b become d and e. Descriptive on-plot titles and explanatory footnotes were removed, with the latter information retained in the central caption. Gene names, sample counts, adjusted effects and FDRs remain visible. The new 9.2 × 5.65 inch figure uses Open Sans and the thesis blue/red condition colours. This reduces the active DepMap set to four main and seven supplementary figures. Prose references and both Results figure blocks were updated; no statistical analysis was refitted. The original exports remain preserved as superseded records. Reproduction: [plotting script](figures/chapter3/code/{STEM}.py); [saved source data](figures/chapter3/source_data/depmap_y_status_dependencies/README.md).
'''
text['Chapter3_source_notes.md'] = text['Chapter3_source_notes.md'].rstrip() + note

archive = ROOT / 'archive' / ('depmap_combined_f01_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
archive.mkdir(parents=True)
for name, value in raw.items():
    assert (ROOT/name).read_bytes() == value, f'Intervening edit: {name}'
    shutil.copy2(ROOT/name, archive/name)
for old_stem in old_stems:
    shutil.copytree(FIG/'selected'/old_stem, archive/old_stem)

selected_dir = FIG/'selected'/SELECTED
selected_dir.mkdir(exist_ok=True)
for ext in ['pdf', 'svg', 'png']:
    shutil.copy2(FIG/f'{STEM}.{ext}', selected_dir/f'{SELECTED}.{ext}')
    assert hashlib.sha256((FIG/f'{STEM}.{ext}').read_bytes()).digest() == hashlib.sha256((selected_dir/f'{SELECTED}.{ext}').read_bytes()).digest()
(selected_dir/'selection_notes.md').write_text(f'''# Combined DepMap main figure

Active working identifier: {CODE}. Combines the former CH3-S3.4.1-F01 (a–c) and CH3-S3.4.2-F01 (now d–e). Canonical stem: `{STEM}`.

The source measurements and model statistics are unchanged. The layout was redrawn from the saved panel CSVs, using Open Sans and the thesis condition colours; interpretive plot titles and in-figure footnotes were removed. Gene labels, sample counts and primary-model effects/FDRs were retained. Display jitter is deterministic. All points are retained. Boxplots use the original median, quartiles and 1.5-IQR convention, with separate gene y-axis scales.

- [Active caption](../../../../Chapter3_figure_captions.md#{ANCHOR})
- [Reproduction code](../../code/{STEM}.py)
- [Source data and hashes](../../source_data/depmap_y_status_dependencies/README.md)

Superseded figures and captions were archived in `{archive.relative_to(ROOT).as_posix()}`. Statistical models were not refitted.
''', encoding='utf-8')
for old_stem in old_stems:
    path = FIG/'selected'/old_stem/'selection_notes.md'
    existing = path.read_text(encoding='utf-8')
    path.write_text(f'> Superseded on 5 October 2026 by [the combined {CODE}](../{SELECTED}/{SELECTED}.pdf). This file preserves the original export provenance.\n\n' + existing, encoding='utf-8')

(FIG/'source_data/depmap_y_status_dependencies/README.md').write_text(f'''# Source data — combined DepMap Y-status and dependencies figure

Active figure: [PDF](../../{STEM}.pdf). [Caption](../../../../Chapter3_figure_captions.md#{ANCHOR}). [Reproduction script](../../code/{STEM}.py).

The five CSVs are exact copies from the 5 October 2026 DepMap analysis export. Their original names indicate the source panels; `provenance.json` records SHA-256 checksums and the final panel mapping. a: RNA histograms (1,616 profiles); b: RNA–DNA comparison (882 male lines); c: call tiers (1,143 male lines); d: gene associations (17,205 tests); e: four gene-effect distributions (568 lines per gene). No analysis models or P values were recalculated.

Reproduce with Python (NumPy, pandas, Matplotlib) using the script above. It loads the bundled thesis Open Sans font and exports matching PDF, SVG and 300-dpi PNG files at 9.2 × 5.65 inches. Original figure titles and footnotes are moved into the active caption; no panel was discarded. Gene names in e are identifiers, not descriptive plot titles.
''', encoding='utf-8')

for name, value in text.items():
    newline = '\r\n' if b'\r\n' in raw[name] else '\n'
    output = value.replace('\n', newline).encode('utf-8')
    if name in ['Chapter3_results.md', 'Results_edit.md']:
        assert output.split(BOUNDARY.encode(), 1)[1] == raw[name].split(BOUNDARY.encode(), 1)[1]
    (ROOT/name).write_bytes(output)

report = {'archive': str(archive), 'active_id': CODE, 'retired_id': 'CH3-S3.4.2-F01',
          'active_depmap_figures': 11, 'active_main': 4, 'active_supplementary': 7,
          'caption_and_references_updated': True, 'later_results_preserved': True,
          'statistical_models_refitted': False,
          'documents': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in doc_names}}
(FIG/'source_data/depmap_y_status_dependencies/integration_checks.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
