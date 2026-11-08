"""Validate adopted figure links, numerical provenance and caption synchrony."""
from pathlib import Path
import csv, difflib, hashlib, json, re, shutil
import xml.etree.ElementTree as ET
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]; FIG = ROOT/'figures/chapter3'
ARCHIVE = ROOT/'archive/regional_y_ypar_adoption_20261028'
DATA = FIG/'source_data/regional_y_atlases_ypar_20261028'
adopted = json.loads((HERE/'adopted_figures.json').read_text())
captions = json.loads((HERE/'captions.json').read_text(encoding='utf-8'))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def backup(path):
    dest = ARCHIVE/path.relative_to(ROOT)
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(path, dest)

pilot = ROOT.parent/'WGS/reanalysis/Y_PAR_masked_20261027'
dest = DATA/'dosage_check'; dest.mkdir(exist_ok=True)
for name in ['nonPAR_XY_dosage.csv', 'conditional_PAR_predictions.csv', 'provenance.json']:
    shutil.copy2(pilot/'atlas_sex_dosage_20261028'/name, dest/name)
shutil.copy2(pilot/'atlas_sex_dosage_20261028.py', dest/'atlas_sex_dosage_20261028.py')

# Put adoption status before historical records so older claims are not read as current.
for path, target in [
    (FIG/'source_data/regional_y_atlases_normalised/README.md', '../regional_y_atlases_ypar_20261028/README.md'),
    (ROOT/'research/model_development/ht1080_du145_ypar_20261027/README.md', '../../../figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md'),
    (ROOT/'research/model_development/a549_dense_bins_20261027/README.md', '../../../figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md')]:
    backup(path); text = path.read_text(encoding='utf-8')
    notice = '> **Adoption update — 28 October 2026:** The three realignment atlases are now adopted into draft 2 as Figs 3.4–3.6. See the [current source and reproduction record]('+target+'). The text and checks below describe the earlier analysis/preview state. The separate-track supplement has not been realigned.\n\n'
    if not text.startswith(notice): path.write_text(notice+text, encoding='utf-8')

agents = ROOT/'AGENTS.md'
note = '''\n\n# Regional atlas realignment adoption in draft 2\n\nAuthor approval, 28 October 2026: adopt the three realigned regional Y atlases into draftv2/chapter3/chapter3_results.md and add sex-chromosome karyotype/dosage context and conditional expected PAR retention to their captions. Current draft 2 numbering is HT1080 Fig 3.4, A549 Fig 3.5 and DU145 Fig 3.6; preserve it. Central captions retain historical numbers 3.10/3.12/3.14 with an explicit current-draft mapping. Use figures/chapter3/code/make_regional_y_atlas.py and frozen source_data/regional_y_atlases_ypar_20261028. These use Y-PAR-masked realignments, combined X-mapped PAR coverage and the approved broader 100-kb filter. RNA values and gene maps are unchanged. Local XY/XXY-equivalent dosage models are not confirmed cytogenetic karyotypes; DU145's rearranged Y and changed X dosage require a qualified expectation. Do not regenerate old atlases from the historical unmasked-reference tables. The separate-track supplement remains historical and is labelled accordingly. Source/reproduction notes and archived predecessors are linked in research/model_development/atlas_adoption_20261028/README.md.\n'''
if '# Regional atlas realignment adoption in draft 2' not in agents.read_text(encoding='utf-8'):
    with agents.open('a', encoding='utf-8') as f: f.write(note)

results_path = ROOT/'draftv2/chapter3/chapter3_results.md'
master_path = ROOT/'Chapter3_figure_captions.md'
results = results_path.read_text(encoding='utf-8'); master = master_path.read_text(encoding='utf-8')
checks = {}
for cell, info in adopted.items():
    stem=info['stem']; selected=info['selected']; folder=FIG/'selected'/selected
    assert results.count(captions[cell]) == master.count(captions[cell]) == 1
    assert results.count(info['inline_png']) == 1
    assert (folder/info['inline_png']).exists()
    marker = f'<!-- Caption source: Chapter3_figure_captions.md; figure: {stem}.'
    end = results.index(marker); start = results.rfind('[![', 0, end)
    block = results[start:end]
    assert f'Fig {info["draft_number"]}.' in block
    for ext in ['pdf', 'svg', 'png']:
        assert sha(folder/f'{selected}.{ext}') == sha(FIG/f'{stem}.{ext}') == info['exports'][ext]
    assert sha(folder/info['inline_png']) == info['exports']['png']
    with Image.open(folder/f'{selected}.png') as im:
        assert im.size[0] > 4000 and im.size[1] > 3900
        png_size = list(im.size)
    svg = (folder/f'{selected}.svg').read_text(encoding='utf-8')
    ET.fromstring(svg)
    assert 'Open Sans' in svg and '<text' in svg
    assert 'PAR1 (X)' in svg and 'PAR2 (X)' in svg
    assert 'PAR expression is not Y-specific' not in svg
    assert 'X-PAR measures combined' not in svg
    val = json.loads((FIG/'source_data/regional_y_atlases_ypar_20261028_rendered'/cell.lower()/'validation.json').read_text())
    assert val['RNA_values_and_means_unchanged'] and val['no_points_clipped']
    for name, expected in info['frozen_sources'].items():
        assert sha(DATA/cell.lower()/name) == expected
    checks[cell] = dict(caption_synchronised=True, links_exist=True, matching_exports=True,
                        editable_svg=True, png_pixels=png_size, sources_unchanged=True,
                        DNA_regions=val['regions'], PAR_regions=val['X_PAR_regions'])
for file in [results_path, master_path]:
    old=(ARCHIVE/file.relative_to(ROOT)).read_text(encoding='utf-8')
    diff=''.join(difflib.unified_diff(old.splitlines(keepends=True), file.read_text(encoding='utf-8').splitlines(keepends=True),
                                   fromfile='before/'+file.name, tofile='after/'+file.name))
    (HERE/(file.stem+'.diff')).write_text(diff, encoding='utf-8')
    # Figure labels and figure count must be unchanged by adoption.
    titles=r'^\*\*(?:Supp )?Fig [0-9.]+\. [^\n]+'
    assert re.findall(titles, old, re.M) == re.findall(titles, file.read_text(encoding='utf-8'), re.M)

# Index entries must agree with the actual selected PDF hashes.
with (FIG/'selected/selected_figure_index.tsv').open(newline='', encoding='utf-8') as f:
    index=list(csv.DictReader(f, delimiter='\t'))
for info in adopted.values():
    row=next(r for r in index if r['source_stem']==info['stem'])
    assert row['sha256_pdf']==info['exports']['pdf']

(HERE/'validation.json').write_text(json.dumps(dict(status='passed', figures=checks,
    figure_titles_and_counts_preserved=True, numerical_equivalence_to_previews='checked in adopt.py',
    manual_visual_review='All three complete previews inspected; no clipping or label overlap found',
    untouched_other_figures='Separate-track caption provenance clarified; no other figure regenerated',
    current_results_sha256=sha(results_path), current_caption_master_sha256=sha(master_path),
    renderer_sha256=sha(FIG/'code/make_regional_y_atlas.py')), indent=2), encoding='utf-8')
print('PASS: three figures, nine exports, content-addressed previews, source hashes, index hashes, caption synchrony and unchanged figure numbering.')
