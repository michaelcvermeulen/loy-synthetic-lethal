"""Archive the approved predecessors and freeze the new atlas render inputs."""
from pathlib import Path
import csv, hashlib, json, shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FIG = ROOT / 'figures/chapter3'
ARCHIVE = ROOT / 'archive/regional_y_ypar_adoption_20261028'
assert not ARCHIVE.exists(), 'Do not replace the pre-adoption archive.'
ARCHIVE.mkdir(parents=True)
DATA = FIG / 'source_data/regional_y_atlases_ypar_20261028'
STAGING = FIG / 'drafts/regional_y_ypar_adoption_20261028'
CONFIG = {
    'HT1080': ('CH3_S3_4_5_1_F04_ht1080_regional_y', '3.4', '3.10',
               FIG/'drafts/ht1080_du145_ypar_20261027/ht1080'),
    'A549': ('CH3_S3_4_5_2_F02_a549_regional_y', '3.5', '3.12',
             FIG/'drafts/a549_dense_bins_20261027/100kb'),
    'DU145': ('CH3_S3_4_5_3_F02_du145_regional_y', '3.6', '3.14',
              FIG/'drafts/ht1080_du145_ypar_20261027/du145')}
PAIRS = {'HT1080': ('HT1080-sgCTRL-D0', 'HT1080-LOY-D0'),
         'A549': ('A549-CTRL2', 'A549-LOY'), 'DU145': ('LNCaP-sgCTRL', 'LNCaP-LOY-2X')}

def archive(path):
    if not path.exists(): return
    dest = ARCHIVE / path.relative_to(ROOT)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if path.is_dir(): shutil.copytree(path, dest)
    else: shutil.copy2(path, dest)

for rel in ['Chapter3_figure_captions.md', 'draftv2/chapter3/chapter3_results.md',
            'draftv2/chapter3/chapter3_figure_numbering.md', 'AGENTS.md',
            'figures/chapter3/code/make_regional_y_atlas.py',
            'figures/chapter3/selected/selected_figure_index.tsv',
            'figures/chapter3/source_data/regional_y_atlases_normalised',
            'research/model_development/ht1080_du145_ypar_20261027/build_renderer.py']:
    archive(ROOT / rel)

inputs = ['qdnaseq_bins.csv', 'ddpcr_target_loci.csv', 'T2T_chrY_cytobands.csv',
          'gene_mapping_plotted.csv', 'relative_expression_per_sample.csv',
          'relative_expression_group_means.csv', 'sample_metadata.csv',
          'binned_coverage.csv', 'autosomal_normalization.csv']
manifest = {}
for cell, (selected, draft_num, master_num, preview) in CONFIG.items():
    slug = cell.lower(); stem = f'ch03_main_{slug}_lpwgs_chrY_atlas'
    archive(FIG/'selected'/selected)
    for ext in ['pdf', 'svg', 'png']: archive(FIG/f'{stem}.{ext}')
    archive(FIG/f'{stem}_notes.md')
    dest = DATA/slug; dest.mkdir(parents=True)
    source = preview/'source_data/regional_y_atlases_normalised'/slug
    for name in inputs: shutil.copy2(source/name, dest/name)
    annotation = source/'bin_annotation.csv'
    if not annotation.exists():
        annotation = ROOT/'research/model_development/wgs_filter_rerun_20261027/bin_annotation.csv'
    shutil.copy2(annotation, dest/'bin_annotation.csv')
    (dest/'sample_map.json').write_text(json.dumps(dict(control=PAIRS[cell][0], LOY=PAIRS[cell][1]), indent=2))
    if cell == 'A549':
        path = dest/'autosomal_normalization.csv'
        with path.open(newline='') as f: rows = list(csv.DictReader(f))
        rows = [r for r in rows if r['alignment'] == 'Y_PAR_masked']
        assert len(rows) == 2
        with path.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    manifest[cell] = dict(selected=selected, stem=stem, draft_number=draft_num,
                          historical_master_number=master_num, preview=str(preview),
                          frozen_sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.iterdir()})

source = ROOT/'research/model_development/ht1080_du145_ypar_20261027/make_atlas_preview.py'
text = source.read_text(encoding='utf-8')
def replace(a, b):
    global text
    assert text.count(a) == 1, (a, text.count(a))
    text = text.replace(a, b)
replace("ROOT = Path('A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS')", 'ROOT = Path(__file__).resolve().parents[3]')
replace("assert CELL in ['HT1080','DU145']\nW = Path('A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/research/model_development/ht1080_du145_ypar_20261027') / CELL / 'plot_inputs'", "W = FIG / 'source_data/regional_y_atlases_ypar_20261028' / SLUG")
replace("R = FIG / f'source_data/{SLUG}_y_expression_confirmation'", 'R = W')
replace("SHARED = FIG / 'source_data/ht1080_y_expression_confirmation'", 'SHARED = W')
replace("SNAP = OUT / 'source_data/regional_y_atlases_normalised' / SLUG", "SNAP = OUT / 'source_data/regional_y_atlases_ypar_20261028_rendered' / SLUG")
replace("FIG / 'source_data/ht1080_lpwgs_confirmation/ddpcr_target_loci.csv'", "W / 'ddpcr_target_loci.csv'")
replace("if CELL=='HT1080':", "if CELL in ['HT1080', 'A549']:")
replace("'status': 'New-alignment preview on the cell-specific template; broad 100-kb filter; not adopted into draft 2'", "'status': 'Adopted 2026-10-28: Y-PAR-masked realignment, X-PAR coverage and broad 100-kb filter'")
replace("'cell_line': CELL,", "'cell_line': CELL,\n    'current_draft_number': {'HT1080': '3.4', 'A549': '3.5', 'DU145': '3.6'}[CELL],\n    'removed_annotations': ['blue explanatory text below a', 'blue explanatory text below c'],")
(FIG/'code/make_regional_y_atlas.py').write_text(text, encoding='utf-8')
# Historical preview builders continue to use their original frozen template.
builder = ROOT/'research/model_development/ht1080_du145_ypar_20261027/build_renderer.py'
bt = builder.read_text(encoding='utf-8')
old = "source=ROOT/'figures/chapter3/code/make_regional_y_atlas.py'"
assert bt.count(old) == 1
builder.write_text(bt.replace(old, "source=HERE/'template_snapshot.py'"), encoding='utf-8')
(HERE/'adoption_manifest.json').write_text(json.dumps(dict(archive=str(ARCHIVE), staging=str(STAGING), figures=manifest), indent=2))
print('Archived predecessors and prepared the unified renderer for all three models.')
