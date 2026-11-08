"""Adopt the three reviewed realignment atlases and synchronise their captions."""
from pathlib import Path
import csv, hashlib, json, re, shutil
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FIG = ROOT/'figures/chapter3'
manifest = json.loads((HERE/'adoption_manifest.json').read_text())
ARCHIVE = Path(manifest['archive']); STAGING = Path(manifest['staging'])
DATA = FIG/'source_data/regional_y_atlases_ypar_20261028'
RESULTS = ROOT/'draftv2/chapter3/chapter3_results.md'
MASTER = ROOT/'Chapter3_figure_captions.md'
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def backup(path):
    dest = ARCHIVE/path.relative_to(ROOT)
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(path, dest)
def write_tsv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

# Refuse to overwrite author changes made since this adoption began.
for path in [RESULTS, MASTER]:
    assert sha(path) == sha(ARCHIVE/path.relative_to(ROOT)), f'Author edited {path}; merge first.'
results = RESULTS.read_text(encoding='utf-8')
master = MASTER.read_text(encoding='utf-8')
original_results = results
original_master = master

KARYOTYPE = {
    'HT1080': 'Local sgCTRL dosage is consistent with an XY-equivalent model. Complete loss of Y-derived PAR with unchanged X dosage would retain approximately 50% of combined PAR signal (XY → X); partial Y retention can leave signal above this baseline. This dosage inference is not a cytogenetic karyotype.',
    'A549': 'Local sgCTRL dosage is consistent with an approximately XXY-equivalent model. Complete loss of Y-derived PAR with unchanged X dosage would retain approximately 67% of combined PAR signal (XXY → XX; 2/3), consistent with the observed PAR medians of 67% and 64%. This dosage inference is not a cytogenetic karyotype.',
    'DU145': 'The reference-stock karyotype comprises one X and a rearranged Y ([ATCC](https://www.atcc.org/products/htb-81)). An intact XY → X model would predict 50% PAR retention. Local X:Yq dosage (approximately 1.8:1) is closer to an XXY-equivalent model, which predicts approximately 67% PAR retention if X dosage is unchanged. Using Yq as the Y-dosage baseline and allowing for the higher X dosage in the depleted pool gives a conditional expectation of approximately 73–74%. Both calculations assume complete loss of Y-derived PAR; bulk lpWGS does not establish the local structural karyotype.'}

COMMON_B = ('**b,** HG002 T2T-Y ideogram with RNA-seq genes above and ddPCR target-gene spans below, including AMELY, DDX3Y and KDM5D. Panels a and b share coordinates and Yq12 compression; X-PAR bins in a are projected within the corresponding labelled PAR spans.')
COMMON_A = ('Grey points show corrected Y-bin copy number as a percentage of sgCTRL; red lines show regional medians. Teal points and lines show combined X/Y PAR coverage mapped to X, normalised to autosomal depth and then sgCTRL. The dashed line marks 100% of sgCTRL. Y bins require ≥60% non-N sequence, ≥30% mappability and ≥50 primary MAPQ ≥37 control reads; PAR bins require ≥20 reads and ≥50 reads per 100 kb. Boundary-crossing Y bins, centromere and Yq12 are excluded; hatching marks the latter two and Yq12 is compressed. PAR1 and PAR2 contain 24 and four bins, respectively. One DNA sample represents each population.')

captions = {}; adopted = {}
for cell, info in manifest['figures'].items():
    slug = cell.lower(); stem = info['stem']; selected = info['selected']
    snap = STAGING/'source_data/regional_y_atlases_ypar_20261028_rendered'/slug
    prior = Path(info['preview'])/'source_data/regional_y_atlases_normalised'/slug
    # Numerical adoption check: reproduce the approved previews exactly.
    for name, keys, cols in [
        ('regional_median_paired_ratios.csv', ['region'], ['n_bins', 'median_percent_of_sgCTRL']),
        ('X_PAR_regional_medians.csv', ['track'], ['n_bins', 'median_percent_of_sgCTRL']),
        ('paired_Y_bins_and_exclusions.csv', ['bin_id'], ['included', 'percent_of_sgCTRL']),
        ('X_PAR_bin_ratios.csv', ['track', 'start0', 'end'], ['included', 'percent_of_sgCTRL'])]:
        old = pd.read_csv(prior/name)
        if name == 'X_PAR_regional_medians.csv' and 'track' not in old.columns:
            old = old.rename(columns={old.columns[0]: 'track'})
        old = old.set_index(keys).sort_index()
        new = pd.read_csv(snap/name).set_index(keys).sort_index()
        pd.testing.assert_index_equal(old.index, new.index)
        for col in cols:
            assert np.allclose(old[col], new[col], equal_nan=True, rtol=1e-12, atol=1e-12), (cell, name, col)
    for name in ['relative_expression_per_sample.csv', 'relative_expression_group_means.csv',
                 'gene_mapping_plotted.csv', 'T2T_chrY_cytobands.csv', 'ddpcr_target_loci.csv']:
        assert sha(prior/name) == sha(snap/name), (cell, name)
    validation = json.loads((snap/'validation.json').read_text())
    assert validation['no_points_clipped'] and validation['RNA_values_and_means_unchanged']
    assert validation['DNA_display_limits'] == [-6, 156]
    assert validation['X_PAR_bins_displayed'] == 28
    dest = FIG/'selected'/selected
    for ext in ['pdf', 'svg', 'png']:
        source = STAGING/f'{stem}.{ext}'
        shutil.copy2(source, FIG/source.name)
        shutil.copy2(source, dest/f'{selected}.{ext}')
    # A content-addressed PNG refreshes the inline preview while preserving old links.
    preview_name = f'{selected}_ypar_{sha(dest/f"{selected}.png")[:12]}.png'
    shutil.copy2(dest/f'{selected}.png', dest/preview_name)
    root_snap = FIG/'source_data/regional_y_atlases_ypar_20261028_rendered'/slug
    if root_snap.exists():
        assert all(sha(p) == sha(root_snap/p.name) for p in snap.iterdir())
    else:
        shutil.copytree(snap, root_snap)

    guides = 'CP + Q53 + P47' if cell == 'DU145' else 'CP + Q53'
    a = f'**a,** lpWGS after two successive {guides} nucleofections, realigned with the Y copies of PAR masked. ' + COMMON_A
    if cell == 'DU145':
        a += ' Distal (~2.46–4.80 Mb) and proximal (~4.80–10.57 Mb) male-specific Yp are summarised separately; the division is not a precise structural breakpoint.'
    if cell == 'HT1080':
        a += ' The near-zero Yq median is labelled <1%; it does not establish complete absence of Yq DNA.'
    # Preserve all original panel-c sampling and normalisation language.
    marker = f'<!-- Caption source: Chapter3_figure_captions.md; figure: {stem}.'
    end = results.index(marker)
    start = results.rfind('**a,**', 0, end)
    old_body = results[start:end].strip()
    c = old_body[old_body.index('**c,**'):]
    body = '\n\n'.join([a, '**Sex-chromosome dosage and expected PAR retention.** ' + KARYOTYPE[cell], COMMON_B, c])
    captions[cell] = body
    results = results[:start] + body + '\n\n' + results[end:]
    image_pattern = rf'{re.escape(selected)}/{re.escape(selected)}_normalised_[0-9a-f]+\.png'
    results, n = re.subn(image_pattern, f'{selected}/{preview_name}', results)
    assert n == 1, (cell, n)

    # Master keeps its historical numbering; draft numbering is recorded explicitly.
    block_start = master.index(f'**Figure:** `{stem}`')
    block_end = master.index('\n## ', block_start)
    block = master[block_start:block_end]
    body_start = block.index('**a,**')
    body_end = block.index('[Source and reproduction notes]', body_start)
    block = block[:body_start] + body + '\n\n' + block[body_end:]
    block = re.sub(r'\[Source and reproduction notes\]\([^\n]+\)',
                   '[Source and reproduction notes](figures/chapter3/source_data/regional_y_atlases_ypar_20261028/README.md)', block)
    block = block.replace(f'**Chapter number:** Fig {info["historical_master_number"]}',
                          f'**Chapter number:** Fig {info["historical_master_number"]} | **Current draft 2 number:** Fig {info["draft_number"]}')
    master = master[:block_start] + block + master[block_end:]
    (dest/'selection_notes.md').write_text(
        f'# {selected}\n\nAdopted 28 October 2026 at the author\'s request: new Y-PAR-masked alignments, broader 100-kb filter and X-mapped combined PAR coverage.\n\n'
        f'Working source: `{stem}`. Current draft 2: **Fig {info["draft_number"]}**; central historical number: Fig {info["historical_master_number"]}.\n\n'
        'Captions include conditional sex-chromosome dosage models, not claims of cytogenetically confirmed local karyotypes. RNA values and gene maps are unchanged.\n\n'
        '[Source, filters and reproduction](../../source_data/regional_y_atlases_ypar_20261028/README.md). '
        '[Central caption master](../../../../Chapter3_figure_captions.md).\n\n'
        f'Predecessors archived at `{ARCHIVE.relative_to(ROOT).as_posix()}`.\n', encoding='utf-8')
    (FIG/f'{stem}_notes.md').write_text(
        f'# {cell} regional Y atlas\n\nCurrent draft 2: Fig {info["draft_number"]}. Adopted 28 October 2026.\n\n'
        '[Authoritative source and reproduction notes](source_data/regional_y_atlases_ypar_20261028/README.md).\n\n'
        'Build with `python figures/chapter3/code/make_regional_y_atlas.py '+cell+'` from the thesis root. '
        'The script reads frozen input tables; it does not rerun alignment or statistical fitting.\n', encoding='utf-8')
    adopted[cell] = dict(info, inline_png=preview_name,
        exports={ext:sha(dest/f'{selected}.{ext}') for ext in ['pdf', 'svg', 'png']},
        Y_regions=validation['regions'], PAR_regions=validation['X_PAR_regions'])

# Only update narrative claims directly tied to the adopted atlas measurements.
changes = {
    'and below 1% in euchromatic Yq, where estimates were limited by the numerical analysis floor':
        'and below 1% in euchromatic Yq',
    'hypotriploid (XXY) A549 lung adenocarcinoma cells':
        'hypotriploid A549 lung adenocarcinoma cells with approximately XXY-equivalent sex-chromosome dosage',
    'retention reduced to 2.0% of sgCTRL across male-specific Yp and 1.3% across euchromatic Yq':
        'retention reduced to 2.1% of sgCTRL across male-specific Yp and 2.7% across euchromatic Yq',
    'proximal Yp (5.3% median DNA retention relative to sgCTRL':
        'proximal Yp (6.1% median DNA retention relative to sgCTRL',
    'euchromatic Yq (1.8% DNA retention and <0.6% RNA expression)':
        'euchromatic Yq (2.2% DNA retention and <0.6% RNA expression)',
    '~2.46–4.80 Mb; 51% median DNA retention':
        '~2.46–4.80 Mb; 54% median DNA retention'}
for old, new in changes.items():
    assert results.count(old) == 1, old
    results = results.replace(old, new)

# The separate-track supplement has not been realigned; describe its provenance
# accurately instead of claiming its old bins are identical to the new atlases.
old = 'in the same usable 100-kb bins included in the corresponding main atlases; coloured lines show regional medians within each population.'
new = 'in the earlier unmasked-reference analysis and its original bin selection; coloured lines show regional medians within each population. These historical tracks differ from the updated main atlases, which use Y-PAR-masked realignments, revised filtering and X-mapped PAR coverage.'
assert results.count(old) == master.count(old) == 1
results = results.replace(old, new); master = master.replace(old, new)

RESULTS.write_text(results, encoding='utf-8')
MASTER.write_text(master, encoding='utf-8')
for path, key in [(FIG/'selected/selected_figure_index.tsv', 'source_stem'), (FIG/'figure_manifest.tsv', 'file')]:
    backup(path)
    with path.open(newline='', encoding='utf-8-sig') as f: rows = list(csv.DictReader(f, delimiter='\t'))
    for row in rows:
        for cell, info in adopted.items():
            if (row[key] == info['stem']) if key == 'source_stem' else (Path(row[key]).stem == info['stem']):
                row['status'] = f'selected; Y-PAR-masked realignment atlas adopted 2026-10-28; draft 2 Fig {info["draft_number"]}'
                if key == 'source_stem': row['sha256_pdf'] = info['exports']['pdf']
                else:
                    row['sha256'] = info['exports'][Path(row[key]).suffix[1:]]
                    row['original_source'] = 'figures/chapter3/code/make_regional_y_atlas.py'
    write_tsv(path, rows)

(HERE/'adopted_figures.json').write_text(json.dumps(adopted, indent=2), encoding='utf-8')
(HERE/'captions.json').write_text(json.dumps(captions, indent=2, ensure_ascii=False), encoding='utf-8')
print('Adopted HT1080, A549 and DU145 as draft 2 Figs 3.4–3.6; caption master and targeted Results values synchronised.')
