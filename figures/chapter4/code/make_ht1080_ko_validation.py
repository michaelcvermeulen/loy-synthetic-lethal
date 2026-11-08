"""Grant-derived PCR panel and workbook-verified clone ddPCR measurements."""
from pathlib import Path
import hashlib
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
PROJECTS = ROOT.parent
sys.path.insert(0, str(PROJECTS / '.figure_python_deps'))
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import openpyxl
from PIL import Image
import pymupdf
from thesis_style import apply_style
from import_ht1080_ko_pcr import BOX, EXPECTED

FIG = ROOT / 'figures/chapter4'
DATA = FIG / 'source_data/ht1080_ko_validation_20261022'
STEM = 'ch04_supp_ht1080_ko_validation'
CLONES = ['LOY1', 'LOY3', 'LOY4', 'LOY5', 'LOY7', 'LOY8',
          'CTRL1', 'CTRL3', 'CTRL4', 'CTRL6', 'CTRL7', 'CTRL8']
ASSAYS = ['RPS4Y1 (Yp)', 'AMELY/AMELX (Yp)', 'KDM5D (Yq)']
COLOURS = ['#2a78d6', '#eb6834', '#1baf7a']
MARKERS = ['o', 's', '^']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    csv = PROJECTS / 'ddPCR/figures/DDX3X_clones_grant_source.csv'
    workbook = PROJECTS / 'ddPCR/data/LOY_clones.xlsx'
    gel = FIG / 'source_data/ht1080_ko_genomic_pcr_20261022/grant_figure5_original.tiff'
    inputs = {str(p): sha(p) for p in (csv, workbook, gel)}
    assert sha(gel) == EXPECTED
    for p in (csv, workbook, gel):
        shutil.copy2(p, DATA / p.name)
    d = pd.read_csv(csv)
    d['label'] = d.clone.str.removeprefix('DDX3X_')
    assert len(d) == 36 and set(d.label) == set(CLONES)
    assert not d.duplicated(['label', 'assay']).any()
    assert set(d.assay) == set(ASSAYS)
    values = openpyxl.load_workbook(workbook, data_only=True)['DDX3X_LOY']
    rows = {values.cell(r, 1).value: r for r in range(3, 22)}
    checks = []
    for x in d.itertuples():
        column = dict(zip(ASSAYS, [5, 6, 7]))[x.assay]
        cell = values.cell(rows[x.clone], column)
        assert np.isclose(float(cell.value), x.ratio, rtol=0, atol=1e-12)
        assert np.isclose(x.dosage, x.ratio / x.exp_ratio, rtol=0, atol=1e-12)
        checks.append({'clone': x.clone, 'assay': x.assay, 'cell': cell.coordinate,
                       'workbook_ratio': cell.value, 'saved_ratio': x.ratio,
                       'displayed_dosage': x.dosage})
    pd.DataFrame(checks).to_csv(DATA / 'workbook_verification.csv', index=False)
    d.to_csv(DATA / 'display_values.csv', index=False)

    apply_style(8)
    fig = plt.figure(figsize=(7.8, 4.25))
    panel = Image.open(gel).crop(BOX)
    panel.save(DATA / 'genomic_pcr_panel_native.png')
    a = fig.add_axes([.006, .24, .405, .66])
    a.imshow(panel, interpolation='none')
    a.axis('off')
    ax = fig.add_axes([.51, .19, .47, .72])
    positions = list(range(13, 7, -1)) + list(range(6, 0, -1))
    y = dict(zip(CLONES, positions))
    ax.axvline(1, color='#777777', lw=.7, zorder=0)
    ax.axvline(0, color='#BBBBBB', lw=.6, zorder=0)
    for yy in positions:
        ax.axhline(yy, color='#EEEEEE', lw=.5, zorder=0)
    for j, assay in enumerate(ASSAYS):
        s = d.loc[d.assay.eq(assay)]
        ax.scatter(s.dosage, s.label.map(y) + [.18, 0, -.18][j],
                   marker=MARKERS[j], s=20, facecolors='white',
                   edgecolors=COLOURS[j], linewidths=.9, zorder=3)
    ax.set_yticks(positions, CLONES)
    ax.tick_params(axis='y', length=0, pad=4)
    ax.set_xlim(-.035, 1.56)
    ax.set_ylim(.35, 14.3)
    ax.set_xticks([0, .5, 1, 1.5])
    ax.set_xlabel('Relative Y dosage', labelpad=4)
    ax.spines['left'].set_visible(False)
    ax.text(0, 13.9, 'LOY clones (n = 6)', fontsize=8, weight='semibold')
    ax.text(0, 7.05, 'Y-retaining clones (n = 6)', fontsize=8, weight='semibold')
    fig.text(.435, .91, 'b', fontsize=12, weight='semibold')
    handles = [Line2D([], [], marker=m, color=c, markerfacecolor='white',
                      markersize=4, linestyle='none', label=label)
               for m, c, label in zip(MARKERS, COLOURS, ASSAYS)]
    fig.legend(handles=handles, loc='lower right', bbox_to_anchor=(.993, .03),
               ncol=3, fontsize=7, columnspacing=.65, handletextpad=.3)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(FIG / f'{STEM}.{ext}', dpi=600, facecolor='white')
    plt.close(fig)
    doc = pymupdf.open(FIG / f'{STEM}.pdf')
    images = doc[0].get_images(full=True)
    intact = False
    import io
    for item in images:
        stored = Image.open(io.BytesIO(doc.extract_image(item[0])['image']))
        if stored.size == panel.size:
            expected = np.array(panel.convert('RGB'))
            actual = np.array(stored.convert('RGB'))
            # Matplotlib stores PDF image rows bottom-up and compensates with
            # the placement transform. The displayed gel remains upright.
            intact = np.array_equal(actual, expected) or np.array_equal(actual[::-1], expected)
    assert intact, 'PDF must retain the original gel pixels'
    doc.close()
    assert all(sha(Path(p)) == h for p, h in inputs.items())
    manifest = {
        'inputs_sha256': inputs, 'workbook_sheet': 'DDX3X_LOY',
        'clones': CLONES, 'measurements': 36, 'all_ratios_match_workbook': True,
        'all_dosages_match_saved_source': True, 'pdf_gel_pixels_unchanged': intact,
        'scaling': {'RPS4Y1 (Yp)': 'ratio / 0.5', 'AMELY/AMELX (Yp)': 'ratio / 1', 'KDM5D (Yq)': 'ratio / 0.5'},
        'interpretation': 'Scaled Y dosage; no percentage of LOY cells inferred.',
        'uncertainty': 'No confidence intervals inferred from concentration values.',
        'selection': 'The same 12 clones retained in the saved grant figure.',
        'ancestry': 'Source identifies DDX3X-targeted clones; chronology relative to WT pools remains unresolved.',
        'source_typography': 'Original gel labels retained; ddPCR panel uses Open Sans.',
        'exports_sha256': {ext: sha(FIG / f'{STEM}.{ext}') for ext in ['pdf', 'svg', 'png']},
    }
    (DATA / 'provenance_and_validation.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'png': str(FIG / f'{STEM}.png'), 'measurements_verified': 36, 'gel_pixels_unchanged': intact}))


if __name__ == '__main__':
    main()
