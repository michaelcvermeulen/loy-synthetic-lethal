"""HT1080 passage ddPCR: exported concentration ratios scaled by reference copies.

The original droplet-ratio analysis is independently reconstructed as provenance.
Current plot values use the concentration-based Chapter 3 copy-number convention.
"""
from pathlib import Path
from html.parser import HTMLParser
import hashlib
import json
import re
import sys

THESIS = Path(__file__).resolve().parents[3]
PROJECTS = THESIS.parent
sys.path.insert(0, str(PROJECTS / '.figure_python_deps'))
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image
import fitz
from thesis_style import apply_style, condition_colours

OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_ddpcr_passage'
REVIEW = OUT / 'review'
RAW = PROJECTS.parent / 'data/ddPCR'
HTML = PROJECTS.parent / 'code/ddPCR_HT1080_passage.html'
RMD = PROJECTS / 'ddPCR/analysis/ddPCR_HT1080_passage.Rmd'
STEM = 'ch04_supp_ht1080_ddpcr_passage'
RUNS = [
    ('AMELX Craig_NOV_13_2025-11-12-22-50', '2025-11-13'),
    ('AMELX Craig_NOV7_2025_2025-11-07-15-53', '2025-11-07'),
    ('AMELX Craig_FEB_24_2026_2026-02-25-17-16', '2026-02-25'),
]
PASSAGES = ['D0', 'P1', 'P2', 'P3', 'P4']
ARM = {'Yp': 'Yp (AMELY/AMELX)', 'Yq': 'Yq (KDM5D or DDX3Y/RNaseP)'}
VALUE = 'estimated_target_copies_per_cell'
REF_COPIES = {'AMELX': 1, 'RNaseP': 2}


def concentration_measurements(wells, sources):
    """Join each measured target and its reference within the same run/well."""
    records = []
    for dirname, run in RUNS:
        path = RAW / dirname / 'data_table.csv'
        sources.add(path)
        raw = pd.read_csv(path)
        col = next(c for c in raw if c.startswith('Conc(copies/'))
        for row in wells.loc[wells.run.eq(run)].itertuples():
            y = raw.loc[raw.Well.eq(row.well) & raw.Target.eq(row.target)]
            ref = raw.loc[raw.Well.eq(row.well) & raw.Target.eq(row.reference)]
            assert len(y) == len(ref) == 1, (run, row.well, row.target)
            assert y.iloc[0]['Sample description 1'] == ref.iloc[0]['Sample description 1'] == row.sample_desc
            yc, rc = float(y.iloc[0][col]), float(ref.iloc[0][col])
            assert np.isfinite(yc) and np.isfinite(rc) and yc >= 0 and rc > 0
            ratio = yc / rc
            copies = REF_COPIES[row.reference]
            records.append(dict(name=row.name, run=run, concentration_source=str(path),
                                target_concentration_copies_per_ul=yc,
                                reference_concentration_copies_per_ul=rc,
                                concentration_ratio=ratio, reference_copy_number=copies,
                                estimated_target_copies_per_cell=ratio * copies))
    result = wells.merge(pd.DataFrame(records), on=['name', 'run'], validate='one_to_one')
    assert len(result) == len(wells) == 29
    assert np.allclose(result[VALUE], result.concentration_ratio * result.reference_copy_number,
                       rtol=0, atol=1e-14)
    return result


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Tables(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tables, self.t, self.row, self.cell = [], None, None, None

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self.t = []
        elif tag == 'tr' and self.t is not None:
            self.row = []
        elif tag in ('td', 'th') and self.row is not None:
            self.cell = ''

    def handle_data(self, s):
        if self.cell is not None:
            self.cell += s

    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.cell is not None:
            self.row.append(self.cell.strip())
            self.cell = None
        elif tag == 'tr' and self.row is not None:
            self.t.append(self.row)
            self.row = None
        elif tag == 'table' and self.t is not None:
            self.tables.append(self.t)
            self.t = None


def prepare():
    DATA.mkdir(parents=True, exist_ok=True)
    sources = {HTML, RMD, Path(__file__), THESIS / 'AGENTS.md',
               THESIS / 'figures/chapter3/code/make_ht1080_ddpcr_rounds.py'}
    rows = []
    for dirname, run in RUNS:
        folder = RAW / dirname
        plates = []
        for f in sorted(folder.glob('*.csv')):
            with f.open(encoding='utf-8-sig') as handle:
                if handle.readline().startswith('ddplate'):
                    plates.append(f)
        if plates:
            meta_path = plates[0]
            plate = pd.read_csv(meta_path, skiprows=4).drop_duplicates('Well')
        else:
            meta_path = folder / 'data_table.csv'
            plate = pd.read_csv(meta_path).rename(columns={'Target': 'TargetName'}).drop_duplicates('Well')
        sources.add(meta_path)
        plate = plate.set_index('Well')
        for f in sorted(folder.glob('*_Amplitude.csv')):
            if not f.stat().st_size:
                continue
            well = re.findall(r'(?:^|_)([A-H][0-9]{2})(?=_|$)', f.stem)[-1]
            if well not in plate.index:
                continue
            meta = plate.loc[well]
            desc, target = str(meta['Sample description 1']), str(meta['TargetName'])
            name = well + '_' + desc
            if not re.search(r'HT1080|^[A-H][0-9]{2}_(D0|P[1-4])', name, re.I):
                continue
            if target not in ['AMELY', 'AMELX', 'KDM5D', 'DDX3Y']:
                continue
            sources.add(f)
            amp = pd.read_csv(f, skiprows=4, usecols=[0, 1, 2, 3])
            y = pd.to_numeric(amp.iloc[:, 2], errors='coerce').to_numpy()
            x = pd.to_numeric(amp.iloc[:, 3], errors='coerce').to_numpy()
            valid = np.isfinite(x) & np.isfinite(y)
            assert np.isin(x[valid], [0, 1]).all() and np.isin(y[valid], [0, 1]).all()
            y, x = y[valid], x[valid]
            calls = {'Y+': int(((y == 1) & (x == 0)).sum()),
                     'Y-': int(((y == 0) & (x == 1)).sum()),
                     'X+Y': int(((y == 1) & (x == 1)).sum()),
                     'empty': int(((y == 0) & (x == 0)).sum())}
            m = re.search(r'D0|P[1-4]|Parental', desc)
            passage = m.group() if m else ('Parental' if 'HT1080' in desc else None)
            condition = next((s for s in ['LOY', 'CTRL', 'Parental'] if s in desc), 'Other')
            rows.append(dict(name=name, well=well, run=run, sample_desc=desc,
                             target=target, reference='AMELX' if target in ['AMELY', 'AMELX'] else 'RNaseP',
                             arm='Yp' if target in ['AMELY', 'AMELX'] else 'Yq',
                             passage=passage, condition=condition,
                             unclassified_excluded=int((~valid).sum()), **calls))
    wells = pd.DataFrame(rows)
    wells['probe_arm'] = wells.arm.map(ARM)
    wells['total'] = wells[['Y+', 'Y-', 'X+Y', 'empty']].sum(axis=1)
    wells['ratio'] = (wells['Y+'] + wells['X+Y']) / (wells['Y-'] + wells['X+Y'])
    pooled = wells.groupby(['passage', 'condition', 'arm', 'probe_arm'], sort=False).agg(
        **{k: (k, 'sum') for k in ['Y+', 'Y-', 'X+Y', 'empty', 'total']},
        n_wells=('name', 'size'), targets=('target', lambda x: ';'.join(sorted(set(x)))),
        wells=('name', lambda x: ';'.join(x))).reset_index()
    pooled['ratio'] = (pooled['Y+'] + pooled['X+Y']) / (pooled['Y-'] + pooled['X+Y'])
    parser = Tables()
    parser.feed(HTML.read_text(encoding='utf-8'))
    checks = {}
    for i, (frame, keys) in enumerate([(wells, ['name', 'run']), (pooled, ['passage', 'condition', 'probe_arm'])]):
        t = parser.tables[i]
        saved = pd.DataFrame(t[1:], columns=t[0])
        assert len(saved) == len(frame), (len(saved), len(frame))
        compared = frame.merge(saved, on=keys, how='outer', suffixes=('', '_saved'), indicator=True, validate='one_to_one')
        assert compared['_merge'].eq('both').all()
        diff = (compared.ratio - pd.to_numeric(compared.ratio_saved)).abs()
        assert diff.max() <= 0.0000500001, diff.max()
        assert (compared.total == pd.to_numeric(compared.total_saved)).all()
        checks['per_well' if i == 0 else 'pooled'] = {'n': len(compared), 'max_difference_from_4dp_saved_ratio': float(diff.max()), 'total_droplets_exact_match': True}
        compared.to_csv(DATA / ('per_well_original_comparison.csv' if i == 0 else 'pooled_original_comparison.csv'), index=False)
    wells.to_csv(DATA / 'reconstructed_all_wells.csv', index=False)
    pooled.to_csv(DATA / 'reconstructed_all_pooled.csv', index=False)
    keep = lambda d: d.loc[d.condition.isin(['CTRL', 'LOY']) & d.passage.isin(PASSAGES)].copy()
    wells, pooled = keep(wells), keep(pooled)
    assert not wells.condition.eq('Parental').any()
    wells.to_csv(DATA / 'legacy_droplet_wells.csv', index=False)
    pooled.to_csv(DATA / 'legacy_droplet_pooled.csv', index=False)
    wells = concentration_measurements(wells, sources)
    summary = wells.groupby(['passage', 'condition', 'arm', 'probe_arm'], sort=False).agg(
        estimated_target_copies_per_cell=(VALUE, 'mean'),
        mean_concentration_ratio=('concentration_ratio', 'mean'),
        reference_copy_number=('reference_copy_number', 'first'),
        n_wells=('name', 'size'), targets=('target', lambda x: ';'.join(sorted(set(x)))),
        wells=('name', lambda x: ';'.join(x))).reset_index()
    assert np.allclose(summary[VALUE], summary.mean_concentration_ratio * summary.reference_copy_number,
                       rtol=0, atol=1e-14)
    assert len(summary) == len(pooled) == 15
    wells.to_csv(DATA / 'plotted_wells.csv', index=False)
    summary.to_csv(DATA / 'plotted_summary.csv', index=False)
    transition = pooled[['passage', 'condition', 'arm', 'ratio']].rename(columns={'ratio': 'previous_pooled_positive_droplet_ratio'}).merge(
        summary, on=['passage', 'condition', 'arm'], validate='one_to_one')
    transition.to_csv(DATA / 'metric_transition.csv', index=False)
    missing = []
    for condition in ['CTRL', 'LOY']:
        for arm in ['Yp', 'Yq']:
            observed = set(pooled.loc[(pooled.condition == condition) & (pooled.arm == arm), 'passage'])
            missing.extend({'condition': condition, 'arm': arm, 'passage': p} for p in PASSAGES if p not in observed)
    manifest = {
        'request': 'Normalise Yp/Yq assays; use Chapter 3 reference-copy scaling. Retain horizontal condition facets and omit Parental.',
        'metric': 'Instrument target concentration / same-well reference concentration, multiplied by reference-locus copies per cell.',
        'reference_copy_numbers': REF_COPIES,
        'reference_assumption': 'Reference-locus copy numbers are assumed unchanged after Y depletion.',
        'summary': 'Arithmetic mean of per-well copy estimates within passage, condition and Y arm, matching the mean used by the original concentration-ratio variant.',
        'metric_change': 'Replaces pooled positive-droplet ratios with instrument concentration-based copy estimates; raw droplet data and summaries retained separately.',
        'not_applied': ['Normalisation to measured sgCTRL', 'Inferential tests', 'Missing-value interpolation'],
        'assay_caveat': 'Original Yq summaries combine DDX3Y and KDM5D. D0 CTRL is DDX3Y; D0 LOY and P4 CTRL combine DDX3Y and KDM5D; P3 and P4 LOY use KDM5D. Arm summaries must not be described as repeated measurements of one identical assay.',
        'yq_assays_by_group': pooled.loc[pooled.arm.eq('Yq'), ['passage', 'condition', 'targets', 'n_wells']].to_dict('records'),
        'plotted_wells': len(wells), 'plotted_summary_points': len(summary),
        'concentration_checks': {'matched_run_well_target_reference_pairs': len(wells), 'all_values_finite': True,
                                 'scaling_identity_verified': True, 'all_original_wells_retained': True},
        'missing_measurements': missing, 'original_report_validation': checks,
        'sources': [{'path': str(p), 'sha256': sha(p)} for p in sorted(sources)],
    }
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    return wells, summary, manifest


def draw_panel(ax, wells, summary, condition, font_size=8.5, show_y=True):
    """One condition facet; identical limits, no markers for missing passages."""
    colours = condition_colours()
    colour = {'CTRL': colours['Y-retaining'], 'LOY': colours['Y-depleted']}
    xmap = {p: i for i, p in enumerate(PASSAGES)}
    ax.set(xlim=(-.2, 4.3), ylim=(-.055, 1.16), xticks=range(5), xticklabels=PASSAGES,
           yticks=[0, .25, .5, .75, 1],
           ylabel='Estimated Y-target\ncopies per cell' if show_y else '')
    ax.set_title('sgCTRL' if condition == 'CTRL' else 'Y-depleted', color=colour[condition],
                 fontsize=9, fontweight='semibold', pad=7)
    ax.spines['bottom'].set_bounds(0, 4)
    ax.spines['left'].set_bounds(0, 1)
    ax.tick_params(length=2.5, pad=3)
    if not show_y:
        ax.tick_params(axis='y', left=False, labelleft=False)
        ax.spines['left'].set_visible(False)
    ax.xaxis.labelpad = 5
    ax.yaxis.labelpad = 6
    for y in [0, .5]:
        ax.hlines(y, 0, 4.08, color='#E5E7EA', linewidth=.5, zorder=0)
    ax.hlines(1, 0, 4.08, color='#ADB3BB', linewidth=.65, linestyle=(0, (3, 3)), zorder=0)
    for arm in ['Yp', 'Yq']:
            d = summary.loc[(summary.condition == condition) & (summary.arm == arm)].copy()
            d['x'] = d.passage.map(xmap)
            d = d.sort_values('x')
            marker = 'o' if arm == 'Yp' else '^'
            style = '-' if arm == 'Yp' else (0, (3.2, 2))
            ax.plot(d.x, d[VALUE], color=colour[condition], linestyle=style, linewidth=1.15, zorder=2)
            for passage, g in wells.loc[(wells.condition == condition) & (wells.arm == arm)].groupby('passage'):
                # Small, deterministic horizontal offsets reveal coincident well observations.
                offsets = np.linspace(-.055, .055, len(g)) if len(g) > 1 else np.array([0.])
                ax.scatter(xmap[passage] + offsets, g[VALUE], marker=marker, s=20,
                           facecolors='white', edgecolors=colour[condition], linewidths=.65, alpha=.6, zorder=3)
            ax.scatter(d.x, d[VALUE], marker=marker, s=17 if arm == 'Yp' else 22,
                       color=colour[condition], edgecolors='white', linewidths=.3, zorder=4)
            last = d.iloc[-1]
            label = 'Yp (AMELY)' if arm == 'Yp' else 'Yq (KDM5D/DDX3Y)'
            below = condition == 'CTRL' and arm == 'Yq'
            ax.annotate(label, (4, last[VALUE]), xytext=(0, -8 if below else 7),
                        textcoords='offset points', color=colour[condition],
                        va='top' if below else 'bottom', ha='right', fontsize=8)
    return ax


def render(wells, summary, manifest):
    REVIEW.mkdir(exist_ok=True)
    family = apply_style(8.5)
    fig = plt.figure(figsize=(5.3, 2.8))
    left = fig.add_axes([.17, .20, .36, .70])
    right = fig.add_axes([.60, .20, .36, .70], sharey=left)
    draw_panel(left, wells, summary, 'CTRL')
    draw_panel(right, wells, summary, 'LOY', show_y=False)
    fig.text(.56, .035, 'Screening passage', ha='center', va='bottom', fontsize=8.5)
    assert left.get_ylim() == right.get_ylim()
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / (STEM + '.' + ext), dpi=600)
    fig.savefig(REVIEW / (STEM + '_preview.png'), dpi=220)
    plt.close(fig)
    # Read the actual vector export for a small, independent bounds/font check.
    doc = fitz.open(OUT / (STEM + '.pdf'))
    page = doc[0]
    spans = [span for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans']]
    outside = [s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
    assert not outside, outside
    assert all('OpenSans' in s['font'] for s in spans)
    assert 'Parental' not in page.get_text()
    assert len(doc) == 1
    pix = page.get_pixmap(matrix=fitz.Matrix(2.5, 2.5), colorspace=fitz.csGRAY)
    pix.save(str(REVIEW / (STEM + '_greyscale.png')))
    doc.close()
    svg = (OUT / (STEM + '.svg')).read_text(encoding='utf-8')
    assert '<text' in svg and 'Parental' not in svg
    for source in manifest['sources']:
        assert sha(Path(source['path'])) == source['sha256']
    result = {
        'figure_role': 'Standalone supplementary figure; numbering pending; later panel integration deferred',
        'stem': STEM, 'size_inches': [5.3, 2.8], 'font_family': family,
        'on_canvas_caption': False,
        'facet_order': ['sgCTRL', 'Y-depleted'], 'shared_y_limits': [-.055, 1.16],
        'curve_labels': ['Yp (AMELY)', 'Yq (KDM5D/DDX3Y)'],
        'plotted_value_column': VALUE, 'reference_copy_numbers': REF_COPIES,
        'regular_font': font_manager.findfont(font_manager.FontProperties(family=family)),
        'main_and_supp1_unchanged': True, 'out_of_bounds_text': outside,
        'svg_editable_text': True, 'png_pixels': list(Image.open(OUT / (STEM + '.png')).size),
        'source_hashes_unchanged': True,
        'exports': [{'path': str(OUT / (STEM + '.' + ext)), 'sha256': sha(OUT / (STEM + '.' + ext))} for ext in ['pdf', 'svg', 'png']],
    }
    (DATA / 'render_validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({'figure': str(OUT / (STEM + '.pdf')), 'plotted_wells': len(wells), 'summary_points': len(summary), 'validation': manifest['concentration_checks']}, indent=2))


if __name__ == '__main__':
    wells, summary, manifest = prepare()
    render(wells, summary, manifest)
