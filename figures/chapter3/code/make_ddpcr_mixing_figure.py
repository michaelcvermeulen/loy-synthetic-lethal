"""Plot the supplied attempt3 cell-mixture table with the thesis style."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager
import numpy as np

THESIS = Path(__file__).resolve().parents[3]
CHAPTER = THESIS / 'figures/chapter3'
STYLE = THESIS / 'skills/thesis-figure-style'
sys.path.insert(0, str(STYLE / 'scripts'))
from thesis_style import apply_style

STEM = 'ch03_main_ddpcr_mixing_validation'
COLOUR = '#00A087'
WIDTH_MM, HEIGHT_MM = 183, 94
DNA_INPUT_NG = {'50uL': 50, '100uL': 100}
SERIES_LABELS = {key: f'{value} ng DNA' for key, value in DNA_INPUT_NG.items()}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=CHAPTER)
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    source = CHAPTER / 'source_data/ddpcr_mixing_attempt3'
    input_path = source / 'ratios.csv'
    with input_path.open(encoding='utf-8-sig', newline='') as handle:
        wide = list(csv.DictReader(handle))
    records = []
    for row in wide:
        for replicate in ('Rep1', 'Rep2', 'Rep3'):
            raw = row[replicate]
            records.append({'series': row['Conc'],
                            'dna_input_ng': DNA_INPUT_NG[row['Conc']],
                            'expected': float(row['Expected']),
                            'replicate': replicate,
                            'observed': '' if raw in ('NA', '') else float(raw),
                            'status': 'missing' if raw in ('NA', '') else 'recorded'})
    assert len(records) == 30 and sum(r['status'] == 'recorded' for r in records) == 24
    write_csv(source / 'measurements_long.csv', records)

    stats, summaries = {}, []
    for series in ('50uL', '100uL'):
        valid = [r for r in records if r['series'] == series and r['status'] == 'recorded']
        x = np.array([r['expected'] for r in valid])
        y = np.array([r['observed'] for r in valid])
        slope, intercept = np.polyfit(x, y, 1)
        correlation = np.corrcoef(x, y)[0, 1]
        stats[series] = {'dna_input_ng': DNA_INPUT_NG[series],
                         'n_measurements': len(valid), 'n_nominal_mixtures': 5,
                         'pearson_r': float(correlation), 'ols_slope': float(slope),
                         'ols_intercept': float(intercept), 'inferential_test': None}
        # Match the descriptive calculations retained during source review.
        retained_r = {'50uL': 0.9929547693942286, '100uL': 0.9882356056062749}
        assert abs(correlation - retained_r[series]) < 1e-10
        for expected in (0, .25, .5, .75, 1):
            vals = y[x == expected]
            summaries.append({'series': series,
                              'dna_input_ng': DNA_INPUT_NG[series],
                              'expected': expected, 'n': len(vals),
                              'mean': float(vals.mean()),
                              'sd': float(vals.std(ddof=1)) if len(vals) > 1 else ''})
    write_csv(source / 'mixture_summary.csv', summaries)
    (source / 'descriptive_statistics.json').write_text(json.dumps(stats, indent=2) + '\n')

    family = apply_style(base_size=9)
    font_path = font_manager.findfont(font_manager.FontProperties(family=family))
    assert Path(font_path).resolve() == (STYLE / 'assets/fonts/OpenSans-Regular.ttf').resolve()
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4),
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=.09, right=.98, bottom=.245, top=.89, wspace=.23)
    offsets = {'Rep1': -.012, 'Rep2': 0., 'Rep3': .012}
    for panel, ax, series in zip(('a', 'b'), axes, ('50uL', '100uL')):
        valid = [r for r in records if r['series'] == series and r['status'] == 'recorded']
        s = stats[series]
        grid = np.linspace(0, 1, 200)
        ax.plot(grid, grid, color='#858585', linestyle=(0, (3, 3)), linewidth=.85,
                zorder=1)
        ax.plot(grid, s['ols_slope'] * grid + s['ols_intercept'],
                color=COLOUR, linewidth=1.1, zorder=2)
        ax.scatter([r['expected'] + offsets[r['replicate']] for r in valid],
                   [r['observed'] for r in valid], s=25, color=COLOUR,
                   edgecolors='white', linewidths=.45, zorder=3)
        ax.set(xlim=(-.055, 1.06), ylim=(-.055, 1.10),
               xticks=[0, .25, .5, .75, 1], yticks=[0, .25, .5, .75, 1])
        ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'])
        ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '1'])
        ax.set_box_aspect(1)
        ax.set_title(SERIES_LABELS[series], pad=9, fontweight='semibold')
        ax.set_xlabel('Expected Y-positive H2452 fraction', labelpad=7)
        ax.text(-.16, 1.075, panel, transform=ax.transAxes, fontsize=12,
                fontweight='semibold', va='center')
        ax.text(.055, .96, f"Pearson r = {s['pearson_r']:.3f}\n"
                f"n = {s['n_measurements']} measurements",
                transform=ax.transAxes, va='top', fontsize=8.5, linespacing=1.5)
    axes[0].set_ylabel('Observed AMELY/AMELX ratio', labelpad=8)
    handles = [Line2D([], [], color=COLOUR, marker='o', linestyle='None', markersize=4,
                       label='Technical measurement'),
               Line2D([], [], color=COLOUR, linewidth=1.1, label='Linear fit'),
               Line2D([], [], color='#858585', linestyle=(0, (3, 3)), linewidth=.85,
                       label='Expected = observed')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .015),
               ncol=3, columnspacing=1.8, handlelength=2, fontsize=8.2)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for ext in ('pdf', 'svg', 'png'):
        path = args.output_dir / f'{STEM}.{ext}'
        fig.savefig(path, dpi=args.dpi)
        outputs.append(path)
    preview = args.output_dir / 'previews' / f'{STEM}.png'
    preview.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(preview, dpi=180)
    outputs.append(preview)
    plt.close(fig)
    provenance = {
        'source_original': 'A:/Dropbox/LOY_paralogs/data/ddPCR/control_test/attempt3/ratios.csv',
        'input_sha256': sha256(input_path), 'script_sha256': sha256(Path(__file__)),
        'style_helper_sha256': sha256(STYLE / 'scripts/thesis_style.py'),
        'font_family': family, 'font_file': str(font_path), 'font_sha256': sha256(Path(font_path)),
        'width_mm': WIDTH_MM, 'height_mm': HEIGHT_MM, 'png_dpi': args.dpi,
        'series_order': ['50uL', '100uL'], 'series_display_labels': SERIES_LABELS,
        'dna_input_ng': DNA_INPUT_NG,
        'dna_input_provenance': 'Author correction, 23 September 2026: total DNA input per reaction, not concentration',
        'measurement_colour': COLOUR,
        'colour_role': 'Measurements in both panels; not Y-status encoding',
        'horizontal_offsets_display_only': offsets,
        'statistics': stats, 'missing_values': 'Six NA entries omitted; never imputed as zero',
        'statistical_scope': 'Descriptive only; technical wells are not independent mixture preparations',
        'axis_definitions': {'expected': 'Y-positive H2452 fraction', 'observed': 'AMELY/AMELX concentration ratio'},
        'axis_confirmation': 'Author confirmed 23 September 2026',
        'unresolved': [],
        'software': {'python': sys.version, 'numpy': np.__version__, 'matplotlib': matplotlib.__version__},
        'outputs': {str(p.relative_to(args.output_dir)): sha256(p) for p in outputs}}
    (source / 'render_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'outputs': [str(p) for p in outputs], 'statistics': stats}, indent=2))


if __name__ == '__main__':
    main()
