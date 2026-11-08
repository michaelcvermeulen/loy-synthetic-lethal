"""DU145 ddPCR by assay/run; preserve historical LNCaP labels in provenance."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

THESIS = Path(__file__).resolve().parents[3]
CHAPTER = THESIS / 'figures/chapter3'
STYLE = THESIS / 'skills/thesis-figure-style'
sys.path.insert(0, str(STYLE / 'scripts'))
from thesis_style import apply_style

STEM = 'ch03_main_du145_ddpcr_assay_comparison'
RUNS = [
    ('2025-11-07', '7 November 2025', 'AMELX Craig_NOV7_2025_2025-11-07-15-53', 'AMELY', 'AMELX', 'AMELY', 'AMELX'),
    ('2025-11-13', '13 November 2025', 'AMELX Craig_NOV_13_2025-11-12-22-50', 'DDX3Y', 'RNaseP', 'DDX3Y', 'RNaseP'),
    # The analysis Rmd explicitly maps raw AMELX to AMELY and AMELX_2 to AMELX.
    ('2026-02-10', '10 February 2026', 'AMELX Craig_FEB_10_2026_2026-02-10-16-54', 'AMELX', 'AMELX_2', 'AMELY', 'AMELX'),
    ('2026-02-25', '24\u201325 February 2026', 'AMELX Craig_FEB_24_2026_2026-02-25-17-16', 'KDM5D', 'RNaseP', 'KDM5D', 'RNaseP'),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-root', type=Path, default=Path('A:/Dropbox/LOY_paralogs/data/ddPCR'))
    parser.add_argument('--output-dir', type=Path, default=CHAPTER)
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = CHAPTER / 'source_data/du145_ddpcr_assay_comparison'
    snapshots = source / 'raw_exports'
    snapshots.mkdir(parents=True, exist_ok=True)
    family = apply_style(8.5)
    font = Path(font_manager.findfont(font_manager.FontProperties(family=family)))
    assert font.resolve() == (STYLE / 'assets/fonts/OpenSans-Regular.ttf').resolve()
    records, manifest = [], []
    for date, display_date, folder, raw_target, raw_ref, target, reference in RUNS:
        path = args.raw_root / folder / 'data_table.csv'
        data = pd.read_csv(path)
        concentration = next(c for c in data if c.startswith('Conc(copies/'))
        mask = data['Sample description 1'].str.contains('LNCAP', case=False, na=False)
        for _, row in data[mask & data.Target.eq(raw_target)].iterrows():
            label = row['Sample description 1']
            upper = label.upper()
            condition = ('Parental' if 'PARENTAL' in upper else 'sgCTRL' if 'CTRL' in upper
                         else 'CP + Q53 + P47' if upper.endswith('_MIX') else 'Y-depleted')
            round_count = 1 if condition == 'CP + Q53 + P47' else None
            if date in ('2026-02-10', '2026-02-25') and condition == 'Y-depleted':
                condition = 'CP + Q53 + P47 2x'
                round_count = 2
            ref_rows = data[data.Well.eq(row.Well) & data.Target.eq(raw_ref)]
            assert len(ref_rows) == 1
            ref = ref_rows.iloc[0]
            yc = pd.to_numeric(row[concentration], errors='coerce')
            xc = pd.to_numeric(ref[concentration], errors='coerce')
            valid = pd.notna(yc) and pd.notna(xc) and xc > 0
            ratio = yc / xc if valid else np.nan
            if valid:
                assert np.isclose(ratio, float(row['Ratio']))
            records.append(dict(cell_line='DU145', historical_cell_line='LNCaP', date=date,
                display_date=display_date, assay=target+'/'+reference, condition=condition,
                source_label=label, export_well=row.Well, raw_target=raw_target,
                raw_reference=raw_ref, target=target, reference=reference,
                target_concentration=yc, reference_concentration=xc,
                concentration_ratio=ratio, instrument_ratio=row['Ratio'],
                accepted_droplets=row['Accepted Droplets'], target_status=row['Status'],
                target_status_reason=row['Status Reason'], reference_status=ref['Status'],
                included=valid, exclusion_reason='' if valid else 'Instrument No Call',
                nucleofection_count=round_count, DNA_sample_id=None, source_path=str(path)))
        snapshot = snapshots / (date+'_data_table.csv')
        shutil.copy2(path, snapshot)
        manifest.append(dict(source_path=str(path), snapshot=str(snapshot.relative_to(source)),
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    all_data = pd.DataFrame(records)
    measurements = all_data[all_data.included].copy()
    assert len(measurements) == 11 and len(all_data) == 12
    assert not measurements.duplicated(['date', 'export_well', 'assay']).any()
    all_data.to_csv(source / 'all_assessed_measurements.csv', index=False)
    measurements.to_csv(source / 'main_measurements.csv', index=False)
    all_data[~all_data.included].to_csv(source / 'excluded_measurements.csv', index=False)
    (source/'raw_source_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')

    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 110/25.4), sharey=True)
    fig.subplots_adjust(left=.095, right=.988, bottom=.255, top=.79, wspace=.26)
    summaries = []
    colours = {'CP + Q53 + P47 2x': '#E64B35', 'CP + Q53 + P47': '#E64B35',
               'sgCTRL': '#3C5488', 'Parental': '#888888'}
    panels = [
        ('AMELY / AMELX', '1\u00d7 population', '2\u00d7 population', '2025-11-07', '2026-02-10'),
        ('Y-linked target / RNaseP', 'DDX3Y', 'KDM5D', '2025-11-13', '2026-02-25'),
    ]
    for ax, panel, letter in zip(axes, panels, 'ab'):
        title, first_heading, second_heading, first_date, second_date = panel
        groups = [(first_date,'CP + Q53 + P47'), (first_date,'sgCTRL'), (first_date,'Parental'),
                  (second_date,'CP + Q53 + P47 2x'), (second_date,'sgCTRL')]
        positions = [0, 1, 2, 3.35, 4.35]
        for x, (date, condition) in zip(positions, groups):
            subset = measurements[measurements.date.eq(date) & measurements.condition.eq(condition)]
            values = subset.concentration_ratio.to_numpy()
            assert len(values)
            med = float(np.median(values))
            q1, q3 = np.quantile(values, [.25, .75])
            ax.bar(x, med, width=.65, color=colours[condition], edgecolor=colours[condition], zorder=2)
            if condition == 'CP + Q53 + P47 2x':
                ax.bar(x, med, width=.65, facecolor='none', edgecolor='#87301E', linewidth=.7, hatch='///', zorder=3)
            if len(values) > 1:
                ax.errorbar(x, med, yerr=[[med-q1], [q3-med]], fmt='none', ecolor='#333333',
                            elinewidth=.8, capsize=2, zorder=5)
            offsets = np.linspace(-.13, .13, len(values)) if len(values)>1 else [0]
            ax.scatter(x+np.asarray(offsets), values, s=15, facecolor='white', edgecolor='#222222', linewidth=.65, zorder=6)
            y = values.max()+.029
            ax.text(x, y, f'{med:.3f}', ha='center', va='bottom', fontsize=7.5, fontweight='semibold')
            ax.text(x, y+.039, f'n={len(values)}', ha='center', va='bottom', fontsize=6.8, color='#666666')
            summaries.append(dict(panel=letter, date=date, assay=subset.assay.iloc[0], condition=condition,
                                  n_measurements=len(values), median=med, q1=q1, q3=q3))
        labels = ['CP + Q53\n+ P47\n1\u00d7', 'sgCTRL', 'Parental', 'CP + Q53\n+ P47\n2\u00d7', 'sgCTRL']
        ax.set_xticks(positions, labels)
        ax.tick_params(axis='x', length=0, pad=7, labelsize=7.2)
        ax.tick_params(axis='y', labelleft=True)
        ax.set_xlim(-.6, 4.95)
        ax.set_ylim(0, .68)
        ax.set_yticks([0, .2, .4, .6])
        ax.set_title(title, fontsize=10, pad=34)
        ax.text(1, 1.035, first_heading, transform=ax.get_xaxis_transform(), ha='center', fontsize=8, color='#555555')
        ax.text(3.85, 1.035, second_heading, transform=ax.get_xaxis_transform(), ha='center', fontsize=8, color='#555555')
        ax.axvline(2.675, color='#DDDDDD', linewidth=.7)
        ax.text(-.14, 1.19, letter, transform=ax.transAxes, fontsize=12, fontweight='semibold')
        ax.grid(axis='y', color='#ECECEC', linewidth=.45)
        ax.set_axisbelow(True)
    axes[0].set_ylabel('Target/reference concentration ratio', labelpad=9)
    fig.text(.545, .045, '1\u00d7 / 2\u00d7 = one / two successive nucleofections', ha='center', fontsize=8, color='#444444')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(args.output_dir/(STEM+'.'+ext), dpi=args.dpi, facecolor='white')
    fig.savefig(args.output_dir/(STEM+'_preview.png'), dpi=170, facecolor='white')
    plt.close(fig)
    pd.DataFrame(summaries).to_csv(source/'main_group_summary.csv', index=False)
    validation = dict(font=str(font), dimensions_mm=[183,110], dpi=args.dpi,
                      included_measurements=11, excluded_no_call=1,
                      source_concentration_recheck='passed', normalization='unscaled concentration ratios',
                      reference_copy_numbers='not confirmed for experimental stock',
                      nucleofection_counts='February AMELY and KDM5D edited samples: CP+Q53+P47 twice, author-confirmed; KDM5D used instead of DDX3Y; November AMELY and DDX3Y edited samples: one nucleofection, author-confirmed 23 September 2026',
                      historical_alias='LNCaP displayed as DU145 per author and chapter notes')
    (source/'figure_validation.json').write_text(json.dumps(validation, indent=2)+'\n', encoding='utf-8')
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == '__main__':
    main()
