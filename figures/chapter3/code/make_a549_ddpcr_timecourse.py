"""Restyle the original A549 DNA-collection time course on the adopted dosage scale."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import dates as mdates, font_manager
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'figures/chapter3'
DATA = OUT / 'source_data/a549_ddpcr_timecourse'
INPUT = OUT / 'source_data/a549_ddpcr_nucleofection_rounds'
STYLE = ROOT / 'skills/thesis-figure-style'
STEM = 'ch03_main_a549_ddpcr_timecourse'
sys.path.insert(0, str(STYLE / 'scripts'))
from thesis_style import apply_style

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    DATA.mkdir(parents=True, exist_ok=True)
    family = apply_style(9)
    actual_font = Path(font_manager.findfont(font_manager.FontProperties(family=family)))
    assert actual_font.resolve() == (STYLE / 'assets/fonts/OpenSans-Regular.ttf').resolve()
    inputs = ['main_measurements.csv', 'a549_master_metadata.csv', 'excluded_measurements.csv']
    for name in inputs:
        if not (DATA / name).exists():
            shutil.copy2(INPUT / name, DATA / name)
    m = pd.read_csv(DATA / 'main_measurements.csv')
    meta = pd.read_csv(DATA / 'a549_master_metadata.csv')
    dates = meta[['Sample', 'DNA_date']].drop_duplicates()
    assert not dates.Sample.duplicated().any()
    m = m.merge(dates, left_on='sample_id', right_on='Sample', validate='many_to_one')
    m['DNA_date'] = pd.to_datetime(m.DNA_date)
    # Preserve the original January-February time-course scope, using collection
    # dates rather than repeat assay dates. March samples lack a parallel series.
    m = m[m.condition.isin(['CTRL2', 'CP + Q53']) &
          m.DNA_date.between('2025-01-24', '2025-02-19')].copy()
    assert m.nucleofection_count.eq(1).all()
    assert not ((m.sample_id == 'I1') & (m.target == 'AMELY')).any()
    assert len(m) == 37
    assert np.allclose(m.estimated_target_copies_per_cell,
                       m.concentration_ratio * m.reference_copy_number)
    keys = ['target', 'condition', 'DNA_date']
    stats = m.groupby(keys).estimated_target_copies_per_cell.agg(
        median='median', q1=lambda x: x.quantile(.25),
        q3=lambda x: x.quantile(.75), n='size').reset_index()
    m.to_csv(DATA / 'timecourse_measurements.csv', index=False)
    stats.to_csv(DATA / 'timecourse_summary.csv', index=False)
    dates.to_csv(DATA / 'dna_collection_date_mapping.csv', index=False)

    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 105/25.4), sharey=True)
    fig.subplots_adjust(left=.105, right=.985, bottom=.28, top=.80, wspace=.17)
    colours = {'CTRL2': '#3C5488', 'CP + Q53': '#E64B35'}
    markers = {'CTRL2': 'o', 'CP + Q53': 's'}
    ticks = pd.to_datetime(['2025-01-24', '2025-01-31', '2025-02-03', '2025-02-05', '2025-02-19'])
    for ax, target, label in zip(axes, ['AMELY', 'DDX3Y'], ['a', 'b']):
        ax.axhline(1, color='#b1b1b1', lw=.65, ls=(0, (4, 3)), zorder=0)
        for condition in ['CTRL2', 'CP + Q53']:
            s = stats[(stats.target == target) & (stats.condition == condition)].sort_values('DNA_date')
            x = mdates.date2num(s.DNA_date)
            ax.errorbar(x, s['median'], yerr=[s['median']-s.q1, s.q3-s['median']],
                        color=colours[condition], marker=markers[condition], ms=4,
                        lw=1.15, elinewidth=.85, capsize=2.2, zorder=3)
            for date, group in m[(m.target == target) & (m.condition == condition)].groupby('DNA_date'):
                dx = np.linspace(-.32, .32, len(group)) if len(group) > 1 else np.array([0.])
                ax.scatter(mdates.date2num(date)+dx, group.estimated_target_copies_per_cell,
                           marker=markers[condition], s=12, facecolor='white',
                           edgecolor=colours[condition], lw=.65, alpha=.8, zorder=4)
        ax.set_xlim(mdates.date2num(ticks[0])-1.4, mdates.date2num(ticks[-1])+1.4)
        ax.set_ylim(-.025, 1.25)
        ax.set_yticks([0, .25, .5, .75, 1])
        ax.set_xticks(mdates.date2num(ticks))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        plt.setp(ax.get_xticklabels(), rotation=50, ha='right', fontsize=8)
        ax.tick_params(length=3, width=.6)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_title(target, fontweight='semibold', fontsize=10, pad=9)
        ax.text(-.15, 1.065, label, transform=ax.transAxes,
                fontweight='semibold', fontsize=12)
        ax.set_xlabel('DNA collection date (2025)', labelpad=8)
    axes[0].set_ylabel('Estimated Y-linked target copies per cell', labelpad=8)
    handles = [Line2D([0], [0], color=colours[c], marker=markers[c], lw=1.15, ms=4,
                      label=l) for c, l in [('CTRL2', 'sgCTRL'), ('CP + Q53', 'CP + Q53 (1×)')]]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.54, .998),
               ncol=2, frameon=False, columnspacing=2, handlelength=2.1)
    fig.text(.105, .022, 'Source: A549 ddPCR measurements; DNA collected 24 January–19 February 2025.',
             color='#555555', fontsize=7)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{STEM}.{ext}', dpi=600, facecolor='white')
    preview = OUT / 'previews'
    preview.mkdir(exist_ok=True)
    fig.savefig(preview / f'{STEM}.png', dpi=180, facecolor='white')
    plt.close(fig)
    manifest = dict(script_sha256=sha(Path(__file__)),
        inputs={n: sha(DATA/n) for n in inputs}, font=str(actual_font),
        n_measurements=len(m), n_dates=m.DNA_date.nunique(),
        sample_ids=sorted(m.sample_id.unique()),
        x_axis='DNA collection date, not assay date or time since nucleofection',
        summary='median and IQR; individual assay observations; no inferential tests',
        reference_copies={'AMELX': 2, 'RNaseP': 3},
        changes_from_legacy='Use validated instrument concentration ratios and full reference-copy scaling; legacy used positive-droplet ratios and RNaseP copies/2.',
        dimensions_mm=[183, 105], dpi=600,
        exports={ext: sha(OUT/f'{STEM}.{ext}') for ext in ['pdf','svg','png']})
    (DATA/'render_provenance.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(stats.to_string(index=False))

if __name__ == '__main__':
    main()
