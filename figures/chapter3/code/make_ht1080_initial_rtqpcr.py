"""Rebuild initial HT1080 RT-qPCR figures from retained instrument tables.

Attempts are repeat assays of the same edits (author confirmed 2026-09-17).
Normalise within each run; never pool wells as biological replicates.
"""
from pathlib import Path
import json
import sys
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'figures/chapter3'
DATA = OUT / 'source_data/ht1080_initial_rtqpcr'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
FONT = apply_style(9.5)
ORDER = ['WT', 'CN', 'CP', 'Q53', 'P47', 'P51']
COLOURS = {1: '#008B8B', 2: '#8064A2'}
MARKERS = {1: 'o', 2: '^'}
records, summaries, qc = [], [], []
for run in (1, 2):
    frame = pd.read_csv(DATA / f'RT_QPCR_HT1080_ATTEMPT{run}.csv')
    frame['source_row'] = np.arange(len(frame)) + 2
    frame['numeric_ct'] = pd.to_numeric(frame['CT'], errors='coerce')
    frame = frame[frame['Sample Name'].str.startswith('HT1080 ', na=False)].copy()
    means = frame.groupby(['Sample Name', 'Target Name'])['Ct Mean'].first()
    for (sample, gene), group in frame.groupby(['Sample Name', 'Target Name']):
        assert group['Ct Mean'].dropna().nunique() <= 1
        gapdh = float(means[(sample, 'GAPDH')])
        wt_delta = float(means[('HT1080 WT', gene)] - means[('HT1080 WT', 'GAPDH')])
        reported_mean = float(means[(sample, gene)])
        qc.append(dict(run=run, sample=sample, gene=gene,
                       n_wells=len(group), n_numeric_ct=int(group.numeric_ct.notna().sum()),
                       reported_mean_ct=reported_mean, recomputed_mean_ct=group.numeric_ct.mean(),
                       mean_difference=reported_mean-group.numeric_ct.mean()))
        for _, row in group.iterrows():
            ddct = row.numeric_ct - gapdh - wt_delta
            records.append(dict(run=run, sample=sample, guide=sample.replace('HT1080 ', ''),
                                gene=gene, source_row=int(row.source_row), ct=row.numeric_ct,
                                reported_mean_ct=reported_mean, gapdh_mean_ct=gapdh,
                                wt_delta_ct=wt_delta, delta_delta_ct=ddct,
                                expression_ratio=2.0**(-ddct)))
        summaries.append(dict(run=run, sample=sample, guide=sample.replace('HT1080 ', ''),
                              gene=gene, n_wells=len(group), n_numeric_ct=int(group.numeric_ct.notna().sum()),
                              expression_ratio=2.0**(-(reported_mean-gapdh-wt_delta))))

WELLS = pd.DataFrame(records)
SUMMARY = pd.DataFrame(summaries)
WELLS.to_csv(DATA / 'normalised_wells.csv', index=False)
SUMMARY.to_csv(DATA / 'normalised_run_summaries.csv', index=False)
pd.DataFrame(qc).to_csv(DATA / 'ct_summary_checks.csv', index=False)
assert len(WELLS) == 252
assert np.allclose(SUMMARY[(SUMMARY.guide == 'WT') & (SUMMARY.gene != 'NLGN4Y')].expression_ratio, 1)
assert SUMMARY[(SUMMARY.guide == 'CN') & (SUMMARY.run == 2)].empty
assert WELLS[WELLS.gene == 'NLGN4Y'].expression_ratio.isna().all()
assert pd.DataFrame(qc).mean_difference.abs().max() < .0011

def draw(genes, stem, width_mm, supplementary=False):
    height_mm = 100 if supplementary else 140
    fig, axes = plt.subplots(1 if supplementary else 2, 2,
                             figsize=(width_mm/25.4, height_mm/25.4), sharey=True)
    axes = np.atleast_1d(axes).ravel()
    fig.subplots_adjust(left=.10, right=.98, bottom=.25 if supplementary else .22,
                        top=.87 if supplementary else .91, wspace=.20, hspace=.48)
    for panel, (ax, gene) in enumerate(zip(axes, genes)):
        for level in [0.5, 1.5]:
            ax.axhline(level, color='#ECECEC', lw=.6, zorder=0)
        ax.axhline(1, color='#919191', ls=(0, (3, 3)), lw=.85, zorder=0)
        for x, guide in enumerate(ORDER):
            for run in (1, 2):
                points = WELLS[(WELLS.gene == gene) & (WELLS.guide == guide) & (WELLS.run == run)]
                if points.empty:
                    continue
                offset = 0 if supplementary else (-.16 if run == 1 else .16)
                values = points.expression_ratio.to_numpy()
                xx = x + offset + np.linspace(-.075, .075, len(values))
                ax.scatter(xx, values, s=22, marker=MARKERS[run],
                           facecolor='white', edgecolor=COLOURS[run], lw=1, zorder=3)
                value = SUMMARY[(SUMMARY.gene == gene) & (SUMMARY.guide == guide) & (SUMMARY.run == run)].expression_ratio.iloc[0]
                ax.plot([x+offset-.12, x+offset+.12], [value, value],
                        color=COLOURS[run], lw=2.0, zorder=4)
        ax.set_title(gene + (' (X/Y)' if gene == 'CD99' else ''), pad=11, fontsize=10.5)
        ax.text(-.13, 1.07, chr(97+panel), transform=ax.transAxes,
                weight='semibold', fontsize=12, ha='left')
        ax.set_xlim(-.55, 5.55)
        ax.set_xticks(range(6), ORDER)
        ax.set_ylim(-.035, 1.65)
        ax.set_yticks([0, .5, 1, 1.5], ['0', '0.5', '1', '1.5'])
        ax.tick_params(axis='y', length=3, pad=4)
        ax.tick_params(axis='x', labelsize=9, length=0, pad=7)
    fig.text(.024, .53 if not supplementary else .56, 'Relative expression (WT = 1)',
             rotation=90, ha='center', va='center', fontsize=10)
    fig.text(.54, .155 if not supplementary else .15, 'HT1080 population', ha='center', fontsize=10)
    handles = [Line2D([], [], marker=MARKERS[r], color=COLOURS[r], markerfacecolor='white',
                      linestyle='none', markersize=5, label=f'Run {r}') for r in ([1] if supplementary else [1, 2])]
    handles += [Line2D([], [], color='#444444', lw=2, label='Ratio from mean Ct')]
    if supplementary:
        fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.53,.053), ncol=2,
                   handletextpad=.5, columnspacing=1.5, fontsize=9)
        fig.text(.53, .023, 'Three technical wells per target/population', ha='center', fontsize=8)
    else:
        handles += [Line2D([], [], color='#919191', lw=.85, linestyle='--', label='WT reference')]
        fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.54,.072), ncol=4,
                   handlelength=1.5, handletextpad=.5, columnspacing=1.4, fontsize=9)
        fig.text(.54, .035, 'Repeat assays of the same edits; 3 wells per run; CN available in Run 1 only',
                 ha='center', fontsize=8.5, color='#555555')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / f'{stem}.{ext}', dpi=600)
    (OUT / 'previews').mkdir(exist_ok=True)
    fig.savefig(OUT / 'previews' / f'{stem}.png', dpi=180)
    plt.close(fig)

main = ['EIF1AY', 'KDM5D', 'RPS4Y1', 'UTY']
assert len(WELLS[WELLS.gene.isin(main)]) == 132
draw(main, 'ch03_main_ht1080_initial_guides_rtqpcr', 183)
draw(['USP9Y', 'ZFY'], 'ch03_supp_ht1080_initial_rtqpcr_additional_targets', 120, True)
provenance = dict(font=FONT, main_size_mm=[183,140], supplementary_size_mm=[120,100], dpi=600,
                  run_colours=COLOURS, run_markers=MARKERS, guide_order=ORDER,
                  main_genes=main, technical_wells_main=132, technical_wells_supplement=36,
                  formula='2^-[(Ct_well - reported_mean_Ct_GAPDH_same_sample_run) - (reported_mean_Ct_target_WT_same_run - reported_mean_Ct_GAPDH_WT_same_run)]',
                  summary='Same formula substituting reported target mean Ct for well Ct; no averaging across runs.',
                  independence='Runs are repeat assays of the same edited populations; author confirmed 2026-09-17. No inferential tests.',
                  omitted='CD99 removed from figure at author request; retained in source and derived tables. NLGN4Y cannot be normalised because all WT wells are undetermined; GAPDH is the reference, not an independent endpoint. NTC retained in copied input only.',
                  historical_change='Correct explicit within-run WT matching instead of positional cross-run assignment in historical combine script. Existing PDFs also used a third missing input; these figures use only the two available HT1080 tables.',
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(DATA / 'render_provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
print('Rendered main and supplementary RT-qPCR figures; numerical checks passed.')
