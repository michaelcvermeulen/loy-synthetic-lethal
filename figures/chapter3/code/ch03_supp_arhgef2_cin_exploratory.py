"""Exploratory display of the saved ARHGEF2/CIN follow-up models.

Run the companion R script with Rscript --vanilla first. The two reconstructed
coefficients and P values are checked against the saved HPC results in R.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

root = Path(__file__).resolve().parents[3]
out = root / 'figures/chapter3'
data = out / 'source_data/arhgef2_cin_exploratory'
points = pd.read_csv(data / 'points.csv')
curves = pd.read_csv(data / 'curves.csv')
contrasts = pd.read_csv(data / 'y_status_contrast.csv')
fits = pd.read_csv(data / 'fits.csv').set_index('assay')
assert fits.verified_against_saved.all()
assert dict(points.groupby('assay').size()) == {'CRISPR': 473, 'RNAi': 244}
for font in ['Regular', 'SemiBold']:
    font_manager.fontManager.addfont(str(root / f'skills/thesis-figure-style/assets/fonts/OpenSans-{font}.ttf'))
plt.rcParams.update({
    'font.family': 'Open Sans', 'font.size': 9, 'axes.labelsize': 9,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .65, 'xtick.major.width': .65, 'ytick.major.width': .65,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'text.color': '#222222',
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Open Sans',
    'mathtext.it': 'Open Sans', 'mathtext.bf': 'Open Sans',
    'axes.labelcolor': '#222222', 'xtick.color': '#333333', 'ytick.color': '#333333',
})
fig = plt.figure(figsize=(9.2, 6.6))
colours = {'Y-positive': '#3C5488', 'Y-negative': '#E64B35'}
for i, assay in enumerate(['CRISPR', 'RNAi']):
    ax = fig.add_axes([.10 + .465 * i, .505, .395, .395])
    stats = fits.loc[assay]
    for group, colour in colours.items():
        p = points[(points.assay == assay) & (points.y_status == group)]
        c = curves[(curves.assay == assay) & (curves.y_status == group)]
        ax.fill_between(c.cin, c.lower, c.upper, color=colour, alpha=.12, lw=0)
        ax.scatter(p.cin, p.adjusted_dependency, s=12, c=colour, alpha=.42,
                   edgecolors='none', label=f'{group} (n = {len(p)})', rasterized=False)
        ax.plot(c.cin, c.estimate, c=colour, lw=1.8,
                ls='-' if group == 'Y-positive' else (0, (4, 1.5)))
    ax.set(xlim=(0, .88), ylim=(-3.15, 3.3), xticks=[0, .2, .4, .6, .8],
           yticks=[-3, -2, -1, 0, 1, 2, 3], xlabel='Chromosomal instability (CIN score)')
    if i == 0:
        ax.set_ylabel('Adjusted ARHGEF2 gene effect\n(rank-normalised)')
    ax.text(0, 1.045, assay, transform=ax.transAxes, weight='semibold', fontsize=10)
    p_label = r'$1.03\times10^{-6}$' if assay == 'CRISPR' else '0.89'
    ax.text(.98, .97, f'LOY × CIN interaction P = {p_label}',
            transform=ax.transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=2))
    ax.legend(loc='lower right', frameon=True, facecolor='white', edgecolor='none',
              framealpha=.9, fontsize=8, markerscale=1.5,
              handletextpad=.25, borderpad=.2, labelspacing=.25)
    fig.text(.025 + .465 * i, .932, chr(97 + i), fontsize=13, weight='semibold')

# The contrast makes the interaction directly visible without dichotomising CIN.
ax = fig.add_axes([.10, .12, .86, .255])
ax.axhline(0, color='#777777', lw=.8, ls=(0, (3, 3)))
assay_colours = {'CRISPR': '#6B4C9A', 'RNAi': '#555555'}
for assay, colour in assay_colours.items():
    c = contrasts[contrasts.assay == assay]
    ax.fill_between(c.cin, c.lower, c.upper, color=colour, alpha=.13, lw=0)
    ax.plot(c.cin, c.estimate, color=colour, lw=2,
            ls='-' if assay == 'CRISPR' else (0, (4, 2)), label=assay)
ax.set(xlim=(0, .88), ylim=(-1.85, 1.55), xticks=[0, .2, .4, .6, .8],
       yticks=[-1.5, -1, -.5, 0, .5, 1, 1.5],
       xlabel='Chromosomal instability (CIN score)',
       ylabel='Predicted gene-effect difference\n(Y-negative − Y-positive)')
ax.text(.012, .97, 'Weaker dependency in LOY', transform=ax.transAxes,
        va='top', fontsize=8, color='#555555')
ax.text(.988, .035, 'Stronger dependency in LOY', transform=ax.transAxes,
        ha='right', va='bottom', fontsize=8, color='#555555')
ax.legend(loc='upper center', bbox_to_anchor=(.51, 1.22), ncol=2, frameon=False,
          handlelength=2.5, columnspacing=2, fontsize=9)
fig.text(.025, .407, 'c', fontsize=13, weight='semibold')
fig.text(.10, .025, 'Points: individual models after covariate adjustment. Bands: pointwise 95% confidence intervals.',
         fontsize=8, color='#444444')
stem = 'ch03_supp_arhgef2_cin_exploratory'
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(out / f'{stem}.{ext}', dpi=250, facecolor='white')
plt.close(fig)

export = root / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07'
source_names = ['cohort.csv', 'context_matrix.csv', 'robustness_context_gene_hits.csv',
                'matrices/crispr_effect.rds', 'matrices/rnai.rds']
source_hashes = {name: hashlib.sha256((export / 'derived/loy_reanalysis_v1' / name).read_bytes()).hexdigest()
                 for name in source_names}
metadata = {
    'status': 'Exploratory visualization requested for interest; not selected or inserted into Results',
    'source_export': str(export), 'source_sha256': source_hashes,
    'model': 'Stage 19 rank-normalised ARHGEF2 outcome ~ lineage + Y status + standardised CIN + Y:CIN; CRISPR also adjusts LibraryKY and LibraryHumagneCD',
    'verification': 'Reconstructed OLS interaction estimates and P values match saved stage-19 coefficients to 1e-9 and log(P ratio) to 1e-7',
    'points': 'Rank-normalised measurements minus centred fitted lineage/library contributions; retain CIN, Y status and interaction',
    'curves': 'Predictions at mean nuisance covariates within each assay; group curves limited to observed CIN ranges, difference limited to overlapping support',
    'uncertainty': 'Classical model-based pointwise 95% t intervals; not simultaneous bands',
    'p_values': 'Nominal OLS follow-up P values. Original moderated discovery test P=1.56361807452549e-6, class FDR=0.0380834818231428, global FDR=0.303762580874048',
    'limits': 'Assay cohorts may overlap and do not establish independent biological replication; ranking scales and covariate distributions are assay-specific. CIN remains continuous.',
    'font': 'Bundled Open Sans',
}
(data / 'provenance.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
(data / 'README.md').write_text('''# Exploratory ARHGEF2–CIN visualization

Requested for interest on 7 October 2026. No figure number assigned and no Results text changed.

The two upper panels show individual cell-line ARHGEF2 measurements after the rank transformation and covariate adjustment used in the original follow-up. Lower values mean stronger dependency. Blue denotes Y-positive models and red Y-negative models. Lines and shaded bands show fitted means and pointwise 95% confidence intervals. CRISPR adjusts for lineage and screening library; RNAi adjusts for lineage. The lower panel shows the predicted Y-negative minus Y-positive difference over the CIN range observed in both Y groups for each assay. Negative differences indicate stronger dependency in LOY.

The R script reconstructs just the two saved stage-19 models and checks both interaction estimates and P values before writing display tables. These nominal follow-up P values differ slightly from the original moderated discovery test (P = 1.56e-6, within-class FDR = 0.038, global FDR = 0.304). RNAi reproduces the saved P = 0.887993. The CRISPR association is therefore a context-screen lead, not a genome-wide confirmed dependency. The plot does not dichotomise CIN or establish a causal effect of LOY or CIN.

Reproduce from the thesis root:

    Rscript --vanilla figures/chapter3/code/ch03_supp_arhgef2_cin_exploratory.R .
    python figures/chapter3/code/ch03_supp_arhgef2_cin_exploratory.py

Source hashes, exact transformations, interval definitions and limitations are recorded in provenance.json. Original source tables and matrices are unchanged.
''', encoding='utf-8')
print(out / f'{stem}.png')
