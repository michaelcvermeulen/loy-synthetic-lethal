"""Descriptive mean expression reduction, linked to a GRCh38 Y ideogram."""
from pathlib import Path
import json
import shutil
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'figures/chapter3'
DATA = OUT / 'source_data/ht1080_initial_rtqpcr'
sys.path.insert(0, str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(10)
REFERENCE = OUT/'source_data/cell_line_profiles/reference_material'
for name in ['y_gene_coordinates_grch38.csv', 'chrY_cytobands_grch38.csv']:
    shutil.copy2(REFERENCE/name, DATA/name)
genes = ['RPS4Y1', 'UTY', 'KDM5D', 'EIF1AY']
guides = ['CN', 'CP', 'Q53', 'P47', 'P51']
d = pd.read_csv(DATA/'normalised_run_summaries.csv')
d = d[d.gene.isin(genes) & d.guide.isin(guides)]
summary = d.groupby(['guide','gene'], sort=False).agg(
    mean_relative_expression=('expression_ratio','mean'), n_runs=('run','nunique'),
    minimum_run_expression=('expression_ratio','min'), maximum_run_expression=('expression_ratio','max')).reset_index()
summary['reduction_percent'] = 100*(1-summary.mean_relative_expression)
summary.to_csv(DATA/'mean_expression_reduction.csv', index=False)
values = summary.pivot(index='guide', columns='gene', values='reduction_percent').loc[guides,genes].to_numpy()
assert values.shape == (5,4) and np.isfinite(values).all()
assert summary[summary.guide.eq('CN')].n_runs.eq(1).all()
assert summary[~summary.guide.eq('CN')].n_runs.eq(2).all()
coords = pd.read_csv(DATA/'y_gene_coordinates_grch38.csv').set_index('gene').loc[genes]
bands = pd.read_csv(DATA/'chrY_cytobands_grch38.csv')
length = bands.end.max()

fig = plt.figure(figsize=(183/25.4, 123/25.4))
left, width = .18, .64
axmap = fig.add_axes([left,.70,width,.23])
axmap.set_xlim(0,1)
axmap.set_ylim(0,1)
axmap.set_axis_off()
shades = {'gneg':'#FFFFFF','gpos50':'#B9B9B9','gvar':'#EAEAEA'}
for band in bands.itertuples():
    x0, x1 = band.start/length, band.end/length
    if band.stain == 'acen':
        verts = [(x0,.56),(x0,.72),(x1,.64)] if band.band.startswith('p') else [(x0,.64),(x1,.72),(x1,.56)]
        axmap.add_patch(Polygon(verts, fc='#777777', ec='#777777', lw=.6))
    else:
        axmap.add_patch(Rectangle((x0,.56),x1-x0,.16,fc=shades[band.stain],ec='#777777',lw=.6))
axmap.text(0,1.02,'Y chromosome',fontsize=10,weight='semibold',ha='left')
axmap.text(1,1.02,'GRCh38 | 57.2 Mb',fontsize=8.5,color='#555555',ha='right')
axmap.text(.08,.82,'Yp',fontsize=8.5,ha='center')
axmap.text(.60,.82,'Yq',fontsize=8.5,ha='center')
for i,(gene,row) in enumerate(coords.iterrows()):
    gx=((row.start-1+row.end)/2)/length
    target=(i+.5)/4
    axmap.plot([gx,gx],[.55,.74],color='#353535',lw=1.3)
    axmap.plot([gx,gx,target,target],[.55,.49,.32,.26],color='#999999',lw=.8)

ax = fig.add_axes([left,.20,width,.48])
cmap=LinearSegmentedColormap.from_list('expression_reduction',['#3B6FB6','#FAFAFA','#C75A2B'])
im=ax.imshow(values,cmap=cmap,norm=TwoSlopeNorm(vmin=-100,vcenter=0,vmax=100),aspect='auto')
ax.set_xticks(range(4),genes,fontsize=10,weight='semibold')
ax.xaxis.tick_top()
ax.set_yticks(range(5),guides,fontsize=10)
ax.tick_params(axis='both',length=0,pad=9)
ax.set_ylabel('Guide',labelpad=16,fontsize=10)
for spine in ax.spines.values(): spine.set_visible(False)
for x in np.arange(.5,4,1): ax.axvline(x,color='white',lw=2)
for y in np.arange(.5,5,1): ax.axhline(y,color='white',lw=2)
for i in range(5):
    for j in range(4):
        v=values[i,j]
        ax.text(j,i,f'{v:.0f}%',ha='center',va='center',fontsize=12,
                color='white' if abs(v)>65 else '#222222',weight='semibold')
cax=fig.add_axes([.86,.27,.023,.34])
cb=fig.colorbar(im,cax=cax,ticks=[-100,-50,0,50,100])
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8,length=2)
cb.set_label('Expression reduction (%)',labelpad=9,fontsize=9)
fig.text(.5,.12,'Reduction = 100 \u00d7 (1 \u2212 mean relative expression)',ha='center',fontsize=10)
fig.text(.5,.068,'Mean of two repeat assays of the same edits; CN: one assay',ha='center',fontsize=8.5,color='#555555')
fig.text(.5,.03,'Positive: lower expression than WT     Negative: higher expression than WT',ha='center',fontsize=8.5,color='#555555')
stem='ch03_main_ht1080_rtqpcr_reduction_ideogram'
for ext in ['pdf','svg','png']: fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
fig.savefig(OUT/'previews'/f'{stem}.png',dpi=180)
plt.close(fig)
(DATA/'reduction_ideogram_provenance.json').write_text(json.dumps(dict(
    figure=stem,size_mm=[183,123],genes=genes,guides=guides,assembly='GRCh38',
    summary='Equal-weight arithmetic mean of within-run mean-Ct expression ratios, then 100*(1-mean).',
    biological_independence='Repeat assays of the same edited populations; not replicate transfections.',
    colour_limits=[-100,100],negative_values='Expression above WT; not clipped.',
    ideogram='Full chromosome to scale; gene midpoints linked to equally spaced columns. Coordinates from frozen Ensembl/UCSC tables.',
    omitted='CD99 excluded at author request; all source measurements retained.'),indent=2),encoding='utf-8')
print('Rendered mean-reduction ideogram with 20 guide/gene summaries.')
