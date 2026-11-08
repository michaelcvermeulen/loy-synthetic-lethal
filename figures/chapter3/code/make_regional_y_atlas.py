"""Render the adopted three-panel regional Y atlas without modifying input tables."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('cell', choices=['HT1080', 'A549', 'DU145'])
parser.add_argument('--output-dir', type=Path, default=FIG)
args = parser.parse_args()
CELL, SLUG, OUT = args.cell, args.cell.lower(), args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
(OUT / 'previews').mkdir(exist_ok=True)
W = FIG / 'source_data/regional_y_atlases_ypar_20261028' / SLUG
R = W
SHARED = W
SNAP = OUT / 'source_data/regional_y_atlases_ypar_20261028_rendered' / SLUG
SNAP.mkdir(parents=True, exist_ok=True)
DNA_TOP = 156  # Shared starting scale; expand below if any observed point requires it.
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)

STEM = f'ch03_main_{SLUG}_lpwgs_chrY_atlas'
COL = {'sgCTRL': '#3C5488', 'LOY': '#E64B35'}
TEAL = '#007786'
SHADE = '#EAF6F7'
GREY = '#858585'
FLOOR = .001  # Saved numerical log transformation floor, not a detection limit.

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

source_paths = [W / 'qdnaseq_bins.csv',
                W / 'ddpcr_target_loci.csv',
                SHARED / 'T2T_chrY_cytobands.csv', R / 'gene_mapping_plotted.csv',
                R / 'relative_expression_per_sample.csv', R / 'relative_expression_group_means.csv',
                R / 'sample_metadata.csv', W / 'binned_coverage.csv', W / 'autosomal_normalization.csv', W / 'bin_annotation.csv', W / 'sample_map.json']
before = {str(p.relative_to(ROOT)): sha(p) for p in source_paths}
for p in source_paths:
    dest = SNAP / p.name
    if dest.exists():
        assert sha(dest) == sha(p), f'Source changed: archive existing derived data before rebuilding: {p}'
    else:
        shutil.copy2(p, dest)

cyto = pd.read_csv(SNAP / 'T2T_chrY_cytobands.csv').query("chrom == 'chrY'")
length_bp = int(cyto.end.max())
c0_bp = int(cyto.loc[cyto.stain.eq('acen'), 'start'].min())
c1_bp = int(cyto.loc[cyto.stain.eq('acen'), 'end'].max())
p1_bp = 2458320
q0_bp = int(cyto.loc[cyto.band.eq('q12'), 'start'].iloc[0])
p2_bp = 62122809
p1, c0, c1, q0, p2, length = np.array(
    [p1_bp, c0_bp, c1_bp, q0_bp, p2_bp, length_bp]) / 1e6
COMPRESSED_WIDTH = 3.8

def gx(x):
    """Genomic Mb -> shared display coordinate; only Yq12 is compressed."""
    a = np.asarray(x, dtype=float)
    result = np.where(a <= q0, a, np.where(
        a < p2, q0 + (a - q0) * COMPRESSED_WIDTH / (p2 - q0),
        q0 + COMPRESSED_WIDTH + a - p2))
    return float(result) if result.ndim == 0 else result

regions = [('PAR1', 0, p1_bp), ('MSY Yp', p1_bp, c0_bp),
           ('centromere', c0_bp, c1_bp), ('euchromatic Yq', c1_bp, q0_bp),
           ('Yq12 heterochromatin', q0_bp, p2_bp), ('PAR2', p2_bp, length_bp)]

if CELL == 'DU145':
    regions[1:2] = [('distal MSY Yp', p1_bp, 4800000), ('proximal MSY Yp', 4800000, c0_bp)]
d = pd.read_csv(SNAP / 'qdnaseq_bins.csv')
d['label'] = d.label.replace({'LOY-2X': 'LOY'})
expected_samples = {'HT1080': {'HT1080-sgCTRL-D0', 'HT1080-LOY-D0'},
                    'A549': {'A549-CTRL2', 'A549-LOY'},
                    'DU145': {'LNCaP-sgCTRL', 'LNCaP-LOY-2X'}}
assert set(d['sample']) == expected_samples[CELL]
assert np.allclose(d.log2, np.log2(np.maximum(d.copynumber, FLOOR)), atol=1e-12, equal_nan=True)
y = d.loc[d.chromosome.eq('chrY')].copy()
keys = ['bin_id', 'start', 'end', 'midpoint']
fields = ['copynumber', 'log2', 'log2_seg', 'use']
assert not y.duplicated(['bin_id', 'label']).any()
paired = y.loc[y.label.eq('sgCTRL'), keys + fields].merge(
    y.loc[y.label.eq('LOY'), keys + fields], on=keys, how='outer',
    suffixes=('_sgCTRL', '_LOY'), validate='one_to_one', indicator=True)
assert paired['_merge'].eq('both').all()
paired = paired.drop(columns='_merge').sort_values('start').reset_index(drop=True)
assert len(paired) == 625

# Original start/end are one-based inclusive; regions are BED-style boundaries.
# Reject bins crossing regional boundaries rather than averaging across a mask edge.
paired['region'] = 'boundary-crossing'
for name, left, right in regions:
    inside = paired.start.gt(left) & paired.end.le(right)
    paired.loc[inside, 'region'] = name
paired['excluded_centromere_or_Yq12'] = (
    (paired.start.le(c1_bp) & paired.end.gt(c0_bp)) |
    (paired.start.le(p2_bp) & paired.end.gt(q0_bp)))
paired['valid_pair'] = (
    paired.use_sgCTRL & paired.use_LOY &
    np.isfinite(paired[['log2_sgCTRL', 'log2_LOY', 'log2_seg_sgCTRL',
                       'log2_seg_LOY', 'copynumber_sgCTRL', 'copynumber_LOY']]).all(axis=1))
paired['control_above_floor'] = paired.copynumber_sgCTRL.gt(FLOOR)
paired['depleted_at_floor'] = paired.copynumber_LOY.le(FLOOR)
sample_map=json.loads((SNAP/'sample_map.json').read_text())
coverage=pd.read_csv(SNAP/'binned_coverage.csv')
control_counts=coverage.loc[coverage.track.eq('Y') & coverage['sample'].eq(sample_map['control'])].set_index('end').read_count
paired['control_primary_reads']=paired.end.map(control_counts)
assert paired.control_primary_reads.notna().all()
annotation=pd.read_csv(SNAP/'bin_annotation.csv').set_index('bin_id')
for field in ['bases','mappability','blacklist']:
    paired['reference_'+field]=paired.bin_id.map(annotation[field])
paired['original_strict_mask']=paired.bin_id.map(annotation['use'])
paired['included'] = (paired.valid_pair & paired.control_above_floor & paired.control_primary_reads.ge(50) &
                      ~paired.excluded_centromere_or_Yq12 &
                      paired.region.ne('boundary-crossing'))
paired['percent_of_sgCTRL'] = np.nan
ok = paired.included
paired.loc[ok, 'percent_of_sgCTRL'] = (
    100 * paired.loc[ok,'copynumber_LOY'].clip(lower=0) / paired.loc[ok,'copynumber_sgCTRL'])
assert paired.loc[ok,'copynumber_LOY'].ge(-1e-12).all()
paired['display_x'] = gx(paired.midpoint)
paired.to_csv(SNAP / 'paired_Y_bins_and_exclusions.csv', index=False)
plot_data = paired.loc[ok].copy()
assert len(plot_data)>0 and not plot_data.region.str.startswith('PAR').any()
assert np.isfinite(plot_data.percent_of_sgCTRL).all()
assert not plot_data.excluded_centromere_or_Yq12.any()
assert plot_data.control_above_floor.all()
summary = plot_data.groupby('region', sort=False).agg(
    n_bins=('bin_id', 'size'),
    median_percent_of_sgCTRL=('percent_of_sgCTRL', 'median'),
    n_depleted_bins_at_floor=('depleted_at_floor', 'sum'),
    min_control_relative_CN=('copynumber_sgCTRL', 'min'))
summary.to_csv(SNAP / 'regional_median_paired_ratios.csv')

mapping = pd.read_csv(SNAP / 'gene_mapping_plotted.csv').sort_values('t2t_start').reset_index(drop=True)
probes = pd.read_csv(SNAP / 'ddpcr_target_loci.csv')
relative = pd.read_csv(SNAP / 'relative_expression_per_sample.csv')
genes = mapping.gene_name.tolist()
n = len(genes)
expected_n = {'HT1080': 6, 'A549': 4, 'DU145': 3}[CELL]
assert n == 11 and len(relative) == 22*expected_n
assert relative.groupby(['gene_name', 'condition']).size().eq(expected_n).all()
meta = pd.read_csv(SNAP / 'sample_metadata.csv')
assert np.allclose(relative.relative_expression_percent,
                   100 * relative.cpm / relative.control_mean_cpm.where(relative.control_mean_cpm.gt(0)), equal_nan=True)
rna_summary = relative.groupby(['condition', 'gene_name']).relative_expression_percent.mean().unstack(0)
assert np.allclose(rna_summary.sgCTRL.dropna(), 100)
assert relative.relative_expression_percent.dropna().between(0, 185).all()
undefined = sorted(relative.loc[relative.relative_expression_percent.isna(), 'gene_name'].unique())
assert undefined == (['NLGN4Y'] if CELL == 'DU145' else [])

norm=pd.read_csv(SNAP/'autosomal_normalization.csv').set_index('sample').median_autosomal_depth
px=coverage.loc[coverage.track.ne('Y')]
keys=['chromosome','start0','end','track']
px=px.loc[px['sample'].eq(sample_map['control'])].merge(px.loc[px['sample'].eq(sample_map['LOY'])],on=keys,
    suffixes=('_control','_LOY'),validate='one_to_one')
px['width']=px.end-px.start0
px['included']=(px.read_count_control>=20)&(px.read_count_control/px.width*100000>=50)&px.depth_control.gt(0)
px['percent_of_sgCTRL']=100*(px.depth_LOY/norm[sample_map['LOY']])/(px.depth_control/norm[sample_map['control']])
assert px.loc[px.included,'percent_of_sgCTRL'].notna().all()
par_summary=px.loc[px.included].groupby('track').agg(n_bins=('start0','size'),
    median_percent_of_sgCTRL=('percent_of_sgCTRL','median'),minimum_control_reads=('read_count_control','min'))
assert set(par_summary.index)=={'PAR1_X','PAR2_X'}
par_summary.to_csv(SNAP/'X_PAR_regional_medians.csv')
maximum=max(plot_data.percent_of_sgCTRL.max(),px.loc[px.included,'percent_of_sgCTRL'].max())
if maximum>DNA_TOP-6: DNA_TOP=float(np.ceil((maximum+12)/50)*50)
assert plot_data.percent_of_sgCTRL.between(0,DNA_TOP).all()
assert px.loc[px.included,'percent_of_sgCTRL'].between(0,DNA_TOP).all()
fig = plt.figure(figsize=(183 / 25.4, 169 / 25.4))
left, right = .115, .955
ax_a = fig.add_axes([left, .665, right-left, .275])
ax_b = fig.add_axes([left, .415, right-left, .215])
ax_c = fig.add_axes([left, .100, right-left, .245])
for ax, letter, y_tag in [(ax_a, 'a', .956), (ax_b, 'b', .618), (ax_c, 'c', .360)]:
    fig.text(.032, y_tag, letter, fontsize=12, fontweight='semibold', va='top')

xlim = (-.15, gx(length) + .15)
for ax in [ax_a, ax_b]:
    ax.set_xlim(xlim)

# a: Current cell-specific Y regions, with X-mapped PAR coverage in labelled spans.
for name,a,b in regions:
    a,b=a/1e6,b/1e6
    if name.startswith('PAR'):
        ax_a.axvspan(gx(a),gx(b),color=SHADE,zorder=0)
    if name in ['centromere','Yq12 heterochromatin']:
        ax_a.axvspan(gx(a),gx(b),facecolor='#F0F0F0',edgecolor='#D5D5D5',hatch='///',lw=.4,zorder=1)
for a,b in [(0,q0),(p2,length)]:
    ax_a.plot([gx(a),gx(b)],[100,100],color=COL['sgCTRL'],lw=.8,ls=(0,(4,3)),zorder=2)
ax_a.scatter(plot_data.display_x,plot_data.percent_of_sgCTRL,s=8,color=GREY,alpha=.65,linewidth=0,zorder=3)
for name,a,b in regions:
    if name in summary.index:
        value=float(summary.loc[name,'median_percent_of_sgCTRL'])
        ax_a.plot([gx(a/1e6),gx(b/1e6)],[value]*2,color=COL['LOY'],lw=2,solid_capstyle='butt',zorder=4)
ax_a.set(ylim=(-6,DNA_TOP),yticks=list(range(0,int(DNA_TOP),50)),ylabel='Y-depleted DNA\n(% of sgCTRL)')
ax_a.spines[['right','top','bottom']].set_visible(False)
ax_a.tick_params(axis='x',bottom=False,labelbottom=False)
def label_at(x,y,label,color,ha='center'):
    ax_a.text(x,y,label,color=color,fontsize=7.6,ha=ha,va='center',fontweight='semibold',
        bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1.1),zorder=7)
def format_median(region):
    value=float(summary.loc[region,'median_percent_of_sgCTRL']);n=int(summary.loc[region,'n_bins'])
    return ('<1%' if value<1 else f'{value:.1f}%'),n,value
if CELL in ['HT1080', 'A549']:
    for region,x in [('MSY Yp',6.5),('euchromatic Yq',19.0)]:
        value,bin_count,v=format_median(region)
        label_at(x,max(v+16,DNA_TOP*.14),f'{value} ({bin_count} bins)',COL['LOY'])
else:
    for region,x in [('distal MSY Yp',3.63),('proximal MSY Yp',7.6),('euchromatic Yq',19.0)]:
        value,bin_count,v=format_median(region)
        sep='\n' if region=='distal MSY Yp' else ' '
        local_max=plot_data.loc[plot_data.region.eq(region),'percent_of_sgCTRL'].max()
        label_at(x,min(DNA_TOP-25,max(v+20,local_max+17)),f'{value}{sep}({bin_count} bins)',COL['LOY'])
    ax_a.axvline(4.8,color='#BBBBBB',lw=.6,ls=':',zorder=1)
    ax_a.text(3.63,DNA_TOP*.85,'Distal Yp',ha='center',fontsize=6.8)
    ax_a.text(7.6,DNA_TOP*.85,'Proximal Yp',ha='center',fontsize=6.8)
for region,x0,x1,y0,y1 in [('PAR1_X',0,2394410,0,p1_bp),('PAR2_X',153925834,154259566,p2_bp,length_bp)]:
    keep=px.track.eq(region)&px.included
    xx=px.loc[keep]
    assert len(xx)>0
    px.loc[keep,'projected_Y_span_midpoint']=(y0+((xx.start0+xx.end)/2-x0)/(x1-x0)*(y1-y0))/1e6
    xx=px.loc[keep];v=float(xx.percent_of_sgCTRL.median())
    ax_a.scatter(gx(xx.projected_Y_span_midpoint),xx.percent_of_sgCTRL,s=8,color=TEAL,alpha=.65,linewidth=0,zorder=3)
    ax_a.plot([gx(y0/1e6),gx(y1/1e6)],[v]*2,color=TEAL,lw=2,zorder=4)
    label_y=min(DNA_TOP-15,max(v+28,DNA_TOP*.718))
    if region=='PAR1_X':label_at(p1/2,label_y,f'{v:.0f}%',TEAL)
    else:
        ax_a.annotate(f'{v:.0f}%',xy=(gx((y0+y1)/2e6),v),xytext=(gx(length)+.08,label_y),
            ha='right',va='center',fontsize=7.6,fontweight='semibold',color=TEAL,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1),
            arrowprops=dict(arrowstyle='-',lw=.6,color=TEAL),zorder=7)
px.to_csv(SNAP/'X_PAR_bin_ratios.csv',index=False)
ax_a.text(gx((q0+p2)/2),DNA_TOP*.48718,'Yq12 heterochromatin\n(unmappable, excluded)',rotation=90,
    ha='center',va='center',fontsize=7.4,color='#666666',zorder=6)
ax_a.annotate('Centromere\nexcluded',xy=(gx((c0+c1)/2),DNA_TOP*.75641),xytext=(12.4,DNA_TOP*.83333),
    fontsize=6.7,ha='left',va='center',color='#666666',arrowprops=dict(arrowstyle='-',color='#999999',lw=.5),zorder=6)
for x,label,color,align in [(p1/2,'PAR1 (X)',TEAL,'center'),(6.5,'Male-specific Yp','#222222','center'),
                           (19.0,'Euchromatic Yq','#222222','center'),(gx(length)+.13,'PAR2 (X)',TEAL,'right')]:
    ax_a.text(x,DNA_TOP+4,label,fontsize=7.4,ha=align,va='bottom',color=color,clip_on=False)
ax_a.text(gx((q0+p2)/2),DNA_TOP-6,'compressed',fontsize=6.5,ha='center',va='center',color='#666666')
handles=[Line2D([],[],marker='o',ls='none',color=GREY,markersize=3,label='100-kb Y bin'),
         Line2D([],[],color=COL['LOY'],lw=2,label='Regional median'),
         Line2D([],[],color=COL['sgCTRL'],lw=.8,ls='--',label='100% of sgCTRL'),
         Line2D([],[],marker='o',ls='none',color=TEAL,markersize=3,label='X-PAR bin')]
ax_a.legend(handles=handles,loc='upper left',bbox_to_anchor=(0,-.025),ncol=4,frameon=False,
    fontsize=6.8,borderaxespad=0,handlelength=1.6,columnspacing=1.1)


# b: The identical coordinate transform is applied to cytobands and every locus.
for band in cyto.itertuples():
    x0, x1 = gx(band.start/1e6), gx(band.end/1e6)
    if band.stain == 'acen':
        points = [(x0, .15), (x0, .65), (x1, .40)] if band.band.startswith('p') else [
            (x0, .40), (x1, .15), (x1, .65)]
        ax_b.add_patch(Polygon(points, facecolor='#777777', edgecolor='#555555', lw=.45))
    else:
        fill = {'gneg': 'white', 'gpos50': '#AAAAAA', 'gvar': '#E5E5E5'}.get(band.stain, '#CCCCCC')
        ax_b.add_patch(Rectangle((x0, .15), x1-x0, .50, facecolor=fill,
                                edgecolor='#555555', lw=.5,
                                hatch='///' if band.stain == 'gvar' else None))
for a, b in [(0, p1), (p2, length)]:
    ax_b.add_patch(Rectangle((gx(a), .15), gx(b)-gx(a), .5,
                            facecolor='#41B6C4', edgecolor='#555555', lw=.5, zorder=4))

label_positions = [.10, 2.1, 3.8, 5.5, 11.9, 13.6, 15.3, 17.0, 20.0, 22.0, gx(length)-1.75]
for row, label_x in zip(mapping.itertuples(), label_positions):
    actual_x = gx((row.t2t_start + row.t2t_end)/2e6)
    color = TEAL if row.region != 'MSY' else '#666666'
    ax_b.plot([actual_x, actual_x, label_x], [.70, .94, 1.27], color='#A5A5A5', lw=.5)
    ax_b.text(label_x, 1.31, row.gene_name, rotation=55, ha='left', va='bottom',
              fontsize=7.1, fontstyle='italic', color=color, clip_on=False)
for row in probes.itertuples():
    a, b = row.t2t_start/1e6, row.t2t_end/1e6
    middle = gx((a+b)/2)
    ax_b.add_patch(Rectangle((gx(a), .15), gx(b)-gx(a), .5,
                            facecolor='#222222', edgecolor='#222222', lw=.55, zorder=5))
    label_x = {'RPS4Y1': 2.7, 'AMELY': 6.65}.get(row.gene_name, middle)
    ax_b.plot([middle, middle, label_x], [.10, -.06, -.23], color='#444444', lw=.6)
    ax_b.text(label_x, -.29, row.gene_name, ha='center', va='top',
              fontsize=7.2, fontweight='semibold')
ax_b.text(25.5, -.39, 'ddPCR targets', ha='center', va='center', fontsize=6.6, color='#666666')

# Slashes mark the two ends of the compressed region without moving PAR2 loci.
for a in [q0, p2]:
    x = gx(a)
    for dx in [-.11, .11]:
        ax_b.plot([x+dx-.07, x+dx+.07], [.10, .24], color='#555555', lw=.65, zorder=7)
ax_b.text(gx((q0+p2)/2), -.26, 'Yq12 compressed', ha='center', va='top',
          fontsize=6.3, color='#666666')
ax_b.set(ylim=(-1.05, 2.8), yticks=[])
ax_b.set_xticks(gx(np.array([0, 5, 10, 15, 20, 25, 62.3])), ['0', '5', '10', '15', '20', '25', '62.3'])
ax_b.tick_params(axis='x', labelsize=7.3, length=2.5, pad=2)
ax_b.set_xlabel('HG002 T2T-Y position (Mb; Yq12 span compressed)', fontsize=7.4, labelpad=4)
ax_b.spines[['left', 'right', 'top']].set_visible(False)
ax_b.spines['bottom'].set_linewidth(.5)
for boundary in [q0, p2]:
    x = gx(boundary)
    ax_b.plot([x-.18, x+.18], [0, 0], color='white', lw=2.2,
              transform=ax_b.get_xaxis_transform(), clip_on=False, zorder=9)
    for dx in [-.09, .09]:
        ax_b.plot([x+dx-.06, x+dx+.06], [-.015, .015], color='#555555', lw=.65,
                  transform=ax_b.get_xaxis_transform(), clip_on=False, zorder=10)

# c: Preserve every current RNA value, group mean and sampling-window symbol.
for i, row in mapping.iterrows():
    if row.region != 'MSY':
        ax_c.axvspan(i-.5, i+.5, color=SHADE, zorder=0)
for condition, off in [('sgCTRL', -.19), ('LOY', .19)]:
    ax_c.scatter(np.arange(n)+off, rna_summary.loc[genes, condition], marker='_',
                 s=95, color=COL[condition], linewidths=2, zorder=4)
    if CELL == 'HT1080':
        for batch, marker, dx in [('August 2025', 'o', -.055), ('May 2026', '^', .055)]:
            for passage in [1, 2, 3]:
                values = relative.loc[relative.condition.eq(condition) & relative.batch.eq(batch) &
                                      relative.passage_label.eq(passage)].set_index('gene_name').loc[genes].relative_expression_percent
                ax_c.scatter(np.arange(n)+off+dx+(passage-2)*.045, values, s=10,
                             facecolor='white', edgecolor=COL[condition], marker=marker, lw=.7, zorder=3)
    else:
        samples = meta.loc[meta.condition.eq(condition), 'sample'].tolist()
        for j, sample in enumerate(samples):
            values = relative.loc[relative['sample'].eq(sample)].set_index('gene_name').loc[genes].relative_expression_percent
            ax_c.scatter(np.arange(n)+off+(j-(len(samples)-1)/2)*.055, values, s=10,
                         facecolor='white', edgecolor=COL[condition], marker='o', lw=.7, zorder=3)
for j, gene in enumerate(genes):
    if gene in undefined:
        ax_c.text(j, 10, 'NA', ha='center', fontsize=7, color='#666666')
ax_c.axhline(100, color='#888888', lw=.7, ls='--', zorder=1)
ax_c.set(xlim=(-.5, n-.5), ylim=(-5, 185), yticks=[0, 50, 100, 150],
         ylabel='RNA expression\n(% of sgCTRL)')
ax_c.set_xticks(range(n), genes, rotation=55, ha='right', fontsize=7.3, fontstyle='italic')
ax_c.tick_params(axis='x', length=0)
ax_c.spines[['right', 'top']].set_visible(False)
for tick, row in zip(ax_c.get_xticklabels(), mapping.itertuples()):
    if row.region != 'MSY':
        tick.set_color(TEAL)
handles = [Line2D([], [], color=COL[c], lw=3, label=label)
           for c, label in [('sgCTRL', 'sgCTRL'), ('LOY', 'Y-depleted')]]
if CELL == 'HT1080':
    handles += [Line2D([], [], marker=m, ls='none', color='#444444', markerfacecolor='white',
                       markersize=3, label=b) for b, m in [('Aug 2025', 'o'), ('May 2026', '^')]]
else:
    handles += [Line2D([], [], marker='o', ls='none', color='#444444', markerfacecolor='white',
                       markersize=3, label=f'n = {expected_n} per group')]
ax_c.legend(handles=handles, ncol=4, frameon=False, fontsize=6.8,
             loc='upper right', bbox_to_anchor=(1, 1.015), columnspacing=1.1,
             handlelength=1.2, borderaxespad=0)


assert np.allclose(ax_a.get_xlim(), ax_b.get_xlim())
assert np.isclose(ax_a.get_position().x0, ax_b.get_position().x0)
assert np.isclose(ax_a.get_position().width, ax_b.get_position().width)

for ext in ['pdf', 'svg', 'png']:
    fig.savefig(OUT / f'{STEM}.{ext}', dpi=600, facecolor='white')
fig.savefig(OUT / 'previews' / f'{STEM}.png', dpi=180, facecolor='white')
plt.close(fig)

after = {str(p.relative_to(ROOT)): sha(p) for p in source_paths}
assert before == after, 'An active source or protected output changed during draft generation.'
validation = {
    'status': 'Adopted 2026-10-28: Y-PAR-masked realignment, X-PAR coverage and broad 100-kb filter',
    'cell_line': CELL,
    'current_draft_number': {'HT1080': '3.4', 'A549': '3.5', 'DU145': '3.6'}[CELL],
    'removed_annotations': ['blue explanatory text below a', 'blue explanatory text below c'],
    'working_id': STEM,
    'panel_map': {'a': 'new paired-bin sgCTRL-normalised lpWGS',
                  'b': 'old d, same locus coordinates with compressed Yq12',
                  'c': 'old e, same RNA values and normalisation'},
    'DNA_formula': '100 * max(corrected_CN_LOY,0)/corrected_CN_control; no positive numerator floor',
    'X_PAR_formula': '100 * (LOY primary-read depth / LOY autosomal median) / (control primary-read depth / control autosomal median)',
    'saved_log2_floor': FLOOR,
    'floor_is_numerical_not_detection_limit': True,
    'bin_rule': 'Reference non-N>=60%, mappability>=30%; no blanket difficulty-overlap exclusion; >=50 primary MAPQ>=37 control reads per full Y bin; finite CN/segmentation, control CN>0.001; whole bins within MSY regions; centromere/Yq12 excluded',
    'separate_mappability_mask_available': True,
    'mappability_operational_definition': 'Original reference annotation with a >=30% threshold; explicit sample control-read support; difficulty overlap retained as metadata',
    'PAR_coordinate_note': 'X-relative positions projected within labelled PAR spans, not a coordinate liftover',
    'X_PAR_bins_displayed': int(px.included.sum()),
    'X_PAR_regions': json.loads(par_summary.reset_index().to_json(orient='records')),
    'DU145_Yp_split_bp': 4800000 if CELL=='DU145' else None,
    'Y_bins_available': len(paired), 'Y_bins_displayed': len(plot_data),
    'excluded_bins': int((~paired.included).sum()),
    'Yq12_interval_bp': [q0_bp, p2_bp], 'centromere_interval_bp': [c0_bp, c1_bp],
    'regions': json.loads(summary.reset_index().to_json(orient='records')),
    'x_axis': '0–27.449937 Mb and PAR2 retain equal scale; only 27.449937–62.122809 Mb compressed to 3.8 display units; a and b identical transform',
    'n_RNA_observations': len(relative), 'n_genes': n,
    'RNA_values_and_means_unchanged': True,
    'RNA_display_limits': [-5, 185],
    'RNA_display_note': 'Small space below zero prevents near-zero point symbols from being cut by the axis; values and ticks unchanged.',
    'no_points_clipped': True,
    'source_files_unchanged': True,
    'DNA_display_limits': [-6, DNA_TOP],
    'RNA_undefined_genes': undefined,
    'active_files_and_source_sha256': before,
    'export_sha256': {ext: sha(OUT / f'{STEM}.{ext}') for ext in ['pdf', 'svg', 'png']},
}
(SNAP / 'validation.json').write_text(json.dumps(validation, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
print(summary.to_string())
print(f'Saved {OUT / (STEM + ".pdf")}')
print(f'Input tables unchanged; {len(plot_data)} paired DNA bins and {len(relative)} RNA observations verified.')
