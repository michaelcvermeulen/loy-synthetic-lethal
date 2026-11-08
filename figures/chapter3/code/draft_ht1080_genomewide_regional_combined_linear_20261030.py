"""Compact draft with unlogged linear DNA/RNA fold-change axes in panels b/c.

Uses frozen display tables; recalculates descriptive Pearson correlations for
the back-transformed displayed values, without refitting the underlying models.
Does not replace selected figures. Regional plotting code is adapted from the adopted
Y-PAR-masked atlas renderer, with the same values and coordinate transform.
"""
from pathlib import Path
import hashlib
import json
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
OUT = FIG / 'review/ht1080_genomewide_regional_combined_20261030_v2_linear'
SOURCE = FIG / 'source_data/ht1080_genomewide_regional_combined_linear_draft_20261030'
SNAP = SOURCE / 'inputs'
STEM = 'ch03_main_ht1080_genomewide_regional_combined_draft'
for folder in [OUT, SNAP]:
    folder.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(7.5)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

global_source = FIG / 'source_data/ht1080_lpwgs_difference_draft'
regional_source = FIG / 'source_data/regional_y_atlases_ypar_20261028_rendered/ht1080'
input_sources = {}
for filename in ['genome_5Mb.csv', 'DNA_RNA_window_1Mb_plotted.csv',
                 'DNA_RNA_arm_plotted.csv', 'DNA_RNA_correlations.csv',
                 'RNA_scatter_arm_labels.csv']:
    input_sources[filename] = global_source / filename
for filename in ['paired_Y_bins_and_exclusions.csv', 'regional_median_paired_ratios.csv',
                 'X_PAR_bin_ratios.csv', 'X_PAR_regional_medians.csv',
                 'T2T_chrY_cytobands.csv', 'ddpcr_target_loci.csv', 'gene_mapping_plotted.csv',
                 'relative_expression_per_sample.csv', 'relative_expression_group_means.csv',
                 'sample_metadata.csv']:
    input_sources[filename] = regional_source / filename
input_sources['all_chromosome_cytobands.bed'] = ROOT.parent / 'WGS/chm13v2.0_cytobands_allchrs.bed'
input_manifest = {}
for filename, original in input_sources.items():
    target = SNAP / filename
    if target.exists():
        assert sha(target) == sha(original), f'Frozen draft input differs: {filename}'
    else:
        shutil.copy2(original, target)
    input_manifest[filename] = {'source': str(original), 'sha256': sha(original)}

protected_paths = [ROOT / 'draftv2/chapter3/chapter3_results.md',
                   ROOT / 'draftv2/chapter3/chapter3_figure_numbering.md',
                   ROOT / 'Chapter3_figure_captions.md',
                   FIG / 'code/make_regional_y_atlas.py',
                   FIG / 'code/ht1080_lpwgs_revision_panels.py']
protected_paths.append(FIG / 'code/draft_ht1080_genomewide_regional_combined_20261030.py')
for prior_folder in [FIG / 'review/ht1080_genomewide_regional_combined_20261030_v1',
                     FIG / 'source_data/ht1080_genomewide_regional_combined_draft_20261030']:
    protected_paths.extend(p for p in prior_folder.rglob('*') if p.is_file())
for folder in ['CH3_S3_4_5_1_F03_ht1080_genomewide', 'CH3_S3_4_5_1_F04_ht1080_regional_y']:
    protected_paths.extend(p for p in (FIG / 'selected' / folder).iterdir() if p.is_file())
protected = {str(p.relative_to(ROOT)): sha(p) for p in protected_paths}

COL = {'sgCTRL': '#3C5488', 'LOY': '#E64B35'}
RED, BLUE, GREEN = COL['LOY'], COL['sgCTRL'], '#00A087'
TEAL, SHADE, GREY = '#007786', '#EAF6F7', '#858585'
CELL = 'HT1080'
DNA_TOP = 156
RNA_YLIM = (-.10, 3.05)
RNA_YTICKS = [0, 1, 2, 3]
DNA_XLIM = (-.05, 1.75)
WIDTH, HEIGHT = 7.2, 8.6
fig = plt.figure(figsize=(WIDTH, HEIGHT))
genome_ax = fig.add_axes([.090, .855, .885, .115])
window_ax = fig.add_axes([.090, .610, .385, .180])
arm_ax = fig.add_axes([.590, .610, .385, .180], sharey=window_ax)
ax_a = fig.add_axes([.090, .405, .885, .130])
ax_b = fig.add_axes([.090, .240, .885, .130])
ax_c = fig.add_axes([.090, .065, .885, .120])
for letter, x, y in [('a', .018, .982), ('b', .018, .802), ('c', .518, .802),
                     ('d', .018, .551), ('e', .018, .380), ('f', .018, .197)]:
    fig.text(x, y, letter, fontsize=12, fontweight='semibold', va='top')

def clean(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=7, pad=2, length=2.5)

def chromosome_key():
    return [Line2D([], [], marker='o', ls='', color='#777777', ms=3, label='Autosomes / X'),
            Line2D([], [], marker='^', ls='', color=RED, ms=4, label='Chromosome Y')]

bins = pd.read_csv(SNAP / 'genome_5Mb.csv')
windows = pd.read_csv(SNAP / 'DNA_RNA_window_1Mb_plotted.csv')
arms = pd.read_csv(SNAP / 'DNA_RNA_arm_plotted.csv')
correlations = pd.read_csv(SNAP / 'DNA_RNA_correlations.csv').set_index('level')
assert [len(bins), len(windows), len(arms)] == [630, 1597, 43]
# Exponentiate the saved log2 summaries; these are geometric fold changes,
# not recomputed arithmetic ratios of raw RNA means. Preserve source columns.
linear_correlations = []
for table, level, chrom_col in [(windows, 'window_1Mb', 'chrom'), (arms, 'arm', 'chromosome')]:
    table['DNA_fold_change'] = np.exp2(table.DNA)
    table['RNA_fold_change'] = np.exp2(table.RNA)
    assert np.allclose(np.log2(table.DNA_fold_change), table.DNA)
    assert np.allclose(np.log2(table.RNA_fold_change), table.RNA)
    nony = table[chrom_col].ne('chrY')
    linear_correlations.append(dict(level=level, n=len(table),
        Pearson=table.DNA_fold_change.corr(table.RNA_fold_change),
        nonY_n=int(nony.sum()),
        nonY_Pearson=table.loc[nony, 'DNA_fold_change'].corr(table.loc[nony, 'RNA_fold_change'])))
    table.to_csv(SOURCE / f'{level}_linear_display.csv', index=False)
linear_correlations = pd.DataFrame(linear_correlations).set_index('level')
linear_correlations.to_csv(SOURCE / 'DNA_RNA_linear_correlations.csv')
cb = pd.read_csv(SNAP / 'all_chromosome_cytobands.bed', sep='\t', header=None,
                 names=['chrom', 'start', 'end', 'band', 'stain'])
chroms = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY']
lengths = cb.groupby('chrom').end.max()
offsets, end = {}, 0
for ch in chroms:
    offsets[ch] = end
    end += lengths[ch] / 1e6 + 20
yq = cb[cb.chrom.eq('chrY') & cb.band.eq('q12')].iloc[0]
for i, ch in enumerate(chroms):
    z = bins[bins.chromosome.eq(ch)]
    genome_ax.scatter(z.mid/1e6 + offsets[ch], z.difference,
                      s=17 if ch == 'chrY' else 4.8,
                      marker='^' if ch == 'chrY' else 'o',
                      color=RED if ch == 'chrY' else ['#555555', '#AAAAAA'][i % 2],
                      lw=0, zorder=4 if ch == 'chrY' else 2)
genome_ax.axhline(0, color='#BBBBBB', lw=.5, zorder=0)
genome_ax.axvspan(offsets['chrY']+yq.start/1e6, offsets['chrY']+yq.end/1e6,
                  color='#EEEEEE', zorder=-1)
genome_ax.annotate('Yq12 excluded', xy=(offsets['chrY']+(yq.start+yq.end)/2e6, 0),
                   xytext=(.988, .96), textcoords='axes fraction', ha='right', va='top',
                   fontsize=6.5, color='#666666',
                   arrowprops=dict(arrowstyle='-', lw=.45, color='#888888'))
extent = max(10, np.ceil(bins.difference.abs().max()))
genome_ax.set(xlim=(-10, end), ylim=(-extent-.5, extent+.5), yticks=[-10, -5, 0, 5, 10],
              ylabel='Δ log$_2$(relative CN)', xlabel='Chromosome')
genome_ax.set_xticks([offsets[ch]+lengths[ch]/2e6 for ch in chroms], [ch[3:] for ch in chroms])
clean(genome_ax)
genome_ax.tick_params(axis='x', labelsize=6.8, pad=2, rotation=45)
genome_ax.get_xticklabels()[-1].set_color(RED)
genome_ax.legend(handles=chromosome_key(), loc='upper left', ncol=2, frameon=False,
                 fontsize=6.8, borderaxespad=.1, handletextpad=.35, columnspacing=1)
genome_ax.text(.01, .045, '5-Mb windows; Y-depleted − sgCTRL', transform=genome_ax.transAxes,
               fontsize=6.6, va='bottom')

for ax in [window_ax, arm_ax]:
    ax.axhline(1, color='#DDDDDD', lw=.5, zorder=0)
    ax.axvline(1, color='#DDDDDD', lw=.5, zorder=0)
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.set(xlim=DNA_XLIM, xticks=[0, .5, 1, 1.5], ylim=RNA_YLIM, yticks=RNA_YTICKS,
           xlabel='DNA fold change', ylabel='RNA fold change')
    clean(ax)
isy = windows.chrom.eq('chrY')
window_ax.scatter(windows.loc[~isy, 'DNA_fold_change'], windows.loc[~isy, 'RNA_fold_change'], s=5.5, color='#777777', alpha=.55, lw=0)
window_ax.scatter(windows.loc[isy, 'DNA_fold_change'], windows.loc[isy, 'RNA_fold_change'], s=22, marker='^', color=RED, lw=0, zorder=4)
window_ax.legend(handles=chromosome_key(), loc='upper left', frameon=False, fontsize=6.8,
                  handletextpad=.35, borderpad=.1, labelspacing=.2)
isy = arms.chromosome.eq('chrY')
arm_ax.scatter(arms.loc[~isy, 'DNA_fold_change'], arms.loc[~isy, 'RNA_fold_change'], s=14, color='#777777', alpha=.85, lw=0)
arm_ax.scatter(arms.loc[isy, 'DNA_fold_change'], arms.loc[isy, 'RNA_fold_change'], s=23, marker='^', color=RED, lw=0)
label_offsets = {'4p': (0, -11), '4q': (15, 19), '18p': (-18, 11),
                 '13q': (4, 6), 'Yp': (3, 6), 'Yq': (3, 6)}
labels = pd.read_csv(SNAP / 'RNA_scatter_arm_labels.csv')
for row in labels.itertuples():
    color = RED if row.arm.startswith('Y') else GREEN if row.arm == '13q' else BLUE
    arm_ax.scatter([np.exp2(row.DNA)], [np.exp2(row.RNA)], s=24, marker='^' if row.arm.startswith('Y') else 'o', color=color, lw=0, zorder=4)
    arm_ax.annotate(row.arm, (np.exp2(row.DNA), np.exp2(row.RNA)), xytext=label_offsets[row.arm],
                    textcoords='offset points', fontsize=6.8, color=color,
                    arrowprops=dict(arrowstyle='-', lw=.4, color='#999999'), zorder=6)
for ax, level, label in [(window_ax, 'window_1Mb', '1-Mb windows'), (arm_ax, 'arm', 'Chromosome arms')]:
    row = linear_correlations.loc[level]
    ax.text(.98, .965, f'{label}\nr (all) = {row.Pearson:.3f}; n = {int(row.n):,}\n'
            f'r (non-Y) = {row.nonY_Pearson:.3f}; n = {int(row.nonY_n):,}',
            transform=ax.transAxes, ha='right', va='top', fontsize=6.7,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=.6))
for values, limits in [(bins.difference, genome_ax.get_ylim()),
                       (windows.DNA_fold_change, window_ax.get_xlim()), (windows.RNA_fold_change, window_ax.get_ylim()),
                       (arms.DNA_fold_change, arm_ax.get_xlim()), (arms.RNA_fold_change, arm_ax.get_ylim())]:
    assert values.between(*limits).all(), 'A retained value would be clipped'
assert arm_ax.get_ylim() == window_ax.get_ylim()
assert arm_ax.get_xlim() == window_ax.get_xlim()
assert all(ax.get_xscale() == ax.get_yscale() == 'linear' for ax in [window_ax, arm_ax])
# Verify scatter coordinates directly against the transformed saved table.
for ax, table, chrom_col in [(window_ax, windows, 'chrom'), (arm_ax, arms, 'chromosome')]:
    is_y = table[chrom_col].eq('chrY')
    for artist, keep in zip(ax.collections[:2], [~is_y, is_y]):
        assert np.allclose(artist.get_offsets(), table.loc[keep, ['DNA_fold_change', 'RNA_fold_change']])

# Read already-verified adopted regional display values; do not re-filter/re-fit.
cyto = pd.read_csv(SNAP / 'T2T_chrY_cytobands.csv').query("chrom == 'chrY'")
length_bp = int(cyto.end.max())
c0_bp = int(cyto.loc[cyto.stain.eq('acen'), 'start'].min())
c1_bp = int(cyto.loc[cyto.stain.eq('acen'), 'end'].max())
p1_bp, p2_bp = 2458320, 62122809
q0_bp = int(cyto.loc[cyto.band.eq('q12'), 'start'].iloc[0])
p1, c0, c1, q0, p2, length = np.array([p1_bp, c0_bp, c1_bp, q0_bp, p2_bp, length_bp])/1e6
COMPRESSED_WIDTH = 3.8
def gx(x):
    a = np.asarray(x, dtype=float)
    result = np.where(a <= q0, a, np.where(a < p2,
        q0+(a-q0)*COMPRESSED_WIDTH/(p2-q0), q0+COMPRESSED_WIDTH+a-p2))
    return float(result) if result.ndim == 0 else result
regions = [('PAR1', 0, p1_bp), ('MSY Yp', p1_bp, c0_bp), ('centromere', c0_bp, c1_bp),
           ('euchromatic Yq', c1_bp, q0_bp), ('Yq12 heterochromatin', q0_bp, p2_bp),
           ('PAR2', p2_bp, length_bp)]
paired = pd.read_csv(SNAP / 'paired_Y_bins_and_exclusions.csv')
plot_data = paired.loc[paired.included].copy()
summary = pd.read_csv(SNAP / 'regional_median_paired_ratios.csv', index_col='region')
px = pd.read_csv(SNAP / 'X_PAR_bin_ratios.csv')
par_summary = pd.read_csv(SNAP / 'X_PAR_regional_medians.csv', index_col='track')
mapping = pd.read_csv(SNAP / 'gene_mapping_plotted.csv').sort_values('t2t_start').reset_index(drop=True)
probes = pd.read_csv(SNAP / 'ddpcr_target_loci.csv')
relative = pd.read_csv(SNAP / 'relative_expression_per_sample.csv')
rna_summary = pd.read_csv(SNAP / 'relative_expression_group_means.csv', index_col=0)
meta = pd.read_csv(SNAP / 'sample_metadata.csv')
genes, n, undefined = mapping.gene_name.tolist(), len(mapping), []
assert len(plot_data) == 168 and int(px.included.sum()) == 28
assert len(relative) == 132 and n == 11
assert np.allclose(plot_data.display_x, gx(plot_data.midpoint))
assert plot_data.percent_of_sgCTRL.between(0, DNA_TOP).all()
assert px.loc[px.included, 'percent_of_sgCTRL'].between(0, DNA_TOP).all()
assert relative.relative_expression_percent.between(0, 185).all()
assert np.allclose(relative.groupby(['gene_name', 'condition']).relative_expression_percent.mean().unstack().loc[genes], rna_summary.loc[genes, ['LOY', 'sgCTRL']])
xlim = (-.15, gx(length)+.15)
ax_a.set_xlim(xlim)
ax_b.set_xlim(xlim)

# Regional drawing block below preserves the adopted atlas geometry and marks.
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
# No source tables are rewritten for this draft.
ax_a.text(gx((q0+p2)/2),DNA_TOP*.48718,'Yq12\nexcluded',rotation=90,
    ha='center',va='center',fontsize=6.8,color='#666666',zorder=6)
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
ax_b.set_xlabel('HG002 T2T-Y position (Mb; Yq12 compressed)', fontsize=7.4, labelpad=4)
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

# Verify coordinate alignment and visible source measurements before exporting.
assert np.allclose(ax_a.get_xlim(), ax_b.get_xlim())
assert np.isclose(ax_a.get_position().x0, ax_b.get_position().x0)
assert np.isclose(ax_a.get_position().width, ax_b.get_position().width)
assert np.isclose(window_ax.get_position().height, arm_ax.get_position().height)
for ax in [ax_a, ax_b, ax_c]:
    ax.tick_params(labelsize=7, pad=2)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
clipped_text = []
from matplotlib.text import Text
for item in fig.findobj(Text):
    if not item.get_visible() or not item.get_text().strip():
        continue
    bounds = item.get_window_extent(renderer)
    if bounds.x0 < -1 or bounds.y0 < -1 or bounds.x1 > fig.bbox.width+1 or bounds.y1 > fig.bbox.height+1:
        clipped_text.append(item.get_text())
assert not clipped_text, clipped_text
for ext in ['pdf', 'svg', 'png']:
    fig.savefig(OUT / f'{STEM}.{ext}', dpi=600, facecolor='white',
                metadata={'Title': 'main draft - HT1080 genome-wide and regional Y characterisation'})
preview = OUT / f'{STEM}_preview.png'
fig.savefig(preview, dpi=220, facecolor='white')
hashed_preview = OUT / f'{STEM}_preview_{sha(preview)[:12]}.png'
shutil.copy2(preview, hashed_preview)
plt.close(fig)
assert protected == {str(p.relative_to(ROOT)): sha(p) for p in protected_paths}
assert all(sha(input_sources[name]) == details['sha256'] == sha(SNAP / name)
           for name, details in input_manifest.items())
validation = {
    'status': 'Draft only; no manuscript or selected-figure changes',
    'size_inches': [WIDTH, HEIGHT], 'font': 'Open Sans', 'dpi': 600,
    'panel_map': {'a': 'Fig 3.3a, inset removed', 'b': 'Fig 3.3b, unlogged fold changes', 'c': 'Fig 3.3c, unlogged fold changes',
                  'd': 'Fig 3.4a', 'e': 'Fig 3.4b', 'f': 'Fig 3.4c'},
    'counts': {'genome_5Mb_windows': len(bins), 'DNA_RNA_1Mb_windows': len(windows),
               'DNA_RNA_arms': len(arms), 'regional_Y_bins': len(plot_data),
               'X_PAR_bins': int(px.included.sum()), 'RNA_observations': len(relative),
               'RNA_genes': n},
    'correlations_from_saved_table': correlations.reset_index().to_dict(orient='records'),
    'correlations_displayed_on_unlogged_values': linear_correlations.reset_index().to_dict(orient='records'),
    'originals_and_manuscript_unchanged': True,
    'prior_combined_draft_preserved': True,
    'shared_scatter_RNA_axis': list(RNA_YLIM),
    'shared_scatter_DNA_axis': list(DNA_XLIM),
    'scatter_axes': 'Both axes linear in b and c; 1 denotes no change.',
    'display_transform': 'DNA_fold_change = 2**DNA; RNA_fold_change = 2**RNA; source log2 summaries retained.',
    'descriptive_correlation_change': 'Pearson correlations recalculated from the displayed unlogged fold changes; no model refitting or new inferential test.',
    'scatter_coordinates_verified_against_transformed_source': True,
    'no_values_clipped': True, 'no_text_outside_canvas': True,
    'no_refitting_or_re_filtering': True,
    'input_manifest': input_manifest, 'protected_sha256': protected,
    'renderer_sha256': sha(Path(__file__)),
    'export_sha256': {ext: sha(OUT / f'{STEM}.{ext}') for ext in ['pdf', 'svg', 'png']},
    'preview': str(hashed_preview)
}
(SOURCE / 'validation.json').write_text(json.dumps(validation, indent=2)+'\n', encoding='utf-8')
(SOURCE / 'README.md').write_text('''# HT1080 combined genome-wide/regional figure draft

Separate review draft combining current Figures 3.3 and 3.4. The Y inset from 3.3a is omitted. At the author's request, panels b and c now plot 2 raised to each saved log2 summary on ordinary linear axes; 1 is unchanged relative to sgCTRL. These are back-transformed geometric fold-change summaries, not arithmetic ratios newly estimated from raw counts. The displayed Pearson correlations are recalculated from these plotted values, including and excluding Y, because Pearson correlation is not invariant to exponentiation. No underlying model, gene set, filtering or measurement is changed. Other panels retain their previous presentation. No manuscript, selected figure or figure number has been changed; the first combined draft is preserved.

The inputs folder freezes the existing display tables and cytobands; validation.json records source hashes, retained counts and protected-file hashes. The regional data are the adopted Y-PAR-masked October 28 atlas, including combined X-mapped PAR coverage. Genome-wide summaries retain their existing source and filtering. The two data preparations must remain distinguished in the caption/Methods.

Reproduce with figures/chapter3/code/draft_ht1080_genomewide_regional_combined_linear_20261030.py using the existing local plotting environment. Review exports and the draft caption are in figures/chapter3/review/ht1080_genomewide_regional_combined_20261030_v2_linear.
''', encoding='utf-8')
print(json.dumps({'pdf': str(OUT/f'{STEM}.pdf'), 'preview': str(hashed_preview),
                  'counts': validation['counts'], 'size_inches': [WIDTH, HEIGHT]}))

