"""ch03_main_du145_rna_sex_chromosome: Y expression against lpWGS segments, regions, PAR and X paralogs."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
from du145_rna_common import (read, T, LPWGS, COND, COND_LABEL, INTERVAL, DARK, GREY, BLUE, RED, MM, panel, export)

STEM = 'ch03_main_du145_rna_sex_chromosome'
y = read('Y_PAR_gene_table')
reg = read('Y_PAR_region_abundance')
par = read('XY_paralog_effects')
FLOOR = 0.001
PAR1_END, YQ_START, YQ_END = 2.45832, 10.883085, 27.449937

fig = plt.figure(figsize=(183 * MM, 180 * MM))
a = fig.add_axes([0.075, 0.745, 0.885, 0.195])
b = fig.add_axes([0.075, 0.435, 0.23, 0.135])
c = fig.add_axes([0.375, 0.435, 0.255, 0.135])
d = fig.add_axes([0.735, 0.435, 0.225, 0.135])
e = fig.add_axes([0.075, 0.115, 0.885, 0.175])

# ---- a: positional RNA (percent of sgCTRL) with the lpWGS segment ratio underneath.
shade = [('PAR1', 0, PAR1_END), ('retained distal Yp', PAR1_END, 4.8), ('depleted MSY', 4.8, YQ_END)]
for name, x0, x1 in shade:
    a.axvspan(x0, x1, color=INTERVAL[name], alpha=.10, lw=0)
    # Region names sit just above the plot area so they never compete with gene labels.
    a.text((x0 + x1) / 2, 1.012, name, transform=a.get_xaxis_transform(), ha='center', va='bottom', fontsize=6.1, color=DARK)
a.text(PAR1_END / 2, 1.0, 'PAR genes:\npanel d', ha='center', va='center', fontsize=5.8, color=DARK)
# lpWGS: saved QDNAseq segment values along chrY (no resegmentation), intersected between specimens.
bins = pd.read_csv(LPWGS / 'qdnaseq_bins.csv')
CTRL_MIN = 0.05   # control segments below this are floor or low-mappability ampliconic intervals; ratio not meaningful


def saved_segments(sample):
    d = bins[(bins['sample'] == sample) & bins['chromosome'].eq('chrY') & bins['use'] & bins['log2_seg'].notna()].sort_values('start')
    d = d.assign(seg=2 ** d.log2_seg)
    run = (d.seg.round(8) != d.seg.round(8).shift()).cumsum()
    return d.groupby(run).agg(start_bp=('start', 'min'), end_bp=('end', 'max'), rel_cn=('seg', 'first')).reset_index(drop=True)


loy, ctl = saved_segments('LNCaP-LOY-2X'), saved_segments('LNCaP-sgCTRL')
dna_rows = []
for _, L in loy.iterrows():
    for _, C in ctl.iterrows():
        s, en = max(L.start_bp, C.start_bp) / 1e6, min(L.end_bp, C.end_bp, YQ_END * 1e6) / 1e6
        if s < en and C.rel_cn >= CTRL_MIN:
            dna_rows.append(dict(start_Mb=s, end_Mb=en, loy_rel_cn=L.rel_cn, ctrl_rel_cn=C.rel_cn,
                                 ratio=100 * max(L.rel_cn, FLOOR) / C.rel_cn, loy_at_floor=L.rel_cn <= FLOOR * 1.1))
dna = pd.DataFrame(dna_rows)
# Ratios below 1% (floor-limited LOY segments) are drawn at the <=1% line, dashed.
dna['plotted_ratio'] = dna.ratio.clip(lower=1)
dna['drawn_at_floor_line'] = dna.ratio < 1
for r in dna.itertuples():
    a.plot([r.start_Mb, r.end_Mb], [r.plotted_ratio] * 2, color='#333333', lw=1.6, solid_capstyle='butt', zorder=4,
           ls=(0, (2.2, 1.2)) if r.drawn_at_floor_line else '-')
# Male-specific genes only: PAR transcripts are chrX-assigned, cannot be attributed to Y, and are shown in panel d.
shown = y[y.expressed & y.t2t_mid_Mb.notna() & (y.CTRL_mean_CPM >= 1) & (y.t2t_mid_Mb < YQ_END) & y.region.eq('MSY')].copy()
pair_cols = ['percent_of_sgCTRL_P12', 'percent_of_sgCTRL_P13', 'percent_of_sgCTRL_P14']
for col in pair_cols:
    a.scatter(shown.t2t_mid_Mb, shown[col].clip(lower=1), s=5, color=RED, alpha=.40, lw=0, zorder=5)
a.scatter(shown.t2t_mid_Mb, shown.percent_of_sgCTRL.clip(lower=1), s=11, facecolors='white',
          edgecolors=RED, linewidths=.7, zorder=6)
a.axhline(100, ls='--', lw=.7, color='#999999', zorder=2)
# Set the scale and limits before adjusting labels so they stay inside the axes.
a.set_yscale('log'); a.set_ylim(.6, 180); a.set_xlim(-.6, YQ_END + .6)
texts = [a.text(r.t2t_mid_Mb, max(r.percent_of_sgCTRL, 1), r.gene_name, fontsize=6.4, color=DARK)
         for r in shown.itertuples() if r.region == 'MSY' and not r.gene_name.startswith('ENSG')]
adjust_text(texts, ax=a, only_move={'text': 'xy'}, arrowprops=dict(arrowstyle='-', color='#AAAAAA', lw=.45), force_text=(.3, .5))
a.set_yticks([1, 3, 10, 30, 100]); a.set_yticklabels(['≤1', '3', '10', '30', '100'])
from matplotlib.ticker import LogLocator, NullFormatter
a.yaxis.set_minor_locator(LogLocator(base=10, subs=range(2, 10)))
a.yaxis.set_minor_formatter(NullFormatter())
a.set_xlabel('CHM13v2.0 chromosome Y position (Mb)')
a.set_ylabel('Y-depleted, % of sgCTRL')
a.legend(handles=[Line2D([], [], ls='', marker='o', mfc='white', mec=RED, mew=.7, ms=3.4, label='RNA, group mean'),
                  Line2D([], [], ls='', marker='o', color=RED, alpha=.4, ms=2.5, label='RNA, each passage pair'),
                  Line2D([], [], color='#333333', lw=1.6, label='lpWGS segment ratio'),
                  Line2D([], [], color='#333333', lw=1.6, ls=(0, (2.2, 1.2)), label='lpWGS, ≤1% (drawn at 1%)')],
         loc='upper center', bbox_to_anchor=(.5, .672), bbox_transform=fig.transFigure, ncol=4,
         fontsize=6.3, handletextpad=.4, columnspacing=1.6, borderaxespad=0, borderpad=.1)
panel(fig, a, 'a', 'Expression along chromosome Y follows the retained and deleted segments', dx=-.058, dy=.03)

# ---- b: regional summed abundance.
order = ['PAR1', 'MSY retained distal Yp', 'MSY depleted', 'PAR2']
rs = reg[reg.region.isin(order)].groupby(['region', 'treatment']).sum_CPM.mean().unstack()
rs = rs.loc[order]
xb = np.arange(len(order))
b.bar(xb - .19, rs['CTRL'], width=.36, color=BLUE, label=COND_LABEL['CTRL'])
b.bar(xb + .19, rs['LOY'], width=.36, color=RED, label=COND_LABEL['LOY'])
for i, r in enumerate(order):
    b.text(i, max(rs.loc[r]) * 1.9, f"{100 * rs.loc[r, 'LOY'] / rs.loc[r, 'CTRL']:.0f}%", ha='center', fontsize=6.6, color=DARK)
b.set_yscale('log'); b.set_ylim(3, 1e4)
b.set_xticks(xb, ['PAR1', 'Retained\nYp', 'Depleted\nMSY', 'PAR2'], fontsize=6.0)
b.set_ylabel('Summed CPM'); b.legend(fontsize=6.4, loc='upper right', handlelength=1, borderaxespad=.2)
panel(fig, b, 'b', 'Regional abundance', dx=-.062, y=.605)

# ---- c: per-gene percent for MSY genes, split by interval.
msy = y[y.expressed & (y.region == 'MSY') & (y.CTRL_mean_CPM >= 1) & y.gene_name.notna()].copy()
msy = msy[~msy.gene_name.str.startswith('ENSG')].sort_values('t2t_mid_Mb')
msy['grp'] = np.where(msy.lpWGS_interval.str.startswith('retained'), 'retained distal Yp', 'depleted MSY')
xc = np.arange(len(msy))
c.bar(xc, msy.percent_of_sgCTRL.clip(lower=.4), color=[INTERVAL[g] for g in msy.grp], width=.7)
for i, r in enumerate(msy.itertuples()):
    c.text(i, max(r.percent_of_sgCTRL, .4) * 1.12, f'{r.percent_of_sgCTRL:.0f}', ha='center', fontsize=6, color=DARK)
c.set_yscale('log'); c.set_ylim(.35, 260); c.set_yticks([1, 10, 100]); c.set_yticklabels(['1', '10', '100'])
c.set_xticks(xc, msy.gene_name, rotation=90, fontsize=6.4)
c.set_ylabel('% of sgCTRL')
c.axhline(100, ls='--', lw=.7, color='#999999')
panel(fig, c, 'c', 'Male-specific Y genes', dx=-.058, y=.605)

# ---- d: PAR1 vs PAR2 per gene.
pr = y[y.expressed & y.PAR & (y.CTRL_mean_CPM >= 5)].copy().sort_values(['region', 'percent_of_sgCTRL'])
xd = np.arange(len(pr))
d.bar(xd, pr.percent_of_sgCTRL, color=[INTERVAL[r] for r in pr.region], width=.7)
d.axhline(100, ls='--', lw=.7, color='#999999')
d.set_xticks(xd, pr.gene_name, rotation=90, fontsize=6.4)
d.set_ylabel('% of sgCTRL'); d.set_ylim(0, 145)
d.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=INTERVAL[k]) for k in ['PAR1', 'PAR2']],
         labels=['PAR1', 'PAR2'], fontsize=6.4, loc='upper left', handlelength=.9, borderaxespad=.2)
panel(fig, d, 'd', 'Pseudoautosomal genes', dx=-.062, y=.605)

# ---- e: X/Y paralog pairs.
pp = par[par.Y_expressed_in_sgCTRL].copy()
pp['label'] = pp.X + ' / ' + pp.Y
pp = pp.sort_values('Y_mean_difference')
xe = np.arange(len(pp))
e.axhline(0, color='#999999', lw=.7)
e.vlines(xe, pp.Y_mean_difference, pp.X_mean_difference, color=GREY, lw=.8, zorder=1)
e.scatter(xe, pp.Y_mean_difference, s=34, color=RED, label='Y-linked', zorder=3)
e.scatter(xe, pp.X_mean_difference, s=34, color=BLUE, label='X-linked paralog', zorder=3)
e.axhspan(0.138, 0.165, color=BLUE, alpha=.12, lw=0, zorder=0)
e.text(len(pp) - .55, .32, 'chrX arm-level RNA shift', fontsize=6.2, color=DARK, va='center', ha='right')
e.set_xticks(xe, pp.label, rotation=28, ha='right', fontsize=6.6)
e.set_ylabel('Mean log$_2$ CPM difference\n(Y-depleted − sgCTRL)')
e.set_xlim(-.6, len(pp) - .4)
e.legend(fontsize=6.5, loc='lower right', borderaxespad=.3)
panel(fig, e, 'e', 'X and Y paralog pairs', dx=-.058, y=.315)

export(fig, STEM,
       {'panel_a_rna_positions': shown[['gene_name', 'region', 't2t_mid_Mb', 'lpWGS_interval', 'CTRL_mean_CPM', 'percent_of_sgCTRL'] + pair_cols],
        'panel_a_lpwgs_segment_ratio': dna,
        'panel_b_region_sums': rs.reset_index(),
        'panel_c_msy_genes': msy[['gene_name', 't2t_mid_Mb', 'grp', 'CTRL_mean_CPM', 'percent_of_sgCTRL', 'logFC', 'BH_Y_PAR_family']],
        'panel_d_par_genes': pr[['gene_name', 'region', 'CTRL_mean_CPM', 'percent_of_sgCTRL', 'logFC', 'BH_Y_PAR_family']],
        'panel_e_paralogs': pp[['X', 'Y', 'X_mean_difference', 'X_BH_inclusive', 'Y_mean_difference', 'Y_BH_inclusive', 'summed_XY_percent_of_sgCTRL']]},
       {'inputs': [T / 'Y_PAR_gene_table.csv', T / 'Y_PAR_region_abundance.csv', T / 'XY_paralog_effects.csv',
                   LPWGS / 'qdnaseq_bins.csv'],
        'notes': ('Percent = 100 * LOY mean CPM / sgCTRL mean CPM on the primary TMM scaling. Panel a shows genes with sgCTRL '
                  'mean CPM >= 1 (male-specific genes only); values at or below 1% are drawn at 1%. lpWGS ratios intersect the saved '
                  'QDNAseq segment values of LNCaP-LOY-2X and LNCaP-sgCTRL along chrY to 27.45 Mb, with the upstream 0.001 '
                  'floor; intervals whose control segment is below 0.05 relative copy number (floor or low-mappability '
                  'ampliconic) are omitted, and ratios below 1% are drawn dashed at 1%. RNA and lpWGS specimens were not matched. PAR features are '
                  'chrX-assigned and cannot be apportioned between X and Y. Panel e shading is the chrX arm-level RNA range '
                  '(Xq +0.138, Xp +0.165) from CN_arm_level_RNA_vs_DNA.csv.')})
