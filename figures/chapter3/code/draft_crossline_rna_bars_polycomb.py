"""ch03_main_crossline_rna_comparison: shared and background-specific responses across HT1080, A549, DU145."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from crossline_common import (ct, CT, LINES, LINE_COLOUR, THEME_LABEL, GROUP_LABEL, BLUE, RED, GREY, DARK, MM, panel, export)

STEM = 'ch03_draft_crossline_rna_comparison_bars_polycomb'
SIGNED = LinearSegmentedColormap.from_list('signed', [BLUE, '#F7F7F7', RED])
y = ct('Y_inventory_by_line'); reg = ct('Y_PAR_region_sums_by_line'); xp = ct('XY_paralogs_by_line')
agr = ct('gene_level_three_way_agreement').iloc[0]; cor = ct('gene_level_correlations')
gs = ct('gene_sets_shared_divergent'); ext = ct('external_replication_sets'); reg_used = ct('external_dataset_registry_used')

fig = plt.figure(figsize=(183 * MM, 232 * MM))
a = fig.add_axes([0.075, 0.7761, 0.5, 0.1691])
b = fig.add_axes([0.715, 0.7761, 0.245, 0.1691])
c1 = fig.add_axes([0.075, 0.502, 0.12, 0.1691])
c2 = fig.add_axes([0.285, 0.502, 0.12, 0.1691])
d = fig.add_axes([0.66, 0.502, 0.2, 0.1691])
e = fig.add_axes([0.215, 0.0457, 0.62, 0.3284])

# ---- a: Y inventory, percent of control on a log axis.
genes = ['RPS4Y1', 'ZFY', 'PRKY', 'USP9Y', 'DDX3Y', 'UTY', 'TXLNGY', 'KDM5D', 'EIF1AY']
FLOOR = 0.1
xa = np.arange(len(genes))
for k, L in enumerate(LINES):
    v = y[y.line.eq(L)].set_index('gene').reindex(genes).percent_of_control
    a.bar(xa + (k - 1) * .25, v, width=.23, color=LINE_COLOUR[L], label=L, zorder=3)
    for xx, vv in zip(xa + (k - 1) * .25, v):
        if np.isfinite(vv):
            a.text(xx, vv + 2, f'{vv:.1f}' if vv < 10 else f'{vv:.0f}', ha='center', va='bottom', rotation=90, fontsize=5)
a.axhline(100, ls='--', lw=.7, color='#999999')
a.axvspan(-.5, 1.5, color='#00A087', alpha=.06, lw=0)
a.text(.5, 117, 'distal Yp', ha='center', fontsize=6.2, color=DARK)
a.set_ylim(0, 150)
a.set_yticks([0,25,50,75,100])
a.set_xticks(xa, genes, rotation=45, ha='right', fontsize=6.6); a.set_xlim(-.6, len(genes) - .4)
a.set_ylabel('Y-depleted, % of sgCTRL')
par = reg[reg.region.eq('PAR1')].set_index('line').percent_of_control
a.text(len(genes) - .45, 114, 'PAR1 summed:  ' + '   '.join(f'{L} {par[L]:.0f}%' for L in LINES), ha='right', va='top', fontsize=6.2, color=DARK)
a.legend(fontsize=6.4, loc='upper right', ncol=3, borderaxespad=.2, handletextpad=.2, columnspacing=.9)
panel(fig, a, 'a', '', dx=-.058, y=0.9589)

# ---- b: X paralogs (expressed in every line), with each line's overall chrX shift.
ok = xp.groupby('X').X_CTRL_CPM.min() >= 1
pairs = [p for p in ['DDX3X', 'EIF1AX', 'KDM5C', 'KDM6A', 'USP9X', 'ZFX', 'RPS4X', 'PRKX', 'TMSB4X'] if ok.get(p, False)]
yb = np.arange(len(pairs))[::-1]
for k, L in enumerate(LINES):
    v = xp[xp.line.eq(L)].set_index('X').reindex(pairs).X_mean_difference
    b.barh(yb + (1 - k) * .25, v, height=.23, color=LINE_COLOUR[L], zorder=3)
    shift = xp[xp.line.eq(L)].chrX_median_shift.iloc[0]
    b.axvline(shift, color=LINE_COLOUR[L], lw=.9, ls=':', zorder=1)
b.axvline(0, color='#999999', lw=.6)
b.set_yticks(yb, pairs, fontsize=6.4); b.set_xlim(-.45, .45)
b.set_xlabel('X paralog, mean log$_2$ CPM difference', fontsize=7)
panel(fig, b, 'b', '', dx=-.075, y=0.9589)

# ---- c: gene-level concordance, magnitude and direction separately.
c1.bar([0, 1], [agr.genes_large_in_all_three, agr.null_large_mean], color=[DARK, GREY], width=.6)
c1.errorbar(1, agr.null_large_mean, yerr=[[agr.null_large_mean - agr['null_large_2.5']], [agr['null_large_97.5'] - agr.null_large_mean]],
            color='#555555', capsize=2, lw=.8)
c1.set_xticks([0, 1], ['Observed', 'Shuffled'], fontsize=6.5); c1.set_ylabel('Genes with |effect| ≥ 0.5\nin all three lines', fontsize=7)
c1.text(.76, .95, f'P = {agr.empirical_P_large_cooccurrence:.3f}', transform=c1.transAxes, ha='center', va='top', fontsize=6.3, color=DARK)
frac = 100 * agr.fraction_same_sign_given_large
c2.bar([0], [frac], color=DARK, width=.55)
c2.errorbar(0, frac, yerr=[[frac - 100 * agr.binomial_CI_low], [100 * agr.binomial_CI_high - frac]], color='#555555', capsize=2, lw=.8)
c2.axhline(25, ls='--', lw=.8, color=RED); c2.text(.42, 25.8, 'chance', color=RED, fontsize=6.1)
c2.set_xlim(-.6, .75); c2.set_ylim(0, 34); c2.set_xticks([0], [f'{int(agr.same_sign_all_three)} of {int(agr.genes_large_in_all_three)}'], fontsize=6.5)
c2.set_ylabel('Same sign in all three (%)', fontsize=7)
c2.text(.5, .95, f'P = {agr.binomial_P_sign_agreement:.3f}', transform=c2.transAxes, ha='center', va='top', fontsize=6.3, color=DARK)
rho = '  '.join(f"{r.line_1}–{r.line_2} {r.Spearman:+.2f}" for r in cor.itertuples())
fig.text(.075, 0.4627, 'Genome-wide Spearman:  ' + rho, fontsize=6.2, color=DARK)
panel(fig, c1, 'c', '', dx=-.058, y=0.7003)

# ---- d: named sets instead of medians over mixed themes (discussion draft).
selected = {
 'HALLMARK_INTERFERON_ALPHA_RESPONSE': 'IFN-alpha response',
 'HALLMARK_INTERFERON_GAMMA_RESPONSE': 'IFN-gamma response',
 'HALLMARK_INFLAMMATORY_RESPONSE': 'Inflammatory response',
 'ONDER_CDH1_TARGETS_2_DN': 'CDH1-loss: down genes',
 'REACTOME_RRNA_PROCESSING': 'rRNA processing',
 'REACTOME_GPCR_LIGAND_BINDING': 'GPCR ligand binding',
 'REACTOME_M_PHASE': 'M phase',
 'REACTOME_MITOTIC_PROPHASE': 'Mitotic prophase',
 'REACTOME_RHO_GTPASE_EFFECTORS': 'RHO GTPase effectors',
 'MIKKELSEN_MCV6_HCP_WITH_H3K27ME3': 'H3K27me3: MCV6',
 'MIKKELSEN_MEF_HCP_WITH_H3K27ME3': 'H3K27me3: MEF',
}
dat = gs.set_index('pathway').loc[list(selected)]
med = dat[[f'NES_{L}' for L in LINES]].to_numpy()
im = d.imshow(med, cmap=SIGNED, norm=TwoSlopeNorm(0, -2.2, 2.2), aspect='auto')
for i in range(med.shape[0]):
    for j,L in enumerate(LINES):
        sig = dat.iloc[i][f'BH_{L}'] < .05
        d.text(j, i, f'{med[i,j]:+.1f}' + ('*' if sig else ''), ha='center', va='center', fontsize=5.6, color='#222222')
d.set_xticks(range(3), LINES, fontsize=6.5); d.xaxis.tick_top()
d.set_yticks(range(len(dat)), list(selected.values()), fontsize=6.2)
for spine in d.spines.values(): spine.set_visible(False)
d.tick_params(length=0)
cb = fig.colorbar(im, ax=d, fraction=.07, pad=.04); cb.set_label('NES (* BH < 0.05)', fontsize=6.5); cb.ax.tick_params(labelsize=6)
panel(fig, d, 'd', '', dx=-.27, y=0.7003)

# ---- e: replication across in-house and external datasets.
from pathlib import Path
EXTRA = Path(__file__).resolve().parents[3] / 'research/molecular_characterisation/crossline_figure_review_2026-09-28/external_replication_sets.csv'
ext = pd.concat([ext, pd.read_csv(EXTRA)], ignore_index=True)
cols = ['HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INFLAMMATORY_RESPONSE',
        'ONDER_CDH1_TARGETS_2_DN', 'BMI1_DN.V1_UP', 'MEL18_DN.V1_UP', 'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_E2F_TARGETS', 'MIKKELSEN_MCV6_HCP_WITH_H3K27ME3', 'MIKKELSEN_MEF_HCP_WITH_H3K27ME3', 'BENPORATH_PRC2_TARGETS', 'BENPORATH_ES_WITH_H3K27ME3']
col_lab = ['IFN-α', 'IFN-γ', 'Inflammatory', 'CDH1 targets', 'BMI1 KD: up', 'MEL18 KD: up', 'G2/M', 'E2F', 'H3K27me3: MCV6', 'H3K27me3: MEF', 'PRC2: hESC', 'H3K27me3: hESC']
short = {'HT1080': 'HT1080', 'A549': 'A549', 'DU145': 'DU145', 'A549_clones_LOYvROY': 'A549 LOY clones',
         'A549_multiome_pb_LOYvROY': 'A549 clones, multiome*', 'GSE313572_HT1080_ctrl_YnegvYpos': 'HT1080 spheroids', 'Tahoe_A549_DMSO_YnegvYpos': 'A549 (Tahoe)', 'Tahoe_J82_DMSO_YnegvYpos': 'J82 (Tahoe)', 'Tahoe_SKMEL2_DMSO_YnegvYpos': 'SK-MEL-2 (Tahoe)', 'Tahoe_HCT15_DMSO_YnegvYpos': 'HCT15 (Tahoe, 6 wells)', 'MB49_cellline_YnegvYpos': 'MB49 cells (mouse)',
         'MB49_tumor_YnegvYpos': 'MB49 tumours (mouse)', 'GSE126712_lymph_XOvXY': 'Lymphocytes', 'GSE207112_troph_XOvXY': 'iPSC trophoblast',
         'GSE264742_NCC_XOvXY': 'iPSC neural crest', 'ERP156313_fetal_XOvXY': 'Fetal tissues', 'PageLab_Fib_X0vXY': 'Fibroblasts',
         'PageLab_LCL_X0vXY': 'Lymphoblastoid lines', 'GSE301297_undiff_UTYKOvWT': 'hESC UTY-KO', 'GSE301297_undiff_DKOvWT': 'hESC UTX/UTY-DKO',
         'XYY_integrated_flipped': 'XYY vs XY (flipped)'}
rows = [r for r in reg_used.dataset if r != 'XYY_integrated_flipped']
A = ext.pivot(index='dataset', columns='gene_set', values='AUC').reindex(index=rows, columns=cols)
P = ext.pivot(index='dataset', columns='gene_set', values='empirical_P').reindex(index=rows, columns=cols)
im2 = e.imshow(A.to_numpy() - .5, cmap=SIGNED, norm=TwoSlopeNorm(0, -.3, .3), aspect='auto')
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        if P.iat[i, j] < .05:
            e.text(j, i, '•', ha='center', va='center', fontsize=7, color='#222222')
e.set_xticks(range(len(cols)), col_lab, fontsize=6.2, rotation=55, ha='left'); e.xaxis.tick_top()
e.set_yticks(range(len(rows)), [short[r] for r in rows], fontsize=6.3)
e.tick_params(length=0)
for s in e.spines.values():
    s.set_visible(False)
grp_of = dict(zip(reg_used.dataset, reg_used.group))
bounds = [i for i in range(1, len(rows)) if grp_of[rows[i]] != grp_of[rows[i - 1]]]
for bnd in bounds:
    e.axhline(bnd - .5, color='white', lw=2.2)
starts = [0] + bounds; ends = bounds + [len(rows)]
for s0, s1 in zip(starts, ends):
    e.text(len(cols) - .4, (s0 + s1 - 1) / 2, GROUP_LABEL[grp_of[rows[s0]]], va='center', ha='left', fontsize=6.3, color=DARK)
e.axvline(5.5, color='white', lw=2.2)
e.axvline(7.5, color='white', lw=2.2)
cb2 = fig.colorbar(im2, ax=e, fraction=.035, pad=.2); cb2.set_label('AUC − 0.5  (blue: lower with Y absent)', fontsize=6.3)
cb2.ax.tick_params(labelsize=6)
panel(fig, e, 'e', '', dx=-.14, y=0.4261)

export(fig, STEM,
       {'panel_a_Y_inventory': y, 'panel_a_region_sums': reg, 'panel_b_XY_paralogs': xp[xp.X.isin(pairs)],
        'panel_c_agreement': ct('gene_level_three_way_agreement'), 'panel_c_correlations': cor,
        'panel_d_named_sets': dat.reset_index(),
        'panel_e_replication_AUC': A.reset_index(), 'panel_e_replication_P': P.reset_index()},
       {'inputs': [EXTRA] + [CT / f'{n}.csv' for n in ['Y_inventory_by_line', 'Y_PAR_region_sums_by_line', 'XY_paralogs_by_line', 'gene_level_three_way_agreement',
                                           'gene_level_correlations', 'gene_sets_shared_divergent', 'external_replication_sets', 'external_dataset_registry_used']],
        'notes': ('Percent of control is the arithmetic ratio of mean CPM per sampling period, averaged across periods. Panel b shows X paralogs with '
                  'sgCTRL CPM >= 1 in every line. Panel c: magnitude co-occurrence against 1,000 within-line shuffles; sign agreement by exact '
                  'binomial against 1/4. Draft: panels a and b use bars, panel a linear percent scale; panel d displays 11 named illustrative sets with BH < 0.05 marked by asterisks; the complete set list is retained in the supplementary atlas and review table. Four additional Polycomb/H3K27me3 signatures in panel e are exploratory, using the original 500-permutation procedure. Panel e: AUC of member-gene effects against non-members, '
                  'dots = empirical two-sided P < 0.05 from 500 size-matched random sets; * the multiome uses the same clones as the A549 LOY clones.')})
