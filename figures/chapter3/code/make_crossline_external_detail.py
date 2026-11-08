"""ch03_supp_crossline_external_detail: PXDN across datasets, and how each dataset's polarity was established."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from crossline_common import (ct, CT, CROSS, LINE_COLOUR, GROUP_LABEL, BLUE, RED, GREY, DARK, MM, panel, export)

STEM = 'ch03_supp_crossline_external_detail'
focus = ct('external_focus_genes'); pol = ct('external_polarity_check'); used = ct('external_dataset_registry_used')
short = {'HT1080': 'HT1080', 'A549': 'A549', 'DU145': 'DU145', 'A549_clones_LOYvROY': 'A549 LOY clones',
         'A549_multiome_pb_LOYvROY': 'A549 clones, multiome', 'GSE313572_HT1080_ctrl_YnegvYpos': 'HT1080 spheroids', 'Tahoe_A549_DMSO_YnegvYpos': 'A549 (Tahoe)', 'Tahoe_J82_DMSO_YnegvYpos': 'J82 (Tahoe)', 'Tahoe_SKMEL2_DMSO_YnegvYpos': 'SK-MEL-2 (Tahoe)', 'Tahoe_HCT15_DMSO_YnegvYpos': 'HCT15 (Tahoe, 6 wells)', 'MB49_cellline_YnegvYpos': 'MB49 cells (mouse)',
         'MB49_tumor_YnegvYpos': 'MB49 tumours (mouse)', 'GSE126712_lymph_XOvXY': 'Lymphocytes', 'GSE207112_troph_XOvXY': 'iPSC trophoblast',
         'GSE264742_NCC_XOvXY': 'iPSC neural crest', 'ERP156313_fetal_XOvXY': 'Fetal tissues', 'PageLab_Fib_X0vXY': 'Fibroblasts',
         'PageLab_LCL_X0vXY': 'Lymphoblastoid lines', 'GSE301297_undiff_UTYKOvWT': 'hESC UTY-KO', 'GSE301297_undiff_DKOvWT': 'hESC UTX/UTY-DKO',
         'XYY_integrated_flipped': 'XYY vs XY (flipped)'}
rows = list(used.dataset)
grp = dict(zip(used.dataset, used.group))
GCOL = {'in_house': '#555555', 'A_somatic_LOY': '#E64B35', 'A2_single_cell_mosaic': '#F39B7F', 'B_45X_vs_46XY': '#3C5488', 'C_not_LOY': '#999999'}

fig = plt.figure(figsize=(183 * MM, 122 * MM))
a = fig.add_axes([0.20, 0.12, 0.33, 0.78])
b = fig.add_axes([0.70, 0.12, 0.27, 0.78], sharey=a)

f = focus[focus.gene.eq('PXDN')].set_index('dataset').reindex(rows)
ya = np.arange(len(rows))[::-1]
a.barh(ya, f.effect_Yloss_polarity.clip(-3, 3), color=[GCOL[grp[r]] for r in rows], height=.66)
for yv, r in zip(ya, rows):
    v, q = f.loc[r, 'effect_Yloss_polarity'], f.loc[r, 'FDR']
    if pd.isna(v):
        a.text(0.05, yv, 'not tested', va='center', fontsize=5.8, color=DARK)
    else:
        lab = f'{v:+.2f}' + (f'  (FDR {q:.2g})' if pd.notna(q) else '')
        if abs(v) > 3:   # clipped bar: label inside, in white
            a.text(np.sign(v) * 2.92, yv, lab + ', clipped', va='center', ha='left' if v < 0 else 'right', fontsize=5.4, color='white')
        else:
            a.text(v + (.08 if v >= 0 else -.08), yv, lab, va='center', ha='left' if v >= 0 else 'right', fontsize=5.6, color=DARK)
a.axvline(0, color='#999999', lw=.6)
a.set_yticks(ya, [short[r] for r in rows], fontsize=6.4); a.set_xlim(-4.2, 2.2)
a.set_xlabel('PXDN effect, Y-loss polarity (log$_2$; clipped at ±3)', fontsize=7)
handles = [plt.Rectangle((0, 0), 1, 1, color=GCOL[g]) for g in GCOL]
a.legend(handles, [GROUP_LABEL[g] for g in GCOL], fontsize=5.9, loc='lower left', borderaxespad=.3, handlelength=.8)
panel(fig, a, 'a', '', dx=-.17)

p = pol.set_index('dataset').reindex(rows)
b.set_xlim(0, 1); b.set_ylim(-.6, len(rows) - .4); b.axis('off')
for yv, r in zip(ya, rows):
    ok = p.loc[r, 'status'] == 'verified_by_Y_gene'
    b.scatter(.03, yv, s=22, color='#00A087' if ok else '#DDDDDD', edgecolors='#777777', linewidths=.4)
    text = (f"{p.loc[r, 'check_gene']} {p.loc[r, 'effect_Yloss_polarity']:+.2f}" if ok else 'Y genes absent; contrast label')
    b.text(.1, yv, text, va='center', fontsize=6, color=DARK)
panel(fig, b, 'b', '', dx=-.05)
export(fig, STEM,
       {'panel_a_PXDN': f.reset_index(), 'panel_b_polarity': p.reset_index()},
       {'inputs': [CT / 'external_focus_genes.csv', CT / 'external_polarity_check.csv', CT / 'external_dataset_registry_used.csv',
                   CROSS / 'config/external_datasets.csv'],
        'notes': ('In-house values are descriptive mean log2 CPM differences; external values are each dataset\'s reported log2 fold change '
                  'oriented so that negative means lower with Y absent (the XYY consensus is sign-flipped). Polarity is verified with a '
                  'male-specific Y gene where one is present in the table; otherwise the direction rests on the dataset\'s contrast label.')})
