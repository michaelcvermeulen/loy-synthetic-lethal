"""Plot-only final add-back displays, 8 October 2026. No model fitting.

Uses frozen library scores, confidence intervals, expression differences and fry
results. Main F02c contains the 150-gene directional tests. Supplement S02
contains the 150-gene heatmap (a), frozen-100 scores/contrasts (b,c) and the
39-gene RNA/ATAC overlap heatmap (d), formerly S04. Use --supp-only to update
the supplement without rewriting the main signature-recovery exports.
Use --main-only for layout-only updates to F02.
"""
from pathlib import Path
import sys
import hashlib
import shutil
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
DATA = FIG / 'source_data/addback_signature_150'
SUP = FIG / 'source_data/addback_signature_final'
UPSTREAM = ROOT.parent / 'RNAseq/analysis/HT1080_ADDBACK/THESIS_ADDBACK_2026/outputs/signature_150'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
SUP.mkdir(exist_ok=True)
inputs = []
for name in ['distance_to_control_by_gene.csv', 'signature_mapping.csv', 'fry_directional.csv']:
    source, target = UPSTREAM / name, SUP / name
    if not target.exists():
        shutil.copy2(source, target)
    assert source.read_bytes() == target.read_bytes(), f'Frozen source differs: {name}'
    inputs.append(dict(file=str(target.relative_to(ROOT)), source=str(source), sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
for name in ['panel_a_library_scores_sig150.csv', 'panel_a_group_means_sig150.csv', 'panel_b_percent_reversal_sig150.csv']:
    p = DATA / name
    inputs.append(dict(file=str(p.relative_to(ROOT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
FROZEN100 = FIG / 'source_data/ht1080_addback_2026'
for name in ['outputs__signature__frozen_signature_scores.csv', 'outputs__signature__score_contrasts.csv']:
    p = FROZEN100 / name
    inputs.append(dict(file=str(p.relative_to(ROOT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
overlap_source = FIG / 'source_data/addback_atac_overlap/gene_table.csv'
inputs.append(dict(file=str(overlap_source.relative_to(ROOT)), sha256=hashlib.sha256(overlap_source.read_bytes()).hexdigest()))
(SUP / 'input_manifest.json').write_text(json.dumps(inputs, indent=2), encoding='utf-8')

GROUPS = ['CTRL', 'LOY', 'LOY_LXHA', 'LOY_KDM5D_AB', 'LOY_UTY_AB', 'LOY_KDM5D_UTY_AB']
COLOUR = dict(zip(GROUPS, ['#3C5488', '#E64B35', '#8C7B6B', '#F0A38A', '#009E73', '#56B4E9']))
SHORT = {'CTRL':'sgCTRL', 'LOY':'Untransduced', 'LOY_LXHA':'Empty vector', 'LOY_KDM5D_AB':'KDM5D', 'LOY_UTY_AB':'UTY', 'LOY_KDM5D_UTY_AB':'UTY + KDM5D'}

def save(fig, stem):
    if '--supp-only' in sys.argv and stem.startswith('ch03_main_'):
        plt.close(fig)
        return
    for ext in ['pdf', 'svg', 'png']:
        target = FIG / f'{stem}.{ext}'
        stage = FIG / f'{stem}.staging.{ext}'
        fig.savefig(stage, dpi=400, facecolor='white')
        stage.replace(target)
    fig.savefig(FIG / 'previews' / f'{stem}_preview.png', dpi=160)
    plt.close(fig)
    print('written', stem)

def fmt_p(p):
    return f'{p:.3f}' if p >= .001 else f'{p:.1e}'.replace('e-0', 'e−')

def draw_directional_tests(ax):
    fry = pd.read_csv(SUP / 'fry_directional.csv')
    assert len(fry) == 12
    contrasts = ['UTY_v_LXHA','KDM5D_v_LXHA','DUAL_v_LXHA','DUAL_v_UTY','DUAL_v_KDM5D','LOY_v_LXHA']
    labels = ['UTY vs vector','KDM5D vs vector','UTY + KDM5D vs vector','UTY + KDM5D vs UTY','UTY + KDM5D vs KDM5D','Untransduced vs vector']
    for i, contrast in enumerate(contrasts):
        for j, arm in enumerate(['August_up','August_down']):
            row = fry[(fry.contrast==contrast)&(fry.signature_arm==arm)].iloc[0]
            increase = row.Direction == 'Up'
            ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,facecolor='#F7DEC5' if increase else '#D6E5F1',edgecolor='white',lw=1.5))
            ax.text(j,i,('Increase' if increase else 'Decrease')+f'\nFDR {fmt_p(row.FDR_12)}',ha='center',va='center',fontsize=8,weight='semibold' if row.FDR_12<.05 else 'normal')
    ax.set(xlim=(-.5,1.5),ylim=(5.5,-.5))
    ax.set_yticks(range(6), labels, fontsize=8)
    ax.set_xticks([0,1], ['Higher in LOY\nExpected: decrease','Lower in LOY\nExpected: increase'],fontsize=8)
    ax.xaxis.tick_top()
    ax.tick_params(length=0,pad=7)
    for spine in ax.spines.values(): spine.set_visible(False)
    fry.to_csv(SUP / 'main_panel_c_fry_display.csv', index=False)

scores = pd.read_csv(DATA / 'panel_a_library_scores_sig150.csv')
means = pd.read_csv(DATA / 'panel_a_group_means_sig150.csv').set_index('group').adjusted_mean
rev = pd.read_csv(DATA / 'panel_b_percent_reversal_sig150.csv')
assert scores.shape[0] == 18 and len(rev) == 4
fig = plt.figure(figsize=(9.4, 6.5))
a = fig.add_axes([.075, .605, .905, .35])
b = fig.add_axes([.13, .095, .285, .36])
a.axhspan(means.CTRL, means.LOY_LXHA, color='#eeeeee', zorder=0)
for group, label in [('CTRL','sgCTRL mean'), ('LOY_LXHA','Empty-vector mean')]:
    a.axhline(means[group], color='#777777', lw=.65, ls=(0,(4,3)), zorder=1)
    a.text(5.45, means[group]+.004, label, ha='right', va='bottom', fontsize=7, color='#555555')
for i, group in enumerate(GROUPS):
    a.plot([i-.28,i+.28], [means[group]]*2, color=COLOUR[group], lw=1.8)
    for _, row in scores[scores.group == group].iterrows():
        occasion = int(row.occasion)
        a.scatter(i+(occasion-2)*.12, row.score, marker=['o','^','s'][occasion-1], s=32, color=COLOUR[group], edgecolor='white', lw=.4, zorder=3)
a.set(xlim=(-.5,5.5), ylim=(-.08,.42), ylabel='LOY signature score')
a.set_xticks(range(6), ['sgCTRL','Untransduced','Empty\nvector','KDM5D','UTY','UTY +\nKDM5D'], fontsize=8)
a.tick_params(axis='x', length=0)
a.legend(handles=[Line2D([],[],linestyle='',marker=m,color='#555555',markersize=4,label=str(i)) for i,m in enumerate(['o','^','s'],1)], title='Collection occasion', ncol=3, loc='upper right', frameon=False, fontsize=7, title_fontsize=7, handletextpad=.3, columnspacing=1)
b.axvspan(0,100,color='#eeeeee',zorder=0)
b.axvline(0,color='#777777',lw=.65,ls=(0,(4,3)))
b.axvline(100,color='#777777',lw=.65,ls=(0,(4,3)))
for i, row in rev.iterrows():
    colour = COLOUR[row.comparison]
    b.plot([row.CI_low,row.CI_high],[i,i],color=colour,lw=1.5)
    b.scatter(row.percent_reversal,i,s=30,color=colour,zorder=3)
    b.text(row.percent_reversal,i-.17,f'{row.percent_reversal:.0f}%  (P = {fmt_p(row.P_difference)})',ha='center',va='bottom',fontsize=7)
b.set_yticks(range(4),[SHORT[g] for g in rev.comparison],fontsize=8)
b.set(xlim=(-165,115),ylim=(3.5,-.65),xlabel='Recovery of signature-score gap (%)')
b.tick_params(axis='y', length=0)
b.spines['left'].set_visible(False)
b.set_xticks([-150,-100,-50,0,50,100])
fig.text(.01,.98,'a',fontsize=12,weight='semibold')
fig.text(.01,.50,'b',fontsize=12,weight='semibold')
c = fig.add_axes([.72,.095,.26,.36])
draw_directional_tests(c)
fig.text(.48,.50,'c',fontsize=12,weight='semibold')
save(fig,'ch03_main_ht1080_addback_signature_recovery')
if '--main-only' in sys.argv:
    sys.exit(0)

mapping = pd.read_csv(SUP / 'signature_mapping.csv')
dist = pd.read_csv(SUP / 'distance_to_control_by_gene.csv')
fry = pd.read_csv(SUP / 'fry_directional.csv')
assert mapping.shape[0] == 150 and fry.shape[0] == 12
order = mapping.assign(abs_effect=mapping.AUG_difference.abs()).sort_values(['direction','abs_effect'],ascending=[False,False]).gene_id_stable.iloc[::-1].tolist()
heat_groups = ['LOY','LOY_LXHA','LOY_UTY_AB','LOY_KDM5D_AB','LOY_KDM5D_UTY_AB']
heat = dist.pivot(index='gene_id', columns='arm', values='arm_vs_CTRL')
heat['LOY_LXHA'] = dist.drop_duplicates('gene_id').set_index('gene_id').vector_vs_CTRL
heat = heat.loc[order,heat_groups]
assert heat.shape == (150,5) and not heat.isna().any().any()
heat.to_csv(SUP / 'panel_a_expression_display.csv')
fig = plt.figure(figsize=(9.6,9.2))
a = fig.add_axes([.135,.67,.31,.285])
cmap = LinearSegmentedColormap.from_list('expression',['#2C6DA4','white','#D9822B'])
im = a.imshow(heat.to_numpy(),aspect='auto',vmin=-4,vmax=4,cmap=cmap,interpolation='nearest',rasterized=True)
a.set_xticks(range(5),[SHORT[g].replace('Empty vector','Empty\nvector').replace('UTY + KDM5D','UTY +\nKDM5D') for g in heat_groups],rotation=40,ha='right',fontsize=8)
a.set_yticks([])
a.axhline(74.5,color='white',lw=1)
a.text(-.09,.75,'Lower in LOY\n(75 genes)',transform=a.transAxes,ha='right',va='center',fontsize=8)
a.text(-.09,.25,'Higher in LOY\n(75 genes)',transform=a.transAxes,ha='right',va='center',fontsize=8)
a.tick_params(length=0)
for spine in a.spines.values(): spine.set_visible(False)
cax = fig.add_axes([.175,.591,.23,.008])
cb = fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-4,0,4])
cb.set_label('Expression difference from sgCTRL (log2)',fontsize=8)
fig.text(.01,.975,'a',fontsize=12,weight='semibold')

# Added panels reproduce the saved frozen-100-gene score panels, not the
# 150-gene score data used in the main figure. No fitting or testing is done.
scores100 = pd.read_csv(FROZEN100 / 'outputs__signature__frozen_signature_scores.csv')
contrasts100 = pd.read_csv(FROZEN100 / 'outputs__signature__score_contrasts.csv').set_index('contrast')
assert len(scores100) == 18 and len(contrasts100) == 15
groups100 = ['CTRL','LOY','LOY_LXHA','LOY_UTY_AB','LOY_KDM5D_AB','LOY_KDM5D_UTY_AB']
colours100 = dict(zip(groups100,['#3C5488','#E64B35','#7E6148','#00A087','#F39B7F','#4DBBD5']))
c = fig.add_axes([.085,.328,.36,.19])
d = fig.add_axes([.245,.065,.20,.17])
for i, group in enumerate(groups100):
    rows = scores100[scores100.final_group == group]
    assert len(rows) == 3
    c.plot([i-.24,i+.24],[rows.signed_percentile_score.mean()]*2,color='#222222',lw=1.5,zorder=4)
    for _,row in rows.iterrows():
        occasion = int(row.collection_occasion)
        c.scatter(i+(occasion-2)*.11,row.signed_percentile_score,marker=['o','^','s'][occasion-1],s=28,color=colours100[group],zorder=3)
c.set(xlim=(-.5,5.5),ylim=(-.28,.235),ylabel='Frozen 100-gene signature score')
c.set_xticks(range(6),['sgCTRL','Untransduced','Empty vector','UTY','KDM5D','UTY + KDM5D'],rotation=35,ha='right',fontsize=8)
c.tick_params(axis='x',length=0)
c.legend(handles=[Line2D([],[],linestyle='',marker=m,color='#555555',markersize=4,label=str(i)) for i,m in enumerate(['o','^','s'],1)],title='Collection occasion',ncol=3,loc='upper right',frameon=False,fontsize=7,title_fontsize=7,handletextpad=.3,columnspacing=1)
forest100 = ['UTY_v_LXHA','KDM5D_v_LXHA','DUAL_v_LXHA','DUAL_v_UTY','DUAL_v_KDM5D']
labels100 = ['UTY vs vector','KDM5D vs vector','UTY + KDM5D vs vector','UTY + KDM5D vs UTY','UTY + KDM5D vs KDM5D']
for i,key in enumerate(forest100):
    row = contrasts100.loc[key]
    d.plot([row.CI_L,row.CI_R],[i,i],color='#222222',lw=1)
    d.plot([row.CI_L,row.CI_L],[i-.08,i+.08],color='#222222',lw=.8)
    d.plot([row.CI_R,row.CI_R],[i-.08,i+.08],color='#222222',lw=.8)
    d.scatter(row.difference,i,color='#222222',s=25,zorder=3)
    d.text(.19,i,fmt_p(row.PValue),ha='right',va='center',fontsize=8)
d.axvline(0,color='#888888',lw=.65,ls=(0,(4,3)))
d.set(xlim=(-.19,.195),ylim=(4.6,-.6),xlabel='Difference in frozen 100-gene score')
d.set_xticks([-.15,-.075,0,.075],['−0.15','−0.075','0','0.075'])
d.set_yticks(range(5),labels100,fontsize=8)
d.tick_params(axis='y',length=0)
d.spines['left'].set_visible(False)
d.text(.19,-.8,'P',ha='right',va='bottom',fontsize=8)
fig.text(.01,.545,'b',fontsize=12,weight='semibold')
fig.text(.01,.268,'c',fontsize=12,weight='semibold')
scores100.to_csv(SUP / 'panel_b_frozen100_scores_display.csv',index=False)
contrasts100.loc[forest100].to_csv(SUP / 'panel_c_frozen100_contrasts_display.csv')
# Former S04: preserve all 39 genes, original order, scale and FDR markers.
overlap = pd.read_csv(overlap_source)
overlap = overlap[(overlap['set'] == 'joint_FDR005') & overlap.tested_in_addback].sort_values('baseline_logFC')
assert len(overlap) == 39
keys = ['UTY_v_LXHA','KDM5D_v_LXHA','DUAL_v_LXHA','LOY_v_LXHA']
cols = ['baseline_logFC'] + [k + '_logFC' for k in keys]
overlap_values = overlap[cols].to_numpy()
lim = np.abs(overlap_values).max()
h = fig.add_axes([.61,.045,.285,.71])
im = h.imshow(overlap_values,aspect='auto',cmap=LinearSegmentedColormap.from_list('effects',['#2166ac','white','#d6604d']),norm=TwoSlopeNorm(0,-lim,lim),interpolation='none')
h.set_yticks(range(len(overlap)),overlap.gene_name,fontstyle='italic',fontsize=7)
h.set_xticks(range(5),['Y-depleted vs sgCTRL','UTY vs vector','KDM5D vs vector','UTY + KDM5D vs vector','Untransduced Y-depleted vs vector'],rotation=90,ha='center',va='bottom',fontsize=7)
h.xaxis.tick_top(); h.tick_params(length=0,pad=3)
for j,key in enumerate(keys,1):
    rows = np.flatnonzero(overlap[key + '_FDR'].to_numpy() < .05)
    h.scatter(np.full(len(rows),j),rows,s=6,c='black',linewidths=0)
for spine in h.spines.values(): spine.set_visible(False)
cax = fig.add_axes([.918,.34,.012,.19])
cb = fig.colorbar(im,cax=cax)
cb.set_label('log2 fold change',fontsize=8)
cb.outline.set_visible(False)
fig.text(.50,.975,'d',fontsize=12,weight='semibold')
overlap.to_csv(SUP / 'panel_d_atac_overlap_display.csv',index=False)
save(fig,'ch03_supp_ht1080_addback_baseline_reversal')
