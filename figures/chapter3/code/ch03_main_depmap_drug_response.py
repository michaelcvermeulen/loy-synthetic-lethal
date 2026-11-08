from pathlib import Path
import json, hashlib, shutil, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
root=Path(__file__).resolve().parents[3]
out=root/'figures/chapter3'
export=root/'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-05'
source=export/'derived/loy_reanalysis_v1/drug_pan_results.csv'
qqsource=export/'results/loy_reanalysis_v1/figure_source_data/fig4a_drug_qq.csv'
d=pd.read_csv(source)
qq=pd.read_csv(qqsource)
data=out/'source_data/depmap_drug_response_main';data.mkdir(parents=True,exist_ok=True)
shutil.copy2(source,data/'drug_pan_results.csv');shutil.copy2(qqsource,data/'fig4a_drug_qq.csv')
assert len(d)==6413 and (d.fdr<.05).sum()==0 and len(qq)==6413
for font in ['Regular','SemiBold']:font_manager.fontManager.addfont(str(root/f'skills/thesis-figure-style/assets/fonts/OpenSans-{font}.ttf'))
plt.rcParams.update({'font.family':'Open Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none'})
sys.path.insert(0,str(root/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
plt.rcParams.update({'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(9.6,10.2))
gs=fig.add_gridspec(1,7,wspace=.22,left=.063,right=.984,bottom=.871,top=.974)
order=['REP.PRIMARY','REP.1M','REP.300','PRISM secondary','GDSC1','GDSC2','CTRP v2']
for i,name in enumerate(order):
    ax=fig.add_subplot(gs[0,i]);q=qq[qq.dataset_label==name]
    ax.scatter(q.expected,q.observed,s=4,c='#555555',alpha=.7,edgecolors='none')
    ax.plot([0,4.6],[0,4.6],lw=.65,ls=(0,(2,3)),c='#666666')
    ax.set(xlim=(-.1,4.6),ylim=(-.1,4.6),xticks=[0,2,4],yticks=[0,2,4])
    ax.tick_params(labelsize=7,length=2)
    ax.text(.5,1.06,name,transform=ax.transAxes,ha='center',fontsize=6.6)
    if i==0:ax.set_ylabel('Observed −log₁₀ P',fontsize=8)
    else:ax.set_yticklabels([])
    if i==3:ax.set_xlabel('Expected −log₁₀ P',fontsize=8)
fig.text(.012,.996,'a',weight='semibold',fontsize=11,va='top')
column_x=[.063,.389,.715]
volcano_axes=[fig.add_axes([x,.629,.269,.178]) for x in column_x]
for i,(key,name) in enumerate([('prism_REP1M','PRISM REP.1M'),('prism_secondary','PRISM secondary'),('gdsc2','GDSC2')]):
    z=d[d.dataset==key].copy();z['neg_log10_p']=-np.log10(z.p);z['bh_first_rank_p']=.05/len(z)
    z.to_csv(data/f'{key}_volcano.csv',index=False)
    ax=volcano_axes[i]
    ax.scatter(z.standardized,z.neg_log10_p,s=11,c='#929292',alpha=.65,edgecolors='none')
    ax.axhline(-np.log10(.05/len(z)),lw=.8,ls=(0,(4,3)),c='#555555')
    ax.axvline(0,lw=.6,c='#CCCCCC',zorder=0)
    ax.text(-.76,-np.log10(.05/len(z))+.08,'BH first-rank reference',fontsize=6.2,va='bottom')
    ax.set(xlim=(-.8,.8),ylim=(-.06,5.15),xticks=[-.5,0,.5],yticks=[0,1,2,3,4,5])
    ax.text(.5,1.09,name,ha='center',transform=ax.transAxes,fontsize=8)
    ax.text(.5,1.015,f'{len(z):,} compounds; minimum FDR = {z.fdr.min():.3f}',ha='center',transform=ax.transAxes,fontsize=6.5)
    ax.text(-.18,1.12,chr(98+i),transform=ax.transAxes,weight='semibold',fontsize=11)
    ax.set_ylabel('−log₁₀ nominal P',fontsize=7)
    ax.set_xlabel('LOY-associated response difference (SD)',fontsize=6.8)
    if key=='prism_secondary':
        selected=z.nsmallest(2,'p')
        positions=[(-.68,3.55),(-.68,2.45)]
    else:
        selected=pd.concat([z[z.standardized<0].nsmallest(1,'p'),z[z.standardized>0].nsmallest(1,'p')])
        positions=[(-.73,2.8),(.06,2.75)] if key=='gdsc2' else [(-.73,3.28),(.06,3.25)]
    for (_,r),(x,y) in zip(selected.iterrows(),positions):
        label=r.compound_name.strip().capitalize().replace('Mycophenolate-mofetil','Mycophenolate\nmofetil')
        annotations = {
            'ENTASOBULIN': 'Tubulin polymerization\ninhibitor',
            'VIXOTRIGINE': 'Sodium channel\nblocker',
            'MYCOPHENOLATE-MOFETIL': 'IMPDH1 / IMPDH2',
            'RIFAMPIN': 'NR1I2 / SLCO1A2\nSLCO1B1 / SLCO1B3',
            'ACETALAX': 'Target not annotated',
            'PYRIDOSTATIN': 'Target not annotated',
        }
        if key != 'gdsc2':
            label += '\n' + annotations[r.compound_name.strip()]
        ax.scatter([r.standardized],[-np.log10(r.p)],s=17,c='#343434',zorder=4)
        ax.annotate(label,(r.standardized,-np.log10(r.p)),xytext=(x,y),fontsize=5.8,linespacing=1.15,ha='left',arrowprops={'arrowstyle':'-','lw':.6,'color':'#555555'},zorder=5)
    ax.tick_params(length=3,labelsize=7)
# Panel f (formerly e) uses the saved SR-1664 within-subset comparison; no model is refitted.
box_data=out/'source_data/sr1664_bcell_response'
points=pd.read_csv(box_data/'plotted_points.csv')
hit=pd.read_csv(box_data/'saved_test.csv').iloc[0]
assert len(points)==24 and not points.ModelID.duplicated().any()
assert hit.compound_name=='SR-1664' and hit.dataset=='prism_REP1M'
shutil.copy2(box_data/'plotted_points.csv',data/'panel_f_sr1664_points.csv')
shutil.copy2(box_data/'saved_test.csv',data/'panel_f_sr1664_test.csv')
ax=fig.add_axes([.715,.330,.269,.235])
rng=np.random.default_rng(20261008)
for x,(group,colour,n) in enumerate([('Y-positive','#3C5488',16),('LOY','#E64B35',8)]):
    vals=points.loc[points.y_status.eq(group),'response'].to_numpy()
    assert len(vals)==n
    ax.boxplot([vals],positions=[x],widths=.46,manage_ticks=False,showfliers=False,patch_artist=True,
        boxprops={'facecolor':colour,'edgecolor':colour,'alpha':.20,'lw':.8},
        medianprops={'color':colour,'lw':1.5},whiskerprops={'color':colour,'lw':.8},
        capprops={'color':colour,'lw':.8})
    ax.scatter(x+rng.uniform(-.12,.12,len(vals)),vals,s=24,c=colour,edgecolors='white',linewidths=.4,alpha=.9,zorder=3)
low,high=min(points.response.min(),0),max(points.response.max(),0)
spread=high-low
ax.set_ylim(low-.08*spread,high+.53*spread)
ax.set_xlim(-.55,1.55)
ax.axhline(0,color='#999999',linestyle=(0,(3,3)),linewidth=.7,zorder=0)
ax.set_xticks([0,1],['Y-positive\nn = 16','LOY\nn = 8'])
ax.set_ylabel('SR-1664 response (log₂ fold change)',fontsize=7)
ax.tick_params(labelsize=7)
ax.text(-.18,1.12,'f',transform=ax.transAxes,weight='semibold',fontsize=11)
ax.text(.5,1.09,'SR-1664 · PRISM REP.1M',transform=ax.transAxes,ha='center',fontsize=8)
ax.text(.5,1.015,'Mature B-cell neoplasms',transform=ax.transAxes,ha='center',fontsize=7)
trans=ax.get_xaxis_transform()
ax.plot([0,0,1,1],[.76,.78,.78,.76],transform=trans,color='#222222',linewidth=.8)
ax.text(.5,.985,f'Δ = {hit.inside:+.3f}',transform=trans,ha='center',va='top',fontsize=8)
exponent=int(np.floor(np.log10(hit.inside_p)))
ax.text(.5,.91,fr'P = {hit.inside_p/10**exponent:.2f} × 10$^{{{exponent}}}$',transform=trans,ha='center',va='top',fontsize=8)
ax.text(.5,.835,f'FDR = {hit.inside_fdr_class:.3f}',transform=trans,ha='center',va='top',fontsize=8)


from matplotlib.colors import Normalize
ROOT=root
FIG=out
EXPORT = ROOT / 'research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-07'
SOURCE = EXPORT / 'derived/loy_reanalysis_v1/final'
DATA = data / 'controls_pooling'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(8)
DATA.mkdir(parents=True, exist_ok=True)
input_names = ['pipeline_controls.csv', 'drug_feature_annotation.csv', 'drug_pan_consensus.csv',
               'drug_pan_mechanism.csv', 'forest_y_status_summary.csv', 'drug_pan_results.csv']
tables, manifest = {}, {}
for name in input_names:
    path = SOURCE / name
    tables[name] = pd.read_csv(path)
    shutil.copy2(path, DATA / name)
    manifest[name] = {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}

controls = tables['pipeline_controls.csv'].query("dataset != 'crispr_effect'").copy()
controls = controls.merge(tables['drug_feature_annotation.csv'][['dataset', 'feature', 'compound_name']],
                          on=['dataset', 'feature'], how='left', validate='one_to_one')
assert len(controls) == 12 and controls.compound_name.notna().all()
assert controls['pass'].all() and controls.p.lt(.01).all()
assert np.array_equal(np.sign(controls.estimate), controls.expected_sign)
pooled = tables['drug_pan_consensus.csv']
classes = tables['drug_pan_mechanism.csv']
rf = tables['forest_y_status_summary.csv']
assert len(pooled) == 352 and pooled.fdr.ge(.05).all()
assert pooled.families.ge(2).all()
assert classes.groupby('class_type').size().to_dict() == {'mechanism': 218, 'target': 385}
assert classes.fdr.ge(.05).all()
assert rf.query("dataset != 'crispr_effect'").stable.eq(0).all()
controls.to_csv(DATA / 'plotted_controls.csv', index=False)

order = ['prism_REPPRIMARY', 'prism_REP1M', 'prism_REP300', 'prism_secondary', 'gdsc1', 'gdsc2', 'ctrp']
labels = ['PRISM REP.PRIMARY', 'PRISM REP.1M', 'PRISM REP.300', 'PRISM secondary', 'GDSC1', 'GDSC2', 'CTRP v2']
types = ['BRAF altered -> BRAF inhibitor sensitivity', 'TP53 altered -> MDM2 inhibitor resistance']
matrix = np.full((7, 2), np.nan)
for row in controls.itertuples():
    matrix[order.index(row.dataset), types.index(row.control)] = -np.log10(row.p)

a = fig.add_axes([.157, .330, .443, .235])
colourmap = plt.get_cmap('magma_r').copy()
colourmap.set_bad('#eeeeee')
im = a.imshow(matrix, cmap=colourmap, norm=Normalize(0, 32), aspect='auto')
a.set_yticks(np.arange(7), labels, fontsize=7)
a.set_xticks([0, 1], ['BRAF alteration\nBRAF inhibitors', 'TP53 alteration\nMDM2 inhibitors'], fontsize=7)
a.xaxis.tick_top()
a.tick_params(length=0, pad=5)
for spine in a.spines.values():
    spine.set_visible(False)
a.set_xticks(np.arange(-.5, 2, 1), minor=True)
a.set_yticks(np.arange(-.5, 7, 1), minor=True)
a.grid(which='minor', color='white', linewidth=1.2)
a.tick_params(which='minor', length=0)
for i, dataset in enumerate(order):
    for j, control in enumerate(types):
        rows = controls[(controls.dataset == dataset) & (controls.control == control)]
        if rows.empty:
            a.text(j, i, 'Not tested', ha='center', va='center', fontsize=7, color='#777777')
            continue
        row = rows.iloc[0]
        colour = 'white' if matrix[i, j] > 17 else '#222222'
        name = row.compound_name.title().replace('Plx', 'PLX').replace('Nutlin-3A', 'Nutlin-3a')
        a.text(j, i, f'{name}\nP = {row.p:.1e}\nn = {row.n_exposed}/{row.n_unexposed}',
               ha='center', va='center', fontsize=6.5, linespacing=1.05, color=colour)
a.text(.25, -.058, 'Altered: more sensitive', transform=a.transAxes, ha='center', va='top', fontsize=7)
a.text(.75, -.058, 'Altered: less sensitive', transform=a.transAxes, ha='center', va='top', fontsize=7)
cax = fig.add_axes([.613, .365, .010, .165])
bar = fig.colorbar(im, cax=cax, ticks=[0, 10, 20, 30])
bar.set_label('−log₁₀ nominal P', fontsize=7, labelpad=4)
bar.ax.tick_params(labelsize=6.5, length=2)
bar.outline.set_visible(False)

b,c,d = [fig.add_axes([x,.063,.269,.187]) for x in column_x]

def association_panel(ax, frame, x, xlim, ticks, xlabel, family_label, colour):
    ax.scatter(frame[x], -np.log10(frame.p), s=13, color=colour, alpha=.65, linewidths=0)
    ax.axvline(0, color='#c9c9c9', lw=.6, zorder=0)
    nominal = -np.log10(.05)
    ax.axhline(nominal, color='#858585', lw=.7, ls=(0, (3, 3)), zorder=0)
    ax.text(.97, nominal + .08, 'Nominal P = 0.05', transform=ax.get_yaxis_transform(),
            fontsize=6.5, ha='right', color='#666666')
    ax.set(xlim=xlim, ylim=(-.05, 3.35), xticks=ticks, yticks=[0, 1, 2, 3], xlabel=xlabel, ylabel='−log₁₀ nominal P')
    ax.tick_params(length=3, pad=2, labelsize=7.5)
    ax.text(.5, 1.10, family_label, transform=ax.transAxes, ha='center', fontsize=7)
    ax.text(.5, 1.025, f'0/{len(frame):,} at FDR <0.05; minimum FDR = {frame.fdr.min():.3f}',
            transform=ax.transAxes, ha='center', fontsize=7)

association_panel(b, pooled, 'estimate', (-.43, .43), [-.4, -.2, 0, .2, .4],
                  'Pooled response difference (SD)', 'Shared compounds', '#555555')
association_panel(c, classes[classes.class_type.eq('mechanism')], 'mean_t', (-1.1, 1.75), [-1, 0, 1],
                  'Mean signed response t statistic', 'Mechanism-of-action classes', '#517491')
association_panel(d, classes[classes.class_type.eq('target')], 'mean_t', (-1.1, 1.75), [-1, 0, 1],
                  'Mean signed response t statistic', 'Target classes', '#517491')

# Label the four lowest nominal P values per panel, retaining cutoff ties.
# Coincident mechanism classes share a callout so their labels do not overlap.
label_layout = {
    'b': {
        'TOSEDOSTAT': ('Tosedostat', -.40, 3.16, 'left'),
        'ABT-737': ('ABT-737', .27, 2.86, 'center'),
        'GSK429286A': ('GSK429286A', .035, 2.04, 'center'),
        'LENALIDOMIDE': ('Lenalidomide', -.405, 2.26, 'left'),
        'BERZOSERTIB': ('Berzosertib', -.39, 1.08, 'left'),
        'ETP-46464': ('ETP-46464', -.40, .79, 'left'),
        'BEZ235': ('BEZ235', -.37, .32, 'left'),
        'VE 821': ('VE-821', .29, .54, 'center'),
        'CERALASERTIB': ('Ceralasertib', .29, .17, 'center'),
    },
    'c': {
        'TACHYKININ ANTAGONIST': ('Tachykinin\nantagonist', .35, 2.78, 'center'),
        'SODIUM CHANNEL BLOCKER': ('Sodium channel\nblocker', -.38, 2.12, 'center'),
        'FGFR2': ('FGFR2; inhibitor\nof FGFR1', 1.23, 2.25, 'center'),
        'INHIBITOR OF FGFR1': None,
    },
    'd': {
        'SLC12A1': ('SLC12A1', .98, 3.10, 'center'),
        'SLC12A2': ('SLC12A2', 1.49, 3.10, 'center'),
        'TRPV1': ('TRPV1', .37, 2.78, 'center'),
        'GPR35': ('GPR35', .54, 2.05, 'center'),
        'ATR': ('ATR', -.91, 1.13, 'center'),
    },
}
label_records = []
for panel, ax, frame, xcolumn, namecolumn in [
    ('b', b, pooled, 'estimate', 'compound_name'),
    ('c', c, classes[classes.class_type.eq('mechanism')], 'mean_t', 'class'),
    ('d', d, classes[classes.class_type.eq('target')], 'mean_t', 'class'),
]:
    ordered = frame.sort_values(['p', namecolumn])
    selected_labels = ordered[ordered.p.le(ordered.iloc[3].p)].copy()
    selected_labels['selection_reason'] = 'top nominal P'
    if panel == 'b':
        atr = ordered[ordered.targets.fillna('').str.contains(r'(?:^|[;,])\s*ATR\s*(?:$|[;,])')].copy()
        assert set(atr.compound_name) == {'BERZOSERTIB', 'ETP-46464', 'BEZ235', 'VE 821', 'CERALASERTIB'}
        assert not selected_labels[namecolumn].isin(atr[namecolumn]).any()
        atr['selection_reason'] = 'ATR target annotation; author-requested label'
        selected_labels = pd.concat([selected_labels, atr], ignore_index=True)
    if panel == 'd':
        atr = ordered[ordered[namecolumn].eq('ATR')].copy()
        assert len(atr) == 1 and not selected_labels[namecolumn].eq('ATR').any()
        atr['selection_reason'] = 'literature-motivated candidate; author-requested label'
        selected_labels = pd.concat([selected_labels, atr], ignore_index=True)
    assert set(selected_labels[namecolumn]) == set(label_layout[panel])
    if panel == 'c':
        overlapping = selected_labels[selected_labels['class'].isin(['FGFR2', 'INHIBITOR OF FGFR1'])]
        assert overlapping[['mean_t', 'p']].drop_duplicates().shape[0] == 1
    for _, row in selected_labels.iterrows():
        source_label = row[namecolumn]
        label_records.append({'panel': {'b':'g','c':'h','d':'i'}[panel], 'source_label': source_label,
                              'x': row[xcolumn], 'p': row.p, 'fdr': row.fdr,
                              'selection_reason': row['selection_reason']})
        layout = label_layout[panel][source_label]
        if layout is None:
            continue
        text, tx, ty, ha = layout
        point = (row[xcolumn], -np.log10(row.p))
        ax.scatter(*point, s=20, color='#333333', linewidths=0, zorder=4)
        ax.annotate(text, point, xytext=(tx, ty), ha=ha, va='center', fontsize=6.3,
                    linespacing=1.12, color='#292929',
                    bbox=dict(facecolor='white', edgecolor='none', pad=.6),
                    arrowprops=dict(arrowstyle='-', lw=.6, color='#666666', shrinkA=2, shrinkB=3),
                    zorder=5)
pd.DataFrame(label_records).to_csv(DATA / 'labelled_associations.csv', index=False)

# ATR is absent from the tested mechanism-class table: preserve that absence.
# Reproduce only its eligibility counts, without combining classes or fitting tests.
drug_results = tables['drug_pan_results.csv']
atr_mechanism_rows = []
for _, row in drug_results.iterrows():
    for term in str(row.mechanism).split(','):
        term = term.strip()
        if 'ATR' in term.split():
            atr_mechanism_rows.append({'class': term,
                'compound': row.compound_id if pd.notna(row.compound_id) else row.compound_name})
atr_mechanism_counts = pd.DataFrame(atr_mechanism_rows).groupby('class').agg(
    compounds=('compound', 'nunique'), tests=('compound', 'size')).reset_index()
assert atr_mechanism_counts.set_index('class').compounds.to_dict() == {
    'ATR KINASE INHIBITOR': 3, 'INHIBITOR OF ATR': 2, 'INHIBITOR OF ATR AND MTOR': 1}
assert not classes[classes.class_type.eq('mechanism')]['class'].isin(atr_mechanism_counts['class']).any()
atr_mechanism_counts.to_csv(DATA / 'atr_mechanism_eligibility.csv', index=False)
c.text(-.98, 3.17, 'ATR groups: 1–3 drugs;\nnot tested (minimum 5)',
       ha='left', va='top', fontsize=6.2, linespacing=1.2, color='#555555')


for x,letter in zip(column_x,['g','h','i']):
    fig.text(x-.048,.279,letter,fontsize=11,weight='semibold')
fig.text(.015,.593,'e',fontsize=11,weight='semibold')
provenance = {
    'inputs': manifest, 'font': font, 'statistics_refitted': False, 'figure_inches': [9.6, 10.2], 'current_panel_mapping': {'a':'e','b':'g','c':'h','d':'i'},
    'control_selection': 'Saved minimum-P eligible matching inhibitor per genotype and dataset; nominal P, not FDR.',
    'control_n_order': 'genotype-altered / unaltered; no control in REP.300',
    'pooled_family': '352 compounds present in at least two screening families; PRISM variants form one family.',
    'pooled_test': 'Saved mixed-model likelihood-ratio P; standardised response, lineage and source fixed effects, ModelID random intercept.',
    'class_test': 'Saved two-sided within-lineage permutation P, 2,000 permutations; BH separately for mechanisms and targets; >=5 compounds/class.',
    'threshold_line': 'Nominal P=0.05 only. FDR results are reported separately above each panel.',
    'labels': 'Four lowest nominal P values per panel b-d, including ties at the fourth value, plus five ATR-annotated pooled compounds in b and ATR in d at the author request; coincident FGFR2 and INHIBITOR OF FGFR1 classes share a label. Highlighting denotes labelled points, not FDR significance.',
    'atr_mechanism_note': 'No ATR mechanism class was tested: the three distinct source annotations have 3, 2 and 1 compounds, below the existing minimum of 5. Panel c displays an explanatory note, not a new data point or pooled class test.',
    'counts': {'drug_controls': 12, 'pooled_compounds': 352, 'mechanism_classes': 218, 'target_classes': 385},
    'forest_reference': 'Combined CH3-S3.4.1-S01 panel d; drug stable counts all zero in current export.'
}
(DATA / 'provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')

d=pd.read_csv(source)
stem='ch03_main_depmap_drug_response'
for ext in ['pdf','svg','png']:
    dest=out/f'{stem}.{ext}';fig.savefig(dest,dpi=300,facecolor='white')
plt.close(fig)
(data/'provenance.json').write_text(json.dumps({'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'qq_source':str(qqsource),'qq_sha256':hashlib.sha256(qqsource.read_bytes()).hexdigest(),'n_tests':len(d),'fdr_hits':int((d.fdr<.05).sum()),'volcano_y':'-log10 nominal P','threshold':'BH first-rank reference P=0.05/m, not an observed BH rejection cutoff; no rejections in any dataset','volcano_x':'saved standardized lineage-adjusted effect','labels':'Most nominally associated compound per direction for REP.1M/GDSC2; two smallest nominal P values for secondary. Labels are not significant hits; 5.8 pt compound plus saved targets or mechanism, GDSC2 labels show compound names only.','font':'bundled Open Sans','code':str(Path(__file__).resolve())},indent=2))
provenance=json.loads((data/'provenance.json').read_text())
provenance.update({'panel_f':'SR-1664 in mature B-cell neoplasms, PRISM REP.1M; 16 Y-positive and 8 LOY male models',
    'panel_f_source':str(box_data),'panel_f_effect':float(hit.inside),'panel_f_p':float(hit.inside_p),
    'panel_f_fdr':float(hit.inside_fdr_class),'panel_f_inference':'Saved within-subset OLS on rank-normalised responses; original BH FDR across the disease class. No interaction displayed.',
    'models_refitted':False,'figure_inches':[9.6,10.2],
    'merged_from':'CH3-S3.3.2-S06; controls e, pooled compounds g, mechanisms h, targets i',
    'controls_pooling_provenance':'controls_pooling/provenance.json',
    'panel_mapping':{'F02 a-d':'a-d','F02 e':'f','S06 a':'e','S06 b-d':'g-i'},'panel_f_points_sha256':hashlib.sha256((box_data/'plotted_points.csv').read_bytes()).hexdigest(),
    'panel_f_test_sha256':hashlib.sha256((box_data/'saved_test.csv').read_bytes()).hexdigest(),
    'exports':{ext:hashlib.sha256((out/f'{stem}.{ext}').read_bytes()).hexdigest() for ext in ['pdf','svg','png']}})
(data/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
print(out/f'{stem}.png')
