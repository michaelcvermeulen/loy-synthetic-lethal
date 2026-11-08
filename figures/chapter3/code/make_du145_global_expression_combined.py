"""Approved DU145 F01/F03 merge and S02 extension; frozen inputs, no model refits.

Main a = former F01a; main b-f = former F03a-e.
S02 a-e retained; f = former F01b top; g = former F01b bottom; h = former F01c.
"""
from pathlib import Path
import argparse, hashlib, json, shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import squareform
from adjustText import adjust_text
from du145_rna_common import COND, COND_LABEL, PASSAGE_MARKER, DARK, GREY, BLUE, RED, MM, short, panel
from du145_reporting import annotate, LABEL
FIG=Path(__file__).resolve().parents[1]
SRC=FIG/'source_data/du145_global_expression_combined'
T=SRC/'inputs'
def read(name):
    return pd.read_csv(T/(name+'.csv'))

def export_current(fig, stem, folder, tables):
    out=SRC/('main' if '_main_' in stem else 'supp')
    out.mkdir(exist_ok=True)
    for name, table in tables.items(): table.to_csv(out/(name+'.csv'),index=False)
    selected=FIG/'selected'/folder
    selected.mkdir(exist_ok=True)
    outputs=[]
    for ext in ['pdf','svg','png']:
        path=FIG/(stem+'.'+ext)
        fig.savefig(path,dpi=600)
        shutil.copy2(path,selected/(folder+'.'+ext))
        outputs.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    fig.savefig(FIG/'previews'/(stem+'_preview.png'),dpi=170)
    plt.close(fig)
    (out/'provenance.json').write_text(json.dumps({'figure':stem,'font':'Open Sans','input_manifest':'../inputs.json',
        'scope':('Added TAF9 and PXDN logCPM displays from saved inputs and compacted layout; no model refits or threshold changes.' if '_main_' in stem else 'Layout only; retained numerical values, thresholds, statistical families and selected genes.'),
        'exports':outputs},indent=2),encoding='utf-8')
    print(stem,flush=True)

def draw_main():
    de = annotate(read('DE_paired_voom_primary'))
    eff = read('DU145_LOY_vs_sgCTRL_primary_effects')
    cn = read('CN_gene_level')
    arm = read('CN_arm_level_RNA_vs_DNA')
    conc = read('CN_RNA_DNA_concordance').set_index('level')

    de = de.merge(cn[['gene_id', 'CN_class', 't2t_arm']], on='gene_id', how='left')
    pair_cols = ['pair_P12', 'pair_P13', 'pair_P14']
    eff = eff.set_index('gene_id')

    fig = plt.figure(figsize=(7.2, 8.1))
    p = fig.add_axes([.10, .755, .34, .21])
    a = fig.add_axes([.60, .755, .37, .21])
    b = fig.add_axes([.10, .425, .34, .225])
    c = fig.add_axes([.60, .425, .37, .225])
    d = fig.add_axes([.10, .08, .245, .22])
    e = fig.add_axes([.445, .08, .12, .22])
    pca, var, assoc, qc, cor = (read(n) for n in ['PCA_core', 'PCA_variance_core', 'PCA_design_associations', 'library_QC', 'sample_correlations_core'])
    assert len(pca) == 6 and len(qc) == 6 and len(cor) == 36
    v = var.set_index('PC').core_variance
    a1 = assoc.set_index('axis')

    # a: PCA with passage-matched libraries joined.
    for passage, grp in pca.groupby('passage_block'):
        grp = grp.set_index('treatment')
        p.plot([grp.loc['CTRL', 'PC1'], grp.loc['LOY', 'PC1']], [grp.loc['CTRL', 'PC2'], grp.loc['LOY', 'PC2']],
               color='#AAAAAA', lw=.7, ls='--', zorder=1)
    for r in pca.itertuples():
        p.scatter(r.PC1, r.PC2, s=34, c=COND[r.treatment], marker=PASSAGE_MARKER[r.passage_block],
                  edgecolors='white', linewidths=.4, zorder=3)
    p.set_xlabel(f'PC1 ({100 * v[1]:.1f}%)'); p.set_ylabel(f'PC2 ({100 * v[2]:.1f}%)')
    p.margins(.18)
    p.text(.03, .97, f"PC1 ~ passage: R² = {a1.loc['PC1', 'passage_R2']:.2f}\nPC2 ~ condition: R² = {a1.loc['PC2', 'treatment_R2']:.2f}",
           transform=p.transAxes, va='top', fontsize=6.8, color=DARK)
    handles = [Line2D([], [], ls='', marker='o', color=COND[k], label=COND_LABEL[k], ms=3.6) for k in ['CTRL', 'LOY']]
    handles += [Line2D([], [], ls='', marker=m, color='#777777', label=p, ms=3.2) for p, m in PASSAGE_MARKER.items()]
    # Off the axes, in the free strip below the x-axis label.
    p.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.0, -0.34), bbox_transform=p.transAxes,
             fontsize=5.8, ncol=5, handlelength=.7, handletextpad=.15, columnspacing=.8, borderpad=.1, borderaxespad=0)
    panel(fig, p, 'a', '', dx=-.062, dy=.012)

    # ---- a: volcano; colour by thesis call and copy-number class.
    sig = de.FDR < .05
    thesis = de.reporting_subset
    cnc = de.CN_class.eq('CN_changed_concordant')
    cats = ['CN_stable', 'CN_changed_concordant', 'CN_changed_discordant']
    cnde = pd.DataFrame({key: de.loc[mask,'CN_class'].value_counts().reindex(cats,fill_value=0) for key,mask in [('FDR<0.05',sig),('reporting_subset',thesis)]}).T
    cnde.index.name='set'

    a.scatter(de.logFC[~sig], -np.log10(de.PValue[~sig]), s=3, color='#DDDDDD', lw=0, rasterized=True)
    a.scatter(de.logFC[sig & ~thesis & ~cnc], -np.log10(de.PValue[sig & ~thesis & ~cnc]), s=9, color='#F0B27A', lw=0, label='FDR < 0.05')
    a.scatter(de.logFC[sig & cnc & ~thesis], -np.log10(de.PValue[sig & cnc & ~thesis]), s=9, color=BLUE, lw=0, alpha=.75, label='FDR < 0.05, CN-concordant')
    a.scatter(de.logFC[thesis], -np.log10(de.PValue[thesis]), s=8, color=RED, lw=0, label=LABEL)
    for x in (-.5, .5):
        a.axvline(x, ls='--', lw=.6, color='#AAAAAA')
    lab = de[thesis & ~de.gene_name.eq('CCR9')].nsmallest(14, 'FDR')

    a.set_xlim(-5.4, 3.4); a.set_ylim(0, -np.log10(de.PValue.min()) * 1.95)
    volcano_labels = {
        'TACSTD2':(-2.8,10.7),'AP1M2':(-2.8,9.7),'THBD':(-2.8,8.7),'CLDN7':(-2.8,7.7),
        'SCEL':(-2.8,5.7),'EBF4':(-2.8,4.7),'GRHL2':(-4.95,6.8),
        'TAF9':(1.65,10.2),'F13A1':(1.65,9.2),'COL4A2':(1.65,8.2),'CREB3L1':(1.65,7.2),
        'SLC16A2':(1.65,6.2),'ST14':(1.65,5.2),'SPOCK2':(1.65,4.2)}
    assert set(volcano_labels) == set(lab.gene_name)
    for row in lab.itertuples():
        a.annotate(row.gene_name, (row.logFC,-np.log10(row.PValue)), volcano_labels[row.gene_name],
                   fontsize=6.3,color=DARK,ha='left',va='center',
                   bbox=dict(boxstyle='square,pad=.05',fc='white',ec='none',alpha=.95),
                   arrowprops=dict(arrowstyle='-',color='#AAAAAA',lw=.45,shrinkA=2,shrinkB=2))
    a.set_xlabel('log$_2$FC (Y-depleted / sgCTRL)'); a.set_ylabel('−log$_{10}$ P')
    a.legend(fontsize=5.7, loc='upper left', borderaxespad=.3, handletextpad=.3)
    n = cnde.loc['reporting_subset']
    a.text(.985, .97, f"147 at FDR < 0.05\n{int(thesis.sum())} meet additional thresholds\n({int(n.CN_stable)} copy-number stable)",
           transform=a.transAxes, ha='right', va='top', fontsize=6.0, color=DARK)
    panel(fig, a, 'b', '', dx=-.062, dy=.012)

    # ---- b: agreement of the three passage-pair differences.
    pd_all = eff.loc[de.gene_id, pair_cols].to_numpy()
    rng = np.random.default_rng(42)
    b.scatter(pd_all[~sig.to_numpy()].min(1), pd_all[~sig.to_numpy()].max(1), s=3, color='#DDDDDD', lw=0, rasterized=True)
    b.scatter(pd_all[(sig & ~thesis).to_numpy()].min(1), pd_all[(sig & ~thesis).to_numpy()].max(1), s=9, color='#F0B27A', lw=0)
    b.scatter(pd_all[thesis.to_numpy()].min(1), pd_all[thesis.to_numpy()].max(1), s=8, color=RED, lw=0)
    lim = (-5.6, 2.8)
    b.axhline(0, lw=.6, color='#CCCCCC'); b.axvline(0, lw=.6, color='#CCCCCC')
    b.fill_between([0, lim[1]], [0, 0], [lim[1], lim[1]], color=RED, alpha=.06, lw=0)
    b.fill_between([lim[0], 0], [lim[0], lim[0]], [0, 0], color=RED, alpha=.06, lw=0)
    b.set_xlim(lim); b.set_ylim(lim)
    b.set_xlabel('Lowest change across passages\n(log$_2$ CPM, Y-depleted − sgCTRL)', fontsize=7.5)
    b.set_ylabel('Highest change across passages\n(log$_2$ CPM, Y-depleted − sgCTRL)', fontsize=7.5)
    b.text(-2.8, -4.75, 'Lower in all\n3 passages', fontsize=7, color=DARK, ha='center', va='center')
    b.text(1.65, .28, 'Higher in all\n3 passages', fontsize=6.3, color=DARK, ha='center')
    b.text(-2.75, 2.15, 'Higher in some passages,\nlower in others', fontsize=6.3, color=DARK, ha='center', va='center')
    b.text(1.42, -3.4, 'Each dot: one gene\nPassages: P12–P14', fontsize=6.3,
           color=DARK, ha='center', va='center')
    # Direct labels identify isolated points and the strongest consistent changes.
    # Positions are fixed for reproducible spacing; leaders end at the measured point.
    label_positions = {
        'LINC02967': (-5.15, -1.35),
        'GRHL2': (-4.95, -3.95),
        'CCR9': (-3.45, -3.35),
        'SFTPA1': (-4.15, -.85),
        'PXDN': (-2.90, -2.05),
        'RSPO3': (1.25, 2.58),
        'ADAM12': (1.50, 2.04),
        'TAF9': (1.80, 1.39),
        'IFITM1': (-3.25, 1.08),
        'S100A14': (-3.40, .30),
    }
    label_data = de.loc[de.gene_name.isin(label_positions),
                        ['gene_id', 'gene_name', 'FDR', 'reporting_subset']].copy()
    assert len(label_data) == len(label_positions)
    label_data['smallest_pair_difference'] = eff.loc[label_data.gene_id, pair_cols].min(axis=1).to_numpy()
    label_data['largest_pair_difference'] = eff.loc[label_data.gene_id, pair_cols].max(axis=1).to_numpy()
    label_data['label_x'] = label_data.gene_name.map(lambda gene: label_positions[gene][0])
    label_data['label_y'] = label_data.gene_name.map(lambda gene: label_positions[gene][1])
    for row in label_data.itertuples():
        b.annotate(row.gene_name,
                   xy=(row.smallest_pair_difference, row.largest_pair_difference),
                   xytext=(row.label_x, row.label_y),
                   fontsize=6.3, color=DARK, ha='left', va='center',
                   bbox=dict(boxstyle='square,pad=0.08', fc='white', ec='none', alpha=.9),
                   arrowprops=dict(arrowstyle='-', color='#888888', lw=.55,
                                   shrinkA=2, shrinkB=2), zorder=5)
    panel(fig, b, 'c', '', dx=-.062, dy=.012)

    # ---- c: arm-level RNA versus lpWGS.
    au = arm[~arm.unit.str.startswith('chrY')].copy()
    c.axhline(0, lw=.6, color='#DDDDDD'); c.axvline(0, lw=.6, color='#DDDDDD')
    fit = np.polyfit(au.DNA_delta_LOY2X, au.RNA_mean_difference, 1)
    xs = np.linspace(au.DNA_delta_LOY2X.min(), au.DNA_delta_LOY2X.max(), 10)
    c.plot(xs, np.polyval(fit, xs), color='#888888', lw=.8)
    c.scatter(au.DNA_delta_LOY2X, au.RNA_mean_difference, s=au.genes / 12, color=BLUE, alpha=.72, lw=0)
    key = au[(au.DNA_delta_LOY2X.abs() > .08) | (au.RNA_mean_difference.abs() > .07)]
    texts = [c.text(r.DNA_delta_LOY2X, r.RNA_mean_difference, r.unit, fontsize=6.2, color=DARK) for r in key.itertuples()]
    adjust_text(texts, ax=c, only_move={'text': 'xy'}, arrowprops=dict(arrowstyle='-', color='#BBBBBB', lw=.45))
    c.set_xlabel('lpWGS arm median log$_2$ difference'); c.set_ylabel('RNA arm mean log$_2$ CPM difference')
    c.text(.03, .97, 'Pearson r\n  gene {:.2f}\n  10 Mb {:.2f}\n  arm {:.2f}'.format(
        conc.loc['gene', 'Pearson'], conc.loc['window_10Mb', 'Pearson'], conc.loc['arm', 'Pearson']),
        transform=c.transAxes, va='top', fontsize=6.4, color=DARK)
    panel(fig, c, 'd', '', dx=-.062, dy=.012)

    # ---- d: composition of the significant sets.
    sets = ['FDR<0.05', 'reporting_subset']
    cats = ['CN_stable', 'CN_changed_concordant', 'CN_changed_discordant']
    colours = {'CN_stable': '#BBBBBB', 'CN_changed_concordant': BLUE, 'CN_changed_discordant': '#F0B27A'}
    bottom = np.zeros(2)
    for cat in cats:
        vals = np.array([cnde.loc[s, cat] for s in sets], dtype=float)
        d.bar([0, 1], vals, bottom=bottom, width=.6, color=colours[cat],
              label={'CN_stable': 'Stable', 'CN_changed_concordant': 'Concordant', 'CN_changed_discordant': 'Discordant'}[cat])
        for i, (v, bm) in enumerate(zip(vals, bottom)):
            if v >= 5:
                d.text(i, bm + v / 2, int(v), ha='center', va='center', fontsize=6.6,
                       color='white' if cat != 'CN_stable' else DARK)
        for i, (value, base) in enumerate(zip(vals, bottom)):
            if 0 < value < 5:
                d.annotate(str(int(value)), (i, base+value/2), (i+.35, base+8), fontsize=6.2, ha='center', arrowprops=dict(arrowstyle='-', lw=.4, color=DARK))
        bottom += vals
    d.set_xticks([0, 1], ['FDR < 0.05\n(147)', f'|log2FC| \u22650.5\n+ consistency\n({int(thesis.sum())})'], fontsize=6.8)
    d.set_ylabel('Genes'); d.set_xlim(-.6, 1.6)
    d.legend(fontsize=5.8, loc='upper right', bbox_to_anchor=(1.08, 1.25), ncol=2, borderaxespad=0, handlelength=.8, columnspacing=.8)
    panel(fig, d, 'e', '', dx=-.062, dy=.012)

    # ---- e: experimental GRHL2 CPM; external DepMap context is supplementary.
    cpm = read('core_CPM_all_genes')
    gid = de.gene_id[de.gene_name.eq('GRHL2')].iloc[0]
    gc = cpm[cpm.gene_id.eq(gid)].copy()
    gc['treatment'] = np.where(gc.sample_key.str.contains('sgCTRL'), 'CTRL', 'LOY')
    for i, (t, sub) in enumerate(gc.groupby('treatment')):
        e.scatter(np.full(len(sub), i) + rng.uniform(-.08, .08, len(sub)), sub.CPM, s=26, color=COND[t], zorder=3)
        e.plot([i - .22, i + .22], [sub.CPM.mean()] * 2, color=COND[t], lw=1.4)
    e.set_xticks([0, 1], ['sgCTRL', 'Y-depleted'], fontsize=5.8, rotation=30)
    e.set_xlim(-.55, 1.55); e.set_ylim(0, 4.6)
    e.set_ylabel('GRHL2 (CPM)', fontsize=7)
    e.tick_params(axis='y', labelsize=6.5)
    e.text(.5, .97, f"FDR = {float(de.FDR[de.gene_name.eq('GRHL2')].iloc[0]):.1e}" + f"\nlog$_2$FC = {float(de.logFC[de.gene_name.eq('GRHL2')].iloc[0]):.2f}",
           transform=e.transAxes, ha='center', va='top', fontsize=6.5, color=DARK)
    panel(fig, e, 'f', '', dx=-.042, dy=.012)

    # g,h: saved TMM log2 CPM (edgeR prior count 2); no model refitting.
    logcpm = read('core_logCPM_all_genes')
    gene_tables = []
    for letter, gene, xpos in [('g', 'TAF9', .65), ('h', 'PXDN', .855)]:
        ax = fig.add_axes([xpos, .08, .12, .22])
        result = de.loc[de.gene_name.eq(gene)].iloc[0]
        values = logcpm.loc[logcpm.gene_id.eq(result.gene_id)].copy()
        assert len(values) == 6
        values['gene_name'] = gene
        values['treatment'] = np.where(values.sample_key.str.contains('sgCTRL'), 'CTRL', 'LOY')
        values['passage'] = values.sample_key.str.extract(r'(P1[234])$')[0]
        for i, treatment in enumerate(['CTRL', 'LOY']):
            sub = values.loc[values.treatment.eq(treatment)].sort_values('passage')
            assert len(sub) == 3
            colour = COND[treatment]
            ax.boxplot([sub.log2_CPM], positions=[i], widths=.55, patch_artist=True,
                       showfliers=False, whis=1.5,
                       boxprops=dict(facecolor=colour, alpha=.18, edgecolor=colour, linewidth=.8),
                       medianprops=dict(color=colour, linewidth=1.2),
                       whiskerprops=dict(color=colour, linewidth=.8),
                       capprops=dict(color=colour, linewidth=.8))
            ax.scatter(i + np.array([-.09, 0, .09]), sub.log2_CPM, s=18,
                       color=colour, edgecolor='white', linewidth=.3, zorder=3)
        span = max(values.log2_CPM.max() - values.log2_CPM.min(), .5)
        ax.set_ylim(values.log2_CPM.min() - .15*span, values.log2_CPM.max() + .5*span)
        ax.set_xlim(-.55, 1.55)
        ax.set_xticks([0, 1], ['sgCTRL', 'Y-depleted'], fontsize=5.8, rotation=30)
        ax.tick_params(axis='y', labelsize=6.5)
        ax.set_ylabel(gene + ' (log$_2$ CPM)', fontsize=7)
        ax.text(.5, .97, f'FDR = {result.FDR:.1e}' + f'\nlog$_2$FC = {result.logFC:.2f}',
                transform=ax.transAxes, ha='center', va='top', fontsize=6.2, color=DARK)
        panel(fig, ax, letter, '', dx=-.042, dy=.012)
        gene_tables.append(values)

    export_current(fig, 'ch03_main_du145_rna_global_structure', 'CH3_S3_4_6_3_F01_DU145_RNA_global', {
        'panel_a_pca':pca[['sample_key','treatment','passage_block','PC1','PC2']],
        'panel_a_variance':var, 'panel_a_design_associations':assoc,
        'panel_b_c_gene_level':de[['gene_id','gene_name','chrom','t2t_arm','logFC','PValue','FDR','reporting_subset','pairs_agreeing_with_logFC','CN_class']],
        'panel_c_pair_differences':eff.loc[de.gene_id,pair_cols].reset_index(),
        'panel_c_gene_labels':label_data, 'panel_d_arm_level':au,
        'panel_e_cn_summary':cnde.reset_index(), 'panel_f_GRHL2_CPM':gc,
        'panel_g_h_TAF9_PXDN_logCPM':pd.concat(gene_tables, ignore_index=True)})


def draw_supp():
    de = annotate(read('DE_paired_voom_primary'))
    eff = read('DU145_LOY_vs_sgCTRL_primary_effects')
    unp = read('DE_unpaired_voom_sensitivity')[['gene_id', 'logFC', 'FDR']].rename(columns={'logFC': 'logFC_unpaired', 'FDR': 'FDR_unpaired'})
    cn = read('CN_gene_level')
    win = read('CN_window_10Mb_RNA_vs_DNA')
    arm = read('CN_arm_level_RNA_vs_DNA')
    byarm = read('CN_FDR05_genes_by_arm')
    cpm = read('core_CPM_all_genes')
    PAIRS = ['pair_P12', 'pair_P13', 'pair_P14']

    g = cn.merge(de[['gene_id', 'gene_name']], on='gene_id', how='left', suffixes=('', '_de'))
    fig = plt.figure(figsize=(183 * MM, 150 * MM))
    a = fig.add_axes([0.075, 0.585, 0.245, 0.335])
    b = fig.add_axes([0.415, 0.585, 0.245, 0.335])
    c = fig.add_axes([0.755, 0.585, 0.205, 0.335])
    d = fig.add_axes([0.075, 0.085, 0.245, 0.335])
    e = fig.add_axes([0.415, 0.085, 0.545, 0.335])

    # a: 10-Mb windows.
    wa = win[~win.unit.str.startswith('chrY')]
    a.axhline(0, lw=.6, color='#DDDDDD'); a.axvline(0, lw=.6, color='#DDDDDD')
    a.scatter(wa.DNA_delta_LOY2X, wa.RNA_mean_difference, s=8, color=BLUE, alpha=.6, lw=0)
    fit = np.polyfit(wa.DNA_delta_LOY2X, wa.RNA_mean_difference, 1)
    xs = np.linspace(wa.DNA_delta_LOY2X.min(), wa.DNA_delta_LOY2X.max(), 10)
    a.plot(xs, np.polyval(fit, xs), color='#888888', lw=.8)
    r = np.corrcoef(wa.DNA_delta_LOY2X, wa.RNA_mean_difference)[0, 1]
    a.text(.97, .04, f'r = {r:.2f}, slope {fit[0]:.2f}\n{len(wa)} windows', transform=a.transAxes, ha='right', va='bottom', fontsize=6.4, color=DARK)
    # Label the ten most extreme windows (combined standardised distance from the origin), one label per
    # chromosome: anchored at that chromosome's most extreme window, naming the span when windows are contiguous.
    wa = wa.assign(extremity=np.hypot(wa.DNA_delta_LOY2X / wa.DNA_delta_LOY2X.std(), wa.RNA_mean_difference / wa.RNA_mean_difference.std()))
    ext = wa.nlargest(10, 'extremity').assign(chrom=lambda d: d.unit.str.split(':').str[0],
                                              start=lambda d: d.unit.str.split(':').str[1].str.removesuffix('Mb').astype(int))
    labels = []
    for ch, grp in ext.groupby('chrom', sort=False):
        starts = sorted(grp.start)
        contiguous = all(b_ - a_ == 10 for a_, b_ in zip(starts, starts[1:]))
        text = f'{ch} {starts[0]}–{starts[-1] + 10} Mb' if contiguous else f'{ch}, {len(starts)} windows'
        anchor = grp.loc[grp.extremity.idxmax()]
        labels.append(dict(chrom=ch, label=text, x=anchor.DNA_delta_LOY2X, y=anchor.RNA_mean_difference, windows=';'.join(grp.unit)))
    labels = pd.DataFrame(labels)
    texts = [a.text(r.x, r.y, r.label, fontsize=5.3, color=DARK, zorder=5,
                    path_effects=[pe.withStroke(linewidth=1.4, foreground='white')]) for r in labels.itertuples()]
    adjust_text(texts, ax=a, only_move={'text': 'xy'}, force_text=(.4, .7), arrowprops=dict(arrowstyle='-', color='#BBBBBB', lw=.45))
    a.set_xlabel('lpWGS 10-Mb median log$_2$ difference'); a.set_ylabel('RNA 10-Mb mean log$_2$ CPM difference')
    panel(fig, a, 'a', '', dx=-.062)

    # b: gene level, the noisiest unit.
    gl = g[np.isfinite(g.CN_delta_LOY2X)]
    b.axhline(0, lw=.6, color='#DDDDDD'); b.axvline(0, lw=.6, color='#DDDDDD')
    b.scatter(gl.CN_delta_LOY2X, gl.mean_difference, s=2, color=GREY, alpha=.4, lw=0, rasterized=True)
    rb = np.corrcoef(gl.CN_delta_LOY2X, gl.mean_difference)[0, 1]
    b.text(.03, .96, f'r = {rb:.2f}\n{len(gl):,} genes', transform=b.transAxes, va='top', fontsize=6.4, color=DARK)
    b.set_xlabel('lpWGS gene-level log$_2$ difference'); b.set_ylabel('RNA mean log$_2$ CPM difference')
    b.set_ylim(-3, 3); b.set_xlim(-1.6, 1.6)
    outside = int((gl.CN_delta_LOY2X.abs() > 1.6).sum() + (gl.mean_difference.abs() > 3).sum())
    b.text(.97, .04, f'{outside} points outside axes', transform=b.transAxes, ha='right', fontsize=6.0, color=DARK)
    panel(fig, b, 'b', '', dx=-.062)

    # c: which arms carry the significant genes.
    top = byarm.head(8).iloc[::-1]
    yc = np.arange(len(top))
    c.barh(yc, top.Freq, color=np.where(top.arm.eq('chr20p'), BLUE, GREY), height=.7)
    c.set_yticks(yc, top.arm, fontsize=6.6)
    c.set_xlabel('Genes at FDR < 0.05')
    c.text(.95, .06, '20p deletion', transform=c.transAxes, ha='right', fontsize=6.4, color=BLUE)
    panel(fig, c, 'c', '', dx=-.075)

    # d: the second Y-depleted lpWGS specimen.
    au = arm[~arm.unit.str.startswith('chrY')]
    d.axhline(0, lw=.6, color='#DDDDDD'); d.axvline(0, lw=.6, color='#DDDDDD')
    d.plot([-.9, .3], [-.9, .3], ls='--', lw=.6, color='#CCCCCC')
    d.scatter(au.DNA_delta_LOY2X, au.DNA_delta_DMSO_LOY1X, s=au.genes / 12, color=BLUE, alpha=.7, lw=0)
    key = au[(au.DNA_delta_LOY2X.abs() > .08) | (au.DNA_delta_DMSO_LOY1X.abs() > .08)]
    texts = [d.text(r.DNA_delta_LOY2X, r.DNA_delta_DMSO_LOY1X, r.unit, fontsize=6.2, color=DARK) for r in key.itertuples()]
    adjust_text(texts, ax=d, only_move={'text': 'xy'}, arrowprops=dict(arrowstyle='-', color='#BBBBBB', lw=.45))
    d.set_xlabel('LOY-2X minus sgCTRL (arm median log$_2$)')
    d.set_ylabel('DMSO-LOY-1X minus sgCTRL')
    panel(fig, d, 'd', '', dx=-.062)

    # e: paired versus unpaired model.
    m = de[['gene_id', 'gene_name', 'logFC', 'FDR', 'reporting_subset']].merge(unp, on='gene_id')
    sig_p = m.FDR < .05
    e.scatter(m.logFC[~sig_p], m.logFC_unpaired[~sig_p], s=3, color='#DDDDDD', lw=0, rasterized=True)
    e.scatter(m.logFC[sig_p], m.logFC_unpaired[sig_p], s=10, color='#F0B27A', lw=0, label='Paired FDR < 0.05')
    th = m.reporting_subset
    e.scatter(m.logFC[th], m.logFC_unpaired[th], s=8, color=RED, lw=0, label=LABEL)
    lim = (-5.6, 3.0)
    e.plot(lim, lim, ls='--', lw=.6, color='#AAAAAA')
    e.set_xlim(lim); e.set_ylim(lim)
    e.set_xlabel('Paired model log$_2$FC (~ passage + condition)')
    e.set_ylabel('Unpaired model log$_2$FC')
    e.text(.03, .96, 'Effect sizes agree closely (Pearson r = {:.3f});\nno gene reaches FDR < 0.05 without passage blocking'.format(
        np.corrcoef(m.logFC, m.logFC_unpaired)[0, 1]), transform=e.transAxes, va='top', fontsize=6.5, color=DARK)
    e.legend(fontsize=6.3, loc='lower right', borderaxespad=.3)
    panel(fig, e, 'e', '', dx=-.045)

    labels_windows=labels.copy()
    model_comparison=m.copy()
    old_height = fig.get_figheight()
    new_height = old_height + 3.05
    fig.set_size_inches(fig.get_figwidth(), new_height)
    for ax in fig.axes:
        pos=ax.get_position()
        ax.set_position([pos.x0, (pos.y0*old_height+3.05)/new_height, pos.width, pos.height*old_height/new_height])
    for text in fig.texts:
        x,y=text.get_position(); text.set_position((x,(y*old_height+3.05)/new_height))
    b1=fig.add_axes([.075,.105,.23,.195])
    b2=fig.add_axes([.42,.105,.21,.195])
    c=fig.add_axes([.755,.105,.17,.195])
    pca, var, assoc, qc, cor = (read(n) for n in ['PCA_core', 'PCA_variance_core', 'PCA_design_associations', 'library_QC', 'sample_correlations_core'])
    assert len(pca) == 6 and len(qc) == 6 and len(cor) == 36
    v = var.set_index('PC').core_variance
    a1 = assoc.set_index('axis')

    # b: library composition, fixed order sgCTRL then Y-depleted.
    qc['label'] = qc.sample_key.map(short)
    qc = qc.sort_values(['treatment', 'passage_block'])
    x = np.arange(len(qc)); cols = [COND[t] for t in qc.treatment]
    for condition, group in pca.groupby('treatment'):
        group = group.sort_values('passage_block')
        b1.plot(range(3), group.PC1, color=COND[condition], lw=.9)
        for i, row in enumerate(group.itertuples()):
            b1.scatter(i, row.PC1, color=COND[condition], marker=PASSAGE_MARKER[row.passage_block],
                       s=22, edgecolors='white', linewidths=.35, zorder=3)
    b1.set_ylabel('PC1 score')
    b1.set_xticks(range(3), ['P12', 'P13', 'P14'], fontsize=6.5)
    b1.set_xlim(-.2, 2.2); b1.set_ylim(-45, 90)
    b1.set_yticks([0, 50])
    b1.text(.02, .98, 'Passage: R\u00b2 = 0.93', transform=b1.transAxes,
            va='top', fontsize=6.5, color=DARK)
    b2.bar(x, 100 * qc.MSY_fraction, color=cols, width=.72)
    b2.set_ylabel('MSY reads\n(% assigned)'); b2.set_xticks(x, qc.label, rotation=90, fontsize=6.5)
    b2.set_xlim(-.6, len(qc) - .4)
    b2.text(.98, .97, 'Y-depleted ÷\nsgCTRL = {:.2f}'.format(qc[qc.treatment == 'LOY'].MSY_fraction.mean() / qc[qc.treatment == 'CTRL'].MSY_fraction.mean()),
            transform=b2.transAxes, ha='right', va='top', fontsize=6.0, color=DARK)
    panel(fig, b1, 'f', '', dx=-.062, y=.332)
    panel(fig, b2, 'g', '', dx=-.062, y=.332)

    # c: Pearson correlation of TMM log2 CPM, average-linkage order (order saved).
    m = cor.pivot(index='sample_1', columns='sample_key', values='correlation')
    dist = 1 - m.to_numpy(); dist = (dist + dist.T) / 2; np.fill_diagonal(dist, 0)
    tree = linkage(squareform(dist), 'average', optimal_ordering=True)
    order = list(m.index[leaves_list(tree)])
    m = m.loc[order, order]
    im = c.imshow(m.to_numpy(), cmap='magma', vmin=float(m.min().min()), vmax=1, interpolation='nearest')
    labels = [('Y' if 'LOY' in s else 'C') + s.split('-')[-1][1:] for s in order]
    c.set_xticks(range(6), labels, rotation=90, fontsize=6.5); c.set_yticks(range(6), labels, fontsize=6.5)
    for tick, s in zip(c.get_yticklabels(), order):
        tick.set_color(COND['LOY' if 'LOY' in s else 'CTRL'])
    for tick, s in zip(c.get_xticklabels(), order):
        tick.set_color(COND['LOY' if 'LOY' in s else 'CTRL'])
    c.tick_params(length=0)
    for s in c.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, cax=fig.add_axes([.937,.119,.008,.145]))
    cb.set_label(''); cb.ax.set_title('r', fontsize=6, pad=4); cb.ax.tick_params(labelsize=6.5)
    # Align dendrogram leaves precisely with the square heatmap columns.
    fig.canvas.draw()
    box = c.get_position()
    dend_ax = fig.add_axes([box.x0, box.y1 + .004, box.width, .047])
    dendrogram(tree, ax=dend_ax, no_labels=True, color_threshold=0,
               above_threshold_color='#777777', link_color_func=lambda _: '#777777')
    dend_ax.set_axis_off()
    panel(fig, c, 'h', '', dx=-.062, y=.332)


    c.text(.5, -.57, 'C, sgCTRL; Y, Y-depleted', transform=c.transAxes, ha='center', fontsize=5.7)
    export_current(fig, 'ch03_supp_du145_rna_copy_number', 'CH3_S3_4_6_3_S02_DU145_RNA_copy_number', {
        'panel_a_windows':wa,'panel_a_window_labels':labels_windows,
        'panel_b_gene_level':gl[['gene_id','gene_name','t2t_arm','CN_delta_LOY2X','mean_difference','CN_class']],
        'panel_c_genes_by_arm':byarm,'panel_d_arm_two_specimens':au,
        'panel_e_paired_vs_unpaired':model_comparison,
        'panel_f_pc1_passage':pca[['sample_key','treatment','passage_block','PC1']],
        'panel_g_library_qc':qc[['sample_key','treatment','passage_block','assigned_reads','MSY_fraction']],
        'panel_h_correlation_ordered':m.reset_index(),
        'panel_h_linkage':pd.DataFrame(tree,columns=['left','right','distance','n_leaves'])})


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only',choices=['main','supp','both'],default='both')
    args=parser.parse_args()
    if args.only in ['main','both']: draw_main()
    if args.only in ['supp','both']: draw_supp()
