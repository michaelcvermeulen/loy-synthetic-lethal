"""One-off archival and renderer construction for the approved F01/F03 merge."""
from pathlib import Path
import hashlib, json, shutil, textwrap, re

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
CODE = FIG / 'code'
ARC = ROOT / 'archive/DU145_F01_F03_combined_20261009'
SRC = FIG / 'source_data/du145_global_expression_combined'
ARC.mkdir(parents=True, exist_ok=True)
stems = ['ch03_main_du145_rna_global_structure', 'ch03_main_du145_rna_expression_changes', 'ch03_supp_du145_rna_copy_number']
folders = ['CH3_S3_4_6_3_F01_DU145_RNA_global', 'CH3_S3_4_6_3_F03_DU145_RNA_expression', 'CH3_S3_4_6_3_S02_DU145_RNA_copy_number']
paths = [ROOT/p for p in ['Results_edit.md', 'Chapter3_figure_captions.md', 'AGENTS.md', 'figures/chapter3/selected/selected_figure_index.tsv']]
paths += [CODE/p for p in ['make_du145_rna_global_structure.py','make_du145_rna_expression_changes.py','make_du145_rna_supp.py']]
for stem in stems:
    paths += list(FIG.glob(stem + '.*'))
    src = FIG / 'source_data' / stem.replace('ch03_main_', '').replace('ch03_supp_', '')
    paths += [p for p in src.rglob('*') if p.is_file()]
for folder in folders:
    paths += [p for p in (FIG/'selected'/folder).rglob('*') if p.is_file()]
records=[]
for p in paths:
    dst=ARC/p.relative_to(ROOT)
    if not dst.exists():
        dst.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(p,dst)
    records.append({'path':p.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(dst.read_bytes()).hexdigest()})
(ARC/'before.json').write_text(json.dumps(records,indent=2),encoding='utf-8')

upstream=ROOT.parent/'RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE/outputs/tables'
names=['PCA_core','PCA_variance_core','PCA_design_associations','library_QC','sample_correlations_core',
       'DE_paired_voom_primary','DU145_LOY_vs_sgCTRL_primary_effects','CN_gene_level','CN_arm_level_RNA_vs_DNA',
       'CN_RNA_DNA_concordance','CN_DE_summary','core_CPM_all_genes','DE_unpaired_voom_sensitivity',
       'CN_window_10Mb_RNA_vs_DNA','CN_FDR05_genes_by_arm']
(SRC/'inputs').mkdir(parents=True,exist_ok=True)
inputs=[]
for name in names:
    p=upstream/(name+'.csv'); dst=SRC/'inputs'/p.name
    if not dst.exists(): shutil.copy2(p,dst)
    assert dst.read_bytes()==p.read_bytes()
    inputs.append({'original':str(p),'local':dst.relative_to(SRC).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(SRC/'inputs.json').write_text(json.dumps(inputs,indent=2),encoding='utf-8')

old=lambda name:(ARC/'figures/chapter3/code'/name).read_text(encoding='utf-8')
global_src=old('make_du145_rna_global_structure.py')
expr=old('make_du145_rna_expression_changes.py')
supp=old('make_du145_rna_supp.py')
global_load=global_src[global_src.index('pca, var,'):global_src.index('fig = plt.figure')]
pca_draw=global_src[global_src.index('# a: PCA'):global_src.index('# b: library composition')]
pca_draw=pca_draw.replace('for p, d in pca.groupby', 'for passage, grp in pca.groupby').replace("d = d.set_index('treatment')", "grp = grp.set_index('treatment')").replace('d.loc[', 'grp.loc[')
pca_draw=re.sub(r'\ba\.', 'p.', pca_draw).replace('panel(fig, a,', 'panel(fig, p,')
pca_draw=pca_draw.replace("bbox_to_anchor=(0.02, 0.012), bbox_transform=fig.transFigure", "bbox_to_anchor=(0.0, -0.24), bbox_transform=p.transAxes")
pca_draw=pca_draw.replace("panel(fig, p, 'a', '', y=.87)", "panel(fig, p, 'a', '', dx=-.062, dy=.012)")
body=expr[expr.index('de = annotate'):expr.index('export(fig, STEM,')]
body=body.replace("panels_csv = pd.read_csv(FOLLOW / 'robustness/marker_panels.csv')\n",'')
body=body.replace("g = panels_csv[panels_csv.gene.eq('GRHL2') & panels_csv.panel.eq('GRHL_family')].iloc[0]\n",'')
start=body.index('fig = plt.figure'); end=body.index('# ---- a: volcano')
layout='''fig = plt.figure(figsize=(7.2, 9.0))
p = fig.add_axes([.10, .745, .34, .22])
a = fig.add_axes([.60, .745, .37, .22])
b = fig.add_axes([.10, .400, .34, .24])
c = fig.add_axes([.60, .400, .37, .24])
d = fig.add_axes([.10, .075, .34, .215])
e = fig.add_axes([.65, .075, .27, .215])
'''
body=body[:start]+layout+global_load+pca_draw+body[end:]
for old_letter,new_letter in zip('abcde','bcdef'):
    body=re.sub(r"panel\(fig, "+old_letter+r", '"+old_letter+r"', '', dx=[^)]+\)",
                f"panel(fig, {old_letter}, '{new_letter}', '', dx=-.062, dy=.012)",body)
body=body.replace("fontsize=6.3, loc='upper left'", "fontsize=5.7, loc='upper left'")
body=body.replace("label='FDR < 0.05, copy-number concordant'", "label='FDR < 0.05, CN-concordant'")
body=body.replace("a.set_ylim(0, -np.log10(de.PValue.min()) * 1.55)", "a.set_ylim(0, -np.log10(de.PValue.min()) * 1.72)")
body=body.replace("d.legend(fontsize=6.2, loc='upper right', borderaxespad=.2, handlelength=.8, title='Copy number', title_fontsize=6.2)",
                  "d.legend(fontsize=6.2, loc='upper right', bbox_to_anchor=(1.03, 1.22), ncol=3, borderaxespad=0, handlelength=.8, columnspacing=.8)")
body=body.replace("if v:\n            d.text", "if v >= 5:\n            d.text")
body=body.replace("    bottom += vals", "    for i, (value, base) in enumerate(zip(vals, bottom)):\n        if 0 < value < 5:\n            d.annotate(str(int(value)), (i, base+value/2), (i+.35, base+8), fontsize=6.2, ha='center', arrowprops=dict(arrowstyle='-', lw=.4, color=DARK))\n    bottom += vals")
body=body.replace("texts = [a.text(r.logFC, -np.log10(r.PValue), r.gene_name, fontsize=6.3, color=DARK) for r in lab.itertuples()]", "")
pos_start=body.index('adjust_text(texts, ax=a,')
pos_end=body.index("a.set_xlabel(",pos_start)
body=body[:pos_start]+'''volcano_labels = {
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
'''+body[pos_end:]
# Preserve the same numerical sources and displayed gene selection.
main_export='''export_current(fig, 'ch03_main_du145_rna_global_structure', 'CH3_S3_4_6_3_F01_DU145_RNA_global', {
    'panel_a_pca':pca[['sample_key','treatment','passage_block','PC1','PC2']],
    'panel_a_variance':var, 'panel_a_design_associations':assoc,
    'panel_b_c_gene_level':de[['gene_id','gene_name','chrom','t2t_arm','logFC','PValue','FDR','reporting_subset','pairs_agreeing_with_logFC','CN_class']],
    'panel_c_pair_differences':eff.loc[de.gene_id,pair_cols].reset_index(),
    'panel_c_gene_labels':label_data, 'panel_d_arm_level':au,
    'panel_e_cn_summary':cnde.reset_index(), 'panel_f_GRHL2_CPM':gc})
'''

supp_load=supp[supp.index('de = annotate'):supp.index('# ================= Supplementary 1:')]
supp_body=supp[supp.index('g = cn.merge'):supp.index("export(fig, 'ch03_supp_du145_rna_copy_number',")]
# Retain the first five panels at their original physical dimensions; append one row.
supp_extra='''old_height = fig.get_figheight()
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
'''+global_load
extra=global_src[global_src.index('# b: library composition'):global_src.index('plotted = pca')]
extra=extra.replace("panel(fig, b1, 'b', '', dx=-.085, y=.87)", "panel(fig, b1, 'f', '', dx=-.062, y=.332)\npanel(fig, b2, 'g', '', dx=-.062, y=.332)")
extra=extra.replace("panel(fig, c, 'c', '', dx=-.09, y=.87)", "panel(fig, c, 'h', '', dx=-.062, y=.332)")
extra=extra.replace("fig.colorbar(im, ax=c, fraction=.05, pad=.03)", "fig.colorbar(im, cax=fig.add_axes([.937,.119,.008,.145]))")
extra=extra.replace("[box.x0, box.y1 + .005, box.width, .14]", "[box.x0, box.y1 + .004, box.width, .047]")
extra=extra.replace("labels = [short(s) for s in order]", "labels = [('Y' if 'LOY' in s else 'C') + s.split('-')[-1][1:] for s in order]")
extra=extra.replace("cb.set_label('Pearson r', fontsize=7)", "cb.set_label(''); cb.ax.set_title('r', fontsize=6, pad=4)")
extra=extra.replace("b2.text(.99, .62,", "b2.text(.98, .97,")
extra=extra.replace("fontsize=6.5, color=DARK)\npanel(fig, b1", "fontsize=6.0, color=DARK)\npanel(fig, b1")
extra += "\nc.text(.5, -.57, 'C, sgCTRL; Y, Y-depleted', transform=c.transAxes, ha='center', fontsize=5.7)\n"
# m is the model comparison above; retain before correlation drawing reuses its name.
supp_extra="model_comparison=m.copy()\n"+supp_extra+extra
supp_export='''export_current(fig, 'ch03_supp_du145_rna_copy_number', 'CH3_S3_4_6_3_S02_DU145_RNA_copy_number', {
    'panel_a_windows':wa,'panel_a_window_labels':labels_windows,
    'panel_b_gene_level':gl[['gene_id','gene_name','t2t_arm','CN_delta_LOY2X','mean_difference','CN_class']],
    'panel_c_genes_by_arm':byarm,'panel_d_arm_two_specimens':au,
    'panel_e_paired_vs_unpaired':model_comparison,
    'panel_f_pc1_passage':pca[['sample_key','treatment','passage_block','PC1']],
    'panel_g_library_qc':qc[['sample_key','treatment','passage_block','assigned_reads','MSY_fraction']],
    'panel_h_correlation_ordered':m.reset_index(),
    'panel_h_linkage':pd.DataFrame(tree,columns=['left','right','distance','n_leaves'])})
'''
header='''"""Approved DU145 F01/F03 merge and S02 extension; frozen inputs, no model refits.

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
        'scope':'Layout only; retained numerical values, thresholds, statistical families and selected genes.',
        'exports':outputs},indent=2),encoding='utf-8')
    print(stem,flush=True)

'''
renderer=header+'def draw_main():\n'+textwrap.indent(body+main_export,'    ')+'\n\ndef draw_supp():\n'+textwrap.indent(supp_load+supp_body+'labels_windows=labels.copy()\n'+supp_extra+supp_export,'    ')
renderer+='''\n\nif __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only',choices=['main','supp','both'],default='both')
    args=parser.parse_args()
    if args.only in ['main','both']: draw_main()
    if args.only in ['supp','both']: draw_supp()
'''
(CODE/'make_du145_global_expression_combined.py').write_text(renderer,encoding='utf-8')
print(f'Archived {len(paths)} files; froze {len(inputs)} inputs; created renderer.')
