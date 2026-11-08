"""Check numerical preservation and document references for the DU145 figure merge."""
from pathlib import Path
import hashlib,json,re
import pandas as pd
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
SRC=FIG/'source_data/du145_global_expression_combined'
ARC=ROOT/'archive/DU145_F01_F03_combined_20261009'
old=ARC/'figures/chapter3/source_data'
pairs=[
('du145_rna_global_structure/panel_a_pca.csv','main/panel_a_pca.csv'),
('du145_rna_global_structure/panel_a_design_associations.csv','main/panel_a_design_associations.csv'),
('du145_rna_expression_changes/panel_a_b_gene_level.csv','main/panel_b_c_gene_level.csv'),
('du145_rna_expression_changes/panel_b_pair_differences.csv','main/panel_c_pair_differences.csv'),
('du145_rna_expression_changes/panel_b_gene_labels.csv','main/panel_c_gene_labels.csv'),
('du145_rna_expression_changes/panel_c_arm_level.csv','main/panel_d_arm_level.csv'),
('du145_rna_expression_changes/panel_d_cn_summary.csv','main/panel_e_cn_summary.csv'),
('du145_rna_expression_changes/panel_e_GRHL2_CPM.csv','main/panel_f_GRHL2_CPM.csv'),
('du145_rna_global_structure/panel_b_pc1_passage.csv','supp/panel_f_pc1_passage.csv'),
('du145_rna_global_structure/panel_b_library_qc.csv','supp/panel_g_library_qc.csv'),
('du145_rna_global_structure/panel_c_correlation_ordered.csv','supp/panel_h_correlation_ordered.csv'),
('du145_rna_global_structure/panel_c_linkage.csv','supp/panel_h_linkage.csv')]
for name in ['panel_a_windows','panel_a_window_labels','panel_b_gene_level','panel_c_genes_by_arm','panel_d_arm_two_specimens','panel_e_paired_vs_unpaired']:
    pairs.append((f'du145_rna_copy_number/{name}.csv',f'supp/{name}.csv'))
for a,b in pairs:
    pd.testing.assert_frame_equal(pd.read_csv(old/a),pd.read_csv(SRC/b),check_exact=True)
for row in json.loads((SRC/'inputs.json').read_text(encoding='utf-8')):
    assert hashlib.sha256((SRC/row['local']).read_bytes()).hexdigest()==row['sha256']
    assert hashlib.sha256(Path(row['original']).read_bytes()).hexdigest()==row['sha256']
pdf=[]
for stem in ['ch03_main_du145_rna_global_structure','ch03_supp_du145_rna_copy_number']:
    reader=PdfReader(FIG/(stem+'.pdf'));assert len(reader.pages)==1
    page=reader.pages[0]
    pdf.append({'file':stem,'inches':[float(page.mediabox.width)/72,float(page.mediabox.height)/72]})
# The retired F03 files remain exact historical copies.
folder='CH3_S3_4_6_3_F03_DU145_RNA_expression'
for ext in ['pdf','svg','png']:
    rel=Path('figures/chapter3/selected')/folder/(folder+'.'+ext)
    assert (ROOT/rel).read_bytes()==(ARC/rel).read_bytes()
report={'matching_display_tables':len(pairs),'unchanged_source_inputs':15,'pdf':pdf,'retired_F03_preserved':True}
if (ARC/'adoption.json').exists():
    text=(ROOT/'Results_edit.md').read_text(encoding='utf-8')
    before=(ARC/'pre_adoption/Results_edit.md').read_text(encoding='utf-8')
    start='#### 3.3.4.3 ';end='#### 3.3.4.4 '
    assert text.split(start)[0]==before.split(start)[0]
    assert text.split(end,1)[1]==before.split(end,1)[1]
    section=text.split(start,1)[1].split(end,1)[0]
    assert 'CH3-S3.3.6.3-F03' not in section
    assert len(re.findall(r'\*\*Figure [^\n]+\(main\)',section))==3
    assert len(re.findall(r'\*\*Figure [^\n]+\(supp\)',section))==8
    for target in re.findall(r'\]\(([^)]+)\)',section):
        if not target.startswith('http'):
            path=Path(target) if re.match('[A-Za-z]:',target) else ROOT/target
            assert path.exists(),str(path)
    central=(ROOT/'Chapter3_figure_captions.md').read_text(encoding='utf-8')
    for code in ['F01','S02']:
        block=re.search(r'\*\*Figure CH3-S3.3.6.3-'+code+r'[^\n]+\n\n[^\n]+\n\n([\s\S]*?)(?=\n---)',section).group(1).strip()
        assert block in central,(code,'caption mismatch')
    assert '\ufffd' not in section
    report.update({'main_figure_blocks':3,'supp_figure_blocks':8,'unrelated_results_unchanged':True,'captions_match':True,'all_section_links_resolve':True})
(SRC/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
