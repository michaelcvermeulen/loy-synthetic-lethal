"""Per-gene robustness for the top DU145 DE genes.

For each focus gene: raw counts and CPM per library, passage-pair differences,
leave-one-pair-out direction, paired/unpaired voom statistics, lpWGS copy-number
class, and agreement with the earlier (superseded) LNCAP analysis branches.
Reads only; writes gene_robustness.csv next to this script.
"""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                       # SEPT_23_2026_REANALYZE
T = ROOT / 'outputs/tables'
LNCAP = ROOT.parent                          # analysis/LNCAP

FOCUS = ['GRHL2', 'TACSTD2', 'AP1M2', 'CLDN7', 'EPS8L1', 'SCEL', 'MAL2', 'TMC4', 'ST14', 'OVOL1', 'S100A14',
         'LCN2', 'CXCL1', 'BST2', 'LTB', 'SFTPA1', 'NAPSA', 'HOXB9', 'PKIA', 'CCR9', 'PXDN', 'IVL', 'KDR',
         'ACSL5', 'MYBPH', 'BICDL2', 'GOLGA7B', 'STEAP4', 'THBD',
         'TAF9', 'F13A1', 'COL4A2', 'CREB3L1', 'RSPO3', 'ADAM12', 'CNTN1', 'TMEM45A', 'SLC16A2',
         'ZEB1', 'ZEB2', 'CDH1', 'CDH2', 'VIM', 'ESRP1', 'ESRP2', 'OVOL2', 'SNAI2', 'STAT1']

de = pd.read_csv(T / 'DE_paired_voom_primary.csv')
unp = pd.read_csv(T / 'DE_unpaired_voom_sensitivity.csv', usecols=['gene_id', 'logFC', 'FDR'])
eff = pd.read_csv(T / 'DU145_LOY_vs_sgCTRL_all_gene_effects.csv')
cn = pd.read_csv(T / 'CN_gene_level.csv', usecols=['gene_id', 'CN_class', 'CN_delta_LOY2X', 't2t_arm'])
lc = pd.read_csv(T / 'core_logCPM_all_genes.csv')
cpm = pd.read_csv(T / 'core_CPM_all_genes.csv')
x = pd.read_rds if False else None  # raw counts come from the source matrix below
raw = pd.read_csv(LNCAP / 'count_matrix/LNCAP_counts_raw.csv').set_index('gene_id')

eff = eff[eff['gene_name'].isin(FOCUS)]
# Keep the representative locus per symbol (highest mean CPM), matching the analysis mapping rule.
eff = eff.sort_values('mean_CPM', ascending=False).drop_duplicates('gene_name')
samples = ['LNCaP-sgCTRL-P12', 'LNCaP-sgCTRL-P13', 'LNCaP-sgCTRL-P14', 'LNCaP-LOY-P12', 'LNCaP-LOY-P13', 'LNCaP-LOY-P14']
cpm_w = cpm.pivot(index='gene_id', columns='sample_key', values='CPM')

rows = []
for _, g in eff.iterrows():
    gid = g['gene_id']
    r = {'gene': g['gene_name'], 'chrom': g['chrom'], 'primary_eligible': g['primary_eligible'],
         'mean_CPM': g['mean_CPM'], 'mean_difference': g['mean_difference']}
    rc = raw.loc[gid, samples] if gid in raw.index else pd.Series(index=samples, dtype=float)
    r['raw_sgCTRL'] = '/'.join(str(int(v)) for v in rc[samples[:3]])
    r['raw_LOY'] = '/'.join(str(int(v)) for v in rc[samples[3:]])
    r['min_raw_any'] = int(rc.min())
    r['CPM_sgCTRL'] = '/'.join(f'{cpm_w.loc[gid, s]:.1f}' for s in samples[:3])
    r['CPM_LOY'] = '/'.join(f'{cpm_w.loc[gid, s]:.1f}' for s in samples[3:])
    for p in ['P12', 'P13', 'P14']:
        r[f'pair_{p}'] = g[f'pair_{p}']
    r['pairs_same_sign'] = int(g['pairs_same_direction_as_mean'])
    r['LOPO_all_same_direction'] = g['LOPO_all_same_direction']
    r['LOPO_min'] = g['LOPO_min']; r['LOPO_max'] = g['LOPO_max']
    d = de[de['gene_id'] == gid]
    r['paired_logFC'] = d['logFC'].iloc[0] if len(d) else None
    r['paired_FDR'] = d['FDR'].iloc[0] if len(d) else None
    r['thesis_DE'] = d['thesis_DE'].iloc[0] if len(d) else None
    u = unp[unp['gene_id'] == gid]
    r['unpaired_FDR'] = u['FDR'].iloc[0] if len(u) else None
    c = cn[cn['gene_id'] == gid]
    r['CN_class'] = c['CN_class'].iloc[0] if len(c) else None
    r['CN_delta'] = c['CN_delta_LOY2X'].iloc[0] if len(c) else None
    r['arm'] = c['t2t_arm'].iloc[0] if len(c) else None
    rows.append(r)
out = pd.DataFrame(rows)

# Earlier branches (superseded; used only as sensitivity). Join by unversioned ENSG, then symbol.
def stable(s): return s.str.replace(r'\.\d+$', '', regex=True)
out = out.merge(eff[['gene_name', 'gene_id']].rename(columns={'gene_name': 'gene'}), on='gene')
out['ensg'] = stable(out['gene_id'])
qc = pd.read_csv(LNCAP / 'QC_reprocess_2026/tables/DE_LOY_vs_CTRL_paired_DESeq2.csv')
idcol = [c for c in qc.columns if c.lower() in ('gene_id', 'geneid', 'ensembl', 'id')]
qc['ensg'] = stable(qc[idcol[0]].astype(str)) if idcol else None
qc = qc[['ensg', 'log2FoldChange', 'padj']].rename(columns={'log2FoldChange': 'QCreprocess_DESeq2_LFC', 'padj': 'QCreprocess_DESeq2_padj'})
rob = pd.read_csv(LNCAP / 'DE_extended/lncap_LOY_v_sgCTRL_DE_robust.csv')
rob['ensg'] = stable(rob['gene_id'].astype(str))
rob = rob[['ensg', 'n_pass', 'FDR_edgeRQL', 'FDR_treat', 'confident_DEG']].rename(columns={'n_pass': 'v86_frameworks_passed_of4'})
out = out.merge(qc, on='ensg', how='left').merge(rob, on='ensg', how='left')

# Flags
out['low_count_flag'] = out['min_raw_any'] < 10
out['CN_flag'] = out['CN_class'].isin(['CN_changed_concordant', 'CN_changed_discordant'])
out = out.drop(columns=['ensg']).sort_values('paired_FDR')
out.to_csv(HERE / 'gene_robustness.csv', index=False)
print(out.to_string(index=False, max_colwidth=22))
