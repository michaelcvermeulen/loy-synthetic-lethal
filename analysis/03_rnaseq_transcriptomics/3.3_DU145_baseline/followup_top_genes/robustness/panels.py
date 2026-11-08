"""Marker panels: interferon/MHC-I, neuroendocrine, lung lineage, GRHL family. Reads only."""
from pathlib import Path
import pandas as pd
T = Path(__file__).resolve().parents[2] / 'outputs/tables'
eff = pd.read_csv(T / 'DU145_LOY_vs_sgCTRL_all_gene_effects.csv')
de = pd.read_csv(T / 'DE_paired_voom_primary.csv', usecols=['gene_id', 'logFC', 'FDR'])
eff = eff.merge(de, on='gene_id', how='left').sort_values('mean_CPM', ascending=False).drop_duplicates('gene_name')
P = {
 'interferon_ISG': ['STAT1','STAT2','IRF1','IRF9','IFI6','IFI27','IFI44','IFI44L','ISG15','IFIT1','IFIT3','MX1','OAS1','OAS2','CXCL10','CXCL11','IDO1','GBP1','GBP4','BST2','IFITM1'],
 'MHC_I_antigen': ['B2M','HLA-A','HLA-B','HLA-C','HLA-E','TAP1','TAP2','TAPBP','NLRC5','PSMB8','PSMB9','CD274'],
 'neuroendocrine': ['CHGA','CHGB','SYP','ENO2','NCAM1','ASCL1','INSM1','SOX2','POU3F2','NEUROD1','SST','UCHL1'],
 'lung_lineage': ['NKX2-1','NAPSA','SFTPA1','SFTPA2','SFTPB','SFTPC','SCGB3A2','HOPX','KRT7'],
 'GRHL_family': ['GRHL1','GRHL2','GRHL3','OVOL1','OVOL2','ESRP1','ESRP2','RAB25','CLDN4','CLDN7','EPCAM','CDH1'],
}
rows = []
for panel, genes in P.items():
    for g in genes:
        r = eff[eff['gene_name'] == g]
        if not len(r):
            rows.append(dict(panel=panel, gene=g, status='absent_from_annotation')); continue
        r = r.iloc[0]
        rows.append(dict(panel=panel, gene=g, CTRL_CPM=round(r['CTRL_mean_CPM'], 2), LOY_CPM=round(r['LOY_mean_CPM'], 2),
                         mean_diff=round(r['mean_difference'], 3), pairs_same=int(r['pairs_same_direction_as_mean']),
                         logFC=None if pd.isna(r['logFC']) else round(r['logFC'], 3),
                         FDR=None if pd.isna(r['FDR']) else float(f"{r['FDR']:.2g}"),
                         status='tested' if not pd.isna(r['FDR']) else ('below_filter' if r['CTRL_mean_CPM'] + r['LOY_mean_CPM'] < 2 else 'not_in_primary')))
out = pd.DataFrame(rows); out.to_csv(Path(__file__).parent / 'marker_panels.csv', index=False)
for p, d in out.groupby('panel', sort=False):
    print(f'== {p} ==  median mean_diff (tested) = {d.loc[d.status=="tested","mean_diff"].median():.3f}; FDR<0.05: {int((d.FDR<0.05).sum())}/{int((d.status=="tested").sum())}')
    print(d.drop(columns='panel').to_string(index=False))
