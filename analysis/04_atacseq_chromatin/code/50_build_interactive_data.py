from pathlib import Path
import hashlib
import json
import shutil
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'integration/interactive_atac_rna'
ASSETS = ROOT / 'report/interactive_assets'
ASSETS.mkdir(parents=True, exist_ok=True)
combined = pd.read_csv(ROOT / 'integration/feature_scatters/feature_gene_effects.csv', dtype={'entrez_id': str}).assign(dataset='Combined')
separate = pd.read_csv(ROOT / 'integration/batch_RNA_overlaps/all_gene_results.csv', dtype={'entrez_id': str}).rename(columns={'batch':'dataset'})
data = pd.concat([combined, separate], ignore_index=True)
assert not data.duplicated(['dataset','feature','entrez_id']).any()
panels = {}
for (dataset,feature), frame in data.groupby(['dataset','feature'], sort=False):
    records = [[r.entrez_id, r.atac_median_logFC, r.rna_logFC, r.minimum_peak_FDR, r.rna_FDR, int(r.tested_peaks)] for r in frame.itertuples()]
    panels[dataset + ':' + feature] = {'rows':records,'rho':frame.atac_median_logFC.corr(frame.rna_logFC,method='spearman')}
genes = {}
for gene_id, group in data.groupby('entrez_id', sort=False):
    names = group.rna_symbol.dropna()
    ensembl = group.rna_ensembl_id.dropna().unique()
    assert len(ensembl)==1
    genes[gene_id] = {'symbol':str(names.iloc[0]) if len(names) else gene_id, 'ensembl':str(ensembl[0])}
models = pd.read_csv(ROOT / 'integration/gene_atlas/data/gene_models.csv',dtype={'gene_id':str}).set_index('gene_id')
for gene_id, record in genes.items():
    g = models.loc[gene_id]
    record.update(chrom=g.chrom, start=int(g.start)-1, end=int(g.end), strand=g.strand)
atlas = pd.read_csv(ROOT / 'integration/gene_atlas/data/selected_genes.csv',dtype={'entrez_id':str}).set_index('entrez_id')
windows = pd.read_csv(ROOT / 'integration/gene_atlas/data/plot_windows.csv',dtype={'entrez_id':str})
for gene_id in genes:
    genes[gene_id]['atlas'] = gene_id in atlas.index
    if gene_id in atlas.index:
        genes[gene_id]['diagrams'] = [{'id':r.window_id, 'name':r.view} for r in windows[windows.entrez_id==gene_id].itertuples()]
stats = {}
paths = {'Combined':ROOT/'config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv',
         **{b:ROOT/f'integration/batch_RNA_overlaps/sources/ht1080_{b}_LOY_v_CTRL_DESeq2.csv' for b in ['AUG_2025','MAY_2026']}}
for dataset,path in paths.items():
    r = pd.read_csv(path,dtype={'entrezid':str})
    unique = r[r.entrezid.notna() & ~r.entrezid.duplicated(keep=False)]
    stats[dataset] = {row.entrezid:[row.logFC if np.isfinite(row.logFC) else None,row.FDR if np.isfinite(row.FDR) else None] for row in unique.itertuples() if row.entrezid in genes}
meta = pd.read_csv(OUT/'RNA_sample_metadata.csv')
meta = meta[meta.treatment.isin(['CTRL','LOY'])].sort_values(['batch','treatment','sample_id'])
sample_names = meta.colname.tolist()
norm = pd.read_csv(OUT/'RNA_normalized_counts.csv').set_index('gene_id_full')
raw = pd.read_csv(OUT/'RNA_raw_counts.csv').set_index('gene_id')
source = pd.read_csv(paths['Combined'],dtype={'entrezid':str}).set_index('gene_id')
expr = {}
for gene_id, gene in genes.items():
    if gene['ensembl'] not in source.index:
        continue
    full = source.loc[gene['ensembl'],'gene_id_full']
    if full in norm.index:
        expr[gene_id] = {'norm':[round(float(x),6) for x in norm.loc[full,sample_names]],
                         'raw':[int(x) for x in raw.loc[full,sample_names]]}
payload = {'genes':genes,'panels':panels,'rna':stats,'expression':expr,
           'samples':meta[['colname','sample_id','batch','treatment','passage']].to_dict('records'),
           'defaults':{'dataset':'Combined','feature':'promoter','atacFDR':.1,'rnaFDR':.05,'gene':'25850'}}
(ASSETS/'data.js').write_text('window.ATAC_RNA_DATA='+json.dumps(payload,separators=(',',':'),allow_nan=False)+';',encoding='utf-8')
plotly = Path('C:/Users/micha/AppData/Local/R/win-library/4.5/plotly/htmlwidgets/lib/plotlyjs/plotly-latest.min.js')
shutil.copy2(plotly,ASSETS/'plotly.min.js')
summary = {'genes':len(genes),'panels':{name:len(v['rows']) for name,v in panels.items()},'genes_with_expression':len(expr),
           'plotly_source':str(plotly),'plotly_sha256':hashlib.sha256(plotly.read_bytes()).hexdigest()}
(OUT/'data_manifest.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
