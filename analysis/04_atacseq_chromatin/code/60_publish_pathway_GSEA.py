"""Publish the verified pathway report after final figure layout corrections."""
from pathlib import Path
import hashlib
import importlib.util
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'integration/pathway_GSEA'
verification=json.loads((OUT/'verification.json').read_text())
assert verification['status']=='passed'
completion=json.loads((OUT/'completion.json').read_text())
source_hash=completion['source_checksums']['integration\\pathway_GSEA\\fgsea_all_results.csv']
assert source_hash==hashlib.sha256((OUT/'fgsea_all_results.csv').read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('plots',ROOT/'code/58_plot_pathway_GSEA.py')
plots=importlib.util.module_from_spec(spec);spec.loader.exec_module(plots)
results=pd.read_csv(OUT/'fgsea_all_results.csv');pairs=pd.read_csv(OUT/'pathway_concordance.csv')
for collection in ('Hallmark','Reactome'):plots.scatter(results,collection,pairs)
for suffix,column in [('NES','NES'),('FDR','padj')]:
    for stem,collection in [('hallmark_pathway_heatmap','Hallmark'),('reactome_pathway_heatmap','Reactome')]:
        matrix=pd.read_csv(OUT/f'{stem}_{suffix}.csv',index_col=0)
        expected=results[results.collection.eq(collection)].pivot(index='pathway',columns='analysis',values=column).reindex(index=matrix.index,columns=matrix.columns)
        pd.testing.assert_frame_equal(matrix,expected,check_names=False,check_exact=False,rtol=1e-12)
main=ROOT/'report/ATACseq_report.html';html=main.read_text(encoding='utf-8')
section='<section id="pathway-gsea"><h2>ATAC and RNA pathway enrichment</h2><p>Hallmark and Reactome GSEA on 11,891 common promoter-linked genes. <a href="pathway_GSEA.html">Pathway heatmaps, ATAC–RNA enrichment comparisons, shared leading-edge genes and sensitivity analyses</a>.</p></section>'
if 'id="pathway-gsea"' not in html:main.write_text(html.replace('</body>',section+'\n</body>'),encoding='utf-8')
completion['status']='complete';completion['verification']=verification
completion['figure_layout_review']='Reviewed PNG exports; corrected Reactome label spacing and scatter legends/labels.'
completion['source_checksums']['code\\58_plot_pathway_GSEA.py']=hashlib.sha256((ROOT/'code/58_plot_pathway_GSEA.py').read_bytes()).hexdigest()
completion['figure_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((OUT/'figures').glob('*'))}
for filename in ('completion.json','ANALYSIS_COMPLETE.json'):(OUT/filename).write_text(json.dumps(completion,indent=2),encoding='utf-8')
print('Published verified pathway report: '+str(ROOT/'report/pathway_GSEA.html'))
