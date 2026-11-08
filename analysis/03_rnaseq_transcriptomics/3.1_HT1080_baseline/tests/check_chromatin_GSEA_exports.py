"""Check GSEA plot coordinates, original annotations and single-panel publication exports."""
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlparse
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from PIL import Image
from pypdf import PdfReader

root = Path(__file__).resolve().parents[1]
read = lambda p: pd.read_csv(root / p)
stats = read('outputs/tables/GSEA_CURVES_selected_signatures.csv')
exports = read('outputs/figures/chromatin_GSEA_individual/manifest.csv')
original = read('outputs/tables/fgsea_all_collections.csv')
assert len(stats) == 12 and len(exports) == 12
assert stats.groupby('figure').panel.nunique().eq(6).all()
assert stats.ES_agreement_error.lt(1e-10).all() and stats.leading_edge_exact_agreement.all()
rank = read('outputs/tables/GSEA_CURVES_frozen_ranking.csv')
assert len(rank) == 25461 and rank.gene_id.nunique() == len(rank)
font_checks = []
for row in exports.itertuples():
    saved = stats[stats.pathway.eq(row.pathway)].iloc[0]
    a = original[original.pathway.eq(row.pathway) & original.contrast.eq('combined')].iloc[0]
    assert abs(saved.NES - a.NES) < 1e-12
    assert abs(saved.original_global_BH - a.BH_across_collections) < 1e-12
    curve = read(row.source_data)
    assert len(curve) == len(rank) + 1
    assert curve['rank'].tolist() == list(range(25462))
    assert curve.iloc[1:].gene_id.tolist() == rank.gene_id.tolist()
    assert np.allclose(curve.iloc[1:].ranking_metric, rank.combined_log2_difference)
    assert curve.hit.sum() == saved.genes and curve.leading_edge.sum() == saved.leading_genes
    assert set(curve[curve.leading_edge].gene_name) == set(a.leadingEdge.split(';'))
    assert abs(curve.iloc[0].running_ES) < 1e-12 and abs(curve.iloc[-1].running_ES) < 1e-10
    assert abs(curve.iloc[int(saved.peak_rank)].running_ES - saved.ES) < 1e-10
    assert curve.running_ES.between(-0.70, 0.50).all(), 'Shared plotting scale clips a curve'
    assert np.allclose(curve.running_ES, curve.increment.cumsum(), atol=1e-10)
    pdf = PdfReader(root / row.pdf)
    assert len(pdf.pages) == 1
    page = pdf.pages[0]
    text = ' '.join(page.extract_text().split())
    for phrase in [saved.title, 'Running ES', 'global BH', 'Hits', 'Gene rank']:
        assert phrase in text, (row.pathway, phrase)
    assert abs(float(page.mediabox.width)*25.4/72 - 89) < 25.4/72
    assert abs(float(page.mediabox.height)*25.4/72 - 78) < 25.4/72
    fonts = []
    for ref in page['/Resources']['/Font'].values():
        font = ref.get_object()
        for item in font.get('/DescendantFonts', [font]):
            f = item.get_object(); desc = f['/FontDescriptor'].get_object()
            assert any(k in desc for k in ['/FontFile', '/FontFile2', '/FontFile3'])
            assert 'Arial' in str(f.get('/BaseFont'))
            fonts.append(str(f.get('/BaseFont')))
    svg = (root / row.svg).read_text(encoding='utf-8')
    assert '<text' in svg and 'Arial' in svg
    sizes = [float(s) for s in re.findall(r'font-size:\s*([0-9.]+)px',svg)]
    assert min(sizes) >= 5
    with Image.open(root / row.png) as image:
        assert abs(image.width - 89/25.4*600) < 1.1 and abs(image.height - 78/25.4*600) < 1.1
        assert all(abs(d-600) < 1 for d in image.info['dpi'])
    for kind in ['pdf','svg','png']:
        assert hashlib.sha256((root / getattr(row,kind)).read_bytes()).hexdigest() == getattr(row,kind+'_sha256')
    font_checks.append(dict(pathway=row.pathway,embedded_Arial_fonts=fonts))
assert stats[stats.pathway.eq('PRC2_EED_DN.V1_DN')].original_global_BH.eq(1).all()
assert stats[stats.pathway.eq('PRC2_SUZ12_UP.V1_UP')].original_global_BH.gt(.05).all()
for file in ['provenance/GSEA_CURVES_frozen_inputs.csv']:
    for row in read(file).itertuples():
        assert hashlib.sha256((root / row.path).read_bytes()).hexdigest() == row.sha256
prior = read('provenance/GSEA_CURVES_prior_figure_hashes.csv')
assert len(prior) == 17
for row in prior.to_dict('records'):
    for kind in ['pdf','svg','png']:
        assert hashlib.sha256((root / row[kind]).read_bytes()).hexdigest() == row[kind+'_sha256']
report = root / 'outputs/report/chromatin_GSEA_curves.html'
soup = BeautifulSoup(report.read_text(encoding='utf-8'),'html.parser')
images = soup.find_all('img')
assert len(images) == 2 and all(i.get('src','').startswith('data:') for i in images)
assert len(soup.find_all('table')) == 1
for a in soup.find_all('a',href=True):
    url = a['href']
    if url.startswith(('#','http:','https:','data:','javascript:','mailto:')):
        continue
    path = unquote(urlparse(url).path)
    assert not path or (report.parent / path).exists(), url
result = dict(curves=12,curve_rows=25462,exact_ES_and_leading_edges=True,original_global_BH_preserved=True,
              individual_exports=36,composite_figures=2,composite_source_panels=12,
              prior_figure_exports_unchanged=51,report_figures=2,broken_local_links=0,fonts=font_checks)
(root / 'provenance/GSEA_CURVES_export_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('Passed 12 exact GSEA curves, 36 individual exports, original statistics, report links and 51 unchanged prior figure exports.')
