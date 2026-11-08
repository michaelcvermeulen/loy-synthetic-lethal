"""Verify the rendered STRING draft against frozen inputs and vector exports."""
from pathlib import Path
import sys,json,re,hashlib,xml.etree.ElementTree as ET
THESIS=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(THESIS.parent/'.figure_python_deps'))
import pandas as pd
import pymupdf
from PIL import Image
OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/ht1080_all_sl_string'
STEM='ch04_main_ht1080_all_sl_string_draft'
nodes=pd.read_csv(DATA/'all_nodes.csv')
source=pd.read_csv(THESIS.parent/'crispr_screen/ht1080_loy_screen_2/results/drugZ/LOY_vs_CTRL_unfiltered.txt',sep='\t')
eligible=source.loc[(source.normZ<0)&(source.fdr_synth<.25)].set_index('GENE')
assert set(nodes.GENE)==set(eligible.index) and len(nodes)==48
for r in nodes.itertuples():
    assert r.normZ==eligible.loc[r.GENE,'normZ']
    assert r.fdr_synth==eligible.loc[r.GENE,'fdr_synth']
fig=pymupdf.open(OUT/f'{STEM}.pdf');assert len(fig)==1
page=fig[0];text=page.get_text()
for gene in nodes.GENE:
    assert len(re.findall(r'(?<![A-Za-z0-9])'+re.escape(gene)+r'(?![A-Za-z0-9])',text))==1,gene
fonts=sorted({x[3] for x in page.get_fonts()});assert all('OpenSans' in x for x in fonts),fonts
spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
outside=[s['text'] for s in spans if not page.rect.contains(pymupdf.Rect(s['bbox']))]
assert not outside,outside
overlaps=[]
for i,a in enumerate(spans):
    for b in spans[i+1:]:
        rect=pymupdf.Rect(a['bbox'])&pymupdf.Rect(b['bbox'])
        if not rect.is_empty and rect.width>.5 and rect.height>.5: overlaps.append([a['text'],b['text']])
assert not overlaps,overlaps
svg=ET.parse(OUT/f'{STEM}.svg')
text_nodes=svg.findall('.//{http://www.w3.org/2000/svg}text');assert len(text_nodes)>=48
png=Image.open(OUT/f'{STEM}.png');assert png.size==(5760,3810)
assert abs(png.info['dpi'][0]-600)<.1
preview=OUT/'review'/f'{STEM}_preview.png'
Image.open(preview).convert('L').save(OUT/'review'/f'{STEM}_greyscale.png')
# Local links in the new notes must resolve.
notes=OUT/'ch04_ht1080_all_sl_string_notes.md'
links=re.findall(r'\]\(([^)]+)\)',notes.read_text(encoding='utf-8'))
missing=[s for s in links if not s.startswith('https:') and not (notes.parent/s).exists()]
assert not missing,missing
manifest=json.loads((DATA/'render_manifest.json').read_text())
for name,digest in manifest['input_sha256'].items():
    assert hashlib.sha256((DATA/name).read_bytes()).hexdigest()==digest
result={'original_drugz_gene_set_and_values_match':True,'all_48_genes_labelled_once':True,
    'pdf_pages':len(fig),'pdf_font_names':fonts,'text_outside_page':outside,
    'overlapping_pdf_text_spans':overlaps,'svg_text_elements':len(text_nodes),
    'png_pixels':png.size,'png_dpi':png.info['dpi'],'notes_links_verified':len(links),
    'frozen_network_data_unchanged':True}
(DATA/'export_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
