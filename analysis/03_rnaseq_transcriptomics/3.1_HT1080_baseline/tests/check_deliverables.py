"""Inspect actual PDF fonts, SVG text, PNG physical size and local report links."""
import csv
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlparse

from bs4 import BeautifulSoup
from PIL import Image
from pypdf import PdfReader

root = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((root / 'outputs/figures/figure_manifest.csv').open(encoding='utf-8')))
figures = {row['figure']: row for row in rows}
checks = []

for name, row in figures.items():
    reader = PdfReader(root / row['pdf'])
    assert len(reader.pages) == 1, name
    page = reader.pages[0]
    width_mm = float(page.mediabox.width) * 25.4 / 72
    height_mm = float(page.mediabox.height) * 25.4 / 72
    # R's Cairo device quantizes its canvas to integer PostScript points (<0.353 mm).
    assert abs(width_mm - float(row['width_mm'])) < 25.4 / 72, name
    assert abs(height_mm - float(row['height_mm'])) < 25.4 / 72, name
    assert page.extract_text().strip(), f'No editable PDF text: {name}'
    fonts = []
    for ref in page['/Resources']['/Font'].values():
        font = ref.get_object()
        descendants = font.get('/DescendantFonts', [font])
        for ref2 in descendants:
            f = ref2.get_object()
            desc = f.get('/FontDescriptor')
            assert desc is not None, f'No font descriptor: {name}'
            desc = desc.get_object()
            embedded = any(key in desc for key in ('/FontFile', '/FontFile2', '/FontFile3'))
            assert embedded, f'Unembedded font in {name}'
            fonts.append(str(f.get('/BaseFont')))
    assert all('Arial' in f for f in fonts), (name, fonts)
    svg = (root / row['svg']).read_text(encoding='utf-8')
    assert '<text' in svg and 'Arial' in svg, f'Missing editable SVG text: {name}'
    sizes = [float(s) for s in re.findall(r'font-size:\s*([0-9.]+)px', svg)]
    # svglite uses CSS pixel coordinates at 72 coordinate units per inch; inspect relative point-scale values.
    assert not sizes or min(sizes) >= 5, (name, min(sizes))
    with Image.open(root / row['png']) as im:
        dpi = im.info.get('dpi')
        assert dpi and all(abs(v - 600) < 1 for v in dpi), (name, dpi)
        assert abs(im.width - float(row['width_mm']) / 25.4 * 600) <= 1.1
        assert abs(im.height - float(row['height_mm']) / 25.4 * 600) <= 1.1
        pixels = im.size
    caption = (root / row['caption']).read_text(encoding='utf-8')
    assert len(caption.split()) < 250 and '{{' not in caption, name
    checks.append(dict(figure=name, width_mm=width_mm, height_mm=height_mm,
                       embedded_fonts=fonts, editable_pdf_text=True, editable_svg_text=True,
                       png_pixels=pixels, png_dpi=dpi, caption_words=len(caption.split())))

for row in rows:
    assert (root / row['source_data']).is_file(), row['source_data']
report = root / 'outputs/report/analysis.html'
soup = BeautifulSoup(report.read_text(encoding='utf-8'), 'html.parser')
broken = []
for a in soup.find_all('a', href=True):
    url = a['href']
    if not url or url.startswith(('#', 'http:', 'https:', 'mailto:', 'data:', 'javascript:')):
        continue
    path = unquote(urlparse(url).path)
    if path and not (report.parent / path).exists():
        broken.append(url)
assert not broken, broken
images = soup.find_all('img')
assert len(images) == len(figures), 'Report figure count differs from the figure manifest'
assert all(im.get('src', '').startswith('data:') for im in images), 'Report images are not self-contained'
assert (root / 'outputs/captions/figure_captions.docx').is_file()
assert (root / 'renv.lock').is_file()
result = dict(figures=checks, source_panels=len(rows), broken_local_links=broken,
              embedded_report_figures=len(images), caption_docx=True, renv_lock=True)
(root / 'provenance/deliverable_checks.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
print(f'Passed: {len(figures)} figures, {len(rows)} source panels, embedded Arial fonts, PDF/SVG text, '
      f'600 dpi PNGs, captions, self-contained report and all local links.')
