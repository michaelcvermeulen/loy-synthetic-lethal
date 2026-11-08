"""Import the original grant Fig. 5a gel without changing its image pixels."""
from pathlib import Path
import base64
import hashlib
import io
import json
import sys
import zipfile

from PIL import Image

FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
PROJECTS = THESIS.parent
sys.path.insert(0, str(PROJECTS / '.figure_python_deps'))
import pymupdf

STEM = 'ch04_supp_ht1080_ko_genomic_pcr'
SOURCE = FIG / 'source_data' / 'ht1080_ko_genomic_pcr_20261022'
DOCX = PROJECTS / 'grant/CIHR_2026/LOY project grant full app_round2.docx'
MEDIA = 'word/media/image5.tiff'
EXPECTED = '089639ab30f220fe03ee81b8044692d6e2a2a8c4005449e46076b4bd57eea74d'
# Whole panel a, excluding the outside border and adjacent grant panels b/c.
# Original gene labels, lane order, ladder and panel letter are retained.
BOX = (1, 1, 312, 279)


def main():
    SOURCE.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(DOCX) as archive:
        original = archive.read(MEDIA)
    assert hashlib.sha256(original).hexdigest() == EXPECTED
    (SOURCE / 'grant_figure5_original.tiff').write_bytes(original)
    full = Image.open(io.BytesIO(original))
    panel = full.crop(BOX)
    png = FIG / f'{STEM}.png'
    panel.save(png)
    assert Image.open(png).tobytes() == panel.tobytes()

    # These are lossless containers for the same raster, not new image data.
    width_pt = 3.36 * 72
    height_pt = width_pt * panel.height / panel.width
    pdf = pymupdf.open()
    page = pdf.new_page(width=width_pt, height=height_pt)
    page.insert_image(page.rect, filename=str(png))
    pdf.save(FIG / f'{STEM}.pdf', deflate=True)
    pdf.close()
    payload = base64.b64encode(png.read_bytes()).decode('ascii')
    (FIG / f'{STEM}.svg').write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="3.36in" height="{height_pt / 72:.8f}in" viewBox="0 0 {panel.width} {panel.height}">'
        '<title>Genomic PCR of DDX3X, TADA2B and TAF6L edited HT1080 populations</title>'
        f'<image width="{panel.width}" height="{panel.height}" xlink:href="data:image/png;base64,{payload}"/>'
        '</svg>\n', encoding='utf-8')
    provenance = {
        'source_docx': str(DOCX), 'source_media': MEDIA,
        'source_sha256': EXPECTED, 'source_figure': '5a', 'grant_pdf_page': 7,
        'source_size_px': list(full.size), 'crop_box_ltrb_px': list(BOX),
        'output_size_px': list(panel.size), 'native_pixels_preserved': True,
        'operations': ['rectangular extraction of original panel a', 'lossless PNG/PDF/SVG containers'],
        'no_adjustment': ['brightness', 'contrast', 'bands', 'lane order', 'labels', 'resampling'],
        'genes_displayed': ['DDX3X', 'TADA2B', 'TAF6L'],
        'missing_image': 'SREBF1, to be supplied by author',
        'assay': 'Genomic PCR; not a western blot or confirmation of protein loss',
        'limitations': 'Grant-resolution image; no size annotations or raw uncropped acquisition supplied.',
        'exports': {ext: hashlib.sha256((FIG / f'{STEM}.{ext}').read_bytes()).hexdigest()
                    for ext in ('png', 'pdf', 'svg')},
    }
    (SOURCE / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'figure': str(png), 'size_px': list(panel.size), 'native_pixels_preserved': True}))


if __name__ == '__main__':
    main()
