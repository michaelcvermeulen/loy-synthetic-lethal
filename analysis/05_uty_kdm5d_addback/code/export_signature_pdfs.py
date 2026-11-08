"""Print editable SVG plots to PDF with explicitly embedded Open Sans fonts."""
from pathlib import Path
import base64
import re
import subprocess
import tempfile
import sys

root = Path(__file__).resolve().parents[1]
fonts = Path('A:/Dropbox/LOY_paralogs/projects/SEPT_11_2026_THESIS/skills/thesis-figure-style/assets/fonts')
browser = Path('C:/Program Files/Google/Chrome/Application/chrome.exe')
css = ''.join(
    '@font-face{font-family:"Open Sans";font-weight:' + weight +
    ';src:url(data:font/ttf;base64,' + base64.b64encode((fonts / filename).read_bytes()).decode() + ') format("truetype");}'
    for weight, filename in [('400', 'OpenSans-Regular.ttf'), ('600 900', 'OpenSans-SemiBold.ttf')]
)
stems = sys.argv[1:] or ['ch03_main_ht1080_addback_baseline_reversal',
         'student_panels/panel_f_frozen_signature_scores',
         'student_panels/panel_g_signature_score_contrasts',
         'student_panels/panel_h_signature_gene_heatmap',
         'student_panels/panel_i_directional_fry_tests']
for stem in stems:
    svg_path = root / 'figures' / (stem + '.svg')
    svg = svg_path.read_text(encoding='utf-8')
    width, height = re.search(r"viewBox=['\"]0 0 ([\d.]+) ([\d.]+)", svg).groups()
    svg = re.sub(r'<\?xml.*?\?>', '', svg)
    svg = re.sub(r'<svg\b', '<svg style="width:100%;height:100%;display:block"', svg, count=1)
    with tempfile.TemporaryDirectory(prefix='signature_pdf_') as tmp:
        html = Path(tmp) / 'figure.html'
        html.write_text('<!doctype html><meta charset="utf-8"><style>' + css +
            f'@page{{size:{width}pt {height}pt;margin:0}}html,body{{margin:0;width:{width}pt;height:{height}pt}}' +
            '</style>' + svg, encoding='utf-8')
        output = svg_path.with_suffix('.pdf')
        subprocess.run([str(browser), '--headless', '--no-pdf-header-footer', '--disable-gpu',
                        '--no-first-run', '--no-default-browser-check', '--disable-extensions',
                        '--disable-background-networking', '--virtual-time-budget=2500',
                        '--user-data-dir=' + str(Path(tmp) / 'profile'),
                        '--print-to-pdf=' + str(output), html.as_uri()],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0), timeout=60)
        data = output.read_bytes()
        names = re.findall(rb'/FontName\s*/([^\s/]+)', data)
        assert names and all(b'OpenSans' in n for n in names), (stem, names)
        assert b'/FontFile2' in data, stem
        print(stem, ': Open Sans embedded; vector text preserved.')
