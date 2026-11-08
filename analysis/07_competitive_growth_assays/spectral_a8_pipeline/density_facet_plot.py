from pathlib import Path
from html import escape
import numpy as np
from spectral.a8_pipeline import competition_report as c

repo = Path(__file__).parents[2]
d4 = repo / "spectral/20260901_125230_Mike_Run1/Mike_Run1"
g = c.load_cef_gates(d4 / "20260901_125230_Mike_Run1.cef")
files = c._file_map(d4 / "1")
samples = [(name, c._primary_features(files[name], g)[0]) for name in sorted(files)]
ranges = [(-4.0, 11.0), (-8.0, 9.0)]
bins = 80
histograms = []
maxima = [0.0, 0.0]
for name, x in samples:
    hs = []
    for channel in range(2):
        h, edges = np.histogram(x[:, channel], bins=bins, range=ranges[channel], density=True)
        hs.append((h, edges)); maxima[channel] = max(maxima[channel], float(h.max()))
    histograms.append((name, hs))
cols, cell_w, cell_h = 4, 320, 290
rows = int(np.ceil(len(samples) / cols))
W, H = cols * cell_w + 20, 65 + rows * cell_h
parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}"><title>Per-sample transformed reporter densities</title><desc>Each sample has eGFP-A arcsinh density above mCherry-A arcsinh density. Primary LDA analysis is unchanged.</desc><rect width="100%" height="100%" fill="white"/><text x="20" y="26" font-family="sans-serif" font-size="18">Per-sample transformed reporter densities</text><text x="20" y="46" font-family="sans-serif" font-size="12">eGFP-A arcsinh(value/1000) above mCherry-A arcsinh(value/1000); shared density scales within channel.</text>']
def sx(v, left, width, r): return left + (v-r[0])/(r[1]-r[0])*width
for i, (name, hs) in enumerate(histograms):
    col, row = i % cols, i // cols; left, top = 10 + col*cell_w, 60 + row*cell_h
    parts.append(f'<text x="{left+cell_w/2}" y="{top+14}" text-anchor="middle" font-family="sans-serif" font-size="11">{escape(name)}</text>')
    for channel, (h, edges) in enumerate(hs):
        pl, pt, pw, ph = left+42, top+24+channel*112, 250, 82; r = ranges[channel]
        points = [f'{sx((edges[j]+edges[j+1])/2,pl,pw,r):.2f},{pt+ph-(val/maxima[channel])*ph:.2f}' for j,val in enumerate(h)]
        color = '#247a3b' if channel==0 else '#b33a63'; label = 'eGFP-A arcsinh' if channel==0 else 'mCherry-A arcsinh'
        parts.append(f'<rect x="{pl}" y="{pt}" width="{pw}" height="{ph}" fill="#fbfbfb" stroke="#aaa"/><path d="M {" L ".join(points)}" fill="none" stroke="{color}" stroke-width="1.6"/><text x="{pl-7}" y="{pt+11}" text-anchor="end" font-family="sans-serif" font-size="10">{label}</text><text x="{pl}" y="{pt+ph+14}" font-family="sans-serif" font-size="9">{r[0]:g}</text><text x="{pl+pw}" y="{pt+ph+14}" text-anchor="end" font-family="sans-serif" font-size="9">{r[1]:g}</text>')
out = repo / "analysis/spectral_competition_20260901/figures/all_samples_density_facets.svg"
out.write_text(''.join(parts)+'</svg>', encoding='utf-8'); print(out)
