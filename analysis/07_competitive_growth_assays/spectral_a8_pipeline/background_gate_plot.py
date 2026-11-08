from pathlib import Path
from html import escape
import numpy as np
from spectral.a8_pipeline import competition_report as c

repo = Path(__file__).parents[2]
d4 = repo / "spectral/20260901_125230_Mike_Run1/Mike_Run1"
g = c.load_cef_gates(d4 / "20260901_125230_Mike_Run1.cef")
files = c._file_map(d4 / "1")
unst = np.vstack([c._primary_features(files[n], g)[0] for n in ("A01-unstained_CTRL", "A02-unstained_LOY")])
thresholds = c.background_gate_thresholds(unst, 99)
names = ["A01-unstained_CTRL", "B05-50-50-CTRL_GFP1", "D01-LOY_mCherry_CTRL_GFP_exp_rep1_001", "E01-LOY_GFP_CTRL_mCherry_exp_rep1_003"]
W,H,CW,CH = 1200, 520, 290, 225
x0,x1,y0,y1 = -3.5,6.5,-8,8
colors={0:'#888888',1:'#247a3b',2:'#b33a63',3:'#d07a00'}
parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"><title>Exploratory unstained background gate</title><desc>Representative transformed reporter clouds with 99th-percentile thresholds. This does not replace the primary LDA analysis.</desc><rect width="100%" height="100%" fill="white"/><text x="20" y="28" font-family="sans-serif" font-size="18">Exploratory 99th-percentile unstained gate</text><text x="20" y="49" font-family="sans-serif" font-size="12">eGFP threshold {thresholds[0]:.3f}; mCherry threshold {thresholds[1]:.3f}. Gray background; green GFP+; magenta mCherry+; orange double+.</text>']
for i,name in enumerate(names):
    f = c._primary_features(files[name], g)[0]; calls = c.apply_background_gate(f, thresholds); left,top=20+i*CW,70; pl,pt=left+35,top+28; pw,ph=235,160
    px=lambda v: pl+(v-x0)/(x1-x0)*pw; py=lambda v: pt+ph-(v-y0)/(y1-y0)*ph
    parts.append(f'<text x="{left+CW/2}" y="{top+13}" text-anchor="middle" font-family="sans-serif" font-size="10">{escape(name)}</text><rect x="{pl}" y="{pt}" width="{pw}" height="{ph}" fill="#fbfbfb" stroke="#777"/><line x1="{px(thresholds[0])}" y1="{pt}" x2="{px(thresholds[0])}" y2="{pt+ph}" stroke="#247a3b" stroke-dasharray="4,3"/><line x1="{pl}" y1="{py(thresholds[1])}" x2="{pl+pw}" y2="{py(thresholds[1])}" stroke="#b33a63" stroke-dasharray="4,3"/>')
    for j in np.linspace(0,len(f)-1,min(800,len(f)),dtype=int):
        parts.append(f'<circle cx="{px(f[j,0]):.2f}" cy="{py(f[j,1]):.2f}" r="1.4" fill="{colors[int(calls[j])]}" fill-opacity=".42"/>')
    parts.append(f'<text x="{pl+pw/2}" y="{pt+ph+18}" text-anchor="middle" font-family="sans-serif" font-size="10">arcsinh(eGFP-A/1000)</text><text x="{pl-13}" y="{pt+ph/2}" transform="rotate(-90 {pl-13} {pt+ph/2})" text-anchor="middle" font-family="sans-serif" font-size="10">arcsinh(mCherry-A/1000)</text>')
parts.append('</svg>')
out=repo/'analysis/spectral_competition_20260901/figures/background_gate_representative.svg'; out.write_text(''.join(parts),encoding='utf-8'); print(out)
