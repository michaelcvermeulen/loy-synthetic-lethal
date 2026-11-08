from pathlib import Path
import json
import re
import numpy as np

root = Path(__file__).parents[2]
source = root / "analysis/spectral_competition_20260901/figures/ratio_ladder.svg"
target = root / "analysis/spectral_competition_20260901/figures/ratio_ladder_publication.svg"
svg = source.read_text(encoding="utf-8")

# Drop the source figure's legend group; it is redrawn below as a single panel
# that also carries the correlation statistics.
svg, removed = re.subn(r'\n?<g data-legend="ratio-ladder">.*?</g>', "", svg, flags=re.DOTALL)
if removed != 1:
    raise SystemExit("source legend group not found; regenerate ratio_ladder.svg first")

svg = svg.replace('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 760 620"', '<svg data-publication-ready="true" xmlns="http://www.w3.org/2000/svg" width="608" height="496" viewBox="0 0 608 496"')
svg = svg.replace('role="img">', 'role="img"><g transform="scale(0.8)">')
svg = svg.replace('font-family="sans-serif"', 'font-family="Arial, Helvetica, sans-serif"')
svg = svg.replace('font-size="11"', 'font-size="12"')
svg = svg.replace('stroke="#e5e5e5"', 'stroke="#dddddd"')
svg = svg.replace('stroke="#777"/>', 'stroke="#555" stroke-width="1.2"/>')
svg = svg.replace('September 1 known-ratio samples', 'Known-ratio calibration samples')
svg = svg.replace('Each panel uses identical transformed reporter axes and the same LDA boundary.', 'Publication export: identical transformed reporter axes and shared LDA boundary.')
svg = svg.replace('<text x="380" y="30" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="18">Known-ratio recovery by dye orientation</text>', '')

rows = json.loads((root / "analysis/spectral_competition_20260901/results.json").read_text(encoding="utf-8"))["ratio_ladder"]
def stats(values):
    x = np.asarray([row["nominal_loy"] for row in values], dtype=float)
    y = np.asarray([row["observed_loy"] for row in values], dtype=float)
    r = float(np.corrcoef(x, y)[0, 1])
    return r, r * r, float(np.sqrt(np.mean((y - x) ** 2)))
overall = stats(rows)
ctrl = stats([row for row in rows if row["orientation"] == "CTRL_GFP"])
loyg = stats([row for row in rows if row["orientation"] == "LOY_GFP"])

# Legend + correlation statistics, boxed and placed in the empty upper-left
# corner of the axes (data hug the 1:1 diagonal, so this region stays clear).
FONT = 'font-family="Arial, Helvetica, sans-serif"'
box_x, box_y, box_w, box_h = 97, 72, 252, 139
marker_x, label_x, text_x = 115, 130, 111
annotation = (
    f'<g data-correlation-statistics="ratio-ladder" data-legend="ratio-ladder-publication">'
    f'<rect x="{box_x}" y="{box_y}" width="{box_w}" height="{box_h}" rx="4" fill="#ffffff" fill-opacity="0.92" stroke="#555" stroke-width="1"/>'
    f'<circle cx="{marker_x}" cy="96" r="5" fill="#247a3b"/>'
    f'<text x="{label_x}" y="100" {FONT} font-size="12">LOY-GFP / CTRL-mCherry</text>'
    f'<polygon points="{marker_x},112 {marker_x - 6},123 {marker_x + 6},123" fill="#b33a3a"/>'
    f'<text x="{label_x}" y="122" {FONT} font-size="12">LOY-mCherry / CTRL-GFP</text>'
    f'<line x1="{box_x + 8}" y1="133" x2="{box_x + box_w - 8}" y2="133" stroke="#cccccc" stroke-width="1"/>'
    f'<text x="{text_x}" y="150" {FONT} font-size="11" fill="#333">Correlation with nominal LOY</text>'
    f'<text x="{text_x}" y="167" {FONT} font-size="10" fill="#555">Overall: r={overall[0]:.3f}, R²={overall[1]:.3f}, RMSE={100*overall[2]:.2f} pp</text>'
    f'<text x="{text_x}" y="183" {FONT} font-size="10" fill="#555">CTRL-GFP: r={ctrl[0]:.3f}, R²={ctrl[1]:.3f}, RMSE={100*ctrl[2]:.2f} pp</text>'
    f'<text x="{text_x}" y="199" {FONT} font-size="10" fill="#555">LOY-GFP: r={loyg[0]:.3f}, R²={loyg[1]:.3f}, RMSE={100*loyg[2]:.2f} pp</text>'
    f'</g>'
)
svg = svg.replace('</svg>', annotation + '</g></svg>')
target.write_text(svg, encoding="utf-8")
print(target)
