"""Shared paths, style and export for the cross-line RNA-seq figures (Results 3.4.6.4).

Reads only the validated cross-line analysis (RNAseq/analysis/crossline_SEPT_2026). No statistic is
recomputed here; plotted tables are written to source_data with hashed provenance.
"""
from pathlib import Path
import json, hashlib
import matplotlib.pyplot as plt
import pandas as pd
from du145_rna_common import FIG, PROJECTS, BLUE, RED, GREY, DARK, MM, panel  # noqa: F401  (applies thesis style)

CROSS = PROJECTS / 'RNAseq/analysis/crossline_SEPT_2026'
CT = CROSS / 'outputs/tables'
LINES = ['HT1080', 'A549', 'DU145']
# Cell-line identity colours: an NPG subset kept distinct from the sgCTRL/Y-depleted condition colours.
LINE_COLOUR = {'HT1080': '#00A087', 'A549': '#7E6148', 'DU145': '#8491B4'}
THEME_LABEL = {'interferon_inflammatory': 'Interferon / inflammatory', 'GRHL2_epithelial': 'GRHL2 epithelial',
               'mitotic': 'G2/M', 'ribosomal/translation': 'Ribosomal', 'polycomb_targets': 'PRC2 / H3K27me3 targets',
               'histone cluster (6p22)': 'Histone (6p22)', 'positional': 'Positional', 'other': 'Other'}
GROUP_LABEL = {'in_house': 'Thesis lines', 'A_somatic_LOY': 'Somatic LOY', 'A2_single_cell_mosaic': 'Natural Y-mosaic,\nsingle-cell', 'B_45X_vs_46XY': '45,X vs 46,XY',
               'C_not_LOY': 'Not LOY'}


def ct(name):
    return pd.read_csv(CT / f'{name}.csv')


def export(fig, stem, tables, notes):
    """Stage then replace PDF/SVG/PNG; write plotted tables and hashed provenance."""
    src = FIG / 'source_data' / stem.removeprefix('ch03_').replace('main_', '').replace('supp_', '')
    src.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(src / f'{name}.csv', index=False)
    records = []
    for ext in ['pdf', 'svg', 'png']:
        path = FIG / f'{stem}.{ext}'
        stage = FIG / f'{stem}.staging.{ext}'
        fig.savefig(stage, dpi=600)
        stage.replace(path)
        records.append({'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    (FIG / 'previews').mkdir(exist_ok=True)
    fig.savefig(FIG / 'previews' / f'{stem}_preview.png', dpi=160)
    plt.close(fig)
    inputs = sorted({str(p) for p in notes.pop('inputs', [])})
    validation = json.loads((CROSS / 'provenance/validation.json').read_text())
    (src / 'provenance.json').write_text(json.dumps({
        'figure': stem, 'font': 'Open Sans', 'analysis_root': str(CROSS), 'analysis_validation': validation,
        'inputs': [{'path': p, 'sha256': hashlib.sha256(Path(p).read_bytes()).hexdigest()} for p in inputs],
        **notes, 'exports': records}, indent=2), encoding='utf-8')
    print(stem, flush=True)
