"""Shared paths, style and export for the DU145 RNA-seq thesis figures (Results 3.4.6.3).

All numerical inputs are read from the validated DU145 reanalysis
(RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE) and its top-gene follow-up.
No statistic is recomputed here; plotted tables are written to source_data.
"""
from pathlib import Path
import sys, json, hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
PROJECTS = THESIS.parent
ANALYSIS = PROJECTS / 'RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE'
T = ANALYSIS / 'outputs/tables'
FOLLOW = ANALYSIS / 'followup_top_genes'
LPWGS = FIG / 'source_data/du145_lpwgs_confirmation'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style  # noqa: E402

apply_style(8)
BLUE, RED, GREY, DARK = '#3C5488', '#E64B35', '#C8C8C8', '#555555'
COND = {'CTRL': BLUE, 'LOY': RED}
COND_LABEL = {'CTRL': 'sgCTRL', 'LOY': 'Y-depleted'}
PASSAGE_MARKER = {'P12': 'o', 'P13': 's', 'P14': '^'}
# Genomic-interval categories (not conditions); NPG subset chosen to stay distinct from blue/vermilion.
INTERVAL = {'PAR1': '#4DBBD5', 'retained distal Yp': '#00A087', 'depleted MSY': '#7E6148', 'PAR2': '#8491B4'}
MM = 1 / 25.4
# Gene-set themes assigned by leading-edge composition (followup_top_genes/gene_set_themes/classify_themes.R).
THEMES = FOLLOW / 'gene_set_themes/significant_sets_by_theme.csv'
THEME_ORDER = ['ribosomal/translation', 'interferon_inflammatory', 'GRHL2_epithelial', 'mitotic', 'other',
               'histone cluster (6p22)', 'positional']
THEME_LABEL = {'ribosomal/translation': 'Ribosomal', 'interferon_inflammatory': 'Interferon',
               'GRHL2_epithelial': 'GRHL2 epithelial', 'mitotic': 'G2/M', 'other': 'Other',
               'histone cluster (6p22)': 'Histone (6p22)', 'positional': 'Positional'}
STATUS_ORDER = ['survives both screens', 'positionally concentrated', 'copy-number sensitive']
STATUS_COLOUR = {'survives both screens': RED, 'positionally concentrated': '#C98A63', 'copy-number sensitive': '#F0B27A'}
STATUS_LABEL = {'survives both screens': 'Survives both screens', 'positionally concentrated': 'Positionally concentrated',
                'copy-number sensitive': 'Copy-number sensitive'}


def short(sample_key):
    """LNCaP-sgCTRL-P12 -> 'sgCTRL P12'; source labels are preserved in the tables."""
    arm, passage = sample_key.split('-')[1], sample_key.split('-')[2]
    return f"{'sgCTRL' if arm == 'sgCTRL' else 'Y-dep.'} {passage}"


def read(name):
    return pd.read_csv(T / f'{name}.csv')


def panel(fig, ax, letter, title, dx=-0.07, dy=0.03, y=None):
    pos = ax.get_position()
    y = pos.y1 + dy if y is None else y
    fig.text(pos.x0 + dx, y, letter, fontsize=12, weight='semibold', va='bottom')
    if title:
        fig.text(pos.x0, y, title, fontsize=9, weight='semibold', va='bottom')


def fmt_p(p, label='adj. P'):
    if p < 1e-3:
        m, e = f'{p:.1e}'.split('e')
        return f'{label} = {m}×10$^{{{int(e)}}}$'
    return f'{label} = {p:.2g}' if p < 0.1 else f'{label} = {p:.2f}'


def export(fig, stem, tables, notes):
    """Stage then replace PDF/SVG/PNG; write plotted tables and provenance."""
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
    (src / 'provenance.json').write_text(json.dumps({
        'figure': stem, 'font': 'Open Sans', 'analysis_root': str(ANALYSIS),
        'inputs': [{'path': p, 'sha256': hashlib.sha256(Path(p).read_bytes()).hexdigest()} for p in inputs],
        **notes, 'exports': records}, indent=2), encoding='utf-8')
    print(stem, flush=True)
