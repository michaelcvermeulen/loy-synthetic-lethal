"""Compose two merged Chapter 4 supplements from existing vector exports (22 Oct 2026).

A549 sensitivity: ch04_supp_a549_p3_gene_sensitivity_review (a-f) stacked above
  ch04_supp_a549_isogenicz_drugz_by_passage_review (relettered g,h).
HT1080 screen QC: ch04_supp_ht1080_crispr_screen_qc (a-h) stacked above a landscape
  redraw of the chromosome-shift heatmap (panel i) from its saved plotted values.

Re-run after either source figure is re-rendered. No values are recomputed.
Usage: python compose_merged_supps.py [a549|qc|all]
"""
from pathlib import Path
import hashlib
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent / '.figure_python_deps'))
OUT = ROOT / 'figures/chapter4'
DATA = OUT / 'source_data/merged_supps_20261022'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd
import pymupdf

sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
SEMIBOLD = ROOT / 'skills/thesis-figure-style/assets/fonts/OpenSans-SemiBold.ttf'
GAP_PT = 6  # vertical gap between stacked sources


def reletter(page, mapping):
    """Replace 13-pt semibold panel letters in place (text only; graphics untouched)."""
    found = []
    for block in page.get_text('dict')['blocks']:
        for line in block.get('lines', []):
            for span in line['spans']:
                t = span['text'].strip()
                if t in mapping and span['size'] > 11 and 'SemiBold' in span['font']:
                    found.append((t, pymupdf.Rect(span['bbox']), span['origin'], span['size']))
    assert sorted(t for t, *_ in found) == sorted(mapping), found
    for _, rect, _, _ in found:
        page.add_redact_annot(rect)
    page.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE, graphics=pymupdf.PDF_REDACT_LINE_ART_NONE)
    for t, _, origin, size in found:
        page.insert_text(origin, mapping[t], fontsize=size, fontname='OpenSansSemiBold', fontfile=str(SEMIBOLD))


def stack(sources, stem, letters=None):
    """sources: list of PDF paths (top to bottom); letters: optional {index: mapping}."""
    docs = [pymupdf.open(p) for p in sources]
    width = max(d[0].rect.width for d in docs)
    assert all(abs(d[0].rect.width - width) < 1 for d in docs), 'sources must share a width'
    height = sum(d[0].rect.height for d in docs) + GAP_PT * (len(docs) - 1)
    out = pymupdf.open()
    page = out.new_page(width=width, height=height)
    y = 0
    for i, d in enumerate(docs):
        src = d
        if letters and i in letters:
            src = pymupdf.open()
            src.insert_pdf(d)
            reletter(src[0], letters[i])
        r = src[0].rect
        page.show_pdf_page(pymupdf.Rect(0, y, width, y + r.height), src, 0)
        y += r.height + GAP_PT
    out.save(OUT / f'{stem}.pdf', garbage=3, deflate=True)
    out[0].get_pixmap(dpi=600).save(OUT / f'{stem}.png')
    (OUT / f'{stem}.svg').write_text(out[0].get_svg_image(text_as_path=False), encoding='utf-8')
    preview = OUT / 'review' / f'{stem}_preview.png'
    out[0].get_pixmap(dpi=170).save(preview)
    hashed = preview.with_name(f'{stem}_{hashlib.sha256(preview.read_bytes()).hexdigest()[:12]}.png')
    preview.replace(hashed)
    text = out[0].get_text()
    return {'stem': stem, 'size_inches': [round(width / 72, 2), round(height / 72, 2)], 'preview': str(hashed),
            'sources': {str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources}}, text


def chromosome_strip(path):
    """Landscape panel i: the chromosome-shift heatmap with chromosomes as columns."""
    v = pd.read_csv(OUT / 'source_data/ht1080_chromosome_shift_heatmap/plotted_values.csv', dtype={'chromosome': str})
    assert len(v) == 25 and v.chromosome.iloc[-1] == 'NTC'
    values = v[['median_delta_relative_NTC', 'median_delta_D0_ge50_relative_NTC']].to_numpy().T
    assert np.isfinite(values).all() and np.abs(values).max() <= .4
    apply_style(8)
    fig = plt.figure(figsize=(9.6, 2.25))
    ax = fig.add_axes([.085, .30, .80, .40])
    cmap = LinearSegmentedColormap.from_list('chromosome_signed_effect', ['#3B6FB6', '#F7F7F7', '#D88036'])
    mesh = ax.pcolormesh(np.arange(26), np.arange(3), values, cmap=cmap, norm=Normalize(-.4, .4),
                         edgecolors='white', linewidth=.8)
    ax.set_xlim(0, 25); ax.set_ylim(2, 0)
    ax.set_xticks(np.arange(25) + .5, v.chromosome.tolist(), fontsize=7.5)
    ax.xaxis.tick_top()
    ax.set_yticks([.5, 1.5], ['All guides', 'D0 ≥50'], fontsize=7.5)
    ax.tick_params(length=0, pad=3)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.axvline(24, color='#777777', lw=.65)
    for r in range(2):
        for c in range(25):
            val = values[r, c]
            label = f'{val:+.3f}'.replace('-', '−') if abs(val) >= .0005 else '0.000'
            ax.text(c + .5, r + .5, label, ha='center', va='center', fontsize=5.6, color='#222222')
    for c, row in v.iterrows():
        ax.text(c + .5, 2.32, f'{row.n_paired_guides:,}\n{row.n_D0_ge50_both:,}', ha='center', va='top',
                fontsize=5.3, color='#666666', linespacing=1.1)
    ax.text(-.25, 2.32, 'n guides', ha='right', va='top', fontsize=6.3, color='#666666')
    fig.text(.012, .97, 'i', fontsize=13, fontweight='semibold', va='top')
    fig.text(.085, .95, 'Median D0-to-endpoint ΔLFC by target chromosome, relative to non-targeting controls; '
             'blue, relative guide depletion in LOY vs CTRL; orange, enrichment', fontsize=7.8, va='top')
    cb = fig.colorbar(mesh, cax=fig.add_axes([.905, .30, .012, .40]), ticks=[-.4, 0, .4])
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6.5, length=2, pad=2)
    fig.savefig(path)
    plt.close(fig)


def main(which):
    records = {}
    if which in ('a549', 'all'):
        rec, text = stack([OUT / 'ch04_supp_a549_p3_gene_sensitivity_review.pdf',
                           OUT / 'ch04_supp_a549_isogenicz_drugz_by_passage_review.pdf'],
                          'ch04_supp_a549_sensitivity_combined', letters={1: {'a': 'g', 'b': 'h'}})
        assert 'Spearman ρ = 0.655' in text and 'Spearman ρ = 0.853' in text
        records['a549'] = rec
    if which in ('qc', 'all'):
        strip = DATA / 'ht1080_chromosome_strip_panel_i.pdf'
        chromosome_strip(strip)
        rec, text = stack([OUT / 'ch04_supp_ht1080_crispr_screen_qc.pdf', strip], 'ch04_supp_ht1080_screen_qc_combined')
        assert '+0.240' in text and '+0.165' in text
        records['qc'] = rec
    prev = json.loads((DATA / 'manifest.json').read_text()) if (DATA / 'manifest.json').exists() else {}
    prev.update(records)
    (DATA / 'manifest.json').write_text(json.dumps(prev, indent=2), encoding='utf-8')
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'all')
