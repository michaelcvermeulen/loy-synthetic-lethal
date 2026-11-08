"""Working comparison: chrY lpWGS tracks for every sequenced sample, grouped by cell line.

Drawn like panel c of ch03_supp_du145_drug_treatment_comparison
(SEPT_11_2026_THESIS/figures/chapter3/code/make_du145_panobinostat.py).
QDNAseq samples use saved 100 kb T2T bins and saved segments (no refitting).
A549-BC14 never went through QDNAseq: it is shown from its mosdepth 1 kb
MAPQ>=20 track, averaged into the same 100 kb bins and divided by the
autosomal median (no GC correction, no segments).
"""
from pathlib import Path
import sys, gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

WGS = Path(__file__).resolve().parents[1]
THESIS = WGS.parent / 'SEPT_11_2026_THESIS'
SD = THESIS / 'figures/chapter3/source_data'
OUT = WGS / 'figures'
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)

GROUPS = [
    ('HT1080', [('HT1080-sgCTRL-D0', 'sgCTRL', '#3C5488'),
                ('HT1080-LOY-D0', 'CP + Q53 2×', '#E64B35')]),
    ('A549', [('A549-CTRL2', 'sgCTRL (CTRL2)', '#3C5488'),
              ('A549-LOY', 'CP + Q53 2×', '#E64B35'),
              ('A549-BC14', 'BC14 (undocumented; mosdepth, not QDNAseq)', '#7E6148')]),
    ('DU145', [('LNCaP-sgCTRL', 'sgCTRL', '#3C5488'),
               ('LNCaP-LOY-2X', '2×, no drug', '#705A8C'),
               ('LNCaP-DMSO-LOY-1X', '1× DMSO', '#777777'),
               ('LNCaP-30', '1× pano 30 nM', '#E99B57'),
               ('LNCaP-70', '1× pano 75 nM', '#E64B35'),
               ('LNCaP-KU', '1× KU', '#4DBBD5'),
               ('LNCaP-KU-30', '1× KU + pano 30 nM', '#00A087'),
               ('LNCaP-KU-70', '1× KU + pano 75 nM', '#8491B4')]),
]
LABEL = {s: l for _, ss in GROUPS for s, l, _ in ss}
COLOR = {s: c for _, ss in GROUPS for s, _, c in ss}
LINE = {s: g for g, ss in GROUPS for s, _, _ in ss}

# --- QDNAseq bins (frozen thesis exports of WGS/*_QDNAseq_cn_df.rds) ---
frames = []
for f in ['ht1080_lpwgs_confirmation/qdnaseq_bins.csv', 'a549_lpwgs_confirmation/qdnaseq_bins.csv',
          'du145_drug_optimisation/lpwgs/chrY_qdnaseq_bins.csv']:
    d = pd.read_csv(SD / f)
    frames.append(d[d.chromosome.astype(str).isin(['chrY', 'Y'])])
w = pd.concat(frames, ignore_index=True)
w = w[w.use.eq(True) & w.log2_seg.notna()].copy()
w['relative_cn'] = 2 ** w.log2
w['segment_cn'] = 2 ** w.log2_seg
w['method'] = 'QDNAseq'
grid = w[w['sample'].eq('A549-CTRL2')][['start', 'end', 'midpoint']].sort_values('start').reset_index(drop=True)


def mosdepth_chrY(sample):
    """Mean 1 kb MAPQ>=20 depth per 100 kb QDNAseq bin, divided by the autosomal median."""
    rows = []
    with gzip.open(WGS / f'tracks/{sample}.regions.bed.gz', 'rt') as fh:
        for line in fh:
            if line.startswith('chrY\t'):
                c, a, b, v = line.split('\t'); rows.append((int(a), float(v)))
            elif rows:
                break
    t = pd.DataFrame(rows, columns=['start0', 'depth'])
    med = float((WGS / f'tracks/{sample}.median.txt').read_text().strip())
    t['bin'] = t.start0 // 100_000
    agg = t.groupby('bin').depth.mean() / med
    g = grid.copy(); g['relative_cn'] = ((g.start - 1) // 100_000).map(agg).to_numpy()
    g['sample'] = sample; g['method'] = 'mosdepth'; g['segment_cn'] = np.nan
    return g, med


bc14, bc14_med = mosdepth_chrY('A549-BC14')
ctrl2_md, _ = mosdepth_chrY('A549-CTRL2')
w = pd.concat([w, bc14], ignore_index=True)

# Scale check: mosdepth vs QDNAseq on A549-CTRL2, male-specific euchromatin (PAR1 end to Yq12 start).
cyto = pd.read_csv(SD / 'ht1080_y_expression_confirmation/T2T_chrY_cytobands.csv')
cyto = cyto[cyto.chrom.eq('chrY')]
genes = pd.read_csv(SD / 'ht1080_lpwgs_confirmation/ddpcr_target_loci.csv')
LENGTH = cyto.end.max() / 1e6
C0 = cyto.loc[cyto.stain.eq('acen'), 'start'].min() / 1e6
C1 = cyto.loc[cyto.stain.eq('acen'), 'end'].max() / 1e6
Q0 = cyto.loc[cyto.band.eq('q12'), 'start'].iloc[0] / 1e6
q = w[w['sample'].eq('A549-CTRL2')].set_index('start').relative_cn
m = ctrl2_md.set_index('start').relative_cn
msy = (ctrl2_md.set_index('start').midpoint.between(2.46, Q0))
ratio = (m[msy] / q.reindex(m[msy].index)).replace([np.inf, -np.inf], np.nan).dropna()
scale_check = dict(n_bins=int(len(ratio)), median_mosdepth_over_qdnaseq=float(ratio.median()),
                   ctrl2_mosdepth_msy_median=float(m[msy].median()), ctrl2_qdnaseq_msy_median=float(q.reindex(m[msy].index).median()),
                   bc14_mosdepth_msy_median=float(bc14.set_index('start').relative_cn[msy.values].median()))

# --- Saved segments: merge runs of identical log2_seg (as in the thesis script) ---
SEG = []
for s, g in w[w.method.eq('QDNAseq')].groupby('sample', sort=False):
    g = g.sort_values('start').copy()
    g['run'] = (g.log2_seg.ne(g.log2_seg.shift()) | g.start.gt(g.end.shift() + 1)).cumsum()
    for _, h in g.groupby('run'):
        SEG.append(dict(sample=s, start=(h.start.min() - 1) / 1e6, end=h.end.max() / 1e6, value=h.segment_cn.iloc[0]))
seg = pd.DataFrame(SEG)

counts = w.groupby('sample').size()
expected = [s for _, ss in GROUPS for s, _, _ in ss]
missing = [s for s in expected if s not in counts.index]
assert not missing, f'missing samples: {missing}'
print('Bins per sample:\n' + counts.reindex(expected).to_string())


def context(ax):
    for a, b in [(0, 2.458320), (62.122809, LENGTH)]: ax.axvspan(a, b, color='#EAF6F7', zorder=0)
    for a, b in [(C0, C1), (Q0, 62.122809)]: ax.axvspan(a, b, color='#ECECEC', zorder=0)


def track(ax, s, first):
    context(ax)
    g = w[w['sample'].eq(s)]
    ax.scatter(g.midpoint, np.minimum(g.relative_cn, 1.27), s=2.5, color='#999999', alpha=.65, lw=0)
    high = g[g.relative_cn > 1.27]
    ax.scatter(high.midpoint, [1.27] * len(high), marker='^', s=7, color='#555555', lw=0)
    for r in seg[seg['sample'].eq(s)].itertuples():
        ax.plot([r.start, r.end], [r.value] * 2, color=COLOR[s], lw=1.2)
    ax.set(xlim=(0, LENGTH), ylim=(-.025, 1.35), yticks=[0, .5, 1])
    ax.tick_params(axis='x', labelbottom=False, bottom=False)
    ax.tick_params(axis='y', labelsize=6.7, pad=2)
    ax.text(.99, .81, LABEL[s], transform=ax.transAxes, ha='right', color=COLOR[s], fontsize=8, fontweight='semibold', zorder=5,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=.6))
    ax.text(.99, .60, s, transform=ax.transAxes, ha='right', color='#888888', fontsize=6, zorder=5,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=.4))
    if first: ax.text(.66, .82, 'Yq12: limited mappability', transform=ax.transAxes, ha='center', fontsize=7, color='#777777')


def ideogram(ax):
    for b in cyto.itertuples():
        lo, hi = b.start / 1e6, b.end / 1e6
        if b.stain == 'acen':
            pts = [(lo, .2), (lo, .65), (hi, .425)] if b.band.startswith('p') else [(lo, .425), (hi, .2), (hi, .65)]
            ax.add_patch(Polygon(pts, facecolor='#888888', edgecolor='#555555', lw=.4))
        else:
            ax.add_patch(Rectangle((lo, .2), hi - lo, .45, facecolor={'gneg': 'white', 'gpos50': '#999999', 'gvar': '#AFC0CD'}.get(b.stain, '#BBBBBB'), edgecolor='#555555', lw=.4))
    for lo, hi in [(0, 2.458320), (62.122809, LENGTH)]: ax.add_patch(Rectangle((lo, .2), hi - lo, .45, color='#41B6C4', lw=0))
    for r in genes.itertuples():
        pos = (r.t2t_start + r.t2t_end) / 2e6
        ax.plot([pos, pos], [.68, 1.0], color='#555555', lw=.5)
        ax.text(pos, 1.05, r.gene_name, ha='left', rotation=35, fontsize=6.5, fontstyle='italic')
    ax.set(xlim=(0, LENGTH), ylim=(0, 2.2), yticks=[], xticks=np.arange(0, 61, 10), xlabel='HG002 T2T-Y position (Mb)')
    ax.spines[['left', 'right', 'top']].set_visible(False); ax.tick_params(axis='x', labelsize=7)


# --- Layout in mm, top-down ---
W, TOP, HEAD, TRACK, GAP, IDEO, AFTER, FOOT = 183, 12, 8, 17, 2.5, 10, 16, 16
L, R = .10, .95
H = TOP + sum(HEAD + len(ss) * (TRACK + GAP) + IDEO + AFTER for _, ss in GROUPS) + FOOT
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
fy = lambda mm: 1 - mm / H
fig.text(L, fy(5), 'All lpWGS runs | chromosome Y relative copy number, grouped by cell line',
         fontsize=9, fontweight='semibold', va='top')
y = TOP
for k, (line, ss) in enumerate(GROUPS):
    fig.text(.025, fy(y), 'abc'[k], fontsize=12, fontweight='semibold', va='top')
    fig.text(L, fy(y), line + (' (sample IDs carry legacy LNCaP names)' if line == 'DU145' else ''),
             fontsize=9, fontweight='semibold', va='top')
    y += HEAD; top_tracks = y
    for i, (s, _, _) in enumerate(ss):
        ax = fig.add_axes([L, fy(y + TRACK), R - L, TRACK / H]); track(ax, s, i == 0)
        y += TRACK + GAP
    fig.text(.036, fy((top_tracks + y) / 2), 'Relative copy number', rotation=90, va='center', fontsize=8)
    ideogram(fig.add_axes([L, fy(y + IDEO), R - L, IDEO / H]))
    y += IDEO + AFTER
fig.text(L, fy(y - 4), 'Relative copy number = 2^log2 of the QDNAseq ratio to the genome median (T2T-CHM13v2.0, 100 kb bins); lines are saved segments.\n'
         'Each cell line was normalised and segmented separately, so 1.0 does not mean the same absolute copy number across lines (A549 near-triploid).\n'
         f'A549-BC14: mosdepth 1 kb MAPQ≥20 depth averaged per 100 kb bin ÷ autosomal median ({bc14_med:g}); no GC correction or segments. '
         f'On A549-CTRL2 this\nmeasure is {scale_check["median_mosdepth_over_qdnaseq"]:.2f}× the QDNAseq value across male-specific euchromatin '
         f'({scale_check["n_bins"]} bins). DU145 "pano 75 nM" = raw sample ID LNCaP-70.',
         fontsize=6.3, color='#666666', va='top', linespacing=1.4)

OUT.mkdir(exist_ok=True)
stem = OUT / 'all_lpwgs_chrY_comparison'
fig.savefig(stem.with_suffix('.pdf'), facecolor='white')
fig.savefig(stem.with_suffix('.png'), dpi=160, facecolor='white')
plt.close(fig)
w.assign(cell_line=w['sample'].map(LINE), display_label=w['sample'].map(LABEL))[
    ['cell_line', 'sample', 'display_label', 'method', 'start', 'end', 'midpoint', 'relative_cn', 'segment_cn']
].to_csv(OUT / 'all_lpwgs_chrY_comparison_bins.csv', index=False)
seg.to_csv(OUT / 'all_lpwgs_chrY_comparison_segments.csv', index=False)
print('Scale check:', scale_check)
print(f'Wrote {stem}.pdf ({W} x {H:.0f} mm)')
