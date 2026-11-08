"""Render ch04_supp_external_isogenicz_comparison: IsogenicZ (T0-adjusted) versus endpoint drugZ, screen by screen.

Inputs are saved scores only; no screen is rescored here.
- External (27): IsogenicZ from research/external_isogenicz_20261022 and the 4.3.3 endpoint drugZ (Supp Fig 4.8).
- HT1080: the mask173_detected IsogenicZ and matched endpoint drugZ fits shown in Fig 4.3b.
- A549 P4: the mask342 + D0-covered IsogenicZ and matched endpoint drugZ fits shown in Supp Fig 4.7h.
"""
from pathlib import Path
import hashlib
import json
import sys

THESIS = Path(__file__).resolve().parents[3]
DEPS = THESIS.parent / '.figure_python_deps'
if DEPS.is_dir():
    sys.path.insert(0, str(DEPS))
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pymupdf
from thesis_style import apply_style
apply_style(7)

FIGS = THESIS / 'figures/chapter4'
SD = FIGS / 'source_data'
EXT = SD / 'external_isogenicz_20261022'
OUT = EXT / 'figure'
REVIEW = FIGS / 'review'
STEM = 'ch04_supp_external_isogenicz_comparison'
FDR = .25
BOTH, DZ_ONLY, ISO_ONLY, BG = '#7E3F98', '#7F8790', '#E9A23B', '#DDE0E4'
STUDY = {'Desjardins': '#00A087', 'Feng': '#8491B4', 'ours': '#222222'}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def label(cid):
    if cid in ('HT1080', 'A549 P4'):
        return cid
    study, alt, line = cid.split('_')
    return f'{alt} {line}' if study.startswith('Desj') else f'{alt} HEK293A'


def load():
    summ = pd.read_csv(EXT / 'contrast_summary.csv').set_index('id')
    pairs, meta, src = {}, {}, []
    for cid in summ.index:
        f1, f2 = EXT / f'{cid}_IsogenicZ.csv', SD / f'external_isogenic_screens_20261021/drugz/{cid}_drugz.txt'
        pairs[cid] = (pd.read_csv(f1).set_index('Gene'), pd.read_csv(f2, sep='\t').set_index('GENE'))
        meta[cid] = ('Desjardins' if cid.startswith('Desj') else 'Feng', summ.at[cid, 't0_spearman'])
        src += [f1, f2]
    for cid, folder, iso, dz, raw, d0 in [
            ('HT1080', 'ht1080_isogenicz_20261020', 'mask173_detected_IsogenicZ_genes.tsv',
             'mask173_detected_matched_endpoint_drugZ_genes.tsv', 'mask173_detected_raw_counts.tsv', ('D0_CTRL', 'D0_LOY')),
            ('A549 P4', 'a549_isogenicz_20261020', 'P4_mask342_D0covered_IsogenicZ_genes.tsv',
             'P4_mask342_D0covered_matched_endpoint_drugZ_genes.tsv', 'P4_mask342_D0covered_raw_counts.tsv',
             ('CTRL_D0_pool', 'LOY_D0_pool'))]:
        a, b, r = (SD / folder / x for x in (iso, dz, raw))
        counts = pd.read_csv(r, sep='\t')
        pairs[cid] = (pd.read_csv(a, sep='\t').set_index('Gene'), pd.read_csv(b, sep='\t').set_index('Gene'))
        meta[cid] = ('ours', float(spearmanr(counts[d0[0]], counts[d0[1]]).statistic))
        src += [a, b, r]
    return pairs, meta, src


def summarise(pairs, meta):
    rows, points = [], []
    for cid, (iso, dz) in pairs.items():
        sh = iso.index.intersection(dz.index)
        i, d = iso.loc[sh], dz.loc[sh]
        si, sd_ = i.fdr_synth < FDR, d.fdr_synth < FDR
        rows.append(dict(screen=cid, label=label(cid), study=meta[cid][0], genes=len(sh),
                         spearman_normZ=spearmanr(i.normZ, d.normZ).statistic, t0_spearman=meta[cid][1],
                         both=int((si & sd_).sum()), drugz_only=int((sd_ & ~si).sum()), isogenicz_only=int((si & ~sd_).sum())))
        points.append(pd.DataFrame({'screen': cid, 'Gene': sh, 'drugZ_normZ': d.normZ.values, 'IsogenicZ_normZ': i.normZ.values,
                                    'drugZ_fdr_dep': d.fdr_synth.values, 'IsogenicZ_fdr_dep': i.fdr_synth.values}))
    t = pd.DataFrame(rows)
    ext = t[t.study != 'ours']
    test = spearmanr(ext.t0_spearman, ext.spearman_normZ)
    return t, pd.concat(points, ignore_index=True), dict(rho=float(test.statistic), p=float(test.pvalue), n=len(ext))


def panel_a(ax, t, test):
    for st, col in STUDY.items():
        d = t[t.study == st]
        ax.scatter(d.t0_spearman, d.spearman_normZ, s=16 if st == 'ours' else 13, color=col, lw=0, zorder=3,
                   marker='D' if st == 'ours' else 'o')
    offs = {'HT1080': (4, -1), 'A549 P4': (4, -1), 'TET2 RPE1': (4, -1), 'IDH1 RPE1': (4, -1), 'STK11 RPE1': (4, -1),
            'CDK12 RPE1': (-4, 0), 'Chr18q RPE1': (4, 2), 'TP53 HEK293A': (4, -4)}
    for _, r in t.iterrows():
        if r.label in offs:
            dx, dy = offs[r.label]
            ax.annotate(r.label, (r.t0_spearman, r.spearman_normZ), xytext=(dx, dy), textcoords='offset points',
                        fontsize=5.8, ha='right' if dx < 0 else 'left', va='center',
                        fontweight='semibold' if r.study == 'ours' else 'normal')
    ax.set(xlabel='Spearman ρ between the two arms’ T0 libraries',
           ylabel='Spearman ρ, IsogenicZ vs drugZ normZ', xlim=(.82, 1.0), ylim=(.66, 1.0))
    ax.text(.97, .04, f"External contrasts: ρ = {test['rho']:.2f}, P = {test['p']:.4f} (n = {test['n']})",
            transform=ax.transAxes, ha='right', fontsize=6)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=STUDY['Desjardins'], ms=3.5, label='Desjardins et al. (2026)'),
                       Line2D([], [], marker='o', ls='', color=STUDY['Feng'], ms=3.5, label='Feng et al. (2022)'),
                       Line2D([], [], marker='D', ls='', color=STUDY['ours'], ms=3.5, label='This study')],
              loc='upper left', frameon=False, fontsize=6, handletextpad=.2)


def panel_b(ax, t):
    tt = t.sort_values('drugz_only', key=lambda s: t.loc[s.index, 'both'] + s)
    y = np.arange(len(tt))
    ax.barh(y, tt.both, color=BOTH, height=.72)
    ax.barh(y, tt.drugz_only, left=tt.both, color=DZ_ONLY, height=.72)
    ax.barh(y, tt.isogenicz_only, left=tt.both + tt.drugz_only, color=ISO_ONLY, height=.72)
    ax.set_yticks(y, tt.label, fontsize=5.4)
    for lab, st in zip(ax.get_yticklabels(), tt.study):
        lab.set_fontweight('semibold' if st == 'ours' else 'normal')
    ax.tick_params(axis='y', length=0, pad=1.5)
    ax.set(xlabel='Depleted genes (FDR <0.25)', ylim=(-.6, len(tt) - .4))
    ax.legend(handles=[Patch(color=BOTH, label='Both methods'), Patch(color=DZ_ONLY, label='drugZ only'),
                       Patch(color=ISO_ONLY, label='IsogenicZ only')], loc='lower right', frameon=False, fontsize=6)


def grid(fig, gs, t, points):
    order = t.sort_values('spearman_normZ').screen.tolist()
    axs = []
    for k, cid in enumerate(order):
        ax = fig.add_subplot(gs[k // 6, k % 6])
        p = points[points.screen == cid]
        dz, iz = p.drugZ_fdr_dep < FDR, p.IsogenicZ_fdr_dep < FDR
        ax.scatter(p.drugZ_normZ[~dz & ~iz], p.IsogenicZ_normZ[~dz & ~iz], s=.25, color=BG, lw=0, rasterized=True)
        for m, col in [(dz & ~iz, DZ_ONLY), (iz & ~dz, ISO_ONLY), (dz & iz, BOTH)]:
            ax.scatter(p.drugZ_normZ[m], p.IsogenicZ_normZ[m], s=1.6, color=col, lw=0, rasterized=True)
        lo = min(p.drugZ_normZ.min(), p.IsogenicZ_normZ.min()); hi = max(p.drugZ_normZ.max(), p.IsogenicZ_normZ.max())
        ax.plot([lo, hi], [lo, hi], color='.55', lw=.4, ls=(0, (2, 2)))
        r = t.set_index('screen').loc[cid]
        ax.set_title(r.label, fontsize=6, pad=1.5, fontweight='semibold' if r.study == 'ours' else 'normal')
        ax.text(.04, .96, f'ρ {r.spearman_normZ:.2f}\n{r.both}/{r.drugz_only}/{r.isogenicz_only}',
                transform=ax.transAxes, va='top', fontsize=5.2, linespacing=1.1)
        ax.tick_params(labelsize=5, length=1.5, pad=1)
        axs.append(ax)
    return axs


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pairs, meta, src = load()
    t, points, test = summarise(pairs, meta)
    assert len(t) == 29
    ht = t.set_index('screen').loc['HT1080']
    assert ht.both + ht.drugz_only == 48  # matches the 48 primary HT1080 sensitizers at FDR <0.25
    t.to_csv(OUT / 'per_screen_agreement.csv', index=False)
    points.to_csv(OUT / 'plotted_gene_scores.csv.gz', index=False)

    fig = plt.figure(figsize=(7.2, 9.6))
    top = fig.add_gridspec(1, 2, left=.085, right=.985, top=.975, bottom=.70, wspace=.42, width_ratios=[1, .9])
    a, b = fig.add_subplot(top[0]), fig.add_subplot(top[1])
    panel_a(a, t, test)
    panel_b(b, t)
    gs = fig.add_gridspec(5, 6, left=.065, right=.985, top=.635, bottom=.055, hspace=.55, wspace=.38)
    axs = grid(fig, gs, t, points)
    fig.text(.525, .025, 'Endpoint drugZ normZ', ha='center', fontsize=7)
    fig.text(.012, .345, 'IsogenicZ normZ (T0-adjusted)', rotation=90, va='center', fontsize=7)
    for ax, letter, x in [(a, 'a', -.2), (b, 'b', -.38)]:
        ax.text(x, 1.0, letter, transform=ax.transAxes, fontsize=10, fontweight='semibold', va='bottom')
    axs[0].text(-.5, 1.18, 'c', transform=axs[0].transAxes, fontsize=10, fontweight='semibold', va='bottom')
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(FIGS / f'{STEM}.{ext}', dpi=600)
    preview = REVIEW / f'{STEM}_preview.png'
    fig.savefig(preview, dpi=170)
    plt.close(fig)
    hashed = REVIEW / f'{STEM}_{sha(preview)[:12]}.png'
    hashed.write_bytes(preview.read_bytes())
    with pymupdf.open(FIGS / f'{STEM}.pdf') as pdf:
        page = pdf[0]
        fonts = [f[3] for f in page.get_fonts()]
        assert all('OpenSans' in f.replace(' ', '') for f in fonts), fonts
        outside = [w[4] for w in page.get_text('words')
                   if w[0] < 0 or w[1] < 0 or w[2] > page.rect.width or w[3] > page.rect.height]
        assert not outside, outside
    manifest = dict(stem=STEM, preview=str(hashed), fdr=FDR, t0_vs_agreement_external=test,
                    hit_overlap_definition='depleted at FDR <0.25 (fdr_synth) in each method',
                    sources={str(p): sha(p) for p in src}, script_sha256=sha(__file__), fonts=fonts)
    (OUT / 'render_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(t.sort_values('spearman_normZ')[['label', 'spearman_normZ', 't0_spearman', 'both', 'drugz_only', 'isogenicz_only']]
          .round(3).to_string(index=False))
    print(json.dumps(test), hashed)


if __name__ == '__main__':
    main()
