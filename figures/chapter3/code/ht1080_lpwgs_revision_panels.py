"""Render the October 2026 main/supplement split from unchanged F03 summaries."""
import io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch
from adjustText import adjust_text

RED, BLUE, GREEN = '#E64B35', '#3C5488', '#00A087'
SUPP_STEM = 'ch03_supp_ht1080_lpwgs_copy_number_detail'
RNA_YLIM = (-8.3, 2.5)
RNA_YTICKS = [-8, -6, -4, -2, 0, 2]


def export(fig, root, stem):
    for ext in ('pdf', 'svg', 'png'):
        buf = io.BytesIO()
        fig.savefig(buf, format=ext, dpi=600, facecolor='white')
        (root/f'{stem}.{ext}').write_bytes(buf.getvalue())
    # A smaller review preview avoids loading the full 600-dpi raster.
    fig.savefig(root/'source_data/ht1080_lpwgs_difference_draft'/f'{stem}_preview.png',
                dpi=160, facecolor='white')


def clean(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=8)


def chromosome_key():
    return [Line2D([], [], marker='o', ls='', color='#777777', ms=4, label='Autosomes / X'),
            Line2D([], [], marker='^', ls='', color=RED, ms=5, label='Chromosome Y')]


def render_arm_scatter(ax, z):
    isy = z.chromosome.eq('chrY'); nony = z.loc[~isy]
    ax.axhline(0, color='#DDDDDD', lw=.6, zorder=0)
    ax.axvline(0, color='#DDDDDD', lw=.6, zorder=0)
    ax.scatter(z.DNA, z.RNA, s=20, color='#777777', alpha=.85, lw=0)
    ax.scatter(z.loc[isy, 'DNA'], z.loc[isy, 'RNA'], s=30, marker='^', color=RED, lw=0)
    score = (((nony[['DNA', 'RNA']]-nony[['DNA', 'RNA']].median())/nony[['DNA', 'RNA']].std())**2).sum(axis=1)
    selected = list(dict.fromkeys(score.nlargest(4).index.tolist()+z.index[z.arm.eq('13q') | isy].tolist()))
    texts, scatter_labels = [], []
    for i in selected:
        r = z.loc[i]
        color = RED if r.chromosome == 'chrY' else GREEN if r.arm == '13q' else BLUE
        ax.scatter([r.DNA], [r.RNA], s=30, marker='^' if r.chromosome == 'chrY' else 'o', color=color, lw=0)
        texts.append(ax.text(r.DNA, r.RNA, r.arm, fontsize=7, color=color))
        scatter_labels.append(dict(arm=r.arm, DNA=r.DNA, RNA=r.RNA))
    ax.set_xscale('symlog', linthresh=.3)
    ax.set_yscale('linear')
    ax.set(xlim=(-12, 1.05), ylim=RNA_YLIM, xlabel='DNA log$_2$ difference',
          ylabel='RNA mean log$_2$ CPM difference')
    ax.set_xticks([-10, -1, -.3, 0, .3, 1], ['−10', '−1', '−0.3', '0', '0.3', '1'])
    ax.set_yticks(RNA_YTICKS)
    np.random.seed(42)
    adjust_text(texts, x=z.DNA.to_numpy(), y=z.RNA.to_numpy(), ax=ax, expand=(1.3, 1.4),
                arrowprops=dict(arrowstyle='-', color='#888888', lw=.5))
    ax.text(.97, .035, f'Chromosome arms\nr (all) = {z.DNA.corr(z.RNA):.3f}; n = {len(z)}\n'
           f'r (non-Y) = {nony.DNA.corr(nony.RNA):.3f}; n = {len(nony)}',
           transform=ax.transAxes, ha='right', va='bottom', fontsize=7,
           bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=1))

    clean(ax)
    ax.tick_params(labelsize=7)
    return scatter_labels


def render_main(root, stem, bins, windows, ydetail, cb, chroms, lengths, offsets, end, yq, arms_rna):
    fig = plt.figure(figsize=(9.4, 6.2))
    # The small Y view is attached to a and sits above it, outside the data area.
    ax = fig.add_axes([.075, .17, .535, .435])
    inset = fig.add_axes([.19, .715, .40, .225])
    scatter = fig.add_axes([.735, .57, .245, .37])
    arm_scatter = fig.add_axes([.735, .075, .245, .37], sharey=scatter)
    for i, ch in enumerate(chroms):
        z = bins[bins.chromosome.eq(ch)]
        ax.scatter(z.mid/1e6+offsets[ch], z.difference,
                   s=20 if ch == 'chrY' else 8, marker='^' if ch == 'chrY' else 'o',
                   color=RED if ch == 'chrY' else ['#555555', '#AAAAAA'][i % 2], lw=0,
                   zorder=4 if ch == 'chrY' else 2)
    ax.axhline(0, color='#BBBBBB', lw=.5, zorder=0)
    ax.axvspan(offsets['chrY']+yq.start/1e6, offsets['chrY']+yq.end/1e6,
               color='#EEEEEE', zorder=-1)
    ax.annotate('Yq12 excluded', xy=(offsets['chrY']+(yq.start+yq.end)/2e6, 0),
                xytext=(.985, .95), textcoords='axes fraction', ha='right', va='top',
                fontsize=7, color='#666666',
                arrowprops=dict(arrowstyle='-', lw=.5, color='#888888'))
    extent = max(10, np.ceil(bins.difference.abs().max()))
    ax.set(xlim=(-10, end), ylim=(-extent-.5, extent+.5), yticks=[-10, -5, 0, 5, 10],
           ylabel='$\u0394$ log$_2$(relative CN)', xlabel='Chromosome')
    ax.set_xticks([offsets[ch]+lengths[ch]/2e6 for ch in chroms], [ch[3:] for ch in chroms])
    clean(ax)
    ax.tick_params(axis='x', labelsize=7, pad=2, rotation=45)
    ax.get_xticklabels()[-1].set_color(RED)
    ax.legend(handles=chromosome_key(), loc='upper left',
              ncol=2, frameon=False, fontsize=7.5, borderaxespad=.5,
              handletextpad=.4, columnspacing=1.3)
    ax.text(.02, .035, '5-Mb windows; Y-depleted − sgCTRL', transform=ax.transAxes,
            fontsize=7.3, va='bottom')
    fig.text(.018, .948, 'a', fontsize=12, weight='semibold')

    cen = cb[cb.chrom.eq('chrY') & cb.stain.eq('acen')]
    c0, c1 = cen.start.min()/1e6, cen.end.max()/1e6
    inset.axvspan(c0, c1, color='#DDDDDD', zorder=0)
    inset.axhline(0, color='#CCCCCC', lw=.4, zorder=0)
    # Exponentiating the existing per-sample log2 medians puts the inset on
    # the original relative-CN scale. P and q are drawn separately across the gap.
    ymax = 0.
    for label, color, marker in [('sgCTRL', BLUE, 'o'), ('LOY', RED, '^')]:
        ymax = max(ymax, float((2**ydetail[label]).max()))
        for _, z in ydetail.groupby('arm', sort=False):
            inset.plot(z.mid/1e6, 2**z[label], color=color, marker=marker,
                       ms=2.5, lw=.65, mew=0, clip_on=False)
    inset.set(xlim=(0, yq.start/1e6), ylim=(-.035, max(1.0, ymax*1.1)),
              xticks=[0, 5, 10, 15, 20, 25], yticks=[0, .5, 1],
              xlabel='Chromosome Y position (Mb)', ylabel='Relative CN')
    inset.set_xlabel('Chromosome Y position (Mb)', fontsize=7.5, labelpad=1)
    inset.set_ylabel('Relative CN', fontsize=7.5, labelpad=3)
    clean(inset)
    inset.tick_params(labelsize=7, pad=1, length=2)
    inset.text(c0/2, 1.01, 'Yp', ha='center', fontsize=7,
               transform=inset.get_xaxis_transform())
    inset.text((c1+yq.start/1e6)/2, 1.01, 'Yq (euchromatic)', ha='center', fontsize=7,
               transform=inset.get_xaxis_transform())
    inset.legend(handles=[Line2D([], [], color=BLUE, marker='o', ms=3, lw=.6, label='sgCTRL'),
                          Line2D([], [], color=RED, marker='^', ms=3, lw=.6, label='Y-depleted')],
                 loc='upper right', fontsize=6.7, frameon=False, ncol=2,
                 columnspacing=.8, handlelength=1.1, handletextpad=.4, borderpad=.2)
    fig.add_artist(ConnectionPatch(xyA=(offsets['chrY']+yq.start/1e6, extent+.5),
                                   coordsA=ax.transData, xyB=(1, 0), coordsB=inset.transAxes,
                                   color='#AAAAAA', lw=.5, zorder=0))

    isy = windows.chrom.eq('chrY')
    nony = windows.loc[~isy]
    scatter.axhline(0, color='#DDDDDD', lw=.5, zorder=0)
    scatter.axvline(0, color='#DDDDDD', lw=.5, zorder=0)
    scatter.scatter(nony.DNA, nony.RNA, s=8, color='#777777', alpha=.55, lw=0)
    scatter.scatter(windows.loc[isy, 'DNA'], windows.loc[isy, 'RNA'],
                    s=27, marker='^', color=RED, lw=0, zorder=4)
    scatter.set(xlim=(-10, 1.2), ylim=RNA_YLIM, xticks=[-10, -8, -6, -4, -2, 0],
                yticks=RNA_YTICKS, xlabel='DNA log$_2$ difference',
                ylabel='RNA mean log$_2$ CPM difference')
    scatter.legend(handles=chromosome_key(), loc='upper left', frameon=False, fontsize=7.2,
                   handletextpad=.4, borderpad=.2, labelspacing=.3)
    scatter.text(.97, .035,
                 f'1-Mb windows\nr (all) = {windows.DNA.corr(windows.RNA):.3f}; n = {len(windows):,}'
                 f'\nr (non-Y) = {nony.DNA.corr(nony.RNA):.3f}; n = {len(nony):,}',
                 transform=scatter.transAxes, ha='right', va='bottom', fontsize=7.2)
    clean(scatter)
    fig.text(.67, .948, 'b', fontsize=12, weight='semibold')
    assert ax.get_yscale() == scatter.get_xscale() == scatter.get_yscale() == 'linear'
    for values, limits in [(bins.difference, ax.get_ylim()), (windows.DNA, scatter.get_xlim()),
                           (windows.RNA, scatter.get_ylim())]:
        assert values.between(*limits).all(), 'A main-panel value would be clipped'
    scatter_labels = render_arm_scatter(arm_scatter, arms_rna)
    assert scatter.get_yscale() == arm_scatter.get_yscale() == 'linear'
    assert scatter.get_ylim() == arm_scatter.get_ylim() == RNA_YLIM
    assert np.array_equal(scatter.get_yticks(), arm_scatter.get_yticks())
    assert np.isclose(scatter.get_position().height, arm_scatter.get_position().height)
    assert arms_rna.RNA.between(*RNA_YLIM).all()
    fig.text(.67, .448, 'c', fontsize=12, weight='semibold')
    export(fig, root, stem)
    plt.close(fig)
    return scatter_labels


def render_supp(root, arms, detail, chr4_detail, cb, lengths, segment_delta):
    fig = plt.figure(figsize=(9.4, 5.4))
    gs = fig.add_gridspec(2, 2, left=.075, right=.985, top=.94, bottom=.11, hspace=.58, wspace=.25)
    a = fig.add_subplot(gs[0, :])
    c = fig.add_subplot(gs[1, 0]); d = fig.add_subplot(gs[1, 1])
    x = np.arange(len(arms))
    bars = a.bar(x, arms.difference, width=.78,
                 color=[GREEN if arm == '13q' else '#666666' for arm in arms.arm])
    for bar, low in zip(bars, arms.low_mappability_arm):
        if low:
            bar.set_facecolor('#DDDDDD'); bar.set_edgecolor('#999999'); bar.set_hatch('///')
    for i in range(2, len(arms), 2):
        a.axvline(i-.5, color='#EEEEEE', lw=.5, zorder=0)
    a.axhline(0, color='#333333', lw=.5)
    a.set(xlim=(-.8, len(arms)-.2), ylim=(-.41, .43),
          ylabel='Median $\u0394$ log$_2$(relative CN)', xlabel='Chromosome arm')
    a.set_xticks(x, arms.arm, rotation=90)
    a.text(.015, .95, 'Autosomes and X', transform=a.transAxes, va='top', fontsize=8)
    ix = arms.index[arms.arm.eq('13q')][0]
    a.text(ix, arms.loc[ix, 'difference']+.025, '13q', ha='center', color='#007C68', fontsize=8)

    for ax, chrom, table in [(c, 'chr13', detail), (d, 'chr4', chr4_detail)]:
        cen = cb[cb.chrom.eq(chrom) & cb.stain.eq('acen')]
        ax.axvspan(cen.start.min()/1e6, cen.end.max()/1e6, color='#EEEEEE', zorder=0)
        ax.axhline(0, color='#888888', lw=.5)
        ax.scatter(table.mid/1e6, table.difference, s=13, color='#555555', lw=0, zorder=3)
        ax.set(xlim=(-2, lengths[chrom]/1e6+2), ylim=(-1, 1.1), yticks=[-1, -.5, 0, .5, 1],
               xlabel=f'Chromosome {chrom[3:]} position (Mb)', ylabel='$\u0394$ log$_2$(relative CN)')
        assert table.difference.between(*ax.get_ylim()).all()
    c.axvline(34.7, color='#AAAAAA', lw=.6, ls=':')
    c.plot([34.7, 113.566686], [segment_delta]*2, color=GREEN, lw=1.5)
    c.text(77, .72, f'$\u0394$ log$_2$ = {segment_delta:.3f}  (1.29-fold)',
           color='#007C68', fontsize=8, ha='center')
    c.text(.025, .04, '1-Mb windows; green line: segment difference', transform=c.transAxes, fontsize=7)
    cen4 = cb[cb.chrom.eq('chr4') & cb.stain.eq('acen')]
    for arm, left, right in [('4p', 0, cen4.start.min()/1e6),
                             ('4q', cen4.end.max()/1e6, lengths['chr4']/1e6)]:
        value = float(arms.loc[arms.arm.eq(arm), 'difference'].iloc[0])
        d.plot([left, right], [value, value], color=BLUE, lw=1.3)
        d.text((left+right)/2, .65, f'{arm}: {value:+.3f}\n({2**value:.2f}-fold)',
               ha='center', color=BLUE, fontsize=8)
    d.text(.025, .04, '1-Mb windows; blue lines: arm medians', transform=d.transAxes, fontsize=7)
    for ax, letter in zip((a, c, d), 'abc'):
        clean(ax)
        ax.text(-.055/ax.get_position().width, 1.045, letter, transform=ax.transAxes,
                fontsize=12, weight='semibold')
    a.tick_params(axis='x', labelsize=7, pad=2)
    for tick, low in zip(a.get_xticklabels(), arms.low_mappability_arm):
        if low: tick.set_color('#999999')
    export(fig, root, SUPP_STEM)
    plt.close(fig)
