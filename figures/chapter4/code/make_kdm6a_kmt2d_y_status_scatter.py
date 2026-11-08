"""Show KDM6A/KMT2D Chronos co-dependencies in the saved XY and LOY cohorts."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS/'figures/chapter4'
PREVIOUS = OUT/'source_data/xy_paralog_loy_ranked_partners'
DATA = OUT/'source_data/kdm6a_kmt2d_y_status_scatter'
SOURCE = Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv')
STEM = 'ch04_supp_kdm6a_kmt2d_xy_loy_scatter_draft'
sys.path.insert(0, str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
from prepare_xy_ranked_partners_by_y_status import sha256


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    previous = json.loads((PREVIOUS/'analysis_manifest.json').read_text())
    digest = sha256(SOURCE)
    assert digest == previous['sources'][str(SOURCE)]
    cohort_path = PREVIOUS/'model_cohort_and_exclusions.csv'
    cohort = pd.read_csv(cohort_path).set_index('ModelID')
    ge = pd.read_csv(SOURCE, usecols=lambda c: c.startswith('Unnamed:') or
                     c.rsplit(' (', 1)[0] in ['KDM6A', 'KMT2D'], index_col=0)
    ge.columns = ge.columns.str.replace(r' \(\d+\)$', '', regex=True)
    assert ge.index.is_unique and cohort.index.is_unique
    assert ge.index.tolist() == cohort.index.tolist()
    assert set(ge.columns) == {'KDM6A', 'KMT2D'}
    selected = cohort.loc[cohort.included].join(ge, validate='one_to_one')
    assert selected.Sex.eq('Male').all() and selected.Sex_calls.eq('Male').all()
    selected['paired_measurements'] = selected[['KDM6A','KMT2D']].notna().all(axis=1)
    selected.to_csv(DATA/'model_values.csv', index_label='ModelID')
    values = selected.loc[selected.paired_measurements].copy()
    old = pd.read_csv(PREVIOUS/'loy_top10_partners.csv')
    expected = old.loc[old.seed.eq('KDM6A') & old.partner.eq('KMT2D')].iloc[0]
    stats = []
    for group in ['XY','LOY']:
        subset = values.loc[values.group.eq(group)]
        x, y = subset.KDM6A.to_numpy(), subset.KMT2D.to_numpy()
        r = float(np.corrcoef(x,y)[0,1])
        assert abs(r-expected[f'{group.lower()}_pearson_r']) < 1e-9
        assert len(subset) == expected[f'{group.lower()}_n_models']
        slope, intercept = np.polyfit(x,y,1)
        stats.append(dict(group=group, label='XY' if group=='XY' else 'X0 / LOY',
                          n_models=len(subset), n_patients=subset.PatientID.nunique(),
                          pearson_r=r, slope=float(slope), intercept=float(intercept),
                          missing_pairs=int((selected.group.eq(group) & ~selected.paired_measurements).sum())))
    stats = pd.DataFrame(stats).set_index('group')
    stats.to_csv(DATA/'correlations.csv')

    apply_style(9)
    colours = condition_colours()
    palette = {'XY': colours['Y-retaining'], 'LOY': colours['Y-depleted']}
    width, height = 7.4, 4.0
    fig, axes = plt.subplots(1,2,figsize=(width,height),sharex=True,sharey=True)
    fig.subplots_adjust(left=.10,right=.975,bottom=.20,top=.88,wspace=.12)
    def padded_bounds(series):
        low, high = float(series.min()), float(series.max())
        padding = .07*(high-low)
        return low-padding, high+padding
    xlim, ylim = padded_bounds(values.KDM6A), padded_bounds(values.KMT2D)
    for letter, ax, group in zip(['a','b'],axes,['XY','LOY']):
        subset = values.loc[values.group.eq(group)]
        colour = palette[group]
        ax.axhline(0,color='#d8d8d8',lw=.6,ls=(0,(3,3)),zorder=0)
        ax.axvline(0,color='#d8d8d8',lw=.6,ls=(0,(3,3)),zorder=0)
        dots = ax.scatter(subset.KDM6A,subset.KMT2D,s=16,alpha=.55,c=colour,
                          linewidths=0,clip_on=True,zorder=2)
        line_x = np.array([subset.KDM6A.min(),subset.KDM6A.max()])
        ax.plot(line_x,stats.loc[group,'intercept']+stats.loc[group,'slope']*line_x,
                color=colour,lw=1.25,zorder=3)
        ax.set(xlim=xlim,ylim=ylim,xlabel='KDM6A gene effect (Chronos)')
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.set_title('XY (Y-retaining)' if group=='XY' else 'X0 / LOY',pad=10,color=colour)
        ax.text(-.15,1.05,letter,transform=ax.transAxes,fontsize=12,fontweight='semibold',va='bottom')
        ax.text(.04,.96,f"Pearson r = {stats.loc[group,'pearson_r']:.2f}\nn = {int(stats.loc[group,'n_models'])}",
                transform=ax.transAxes,ha='left',va='top',fontsize=9,
                bbox=dict(boxstyle='square,pad=.2',facecolor='white',alpha=.8,edgecolor='none'),zorder=4)
        assert len(dots.get_offsets()) == len(subset)
        assert np.allclose(dots.get_offsets(),subset[['KDM6A','KMT2D']].to_numpy())
        assert subset.KDM6A.between(*xlim).all() and subset.KMT2D.between(*ylim).all()
    axes[0].set_ylabel('KMT2D gene effect (Chronos)')
    fig.text(.10,.035,'More negative gene effects indicate stronger dependencies.',fontsize=8,color='#555555')
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=180)
    plt.close(fig)

    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in font[3] for font in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside = [s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps = []
        for i,a in enumerate(spans):
            for b in spans[i+1:]:
                ra,rb = fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
                area = (ra&rb).get_area()
                if area > 1 and area/min(ra.get_area(),rb.get_area()) > .1:
                    overlaps.append([a['text'],b['text']])
        assert not outside, outside
        assert not overlaps, overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.8,1.8),colorspace=fitz.csGRAY).save(
            OUT/'review'/f'{STEM}_greyscale.png')
    manifest = dict(stem=STEM,dependency_release=previous['dependency_release'],
                    y_status_release=previous['y_status_release'],
                    cohort='Same eligible male cohorts as the preceding XY/LOY heatmaps.',
                    displayed_model_counts=stats.n_models.to_dict(),
                    interpretation='X0/LOY is the primary Y-negative call; it does not establish a 45,X karyotype. Unadjusted descriptive correlations; no between-group correlation test or donor adjustment.',
                    regression='Ordinary least squares of KMT2D on KDM6A, separately by group, shown only over observed x ranges; no confidence bands.',
                    axes={'x':'KDM6A Chronos gene effect','y':'KMT2D Chronos gene effect',
                          'shared_xlim':xlim,'shared_ylim':ylim,'no_outlier_exclusion':True},
                    dimensions_inches=[width,height],font='Open Sans',colours=palette,
                    validation={'matches_previous_correlations':True,'all_points_plotted':True,
                                'text_outside_page':outside,'text_overlaps':overlaps},
                    source_sha256={str(SOURCE):digest,**{str(p):sha256(p) for p in [cohort_path,
                        PREVIOUS/'loy_top10_partners.csv',PREVIOUS/'analysis_manifest.json']}})
    (DATA/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(stats.to_string())
    print(json.dumps({k:v for k,v in manifest.items() if k!='source_sha256'},indent=2))


if __name__ == '__main__':
    main()
