"""DDX3X versus IGBP1, GBF1 and IK dependencies in the saved XY/LOY cohorts."""
from pathlib import Path
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
DATA = OUT/'source_data/ddx3x_partner_y_status_scatter'
SOURCE = Path('A:/DEPMAP/DepMap/CRISPRGeneEffect.csv')
STEM = 'ch04_supp_ddx3x_igbp1_gbf1_ik_xy_loy_scatter_draft'
PARTNERS = ['IGBP1','GBF1','IK']
sys.path.insert(0, str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
from prepare_xy_ranked_partners_by_y_status import sha256


def padded_bounds(series):
    low, high = float(series.min()), float(series.max())
    pad = .07*(high-low)
    return low-pad, high+pad


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    previous = json.loads((PREVIOUS/'analysis_manifest.json').read_text())
    digest = sha256(SOURCE)
    assert digest == previous['sources'][str(SOURCE)]
    cohort_path = PREVIOUS/'model_cohort_and_exclusions.csv'
    cohort = pd.read_csv(cohort_path).set_index('ModelID')
    ge = pd.read_csv(SOURCE, usecols=lambda c: c.startswith('Unnamed:') or
                     c.rsplit(' (', 1)[0] in ['DDX3X']+PARTNERS, index_col=0)
    ge.columns = ge.columns.str.replace(r' \(\d+\)$', '', regex=True)
    assert ge.index.is_unique and cohort.index.is_unique
    assert ge.index.tolist() == cohort.index.tolist()
    assert set(ge.columns) == set(['DDX3X']+PARTNERS)
    selected = cohort.loc[cohort.included].join(ge, validate='one_to_one')
    assert selected.Sex.eq('Male').all() and selected.Sex_calls.eq('Male').all()
    for partner in PARTNERS:
        selected[f'{partner}_paired_measurements'] = selected[['DDX3X',partner]].notna().all(axis=1)
    selected.to_csv(DATA/'model_values.csv', index_label='ModelID')
    old = pd.read_csv(PREVIOUS/'loy_top10_partners.csv')
    old = old.loc[old.seed.eq('DDX3X') & old.partner.isin(PARTNERS)].set_index('partner')
    stats, subsets = [], {}
    for partner in PARTNERS:
        for group in ['XY','LOY']:
            subset = selected.loc[selected.group.eq(group) & selected[f'{partner}_paired_measurements']]
            subsets[partner,group] = subset
            x, y = subset.DDX3X.to_numpy(), subset[partner].to_numpy()
            r = float(np.corrcoef(x,y)[0,1])
            assert abs(r-old.loc[partner,f'{group.lower()}_pearson_r']) < 1e-9
            assert len(subset) == old.loc[partner,f'{group.lower()}_n_models']
            slope, intercept = np.polyfit(x,y,1)
            stats.append(dict(partner=partner, group=group, n_models=len(subset),
                              n_patients=subset.PatientID.nunique(), pearson_r=r,
                              slope=float(slope), intercept=float(intercept),
                              missing_pairs=int((selected.group.eq(group) & ~selected[f'{partner}_paired_measurements']).sum())))
    stats = pd.DataFrame(stats).set_index(['partner','group'])
    stats.to_csv(DATA/'correlations.csv')

    apply_style(9)
    colours = condition_colours()
    palette = {'XY':colours['Y-retaining'], 'LOY':colours['Y-depleted']}
    width,height = 7.4,8.5
    fig,axes = plt.subplots(3,2,figsize=(width,height),sharex=True,sharey='row')
    fig.subplots_adjust(left=.105,right=.975,bottom=.095,top=.94,wspace=.13,hspace=.20)
    xlim = padded_bounds(selected.DDX3X.dropna())
    ylimits = {}
    for row,partner in enumerate(PARTNERS):
        ylim = padded_bounds(selected.loc[selected[f'{partner}_paired_measurements'],partner])
        ylimits[partner] = ylim
        for col,group in enumerate(['XY','LOY']):
            ax,subset = axes[row,col],subsets[partner,group]
            stat,colour = stats.loc[partner,group],palette[group]
            ax.axhline(0,color='#d8d8d8',lw=.6,ls=(0,(3,3)),zorder=0)
            ax.axvline(0,color='#d8d8d8',lw=.6,ls=(0,(3,3)),zorder=0)
            dots = ax.scatter(subset.DDX3X,subset[partner],s=13,alpha=.55,c=colour,
                              linewidths=0,clip_on=True,zorder=2)
            line_x = np.array([subset.DDX3X.min(),subset.DDX3X.max()])
            ax.plot(line_x,stat.intercept+stat.slope*line_x,color=colour,lw=1.2,zorder=3)
            ax.set(xlim=xlim,ylim=ylim)
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            if row == 0:
                ax.set_title('XY (Y-retaining)' if group=='XY' else 'LOY (Y-depleted)',pad=12,color=colour)
            ax.text(-.155,1.025,chr(97+row*2+col),transform=ax.transAxes,
                    fontsize=12,fontweight='semibold',va='bottom')
            ax.text(.04,.96,f'Pearson r = {stat.pearson_r:.2f}\nn = {int(stat.n_models)}',
                    transform=ax.transAxes,ha='left',va='top',fontsize=8.5,
                    bbox=dict(boxstyle='square,pad=.2',facecolor='white',alpha=.8,edgecolor='none'),zorder=4)
            assert len(dots.get_offsets()) == len(subset)
            assert np.allclose(dots.get_offsets(),subset[['DDX3X',partner]].to_numpy())
            assert subset.DDX3X.between(*xlim).all() and subset[partner].between(*ylim).all()
        axes[row,0].set_ylabel(f'{partner} gene effect (Chronos)')
    for ax in axes[-1]:
        ax.set_xlabel('DDX3X gene effect (Chronos)')
    fig.text(.105,.02,'More negative gene effects indicate stronger dependencies.',fontsize=8,color='#555555')
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=170)
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
                if area>1 and area/min(ra.get_area(),rb.get_area())>.1:
                    overlaps.append([a['text'],b['text']])
        assert not outside, outside
        assert not overlaps, overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.5,1.5),colorspace=fitz.csGRAY).save(
            OUT/'review'/f'{STEM}_greyscale.png')
    manifest = dict(stem=STEM,partners=PARTNERS,dependency_release=previous['dependency_release'],
                    y_status_release=previous['y_status_release'],
                    cohort='Same eligible male cohorts as preceding XY/LOY heatmaps; no new selection of models.',
                    selection='User-selected IGBP1, GBF1 and IK; the first three partners in the preceding within-LOY DDX3X ranking.',
                    interpretation='Descriptive model-level correlations. Partners selected after the LOY ranking; no difference-in-correlation test or adjustment for lineage or donor relationships.',
                    regression='Ordinary least squares of partner gene effect on DDX3X, separately by group, over observed x ranges; no confidence bands.',
                    axes={'x':'DDX3X Chronos gene effect','shared_xlim':xlim,'shared_y_limits_within_partner':ylimits,
                          'no_outlier_exclusion':True},
                    dimensions_inches=[width,height],font='Open Sans',colours=palette,
                    validation={'matches_all_six_previous_correlations':True,'all_points_plotted':True,
                                'text_outside_page':outside,'text_overlaps':overlaps},
                    source_sha256={str(SOURCE):digest,**{str(p):sha256(p) for p in [cohort_path,
                        PREVIOUS/'loy_top10_partners.csv',PREVIOUS/'analysis_manifest.json']}})
    (DATA/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(stats.to_string())
    print(json.dumps(manifest['validation'],indent=2))


if __name__ == '__main__':
    main()
