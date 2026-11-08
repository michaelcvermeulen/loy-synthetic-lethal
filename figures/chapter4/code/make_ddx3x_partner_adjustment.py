"""Raw and adjusted DDX3X partner scatters, with nuisance-check intervals."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS/'figures/chapter4'
DATA = OUT/'source_data/ddx3x_partner_adjustment'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style

PARTNERS = ['IGBP1','GBF1','IK']
GROUPS = ['XY','LOY','Female']
PALETTE = {'XY':'#3C5488','LOY':'#E64B35','Female':'#00A087'}
MARKERS = {'XY':'o','LOY':'s','Female':'^'}
LABELS = {'XY':'XY (Y-retaining)','LOY':'LOY (Y-depleted)','Female':'Female'}
ANALYSES = ['raw','lineage','lineage_efficacy','lineage_essential_mean']
METHOD_LABELS = ['Raw','Lineage','Lineage +\nefficacy','Lineage +\nessential mean']


def bounds(s):
    lo,hi = float(s.min()),float(s.max())
    pad = .075*(hi-lo)
    return [lo-pad,hi+pad]


def save_and_review(fig,stem):
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{stem}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{stem}_preview.png',dpi=150)
    plt.close(fig)
    with fitz.open(OUT/f'{stem}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside = [s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps = []
        for i,a in enumerate(spans):
            for b in spans[i+1:]:
                ra,rb = fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
                overlap = (ra&rb).get_area()
                if overlap > 1 and overlap/min(ra.get_area(),rb.get_area()) > .1:
                    overlaps.append([a['text'],b['text']])
        assert not outside,outside
        assert not overlaps,overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.5,1.5),colorspace=fitz.csGRAY).save(
            OUT/'review'/f'{stem}_greyscale.png')
    return {'stem':stem,'font':'Open Sans','text_outside_page':outside,'text_overlaps':overlaps}


def scatter_panels(fig,axes,values,stats,adjusted):
    xlim = bounds(values.DDX3X)
    audits = []
    for row,partner in enumerate(PARTNERS):
        ylim = bounds(values[partner])
        for col,group in enumerate(GROUPS):
            ax = axes[row,col]
            frame = values.loc[values.group.eq(group)]
            st = stats.loc[(partner,group)]
            ax.axhline(0,color='#d8d8d8',lw=.6,ls=(0,(3,3)),zorder=0)
            ax.axvline(0,color='#d8d8d8',lw=.6,ls=(0,(3,3)),zorder=0)
            dots = ax.scatter(frame.DDX3X,frame[partner],s=11,alpha=.5,
                              color=PALETTE[group],linewidths=0,zorder=2)
            slope,intercept = np.polyfit(frame.DDX3X,frame[partner],1)
            xs = np.array([frame.DDX3X.min(),frame.DDX3X.max()])
            ax.plot(xs,intercept+slope*xs,lw=1.15,color=PALETTE[group],zorder=3)
            ax.set(xlim=xlim,ylim=ylim)
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            if row == 0:
                ax.set_title(LABELS[group],pad=12,color=PALETTE[group],fontsize=10)
            ax.text(-.17,1.025,chr(97+row*3+col),transform=ax.transAxes,
                    fontsize=12,fontweight='semibold',va='bottom')
            coefficient = 'Partial r' if adjusted else 'Pearson r'
            ax.text(.035,.965,f'{coefficient} = {st.pearson_r:.2f}\nn = {len(frame)}',
                    transform=ax.transAxes,ha='left',va='top',fontsize=8,
                    bbox=dict(boxstyle='square,pad=.18',facecolor='white',alpha=.85,edgecolor='none'),zorder=4)
            if col == 0:
                ax.set_ylabel(f'{partner} '+('residual gene effect' if adjusted else 'gene effect (Chronos)'))
            if row == 2:
                ax.set_xlabel('DDX3X '+('residual gene effect' if adjusted else 'gene effect (Chronos)'))
            assert len(dots.get_offsets()) == len(frame)
            assert np.allclose(dots.get_offsets(),frame[['DDX3X',partner]].to_numpy())
            assert abs(np.corrcoef(frame.DDX3X,frame[partner])[0,1]-st.pearson_r) < 1e-12
            assert frame.DDX3X.between(*xlim).all() and frame[partner].between(*ylim).all()
            audits.append(dict(partner=partner,group=group,n_points=len(frame),xlim=xlim,ylim=ylim))
    return audits


def main():
    apply_style(9)
    values = pd.read_csv(DATA/'model_values.csv')
    residuals = pd.read_csv(DATA/'residuals.csv')
    stats = pd.read_csv(DATA/'correlations.csv')
    export_records = []
    raw_stem = 'ch04_supp_ddx3x_partners_xy_loy_female_raw_draft'
    fig,axes = plt.subplots(3,3,figsize=(9.6,8.4),sharex=True,sharey='row')
    fig.subplots_adjust(left=.085,right=.985,bottom=.075,top=.94,wspace=.17,hspace=.20)
    audit = scatter_panels(fig,axes,values,stats[stats.analysis.eq('raw')].set_index(['partner','group']),False)
    record = save_and_review(fig,raw_stem)
    record.update(dimensions_inches=[9.6,8.4],points=audit)
    export_records.append(record)

    adjusted_stem = 'ch04_supp_ddx3x_partners_xy_loy_female_adjusted_draft'
    fig = plt.figure(figsize=(10.2,11.8))
    grid = fig.add_gridspec(4,3,left=.09,right=.985,bottom=.115,top=.95,
                           height_ratios=[1,1,1,.95],wspace=.22,hspace=.50)
    axes = np.empty((3,3),dtype=object)
    for row in range(3):
        for col in range(3):
            axes[row,col] = fig.add_subplot(grid[row,col],sharex=axes[0,0] if row+col else None,
                                           sharey=axes[row,0] if col else None)
            if row < 2:
                axes[row,col].tick_params(labelbottom=False)
            if col:
                axes[row,col].tick_params(labelleft=False)
    audit = scatter_panels(fig,axes,residuals[residuals.analysis.eq('lineage_efficacy')],
                           stats[stats.analysis.eq('lineage_efficacy')].set_index(['partner','group']),True)
    method_stats = stats.set_index(['partner','group','analysis'])
    xmin = min(-.1,float(stats.bootstrap_ci_low.min())-.04)
    xmax = max(.6,float(stats.bootstrap_ci_high.max())+.04)
    summary_checks = []
    for col,partner in enumerate(PARTNERS):
        ax = fig.add_subplot(grid[3,col])
        ax.axvline(0,lw=.7,color='#c7c7c7',ls=(0,(3,3)),zorder=0)
        for m,analysis in enumerate(ANALYSES):
            for j,group in enumerate(GROUPS):
                st = method_stats.loc[(partner,group,analysis)]
                y = 3-m+(1-j)*.19
                ax.plot([st.bootstrap_ci_low,st.bootstrap_ci_high],[y,y],color=PALETTE[group],lw=1.1)
                ax.scatter([st.pearson_r],[y],color=PALETTE[group],s=22,marker=MARKERS[group],linewidths=0,zorder=3)
                summary_checks.append(dict(partner=partner,group=group,analysis=analysis,
                                           r=st.pearson_r,ci_low=st.bootstrap_ci_low,ci_high=st.bootstrap_ci_high))
        ax.set(xlim=(xmin,xmax),ylim=(-.5,3.5),yticks=[3,2,1,0],yticklabels=METHOD_LABELS,
               xlabel='Pearson / partial r')
        ax.tick_params(axis='y',labelsize=7.5,length=0,pad=3)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.spines['left'].set_visible(False)
        ax.set_title(partner,fontsize=9.5,pad=8)
        ax.text(-.20,1.025,chr(106+col),transform=ax.transAxes,
                fontsize=12,fontweight='semibold',va='bottom')
    handles = [Line2D([],[],linestyle='none',marker=MARKERS[g],color=PALETTE[g],label=LABELS[g],markersize=5) for g in GROUPS]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.55,.035),ncol=3,frameon=False,fontsize=8)
    fig.text(.55,.016,'Intervals: 95% donor-bootstrap estimates; selected partners retained.',ha='center',fontsize=8,color='#555555')
    record = save_and_review(fig,adjusted_stem)
    record.update(dimensions_inches=[10.2,11.8],points=audit,summary_points=summary_checks)
    export_records.append(record)
    (DATA/'render_manifest.json').write_text(json.dumps({'colours':PALETTE,'exports':export_records},indent=2),encoding='utf-8')
    print(json.dumps([{k:v for k,v in r.items() if k not in ['points','summary_points']} for r in export_records],indent=2))


if __name__ == '__main__':
    main()
