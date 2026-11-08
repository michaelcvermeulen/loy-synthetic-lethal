"""Render A549 counterparts of the HT1080 main and supplementary QC panels."""
from pathlib import Path
import hashlib
import json
import os
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent / '.figure_python_deps'))
OUT = ROOT / 'figures/chapter4'
DATA = OUT / 'source_data/a549_matched_qc'
REVIEW = OUT / 'review'
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np
import pandas as pd
import pymupdf
from PIL import Image
from make_a549_screen_schematic import draw_schematic
from draw_a549_molecular_confirmation import draw_molecular_confirmation
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
apply_style(9)
COL = condition_colours()
BLUE, RED, DARK = COL['Y-retaining'], COL['Y-depleted'], '#222222'
SUMMARY = json.loads((DATA / 'summary.json').read_text())
SAMPLES = SUMMARY['samples']
LABELS = ['CTRL D0 1', 'CTRL D0 2', 'LOY D0 1', 'LOY D0 2', 'CTRL P4 R2', 'LOY P4 R1', 'LOY P4 R2']
COLOURS = [BLUE, BLUE, RED, RED, BLUE, RED, RED]
LOG = pd.read_csv(DATA / 'guide_log2_cpm.csv', index_col=0)
QC = pd.read_csv(DATA / 'sample_qc.csv').set_index('sample').loc[SAMPLES]

def letter(fig, x, y, value):
    fig.text(x, y, value, fontsize=13, fontweight='semibold', va='top')

def sample_ticks(ax):
    ax.set_xticks(range(7), LABELS, rotation=55, ha='right', fontsize=7.3)
    ax.set_xlim(-.6, 6.6)
    ax.axvline(3.5, color='#BCC2C9', lw=.8, ls='--')

def depletion(ax):
    genes = pd.read_csv(DATA / 'baseline_gene_qc.csv')
    values = [genes.loc[genes.genotype.eq(g) & genes.category.eq(c), 'lfc'].to_numpy()
              for g, c in [('CTRL', 'Essential'), ('CTRL', 'Non-essential'), ('LOY', 'Essential'), ('LOY', 'Non-essential')]]
    positions = [0, 1, 2.4, 3.4]
    v = ax.violinplot(values, positions=positions, widths=.82, showmedians=True, showextrema=False, points=120)
    for body, color in zip(v['bodies'], [BLUE, BLUE, RED, RED]):
        body.set(facecolor=color, edgecolor=color, alpha=.58, lw=.7)
    v['cmedians'].set(color=DARK, linewidth=1.1)
    ax.axhline(0, color='#9AA1A8', ls='--', lw=.8)
    lo, hi = min(map(np.min, values)), max(map(np.max, values))
    ax.set(xlim=(-.65, 4.05), ylim=(np.floor(lo)-.2, np.ceil(hi)+1.2), ylabel='Gene log₂(P4 / D0)')
    ax.set_xticks(positions, ['Ess.', 'Non-ess.', 'Ess.', 'Non-ess.'], fontsize=8)
    for x, values in zip(positions, values):
        ax.text(x, np.ceil(hi)+.2, f'{len(values):,}', ha='center', fontsize=8)
    ax.text(.02, .98, 'n genes', transform=ax.transAxes, va='top', fontsize=8)
    for x, group, color in [(.5, 'CTRL', BLUE), (2.9, 'LOY', RED)]:
        ax.text(x, .065, f'ROC AUC {SUMMARY["roc"][group]["auc"]:.3f}',
                transform=ax.get_xaxis_transform(), ha='center', va='bottom', color=color, fontsize=8)
        ax.text(x, -.20, group, transform=ax.get_xaxis_transform(), ha='center', color=color, fontweight='semibold')

def matrix(ax, fig):
    data = pd.read_csv(DATA / 'spearman_guide_counts.csv', index_col=0).loc[SAMPLES, SAMPLES]
    minimum = np.floor(data.to_numpy().min()*20)/20
    norm = Normalize(vmin=minimum, vmax=1)
    im = ax.imshow(data, cmap='magma', norm=norm, aspect='equal')
    for i in range(7):
        for j in range(7):
            ax.text(j, i, f'{data.iloc[i,j]:.3f}', ha='center', va='center', fontsize=7,
                    color=DARK if norm(data.iloc[i,j]) > .68 else 'white')
    ax.set_xticks(range(7), LABELS, rotation=55, ha='right', fontsize=7)
    ax.set_yticks(range(7), LABELS, fontsize=7)
    for labels in [ax.get_xticklabels(), ax.get_yticklabels()]:
        for lab, color in zip(labels, COLOURS):
            lab.set_color(color)
    ax.tick_params(length=0, pad=3)
    ax.axhline(3.5, color='white', lw=1)
    ax.axvline(3.5, color='white', lw=1)
    for sp in ax.spines.values(): sp.set_visible(False)
    cb = fig.colorbar(im, cax=ax.inset_axes([0, 1.075, 1, .045]), orientation='horizontal')
    cb.ax.tick_params(labelsize=8, length=2, pad=2)
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
    cb.set_label('Spearman correlation', fontsize=9, labelpad=4)

def pca(ax):
    d = pd.read_csv(DATA / 'pca_sample_scores.csv').set_index('sample').loc[SAMPLES]
    offsets = [(5, 12), (5, -13), (5, 12), (5, -13), (5, 10), (5, 13), (5, -13)]
    for sample, label, color, offset in zip(SAMPLES, LABELS, COLOURS, offsets):
        r = d.loc[sample]
        ax.scatter(r.PC1, r.PC2, s=30, marker='D' if r.stage == 'D0' else 'o', color=color, edgecolor='white', lw=.5, zorder=3)
        ax.annotate(label, (r.PC1, r.PC2), xytext=offset, textcoords='offset points', color=color, fontsize=7,
                    va='center', arrowprops={'arrowstyle': '-', 'color': color, 'lw': .5, 'alpha': .55})
    pct = SUMMARY['pca']['variance_percent']
    dx, dy = np.ptp(d.PC1), np.ptp(d.PC2)
    ax.set(xlabel=f'PC1 ({pct[0]:.1f}%)', ylabel=f'PC2 ({pct[1]:.1f}%)',
           xlim=(d.PC1.min()-.12*dx, d.PC1.max()+.40*dx), ylim=(d.PC2.min()-.22*dy, d.PC2.max()+.22*dy))

def densities(ax):
    table = pd.read_csv(DATA / 'endpoint_reference_density_guides.csv')
    curves = pd.read_csv(DATA / 'endpoint_reference_density_curves.csv')
    meta = SUMMARY['density']['groups']
    ntc = table.loc[table.included & table.category.eq('Non-targeting'), 'LFC']
    for cat, color, ls in [('Essential', '#8473AC', '-'), ('Non-essential', '#B79549', '--')]:
        ax.fill_between(curves.LFC, curves[cat], color=color, alpha=.14, zorder=1, lw=0)
        ax.plot(curves.LFC, curves[cat], color=color, alpha=.55, lw=1.1, ls=ls, zorder=2,
                label=f'{cat} (n = {meta[cat]["included_guides"]:,})')
    edges = np.arange(np.floor(ntc.min()*10)/10, np.ceil(ntc.max()*10)/10+.1001, .1)
    ax.hist(ntc, bins=edges, density=True, color='#858E9D', alpha=.22, edgecolor='white', lw=.25, zorder=3)
    ax.plot(curves.LFC, curves['Non-targeting'], color=DARK, lw=1.5, zorder=4, label=f'Non-targeting (n = {len(ntc):,})')
    ax.axvline(0, color='#9AA1A8', ls=':', lw=.7, zorder=0)
    ax.set(xlabel='Guide log₂ fold change (LOY / CTRL P4)', ylabel='Density',
           xlim=(curves.LFC.min(), curves.LFC.max()), ylim=(0, curves.drop(columns='LFC').to_numpy().max()*1.32))
    handles, labels = ax.get_legend_handles_labels()
    ax.legend([handles[i] for i in [2,0,1]], [labels[i] for i in [2,0,1]], loc='upper left', fontsize=7, handlelength=1.8, labelspacing=.3)
    s = meta['Non-targeting']
    ax.text(.985, .96, f'NTC median = {s["median"]:.3f}\nSD = {s["sd"]:.3f}\nCentral 95%: {s["q025"]:.2f} to {s["q975"]:.2f}',
            ha='right', va='top', transform=ax.transAxes, fontsize=7, linespacing=1.15)

def main_figure():
    width, height = 243.84, 224
    fig = plt.figure(figsize=(width/25.4, height/25.4))
    def bounds(x, y, w, h): return [x/width, y/height, w/width, h/height]
    # Place the Chapter 3 inset beside the two middle-row QC panels. All
    # quantitative screen inputs and panel letters b-e remain unchanged.
    draw_schematic(fig, bounds(10, 156, 224, 64))
    draw_molecular_confirmation(fig, bounds(18, 96, 78, 50))
    depletion(fig.add_axes(bounds(111, 88, 47, 48)))
    matrix(fig.add_axes(bounds(183, 85, 52, 52)), fig)
    pca(fig.add_axes(bounds(18, 14, 87, 48)))
    densities(fig.add_axes(bounds(130, 14, 106, 48)))
    for x, y, label in [(3,222,'a'),(103,149,'b'),(166,149,'c'),(3,68,'d'),(116,68,'e')]:
        letter(fig,x/width,y/height,label)
    return fig


def mapping(fig):
    a = fig.add_axes([.085,.835,.225,.10])
    b = fig.add_axes([.085,.735,.225,.072])
    a.bar(range(7), QC.mapping_percent, color=COLOURS, width=.65, alpha=.9)
    a.set(ylim=(0,105), ylabel='Mapped reads (%)', xlim=(-.6,6.6))
    a.set_yticks([0,50,100]); a.set_xticks([])
    a.axvline(3.5,color='#BCC2C9',lw=.8,ls='--')
    for x,v in enumerate(QC.mapping_percent): a.text(x,v+2,f'{v:.1f}',ha='center',fontsize=6.8)
    b.scatter(range(7), QC.detected_percent, c=COLOURS,s=22,zorder=3)
    b.set(ylim=(98.8,100),ylabel='Guides\ndetected (%)')
    b.set_yticks([98.8,99.4,100]); b.yaxis.label.set_size(8)
    sample_ticks(b)
    for x,v in enumerate(QC.detected_percent): b.text(x,v+.06,f'{v:.2f}',ha='center',fontsize=6.4)

def abundance(ax):
    violin=ax.violinplot([LOG[s].to_numpy() for s in SAMPLES],positions=range(7),widths=.78,showmedians=True,showextrema=False,points=120)
    for body,col in zip(violin['bodies'],COLOURS): body.set(facecolor=col,edgecolor=col,alpha=.58,lw=.7)
    violin['cmedians'].set(color=DARK,lw=1)
    sample_ticks(ax)
    ax.set(ylabel='Guide abundance, log₂(CPM + 1)',ylim=(-.25,10))
    for x,s in enumerate(SAMPLES): ax.text(x,1.015,str(QC.loc[s,'zero_guides']),ha='center',fontsize=7,transform=ax.get_xaxis_transform())
    ax.text(.02,1.09,'Zero-count guides',transform=ax.transAxes,va='bottom',fontsize=8)
    ax.text(.98,.92,'77,441 guides / sample',transform=ax.transAxes,ha='right',fontsize=7.5)

def gini(ax):
    d=QC.gini_mageck_export
    ax.bar(range(7),d,color=COLOURS,width=.65,alpha=.9)
    ax.set(ylim=(0,.14),ylabel='MAGeCK Gini index (log counts)')
    ax.set_yticks([0,.04,.08,.12]); sample_ticks(ax)
    for x,v in enumerate(d): ax.text(x,v+.003,f'{v:.5f}',ha='center',va='bottom',fontsize=7,rotation=55)

def scatter(ax,key):
    if key=='D0_pooled':
        d=pd.read_csv(DATA/'pooled_d0_log2_cpm.csv'); x,y=d.CTRL,d.LOY
        title,color='D0: pooled CTRL vs LOY',DARK
        xlabel,ylabel='CTRL D0','LOY D0'
    else:
        pairs={'CTRL_D0_technical':('CTRL_D0_1','CTRL_D0_2','CTRL D0: technical libraries',BLUE,'Library 1','Library 2'),
               'LOY_D0_technical':('LOY_D0_1','LOY_D0_2','LOY D0: technical libraries',RED,'Library 1','Library 2'),
               'LOY_P4':('LOY_P4_R1','LOY_P4_R2','LOY P4: outgrowth cultures',RED,'R1','R2')}
        xs,ys,title,color,xlabel,ylabel=pairs[key]; x,y=LOG[xs],LOG[ys]
    limit=9 if 'D0' in key else 8.25
    assert max(x.max(),y.max())<limit
    norm=LogNorm(vmin=1,vmax=10000 if 'D0' in key else 3000)
    hb=ax.hexbin(x,y,gridsize=95,extent=(0,limit,0,limit),mincnt=1,cmap='magma',norm=norm,lw=0,rasterized=True)
    assert hb.get_array().max() <= norm.vmax
    s=SUMMARY['regression'][key]
    ax.plot([0,limit],[0,limit],color='#9AA1A8',ls='--',lw=.9,label='Identity')
    ax.plot([0,limit],[s['intercept'],s['intercept']+limit*s['slope']],color=color,lw=1.35,label='Linear fit')
    ax.set(xlim=(0,limit),ylim=(0,limit),xlabel=xlabel+', log₂(CPM + 1)',ylabel=ylabel+', log₂(CPM + 1)')
    ticks=[0,3,6,9] if 'D0' in key else [0,2,4,6,8]
    ax.set_xticks(ticks);ax.set_yticks(ticks)
    ax.text(.025,.975,title,transform=ax.transAxes,va='top',fontweight='semibold',color=color,fontsize=7.6)
    ax.text(.025,.86,f'Pearson r = {s["pearson_r"]:.3f}\nn = 77,441 guides',transform=ax.transAxes,va='top',fontsize=8)
    ax.legend(loc='lower right',fontsize=7,handlelength=2)
    return hb

def roc(ax):
    for group,color in [('CTRL',BLUE),('LOY',RED)]:
        d=pd.read_csv(DATA/f'roc_curve_{group}.csv')
        ax.plot(d.false_positive_rate,d.true_positive_rate,color=color,lw=1.6,label=f'{group}  AUC {SUMMARY["roc"][group]["auc"]:.3f}')
    ax.plot([0,1],[0,1],color='#B4B8BF',lw=.8,ls='--',zorder=0)
    ax.set(xlim=(0,1),ylim=(0,1.02),xlabel='False-positive rate',ylabel='True-positive rate')
    ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1]);ax.legend(loc='lower right',fontsize=8)
    ax.text(.97,.39,'Essential / non-essential genes\n1,206 / 749 in each genotype',ha='right',transform=ax.transAxes,fontsize=7.5)

def concentration(ax):
    values=QC.top100_read_fraction*100
    ax.bar(range(7),values,color=COLOURS,width=.65,alpha=.9)
    ax.set(ylim=(0,1.2),ylabel='Reads in top 100 guides (%)')
    sample_ticks(ax)
    for x,v in enumerate(values): ax.text(x,v+.025,f'{v:.2f}',ha='center',fontsize=7)

def supplement():
    fig=plt.figure(figsize=(9.6,9.8))
    mapping(fig)
    abundance(fig.add_axes([.405,.735,.235,.20]))
    gini(fig.add_axes([.735,.735,.235,.20]))
    hbs=[]
    for x,key in zip([.075,.405,.735],['D0_pooled','CTRL_D0_technical','LOY_D0_technical']):
        hbs.append(scatter(fig.add_axes([x,.395,.235,.21]),key))
    cb=fig.colorbar(hbs[0],cax=fig.add_axes([.075,.645,.895,.008]),orientation='horizontal',ticks=[1,10,100,1000,10000])
    cb.set_label('D0 guides / hexagon',fontsize=8,labelpad=3);cb.ax.xaxis.set_label_position('top');cb.ax.tick_params(labelsize=7,length=2,pad=2)
    hb=scatter(fig.add_axes([.075,.075,.235,.21]),'LOY_P4')
    cb=fig.colorbar(hb,cax=fig.add_axes([.075,.325,.235,.008]),orientation='horizontal',ticks=[1,10,100,1000])
    cb.set_label('P4 guides / hexagon',fontsize=8,labelpad=3);cb.ax.xaxis.set_label_position('top');cb.ax.tick_params(labelsize=7,length=2,pad=2)
    roc(fig.add_axes([.405,.075,.235,.21]))
    concentration(fig.add_axes([.735,.075,.235,.21]))
    for y,letters in [(.97,'abc'),(.63,'def'),(.31,'ghi')]:
        for x,l in zip([.02,.35,.68],letters):letter(fig,x,y,l)
    return fig

def save(fig,stem):
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
    preview=REVIEW/f'{stem}_preview.png'
    fig.savefig(preview,dpi=160)
    image=Image.open(preview); image.convert('L').save(REVIEW/f'{stem}_greyscale.png')
    size=(fig.get_size_inches()*25.4).tolist();plt.close(fig)
    with pymupdf.open(OUT/f'{stem}.pdf') as doc:
        page=doc[0]
        outside=[w for w in page.get_text('words') if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
        assert not outside,outside
        has_workflow=stem.startswith('ch04_main_')
        assert ('P3 collection' in page.get_text()) if has_workflow else ('P3' not in page.get_text())
        assert all('OpenSans' in f[3].replace(' ','') for f in page.get_fonts(full=True))
    assert '<text' in (OUT/f'{stem}.svg').read_text(encoding='utf-8')
    with Image.open(OUT/f'{stem}.png') as im: assert min(im.info['dpi'])>599
    return {'stem':stem,'figure_mm':size,'pdf_text_in_bounds':True,'open_sans':True,'editable_svg':True,'png_dpi':600,
            'p3_shown_in_workflow_only':has_workflow,'qc_uses_d0_p4_only':True}

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--main-only',action='store_true')
    args=parser.parse_args()
    hashes={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob('*.csv')}
    results=[save(main_figure(),'ch04_main_a549_crispr_screen_qc_draft')]
    if not args.main_only:results.append(save(supplement(),'ch04_supp_a549_crispr_screen_qc_draft'))
    elif (DATA/'figure_validation.json').exists():
        results.extend(r for r in json.loads((DATA/'figure_validation.json').read_text())['exports'] if r['stem'].startswith('ch04_supp_'))
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==h for p,h in hashes.items())
    (DATA/'figure_validation.json').write_text(json.dumps({'exports':results,'source_tables_unchanged':True},indent=2))
    print(json.dumps(results,indent=2))

if __name__=='__main__':
    main()
