"""Supp Fig 4.7: display-only refinement of frozen molecular-integration tables."""
from pathlib import Path
import hashlib
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent / '.figure_python_deps'))
OUT = ROOT / 'figures/chapter4'
SOURCE = OUT / 'source_data/ch4_molecular_integration_20261021'
DATA = OUT / 'source_data/ht1080_molecular_integration_supp_20261028'
DATA.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
import pandas as pd
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)

STEM = 'ch04_supp_ht1080_molecular_integration'
SCREEN = LinearSegmentedColormap.from_list('screen_teal_purple', ['#218C8D', '#F7F7F7', '#7A3E85'])
MOLECULAR = LinearSegmentedColormap.from_list('molecular_blue_red', ['#2166AC', '#F7F7F7', '#B2182B'])
WIDTH, HEIGHT = 9.6, 10.5
CELL_W = .79
COL_X = [3.65, 4.44, 5.35, 6.14, 7.05, 7.84, 8.63]
COL_LABELS = ['drugZ', 'IsogenicZ', 'RNA', 'ATAC', 'UTY', 'KDM5D', 'UTY+\nKDM5D']
A_COLUMNS = ['Screen drugZ', 'Screen IsogenicZ', 'RNA baseline', 'ATAC promoter',
             'Add-back UTY', 'Add-back KDM5D', 'Add-back UTY+KDM5D']
A_GROUPS = [
    ('Screen sensitizer complexes', [
        ('SAGA', ['SAGA complex']),
        ('SWI/SNF', ['SWI/SNF complex', 'ATP-dependent remodelers'])]),
    ('Downstream regulated target programmes', [
        ('Sterol', ['SREBF target activation', 'Cholesterol homeostasis']),
        ('MYC', ['MYC targets V2', 'MYC targets V1']),
        ('Prolif.', ['G2M checkpoint', 'E2F targets'])]),
    ('Screen resistor pathways', [
        ('Mito. transl.', ['Mitochondrial translation', 'Mito. large ribosomal subunit']),
        ('Resp. chain', ['Oxidative phosphorylation', 'Complex I biogenesis']),
        ('Cyto. transl.', ['Translation', 'Ribosome'])])]
B_GROUPS = [
    ('SAGA', ['TAF6L', 'TADA2B', 'SUPT20H', 'ATXN7L3', 'TADA1', 'TAF5L']),
    ('SWI/SNF', ['SMARCA4', 'SMARCE1', 'ARID1A']),
    ('Co-act.', ['MED12', 'KMT2D', 'CREBBP']),
    ('Sterol', ['SREBF1', 'SREBF2']),
    ('Other', ['DDX3X', 'TIMM8A', 'MYC'])]


def hashes():
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in SOURCE.iterdir() if p.is_file()}


def column_headers(ax, y):
    for text, start, end in [('Screen', 0, 1), ('Baseline', 2, 3), ('Add-back vs EV', 4, 6)]:
        ax.text((COL_X[start] + COL_X[end] + CELL_W) / 2, y, text,
                ha='center', va='center', fontsize=9, fontweight='semibold')


def column_ticks(ax, y):
    for x, label in zip(COL_X, COL_LABELS):
        ax.text(x + CELL_W/2, y, label, ha='center', va='top', fontsize=8)


def cell(ax, x, y, height, value, label, bold, cmap, limit):
    colour = '#E6E6E6' if pd.isna(value) else cmap(Normalize(-limit, limit, clip=True)(value))
    ax.add_patch(Rectangle((x, y), CELL_W, height, facecolor=colour, edgecolor='white', linewidth=.65))
    # Keep labels readable in the saturated tails without changing significance encoding.
    colour_rgb = matplotlib.colors.to_rgb(colour)
    linear = [c/12.92 if c <= .04045 else ((c+.055)/1.055)**2.4 for c in colour_rgb]
    luminance = .2126*linear[0] + .7152*linear[1] + .0722*linear[2]
    ink = 'white' if luminance < .19 else '#222222'
    label = str(label).replace('×10', 'e').translate(str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻', '0123456789−'))
    ax.text(x + CELL_W/2, y + height/2, label, ha='center', va='center',
            fontsize=7.6, fontweight='semibold' if bold else 'normal', color=ink)


def colourbars(fig, ax, top, panel):
    for x, cmap, limit, title, unit in [
        (3.65, SCREEN, 8, 'Screen: resistor <-> sensitizer',
         'Signed −log₁₀ FDR' if panel == 'a' else '−normZ'),
        (6.85, MOLECULAR, 8 if panel == 'a' else 1.5, 'RNA / ATAC: lower <-> higher',
         'Signed −log₁₀ FDR' if panel == 'a' else 'log₂ fold change')]:
        width = 2.55
        ax.text(x+width/2, top, title, ha='center', va='bottom', fontsize=8)
        cax = fig.add_axes([x/WIDTH, (HEIGHT-top-.17)/HEIGHT, width/WIDTH, .075/HEIGHT])
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=Normalize(-limit, limit), cmap=cmap),
                          cax=cax, orientation='horizontal', ticks=[-limit, 0, limit])
        cb.outline.set_visible(False)
        cb.ax.tick_params(length=2, pad=2, labelsize=7)
        cb.set_ticklabels([f'≤−{limit:g}', '0', f'≥{limit:g}'])
        ax.text(x+width/2, top+.44, unit, ha='center', va='top', fontsize=7.5)


def text_layout(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = list(fig.texts)
    for ax in fig.axes:
        texts.extend(ax.texts)
        if ax.axison:
            texts.extend([ax.title, ax.xaxis.label, ax.yaxis.label])
            texts.extend(ax.get_xticklabels() + ax.get_yticklabels())
    texts = [t for t in texts if t.get_visible() and t.get_text().strip()]
    boxes = [(t.get_text(), Text.get_window_extent(t, renderer=renderer)) for t in texts]
    overlaps = []
    for i, (text, box) in enumerate(boxes):
        for other, obox in boxes[i+1:]:
            if min(box.x1, obox.x1)-max(box.x0, obox.x0)>.4 and min(box.y1, obox.y1)-max(box.y0, obox.y0)>.4:
                overlaps.append([text, other])
    outside = [text for text,b in boxes if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5]
    report = {'text_items':len(texts), 'overlaps':overlaps, 'outside_page':outside}
    (DATA/'text_layout.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    assert not overlaps and not outside, report
    return report


def main():
    original = hashes()
    a = pd.read_csv(SOURCE/'figure_panel_a_values.csv')
    a = a.loc[a.column.isin(A_COLUMNS)].copy()
    b = pd.read_csv(SOURCE/'figure_panel_b_values.csv')
    alook = a.set_index(['set','column'])
    blook = b.set_index(['symbol','x'])
    assert len(a)==15*7 and len(b)==17*7
    fig = plt.figure(figsize=(WIDTH,HEIGHT))
    ax = fig.add_axes([0,0,1,1]);ax.set_xlim(0,WIDTH);ax.set_ylim(HEIGHT,0);ax.axis('off')
    ax.text(.12,.12,'a',fontsize=14,fontweight='semibold',va='top')
    column_headers(ax,.26)
    row_h = .18
    y = .48
    for category, themes in A_GROUPS:
        ax.plot([.20,9.42],[y-.055,y-.055],color='#D5D5D5',lw=.6,zorder=0)
        ax.text(.22,y,category,fontsize=8.5,fontweight='semibold',va='top')
        y += .22
        for theme, sets in themes:
            ax.text(.22,y+row_h*len(sets)/2,theme,fontsize=7.8,va='center',fontweight='semibold')
            for name in sets:
                ax.text(3.50,y+row_h/2,name,fontsize=7.8,ha='right',va='center')
                for j,column in enumerate(A_COLUMNS):
                    row=alook.loc[(name,column)]
                    cell(ax,COL_X[j],y,row_h,row.signed_log10_FDR,row.lab,row.face=='bold',SCREEN if j<2 else MOLECULAR,8)
                y += row_h
            y += .045
        y += .09
    column_ticks(ax,y-.08)
    colourbars(fig,ax,4.90,'a')

    ax.text(.12,5.65,'b',fontsize=14,fontweight='semibold',va='top')
    column_headers(ax,5.79)
    y = 6.01
    for group, genes in B_GROUPS:
        ax.text(1.48,y+row_h*len(genes)/2,group,fontsize=8,va='center',fontweight='semibold')
        for gene in genes:
            ax.text(3.50,y+row_h/2,gene,fontsize=8.3,ha='right',va='center',fontstyle='italic')
            for j,column in enumerate(COL_LABELS):
                row=blook.loc[(gene,column)]
                cell(ax,COL_X[j],y,row_h,row.v,row.lab,row.face=='bold',SCREEN if j<2 else MOLECULAR,8 if j<2 else 1.5)
            y+=row_h
        y+=.055
    column_ticks(ax,y-.025)
    colourbars(fig,ax,9.85,'b')

    layout=text_layout(fig)
    a.to_csv(DATA/'panel_a_values.csv',index=False)
    b.to_csv(DATA/'panel_b_values.csv',index=False)
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    preview=OUT/'review'/f'{STEM}_preview.png'
    fig.savefig(preview,dpi=160)
    plt.close(fig)
    assert hashes()==original
    manifest={'figure_number':'Supp Fig 4.7','stem':STEM,'inches':[WIDTH,HEIGHT],
              'png_dpi':600,'font':'Open Sans','source_hashes':original,'source_hashes_unchanged':True,
              'panel_a_cells':len(a),'panel_b_cells':len(b),'columns':A_COLUMNS,
              'excluded_display_column':'Untransduced LOY vs EV','panel_a_groups':A_GROUPS,
              'screen_palette':{'negative_resistor':'#218C8D','zero':'#F7F7F7','positive_sensitizer':'#7A3E85'},
              'molecular_palette':{'negative_lower':'#2166AC','zero':'#F7F7F7','positive_higher':'#B2182B'},
              'panel_a_limits':[-8,8],'panel_b_screen_limits':[-8,8],'panel_b_molecular_limits':[-1.5,1.5],
              'text_layout':layout,'preview':preview.as_posix()}
    (DATA/'render_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'figure':STEM,'text_layout':layout,'panel_a_cells':len(a),'panel_b_cells':len(b)},indent=2))


if __name__=='__main__':main()
