"""Compact workflow, illustrative trajectories and five measured flow panels."""
from pathlib import Path
import argparse
import json
import os
import sys

THESIS = Path(__file__).resolve().parents[3]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--unstained-overlay-draft', action='store_true',
                    help='Export a separate review draft with October boundaries over the historical unstained panel.')
DRAFT_OVERLAY = parser.parse_args().unstained_overlay_draft
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/competition_workflow'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc, FancyArrowPatch, Rectangle
import numpy as np
import pymupdf

sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours
apply_style(8)
COL = condition_colours()
BLUE, RED = COL['Y-retaining'], COL['Y-depleted']
GREEN, CHERRY = '#32965B', '#D34F4F'
INK, GREY, LIGHT = '#222222', '#777E86', '#D9DDE2'
STEM = 'ch04_main_competition_assay_workflow_draft'
if DRAFT_OVERLAY:
    OUT = OUT / 'review/unstained_overlay_draft'
    OUT.mkdir(parents=True, exist_ok=True)
    STEM = 'ch04_main_competition_assay_workflow_unstained_overlay_draft'
FIG_W, FIG_H, Y_MIN = 183, 141.5, -51.5
fig = plt.figure(figsize=(FIG_W / 25.4, FIG_H / 25.4))
ax = fig.add_axes([0, 0, 1, 1], xlim=(0, FIG_W), ylim=(Y_MIN, Y_MIN+FIG_H))
ax.axis('off')
texts = []

def text(x, y, label, size=8, color=INK, weight='normal', ha='center', va='center', linespacing=1.23):
    artist = ax.text(x, y, label, fontsize=size, color=color, fontweight=weight,
                     ha=ha, va=va, linespacing=linespacing)
    texts.append(artist)
    return artist

def arrow(x1, y1, x2, y2, width=.85):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',
                 mutation_scale=8, lw=width, color=GREY, shrinkA=0, shrinkB=0))

def matrix_icon(x, y):
    # Schematic grid: the drawn rows and shading are not matrix dimensions
    # or coefficients. Colours echo GFP, mCherry and autofluorescence.
    for column, color in enumerate((GREEN, CHERRY, GREY)):
        for row in range(4):
            ax.add_patch(Rectangle((x-3.4+column*2.5, y-3.1+row*1.7),
                                   1.8, 1.1, fc=color, ec='none', alpha=.65))
    for side in (-1, 1):
        edge = x+side*4.6
        ax.plot([edge-side*.9, edge, edge, edge-side*.9],
                [y+3.7, y+3.7, y-3.7, y-3.7], color=INK, lw=.75)


def plasmid(x, y, color):
    ax.add_patch(Circle((x,y),4.4,fill=False,ec=GREY,lw=1.05))
    ax.add_patch(Arc((x,y),8.8,8.8,theta1=20,theta2=145,ec=color,lw=3))
    t=np.deg2rad(145)
    arrow(x+4.4*np.cos(t-0.07),y+4.4*np.sin(t-0.07),
          x+4.4*np.cos(t+0.08),y+4.4*np.sin(t+0.08),width=1.1)

def cell(x,y,genotype,reporter):
    border=BLUE if genotype=='CTRL' else RED
    fill=GREEN if reporter=='GFP' else CHERRY
    ax.add_patch(Circle((x,y),4.6,fc='white',ec=border,lw=1.4))
    ax.add_patch(Circle((x,y),3.3,fc=fill,ec='none',alpha=.18))
    ax.add_patch(Circle((x-.8,y+.6),1.25,fc=fill,ec='none',alpha=.75))
    text(x,y-7.4,genotype,8,color=border,weight='semibold')
    text(x,y-11.0,reporter,7.4,color=fill)

dish_icon_path = DATA / 'cartoon_green_red_dish_icon.png'
dish_icon = plt.imread(dish_icon_path)

def dish(x,y):
    # Both reciprocal cultures contain the same two reporters. Their genotype
    # assignment remains explicit in the upstream cell labels, not icon counts.
    width = 22
    height = width * dish_icon.shape[0] / dish_icon.shape[1]
    ax.imshow(dish_icon,extent=(x-width/2,x+width/2,y-height/2,y+height/2),zorder=3)

# Reporter preparation, reciprocal mixtures and cytometric readout.
text(2,87,'a',10,weight='semibold',ha='left')
text(24,86,'Reporter plasmids',8.5,weight='semibold')
text(76,86,'Dye-swap transduction',8.5,weight='semibold')
text(113,86,'Mixing',8.5,weight='semibold')
text(156,86,'Flow cytometry',8.5,weight='semibold')

plasmid(8,73,GREEN)
text(15,75,'pGIPZ',8.3,ha='left',weight='semibold')
text(15,70.5,'TurboGFP',7.8,color=GREEN,ha='left')
plasmid(8,53,CHERRY)
text(15,56,'pCDH-CMV-\nmCherry-T2A-Puro',7.8,ha='left',weight='semibold')
text(15,48.5,'mCherry',7.8,color=CHERRY,ha='left')
# A shared connector indicates both reporters are used in the labelling stage.
ax.plot([44.5,47,47,44.5],[76,76,50,50],color=LIGHT,lw=.8)
arrow(47,63,54,63)

for y,swap,number in [(74,False,1),(54,True,2)]:
    cell(65,y,'CTRL','mCherry' if swap else 'GFP')
    cell(87,y,'LOY','GFP' if swap else 'mCherry')
    text(76,y,'+',10,color=GREY)
    arrow(94,y,103,y)
    dish(113,y)
    text(113,y-8.5,f'Mixture {number}',7.2,color=GREY)

ax.plot([123,126,126,123],[74,74,54,54],color=LIGHT,lw=.8)
if DRAFT_OVERLAY:
    # Keep the flow block on one centreline, with a compact icon/label row.
    arrow(126,68.5,131.8,68.5)
    text(156,79.5,'BD FACSDiscover A8',8.1,weight='semibold')
    arrow(156,76.5,156,73.2,width=.65)
    matrix_icon(138,68.5)
    text(146.2,68.5,'Spectral unmixing',7.2,weight='semibold',ha='left')
    text(156,60.5,'Pre-existing GFP/mCherry matrix\nwith autofluorescence reference',
         6.9,color=GREY,linespacing=1.06)
    arrow(156,55,156,52.5,width=.65)
    text(156,47.5,'Background exclusion\nand GFP:mCherry quantification',
         7,linespacing=1.08)
else:
    arrow(126,68.5,133,68.5)
    text(157,79.5,'BD FACSDiscover A8',8.1,weight='semibold')
    arrow(157,76.7,157,75,width=.65)
    matrix_icon(140,68.5)
    text(164,70.5,'Pre-existing\nunmixing matrix',7.2,weight='semibold')
    text(164,64.5,'mcherry_GFP_matrix',6.5,color=GREY)
    text(157,58,'GFP · mCherry\n+ autofluorescence reference',6.9,color=GREY)
    arrow(157,53,157,50.5,width=.65)
    text(157,46.5,'Background exclusion\n+ reporter assignment',7)

# One time-proportional schedule is shared by the reciprocal cultures.
ax.plot([3,180],[40.5,40.5],color=LIGHT,lw=.65)
text(4,36.5,'28-day follow-up',8.1,weight='semibold',ha='left')
day_x=lambda day: 11 + day/28*83
arrow(day_x(0),22,day_x(28)+4,22)
READ_DAYS=[0,4,8,14,21,28]
for day in READ_DAYS:
    x=day_x(day)
    ax.add_patch(Circle((x,22),1.1,fc='white',ec=INK,lw=.9))
    text(x,27.5,f'D{day}',7.7,weight='semibold')
text(day_x(0),16.5,'Baseline',7.1,color=GREY)

# Deliberately schematic outcomes, independent of the acquired D0 events.
# The ordinate communicates direction only; these are not fitted effects.
if not DRAFT_OVERLAY:
    text(143,37.5,'Illustrative hit outcomes',7.6,weight='semibold')
trajectory_ax = fig.add_axes([115/FIG_W,(16-Y_MIN)/FIG_H,47/FIG_W,19/FIG_H])
example_days = np.asarray(READ_DAYS,dtype=float)
example_change = example_days/28
enrichment_color, depletion_color = '#7954A1', '#A47720'
trajectory_ax.axhline(0,color=GREY,lw=.6,ls=(0,(2,2)),zorder=0)
trajectory_ax.plot(example_days,example_change,color=enrichment_color,lw=1.1,
                   marker='o',ms=2.5,mfc='white',mew=.8)
trajectory_ax.plot(example_days,-example_change,color=depletion_color,lw=1.1,
                   ls=(0,(4,2)),marker='s',ms=2.4,mfc='white',mew=.8)
trajectory_ax.set(xlim=(-.6,28.6),ylim=(-1.17,1.17),xticks=READ_DAYS,
                  xticklabels=[f'D{day}' for day in READ_DAYS],
                  yticks=[-1,0,1],yticklabels=['−','0','+'],
                  ylabel='log(LOY:CTRL)\nchange from D0')
trajectory_ax.tick_params(axis='both',labelsize=6.5,length=2,pad=1.5)
trajectory_ax.xaxis.label.set_size(6.8)
trajectory_ax.yaxis.label.set_size(6.5)
trajectory_ax.xaxis.labelpad=1
trajectory_ax.yaxis.labelpad=2
for spine in trajectory_ax.spines.values():
    spine.set_linewidth(.6)
for y,label,color in [(1,'Enrichment',enrichment_color),(0,'No change',GREY),(-1,'Depletion',depletion_color)]:
    texts.append(trajectory_ax.text(30,y,label,color=color,fontsize=6.5,va='center',clip_on=False))
illustrative_record = {
    'status':'ILLUSTRATIVE ONLY — not experimental measurements or fitted estimates',
    'days':READ_DAYS,
    'schematic_enrichment':example_change.tolist(),
    'schematic_depletion':(-example_change).tolist(),
    'ordinate':'Change from D0 in natural log(LOY:CTRL); only direction is labelled on the figure.',
    'interpretation':'Alternative hypothetical enrichment or depletion outcomes of a hit. These curves do not denote the two reciprocal reporter orientations.',
    'magnitude':'Arbitrary drawing coordinates; not a forecast or effect-size estimate.',
}
if not DRAFT_OVERLAY:
    (DATA/'illustrative_hit_trajectories.json').write_text(json.dumps(illustrative_record,indent=2)+'\n',encoding='utf-8')

# Actual events: October single-reporter controls and mixtures, plus a
# September unstained reference with its own acquisition's singlet gates.
manifest = json.loads((DATA / 'flow_examples_manifest.json').read_text(encoding='utf-8'))
control_manifest = json.loads((DATA / 'flow_controls_manifest.json').read_text(encoding='utf-8'))
weights = np.asarray(manifest['classifier']['weights'])
intercept = manifest['classifier']['intercept']
gfp_limit, mcherry_limit = manifest['negative_rule']['limits']
ax.plot([3,180],[9.5,9.5],color=LIGHT,lw=.65)
# Extend the shared lower limit slightly to retain every unstained event.
plot_limits = (-8.5,12.5)
controls_by_label = {well['label']: well for well in control_manifest['controls']}
panel_records = [controls_by_label[label] for label in ('Unstained', 'GFP only', 'mCherry only')] + manifest['wells']
for index, well in enumerate(panel_records):
    left = 12 + index*34.3
    centre = left+15.9
    historical = well.get('historical_unstained', False)
    mixture = well.get('mixture')
    title = f'Mixture {mixture}' if mixture else well['label']
    subtitle = ('WT · Oct D0' if mixture else
                'CTRL · 1 Sep 2026' if historical else 'CTRL · Oct D0')
    if DRAFT_OVERLAY and index < 4:
        subtitle = 'WT' if mixture else 'CTRL'
    text(left,6,chr(ord('b')+index),9,weight='semibold',ha='left')
    text(centre+1.5,6,title,7.5,weight='semibold')
    text(centre,1.5,subtitle,6.4,color=GREY)
    plot_ax = fig.add_axes([left/FIG_W, (-33.3-Y_MIN)/FIG_H, 31.8/FIG_W, 31.8/FIG_H])
    events = np.genfromtxt(DATA / well['events_csv'], delimiter=',', names=True, dtype=None, encoding='utf-8')
    values = np.column_stack([events['eGFP_arcsinh'], events['mCherry_arcsinh']])
    assert len(values) == well['primary_singlets']
    assert ((values >= plot_limits[0]) & (values <= plot_limits[1])).all()
    if not historical or DRAFT_OVERLAY:
        plot_ax.add_patch(Rectangle((plot_limits[0], plot_limits[0]),
            gfp_limit-plot_limits[0], mcherry_limit-plot_limits[0],
            facecolor='#F2F3F5', edgecolor='#A3A9B0', linewidth=.5,
            linestyle=(0,(3,2)), zorder=0))
        line_x = np.linspace(*plot_limits,500)
        line_y = -(weights[0]*line_x+intercept)/weights[1]
        plot_ax.plot(line_x,line_y,color=INK,lw=.55,ls=(0,(4,3)),zorder=3)
    if historical:
        plot_ax.scatter(values[:,0],values[:,1],s=.75,c=GREY,alpha=.23,
                        edgecolors='none',rasterized=True,zorder=2)
    else:
        excluded = events['excluded'].astype(bool)
        gfp = events['reporter_class'] == 'GFP'
        for mask, color in ((excluded, '#A3A9B0'), (~excluded & ~gfp, CHERRY), (~excluded & gfp, GREEN)):
            plot_ax.scatter(values[mask,0],values[mask,1],s=.85,c=color,alpha=.32,
                            edgecolors='none',rasterized=True,zorder=2)
    plot_ax.set(xlim=plot_limits, ylim=plot_limits, xticks=[-5,0,5,10], yticks=[-5,0,5,10],
                ylabel='mCherry' if index == 0 else '')
    plot_ax.tick_params(axis='both',labelsize=6.4,length=1.6,pad=1.3)
    plot_ax.tick_params(axis='y',labelleft=index == 0)
    plot_ax.yaxis.label.set_size(7.2)
    plot_ax.yaxis.labelpad=2
    for spine in plot_ax.spines.values():
        spine.set_linewidth(.6)
    if historical:
        text(centre,-38.5,'Historical reference',6.5,color=GREY)
        if DRAFT_OVERLAY:
            text(centre,-41.9,'Draft gate overlay',6.3,color='#7954A1')
    elif mixture:
        gfp_genotype, mcherry_genotype = ('CTRL','LOY') if mixture == 1 else ('LOY','CTRL')
        text(centre,-38.5,f"{gfp_genotype}–GFP {100*well['gfp_fraction']:.1f}%",color=GREEN,size=6.5)
        text(centre,-41.9,f"{mcherry_genotype}–mCherry {100*(1-well['gfp_fraction']):.1f}%",color=CHERRY,size=6.5)
    else:
        expected_gfp = well['label'] == 'GFP only'
        percent = 100*(well['gfp_fraction'] if expected_gfp else 1-well['gfp_fraction'])
        text(centre,-38.5,f"{'GFP' if expected_gfp else 'mCherry'} {percent:.1f}%",6.5,color=GREEN if expected_gfp else CHERRY)
        text(centre,-41.9,f"Excluded {100*well['negative_broad_fraction']:.1f}%",6.3,color=GREY)
text(96,-48,'GFP   ·   Both axes: arcsinh(unmixed signal / 1000)',7)

fig.canvas.draw()
renderer=fig.canvas.get_renderer()
bounds=[]
for artist in texts:
    box=artist.get_window_extent(renderer)
    if not fig.bbox.contains(box.x0,box.y0) or not fig.bbox.contains(box.x1,box.y1):
        bounds.append(artist.get_text())
assert not bounds, f'Text outside figure: {bounds}'
for ext in ('pdf','svg','png'):
    fig.savefig(OUT/f'{STEM}.{ext}',dpi=600)
review=OUT if DRAFT_OVERLAY else OUT/'review'
review.mkdir(exist_ok=True)
if DRAFT_OVERLAY:
    fig.savefig(review/f'{STEM}_preview.png',dpi=190)
else:
    fig.savefig(review/f'{STEM}.png',dpi=190)
    fig.savefig(review/f'{STEM}_unstained_first.png',dpi=190)
plt.close(fig)
doc=pymupdf.open(OUT/f'{STEM}.pdf')
page=doc[0]
page.get_pixmap(matrix=pymupdf.Matrix(1.8,1.8),colorspace=pymupdf.csGRAY).save(review/f'{STEM}_greyscale.png')
svg=(OUT/f'{STEM}.svg').read_text(encoding='utf-8')
validation={
    'figure_mm':[FIG_W,FIG_H], 'png_dpi':600, 'svg_editable_text': '<text' in svg,
    'instrument_icon_displayed':False, 'instrument_label':'BD FACSDiscover A8',
    'unmixing_matrix':'mcherry_GFP_matrix',
    'unmixing_matrix_icon':'Editable vector bracketed grid; schematic dimensions and shading, no coefficients shown',
    'unmixing_components':['eGFP','mCherry','Autofluorescence_001'],
    'unmixing_distinct_from_classifier':True,
    'dish_icon':str(dish_icon_path), 'dish_icon_format':'embedded transparent PNG',
    'dish_icon_semantics':'Simple cartoon green/red cells in retained dish; counts and proportions are not data',
    'reading_days':READ_DAYS, 'drug_content_displayed':False,
    'illustrative_trajectories':illustrative_record,
    'illustrative_trajectory_heading_displayed':not DRAFT_OVERLAY,
    'pdf_pages':len(doc), 'text_outside_canvas':bounds,
    'fonts':[f[3] for f in page.get_fonts()],
    'schematic_only':False,
    'flow_samples':[w['sample'] for w in panel_records],
    'flow_points_displayed':[w['primary_singlets'] for w in panel_records],
    'historical_unstained_sample':'A01-unstained_CTRL; 1 September 2026',
    'historical_unstained_october_gates_or_classifier_applied':False,
    'historical_unstained_october_boundaries_displayed':DRAFT_OVERLAY,
    'date_labels_omitted_from_panels':['b','c','d','e'] if DRAFT_OVERLAY else [],
    'flow_classifier_refitted':False, 'flow_axes_limits':list(plot_limits),
    'flow_points_clipped_or_subsampled':False,
    'notes':'Panel a and its enrichment/depletion trajectories are illustrative. Panels c-f show October D0 singlets. Panel b is a historical September unstained CTRL reference with its own CEF gates, not a same-day October control. Percentages use saved reporter-positive denominators. No gates or classifier parameters were refitted.'
}
if DRAFT_OVERLAY:
    validation['notes'] += ' Review-only overlay: the saved October classifier boundary and exclusion rectangle are displayed on panel b without classifying or deleting its grey events. Main figure exports are unchanged.'
(OUT/'export_validation.json' if DRAFT_OVERLAY else DATA/'export_validation.json').write_text(json.dumps(validation,indent=2)+'\n',encoding='utf-8')
print(json.dumps(validation))
