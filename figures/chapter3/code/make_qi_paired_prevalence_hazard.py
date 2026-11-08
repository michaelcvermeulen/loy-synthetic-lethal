"""Paired rank display of existing Qi prevalence and published PFS estimates."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'figures/chapter3'
DATA = OUT / 'source_data/qi_paired_prevalence_hazard'
DATA.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
font = apply_style(9)
src_a = OUT / 'source_data/qi_tcga_loy_prevalence'
src_b = OUT / 'source_data/qi_tcga_pfi'
ref = pd.read_csv(src_a / 'prevalence_by_tumour_type.csv')
long = pd.read_csv(src_a / 'y_status_distribution_long.csv')
hazard = pd.read_csv(src_b / 'forest_display_table.csv')
palette = json.loads((src_a / 'y_status_distribution_provenance.json').read_text())
overall = json.loads((src_b / 'published_overall_estimate.json').read_text())
statuses = palette['status_order']
colours = palette['colours']
order_a = ref.Tumor_type.tolist()
ranked = hazard[hazard.LOY_HR.notna()].sort_values('LOY_HR', ascending=False)
missing = hazard[hazard.LOY_HR.isna()]
order_b = ranked.Tumor_type.tolist() + missing.Tumor_type.tolist()
assert len(order_a) == len(order_b) == len(set(order_b)) == 28
assert set(order_a) == set(order_b)
assert ranked.LOY_HR.is_monotonic_decreasing
assert set(missing.Tumor_type) == {'LAML', 'PCPG'}
a_rank = {t:i for i,t in enumerate(order_a)}
b_row = {t:i for i,t in enumerate(order_b)}
# Unreported estimates occupy a separate bottom block and have no hazard rank.
b_y = {t:float(i) if i<26 else i+.6 for t,i in b_row.items()}
joined = ref[['Tumor_type','tumour_name','n_called','n_complete_loy','percentage_complete_loy']].merge(
    hazard[['Tumor_type','LOY_HR','ci_lower','ci_upper','LOY_Cox_p','LOY_q']],
    on='Tumor_type', validate='one_to_one')
joined['prevalence_rank'] = joined.Tumor_type.map(a_rank)+1
joined['hazard_rank'] = joined.Tumor_type.map({t:i+1 for i,t in enumerate(ranked.Tumor_type)})
joined['right_display_y'] = joined.Tumor_type.map(b_y)
joined.to_csv(DATA/'paired_rank_table.csv', index=False)
lookup = joined.set_index('Tumor_type')
percent = long.pivot(index='cohort', columns='status', values='percentage').reindex(columns=statuses)
assert np.allclose(percent.loc[order_a].sum(axis=1), 100)
assert np.allclose(percent.loc[order_a,'LOY'], ref.percentage_complete_loy)

fig = plt.figure(figsize=(280/25.4, 215/25.4))
bottom, height = .14, .69
ax_a = fig.add_axes([.060,bottom,.235,height])
ax_b = fig.add_axes([.565,bottom,.205,height])
ax_link = fig.add_axes([.345,bottom,.160,height])
for ax in [ax_a, ax_b, ax_link]:
    ax.set_ylim(30.3,-.7)
    ax.set_yticks([])
ax_link.set_xlim(0,1)
ax_link.axis('off')
ax_a.set_xlim(100,0)
ax_a.set_xticks([100,80,60,40,20,0])
ax_a.set_xlabel('Tumours by Y-chromosome status (%)', fontsize=9, labelpad=8)
ax_a.spines['left'].set_visible(False)
ax_b.spines['left'].set_visible(False)
ax_b.set_xscale('log')
ax_b.set_xlim(.025,64)
ticks=[.05,.25,1,4,16,64]
ax_b.xaxis.set_major_locator(FixedLocator(ticks))
ax_b.xaxis.set_major_formatter(FixedFormatter([str(t) for t in ticks]))
ax_b.xaxis.set_minor_locator(NullLocator())
ax_b.axvline(1, color='#999999', lw=.7, ls=(0,(3,3)))
ax_b.set_xlabel('Hazard ratio (LOY versus WT)', fontsize=9, labelpad=8)
for ax in [ax_a,ax_b]:ax.tick_params(axis='x',labelsize=8,length=3)

counter={}
def tag(artist, code, role):
    key=f'cancer-{code}-{role}'
    counter[key]=counter.get(key,0)+1
    artist.set_gid(f'{key}-{counter[key]}')
    return artist
def label(ax, x, y, txt, code=None, **kwargs):
    artist=ax.text(x,y,txt,transform=ax.get_yaxis_transform(),va='center',clip_on=False,**kwargs)
    if code:tag(artist,code,'text')
    return artist

# A is mirrored horizontally: the complete-LOY segment starts at the inner edge.
for code in order_a+['All tumours']:
    y=a_rank[code] if code in a_rank else 29.3
    left=0
    for status in statuses:
        value=percent.loc[code,status]
        rect=ax_a.barh(y,value,left=left,height=.66,color=colours[status],edgecolor='white',
                       linewidth=.35,hatch=palette['hatches'].get(status)).patches[0]
        if code in a_rank:tag(rect,code,'bar')
        if value>=8:
            t=ax_a.text(left+value/2,y,f'{value:.0f}',ha='center',va='center',fontsize=7.3,
                          color='white' if status in ['LOY','Gain','WT'] else '#222222')
            if code in a_rank:tag(t,code,'percent')
        left+=value
    significant=code in a_rank and pd.notna(lookup.loc[code,'LOY_q']) and lookup.loc[code,'LOY_q']<.05
    label(ax_a,1.045,y,code if code in a_rank else 'Overall',code if code in a_rank else None,
          fontsize=8.5,weight='semibold' if significant or code=='All tumours' else 'normal')
    n=int(lookup.loc[code,'n_called']) if code in a_rank else 5014
    label(ax_a,-.065,y,f'{n:,}',ha='right',fontsize=8,color='#555555')
label(ax_a,-.065,-1.6,'n',ha='right',fontsize=8.5,weight='semibold')

for code in order_b:
    row=lookup.loc[code];y=b_y[code];sig=pd.notna(row.LOY_q) and row.LOY_q<.05
    label(ax_b,-.060,y,code,code,ha='right',fontsize=8.5,weight='semibold' if sig else 'normal')
    if pd.isna(row.LOY_HR):
        label(ax_b,.5,y,'Not reported',code,ha='center',fontsize=8,color='#777777')
        continue
    artists=ax_b.errorbar(row.LOY_HR,y,xerr=[[row.LOY_HR-row.ci_lower],[row.ci_upper-row.LOY_HR]],
                  fmt='o',ms=4,color=colours['LOY'],elinewidth=.85,capsize=2,markeredgewidth=.7)
    tag(artists.lines[0],code,'point')
    for artist in list(artists.lines[1])+list(artists.lines[2]):tag(artist,code,'interval')
    label(ax_b,1.415,y,f'{row.LOY_HR:.2f} ({row.ci_lower:.2f}\u2013{row.ci_upper:.2f})',code,
          fontsize=8,ha='center',weight='semibold' if sig else 'normal')
    label(ax_b,2.025,y,f'{row.LOY_q:.3f}',code,fontsize=8,ha='right',
          color=colours['LOY'] if sig else '#444444',weight='semibold' if sig else 'normal')
label(ax_b,1.415,-1.6,'HR (95% CI)',ha='center',fontsize=8.5,weight='semibold')
label(ax_b,2.025,-1.6,'q',ha='right',fontsize=8.5,weight='semibold')
# Lines show matching identities and changes in order, not effect magnitudes.
for code in sorted(order_a, key=lambda t: pd.notna(lookup.loc[t,'LOY_q']) and lookup.loc[t,'LOY_q']<.05):
    sig=pd.notna(lookup.loc[code,'LOY_q']) and lookup.loc[code,'LOY_q']<.05
    unavailable=pd.isna(lookup.loc[code,'LOY_HR'])
    ya,yb=a_rank[code],b_y[code]
    path=MplPath([(0,ya),(.35,ya),(.65,yb),(1,yb)],
                 [MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4])
    line=PathPatch(path,facecolor='none',edgecolor=colours['LOY'] if sig else '#A4A8AC',
                   lw=1.15 if sig else .65,alpha=.92 if sig else .60,
                   ls=(0,(3,3)) if unavailable else 'solid')
    tag(line,code,'connector');ax_link.add_patch(line)

ax_b.axhline(25.6,color='#DDDDDD',lw=.65)
for ax in [ax_a,ax_b]:ax.axhline(28.45,color='#AAAAAA',lw=.65)
ax_b.errorbar(overall['HR'],29.3,xerr=[[overall['HR']-overall['ci_lower']],
    [overall['ci_upper']-overall['HR']]],fmt='D',ms=4,color='#444444',capsize=2,elinewidth=.85)
label(ax_b,-.06,29.3,'Overall',ha='right',fontsize=8.5,weight='semibold')
label(ax_b,1.415,29.3,'1.17 (1.06\u20131.28)',ha='center',fontsize=8,weight='semibold')
label(ax_b,.03,29.3,'Unadjusted',ha='left',fontsize=7,color='#555555')

fig.text(.025,.97,'a',fontsize=13,weight='semibold',va='top')
fig.text(.060,.97,'Y-chromosome status',fontsize=11,weight='semibold',va='top')
fig.text(.060,.936,'Ranked by complete-LOY prevalence',fontsize=9,va='top')
fig.text(.525,.97,'b',fontsize=13,weight='semibold',va='top')
fig.text(.565,.97,'Complete LOY and PFS',fontsize=11,weight='semibold',va='top')
fig.text(.565,.936,'Ranked by adjusted hazard ratio',fontsize=9,va='top')
fig.text(.565,.904,'Adjusted for age and chromosome-arm losses',fontsize=8.5,va='top')
handles=[Patch(facecolor=colours[s],edgecolor='white',hatch=palette['hatches'].get(s),label=palette['labels'][s]) for s in statuses]
fig.legend(handles=handles,ncol=3,loc='upper left',bbox_to_anchor=(.052,.912),fontsize=8,
           columnspacing=1.25,handlelength=1.5,labelspacing=.5)
fig.text(.425,.862,'Matched cancer types',ha='center',fontsize=8,color='#666666')
fig.text(.565,.071,'Lower hazard',fontsize=8,color='#555555')
fig.text(.770,.071,'Higher hazard',ha='right',fontsize=8,color='#555555')
fig.text(.025,.027,'Source: Qi et al. (2023), Tables S1A and S4; published analysis notebook.',fontsize=7.5,color='#555555')
stem='ch03_main_intro_tcga_prevalence_hazard_paired'
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{stem}.{ext}',dpi=600)
fig.savefig(OUT/'previews'/f'{stem}.png',dpi=160)
plt.close(fig)
inputs=[src_a/'prevalence_by_tumour_type.csv',src_a/'y_status_distribution_long.csv',
        src_a/'y_status_distribution_provenance.json',src_b/'forest_display_table.csv',
        src_b/'published_overall_estimate.json']
meta=dict(figure=stem,size_mm=[280,215],font=font,left_order=order_a,right_order=order_b,
          right_unranked=missing.Tumor_type.tolist(),connectors=28,
          highlight='Published LOY q < 0.05 (UVM, PRAD, MESO)',
          ranks='Prevalence: decreasing complete LOY. Hazard: decreasing adjusted HR; unreported rows unranked.',
          source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
          script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          analysis='Layout only; no new survival models, confidence intervals, or association tests.',
          overall='Shown separately on both sides without a connector; different source cohorts/denominators.')
(DATA/'render_provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print('Rendered paired figure: 28 identity links, 26 ranked adjusted HRs, 2 unranked missing estimates.')

# Standalone interactive companion: the scientific marks are the exported SVG.
# Hover or focus any cancer label/mark to trace the same cancer on both sides.
import base64
import xml.etree.ElementTree as ET
svg_ns='http://www.w3.org/2000/svg'
ET.register_namespace('',svg_ns)
ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
tree=ET.parse(OUT/f'{stem}.svg')
svg=tree.getroot()
details={}
for code in order_a:
    r=lookup.loc[code]
    h=f'Adjusted HR {r.LOY_HR:.2f} (95% CI {r.ci_lower:.2f}\u2013{r.ci_upper:.2f}); q = {r.LOY_q:.3f}' if pd.notna(r.LOY_HR) else 'Adjusted HR not reported'
    details[code]=f'{r.tumour_name} ({code}) | Complete LOY: {r.percentage_complete_loy:.1f}% ({int(r.n_complete_loy)}/{int(r.n_called)}) | {h}'
for el in svg.iter():
    gid=el.attrib.get('id','')
    if gid.startswith('cancer-'):
        code=gid.split('-')[1]
        el.set('class','cancer')
        el.set('data-cancer',code)
        title=ET.SubElement(el,f'{{{svg_ns}}}title');title.text=details[code]
        if '-connector-' in gid:
            el.set('class','cancer connector')
        if '-text-1' in gid:
            el.set('tabindex','0')
            el.set('role','button')
            el.set('aria-label',details[code])
style=ET.SubElement(svg,f'{{{svg_ns}}}style')
style.text='.cancer{cursor:pointer;transition:opacity .12s} .dim{opacity:.10} .active.connector path{stroke-width:2.5!important;stroke-opacity:1!important} .active text{font-weight:600!important}'
svg_text=ET.tostring(svg,encoding='unicode')
font_root=ROOT/'skills/thesis-figure-style/assets/fonts'
font_css=''
for pattern,weight in [('OpenSans-Regular.ttf',400),('OpenSans-SemiBold.ttf',600)]:
    candidates=list(font_root.rglob(pattern))
    if candidates:
        encoded=base64.b64encode(candidates[0].read_bytes()).decode()
        font_css+=f"@font-face{{font-family:'Open Sans';font-weight:{weight};src:url(data:font/ttf;base64,{encoded}) format('truetype');}}"
html='''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LOY prevalence and PFS: paired cancer types</title>
<style>FONT_CSS
body{margin:0;background:#f4f5f6;color:#222;font-family:'Open Sans',Arial,sans-serif}
header{padding:14px 24px;background:white;border-bottom:1px solid #ddd;position:sticky;top:0;z-index:2}
label{font-size:14px} select,button{font:inherit;margin-left:10px;padding:5px 9px;border:1px solid #aaa;border-radius:4px;background:white}
#detail{font-size:13px;min-height:36px;margin-top:8px;max-width:1400px}
main{background:white;max-width:1600px;margin:12px auto}svg{width:100%;height:auto;display:block}
</style><header><label for="cancer">Trace a cancer type</label><select id="cancer"><option value="">All cancer types</option></select><button id="reset">Reset</button>
<div id="detail" aria-live="polite">Hover over a label, bar or connecting line, or select a cancer type. Click to keep it highlighted.</div></header>
<main>SVG_CONTENT</main><script>
const details=DETAILS;
const menu=document.querySelector('#cancer'),info=document.querySelector('#detail');
const marks=[...document.querySelectorAll('[data-cancer]')];let pinned='';
Object.keys(details).sort().forEach(c=>{const o=document.createElement('option');o.value=c;o.textContent=details[c].split(' | ')[0];menu.append(o)});
function show(c){marks.forEach(e=>{e.classList.toggle('dim',!!c&&e.dataset.cancer!==c);e.classList.toggle('active',!!c&&e.dataset.cancer===c)});
info.textContent=c?details[c]:'Hover over a label, bar or connecting line, or select a cancer type. Click to keep it highlighted.';}
marks.forEach(e=>{const c=e.dataset.cancer;e.addEventListener('mouseenter',()=>show(c));e.addEventListener('mouseleave',()=>show(pinned));e.addEventListener('focus',()=>show(c));e.addEventListener('blur',()=>show(pinned));e.addEventListener('click',()=>{pinned=c;menu.value=c;show(c)});e.addEventListener('keydown',ev=>{if(ev.key==='Enter'||ev.key===' '){ev.preventDefault();pinned=c;menu.value=c;show(c)}})});
menu.addEventListener('change',()=>{pinned=menu.value;show(pinned)});
function reset(){pinned='';menu.value='';show('')}
document.querySelector('#reset').addEventListener('click',reset);document.addEventListener('keydown',e=>{if(e.key==='Escape')reset()});
</script></html>'''.replace('FONT_CSS',font_css).replace('SVG_CONTENT',svg_text).replace('DETAILS',json.dumps(details))
(OUT/f'{stem}.html').write_text(html,encoding='utf-8')
print('Saved standalone hoverable HTML companion.')
