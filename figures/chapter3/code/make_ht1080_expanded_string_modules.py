"""Display a thematic selection from the frozen broader STRING analysis."""
from pathlib import Path
import sys,json,shutil,hashlib,textwrap
import pandas as pd
import networkx as nx
FIG=Path(__file__).resolve().parents[1]
THESIS=FIG.parents[1]
ANALYSIS=THESIS.parent/'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE'
sys.path.insert(0,str(ANALYSIS))
from plot_network_modules import plot_module,symbol_size
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(7)
SOURCE=FIG/'source_data/ht1080_expanded_string_modules'
SOURCE.mkdir(parents=True,exist_ok=True)
# Use the frozen thesis tables; this renderer never refits or refreshes analyses.
input_files=[SOURCE/name for name in ['modules.csv','candidate_node_evidence.csv','annotation_significance_plot_data.csv','module_top_annotations.csv','module_labels_reviewed.json','STRING_association_DOWN.graphml','STRING_association_UP.graphml']]
input_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in input_files}
modules=pd.read_csv(SOURCE/'modules.csv').set_index('module_id',drop=False)
genes=pd.read_csv(SOURCE/'candidate_node_evidence.csv').set_index('gene_name')
labels=json.loads((SOURCE/'module_labels_reviewed.json').read_text(encoding='utf-8'))
ids={'DOWN':[1,2,3,4,8,10,11,15,18],'UP':[7,3,16,5,6,8,9,11,18]}
fig=plt.figure(figsize=(183/25.4,250/25.4))
positions=[]; chosen=[]; texts=[]
for block,direction in enumerate(['DOWN','UP']):
 color={'DOWN':'#0072B2','UP':'#D55E00'}[direction]
 fig.text(.02,.983-block*.475,'ab'[block],fontsize=10,fontweight='bold',va='top')
 fig.text(.06,.983-block*.475,('Lower' if direction=='DOWN' else 'Higher')+' expression in Y-depleted HT1080',fontsize=9,fontweight='semibold',color=color,va='top')
 graph=nx.read_graphml(SOURCE/f'STRING_association_{direction}.graphml')
 for i,num in enumerate(ids[direction]):
  mid=f'{direction}_{num:02d}';row=modules.loc[mid];chosen.append(mid)
  x=.045+(i%3)*.325;y=.829-block*.475-(i//3)*.147
  ax=fig.add_axes([x,y,.285,.104])
  label='\n'.join(textwrap.wrap(labels[mid],width=34))
  if mid=='UP_07':label='Eicosanoid metabolism,\nglutathione enzymes and cytokines'
  if mid=='UP_03':label='Chemokine receptors and\nimmune co-signalling'
  if mid=='UP_16':label='TREM receptors and\nCKLF/CMTM family'
  drawn=plot_module(ax,row,graph,genes,label,color,positions,'main')
  texts.extend(drawn)
  ax.text(0,-.055,f'n = {row.nodes}; annotation q = {row.BH_all_modules_and_terms:.1e}',transform=ax.transAxes,fontsize=5.5,color='#59636A')
handles=[Line2D([0],[0],color='#596A72',lw=.7),Line2D([0],[0],color='#B8BEC2',lw=.7,ls='--')]
fig.legend(handles,['Physical support '+chr(8805)+' 0.7','Other STRING association'],loc='lower left',bbox_to_anchor=(.035,.006),ncol=2,frameon=False,fontsize=6)
fig.text(.97,.020,'Node area scales with |mean log'+chr(8322)+' RNA difference|',ha='right',fontsize=6)
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
overlap=[]
for ax in fig.axes:
 labs=[t for t in ax.texts if t.get_fontstyle()=='italic']
 for i,t in enumerate(labs):
  for u in labs[i+1:]:
   if t.get_window_extent(renderer).overlaps(u.get_window_extent(renderer)):overlap.append([t.get_text(),u.get_text()])
assert not overlap,overlap
stem='ch03_main_ht1080_expanded_string_modules'
from io import BytesIO
for ext in ['pdf','svg','png']:
 buffer=BytesIO();fig.savefig(buffer,format=ext,dpi=400);(FIG/f'{stem}.{ext}').write_bytes(buffer.getvalue())
modules.loc[chosen].to_csv(SOURCE/'displayed_modules.csv',index=False)
pd.DataFrame(positions).to_csv(SOURCE/'node_and_label_positions.csv',index=False)
assert all(hashlib.sha256(p.read_bytes()).hexdigest()==input_hashes[p.name] for p in input_files)
(SOURCE/'provenance.json').write_text(json.dumps({'analysis':str(ANALYSIS),'selection':ids,'rationale':'Thematic selection unchanged: 9 lower- and 9 higher-expression modules. U07 and U16 retain the author-selected positions. Not a significance-ranked subset.','candidate_rule':'absolute mean log2 expression difference >=0.5, period_direction_agrees, LOO_all_same_direction; no gene-level P-value threshold','label_review':'2026-10-09; descriptive labels broadened for mixed modules; NLGN4X dagger and explanatory note removed at author request; British signalling. No membership, topology, effects or annotation statistics changed. D15 is a module label, not evidence of an assembled complex in HT1080.','label_overlap_pairs':overlap,'input_sha256':input_hashes,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCE.glob('*.csv')}},indent=2),encoding='utf-8')
fig.savefig(FIG/'previews'/f'{stem}_preview.png',dpi=180)
print('Saved 18 modules:',int(modules.loc[chosen].nodes.sum()),'displayed nodes')

