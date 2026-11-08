"""All-collection GSEA report, representative main figure and complete hit atlas."""
from pathlib import Path
import sys,json,hashlib,html,textwrap,base64,os
import networkx as nx
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parent;BASE=ROOT/'outputs/programme_analysis_v1';OUT=BASE/'all_collections'
THESIS=ROOT.parents[3]/'SEPT_11_2026_THESIS';FIG=THESIS/'figures/chapter3'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
cmap=LinearSegmentedColormap.from_list('signed',['#3C5488','#F7F7F7','#E64B35'])
g=pd.read_csv(OUT/'programme_consistency.csv');allg=pd.read_csv(OUT/'GSEA_all_contrasts.csv')
tf=pd.read_csv(BASE/'tables/TF_consistency.csv')
for contrast,q in allg.groupby('contrast'):
 assert q.pathway.is_unique
 valid=q.pval.notna();p=q.loc[valid,'pval'].to_numpy()
 # Adjust with the full tested-family size, including any unavailable results.
 order=np.argsort(p);adjusted=np.minimum.accumulate((p[order]*len(q)/np.arange(1,len(p)+1))[::-1])[::-1]
 expected=np.empty_like(p);expected[order]=np.minimum(adjusted,1)
 assert np.allclose(expected,q.loc[valid,'FDR'],rtol=1e-9,atol=1e-12)
for r in pd.read_csv(OUT/'input_hashes.csv').itertuples():
 assert hashlib.sha256((ROOT/r.path).read_bytes()).hexdigest()==r.sha256
labels={'H':'Hallmark','C1':'Positional','C2':'Curated','C3':'Regulatory targets','C4':'Computational','C5':'Ontology','C6':'Oncogenic','C7':'Immunologic','C8':'Cell type','C9':'Cell perturbation'}
order=['H']+[f'C{i}' for i in range(1,10)]
summary=[]
for (contrast,col),q in allg.groupby(['contrast','collection']):
 s=q.FDR<.05
 summary.append(dict(contrast=contrast,collection=col,description=labels.get(col,col),tested=len(q),
  finite_P=int(q.pval.notna().sum()),FDR05=int(s.sum()),higher=int((s&(q.NES>0)).sum()),lower=int((s&(q.NES<0)).sum())))
pd.DataFrame(summary).to_csv(OUT/'collection_summary.csv',index=False)
g['direction']=np.where(g.NES>0,'Higher','Lower')
hits=g[g.FDR<.05].copy();consistent=g[g.directionally_consistent].copy()
hits.to_csv(OUT/'all_overall_significant_pathways.csv',index=False)
# Author-approved biological themes; explicit representative sets, not an exhaustive module assignment.
themes=[
 ('Up','Axon guidance and synaptic organisation',[
 ('GOBP_NEURON_PROJECTION_GUIDANCE','Neuron projection guidance'),
 ('GOBP_AXON_DEVELOPMENT','Axon development'),
 ('GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY','Regulation of synapse structure or activity')]),
 ('Up','Polycomb-associated developmental signatures',[
 ('BENPORATH_PRC2_TARGETS','PRC2 target signature'),
 ('MIKKELSEN_MCV6_HCP_WITH_H3K27ME3','H3K27me3-associated promoters (MCV6)'),
 ('MIKKELSEN_NPC_HCP_WITH_H3K27ME3','H3K27me3-associated promoters (NPC)')]),
 ('Up','N-glycan processing',[
 ('GOBP_N_GLYCAN_PROCESSING','N-glycan processing')]),
 ('Up','Sterol / cholesterol responses',[
 ('GOBP_RESPONSE_TO_STEROL','Response to sterol'),
 ('GOBP_CELLULAR_RESPONSE_TO_CHOLESTEROL','Cellular response to cholesterol')]),
 ('Down','Interferon and cytokine responses',[
 ('HALLMARK_INTERFERON_ALPHA_RESPONSE','Interferon-alpha response'),
 ('HALLMARK_INTERFERON_GAMMA_RESPONSE','Interferon-gamma response'),
 ('REACTOME_SIGNALING_BY_INTERLEUKINS','Signalling by interleukins')]),
 ('Down','Extracellular matrix and adhesion',[
 ('REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','Extracellular matrix organisation'),
 ('REACTOME_COLLAGEN_DEGRADATION','Collagen degradation'),
 ('REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS','Integrin cell-surface interactions')]),
 ('Down','Migration, cytoskeleton and wound response',[
 ('GOBP_ACTIN_FILAMENT_ORGANIZATION','Actin filament organisation'),
 ('GOBP_POSITIVE_REGULATION_OF_LOCOMOTION','Positive regulation of locomotion'),
 ('GOBP_WOUND_HEALING','Wound healing')]),
 ('Down','Growth-factor signalling',[
 ('GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE','Positive regulation of MAPK cascade'),
 ('REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES','Receptor tyrosine kinase signalling'),
 ('GOBP_POSITIVE_REGULATION_OF_PHOSPHATIDYLINOSITOL_3_KINASE_PROTEIN_KINASE_B_SIGNAL_TRANSDUCTION','Positive regulation of PI3K-AKT signalling')])]
rows=[]
for direction,theme,terms in themes:
 for pathway,label in terms:
  q=g.loc[g.pathway==pathway];assert len(q)==1
  row=q.iloc[0].to_dict();assert row['FDR']<.05 and row['directionally_consistent']
  assert (row['NES']>0)==(direction=='Up')
  row.update(theme=theme,display_direction=direction,display_label=label);rows.append(row)
chosen=pd.DataFrame(rows)
chosen["node_id"]=np.arange(1,len(chosen)+1)
ft=tf[tf.directionally_consistent].sort_values('FDR').head(12).sort_values('score')
def nice(name,width=47):
 for prefix in ['HALLMARK_','REACTOME_','GOBP_','GOCC_','GOMF_']:name=name.removeprefix(prefix)
 name=name.replace('_',' ')
 if len(name)>2*width-5:name=name[:2*width-8]+'...'
 lines=textwrap.wrap(name,width)
 if len(lines)>2:lines=[lines[0],lines[1][:width-3]+'...']
 return '\n'.join(lines)
def heat(ax,frame,vals,ps,names,vmax=None):
 mat=frame[vals].to_numpy();p=frame[ps].to_numpy();mx=vmax or max(np.nanmax(abs(mat)),1)
 im=ax.imshow(mat,aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=-mx,vcenter=0,vmax=mx))
 ax.set_yticks(range(len(frame)),names,fontsize=6.4);ax.set_xticks(range(3),['Overall','August','May'],fontsize=7);ax.tick_params(length=0)
 for i,j in zip(*np.where(p<.05)):
  ax.text(j,i,'*',ha='center',va='center',fontsize=9,color='white' if abs(mat[i,j])>.6*mx else '#222222')
 for i,v in enumerate(p[:,0]):
  if np.isfinite(v):ax.text(3.05,i,f'{v:.1e}',fontsize=6,va='center',clip_on=False)
 ax.text(3.05,-1,'Overall\nFDR',fontsize=6,clip_on=False)
 return im
stem='ch03_main_ht1080_gsea_themes_network'
fig=plt.figure(figsize=(183/25.4,320/25.4))
def grouped(direction):
 rows=[]
 for theme,q in chosen[chosen.display_direction==direction].groupby('theme',sort=False):
  row={c:np.nan for c in chosen.columns};row['display_label']=theme;row['is_header']=True;rows.append(row)
  for record in q.to_dict('records'):record['is_header']=False;rows.append(record)
 return pd.DataFrame(rows)
up=grouped('Up');down=grouped('Down')
a_up=fig.add_axes([.56,.675,.265,.235]);a_down=fig.add_axes([.56,.345,.265,.285])
fig.text(.11,.965,'GSEA across all MSigDB collections',fontsize=9,weight='semibold')
fig.text(.11,.925,'Up in Y-depleted',fontsize=8,weight='semibold',color='#E64B35')
fig.text(.11,.648,'Down in Y-depleted',fontsize=8,weight='semibold',color='#3C5488')
shared_max=max(chosen[['NES','AUG_2025_NES','MAY_2026_NES']].abs().max().max(),1)
for ax,q in [(a_up,up),(a_down,down)]:
 im=heat(ax,q,['NES','AUG_2025_NES','MAY_2026_NES'],['FDR','AUG_2025_FDR','MAY_2026_FDR'],[label if header else f'{int(node)}. {label}' for label,header,node in zip(q.display_label,q.is_header,q.node_id)],shared_max)
 for tick,header in zip(ax.get_yticklabels(),q.is_header):
  if header:tick.set_weight('semibold');tick.set_color('#222222')
cb=fig.colorbar(im,cax=fig.add_axes([.95,.345,.012,.565]));cb.ax.set_title('NES',fontsize=6.5);cb.ax.tick_params(labelsize=6)
# Retain only the GSEA heatmaps and crop the former network area.
for ax in fig.axes:
 pos=ax.get_position();ax.set_position([pos.x0,(pos.y0-.30)/.70,pos.width,pos.height/.70])
for label in fig.texts:
 x,y=label.get_position();label.set_position((x,(y-.30)/.70))
fig.set_size_inches(183/25.4,224/25.4)
fig.text(.11,.015,'* Global FDR < 0.05 in the indicated comparison',fontsize=6.5)
exports=[]
for ext in ['pdf','svg','png']:
 p=FIG/f'{stem}.{ext}';stage=FIG/f'{stem}.staging.{ext}';fig.savefig(stage,dpi=600);stage.replace(p)
 exports.append(dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
fig.savefig(FIG/'previews'/f'{stem}_preview.png',dpi=150);plt.close(fig)
source=FIG/'source_data/ht1080_rna_all_collections';source.mkdir(parents=True,exist_ok=True)
chosen.to_csv(source/'main_pathways.csv',index=False);ft.to_csv(source/'main_TFs.csv',index=False)
g.to_csv(source/'all_overall_pathways.csv',index=False)

atlas='ch03_supp_ht1080_gsea_all_collections_atlas'
hits['collection_order']=hits.collection.map({c:i for i,c in enumerate(order)})
hits=hits.sort_values(['collection_order','FDR','pathway'])
pages=[];vmax=max(hits[['NES','AUG_2025_NES','MAY_2026_NES']].abs().max().max(),1)
if '--reuse-atlas' in sys.argv:
 pages=pd.read_csv(source/'atlas_page_index.csv').to_dict('records')
else:
 with PdfPages(FIG/f'{atlas}.pdf') as pdf:
  for page,start in enumerate(range(0,len(hits),24),1):
   q=hits.iloc[start:start+24]
   fig=plt.figure(figsize=(210/25.4,260/25.4));ax=fig.add_axes([.59,.10,.245,.79])
   im=heat(ax,q,['NES','AUG_2025_NES','MAY_2026_NES'],['FDR','AUG_2025_FDR','MAY_2026_FDR'],[f'{r.collection} | '+nice(r.pathway,56) for r in q.itertuples()],vmax)
   fig.text(.06,.97,f'All-collection GSEA | significant-set atlas | page {page}',fontsize=10,weight='semibold',va='top')
   fig.text(.06,.935,f'Overall global FDR < 0.05; rows {start+1}–{start+len(q)} of {len(hits)}',fontsize=8)
   cb=fig.colorbar(im,cax=fig.add_axes([.95,.10,.012,.79]));cb.ax.set_title('NES',fontsize=7);cb.ax.tick_params(labelsize=6)
   fig.text(.06,.04,'* Global FDR < 0.05 in the indicated contrast. Shared colour scale across all pages.',fontsize=7)
   pdf.savefig(fig)
   if page==1:fig.savefig(FIG/'previews'/f'{atlas}_preview.png',dpi=150)
   plt.close(fig)
   pages.extend(dict(page=page,pathway=r.pathway,collection=r.collection,FDR=r.FDR) for r in q.itertuples())
 pd.DataFrame(pages).to_csv(source/'atlas_page_index.csv',index=False)
exports.append(dict(file=f'{atlas}.pdf',sha256=hashlib.sha256((FIG/f'{atlas}.pdf').read_bytes()).hexdigest()))
# Original three-collection results remain intact; compare overlapping terms explicitly.
old=pd.read_csv(BASE/'tables/programme_consistency.csv')
compare=old[['pathway','FDR','NES','directionally_consistent']].merge(g[['pathway','FDR','NES','directionally_consistent']],on='pathway',suffixes=('_three_collections','_all_collections'))
compare.to_csv(OUT/'three_vs_all_collections.csv',index=False)
stats={'tested_overall':len(g),'significant_overall':len(hits),'directionally_consistent':len(consistent),
 'both_periods_significant':int(g.both_periods_FDR05.sum()),'higher':int((hits.NES>0).sum()),'lower':int((hits.NES<0).sum()),
 'atlas_pages':len(set(p['page'] for p in pages)),'displayed_main_pathways':len(chosen),
 'selection':'Explicit representative pathways within eight author-approved biological themes; overall global FDR < 0.05 and matching NES directions across periods. Themes are interpretive groupings, not recomputed network modules. Heatmap-only display; the prior gene-set network has been removed.',
 'exports':exports}
(source/'provenance.json').write_text(json.dumps(stats,indent=2));(OUT/'summary.json').write_text(json.dumps(stats,indent=2))
def table(d):return '<div class="scroll">'+d.to_html(index=False,border=0,float_format=lambda v:f'{v:.4g}')+'</div>'
overview=pd.DataFrame(summary);main=overview[overview.contrast=='overall']
preview=FIG/'previews'/f'{stem}_preview.png'
body=f'''<h1>HT1080: GSEA across all MSigDB collections</h1><p>{len(g):,} eligible sets; {len(hits):,} overall global-FDR discoveries; {len(consistent):,} with matching directions across periods; {int(g.both_periods_FDR05.sum()):,} also significant in both periods.</p>
<p>All available frozen MSigDB 2026.1.Hs H and C1–C9 collections were screened using current moderated t rankings and the existing autosomal/non-PAR X background. Set size: 15–500 represented genes. BH correction spans every eligible set across collections separately per contrast. Within-collection FDR is secondary. Current rankings differ from the older descriptive all-collection workflow.</p>
<p>The figure displays representative pathways in eight agreed biological themes. Source signatures describe their originating experiments or cell types; enrichment does not establish the same mechanism or cell identity in HT1080. The prior three-collection CAMERA sensitivity does not test this expanded family. GSEA remains exploratory.</p>
<p>Numerical accounting: {int(allg.p_at_floor.sum()):,} tests across all contrasts reached the 1e-10 P-value floor; {int(allg.log2err.isna().sum()):,} have no Monte Carlo error estimate. Floor values are computational bounds, not exact probabilities. Warnings and full precision are retained with each contrast.</p>
<p><a href="GSEA_all_contrasts.csv">All tests, all five contrasts</a> | <a href="all_overall_significant_pathways.csv">All overall discoveries</a> | <a href="three_vs_all_collections.csv">Comparison with three-collection screen</a> | <a href="../report.html">Previous programme report</a></p>
<img style="max-width:850px;width:100%" src="data:image/png;base64,{base64.b64encode(preview.read_bytes()).decode()}">
<h2>Collection coverage</h2>{table(main)}<h2>All contrasts</h2>{table(overview)}
<h2>All overall significant sets</h2><input id="search" placeholder="Filter pathway or collection"><p>FDR columns below use the full cross-collection family.</p>'''
body+='<p><a href="'+Path(os.path.relpath(FIG/f'{stem}.pdf',OUT)).as_posix()+'">Representative programme/regulator figure</a> | <a href="'+Path(os.path.relpath(FIG/f'{atlas}.pdf',OUT)).as_posix()+'">Complete significant-set atlas</a></p>'
body+=table(hits[['pathway','collection','subcollection','NES','FDR','AUG_2025_NES','AUG_2025_FDR','MAY_2026_NES','MAY_2026_FDR','same_direction','both_periods_FDR05']])
body+='''<script>document.getElementById('search').addEventListener('input',function(){const q=this.value.toLowerCase();document.querySelectorAll('tbody tr').forEach(r=>r.style.display=r.textContent.toLowerCase().includes(q)?'':'none');});</script>'''
style='body{font-family:Arial,sans-serif;max-width:1250px;margin:30px auto;padding:0 20px;line-height:1.5}.scroll{overflow:auto;max-height:600px}td,th{padding:6px;border-bottom:1px solid #ddd;font-size:12px}th{background:#edf1f7}input{padding:10px;width:95%}h1,h2{color:#3C5488}'
(OUT/'report.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>All-collection GSEA</title><style>'+style+'</style></head><body>'+body+'</body></html>',encoding='utf-8')
(OUT/'validation.json').write_text(json.dumps({'all_contrast_BH_verified':True,'original_input_hashes_verified':True,'finite_tests':int(allg.pval.notna().sum()),'tests_total':len(allg),'pvalues_at_floor':int(allg.p_at_floor.sum())},indent=2))
print(json.dumps(stats,indent=2))
