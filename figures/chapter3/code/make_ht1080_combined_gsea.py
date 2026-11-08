"""Combine saved primary and all-MSigDB HT1080 GSEA displays without reanalysis."""
from pathlib import Path
import sys,json,hashlib,textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
import pandas as pd
import numpy as np
FIG=Path(__file__).resolve().parents[1]; ROOT=FIG.parents[1]
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
paths=[FIG/'source_data/ht1080_primary_gsea/displayed_pathways.csv',FIG/'source_data/ht1080_rna_all_collections/main_pathways.csv',FIG/'source_data/ht1080_primary_gsea/discovery_counts.csv',FIG/'source_data/ht1080_rna_all_collections/network_nodes.csv',FIG/'source_data/ht1080_combined_gsea/CAMERA_overall_sensitivity.csv']
hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
a=pd.read_csv(paths[0]);allb=pd.read_csv(paths[1]);nodes=pd.read_csv(paths[3])
camera=pd.read_csv(paths[4]);camera=camera[camera.correlation_mode.eq('fixed_0.01')].copy()
assert len(camera)==4868 and camera.pathway.is_unique and camera.FDR.lt(.05).sum()==17
camera=camera[['pathway','Direction','PValue','FDR']].rename(columns={'Direction':'CAMERA_direction','PValue':'CAMERA_PValue','FDR':'CAMERA_FDR'})
a=a.merge(camera,on='pathway',how='left',sort=False,validate='one_to_one')
assert a.CAMERA_FDR.notna().all()
assert ((a.NES.gt(0)&a.CAMERA_direction.eq('Up'))|(a.NES.lt(0)&a.CAMERA_direction.eq('Down'))).all()
existing=allb[~allb.pathway.isin(a.pathway)]
eligible=nodes[~nodes.pathway.isin(a.pathway)&nodes.FDR.lt(.05)&nodes.same_direction.eq(True)]
expanded=[]
for theme,g in eligible.groupby('theme',sort=False):
 keep=g[g.pathway.isin(existing.pathway)]
 extra=g[~g.pathway.isin(keep.pathway)].sort_values(['FDR','pathway']).head(max(0,4-len(keep)))
 expanded.append(pd.concat([keep,extra]))
b=pd.concat(expanded,ignore_index=True)
assert len(a)==15 and len(b)==28 and not set(a.pathway)&set(b.pathway)
assert set(existing.pathway)<=set(b.pathway)
cols=['NES','AUG_2025_NES','MAY_2026_NES'];fcols=['FDR','AUG_2025_FDR','MAY_2026_FDR']
cmap=LinearSegmentedColormap.from_list('NES',['#3C5488','#FAFAFA','#E64B35']);limit=np.ceil(max(abs(a[cols].to_numpy()).max(),abs(b[cols].to_numpy()).max())*10)/10
fig=plt.figure(figsize=(183/25.4,270/25.4))
short={'GOBP_POSITIVE_REGULATION_OF_TYPE_II_INTERFERON_PRODUCTION':'Positive regulation of IFN-gamma production','GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS':'CNS projection neuron axonogenesis'}
def label(row,primary):
 name=short.get(row.pathway,row.pathway.removeprefix('HALLMARK_').removeprefix('REACTOME_').removeprefix('GOBP_').replace('_',' ').capitalize()) if primary else row.display_label
 name=name.replace('Myc targets v1','MYC targets V1').replace('Esr mediated signaling','ESR-mediated signalling').replace('mapk cascade','MAPK cascade')
 code={'Hallmark':'H','Reactome':'R','GO Biological Process':'GO'}.get(getattr(row,'analysis_collection',''),row.collection)
 if not primary:
  name=name.replace('\\n',' ').replace('\n',' ')
  if name.isupper():name=name.capitalize()
  name=name.replace('Central nervous system neuron axonogenesis','CNS neuron axonogenesis').replace('Mikkelsen mef hcp with h3k27me3','H3K27me3-associated promoters (MEF)').replace('Cell surface receptor signaling pathway via jak stat','JAK-STAT receptor signalling').replace('Phosphatidylinositol 3 kinase protein kinase b signal transduction','PI3K-AKT signal transduction')
 name=name.replace('Type ii','Type II').replace('type ii','type II')
 return '['+code+'] '+textwrap.fill(name,width=62)
# Preserve the original interpretive themes, with spacing and labelled headings.
theme_map=allb.set_index('pathway').theme.to_dict()
a['theme']=a.pathway.map(theme_map)
a.loc[a.pathway.eq('GOBP_POSITIVE_REGULATION_OF_TYPE_II_INTERFERON_PRODUCTION'),'theme']='Interferon and cytokine responses'
a.loc[a.pathway.eq('HALLMARK_INFLAMMATORY_RESPONSE'),'theme']='Interferon and cytokine responses'
a.loc[a.NES.gt(0),'theme']='Axon guidance and synaptic organisation'
a['theme']=a['theme'].fillna('Other negatively enriched sets')
theme_order=list(dict.fromkeys(allb.theme))+['Other negatively enriched sets']
for data,rect,letter,title,subtitle,primary in [(a,[.54,.635,.255,.255],'a','Primary gene-set analysis','Hallmark, Reactome and GO Biological Process | 4,868 sets tested',True),(b,[.54,.105,.255,.425],'b','Additional sets from the all-MSigDB analysis','All collections | 23,696 sets tested | sets shown in a omitted',False)]:
 top=rect[1]+rect[3];fig.text(.035,top+.063,letter,fontsize=12,weight='semibold');fig.text(.08,top+.064,title,fontsize=10,weight='semibold');fig.text(.08,top+.042,subtitle,fontsize=7)
 if primary:fig.text(.08,top+.023,'220 directionally consistent: 209 lower, 11 higher; 55 significant in both periods',fontsize=6.8,color='#555555')
 rows=[];heads=[];matrix=[];cursor=0
 for theme in theme_order:
  sub=data[data.theme.eq(theme)]
  if sub.empty:continue
  heads.append((cursor,theme,sub.NES.mean()>0));matrix.append([np.nan]*3);cursor+=1
  for row in sub.itertuples():
   rows.append((cursor,row));matrix.append([getattr(row,c) for c in cols]);cursor+=1
 ax=fig.add_axes(rect);vals=np.array(matrix);im=ax.imshow(vals,aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit))
 ax.set_yticks([i for i,r in rows],[label(r,primary) for i,r in rows],fontsize=6.1)
 ax.set_xticks([0,1,2],['Overall','August','May'],fontsize=7);ax.tick_params(length=0);ax.spines[:].set_visible(False)
 for i,theme,positive in heads:
  ax.text(-.57,i,theme,ha='right',va='center',fontsize=6.5,weight='semibold',color='#B5422B' if positive else '#3C5488',clip_on=False)
 for i,row in rows:
  for j,c in enumerate(fcols):
   if getattr(row,c)<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=9,color='white' if abs(vals[i,j])>limit*.55 else '#242424')
  if primary:
   ax.text(2.97,i,f'{row.FDR:.1e}',ha='center',va='center',fontsize=6.5)
   ax.text(4.02,i,f'{row.CAMERA_FDR:.2g}',ha='center',va='center',fontsize=6.5,
           weight='semibold' if row.CAMERA_FDR<.05 else 'normal')
  else:
   ax.text(2.7,i,f'{row.FDR:.1e}',va='center',fontsize=6)
 if primary:
  ax.text(2.97,-.9,'fgsea\nFDR',ha='center',fontsize=6.5,va='bottom')
  ax.text(4.02,-.9,'CAMERA\nFDR',ha='center',fontsize=6.5,va='bottom')
 else:
  ax.text(2.7,-.9,'Overall fgsea FDR',fontsize=6.5,va='bottom')
cax=fig.add_axes([.34,.055,.35,.012]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.ax.tick_params(labelsize=7);cb.ax.set_title('Normalised enrichment score (NES)',fontsize=7,pad=6)
fig.text(.08,.025,'* fgsea FDR < 0.05 in the indicated contrast. Bold CAMERA values: overall FDR < 0.05.',fontsize=6.5)
fig.text(.08,.01,'H, Hallmark; R, Reactome; GO, GO Biological Process; C2, curated; C5, ontology.',fontsize=6.5)
stem='ch03_main_ht1080_combined_gsea'
for ext in ['pdf','svg','png']:fig.savefig(FIG/(stem+'.'+ext),dpi=600)
fig.savefig(FIG/'previews'/(stem+'_preview.png'),dpi=170)
out=FIG/'source_data/ht1080_combined_gsea';out.mkdir(exist_ok=True)
a.to_csv(out/'panel_a_primary.csv',index=False);b.to_csv(out/'panel_b_additional.csv',index=False)
assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in hashes.items())
(out/'provenance.json').write_text(json.dumps({'source_sha256':hashes,'primary_rows':len(a),'additional_rows':len(b),'duplicate_sets_removed':sorted(set(a.pathway)&set(allb.pathway)),'values':'Original saved NES and FDR values retained; no analysis rerun','shared_NES_limits':[-float(limit),float(limit)],'font':'Open Sans','size_mm':[183,270], 'expanded_selection':'Retain previous panel b rows; add lowest-overall-FDR eligible saved themed sets up to four per theme, excluding primary display', 'grouping':'Original all-MSigDB biological themes restored; primary-only unmatched sets labelled Other negatively enriched sets', 'CAMERA':{'added':'2026-10-09','source':'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/programme_analysis_v1/tables/CAMERA_overall_sensitivity.csv, relative to projects','correlation':'fixed 0.01','family':4868,'full_family_FDR05':17,'displayed_FDR05':int(a.CAMERA_FDR.lt(.05).sum()),'display':'Overall CAMERA FDR beside overall fgsea FDR in panel a; bold denotes CAMERA FDR <0.05; all displayed directions agree','join':'exact pathway identifier; one-to-one; no missing values','panel_b':'No CAMERA added: its larger all-MSigDB family is distinct','analysis_refitted':False}},indent=2),encoding='utf-8')
print('Combined figure rendered: 15 primary and 28 additional sets; saved inputs unchanged.')
