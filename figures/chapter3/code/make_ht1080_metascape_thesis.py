from pathlib import Path
import sys,json,hashlib,textwrap,xml.etree.ElementTree as ET
from collections import Counter
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from adjustText import adjust_text
FIG=Path(__file__).resolve().parents[1];THESIS=FIG.parents[1];PROJECT=THESIS.parent
MS=PROJECT/'RNAseq/analysis/HT1080/reannotate_v105/metascape'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
OUT=FIG/'source_data/ht1080_metascape_latest';OUT.mkdir(parents=True,exist_ok=True)
D=pd.read_csv(MS/'pathway_network_combinedrun_nodes.csv')
NAMES={'up':{1:'Polycomb-associated promoter signatures',2:'Neurite and synapse morphogenesis',3:'CDH1-associated signatures',4:'Embryonic patterning and organogenesis',5:'Ion channels and neuronal excitability',6:'GPCR / neuroactive ligand signalling',7:'Secreted ligands and growth factors',8:'Nerve development and chemotaxis',9:'Synaptic transmission and cognition'},'down':{1:'Stromal ECM and collagen turnover',2:'Lymphoma / KRAS-associated signatures',3:'RTK / PI3K / MAPK signalling',4:'Angiogenesis and vascular sprouting',5:'CDH1 / beta-catenin signatures',6:'TGF-beta1 response signatures',7:'HDAC1/2 perturbation signatures'}}
PALETTE=['#3C5488','#E64B35','#00A087','#F39B7F','#4DBBD5','#8491B4','#91D1C2','#7E6148']
ns={'x':'http://www.cs.rpi.edu/XGMML'}
exports=[];audit=[];figures=[]
def wrap(x,n):return '\n'.join(textwrap.wrap(str(x),n))
def short(s):
 s=s.replace('MIKKELSEN ','').replace(' HCP WITH H3K27ME3',' HCP H3K27me3').replace('_',' ')
 return s if len(s)<29 else s[:26]+'...'
for direction in ['up','down']:
 folder=MS/f'BOTH_BATCH_LOY_sgCTRL_background_0.05_0.5LFC_{direction.upper()}'
 xp=folder/'Enrichment_GO/GONetwork.xgmml';root=ET.parse(xp).getroot()
 lookup={};raw=[]
 for node in root.findall('x:node',ns):
  att={a.attrib['name']:a.attrib.get('value') for a in node.findall('x:att',ns)}
  lookup[node.attrib['id']]=att['canonicalName'];raw.append(att)
 q=D[D.direction==direction].copy();G=nx.Graph();G.add_nodes_from(q.name)
 edges=[]
 for edge in root.findall('x:edge',ns):
  u=lookup[edge.attrib['source']];v=lookup[edge.attrib['target']]
  if u in G and v in G and u!=v:G.add_edge(u,v)
 assert set(G)==set(q.name)
 rawdf=pd.DataFrame(raw).set_index('canonicalName')
 assert np.allclose(q.mlogp.to_numpy(),-rawdf.loc[q.name,'LogP'].astype(float).to_numpy())
 top5=set(q.nlargest(5,'mlogp').name)
 groups=[]
 for c,g in q.groupby('community'):
  best=g.nlargest(2,'mlogp');cnt=Counter(g.hits.str.split('|').explode());genes=[x for x,n in sorted(cnt.items(),key=lambda z:(-z[1],z[0]))[:4]]
  groups.append(dict(community=int(c),theme=NAMES[direction][c],terms=len(g),best_logP=g.mlogp.max(),representative_pathways='; '.join(short(s) for s in best.desc),recurring_hit_genes=', '.join(genes),priority=bool(set(g.name)&top5)))
 summary=pd.DataFrame(groups).sort_values(['priority','terms','community'],ascending=[False,False,True])
 major=summary[(summary.terms>=4)|summary.priority].head(8).community.tolist();colors={c:PALETTE[i] for i,c in enumerate(major)}
 summary['color']=summary.community.map(colors).fillna('#BBBBBB');q['color']=q.community.map(colors).fillna('#BBBBBB');q['theme']=q.community.map(NAMES[direction])
 comps=sorted(nx.connected_components(G),key=lambda c:(-len(c),sorted(c)[0]));pos={}
 for i,comp in enumerate(comps):
  sub=G.subgraph(sorted(comp));xy=nx.kamada_kawai_layout(sub,weight=None) if len(comp)>1 else {next(iter(comp)):np.zeros(2)}
  if i==0:center=np.array([0.,0.]);scale=1.
  else:center=np.array([1.65,1.0-(i-1)*.62]);scale=.19 if len(comp)<6 else .28
  for node,pt in xy.items():pos[node]=pt*scale+center
 q['x']=q.name.map(lambda n:pos[n][0]);q['y']=q.name.map(lambda n:pos[n][1])
 fig=plt.figure(figsize=(183/25.4,230/25.4));ax=fig.add_axes([.05,.34,.90,.57]);ax.set_axis_off()
 nx.draw_networkx_edges(G,pos,ax=ax,width=.5,edge_color='#B9B9B9',alpha=.55)
 sizes=12+q.mlogp.to_numpy()*4.0
 nx.draw_networkx_nodes(G,pos,ax=ax,nodelist=q.name.tolist(),node_size=sizes,node_color=q.color.tolist(),edgecolors='#555555',linewidths=.4)
 texts=[]
 for c in major:
  z=q[q.community==c];x=z.x.mean();y=z.y.mean();texts.append(ax.text(x,y,wrap(NAMES[direction][c],23),fontsize=7,weight='semibold',ha='center',va='center',bbox=dict(facecolor='white',alpha=.72,edgecolor='none',pad=1)))
 np.random.seed(42)
 adjust_text(texts,x=q.x.to_numpy(),y=q.y.to_numpy(),ax=ax,arrowprops=dict(arrowstyle='-',color='#777777',lw=.5),iter_lim=250)
 ax.margins(.12)
 fig.text(.055,.957,('Higher' if direction=='up' else 'Lower')+' expression in Y-depleted HT1080',fontsize=10,weight='semibold')
 fig.text(.055,.93,f'Metascape | {len(q)} terms, {G.number_of_edges()} connections',fontsize=8,color='#555555')
 handles=[ax.scatter([],[],s=12+v*4,c='#BBBBBB',edgecolors='#555555',linewidths=.4,label=str(v)) for v in [5,10,15]]
 fig.legend(handles=handles,title=r'$-\log_{10}P$',loc='lower center',bbox_to_anchor=(.5,.30),ncol=3,frameon=False,fontsize=7,title_fontsize=7)
 tx=fig.add_axes([.04,.045,.93,.24]);tx.set_axis_off();tx.set_xlim(0,1);tx.set_ylim(0,1)
 xs=[.025,.305,.385,.45,.775];headers=['Theme','Terms','Best\n-log10 P','Representative pathways','Recurring hit genes']
 for x,h in zip(xs,headers):tx.text(x,.99,h,fontsize=6.2,weight='semibold',va='top')
 tx.plot([0,1],[.87,.87],color='#333333',lw=.5)
 heights=.80/len(summary)
 for i,row in enumerate(summary.itertuples()):
  y=.83-i*heights;tx.scatter(.008,y-.018,s=19,color=row.color)
  vals=[wrap(row.theme,26),str(row.terms),f'{row.best_logP:.1f}',row.representative_pathways.replace('; ','\n'),wrap(row.recurring_hit_genes,23)]
  for x,v in zip(xs,vals):tx.text(x,y,v,fontsize=6.0,va='top')
 fig.text(.055,.019,'Grey: smaller communities. Gene lists show recurring enrichment hits, not GSEA leading edges.',fontsize=5.8,color='#555555')
 stem=f'ch03_main_ht1080_metascape_{direction}_network'
 for ext in ['pdf','svg','png']:
  target=FIG/f'{stem}.{ext}';fig.savefig(target,dpi=600);exports.append(dict(file=target.name,sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
 fig.savefig(FIG/'previews'/f'{stem}_preview.png',dpi=150);figures.append(fig)
 q.to_csv(OUT/f'{direction}_nodes.csv',index=False);pd.DataFrame(G.edges,columns=['source','target']).to_csv(OUT/f'{direction}_edges.csv',index=False);summary.to_csv(OUT/f'{direction}_themes.csv',index=False)
 audit.append(dict(direction=direction,nodes=len(G),edges=G.number_of_edges(),omitted_isolates=len(raw)-len(G),source=str(xp),sha256=hashlib.sha256(xp.read_bytes()).hexdigest(),export_date=root.find('.//{http://purl.org/dc/elements/1.1/}date').text))
with PdfPages(FIG/'ch03_supp_ht1080_metascape_networks_combined.pdf') as pdf:
 for fig in figures:pdf.savefig(fig);plt.close(fig)
(OUT/'provenance.json').write_text(json.dumps(dict(runs=audit,exports=exports,community_source=str(MS/'pathway_network_combinedrun_nodes.csv'),selection='Existing September 8 combined DESeq2 batch+treatment Metascape export; FDR<0.05, absolute log2FC>0.5. Not current GSEA.',layout='Original connected nodes and communities retained; exported edges unchanged; Kamada-Kawai within components, separated component packing.'),indent=2))
print(json.dumps(audit,indent=2))
