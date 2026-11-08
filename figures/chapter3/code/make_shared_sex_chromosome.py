"""Matched descriptive RNA/lpWGS panels; frozen CPM and effects, no DE refitting."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text
F=Path(__file__).resolve().parents[1]; ROOT=F.parents[1]; P=ROOT.parent
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
S=F/'source_data/shared_sex_chromosome'
OUT=F/'drafts/matched_sex_chromosomes_20260926'; OUT.mkdir(parents=True,exist_ok=True)
BLUE,RED,TEAL,BROWN='#3C5488','#E64B35','#00A087','#7E6148'
COL={'PAR1':'#4DBBD5','Distal Yp':TEAL,'Remaining MSY':BROWN,'PAR2':'#8491B4'}
REG=list(COL); END=27.449937; BREAK=4.8
posfile=P/'RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE/outputs/tables/Y_PAR_gene_table.csv'
pos=pd.read_csv(posfile).drop_duplicates('gene_id').set_index('gene_id').t2t_mid_Mb
pairs=[('RPS4X','RPS4Y1'),('ZFX','ZFY'),('PRKX','PRKY'),('USP9X','USP9Y'),('DDX3X','DDX3Y'),('KDM6A','UTY'),('TMSB4X','TMSB4Y'),('NLGN4X','NLGN4Y'),('KDM5C','KDM5D'),('EIF1AX','EIF1AY')]
WGS={'HT1080':('HT1080-sgCTRL-D0','HT1080-LOY-D0'),'A549':('A549-CTRL2','A549-LOY'),'DU145':('LNCaP-sgCTRL','LNCaP-LOY-2X')}
all_data={}
for line in WGS:
 e=pd.read_csv(S/f'{line}_effects.csv').set_index('gene_id'); cp=pd.read_csv(S/f'{line}_CPM.csv').set_index('gene_id'); meta=pd.read_csv(S/f'{line}_samples.csv')
 assert cp.columns.tolist()==meta.sample_key.tolist()
 e['pos']=pos.reindex(e.index)
 e['region']=np.where(e.PAR & e.chrom.eq('chrX'),np.where(e.start<5e6,'PAR1','PAR2'),np.where(e.chrom.eq('chrY') & ~e.PAR,np.where(e.pos<BREAK,'Distal Yp','Remaining MSY'),'X'))
 periods=sorted(meta.batch.unique()); rows=[]; sums=[]
 for batch in periods:
  m=meta[meta.batch==batch]; ctrl=m[m.treatment=='CTRL'].sample_key; loy=m[m.treatment=='LOY'].sample_key
  d=e.copy(); d['CTRL']=cp[ctrl].mean(axis=1);d['LOY']=cp[loy].mean(axis=1);d['percent']=100*d.LOY/d.CTRL;d['batch']=batch
  if line=='HT1080':
   ef=pd.read_csv(P/'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE/outputs/tables/HT1080_LOY_vs_sgCTRL_all_gene_effects.csv').set_index('gene_id')
   d['effect']=ef[batch].reindex(d.index)
  else: d['effect']=d.mean_difference
  rows.append(d.reset_index())
  for r in REG:
   ids=d.index[d.region.eq(r)]
   for sample in m.itertuples(): sums.append(dict(batch=batch,region=r,sample=sample.sample_key,treatment=sample.treatment,CPM=cp.loc[ids,sample.sample_key].sum()))
 data=pd.concat(rows); sums=pd.DataFrame(sums)
 all_data[line]=(data,sums,meta,cp)
# Same named Y/PAR gene order across lines; omit genes below display thresholds in every line.
y_order={}; par_order={}
for line,(d,_,_,_) in all_data.items():
 for r in d.itertuples():
  if str(r.gene_name).startswith('ENSG'): continue
  if r.region in ['Distal Yp','Remaining MSY'] and r.CTRL>=1 and np.isfinite(r.pos) and r.pos<END: y_order[r.gene_name]=r.pos
  if r.region in ['PAR1','PAR2'] and r.CTRL>=5: par_order[r.gene_name]=(0 if r.region=='PAR1' else 1,r.start)
y_order=sorted(y_order,key=y_order.get); par_order=sorted(par_order,key=par_order.get)
for line,(data,sums,meta,cp) in all_data.items():
 if '--line' in sys.argv and line != sys.argv[sys.argv.index('--line')+1]: continue
 if line in ['DU145','A549']:
  fig=plt.figure(figsize=(7.2,6.5))
  positions=[[.09,.73,.87,.22],[.09,.395,.24,.19],[.41,.395,.25,.19],[.75,.395,.21,.19],[.09,.09,.40,.17]]
 else:
  fig=plt.figure(figsize=(7.2,7.4))
  positions=[[.09,.74,.87,.20],[.09,.415,.24,.18],[.41,.415,.25,.18],[.75,.415,.21,.18],[.09,.095,.87,.17]]
 axes=[fig.add_axes(v) for v in positions]; a,b,c,d,e=axes
 for ax,letter in zip(axes,'abcde'): ax.text(-.085,1.07,letter,transform=ax.transAxes,weight='semibold',fontsize=12)
 for name,x0,x1 in [('PAR1',0,2.45832),('Distal Yp',2.45832,BREAK),('Remaining MSY',BREAK,END)]:
  a.axvspan(x0,x1,color=COL[name],alpha=.1,lw=0);a.text((x0+x1)/2,1.02,name,transform=a.get_xaxis_transform(),ha='center',fontsize=7,color='#555555')
 binsfile=F/f'source_data/{line.lower()}_lpwgs_confirmation/qdnaseq_bins.csv'; bins=pd.read_csv(binsfile)
 def segments(sample):
  q=bins[bins['sample'].eq(sample)&bins.chromosome.eq('chrY')&bins['use']&bins.log2_seg.notna()].sort_values('start').copy();q['cn']=2**q.log2_seg;q['run']=(q.cn.round(8)!=q.cn.round(8).shift()).cumsum()
  return q.groupby('run').agg(start=('start','min'),end=('end','max'),cn=('cn','first'))
 ctrl,loy=map(segments,WGS[line]); dna=[]
 for C in ctrl.itertuples():
  for L in loy.itertuples():
   st=max(C.start,L.start)/1e6; en=min(C.end,L.end,END*1e6)/1e6
   if line=='A549': st=max(st,2.45832)
   if st<en and C.cn>=.05:
    pct=100*max(L.cn,.001)/C.cn;dna.append(dict(start=st,end=en,percent=pct));a.plot([st,en],[max(pct,1)]*2,color='#333333',lw=1.3,ls='--' if pct<1 else '-')
 if line=='A549':
  # PAR1 segmentation disagrees with corrected bins; use the regional bin median.
  par=bins[bins.chromosome.eq('chrY') & bins['use'] & (bins.midpoint<2.45832)]
  med=par.groupby('sample').copynumber.median()
  pct=100*med[WGS[line][1]]/med[WGS[line][0]]
  a.plot([0,2.45832],[pct,pct],color='#333333',lw=1.3,ls=':')
  a.annotate(f'{100-pct:.0f}% lower\nDNA signal',xy=(1.22916,pct),xytext=(1.22916,12),ha='center',va='center',fontsize=6,color='#333333',arrowprops=dict(arrowstyle='-',color='#555555',lw=.6))
  dna.append(dict(start=0,end=2.45832,percent=pct))
  pd.DataFrame([dict(region='PAR1',CTRL=med[WGS[line][0]],LOY=med[WGS[line][1]],percent=pct,n_bins=len(par)//2,method='Median corrected bin signal; bins assigned by midpoint')]).to_csv(S/'A549/PAR1_bin_check.csv',index=False)
 periods=sorted(data.batch.unique()); handles=[]; pp=[]
 for j,batch in enumerate(periods):
  marker=['o','^'][j]; label=('Batch '+str(j+1)) if line=='HT1080' else 'RNA mean'
  q=data[data.batch==batch].set_index('gene_name'); q=q[~q.index.duplicated()]
  yy=q.reindex(y_order); yy=yy[yy.CTRL>=1]
  a.scatter(yy.pos,yy.percent.clip(lower=1),marker=marker,s=20,facecolor='white' if j==0 else RED,edgecolors=RED,lw=.7,zorder=4)
  handles.append(Line2D([],[],marker=marker,color=RED,mfc='white' if j==0 else RED,ls='',label=label,ms=4))
  # Regional mean abundances, with sampling periods offset and differentiated by hatching.
  for k,group in enumerate(['CTRL','LOY']):
   ss=sums[(sums.batch==batch)&sums.treatment.eq(group)].groupby('region').CPM.mean().reindex(REG)
   width=.36/len(periods); offset=(j*2+k-(2*len(periods)-1)/2)*width
   b.bar(np.arange(4)+offset,ss,width=width,color=BLUE if group=='CTRL' else RED,hatch='///' if j else None,edgecolor='white',lw=.2)
  shift=(j-(len(periods)-1)/2)*.25
  for ax,order,threshold in [(c,y_order,1),(d,par_order,5)]:
   qq=q.reindex(order); valid=qq.CTRL>=threshold; xx=np.arange(len(order))+shift
   col=[COL.get(r,BROWN) for r in qq.region]
   ax.scatter(xx[valid],qq.percent[valid].clip(lower=1),c=np.array(col)[valid],s=16,marker=marker,edgecolors='white',lw=.3,zorder=3)
  for idx,(X,Y) in enumerate(pairs):
   if X not in q.index or Y not in q.index or q.loc[Y,'CTRL']<5: continue
   xx=idx+shift; xv=q.loc[X,'effect']; yv=q.loc[Y,'effect'];e.plot([xx,xx],[xv,yv],c='#cccccc',lw=.7)
   e.scatter(xx,xv,c=BLUE,s=22,marker=marker,zorder=3);e.scatter(xx,yv,c=RED,s=22,marker=marker,zorder=3)
   pp.append(dict(batch=batch,X=X,Y=Y,X_effect=xv,Y_effect=yv))
 # Label positional means once, keeping each period's marks visible.
 avg=data[data.gene_name.isin(y_order)].groupby('gene_name').agg(pos=('pos','first'),percent=('percent','mean'),CTRL=('CTRL','mean'))
 texts=[a.text(r.pos,max(r.percent,1),name,fontsize=6) for name,r in avg.iterrows() if r.CTRL>=1]
 a.set(xlim=(-.5,END+.5),yscale='log',ylim=(.65,200),xlabel='CHM13v2.0 chromosome Y position (Mb)',ylabel='Expression / DNA (% sgCTRL)')
 a.set_yticks([1,10,100],['\u22641','10','100']);a.axhline(100,c='#999999',ls='--',lw=.6)
 np.random.seed(42);adjust_text(texts,ax=a,arrowprops=dict(arrowstyle='-',color='#aaaaaa',lw=.4))
 handles += [Line2D([],[],color='#333333',lw=1.3,label='lpWGS'),Line2D([],[],color='#333333',lw=1.3,ls='--',label='lpWGS \u22641%')]
 a.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.30),ncol=4,fontsize=6.5,frameon=False)
 for i,region in enumerate(REG):
  vals=[]
  for batch in periods:
   rr=sums[(sums.batch==batch)&sums.region.eq(region)].groupby('treatment').CPM.mean()
   vals.append(f'{100*rr.LOY/rr.CTRL:.0f}%')
  height=sums[sums.region.eq(region)].groupby(['batch','treatment']).CPM.mean().max()
  b.text(i,height*1.8,' / '.join(vals),ha='center',fontsize=5.6)
 b.set(yscale='log',ylim=(.1,1e5),ylabel='Summed CPM');b.set_xticks(range(4),['PAR1','Distal\nYp','Remaining\nMSY','PAR2'],fontsize=6)
 b.legend(handles=[plt.Rectangle((0,0),1,1,color=BLUE,label='sgCTRL'),plt.Rectangle((0,0),1,1,color=RED,label='Y-depleted')],fontsize=6,frameon=False)
 for ax,order in [(c,y_order),(d,par_order)]:
  ax.set(yscale='log',ylim=(.65,250),ylabel='% of sgCTRL',xlim=(-.6,len(order)-.4));ax.set_yticks([1,10,100],['\u22641','10','100']);ax.axhline(100,c='#999999',ls='--',lw=.6);ax.set_xticks(range(len(order)),order,rotation=90,fontsize=5.6)
 # Matched DU145/A549 layout: regional percentages and per-gene bars.
 if line in ['DU145','A549']:
  for ax in [b,c,d]: ax.clear()
  for ax,letter in [(b,'b'),(c,'c'),(d,'d')]:
   ax.text(-.085,1.07,letter,transform=ax.transAxes,weight='semibold',fontsize=12)
  regional=sums.groupby(['region','treatment']).CPM.mean().unstack().reindex(REG)
  regional['percent']=100*regional.LOY/regional.CTRL
  xx=np.arange(4)
  b.bar(xx-.18,[100]*4,width=.34,color=BLUE,label='sgCTRL')
  b.bar(xx+.18,regional.percent,width=.34,color=RED,label='Y-depleted')
  for i,value in enumerate(regional.percent):
   label=f'{value:.1f}%' if line=='A549' and value<10 else f'{value:.0f}%'
   b.text(i+.18,value+4,label,ha='center',fontsize=6)
  b.axhline(100,c='#999999',ls='--',lw=.6)
  b.set(ylim=(0,140),ylabel='Expression (% of sgCTRL)')
  b.set_yticks([0,50,100]);b.set_xticks(xx,['PAR1','Distal\nYp','Remaining\nMSY','PAR2'],fontsize=6)
  b.legend(fontsize=5.5,frameon=False,ncol=2,loc='upper left',bbox_to_anchor=(-.03,1.035),handlelength=.8,columnspacing=.6)
  for ax,order,threshold in [(c,y_order,1),(d,par_order,5)]:
   qq=q.reindex(order); qq=qq[qq.CTRL>=threshold]
   xx=np.arange(len(qq))
   ax.bar(xx,qq.percent,color=[COL[r] for r in qq.region],width=.72)
   ax.axhline(100,c='#999999',ls='--',lw=.6)
   ax.set(ylabel='% of sgCTRL',xlim=(-.6,len(qq)-.4))
   ax.set_xticks(xx,qq.index,rotation=90,fontsize=6)
   if ax is c:
    if line=='A549':
     ax.set(yscale='log',ylim=(.03,200));ax.set_yticks([.1,1,10,100],['0.1','1','10','100'])
     for i,value in enumerate(qq.percent):
      label='0%' if value==0 else (f'{value:.2f}%' if value<1 else f'{value:.1f}%')
      y=max(value*(1.85 if i%2 else 1.25),.047)
      if value==0: ax.scatter(i,.035,marker='v',s=15,facecolor='white',edgecolor=COL[qq.region.iloc[i]],lw=.6,zorder=4)
      ax.text(i,y,label,ha='center',fontsize=5.1,rotation=0)
    else:
     ax.set(yscale='log',ylim=(.65,200));ax.set_yticks([1,10,100],['1','10','100'])
     for i,value in enumerate(qq.percent): ax.text(i,value*1.16,f'{value:.0f}%',ha='center',fontsize=5.7)
   else:
    # Pool all annotated PAR1 + PAR2 features before taking the condition ratio.
    pooled=sums[sums.region.isin(['PAR1','PAR2'])].groupby(['sample','treatment']).CPM.sum().reset_index().groupby('treatment').CPM.mean()
    combined=100*pooled.LOY/pooled.CTRL
    combined_x=len(qq)+.45
    ax.bar(combined_x,combined,width=.8,color='#555555')
    ax.text(combined_x,combined+4,f'{combined:.0f}%',ha='center',fontsize=6)
    ax.axvline(len(qq)-.25,color='#cccccc',lw=.6)
    ax.set_xlim(-.6,combined_x+.65)
    ax.set_xticks(list(xx)+[combined_x],list(qq.index)+['Combined PAR'],rotation=90,fontsize=6)
    ax.set(ylim=(0,140));ax.set_yticks([0,50,100])
    pd.DataFrame([dict(region='PAR1 + PAR2',CTRL_CPM=pooled.CTRL,LOY_CPM=pooled.LOY,percent=combined)]).to_csv(S/line/'panel_d_combined_PAR.csv',index=False)
  regional.to_csv(S/line/'panel_b_regional_percentages.csv')
 if line in ['DU145','A549']:
  e.clear()
  e.text(-.085,1.07,'e',transform=e.transAxes,weight='semibold',fontsize=12)
  ep=pd.DataFrame(pp)
  for i,row in enumerate(ep.itertuples()):
   e.plot([i,i],[row.X_effect,row.Y_effect],c='#cccccc',lw=.7)
   e.scatter(i,row.X_effect,c=BLUE,s=20,zorder=3)
   e.scatter(i,row.Y_effect,c=RED,s=20,zorder=3)
 e.axhline(0,c='#aaaaaa',lw=.7);e.set(ylim=(-11,5),ylabel='Mean log$_2$CPM difference',xlim=(-.6,len(pairs)-.4));e.set_xticks(range(len(pairs)),[X+' / '+Y for X,Y in pairs],rotation=30,ha='right',fontsize=6.4)
 e.yaxis.label.set_size(7)
 e.legend(handles=[Line2D([],[],color=BLUE,marker='o',ls='',label='X-linked'),Line2D([],[],color=RED,marker='o',ls='',label='Y-linked')],fontsize=6,loc='upper right')
 if line=='A549':
  e.set_xlim(-.6,len(ep)-.4)
  e.set_xticks(range(len(ep)),[r.X+' / '+r.Y for r in ep.itertuples()],rotation=30,ha='right',fontsize=6.4)
 if line in ['DU145','A549']:
  e.set_xlim(-.6,len(ep)-.4)
  e.set_xticks(range(len(ep)),[r.X+' / '+r.Y for r in ep.itertuples()],rotation=45,ha='right',fontsize=5.7)
  e.legend(handles=[Line2D([],[],color=BLUE,marker='o',ls='',label='X-linked'),Line2D([],[],color=RED,marker='o',ls='',label='Y-linked')],fontsize=5.5,loc='lower right',bbox_to_anchor=(1,1.01),ncol=2,handletextpad=.2,columnspacing=.6)
 if line in ['DU145','A549']:
  e.set_ylim(-5 if line=='DU145' else -10,.9)
  fax=fig.add_axes([.60,.09,.36,.17])
  fax.text(-.10,1.07,'f',transform=fax.transAxes,weight='semibold',fontsize=12)
  totals=[]
  for i,row in enumerate(ep.itertuples()):
   ctrl_sum=q.loc[row.X,'CTRL']+q.loc[row.Y,'CTRL']
   loy_sum=q.loc[row.X,'LOY']+q.loc[row.Y,'LOY']
   totals.append(dict(X=row.X,Y=row.Y,CTRL_CPM=ctrl_sum,LOY_CPM=loy_sum,percent=100*loy_sum/ctrl_sum))
  totals=pd.DataFrame(totals); xx=np.arange(len(totals))
  fax.bar(xx-.17,[100]*len(totals),width=.32,color=BLUE)
  fax.bar(xx+.17,totals.percent,width=.32,color=RED)
  for i,value in enumerate(totals.percent): fax.text(i+.17,value+3,f'{value:.0f}',ha='center',fontsize=5.5)
  fax.axhline(100,c='#999999',lw=.6,ls='--')
  fax.set(ylim=(0,125),ylabel='Combined X + Y (% sgCTRL)',xlim=(-.6,len(totals)-.4))
  fax.yaxis.label.set_size(7)
  fax.set_xticks(xx,[r.X+' / '+r.Y for r in totals.itertuples()],rotation=45,ha='right',fontsize=5.7)
  fax.set_yticks([0,50,100])
  totals.to_csv(S/line/'panel_f_combined_XY.csv',index=False)
 stem=f'ch03_main_{line.lower()}_rna_sex_chromosomes_matched'
 for ext in ['pdf','svg','png']: fig.savefig(OUT/f'{stem}.{ext}',dpi=240)
 plt.close(fig)
 dst=S/line;dst.mkdir(exist_ok=True);data.to_csv(dst/'gene_period_values.csv',index=False);sums.to_csv(dst/'region_sample_sums.csv',index=False);pd.DataFrame(dna).to_csv(dst/'lpwgs_ratios.csv',index=False);pd.DataFrame(pp).to_csv(dst/'paralog_period_effects.csv',index=False)
 inputs=[S/f'{line}_{n}.csv' for n in ['CPM','effects','samples']]+[posfile,binsfile]
 (dst/'provenance.json').write_text(json.dumps({'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'region_rule':'PAR1/PAR2 from chrX annotation; distal MSY midpoint <4.8 Mb on CHM13v2; remaining MSY otherwise. Common comparison windows, not inferred deletion boundaries.','thresholds':'Y display CTRL mean CPM >=1; PAR >=5; pairs Y >=5. Missing marks are below display threshold.','RNA':'Ratios of condition mean CPM within each sampling period. No sample pairing inferred. Gene sums include all annotated regional features.','WGS':'Saved segments; omit ctrl CN <0.05; floor 0.001. Positional panel a displays percentages <=1 at 1. A549 PAR1 uses the ratio of median corrected bin signals, assigned by midpoint, rather than the saved segmentation.', 'display':'A549 and DU145 use six matched panels. Regional and PAR bars use linear percentage axes. A549 Y-gene bars retain positive percentages below 1; zero expression is marked with an open downward triangle at the lower display boundary. A549 panel e uses a wider range to retain its larger Y decreases.','periods':periods},indent=2),encoding='utf-8')
 print(stem)
