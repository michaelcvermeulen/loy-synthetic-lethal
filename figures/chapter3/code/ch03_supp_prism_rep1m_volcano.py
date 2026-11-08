from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[3]
source=ROOT/"research/depmap_loy_reanalysis_v1/loy_reanalysis_v1_export_2026-10-05/derived/loy_reanalysis_v1/drug_pan_results.csv"
out=ROOT/"figures/chapter3"
stem="ch03_supp_prism_rep1m_volcano"
data=out/"source_data/prism_rep1m_volcano"
data.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(source).query("dataset == 'prism_REP1M'").copy()
assert len(d)==1265 and d[['standardized','p','fdr']].notna().all().all()
assert ((d.p>0)&(d.p<=1)).all()
d['neg_log10_p']=-np.log10(d.p)
d['labelled']=False
for sign in [-1,1]:
    idx=d[d.standardized*sign>0].nsmallest(4,'p').index
    d.loc[idx,'labelled']=True
mechanisms = {
    'ENTASOBULIN': 'Tubulin polymerization inhibitor',
    'VIXOTRIGINE': 'Sodium channel blocker',
    'LORACARBEF': 'Bacterial cell-wall inhibitor',
    'LARAZOTIDE-ACETATE': 'Tight-junction regulator',
    'CIPARGAMIN': 'Antimalarial; target unspecified',
    'OMADACYCLINE': 'Bacterial 30S ribosome inhibitor',
    'EVOCALCET': 'Calcium sensitizer',
}
d['label_annotation'] = d.apply(lambda r: r.targets if pd.notna(r.targets) and str(r.targets).strip() else mechanisms.get(r.compound_name.strip(), ''), axis=1)
d['label_annotation_type'] = np.where(d.targets.notna(), 'recorded target', 'recorded mechanism / target unspecified')
d.to_csv(data/"plotted_results.csv",index=False)
font_manager.fontManager.addfont(str(ROOT/"skills/thesis-figure-style/assets/fonts/OpenSans-Regular.ttf"))
plt.rcParams.update({'font.family':'Open Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(6.3,4.3))
ax.scatter(d.standardized,d.neg_log10_p,s=13,c='#9B9B9B',alpha=.65,edgecolors='none')
sel=d[d.labelled]
ax.scatter(sel.standardized,sel.neg_log10_p,s=20,c='#343434',edgecolors='white',linewidths=.35,zorder=3)
ax.axvline(0,color='#BBBBBB',lw=.7,zorder=0)
ax.axhline(-np.log10(.05),color='#777777',lw=.7,ls=(0,(4,3)),zorder=0)
ax.text(.505,-np.log10(.05)+.045,'Nominal P = 0.05',ha='right',va='bottom',fontsize=8,color='#555555')
for sign in [-1,1]:
    rows=sel[sel.standardized*sign>0].sort_values('p')
    for (_,r),y in zip(rows.iterrows(),[2.92,2.54,2.16,1.78]):
        ax.annotate(('IT1t' if r.compound_name.strip()=='IT1T' else r.compound_name.strip().capitalize())+'\n'+r.label_annotation,xy=(r.standardized,r.neg_log10_p),xytext=(sign*.64,y),ha='left' if sign<0 else 'right',va='center',fontsize=6.5,linespacing=1.2,arrowprops={'arrowstyle':'-','lw':.65,'color':'#666666','shrinkA':2,'shrinkB':3},zorder=4)
ax.set(xlim=(-.66,.66),ylim=(0,3.38),xlabel='LOY-associated drug-response difference (SD units)',ylabel=r'$-\log_{10}$ P value')
ax.set_xticks([-.4,-.2,0,.2,.4])
ax.text(.02,.975,f'{len(d):,} compounds; 0 at FDR <0.05\nMinimum FDR = {d.fdr.min():.3f}',transform=ax.transAxes,va='top',fontsize=8.5)
fig.subplots_adjust(left=.12,right=.985,bottom=.23,top=.98)
ax.text(0,-.205,'Greater sensitivity in LOY',transform=ax.transAxes,ha='left',fontsize=8)
ax.text(1,-.205,'Lower sensitivity in LOY',transform=ax.transAxes,ha='right',fontsize=8)
for ext in ['pdf','svg','png']:
    fig.savefig(out/f'{stem}.{ext}',dpi=300,facecolor='white')
plt.close(fig)
(data/'provenance.json').write_text(json.dumps({'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'dataset':'prism_REP1M','n':len(d),'fdr_hits':int((d.fdr<.05).sum()),'x':'Saved standardized regression effect; negative denotes greater sensitivity in LOY','y':'-log10 nominal P','labels':'Four smallest nominal P values per direction; recorded target if available, otherwise saved mechanism. Font 6.5 pt. Illustrative, not FDR significant','correction':'Saved BH correction within all 1265 REP.1M tests, unchanged','font':'Bundled Open Sans Regular','code':str(Path(__file__).resolve())},indent=2))
(data/'README.md').write_text('Exploratory PRISM REP.1M volcano. No models were refitted. Full plotted statistics are in plotted_results.csv; source and conventions are in provenance.json. Active caption: [Chapter3_figure_captions.md](../../../../Chapter3_figure_captions.md). Not assigned a coded figure number or inserted into Results.\n')
print(out/f'{stem}.png')
