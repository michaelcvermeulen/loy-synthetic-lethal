"""Extend S02 with saved sgCTRL-versus-parental GSEA; no analysis refit."""
from pathlib import Path
import hashlib,json,textwrap
import pandas as pd
import numpy as np
FIG=Path(__file__).resolve().parents[1]
THESIS=FIG.parents[1]
ROOT=THESIS.parent/'RNAseq/analysis/HT1080/SEPT_11_2026_REANALYZE'
source_path=ROOT/'plot_programme_analysis_v1.py'
source=source_path.read_text(encoding='utf-8')
source=source[:source.index("src=FIG/'source_data/ht1080_rna_programme_overview'")]
source=source.replace("ROOT=Path(__file__).resolve().parent;",f"ROOT=Path({str(ROOT)!r});")
source=source.replace('Leading pathway patterns','Y-depleted versus sgCTRL: pathways').replace('Inferred transcription-factor activity','Y-depleted versus sgCTRL: TF activity')
exec(compile(source,str(source_path),'exec'))
# Preserve panels a-d and their data; use the new right half for the distinct contrast.
fig.set_size_inches(14,9.25)
for ax in fig.axes:
    pos=ax.get_position();ax.set_position([pos.x0*.52,pos.y0,pos.width*.52,pos.height])
for text in fig.texts:
    x,y=text.get_position();text.set_position((x*.52,y))
fig.texts[-1].set_text('HT1080 expression comparisons')
path=ROOT/'outputs/programme_analysis_v1/all_collections/GSEA_CTRL_vs_parental.csv'
allsets=pd.read_csv(path)
groups=[('e','Hallmark',allsets.collection.eq('H')),
        ('f','GO Biological Process',allsets.subcollection.eq('GO:BP')),
        ('g','Reactome',allsets.subcollection.eq('CP:REACTOME'))]
chosen=[]
fig.text(.55,.985,'sgCTRL versus parental (August): GSEA',fontsize=10,weight='semibold',va='top')
for (letter,name,mask),bottom in zip(groups,[.72,.405,.09]):
    eligible=allsets[mask & allsets.FDR.lt(.05)].copy()
    top=pd.concat([eligible[np.sign(eligible.NES)==sign].sort_values(['FDR','pval','pathway']).head(3) for sign in [-1,1]])
    top['display_collection']=name;chosen.append(top)
    ax=fig.add_axes([.77,bottom,.12,.205])
    y=np.arange(len(top));ax.barh(y,top.NES,color=np.where(top.NES>0,red,blue),height=.55)
    names=[]
    for p in top.pathway:
        for prefix in ['HALLMARK_','GOBP_','REACTOME_']:p=p.removeprefix(prefix)
        names.append('\n'.join(textwrap.wrap(p.replace('_',' ').lower(),34)))
    ax.set_yticks(y,names,fontsize=7);ax.invert_yaxis();ax.tick_params(axis='y',length=0,pad=4)
    ax.axvline(0,c='#777777',lw=.6);ax.set_xlim(-3,3);ax.set_xticks([-3,0,3]);ax.set_xlabel('NES',fontsize=8)
    ax.spines[['top','right','left']].set_visible(False)
    ax.text(1.10,1.05,'NES / FDR',transform=ax.transAxes,fontsize=7,ha='left')
    for j,q in enumerate(top.itertuples()):
        ax.text(1.10,j,f'{q.NES:+.2f} / {q.FDR:.2g}',transform=ax.get_yaxis_transform(),fontsize=7,va='center',ha='left')
    fig.text(.535,bottom+.235,letter,fontsize=12,weight='semibold')
    fig.text(.565,bottom+.235,name,fontsize=8,weight='semibold')
fig.text(.56,.028,'Negative NES: lower in sgCTRL; positive NES: higher in sgCTRL\nTop three significant sets per direction; all-collection BH FDR < 0.05',fontsize=7,color='#555555')
src=FIG/'source_data/ht1080_rna_programme_overview'
top=pd.concat(chosen);top.to_csv(src/'parental_GSEA_display.csv',index=False)
for ext in ['pdf','svg','png']:fig.savefig(FIG/f'{STEM}.{ext}',dpi=300)
fig.savefig(FIG/'previews'/f'{STEM}_preview.png',dpi=160)
(src/'parental_GSEA_provenance.json').write_text(json.dumps({'source':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'contrast':'sgCTRL minus parental, August only','selection':'Up to three FDR <0.05 sets per direction per collection; ordered by FDR, pval, pathway','FDR_family':'all 23526 sets in the saved contrast','ranking':'saved moderated t','displayed':len(top),'analysis_refitted':False,'original_panels':'a-d unchanged; new e Hallmark, f GO BP, g Reactome'},indent=2))
print('Saved S02 with',len(top),'parental GSEA hits.')
