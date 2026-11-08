from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
root=Path(__file__).resolve().parents[3]
out=root/'figures/chapter3'
data=out/'source_data/prism_rep1m_boxplots'
d=pd.read_csv(data/'cell_line_responses.csv')
s=pd.read_csv(data/'compound_statistics.csv').sort_values('p')
font_manager.fontManager.addfont(str(root/'skills/thesis-figure-style/assets/fonts/OpenSans-Regular.ttf'))
plt.rcParams.update({'font.family':'Open Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,4,figsize=(9.1,5.7))
rng=np.random.default_rng(121)
for ax,(_,r) in zip(axs.flat,s.iterrows()):
    z=d[d.feature==r.feature]
    arrays=[]
    counts=[]
    for x,(group,col) in enumerate([('Y_positive','#3C5488'),('Y_negative','#E64B35')]):
        y=z.loc[z.y_call==group,'response'].to_numpy()
        arrays.append(y);counts.append(len(y))
        ax.scatter(x+rng.uniform(-.22,.22,len(y)),y,s=7,color=col,alpha=.32,edgecolors='none',zorder=2)
    bp=ax.boxplot(arrays,positions=[0,1],widths=.46,showfliers=False,patch_artist=True,medianprops={'color':'#202020','linewidth':1},whiskerprops={'linewidth':.8},capprops={'linewidth':.8},boxprops={'linewidth':.8},zorder=3)
    for patch,col in zip(bp['boxes'],['#3C5488','#E64B35']):
        patch.set_facecolor(col);patch.set_alpha(.22)
    name='IT1t' if r.compound_name.strip()=='IT1T' else r.compound_name.strip().capitalize()
    ax.text(.5,1.26,name,ha='center',transform=ax.transAxes,fontsize=9)
    annotation=r.label_annotation.replace('Bacterial cell-wall inhibitor','Bacterial cell-wall\ninhibitor').replace('Tubulin polymerization inhibitor','Tubulin polymerization\ninhibitor').replace('Bacterial 30S ribosome inhibitor','Bacterial 30S ribosome\ninhibitor').replace('Antimalarial; target unspecified','Antimalarial; target\nunspecified')
    ax.text(.5,1.14,annotation,ha='center',va='center',transform=ax.transAxes,fontsize=6.5)
    ax.text(.5,1.015,f'P = {r.p:.3g}; FDR = {r.fdr:.3f}',ha='center',va='bottom',transform=ax.transAxes,fontsize=7)
    ax.set_xticks([0,1],[f'Y-positive\n(n = {counts[0]})',f'LOY\n(n = {counts[1]})'])
    ax.set_xlim(-.5,1.5)
    ax.axhline(0,lw=.5,color='#CCCCCC',zorder=0)
    ax.tick_params(length=3)
    ax.margins(y=.07)
for ax in axs[:,0]: ax.set_ylabel('PRISM response (log₂ fold change)')
fig.subplots_adjust(left=.075,right=.99,bottom=.12,top=.875,hspace=.78,wspace=.32)
fig.text(.5,.018,'Lower response indicates greater drug sensitivity. Points: cell lines; boxes: median and IQR; whiskers: 1.5 × IQR.',ha='center',fontsize=7)
stem='ch03_supp_prism_rep1m_boxplots'
for ext in ['png','pdf','svg']:fig.savefig(out/f'{stem}.{ext}',dpi=300,facecolor='white')
plt.close(fig)
(data/'README.md').write_text('Exploratory raw PRISM REP.1M responses for the eight compounds labelled on the volcano (four smallest nominal P values per direction). Cell-line group counts and means verified against drug_pan_results.csv. Displayed P and FDR are the original lineage-adjusted, donor-clustered results, not an unadjusted box-plot comparison. All FDR >0.05. Independent y scales retain all observations. Extraction and rendering scripts: ../../code/ch03_supp_prism_rep1m_boxplots_extract.R and ../../code/ch03_supp_prism_rep1m_boxplots.py. Central caption: ../../../../Chapter3_figure_captions.md. Not inserted into Results.\n')
print(out/f'{stem}.png')
