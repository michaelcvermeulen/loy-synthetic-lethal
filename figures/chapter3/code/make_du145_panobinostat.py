"""Main and supplementary DU145 treatment figures from verified saved outputs."""
from pathlib import Path
import sys, json, hashlib, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / 'figures/chapter3'
SRC = FIG / 'source_data/du145_drug_optimisation'
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
SAMPLES = ['LNCaP-sgCTRL','LNCaP-DMSO-LOY-1X','LNCaP-30','LNCaP-70','LNCaP-LOY-2X','LNCaP-KU','LNCaP-KU-30','LNCaP-KU-70']
LABEL = dict(zip(SAMPLES,['sgCTRL','1× DMSO','1× pano 30 nM','1× pano 75 nM','2×, no drug','1× KU','1× KU + pano 30 nM','1× KU + pano 75 nM']))
COLOR = dict(zip(SAMPLES,['#3C5488','#777777','#E99B57','#E64B35','#705A8C','#4DBBD5','#00A087','#8491B4']))
MARKER = dict(zip(SAMPLES,['o','o','s','D','^','v','P','X']))
COND = dict(zip(SAMPLES,['sgCTRL','DMSO vehicle','Panobinostat 30 nM','Panobinostat 75 nM','LOY 2x','KU-55933','KU-55933 + panobinostat 30 nM','KU-55933 + panobinostat 75 nM']))
w = pd.read_csv(SRC / 'lpwgs/chrY_qdnaseq_bins.csv')
w = w[w.use.eq(True) & w.log2_seg.notna()].copy()
w['relative_cn'] = 2 ** w.log2
w['segment_cn'] = 2 ** w.log2_seg
dd = pd.read_csv(SRC / 'drug_vs_vehicle.csv')
all_dd = pd.read_csv(SRC / 'all_target_measurements.csv')
controls = all_dd[all_dd.condition.eq('sgCTRL')].set_index('target').concentration_ratio
dd['percent_reduction_vs_sgCTRL'] = 100 * (1 - dd.concentration_ratio / dd.target.map(controls))
baseline = pd.read_csv(FIG / 'source_data/du145_ddpcr_assay_comparison/main_group_summary.csv')
two_round = {}
for target in ['AMELY', 'KDM5D']:
    rows = baseline[baseline.assay.str.startswith(target + '/')]
    edited = rows[rows.condition.eq('CP + Q53 + P47 2x')].iloc[0]
    control = rows[rows.condition.eq('sgCTRL') & rows.date.eq(edited.date)].iloc[0]
    two_round[target] = dict(ratio=edited['median'], control=control['median'],
                             reduction=100*(1-edited['median']/control['median']), date=edited.date)
dd.to_csv(SRC / 'drug_vs_vehicle_and_sgCTRL.csv', index=False)
pd.DataFrame(two_round).T.to_csv(SRC / 'two_round_ddpcr_comparator.csv', index_label='target')
cyto = pd.read_csv(FIG / 'source_data/ht1080_y_expression_confirmation/T2T_chrY_cytobands.csv')
cyto = cyto[cyto.chrom.eq('chrY')]
genes = pd.read_csv(FIG / 'source_data/ht1080_lpwgs_confirmation/ddpcr_target_loci.csv')
LENGTH = cyto.end.max()/1e6
C0 = cyto.loc[cyto.stain.eq('acen'),'start'].min()/1e6
C1 = cyto.loc[cyto.stain.eq('acen'),'end'].max()/1e6
Q0 = cyto.loc[cyto.band.eq('q12'),'start'].iloc[0]/1e6
SEG = []
for s,g in w.groupby('sample',sort=False):
    g=g.sort_values('start').copy()
    g['run']=(g.log2_seg.ne(g.log2_seg.shift()) | g.start.gt(g.end.shift()+1)).cumsum()
    for _,h in g.groupby('run'):
        SEG.append(dict(sample=s,start=(h.start.min()-1)/1e6,end=h.end.max()/1e6,value=h.segment_cn.iloc[0]))
seg=pd.DataFrame(SEG)
seg.to_csv(SRC/'lpwgs/plotted_saved_segments.csv',index=False)
w.to_csv(SRC/'lpwgs/plotted_chrY_bins.csv',index=False)

def tag(fig,letter,y):
    fig.text(.025,y,letter,fontsize=12,fontweight='semibold',va='top')

def context(ax):
    for a,b in [(0,2.458320),(62.122809,LENGTH)]: ax.axvspan(a,b,color='#EAF6F7',zorder=0)
    for a,b in [(C0,C1),(Q0,62.122809)]: ax.axvspan(a,b,color='#ECECEC',zorder=0)

def segments(ax,s,clip=None,lw=1.2,ls='-'):
    for r in seg[seg['sample'].eq(s)].itertuples():
        lo,hi=r.start,r.end
        if clip:lo,hi=max(lo,clip[0]),min(hi,clip[1])
        if hi>lo:ax.plot([lo,hi],[r.value]*2,color=COLOR[s],lw=lw,ls=ls)

def profiles(fig,rect,samples):
    left,bottom,width,height=rect
    gs=fig.add_gridspec(len(samples),1,left=left,right=left+width,bottom=bottom,top=bottom+height,hspace=.18)
    for i,s in enumerate(samples):
        ax=fig.add_subplot(gs[i]);context(ax)
        g=w[w['sample'].eq(s)]
        ax.scatter(g.midpoint,np.minimum(g.relative_cn,1.27),s=2.5,color='#999999',alpha=.65,lw=0)
        high=g[g.relative_cn>1.27];ax.scatter(high.midpoint,[1.27]*len(high),marker='^',s=7,color='#555555',lw=0)
        segments(ax,s)
        ax.set(xlim=(0,LENGTH),ylim=(-.025,1.35),yticks=[0,.5,1])
        ax.tick_params(axis='x',labelbottom=False,bottom=False)
        ax.tick_params(axis='y',labelsize=6.7,pad=2)
        ax.text(.99,.81,LABEL[s],transform=ax.transAxes,ha='right',color=COLOR[s],fontsize=8,fontweight='semibold')
        if i==0:ax.text(.66,.82,'Yq12: limited mappability',transform=ax.transAxes,ha='center',fontsize=7,color='#777777')
    fig.text(.036,bottom+height/2,'Relative copy number',rotation=90,va='center',fontsize=8)

def ideogram(fig,rect):
    ax=fig.add_axes(rect)
    for b in cyto.itertuples():
        lo,hi=b.start/1e6,b.end/1e6
        if b.stain=='acen':
            pts=[(lo,.2),(lo,.65),(hi,.425)] if b.band.startswith('p') else [(lo,.425),(hi,.2),(hi,.65)]
            ax.add_patch(Polygon(pts,facecolor='#888888',edgecolor='#555555',lw=.4))
        else:ax.add_patch(Rectangle((lo,.2),hi-lo,.45,facecolor={'gneg':'white','gpos50':'#999999','gvar':'#AFC0CD'}.get(b.stain,'#BBBBBB'),edgecolor='#555555',lw=.4))
    for lo,hi in [(0,2.458320),(62.122809,LENGTH)]:ax.add_patch(Rectangle((lo,.2),hi-lo,.45,color='#41B6C4',lw=0))
    for r in genes.itertuples():
        pos=(r.t2t_start+r.t2t_end)/2e6
        ax.plot([pos,pos],[.68,1.0],color='#555555',lw=.5)
        ax.text(pos,1.05,r.gene_name,ha='left',rotation=35,fontsize=6.5,fontstyle='italic')
    ax.set(xlim=(0,LENGTH),ylim=(0,2.2),yticks=[],xticks=np.arange(0,61,10),xlabel='HG002 T2T-Y position (Mb)')
    ax.spines[['left','right','top']].set_visible(False);ax.tick_params(axis='x',labelsize=7)

def ddpanel(ax,target,samples,full=False):
    ref='AMELX' if target=='AMELY' else 'RNaseP'
    ys=[]
    for i,s in enumerate(samples):
        row=dd[dd.condition.eq(COND[s]) & dd.target.eq(target)]
        assert len(row)==1
        row=row.iloc[0];y=row.concentration_ratio;ys.append(y)
        ax.vlines(i,0,y,color=COLOR[s],lw=1.2,alpha=.8)
        ax.scatter(i,y,s=25,marker=MARKER[s],facecolor='white',edgecolor=COLOR[s],lw=1.2,zorder=3)
        if not full:
            ax.text(i,y+.004,f'{row.percent_reduction_vs_sgCTRL:.1f}%',ha='center',fontsize=7,color=COLOR[s])
    ax.set(title=f'{target} / {ref}',ylim=(0,.09),yticks=[0,.03,.06,.09],xlim=(-.5,len(samples)-.5))
    ax.set_xticks(range(len(samples)))
    if full:ax.set_xticklabels(['1× DMSO','Pano 30','Pano 75','KU','KU + 30','KU + 75'],rotation=55,ha='right',fontsize=7)
    else:ax.set_xticklabels(['1× DMSO','Pano\n30 nM','Pano\n75 nM'],fontsize=7)
    ax.axhline(ys[0],ls=':',lw=.7,color='#999999',zorder=0)
    ax.grid(axis='y',color='#EEEEEE',lw=.4);ax.set_axisbelow(True)

def save(fig,stem):
    destinations=[(FIG/f'{stem}.{ext}',600) for ext in ['pdf','svg','png']]+[(FIG/'previews'/f'{stem}.png',160)]
    for dest,dpi in destinations:
        tmp=dest.with_name(dest.stem+'.staging'+dest.suffix)
        fig.savefig(tmp,dpi=dpi,facecolor='white')
        for attempt in range(10):
            try:tmp.replace(dest);break
            except OSError:
                if attempt==9:raise
                time.sleep(.5)
    plt.close(fig)

main='ch03_main_du145_panobinostat_single_round'
supp='ch03_supp_du145_drug_treatment_comparison'
fig=plt.figure(figsize=(183/25.4,252/25.4))
tag(fig,'a',.979)
ax=fig.add_axes([.10,.872,.85,.103]);ax.axis('off')
ax.text(0,.92,'Testing a single-round route to Y depletion',fontweight='semibold',fontsize=9)
for y,left,right in [(.55,'DMSO or pano (24 h before)','1× CP + Q53 + P47'),(.12,'2× CP + Q53 + P47','Two successive nucleofections; no drug')]:
    ax.text(.01,y,left,va='center',fontsize=8)
    ax.annotate('',xy=(.48,y),xytext=(.43,y),arrowprops=dict(arrowstyle='->',lw=.8,color='#555555'))
    ax.text(.50,y,right,va='center',fontsize=8)
ax.text(1,.92,'DU145',ha='right',fontweight='semibold',fontsize=9)
tag(fig,'c',.659)
profiles(fig,[.10,.29,.85,.359],SAMPLES[:5])
ideogram(fig,[.10,.242,.85,.042])
tag(fig,'d',.186)
# Expanded saved segments and separate retained-Yp summaries.
styles=['-', '--', '-', ':']
for j,(lo,hi,ymax,title) in enumerate([(6,10,.20,'Proximal Yp'),(13,Q0,.04,'Euchromatic Yq')]):
    ax=fig.add_axes([.10+j*.30,.084,.25,.091])
    for s,ls in zip(SAMPLES[1:5],styles):segments(ax,s,clip=(lo,hi),lw=1.25,ls=ls)
    ax.set(xlim=(lo,hi),ylim=(0,ymax),xlabel='T2T-Y position (Mb)')
    if j==0:ax.set_ylabel('Relative CN')
    ax.set_title(title,loc='left',fontsize=8);ax.tick_params(labelsize=7)
ax2=fig.add_axes([.70,.084,.25,.091]); reg=[]
for i,s in enumerate(SAMPLES[1:5]):
    x=w[w['sample'].eq(s)&w.midpoint.ge(2.46)&w.midpoint.lt(4.8)]
    v=x.relative_cn.median();reg.append(dict(sample=s,n_bins=len(x),median_relative_cn=v))
    ax2.scatter(i,v,s=25,marker=MARKER[s],color=COLOR[s]);ax2.text(i,v+.04,f'{v:.2f}',ha='center',fontsize=6.5)
ax2.set(title='Retained distal Yp',ylim=(0,.6),yticks=[0,.3,.6],xticks=range(4),xticklabels=['Veh.','30','75','2×'])
ax2.tick_params(labelsize=7);ax2.set_xlabel('1×: vehicle / pano (nM)',fontsize=6.5)
handles=[Line2D([],[],color=COLOR[s],lw=1.4,ls=ls,marker=MARKER[s],markersize=3,label=LABEL[s]) for s,ls in zip(SAMPLES[1:5],styles)]
fig.legend(handles=handles,loc='center',bbox_to_anchor=(.53,.035),ncol=4,frameon=False,fontsize=7,handlelength=1.8,columnspacing=1.2)
tag(fig,'b',.854)
fig.text(.10,.862,'ddPCR | labels: depletion (%) versus sgCTRL',fontsize=9,fontweight='semibold',va='top')
for i,target in enumerate(['AMELY','DDX3Y','KDM5D']):
    ax=fig.add_axes([.10+i*.30,.717,.25,.108]);ddpanel(ax,target,SAMPLES[1:4])
    if i==0:ax.set_ylabel('Target/reference ratio')
fig.text(.10,.012,'Source: DU145 treatment lpWGS and 25 March 2026 ddPCR.',fontsize=6.7,color='#666666')
save(fig,main)
pd.DataFrame(reg).to_csv(SRC/'lpwgs/plotted_distal_Yp_summary.csv',index=False)

fig=plt.figure(figsize=(183/25.4,265/25.4))
tag(fig,'c',.537)
profiles(fig,[.10,.125,.85,.395],[SAMPLES[i] for i in [0,1,2,3,5,6,7,4]])
ideogram(fig,[.10,.081,.85,.039])
tag(fig,'a',.970)
fig.text(.10,.970,'DU145 ddPCR | all six single-round conditions',fontsize=9,fontweight='semibold',va='top')
ss=[SAMPLES[i] for i in [1,2,3,5,6,7]]
for i,t in enumerate(['AMELY','DDX3Y','KDM5D']):
    ax=fig.add_axes([.10+i*.30,.817,.25,.10]);ddpanel(ax,t,ss,full=True)
    if i==0:ax.set_ylabel('Target/reference ratio')
tag(fig,'b',.737)
for bottom,metric,title in [(.641,'percent_reduction_vs_sgCTRL','Overall target-dosage reduction (%) versus sgCTRL'),
                            (.545,'percent_reduction_vs_vehicle','Additional reduction of remaining signal (%) versus 1X DMSO')]:
    overall = metric == 'percent_reduction_vs_sgCTRL'
    cmap = plt.get_cmap('magma_r' if overall else 'RdBu_r')
    norm = matplotlib.colors.Normalize(vmin=0 if overall else -100, vmax=100)
    ax=fig.add_axes([.10,bottom,.77,.061]);ax.axis('off')
    fig.text(.10,bottom+.074,title,fontsize=7.5,fontweight='semibold')
    table=[]
    for t in ['AMELY','DDX3Y','KDM5D']:
        values=[f'{dd.loc[dd.condition.eq(COND[s]) & dd.target.eq(t),metric].iloc[0]:.1f}' for s in ss]
        if overall:
            values.append(f'{two_round[t]["reduction"]:.1f}' if t in two_round else 'Not measured')
        table.append([t]+values)
    headers=['Target','1X DMSO','Pano 30','Pano 75','KU','KU + 30','KU + 75']
    if overall: headers.append('2X*')
    tab=ax.table(cellText=table,colLabels=headers,cellLoc='center',bbox=[0,0,1,1])
    tab.auto_set_font_size(False);tab.set_fontsize(6.6)
    for (r,c),cell in tab.get_celld().items():
        cell.set_edgecolor('#DDDDDD');cell.set_linewidth(.4)
        if r==0:cell.set_facecolor('#F1F1F1')
        elif c>0:
            try:
                value=float(cell.get_text().get_text())
            except ValueError:
                cell.set_facecolor('#E7E7E7');cell.get_text().set_fontsize(5.8)
                continue
            rgba=cmap(norm(value));cell.set_facecolor(rgba)
            luminance=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
            cell.get_text().set_color('white' if luminance<.45 else '#202020')
    cax=fig.add_axes([.90,bottom,.012,.061])
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,
                    ticks=[0,50,100] if overall else [-100,0,100])
    cb.ax.tick_params(labelsize=6,pad=2,length=2);cb.outline.set_linewidth(.4)
    cb.ax.set_title('%',fontsize=6,pad=3)
fig.text(.10,.744,'All drug conditions: 1X pooled CP + Q53 + P47. *2X: no drug, separate February assays.',fontsize=6.3)
fig.text(.10,.013,'March: same-run sgCTRL / DMSO. 2X: assay-matched sgCTRL; no cross-run comparison with March DMSO.',fontsize=6.3,color='#666666')
save(fig,supp)
manifest={'figures':[main,supp],'font':'Open Sans','main_mm':[183,252],'supp_mm':[183,265],
          'png_dpi':600,'n_WGS_per_condition':1,'n_ddPCR_per_target_condition':1,'inferential_tests':False,
          'color_mapping':COLOR,'segments':'Existing saved QDNAseq segments, no refitting',
          'relative_CN':'2**saved log2; original numerical floor 0.001 retained',
          'provenance':'Input recovery and hashes are in lpwgs/provenance.json and ../provenance.json',
          'outputs':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for stem in [main,supp] for ext in ['pdf','svg','png'] for f in [FIG/f'{stem}.{ext}']}}
(SRC/'figure_render_provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Rendered',main,'and',supp)
