"""Alternative display of verified DNA profiles and paired DNA/RNA arm summaries."""
from pathlib import Path
import sys, io, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

FIG=Path(__file__).resolve().parents[1]
THESIS=FIG.parents[1]
SRC=FIG/'source_data/a549_du145_rna_dna'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
STEM='ch03_supp_a549_du145_dna_rna_arm_comparison'
RNA_SCATTER_YLIM=(-8.3,2.5)
RNA_SCATTER_YTICKS=[-8,-6,-4,-2,0,2]
DNA_SCATTER_XLIM=(-8,1)
DNA_SCATTER_XTICKS=[-8,-6,-4,-2,0]
chroms=['chr'+str(i) for i in range(1,23)]+['chrX','chrY']
cb=pd.read_csv(THESIS.parent/'WGS/chm13v2.0_cytobands_allchrs.bed',sep='\t',header=None)
lengths=cb.groupby(0)[2].max()
yq12=cb[(cb[0]=='chrY') & (cb[3]=='q12')]
assert len(yq12)==1
YQ12_START,YQ12_END=map(int,yq12.iloc[0,[1,2]])  # BED: zero-based, half-open.
excluded_counts={}
offsets={};end=0
for ch in chroms:
    offsets[ch]=end;end+=lengths[ch]/1e6+20
profiles={}
BIN_SIZE=5_000_000
for line,ctrl,loy in [('A549','A549-CTRL2','A549-LOY'),('DU145','LNCaP-sgCTRL','LNCaP-LOY-2X')]:
    # Aggregate unsegmented corrected values, not the saved segment estimates.
    source=FIG/f'source_data/{line.lower()}_lpwgs_confirmation/qdnaseq_bins.csv'
    raw=pd.read_csv(source)
    raw=raw[raw['use'] & raw['sample'].isin([ctrl,loy]) & np.isfinite(raw.log2)]
    paired=raw.pivot(index=['chromosome','start','end'],columns='sample',values='log2').dropna().reset_index()
    # Source QDNAseq coordinates are one-based inclusive. Exclude any overlap,
    # including the 100-kb bin crossing the cytoband boundary; never filter by value.
    excluded=paired.chromosome.eq('chrY') & (paired.start<=YQ12_END) & (paired.end>YQ12_START)
    paired.loc[excluded].to_csv(SRC/f'{line}_Yq12_excluded_source_bins.csv',index=False)
    excluded_counts[line]=int(excluded.sum())
    paired=paired.loc[~excluded].copy()
    paired['window_start']=(((paired.start+paired.end)/2)//BIN_SIZE).astype(int)*BIN_SIZE
    b=paired.groupby(['chromosome','window_start']).agg(control=(ctrl,'median'),LOY=(loy,'median'),n_source_bins=('start','size'),source_start=('start','min'),source_end=('end','max')).reset_index()
    b['window_end']=np.minimum(b.window_start+BIN_SIZE,b.chromosome.map(lengths))
    b['window_end']=np.minimum(b.window_end,b.source_end)
    b['mid']=(b.window_start+b.window_end)/2
    b['DNA']=b.LOY-b.control
    b.to_csv(SRC/f'{line}_DNA_unsegmented_5Mb.csv',index=False)
    profiles[line]=b
arms={line:pd.read_csv(SRC/f'{line}_arm.csv').set_index('unit') for line in profiles}
order=[ch+arm for ch in chroms[:-1] for arm in ['p','q'] if ch+arm in arms['A549'].index]
assert set(order)==set(arms['DU145'].index) and len(order)==41
track_extent=max(b.DNA.abs().max() for b in profiles.values())*1.1
vals=np.concatenate([a[['DNA','RNA']].to_numpy().ravel() for a in arms.values()])
bar_limits=(min(-.1,vals.min()-.08),max(.1,vals.max()+.1))
# Regional MSY summaries: exclude PAR, centromere and Yq12 explicitly.
cen=cb[(cb[0]=='chrY') & (cb[4]=='acen')]
regions={'MSY Yp':(2458320,int(cen[1].min())), 'MSY Yq':(int(cen[2].max()),YQ12_START)}
msy={}; correlations=[]
for line,ctrl,loy in [('A549','A549-CTRL2','A549-LOY'),('DU145','LNCaP-sgCTRL','LNCaP-LOY-2X')]:
    raw=pd.read_csv(FIG/f'source_data/{line.lower()}_lpwgs_confirmation/qdnaseq_bins.csv')
    raw=raw[raw.use & raw.chromosome.eq('chrY') & raw['sample'].isin([ctrl,loy]) & np.isfinite(raw.log2_seg)]
    paired=raw.pivot(index=['start','end'],columns='sample',values='log2_seg').dropna().reset_index()
    genes=pd.read_csv(SRC/f'{line}_8panel_Y_gene_mapping.csv')
    rows=[]
    for region,(start,end_region) in regions.items():
        b=paired[(paired.start>start) & (paired.end<=end_region)]
        g=genes[genes.included & genes.mid.ge(start) & genes.mid.lt(end_region)]
        assert len(b)>0 and len(g)>0
        rows.append(dict(unit=region,DNA=(b[loy]-b[ctrl]).median(),RNA=g.mean_difference.mean(),bins=len(b),genes=len(g),BED_start=start,BED_end=end_region))
    msy[line]=pd.DataFrame(rows).set_index('unit')
    msy[line].to_csv(SRC/f'{line}_S07_MSY_regions.csv')
    combined=pd.concat([arms[line],msy[line]])
    assert len(combined)==43
    correlations.append(dict(cell_line=line,n=43,Pearson=combined.DNA.corr(combined.RNA),nonY_n=41,nonY_Pearson=arms[line].DNA.corr(arms[line].RNA)))
pd.DataFrame(correlations).to_csv(SRC/'S07_MSY_correlations.csv',index=False)
windows={line:pd.read_csv(SRC/f'{line}_window_10Mb.csv') for line in profiles}
window_all=pd.concat(list(windows.values()))
fig=plt.figure(figsize=(15,8))
gs=fig.add_gridspec(4,3,width_ratios=[3.2,1.2,1.5],height_ratios=[.65,1.35,.65,1.35],
                   left=.055,right=.985,bottom=.075,top=.975,hspace=.53,wspace=.32)
labels=[]
rna_scatter_axes=[]
scatter_all=pd.concat([a[['DNA','RNA']] for a in arms.values()])
def limits(v):
    lo,hi=v.min(),v.max();return lo-.12*(hi-lo),hi+.20*(hi-lo)
for row,line in enumerate(profiles):
    ax=fig.add_subplot(gs[row*2,:]);b=profiles[line]
    for i,ch in enumerate(chroms):
        z=b[b.chromosome.eq(ch)]
        ax.scatter(z.mid/1e6+offsets[ch],z.DNA,s=6,c='#E64B35' if ch=='chrY' else ['#555555','#AAAAAA'][i%2],lw=0)
    ax.axhline(0,c='#333333',lw=.5)
    mask_left=offsets['chrY']+YQ12_START/1e6
    mask_right=offsets['chrY']+YQ12_END/1e6
    ax.axvspan(mask_left,mask_right,color='#EEEEEE',zorder=-1)
    ax.annotate('Yq12 excluded',xy=((mask_left+mask_right)/2,0),xytext=(.99,.78),
                textcoords='axes fraction',ha='right',fontsize=7,color='#666666',
                arrowprops=dict(arrowstyle='-',lw=.5,color='#888888'))
    ax.set_yscale('linear');ax.set_ylim(-track_extent,track_extent)
    ax.set_yticks([-10,-5,0,5,10])
    assert b.DNA.between(-track_extent,track_extent).all(), 'DNA-track point would be clipped'
    ax.set_xlim(-10,end)
    ax.set_xticks([offsets[ch]+lengths[ch]/2e6 for ch in chroms],[ch[3:] for ch in chroms],fontsize=7)
    ax.set_ylabel('DNA log$_2$ difference')
    ax.text(.01,.08,line,transform=ax.transAxes,weight='semibold',fontsize=10)
    ax.text(-.055,1.02,'a' if row==0 else 'd',transform=ax.transAxes,weight='semibold',fontsize=12)
    ax.text(.09,.90,'5-Mb bins',transform=ax.transAxes,ha='left',va='top',fontsize=8)
    ax=fig.add_subplot(gs[row*2+1,0]);a=arms[line].loc[order];x=np.arange(len(order))
    ax.bar(x-.2,a.DNA,width=.38,color='#555555',label='DNA: median segment difference')
    ax.bar(x+.2,a.RNA,width=.38,color='#00A087',label='RNA: mean log$_2$ CPM difference')
    ax.axhline(0,c='#333333',lw=.6)
    for i in range(1,len(order)):
        if order[i][:-1]!=order[i-1][:-1]:ax.axvline(i-.5,color='#EEEEEE',lw=.5,zorder=0)
    ax.set_xticks(x,[v[3:] for v in order],rotation=90,fontsize=7)
    ax.set_xlim(-.7,len(order)-.3);ax.set_ylim(*bar_limits)
    ax.set_ylabel('LOY minus sgCTRL\n(log$_2$ units)')
    ax.text(-.075,1.02,'b' if row==0 else 'e',transform=ax.transAxes,weight='semibold',fontsize=12)
    if row==0:ax.legend(loc='lower left',frameon=False,fontsize=7.5,ncol=1)
    a.reset_index().to_csv(SRC/f'{line}_arm_bar_display.csv',index=False)
    ax=fig.add_subplot(gs[row*2+1,1])
    ax.axhline(0,c='#DDDDDD',lw=.6);ax.axvline(0,c='#DDDDDD',lw=.6)
    ax.scatter(a.DNA,a.RNA,s=18,c='#888888',lw=0)
    # Four furthest arms from the within-line median, scaled by each axis SD.
    score=(((a[['DNA','RNA']]-a[['DNA','RNA']].median())/a[['DNA','RNA']].std())**2).sum(axis=1)
    selected=score.nlargest(4).index
    ax.scatter(a.loc[selected,'DNA'],a.loc[selected,'RNA'],s=23,c='#3C5488',lw=0)
    y=msy[line]
    ax.scatter(y.DNA,y.RNA,s=27,c='#E64B35',lw=0)
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.set_xlim(*DNA_SCATTER_XLIM);ax.set_ylim(*RNA_SCATTER_YLIM)
    ax.set_xticks(DNA_SCATTER_XTICKS);ax.set_yticks(RNA_SCATTER_YTICKS)
    texts=[]
    for arm in selected:
        r=a.loc[arm]
        texts.append(ax.text(r.DNA,r.RNA,arm[3:],fontsize=7.5))
        labels.append(dict(cell_line=line,arm=arm,DNA=r.DNA,RNA=r.RNA,label_score=score[arm]))
    for region,r in y.iterrows():texts.append(ax.text(r.DNA,r.RNA,region,fontsize=7,color='#E64B35'))
    np.random.seed(42)
    adjust_text(texts,x=np.r_[a.DNA,y.DNA],y=np.r_[a.RNA,y.RNA],ax=ax,expand=(1.25,1.35),ensure_inside_axes=True,
                arrowprops=dict(arrowstyle='-',color='#666666',lw=.5))
    combined=pd.concat([a,y])
    ax.text(.97,.035,f'All: r = {combined.DNA.corr(combined.RNA):.3f}; n = 43\nNon-Y: r = {a.DNA.corr(a.RNA):.3f}; n = 41',transform=ax.transAxes,ha='right',va='bottom',fontsize=7)
    ax.set_xlabel('DNA log$_2$ difference')
    ax.set_ylabel('RNA mean log$_2$ CPM difference')
    rna_scatter_axes.append(ax)
    assert combined.RNA.between(*RNA_SCATTER_YLIM).all(), 'RNA point would be clipped'
    assert combined.DNA.between(*DNA_SCATTER_XLIM).all(), 'DNA point would be clipped'
    ax.text(-.24,1.02,'c' if row==0 else 'f',transform=ax.transAxes,weight='semibold',fontsize=12)
    ax=fig.add_subplot(gs[row*2+1,2]);w=windows[line]
    ax.axhline(0,color='#DDDDDD',lw=.6);ax.axvline(0,color='#DDDDDD',lw=.6)
    ax.scatter(w.DNA,w.RNA,s=12,color='#777777',alpha=.75,lw=0)
    ax.set_xlim(*limits(window_all.DNA));ax.set_ylim(*limits(window_all.RNA))
    ax.text(.04,.96,f'10-Mb windows\nr = {w.DNA.corr(w.RNA):.3f}; n = {len(w)}',transform=ax.transAxes,va='top',fontsize=8)
    texts=[]
    if line=='DU145':
        label_rows=pd.read_csv(SRC/'DU145_10Mb_cytoband_labels.csv')
        for q in label_rows.itertuples():
            texts.append(ax.text(q.DNA,q.RNA,q.label,fontsize=6.7,color='#3C5488'))
            ax.scatter([q.DNA],[q.RNA],s=17,color='#3C5488',lw=0)
    else:
        q=w.loc[w.DNA.idxmin()]
        texts.append(ax.text(q.DNA,q.RNA,q.unit.replace('chr',''),fontsize=7))
    np.random.seed(42)
    adjust_text(texts,ax=ax,expand=(1.3,1.5),force_text=(.5,.7),arrowprops=dict(arrowstyle='-',color='#888888',lw=.5))
    ax.set_xlabel('DNA log$_2$ difference');ax.set_ylabel('RNA mean log$_2$ CPM difference')
    ax.text(-.24,1.02,'g' if row==0 else 'h',transform=ax.transAxes,weight='semibold',fontsize=12)
for ax in fig.axes:ax.spines[['top','right']].set_visible(False)
assert all(ax.get_xscale()==ax.get_yscale()=='linear' for ax in fig.axes)
assert all(ax.get_yscale()=='linear' and ax.get_ylim()==RNA_SCATTER_YLIM for ax in rna_scatter_axes)
assert all(ax.get_xlim()==DNA_SCATTER_XLIM for ax in rna_scatter_axes)
assert np.array_equal(rna_scatter_axes[0].get_xticks(),rna_scatter_axes[1].get_xticks())
assert np.array_equal(rna_scatter_axes[0].get_yticks(),rna_scatter_axes[1].get_yticks())
pd.DataFrame(labels).to_csv(SRC/'arm_scatter_labels.csv',index=False)
for ext in ['pdf','svg','png']:
    buf=io.BytesIO();fig.savefig(buf,format=ext,dpi=600);(FIG/f'{STEM}.{ext}').write_bytes(buf.getvalue())
fig.savefig(SRC/f'{STEM}_preview.png',dpi=144,facecolor='white')
plt.close(fig)
(SRC/'arm_bar_display_notes.json').write_text(json.dumps({'comparisons':list(profiles),'arms':order,
 'omitted_arms':['chr13p','chr14p','chr15p','chr21p','chr22p'],
 'reason':'fewer than 20 mapped expressed primary genes; DNA remains in genome tracks',
 'statistics':'unchanged from paired source tables; no uncertainty intervals across genomic units',
 'DNA_track_scale':'linear spacing of existing log2 differences; shared symmetric limits; no clipping of 5-Mb values',
 'DNA_track_aggregation':'median unsegmented corrected log2 values per population in fixed 5-Mb windows, then LOY minus control; paired finite usable source bins',
 'Yq12_mask':{'chromosome':'chrY','BED_start':YQ12_START,'BED_end':YQ12_END,
              'source':'WGS/chm13v2.0_cytobands_allchrs.bed; q12',
              'rule':'exclude every overlapping source bin before aggregation; boundary windows truncate to retained source span',
              'excluded_source_bins':excluded_counts},
 'scatter_labels':'four largest sums of squared median-centred, SD-scaled DNA/RNA deviations per cell line; descriptive selection, not outlier test',
 'MSY_scatter':{'regions_BED':regions,'DNA':'median paired segment differences; only fully contained bins','RNA':'mean saved expressed non-PAR gene effects; region midpoint assignment','scale':'linear spacing on both DNA and RNA axes','DNA_limits':DNA_SCATTER_XLIM,'DNA_ticks':DNA_SCATTER_XTICKS,'RNA_limits':RNA_SCATTER_YLIM,'RNA_ticks':RNA_SCATTER_YTICKS,'statistics':correlations},
 'bar_scale':'linear, common limits for both cell lines; DNA and RNA measure different quantities'},indent=2)+'\n',encoding='utf-8')
print('Saved alternative DNA-track / paired-arm-bar display: 41 arms per cell line.')
