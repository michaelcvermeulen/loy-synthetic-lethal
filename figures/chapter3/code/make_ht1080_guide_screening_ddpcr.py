"""Author-adopted five-panel Fig 3.8: HT1080 guide screen and HT1080/A549 ddPCR."""
from pathlib import Path
import sys, json, hashlib, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT/'figures/chapter3'
OUT = FIG
SNAP = FIG/'source_data/ht1080_a549_guide_screening_ddpcr'
SNAP.mkdir(exist_ok=True)
STEM = 'ch03_main_ht1080_guide_screening_ddpcr'
sys.path.insert(0, str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)

inputs = {
    'rtqpcr_values.csv': FIG/'source_data/ht1080_initial_rtqpcr/mean_expression_reduction_three_genes.csv',
    'rna_gene_coordinates.csv': FIG/'source_data/ht1080_initial_rtqpcr/y_gene_coordinates_grch38.csv',
    'ddpcr_gene_coordinates.csv': FIG/'source_data/ht1080_guide_screening_ddpcr/ddpcr_gene_coordinates_grch38.csv',
    'cytobands.csv': FIG/'source_data/ht1080_initial_rtqpcr/chrY_cytobands_grch38.csv',
    'ht1080_measurements.csv': FIG/'source_data/ht1080_ddpcr_nucleofection_rounds/main_measurements.csv',
    'ht1080_summary.csv': FIG/'source_data/ht1080_ddpcr_nucleofection_rounds/main_group_summary.csv',
    'a549_measurements.csv': FIG/'source_data/a549_ddpcr_nucleofection_rounds/main_measurements.csv',
    'a549_summary.csv': FIG/'source_data/a549_ddpcr_nucleofection_rounds/main_group_summary.csv',
}
source_hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs.values()}
for name, p in inputs.items():
    shutil.copy2(p, SNAP/name)

rna = pd.read_csv(SNAP/'rtqpcr_values.csv')
ht = pd.read_csv(SNAP/'ht1080_measurements.csv')
a549 = pd.read_csv(SNAP/'a549_measurements.csv')
assert len(rna) == 15 and len(ht) == 28 and len(a549) == 54
assert not ht.sample_id.isin(['C5', 'B6']).any()
assert ht[ht.plot_group.eq('CP + Q53 2x')].assay.tolist() == ['AMELY/AMELX']
assert not ((a549.sample_id == 'I1') & (a549.assay == 'AMELY/AMELX')).any()
assert a549.included.all()
a549['plot_group'] = a549.condition.replace({'CTRL2': 'sgCTRL', 'MIX/MIX': 'CP + Q53 2x'})

genes = ['RPS4Y1', 'KDM5D', 'EIF1AY']
guides = ['CN', 'CP', 'Q53', 'P47', 'P51']
values = rna.pivot(index='guide', columns='gene', values='reduction_percent').loc[guides, genes].to_numpy()
coords = pd.read_csv(SNAP/'rna_gene_coordinates.csv').set_index('gene').loc[genes]
ddpcr_coords = pd.read_csv(SNAP/'ddpcr_gene_coordinates.csv').set_index('gene')
bands = pd.read_csv(SNAP/'cytobands.csv'); length = bands.end.max()
assert np.isfinite(values).all()

fig = plt.figure(figsize=(183/25.4, 230/25.4))
left, width = .18, .67
fig.text(.025, .978, 'a', fontsize=11, weight='semibold', va='top')
axmap = fig.add_axes([left, .875, width, .09])
axmap.set(xlim=(0, 1), ylim=(0, 1)); axmap.set_axis_off()
shades = {'gneg': '#FFFFFF', 'gpos50': '#B9B9B9', 'gvar': '#EAEAEA'}
for band in bands.itertuples():
    x0, x1 = band.start/length, band.end/length
    if band.stain == 'acen':
        verts = [(x0,.56),(x0,.72),(x1,.64)] if band.band.startswith('p') else [(x0,.64),(x1,.72),(x1,.56)]
        axmap.add_patch(Polygon(verts, fc='#777777', ec='#777777', lw=.6))
    else:
        axmap.add_patch(Rectangle((x0,.56), x1-x0, .16, fc=shades[band.stain], ec='#777777', lw=.6))
axmap.text(0, 1.19, 'HT1080 guide screen', fontsize=8.5, weight='semibold')
axmap.text(1, 1.19, 'Y chromosome | GRCh38', fontsize=7.5, color='#555555', ha='right')
axmap.text(.08, .82, 'Yp', fontsize=7.5, ha='center')
axmap.text(.60, .82, 'Yq', fontsize=7.5, ha='center')
for i, (gene, row) in enumerate(coords.iterrows()):
    gx = ((row.start-1+row.end)/2)/length; target = (i+.5)/3
    axmap.plot([gx,gx], [.55,.74], color='#353535', lw=1.1)
    axmap.plot([gx,gx,target,target], [.55,.49,.32,.26], color='#999999', lw=.7)
for gene, label_x, panels in [('AMELY',.12,'b,d'), ('DDX3Y',.30,'c,e')]:
    row = ddpcr_coords.loc[gene]; gx = ((row.start-1+row.end)/2)/length
    axmap.plot([gx,gx], [.55,.74], color='#353535', lw=1.1)
    axmap.plot([gx,gx,label_x,label_x], [.74,.84,.91,.94], color='#777777', lw=.7)
    axmap.text(label_x,.97,f'{gene} ({panels})',fontsize=7.2,weight='semibold',ha='center',va='bottom')

heat = fig.add_axes([left,.715,width,.14])
cmap = LinearSegmentedColormap.from_list('expression_reduction',['#C75A2B','#FAFAFA','#3B6FB6'])
im = heat.imshow(values,cmap=cmap,norm=TwoSlopeNorm(vmin=-100,vcenter=0,vmax=100),aspect='auto')
heat.set_xticks(range(3),genes,fontsize=8,weight='semibold'); heat.xaxis.tick_top()
heat.set_yticks(range(5),guides,fontsize=8); heat.tick_params(axis='both',length=0,pad=5)
heat.set_ylabel('Guide',labelpad=10)
for spine in heat.spines.values(): spine.set_visible(False)
for x in np.arange(.5,3,1): heat.axvline(x,color='white',lw=2)
for y in np.arange(.5,5,1): heat.axhline(y,color='white',lw=2)
for i in range(5):
    for j in range(3):
        v = values[i,j]
        heat.text(j,i,f'{v:.0f}%',ha='center',va='center',fontsize=8.5,
                  color='white' if abs(v)>65 else '#222222',weight='semibold')
cax = fig.add_axes([.88,.720,.016,.13])
cb = fig.colorbar(im,cax=cax,ticks=[-100,-50,0,50,100])
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6.3,length=2)
cb.set_label('Expression reduction (%)',fontsize=7,labelpad=5)

groups = ['CP','Q53','CP + Q53','CP + Q53 2x','sgCTRL']
markers = {'Transfection experiment 1':'s','Transfection experiment 2':'o','Two successive nucleofections':'D'}
summary = []
axes = []
for line, data, y0, letters in [('HT1080',ht,.420,'bc'), ('A549',a549,.10,'de')]:
    fig.text(.105,y0+.257,line,fontsize=9.5,weight='semibold')
    for x0, assay, letter in zip([.105,.59],['AMELY/AMELX','DDX3Y/RNaseP'],letters):
        ax = fig.add_axes([x0,y0,.38,.215]); axes.append(ax)
        for i, group in enumerate(groups):
            g = data[data.assay.eq(assay) & data.plot_group.eq(group)]
            if line == 'HT1080': g = g.sort_values(['date','export_well'])
            v = g.estimated_target_copies_per_cell.to_numpy()
            if not len(v):
                ax.text(i,.045,'NA',transform=ax.get_xaxis_transform(),ha='center',fontsize=8,color='#777777')
                summary.append(dict(cell_line=line,assay=assay,group=group,n_measurements=0,median=np.nan,q1=np.nan,q3=np.nan))
                continue
            med = float(np.median(v)); q1,q3 = np.quantile(v,[.25,.75])
            color = '#3C5488' if group == 'sgCTRL' else '#E64B35'
            ax.bar(i,med,width=.67,color=color,zorder=2)
            if group == 'CP + Q53 2x':
                ax.bar(i,med,width=.67,facecolor='none',edgecolor='#87301E',linewidth=.7,hatch='///',zorder=3)
            ax.errorbar(i,med,yerr=[[med-q1],[q3-med]],fmt='none',ecolor='#333333',elinewidth=.8,capsize=2,zorder=4)
            offsets = np.linspace(-.18,.18,len(g)) if len(g)>1 else [0]
            for offset, (_, row) in zip(offsets,g.iterrows()):
                marker = markers[row.point_class] if line == 'HT1080' else 'o'
                ax.scatter(i+offset,row.estimated_target_copies_per_cell,s=13,marker=marker,
                           facecolor='white',edgecolor='#222222',linewidth=.65,zorder=5)
            y = max(v.max(),q3)+.058
            if .93 < y < 1.05: y = 1.05
            ax.text(i,y,f'{med:.3f}',ha='center',va='bottom',fontsize=7,weight='semibold')
            ax.text(i,y+.065,f'n={len(v)}',ha='center',va='bottom',fontsize=6.4,color='#666666')
            summary.append(dict(cell_line=line,assay=assay,group=group,n_measurements=len(v),median=med,q1=q1,q3=q3))
            assert np.all((v>=0) & (v<=1.35))
        ax.axhline(1,color='#888888',lw=.8,ls=(0,(4,3)),zorder=1)
        ax.axvline(3.5,color='#DDDDDD',lw=.65)
        ax.set(xlim=(-.55,4.55),ylim=(0,1.35),xticks=range(5),
               xticklabels=['CP\n1×','Q53\n1×','CP + Q53\n1×','CP + Q53\n2×','sgCTRL'],
               yticks=np.arange(0,1.21,.2))
        ax.tick_params(axis='x',length=0,pad=6,labelsize=7)
        ax.tick_params(axis='y',labelsize=7)
        ax.set_title(assay.split('/')[0],pad=9,fontsize=9)
        ax.text(-.19,1.05,letter,transform=ax.transAxes,fontsize=11,weight='semibold')
        ax.grid(axis='y',color='#ECECEC',lw=.45);ax.set_axisbelow(True)
        if x0 == .105: ax.set_ylabel('Estimated Y-linked target copies per cell',fontsize=7.6,labelpad=7)

# Restrict the transfection-experiment shape key explicitly to the HT1080 row.
handles = [Line2D([],[],marker=markers[k],color='none',markerfacecolor='white',markeredgecolor='#222222',
                  markersize=3.6,label=label) for k,label in zip(markers,['Experiment 1','Experiment 2','Two successive nucleofections'])]
fig.text(.105,.037,'HT1080 points:',fontsize=6.7,va='center')
fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.27,.037),ncol=3,frameon=False,
           fontsize=6.4,handletextpad=.35,columnspacing=.9,borderaxespad=0)

summary = pd.DataFrame(summary)
common = ['assay','group','n_measurements','median','q1','q3']
old_ht = pd.read_csv(SNAP/'ht1080_summary.csv')
old_a = pd.read_csv(SNAP/'a549_summary.csv').rename(columns={'condition':'group'})
old_a['group'] = old_a.group.replace({'CTRL2':'sgCTRL','MIX/MIX':'CP + Q53 2x'})
for line, old in [('HT1080',old_ht),('A549',old_a)]:
    new = summary[summary.cell_line.eq(line)][common].reset_index(drop=True)
    pd.testing.assert_frame_equal(new,old[common].reset_index(drop=True),check_dtype=False,rtol=1e-12,atol=1e-12)
summary.to_csv(SNAP/'combined_display_summary.csv',index=False)
assert all(ax.get_ylim() == (0,1.35) for ax in axes)
for ext in ['pdf','svg','png']:
    fig.savefig(OUT/f'{STEM}.{ext}',dpi=600,facecolor='white')
fig.savefig(FIG/'previews'/f'{STEM}.png',dpi=150,facecolor='white')
plt.close(fig)

for p in inputs.values():
    assert hashlib.sha256(p.read_bytes()).hexdigest() == source_hashes[p.relative_to(ROOT).as_posix()]
validation = dict(status='adopted as Fig 3.8, 27 October 2026',dimensions_mm=[183,230],dpi=600,font='Open Sans',
                  RNA_summaries=15,HT1080_ddPCR_measurements=28,A549_ddPCR_measurements=54,
                  both_original_group_summaries_match=True,missing_HT1080_2x_DDX3Y_preserved=True,
                  original_numerical_inputs_unchanged=True,source_sha256=source_hashes,
                  panel_map={'a':'pre-merge Fig 3.8a, compacted','b':'pre-merge Fig 3.8b','c':'pre-merge Fig 3.8c',
                             'd':'former Fig 3.11a','e':'former Fig 3.11b'})
(SNAP/'validation.json').write_text(json.dumps(validation,indent=2)+'\n',encoding='utf-8')
print(json.dumps(validation,indent=2))
