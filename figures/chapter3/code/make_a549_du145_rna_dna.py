"""Current thesis RNA effects and saved lpWGS segments; no DE or segmentation refit."""
from pathlib import Path
import sys, json, hashlib, io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
from a549_du145_y_arm_extension import extend

FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
ROOT = THESIS.parent
SRC = FIG/'source_data/a549_du145_rna_dna'
SRC.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
STEM = 'ch03_supp_a549_du145_rna_dna_copy_number'
AT = ROOT/'RNAseq/analysis/A549/SEPT_18_2026_REANALYZE/outputs/tables'
DT = ROOT/'RNAseq/analysis/LNCAP/SEPT_23_2026_REANALYZE/outputs/tables'
inputs = [AT/'DE_unpaired_voom_exploratory.csv', DT/'CN_gene_level.csv',
          DT/'CN_window_10Mb_RNA_vs_DNA.csv', DT/'CN_arm_level_RNA_vs_DNA.csv',
          ROOT/'WGS/T2T_genes.gz', ROOT/'WGS/chm13v2.0_cytobands_allchrs.bed']
cb = pd.read_csv(inputs[-1], sep='\t', header=None, names=['chrom','start','end','band','stain'])
pends = cb[cb.band.str.startswith('p')].groupby('chrom').end.max()
qstarts = cb[cb.band.str.startswith('q')].groupby('chrom').start.min()
chroms = ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']
lengths = cb.groupby('chrom').end.max()
offsets = {}; offset = 0
for ch in chroms:
    offsets[ch] = offset
    offset += lengths[ch]/1e6 + 20

def units(d, ch, mid):
    d['arm'] = np.where(mid<=ch.map(pends), ch+'p', np.where(mid>=ch.map(qstarts),ch+'q',None))
    d['window'] = ch+':'+(np.floor(mid/1e7)*10).astype('Int64').astype(str)+'Mb'

profiles = {}
for line, ctrl, loy in [('A549','A549-CTRL2','A549-LOY'),('DU145','LNCaP-sgCTRL','LNCaP-LOY-2X')]:
    path = FIG/f'source_data/{line.lower()}_lpwgs_confirmation/qdnaseq_bins.csv'
    inputs.append(path)
    b = pd.read_csv(path)
    b = b[b['use'] & b['sample'].isin([ctrl,loy]) & np.isfinite(b.log2_seg)]
    b = b.pivot(index=['chromosome','start','end'],columns='sample',values='log2_seg').dropna().reset_index()
    b['DNA'] = b[loy]-b[ctrl]
    b['mid'] = (b.start+b.end)/2
    units(b,b.chromosome,b.mid)
    b.to_csv(SRC/f'{line}_DNA_bins.csv',index=False)
    profiles[line] = b

# Match the DU145 mapping rules: stable ENSG on the same chromosome, then unique symbol.
t = pd.read_csv(ROOT/'WGS/T2T_genes.gz',sep='\t',usecols=[0,1,2,12,20],low_memory=False)
t.columns = ['chrom','start','end','symbol','source_gene']
t = t[t.source_gene.str.startswith('ENSG',na=False)].copy()
t['stable'] = t.source_gene.str.split('.').str[0]
t = t.groupby(['chrom','stable'],sort=False).agg(t2t_start=('start','min'),t2t_end=('end','max'),symbol=('symbol','first')).reset_index()
a = pd.read_csv(AT/'DE_unpaired_voom_exploratory.csv')
assert a.chrom.isin(chroms[:-1]).all() and not a.gene_id.duplicated().any()
a['stable'] = a.gene_id.str.split('.').str[0]
a = a.merge(t,on=['chrom','stable'],how='left',validate='many_to_one')
a['mapping'] = np.where(a.t2t_start.notna(),'stable_ENSG_same_chrom','unmapped')
unique = t[~t.duplicated(['chrom','symbol'],keep=False)].set_index(['chrom','symbol'])
for i in a.index[a.t2t_start.isna()]:
    key = (a.at[i,'chrom'],a.at[i,'gene_name'])
    if key in unique.index:
        z = unique.loc[key]
        a.loc[i,['t2t_start','t2t_end','mapping']] = [z.t2t_start,z.t2t_end,'symbol_same_chrom_unique']
a['DNA'] = np.nan
for ch,g in a[a.t2t_start.notna()].groupby('chrom'):
    bb = profiles['A549'].query('chromosome == @ch')
    for i,r in g.iterrows():
        hit = bb[(bb.start<=r.t2t_end)&(bb.end>=r.t2t_start)]
        if len(hit):
            # Difference of per-sample medians, as used by DU145.
            a.at[i,'DNA'] = hit['A549-LOY'].median()-hit['A549-CTRL2'].median()
a['mid'] = (a.t2t_start+a.t2t_end)/2
units(a,a.chrom,a.mid)
a['RNA'] = a.mean_difference
a.to_csv(SRC/'A549_gene_mapping_and_effects.csv',index=False)
tables = {'A549':{'gene':a[['gene_id','DNA','RNA']].dropna()}}
for level,key,min_genes in [('window_10Mb','window',10),('arm','arm',20)]:
    r = a[a.t2t_start.notna()].groupby(key).agg(RNA=('RNA','mean'),genes=('gene_id','size'))
    b = profiles['A549'].groupby(key).agg(DNA=('DNA','median'),bins=('DNA','size'))
    z = r.join(b).dropna().query('genes >= @min_genes').reset_index().rename(columns={key:'unit'})
    tables['A549'][level] = z
g = pd.read_csv(DT/'CN_gene_level.csv')
tables['DU145'] = {'gene':g[['gene_id','CN_delta_LOY2X','mean_difference']].rename(columns={'CN_delta_LOY2X':'DNA','mean_difference':'RNA'}).dropna()}
for level,filename in [('window_10Mb','CN_window_10Mb_RNA_vs_DNA.csv'),('arm','CN_arm_level_RNA_vs_DNA.csv')]:
    tables['DU145'][level] = pd.read_csv(DT/filename).rename(columns={'DNA_delta_LOY2X':'DNA','RNA_mean_difference':'RNA'})
stats = []
for line,ts in tables.items():
    for level,z in ts.items():
        if 'unit' in z: assert not z.unit.str.startswith('chrY').any()
        assert z[['DNA','RNA']].notna().all().all()
        stats.append(dict(cell_line=line,level=level,n=len(z),Pearson=z.DNA.corr(z.RNA)))
        z.to_csv(SRC/f'{line}_{level}.csv',index=False)
stats = pd.DataFrame(stats)
old = pd.read_csv(DT/'CN_RNA_DNA_concordance.csv').set_index('level')
for r in stats[stats.cell_line.eq('DU145')].itertuples():
    assert np.isclose(r.Pearson,old.loc[r.level,'Pearson'],atol=1e-12)
stats.to_csv(SRC/'correlations.csv',index=False)

# Record DNA-only arm summaries too, including gene-poor regions omitted from correlations.
allarms = []
for line,b in profiles.items():
    z = b.groupby('arm').agg(DNA=('DNA','median'),bins=('DNA','size')).reset_index()
    z['cell_line'] = line
    z['RNA_genes'] = z.arm.map(tables[line]['arm'].set_index('unit').genes).fillna(0).astype(int)
    allarms.append(z)
pd.concat(allarms).to_csv(SRC/'DNA_arm_summary.csv',index=False)
tracks,armplots,extension_report=extend(FIG,SRC,AT,DT,cb,t,units,profiles,tables,inputs)

fig = plt.figure(figsize=(11.5,8.0))
gs = fig.add_gridspec(4,3,height_ratios=[.75,1.25,.75,1.25],hspace=.60,wspace=.34,
                      left=.075,right=.985,bottom=.075,top=.97)
for row,line in enumerate(['A549','DU145']):
    ax = fig.add_subplot(gs[2*row,:]); b=tracks[line]
    for i,ch in enumerate(chroms):
        z=b[b.chromosome.eq(ch)]
        ax.scatter(z.mid/1e6+offsets[ch],z.DNA,s=7,c='#E64B35' if ch=='chrY' else ['#666666','#AAAAAA'][i%2],lw=0)
    ax.axhline(0,c='#333333',lw=.6)
    ax.set_xticks([offsets[c]+lengths[c]/2e6 for c in chroms],[c[3:] for c in chroms],fontsize=7)
    ax.set_xlim(-10,offset)
    ax.set_yscale('symlog',linthresh=.3)
    extent=max(abs(z.DNA).max() for z in tracks.values())*1.1
    ax.set_ylim(-extent,extent)
    ax.set_yticks([-10,-1,0,1,10],['-10','-1','0','1','10'])
    ax.set_ylabel('DNA log$_2$ difference',fontsize=8)
    ax.text(-.065,1.02,'a' if row==0 else 'e',transform=ax.transAxes,weight='semibold',fontsize=12)
    ax.text(.01,.08,line,transform=ax.transAxes,weight='semibold',fontsize=9)
    yq=cb[cb.chrom.eq('chrY') & cb.band.eq('q12')].iloc[0]
    ax.axvspan(offsets['chrY']+yq.start/1e6,offsets['chrY']+yq.end/1e6,color='#EEEEEE',zorder=-1)
    ax.text(.01,.91,'5-Mb windows',transform=ax.transAxes,va='top',fontsize=7)
    ax.annotate('Yq12 excluded',xy=(offsets['chrY']+(yq.start+yq.end)/2e6,0),xytext=(.99,.82),textcoords='axes fraction',ha='right',fontsize=7,color='#777777',arrowprops=dict(arrowstyle='-',lw=.5,color='#999999'))
    for j,(level,label) in enumerate([('gene','Genes'),('window_10Mb','10-Mb windows'),('arm','Chromosome arms')]):
        ax=fig.add_subplot(gs[2*row+1,j]); z=armplots[line] if level=='arm' else tables[line][level]
        ax.axhline(0,c='#DDDDDD',lw=.6);ax.axvline(0,c='#DDDDDD',lw=.6)
        ax.scatter(z.DNA,z.RNA,s=2 if j==0 else 12,alpha=.3 if j==0 else .75,c='#555555',lw=0,rasterized=j==0)
        r=z.DNA.corr(z.RNA)
        stat=f'{label}\nr = {r:.3f}; n = {len(z):,}'
        if level=='arm':stat=f'{label}\nr (all) = {r:.3f}; n = {len(z)}\nr (non-Y) = {tables[line][level].DNA.corr(tables[line][level].RNA):.3f}'
        ax.text(.04,.96,stat,transform=ax.transAxes,va='top',fontsize=7.5,bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1))
        ax.text(-.22,1.04,chr(ord('b')+row*4+j),transform=ax.transAxes,fontsize=12,weight='semibold')
        ax.set_xlabel('DNA log$_2$ difference')
        ax.set_ylabel('RNA mean log$_2$ CPM difference')
        both=pd.concat([tables[x][level] for x in tables])
        xmin,xmax=both.DNA.min(),both.DNA.max(); ymin,ymax=both.RNA.min(),both.RNA.max()
        ax.set_xlim(xmin-.08*(xmax-xmin),xmax+.08*(xmax-xmin))
        if level=='gene':
            ax.set_xscale('symlog',linthresh=.25)
            ax.set_xticks([-1,0,1,10],['-1','0','1','10'])
        ax.set_ylim(ymin-.1*(ymax-ymin),ymax+.22*(ymax-ymin))
        if level=='arm':
            ax.set_xscale('symlog',linthresh=.3);ax.set_yscale('symlog',linthresh=.3)
            ax.set_xlim(-12,1);ax.set_ylim(-12,1)
            ticks=[-10,-1,-.3,0,.3,1];ax.set_xticks(ticks,[str(v) for v in ticks]);ax.set_yticks(ticks,[str(v) for v in ticks]);ax.tick_params(labelsize=7)
            ys=z.unit.str.startswith('chrY');ax.scatter(z.loc[ys,'DNA'],z.loc[ys,'RNA'],s=25,color='#E64B35',lw=0,zorder=4)
            texts=[]
            labels=['chrYp','chrYq']+(['chr20p','chrXp','chrXq','chr15q'] if line=='DU145' else [])
            for unit in labels:
                q=z[z.unit.eq(unit)].iloc[0]
                texts.append(ax.text(q.DNA,q.RNA,unit[3:],fontsize=7,color='#E64B35' if unit.startswith('chrY') else '#333333'))
            np.random.seed(42);adjust_text(texts,ax=ax,expand=(1.3,1.4),arrowprops=dict(arrowstyle='-',color='#888888',lw=.5))
        if level=='window_10Mb' and line=='A549':
            q=z.loc[z.DNA.idxmin()]
            ax.annotate(q.unit.replace('chr',''),(q.DNA,q.RNA),xytext=(8,-13),textcoords='offset points',fontsize=7)
        if level=='window_10Mb' and line=='DU145':
            score=(((z[['DNA','RNA']]-z[['DNA','RNA']].median())/z[['DNA','RNA']].std())**2).sum(axis=1)
            selected=list(dict.fromkeys(score.nlargest(4).index.tolist()+[z.DNA.idxmax(),z.RNA.idxmax()]))
            texts=[]; label_rows=[]
            for i in selected:
                q=z.loc[i];chrom,position=q.unit.split(':')
                left=int(position[:-2])*1_000_000;right=min(left+10_000_000,int(lengths[chrom]))
                bands=cb[cb.chrom.eq(chrom)&cb.start.lt(right)&cb.end.gt(left)].sort_values('start')
                first,last=bands.band.iloc[0],bands.band.iloc[-1]
                label=chrom[3:]+first+('–'+last if first!=last else '')
                texts.append(ax.text(q.DNA,q.RNA,label,fontsize=6.5,color='#3C5488'))
                ax.scatter([q.DNA],[q.RNA],s=16,color='#3C5488',lw=0,zorder=4)
                label_rows.append(dict(unit=q.unit,BED_start=left,BED_end=right,label=label,DNA=q.DNA,RNA=q.RNA))
            np.random.seed(42)
            adjust_text(texts,ax=ax,expand=(1.3,1.5),force_text=(.5,.7),arrowprops=dict(arrowstyle='-',color='#888888',lw=.5))
            pd.DataFrame(label_rows).to_csv(SRC/'DU145_10Mb_cytoband_labels.csv',index=False)
for ax in fig.axes:
    ax.spines[['top','right']].set_visible(False)
for ext in ['pdf','svg','png']:
    buf=io.BytesIO();fig.savefig(buf,format=ext,dpi=220)
    (FIG/f'{STEM}.{ext}').write_bytes(buf.getvalue())
plt.close(fig)
provenance={'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
 'RNA':'saved mean difference in TMM log2 CPM, LOY minus control; expressed autosomal/non-PAR X genes',
 'DNA':'saved T2T QDNAseq segment log2 difference; usable bins; no segmentation refit',
 'mapping':'stable ENSG same chromosome then unique symbol; gene spans mapped to overlapping bins',
 'aggregation':'gene: difference of sample medians; windows/arms: median bin DNA difference versus mean gene RNA effect',
 'minimum_genes':{'window_10Mb':10,'arm':20},'DU145_saved_correlations_reproduced':True,
 'scope':'descriptive population comparisons; RNA and lpWGS specimens not established as matched; genomic units not independent replicates'}
(SRC/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
(SRC/'8panel_Y_extension_provenance.json').write_text(json.dumps(extension_report,indent=2)+'\n',encoding='utf-8')
print(stats.to_string(index=False))
