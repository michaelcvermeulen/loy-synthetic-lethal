"""Exploratory period-specific RNA versus fixed ATAC arm summaries; no model refitting."""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from adjustText import adjust_text

R = Path(__file__).resolve().parents[3]
F = R / 'figures/chapter3'
D = F / 'source_data/ht1080_rna_period_atac_arms'
D.mkdir(exist_ok=True)
sys.path.insert(0, str(R / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
out = pd.read_csv(F / 'source_data/ht1080_atac_copy_number/panel_b_arm_correlations.csv').set_index('arm')
g = pd.read_csv(F / 'source_data/ht1080_rna_atac_integration/integration_genes.csv')
ann = pd.read_csv(R.parent / 'ATACseq/SEPT_08_2026_REPROCESS/config/genes_knownGene.tsv', sep='\t').rename(columns={'gene_id':'entrez_id'})
g = g.merge(ann[['entrez_id','start','end']], on='entrez_id', validate='one_to_one')
bands = pd.read_csv(F / 'source_data/ht1080_profile_wgs_cn/reference_material/cytoBand_hg38.txt.gz', sep='\t', header=None, names=['chrom','start','end','band','stain'])
cen = bands[bands.stain.eq('acen')].groupby('chrom').agg(start=('start','min'), end=('end','max'))
mid = (g.start + g.end)/2
g['arm'] = pd.Series(index=g.index, dtype='object')
for chrom in ['chr'+str(i) for i in range(1,23)]:
    a,b = cen.loc[chrom,['start','end']]
    g.loc[g.chrom.eq(chrom) & mid.lt(a),'arm'] = chrom+'p'
    g.loc[g.chrom.eq(chrom) & mid.ge(b),'arm'] = chrom+'q'
for label,col in [('August','AUG_logFC'),('May','MAY_logFC')]:
    out[label+'_RNA'] = g.groupby('arm')[col].median() - g[col].median()
assert out[['August_RNA','May_RNA','ATAC_mean_shift']].notna().all().all()
assert np.allclose(out.RNA_median_logFC, (g.groupby('arm').rna_logFC.median()-g.rna_logFC.median()).reindex(out.index))
out.to_csv(D / 'arm_values.csv')
stats = {}
for label,col in [('August','August_RNA'),('May','May_RNA'),('Combined','RNA_median_logFC')]:
    stats[label] = {}
    for subset,mask in [('all_39',np.ones(len(out),dtype=bool)),('without_chr18',~out.index.str.startswith('chr18'))]:
        x,y = out.loc[mask,col],out.loc[mask,'ATAC_mean_shift']
        s,p = spearmanr(x,y),pearsonr(x,y)
        stats[label][subset] = dict(n=len(x),spearman_rho=float(s.statistic),spearman_p=float(s.pvalue),pearson_r=float(p.statistic),pearson_p=float(p.pvalue))
(D / 'statistics.json').write_text(json.dumps(stats,indent=2))
fig,axes = plt.subplots(1,2,figsize=(7.2,3.3),sharex=True,sharey=True)
fig.subplots_adjust(left=.10,right=.985,bottom=.19,top=.96,wspace=.09)
xs = out[['August_RNA','May_RNA']].to_numpy()
ys = out.ATAC_mean_shift.to_numpy()
xpad = np.ptp(xs)*.14; ypad = np.ptp(ys)*.18
for ax,label,letter in zip(axes,['August','May'],['a','b']):
    x,y = out[label+'_RNA'],out.ATAC_mean_shift
    ax.axhline(0,color='.85',lw=.6,zorder=0); ax.axvline(0,color='.85',lw=.6,zorder=0)
    ax.scatter(x,y,s=19,color='.42',edgecolor='white',linewidth=.35,zorder=3)
    ax.set_xlim(xs.min()-xpad,xs.max()+xpad); ax.set_ylim(ys.min()-ypad,ys.max()+ypad)
    s = stats[label]['all_39']
    ax.text(.04,.96,f"{label} RNA\nSpearman ρ = {s['spearman_rho']:.2f}\nP = {s['spearman_p']:.3f}; n = {s['n']} arms",transform=ax.transAxes,va='top',fontsize=8)
    ax.text(-.10,1.01,letter,transform=ax.transAxes,fontweight='semibold',fontsize=11)
    labels = set(['chr13q','chr18p','chr18q']) | set(x.abs().nlargest(2).index) | set(y.abs().nlargest(2).index)
    texts = [ax.text(x.loc[k],y.loc[k],k.removeprefix('chr'),fontsize=7) for k in sorted(labels)]
    adjust_text(texts,x=x.to_numpy(),y=y.to_numpy(),ax=ax,expand=(1.3,1.5),arrowprops=dict(arrowstyle='-',color='.45',lw=.45))
    ax.set_xlabel('RNA arm median log₂ fold change\n(Y-depleted − sgCTRL; centred)')
    ax.spines[['top','right']].set_visible(False)
axes[0].set_ylabel('ATAC arm mean shift\n(Y-depleted − sgCTRL; centred log₂ CPM)')
stem = 'ch03_supp_ht1080_rna_period_atac_arms_draft'
for ext in ['pdf','svg','png']:
    fig.savefig(F / (stem+'.'+ext),dpi=220)
plt.close(fig)
(D / 'README.md').write_text('''# Exploratory RNA period versus ATAC arm comparison

Reproduce with code/plot_ht1080_rna_period_atac_arms.py. No differential models refitted.
Uses the same 39 complete autosomal arms and saved ATAC summaries as F07 panel b.
RNA input: ht1080_rna_atac_integration/integration_genes.csv, AUG_logFC and MAY_logFC;
same gene universe as existing F07, without differential-expression filtering.
Gene midpoints assign arms, excluding centromeres. Each period's median gene effect
per arm is centred on that period's genome-wide median in this input gene universe.
ATAC is identical in both panels: mean of the three centred passage-level arm medians.
Both panels have identical axes. Correlations use arm summaries, not biological replicates;
P values are two-sided, unadjusted and exploratory. No test of the difference between
correlations is claimed. statistics.json also reports Pearson correlations, the joint RNA
comparison and descriptive sensitivity excluding both chromosome 18 arms.
The comparison does not establish passage/time matching between RNA and ATAC samples.
Font: bundled Open Sans. Draft only; existing F07 and manuscript are unchanged.
''',encoding='utf-8')
print(json.dumps(stats,indent=2))
