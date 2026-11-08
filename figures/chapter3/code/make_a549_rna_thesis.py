"""Restyle frozen A549 outputs; no model fitting or hypothesis testing."""
from pathlib import Path
import sys, json, hashlib, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'figures/chapter3'
SOURCE = ROOT.parent / 'RNAseq/analysis/A549/SEPT_18_2026_REANALYZE/outputs/tables'
DATA = OUT / 'source_data/a549_rna_thesis'
DATA.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
FONT = apply_style(8)
COLOURS = {'CTRL':'#3C5488', 'LOY':'#E64B35', 'parental':'#777777', 'PARENTAL':'#777777'}
LABELS = {'CTRL':'sgCTRL', 'LOY':'Y-depleted', 'parental':'Parental', 'PARENTAL':'Parental'}
used = {}
outputs = []

def read(name):
    p = SOURCE / (name+'.csv')
    shutil.copy2(p, DATA / p.name)
    used[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p)

def panel(ax, letter, title):
    ax.set_title(title, loc='left', pad=10)
    ax.text(-.13, 1.07, letter, transform=ax.transAxes, fontweight='semibold', fontsize=11)

def save(fig, stem):
    for ext in ['pdf','svg','png']:
        p = OUT / (stem+'.'+ext)
        fig.savefig(p, dpi=300, facecolor='white')
        outputs.append({'file':p.name, 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    plt.close(fig)

def legend(ax, groups=('CTRL','LOY'), **kwargs):
    use_shapes=kwargs.pop('abundance',False)
    ax.legend(handles=[Line2D([],[],color=COLOURS[g],marker='^' if g=='LOY' and use_shapes else 'o',linestyle='',label=LABELS[g]) for g in groups], **kwargs)

def pca(ax, d, xlabel, ylabel):
    for group, s in d.groupby('treatment', sort=False):
        col=COLOURS.get(group, '#777777')
        ax.scatter(s.PC1,s.PC2,c=col,s=33,edgecolor='white',linewidth=.4)
        for r in s.itertuples():
            ax.annotate(str(r.sample_index),(r.PC1,r.PC2),xytext=(5,4),textcoords='offset points',fontsize=7,color=col)
    ax.set(xlabel=xlabel,ylabel=ylabel)
    ax.margins(.18)

p = read('PCA_core')
if '--qc-only' not in sys.argv:
    import runpy
    runpy.run_path(str(OUT/'code/make_a549_volcano_PAR_RINe.py'),run_name='__main__')
    xy=read('XY_paralog_sample_abundance'); effects=read('XY_paralog_effects_and_length_proxy')
    par=read('PAR_selected_gene_sample_abundance'); par_de=read('DE_PAR_secondary_voom')
    
    def abundance(ax, data, genes, relative=False):
        for i,gene in enumerate(genes):
            d=data[data.gene_name==gene]
            denominator=d.loc[d.treatment=='CTRL','CPM'].mean() if relative else 1
            for group,offset,marker in [('CTRL',-.16,'o'),('LOY',.16,'^')]:
                s=d[d.treatment==group].sort_values('sample_index')
                vals=s.CPM.to_numpy()/denominator*100 if relative else np.log2(s.CPM.to_numpy()+1)
                xx=i+offset+np.linspace(-.07,.07,len(vals))
                ax.scatter(xx,vals,s=20,c=COLOURS[group],marker=marker,edgecolor='white',lw=.3,zorder=3)
                ax.plot([i+offset-.11,i+offset+.11],[np.mean(vals)]*2,color=COLOURS[group],lw=1.3)
        ax.set_xticks(range(len(genes)),genes,rotation=45,ha='right',fontstyle='italic')
        ax.set_xlim(-.55,len(genes)-.45)
        ax.set_ylabel('Expression (% of mean sgCTRL)' if relative else 'log₂(CPM + 1)')
        if relative: ax.axhline(100,color='.7',lw=.6,ls='--',zorder=0)
        ax.set_ylim(bottom=0)
    
    from a549_par_heatmap import draw as draw_par_heatmap
    fig=draw_par_heatmap(par,xy,effects,par_de)
    save(fig,'ch03_main_a549_rna_sex_chromosomes')
    
    # Pathways are already included in F01 by the delegated combined renderer.
    # The standalone F03 remains historical.
    # Complete Hallmark enrichment supplement retired; rebuild its ssGSEA replacement with make_a549_hallmark_ssgsea.py.
    

# Expanded QC uses provider integrity and the recent sequencing-QC archive.
import runpy
runpy.run_path(str(OUT/'code/make_a549_qc_expanded.py'),run_name='__main__')

# Retired S05 is generated only on explicit archival request.
if '--include-retired-par' in sys.argv and '--qc-only' not in sys.argv:
    genes=par.gene_name.drop_duplicates().tolist()
    fig,axes=plt.subplots(4,4,figsize=(7.1,8.1),layout='constrained')
    for a,gene in zip(axes.flat,genes):
        abundance(a,par,[gene],True);a.set_xticks([]);a.set_title(gene,fontstyle='italic');a.set_ylabel('% of sgCTRL',fontsize=7)
        row=par_de[par_de.gene_name==gene].iloc[0]
        a.text(.98,.98,f'PAR adj. P\n{row.BH_PAR_family:.3f}',transform=a.transAxes,ha='right',va='top',fontsize=6.5)
        a.set_ylim(0,a.get_ylim()[1]*1.22)
    for a in list(axes.flat)[len(genes):]: a.set_visible(False)
    legend(list(axes.flat)[len(genes)-1],loc='upper left',bbox_to_anchor=(1.3,1),abundance=True)
    save(fig,'ch03_supp_a549_rna_par_expression')
    

(DATA/('provenance_qc_wrapper.json' if '--qc-only' in sys.argv else 'provenance.json')).write_text(json.dumps({'source':str(SOURCE),'font':FONT,'input_sha256':used,'exports':outputs,'note':'Frozen analysis restyled without fitting models. Drafts pending author review.'},indent=2),encoding='utf-8')
print(f'Wrote {len(outputs)} exports; font: {FONT}')
