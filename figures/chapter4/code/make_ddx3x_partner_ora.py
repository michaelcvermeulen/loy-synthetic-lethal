"""Show representative ORA terms and the actual contributing DDX3X partners."""
from pathlib import Path
import json
import sys
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import pymupdf as fitz

THESIS = Path(__file__).resolve().parents[3]
OUT = THESIS/'figures/chapter4'
DATA = OUT/'source_data/ddx3x_partner_ora'
STEM = 'ch04_supp_ddx3x_partner_overrepresentation_draft'
sys.path.insert(0,str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
PALETTE = {'Reactome':'#BC8A26','GO cellular component':'#00857A','GO biological process':'#7A3E85'}
LABELS = {
    'REACTOME_SYNTHESIS_OF_DNA':'Synthesis of DNA',
    'GOCC_PROTEASOME_ACCESSORY_COMPLEX':'Proteasome accessory complex',
    'GOCC_PROTEASOME_REGULATORY_PARTICLE_BASE_SUBCOMPLEX':'Proteasome regulatory particle, base subcomplex',
    'REACTOME_DENGUE_VIRUS_INFECTION':'Dengue virus infection',
    'REACTOME_PROCESSING_OF_CAPPED_INTRON_CONTAINING_PRE_MRNA':'Processing of capped intron-containing pre-mRNA',
    'GOCC_CATALYTIC_STEP_2_SPLICEOSOME':'Catalytic step 2 spliceosome',
    'GOCC_U2_TYPE_SPLICEOSOMAL_COMPLEX':'U2-type spliceosomal complex',
    'GOBP_MRNA_PROCESSING':'mRNA processing',
    'GOCC_NUCLEAR_DNA_DIRECTED_RNA_POLYMERASE_COMPLEX':'Nuclear DNA-directed RNA polymerase complex',
    'GOBP_RNA_EXPORT_FROM_NUCLEUS':'RNA export from nucleus',
    'GOBP_NUCLEAR_TRANSPORT':'Nuclear transport',
    'REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX':'Activation of the pre-replicative complex',
}


def main():
    apply_style(9)
    data = pd.read_csv(DATA/'display_terms.csv')
    assert len(data)>0 and data.fdr_bh_global.lt(.05).all()
    width,height = 12.4,6.5
    fig = plt.figure(figsize=(width,height))
    grid = fig.add_gridspec(1,4,left=.018,right=.985,top=.87,bottom=.12,
                           width_ratios=[3.45,1.85,.90,4.45],wspace=.10)
    axes = [fig.add_subplot(grid[0,i]) for i in range(4)]
    for ax in axes:
        ax.set_ylim(len(data)-.5,-.5)
        ax.set_yticks([])
        for row in range(len(data)):
            if row%2==0:
                ax.axhspan(row-.5,row+.5,color='#f4f4f4',zorder=0)
    for i in [0,2,3]:
        axes[i].set_xlim(0,1)
        axes[i].axis('off')
    scores = -np.log10(data.fdr_bh_global.to_numpy())
    ax = axes[1]
    bars = ax.barh(np.arange(len(data)),scores,height=.42,
                   color=data.source.map(PALETTE),edgecolor='none',zorder=3)
    ax.set_xlim(0,max(scores)*1.06)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.set_xlabel('−log₁₀(FDR)',labelpad=6)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y',left=False)
    for row,term in enumerate(data.itertuples()):
        label = LABELS.get(term.pathway,term.pathway.split('_',1)[-1].replace('_',' ').capitalize())
        axes[0].text(.99,row,textwrap.fill(label,43),ha='right',va='center',fontsize=8.5,
                     linespacing=1.15)
        axes[2].text(.5,row,f'{term.fdr_bh_global:.2g}',ha='center',va='center',fontsize=8.5)
        genes = term.overlap_genes.split(';')
        assert len(genes)==term.overlap_n
        text = f'{term.overlap_n}/31: '+', '.join(genes)
        axes[3].text(.02,row,textwrap.fill(text,64),ha='left',va='center',fontsize=8.1,linespacing=1.15)
        assert abs(bars[row].get_width()-scores[row])<1e-12
    for i,title in [(0,'Term'),(1,'Enrichment significance'),(2,'BH FDR'),(3,'Contributing genes')]:
        axes[i].text(.99 if i==0 else (.02 if i==3 else .5),1.025,title,transform=axes[i].transAxes,
                     ha='right' if i==0 else ('left' if i==3 else 'center'),va='bottom',
                     fontsize=9,fontweight='semibold')
    handles = [Line2D([],[],marker='s',linestyle='none',color=c,label=s,markersize=6)
               for s,c in PALETTE.items()]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.995),ncol=3,
               frameon=False,fontsize=9,columnspacing=2)
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=155)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page = doc[0]
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b
                 for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        outside = [s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        overlaps = []
        for i,a in enumerate(spans):
            for b in spans[i+1:]:
                ra,rb = fitz.Rect(a['bbox']),fitz.Rect(b['bbox'])
                area = (ra&rb).get_area()
                if area>1 and area/min(ra.get_area(),rb.get_area())>.1:
                    overlaps.append([a['text'],b['text']])
        assert not outside,outside
        assert not overlaps,overlaps
        page.get_pixmap(matrix=fitz.Matrix(1.4,1.4),colorspace=fitz.csGRAY).save(
            OUT/'review'/f'{STEM}_greyscale.png')
    manifest = dict(stem=STEM,font='Open Sans',dimensions_inches=[width,height],
                    displayed_terms=data.pathway.tolist(),colour_mapping=PALETTE,
                    labels={p:LABELS.get(p,p) for p in data.pathway},
                    validation={'bar_values_match':True,'all_overlap_counts_match_gene_lists':True,
                                'text_outside_page':outside,'text_overlaps':overlaps})
    (DATA/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest['validation'],indent=2))


if __name__=='__main__':
    main()
