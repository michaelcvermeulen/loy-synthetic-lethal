"""Two compact columns using the approved HT1080 profile conventions."""
from pathlib import Path
import json
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from make_cell_line_profiles import THESIS, BLUE, apply_style, chromosome_diagram, ordered_y_genes

STEM='ch03_main_a549_du145_profiles'
OUT=THESIS/'figures/chapter3'
SOURCE=OUT/'source_data/a549_du145_profiles'
BACKGROUND={
 'A549':{'disease':'Lung adenocarcinoma','catalogue':'ATCC CCL-185  |  RRID:CVCL_0023',
   'origin':'Primary lung tumour; 58-year-old man (1972)',
   'growth':'Adherent epithelial cells; doubling time ~20\u201330 h',
   'karyotype':'Hypotriploid; modal count 65, XXY (DSMZ)',
   'genotype':['KRAS G12S; STK11 Q37*','TP53 wild-type (reported in ATCC-derived cells)']},
 'DU145':{'disease':'Prostate carcinoma','catalogue':'ATCC HTB-81  |  RRID:CVCL_0105',
   'origin':'Brain metastasis; 69-year-old man (1975)',
   'growth':'Adherent epithelial cells; doubling time ~30\u201340 h',
   'karyotype':'Hypotriploid; modal count 62 (DSMZ)',
   'genotype':['TP53 P223L / V274F; RB1 K715*','CDKN2A D84Y']}}

def main():
    apply_style(base_size=8)
    coords=ordered_y_genes(SOURCE)
    genes=coords.gene.tolist()
    expression=pd.read_csv(SOURCE/'y_expression.csv')
    copy=pd.read_csv(SOURCE/'y_gene_wgs_cn.csv')
    deps=pd.read_csv(SOURCE/'selected_dependencies.csv')
    segments=pd.read_csv(SOURCE/'genome_wgs_segments.csv')
    lengths=pd.read_csv(SOURCE/'chromosome_lengths_grch38.csv')
    meta=json.loads((SOURCE/'provenance.json').read_text())
    fig=plt.figure(figsize=(183/25.4,205/25.4))
    fig.add_artist(Line2D([.50,.50],[.035,.97],transform=fig.transFigure,color='#DDDDDD',lw=.65))
    dep_lim=(np.floor(min(deps.gene_effect.min(),deps.reference_q25.min())*2)/2-.15,
             max(.2,deps.reference_q75.max()+.15))
    ymax=np.ceil(segments.estimated_copy_number.max())
    for x,line,letter in zip([.045,.555],['A549','DU145'],'ab'):
        info=BACKGROUND[line]
        fig.text(x,.972,letter,fontsize=12,weight='semibold',va='top')
        fig.text(x+.032,.974,line,fontsize=13,weight='semibold',va='top')
        fig.text(x,.935,info['disease'],fontsize=9,weight='semibold',va='top')
        fig.text(x,.911,info['catalogue'],fontsize=6.8,color='#666666',va='top')
        for y,key in zip([.876,.852,.828],['origin','growth','karyotype']):
            fig.text(x,y,info[key],fontsize=7.4,va='top')
        fig.text(x,.791,'Reported genotype',fontsize=8,weight='semibold',va='top')
        for i,txt in enumerate(info['genotype']):
            fig.text(x,.765-i*.023,txt,fontsize=7.5,va='top')
        fig.text(x,.696,'Y-linked expression and copy number',fontsize=8.5,weight='semibold')
        chromosome_diagram(fig,[x-.025,.425,.115,.235],SOURCE,coords)
        exp=expression.loc[expression.cell_line.eq(line)].set_index('gene').loc[genes]
        cn=copy.loc[copy.cell_line.eq(line)].set_index('gene').loc[genes]
        ax=fig.add_axes([x+.16,.425,.19,.235])
        yy=np.arange(len(genes))
        ax.barh(yy,exp.expression,color=BLUE,height=.62)
        ax.set_yticks(yy,genes,fontsize=7.1)
        ax.set_ylim(8.5,-.5)
        ax.set_xlim(0,10.5)
        ax.set_xticks([0,5,10])
        ax.set_xlabel('RNA expression [log\u2082(TPM + 1)]',fontsize=7)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y',length=0,pad=3)
        ax.tick_params(axis='x',labelsize=7)
        for i,value in enumerate(exp.expression):
            ax.text(value+.12,i,f'{value:.2f}',va='center',fontsize=6.6)
        cnax=fig.add_axes([x+.366,.425,.065,.235])
        cnax.set_xlim(0,1);cnax.set_ylim(8.5,-.5);cnax.set_axis_off()
        cnax.axvspan(0,1,color='#F3F3F3')
        cnax.text(.5,1.014,'WGS CN\n(est. copies)',transform=cnax.transAxes,fontsize=6.5,ha='center',va='bottom')
        for i,value in enumerate(cn.estimated_copy_number):
            label='<0.001' if value<.001 else f'{value:.2f}'
            cnax.text(.5,i,label,ha='center',va='center',fontsize=7)
        if line=='A549':
            fig.text(x,.370,'Public RNA/WGS profiles disagree on Y status; sampled separately.',
                     fontsize=6.5,color='#555555')
        fig.text(x,.350,'Selected genetic dependencies',fontsize=8.5,weight='semibold')
        d=deps.loc[deps.cell_line.eq(line)]
        assert d.n_reference.eq(1207).all()
        ax=fig.add_axes([x+.095,.207,.32,.116])
        for i,r in enumerate(d.itertuples()):
            ax.plot([r.reference_q25,r.reference_q75],[i,i],color='#BDBDBD',lw=3,solid_capstyle='butt')
        ax.scatter(d.reference_median,np.arange(6),marker='|',s=85,color='#777777',zorder=3)
        ax.scatter(d.gene_effect,np.arange(6),s=22,color=BLUE,zorder=4)
        for i,value in enumerate(d.gene_effect):
            ax.annotate(f'{value:.2f}',(value,i),xytext=(5,0),textcoords='offset points',
                        fontsize=6.4,va='center')
        ax.set_yticks(np.arange(6),d.gene,fontsize=7.3)
        ax.set_ylim(5.6,-.6);ax.set_xlim(*dep_lim)
        ax.axvline(0,color='#CCCCCC',lw=.6)
        ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
        ax.tick_params(axis='x',labelsize=7)
        ax.set_xlabel('Chronos gene effect',fontsize=7)
        handles=[Line2D([],[],marker='o',color=BLUE,ls='',markersize=4,label=line),
                 Line2D([],[],color='#BDBDBD',lw=3,marker='|',markeredgecolor='#777777',
                        label='Others: median / IQR (n = 1,207)')]
        fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(x,.342),ncol=2,
                   fontsize=6.5,frameon=False,borderaxespad=0,handlelength=1.4,
                   handletextpad=.5,columnspacing=1.2)
        fig.text(x,.139,'Genome-wide WGS copy number',fontsize=8.1,weight='semibold')
        fig.text(x+.431,.120,f"Ploidy: {meta['models'][line]['matched_ploidy']:.2f}",fontsize=6.7,ha='right',color='#555555')
        ax=fig.add_axes([x+.04,.049,.391,.065])
        ss=segments.loc[segments.cell_line.eq(line)]
        offset=0;centres=[];labels=[]
        for i,r in enumerate(lengths.itertuples()):
            s=ss.loc[ss.Chrom.eq(r.chromosome)]
            colour=BLUE if r.chromosome=='chrY' else ['#555555','#999999'][i%2]
            if i%2==0:ax.axvspan(offset,offset+r.length_bp,color='#F5F5F5',lw=0)
            ax.hlines(s.estimated_copy_number,s.Start+offset,s.End+offset,colors=colour,lw=.65)
            ax.scatter((s.Start+s.End)/2+offset,s.estimated_copy_number,s=.8,color=colour)
            centres.append(offset+r.length_bp/2);labels.append(r.chromosome.replace('chr',''))
            offset+=r.length_bp+12e6
        ax.set_xlim(0,offset-12e6);ax.set_ylim(-.2,ymax+.3)
        ax.set_yticks([0,5,10] if ymax>=10 else [0,2,4,6,8])
        ax.set_xticks(centres,labels,fontsize=5.5,rotation=60,ha='right')
        ax.tick_params(axis='x',length=0,pad=4)
        ax.tick_params(axis='y',labelsize=6.3,pad=2)
        ax.set_ylabel('Est. copies',fontsize=6.5,labelpad=2)
        ax.set_xlabel('Chromosome (GRCh38)',fontsize=6.8,labelpad=4)
        ax.axhline(2,color='#BDBDBD',lw=.5,ls='--')
        ax.spines['bottom'].set_visible(False)
    for ext in ['pdf','svg','png']:
        target=OUT/f'{STEM}.{ext}'
        stage=OUT/f'{STEM}.staging.{ext}'
        fig.savefig(stage,dpi=600)
        for attempt in range(10):
            try:
                stage.replace(target)
                break
            except OSError:
                if attempt==9:raise
                time.sleep(.5)
    target=OUT/'previews'/f'{STEM}.png'
    stage=target.with_name(f'{STEM}.staging.png')
    fig.savefig(stage,dpi=180)
    stage.replace(target)
    plt.close(fig)
    (SOURCE/'background.json').write_text(json.dumps(BACKGROUND,indent=2)+'\n',encoding='utf-8')
    print(STEM)

if __name__=='__main__':main()
