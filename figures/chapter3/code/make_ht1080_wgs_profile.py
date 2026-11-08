"""Separate HT1080 profile variant with gene and genome-wide WGS copy number."""
from pathlib import Path
import argparse
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from make_cell_line_profiles import (
    THESIS, BLUE, apply_style, background_with_photo, chromosome_diagram,
    ordered_y_genes, export,
)

STEM = 'ch03_main_ht1080_profile_wgs_cn'


def draw_profile(source, output):
    meta = json.loads((source/'ht1080_background.json').read_text(encoding='utf-8'))
    cnmeta = json.loads((source/'wgs_cn_provenance.json').read_text(encoding='utf-8'))
    coords = ordered_y_genes(source)
    expr = pd.read_csv(source/'ht1080_y_expression.csv').set_index('gene').loc[coords.gene].reset_index()
    copy = pd.read_csv(source/'ht1080_y_gene_wgs_cn.csv').set_index('gene').loc[coords.gene].reset_index()
    deps = pd.read_csv(source/'ht1080_dependency_context.csv')
    assert len(expr)==len(copy)==9 and len(deps)==16
    fig = plt.figure(figsize=(183/25.4,200/25.4))
    background_with_photo(fig,source,meta)
    fig.text(.045,.585,'c',fontsize=11,weight='semibold',va='bottom')
    fig.text(.079,.585,'Y-linked expression and copy number',fontsize=9,weight='semibold',va='bottom')
    chromosome_diagram(fig,[.035,.16,.205,.385],source,coords)
    ax = fig.add_axes([.30,.16,.145,.385])
    yy = np.arange(len(expr))
    ax.barh(yy,expr.expression,color=BLUE,height=.60)
    ax.set_yticks(yy,expr.gene,fontsize=8)
    ax.set_ylim(len(expr)-.5,-.5)
    ax.set_xlim(0,11)
    ax.set_xticks([0,5,10])
    ax.set_xlabel('RNA expression\nlog₂(TPM + 1)',fontsize=7.3)
    for y,val in zip(yy,expr.expression):
        ax.text(val+.10,y,f'{val:.2f}',va='center',fontsize=7)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y',length=0)
    ax.text(.5,1.016,'RNA',transform=ax.transAxes,ha='center',va='bottom',fontsize=7.2,color='#555555')
    cnax = fig.add_axes([.505,.16,.07,.385])
    cnax.set_xlim(0,1)
    cnax.set_ylim(len(expr)-.5,-.5)
    cnax.set_axis_off()
    cnax.axvspan(0,1,color='#F3F3F3',zorder=0)
    for i,val in enumerate(copy.estimated_copy_number):
        cnax.text(.5,i,f'{val:.2f}',fontsize=7.5,ha='center',va='center',color='#333333')
    cnax.text(.5,1.016,'WGS CN\n(est. copies)',transform=cnax.transAxes,
              ha='center',va='bottom',fontsize=6.8,color='#555555')

    ax = fig.add_axes([.69,.16,.27,.385])
    fig.text(.595,.585,'d',fontsize=11,weight='semibold',va='bottom')
    fig.text(.629,.585,'Selected genetic dependencies',fontsize=9,weight='semibold',va='bottom')
    assert deps.n_reference.eq(1207).all()
    ax.text(0,1.008,'Other models: n = 1,207',transform=ax.transAxes,fontsize=7,color='#666666',va='bottom')
    for y,row in enumerate(deps.itertuples()):
        ax.plot([row.reference_q25,row.reference_q75],[y,y],color='#BDBDBD',lw=3,solid_capstyle='butt',zorder=1)
    yy = np.arange(len(deps))
    ax.scatter(deps.reference_median,yy,marker='|',s=100,color='#777777',label='Reference median',zorder=3)
    ax.scatter(deps.gene_effect,yy,s=27,color=BLUE,label='HT1080',zorder=4)
    ax.set_yticks(yy,deps.gene,fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Chronos gene effect')
    ax.axvline(0,lw=.5,color='#BDBDBD',zorder=0)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y',length=0)
    ax.legend(loc='upper left',bbox_to_anchor=(-.05,-.15),ncol=2,
              fontsize=6.7,handletextpad=.3,columnspacing=.8)

    # Keep the approved upper panels at their existing physical sizes and add
    # space underneath for chromosome-wide CN. Fonts are not scaled.
    old_height, new_height, added_mm = 200,250,50
    for ax in fig.axes:
        pos = ax.get_position()
        ax.set_position([pos.x0,(pos.y0*old_height+added_mm)/new_height,
                         pos.width,pos.height*old_height/new_height])
    for text in fig.texts:
        x,y = text.get_position()
        text.set_position((x,(y*old_height+added_mm)/new_height))
    fig.set_size_inches(183/25.4,new_height/25.4)

    fig.text(.045,.265,'DepMap Public 26Q1  |  Matched HT1080 model condition',fontsize=7.1,color='#666666')
    fig.text(.045,.246,'Background and image: repository properties; doubling time depends on culture conditions.',fontsize=7,color='#666666')
    fig.text(.045,.211,'e',fontsize=11,weight='semibold',va='bottom')
    fig.text(.079,.211,'Genome-wide WGS copy number',fontsize=9,weight='semibold',va='bottom')
    fig.text(.96,.212,f"WGS-matched ploidy: {cnmeta['matched_ploidy']:.2f}",fontsize=7.2,color='#555555',ha='right',va='bottom')
    ax = fig.add_axes([.12,.058,.84,.132])
    segments = pd.read_csv(source/'ht1080_genome_wgs_segments.csv')
    lengths = pd.read_csv(source/'chromosome_lengths_grch38.csv')
    gap = 15e6
    offset = 0
    centres, labels = [], []
    for i,row in enumerate(lengths.itertuples()):
        sub = segments.loc[segments.Chrom.eq(row.chromosome)]
        colour = BLUE if row.chromosome=='chrY' else ['#555555','#999999'][i%2]
        if i%2==0:
            ax.axvspan(offset,offset+row.length_bp,color='#F5F5F5',lw=0,zorder=0)
        ax.hlines(sub.estimated_copy_number,sub.Start+offset,sub.End+offset,
                  colors=colour,lw=1.05,zorder=3)
        # Midpoint dots make very short segments visible at this physical size.
        ax.scatter((sub.Start+sub.End)/2+offset,sub.estimated_copy_number,
                   s=3,color=colour,edgecolors='none',zorder=4)
        centres.append(offset+row.length_bp/2)
        labels.append(row.chromosome.replace('chr',''))
        offset += row.length_bp+gap
    ax.set_xlim(0,offset-gap)
    ax.set_ylim(-.2,np.ceil(segments.estimated_copy_number.max())+.25)
    ax.set_yticks([0,2,4,6,8])
    ax.set_xticks(centres,labels,fontsize=6.3)
    ax.tick_params(axis='x',length=0,pad=5)
    ax.tick_params(axis='y',labelsize=7)
    ax.axhline(2,color='#C0C0C0',linestyle='--',lw=.6,zorder=1)
    ax.set_ylabel('WGS CN\n(est. copies)',fontsize=7.4)
    ax.set_xlabel('Chromosome (GRCh38)',fontsize=7.4,labelpad=5)
    ax.spines['bottom'].set_visible(False)
    export(fig,output,STEM)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-dir',type=Path,default=THESIS/'figures/chapter3/source_data/ht1080_profile_wgs_cn')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    apply_style(base_size=8)
    draw_profile(args.source_dir,args.output_dir)


if __name__=='__main__':
    main()
