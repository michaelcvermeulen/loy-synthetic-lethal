"""Display each DDX3X guide's CPM in the seven retained A549 libraries."""
from pathlib import Path
import hashlib
import json
import os
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT.parent/'.figure_python_deps'))
OUT=ROOT/'figures/chapter4'
DATA=OUT/'source_data/a549_ddx3x_cpm'
DATA.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter
import pandas as pd
import numpy as np
import pymupdf
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style,condition_colours
apply_style(8)
COL=condition_colours()
SOURCE=ROOT.parent/'crispr_screen/a549_loy_screen1/CRISPR_screen/Screen1/Mageck/A549_full.count.txt'
COUNTS=OUT/'source_data/a549_initial_qc/retained_counts.tsv'
STEM='ch04_supp_a549_ddx3x_guide_cpm_draft'

def main():
    source_hash=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    raw=pd.read_csv(COUNTS,sep='\t')
    samples=raw.columns[2:].tolist()
    original=pd.read_csv(SOURCE,sep='\t')
    assert len(raw)==len(original)==77441
    assert raw.sgRNA.is_unique and original.sgRNA.is_unique
    pd.testing.assert_frame_equal(raw.set_index('sgRNA').sort_index(),
                                  original.set_index('sgRNA')[raw.columns[1:]].sort_index())
    d=raw.loc[raw.Gene.eq('DDX3X')].sort_values('sgRNA').copy()
    assert len(d)==4 and len(samples)==7 and not any('_P3_' in x for x in samples)
    d.to_csv(DATA/'raw_counts.csv',index=False)
    totals=raw[samples].sum()
    totals.rename_axis('sample').rename('total_guide_reads').to_csv(DATA/'library_totals.csv')
    d[samples]=d[samples].div(totals,axis=1)*1e6
    assert np.isfinite(d[samples]).all().all() and (d[samples]>=0).all().all()
    d.to_csv(DATA/'plotted_cpm.csv',index=False)
    zero_guides=d.loc[d[samples].sum(axis=1).eq(0),'sgRNA'].tolist()
    labels=['CTRL\nD0-1','CTRL\nD0-2','LOY\nD0-1','LOY\nD0-2','CTRL\nP4-R2','LOY\nP4-R1','LOY\nP4-R2']
    colours=[COL['Y-retaining'] if s.startswith('CTRL') else COL['Y-depleted'] for s in samples]
    fig,axes=plt.subplots(2,2,figsize=(183/25.4,131/25.4),sharey=True)
    fig.subplots_adjust(left=.10,right=.985,bottom=.16,top=.91,wspace=.16,hspace=.42)
    for ax,(_,row),letter in zip(axes.flat,d.iterrows(),'abcd'):
        ax.axvline(3.5,color='#B7BDC5',lw=.7,ls='--',zorder=0)
        for i,(sample,colour) in enumerate(zip(samples,colours)):
            value=row[sample]
            ax.scatter(i,value,s=30,marker='D' if '_D0_' in sample else 'o',color=colour,edgecolor='white',lw=.55,zorder=3,clip_on=False)
            ax.annotate(f'{value:.3g}',(i,value),xytext=(0,6),textcoords='offset points',ha='center',va='bottom',fontsize=6.7,color=colour)
        ax.set_yscale('linear')
        ax.set(xlim=(-.55,6.55),ylim=(0,45))
        ax.yaxis.set_major_locator(FixedLocator([0,10,20,30,40]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value,_:f'{value:g}'))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis='y',which='minor',length=0)
        ax.set_xticks(range(7),labels,fontsize=6.7)
        ax.tick_params(axis='x',length=0,pad=4)
        for label,colour in zip(ax.get_xticklabels(),colours):label.set_color(colour)
        if row.sgRNA in zero_guides:
            ax.text(.5,.45,'No reads detected in any retained library',transform=ax.transAxes,ha='center',fontsize=7.5,color='#666666')
        ax.set_title('DDX3X guide '+row.sgRNA.split('_')[1],fontsize=9,pad=8)
        ax.text(-.16,1.08,letter,transform=ax.transAxes,fontweight='semibold',fontsize=11)
        ax.grid(axis='y',which='major',color='#E9EBEE',lw=.4,zorder=0)
    fig.supylabel('Guide abundance (CPM)',x=.012,fontsize=9)
    fig.legend(handles=[Line2D([],[],marker='D',ls='',color='#666666',label='D0: technical libraries'),
                        Line2D([],[],marker='o',ls='',color='#666666',label='P4: outgrowth cultures')],
               loc='lower center',bbox_to_anchor=(.53,.008),ncol=2,fontsize=7.5,handletextpad=.4,columnspacing=2)
    fig.text(.53,.077,'Individual libraries shown; no pseudocount added.',ha='center',fontsize=7,color='#666666')
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600)
    preview=OUT/'review'/f'{STEM}_preview.png';fig.savefig(preview,dpi=180)
    preview_hashed=preview.with_name(STEM+'_linear_'+hashlib.sha256(preview.read_bytes()).hexdigest()[:12]+'.png')
    preview_hashed.write_bytes(preview.read_bytes())
    plt.close(fig)
    with pymupdf.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all('OpenSans' in f[3].replace(' ','') for f in page.get_fonts(full=True))
        outside=[w for w in page.get_text('words') if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
        assert not outside,outside
        assert 'P3' not in page.get_text()
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==source_hash
    result={'source':str(SOURCE),'source_sha256':source_hash,'count_matrix':str(COUNTS),
            'counts_sha256':hashlib.sha256(COUNTS.read_bytes()).hexdigest(),
            'zero_read_guides':zero_guides,'retained_matrix_matches_original':True,'n_guides':4,'n_libraries':7,'n_points':28,'y_axis':'linear CPM; shared limits 0–45; no pseudocount',
            'normalisation':'Raw guide count / total mapped guide counts within each library * 1e6',
            'data_match_counts':True,'fonts_open_sans':True,'text_within_page':True,'source_unchanged':True,
            'figure_mm':[183,131],'png_dpi':600,'preview':str(preview_hashed)}
    (DATA/'manifest.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
