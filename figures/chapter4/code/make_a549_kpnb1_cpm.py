"""Plot raw-library CPM for every KPNB1 guide and all 11 A549 libraries."""
from pathlib import Path
import hashlib
import json
import os
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT.parent/'.figure_python_deps'))
OUT=ROOT/'figures/chapter4'
DATA=OUT/'source_data/a549_kpnb1_cpm_20261020'
DATA.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(DATA/'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pymupdf
sys.path.insert(0,str(ROOT/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style, condition_colours

SOURCE=ROOT.parent/'crispr_screen/a549_loy_screen1/CRISPR_screen/Screen1/Mageck/A549_full.count.txt'
ELIGIBILITY=ROOT/'research/a549_mle_timecourse_20261020/inputs/guide_eligibility.csv'
STEM='ch04_supp_a549_kpnb1_guide_cpm_review'
SAMPLES=['CTRL_D0_1','CTRL_D0_2','LOY_D0_1','LOY_D0_2',
         'CTRL_P3_R2','LOY_P3_R1','LOY_P3_R2','CTRL_P3_R1',
         'CTRL_P4_R2','LOY_P4_R1','LOY_P4_R2']
LABELS=['CTRL\nD0-1','CTRL\nD0-2','LOY\nD0-1','LOY\nD0-2',
        'CTRL\nP3-R2','LOY\nP3-R1','LOY\nP3-R2','LOY*\nP3-R1',
        'CTRL\nP4-R2','LOY\nP4-R1','LOY\nP4-R2']
ORDER=['10271','10272','10269','10270']


def main():
    raw_hash=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    raw=pd.read_csv(SOURCE,sep='\t').set_index('sgRNA')
    assert set(raw.columns)==set(['Gene']+SAMPLES)
    counts=raw.loc[raw.Gene.eq('KPNB1'),SAMPLES].copy()
    guides={sg.split('_')[1]:sg for sg in counts.index}
    assert set(guides)==set(ORDER)
    counts=counts.loc[[guides[g] for g in ORDER]]
    totals=raw[SAMPLES].sum()
    cpm=counts.div(totals,axis=1)*1e6
    assert np.isfinite(cpm.to_numpy()).all()
    np.testing.assert_allclose(raw[SAMPLES].div(totals).sum(),1,rtol=1e-12)
    np.testing.assert_allclose(cpm.mul(totals,axis=1)/1e6,counts,atol=1e-10)
    eligibility=pd.read_csv(ELIGIBILITY).set_index('sgRNA').loc[counts.index]
    assert eligibility.retained_for_MLE.tolist()==[True,True,False,False]
    counts.to_csv(DATA/'KPNB1_raw_counts_all_libraries.csv')
    cpm.to_csv(DATA/'KPNB1_CPM_all_libraries.csv')
    totals.rename('total_mapped_guide_reads').to_csv(DATA/'library_denominators.csv',index_label='sample')
    eligibility.to_csv(DATA/'KPNB1_MLE_eligibility.csv')
    rows=[]
    for sg in counts.index:
        for sample in SAMPLES:
            rows.append({'sgRNA':sg,'guide_id':sg.split('_')[1],'sample':sample,
                'condition':'LOY' if sample.startswith('LOY') or sample=='CTRL_P3_R1' else 'CTRL',
                'timepoint':sample.split('_')[1],'reassigned_as_LOY':sample=='CTRL_P3_R1',
                'raw_count':int(counts.loc[sg,sample]),'total_mapped_guide_reads':int(totals[sample]),
                'CPM':float(cpm.loc[sg,sample]),
                'retained_for_MLE':bool(eligibility.loc[sg,'retained_for_MLE'])})
    pd.DataFrame(rows).to_csv(DATA/'plotted_values.csv',index=False)

    apply_style(8)
    colours=condition_colours()
    cs=[colours['Y-depleted'] if s.startswith('LOY') or s=='CTRL_P3_R1'
        else colours['Y-retaining'] for s in SAMPLES]
    x=np.array([0,1,2,3,4.4,5.4,6.4,7.4,8.8,9.8,10.8])
    fig,axes=plt.subplots(2,2,figsize=(9.6,6.3))
    fig.subplots_adjust(left=.075,right=.99,bottom=.17,top=.94,wspace=.17,hspace=.39)
    for ax,sg,letter in zip(axes.flat,counts.index,'abcd'):
        maximum=float(cpm.loc[sg].max())
        ax.set(xlim=(-.55,11.35),ylim=(0,maximum*1.30),ylabel='CPM')
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=3))
        ax.axvline(3.7,color='#C9CDD2',lw=.65,ls='--')
        ax.axvline(8.1,color='#C9CDD2',lw=.65,ls='--')
        ax.grid(axis='y',color='#ECEEF0',lw=.4)
        for xi,sample,colour in zip(x,SAMPLES,cs):
            value=float(cpm.loc[sg,sample])
            marker='D' if '_D0_' in sample else 's' if sample=='CTRL_P3_R1' else 'o'
            ax.scatter(xi,value,s=28,marker=marker,color=colour,edgecolors='white',
                       linewidths=.45,zorder=4,clip_on=False)
            label='0' if value==0 else f'{value:.2g}'
            ax.annotate(label,(xi,value),xytext=(0,6),textcoords='offset points',
                        ha='center',fontsize=6.6,color=colour)
        ax.set_xticks(x,LABELS,fontsize=6.7)
        ax.tick_params(axis='x',length=0,pad=4)
        for tick,colour in zip(ax.get_xticklabels(),cs):tick.set_color(colour)
        guide=sg.split('_')[1]
        status='retained in MLE' if eligibility.loc[sg,'retained_for_MLE'] else 'excluded: low D0 counts'
        ax.set_title(f'KPNB1 guide {guide}  |  {status}',fontsize=8.8,pad=8)
        ax.text(-.15,1.055,letter,transform=ax.transAxes,fontsize=12,fontweight='semibold')
    fig.text(.075,.092,'*Recorded CTRL_P3_R1; reassigned LOY (square). D0 diamonds are technical libraries.',fontsize=7.5,color='#555555')
    fig.text(.075,.057,'CPM = raw count / all mapped guide reads × 10⁶; no pseudocount. Linear, panel-specific y scales.',fontsize=7.5,color='#555555')
    fig.text(.075,.022,'Unfiltered library totals are used here; P3 CPM therefore retains the effect of PCR-dominated libraries.',fontsize=7.5,color='#555555')
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{STEM}.{ext}',dpi=600)
    preview=OUT/'review'/f'{STEM}_preview.png'
    fig.savefig(preview,dpi=160)
    hashed=preview.with_name(f'{STEM}_{hashlib.sha256(preview.read_bytes()).hexdigest()[:12]}.png')
    hashed.write_bytes(preview.read_bytes())
    plt.close(fig)
    with pymupdf.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all('OpenSans' in f[3].replace(' ','') for f in page.get_fonts(full=True))
        outside=[w for w in page.get_text('words') if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
        assert not outside,outside
        assert all(g in page.get_text() for g in ORDER)
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==raw_hash
    result={'source':str(SOURCE),'source_sha256':raw_hash,'source_unchanged':True,
        'normalization':'Raw count / total mapped guide counts over all 77,441 guides within the same unfiltered library * 1e6',
        'not_MLE_NTC_normalized':True,'no_pseudocount':True,'n_guides':4,'n_samples':11,'n_points':44,
        'reassigned_sample':'CTRL_P3_R1 shown as LOY; original identifier retained in data',
        'y_axis':'linear, separate limits per guide','font':'Open Sans','pdf_text_within_page':True,
        'retained_guides':['10271','10272'],'D0_coverage_excluded_guides':['10269','10270'],
        'figure_inches':[9.6,6.3],'preview':str(hashed),'matplotlib_version':matplotlib.__version__,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (DATA/'manifest.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    (DATA/'README.md').write_text('''# A549 KPNB1 guide CPM, all libraries

Four library guides across all 11 original A549 libraries are shown on linear,
panel-specific y scales. CPM uses the total counts across every guide in the
original unfiltered sample; it is not the NTC-scaled quantity used by MLE.
No pseudocount is added. PCR-dominated P3 libraries retain their original
denominators, so their CPM is affected by the dominance of other guides.

Recorded CTRL_P3_R1 is displayed as LOY with a square and an asterisk; its
source identifier is preserved in the tables. D0 technical libraries are
shown separately. Guides 10271 and 10272 passed the MLE D0 coverage rule;
10269 and 10270 did not. None of these four guides is in the 342-guide mask.

The two retained guides are discordant at P4: 10272 has 1,128 reads in CTRL
versus 3 and 1 in LOY; 10271 has 8 in CTRL versus 25 and 80 in LOY. The strong
gene-level MLE score therefore warrants guide-level review before prioritizing
KPNB1 as a replicated LOY dependency.
''',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
