"""Matching rank views of saved HT1080 comparisons; no screen refitting."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import make_ht1080_candidates as C

ISO = C.OUT / 'source_data/ht1080_isogenicz_20261020'
OUT = C.OUT / 'source_data/ht1080_method_panels_20261022'
SOURCES = [C.DATA / 'method_comparison.csv', C.DATA / 'summary.json',
           ISO / 'mask173_detected_IsogenicZ_genes.tsv',
           ISO / 'mask173_detected_matched_endpoint_drugZ_genes.tsv',
           ISO / 'method_agreement.csv', ISO / 'figure_table_values.csv']


def hashes():
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES}


def rank_panel(ax, frame, ylabel, xlabel, labels, rho, n):
    ax.scatter(frame.x, frame.y, s=2.2, color='#CBCBCB', alpha=.55,
               linewidths=0, rasterized=True, zorder=1)
    hits = frame.loc[frame.primary_drugZ_FDR.lt(.25)]
    ax.scatter(hits.x, hits.y, s=15, c=[C.colour(v) for v in hits.primary_drugZ_FDR],
               edgecolor='white', linewidths=.2, zorder=3)
    for gene, xytext in labels.items():
        row = frame.loc[gene]
        if row.primary_drugZ_FDR >= .25:
            ax.scatter([row.x], [row.y], s=17, facecolor='white', edgecolor='#444444',
                       linewidths=.6, zorder=4)
        ax.annotate(gene, (row.x, row.y), xytext=xytext, fontsize=6.5,
                    color='#444444', arrowprops={'arrowstyle':'-', 'color':'#999999', 'lw':.45})
    ax.plot([1,25000], [1,25000], color='#888888', lw=.6, ls='--', zorder=0)
    ax.set(xscale='log', yscale='log', xlim=(.8,27000), ylim=(.7,27000))
    ticks = [1,10,100,1000,10000]
    ticklabels = ['1','10','100','1,000','10,000']
    ax.set_xticks(ticks,ticklabels); ax.set_yticks(ticks,ticklabels)
    ax.xaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_locator(NullLocator())
    ax.tick_params(labelsize=7,pad=2)
    ax.set_xlabel(xlabel,fontsize=7.5,labelpad=2)
    ax.set_ylabel(ylabel,fontsize=7.5,labelpad=3)
    ax.text(.52,.035,f'Spearman ρ = {rho:.3f}\n{n:,} genes',
            transform=ax.transAxes,fontsize=6.8)
    return {'genes':len(frame), 'highlighted_primary_candidates':len(hits),
            'labels':list(labels), 'label_positions':labels,
            'rho_from_saved_analysis':float(rho), 'xlim':[.8,27000], 'ylim':[.7,27000]}


def draw_comparisons(fig, height, *, axes_bounds=None, notes_positions=None):
    before = hashes()
    comp = pd.read_csv(C.DATA/'method_comparison.csv').set_index('GENE')
    iso = pd.read_csv(ISO/'mask173_detected_IsogenicZ_genes.tsv',sep='\t').set_index('Gene')
    endpoint = pd.read_csv(ISO/'mask173_detected_matched_endpoint_drugZ_genes.tsv',sep='\t').set_index('Gene')
    agreement = pd.read_csv(ISO/'method_agreement.csv').set_index('comparison').loc[
        'mask173_detected: IsogenicZ versus matched endpoint drugZ']
    a = pd.DataFrame({'x':comp.rank_synth,'y':comp['neg|rank'], 'primary_drugZ_FDR':comp.fdr_synth})
    b = pd.DataFrame({'x':endpoint.loc[iso.index,'rank_synth'],'y':iso.rank_synth,
                      'primary_drugZ_FDR':comp.reindex(iso.index).fdr_synth})
    assert len(a)==19112 and len(b)==19082 and b.notna().all().all()
    labels_a={'DDX3X':(1.8,1.1),'TIMM8A':(1.0,7),
              'SREBF2':(1.0,600),'SREBF1':(6,1800),'EIF4A3':(60,5500),
              'TADA2B':(150,28),'TAF6L':(70,7),'SUPT20H':(18,1.1)}
    labels_b={'DDX3X':(1.8,1.0),'TIMM8A':(1.0,7),
              'SREBF2':(1.0,40),'SREBF1':(1.0,250),
              'TAF6L':(30,1.5),'TADA2B':(100,5),'SUPT20H':(110,22),
              'EIF4A3':(80,100),'EIF1AX':(100,14000),'L3MBTL2':(1800,23000)}
    axes_bounds=axes_bounds or [[.064,8.65/height,.255,2.28/height],[.389,8.65/height,.255,2.28/height]]
    notes_positions=notes_positions or [[.064,8.12/height],[.389,8.12/height]]
    axa=fig.add_axes(axes_bounds[0])
    axb=fig.add_axes(axes_bounds[1])
    meta_a=rank_panel(axa,a,'MAGeCK depletion rank','Primary drugZ depletion rank',labels_a,
                      C.SUMMARY['method_rank_spearman'],len(a))
    meta_b=rank_panel(axb,b,'IsogenicZ depletion rank','Primary drugZ depletion rank',labels_b,
                      agreement.spearman,len(b))
    fig.text(*notes_positions[0],'MAGeCK support for 48 primary drugZ hits:\n21 at FDR <0.25; 9 at FDR <0.05',
             fontsize=7,va='top',linespacing=1.25)
    fig.text(*notes_positions[1],f'{int(agreement.shared_top100)} shared top-100 ranks\nSame 76,033 guides in both fits',
             fontsize=7,va='top',linespacing=1.25)
    OUT.mkdir(parents=True,exist_ok=True)
    a.to_csv(OUT/'panel_a_primary_drugz_mageck_ranks.csv',index_label='Gene')
    b.to_csv(OUT/'panel_b_matched_drugz_isogenicz_ranks.csv',index_label='Gene')
    assert hashes()==before
    meta={'a':meta_a,'b':meta_b,'sources':before,'source_files_unchanged':True,
          'same_axes':True,'colour_definition':'Historical primary drugZ depletion FDR in both panels',
          'comparison_a':'Historical primary drugZ versus MAGeCK; 19,112 shared genes',
          'comparison_b':'mask173_detected matched endpoint drugZ versus IsogenicZ; 19,082 genes',
          'display_change':'Panel b now shows saved depletion ranks, replacing the prior normZ hexbin view; no scores or FDRs refitted.'}
    (OUT/'comparison_manifest.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    return meta


def render_coverage_table():
    """Keep former Supp 4.6b as the sole content of Supp 4.6."""
    values=pd.read_csv(ISO/'figure_table_values.csv',index_col=0)
    def q(v):
        if pd.isna(v): return '—'
        return f'{v:.3f}' if v>=.001 else f'{v:.2g}'
    rows=[]
    for gene,r in values.iterrows():
        rows.append([gene,str(int(r.review_retained_guides)),f'{r.review_normZ:.2f}',q(r.review_fdr_synth),
                     '0' if pd.isna(r.D0covered_normZ) else str(int(r.D0covered_retained_guides)),q(r.D0covered_fdr_synth)])
    fig=plt.figure(figsize=(5.4,3.65))
    ax=fig.add_axes([.025,.16,.95,.82]);ax.axis('off')
    tab=ax.table(cellText=rows,colLabels=['Gene','Guides','normZ','FDR','≥50 n','≥50 FDR'],
                 colWidths=[.23,.105,.14,.175,.12,.23],bbox=[0,0,1,1],cellLoc='center')
    tab.auto_set_font_size(False);tab.set_fontsize(8)
    for (i,j),cell in tab.get_celld().items():
        cell.set_edgecolor('white');cell.set_facecolor('#EAEAEA' if i==0 else '#F5F5F5' if i%2 else 'white')
        if j==0:cell.set_text_props(ha='left')
    fig.text(.025,.105,'FDR: depletion, BH within each fit. ≥50: reads in both D0 libraries.',fontsize=7.3,color='#555555')
    fig.text(.025,.045,'—: not tested; all guides fail the coverage filter.',fontsize=7.3,color='#555555')
    C.export(fig,'ch04_supp_ht1080_isogenicz_review')
    plt.close(fig)
    values.to_csv(OUT/'supp_4_6_coverage_table_values.csv')
    return {'inches':[5.4,3.65], 'rows':len(rows), 'values_source':str(ISO/'figure_table_values.csv'),
            'table_values_unchanged':True,'removed_panel':'Former a moved to Supp Fig 4.3b as matching rank plot'}
