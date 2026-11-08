"""Combine frozen hit/SAGA co-dependencies and DepMap–STRING evidence."""
from pathlib import Path
import json,hashlib
import make_ht1080_depmap_codependency as D
import make_ht1080_depmap_string_network as N
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,Rectangle
from matplotlib.colors import LinearSegmentedColormap,Normalize
import pymupdf as fitz

OUT=D.OUT
DATA=OUT/'source_data/ht1080_combined_dependency_context'
DATA.mkdir(exist_ok=True)
STEM='ch04_supp_ht1080_combined_dependency_context_draft'
W,H=12.4,15.0
EXTRA_PARTNER_SEEDS=['TADA2B','CREBBP','TADA1','SUPT20H','ATXN7L3','KIF2A']
PARTNER_BLOCK_SIZE=8


def rect(x,y,w,h):return [x/W,y/H,w/W,h/H]


def top_partner_panel(fig):
    """Rank heatmaps: each seed column has its own ten genome-wide partners."""
    ranked=N.NODES.assign(abs_normZ=N.NODES.normZ.abs()).sort_values('abs_normZ',ascending=False,kind='stable')
    top_ten=ranked.head(10).index.tolist()
    seeds=top_ten+[g for g in EXTRA_PARTNER_SEEDS if g not in top_ten]
    screen_rank={g:i+1 for i,g in enumerate(ranked.index)}
    assert len(seeds)==16
    assert {'DDX3X','SREBF1','SREBF2'}.issubset(seeds)
    all_scores=pd.read_csv(DATA.parent/'ht1080_candidate_discovery/all_gene_scores.csv').set_index('GENE')
    partner_source=pd.concat([D.TOP,pd.read_csv(DATA/'additional_top_partners.csv').fillna('')],ignore_index=True)
    fig.text(7.4/W,7.48/H,'d',fontsize=13,fontweight='semibold',va='top')
    fig.text(7.85/W,7.39/H,'Ranked co-dependency partners',fontsize=8.5)
    fig.text(7.85/W,7.16/H,'16 selected genes · 10 independently ranked partners each',fontsize=7.1,color='#555555')
    tables=[]
    axes=[]
    for i,gene in enumerate(seeds):
        t=partner_source.loc[partner_source.seed.eq(gene)&partner_source.selection.eq('Top positive')].sort_values('rank_positive').head(10).copy()
        assert len(t)==10 and t.rank_positive.tolist()==list(range(1,11))
        assert (t.partner!=gene).all() and (t.n_models>=200).all() and (t.pearson_r>0).all()
        t.insert(0,'seed_screen_abs_normZ_rank',screen_rank.get(gene,np.nan))
        t.insert(0,'seed_primary_candidate',gene in N.NODES.index)
        t.insert(0,'seed_display_order',i+1)
        t.insert(1,'seed_selection','Top ten by absolute drugZ normZ' if gene in top_ten else 'Additional author-selected gene')
        t.insert(2,'seed_normZ',float(all_scores.loc[gene,'normZ']))
        tables.append(t)
    # Retain the positive red half of the Pearson palette used in a and b.
    cmap=LinearSegmentedColormap.from_list('positive_pearson',D.CMAP(np.linspace(.5,1,256)))
    norm=Normalize(0,1)
    for block,y0 in enumerate([4.05,1.10]):
        ax=fig.add_axes(rect(7.65,y0,4.55,2.65))
        block_tables=tables[block*PARTNER_BLOCK_SIZE:(block+1)*PARTNER_BLOCK_SIZE]
        values=np.column_stack([t.pearson_r.to_numpy() for t in block_tables])
        im=ax.imshow(values,cmap=cmap,norm=norm,aspect='auto',interpolation='none')
        ax.set_xticks(range(PARTNER_BLOCK_SIZE),seeds[block*PARTNER_BLOCK_SIZE:(block+1)*PARTNER_BLOCK_SIZE],fontsize=7.4,fontweight='semibold')
        ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=5)
        ax.set_yticks(range(10),range(1,11),fontsize=7.2)
        ax.tick_params(axis='y',length=0,pad=3)
        ax.set_ylabel('Partner rank',fontsize=7.5,labelpad=5)
        ax.set_xticks(np.arange(-.5,PARTNER_BLOCK_SIZE,1),minor=True)
        ax.set_yticks(np.arange(-.5,10,1),minor=True)
        ax.grid(which='minor',color='white',linewidth=.7)
        ax.tick_params(which='minor',length=0)
        for spine in ax.spines.values():spine.set_visible(False)
        for col,t in enumerate(block_tables):
            for row,p in enumerate(t.itertuples()):
                rgb=np.array(cmap(norm(p.pearson_r))[:3])
                linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
                luminance=float(linear@np.array([.2126,.7152,.0722]))
                colour='#222222' if luminance>.179 else 'white'
                is_saga=p.partner in D.A.index and bool(D.A.loc[p.partner,'saga_reference_20'])
                hit=p.partner_screen_direction!='Not nominated'
                ax.text(col,row-.19,p.partner+('*' if hit else ''),ha='center',va='center',
                        fontsize=6.9,fontweight='semibold' if is_saga else 'normal',color=colour)
                ax.text(col,row+.29,f'{p.pearson_r:.2f}',ha='center',va='center',fontsize=6.0,color=colour)
        axes.append(ax)
    cb=fig.colorbar(im,cax=fig.add_axes(rect(8.55,.65,2.75,.10)),orientation='horizontal',ticks=[0,.5,1])
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=7,length=2,pad=2)
    cb.set_label('Pearson r · unadjusted',fontsize=7.8,labelpad=2)
    fig.text(7.7/W,.12/H,'* Depleted screen hit; semibold: SAGA reference member.',fontsize=7.1,color='#555555')
    table=pd.concat(tables,ignore_index=True)
    table.to_csv(DATA/'panel_d_top_partners.csv',index=False)
    return table,axes


def main():
    D.apply_style(8)
    order=[g for g in pd.read_csv(D.DATA/'hit_heatmap_order.csv').gene if g!='ZCCHC16']
    leaves=pd.read_csv(D.DATA/'hit_heatmap_leaf_ids.csv').gene.tolist()
    linkage=pd.read_csv(D.DATA/'hit_heatmap_linkage.csv').to_numpy().tolist()
    saga=pd.read_csv(D.DATA/'saga_heatmap_order.csv').gene.tolist()
    groups,module,core,strip,evidence_order,pos=N.load_saved_structure()
    evidence_order=[g for g in evidence_order if g!='ZCCHC16']
    assert len(order)==56 and len(leaves)==56 and len(saga)==24 and len(evidence_order)==56
    fig=plt.figure(figsize=(W,H))
    fig.text(.16/W,14.90/H,'a',fontsize=13,fontweight='semibold',va='top')
    fig.text(.75/W,14.82/H,'Screen hits · unadjusted · 56 genes · 1,186 models',fontsize=8.5)
    ax=fig.add_axes(rect(.75,8.05,6.0,6.0))
    im=D.matrix(ax,D.R,order,size=7.1)
    ax.yaxis.tick_right();ax.tick_params(axis='y',pad=9,labelright=True,labelleft=False)
    top=fig.add_axes(rect(.75,14.12,6.0,.53))
    left=fig.add_axes(rect(.08,8.05,.59,6.0))
    D.draw_hit_dendrograms(top,left,leaves,order,linkage)
    for label in ax.get_xticklabels()+ax.get_yticklabels():
        g=label.get_text();label.set_color(D.DIR[D.A.loc[g,'screen_direction']])
        if D.A.loc[g,'screen_fdr']<.05:label.set_fontweight('semibold')
        if not bool(D.A.loc[g,'depmap_available']):label.set_color('#777777')
    for i,g in enumerate(order):
        if bool(D.A.loc[g,'saga_reference_20']):
            ax.add_patch(Rectangle((len(order)-.25,i-.37),.22,.74,facecolor='#BAA04B',edgecolor='none',clip_on=False))

    fig.text(7.85/W,14.90/H,'b',fontsize=13,fontweight='semibold',va='top')
    fig.text(8.5/W,14.82/H,'SAGA · lineage/sex-adjusted · 1,186 models',fontsize=8.5)
    bx=fig.add_axes(rect(8.5,11.15,3.5,3.5))
    D.matrix(bx,D.ADJ,saga,size=7.2,annotate=True)
    for label in bx.get_xticklabels()+bx.get_yticklabels():
        g=label.get_text();label.set_color(D.MOD[D.A.loc[g,'saga_module']])
        if bool(D.A.loc[g,'screen_hit']):label.set_fontweight('semibold')
    for i,g in enumerate(saga[1:],1):
        if D.A.loc[g,'saga_module']!=D.A.loc[saga[i-1],'saga_module']:
            bx.axhline(i-.5,color='white',lw=1);bx.axvline(i-.5,color='white',lw=1)
    fig.legend(handles=[Patch(facecolor=c,label=m if m!='Alternatives' else 'Alternatives') for m,c in D.MOD.items()],
               loc='upper left',bbox_to_anchor=(8.0/W,10.42/H),ncol=3,fontsize=7.5,frameon=False,columnspacing=1)
    fig.text(8.05/W,9.85/H,'b: semibold, six screen hits; numbers, |r| ≥0.40.',fontsize=7.5,color='#555555')
    fig.legend(handles=[Patch(facecolor=D.DIR['Depleted'],label='Depleted in LOY (47)'),
                        Patch(facecolor=D.DIR['Enriched'],label='Enriched in LOY (9)'),
                        Patch(facecolor='#BAA04B',label='SAGA reference member')],
               loc='upper left',bbox_to_anchor=(7.95/W,9.48/H),ncol=2,fontsize=7.5,frameon=False,columnspacing=1)
    fig.text(8.05/W,8.82/H,'a: label colour, screen direction; semibold, FDR <0.05.',fontsize=7.4,color='#555555')
    cb=fig.colorbar(im,cax=fig.add_axes(rect(8.5,8.37,3.2,.11)),orientation='horizontal',ticks=[-1,0,1])
    cb.set_label('Pearson r (a–b)',fontsize=8,labelpad=2)
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=7.5,length=2)

    fig.text(.16/W,7.48/H,'c',fontsize=13,fontweight='semibold',va='top')
    # Keep square cells, with partner profiles in the lower-right area.
    cx=N.matrix_panel(fig,evidence_order,groups,rect=rect(.95,1.0,6.0,6.0),bar_rect=rect(.95,7.12,6.0,.13),
                      output_path=DATA/'panel_c_evidence_codes.csv')
    cx.set_aspect('equal')
    for label in cx.texts:
        if label.get_text().startswith('Order: positive DepMap modules'):
            label.set_text('Order: positive DepMap modules,\nthen remaining hits')
            label.set_va('top')
    fig.text(.95/W,.17/H,'Each pair is shown once; the upper triangle and diagonal are omitted.',fontsize=8,color='#555555')
    partners,partner_axes=top_partner_panel(fig)
    fig.canvas.draw()
    matrix_inches=cx.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    assert np.isclose(matrix_inches.width,matrix_inches.height)
    assert np.allclose(ax.images[0].get_array(),D.R.reindex(index=order,columns=order).to_numpy(),equal_nan=True)
    assert np.allclose(bx.images[0].get_array(),D.ADJ.loc[saga,saga].to_numpy(),equal_nan=True)
    old_codes=pd.read_csv(N.DATA/'panel_c_evidence_codes.csv',index_col=0)
    assert np.allclose(cx.images[0].get_array(),old_codes.loc[evidence_order,evidence_order].to_numpy(),equal_nan=True)
    for block,dx in enumerate(partner_axes):
        expected=np.column_stack([partners.loc[partners.seed_display_order.eq(rank),'pearson_r'].to_numpy()
                                  for rank in range(block*PARTNER_BLOCK_SIZE+1,(block+1)*PARTNER_BLOCK_SIZE+1)])
        assert np.allclose(dx.images[0].get_array(),expected)
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{STEM}.{ext}',dpi=600 if ext=='png' else 300)
    fig.savefig(OUT/'review'/f'{STEM}_preview.png',dpi=140)
    plt.close(fig)
    with fitz.open(OUT/f'{STEM}.pdf') as doc:
        page=doc[0]
        assert all(g in page.get_text() for g in set(order+saga))
        assert all(g in page.get_text() for g in partners.partner)
        assert 'ZCCHC16' not in page.get_text()
        assert all('OpenSans' in f[3] for f in page.get_fonts())
        spans=[s for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for s in line['spans'] if s['text'].strip()]
        outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
        page.get_pixmap(matrix=fitz.Matrix(1.2,1.2),colorspace=fitz.csGRAY).save(OUT/'review'/f'{STEM}_greyscale.png')
    source_paths=[D.DATA/f for f in ['seed_pairwise_pearson.csv','seed_pairwise_lineage_sex_residual_pearson.csv',
                  'hit_heatmap_order.csv','hit_heatmap_leaf_ids.csv','hit_heatmap_linkage.csv','saga_heatmap_order.csv','top_partners.csv']]
    source_paths += [N.DATA/f for f in ['pairwise_evidence.csv','node_layout_and_modules.csv','panel_c_evidence_codes.csv']]
    source_paths += [DATA/'additional_top_partners.csv',DATA.parent/'ht1080_candidate_discovery/all_gene_scores.csv']
    record=dict(stem=STEM,page_inches=[W,H],png_dpi=600,panels=['a: all-hit primary Pearson correlations',
        'b: adjusted SAGA correlations (original panel b)','c: complete DepMap/STRING evidence',
        'd: ranked heatmaps of ten genome-wide positive partners for sixteen selected genes'],
        primary_candidates=57,hit_genes=56,hit_observed_genes=56,hit_order=order,saga_genes=24,saga_order=saga,evidence_order=evidence_order,
        display_omissions={'ZCCHC16':'No DepMap profile; omitted from this supplement at author request'},
        unique_evidence_pairs=1540,correlation_scale=[-1,1],source_values_exactly_preserved=True,
        panel_c_size_inches=[matrix_inches.width,matrix_inches.height],
        panel_d_seeds=partners.seed.drop_duplicates().tolist(),panel_d_partners_per_seed=10,
        panel_d_selection='Top ten candidates by absolute primary drugZ normZ, followed by author-selected TADA2B, CREBBP, TADA1, SUPT20H, ATXN7L3 and KIF2A; ten highest positive genome-wide Pearson correlations per seed, self-correlations excluded. KIF2A (a primary candidate) replaced KIF1A at author request on 22 October 2026 and uses its saved genome-wide ranking.',
        panel_d_additional_selected_genes=EXTRA_PARTNER_SEEDS,
        panel_d_layout='Two 10-rank by 8-seed heatmaps; rows are partner ranks, not shared genes; names and r printed in each cell.',
        panel_d_model_count_range=[int(partners.n_models.min()),int(partners.n_models.max())],
        panel_d_pearson_scale=[0,1],panel_d_correlation_adjustment='Unadjusted',
        existing_analysis_refit=False,kif1a_profile_newly_computed=False,text_outside_page=outside,
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (DATA/'render_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k not in ['source_sha256','hit_order','saga_order','evidence_order']},indent=2))
    assert not outside,outside


if __name__=='__main__':main()
