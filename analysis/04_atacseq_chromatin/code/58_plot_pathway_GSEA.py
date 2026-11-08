"""Render GSEA comparison figures, concordance/leading-edge tables, and HTML report."""
from pathlib import Path
from html import escape
import hashlib
import json
import math
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from adjustText import adjust_text

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'integration/pathway_GSEA'
FIG=OUT/'figures'
MAIN=['ATAC','RNA_Combined','RNA_AUG_2025','RNA_MAY_2026']
LABELS={'ATAC':'Promoter\nATAC','RNA_Combined':'RNA\nCombined','RNA_AUG_2025':'RNA\nAugust 2025','RNA_MAY_2026':'RNA\nMay 2026',
        'ATAC_representative':'ATAC\nAbundant peak','ATAC_median_logFC':'ATAC\nMedian log₂FC','ATAC_passage':'ATAC\nPassage model','ATAC_no_chr19':'ATAC\nWithout chr19',
        'ATAC_tie_shuffle_1':'ATAC\nTie order 1','ATAC_tie_shuffle_2':'ATAC\nTie order 2'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
CMAP=LinearSegmentedColormap.from_list('LOY_enrichment',['#2166ac','#f7f7f7','#b2182b'])
CMAP.set_bad('#dedede')


def pretty(name):
    words=name.replace('HALLMARK_','').replace('REACTOME_','').split('_')
    acronyms={'KRAS','MYC','E2F','G2M','G2','M','DNA','RNA','ECM','IL6','IL2','JAK','STAT3','STAT5','TNFA','NFKB','UV','PI3K','AKT','MTOR','MTORC1','TGF','P53','IGF','IGFBPS','PTK6','AURKA','TPX2','PLK1','SARS','COV','ROBO','UP','DN','V1','V2'}
    text=' '.join(w if w in acronyms else w.lower() for w in words)
    return text[0].upper()+text[1:]


def save(fig,stem):
    for ext in ('pdf','svg','png'):
        fig.savefig(FIG/f'{stem}.{ext}',dpi=300,facecolor='white')
    plt.close(fig)


def heatmap(results,collection,columns,paths,stem,title,caption):
    frame=results[results.collection.eq(collection)]
    nes=frame.pivot(index='pathway',columns='analysis',values='NES').reindex(index=paths,columns=columns)
    q=frame.pivot(index='pathway',columns='analysis',values='padj').reindex(index=paths,columns=columns)
    values=nes.to_numpy(float)
    bound=max(2,math.ceil(np.nanmax(np.abs(values))*2)/2)
    height=max(6.5,len(paths)*(.42 if collection=='Reactome' else .235)+2.25)
    width=11.8 if len(columns)>4 else 10.7
    fig,ax=plt.subplots(figsize=(width,height))
    fig.subplots_adjust(left=.48 if collection=='Reactome' else .39,right=.90,top=1-.9/height,bottom=1.0/height)
    im=ax.imshow(values,aspect='auto',cmap=CMAP,norm=TwoSlopeNorm(vmin=-bound,vcenter=0,vmax=bound),interpolation='nearest')
    ax.set_xticks(range(len(columns)),[LABELS.get(c,c) for c in columns],fontsize=9)
    ax.xaxis.tick_top();ax.tick_params(axis='both',length=0,pad=7)
    ax.set_yticks(range(len(paths)),[textwrap.fill(pretty(p),48 if collection=='Reactome' else 42) for p in paths],fontsize=8.5)
    ax.set_xticks(np.arange(-.5,len(columns),1),minor=True);ax.set_yticks(np.arange(-.5,len(paths),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=.7);ax.tick_params(which='minor',length=0)
    for i in range(len(paths)):
        for j in range(len(columns)):
            v=values[i,j];fdr=q.iloc[i,j]
            star='**' if fdr<.01 else '*' if fdr<.05 else ''
            label='NA' if not np.isfinite(v) else f'{v:.2f}{star}'
            ax.text(j,i,label,ha='center',va='center',fontsize=7.4 if len(columns)>4 else 8.1,color='white' if abs(v)>bound*.7 else '#203040')
    cax=fig.add_axes([.923,1.0/height,.018,1-2/height]);cb=fig.colorbar(im,cax=cax);cb.set_label('Normalized enrichment score (NES)',fontsize=9)
    fig.suptitle(title,x=.025,y=.992,ha='left',fontsize=14.5,fontweight='bold')
    fig.text(.025,.025,caption+'\n* pathway FDR < 0.05; ** < 0.01. BH within collection and analysis. Cells show unscaled NES.',fontsize=8.3,va='bottom')
    save(fig,stem)
    nes.to_csv(OUT/f'{stem}_NES.csv');q.to_csv(OUT/f'{stem}_FDR.csv')


def scatter(results,collection,pairs):
    fig,axes=plt.subplots(1,3,figsize=(15,6.5))
    fig.subplots_adjust(left=.065,right=.985,bottom=.29,top=.84,wspace=.25)
    colors={'Neither':'#cbd1d5','ATAC only':'#0072b2','RNA only':'#d55e00','Both':'#7b3294'}
    columns=MAIN[1:]
    allf=pairs[pairs.collection.eq(collection)]
    bx=max(2,math.ceil(allf.ATAC_NES.abs().max()*2)/2+.15)
    by=max(2,math.ceil(allf.RNA_NES.abs().max()*2)/2+.15)
    for ax,column in zip(axes,columns):
        data=allf[allf.RNA_analysis.eq(column)].copy()
        for group,color in colors.items():
            points=data[data.membership.eq(group)]
            ax.scatter(points.ATAC_NES,points.RNA_NES,s=30 if group=='Both' else 15,c=color,alpha=.85 if group=='Both' else .55,label=f'{group} ({len(points)})',linewidths=0)
        ax.axhline(0,color='#919ba4',linewidth=.8,linestyle='--');ax.axvline(0,color='#919ba4',linewidth=.8,linestyle='--')
        ax.set(xlim=(-bx,bx),ylim=(-by,by),xlabel='Promoter ATAC NES')
        if ax is axes[0]:ax.set_ylabel('RNA NES')
        both=int(data.membership.eq('Both').sum());same=int((data.membership.eq('Both')&data.concordant).sum())
        ax.set_title(LABELS[column].replace('\n',' ') + f'\nBoth q < 0.05: {both}; same direction: {same}',fontsize=10)
        chosen=data.sort_values(['worst_FDR','best_FDR','pathway']).head(5)
        labels=[ax.text(r.ATAC_NES-.15,r.RNA_NES+.3+.12*i,textwrap.fill(pretty(r.pathway),24),fontsize=6.5,color='#263949') for i,r in enumerate(chosen.itertuples())]
        np.random.seed(20260909)
        adjust_text(labels,x=data.ATAC_NES.to_numpy(),y=data.RNA_NES.to_numpy(),ax=ax,
                    target_x=chosen.ATAC_NES.to_numpy(),target_y=chosen.RNA_NES.to_numpy(),expand=(1.3,1.4),force_text=(.7,.7),force_static=(.2,.2),force_pull=(.005,.005),max_move=(20,20),
                    arrowprops={'arrowstyle':'-','color':'#84919b','lw':.5},ensure_inside_axes=True,iter_lim=800)
        ax.legend(loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2,fontsize=6.8,frameon=False,handletextpad=.1)
    fig.suptitle(f'{collection}: pathway concordance between ATAC and RNA',x=.025,y=.98,ha='left',fontweight='bold',fontsize=16)
    fig.text(.025,.03,'11,891 common eligible genes; no gene-level FDR cutoff. Each point is a pathway. Colors use collection-wise pathway FDR < 0.05.\nLabels: five pathways with the lowest maximum ATAC/RNA FDR per panel. The three panels reuse the same ATAC experiment.',fontsize=8.5)
    save(fig,collection.lower()+'_ATAC_RNA_NES_scatter')


def main():
    FIG.mkdir(exist_ok=True)
    results=pd.read_csv(OUT/'fgsea_all_results.csv')
    ranks=pd.read_csv(OUT/'gene_ranks.csv',dtype={'entrez_id':str}).set_index('entrez_id')
    settings=json.loads((OUT/'analysis_settings.json').read_text())
    assert len(ranks)==11891 and not results.duplicated(['analysis','pathway']).any()
    pairs=[];edges=[]
    atac=results[results.analysis.eq('ATAC')].set_index('pathway')
    for column in MAIN[1:]:
        rna=results[results.analysis.eq(column)].set_index('pathway')
        for path,a in atac.iterrows():
            r=rna.loc[path]
            aset=set(str(a.leading_edge_entrez).split(';')) if pd.notna(a.leading_edge_entrez) else set()
            rset=set(str(r.leading_edge_entrez).split(';')) if pd.notna(r.leading_edge_entrez) else set()
            shared=sorted(aset&rset)
            both=a.padj<.05 and r.padj<.05
            direction=bool(np.sign(a.NES)==np.sign(r.NES))
            member='Both' if both else 'ATAC only' if a.padj<.05 else 'RNA only' if r.padj<.05 else 'Neither'
            pairs.append({'pathway':path,'collection':a.collection,'RNA_analysis':column,'ATAC_NES':a.NES,'ATAC_FDR':a.padj,'RNA_NES':r.NES,'RNA_FDR':r.padj,
                'membership':member,'concordant':direction,'worst_FDR':max(a.padj,r.padj),'best_FDR':min(a.padj,r.padj),
                'ATAC_leading_edge_genes':len(aset),'RNA_leading_edge_genes':len(rset),'shared_leading_edge_genes':len(shared),
                'leading_edge_Jaccard':len(shared)/len(aset|rset) if aset|rset else 0,'shared_entrez':';'.join(shared),'shared_symbols':';'.join(ranks.loc[shared,'symbol'].fillna('NA'))})
            for gid in shared:
                g=ranks.loc[gid]
                edges.append({'pathway':path,'collection':a.collection,'RNA_analysis':column,'both_pathway_FDR005':both,'concordant_pathway_direction':direction,
                    'entrez_id':gid,'symbol':g.symbol,'chrom':g.chrom,'ATAC_median_logFC':g.ATAC_median_logFC,'ATAC_rank_stat':g.ATAC,
                    'RNA_logFC':g[column+'_logFC'],'RNA_FDR':g[column+'_FDR'],'RNA_rank_stat':g[column]})
    pairs=pd.DataFrame(pairs);edges=pd.DataFrame(edges)
    pairs.to_csv(OUT/'pathway_concordance.csv',index=False)
    edges.to_csv(OUT/'shared_leading_edge_genes.csv',index=False)
    summary=pairs.groupby(['collection','RNA_analysis']).apply(lambda g:pd.Series({'tested_pathways':len(g),'ATAC_FDR005':int((g.ATAC_FDR<.05).sum()),'RNA_FDR005':int((g.RNA_FDR<.05).sum()),'both_FDR005':int(g.membership.eq('Both').sum()),'both_concordant_FDR005':int((g.membership.eq('Both')&g.concordant).sum()),'both_opposing_FDR005':int((g.membership.eq('Both')&~g.concordant).sum())}),include_groups=False).reset_index()
    summary.to_csv(OUT/'pathway_concordance_summary.csv',index=False)
    hallmark=atac[atac.collection.eq('Hallmark')].sort_values('NES').index.tolist()
    heatmap(results,'Hallmark',MAIN,hallmark,'hallmark_pathway_heatmap','Pathway enrichment: promoter accessibility and RNA expression',
        f'All {len(hallmark)} testable Hallmark sets; 11,891 common genes. Negative NES: toward lower LOY signal; positive NES: toward higher LOY signal.\nRows ordered by ATAC NES. One shared ATAC experiment is compared with three related RNA results.')
    score=pairs[pairs.collection.eq('Reactome')].groupby('pathway').worst_FDR.min().sort_values(kind='stable')
    chosen=score.head(30).index.tolist()
    chosen=atac.loc[chosen].sort_values('NES').index.tolist()
    heatmap(results,'Reactome',MAIN,chosen,'reactome_pathway_heatmap','Reactome pathway enrichment: ATAC and RNA',
        '30 pathways selected by the lowest maximum of ATAC FDR and best RNA FDR; all results are supplied in CSV.\n11,891 common genes; rows ordered by ATAC NES. This selection score is descriptive and is not a joint FDR.')
    for collection in ('Hallmark','Reactome'):scatter(results,collection,pairs)
    sensitivities=['ATAC','ATAC_representative','ATAC_median_logFC','ATAC_passage','ATAC_no_chr19','ATAC_tie_shuffle_1','ATAC_tie_shuffle_2']
    heatmap(results,'Hallmark',sensitivities,hallmark,'hallmark_ATAC_sensitivity','Hallmark enrichment: sensitivity to ATAC ranking and model',
        'Identical gene universe except the chromosome 19 exclusion (11,070 genes). Values are raw NES; rows follow primary ATAC NES.\nPrimary: median signed √QL F; abundance-selected peak, median log₂FC, passage model and tied-rank orderings are sensitivity analyses.')
    robust=[]
    for column in sensitivities[1:]:
        s=results[results.analysis.eq(column)].set_index('pathway')
        for path,a in atac.iterrows():
            if path not in s.index:continue
            r=s.loc[path]
            robust.append({'pathway':path,'collection':a.collection,'sensitivity':column,'primary_NES':a.NES,'primary_FDR':a.padj,'sensitivity_NES':r.NES,'sensitivity_FDR':r.padj,
                           'same_direction':bool(np.sign(a.NES)==np.sign(r.NES)),'primary_FDR005':bool(a.padj<.05),'sensitivity_FDR005':bool(r.padj<.05)})
    robust=pd.DataFrame(robust);robust.to_csv(OUT/'ATAC_sensitivity_comparison.csv',index=False)
    robust_summary=robust.groupby(['collection','sensitivity']).apply(lambda g:pd.Series({'shared_testable_pathways':len(g),'primary_significant':int(g.primary_FDR005.sum()),'primary_significant_retained':int((g.primary_FDR005&g.sensitivity_FDR005&g.same_direction).sum()),'NES_spearman':g.primary_NES.corr(g.sensitivity_NES,method='spearman'),'max_abs_NES_change':float((g.primary_NES-g.sensitivity_NES).abs().max())}),include_groups=False).reset_index()
    robust_summary.to_csv(OUT/'ATAC_sensitivity_summary.csv',index=False)
    significant=pairs[pairs.membership.eq('Both')].copy()
    significant[['collection','RNA_analysis','pathway','ATAC_NES','ATAC_FDR','RNA_NES','RNA_FDR','concordant','shared_leading_edge_genes']].to_csv(OUT/'jointly_significant_pathways.csv',index=False)
    settings['status']='complete'
    settings['figures']=[p.name for p in FIG.glob('*.pdf')]
    (OUT/'analysis_settings.json').write_text(json.dumps(settings,indent=2),encoding='utf-8')
    methods='''Ranked pathway enrichment compares promoter ATAC with the combined, August 2025 and May 2026 RNA results on the intersection of 11,891 eligible genes. Existing chromosome/mapping/PAR exclusions are retained. No peak or gene FDR or fold-change filter is used before ranking.
Primary ATAC score: median of sign(logFC) * sqrt(edgeR QL F) over all tested linked promoter peaks. This is a ranking summary, not a calibrated gene-level test statistic. RNA score: original DESeq2 signed Wald statistic (source column t), with LOY / CTRL direction. Combined RNA uses the original batch + treatment fit; separate results use treatment fits. No differential model is refitted.
MSigDB human 2026.1 Hallmark and Reactome sets are obtained through the locally installed msigdbr 26.1.1. Human-native Entrez IDs are matched directly; each pathway is intersected with the rank universe. Duplicate memberships are removed. Set size range 15–500 after intersection; complete size and exclusion audits are provided.
fgsea 1.36.2 multilevel: scoreType std, gseaParam 1, sampleSize 101, nPermSimple 5000, eps 0; seeded serial execution. BH correction is within each collection and analysis. padj_global supplies the alternative correction across both collections within each analysis. log2err is retained in the complete tables. Every main assay tests the same pathways with the same effective sizes.
Heatmaps show raw NES without row or column scaling. Color limits are symmetric around zero. No hierarchical clustering is used: rows are ordered by primary ATAC NES, columns stay in assay/batch order. Hallmark includes all eligible sets. Reactome displays 30 sets chosen by the smallest maximum of ATAC FDR and best RNA FDR, then ordered by ATAC NES. This display-selection score is not a calibrated joint FDR.
Concordance means matching NES signs. Jointly significant means each separate assay passes collection-wise FDR < 0.05; it does not imply a separately calibrated joint FDR. Shared leading edges are intersections of fgsea leading-edge gene IDs. All intersections are supplied, including nonsignificant or opposite-direction pathways, with explicit flags. The combined RNA result overlaps the two constituent batch results; the same ATAC data are reused in all three comparisons.
Sensitivity analyses: one promoter peak per gene selected by highest primary logCPM (ties by peak ID), median promoter logFC, median signed sqrt(F) from the existing passage model, exclusion of chr19 from all four rankings, and two seeded reorderings confined to equal ATAC scores. The latter use the primary ATAC fgsea seed to assess tied-rank sensitivity.
Shared promoters create correlated/identical ATAC ranks; deterministic primary tie ordering uses gene ID. Tied-rank tests assess ordering sensitivity, not inter-gene correlation. These preranked tests do not explicitly model shared-peak or inter-gene dependence. Correlation-aware gene-set testing remains a potential follow-up. No copy-number correction, TF-binding inference or pathway-activity measurement is implied by NES. Removing chr19 is a sensitivity check, not a substitute for copy-number adjustment.
References: fgsea https://bioconductor.org/packages/release/bioc/html/fgsea.html ; MSigDB https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
'''
    (OUT/'METHODS.txt').write_text(methods,encoding='utf-8')
    main_caption='Hallmark pathway enrichment in HT1080 LOY versus CTRL. Columns show promoter ATAC and the existing combined, August 2025 and May 2026 RNA results. All four analyses use the same 11,891 eligible genes without a gene-level significance cutoff. ATAC ranks use the median signed square root of the edgeR QL F statistic across linked promoter peaks; RNA ranks use original signed DESeq2 Wald statistics. Colors and printed values show unscaled normalized enrichment scores (NES); asterisks mark BH FDR within the Hallmark collection (* <0.05; ** <0.01). Rows are ordered by ATAC NES. Only sets containing 15–500 eligible genes are tested. Positive/negative NES indicates relative enrichment toward more positive/negative statistics, respectively. NES does not directly measure pathway activity. The ATAC experiment is shared across comparisons, and the combined RNA model includes the two displayed RNA batches.'
    (OUT/'figure_caption.txt').write_text(main_caption,encoding='utf-8')
    def table(frame):return frame.to_html(index=False,border=0,float_format=lambda v:f'{v:.4g}',escape=True)
    def figure(stem,caption):
        return f'<figure><a href="../integration/pathway_GSEA/figures/{stem}.pdf"><img src="../integration/pathway_GSEA/figures/{stem}.png" alt="{escape(caption)}"></a><figcaption>{escape(caption)} <a href="../integration/pathway_GSEA/figures/{stem}.pdf">Vector PDF</a> · <a href="../integration/pathway_GSEA/figures/{stem}.svg">SVG</a></figcaption></figure>'
    hm=results[(results.analysis.eq('ATAC'))&(results.collection.eq('Hallmark'))&(results.padj<.05)][['pathway','NES','padj']].sort_values('NES')
    report='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HT1080 ATAC–RNA pathway enrichment</title><style>body{font:15px/1.5 Arial,sans-serif;max-width:1350px;margin:25px auto;padding:0 25px;color:#233746}h1,h2{color:#17354b}a{color:#0072b2}figure{margin:24px 0}img{width:100%;height:auto}figcaption,.small{font-size:13px;color:#536878}.note{background:#edf3f7;padding:18px;border-left:4px solid #0072b2}table{border-collapse:collapse;font-size:12px;width:100%}td,th{padding:7px;border-bottom:1px solid #dce3e8;text-align:left}.scroll{overflow:auto}pre{white-space:pre-wrap}</style><body><p><a href="ATACseq_report.html">Main ATAC report</a> · <a href="interactive_ATAC_RNA.html">Interactive gene explorer</a></p><h1>Pathway enrichment: promoter ATAC and RNA</h1>'''
    report+='<p class="note"><strong>Interpretation:</strong> primary ATAC shows enrichment of cell-cycle and MYC-related gene sets toward more positive accessibility statistics. RNA pathway support varies by batch. No tested pathway passes FDR 0.05 in both ATAC and combined RNA, or in both ATAC and May RNA. August RNA shares three concordant pathways with ATAC; one additional jointly significant pathway has opposing directions. This does not support a uniformly concordant pathway response across RNA batches.</p>'
    report+='<p>11,891 common eligible genes; 49 Hallmark and 933 Reactome sets pass the 15–500 gene size filter. No input gene-level FDR cutoff. Pathway FDR is BH-adjusted within each collection and analysis. Results use MSigDB 2026.1 and fgsea 1.36.2.</p>'
    report+='<h2>Hallmark pathway heatmap</h2>'+figure('hallmark_pathway_heatmap',main_caption)
    report+='<h2>Pathway support across assays</h2><div class="scroll">'+table(summary)+'</div>'
    report+=figure('hallmark_ATAC_RNA_NES_scatter','Hallmark pathway NES comparisons. Labels identify the five pathways with the lowest maximum ATAC/RNA FDR per panel; labels do not imply significance in both assays.')
    report+='<h2>Reactome pathways</h2>'+figure('reactome_pathway_heatmap','Thirty Reactome pathways selected for the strongest joint statistical support across ATAC and any RNA result. The complete results include all 933 tested sets.')
    report+=figure('reactome_ATAC_RNA_NES_scatter','Reactome pathway NES comparisons, with all tested pathways shown.')
    report+='<h2>Jointly significant pathways and shared leading-edge genes</h2><div class="scroll">'+table(significant[['RNA_analysis','pathway','ATAC_NES','ATAC_FDR','RNA_NES','RNA_FDR','concordant','shared_leading_edge_genes']])+'</div>'
    report+='<p><a href="../integration/pathway_GSEA/pathway_concordance.csv">All pathway comparisons and shared gene lists</a> · <a href="../integration/pathway_GSEA/shared_leading_edge_genes.csv">Shared leading-edge genes with gene-level ATAC/RNA effects</a></p>'
    report+='<h2>Primary ATAC Hallmark results</h2>'+table(hm)
    report+='<h2>Sensitivity checks</h2>'+figure('hallmark_ATAC_sensitivity','Alternative promoter summaries, the passage-adjusted model, chromosome 19 exclusion and tied-rank ordering. Asterisks use pathway FDR within each analysis.')
    report+='<p>About 14.5% of primary ATAC ranks are tied, often because genes share a promoter peak. Tie reordering tests address ordering, while shared-peak correlation remains a limitation of these preranked enrichment tests. Chr19 exclusion is a sensitivity analysis; copy-number effects have not been adjusted.</p><div class="scroll">'+table(robust_summary)+'</div>'
    report+='<h2>Data and methods</h2><p><a href="../integration/pathway_GSEA/fgsea_all_results.csv">All primary and sensitivity results</a> · <a href="../integration/pathway_GSEA/gene_ranks.csv">Exact gene ranks</a> · <a href="../integration/pathway_GSEA/sources/gene_set_membership.csv">Frozen gene sets</a> · <a href="../integration/pathway_GSEA/METHODS.txt">Methods</a> · <a href="../integration/pathway_GSEA/figure_caption.txt">Figure caption</a></p><details><summary>Complete methods</summary><pre>'+escape(methods)+'</pre></details></body></html>'
    if (ROOT/'report/ribosomal_GSEA_sensitivity.html').exists():
        report=report.replace('</body>','<section id="ribosomal-sensitivity"><h2>Ribosomal protein exclusion sensitivity</h2><p><a href="ribosomal_GSEA_sensitivity.html">Compare the original GSEA with cytosolic-only and cytosolic-plus-mitochondrial ribosomal protein exclusions</a>.</p></section></body>')
    if (ROOT/'report/pathway_gene_drivers.html').exists():
        report=report.replace('</body>','<section id="pathway-gene-drivers"><h2>Inspect the genes driving enrichment</h2><p><a href="pathway_gene_drivers.html">E2F, MYC and interferon: gene effects, individual RNA and ATAC samples, and interactive promoter peak inspection</a>.</p></section></body>')
    if (ROOT/'report/RNA_interaction_CN.html').exists():
        report=report.replace('</body>','<section id="rna-interaction-cn"><h2>Formal RNA batch interaction and D0 copy-number context</h2><p><a href="RNA_interaction_CN.html">Interaction tests, correlation sensitivity, hg38 WGS profiles and joint gene tables</a>.</p></section></body>')
    (ROOT/'report/pathway_GSEA.html').write_text(report,encoding='utf-8')
    completion={'status':'figures_ready_pending_validation','common_genes':len(ranks),'pathway_summary':summary.to_dict('records'),'figures':settings['figures'],
                'source_checksums':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),OUT/'fgsea_all_results.csv',ROOT/'code/57_pathway_GSEA.R']}}
    (OUT/'completion.json').write_text(json.dumps(completion,indent=2),encoding='utf-8')
    print(json.dumps(completion,indent=2))


if __name__=='__main__':main()
