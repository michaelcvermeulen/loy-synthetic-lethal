"""Compare ribosomal-excluded GSEA with the unchanged primary analysis."""
from pathlib import Path
from html import escape
import hashlib
import importlib.util
import json
import math
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'integration/pathway_GSEA';OUT=BASE/'ribosomal_sensitivity';FIG=OUT/'figures'
FIG.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('plot_base',ROOT/'code/58_plot_pathway_GSEA.py')
plot=importlib.util.module_from_spec(spec);spec.loader.exec_module(plot)
plot.OUT=OUT;plot.FIG=FIG
MAIN=plot.MAIN
MODES=['original','no_cytosolic','no_cytosolic_or_mitochondrial']
NAMES={'original':'Original','no_cytosolic':'Exclude cytosolic','no_cytosolic_or_mitochondrial':'Exclude cytosolic + mitochondrial'}
GENES={'original':11891,'no_cytosolic':11822,'no_cytosolic_or_mitochondrial':11749}


def bh(values):
    values=np.asarray(values,float);ix=np.argsort(values);out=np.empty(len(ix));out[ix]=np.minimum(1,np.minimum.accumulate((values[ix]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1]);return out


def comparison_heatmap(data,collection,paths,stem):
    columns=pd.MultiIndex.from_product([MODES,MAIN],names=['exclusion','analysis'])
    frame=data[data.collection.eq(collection)]
    nes=frame.pivot(index='pathway',columns=['exclusion','analysis'],values='NES').reindex(index=paths,columns=columns)
    q=frame.pivot(index='pathway',columns=['exclusion','analysis'],values='padj').reindex(index=paths,columns=columns)
    values=nes.to_numpy();bound=max(2,math.ceil(np.nanmax(abs(values))*2)/2)
    height=len(paths)*(.42 if collection=='Reactome' else .245)+2.25
    fig,ax=plt.subplots(figsize=(17,height));fig.subplots_adjust(left=.35,right=.93,top=1-1.25/height,bottom=1/height)
    im=ax.imshow(values,aspect='auto',interpolation='nearest',cmap=plot.CMAP,norm=TwoSlopeNorm(vmin=-bound,vcenter=0,vmax=bound))
    ax.set_yticks(range(len(paths)),[textwrap.fill(plot.pretty(p),49) for p in paths],fontsize=8)
    ax.set_xticks(range(12),[plot.LABELS[c] for m,c in columns],fontsize=8);ax.xaxis.tick_top();ax.tick_params(length=0,pad=7)
    ax.set_xticks(np.arange(-.5,12,1),minor=True);ax.set_yticks(np.arange(-.5,len(paths),1),minor=True);ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',length=0)
    for x in (3.5,7.5):ax.axvline(x,color='white',lw=6)
    for i in range(len(paths)):
        for j in range(12):
            v=values[i,j];fdr=q.iloc[i,j];stars='**' if fdr<.01 else '*' if fdr<.05 else ''
            ax.text(j,i,'NA' if not np.isfinite(v) else f'{v:.2f}{stars}',ha='center',va='center',fontsize=7,color='white' if abs(v)>bound*.7 else '#233746')
    for k,mode in enumerate(MODES):ax.text((k*4+2)/12,1+(.72/height),f'{NAMES[mode]}\n{GENES[mode]:,} genes',transform=ax.transAxes,ha='center',va='bottom',fontsize=10,fontweight='bold')
    cax=fig.add_axes([.947,1/height,.012,1-2.25/height]);fig.colorbar(im,cax=cax,label='Normalized enrichment score (NES)')
    fig.suptitle(f'{collection}: sensitivity to ribosomal protein exclusion',x=.025,y=.99,ha='left',fontsize=17,fontweight='bold')
    fig.text(.025,.024,'Identical exclusions across ATAC and all three RNA rankings. Remaining gene statistics are unchanged; pathway membership and FDR are recalculated.\nRows and color scale are shared across scenarios. * collection-wise pathway FDR < 0.05; ** < 0.01. Gray / NA: outside the 15–500 gene size filter.',fontsize=9)
    plot.save(fig,stem);nes.to_csv(OUT/f'{stem}_NES.csv');q.to_csv(OUT/f'{stem}_FDR.csv')


def main():
    baseline=pd.read_csv(BASE/'fgsea_all_results.csv');baseline=baseline[baseline.analysis.isin(MAIN)].assign(exclusion='original')
    filtered=pd.read_csv(OUT/'fgsea_exclusion_results.csv')
    all_results=pd.concat([baseline,filtered],ignore_index=True)
    excluded=pd.read_csv(OUT/'excluded_genes.csv',dtype={'entrez_id':str})
    ranks=pd.read_csv(BASE/'gene_ranks.csv',dtype={'entrez_id':str}).set_index('entrez_id')
    sets=pd.read_csv(BASE/'sources/gene_set_membership.csv',dtype={'entrez_id':str}).groupby('pathway').entrez_id.agg(set).to_dict()
    # Re-adjust within the identical collection-specific pathway family across all scenarios.
    all_results['FDR_common_pathway_family']=np.nan
    for collection in all_results.collection.unique():
        family=set.intersection(*[set(all_results.loc[all_results.collection.eq(collection)&all_results.exclusion.eq(m),'pathway']) for m in MODES])
        for mode in MODES:
            for analysis in MAIN:
                ix=all_results.index[all_results.collection.eq(collection)&all_results.exclusion.eq(mode)&all_results.analysis.eq(analysis)&all_results.pathway.isin(family)]
                all_results.loc[ix,'FDR_common_pathway_family']=bh(all_results.loc[ix,'pval'])
    all_results.to_csv(OUT/'all_scenarios_results.csv',index=False)
    summary=all_results.groupby(['exclusion','analysis','collection']).apply(lambda g:pd.Series({'tested_pathways':len(g),'FDR005':int(g.padj.lt(.05).sum()),'positive_FDR005':int((g.padj.lt(.05)&g.NES.gt(0)).sum()),'negative_FDR005':int((g.padj.lt(.05)&g.NES.lt(0)).sum())}),include_groups=False).reset_index()
    summary.to_csv(OUT/'significant_pathway_counts.csv',index=False)
    comparisons=[];concordance=[];concordance_rows=[]
    base=all_results[all_results.exclusion.eq('original')]
    for mode in MODES[1:]:
        removed=set(excluded.loc[excluded['class'].eq('cytosolic'),'entrez_id']) if mode=='no_cytosolic' else set(excluded.entrez_id)
        merged=base.merge(all_results[all_results.exclusion.eq(mode)],on=['analysis','collection','pathway'],how='outer',suffixes=('_original','_excluded'))
        merged['scenario']=mode
        merged['same_direction']=np.sign(merged.NES_original)==np.sign(merged.NES_excluded)
        merged['newly_significant']=merged.padj_original.ge(.05)&merged.padj_excluded.lt(.05)
        merged['lost_significance']=merged.padj_original.lt(.05)&merged.padj_excluded.ge(.05)
        merged['not_tested_after']=merged.NES_excluded.isna()
        merged['newly_testable']=merged.NES_original.isna()
        merged['NES_change']=merged.NES_excluded-merged.NES_original
        fractions=[];counts=[]
        for r in merged.itertuples():
            edge=set(str(r.leading_edge_entrez_original).split(';')) if pd.notna(r.leading_edge_entrez_original) else set()
            weights=ranks.loc[sorted(edge),r.analysis].abs() if edge else pd.Series(dtype=float)
            fractions.append(float(weights[weights.index.isin(removed)].sum()/weights.sum()) if weights.sum()>0 else 0)
            counts.append(len(sets[r.pathway]&removed))
        merged['removed_baseline_leading_edge_weight_fraction']=fractions
        merged['removed_pathway_genes']=counts
        comparisons.append(merged)
    comparisons=pd.concat(comparisons,ignore_index=True)
    comparisons.to_csv(OUT/'pathway_changes.csv',index=False)
    for mode in MODES:
        d=all_results[all_results.exclusion.eq(mode)]
        a=d[d.analysis.eq('ATAC')]
        for assay in MAIN[1:]:
            z=a.merge(d[d.analysis.eq(assay)],on=['collection','pathway'],suffixes=('_ATAC','_RNA'))
            z['scenario']=mode;z['RNA_result']=assay
            z['both_FDR005']=z.padj_ATAC.lt(.05)&z.padj_RNA.lt(.05)
            z['concordant']=np.sign(z.NES_ATAC)==np.sign(z.NES_RNA)
            z['shared_leading_edge_genes']=[len(set(str(x).split(';'))&set(str(y).split(';'))) for x,y in zip(z.leading_edge_entrez_ATAC,z.leading_edge_entrez_RNA)]
            concordance_rows.append(z)
            for collection,g in z.groupby('collection'):
                concordance.append({'scenario':mode,'RNA_result':assay,'collection':collection,'tested_pathways':len(g),'both_FDR005':int(g.both_FDR005.sum()),'same_direction_both_FDR005':int((g.both_FDR005&g.concordant).sum()),'opposing_both_FDR005':int((g.both_FDR005&~g.concordant).sum())})
    concordance=pd.DataFrame(concordance);concordance.to_csv(OUT/'ATAC_RNA_concordance_summary.csv',index=False)
    concordance_rows=pd.concat(concordance_rows,ignore_index=True)
    concordance_rows.to_csv(OUT/'ATAC_RNA_pathway_comparisons.csv',index=False)
    concordance_rows[concordance_rows.both_FDR005].to_csv(OUT/'jointly_significant_pathways.csv',index=False)
    hallmark=baseline[baseline.analysis.eq('ATAC')&baseline.collection.eq('Hallmark')].sort_values('NES').pathway.tolist()
    reactome=pd.read_csv(BASE/'reactome_pathway_heatmap_NES.csv',index_col=0).index.tolist()
    comparison_heatmap(all_results,'Hallmark',hallmark,'hallmark_ribosomal_comparison')
    comparison_heatmap(all_results,'Reactome',reactome,'reactome_ribosomal_comparison')
    for mode in MODES[1:]:
        plot.heatmap(all_results[all_results.exclusion.eq(mode)],'Hallmark',MAIN,hallmark,f'hallmark_{mode}',
            'Pathway enrichment after ribosomal protein exclusion',
            f'{NAMES[mode]}; {GENES[mode]:,} common genes. The primary ATAC and RNA statistics of retained genes are unchanged.\nRows retain the original ATAC NES order; this is a sensitivity analysis with recomputed gene-set membership and FDR.')
    focus=['HALLMARK_MYC_TARGETS_V1','HALLMARK_MYC_TARGETS_V2','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT']
    focus_table=all_results[all_results.pathway.isin(focus)][['exclusion','analysis','pathway','size','NES','padj','FDR_common_pathway_family']]
    focus_table.to_csv(OUT/'MYC_E2F_cell_cycle_comparison.csv',index=False)
    # Summarize test-family changes separately from changes among jointly testable pathways.
    changes_summary=comparisons.groupby(['scenario','analysis','collection']).apply(lambda g:pd.Series({'newly_significant':int(g.newly_significant.sum()),'lost_significance':int(g.lost_significance.sum()),'not_tested_after':int(g.not_tested_after.sum()),'newly_testable':int(g.newly_testable.sum())}),include_groups=False).reset_index()
    changes_summary.to_csv(OUT/'pathway_change_summary.csv',index=False)
    methods='''Ribosomal exclusion sensitivity of the existing promoter ATAC/RNA GSEA. Baseline uses 11,891 common eligible genes. The exclusions are fixed from HGNC groups 728 (S ribosomal proteins), 729 (L ribosomal proteins), and 646 (mitochondrial ribosomal proteins), retrieved and frozen with checksums. Group 728/729 includes FAU and UBA52. Cytosolic-only exclusion removes 69 eligible genes, leaving 11,822. Cytosolic plus mitochondrial exclusion removes 142, leaving 11,749.
Matching uses existing Entrez ID OR exact approved HGNC symbol. One source mapping discrepancy is explicit: the existing rank table calls Entrez 102465978 RPL8, while HGNC assigns RPL8 to Entrez 6132. That row is excluded by its RNA symbol; original gene mappings and scores are otherwise unchanged. The discrepancy is recorded in exclusion_mapping_discrepancies.csv. This is an HGNC protein-family exclusion, not removal of all ribosome-associated genes, rRNA genes, pseudogenes, or ribosome-biogenesis factors. RPS6K kinases and antisense genes are retained. Nuclear-encoded mitochondrial ribosomal proteins are separate from mitochondrial-chromosome transcripts.
Each exclusion is identical across ATAC and all three RNA results. Retained gene scores match the primary gene_ranks.csv table; no differential-expression/accessibility model or count normalization is refitted. The frozen MSigDB 2026.1 Hallmark and Reactome memberships are re-intersected with the reduced rank universe. The same fgseaMultilevel settings and assay-specific seeds are used: gseaParam 1, scoreType std, nPermSimple 5000, sampleSize 101, eps 0, minimum 15 and maximum 500 effective genes, serial execution.
P-values, NES, leading edges and BH FDR are recomputed. Default FDR is within each collection/scenario/assay. FDR_common_pathway_family repeats BH on only pathways testable in every scenario within that collection, to distinguish test-family changes from enrichment changes. padj_global is available across Hallmark plus Reactome within each scenario/assay. A disappearing pathway may have fallen below the size filter and is not automatically a nonsignificant result; eligibility changes are recorded separately.
Heatmaps use raw NES, a shared symmetric scale within each comparison figure, fixed original row order and fixed assay order. Reactome retains the 30 pathways selected for the original figure. Gray/NA denotes an untestable set. Stars mark collection-wise pathway FDR (<0.05, <0.01). Original results remain the primary analysis. Comparisons after exclusion test enrichment among non-ribosomal genes and do not automatically represent improved statistical validity. Filtering changes pathway definitions and the ranked background. Shared promoter/inter-gene correlation remains a limitation of preranked GSEA.
Sources: https://rest.genenames.org/fetch/gene_group_id/728 ; https://rest.genenames.org/fetch/gene_group_id/729 ; https://rest.genenames.org/fetch/gene_group_id/646 ; https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/
'''
    (OUT/'METHODS.txt').write_text(methods,encoding='utf-8')
    def table(frame):return frame.to_html(index=False,border=0,float_format=lambda v:f'{v:.4g}',escape=True)
    def fig(stem,caption):return f'<figure><a href="../integration/pathway_GSEA/ribosomal_sensitivity/figures/{stem}.pdf"><img src="../integration/pathway_GSEA/ribosomal_sensitivity/figures/{stem}.png" alt="{escape(caption)}"></a><figcaption>{escape(caption)} <a href="../integration/pathway_GSEA/ribosomal_sensitivity/figures/{stem}.pdf">Vector PDF</a> · <a href="../integration/pathway_GSEA/ribosomal_sensitivity/figures/{stem}.svg">SVG</a></figcaption></figure>'
    combined_shared=concordance[(concordance.RNA_result=='RNA_Combined')].both_FDR005.sum()
    text=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Ribosomal exclusion sensitivity</title><style>body{{font:15px/1.5 Arial,sans-serif;max-width:1500px;margin:30px auto;padding:0 25px;color:#233746}}a{{color:#0072b2}}img{{width:100%;height:auto}}figure{{margin:25px 0}}table{{font-size:12px;border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid #dbe4ea;padding:7px;text-align:left}}.scroll{{overflow:auto}}.note{{background:#edf3f7;padding:18px}}pre{{white-space:pre-wrap}}</style><body><p><a href="pathway_GSEA.html">Original pathway report</a> · <a href="ATACseq_report.html">Main ATAC report</a></p><h1>Does removing ribosomal proteins change the GSEA conclusions?</h1>
<p class="note">The main ATAC MYC and E2F enrichments persist after both exclusions. August RNA loses many significant Reactome results, while combined RNA retains four significant Hallmark sets. Removing cytosolic genes yields two significant May RNA Reactome sets, but neither remains significant when mitochondrial ribosomal proteins are also excluded. This indicates sensitivity of some pathway results without resolving the overall ATAC–RNA discrepancy.</p>
<p>Original: 11,891 genes. Exclude cytosolic ribosomal proteins: 11,822 genes (69 removed). Exclude cytosolic and mitochondrial ribosomal proteins: 11,749 genes (142 removed). Same exclusions in all four assays; retained gene statistics are unchanged. The original analysis remains available.</p>'''
    text+='<h2>Hallmark comparison</h2>'+fig('hallmark_ribosomal_comparison','All 49 originally testable Hallmark sets, with fixed row order and a common NES scale across the original and two exclusion scenarios.')
    text+='<h2>Reactome comparison</h2>'+fig('reactome_ribosomal_comparison','The same 30 Reactome pathways shown in the original figure. Gray cells are outside the size filter after exclusion.')
    text+='<h2>Number of significant pathways</h2><div class="scroll">'+table(summary)+'</div>'
    text+='<h2>ATAC–RNA agreement</h2><div class="scroll">'+table(concordance)+'</div>'
    text+='<h2>MYC and E2F results</h2><div class="scroll">'+table(focus_table[focus_table.analysis.eq('ATAC')])+'</div>'
    text+='<h2>Standalone figures after exclusion</h2>'+fig('hallmark_no_cytosolic','Hallmark results after cytosolic ribosomal protein exclusion.')+fig('hallmark_no_cytosolic_or_mitochondrial','Hallmark results after cytosolic and mitochondrial ribosomal protein exclusion.')
    text+='<h2>Downloads and methods</h2><p>'+ ' · '.join(f'<a href="../integration/pathway_GSEA/ribosomal_sensitivity/{p}">{label}</a>' for p,label in [('excluded_genes.csv','Exact excluded genes'),('pathway_changes.csv','Pathway changes and removed leading-edge contribution'),('all_scenarios_results.csv','All GSEA results'),('ATAC_RNA_pathway_comparisons.csv','ATAC–RNA comparisons'),('prefix_matches_retained.csv','Similarly named genes retained'),('exclusion_mapping_discrepancies.csv','Mapping discrepancy audit'),('METHODS.txt','Methods')])+'</p>'
    text+='<p>HGNC family definitions and source responses are frozen under sources/. <a href="https://www.genenames.org/help/rest/">HGNC REST documentation</a>. The <a href="https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/">GSEA user guide</a> does not generally recommend filtering genes to improve enrichment; these are explicitly reported as sensitivity analyses.</p><details><summary>Complete methods</summary><pre>'+escape(methods)+'</pre></details></body></html>'
    (ROOT/'report/ribosomal_GSEA_sensitivity.html').write_text(text,encoding='utf-8')
    result={'status':'report_ready_pending_verification','excluded_cytosolic':69,'excluded_mitochondrial':73,'concordance':concordance.to_dict('records'),
            'baseline_results_sha256':hashlib.sha256((BASE/'fgsea_all_results.csv').read_bytes()).hexdigest(),'filtered_results_sha256':hashlib.sha256((OUT/'fgsea_exclusion_results.csv').read_bytes()).hexdigest()}
    (OUT/'completion.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(concordance.to_string(index=False))


if __name__=='__main__':main()
