"""Inspect all leading-edge genes in E2F, MYC and interferon sets without refitting models."""
from pathlib import Path
from html import escape
import hashlib
import json
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'integration/pathway_GSEA';OUT=BASE/'gene_drivers';FIG=OUT/'figures'
ASSETS=ROOT/'report/pathway_driver_assets'
for p in [OUT,FIG,ASSETS]:p.mkdir(parents=True,exist_ok=True)
PATHWAYS={'E2F':'HALLMARK_E2F_TARGETS','MYC_V1':'HALLMARK_MYC_TARGETS_V1','MYC_V2':'HALLMARK_MYC_TARGETS_V2','IFN_alpha':'HALLMARK_INTERFERON_ALPHA_RESPONSE','IFN_gamma':'HALLMARK_INTERFERON_GAMMA_RESPONSE'}
TITLES={'E2F':'E2F targets','MYC_V1':'MYC targets V1','MYC_V2':'MYC targets V2','IFN_alpha':'Interferon alpha response','IFN_gamma':'Interferon gamma response'}
MAIN=['ATAC','RNA_Combined','RNA_AUG_2025','RNA_MAY_2026'];SAMPLES=['CTRL_11','CTRL_12','CTRL_21','LOY_11','LOY_12','LOY_21']
BLUE='#0072b2';ORANGE='#d55e00';PURPLE='#7b3294'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})


def save(fig,stem):
    for ext in ['pdf','png','svg']:fig.savefig(FIG/f'{stem}.{ext}',dpi=300,facecolor='white')
    plt.close(fig)


def main():
    ranks=pd.read_csv(BASE/'gene_ranks.csv',dtype={'entrez_id':str}).set_index('entrez_id')
    results=pd.read_csv(BASE/'fgsea_all_results.csv');results=results[results.analysis.isin(MAIN)&results.pathway.isin(PATHWAYS.values())]
    membership=pd.read_csv(BASE/'sources/gene_set_membership.csv',dtype={'entrez_id':str})
    sets={p:set(membership.loc[membership.pathway.eq(p),'entrez_id'])&set(ranks.index) for p in PATHWAYS.values()}
    leading={(r.pathway,r.analysis):set(str(r.leading_edge_entrez).split(';')) for r in results.itertuples()}
    all_drivers=sorted(set.union(*leading.values()));assert len(all_drivers)==379
    rna_source=pd.read_csv(ROOT/'config/RNAseq/ht1080_batches_LOY_v_CTRL_DESeq2_batch_treatment.csv',dtype={'entrezid':str}).set_index('entrezid')
    rna_sources={'RNA_Combined':rna_source}
    for batch in ['AUG_2025','MAY_2026']:rna_sources['RNA_'+batch]=pd.read_csv(ROOT/f'integration/batch_RNA_overlaps/sources/ht1080_{batch}_LOY_v_CTRL_DESeq2.csv',dtype={'entrezid':str}).set_index('entrezid')
    meta=pd.read_csv(ROOT/'integration/interactive_atac_rna/RNA_sample_metadata.csv')
    meta=meta[meta.treatment.isin(['CTRL','LOY'])].sort_values(['batch','treatment','sample_id']).reset_index(drop=True)
    assert len(meta)==12 and (meta.groupby(['batch','treatment']).size()==3).all()
    norm=pd.read_csv(ROOT/'integration/interactive_atac_rna/RNA_normalized_counts.csv').set_index('gene_id_full')
    raw=pd.read_csv(ROOT/'integration/interactive_atac_rna/RNA_raw_counts.csv').set_index('gene_id')
    sizes=pd.read_csv(ROOT/'integration/interactive_atac_rna/RNA_size_factors.csv').set_index('sample').size_factor
    peaks=pd.read_csv(ROOT/'differential/primary_TMM_all.csv').set_index('peak_id')
    counts=pd.read_csv(ROOT/'counts/peak_counts.tsv',sep='\t').set_index('peak_id')
    nf=pd.read_csv(ROOT/'differential/normalization_factors.csv');effective=nf[nf.arm.eq('primary_TMM')].set_index('sample').loc[SAMPLES,'effective_size']
    cpm=counts[SAMPLES].div(effective,axis=1)*1e6
    links=pd.read_csv(BASE/'promoter_peak_ranking_audit.csv',dtype={'entrez_id':str})
    ribosomal=pd.read_csv(BASE/'ribosomal_sensitivity/excluded_genes.csv',dtype={'entrez_id':str}).set_index('entrez_id')['class'].to_dict()
    gene_data={};expression_rows=[];peak_rows=[];robust_rows=[]
    for gid in all_drivers:
        rank=ranks.loc[gid];source=rna_source.loc[gid]
        assert isinstance(source,pd.Series)
        normalized=norm.loc[source.gene_id_full,meta.colname].to_numpy(float)
        rawcounts=raw.loc[source.gene_id_full,meta.colname].to_numpy(int)
        np.testing.assert_allclose(normalized,rawcounts/sizes.loc[meta.colname].to_numpy(),rtol=1e-12,atol=1e-10)
        rna=[]
        for a in MAIN[1:]:
            r=rna_sources[a].loc[gid]
            rna.append({'analysis':a,'logFC':float(r.logFC),'SE':float(r.lfcSE),'FDR':float(r.FDR)})
        linked=links[links.entrez_id.eq(gid)].sort_values(['logCPM','peak_id'],ascending=[False,True])
        plist=[]
        for p in linked.itertuples():
            record={'id':p.peak_id,'logFC':float(p.logFC),'FDR':float(p.FDR),'raw':counts.loc[p.peak_id,SAMPLES].astype(int).tolist(),'cpm':cpm.loc[p.peak_id,SAMPLES].tolist()}
            plist.append(record)
            peak_rows.append({'entrez_id':gid,'symbol':rank.symbol,'peak_id':p.peak_id,'representative':p.peak_id==rank.representative_peak,'logFC':p.logFC,'FDR':p.FDR,**{s+'_raw':record['raw'][i] for i,s in enumerate(SAMPLES)},**{s+'_TMM_CPM':record['cpm'][i] for i,s in enumerate(SAMPLES)}})
        gene_data[gid]={'symbol':rank.symbol,'chrom':rank.chrom,'atac_median_logFC':rank.ATAC_median_logFC,'representative_peak':rank.representative_peak,'ribosomal':ribosomal.get(gid,''),'rna':rna,'norm':normalized.tolist(),'raw':rawcounts.tolist(),'peaks':plist}
        for i,m in meta.iterrows():expression_rows.append({'entrez_id':gid,'symbol':rank.symbol,'sample':m.colname,'batch':m.batch,'condition':m.treatment,'raw':int(rawcounts[i]),'normalized':normalized[i],'log2_normalized_plus1':float(np.log2(normalized[i]+1))})
        log=np.log2(normalized+1)
        for batch in ['AUG_2025','MAY_2026']:
            ci=np.flatnonzero((meta.batch==batch)&(meta.treatment=='CTRL'));li=np.flatnonzero((meta.batch==batch)&(meta.treatment=='LOY'))
            cv,lv=log[ci],log[li];difference=float(lv.mean()-cv.mean())
            loo=[]
            for j in range(3):loo.extend([float(lv.mean()-np.delete(cv,j).mean()),float(np.delete(lv,j).mean()-cv.mean())])
            robust_rows.append({'entrez_id':gid,'symbol':rank.symbol,'batch':batch,'mean_log_expression_difference_LOY_CTRL':difference,
                'LOY_CTRL_pairwise_fraction_positive':float((lv[:,None]>cv[None,:]).mean()),'all_LOY_above_all_CTRL':bool(lv.min()>cv.max()),'all_LOY_below_all_CTRL':bool(lv.max()<cv.min()),
                'leave_one_out_min_difference':min(loo),'leave_one_out_max_difference':max(loo),'leave_one_out_same_sign':bool(np.all(np.sign(loo)==np.sign(difference)))})
    pd.DataFrame(expression_rows).to_csv(OUT/'RNA_sample_expression.csv',index=False)
    pd.DataFrame(peak_rows).to_csv(OUT/'ATAC_promoter_peak_measurements.csv',index=False)
    robustness=pd.DataFrame(robust_rows);robustness.to_csv(OUT/'RNA_sample_robustness.csv',index=False)
    pathway_data={};driver_rows=[];selection=[];summaries=[];concentration=[]
    for key,p in PATHWAYS.items():
        ids=sorted(sets[p]);frame=ranks.loc[ids]
        union=set.union(*(leading[(p,a)] for a in MAIN));weights={};candidate=[];reasons={}
        for a in MAIN:
            le=leading[(p,a)];w=ranks.loc[sorted(le),a].abs();w=w/w.sum();weights[a]=w.to_dict()
            order=w.sort_values(ascending=False,kind='stable')
            for gid in order.head(3).index:
                if gid not in candidate:candidate.append(gid)
                reasons.setdefault(gid,[]).append('top3_'+a+'_leading_weight')
            concentration.append({'pathway':p,'analysis':a,'leading_edge_genes':len(le),'top1_leading_weight_fraction':float(order.iloc[0]),'top5_leading_weight_fraction':float(order.head(5).sum()),
                'genes_for_50pct_leading_weight':int(np.searchsorted(order.cumsum().to_numpy(),.5)+1),'top5_symbols':';'.join(ranks.loc[order.head(5).index,'symbol'])})
        for gid in sorted(union,key=lambda g:(-max(weights[a].get(g,0) for a in MAIN),g)):
            if len(candidate)>=12:break
            if gid not in candidate:candidate.append(gid);reasons[gid]=['fill_by_max_leading_weight_across_assays']
        candidate=sorted(candidate,key=lambda g:(ranks.loc[g,'ATAC_median_logFC'],g))
        for gid in candidate:selection.append({'pathway':p,'entrez_id':gid,'symbol':ranks.loc[gid,'symbol'],'selection_reason':';'.join(reasons[gid])})
        for gid in ids:
            g=ranks.loc[gid]
            driver_rows.append({'pathway':p,'entrez_id':gid,'symbol':g.symbol,'any_leading_edge':gid in union,'selected_example':gid in candidate,'ribosomal':ribosomal.get(gid,''),
                'ATAC_median_logFC':g.ATAC_median_logFC,**{a+'_logFC':g[a+'_logFC'] for a in MAIN[1:]},**{a+'_FDR':g[a+'_FDR'] for a in MAIN[1:]},
                **{a+'_leading_edge':gid in leading[(p,a)] for a in MAIN},**{a+'_leading_weight_fraction':weights[a].get(gid,0) for a in MAIN}})
        aug=frame.RNA_AUG_2025_logFC;may=frame.RNA_MAY_2026_logFC
        sums={'pathway':p,'eligible_genes':len(ids),'leading_edge_union_genes':len(union),'ATAC_positive_genes':int(frame.ATAC_median_logFC.gt(0).sum()),
            'RNA_Aug_positive_genes':int(aug.gt(0).sum()),'RNA_May_positive_genes':int(may.gt(0).sum()),'RNA_both_negative_genes':int((aug.lt(0)&may.lt(0)).sum()),
            'RNA_both_positive_genes':int((aug.gt(0)&may.gt(0)).sum()),'RNA_opposite_sign_genes':int((np.sign(aug)!=np.sign(may)).sum()),
            'RNA_batch_effect_spearman':float(aug.corr(may,method='spearman')),'median_ATAC_logFC':float(frame.ATAC_median_logFC.median()),
            'median_RNA_Aug_logFC':float(aug.median()),'median_RNA_May_logFC':float(may.median())}
        summaries.append(sums)
        pathway_data[key]={'name':TITLES[key],'pathway':p,'eligible':ids,'drivers':sorted(union,key=lambda g:ranks.loc[g,'symbol']),'examples':candidate,
            'gsea':results[results.pathway.eq(p)][['analysis','NES','padj','size']].to_dict('records'),'summary':sums}
    all_effects=pd.DataFrame(driver_rows);all_effects.to_csv(OUT/'all_pathway_gene_effects.csv',index=False)
    pd.DataFrame(selection).to_csv(OUT/'selected_example_genes.csv',index=False)
    summary=pd.DataFrame(summaries);summary.to_csv(OUT/'pathway_direction_summary.csv',index=False)
    concentration=pd.DataFrame(concentration);concentration.to_csv(OUT/'leading_edge_concentration.csv',index=False)
    payload={'genes':gene_data,'pathways':pathway_data,'rna_samples':meta[['colname','sample_id','batch','treatment']].to_dict('records'),'atac_samples':SAMPLES}
    (ASSETS/'data.js').write_text('window.PATHWAY_DRIVERS='+json.dumps(payload,separators=(',',':'),allow_nan=False)+';',encoding='utf-8')
    draw_overview(ranks,pathway_data)
    for key,pdct in pathway_data.items():
        draw_examples(key,pdct,gene_data,meta)
        print(f'{key}: {len(pdct["drivers"])} leading-edge genes; {len(pdct["examples"])} plotted examples',flush=True)
    methods='''Gene inspection covers HALLMARK_E2F_TARGETS, HALLMARK_MYC_TARGETS_V1/V2, and HALLMARK_INTERFERON_ALPHA/GAMMA_RESPONSE from the original (ribosomal-inclusive) GSEA. All eligible pathway genes are included in the effect-distribution summaries. The interactive inspector contains the union of leading-edge genes across the four primary GSEA rankings: 379 distinct genes. Leading edges from nonsignificant pathways are explicitly retained as descriptive, exploratory gene lists; membership is not evidence of gene-level significance.
For each pathway, the static examples include the top three absolute-rank-statistic contributors within each assay's leading edge. Duplicates are collapsed and the list is filled to 12 by maximum within-leading-edge weight across the assays. This rule is applied regardless of gene significance or direction; exact selection reasons are saved. Static figures use a subset; the complete leading-edge gene union is available in the browser inspector and CSV files.
ATAC effect is the median log2 fold change of all tested linked promoter peaks. The representative peak is the highest-abundance promoter peak (primary logCPM), as fixed in the preceding sensitivity analysis. Its raw counts and primary-TMM CPM are shown for all six independent ATAC replicates. All linked promoter effects and FDR values are available in the inspector. The representative peak does not necessarily describe every promoter peak; the promoter median has no gene-level FDR.
RNA effects and 95% Wald intervals reuse source logFC and lfcSE (logFC +/- 1.96*SE); they are unadjusted intervals, whereas source FDR values are multiple-testing adjusted. No interaction model, new DE model, or count normalization is fitted. RNA dots show log2(normalized count + 1), using the shared combined-model DESeq2 size factors, with every CTRL/LOY sample shown separately by RNA batch. Counts are not batch-corrected. Individual RNA and ATAC samples are not treated as matched across assays.
Example expression heatmaps subtract each gene's mean log2(normalized count + 1) across all 12 RNA samples; no row standard-deviation scaling or batch correction is applied. Colors therefore represent within-gene log-expression deviations, not absolute abundance. Effect heatmaps use actual log2 fold changes and display each gene's ATAC median plus three original RNA effects. RNA stars use source gene FDR < 0.05; no ATAC gene-significance stars are shown.
Direction summaries use all eligible genes in each pathway, rather than only its selected examples. Leading-edge concentration measures the fraction of absolute ranking-statistic weight carried by the top one/five leading-edge genes and the number needed to reach 50% of that weight. This is descriptive, not a test of the pathway or a leave-gene-out GSEA refit. Pathway gene sets overlap.
RNA sample robustness records the mean difference of log2(normalized count + 1) between LOY and CTRL within each batch, all nine between-group sample comparisons, complete group separation, and the range after omitting each of the six samples once. These summaries retain original size factors and are descriptive; their mean log-expression differences need not equal the source DESeq2 coefficient. Opposite signs between batches do not by themselves establish a statistically significant interaction.
The full analysis and source mapping exclusions are inherited from integration/pathway_GSEA. Ribosomal-family membership is marked using the preceding HGNC exclusion audit, while the primary gene inspection retains those genes. No copy-number adjustment, direct TF-binding conclusion, or formal RNA batch-by-treatment interaction test is added.
'''
    (OUT/'METHODS.txt').write_text(methods,encoding='utf-8')
    write_report(summary,concentration,pathway_data,methods)
    sources=[BASE/'gene_ranks.csv',BASE/'fgsea_all_results.csv',ROOT/'integration/interactive_atac_rna/RNA_normalized_counts.csv',ROOT/'integration/interactive_atac_rna/RNA_raw_counts.csv',ROOT/'counts/peak_counts.tsv',ROOT/'differential/normalization_factors.csv',Path(__file__)]
    result={'status':'built_pending_verification','genes':len(all_drivers),'pathways':len(PATHWAYS),'pathway_summaries':summaries,
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (OUT/'completion.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(summary.to_string(index=False),flush=True)


def draw_overview(ranks,pathways):
    fig,axes=plt.subplots(1,5,figsize=(18,4.9));fig.subplots_adjust(left=.05,right=.995,bottom=.24,top=.77,wspace=.3)
    for ax,(key,p) in zip(axes,pathways.items()):
        d=ranks.loc[p['eligible']]
        lim=max(1.5,float(np.max(abs(d[['RNA_AUG_2025_logFC','RNA_MAY_2026_logFC']].to_numpy())))*1.08)
        colors=np.clip(d.ATAC_median_logFC,-1,1)
        im=ax.scatter(d.RNA_AUG_2025_logFC,d.RNA_MAY_2026_logFC,c=colors,cmap='RdBu_r',vmin=-1,vmax=1,s=20,alpha=.8,linewidths=0)
        ax.axhline(0,color='#a1aab2',lw=.6);ax.axvline(0,color='#a1aab2',lw=.6);ax.plot([-lim,lim],[-lim,lim],ls='--',lw=.6,color='#b9c1c8')
        ax.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel='August RNA log₂FC');ax.set_title(TITLES[key]+f'\nn = {len(d)}',fontsize=10)
        if ax is axes[0]:ax.set_ylabel('May RNA log₂FC')
    cax=fig.add_axes([.82,.09,.15,.027]);fig.colorbar(im,cax=cax,orientation='horizontal',label='Median promoter ATAC log₂FC (colors clipped at ±1)')
    fig.suptitle('Gene-level RNA responses across batches',x=.025,y=.98,ha='left',fontweight='bold',fontsize=17)
    fig.text(.025,.035,'All eligible members of each pathway; one point per gene. Effects: LOY / CTRL.\nOpposite RNA effect signs are descriptive; no batch-by-treatment interaction has been tested.',fontsize=9)
    save(fig,'RNA_batch_gene_effects')


def draw_examples(key,pathway,genes,meta):
    ids=pathway['examples'];names=[genes[g]['symbol'] for g in ids]
    effects=np.array([[genes[g]['atac_median_logFC']]+[r['logFC'] for r in genes[g]['rna']] for g in ids])
    expr=np.array([np.log2(np.array(genes[g]['norm'])+1) for g in ids]);expr-=expr.mean(axis=1,keepdims=True)
    fig,axes=plt.subplots(1,2,figsize=(14,6.5),gridspec_kw={'width_ratios':[1,2.4]});fig.subplots_adjust(left=.12,right=.95,bottom=.20,top=.82,wspace=.55)
    for ax,values in zip(axes,[effects,expr]):
        bound=max(.5,math.ceil(float(abs(values).max())*2)/2)
        im=ax.imshow(values,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-bound,vcenter=0,vmax=bound),aspect='auto')
        ax.set_yticks(range(len(ids)),names,fontsize=8);ax.tick_params(length=0);fig.colorbar(im,ax=ax,shrink=.8,pad=.02)
    axes[0].set_xticks(range(4),['ATAC\nmedian','RNA\nCombined','RNA\nAugust','RNA\nMay'],fontsize=8);axes[0].set_title('Original log₂ fold changes',fontsize=11)
    for i,g in enumerate(ids):
        for j,r in enumerate(genes[g]['rna']):
            if r['FDR']<.05:axes[0].text(j+1,i,'*',ha='center',va='center',fontsize=12)
    axes[1].set_xticks(range(12),[f'{"Aug" if m.batch=="AUG_2025" else "May"}\n{m.treatment}\n{str(m.sample_id).split("-")[-1]}' for m in meta.itertuples()],fontsize=7.5)
    for x in [2.5,5.5,8.5]:axes[1].axvline(x,color='white',lw=3)
    axes[1].set_title('RNA by sample: mean-centered log₂(normalized count + 1)',fontsize=10)
    fig.suptitle(TITLES[key]+': selected leading-edge genes',x=.025,y=.98,ha='left',fontweight='bold',fontsize=17)
    fig.text(.025,.035,'Examples: top three leading-edge contributors per assay, deduplicated and filled to 12 by maximum contribution.\nLeft: actual log₂FC; * RNA gene FDR < 0.05. ATAC median has no gene-level FDR. Right: within-gene deviations without SD scaling or batch correction.\nAll leading-edge genes and individual-sample dots are available in the browser report.',fontsize=9)
    save(fig,key+'_gene_effects_expression')
    with PdfPages(FIG/(key+'_sample_dots.pdf')) as book:
        for start in range(0,len(ids),4):
            subset=ids[start:start+4];fig,axes=plt.subplots(len(subset),3,figsize=(14,2.2*len(subset)+1.5),squeeze=False)
            fig.subplots_adjust(left=.07,right=.98,bottom=.09,top=.90,hspace=.7,wspace=.4)
            for i,gid in enumerate(subset):
                g=genes[gid];ax=axes[i,0]
                for j,r in enumerate(g['rna']):ax.errorbar(r['logFC'],j,xerr=1.96*r['SE'],fmt='o',color=PURPLE if r['FDR']<.05 else '#728391',markersize=4,capsize=2)
                ax.axvline(0,color='#a5afb7',lw=.7);ax.set_yticks(range(3),['Combined','August','May'],fontsize=8);ax.invert_yaxis();ax.set_xlabel('RNA log₂FC (95% Wald CI)',fontsize=8);ax.set_title(g['symbol'],loc='left',fontweight='bold',fontsize=11)
                ax=axes[i,1];values=np.log2(np.array(g['norm'])+1)
                for j,(batch,condition) in enumerate([('AUG_2025','CTRL'),('AUG_2025','LOY'),('MAY_2026','CTRL'),('MAY_2026','LOY')]):
                    ix=np.flatnonzero((meta.batch==batch)&(meta.treatment==condition));v=values[ix];color=BLUE if condition=='CTRL' else ORANGE
                    ax.scatter(j+np.array([-.13,0,.13]),v,c=color,s=22);ax.hlines(v.mean(),j-.22,j+.22,color=color,lw=1.5)
                ax.set_xticks(range(4),['Aug\nCTRL','Aug\nLOY','May\nCTRL','May\nLOY'],fontsize=8);ax.set_ylabel('RNA log₂(norm. + 1)',fontsize=8)
                ax=axes[i,2];p=next(p for p in g['peaks'] if p['id']==g['representative_peak'])
                for j in range(2):
                    v=p['cpm'][j*3:j*3+3];color=BLUE if j==0 else ORANGE;ax.scatter(j+np.array([-.13,0,.13]),v,c=color,s=22);ax.hlines(np.mean(v),j-.22,j+.22,color=color,lw=1.5)
                ax.set_xticks([0,1],['CTRL','LOY'],fontsize=8);ax.set_ylim(bottom=0);ax.set_ylabel('ATAC peak TMM CPM',fontsize=8)
                ax.set_title(f'Median ATAC log₂FC {g["atac_median_logFC"]:.2f}\nShown peak log₂FC {p["logFC"]:.2f}; FDR {p["FDR"]:.3g}',fontsize=8)
            fig.suptitle(TITLES[key]+': individual samples',x=.025,y=.98,ha='left',fontweight='bold',fontsize=16)
            fig.text(.025,.025,'Each dot is a sample; horizontal bars show means. RNA shares combined-model size factors. ATAC shows the abundance-selected promoter peak.\nRNA and ATAC samples are not paired across assays. RNA interval colors indicate source FDR < 0.05 (purple).',fontsize=8)
            book.savefig(fig);plt.close(fig)


def write_report(summary,concentration,pathways,methods):
    table=lambda f:f.to_html(index=False,border=0,float_format=lambda v:f'{v:.3g}',escape=True)
    intro='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>E2F, MYC and interferon gene inspection</title><script src="interactive_assets/plotly.min.js"></script><script src="pathway_driver_assets/data.js"></script><style>body{font:15px/1.45 Arial,sans-serif;color:#233746;max-width:1500px;margin:25px auto;padding:0 25px}a{color:#0072b2}h1,h2{color:#15384e}.note{background:#edf3f7;padding:16px}img{width:100%;height:auto}figure{margin:25px 0}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:7px;border-bottom:1px solid #dae3e9;text-align:left}.scroll{overflow:auto;max-height:420px}.plots{display:grid;grid-template-columns:1fr 1fr;gap:18px}.plot{height:350px}.controls{display:flex;gap:15px;flex-wrap:wrap;align-items:center;margin:20px 0}select,input,button{font:inherit;padding:8px}tr[data-gene]{cursor:pointer}tr.selected{background:#f0e6f5}pre{white-space:pre-wrap}.small{font-size:12px;color:#526674}@media(max-width:900px){.plots{grid-template-columns:1fr}}</style><body><p><a href="pathway_GSEA.html">Pathway GSEA report</a> · <a href="ribosomal_GSEA_sensitivity.html">Ribosomal exclusion sensitivity</a></p><h1>Genes driving E2F, MYC and interferon enrichment</h1><p class="note">Explore all 379 distinct leading-edge genes from the original ATAC and three RNA rankings. Leading edges from nonsignificant pathways are descriptive. The figures below use 12 reproducibly selected examples per pathway; full pathway membership is retained in the direction summaries.</p>
<h2>All pathway genes: August versus May RNA effects</h2><a href="../integration/pathway_GSEA/gene_drivers/figures/RNA_batch_gene_effects.pdf"><img src="../integration/pathway_GSEA/gene_drivers/figures/RNA_batch_gene_effects.png" alt="RNA effect comparison for all pathway members"></a>
<h2>Inspect a leading-edge gene</h2><div class="controls"><label>Pathway <select id="driver-pathway"></select></label><label>Gene <select id="driver-gene"></select></label><input id="driver-search" placeholder="Gene symbol or Entrez ID"><button id="driver-find">Find gene</button></div><p id="driver-status" class="small"></p><div id="driver-pathway-stats"></div><h3 id="driver-heading"></h3><p id="driver-info" class="small"></p><div class="plots"><div id="driver-effects" class="plot"></div><div id="driver-expression" class="plot"></div><div id="driver-promoters" class="plot"></div><div><label>Promoter peak <select id="driver-peak"></select></label><div id="driver-atac-counts" class="plot"></div></div></div><p class="small">RNA: original effects and 95% Wald intervals; sample expression uses shared combined-model size factors and is not batch-corrected. ATAC: every linked promoter effect, with per-sample counts for the selected peak. No gene-level ATAC FDR is assigned.</p><div class="scroll" id="driver-table"></div>
'''
    findings='''<section id="gene-level-findings"><h2>What the genes show</h2>
<p><strong>E2F is a broad accessibility pattern, with RNA support that differs between batches.</strong> Positive median promoter effects occur in 147/173 eligible E2F genes. RNA effects are positive in 125/173 genes in August, compared with 55/173 in May. Most effects are modest: the median RNA log2FC is +0.118 in August and −0.073 in May. The top five genes contribute only 9.6% of ATAC leading-edge weight, so that leading edge is spread across many genes. These direction counts do not imply individual gene significance or a formally significant batch-by-treatment interaction.</p>
<p><strong>MYC target accessibility does not establish MYC transcriptional activation.</strong> Median promoter effects are positive in 137/163 MYC V1 genes and 43/49 MYC V2 genes, while RNA effects are generally near zero or negative. MYC itself has negative RNA coefficients in both batches, neither significant at FDR 0.05. Target-set enrichment therefore should not be interpreted as increased MYC expression or direct binding.</p>
<p><strong>Interferon has clearer examples of reproducible RNA decreases, with heterogeneous ATAC effects.</strong> SAMD9L and TNFAIP2 have negative RNA effects with FDR &lt;0.05 in both batches, and negative median promoter ATAC effects. All three LOY RNA samples lie below all three CTRL samples within each batch for both genes. MX1 also has this sample separation, but only May reaches RNA FDR &lt;0.05. BST2 decreases in both RNA batches but reaches FDR &lt;0.05 only in August. These examples do not imply significant ATAC changes at every linked peak.</p>
<p>Across the complete interferon gene sets, ATAC effects are mixed: 33/63 alpha-response genes and 71/118 gamma-response genes have positive promoter medians. A negative GSEA score describes position within the full ranked universe; it does not mean every pathway gene has a negative effect. Interferon ATAC leading edges are more concentrated than E2F: their top five genes account for 46% (alpha) and 34% (gamma) of leading-edge weight. This weight summary is descriptive, not a leave-gene-out enrichment test.</p></section>'''
    intro=intro.replace('<h2>All pathway genes:',findings+'<h2>All pathway genes:')
    for key,p in pathways.items():
        stem=key+'_gene_effects_expression';intro+=f'<h2>{escape(TITLES[key])}: example genes</h2><figure><a href="../integration/pathway_GSEA/gene_drivers/figures/{stem}.pdf"><img src="../integration/pathway_GSEA/gene_drivers/figures/{stem}.png" alt="{escape(TITLES[key])} effects and RNA expression"></a></figure><p><a href="../integration/pathway_GSEA/gene_drivers/figures/{stem}.pdf">Vector heatmap PDF</a> · <a href="../integration/pathway_GSEA/gene_drivers/figures/{key}_sample_dots.pdf">Individual-sample plots and RNA effect intervals</a></p>'
    intro+='<h2>Direction across all eligible pathway genes</h2><div class="scroll">'+table(summary)+'</div><h2>Is the leading edge dominated by a few genes?</h2><p>The following weight summaries are descriptive. They are based on absolute ranking statistics within each leading edge and do not constitute a new enrichment test.</p><div class="scroll">'+table(concentration)+'</div>'
    intro+='<h2>Data and methods</h2><p>'+ ' · '.join(f'<a href="../integration/pathway_GSEA/gene_drivers/{p}">{label}</a>' for p,label in [('all_pathway_gene_effects.csv','All pathway genes and effects'),('RNA_sample_expression.csv','RNA sample counts'),('ATAC_promoter_peak_measurements.csv','ATAC sample counts and peaks'),('RNA_sample_robustness.csv','Within-batch sample robustness'),('selected_example_genes.csv','Example selection audit'),('METHODS.txt','Methods')])+'</p><details><summary>Complete methods</summary><pre>'+escape(methods)+'</pre></details><script src="pathway_driver_assets/viewer.js"></script></body></html>'
    if (ROOT/'report/RNA_interaction_CN.html').exists():
        intro=intro.replace('</body>','<section id="rna-interaction-cn"><h2>Formal RNA batch interaction and D0 copy-number context</h2><p><a href="RNA_interaction_CN.html">Interaction tests, correlation sensitivity, hg38 WGS profiles and joint gene tables</a>.</p></section></body>')
    (ROOT/'report/pathway_gene_drivers.html').write_text(intro,encoding='utf-8')
    section='<section id="pathway-gene-drivers"><h2>Inspect the genes driving enrichment</h2><p><a href="pathway_gene_drivers.html">E2F, MYC and interferon: gene effects, individual RNA and ATAC samples, and interactive promoter peak inspection</a>.</p></section>'
    for name in ['pathway_GSEA.html','ATACseq_report.html']:
        path=ROOT/'report'/name
        if path.exists():
            text=path.read_text(encoding='utf-8')
            if 'id="pathway-gene-drivers"' not in text:path.write_text(text.replace('</body>',section+'</body>'),encoding='utf-8')


if __name__=='__main__':main()
