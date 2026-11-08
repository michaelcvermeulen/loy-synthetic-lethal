"""Join newly fitted RNA interaction results to hg38 D0 dosage context."""
from pathlib import Path
from html import escape
import json,hashlib,math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'integration/RNA_interaction_CN';FIG=OUT/'figures'
FIG.mkdir(parents=True,exist_ok=True)
PATHS={'E2F':'HALLMARK_E2F_TARGETS','MYC V1':'HALLMARK_MYC_TARGETS_V1','MYC V2':'HALLMARK_MYC_TARGETS_V2','IFN alpha':'HALLMARK_INTERFERON_ALPHA_RESPONSE','IFN gamma':'HALLMARK_INTERFERON_GAMMA_RESPONSE'}
FOCUS=['SAMD9L','TNFAIP2','MX1','BST2','TFAP2C','MYC','POLD1','CHEK2']
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
BLUE='#0072b2';ORANGE='#d55e00';PURPLE='#7b3294'
def save(fig,name):
    for ext in ['pdf','png','svg']:fig.savefig(FIG/f'{name}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def rna_figures(rna):
    d=rna[np.isfinite(rna.interaction_padj)].copy();sig=d.interaction_padj<.05
    fig,ax=plt.subplots(figsize=(9,6.5));fig.subplots_adjust(bottom=.2,top=.86)
    ax.scatter(d.interaction_log2FoldChange,-np.log10(d.interaction_padj.clip(lower=1e-300)),s=8,c=np.where(sig,np.where(d.interaction_log2FoldChange>0,ORANGE,BLUE),'#c4ccd2'),alpha=.65,linewidths=0)
    ax.axhline(-np.log10(.05),ls='--',lw=.8,c='#6a7882');ax.axvline(0,lw=.7,c='#adb7bf')
    labels=pd.concat([d.nsmallest(8,'interaction_padj'),d[d.gene_name.isin(['TFAP2C','LENG1'])]]).drop_duplicates('gene_id')
    texts=[ax.text(r.interaction_log2FoldChange,-np.log10(r.interaction_padj),r.gene_name,fontsize=8) for r in labels.itertuples()]
    adjust_text(texts,ax=ax,arrowprops={'arrowstyle':'-','color':'#687985','lw':.6},ensure_inside_axes=True)
    ax.set(xlabel='Interaction log₂FC: May (LOY / CTRL) − August (LOY / CTRL)',ylabel='−log₁₀ interaction FDR')
    fig.suptitle(f'RNA batch × treatment interaction: {sig.sum()} Wald hits at FDR < 0.05',fontweight='bold',fontsize=14)
    fig.text(.02,.02,'DESeq2 on raw counts; 3 CTRL + 3 LOY per batch. BH across all finite gene interaction tests.\nPositive: LOY effect is more positive in May. CTRL versus sgCTRL context differs between batches.\nLENG1 has an unstable effect magnitude due to a zero-count group; see the likelihood-ratio sensitivity in the report.',fontsize=9)
    save(fig,'RNA_interaction_volcano')
    genes=FOCUS+rna.nsmallest(4,'interaction_padj').gene_name.tolist()
    rows=rna.set_index('gene_name').loc[genes]
    fig,axes=plt.subplots(1,2,figsize=(12,7.5),sharey=True);fig.subplots_adjust(left=.14,right=.90,bottom=.16,top=.85,wspace=.28)
    for i,(name,row) in enumerate(rows.iterrows()):
        for j,(batch,color) in enumerate([('August',BLUE),('May',ORANGE)]):
            axes[0].errorbar(row[batch+'_log2FoldChange'],i+(-.13 if j==0 else .13),xerr=1.96*row[batch+'_lfcSE'],fmt='o',c=color,ms=4,capsize=2,label=batch if i==0 else None)
        axes[1].errorbar(row.interaction_log2FoldChange,i,xerr=1.96*row.interaction_lfcSE,fmt='o',c=PURPLE if row.interaction_padj<.05 else '#7b8c99',ms=5,capsize=2)
        axes[1].text(1.02,i,f'{row.interaction_padj:.3g}',transform=axes[1].get_yaxis_transform(),va='center',fontsize=8)
    for ax in axes:ax.axvline(0,c='#aab5bd',lw=.8);ax.set_ylim(len(rows)-.5,-.5)
    axes[0].set_yticks(range(len(rows)),genes);axes[0].set_xlabel('RNA log₂FC (LOY / CTRL)');axes[0].legend(loc='upper left',bbox_to_anchor=(0,1.13),ncol=2,frameon=False)
    axes[1].set_xlabel('Interaction log₂FC (May − August)');axes[1].set_title('Difference between responses',fontsize=11)
    axes[1].text(1.02,1.025,'FDR',transform=axes[1].transAxes,fontsize=9)
    fig.suptitle('Key genes: batch effects and formal interaction',fontweight='bold',fontsize=16)
    fig.text(.02,.025,'Points and 95% Wald intervals from the new joint interaction model; previous separate-fit RNA results are unchanged.\nExamples: eight previously discussed genes plus the four lowest interaction FDRs. Nonsignificance does not establish equivalent effects.',fontsize=9)
    save(fig,'RNA_key_gene_interactions')

def interval_context(chrom,start,end,bychrom):
    """Length-weighted context, without filling filtered bins or assembly gaps."""
    result={'covered_fraction':0.,'delta_mean':np.nan,'delta_min':np.nan,'delta_max':np.nan,'delta_max_abs':np.nan,'bins':0}
    if chrom not in bychrom or pd.isna(start) or pd.isna(end) or end<=start:return result
    b=bychrom[chrom];i=np.searchsorted(b['end0'],start,side='right');j=np.searchsorted(b['start0'],end,side='left')
    if j<=i:return result
    w=np.minimum(b['end0'][i:j],end)-np.maximum(b['start0'][i:j],start);v=b['delta_segment_log2'][i:j]
    good=(w>0)&b['use'][i:j]&np.isfinite(v);w=w[good];v=v[good]
    if not len(v):return result
    result.update(covered_fraction=float(w.sum()/(end-start)),delta_mean=float(np.average(v,weights=w)),delta_min=float(v.min()),delta_max=float(v.max()),delta_max_abs=float(np.max(abs(v))),bins=len(v))
    return result

def cn_join(rna,bins):
    bychrom={c:{k:f[k].to_numpy() for k in ['start0','end0','delta_segment_log2','use']} for c,f in bins.groupby('chrom',sort=False)}
    peaks=pd.read_csv(ROOT/'differential/primary_TMM_all.csv')
    links=pd.read_csv(ROOT/'integration/feature_overlaps/peak_gene_link_audit.csv',dtype={'entrez_id':str})
    links=links[links.feature.eq('promoter')][['entrez_id','peak_id']].drop_duplicates()
    p=peaks[peaks.peak_id.isin(links.peak_id)].copy()
    contexts=[interval_context(r.chrom,r.start,r.end,bychrom) for r in p.itertuples()]
    p=pd.concat([p.reset_index(drop=True),pd.DataFrame(contexts).add_prefix('CN_')],axis=1)
    p.to_csv(OUT/'ATAC_promoter_peak_CN_context.csv',index=False)
    linked=links.merge(p[['peak_id','CN_covered_fraction','CN_delta_mean','CN_delta_min','CN_delta_max','CN_delta_max_abs']],on='peak_id',validate='many_to_one')
    linked.to_csv(OUT/'promoter_gene_peak_CN_links.csv',index=False)
    group=linked.groupby('entrez_id')
    pg=group.agg(CN_promoter_peak_count=('peak_id','size'),CN_promoter_min_coverage=('CN_covered_fraction','min'),CN_promoter_delta_median=('CN_delta_mean','median'),CN_promoter_delta_min=('CN_delta_min','min'),CN_promoter_delta_max=('CN_delta_max','max'),CN_promoter_max_abs_delta=('CN_delta_max_abs','max')).reset_index()
    for threshold in [.1,.2,.3]:
        pg[f'CN_flag_{threshold:.2f}']=np.where(pg.CN_promoter_max_abs_delta>=threshold,'D0 dosage difference',np.where(pg.CN_promoter_min_coverage>=.9,'No large D0 difference','Insufficient coverage'))
    base=pd.read_csv(ROOT/'integration/feature_scatters/feature_gene_effects.csv',dtype={'entrez_id':str})
    base=base[base.feature.eq('promoter')].copy();assert len(base)==11898
    joint=base.merge(pg,on='entrez_id',how='left',validate='one_to_one')
    rnau=rna[~rna.entrezid.duplicated(False)].copy().rename(columns={'entrezid':'entrez_id'})
    joint=joint.merge(rnau,on='entrez_id',how='left',validate='one_to_one',suffixes=('','_interaction'))
    genes=pd.read_csv(ROOT/'config/genes_knownGene.tsv',sep='\t',dtype={'gene_id':str}).set_index('gene_id')
    bodies=[]
    for gid in joint.entrez_id:
        if gid in genes.index:
            g=genes.loc[gid];assert isinstance(g,pd.Series)
            c=interval_context(g.chrom,g.start-1,g.end,bychrom)
            c.update(chrom=g.chrom,start0=int(g.start-1),end0=int(g.end),tss0=int(g.tss-1))
        else:c={}
        bodies.append(c)
    joint=pd.concat([joint,pd.DataFrame(bodies).add_prefix('gene_body_')],axis=1)
    joint['Both_ATAC010_RNA005']=(joint.minimum_peak_FDR<.1)&(joint.rna_FDR<.05)
    joint['down_down']=joint.Both_ATAC010_RNA005&(joint.atac_median_logFC<0)&(joint.rna_logFC<0)
    joint['up_up']=joint.Both_ATAC010_RNA005&(joint.atac_median_logFC>0)&(joint.rna_logFC>0)
    exact=pd.read_csv(ROOT/'integration/down70_ATAC010_RNA005/exact_70_genes.csv',dtype={'entrez_id':str})
    assert joint.Both_ATAC010_RNA005.sum()==82 and set(joint.loc[joint.down_down,'entrez_id'])==set(exact.entrez_id)
    joint.to_csv(OUT/'ATAC_RNA_CN_joint_gene_table.csv',index=False)
    joint[joint.Both_ATAC010_RNA005].to_csv(OUT/'overlap82_CN_context.csv',index=False)
    summary=[]
    membership=pd.read_csv(ROOT/'integration/pathway_GSEA/sources/gene_set_membership.csv',dtype={'entrez_id':str})
    common=set(pd.read_csv(ROOT/'integration/pathway_GSEA/gene_ranks.csv',dtype={'entrez_id':str}).entrez_id)
    groups={'All promoter-eligible':joint,'Overlap 82':joint[joint.Both_ATAC010_RNA005],'Down/down 70':joint[joint.down_down],'Up/up 9':joint[joint.up_up],
      **{k:joint[joint.entrez_id.isin(set(membership.loc[membership.pathway.eq(v),'entrez_id'])&common)] for k,v in PATHS.items()}}
    for name,g in groups.items():
        for threshold in [.1,.2,.3]:
            flag=f'CN_flag_{threshold:.2f}';stable=g[g[flag].eq('No large D0 difference')]
            summary.append({'group':name,'threshold':threshold,'genes':len(g),'D0_dosage_difference':int(g[flag].eq('D0 dosage difference').sum()),'no_large_D0_difference':len(stable),'insufficient_coverage':int(g[flag].eq('Insufficient coverage').sum()),'ATAC_positive_in_no_large_difference':int(stable.atac_median_logFC.gt(0).sum()),'median_ATAC_in_no_large_difference':stable.atac_median_logFC.median(),'median_RNA_in_no_large_difference':stable.rna_logFC.median()})
    summary=pd.DataFrame(summary);summary.to_csv(OUT/'CN_stratification_summary.csv',index=False)
    return joint,summary

def cn_figures(bins,joint,summary):
    use=bins[bins.use&np.isfinite(bins.delta_segment_log2)].copy()
    chromosomes=[f'chr{i}' for i in range(1,23)];sizes=bins.groupby('chrom').end0.max();offset={};n=0
    for c in chromosomes:offset[c]=n;n+=sizes[c]+3000000
    fig,axes=plt.subplots(3,1,figsize=(15,8),sharex=True);fig.subplots_adjust(bottom=.17,top=.87,hspace=.3)
    for c in chromosomes:
        d=use[use.chrom.eq(c)];x=(d.start0+d.end0)/2+offset[c]
        for ax,col in zip(axes,['CTRL_segment_log2','LOY_segment_log2','delta_segment_log2']):
            ax.scatter(x/1e6,d[col],s=2,c=BLUE if int(c[3:])%2 else '#7c9cac',alpha=.75,linewidths=0)
            ax.axvline(offset[c]/1e6,c='#dde3e8',lw=.5)
    for ax in axes:ax.axhline(0,c='#a8b4bd',lw=.7);ax.set_ylabel('Relative log₂ ratio')
    axes[0].set_title('CTRL D0',loc='left');axes[1].set_title('LOY D0',loc='left');axes[2].set_title('LOY − CTRL D0',loc='left')
    axes[2].axhline(.2,c=ORANGE,ls='--',lw=.7);axes[2].axhline(-.2,c=ORANGE,ls='--',lw=.7)
    axes[2].set_xticks([(offset[c]+sizes[c]/2)/1e6 for c in chromosomes],[c[3:] for c in chromosomes]);axes[2].set_xlabel('Autosomal chromosome (hg38)')
    fig.suptitle('D0 WGS copy-number context: 100 kb QDNAseq',fontweight='bold',fontsize=16)
    fig.text(.02,.02,'GC/mappability-corrected, segmented relative depth; one WGS sample per condition. Dashed lines: descriptive |difference| = 0.20.\nFiltered bins remain missing. Median normalization does not measure absolute ploidy. Matching to later RNA/ATAC cultures is unconfirmed.',fontsize=9)
    save(fig,'WGS_genome_dosage_context')
    s=summary[summary.threshold.eq(.2)].iloc[1:]
    fig,ax=plt.subplots(figsize=(11,6.5));fig.subplots_adjust(left=.18,bottom=.26,top=.9);left=np.zeros(len(s))
    for col,color,label in [('no_large_D0_difference','#88b6a0','No large D0 difference'),('D0_dosage_difference',ORANGE,'D0 dosage difference'),('insufficient_coverage','#c9cfd3','Insufficient coverage')]:
        ax.barh(range(len(s)),s[col],left=left,color=color,label=label);left+=s[col].to_numpy()
    ax.set_yticks(range(len(s)),s.group);ax.invert_yaxis();ax.set_xlabel('Genes');ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.17),ncol=3,fontsize=8)
    ax.set_title('Promoter genes stratified by D0 copy-number context',loc='left',fontweight='bold')
    fig.text(.02,.02,'Dosage flag: any covered linked promoter peak overlaps a bin with |segmented LOY−CTRL log₂ ratio| ≥ 0.20.\n“No large difference” requires ≥90% usable coverage of every linked promoter peak; it does not establish copy-number stability in later cultures.',fontsize=8)
    save(fig,'CN_pathway_overlap_stratification')
    focus=joint[joint.rna_symbol.isin(FOCUS)].copy().set_index('rna_symbol').reindex(FOCUS).reset_index()
    focus.to_csv(OUT/'key_gene_CN_and_interaction.csv',index=False)
    fig,axes=plt.subplots(4,2,figsize=(13,10));fig.subplots_adjust(top=.91,bottom=.10,hspace=.7,wspace=.26)
    for ax,r in zip(axes.flat,focus.itertuples()):
        if pd.isna(r.gene_body_start0):ax.text(.5,.5,'No mapped gene coordinates',ha='center');continue
        c=r.gene_body_chrom;start=max(0,r.gene_body_start0-500000);end=r.gene_body_end0+500000;d=use[use.chrom.eq(c)&(use.end0>start)&(use.start0<end)]
        x=(d.start0+d.end0)/2/1e6
        ax.scatter(x,d.delta_bin_log2,s=12,c='#98aab6',label='Normalized bin')
        ax.scatter(x,d.delta_segment_log2,s=12,c=PURPLE,label='Segmented bin')
        ax.axvspan(r.gene_body_start0/1e6,r.gene_body_end0/1e6,color='#e4b65c',alpha=.25)
        for y in [-.2,.2]:ax.axhline(y,c='#adb7bf',lw=.6,ls='--')
        flag=focus.loc[focus.rna_symbol.eq(r.rna_symbol),'CN_flag_0.20'].iloc[0]
        ax.axhline(0,c='#b7c1c8',lw=.6);ax.set_title(r.rna_symbol+' · '+flag,loc='left',fontsize=10)
        ax.set_xlabel(c+' position (Mb)',fontsize=8);ax.set_ylabel('LOY − CTRL log₂ ratio',fontsize=8)
        ax.set_ylim(min(-.5,d.delta_bin_log2.min()-.1) if len(d) else -.5,max(.5,d.delta_bin_log2.max()+.1) if len(d) else .5)
    axes[0,0].legend(frameon=False,fontsize=8)
    fig.suptitle('D0 WGS around key genes',fontweight='bold',fontsize=16)
    fig.text(.02,.02,'Shading: gene body. Each point is a 100 kb bin; excluded bins are omitted. Values are relative dosage context, not ATAC coverage.\nLocal patterns cannot establish copy number in the later assay cultures or exclude sub-bin focal changes.',fontsize=9)
    save(fig,'WGS_key_gene_regions')

def write_report(rna,joint,summary):
    def table(d):return '<div class="scroll">'+d.to_html(index=False,border=0,float_format=lambda x:f'{x:.3g}',escape=True)+'</div>'
    def figure(stem,caption):return f'<figure><a href="../integration/RNA_interaction_CN/figures/{stem}.pdf"><img src="../integration/RNA_interaction_CN/figures/{stem}.png" alt="{escape(caption)}"></a><figcaption>{escape(caption)} <a href="../integration/RNA_interaction_CN/figures/{stem}.pdf">Vector PDF</a></figcaption></figure>'
    model=json.loads((OUT/'RNA_interaction_model_summary.json').read_text())
    cn=json.loads((OUT/'WGS_model_summary.json').read_text())
    fixed=pd.read_csv(OUT/'RNA_interaction_CAMERA_fixed_correlation.csv').rename(columns={'Direction':'Direction_fixed','PValue':'P_fixed','FDR':'FDR_fixed'})
    estimated=pd.read_csv(OUT/'RNA_interaction_CAMERA_estimated_correlation.csv').rename(columns={'Direction':'Direction_estimated','PValue':'P_estimated','FDR':'FDR_estimated'})
    cam=fixed.merge(estimated.drop(columns='NGenes'),on='pathway',validate='one_to_one')
    cam.to_csv(OUT/'RNA_interaction_CAMERA_comparison.csv',index=False)
    key=joint[joint.rna_symbol.isin(FOCUS)].copy().set_index('rna_symbol').reindex(FOCUS).reset_index()
    rows=summary[(summary.threshold==.2)].set_index('group')
    changed=int(rows.loc['Down/down 70','D0_dosage_difference']);stable=int(rows.loc['Down/down 70','no_large_D0_difference']);unknown=int(rows.loc['Down/down 70','insufficient_coverage'])
    html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>RNA interaction and D0 copy-number context</title><style>body{font:15px/1.5 Arial,sans-serif;color:#243745;max-width:1400px;margin:25px auto;padding:0 25px}h1,h2{color:#15384e}a{color:#0072b2}.note{padding:18px;background:#edf3f7}figure{margin:25px 0}img{width:100%;height:auto}table{border-collapse:collapse;font-size:12px;width:100%}td,th{text-align:left;padding:8px;border-bottom:1px solid #d8e1e8}.scroll{overflow:auto;max-height:600px}pre{white-space:pre-wrap}</style></head><body><p><a href="pathway_gene_drivers.html">Gene inspection</a> · <a href="pathway_GSEA.html">Original pathway GSEA</a> · <a href="ATACseq_report.html">Main ATAC report</a></p><h1>RNA batch interaction and D0 copy-number context</h1>'''
    html+=f'<div class="note"><strong>The RNA analysis finds {model["interaction_FDR005"]} Wald interaction hits; 52 also pass the likelihood-ratio sensitivity (FDR &lt;0.05 in both).</strong> LENG1 is the unsupported Wald hit and has an unstable sparse-count effect estimate. The broad E2F/MYC pathway difference is sensitive to the assumed correlation between genes. The hg38 WGS assessment supplies D0 line-level context, with one sample per condition; it does not correct or establish copy number in the later RNA/ATAC cultures.</div>'
    html+='<h2>What this changes in the interpretation</h2><p><strong>Batch dependence is supported for a subset of genes.</strong> The interaction compares May LOY-versus-CTRL log2FC minus August LOY-versus-CTRL log2FC. The sample design is balanced (three per condition per batch), but August CTRL and May sgCTRL differ in experimental context. These data cannot attribute an interaction specifically to passage, technical batch, or control manipulation.</p>'
    html+='<p><strong>SAMD9L and TNFAIP2 remain reproducibly decreased in RNA.</strong> Both have significant negative within-batch contrasts in the new interaction fit, without a significant interaction. This supports their direction across these experiments; a nonsignificant interaction does not establish equivalent effect sizes. TFAP2C has a more negative response in May, but its interaction FDR is 0.0524, just above the prespecified 0.05 cutoff.</p>'
    html+='<p><strong>The E2F/MYC pathway difference is not robust to correlation handling.</strong> CAMERA gives FDR &lt;0.05 for the three target sets when using its default correlation of 0.01. Estimated within-set residual correlations are about 0.49–0.65, and none of the five pathway tests is significant with that adjustment. We therefore retain the observed direction shift as descriptive and avoid claiming a firmly established pathway-wide interaction. Interferon sets have no significant interaction under either setting.</p>'
    html+=f'<p><strong>D0 dosage screening of the 70 down/down genes:</strong> {changed} overlap a promoter bin with |relative dosage difference| ≥0.20, {stable} have no difference this large with sufficient promoter coverage, and {unknown} have insufficient coverage. These categories use 100 kb bins and do not establish later-culture copy-number stability. Results at thresholds 0.10 and 0.30 are also provided.</p>'
    html+='<p><strong>The measured D0 dosage does not explain the four highlighted interferon examples.</strong> SAMD9L, TNFAIP2, MX1 and BST2 all have complete linked-promoter coverage and small D0 relative dosage differences (absolute log2 ratio below 0.03). MYC and CHEK2, in contrast, fall in regions with lower relative dosage in LOY D0 (about −0.21 and −0.26). These findings qualify gene-level interpretation without proving a regulatory mechanism.</p>'
    html+='<p><strong>Chromosome 19 does not show a large broad D0 loss.</strong> Its median segmented LOY−CTRL log2 ratio is −0.022. Among the 19 chromosome-19 genes in the down/down set, 13 have adequate promoter coverage without a large difference and six lack adequate coverage. This argues against a broad D0 chromosome-19 loss explaining the observed cluster, while leaving focal and later-culture changes unresolved.</p>'
    html+='<h2>Formal RNA interaction</h2>'+figure('RNA_interaction_volcano','All finite gene interaction tests; each gene is one point.')+figure('RNA_key_gene_interactions','New joint-model within-batch effects, interaction differences and 95% Wald intervals.')
    html+=f'<p>{model["tested_genes"]:,} genes passed the annotation/count filter; {model["finite_interaction_p"]:,} have finite interaction P-values. {model["interaction_FDR005"]} pass FDR 0.05 and {model["interaction_FDR010"]} pass FDR 0.10. No coefficient fits failed to converge. Nine genes lack a finite interaction P-value under the default Cook outlier handling.</p>'
    lr=pd.read_csv(OUT/'RNA_interaction_LRT_sensitivity.csv')
    both=((lr.Wald_FDR<.05)&(lr.LRT_FDR<.05)).sum()
    leng=lr[lr.symbol.eq('LENG1')].iloc[0]
    annotated=rna.merge(lr.drop(columns='symbol'),on='gene_id',validate='one_to_one')
    annotated.to_csv(OUT/'RNA_interaction_results_with_validation.csv',index=False)
    annotated[(annotated.Wald_FDR<.05)&(annotated.LRT_FDR<.05)].sort_values('Wald_FDR').to_csv(OUT/'RNA_interaction_supported_both_tests.csv',index=False)
    html+=f'<p><strong>Sparse-count check:</strong> reversing the batch reference flags LENG1 as having an unstable Wald effect magnitude because an entire sample group has zero counts. Its Wald FDR is {leng.Wald_FDR:.3g}; likelihood-ratio-test FDR is {leng.LRT_FDR:.3g}. Treat its estimated interaction magnitude cautiously. Across all genes, {both} of the 53 Wald hits also pass LRT FDR 0.05. <a href="../integration/RNA_interaction_CN/RNA_interaction_LRT_sensitivity.csv">Full LRT sensitivity</a> · <a href="../integration/RNA_interaction_CN/RNA_reparameterization_audit.csv">Reference-reversal audit</a>.</p>'
    html+='<p>Prioritize the <a href="../integration/RNA_interaction_CN/RNA_interaction_supported_both_tests.csv">52 interactions supported by both tests</a>; LENG1 is not included in that list. The LRT uses the fitted DESeq2 dispersions and compares the full interaction model with the additive model, ~batch + treatment, with one degree of freedom and genome-wide BH correction. This additional check was triggered by the sparse-count numerical instability.</p>'
    html+='<h2>Pathway interaction sensitivity</h2><p>These five previously selected sets are tested on the RNA count matrix using limma-voom and CAMERA. They include RNA-eligible genes without requiring an ATAC promoter, so set sizes differ from the original common-universe GSEA. Down means a relatively more negative May response; Up means a relatively more positive May response. BH is applied within these five exploratory questions. Fixed and estimated correlation results are both shown.</p>'+table(cam)
    html+='<h2>Key genes across assays</h2><p>Original combined RNA and ATAC effects are retained in the joint table. The August, May and interaction columns below come from the new interaction model. Promoter CN values are summaries of the bins overlapping all linked promoter peaks; they are descriptive relative dosage differences.</p>'
    cols=['rna_symbol','atac_median_logFC','August_log2FoldChange','August_padj','May_log2FoldChange','May_padj','interaction_log2FoldChange','interaction_padj','CN_promoter_delta_median','CN_promoter_min_coverage','CN_flag_0.20']
    html+=table(key[cols])
    html+='<h2>hg38 WGS assessment</h2>'+figure('WGS_genome_dosage_context','GC/mappability-corrected and segmented autosomal relative dosage.')+figure('WGS_key_gene_regions','D0 dosage context around the eight previously discussed genes.')
    html+=f'<p>{cn["usable_autosomal_bins"]:,} of {cn["all_bins"]:,} annotation bins are retained for the autosomal analysis. X/Y and filtered bins are not treated as neutral. We use compatible hg38 BAMs and frozen QDNAseq.hg38 100 kb annotations; the older T2T results are not overlapped.</p>'
    html+='<h2>Overlap and pathway stratification</h2>'+figure('CN_pathway_overlap_stratification','Main CN effect-size threshold: absolute relative log2 difference ≥0.20.')+table(summary)
    html+='<p>A dosage flag means at least one covered linked promoter peak overlaps an altered bin. “No large D0 difference” requires ≥90% coverage of every linked promoter peak and no covered bin reaching the threshold. It is deliberately a relative D0 description: both populations can share an amplification/deletion, median normalization removes global ploidy shifts, and later cultures may differ. A CN overlap is not evidence that dosage caused the RNA or ATAC effect.</p>'
    html+='<h2>Download data and methods</h2><p>'+ ' · '.join(f'<a href="../integration/RNA_interaction_CN/{name}">{label}</a>' for name,label in [('RNA_interaction_all_genes.csv','All gene interaction results'),('RNA_interaction_FDR005.csv','53 significant interactions'),('ATAC_RNA_CN_joint_gene_table.csv','Full joint gene table'),('overlap82_CN_context.csv','82-gene overlap with CN context'),('ATAC_promoter_peak_CN_context.csv','Promoter peak CN context'),('WGS_hg38_100kb_bins.csv','All WGS bins and filter flags'),('WGS_hg38_segments.csv','WGS segments'),('CN_stratification_summary.csv','CN threshold sensitivity'),('RNA_model_sample_metadata.csv','RNA samples'),('RNA_INTERACTION_METHODS.txt','RNA methods'),('WGS_CN_METHODS.txt','WGS methods')])+'</p>'
    html+='<p>Methods references: <a href="https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions">DESeq2 interactions</a>; <a href="https://doi.org/10.1093/nar/gks461">CAMERA</a>; <a href="https://bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf">QDNAseq workflow</a>; <a href="https://doi.org/10.5281/zenodo.17948766">hg38 bin annotation v1.2.0</a>.</p>'
    for name in ['RNA_INTERACTION_METHODS.txt','WGS_CN_METHODS.txt']:html+='<details><summary>'+name+'</summary><pre>'+escape((OUT/name).read_text())+'</pre></details>'
    html+='</body></html>'
    (ROOT/'report/RNA_interaction_CN.html').write_text(html,encoding='utf-8')
    section='<section id="rna-interaction-cn"><h2>Formal RNA batch interaction and D0 copy-number context</h2><p><a href="RNA_interaction_CN.html">Interaction tests, correlation sensitivity, hg38 WGS profiles and joint gene tables</a>.</p></section>'
    for name in ['ATACseq_report.html','pathway_gene_drivers.html','pathway_GSEA.html']:
        path=ROOT/'report'/name;text=path.read_text(encoding='utf-8')
        if 'id="rna-interaction-cn"' not in text:text=text.replace('</body>',section+'</body>')
        if name=='ATACseq_report.html':
            text=text.replace('<p>lpWGS dataset selection and assay matching remain pending. No copy-number adjustment has been applied. Candidate hg38 WGS inputs are recorded in PROCESSING_PLAN.md; T2T WGS coordinates require conversion or recomputation before overlap.</p>','<p>hg38 D0 WGS copy-number context and the formal RNA batch interaction are complete. <a href="RNA_interaction_CN.html">Open the interaction and copy-number report</a>. Direct culture matching remains unconfirmed; no sample-specific copy-number correction has been applied.</p>')
        path.write_text(text,encoding='utf-8')
    return {'status':'built_pending_verification','RNA_model':model,'interaction_genes_supported_both_tests':int(both),'Wald_hit_not_supported_by_LRT':'LENG1','WGS_model':cn,'report':'report/RNA_interaction_CN.html','down70_CN_threshold020':{'dosage_flag':changed,'no_large_difference':stable,'insufficient_coverage':unknown},'copy_number_adjustment_applied':False}

def main():
    rna=pd.read_csv(OUT/'RNA_interaction_all_genes.csv',dtype={'entrezid':str})
    bins=pd.read_csv(OUT/'WGS_hg38_100kb_bins.csv')
    rna_figures(rna)
    joint,summary=cn_join(rna,bins)
    cn_figures(bins,joint,summary)
    completion=write_report(rna,joint,summary)
    completion['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),OUT/'RNA_interaction_all_genes.csv',OUT/'WGS_hg38_100kb_bins.csv',ROOT/'integration/feature_scatters/feature_gene_effects.csv',ROOT/'config/genes_knownGene.tsv']}
    (OUT/'completion.json').write_text(json.dumps(completion,indent=2),encoding='utf-8')
    print(json.dumps(completion,indent=2))
    print(summary[summary.threshold.eq(.2)].to_string(index=False))

if __name__=='__main__':main()
