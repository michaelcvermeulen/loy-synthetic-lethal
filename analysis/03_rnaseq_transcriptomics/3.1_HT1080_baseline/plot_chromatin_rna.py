"""Publication exports and standalone report for the measured-reference RNA extension."""
from pathlib import Path
import argparse
import base64
import hashlib
import html
import json
import string
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'outputs/tables/chromatin_rna'
FIG = ROOT / 'outputs/figures/chromatin_rna'
PROV = ROOT / 'provenance/chromatin_rna'
UP, DOWN, INK = '#D55E00', '#0072B2', '#252525'
LABELS = {
    'H2AK119ub1_PROMOTER': 'H2AK119ub1 promoter',
    'H3K27me3_PROMOTER': 'H3K27me3 promoter',
    'H3K9me3_PROMOTER': 'H3K9me3 promoter',
    'H4K20me3_PROMOTER': 'H4K20me3 promoter',
    'H4K20me1_PROMOTER': 'H4K20me1 promoter',
    'H3K4me1_PROMOTER': 'H3K4me1 promoter',
    'H3K27ac_PROMOTER': 'H3K27ac promoter',
    'H3K4me3_PROMOTER': 'H3K4me3 promoter',
    'LaminB1_PROMOTER': 'Lamin B1 promoter peak',
    'LaminB2_PROMOTER': 'Lamin B2 promoter peak',
    'H3K4me3_H3K27me3_COMARKED_PROMOTER': 'H3K4me3 / H3K27me3 co-marked',
    'H3K27me3_H2AK119ub1_JOINT_PROMOTER': 'H3K27me3 / H2AK119ub1 joint',
    'H3K4me3_H3K27ac_JOINT_PROMOTER': 'H3K4me3 / H3K27ac joint',
    'H3K9me3_H4K20me3_JOINT_PROMOTER': 'H3K9me3 / H4K20me3 joint',
    'LaminB1_LaminB2_JOINT_PROMOTER': 'Lamin B1 / B2 joint peak',
    'H3K27me3_DETECTED_H3K4me3_NOT_DETECTED': 'H3K27me3 detected; H3K4me3 undetected',
    'H3K4me3_DETECTED_H3K27me3_NOT_DETECTED': 'H3K4me3 detected; H3K27me3 undetected',
    'H3K9me3_DETECTED_H3K27me3_NOT_DETECTED': 'H3K9me3 detected; H3K27me3 undetected',
    'H3K27me3_GENE_BODY_PEAK_COVERAGE_GE50PCT': 'H3K27me3 body coverage ≥50%',
    'H3K9me3_GENE_BODY_PEAK_COVERAGE_GE50PCT': 'H3K9me3 body coverage ≥50%',
}


def qtext(q):
    if not np.isfinite(q):
        return 'NA'
    return f'{q:.2g}' if q < .01 else f'{q:.3f}'


def style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(length=2.5, width=.6)


def save(fig, name, caption, manifest):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text().strip():
            continue
        box = text.get_window_extent(renderer)
        assert box.x0 >= -2 and box.y0 >= -2 and box.x1 <= canvas.width + 2 and box.y1 <= canvas.height + 2, (name, text.get_text(), box.bounds)
    entry = {'figure': name, 'width_mm': fig.get_figwidth()*25.4, 'height_mm':fig.get_figheight()*25.4}
    for ext in ['pdf', 'svg', 'png']:
        path = FIG / f'{name}.{ext}'
        fig.savefig(path, dpi=600, facecolor='white')
        entry[ext] = str(path.relative_to(ROOT)).replace('\\', '/')
        entry[f'{ext}_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    entry['caption'] = caption
    manifest.append(entry)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--primary-only',action='store_true')
    args=parser.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family':'Arial', 'font.size':7, 'axes.titlesize':8,
                         'axes.labelsize':7, 'xtick.labelsize':6.5, 'ytick.labelsize':6.5,
                         'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none',
                         'axes.linewidth':.6, 'lines.linewidth':1.05})
    cfg = json.loads((ROOT/'config/chromatin_rna.json').read_text())
    summary = pd.read_csv(OUT/'primary_results_with_effects.csv').set_index('pathway')
    fg = pd.read_csv(OUT/('fgsea_primary_complete.csv' if args.primary_only else 'fgsea_results.csv'))
    ca = pd.read_csv(OUT/('camera_primary_complete.csv' if args.primary_only else 'camera_results.csv'))
    catalog = pd.read_csv(ROOT/'outputs/tables/chromatin_atlas/gene_set_catalogue.csv')
    names = list(catalog.set_name)
    summary = summary.loc[names]
    selected = cfg['main_curve_panels']
    manifest = []
    scope = ('Reference marks were measured in external HT1080 cultures. RNA comprises 12 libraries, six matched passage/day '
             'comparisons from the same established CTRL and LOY populations sampled in August 2025 and May 2026. '
             'These are conditional associations, without independent population replication or a measurement of LOY-associated histone changes.')
    statistics = ('GSEA uses signed, equally period-weighted mean differences in TMM log2 CPM (prior count 2), all 25,461 fixed '
                  'representative genes, weighted exponent 1, minimum set size 15 and maximum universe size minus one. '
                  'CAMERA reuses the saved voom expression/precision weights and six passage-block intercepts plus separate period '
                  'effects (4 residual degrees of freedom); its residual-independence assumption is unverified for serial passages. '
                  'BH q values span the 20 custom sets separately for each method and contrast; original MSigDB global q values are unchanged.')
    # Complete 20-set overview, preserving frozen catalogue order rather than significance sorting.
    fig = plt.figure(figsize=(183/25.4, 190/25.4))
    ax = fig.add_axes([.425, .16, .235, .70])
    qa = fig.add_axes([.69, .16, .295, .70])
    fig.text(.025, .96, 'Measured HT1080 chromatin sets and LOY-associated RNA', fontsize=10, weight='bold')
    fig.text(.025, .925, 'Complete frozen collection • higher in LOY → positive NES', fontsize=7, color=INK)
    for i, name in enumerate(names):
        for axis in [ax, qa]:
            if i % 2 == 0:
                axis.axhspan(i-.5, i+.5, color='#f4f4f4', zorder=0)
        fig.text(.025, .16+.70*(1-(i+.5)/20), f'{LABELS[name]}  ({int(summary.loc[name,"size"]):,})',
                 va='center', fontsize=6.2)
    for contrast, marker, offset, color in [('combined','o',0,INK), ('AUG_2025','^',-.20,DOWN), ('MAY_2026','s',.20,UP)]:
        z = fg[fg.scenario.eq('window_2000') & fg.contrast.eq(contrast)].set_index('pathway').loc[names]
        ax.scatter(z.NES, np.arange(20)+offset, s=14 if contrast=='combined' else 10, marker=marker,
                   color=color, linewidths=0, label={'combined':'Combined','AUG_2025':'August','MAY_2026':'May'}[contrast], zorder=3)
    ax.set_ylim(19.5, -.5); ax.set_yticks([]); ax.axvline(0, color='#aaa', lw=.6)
    ax.set_xlabel('GSEA NES'); style(ax); ax.spines['left'].set_visible(False)
    ax.set_title('a  Enrichment direction', loc='left', pad=16, fontsize=8)
    ax.legend(frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(.5, -.10),
              fontsize=6, handletextpad=.25, columnspacing=.8)
    columns = ['BH_custom_20', 'camera_fixed_BH_custom_20', 'camera_estimated_BH_custom_20']
    for i, name in enumerate(names):
        row = summary.loc[name]
        for j, column in enumerate(columns):
            q = row[column]
            if np.isfinite(q) and q < .05:
                qa.add_patch(Rectangle((j-.46, i-.43), .92, .86, color='#dcebd8', zorder=1))
            qa.text(j, i, qtext(q), va='center', ha='center', fontsize=6.4,
                    weight='bold' if q < .05 else 'normal')
    qa.set_xlim(-.5, 2.5); qa.set_ylim(19.5, -.5); qa.set_yticks([])
    qa.set_xticks([0,1,2], ['GSEA','CAMERA\nρ = 0.01','CAMERA\nestimated ρ'])
    qa.tick_params(axis='x', length=0, pad=7)
    qa.spines[:].set_visible(False)
    qa.set_title('b  Combined custom-set BH q', loc='left', pad=16, fontsize=8)
    fig.text(.025, .052, 'Numbers in parentheses: genes. Shading: q < 0.05 within the 20-set custom family.', fontsize=6.8)
    fig.text(.025, .028, 'Serially sampled populations; reference chromatin association does not measure histone changes.', fontsize=6.8)
    caption = ('**All measured-reference gene sets.** (a) GSEA normalized enrichment scores for the combined and period-specific '
               'RNA rankings in frozen catalogue order. Parentheses give primary eligible gene counts. (b) Combined BH-adjusted '
               'P values from GSEA and CAMERA with fixed ρ=0.01 or estimated nonnegative within-set correlation. Green shading '
               'marks q<0.05, without implying equivalent null hypotheses across methods. '+statistics+' '+scope)
    save(fig, 'HT1080_chromatin_RNA_all_sets', caption, manifest)
    # Six question-led curves, selected in config before testing.
    curves = pd.read_csv(OUT/'selected_running_curves.csv.gz')
    fig = plt.figure(figsize=(183/25.4, 224/25.4))
    outer = fig.add_gridspec(3,2,left=.095,right=.98,bottom=.085,top=.87,hspace=.53,wspace=.28)
    fig.text(.025,.966,'RNA enrichment across measured chromatin contexts',fontsize=10,weight='bold')
    fig.text(.025,.941,'Six prespecified contexts • full primary RNA ranking',fontsize=7)
    for i,name in enumerate(selected):
        sub = outer[i//2,i%2].subgridspec(3,1,height_ratios=[3.4,.40,.8],hspace=.07)
        ax=fig.add_subplot(sub[0]); rug=fig.add_subplot(sub[1],sharex=ax); metric=fig.add_subplot(sub[2],sharex=ax)
        d=curves[curves.pathway.eq(name)]; row=summary.loc[name]
        color=UP if row.NES>0 else DOWN
        x=np.r_[0,d['rank'].to_numpy()]; y=np.r_[0,d.running_ES.to_numpy()]
        ax.plot(x,y,color=color,lw=1.2); ax.axhline(0,color='#aaa',lw=.55)
        ax.set_xlim(0,25461); ax.set_ylim(-.65,.65); ax.set_yticks([-.5,0,.5])
        ax.set_title(f'{string.ascii_lowercase[i]}  {LABELS[name]}',loc='left',fontsize=7.5,pad=24)
        ax.text(0,1.10,f'n = {int(row["size"]):,}; NES = {row.NES:.2f}; BH q = {qtext(row.BH_custom_20)}',
                transform=ax.transAxes,fontsize=6.4)
        ax.set_ylabel('Running ES'); style(ax); ax.tick_params(labelbottom=False,bottom=False)
        hit=d[d.is_hit]; rug.vlines(hit['rank'],0,1,color='#bdbdbd',lw=.22,rasterized=True)
        le=hit[hit.is_leading_edge]; rug.vlines(le['rank'],0,1,color=color,lw=.28,rasterized=True)
        rug.set_ylim(0,1); rug.axis('off')
        metric.fill_between(d['rank'],0,d.rank_metric,color='#bababa',linewidth=0,rasterized=True)
        metric.axhline(0,color='#888',lw=.5); metric.set_xticks([0,10000,20000,25461],['0','10k','20k','25.5k'])
        metric.set_yticks([]); metric.set_xlabel('Rank: higher → lower RNA in LOY',fontsize=6.3,labelpad=3); style(metric)
    fig.text(.025,.026,'Ticks: complete set (grey), leading edge (colour). Lower trace: signed RNA ranking metric.',fontsize=6.8)
    caption=('**Prespecified chromatin-context GSEA curves.** Six reference contexts were selected before testing. Upper traces '
             'show complete weighted running enrichment scores; rugs retain every member and highlight the leading edge. '
             'Lower traces show the signed RNA metric. Curves use identical ES and rank axes; leading edges and ES were independently '
             'reconstructed and checked against fgsea. '+statistics+' '+scope)
    save(fig,'HT1080_chromatin_RNA_GSEA_curves',caption,manifest)
    # ECDFs show entire distributions, including tails, without treating genes as sample replicates.
    effects=pd.read_csv(OUT/'gene_effect_evidence.csv')
    members=pd.read_csv(ROOT/'outputs/tables/chromatin_atlas/gene_set_membership_RNA_universe.csv')
    fig,axes=plt.subplots(3,2,figsize=(183/25.4,204/25.4))
    fig.subplots_adjust(left=.10,right=.98,bottom=.14,top=.89,wspace=.28,hspace=.62)
    fig.text(.025,.967,'RNA effect distributions within reference chromatin sets',fontsize=10,weight='bold')
    fig.text(.025,.937,'Empirical cumulative distributions • every gene retained',fontsize=7)
    source=[]
    for i,(ax,name) in enumerate(zip(axes.flat,selected)):
        genes=set(members.loc[members.set_name.eq(name),'gene_name']); hit=effects.gene_name.isin(genes)
        groups=[('Combined',effects.loc[hit,'mean_difference'],INK,'-'),('August',effects.loc[hit,'AUG_2025'],DOWN,'-'),
                ('May',effects.loc[hit,'MAY_2026'],UP,'-'),('Complement, combined',effects.loc[~hit,'mean_difference'],'#999','--')]
        for label,v,c,ls in groups:
            x=np.sort(v); y=np.arange(1,len(x)+1)/len(x)
            ax.step(x,y,where='post',color=c,lw=1.0,ls=ls,label=label)
            source.append(pd.DataFrame({'pathway':name,'group':label,'effect':x,'cumulative_fraction':y}))
        ax.axvline(0,color='#bbb',lw=.6); ax.set_ylim(0,1)
        ax.set_xlabel('LOY − CTRL (log$_2$ CPM)'); ax.set_ylabel('Fraction of genes')
        ax.set_title(f'{string.ascii_lowercase[i]}  {LABELS[name]}\nSet n = {len(genes):,}; complement n = {len(effects)-len(genes):,}',loc='left',fontsize=7.2,pad=7)
        style(ax)
    handles,labels=axes.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.045),ncol=4,frameon=False,fontsize=6.5)
    fig.text(.025,.019,'Genes describe distributions; they are not independent biological replicates. No gene-level P values are assigned.',fontsize=6.5)
    pd.concat(source).to_csv(OUT/'figure_ECDF_source.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    caption=('**Complete RNA effect distributions.** Empirical cumulative distribution functions show all genes in each of the '
             'six prespecified reference sets for combined, August and May differences. Grey dashed curves show the combined '
             'distribution of the full set complement, including unevaluable promoters in the primary background. The separate '
             'evaluable-promoter sensitivity excludes those genes. Plot limits retain the full observed range; set sizes and '
             'complement sizes are indicated. No confidence intervals or gene-as-replicate hypothesis tests are used. '+scope)
    save(fig,'HT1080_chromatin_RNA_effect_distributions',caption,manifest)
    if args.primary_only:
        pd.DataFrame(manifest).to_csv(FIG/'primary_manifest.csv',index=False)
        (ROOT/'outputs/captions/chromatin_rna_primary.md').write_text(
            '\n\n'.join(m['caption'] for m in manifest),encoding='utf-8')
        print('Three complete primary figures exported; sensitivity figure and final report remain pending.')
        return
    # Sensitivity figure: predetermined scenarios and all sets, no selection on q.
    scenarios=['window_2000','window_1000','window_2500','relaxed_q_primary_window','any_annotated_TSS',
               'evaluable_promoters_only','prior025_mean_difference','RLE_mean_difference',
               'residualized_abundance_length','residualized_abundance_length_arm']
    labels=['Primary','±1 kb','±2.5 kb','Relaxed peak q','Any TSS','Evaluable only','Prior 0.25','RLE',
            'Abundance + length','+ chromosome arm']
    matrix=fg[fg.contrast.eq('combined') & fg.scenario.isin(scenarios)].pivot(index='pathway',columns='scenario',values='NES').reindex(index=names,columns=scenarios)
    qs=fg[fg.contrast.eq('combined') & fg.scenario.isin(scenarios)].pivot(index='pathway',columns='scenario',values='BH_custom_20').reindex(index=names,columns=scenarios)
    fig=plt.figure(figsize=(183/25.4,190/25.4)); ax=fig.add_axes([.43,.25,.54,.60])
    fig.text(.025,.963,'Sensitivity to reference definitions and RNA covariates',fontsize=10,weight='bold')
    fig.text(.025,.927,'All custom sets • GSEA NES • fixed scenario order',fontsize=7)
    cmap=plt.get_cmap('RdBu_r').copy(); cmap.set_bad('#bfbfbf')
    im=ax.imshow(matrix,aspect='auto',vmin=-2.5,vmax=2.5,cmap=cmap,interpolation='nearest')
    ax.set_yticks(np.arange(20),[LABELS[n] for n in names],fontsize=6.2)
    ax.set_xticks(np.arange(len(scenarios)),labels,rotation=55,ha='right',fontsize=6.4)
    ax.tick_params(length=0,pad=5); ax.spines[:].set_visible(False)
    for i in range(20):
        for j in range(len(scenarios)):
            if qs.iloc[i,j]<.05:
                ax.text(j,i,'•',ha='center',va='center',fontsize=8,color='black' if abs(matrix.iloc[i,j])<1.7 else 'white')
    cax=fig.add_axes([.69,.09,.22,.018]); fig.colorbar(im,cax=cax,orientation='horizontal',label='GSEA NES')
    fig.text(.025,.068,'Dot: within-scenario custom BH q < 0.05. Grey: fewer than 15 eligible genes.',fontsize=6.6)
    fig.text(.025,.035,'Sensitivity scenarios do not replace the primary test. Covariate residuals are descriptive; arm adjustment is not CN correction.',fontsize=6.1)
    matrix.to_csv(OUT/'figure_sensitivity_NES.csv')
    caption=('**Prespecified sensitivity analyses.** All 20 sets are shown in frozen order. Colour represents GSEA NES, clipped '
             'at ±2.5 for display; complete numerical results are exported. Dots denote BH q<0.05 across the 20 sets within '
             'each scenario and do not imply a single family across alternative analyses. Grey means the set has fewer than '
             '15 genes; primary gene-body definitions remain unchanged under promoter alternatives. Additional pooled-sensitivity '
             'BH and per-scenario BY values are exported. Residualized rankings remove fitted trends in mean expression and '
             'log10 feature length (four-degree natural splines), then optionally chromosome arm. These descriptive regressions '
             'can remove real biology as well as technical structure and are not matched copy-number correction. '+scope)
    save(fig,'HT1080_chromatin_RNA_sensitivity',caption,manifest)
    pd.DataFrame(manifest).to_csv(FIG/'manifest.csv',index=False)
    captions='\n\n'.join(f'## {m["figure"]}\n\n{m["caption"]}' for m in manifest)+'\n'
    (ROOT/'outputs/captions/chromatin_rna.md').write_text(captions,encoding='utf-8')
    # Read reviewed interpretation written alongside the analysis, followed by complete tables and figures.
    narrative=(ROOT/'CHROMATIN_RNA_ENRICHMENT.md').read_text(encoding='utf-8') if (ROOT/'CHROMATIN_RNA_ENRICHMENT.md').exists() else 'Interpretation pending.'
    import mistune
    from bs4 import BeautifulSoup
    body=BeautifulSoup(mistune.create_markdown(plugins=['table'])(narrative),'html.parser')
    # The markdown lives at the analysis root; the HTML lives two directories below it.
    for a in body.find_all('a',href=True):
        if not a['href'].startswith(('http://','https://','#','mailto:')):
            a['href']='../../'+a['href']
    body=str(body)
    pages=[]
    for m in manifest:
        data=base64.b64encode((ROOT/m['png']).read_bytes()).decode()
        links=' · '.join(f'<a href="../../{m[ext]}">{ext.upper()}</a>' for ext in ['pdf','svg','png'])
        pages.append(f'<section><h2>{html.escape(m["figure"])}</h2><p>{links}</p><img alt="{html.escape(m["figure"])}" src="data:image/png;base64,{data}"><p>{html.escape(m["caption"])}</p></section>')
    table=summary[['size','NES','BH_custom_20','camera_fixed_BH_custom_20','camera_estimated_BH_custom_20',
                   'combined_mean','AUG_2025_mean','MAY_2026_mean','combined_fraction_positive']].rename(index=LABELS)
    table=table.rename(columns={'size':'Genes','NES':'GSEA NES','BH_custom_20':'GSEA BH q',
        'camera_fixed_BH_custom_20':'CAMERA fixed BH q','camera_estimated_BH_custom_20':'CAMERA estimated BH q',
        'combined_mean':'Mean combined','AUG_2025_mean':'Mean August','MAY_2026_mean':'Mean May',
        'combined_fraction_positive':'Fraction RNA positive'})
    doc='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>HT1080 measured chromatin and RNA</title><style>body{font:16px Arial,sans-serif;max-width:1250px;margin:40px auto;padding:0 22px;color:#252525;line-height:1.5}img{width:100%;max-width:950px}table{border-collapse:collapse;font-size:12px}th,td{padding:7px;border-bottom:1px solid #ddd}section{margin-top:50px}a{color:#0072b2}.table-scroll{overflow-x:auto}</style><body>'+body+'<h2>Complete primary results</h2><div class="table-scroll">'+table.to_html(float_format=lambda x:f'{x:.4g}')+'</div>'+''.join(pages)+'</body></html>'
    (ROOT/'outputs/report/chromatin_rna_enrichment.html').write_text(doc,encoding='utf-8')
    print('Exported four figures, source data, captions and standalone HTML.')


if __name__=='__main__':
    main()
