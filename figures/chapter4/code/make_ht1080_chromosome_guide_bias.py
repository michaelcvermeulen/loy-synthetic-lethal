"""Exploratory chromosome-level guide trajectories and baseline sensitivity."""
from pathlib import Path
import sys, os, json, hashlib

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_chromosome_guide_bias'
DATA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(DATA / 'mpl_cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import pymupdf
sys.path.insert(0, str(THESIS / 'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(9)

SCREEN = THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
FILES = {
    'annotation': SCREEN/'brunello/broadgpp-brunello-library-contents.txt',
    'library': SCREEN/'brunello/broadgpp-brunello-library-corrected.txt',
    'counts': SCREEN/'count/ht1080_screen2_R1.count.txt',
    'mageck': SCREEN/'test/POST_LOY_vs_CTRL.sgrna_summary.txt',
}
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']
GROUPS = CHROMS + ['NTC']
STEM = 'ch04_supp_ht1080_chromosome_guide_bias'


def prepare():
    ann = pd.read_csv(FILES['annotation'], sep='\t')
    lib = pd.read_csv(FILES['library'], sep='\t').rename(columns={'sgRNAID':'guide_id', 'Seq':'sequence'})
    counts = pd.read_csv(FILES['counts'], sep='\t')
    samples = counts.columns[2:].tolist()
    totals = counts[samples].sum()
    ann = ann.rename(columns={'sgRNA Target Sequence':'sequence', 'Target Gene Symbol':'original_gene',
        'Genomic Sequence':'reference_accession', 'Position of Base After Cut (1-based)':'cut_position_1based'})
    assert ann.sequence.is_unique and lib.sequence.is_unique and counts.sgRNA.is_unique
    x = lib.merge(ann[['sequence','original_gene','reference_accession','cut_position_1based']],
        on='sequence', how='left', validate='one_to_one', indicator=True)
    assert x['_merge'].eq('both').all()
    x = x.drop(columns='_merge').merge(counts.rename(columns={'sgRNA':'guide_id','Gene':'count_gene'}),
        on='guide_id', validate='one_to_one')
    assert x.gene.eq(x.count_gene).all() and len(x) == 77441
    accession_number = x.reference_accession.str.extract(r'^NC_0000(\d{2})\.')[0]
    mapping = {f'{i:02d}': str(i) for i in range(1, 23)} | {'23':'X','24':'Y'}
    x['chromosome'] = accession_number.map(mapping)
    ntc = x.gene.eq('Non-Targeting-Control')
    assert ntc.sum() == 1000 and x.loc[~ntc,'chromosome'].notna().all()
    assert x.loc[ntc,'reference_accession'].isna().all()
    x.loc[ntc,'chromosome'] = 'NTC'
    saved = pd.read_csv(FILES['mageck'], sep='\t')[['sgrna','LFC']].rename(columns={'sgrna':'guide_id','LFC':'saved_mageck_endpoint_LFC'})
    x = x.merge(saved, on='guide_id', validate='one_to_one')
    for sample in samples:
        x[sample+'_CPM'] = x[sample]/totals[sample]*1e6
    x['endpoint_CTRL_mean_CPM'] = x[[s+'_CPM' for s in samples[2:4]]].mean(axis=1)
    x['endpoint_LOY_mean_CPM'] = x[[s+'_CPM' for s in samples[4:6]]].mean(axis=1)
    x['D0_absent_both'] = x[samples[:2]].sum(axis=1).eq(0)
    x['endpoint_absent_all4'] = x[samples[2:]].sum(axis=1).eq(0)
    x['paired_observed'] = ~(x.D0_absent_both | x.endpoint_absent_all4)
    x['D0_ge50_both'] = x.D0_CTRL.ge(50) & x.D0_LOY.ge(50)
    x['D0_LFC_CPM_plus_1'] = (np.log2(x.D0_LOY_CPM+1)-np.log2(x.D0_CTRL_CPM+1)).mask(x.D0_absent_both)
    x['endpoint_LFC_CPM_plus_1'] = (np.log2(x.endpoint_LOY_mean_CPM+1)-np.log2(x.endpoint_CTRL_mean_CPM+1)).mask(x.endpoint_absent_all4)
    x['endpoint_mageck_display_LFC'] = x.saved_mageck_endpoint_LFC.mask(x.endpoint_absent_all4)
    x['delta_LFC_matched'] = x.endpoint_LFC_CPM_plus_1-x.D0_LFC_CPM_plus_1
    valid = x.loc[x.paired_observed].copy()
    ntc_delta = float(valid.loc[valid.chromosome.eq('NTC'),'delta_LFC_matched'].median())
    ntc_high_delta = float(valid.loc[valid.chromosome.eq('NTC') & valid.D0_ge50_both,'delta_LFC_matched'].median())
    x['delta_relative_NTC_median'] = x.delta_LFC_matched-ntc_delta
    gene = valid.loc[~valid.chromosome.eq('NTC')].groupby(['chromosome','gene'], sort=False).agg(
        n_observed_guides=('guide_id','size'), median_guide_delta=('delta_LFC_matched','median')).reset_index()
    gene['median_guide_delta_relative_NTC'] = gene.median_guide_delta-ntc_delta
    records = []
    for chrom in GROUPS:
        all_rows = x.loc[x.chromosome.eq(chrom)]
        sub = all_rows.loc[all_rows.paired_observed]
        high = sub.loc[sub.D0_ge50_both]
        gd = gene.loc[gene.chromosome.eq(chrom)]
        records.append({
            'chromosome':chrom, 'n_library_guides':len(all_rows), 'n_paired_guides':len(sub),
            'n_missing_either_stage':len(all_rows)-len(sub), 'n_genes':sub.gene.nunique() if chrom!='NTC' else 0,
            'median_D0_LFC':sub.D0_LFC_CPM_plus_1.median(),
            'median_endpoint_matched_LFC':sub.endpoint_LFC_CPM_plus_1.median(),
            'median_endpoint_saved_mageck_LFC':sub.endpoint_mageck_display_LFC.median(),
            'median_delta':sub.delta_LFC_matched.median(),
            'q25_delta':sub.delta_LFC_matched.quantile(.25), 'q75_delta':sub.delta_LFC_matched.quantile(.75),
            'median_delta_relative_NTC':sub.delta_LFC_matched.median()-ntc_delta,
            'q25_delta_relative_NTC':sub.delta_LFC_matched.quantile(.25)-ntc_delta,
            'q75_delta_relative_NTC':sub.delta_LFC_matched.quantile(.75)-ntc_delta,
            'fraction_positive_delta_relative_NTC':sub.delta_LFC_matched.gt(ntc_delta).mean(),
            'n_D0_ge50_both':len(high), 'median_delta_D0_ge50':high.delta_LFC_matched.median(),
            'median_delta_D0_ge50_relative_NTC':high.delta_LFC_matched.median()-ntc_high_delta,
            'median_gene_delta_relative_NTC':gd.median_guide_delta_relative_NTC.median() if len(gd) else np.nan,
        })
    summary = pd.DataFrame(records)
    # This is annotation-based chromosome assignment; record genes with guides on
    # multiple chromosomes instead of silently relocating them by gene symbol.
    cross = x.loc[~ntc].groupby('gene').chromosome.nunique()
    cross_names = cross[cross.gt(1)].index.tolist()
    metadata = {
        'library_guides':len(x), 'chromosome_annotated_guides':int((x.chromosome!='NTC').sum()),
        'non_targeting_guides':int((x.chromosome=='NTC').sum()),
        'paired_guides_all':int(x.paired_observed.sum()),
        'paired_chromosome_guides':int((x.paired_observed & x.chromosome.ne('NTC')).sum()),
        'paired_NTC_guides':int((x.paired_observed & x.chromosome.eq('NTC')).sum()),
        'guides_missing_either_stage':int((~x.paired_observed).sum()),
        'NTC_median_delta':ntc_delta, 'NTC_D0_ge50_median_delta':ntc_high_delta,
        'sample_library_totals':{k:int(v) for k,v in totals.items()},
        'gene_symbols_with_multiple_chromosome_annotations':cross_names,
        'chromosome_order':GROUPS,
        'paired_matched_range':[float(valid[['D0_LFC_CPM_plus_1','endpoint_LFC_CPM_plus_1']].min().min()),
            float(valid[['D0_LFC_CPM_plus_1','endpoint_LFC_CPM_plus_1']].max().max())],
        'paired_original_range':[float(valid[['D0_LFC_CPM_plus_1','endpoint_mageck_display_LFC']].min().min()),
            float(valid[['D0_LFC_CPM_plus_1','endpoint_mageck_display_LFC']].max().max())],
    }
    assert len(summary)==25 and summary.n_paired_guides.sum()==77086
    assert metadata['paired_NTC_guides']==998
    x.to_csv(DATA/'all_guide_chromosome_values.csv',index=False)
    summary.to_csv(DATA/'chromosome_summary.csv',index=False)
    gene.to_csv(DATA/'gene_median_deltas.csv',index=False)
    x.loc[x.gene.isin(cross_names),['guide_id','gene','original_gene','chromosome','reference_accession']].to_csv(DATA/'genes_with_multiple_chromosome_annotations.csv',index=False)
    (DATA/'summary.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    return x, summary, gene, metadata


def trajectories(x, mode):
    fig, axes = plt.subplots(6,4,figsize=(11.6,12.2))
    fig.subplots_adjust(left=.065,right=.985,top=.93,bottom=.10,wspace=.20,hspace=.67)
    endpoint = 'endpoint_LFC_CPM_plus_1' if mode=='matched' else 'endpoint_mageck_display_LFC'
    limits = (-4.5,4.5) if mode=='matched' else (-10,10)
    totals = 0
    for i,(ax,chrom) in enumerate(zip(axes.flat,CHROMS)):
        sub=x.loc[x.chromosome.eq(chrom)&x.paired_observed]
        values=sub[['D0_LFC_CPM_plus_1',endpoint]].to_numpy()
        segments=np.stack([np.column_stack([np.zeros(len(sub)),values[:,0]]),
            np.column_stack([np.ones(len(sub)),values[:,1]])],axis=1)
        ax.add_collection(LineCollection(segments,colors='#555555',linewidths=.35,alpha=.025,rasterized=True,zorder=1))
        for xx in [0,1]:
            ax.scatter(np.full(len(sub),xx),values[:,xx],s=2,color='#555555',alpha=.025,lw=0,rasterized=True,zorder=2)
        ax.plot([0,1],np.median(values,axis=0),color='#222222',lw=1.5,marker='o',ms=3,zorder=4)
        ax.axhline(0,color='#AAAAAA',lw=.6,ls='--',zorder=0)
        ax.set_xlim(-.12,1.12);ax.set_ylim(*limits)
        ax.set_xticks([0,1],['D0','Endpoint'])
        ax.set_yticks([-4,0,4] if mode=='matched' else [-10,0,10])
        ax.tick_params(labelsize=8,pad=2,length=2.5)
        if i%4:ax.tick_params(labelleft=False)
        ax.set_title(f'Chr {chrom}  (n = {len(sub):,})',fontsize=9,pad=4)
        totals+=len(sub)
    fig.supylabel('Guide log₂ fold change (LOY / sgCTRL)',x=.012,fontsize=10)
    heading='Matched CPM + 1 contrasts at D0 and endpoint' if mode=='matched' else 'Original display: D0 CPM + 1 contrast to saved MAGeCK endpoint LFC'
    fig.text(.065,.975,heading,fontsize=11,weight='semibold')
    fig.text(.985,.975,'Exploratory QC',fontsize=9,ha='right',color='#666666')
    fig.text(.065,.055,'Grey: all observed guide pairs (2.5% opacity). Black: stage-wise medians. Density reflects guide count; no downsampling.',fontsize=8)
    if mode=='matched':
        fig.text(.065,.034,'Both stages: log₂[(LOY CPM + 1)/(sgCTRL CPM + 1)]. Endpoint CPM is the within-genotype mean of R1/R2.',fontsize=8)
    else:
        fig.text(.065,.034,'The two estimators use different normalisation. This companion view is not used to quantify chromosome shifts.',fontsize=8)
    fig.text(.065,.013,'All panels share limits and retain all finite values. Guides absent from both genotypes at either stage are omitted.',fontsize=8)
    assert totals==76088
    return fig


def summary_plot(summary, meta):
    fig=plt.figure(figsize=(11.6,9.0))
    a=fig.add_axes([.085,.19,.27,.68])
    b=fig.add_axes([.45,.19,.34,.68],sharey=a)
    yy=np.arange(len(summary))
    labels=[f'Chr {c}' if c!='NTC' else 'NTC' for c in summary.chromosome]
    for ax in [a,b]:
        ax.axvline(0,color='#AAAAAA',lw=.65,ls='--',zorder=0)
        ax.axhline(23.5,color='#CCCCCC',lw=.7)
        ax.set_ylim(24.7,-.9);ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y',length=0)
    a.hlines(yy,summary.median_D0_LFC,summary.median_endpoint_matched_LFC,color='#AAAAAA',lw=1)
    a.scatter(summary.median_D0_LFC,yy,s=26,facecolors='white',edgecolors='#333333',lw=.8,label='D0',zorder=3)
    a.scatter(summary.median_endpoint_matched_LFC,yy,s=20,color='#333333',label='Endpoint',zorder=3)
    a.set_yticks(yy,labels);a.set_xlim(-.20,.42)
    a.set_xticks([-.2,0,.2,.4]);a.set_xlabel('Median guide LFC')
    a.legend(loc='lower left',bbox_to_anchor=(-.02,1.025),ncol=2,fontsize=9,handletextpad=.3,columnspacing=1)
    b.hlines(yy,summary.q25_delta_relative_NTC,summary.q75_delta_relative_NTC,color='#C5C5C5',lw=2,zorder=1)
    b.scatter(summary.median_delta_relative_NTC,yy,s=25,color='#222222',label='All observed guides',zorder=3)
    b.scatter(summary.median_delta_D0_ge50_relative_NTC,yy+.14,s=24,marker='D',facecolor='white',edgecolor='#222222',lw=.8,label='D0 ≥50 reads in both',zorder=4)
    b.set_xlim(-.58,.60);b.set_xticks([-.5,-.25,0,.25,.5]);b.tick_params(labelleft=False)
    b.set_xlabel('Median ΔLFC relative to NTC')
    b.legend(loc='lower left',bbox_to_anchor=(-.02,1.025),fontsize=8.5,handletextpad=.3)
    for i,r in summary.iterrows():
        b.text(1.04,i,f'{r.median_delta_relative_NTC:+.3f}',transform=b.get_yaxis_transform(),va='center',fontsize=8)
        b.text(1.28,i,f'{r.n_paired_guides:,} / {r.n_D0_ge50_both:,}',transform=b.get_yaxis_transform(),va='center',fontsize=8)
    b.text(1.04,1.02,'Shift',transform=b.transAxes,fontsize=8)
    b.text(1.28,1.02,'n all / ≥50',transform=b.transAxes,fontsize=8)
    fig.text(.025,.925,'a',fontsize=14,weight='semibold');fig.text(.40,.925,'b',fontsize=14,weight='semibold')
    fig.text(.085,.965,'Chromosome-level summaries using matched CPM + 1 contrasts',fontsize=11,weight='semibold')
    fig.text(.085,.105,'ΔLFC = endpoint − D0 for each guide. Panel b subtracts the median NTC ΔLFC within each guide set.',fontsize=8.5)
    fig.text(.085,.079,'Grey intervals: interquartile ranges of guide shifts, not confidence intervals. Diamonds: D0 read-threshold sensitivity.',fontsize=8.5)
    fig.text(.085,.053,f'NTC median ΔLFC: {meta["NTC_median_delta"]:+.3f} (n = {meta["paired_NTC_guides"]}); ≥50 subset: {meta["NTC_D0_ge50_median_delta"]:+.3f}.',fontsize=8.5)
    fig.text(.085,.027,'Descriptive chromosome associations; no significance tests. Guide counts are not independent biological replicates.',fontsize=8.5)
    return fig


def export(fig,name,report):
    stem=STEM+'_'+name
    report.savefig(fig,dpi=300)
    fig.savefig(OUT/f'{stem}.pdf',dpi=300)
    fig.savefig(OUT/f'{stem}.svg',dpi=300)
    fig.savefig(OUT/f'{stem}.png',dpi=600)
    fig.savefig(OUT/f'review/{stem}_preview.png',dpi=160)
    plt.close(fig)


def workbook(summary,gene,meta):
    definitions=[
        ('Purpose','Exploratory chromosome-associated guide shifts, not proof of technical bias or new gene hit calling.'),
        ('Chromosomes','Original Brunello NC_000001-24 reference accession assigned per unique guide sequence; 1-22, X and Y. NTC guides shown separately.'),
        ('CPM','Sample raw counts divided by the sum across all 77,441 guide rows, times 1,000,000.'),
        ('Matched LFC','log2[(LOY CPM + 1)/(CTRL CPM + 1)] at both D0 and endpoint; endpoint uses mean CPM across two outgrowth cultures per genotype.'),
        ('Delta','Paired guide endpoint matched LFC minus D0 LFC. This descriptive CPM calculation is separate from saved MAGeCK/drugZ results.'),
        ('NTC reference','Median paired-guide delta among observed NTC guides. Subtracted from each chromosome median; the >=50 sensitivity uses its own eligible NTC median.'),
        ('Eligibility','At least one count across the two D0 samples and at least one count across the four endpoint samples. One-genotype absence is retained.'),
        ('Sensitivity','At least 50 raw reads in EACH D0 sample; descriptive sensitivity using the existing QC read-count threshold.'),
        ('Gene weighted','Median of within-gene guide-median deltas within each chromosome, minus the all-observed NTC guide-median delta.'),
        ('Uncertainty','No P values or confidence intervals. Guide IQRs describe variation among guides; guides and genes are not independent biological replicates.'),
        ('Limitations','Differences can reflect target composition, biology, chromosome copy number, guide specificity, starting representation or technical effects. No cause is assigned.'),
        ('Saved results','Original gene-level drugZ scores and MAGeCK estimates are not rerun or replaced. The original mixed-estimator overlay is a companion view only.'),
        ('Full guide data','all_guide_chromosome_values.csv contains all 77,441 guide rows, counts, CPM, membership, eligibility and both LFC estimators.'),
    ]
    with pd.ExcelWriter(DATA/'ht1080_chromosome_bias_summary.xlsx',engine='openpyxl') as writer:
        for name,frame in [('Chromosome summary',summary),('Gene medians',gene),('Definitions',pd.DataFrame(definitions,columns=['Item','Definition']))]:
            frame.to_excel(writer,sheet_name=name,index=False)
            ws=writer.sheets[name];ws.freeze_panes='B2';ws.auto_filter.ref=ws.dimensions
            from openpyxl.styles import Font,PatternFill
            from openpyxl.utils import get_column_letter
            for cell in ws[1]:cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='555555')
            for i,col in enumerate(frame,1):ws.column_dimensions[get_column_letter(i)].width=min(42,max(15,len(col)+2))
            if name=='Definitions':ws.column_dimensions['B'].width=140


def main():
    x,summary,gene,meta=prepare()
    workbook(summary,gene,meta)
    with PdfPages(OUT/f'{STEM}_report.pdf') as report:
        export(trajectories(x,'matched'),'matched_trajectories',report)
        export(summary_plot(summary,meta),'summary',report)
        export(trajectories(x,'original'),'original_estimators',report)
    doc=pymupdf.open(OUT/f'{STEM}_report.pdf')
    assert len(doc)==3
    clipped=[]
    for p,page in enumerate(doc):
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines',[]):
                for span in line['spans']:
                    x0,y0,x1,y1=span['bbox']
                    if x0<0 or y0<0 or x1>page.rect.width or y1>page.rect.height:clipped.append([p+1,span['text']])
    assert not clipped,clipped
    fonts=sorted({f[3] for page in doc for f in page.get_fonts()})
    assert all('OpenSans' in f for f in fonts)
    # Preserve the original Y-guide values exactly, including its omitted guides.
    y=pd.read_csv(OUT/'source_data/ht1080_y_guides/all_y_guides.csv')
    yy=x.loc[x.chromosome.eq('Y')].merge(y,on='guide_id',validate='one_to_one',suffixes=('','_prior'))
    assert len(yy)==188
    assert np.allclose(yy.D0_LFC_CPM_plus_1,yy.D0_log2FC_CPM_plus_1,equal_nan=True)
    assert np.allclose(yy.saved_mageck_endpoint_LFC,yy.mageck_LFC)
    assert meta['paired_matched_range'][0]>-4.5 and meta['paired_matched_range'][1]<4.5
    assert meta['paired_original_range'][0]>-10 and meta['paired_original_range'][1]<10
    manifest={
        'status':'Exploratory chromosome bias review',
        'inputs':{k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in FILES.items()},
        'matched_metric':'log2[(LOY CPM + 1)/(CTRL CPM + 1)] at both stages; endpoint uses within-genotype mean CPM',
        'NTC_reference':'Median guide-level delta, separately for all-observed and D0>=50-in-both subsets',
        'no_new_significance_tests':True,'no_primary_analysis_changes':True,
        'line_alpha':.025,'same_opacity_all_chromosomes':True,'no_downsampling':True,
        'dense_guide_marks':'Rasterised at 300 dpi in PDF/SVG; labels and median traces remain vector',
        'font':fonts,'metadata':meta,
        'validation':{'all_77441_guides_joined_uniquely':True,'24_chromosomes_plus_NTC':True,
            'prior_Y_values_agree':True,'all_finite_guide_values_within_axes':True,'no_clipped_text':True},
    }
    (DATA/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({'report':str(OUT/f'{STEM}_report.pdf'),'metadata':meta,
        'largest_shifts':summary.loc[summary.chromosome.isin(['Y','18','13','22','9']),
            ['chromosome','n_paired_guides','median_delta_relative_NTC','median_delta_D0_ge50_relative_NTC','median_gene_delta_relative_NTC']].to_dict('records')},indent=2))


if __name__=='__main__':
    main()
