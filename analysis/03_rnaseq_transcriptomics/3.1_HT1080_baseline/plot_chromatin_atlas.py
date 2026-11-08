"""Publication figure and standalone methods report for the measured reference atlas."""
from pathlib import Path
import base64
import html
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs/tables/chromatin_atlas'
FIG=ROOT/'outputs/figures/chromatin_atlas'
PROV=ROOT/'provenance/chromatin_atlas'
ORDER=['H3K4me3','H3K27ac','H3K4me1','H3K27me3','H2AK119ub1','H3K9me3','H4K20me3','H4K20me1','LaminB1','LaminB2']
COMBINATIONS={
 'H3K4me3_H3K27me3_COMARKED_PROMOTER':'H3K4me3 + H3K27me3',
 'H3K27me3_H2AK119ub1_JOINT_PROMOTER':'H3K27me3 + H2AK119ub1',
 'H3K4me3_H3K27ac_JOINT_PROMOTER':'H3K4me3 + H3K27ac',
 'H3K9me3_H4K20me3_JOINT_PROMOTER':'H3K9me3 + H4K20me3'}


def main():
    FIG.mkdir(parents=True,exist_ok=True)
    genes=pd.read_csv(OUT/'HT1080_promoter_marks_RNA_universe.csv')
    marks=pd.read_csv(OUT/'promoter_mark_summary.csv').set_index('mark').loc[ORDER]
    sets=pd.read_csv(OUT/'gene_set_catalogue.csv').set_index('set_name')
    sensitivity=pd.read_csv(OUT/'gene_set_size_sensitivity.csv').set_index('set_name')
    summary=json.loads((PROV/'atlas_summary.json').read_text())
    pqc=pd.read_csv(OUT/'peak_input_qc.csv')
    signal_qc=pd.read_csv(OUT/'bigwig_interval_qc.csv')
    scales=[]
    for mark,sub in signal_qc.groupby('mark'):
        standard=sub[sub.chrom.isin(['chr'+str(i) for i in range(1,23)]+['chrX','chrY'])]
        scales.append({'mark':mark,'signal_units': 'Signed deposited signal; normalization unresolved' if mark.startswith('Lamin') else 'RPKM according to GEO processing metadata',
                       'median_standard_chromosome_interval_bp':float(standard.median_interval_width.median()),
                       'has_negative_signal':bool(sub.minimum_signal.lt(0).any()),
                       'use':'Coarse domain context only; do not interpret as promoter-resolution RPKM' if mark.startswith('Lamin') else 'Reference signal; do not compare absolute levels between antibodies'})
    pd.DataFrame(scales).to_csv(OUT/'signal_scale_and_resolution.csv',index=False)
    font=Path('C:/Windows/Fonts/arial.ttf')
    if font.exists():font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,
                         'xtick.labelsize':6.5,'ytick.labelsize':6.5,'pdf.fonttype':42,
                         'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6,
                         'savefig.facecolor':'white','figure.facecolor':'white'})
    fig=plt.figure(figsize=(183/25.4,215/25.4))
    grid=fig.add_gridspec(2,2,left=.14,right=.98,bottom=.16,top=.88,
                         width_ratios=[1.05,1],height_ratios=[1,1.05],wspace=.72,hspace=.60)
    axes=[fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,0]),fig.add_subplot(grid[1,1])]
    for label,ax in zip('abcd',axes):
        ax.text(-.24,1.08,label,transform=ax.transAxes,fontweight='bold',fontsize=10)
        ax.spines[['top','right']].set_visible(False)
    ax=axes[0]
    colors=['#3475A8']*8+['#858585']*2
    ax.barh(np.arange(10),marks.detected_any_bp,color=colors,height=.63)
    ax.set_yticks(np.arange(10),ORDER);ax.invert_yaxis()
    ax.set_xlim(0,12000);ax.set_xticks([0,5000,10000],['0','5,000','10,000'])
    for y,n in enumerate(marks.detected_any_bp):ax.text(n+130,y,f'{n:,}',va='center',fontsize=6.5)
    ax.set_xlabel('Genes with a promoter peak')
    ax.set_title('Measured reference marks',loc='left',pad=9)

    ax=axes[1]
    ns=[int(sets.loc[x,'n_RNA_universe']) for x in COMBINATIONS]
    ax.barh(range(4),ns,color=['#76559C','#AB7A3A','#3475A8','#397F6D'],height=.56)
    ax.set_yticks(range(4),[x.replace(' + ','\n+ ') for x in COMBINATIONS.values()])
    ax.invert_yaxis();ax.set_xlim(0,1900);ax.set_xticks([0,500,1000,1500])
    for y,n in enumerate(ns):ax.text(n+25,y,f'{n:,}',va='center',fontsize=7)
    ax.set_title('Co-marked promoters',loc='left',pad=9)
    ax.set_xlabel('Genes; same promoter window')

    ax=axes[2]
    histones=ORDER[:8]
    binary=genes[[m+'__detected' for m in histones]].astype('boolean')
    values=binary.astype(float).to_numpy()
    known=np.isfinite(values).all(axis=1)
    pattern=np.nan_to_num(values,nan=0).dot(2**np.arange(7,-1,-1))
    # Reference-mark ordering only; RNA effect values are never loaded here.
    order=np.lexsort((genes.gene_id.to_numpy(),-pattern,~known))
    rows=[];binned=[]
    for i,start in enumerate(range(0,len(order),50)):
        ids=order[start:start+50];v=values[ids]
        denominator=np.isfinite(v).sum(axis=0)
        means=np.divide(np.nansum(v,axis=0),denominator,out=np.full(8,np.nan),where=denominator>0)
        binned.append(means)
        for j,mark in enumerate(histones):
            rows.append({'bin':i,'first_gene_order':start+1,'last_gene_order':start+len(ids),
                         'mark':mark,'fraction_detected':means[j],'evaluable_genes':int(denominator[j])})
    pd.DataFrame(rows).to_csv(OUT/'figure_heatmap_bins.csv',index=False)
    pd.DataFrame({'order':np.arange(1,len(order)+1),'gene_id':genes.iloc[order].gene_id.to_numpy()}).to_csv(OUT/'figure_gene_order.csv',index=False)
    cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#D9D9D9')
    im=ax.imshow(np.array(binned),aspect='auto',interpolation='nearest',vmin=0,vmax=1,cmap=cmap,rasterized=True)
    ax.set_xticks(range(8),histones,rotation=60,ha='right')
    ax.set_yticks([0,199,399,len(binned)-1],['1','10,000','20,000','25,461'])
    ax.set_ylabel('Genes ordered by reference marks')
    ax.set_title('All ranked genes',loc='left',pad=9)
    cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.048,pad=.31)
    cb.set_ticks([0,.5,1]);cb.set_label('Fraction with a peak per 50-gene bin',fontsize=6.5)

    ax=axes[3]
    chosen=[list(COMBINATIONS)[i] for i in [0,1,3]]
    palette=['#76559C','#AB7A3A','#397F6D']
    for name,color,marker in zip(chosen,palette,['o','s','^']):
        y=[sensitivity.loc[name,f'window_{w}'] for w in [1000,2000,2500]]
        ax.plot([1,2,2.5],y,color=color,marker=marker,markersize=4,linewidth=1,
                label=COMBINATIONS[name])
    ax.set_xticks([1,2,2.5]);ax.set_xlabel('Promoter half-width (kb)')
    ax.set_ylabel('Genes with both marks');ax.set_ylim(bottom=0)
    ax.set_title('Promoter-window sensitivity',loc='left',pad=9)
    ax.axvline(2,color='.65',lw=.6,ls=':',zorder=0)
    ax.legend(loc='upper left',bbox_to_anchor=(-.25,-.23),frameon=False,fontsize=6.1,handlelength=1.7)
    fig.suptitle('Measured HT1080 chromatin reference',x=.04,ha='left',fontsize=11,fontweight='bold',y=.985)
    fig.text(.04,.955,'GSE261834  |  25,461 RNA-ranked genes  |  reference annotation; no LOY contrast',fontsize=7,color='.3')
    stem=FIG/'HT1080_measured_chromatin_reference'
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    canvas=fig.bbox
    for text in fig.findobj(matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            box=text.get_window_extent(renderer)
            assert box.x0>=-1 and box.y0>=-1 and box.x1<=canvas.width+1 and box.y1<=canvas.height+1, text.get_text()
    for ext in ['pdf','svg','png']:
        fig.savefig(stem.with_suffix('.'+ext),dpi=600)
    plt.close(fig)
    caption=(
      '**Measured HT1080 chromatin reference from GSE261834.** '
      '(a) Numbers of genes with at least one unmasked base of a filtered reference peak in the selected promoter. '
      'Eight histone marks and two lamin peak tracks are shown. (b) Selected intersections require both marks in the same promoter window. '
      '(c) All 25,461 RNA-ranked genes are ordered by their eight-mark binary pattern and gene ID, without using RNA effect direction. '
      'Colour indicates the fraction with a detected peak within consecutive bins of up to 50 genes; grey denotes bins with no evaluable promoters. '
      '(d) Co-marked gene counts for promoter half-widths of 1, 2 and 2.5 kb; the dotted line marks the primary definition. '
      'Primary promoters include the TSS base plus 2,000 bp on either side (4,001 bp before clipping/masking), selected using MANE_Select then Ensembl_canonical tags. '
      'Peak filters are minus-log10(q) >=5 for the deposited narrow peaks and >=2 for broad peaks; ENCODE hg38 blacklist bases are subtracted. '
      'Promoters with more than 20% masked bases are unevaluable (467 ranked genes); 24,994 are evaluable. '
      'Counts are descriptive and gene sets overlap. No hypothesis tests, error bars or independent replicate counts are implied. '
      'Source tracks may pool replicates. Co-marking does not establish bivalency within the same cell or nucleosome. '
      'The data describe reference HT1080 cultures, not chromatin changes following LOY.\n')
    cap=ROOT/'outputs/captions/chromatin_reference_atlas.md'
    cap.write_text(caption,encoding='utf-8')
    short_table=sets.loc[list(COMBINATIONS),['n_RNA_universe']].copy()
    short_table.index=[COMBINATIONS[x] for x in short_table.index]
    overview=marks[['detected_any_bp','coverage_ge10pct','coverage_ge50pct']].rename(columns={
        'detected_any_bp':'Any peak overlap','coverage_ge10pct':'>=10% coverage','coverage_ge50pct':'>=50% coverage'})
    candidates=pd.read_csv(OUT/'selected_locus_reference_annotations.csv')
    selected=candidates[candidates.gene_name.isin(['CHRNB2','CHRNA7','MMP1','MMP3','ZMYND8'])]
    methods=f'''# Measured HT1080 chromatin-reference atlas

Completed 14 September 2026. Step 1 only: reference annotation and gene sets. No new GSEA, CAMERA, RNA differential analysis or ATAC integration was run.

## Deliverables

- [Publication figure (PDF)](outputs/figures/chromatin_atlas/HT1080_measured_chromatin_reference.pdf), [PNG](outputs/figures/chromatin_atlas/HT1080_measured_chromatin_reference.png), [SVG](outputs/figures/chromatin_atlas/HT1080_measured_chromatin_reference.svg).
- [Standalone illustrated report](outputs/report/chromatin_reference_atlas.html).
- [Compact promoter matrix: all 25,461 ranked genes](outputs/tables/chromatin_atlas/HT1080_promoter_marks_RNA_universe.csv).
- [Full ranked-gene feature matrix](outputs/tables/chromatin_atlas/HT1080_chromatin_features_RNA_universe.csv), including gene-body and signal summaries.
- [All 78,899 annotated genes](outputs/tables/chromatin_atlas/HT1080_chromatin_features_all_genes.csv.gz), including excluded/unknown features and Y genes as context.
- [Twenty gene-set definitions and sizes](outputs/tables/chromatin_atlas/gene_set_catalogue.csv), [symbol GMT](outputs/tables/chromatin_atlas/HT1080_GSE261834_measured_symbols.gmt), [gene-ID GMT](outputs/tables/chromatin_atlas/HT1080_GSE261834_measured_gene_ids.gmt), [membership](outputs/tables/chromatin_atlas/gene_set_membership_RNA_universe.csv).
- [Window, threshold and alternative-TSS sensitivity](outputs/tables/chromatin_atlas/gene_set_size_sensitivity.csv), [uncertain ranked promoters](outputs/tables/chromatin_atlas/RNA_genes_with_uncertain_promoter_annotation.csv).

## What the atlas contains

Eight standard NEED-seq histone tracks and two lamin peak tracks from GSE261834 were annotated against the frozen GENCODE v49 gene table. All 25,461 original representative gene IDs and symbols were preserved. There are 24,994 evaluable ranked promoters and 467 with more than 20% of the promoter covered by the hg38 blacklist. Unknown calls are retained, never converted to absence.

{overview.to_markdown()}

Selected co-marked promoter sets:

{short_table.to_markdown()}

The two 318-gene sets are separate definitions and must not be assumed to contain the same genes. Broad gene-body coverage sets and individual mark sets are also provided. Fourteen of the 20 sets exceed the original fgsea maximum of 500 genes; six fit the existing 15–500 limits. Sets were not truncated to force eligibility. Before step 2, prespecify size handling for the custom collection and assess the 467 unevaluable promoters as a background sensitivity.

## Frozen definitions and provenance

Inputs, peak thresholds and TSS rules were hashed before annotation. Only RNA gene IDs/names were read to define the universe; effect values did not determine genes, windows or thresholds. The input manifest includes the exact peak/bigWig files, blacklist, GTF, gene annotation and source metadata. Positive labels refer to deposited reference peaks under the chosen filter, not replicate-supported differential chromatin.

GENCODE contains 509,650 transcripts at 404,269 distinct gene-associated TSSs. One transcript per gene was selected by MANE_Select, then Ensembl_canonical, with deterministic fallback rules frozen in configuration (no fallback was needed). Transcript-level coordinates and all distinct-TSS peak calls are retained. The primary window is [TSS0−2000, TSS0+2001) in BED coordinates, after chromosome clipping. Minus-strand TSS0 is GTF transcript end minus one. The sensitivity windows use half-widths 1,000 and 2,500 bp. This symmetric convention includes the TSS base; the old `promoters(upstream=2500,downstream=2500)` definition contained 5,000 rather than 5,001 bp.

Narrow-peak minus-log10(q) must be >=5; broad-peak values must be >=2. The prespecified relaxed comparison uses >=3 and >=1.3, respectively. Filtering uses peak q-values within each source analysis; these are not RNA p-values. Blacklist bases are subtracted from both peaks and queried regions. Any nonzero overlap defines the main positive call; exact covered bases, covered fraction and stricter 10%/50% coverage summaries remain available. Gene-body overlap is reported separately from promoter overlap.

Co-marked sets require two marks in the same selected promoter window. The alternative-TSS analysis also requires both marks in the same TSS window, preventing false co-marking by combining marks at distant alternative promoters. This remains bulk co-marking, not proof of simultaneous marks in one cell or nucleosome.

Signal summaries integrate original bigWig intervals exactly over unmasked regions. Both represented-base mean and mean with unrepresented bases treated as zero are exported, together with represented fraction. Chromosome absence remains missing. All ten files have the expected GRCh38 chr1/chrX lengths and agree on shared sequence lengths. Eight histone tracks are described as RPKM in GEO; absolute values are not comparable between antibodies. **The lamin bigWigs are signed, approximately 20-kb signals, inconsistent with treating them as promoter-resolution RPKM.** Their normalization is unresolved at file level; they are retained as separately labelled coarse context. Lamin gene sets use the deposited peak files, not those coarse signal values, and are not claimed to reconstruct the paper's LAD states.

The publication describes replicate merging and GEO associates the processed histone files with the first sample accession. Exact per-file pooling and original assay-control relationships remain unresolved. This atlas therefore makes no independent-replicate support claim. No K562, protocol-comparison/MNase, imputed ENCODE or other-cell-line tracks enter these sets. Primary source: [GSE261834](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261834), [NEED-seq methods paper](https://academic.oup.com/nar/article/53/2/gkae1317/7958221).

## Interpretation and limitations

The H3K4me1 source file contains 1,905 deposited peaks, of which 293 pass the primary threshold; only 16 ranked promoters overlap a retained peak. The relaxed threshold identifies 92 such promoters. This is a limited reference feature, not evidence that other genes lack enhancer regulation. Distal overlapping H3K4me1/H3K27ac intervals outside all annotated ±2-kb TSS windows are exported only as candidate regulatory regions (minimum residual length 100 bp); no unsupported enhancer-to-gene assignments were made.

The deposited broadPeak files mostly contain short intervals (filtered medians approximately 244–335 bp). The broadPeak format alone does not establish broad heterochromatin domains. Coverage and continuous histone signal are retained to expose this limitation.

Candidate co-marked ranked genes number 260, 318 and 336 at ±1, ±2 and ±2.5 kb; the relaxed peak thresholds give 1,002 and the any-annotated-TSS definition gives 356. This demonstrates substantial threshold dependence. The old 432-gene list covered a different annotation/background; its 209 expressed genes are not directly comparable to 318 in the current ranked universe.

CHRNB2, CHRNA7, MMP1, MMP3 and ZMYND8 have no qualifying peak from the ten primary reference peak tracks at the selected promoter. This does not establish absence of chromatin regulation, and reference cultures can differ from ours. It does mean that Mikkelsen leading-edge membership is not independently corroborated as H3K27me3 promoter occupancy for CHRNB2 by this particular primary peak definition. Inspect quantitative signal, threshold/TSS sensitivity and independent ChIP before drawing a stronger conclusion.

Our ATAC samples are successive passages with uncertain RNA partners. No sample matching or population-level validation was inferred here. These are baseline reference features; measuring LOY-associated histone changes requires chromatin assays in the CTRL and LOY populations.

## Validation and rerun

Validation compares frozen hashes and identities, strand-aware TSS coordinates, independently reconstructed peak/blacklist overlaps and exact basewise bigWig signal summaries. See [validation](provenance/chromatin_atlas/validation.json) and [signal checks](provenance/chromatin_atlas/signal_validation.json). Plotting inputs and complete sensitivity membership are retained.

Run `python run_chromatin_atlas.py` from this analysis directory. It uses Windows Python/R for annotation and the existing WSL ATAC Python environment for bigWigs, then validates and regenerates the report. No package installation or original input modification is required. Original RNA/GSEA and source ATAC outputs remain unchanged.

## Figure caption

{caption}
'''
    (ROOT/'CHROMATIN_REFERENCE_ATLAS.md').write_text(methods,encoding='utf-8')
    png=base64.b64encode(stem.with_suffix('.png').read_bytes()).decode()
    report=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>HT1080 measured chromatin reference</title>
<style>body{{font:16px/1.55 Arial,sans-serif;color:#222;max-width:1150px;margin:35px auto;padding:0 25px}}h1{{font-size:28px}}table{{border-collapse:collapse;font-size:14px}}td,th{{padding:7px 12px;border-bottom:1px solid #ddd;text-align:left}}img{{width:100%;max-width:950px}}a{{color:#1c5f92}}.note{{background:#f4f5f6;padding:16px}}</style>
<h1>Measured HT1080 chromatin reference</h1><p>25,461 RNA-ranked genes · 24,994 evaluable promoters · 20 reference gene sets</p>
<p class="note">Step 1 is complete: this report annotates reference chromatin. No RNA enrichment or LOY chromatin comparison was performed. Replicate pooling remains unresolved. Co-marking does not establish bivalency in individual cells.</p>
<p><a href="../../CHROMATIN_REFERENCE_ATLAS.md">Full methods and interpretation</a> · <a href="../tables/chromatin_atlas/HT1080_promoter_marks_RNA_universe.csv">Promoter matrix</a> · <a href="../tables/chromatin_atlas/HT1080_chromatin_features_RNA_universe.csv">Full feature matrix</a> · <a href="../tables/chromatin_atlas/HT1080_GSE261834_measured_symbols.gmt">Gene sets (GMT)</a> · <a href="../figures/chromatin_atlas/HT1080_measured_chromatin_reference.pdf">Figure PDF</a> · <a href="../figures/chromatin_atlas/HT1080_measured_chromatin_reference.svg">Figure SVG</a></p>
<img src="data:image/png;base64,{png}" alt="Four-panel measured HT1080 reference atlas"><p>{html.escape(caption.replace('**',''))}</p>
<h2>Promoter assignments</h2>{overview.to_html()}<h2>All gene-set definitions</h2>{sets.reset_index().to_html(index=False)}
<h2>Interpretive limits</h2><p>467 masked promoters remain unknown. H3K4me1 calls are sparse. BroadPeak format does not imply broad domains. Lamin signal is coarse and signed; its normalization remains unresolved. Fourteen gene sets exceed the previous GSEA size cap of 500, requiring a prespecified decision for the next testing stage. All windows, threshold comparisons and gene memberships are available in the linked methods.</p></html>'''
    (ROOT/'outputs/report/chromatin_reference_atlas.html').write_text(report,encoding='utf-8')
    print('Figure, caption and reference report written.',flush=True)


if __name__=='__main__':main()
