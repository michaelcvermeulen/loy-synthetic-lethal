"""Audit and visualize the nicotinic RNA pattern using existing results only."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
TABLES=ROOT/'outputs/tables/nicotinic_receptors'
FIGURES=ROOT/'outputs/figures/nicotinic_receptors'
PROV=ROOT/'provenance/nicotinic_receptors'
UP=['CHRNA6','CHRNA7','CHRNB2','CHRNB4','CHRNE']
FAMILY=[f'CHRNA{i}' for i in [1,2,3,4,5,6,7,9,10]]+[f'CHRNB{i}' for i in [1,2,3,4]]+['CHRND','CHRNE','CHRNG','CHRFAM7A']
BLUE='#0072B2'; ORANGE='#D55E00'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return pd.read_csv(path,float_precision='round_trip')


def clean_axis(ax):
    for side in ['top','right']:
        ax.spines[side].set_visible(False)
    for side in ['bottom','left']:
        ax.spines[side].set_linewidth(.5)
        ax.spines[side].set_color('#6C7780')
    ax.tick_params(labelsize=5.5,length=2,width=.5)


def main():
    for folder in [FIGURES,PROV]:
        folder.mkdir(parents=True,exist_ok=True)
    inputs=[ROOT/p for p in ['outputs/models/01_ingest.rds','outputs/models/02_expression.rds',
        'outputs/tables/HT1080_LOY_vs_sgCTRL_all_gene_effects.csv','outputs/tables/fgsea_all_collections.csv',
        'outputs/tables/camera_all_collections.csv','outputs/tables/camera_estimated_correlation_sensitivity.csv',
        'outputs/tables/network_modules/UP_association_edges.csv','outputs/tables/pathway_resource_coverage.csv',
        'outputs/tables/nicotinic_receptors/sample_expression.csv','outputs/tables/nicotinic_receptors/gene_effects.csv']]
    frozen={p.relative_to(ROOT).as_posix():sha(p) for p in inputs}
    samples=read(inputs[-2]); effects=read(inputs[-1]).set_index('gene_name')
    original=read(inputs[2]).set_index('gene_id')
    assert len(samples)==26*12 and not samples.duplicated(['gene_id','sample_key']).any()
    assert (samples.groupby('sample_key').size()==26).all()
    assert samples.raw_count.ge(0).all() and samples.raw_count.eq(np.floor(samples.raw_count)).all()
    np.testing.assert_allclose(effects.mean_difference,original.loc[effects.gene_id,'mean_difference'],rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(samples.groupby('gene_name').CPM.mean().loc[effects.index],effects.mean_CPM,rtol=1e-12,atol=1e-12)
    differences=samples.pivot(index=['gene_name','batch','passage'],columns='treatment',values='logCPM')
    differences['difference']=differences.LOY-differences.CTRL
    differences=differences.reset_index()
    for row in differences.itertuples():
        saved=effects.loc[row.gene_name,f'{row.batch}__passage_{row.passage}']
        assert abs(row.difference-saved)<1e-10,(row.gene_name,row.batch,row.passage)
    differences.to_csv(TABLES/'matched_differences.csv',index=False,float_format='%.17g')
    abundance=samples.groupby(['gene_name','batch','treatment']).agg(mean_CPM=('CPM','mean'),
        min_reads=('raw_count','min'),max_reads=('raw_count','max'),total_reads=('raw_count','sum')).reset_index()
    abundance.to_csv(TABLES/'period_abundance.csv',index=False,float_format='%.17g')
    effects['total_assigned_reads']=samples.groupby('gene_name').raw_count.sum()
    effects['zero_count_libraries']=samples.groupby('gene_name').raw_count.agg(lambda v:int(v.eq(0).sum()))
    effects['all_zero_counts']=effects.total_assigned_reads.eq(0)
    effects['module_UP_01']=effects.index.isin(UP)
    effects['family_member']=effects.index.isin(FAMILY)
    effects.to_csv(TABLES/'gene_summary.csv',float_format='%.17g')
    fg=read(inputs[3])
    # Nicotinate/nicotinamide pathways are vitamin metabolism, not this receptor family.
    mask=fg.pathway.str.contains('NICOTINIC|ACETYLCHOLINE|CHOLINERGIC|NICOTINE',regex=True)
    focused=fg[mask].copy()
    for path,prefix in [(inputs[4],'camera_fixed'),(inputs[5],'camera_estimated')]:
        cam=read(path)[['pathway','contrast','PValue','Direction','BH_across_collections']]
        cam=cam.rename(columns={c:f'{prefix}_{c}' for c in ['PValue','Direction','BH_across_collections']})
        focused=focused.merge(cam,on=['pathway','contrast'],how='left',validate='one_to_one')
    focused.to_csv(TABLES/'original_cholinergic_enrichments.csv',index=False,float_format='%.17g')
    coverage=read(inputs[7]); coverage=coverage[coverage.pathway.str.contains('NICOTINIC|ACETYLCHOLINE|CHOLINERGIC|NICOTINE')]
    coverage.to_csv(TABLES/'cholinergic_set_coverage.csv',index=False)
    shared=fg[fg.contrast.eq('combined')&fg.BH_across_collections.lt(.05)].copy()
    shared['UP_module_leading_edge_genes']=shared.leadingEdge.map(lambda s:';'.join(sorted(set(s.split(';'))&set(UP))))
    shared=shared[shared.UP_module_leading_edge_genes.ne('')]
    shared.to_csv(TABLES/'UP_module_links_to_enriched_sets.csv',index=False,float_format='%.17g')
    edges=read(inputs[6]); edges=edges[edges.gene_a.isin(UP)&edges.gene_b.isin(UP)]
    assert len(edges)==10 and edges.physical_supported.all()
    edges.to_csv(TABLES/'UP_module_STRING_edges.csv',index=False)

    plt.rcParams.update({'font.family':'Arial','font.size':6,'pdf.fonttype':42,'ps.fonttype':42,
                         'svg.fonttype':'none','svg.hashsalt':'HT1080_nicotinic_42'})
    fig=plt.figure(figsize=(183/25.4,220/25.4),facecolor='white')
    fig.text(.04,.978,'Nicotinic receptor subunit RNA changes in HT1080 LOY',fontsize=9,fontweight='bold',va='top')
    fig.text(.04,.951,'Selective changes across the family; five subunits form the UP_01 association module',fontsize=6.5,color='#53616A')
    fig.text(.02,.918,'a',fontsize=8,fontweight='bold');fig.text(.07,.918,'Family-wide effects',fontsize=7,fontweight='bold')
    ax=fig.add_axes([.145,.423,.34,.465]); selected=effects.loc[FAMILY]
    ax.axvline(0,color='#929DA4',lw=.55,ls='--'); ax.set_xlim(-1.8,3.35); ax.set_ylim(len(FAMILY)-.5,-.5)
    for i,(gene,row) in enumerate(selected.iterrows()):
        if i%2==0: ax.axhspan(i-.5,i+.5,color='#F5F7F8',zorder=0)
        if row.all_zero_counts:
            ax.text(.05,i,'no assigned reads',fontsize=5.4,color='#777F85',va='center')
        else:
            alpha=1 if row.primary_eligible else .38
            ax.plot([row.AUG_2025,row.MAY_2026],[i,i],color='#C3CBD1',lw=.7,alpha=alpha)
            ax.scatter(row.AUG_2025,i,s=13,c=BLUE,marker='o',alpha=alpha,zorder=3)
            ax.scatter(row.MAY_2026,i,s=15,c=ORANGE,marker='^',alpha=alpha,zorder=3)
            ax.scatter(row.mean_difference,i,s=7,c='#263640',marker='D',alpha=alpha,zorder=4)
        ax.text(3.63,i,f'{row.mean_CPM:.2f}',fontsize=5.7,va='center',ha='right',clip_on=False,color='#53616A')
    ax.set_yticks(np.arange(len(FAMILY)));ax.set_yticklabels([g+('†' if not effects.loc[g,'primary_eligible'] else '') for g in FAMILY],fontsize=6,fontstyle='italic')
    for tick,g in zip(ax.get_yticklabels(),FAMILY):
        if g in UP:tick.set_fontweight('bold')
        if not effects.loc[g,'primary_eligible']:tick.set_color('#8A9298')
    ax.set_xlabel('LOY − sgCTRL (log2 CPM difference)',fontsize=6)
    clean_axis(ax);ax.tick_params(axis='y',length=0)
    fig.text(.503,.897,'Mean CPM',fontsize=5.7,ha='right',fontweight='bold')
    h=[Line2D([0],[0],marker=m,color='none',markerfacecolor=c,markeredgecolor=c,markersize=3,label=l) for m,c,l in [('o',BLUE,'August'),('^',ORANGE,'May'),('D','#263640','Combined')]]
    fig.legend(handles=h,loc='upper left',bbox_to_anchor=(.115,.393),ncol=3,frameon=False,fontsize=5.7,columnspacing=1.1,handletextpad=.25)
    fig.text(.07,.365,'† Below original expression filter. Bold: UP_01 member.',fontsize=5.3,color='#53616A')

    fig.text(.56,.918,'b',fontsize=8,fontweight='bold');fig.text(.595,.918,'All six passage/day matches',fontsize=7,fontweight='bold')
    bx=fig.add_axes([.635,.659,.29,.22]); matrix=effects.loc[UP,[f'{b}__passage_{i}' for b in ['AUG_2025','MAY_2026'] for i in [1,2,3]]].to_numpy()
    image=bx.imshow(matrix,cmap='RdBu_r',vmin=-3.6,vmax=3.6,aspect='auto')
    bx.set_xticks(range(6),['1','2','3','1','2','3'],fontsize=6);bx.set_yticks(range(5),UP,fontsize=6,fontstyle='italic')
    bx.tick_params(length=0);bx.axvline(2.5,color='white',lw=2)
    for i in range(5):
        for j in range(6):bx.text(j,i,f'{matrix[i,j]:+.2f}',ha='center',va='center',fontsize=5.4,color='white' if abs(matrix[i,j])>2.1 else '#25343E')
    for spine in bx.spines.values():spine.set_visible(False)
    bx.text(1,-.85,'August 2025',fontsize=6,ha='center');bx.text(4,-.85,'May 2026',fontsize=6,ha='center')
    bx.set_xlabel('Matched passage number',fontsize=5.7)
    cb=fig.colorbar(image,cax=fig.add_axes([.945,.659,.012,.22]),ticks=[-3,0,3]);cb.ax.tick_params(labelsize=5.3,length=2,width=.4)
    cb.set_label('LOY − sgCTRL (log2 CPM)',fontsize=5.7);cb.outline.set_linewidth(.4)
    fig.text(.575,.593,'Expression pattern',fontsize=6.3,fontweight='bold')
    fig.text(.575,.572,'CHRNA7 and CHRNB2 have stronger read support.\nCHRNA6 and CHRNB4 are low-abundance genes;\ntheir effect magnitudes depend on the prior count.',fontsize=5.7,linespacing=1.5,va='top')
    fig.text(.575,.502,'Pathway context',fontsize=6.3,fontweight='bold')
    fig.text(.575,.48,'Broader cholinergic GSEA: no global-BH < 0.05 hits.\nThe UP_01 annotation is highly coherent (q = 5.51e−12).\nCHRNB2 is in both Mikkelsen H3K27me3 leading edges.\nThese statistics address different questions.',fontsize=5.7,linespacing=1.5,va='top')

    fig.text(.02,.326,'c',fontsize=8,fontweight='bold');fig.text(.055,.326,'Original sample abundance in the five UP-module genes',fontsize=7,fontweight='bold')
    for index,gene in enumerate(UP):
        cx=fig.add_axes([.055+index*.19,.105,.147,.17]); sub=samples[samples.gene_name.eq(gene)]
        for batch,start in [('AUG_2025',0),('MAY_2026',3)]:
            for passage,marker in [(1,'o'),(2,'s'),(3,'^')]:
                pair=sub[sub.batch.eq(batch)&sub.passage.eq(passage)].set_index('treatment').loc[['CTRL','LOY']]
                offset=(passage-2)*.09
                cx.plot([start+offset,start+1+offset],pair.CPM,color='#ACB6BE',lw=.65,zorder=1)
                cx.scatter([start+offset,start+1+offset],pair.CPM,s=14,marker=marker,c=[BLUE,ORANGE],edgecolor='white',linewidth=.25,zorder=3)
        cx.set_title(gene,fontstyle='italic',fontsize=7,fontweight='bold',pad=7)
        cx.set_xticks([0,1,3,4],['C','L','C','L']);cx.set_xlim(-.5,4.5)
        cx.set_ylim(-.05*sub.CPM.max(),1.08*sub.CPM.max())
        cx.text(.5,-.22,'Aug',transform=cx.get_xaxis_transform(),ha='center',fontsize=5.5)
        cx.text(3.5,-.22,'May',transform=cx.get_xaxis_transform(),ha='center',fontsize=5.5)
        clean_axis(cx)
        if index==0:cx.set_ylabel('TMM CPM',fontsize=6)
    fig.text(.055,.052,'C, sgCTRL; L, LOY. Lines link matched passages; circles, squares and triangles denote passages 1, 2 and 3.',fontsize=5.5,color='#53616A')
    fig.text(.055,.033,'Same established populations sampled again. RNA abundance does not measure receptor assembly, surface expression or channel activity.',fontsize=5.5,color='#53616A')
    stem='HT1080_nicotinic_receptor_followup'
    for extension in ['pdf','svg','png']:fig.savefig(FIGURES/f'{stem}.{extension}',dpi=600,facecolor='white')
    plt.close(fig)

    def table(frame,columns,headers):
        lines=['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(columns))+'|']
        for gene,row in frame.iterrows():
            vals=[gene if c=='gene' else (f'{row[c]:.2f}' if isinstance(row[c],(float,np.floating)) else str(row[c])) for c in columns]
            lines.append('| '+' | '.join(vals)+' |')
        return '\n'.join(lines)
    top_table=table(effects.loc[UP],['gene','mean_CPM','AUG_2025','MAY_2026','mean_difference','matches_agreeing_with_mean'],
        ['Gene','Mean CPM','August Δlog2 CPM','May Δlog2 CPM','Combined Δlog2 CPM','Same-direction matches / 6'])
    report='''# Nicotinic receptor follow-up in HT1080 LOY

The current data support selective remodeling of nicotinic-subunit RNA abundance, with clearer increases in CHRNA7 and CHRNB2, low-abundance increases in CHRNA6 and CHRNB4, and decreases in CHRNA1, CHRNA9 and CHRNA10. They do not establish increased receptor activity or one five-subunit receptor complex.

'''+top_table+'''

The five-member UP_01 module has MCODE score 5, ten association edges and ten edges above the STRING physical-score threshold of 700. Its top annotation, Reactome acetylcholine binding and downstream events, has raw P = 8.47e−17 and global annotation BH q = 5.51e−12. This quantifies functional coherence of the selected genes. It is not a test of RNA change or evidence that all five products assemble together in HT1080. STRING physical edges can encode co-complex evidence and transferred evidence; pair-specific experimental provenance is not resolved by the cached combined-score files. [STRING scoring documentation](https://version-12-0.string-db.org/help/scores/).

## Abundance and temporal consistency

CHRNA7 and CHRNB2 have appreciably more assigned reads than CHRNA6 and CHRNB4. In May, CHRNA7 averages 2.12 CPM in sgCTRL and 5.99 in LOY; CHRNB2 averages 3.39 and 10.80 CPM. Both are higher in five of six passage/day comparisons and remain positive after leaving out any one match. Their effect estimates change by about 0.01 log2 units with the smaller prior count or alternative RLE normalization.

CHRNA6 averages 0.41 CPM overall and includes zero-read controls. CHRNB4 averages 0.28 CPM and is higher in all six matches. Its May controls have 0–4 reads, while LOY samples have 12–22. All five genes remain eligible at the saved 50- and 100-total-read sensitivity thresholds. However, reducing the logCPM prior count from 2 to 0.25 changes the combined CHRNA6 estimate from +1.18 to +1.95, and CHRNB4 from +2.12 to +3.29. Direction is more stable than the exact magnitude at these low counts. Prior-count effects are not confidence intervals.

The mean module contrast is +0.60 in August and +1.75 in May. This is descriptive temporal variation. The larger CHRNB4 May contrast partly reflects lower sgCTRL abundance: its mean LOY abundance is 0.52 CPM in August and 0.40 in May, while controls fall from 0.17 to 0.03 CPM. A larger contrast therefore does not imply increasing LOY abundance over time. The same established populations were sampled in both periods; the six matches are not six independently generated populations.

## The rest of the receptor family changes in different directions

CHRNA1 is lower by −1.33 log2 CPM and CHRNA9 by −1.21, both in all six matches. CHRNA10 is lower by −0.62, with a stronger August effect and near-zero May mean. CHRNA3 and CHRNA5 show small positive combined means but opposite period signs. CHRNA4 is below the expression filter, and CHRNB3 has no assigned reads. The full family, including excluded low-expression genes, is retained in the figure and source tables. All-zero genes are displayed as having no assigned reads; tiny floating-point residuals in the original effect table are not interpreted biologically.

## What could the subunits assemble into?

| Candidate assembly | Relevant observations | Interpretation |
|---|---|---|
| α7-containing receptor | CHRNA7 +0.93; CHRFAM7A +0.51; RIC3 and TMEM35A below the expression filter | A testable receptor hypothesis, with assembly and gene-specific measurement unresolved. |
| α6β4-containing receptor | CHRNA6 +1.18 and CHRNB4 +2.12, both low abundance | A plausible pairing; not evidence that a functional channel exists in these cells. |
| β2-containing receptor | CHRNB2 +1.13; CHRNA6 increases, while CHRNA4 is very low and CHRNB3 has no reads | The RNA data do not identify the exact β2-containing receptor subtype. |
| Muscle-type receptor | CHRNE +0.52, CHRNA1 −1.33, CHRNB1 +0.21 and CHRND +0.43 | The muscle-type subunits do not show coordinated induction. |

α7 can form a homopentamer ([structural study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8135066/)); human α6β4 has been structurally characterized as two α6 plus three β4 subunits ([Su et al., 2025](https://www.nature.com/articles/s41467-025-58333-0)). CHRNE is the epsilon component of adult muscle-type receptors containing α1, β1 and δ ([NCBI RefSeq summary](https://www.ncbi.nlm.nih.gov/gene/1145/)). These are distinct known assemblies, so the five-node MCODE clique should not be drawn as one demonstrated heteropentamer. Gene-level CPM is not a measurement of protein abundance or subunit stoichiometry.

CHRFAM7A is a partial duplication of CHRNA7. Coexpression suppressed α7-mediated currents in heterologous experiments ([Araud et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3162115/)). Here CHRFAM7A also increases modestly (+0.51; 0.42 CPM). We cannot infer its net functional effect in HT1080, or resolve shared versus gene-specific reads from this count matrix. Gene-specific exon/junction inspection and orthogonal assays are needed before treating the CHRNA7 signal as a validated α7 change. Low RIC3/TMEM35A RNA is further context for assembly uncertainty, not proof that assembly is impossible; the NACHO/TMEM35A requirement was demonstrated in a different experimental setting ([Gu et al., 2016](https://pubmed.ncbi.nlm.nih.gov/26875622/)).

## Reconciliation with GSEA and the chromatin signatures

The original broader cholinergic sets are not globally significant: acetylcholine-gated channel complex has NES +1.20, global BH q = 0.620; cholinergic synaptic transmission has NES +1.11, q = 0.693; response to nicotine has q = 0.804; response to acetylcholine has q = 0.908. Fixed-correlation CAMERA also fails global q < 0.05 for these sets. The module's highly significant annotation does not override these whole-ranking tests.

The specific Reactome annotation used for UP_01 has only 11 represented ranking genes and was below the original GSEA minimum of 15. It was therefore not tested in that GSEA analysis; it is not a failed enrichment test. The network annotation used a separate, documented size range and annotation-testing family. No thresholds have been relaxed in this follow-up.

CHRNB2 occurs in the leading edges of both enriched Mikkelsen H3K27me3 source signatures and both CNS axonogenesis GO sets. It is the only UP_01 gene in the leading edges of any combined global-BH < 0.05 set. This links one receptor-subunit RNA change to the previously observed developmental/chromatin-associated pattern, but does not establish H3K27me3 loss at CHRNB2 in HT1080 or chromatin regulation of the whole module.

CHAT is very low and below the expression filter, SLC5A7 has no assigned reads, and SLC18A3 is below the filter with opposite period effects. This analysis supplies no positive evidence for an induced endogenous acetylcholine-production/release program. It does not measure ligand availability in the culture medium or exclude alternative sources of acetylcholine.

## Practical follow-up

1. Validate CHRNA7 and CHRNB2 RNA first, with CHRNA9 as a decreasing family comparator. Include CHRFAM7A-specific measurements and gene-specific CHRNA7 coverage. Validate CHRNA6/CHRNB4 with assays that work at their observed low abundance.
2. Establish surface receptor expression and an acetylcholine-evoked response before assigning channel activity. Combine functional readouts with genetic perturbation of CHRNA7 or CHRNB2 and suitable receptor-selectivity controls; a single drug response would not identify the receptor subtype.
3. If the RNA pattern and receptor function are confirmed, test whether the relevant subtype contributes to the LOY-associated phenotype, then examine the CHRNB2 chromatin hypothesis with locus-specific chromatin measurements. Independent LOY derivations or addback are required for stronger causal attribution.

## Deliverables and reproducibility

[Publication PDF](outputs/figures/nicotinic_receptors/HT1080_nicotinic_receptor_followup.pdf) · [600 dpi PNG](outputs/figures/nicotinic_receptors/HT1080_nicotinic_receptor_followup.png) · [Editable SVG](outputs/figures/nicotinic_receptors/HT1080_nicotinic_receptor_followup.svg) · [Figure caption](outputs/captions/nicotinic_receptor_followup.md).

Source data: [gene summary](outputs/tables/nicotinic_receptors/gene_summary.csv), [all sample counts and expression](outputs/tables/nicotinic_receptors/sample_expression.csv), [matched differences](outputs/tables/nicotinic_receptors/matched_differences.csv), [original enrichment results](outputs/tables/nicotinic_receptors/original_cholinergic_enrichments.csv), [gene-set coverage](outputs/tables/nicotinic_receptors/cholinergic_set_coverage.csv), and [links to enriched signatures](outputs/tables/nicotinic_receptors/UP_module_links_to_enriched_sets.csv).

Run `extract_nicotinic_followup.R` and then `nicotinic_followup.py` from any working directory. The frozen 01_ingest/02_expression models, original enrichment results and STRING tables are reused. Counts are assigned reads; no count files, models, gene-level tests, normalization or original enrichment statistics are changed. This follow-up is descriptive and post hoc.
'''
    md=ROOT/'NICOTINIC_RECEPTOR_FOLLOWUP.md';md.write_text(report,encoding='utf-8')
    caption=('**Selective nicotinic receptor subunit RNA changes in HT1080 LOY.** '
        'a, All 16 canonical human nicotinic receptor subunit genes and CHRFAM7A. Blue circles and orange triangles show '
        'August and May mean passage-matched LOY-minus-sgCTRL log2 CPM differences; black diamonds show their equally weighted mean. '
        'Bold labels identify UP_01 members; daggers and faded points mark genes below the original expression filter. '
        'CHRNB3 has no assigned reads. Mean CPM is averaged across the 12 core libraries. '
        'b, The five UP_01 members across all six passage/day matches, with exact log2 CPM differences shown in each cell. '
        'c, Original TMM CPM values in the 12 libraries for each UP_01 gene. Lines join within-period matched passages; '
        'blue denotes sgCTRL and orange LOY, with circles, squares and triangles denoting passages 1–3. Axes use gene-specific linear CPM ranges. '
        'All abundance values come from the saved normalization; logCPM uses prior count 2. CHRNA6 and CHRNB4 are low-abundance genes '
        'whose effect magnitudes depend on that prior. No P values or confidence intervals for RNA change are inferred from these repeated '
        'samples of the same established populations. The annotation q value and GSEA results are unchanged original statistics from distinct '
        'testing families. Source H3K27me3 membership does not measure chromatin marks in HT1080. RNA abundance does not determine receptor assembly or activity.')
    caption_path=ROOT/'outputs/captions/nicotinic_receptor_followup.md'
    caption_path.write_text(caption+'\n',encoding='utf-8')
    pdf=PdfReader(FIGURES/f'{stem}.pdf'); assert len(pdf.pages)==1
    assert abs(float(pdf.pages[0].mediabox.width)*25.4/72-183)<.01
    assert abs(float(pdf.pages[0].mediabox.height)*25.4/72-220)<.01
    text=pdf.pages[0].extract_text(); assert all(g in text for g in FAMILY)
    for ref in pdf.pages[0]['/Resources']['/Font'].values():
        font=ref.get_object()
        for item in font.get('/DescendantFonts',[font]):
            f=item.get_object();desc=f['/FontDescriptor'].get_object()
            assert 'Arial' in str(f.get('/BaseFont')) and any(k in desc for k in ['/FontFile','/FontFile2','/FontFile3'])
    with Image.open(FIGURES/f'{stem}.png') as im:
        assert abs(im.width-183/25.4*600)<1.1 and abs(im.height-220/25.4*600)<1.1
        assert all(abs(v-600)<1 for v in im.info['dpi'])
    assert frozen=={p:sha(ROOT/p) for p in frozen}
    outputs=list(TABLES.glob('*'))+list(FIGURES.glob('*'))+[md,caption_path]
    result={'created_utc':datetime.now(timezone.utc).isoformat(),'genes':26,'family_genes':17,'samples':12,
        'source_inputs_unchanged':frozen,'new_tests_performed':False,'script_sha256':sha(__file__),
        'outputs_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in outputs},'visual_review':'pending'}
    (PROV/'completion.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'genes':26,'sample_observations':len(samples),'module_edges':len(edges),
                     'links_to_enriched_sets':len(shared),'report':str(md),'figure':stem},indent=2))


if __name__=='__main__':
    main()
