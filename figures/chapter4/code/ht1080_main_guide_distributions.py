"""Supporting guide distributions for the primary drugZ candidate figure.

All guide values are the saved MAGeCK endpoint contrast. No significance tests
or residual model are fitted. Guide ticks have equal height and no x jitter.
"""
import hashlib
import json
import make_ht1080_candidates as C
import numpy as np
import pandas as pd

DATA = C.OUT / 'source_data/ht1080_main_guide_distributions'
SCREEN = C.THESIS.parent / 'crispr_screen/ht1080_loy_screen_2'
SAMPLES = ['D0_CTRL', 'D0_LOY', 'POST_CTRL_REP1', 'POST_CTRL_REP2', 'POST_LOY_REP1', 'POST_LOY_REP2']
XLIM = (-8.5, 5.8)


def prepare(named, gametologs):
    DATA.mkdir(exist_ok=True)
    paths = {'counts':SCREEN/'count/ht1080_screen2_R1.count.txt',
             'mageck_guides':SCREEN/'test/POST_LOY_vs_CTRL.sgrna_summary.txt',
             'drugz_genes':SCREEN/'results/drugZ/LOY_vs_CTRL_unfiltered.txt'}
    counts=pd.read_csv(paths['counts'],sep='\t')
    mg=pd.read_csv(paths['mageck_guides'],sep='\t')
    joined=counts.merge(mg[['sgrna','Gene','LFC','FDR']],left_on=['sgRNA','Gene'],
                       right_on=['sgrna','Gene'],how='left',validate='one_to_one').drop(columns='sgrna')
    assert len(joined)==77441 and joined.LFC.notna().all()
    joined['endpoint_read_total']=joined[SAMPLES[2:]].sum(axis=1)
    joined['endpoint_detected_samples']=joined[SAMPLES[2:]].gt(0).sum(axis=1)
    joined['all_samples_zero']=joined[SAMPLES].sum(axis=1).eq(0)
    joined['lfc_displayed']=joined.endpoint_read_total.gt(0)&np.isfinite(joined.LFC)
    genes=set(named.GENE)|set(gametologs)
    selected=joined.loc[joined.Gene.isin(genes)].sort_values(['Gene','sgRNA']).copy()
    selected['guide_number']=selected.groupby('Gene').cumcount()+1
    selected['sequence_5prime_to_3prime']=selected.sgRNA.str.rsplit('_',n=1).str[-1]
    selected['drugz_normZ']=selected.Gene.map(C.INDEXED.normZ)
    enriched=set(named.loc[named.direction.eq('enriched'),'GENE'])
    selected['drugz_fdr_direction']=np.where(selected.Gene.isin(enriched),'enrichment','depletion')
    selected['drugz_directional_fdr']=[C.INDEXED.loc[r.Gene,'fdr_supp' if r.Gene in enriched else 'fdr_synth']
                                       for r in selected.itertuples()]
    selected=selected.rename(columns={'LFC':'mageck_reported_log2FC','FDR':'mageck_guide_FDR'})
    assert len(genes)==34 and len(selected)==136 and selected.lfc_displayed.sum()==134
    assert selected.groupby('Gene').size().eq(4).all()
    selected.to_csv(DATA/'selected_guides.csv',index=False)

    ntc=joined.loc[joined.Gene.eq('Non-Targeting-Control')].copy()
    assert len(ntc)==1000 and ntc.lfc_displayed.sum()==999
    ntc.to_csv(DATA/'non_targeting_guides.csv',index=False)
    values=ntc.loc[ntc.lfc_displayed,'LFC'].to_numpy()
    bandwidth=.9*min(values.std(ddof=1),np.subtract(*np.percentile(values,[75,25]))/1.34)*len(values)**(-.2)
    grid=np.linspace(*XLIM,2861)
    density=np.exp(-.5*((grid[:,None]-values[None,:])/bandwidth)**2).mean(axis=1)/(bandwidth*np.sqrt(2*np.pi))
    pd.DataFrame({'guide_log2FC':grid,'kde_density':density,'display_density_peak_scaled':density/density.max()}).to_csv(DATA/'ntc_density.csv',index=False)
    manifest={'inputs':{k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()},
              'panel_b_genes':gametologs,'panel_c_genes':named.GENE.tolist(),
              'unique_genes':34,'unique_guides':136,'unique_guides_displayed':134,
              'duplicated_gene_across_panels':'DDX3X',
              'ntc_guides':1000,'ntc_with_endpoint_reads':999,'ntc_zero_endpoint_reads_excluded':1,
              'density_method':'Gaussian KDE; 0.9 * min(sample SD, IQR/1.34) * n^(-1/5)',
              'density_bandwidth':float(bandwidth),'x_limits':list(XLIM),
              'tick_height':'equal; no biological meaning','horizontal_jitter':False,
              'interpretation':'Observed guide fold changes and empirical NTC reference, not fitted residuals or new significance tests.',
              'normalisation_caveat':'Saved MAGeCK used non-targeting controls; their centring is not independent QC.'}
    (DATA/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return selected,grid,density/density.max(),manifest


def draw(ax, genes, selected, grid, curve, panel, column, marks, *, saga=(), slots=None, font_size=7.1):
    slots=slots or len(genes)
    for j,gene in enumerate(genes):
        base=slots-1-j
        data=selected.loc[selected.Gene.eq(gene)]
        ax.fill_between(grid,base+.06,base+.06+curve*.76,color='#EAEAEA',lw=0,zorder=0)
        ax.plot(grid,base+.06+curve*.76,color='#A7A7A7',lw=.6,zorder=1)
        ax.hlines(base+.06,*XLIM,color='#DEDEDE',lw=.4,zorder=1)
        ax.vlines(0,base+.06,base+.86,color='#C3C3C3',lw=.45,ls=(0,(2,2)),zorder=1)
        for r in data.itertuples():
            if not r.lfc_displayed: continue
            ax.vlines(r.mageck_reported_log2FC,base+.06,base+.86,color='#333333',lw=.85,zorder=3)
            marks.append({'panel':panel,'column':column,'row':j,'Gene':gene,'sgRNA':r.sgRNA,
                          'x':r.mageck_reported_log2FC,'height':.8})
        label=gene+('*' if not data.lfc_displayed.all() else '')
        ax.text(-.04,base+.47,label,transform=ax.get_yaxis_transform(),va='center',ha='right',
                fontsize=font_size,fontweight='semibold' if gene in saga else 'normal')
        ax.text(1.045,base+.47,C.fdr_text(data.drugz_directional_fdr.iloc[0]),
                transform=ax.get_yaxis_transform(),va='center',fontsize=font_size-.2)
    ax.text(1.045,1.018,'drugZ FDR',transform=ax.transAxes,fontsize=6.8,ha='left',va='bottom')
    ax.set(xlim=XLIM,ylim=(-.08,slots+.02),yticks=[])
    ax.spines['left'].set_visible(False)
    ax.set_xticks([-8,-4,0,4])
    ax.tick_params(axis='x',labelsize=7,length=2,pad=2)
    ax.set_xlabel('Guide log₂FC (LOY / CTRL)',fontsize=7.2,labelpad=3)


def finish(marks,manifest):
    table=pd.DataFrame(marks)
    assert len(table)==137 and len(table.sgRNA.unique())==134
    table.to_csv(DATA/'displayed_tick_positions.csv',index=False)
    manifest['displayed_tick_count']=137
    manifest['render_helper_sha256']=hashlib.sha256(__import__('pathlib').Path(__file__).read_bytes()).hexdigest()
    (DATA/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
