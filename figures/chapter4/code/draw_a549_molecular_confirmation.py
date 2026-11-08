"""A549 Y DNA/RNA inset from frozen Chapter 3 characterisation data."""
from pathlib import Path
import hashlib
import json
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent / '.figure_python_deps'))
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle, Polygon, ConnectionPatch

CH3 = ROOT / 'figures/chapter3/source_data'
DATA = ROOT / 'figures/chapter4/source_data/a549_molecular_inset'
GENES = ['RPS4Y1', 'ZFY', 'USP9Y', 'DDX3Y', 'UTY', 'NLGN4Y', 'KDM5D', 'EIF1AY']
REGIONS = ['MSY Yp', 'MSY Yq (euchromatic)']
SOURCES = {
    'dna': CH3 / 'a549_lpwgs_confirmation/regional_medians_verified.csv',
    'dna_bins': CH3 / 'a549_lpwgs_confirmation/qdnaseq_bins.csv',
    'rna': CH3 / 'shared_sex_chromosome/A549/gene_period_values.csv',
    'rna_cpm': CH3 / 'shared_sex_chromosome/A549_CPM.csv',
    'rna_samples': CH3 / 'shared_sex_chromosome/A549_samples.csv',
    'gene_coordinates': CH3 / 'a549_y_expression_confirmation/gene_mapping_plotted.csv',
    'cytobands': CH3 / 'ht1080_y_expression_confirmation/T2T_chrY_cytobands.csv',
}
BLUE, RED, DARK = '#3C5488', '#E64B35', '#222222'


def prepare():
    """Verify condition means and regional medians, then freeze plotted values."""
    DATA.mkdir(exist_ok=True)
    dna = pd.read_csv(SOURCES['dna'])
    dna = dna.loc[dna.region.isin(REGIONS)].copy()
    dna['analysis_floor'] = dna.median_relative_cn.le(.001000001)
    bands = pd.read_csv(SOURCES['cytobands'])
    cen = bands.loc[bands.stain.eq('acen')]
    windows = [(2.458320, cen.start.min()/1e6),
               (cen.end.max()/1e6, bands.loc[bands.band.eq('q12'), 'start'].iloc[0]/1e6)]
    bins = pd.read_csv(SOURCES['dna_bins'])
    bins = bins.loc[bins.use & bins.log2_seg.notna() & bins.chromosome.eq('chrY')]
    assert set(bins['sample']) == {'A549-CTRL2', 'A549-LOY'}
    for region, (start, end) in zip(REGIONS, windows):
        for condition in ['sgCTRL', 'LOY']:
            selected = bins.loc[bins.label.eq(condition) & bins.midpoint.ge(start) & bins.midpoint.lt(end)]
            row = dna.loc[dna.region.eq(region) & dna.condition.eq(condition)].iloc[0]
            assert len(selected) == row.n_bins
            np.testing.assert_allclose(np.median(2**selected.log2), row.median_relative_cn, atol=1e-12)
    source_rna = pd.read_csv(SOURCES['rna']).set_index('gene_name').loc[GENES]
    assert source_rna.batch.eq('AUG_2025').all()
    cpm = pd.read_csv(SOURCES['rna_cpm']).set_index('gene_id').loc[source_rna.gene_id]
    metadata = pd.read_csv(SOURCES['rna_samples'])
    means = {}
    for group in ['CTRL', 'LOY']:
        samples = metadata.loc[metadata.treatment.eq(group), 'sample_key'].tolist()
        assert len(samples) == 4
        means[group] = cpm[samples].mean(axis=1).to_numpy()
        np.testing.assert_allclose(means[group], source_rna[group], rtol=1e-12)
    ratios = 100 * means['LOY'] / means['CTRL']
    np.testing.assert_allclose(ratios, source_rna.percent, rtol=1e-12)
    rna = pd.DataFrame({'gene_name': GENES, 'sgCTRL': 100., 'LOY': ratios,
                        'CTRL_mean_CPM': means['CTRL'], 'LOY_mean_CPM': means['LOY']})
    dna.to_csv(DATA / 'dna_inset.csv', index=False)
    rna.to_csv(DATA / 'rna_inset.csv', index=False)
    mapping = pd.read_csv(SOURCES['gene_coordinates']).set_index('gene_name').loc[GENES]
    assert mapping.chromosome.eq('chrY').all()
    mapping.to_csv(DATA / 'gene_coordinates.csv')
    bands.to_csv(DATA / 'cytobands.csv', index=False)
    manifest = {
        'inputs': {key: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                   for key, path in SOURCES.items()},
        'dna': 'Median corrected, unsegmented relative copy number in non-PAR Yp and euchromatic Yq; one WGS specimen per condition; analytical floor 0.001.',
        'rna': '100 x LOY / CTRL condition mean CPM, four August 2025 RNA samples per condition; current shared-sex-chromosome Chapter 3 export. No inferred pairing or new DE tests.',
        'context': 'Chapter 3 A549 characterisation; not screen-passage or specimen-matched DNA/RNA measurements.',
        'selection': 'Same eight MSY genes and genomic order as Fig 4.1a.',
        'colour_scales': {'DNA_relative_CN': [0, .55], 'RNA_percent_CTRL': [0, 100]},
        'verification': {'DNA_medians_reproduced_from_bins': True, 'RNA_means_reproduced_from_CPM': True},
    }
    (DATA / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    return manifest


def draw_molecular_confirmation(fig, rect):
    ax = fig.add_axes(rect, xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ideogram = ax.inset_axes([.02, .92, .96, .045])
    bands = pd.read_csv(DATA / 'cytobands.csv')
    for band in bands.itertuples():
        start, end = band.start/1e6, band.end/1e6
        if band.stain == 'acen':
            points = ([(start, .1), (start, .9), (end, .5)] if band.band.startswith('p')
                      else [(start, .5), (end, .1), (end, .9)])
            ideogram.add_patch(Polygon(points, facecolor='#777777', edgecolor='#555555', lw=.4))
        else:
            ideogram.add_patch(Rectangle((start, .1), end-start, .8,
                facecolor={'gneg':'white', 'gpos50':'#999999', 'gvar':'#AFC0CD'}[band.stain],
                edgecolor='#555555', lw=.4))
    length = bands.end.max()/1e6; cen = bands.loc[bands.stain.eq('acen')]
    ideogram.set(xlim=(0,length), ylim=(0,1)); ideogram.axis('off')
    ideogram.text(-.025,.5,'Y',transform=ideogram.transAxes,ha='right',va='center',fontsize=8.5)
    ideogram.text(cen.start.min()/2e6,1.3,'Yp',ha='center',va='bottom',fontsize=7.5)
    ideogram.text((cen.end.max()/1e6+length)/2,1.3,'Yq',ha='center',va='bottom',fontsize=7.5)
    d = pd.read_csv(DATA/'dna_inset.csv').pivot(index='condition',columns='region',values='median_relative_cn')
    values = d.loc[['sgCTRL','LOY'],REGIONS].to_numpy()
    dna_ax = ax.inset_axes([.02,.57,.23,.20])
    di = dna_ax.imshow(values,cmap='viridis',vmin=0,vmax=.55,aspect='auto')
    dna_ax.set_xticks([0,1],['Yp','Yq'],fontsize=8)
    dna_ax.set_yticks([0,1],['CTRL','LOY'],fontsize=8)
    for label,colour in zip(dna_ax.get_yticklabels(),[BLUE,RED]): label.set_color(colour)
    dna_ax.tick_params(length=0,pad=3); dna_ax.set_title('DNA',fontsize=8.5,pad=4)
    for i in range(2):
        for j in range(2):
            value=values[i,j]
            label=f'{value:.3f}' if value < .01 else f'{value:.2f}'
            dna_ax.text(j,i,label,ha='center',va='center',fontsize=7.5,color='white' if value<.3 else DARK)
    rna = pd.read_csv(DATA/'rna_inset.csv').set_index('gene_name').loc[GENES]
    rv = rna[['sgCTRL','LOY']].to_numpy().T
    rna_ax = ax.inset_axes([.35,.57,.63,.20])
    ri = rna_ax.imshow(rv,cmap='magma',vmin=0,vmax=100,aspect='auto')
    for i,row in enumerate(rv):
        for j,value in enumerate(row):
            rna_ax.add_patch(Rectangle((j-.5,i-.5),1,1,fill=False,edgecolor='white',lw=.8))
            label=f'{value:.0f}' if i==0 or value==0 else (f'{value:.2f}' if value<.1 else f'{value:.1f}')
            rna_ax.text(j,i,label,ha='center',va='center',fontsize=7.5,color=DARK if value>=65 else 'white')
    rna_ax.set_xticks(range(len(GENES)),GENES,rotation=55,ha='right',fontsize=8,fontstyle='italic')
    rna_ax.set_yticks([]);rna_ax.tick_params(length=0,pad=2)
    rna_ax.set_title('RNA',fontsize=8.5,pad=4,bbox={'facecolor':'white','edgecolor':'none','pad':1})
    genes=pd.read_csv(DATA/'gene_coordinates.csv').set_index('gene_name').loc[GENES]
    for j,gene in enumerate(genes.itertuples()):
        ax.add_artist(ConnectionPatch(xyA=((gene.t2t_start+gene.t2t_end)/2e6,.1),coordsA=ideogram.transData,
            xyB=(j,-.5),coordsB=rna_ax.transData,color='#808080',alpha=.2,lw=.7,zorder=1,clip_on=False))
    for image, box, ticks, label in [(di,[.02,.19,.23,.025],[0,.5],'Relative CN'),
                                      (ri,[.52,.19,.46,.025],[0,50,100],'% of CTRL expression')]:
        cb=fig.colorbar(image,cax=ax.inset_axes(box),orientation='horizontal',ticks=ticks)
        cb.ax.tick_params(labelsize=7,length=2,pad=1);cb.set_label(label,fontsize=7.5,labelpad=1)
    for panel in [dna_ax,rna_ax]:
        for spine in panel.spines.values():spine.set_visible(False)
    ax.text(.5,-.12,'Chapter 3 characterisation',ha='center',fontsize=7.5,color='#727984')
    return ax


if __name__ == '__main__':
    print(json.dumps(prepare()['verification'],indent=2))
