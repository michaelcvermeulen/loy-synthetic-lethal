"""HT1080 F03 main/supplement split; preserves original measurements."""
from pathlib import Path
import sys, io, json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ht1080_dna_rna_summary import summarize

FIG = Path(__file__).resolve().parents[1]
THESIS = FIG.parents[1]
SOURCE = FIG/'source_data/ht1080_lpwgs_confirmation/qdnaseq_bins.csv'
CYTO = THESIS.parent/'WGS/chm13v2.0_cytobands_allchrs.bed'
OUT = FIG/'source_data/ht1080_lpwgs_difference_draft'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(THESIS/'skills/thesis-figure-style/scripts'))
from thesis_style import apply_style
apply_style(8)
STEM = 'ch03_main_ht1080_lpwgs_genomewide'
chroms = ['chr'+str(i) for i in range(1, 23)]+['chrX', 'chrY']
cb = pd.read_csv(CYTO, sep='\t', header=None, names=['chrom','start','end','band','stain'])
lengths = cb.groupby('chrom').end.max()
yq = cb[cb.chrom.eq('chrY') & cb.band.eq('q12')].iloc[0]
d = pd.read_csv(SOURCE)
d = d[d.use & np.isfinite(d.log2)]
paired = d.pivot(index=['chromosome','start','end'], columns='label', values='log2').dropna().reset_index()
assert len(paired) == 31185
excluded = paired.chromosome.eq('chrY') & paired.start.le(yq.end) & paired.end.gt(yq.start)
paired.loc[excluded].to_csv(OUT/'Yq12_excluded_source_bins.csv', index=False)
paired = paired.loc[~excluded].copy()
paired['mid'] = (paired.start+paired.end)/2
paired['difference'] = paired.LOY-paired.sgCTRL

def aggregate(data, size):
    data = data.copy()
    data['window_start'] = (data.mid//size).astype(int)*size
    b = data.groupby(['chromosome','window_start']).agg(
        sgCTRL=('sgCTRL','median'), LOY=('LOY','median'), n_bins=('start','size'),
        source_start=('start','min'), source_end=('end','max')).reset_index()
    b['window_end'] = np.minimum(b.window_start+size, b.source_end)
    b['mid'] = (b.window_start+b.window_end)/2
    b['difference'] = b.LOY-b.sgCTRL
    return b

bins = aggregate(paired, 5_000_000)
bins.to_csv(OUT/'genome_5Mb.csv', index=False)
arms = []
for ch in chroms[:-1]:
    cen = cb[cb.chrom.eq(ch) & cb.stain.eq('acen')]
    assert len(cen) == 2
    for arm, left, right in [('p',0,cen.start.min()),('q',cen.end.max(),lengths[ch])]:
        v = paired[paired.chromosome.eq(ch) & paired.mid.gt(left) & paired.mid.le(right)]
        arms.append(dict(chromosome=ch, arm=ch[3:]+arm, n_bins=len(v),
                         difference=v.difference.median(),
                         low_mappability_arm=ch in ['chr13','chr14','chr15','chr21','chr22'] and arm=='p'))
arms = pd.DataFrame(arms)
arms.to_csv(OUT/'arm_medians.csv', index=False)
windows_rna,arms_rna,rna_provenance=summarize(THESIS,OUT,paired,arms,cb,aggregate)
# Mask centromere-overlapping bins in the two non-Y details and the new Y inset.
# BED boundaries are zero-based half-open; QDNAseq bins are one-based inclusive.
detail_masks={}
def centromere_masked_detail(chrom):
    cen=cb[cb.chrom.eq(chrom) & cb.stain.eq('acen')]
    left,right=int(cen.start.min()),int(cen.end.max())
    source=paired[paired.chromosome.eq(chrom)].copy()
    excluded=source.start.le(right) & source.end.gt(left)
    source.loc[excluded].to_csv(OUT/f'{chrom}_centromere_excluded_source_bins.csv',index=False)
    retained=source.loc[~excluded].copy()
    parts=[]
    for arm,v in [('p',retained[retained.end.le(left)]),('q',retained[retained.start.gt(right)])]:
        b=aggregate(v,1_000_000)
        b['arm']=arm
        # Position boundary-window points in their retained source span.
        b['mid']=(b.source_start+b.source_end)/2
        parts.append(b)
    result=pd.concat(parts,ignore_index=True)
    assert not ((result.source_start.le(right)) & (result.source_end.gt(left))).any()
    assert result.n_bins.sum()==len(retained)
    detail_masks[chrom]=dict(BED_start=left,BED_end=right,excluded_source_bins=int(excluded.sum()),retained_source_bins=len(retained),plotted_windows=len(result))
    return result
detail = centromere_masked_detail('chr13')
detail.to_csv(OUT/'chr13_1Mb.csv', index=False)
chr4_detail=centromere_masked_detail('chr4')
chr4_detail.to_csv(OUT/'chr4_1Mb.csv',index=False)
segment_delta = 0.479599094060819-0.111904199476798
offsets = {}; end = 0
for ch in chroms:
    offsets[ch] = end
    end += lengths[ch]/1e6+20

# Reorganise the display without changing the retained genomic/RNA summaries.
from ht1080_lpwgs_revision_panels import render_main, render_supp, SUPP_STEM
ydetail = centromere_masked_detail('chrY')
ydetail.to_csv(OUT/'chrY_inset_1Mb.csv', index=False)
scatter_labels = render_main(FIG, STEM, bins, windows_rna, ydetail, cb, chroms, lengths, offsets, end, yq, arms_rna)
if '--main-only' not in sys.argv:
    render_supp(FIG, arms, detail, chr4_detail, cb, lengths, segment_delta)
pd.DataFrame(scatter_labels).to_csv(OUT/'RNA_scatter_arm_labels.csv', index=False)
assert len(bins)==630 and excluded.sum()==351
assert len(arms)==46 and np.isfinite(arms.difference).all()
assert arms.difference.between(-.41,.43).all()
assert np.isclose(2**segment_delta,1.29029,rtol=1e-4)
report=dict(input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE,CYTO]},
  comparison='HT1080-LOY-D0 (CP + Q53 2x) minus HT1080-sgCTRL-D0; one WGS specimen per condition',
  source_bins_before_mask=31185,excluded_Yq12_bins=int(excluded.sum()),genome_windows=len(bins),
  Yq12_BED=[int(yq.start),int(yq.end)],coordinate_rule='exclude every overlapping one-based inclusive source bin',
  window_estimator='median corrected unsegmented log2 per sample, then subtract sgCTRL from LOY',
  arm_estimator='median paired 100-kb log2 difference; midpoints outside centromere; all 46 autosomal/X arms',
  arm_caution='Acrocentric short arms are retained, but their repetitive sequence and floor-level values limit interpretation.',
  genome_axis='linear axis spacing of log2 differences; all 630 window values retained',
  main_RNA_axes='linear axis spacing of log2 differences; all 1597 windows retained',
  Y_inset='1-Mb per-sample log2 medians back-transformed to relative CN; centromere/Yq12 excluded',
  supplementary_stem=SUPP_STEM,
  main_arm_scatter_axes='DNA x: symlog, linear within +/-0.3; RNA y: linear and shared with main b',
  shared_RNA_axis={'panels':['b','c'], 'scale':'linear', 'limits':[-8.3,2.5], 'ticks':[-8,-6,-4,-2,0,2], 'equal_physical_height':True},
  panel_mapping={'old a':'main a', 'old d':'main b', 'old b':'supp a', 'old e':'main c', 'old c':'supp b', 'old f':'supp c'},
  export_dpi=600,
  chr13_original_segment_delta=segment_delta,chr13_relative_ratio=2**segment_delta,
  font='Open Sans',source_measurements_unchanged=True,selected_exports_refreshed_separately=True)
report['RNA_comparison']=rna_provenance
report['detail_centromere_masks']=detail_masks
(OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
