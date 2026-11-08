"""Verify the frozen guide values, eligibility and publication exports."""
from pathlib import Path
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent / '.figure_python_deps'))
import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image

OUT = THESIS / 'figures/chapter4'
DATA = OUT / 'source_data/ht1080_guide_rugs'
STEM = 'ch04_supp_ht1080_guide_consistency_draft'


def main():
    manifest = json.loads((DATA/'source_manifest.json').read_text())
    for entry in manifest['inputs'].values():
        assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest() == entry['sha256']
    original = pd.read_csv(manifest['inputs']['mageck_guides']['path'], sep='\t').set_index('sgrna')
    counts = pd.read_csv(manifest['inputs']['raw_counts']['path'], sep='\t').set_index('sgRNA')
    drugz = pd.read_csv(manifest['inputs']['drugz_genes']['path'], sep='\t').set_index('GENE')
    selected = pd.read_csv(DATA/'selected_guides.csv', float_precision='round_trip')
    ntc = pd.read_csv(DATA/'non_targeting_guides.csv', float_precision='round_trip')
    ticks = pd.read_csv(DATA/'displayed_tick_positions.csv', float_precision='round_trip')
    assert len(selected) == 64 and selected.Gene.nunique() == 16
    assert selected.groupby('Gene').size().eq(4).all()
    assert len(ticks) == 62 and ticks.sgRNA.is_unique
    assert len(ntc) == 1000 and ntc.lfc_displayed.sum() == 999
    assert set(ticks.sgRNA) == set(selected.loc[selected.lfc_displayed,'sgRNA'])
    assert selected.loc[~selected.lfc_displayed,'Gene'].tolist() == ['DDX3X','MED12']
    assert selected.loc[~selected.lfc_displayed,'all_samples_zero'].all()
    samples=['D0_CTRL','D0_LOY','POST_CTRL_REP1','POST_CTRL_REP2','POST_LOY_REP1','POST_LOY_REP2']
    for row in selected.itertuples():
        assert row.mageck_reported_log2FC == original.loc[row.sgRNA,'LFC']
        assert row.mageck_guide_FDR == original.loc[row.sgRNA,'FDR']
        assert row.drugz_normZ == drugz.loc[row.Gene,'normZ']
        field='fdr_supp' if row.drugz_fdr_direction == 'enrichment' else 'fdr_synth'
        assert row.drugz_directional_fdr == drugz.loc[row.Gene,field]
        assert all(getattr(row,s) == counts.loc[row.sgRNA,s] for s in samples)
        assert row.lfc_displayed == (sum(getattr(row,s) for s in samples[2:]) > 0)
    for row in ticks.itertuples():
        assert row.x == original.loc[row.sgRNA,'LFC']
        assert manifest['x_limits'][0] < row.x < manifest['x_limits'][1]
    assert ntc.Gene.eq('Non-Targeting-Control').all()
    assert np.array_equal(ntc.LFC.to_numpy(), original.loc[ntc.sgRNA,'LFC'].to_numpy())
    assert np.array_equal(ntc.lfc_displayed.to_numpy(), counts.loc[ntc.sgRNA,samples[2:]].sum(axis=1).gt(0).to_numpy())
    vals=ntc.loc[ntc.lfc_displayed,'LFC'].to_numpy()
    assert vals.min() > manifest['x_limits'][0] and vals.max() < manifest['x_limits'][1]
    density=pd.read_csv(DATA/'ntc_density.csv')
    assert np.isclose(np.trapezoid(density.kde_density, density.guide_log2FC),1,atol=1e-8)
    assert density.display_density_peak_scaled.max() == 1

    doc=fitz.open(OUT/f'{STEM}.pdf')
    assert len(doc) == 1
    page=doc[0];text=page.get_text()
    assert all(g in text for g in selected.Gene.unique())
    assert all(s in text for s in ['drugZ FDR','MAGeCK','999','staggered','no reads'])
    spans=[s for block in page.get_text('dict')['blocks'] if 'lines' in block
           for line in block['lines'] for s in line['spans'] if s['text'].strip()]
    bounds=[];overlaps=[]
    for i,s in enumerate(spans):
        rect=fitz.Rect(s['bbox'])
        if not page.rect.contains(rect): bounds.append(s['text'])
        for t in spans[i+1:]:
            intersection=rect & fitz.Rect(t['bbox'])
            if intersection.width > .5 and intersection.height > .5:
                overlaps.append([s['text'],t['text']])
    assert not bounds,bounds
    assert not overlaps,overlaps
    assert all('OpenSans' in s['font'] for s in spans),set(s['font'] for s in spans)
    svg=ET.parse(OUT/f'{STEM}.svg')
    editable=len(svg.findall('.//{http://www.w3.org/2000/svg}text'))
    assert editable >= 50
    png=Image.open(OUT/f'{STEM}.png')
    assert png.size == (5760,3690)
    assert np.isclose(png.info['dpi'][0],600,atol=.1)
    review=Image.open(OUT/'review'/f'{STEM}_preview.png')
    review.convert('L').save(OUT/'review'/f'{STEM}_greyscale_check.png')
    result={'original_sources_unchanged':True,'exact_saved_guide_and_gene_statistics':True,
            'selected_guides':64,'displayed_guides':62,'background_guides':999,
            'all_observed_values_within_axes':True,'kde_integrates_to_one':True,
            'pdf_pages':1,'text_outside_page':bounds,'text_overlaps':overlaps,
            'fonts':sorted(set(s['font'] for s in spans)), 'svg_text_nodes':editable,
            'png_pixels':list(png.size),'png_dpi':list(png.info['dpi']),
            'validator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (DATA/'export_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
