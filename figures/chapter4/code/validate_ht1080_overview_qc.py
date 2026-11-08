"""Check export readability metadata, source integrity, and PDF text bounds."""
from pathlib import Path
import hashlib
import json
import sys
import xml.etree.ElementTree as ET

THESIS=Path(__file__).resolve().parents[3]
DEPS=THESIS.parent/'.figure_python_deps'
if DEPS.is_dir(): sys.path.insert(0,str(DEPS))
import pymupdf as fitz
from PIL import Image
import numpy as np
import pandas as pd

OUT=THESIS/'figures/chapter4'
DATA=OUT/'source_data/ht1080_figure1_qc'
STEMS=['ch04_main_ht1080_crispr_screen_workflow','ch04_supp_ht1080_crispr_screen_qc']


def main():
    manifest=json.loads((DATA/'source_manifest.json').read_text())
    for key,entry in manifest.items():
        actual=hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()
        assert actual==entry['sha256'],f'Source changed: {key}'
    result={'original_source_hashes_unchanged':True,'exports':{}}
    pca=json.loads((DATA/'pca_manifest.json').read_text())
    for entry in pca['sources'].values():
        assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()==entry['sha256']
    scores=pd.read_csv(DATA/'pca_sample_scores.csv')
    assert len(scores)==6 and scores['sample'].is_unique
    assert pca['guides_retained']==77205 and pca['zero_variance_guides_excluded']==236
    assert np.isclose(sum(pca['variance_percent']),100)
    assert np.allclose(scores[[f'PC{i}' for i in range(1,6)]].sum(),0,atol=1e-8)
    result['pca']={'retained_guides':pca['guides_retained'],'variance_percent':pca['variance_percent'][:2],
                   'source_hashes_unchanged':True,'independent_gram_eigenvalue_check':pca['independent_gram_eigenvalue_check']}
    density=json.loads((DATA/'endpoint_reference_density_manifest.json').read_text())
    for entry in density['sources'].values():
        assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()==entry['sha256']
    guide_table=pd.read_csv(DATA/'endpoint_reference_density_guides.csv')
    curve_table=pd.read_csv(DATA/'endpoint_reference_density_curves.csv')
    original_guides=pd.read_csv(density['sources']['mageck_guides']['path'],sep='\t').set_index('sgrna')
    raw_counts=pd.read_csv(density['sources']['counts']['path'],sep='\t').set_index('sgRNA')
    endpoint_columns=['POST_CTRL_REP1','POST_CTRL_REP2','POST_LOY_REP1','POST_LOY_REP2']
    assert guide_table.sgrna.is_unique and len(guide_table)==8820
    assert np.allclose(guide_table.LFC,original_guides.loc[guide_table.sgrna,'LFC'],rtol=0,atol=1e-12)
    totals=raw_counts.loc[guide_table.sgrna,endpoint_columns].sum(axis=1).to_numpy()
    assert np.array_equal(guide_table.endpoint_raw_reads,totals)
    assert np.array_equal(guide_table.included,totals>0)
    expected_counts={'Non-targeting':999,'Essential':4773,'Non-essential':2983}
    for category,n in expected_counts.items():
        subset=guide_table.loc[guide_table.category.eq(category) & guide_table.included]
        assert len(subset)==n==density['groups'][category]['included_guides']
        assert curve_table.LFC.min()<subset.LFC.min() and curve_table.LFC.max()>subset.LFC.max()
        assert np.all(curve_table[category]>=0)
        assert np.isclose(np.trapezoid(curve_table[category],curve_table.LFC),1,atol=1e-4)
    result['endpoint_reference_densities']={'original_guide_lfc_preserved':True,
        'eligibility_matches_raw_endpoint_reads':True,'included_guides':expected_counts,
        'common_bandwidth_log2_units':density['bandwidth_log2_units'],
        'full_observed_range_included':True,'density_integrals_equal_one':True}
    ddpcr=json.loads((DATA/'ddpcr_integration_manifest.json').read_text())
    for entry in ddpcr['sources'].values():
        assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()==entry['sha256']
    wells=pd.read_csv(ddpcr['sources']['plotted_wells.csv']['path'])
    summary=pd.read_csv(ddpcr['sources']['plotted_summary.csv']['path'])
    value='estimated_target_copies_per_cell'
    keys=['condition','passage','arm']
    assert len(wells)==29 and len(summary)==15
    observed=wells.groupby(keys)[value].mean().sort_index()
    saved=summary.set_index(keys)[value].sort_index()
    assert observed.index.equals(saved.index)
    assert np.allclose(observed,saved,rtol=0,atol=1e-14)
    assert np.allclose(wells[value],wells.concentration_ratio*wells.reference_copy_number,rtol=0,atol=1e-14)
    result['ddpcr_integration']={'panel':'Supplementary Figure 1h','source_hashes_unchanged':True,
        'wells':len(wells),'summary_points':len(summary),'summary_matches_well_means':True,
        'reference_copy_scaling_preserved':True}
    for stem in STEMS:
        pdf=fitz.open(OUT/f'{stem}.pdf')
        assert len(pdf)==1
        page=pdf[0]
        words=page.get_text('words')
        assert len(words)>80
        outside=[w[:5] for w in words if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
        assert not outside,(stem,outside)
        fonts=page.get_fonts()
        assert fonts and all('OpenSans' in font[3] for font in fonts),fonts
        svg=ET.parse(OUT/f'{stem}.svg')
        text_elements=len(svg.findall('.//{http://www.w3.org/2000/svg}text'))
        assert text_elements>70
        png=Image.open(OUT/f'{stem}.png')
        assert min(png.info.get('dpi',(0,0)))>590
        assert png.width>5000
        pix=page.get_pixmap(matrix=fitz.Matrix(1.4,1.4),colorspace=fitz.csGRAY)
        pix.save(OUT/'review'/f'{stem}_greyscale_check.png')
        result['exports'][stem]={'pdf_pages':len(pdf),'pdf_text_words':len(words),
            'out_of_bounds_text':len(outside),'fonts':[f[3] for f in fonts],
            'svg_editable_text_elements':text_elements,'png_pixels':[png.width,png.height],
            'png_dpi':list(png.info['dpi'])}
        pdf.close()
    script_paths=list(Path(__file__).parent.glob('*ht1080_overview_qc.py'))
    script_paths += [Path(__file__).parent/name for name in ['make_ht1080_screen_schematic.py','prepare_ht1080_reference_densities.py','make_ht1080_ddpcr_passage.py']]
    result['scripts']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in script_paths}
    (DATA/'export_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
