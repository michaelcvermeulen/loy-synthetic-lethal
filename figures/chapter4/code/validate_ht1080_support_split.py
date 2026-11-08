"""Validate both figures from the 28 October Fig 4.3 split."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT.parent/'.figure_python_deps'))
import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image

FIG=ROOT/'figures/chapter4'
DATA=FIG/'source_data/ht1080_support_split_20261028'
DRAFT=ROOT/'draftv2/chapter4'

def main():
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    archive=Path((DATA/'archive_path.txt').read_text(encoding='utf-8'))
    original=json.loads((archive/'original_hashes.json').read_text(encoding='utf-8'))
    original_csvs=[p for p in original if p.endswith('.csv')]
    for p in original_csvs:assert sha(ROOT/p)==original[p],p
    manifest=json.loads((DATA/'render_manifest.json').read_text(encoding='utf-8'))
    profiles=pd.read_csv(FIG/'source_data/ht1080_candidate_discovery/guide_profiles.csv')
    checks={}
    for name,size,letters in [('ch04_main_ht1080_candidate_support',(9.6,8.6),'abc'),
                              ('ch04_supp_ht1080_method_and_d0_support',(9.6,7.0),'abcd')]:
        record=manifest[name]
        for p,digest in record['source_hashes'].items():assert sha(Path(p))==digest,p
        assert record['source_hashes_unchanged']
        assert not record['text_layout']['overlaps'] and not record['text_layout']['outside_page']
        doc=fitz.open(FIG/(name+'.pdf'));page=doc[0]
        assert len(doc)==1
        np.testing.assert_allclose([page.rect.width,page.rect.height],np.array(size)*72,atol=.01)
        fonts=page.get_fonts()
        assert fonts and all('OpenSans' in f[3] for f in fonts),fonts
        for font in fonts:assert len(doc.extract_font(font[0])[3])>0,font
        outside=[w[4] for w in page.get_text('words') if not (page.rect+(-.3,-.3,.3,.3)).contains(fitz.Rect(w[:4]))]
        assert not outside,outside
        text=page.get_text()
        for letter in letters:assert letter in text.splitlines(),(name,letter)
        assert 'Brunello sgRNA' in text
        assert ('Zero reads in all six samples' if letters=='abcd' else 'Zero reads in all six libraries') in text
        svg=(FIG/(name+'.svg')).read_text(encoding='utf-8')
        assert '<text' in svg and 'Open Sans' in svg
        with Image.open(FIG/(name+'.png')) as im:
            assert im.size==tuple(round(x*600) for x in size),im.size
            assert all(abs(x-600)<.1 for x in im.info['dpi'])
        checks[name]={'inches':size,'pdf_page_count':1,'embedded_fonts':[f[3] for f in fonts],
            'all_fonts_embedded':True,'editable_svg_text':True,'png_dpi':600,
            'text_overlaps':0,'text_outside_page':0,'sha256':{ext:sha(FIG/(name+'.'+ext)) for ext in ['pdf','svg','png']}}
        if letters=='abcd':
            for s in ['Spearman ρ = 0.877','Spearman ρ = 0.735','19,112 genes','19,082 genes','14.4','14.2','12.2']:assert s in text,s
            assert text.count('Brunello sgRNA (1–4 per gene)')==2
            for cell in ['14.4†','14.2†','12.2†']:assert text.count(cell)==2,(cell,text.count(cell))
            assert '† D0 PCR-jackpotted guide (masked in IsogenicZ)' in text
        doc.close()

    for name,timepoint,expected in [('main_panel_a_endpoint.csv','Endpoint',160),('main_panel_b_endpoint.csv','Endpoint',96),
                                     ('supp_panel_c_d0.csv','D0',88),('supp_panel_d_d0.csv','D0',48)]:
        actual=pd.read_csv(DATA/name)
        assert len(actual)==expected and actual.timepoint.eq(timepoint).all()
        source=profiles.loc[profiles.timepoint.eq(timepoint)&profiles.gene.isin(actual.gene.unique())]
        pd.testing.assert_frame_equal(actual[profiles.columns].reset_index(drop=True),source.reset_index(drop=True))
        if timepoint=='Endpoint':np.testing.assert_allclose(actual.plotted_cpm_plus_1,actual.cpm+1,rtol=1e-14,atol=1e-14)
        else:np.testing.assert_array_equal(actual.plotted_log2_cpm_plus_1,actual.log2_cpm_plus_1)
    source=FIG/'source_data/ht1080_candidate_discovery'
    method=FIG/'source_data/ht1080_method_panels_20261022'
    iso=FIG/'source_data/ht1080_isogenicz_20261020'
    comp=pd.read_csv(source/'method_comparison.csv').set_index('GENE')
    a=pd.read_csv(method/'panel_a_primary_drugz_mageck_ranks.csv').set_index('Gene')
    b=pd.read_csv(method/'panel_b_matched_drugz_isogenicz_ranks.csv').set_index('Gene')
    ig=pd.read_csv(iso/'mask173_detected_IsogenicZ_genes.tsv',sep='\t').set_index('Gene')
    dz=pd.read_csv(iso/'mask173_detected_matched_endpoint_drugZ_genes.tsv',sep='\t').set_index('Gene')
    for actual,expected in [(a.x,comp.rank_synth),(a.y,comp['neg|rank']),(b.x,dz.loc[b.index,'rank_synth']),(b.y,ig.loc[b.index,'rank_synth'])]:
        np.testing.assert_array_equal(actual,expected)
    for f in [a,b]:np.testing.assert_allclose(f.primary_drugZ_FDR,comp.loc[f.index,'fdr_synth'],atol=0,rtol=1e-14)
    assert (len(a),len(b))==(19112,19082)
    saga=pd.read_csv(source/'saga_reference_scores.csv')
    hits=saga.loc[saga.reference_20 & saga.fdr_synth.lt(.25)]
    assert saga.reference_20.sum()==20 and len(saga)==24 and len(hits)==6 and hits.fdr_synth.lt(.05).sum()==5
    assert set(hits.module)=={'Core','HAT','DUB'}
    source_schematic=pd.read_csv(source/'saga_schematic_members.csv')
    pd.testing.assert_frame_equal(pd.read_csv(DATA/'saga_schematic_members.csv'),source_schematic)
    outliers=pd.read_csv(iso/'mask173_raw_counts.tsv',sep='\t',usecols=[0]).iloc[:,0]
    assert not outliers.isin(['DDX3X_4638_GCACCACCATAAACCACGCA','KDM6A_20072_TCTTTGTATGAACAGCTGGG']).any()

    central=(ROOT/'Chapter4_figure_captions.md').read_text(encoding='utf-8')
    blocks=[b for b in re.split(r'(?=^## )',central,flags=re.M) if re.match(r'## (main|supp) - (?:Supp )?Fig 4\.\d+\.',b)]
    mains=[];supps=[]
    for block in blocks:
        match=re.match(r'## (main|supp) - (?:Supp )?Fig 4\.(\d+)\.',block)
        (mains if match[1]=='main' else supps).append(int(match[2]))
    assert sorted(mains)==list(range(1,9)),mains
    assert sorted(supps)==list(range(1,13)),supps
    anchors=re.findall(r'<a id="([^"]+)"',central)
    assert len(anchors)==len(set(anchors))
    results=(DRAFT/'chapter4_results.md').read_text(encoding='utf-8')
    integration=json.loads((DATA/'integration.json').read_text(encoding='utf-8'))
    for key in ['main','supplement']:
        cap=integration[key]['caption'];assert cap in results and cap in central
        p=Path(integration['previews'][integration[key]['stem']]);assert p.exists() and p.as_posix() in results
    for p in DRAFT.glob('*.md'):
        if 'source_notes' in p.name or 'numbering' in p.name:continue
        text=p.read_text(encoding='utf-8')
        for anchor in re.findall(r'Chapter4_figure_captions\.md#([A-Za-z0-9_-]+)',text):assert anchor in anchors,(p.name,anchor)
        # Only main-figure calls: supplementary panel d remains valid.
        bad=re.findall(r'(?<!Supp )(?:\*\*Fig\*\*|Fig)\s+4\.3[d-g](?:\b|,)',text)
        assert not bad,(p.name,bad)
        assert '\ufffd' not in text
    active_numbers=[int(n) for n in re.findall(r'^\*\*Supp Fig 4\.(\d+)\.',results,flags=re.M)]
    assert active_numbers==[1]+list(range(3,13)),active_numbers
    record={'source_csvs_byte_identical':len(original_csvs),'all_input_hashes_preserved':True,
        'endpoint_points':{'a':160,'b':96},'D0_cells':{'c':88,'d':48},'ranks_scores_fdrs_preserved':True,
        'saga_reference_subunits':20,'saga_alternatives':4,'saga_hits_fdr025':6,'saga_hits_fdr005':5,
        'D0_outlier_guides_confirmed':True,'active_main_figures':8,'active_supplementary_figures':12,
        'central_and_inline_captions_match':True,'caption_anchors_resolve':True,'old_main_panel_calls_removed':True,
        'visual_review':'Both complete previews inspected at 160 dpi; colour scales, guide order, SAGA labels and new panel letters are legible and unclipped.',
        'exports':checks}
    (DATA/'validation.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if k!='exports'},indent=2))

if __name__=='__main__':main()
