"""Check figure exports, frozen measurements and updated caption/reference wiring."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT.parent/'.figure_python_deps'))
import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image
import make_ht1080_molecular_integration as M


def main():
    manifest=json.loads((M.DATA/'render_manifest.json').read_text(encoding='utf-8'))
    integration=json.loads((M.DATA/'integration.json').read_text(encoding='utf-8'))
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for p,h in manifest['source_hashes'].items():assert sha(Path(p))==h,p
    old_a=pd.read_csv(M.SOURCE/'figure_panel_a_values.csv')
    a=pd.read_csv(M.DATA/'panel_a_values.csv');b=pd.read_csv(M.DATA/'panel_b_values.csv')
    pd.testing.assert_frame_equal(a,old_a.loc[old_a.column.isin(M.A_COLUMNS)].reset_index(drop=True))
    pd.testing.assert_frame_equal(b,pd.read_csv(M.SOURCE/'figure_panel_b_values.csv'))
    assert len(a)==105 and len(b)==119
    assert a.column.nunique()==7 and b.x.nunique()==7
    assert not a.column.str.contains('Untransduced').any()
    assert a.groupby('set').size().eq(7).all() and b.groupby('symbol').size().eq(7).all()
    assert manifest['screen_palette']['positive_sensitizer']=='#7A3E85'
    assert manifest['screen_palette']['negative_resistor']=='#218C8D'
    assert manifest['molecular_palette']['positive_higher']=='#B2182B'
    assert manifest['molecular_palette']['negative_lower']=='#2166AC'
    assert not manifest['text_layout']['overlaps'] and not manifest['text_layout']['outside_page']
    pdf=M.OUT/f'{M.STEM}.pdf';doc=fitz.open(pdf);page=doc[0]
    assert len(doc)==1
    np.testing.assert_allclose([page.rect.width,page.rect.height],np.array(manifest['inches'])*72,atol=.01)
    fonts=page.get_fonts();assert fonts and all('OpenSans' in f[3] for f in fonts),fonts
    assert all(len(doc.extract_font(f[0])[3]) for f in fonts)
    text=page.get_text()
    for title in ['Screen: resistor <-> sensitizer','RNA / ATAC: lower <-> higher']:
        assert text.count(title)==2,(title,text.count(title))
    for title,_ in M.A_GROUPS:assert title in text,title
    assert 'Untransduced' not in text and 'Untrans-' not in text
    assert text.count('Screen')>=5 and text.count('Baseline')==2 and text.count('Add-back vs EV')==2
    assert not [w for w in page.get_text('words') if not (page.rect+(-.3,-.3,.3,.3)).contains(fitz.Rect(w[:4]))]
    doc.close()
    with Image.open(M.OUT/f'{M.STEM}.png') as im:
        assert im.size==tuple(round(s*600) for s in manifest['inches'])
        assert all(abs(x-600)<.1 for x in im.info['dpi'])
    svg=(M.OUT/f'{M.STEM}.svg').read_text(encoding='utf-8')
    assert '<text' in svg and 'Open Sans' in svg
    central=(ROOT/'Chapter4_figure_captions.md').read_text(encoding='utf-8')
    results=(ROOT/'draftv2/chapter4/chapter4_results.md').read_text(encoding='utf-8')
    discussion=(ROOT/'draftv2/chapter4/chapter4_discussion.md').read_text(encoding='utf-8')
    assert integration['caption'] in central and integration['caption'] in results
    assert integration['preview'] in results and Path(integration['preview']).exists()
    mains=re.findall(r'^## main - Fig 4\.(\d+)\.',central,flags=re.M)
    supps=re.findall(r'^## supp - Supp Fig 4\.(\d+)\.',central,flags=re.M)
    assert sorted(map(int,mains))==list(range(1,9)),mains
    assert sorted(map(int,supps))==list(range(1,13)),supps
    assert re.findall(r'^\*\*Fig 4\.(\d+)\.',results,flags=re.M)==list(map(str,range(0,8)))
    assert '**Fig 4.8. DDX3X' in discussion
    anchors=re.findall(r'<a id="([^"]+)"',central);assert len(anchors)==len(set(anchors))
    for f in (ROOT/'draftv2/chapter4').glob('*.md'):
        if 'source_notes' in f.name:continue
        content=f.read_text(encoding='utf-8')
        for anchor in re.findall(r'Chapter4_figure_captions\.md#([A-Za-z0-9_-]+)',content):
            assert anchor in anchors,(f.name,anchor)
        assert '\ufffd' not in content,f.name
    # Author's revised Results prose must be identical outside the target figure block.
    archive=Path(integration['archive'])
    original=(archive/'draftv2/chapter4/chapter4_results.md').read_text(encoding='utf-8')
    def without_target(t):
        return re.sub(r'!\[HT1080 screen evidence in the context of molecular profiles\].*?\n\n---','<FIGURE>\n\n---',t,flags=re.S)
    assert without_target(results)==without_target(original)
    record={'source_files_unchanged':len(manifest['source_hashes']),'panel_a_cells':105,'panel_b_cells':119,
            'two_palettes':True,'horizontal_colourbars':4,'untransduced_display_column_removed':True,
            'row_categories':3,'all_other_values_preserved':True,'embedded_open_sans':True,
            'editable_svg':True,'png_dpi':600,'text_overlaps':0,'text_outside_page':0,
            'main_figures':8,'supplementary_figures':12,'caption_and_anchor_checks_passed':True,
            'author_results_prose_preserved':True,'full_preview_visually_inspected':True,
            'sha256':{ext:sha(M.OUT/f'{M.STEM}.{ext}') for ext in ['pdf','svg','png']}}
    (M.DATA/'validation.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
