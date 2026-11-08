"""Workflow, flow QC, known-ratio recovery, genomic PCR and clone ddPCR."""
from pathlib import Path
import hashlib
import io
import json
import re
import sys
import xml.etree.ElementTree as ET

THESIS = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(THESIS.parent/'.figure_python_deps'))
import pymupdf
from PIL import Image
from draw_competition_qc_panels import render_lower, HEIGHT as LOWER_HEIGHT

OUT = THESIS/'figures/chapter4'
DATA = OUT/'source_data/competition_workflow_qc'
REVIEW = OUT/'review'
STEM = 'ch04_main_competition_workflow_qc_draft'
WORKFLOW = OUT/'review/unstained_overlay_draft/ch04_main_competition_assay_workflow_unstained_overlay_draft'
RATIOS = OUT/'ch04_main_competition_known_ratios_draft'
SVG = 'http://www.w3.org/2000/svg'
ET.register_namespace('', SVG)
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def svg_group(path,prefix,y,rename=False):
    root=ET.parse(path).getroot()
    width,height=map(float,root.attrib['viewBox'].split()[2:])
    ids={e.attrib['id']:prefix+e.attrib['id'] for e in root.iter() if 'id' in e.attrib}
    letters=0
    for e in root.iter():
        for key,value in list(e.attrib.items()):
            if key=='id':
                e.set(key,ids[value])
            else:
                value=re.sub(r'url\(#([^)]*)\)',lambda m:'url(#'+ids.get(m[1],m[1])+')',value)
                if value.startswith('#') and value[1:] in ids:
                    value='#'+ids[value[1:]]
                e.set(key,value)
        if rename and e.tag==f'{{{SVG}}}text' and e.text in ('a','b'):
            e.text={'a':'g','b':'h'}[e.text]
            e.set('style',e.attrib['style'].replace('12px','9px'))
            letters+=1
    if rename:
        assert letters==2
    group=ET.Element(f'{{{SVG}}}g',{'transform':f'translate(0 {y})'})
    for e in root:
        group.append(e)
    return group,width,height


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    sources=[p.with_suffix(ext) for p in (WORKFLOW,RATIOS) for ext in ('.pdf','.svg')]
    sources += [OUT/'source_data/competition_known_ratios'/name for name in
                ('measurements.csv','summary.csv','source_manifest.json')]
    sources += [OUT/'source_data/competition_workflow'/name for name in
                ('flow_examples_manifest.json','flow_controls_manifest.json')]
    sources += list((OUT/'source_data/competition_workflow').glob('*events.csv'))
    sources += [OUT/'source_data/ht1080_ko_validation_20261022'/name for name in
                ('grant_figure5_original.tiff','display_values.csv','DDX3X_clones_grant_source.csv',
                 'LOY_clones.xlsx','workbook_verification.csv','provenance_and_validation.json')]
    sources += [OUT/'source_data/ht1080_author_ko_gel_20261022'/name for name in
                ('author_supplied_ko_gel.png','provenance.json')]
    source_hashes={str(p):sha(p) for p in sources}
    lower_stem=render_lower()
    top=pymupdf.open(WORKFLOW.with_suffix('.pdf'))
    bottom=pymupdf.open(lower_stem.with_suffix('.pdf'))
    assert len(top)==len(bottom)==1
    assert abs(top[0].rect.width-bottom[0].rect.width)<.001
    width=top[0].rect.width
    top_h=80.5*72/25.4
    bottom_h=bottom[0].rect.height
    # Retain panel a at original size; flow panels are redrawn in the lower grid.
    top[0].add_redact_annot(pymupdf.Rect(0,top_h,width,top[0].rect.height),fill=(1,1,1))
    top[0].apply_redactions(images=0,graphics=0,text=0)
    top[0].set_cropbox(pymupdf.Rect(0,0,width,top_h))
    combined=pymupdf.open()
    page=combined.new_page(width=width,height=top_h+bottom_h)
    page.show_pdf_page(pymupdf.Rect(0,0,width,top_h),top,0)
    page.show_pdf_page(pymupdf.Rect(0,top_h,width,top_h+bottom_h),bottom,0)
    combined.set_metadata({'title':'Competition-assay workflow, measurement QC and population characterisation',
                           'creator':'Vector composition of approved thesis panels'})
    combined.save(OUT/f'{STEM}.pdf',garbage=4,deflate=True)
    for dpi,path in [(600,OUT/f'{STEM}.png'),(190,REVIEW/f'{STEM}_preview.png')]:
        pix=page.get_pixmap(dpi=dpi,alpha=False)
        Image.frombytes('RGB',(pix.width,pix.height),pix.samples).save(path,dpi=(dpi,dpi))
    page.get_pixmap(dpi=140,colorspace=pymupdf.csGRAY).save(REVIEW/f'{STEM}_greyscale_check.png')
    group,w,h=svg_group(WORKFLOW.with_suffix('.svg'),'workflow_',0)
    h=top_h
    # Nested SVG viewport clips the source to panel a without scaling it.
    viewport=ET.Element(f'{{{SVG}}}svg',{'width':str(w),'height':str(h),
                        'viewBox':f'0 0 {w} {h}','overflow':'hidden'})
    viewport.append(group)
    lower,w2,h2=svg_group(lower_stem.with_suffix('.svg'),'compact_',h)
    assert abs(w-w2)<.001
    svg=ET.Element(f'{{{SVG}}}svg',{'width':f'{w}pt','height':f'{h+h2}pt',
                                'viewBox':f'0 0 {w} {h+h2}','version':'1.1'})
    desc=ET.SubElement(svg,f'{{{SVG}}}desc')
    desc.text='Panel a: approved workflow and schematic follow-up. Panels b-f: frozen control and baseline flow examples. Panel g: September known-ratio recovery. Panel h: original genomic PCR gel with retypeset labels. Panel i: 36 ddPCR measurements for 12 DDX3X-targeted clones.'
    svg.extend([viewport,lower])
    ET.ElementTree(svg).write(OUT/f'{STEM}.svg',encoding='utf-8',xml_declaration=True)
    text=page.get_text()
    assert 'Illustrative hit outcomes' not in text
    assert 'Pre-existing GFP/mCherry matrix' in text
    assert all(label in text.splitlines() for label in 'abcdefghi')
    words=page.get_text('words')
    outside=[w for w in words if w[0]<-.5 or w[1]<-.5 or w[2]>page.rect.width+.5 or w[3]>page.rect.height+.5]
    assert not outside,outside
    fonts=page.get_fonts(full=True)
    assert fonts and all('OpenSans' in f[3].replace(' ','') for f in fonts),fonts
    # Check that composition preserves visible workflow and redrawn-panel text.
    expected=top[0].get_text('words')+bottom[0].get_text('words')
    assert sorted(w[4] for w in words)==sorted(w[4] for w in expected)
    for path,digest in source_hashes.items():
        assert sha(Path(path))==digest
    ids=[e.attrib['id'] for e in svg.iter() if 'id' in e.attrib]
    assert len(ids)==len(set(ids))
    png=Image.open(OUT/f'{STEM}.png')
    assert min(png.info['dpi'])>590
    gel=Image.open(DATA/'panel_h_gel_native.png').convert('RGB')
    intact=False
    for entry in page.get_images(full=True):
        stored=Image.open(io.BytesIO(combined.extract_image(entry[0])['image'])).convert('RGB')
        if stored.size==gel.size:
            intact=(stored.tobytes()==gel.tobytes() or
                    stored.transpose(Image.Transpose.FLIP_TOP_BOTTOM).tobytes()==gel.tobytes())
    assert intact,'Gel pixels must be retained in the final PDF'
    result={'figure_mm':[183,80.5+LOWER_HEIGHT],'panel_map':{'a':'workflow and schematic follow-up',
        'b':'historical unstained, visual gate overlay','c':'GFP-only','d':'mCherry-only',
        'e':'mixture 1 baseline','f':'mixture 2 baseline','g':'known-ratio recovery',
        'h':'genomic PCR','i':'Y status of 12 DDX3X-targeted clones'},
        'source_hashes':source_hashes,'source_files_unchanged':True,'workflow_panel_a_unscaled':True,
        'layout':'Five flow plots across one row; g, h and i share the bottom row. Panels a-f retain their physical size and position.',
        'gel_pixels_preserved_in_pdf':intact,'clone_measurements':36,
        'flow_data_limits_unchanged':True,'known_ratio_means_and_ranges_unchanged':True,
        'workflow_variant':'Latest approved unstained-overlay draft with shortened flow wording; dates omitted from b-e and illustrative-heading omitted.',
        'historical_overlay':'Visual only; unstained events remain grey and unclassified.',
        'pdf_pages':1,'pdf_text_outside_page':0,'pdf_composition_text_preserved':True,
        'svg_text_elements':sum(e.tag==f'{{{SVG}}}text' for e in svg.iter()),
        'svg_ids_unique':True,'png_pixels':list(png.size),'png_dpi':list(png.info['dpi']),
        'renderer_sha256':sha(Path(__file__)),
        'lower_panel_renderer_sha256':sha(Path(__file__).with_name('draw_competition_qc_panels.py'))}
    (DATA/'composition_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))
    combined.close();top.close();bottom.close()


if __name__=='__main__':
    main()
