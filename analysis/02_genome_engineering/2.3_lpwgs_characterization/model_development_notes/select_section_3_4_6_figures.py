"""Freeze reviewed section 3.4.6 figures and assign citations in first-mention order."""
from pathlib import Path
import csv, re, json, shutil, hashlib, io
from pypdf import PdfReader
from PIL import Image
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
FIG=ROOT/'figures/chapter3'
SELECTED=FIG/'selected'
entries=[
('3.4.6.1','S01','ht1080_profile','ch03_main_ht1080_profile_wgs_cn','HT1080 reference profile'),
('3.4.6.1','F01','ht1080_guide_rtqpcr','ch03_main_ht1080_rtqpcr_three_genes_ideogram','Initial guide testing in HT1080'),
('3.4.6.1','S02','ht1080_all_rtqpcr','ch03_supp_ht1080_rtqpcr_all_genes_ideogram','All transcripts in the initial HT1080 guide assays'),
('3.4.6.1','S03','ddpcr_mixing_check','ch03_main_ddpcr_mixing_validation','ddPCR measurements across defined cell mixtures'),
('3.4.6.1','F02','ht1080_ddpcr','ch03_main_ht1080_ddpcr_nucleofection_rounds','HT1080 guide and nucleofection comparisons'),
('3.4.6.1','F03','ht1080_genomewide','ch03_main_ht1080_lpwgs_genomewide','Genome-wide copy-number profiles in HT1080'),
('3.4.6.1','F04','ht1080_regional_y','ch03_main_ht1080_lpwgs_chrY_atlas','Regional Y depletion and expression in HT1080'),
('3.4.6.2','S01','a549_du145_profiles','ch03_main_a549_du145_profiles','A549 and DU145 reference profiles'),
('3.4.6.2','F01','a549_ddpcr','ch03_main_a549_ddpcr_nucleofection_rounds','A549 guide and nucleofection comparisons'),
('3.4.6.2','F02','a549_regional_y','ch03_main_a549_lpwgs_chrY_atlas','Regional Y depletion and expression in A549'),
('3.4.6.3','F01','du145_ddpcr','ch03_main_du145_ddpcr_assay_comparison','DU145 ddPCR after one and two nucleofections'),
('3.4.6.3','F02','du145_regional_y','ch03_main_du145_lpwgs_chrY_atlas','Regional Y depletion and expression in DU145'),
('3.4.6.3','S01','du145_distal_yp','ch03_supp_du145_lpwgs_PAR1_Yp_zoom','Retained distal Yp in DU145'),
('3.4.6.4','S01','h2452_initial_guides','ch03_supp_h2452_ddpcr_verified_initial_guides','H2452 initial guides: verified ddPCR assay'),
('3.4.6.5','F01','du145_panobinostat','ch03_main_du145_panobinostat_single_round','Panobinostat during single-round Y targeting in DU145'),
('3.4.6.5','S01','du145_y_content','ch03_main_du145_lpwgs_Y_content_summary','Whole-Y and regional DNA content across DU145 treatments'),
('3.4.6.5','S02','du145_all_treatments','ch03_supp_du145_drug_treatment_comparison','Complete DU145 pharmacological comparison'),
('3.4.6.5','F02','du145_pano_grouped','ch03_supp_du145_panobinostat_grouped','DU145 comparison grouped by panobinostat exposure'),
]
configs=[]
for section,number,slug,stem,title in entries:
    code=f'CH3-S{section}-{number}'
    folder=f'CH3_S{section.replace(".","_")}_{number}_{slug}'
    configs.append(dict(section=section,number=number,stem=stem,title=title,code=code,folder=folder,
                        role='main' if number.startswith('F') else 'supp',
                        path=f'figures/chapter3/selected/{folder}/{folder}'))

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def slugify(h):return re.sub(r'[^\w\- ]','',h.lower()).replace(' ','-')
def table_read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t');return reader.fieldnames,list(reader)
def table_write(path,fields,rows):
    out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
    w.writeheader();w.writerows(rows);path.write_bytes(out.getvalue().encode('utf-8'))

rp=ROOT/'Chapter3_results.md';cp=ROOT/'Chapter3_figure_captions.md'
original=rp.read_bytes();results=original.decode('utf-8')
cap_original=cp.read_bytes();captions=cap_original.decode('utf-8')
start=results.index('### 3.4.6 ');end=results.index('### 3.4.7 ',start)
body=results[start:end].replace('\r\n','\n').replace('\r','')
body=body.replace('ch03_main_h2452_ddpcr_initial_guides','ch03_supp_h2452_ddpcr_verified_initial_guides')
# Remove only existing figure display blocks; rebuild main displays at first citation.
body=re.sub(r'^\[!\[.*\n?', '',body,flags=re.M)
body=re.sub(r'^\*\*(?:main|supp) - .*\n?', '',body,flags=re.M)
body=re.sub(r'^\[Figure caption\].*\n?', '',body,flags=re.M)
body=re.sub(r'^\*\[A549 and DU145 baseline profiles\].*\n?',
            '[A549 and DU145 baseline profiles](figures/chapter3/ch03_main_a549_du145_profiles.pdf).\n',body,flags=re.M)
old='We assessed Y depletion by lpWGS and ddPCR targeting AMELY, DDX3Y and KDM5D.'
assert old in body
body=body.replace(old,old[:-1]+' ([treatment figure](figures/chapter3/ch03_main_du145_panobinostat_single_round.pdf)).')
old='These whole-Y summaries therefore showed modest additional reductions alongside persistent regional retention.'
assert old in body
body=body.replace(old,old[:-1]+' ([Y-content summary](figures/chapter3/ch03_main_du145_lpwgs_Y_content_summary.pdf)).')

# Add caption for the audited H2452 batch; preserve legacy combined plot as unselected provenance.
assert '`ch03_supp_h2452_ddpcr_verified_initial_guides`' not in captions
captions+='''\n\n## supp - H2452 initial guides in the verified ddPCR assay

**Figure:** `ch03_supp_h2452_ddpcr_verified_initial_guides` | **Results:** 3.4.6.4

[PDF](figures/chapter3/ch03_supp_h2452_ddpcr_verified_initial_guides.pdf) | [Editable SVG](figures/chapter3/ch03_supp_h2452_ddpcr_verified_initial_guides.svg) | [PNG](figures/chapter3/ch03_supp_h2452_ddpcr_verified_initial_guides.png)

AMELY/AMELX concentration ratios in sgCTRL (source CTRL2), CP, CN and Q53 samples from the verified 3 December 2024 assay. Points show individual DNA measurements; circles identify DNA collected on 21 November and triangles 26 November. Horizontal lines and values show medians. Measurement counts are two for sgCTRL, CP and CN, and one for Q53; these counts do not establish independent nucleofection replication. Ratios were calculated from instrument target and reference concentrations and checked against the exported ratio. Reference-copy-number scaling and inferential testing were not applied. This selected display excludes later batches with unresolved channel orientation and legacy ratio conventions. The initial-guide conclusion remains separate from the later five-guide experiment awaiting data transfer.

[Source and reproduction notes](figures/chapter3/ch03_supp_h2452_ddpcr_verified_initial_guides_notes.md)
'''
captions=captions.replace('The label 2× denotes two successive nucleofections;',
                          'The labels 1× and 2× denote one and two successive nucleofections, respectively;')
captions=captions.replace('The second-nucleofection assignments are author-confirmed.',
                          'The November one-nucleofection and February two-nucleofection assignments are author-confirmed.')
captions=captions.replace('KP4 and H2452 cells were mixed in defined proportions before DNA extraction,',
                          'Y-negative KP4 and Y-positive H2452 cells were mixed in defined proportions before DNA extraction,')

export_checks=[]
blocks=re.split(r'(?m)(?=^## )',captions)
for c in configs:
    assert not (SELECTED/c['folder']).exists(),c['folder']
    for ext in ['pdf','svg','png']:assert (FIG/(c['stem']+'.'+ext)).is_file()
    page=PdfReader(FIG/(c['stem']+'.pdf')).pages[0]
    fonts=[str(v.get_object().get('/BaseFont')) for v in page['/Resources']['/Font'].get_object().values()]
    assert fonts and all('OpenSans' in f for f in fonts),fonts
    for v in page['/Resources']['/Font'].get_object().values():
        font=v.get_object()
        descendants=font.get('/DescendantFonts',[font])
        for child in descendants:
            desc=child.get_object().get('/FontDescriptor')
            if desc: assert any(k in desc.get_object() for k in ['/FontFile','/FontFile2','/FontFile3'])
    with Image.open(FIG/(c['stem']+'.png')) as im:
        assert min(im.info.get('dpi',(0,0)))>590
        dimensions=im.size
    assert '<text' in (FIG/(c['stem']+'.svg')).read_text(encoding='utf-8')
    export_checks.append(dict(code=c['code'],fonts=fonts,png_dimensions=dimensions,
        dimensions_mm=[round(float(v)*25.4/72,2) for v in [page.mediabox.width,page.mediabox.height]],
        sha256={ext:digest(FIG/(c['stem']+'.'+ext)) for ext in ['pdf','svg','png']}))
    hits=[i for i,b in enumerate(blocks) if f'**Figure:** `{c["stem"]}`' in b]
    assert len(hits)==1,(c['stem'],hits)
    i=hits[0];b=blocks[i];heading=b.splitlines()[0]
    heading=re.sub(r'^## (main|supp)',f'## {c["role"]}',heading)
    b=heading+'\n'+b.split('\n',1)[1]
    c['anchor']=slugify(heading.removeprefix('## '))
    b=re.sub(r'\*\*Figure:\*\* `'+re.escape(c['stem'])+r'`[^\n]*',
        f'**Figure:** `{c["stem"]}` | **Working ID:** {c["code"]} | **Results:** {c["section"]}\n\n'
        f'[Selected PDF]({c["path"]}.pdf) | [Selected SVG]({c["path"]}.svg) | [Selected PNG]({c["path"]}.png)',b,count=1)
    b=re.sub(r'^\*\*Status:\*\*[^\n]*\n','',b,flags=re.M)
    blocks[i]=b
    pattern=r'\[([^\]]*)\]\(figures/chapter3/'+re.escape(c['stem'])+r'\.pdf\)'
    def ref(m):
        panel=re.search(r'panel\s+([a-z])',m.group(1))
        return f'[Figure {c["code"]}'+(f', panel {panel.group(1)}' if panel else '')+f']({c["path"]}.pdf)'
    body,count=re.subn(pattern,ref,body)
    assert count>=1,(c['stem'],'not cited')
captions=''.join(blocks)

# Insert main figure displays once, immediately after first citing paragraph.
seen=set();paragraphs=[]
for para in re.split(r'\n\s*\n',body):
    paragraphs.append(para)
    for c in configs:
        if c['role']=='main' and c['code'] not in seen and c['path']+'.pdf' in para:
            paragraphs.append(f'[![Figure {c["code"]}. {c["title"]}.]({c["path"]}.png)]({c["path"]}.pdf)\n\n'
              f'**Figure {c["code"]} (main).** {c["title"]}. [Full caption](Chapter3_figure_captions.md#{c["anchor"]}).')
            seen.add(c['code'])
body='\n\n'.join(p for p in paragraphs if p.strip())+'\n\n'
assert len(seen)==sum(c['role']=='main' for c in configs)
assert not re.search(r'figures/chapter3/(?:previews/)?ch03_',body),'Uncoded reference remains'
for section in sorted({c['section'] for c in configs}):
    for role in ['main','supp']:
        group=[c for c in configs if c['section']==section and c['role']==role]
        actual=sorted(group,key=lambda c:body.index(c['code']))
        assert [c['number'] for c in actual]==sorted(c['number'] for c in group)
updated=results[:start]+body+results[end:]
assert updated[:start]==results[:start] and updated[updated.index('### 3.4.7 '):]==results[end:]
assert not any(x in updated+captions for x in ['\u00c3','\u00c2','\ufffd'])

# Numerical checks: frozen DU145 measurements retain values after the label correction.
backup=json.loads((ROOT/'research/model_development/selection_3_4_6_backup.json').read_text())['archive']
before=pd.read_csv(Path(backup)/'du145_ddpcr_source/main_group_summary.csv')
after=pd.read_csv(FIG/'source_data/du145_ddpcr_assay_comparison/main_group_summary.csv')
pd.testing.assert_frame_equal(before,after)
index_path=SELECTED/'selected_figure_index.tsv';fields,index=table_read(index_path)
assert not any(row['section'].startswith('3.4.6') for row in index)
assert rp.read_bytes()==original and cp.read_bytes()==cap_original,'Concurrent document changes'
for c,check in zip(configs,export_checks):
    dest=SELECTED/c['folder'];dest.mkdir()
    for ext in ['pdf','svg','png']:
        shutil.copy2(FIG/(c['stem']+'.'+ext),dest/(c['folder']+'.'+ext))
        assert digest(dest/(c['folder']+'.'+ext))==check['sha256'][ext]
    (dest/'selection_notes.md').write_text(f'# {c["code"]}: {c["title"]}\n\n'
        f'Selected 23 September 2026 as {c["role"]}, section {c["section"]}.\n\n'
        f'Working source: `{c["stem"]}`. Active caption: [full caption](../../../../Chapter3_figure_captions.md#{c["anchor"]}).\n\n'
        'PDF/SVG/PNG are frozen copies of the reviewed working exports. Reproduction and input records remain with the working figure and source_data. See the section review report in research/model_development/figure_selection_review.\n',encoding='utf-8')
    index.append(dict(working_id=c['code'],chapter='3',section=c['section'],role=c['role'],
        content=c['title'],status='selected',selected_pdf=f'{c["folder"]}/{c["folder"]}.pdf',
        source_stem=c['stem'],sha256_pdf=check['sha256']['pdf']))
# Preserve existing row order, inserting this section before 3.4.7.
newrows=index[-len(configs):];oldrows=index[:-len(configs)]
pos=next((i for i,row in enumerate(oldrows) if tuple(map(int,row['section'].split('.')))>(3,4,6,9)),len(oldrows))
index=oldrows[:pos]+newrows+oldrows[pos:]
table_write(index_path,fields,index)
mp=FIG/'figure_manifest.tsv';mf,mrows=table_read(mp)
for c in configs:
    for ext in ['pdf','svg','png']:
        file=c['stem']+'.'+ext
        found=[row for row in mrows if row['file']==file]
        if not found:
            row={k:'' for k in mf};row.update(file=file,title=c['title'],original_source='Reviewed retained source data');mrows.append(row);found=[row]
        for row in found:row.update(role=c['role'],sha256=digest(FIG/file),status='reviewed and selected 2026-09-23: '+c['code'])
table_write(mp,mf,mrows)
rp.write_bytes(updated.encode('utf-8'));cp.write_bytes(captions.encode('utf-8'))
report=dict(date='2026-09-23',section='3.4.6',figures=configs,exports=export_checks,
    checks=['All 18 export sets exist; embedded Open Sans, editable SVG text, 600-dpi PNG',
            'Visual inspection of all reviewed figure exports via six contact sheets',
            'DU145 ddPCR source medians unchanged after 1x/2x label correction',
            'H2452 uses seven verified concentration-ratio measurements; unresolved combined batches excluded',
            'Main and supplementary codes ascend independently by first mention within each subsection',
            'All 3.4.6 figure links use coded selected copies; other Results sections preserved'],
    limits=['No alignment, segmentation, RNA-seq model or raw-droplet gating rerun',
            'Atlas RNA panels retain explicitly captioned earlier PAR-fix CPM inputs; do not substitute later DE pipelines',
            'Mixing plot is descriptive tracking of recorded ratios; exact ratio convention and fraction direction remain qualified in caption',
            'H2452 five-guide results and DU145 RPS4Y1 testing remain pending; HepG2 is text-only',
            'Selection does not resolve remaining sample/passage links or low-mappability intervals'],archive=backup)
review=ROOT/'research/model_development/figure_selection_review'
(review/'selection_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
(review/'selection_report.md').write_text('# Section 3.4.6 figure selection\n\n'+
    '\n'.join(f'- {c["code"]} ({c["role"]}): {c["title"]}' for c in configs)+
    '\n\n## Checks\n\n'+'\n'.join('- '+x for x in report['checks'])+
    '\n\n## Scope and remaining work\n\n'+'\n'.join('- '+x for x in report['limits'])+'\n',encoding='utf-8')
print(f'Selected {len(configs)} figures: {sum(c["role"]=="main" for c in configs)} main, {sum(c["role"]=="supp" for c in configs)} supplementary.')
