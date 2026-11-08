"""Integrate LOY-frequency panel c and preserve the former panel as S01."""
from pathlib import Path
import csv, hashlib, io, json, re, shutil

ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'figures/chapter3'
DATA=FIG/'source_data/depmap_y_status_dependencies'
revision=json.loads((DATA/'current_revision.json').read_text(encoding='utf-8'))
archive=Path(revision['archive'])
MAIN='CH3-S3.4.1-F01'
SUP='CH3-S3.4.1-S01'
MAIN_STEM='ch03_main_depmap_y_status_dependencies'
SUP_STEM='ch03_supp_depmap_y_call_evidence'
MAIN_SELECTED='CH3_S3_4_1_F01_depmap_y_status_dependencies'
SUP_SELECTED='CH3_S3_4_1_S01_depmap_y_call_evidence'
MAIN_PATH=f'figures/chapter3/selected/{MAIN_SELECTED}/{MAIN_SELECTED}'
SUP_PATH=f'figures/chapter3/selected/{SUP_SELECTED}/{SUP_SELECTED}'
SUP_TITLE='Evidence supporting integrated Y-status calls in DepMap'
SUP_ANCHOR='supp--evidence-supporting-integrated-y-status-calls-in-depmap'
files=['Chapter3_results.md','Results_edit.md','Chapter3_figure_captions.md','Chapter3_source_notes.md']
original={n:(ROOT/n).read_bytes() for n in files}
texts={n:b.decode('utf-8').replace('\r\n','\n') for n,b in original.items()}
assert f'**Working ID:** {SUP}' not in texts['Chapter3_figure_captions.md']
blocks=re.split(r'(?=^## )',texts['Chapter3_figure_captions.md'],flags=re.M)
main_block=next(b for b in blocks if f'**Working ID:** {MAIN} ' in b)
caption_match=re.search(r'^\*\*Y-status classification[^\n]+$',main_block,re.M)
assert caption_match
old_caption=caption_match[0]
old_c=re.search(r'\*\*c,\*\* .*?(?=\*\*d,\*\*)',old_caption).group(0)
new_c=('**c,** Percentage of Y-negative cell lines by cancer lineage (OncoTree lineage), calculated as '
       '100 × Y-negative/(Y-negative + Y-positive). This summary uses non-engineered, male-derived cancer models '
       'before restriction to a particular assay. Of 1,076 eligible models, 910 had definitive integrated calls '
       '(377 Y-negative and 533 Y-positive); 166 set-aside models were excluded from the denominator. '
       'The 17 lineages with at least 10 classified models are shown, comprising 875 models. '
       'Labels give Y-negative/classified counts. ')
new_caption=old_caption.replace(old_c,new_c)
new_main_block=main_block.replace(old_caption,new_caption)
supp_body=('Counts of all 1,143 annotated male models by the evidence supporting their integrated Y-status call. '
           'Blue, red and grey denote Y-positive (588), Y-negative (383) and set-aside (172) models, respectively; '
           'counts are printed for segments containing at least 40 models. Expression determines the integrated call, '
           'DNA supplies a call where RNA is unavailable, and definite disagreement sets a model aside. '
           'Intermediate RNA profiles, ambiguous DNA without RNA, and models lacking both measurements are also set aside. '
           'DNA copy ratios derive from whole-genome sequencing where available and otherwise from exome sequencing. '
           'This inventory precedes the cancer-model and engineering exclusions used for the lineage frequencies in '
           f'[Figure {MAIN}, panel c]({MAIN_PATH}.pdf).')
supp_block=(f'## supp — {SUP_TITLE}\n\n'
            f'**Figure:** `{SUP_STEM}` | **Working ID:** {SUP} | **Results:** 3.4.1 | **Title:** {SUP_TITLE}.\n\n'
            f'[Selected PDF]({SUP_PATH}.pdf) | [Selected SVG]({SUP_PATH}.svg) | [Selected PNG]({SUP_PATH}.png)\n\n'
            f'[PDF](figures/chapter3/{SUP_STEM}.pdf) | [Editable SVG](figures/chapter3/{SUP_STEM}.svg) | [PNG](figures/chapter3/{SUP_STEM}.png)\n\n'
            f'**{SUP_TITLE}.** {supp_body}\n\n'
            f'[Source and reproduction notes](figures/chapter3/selected/{SUP_SELECTED}/selection_notes.md)\n\n')
texts['Chapter3_figure_captions.md']=texts['Chapter3_figure_captions.md'].replace(main_block,new_main_block+supp_block,1)

for name in ['Chapter3_results.md','Results_edit.md']:
    t=texts[name]
    old_ref=f'([Figure {MAIN}, panels a–c]({MAIN_PATH}.pdf)).'
    new_ref=(f'([Figure {MAIN}, panels a and b]({MAIN_PATH}.pdf); [Figure {SUP}]({SUP_PATH}.pdf)). '
             f'Among eligible cancer models, Y-negative frequencies varied across cancer lineages '
             f'([Figure {MAIN}, panel c]({MAIN_PATH}.pdf)).')
    assert t.count(old_ref)==1,name
    t=t.replace(old_ref,new_ref,1)
    if name=='Results_edit.md':
        old_body=re.sub(r'^\*\*[^*]+\*\*\s*','',old_caption,count=1)
        new_body=re.sub(r'^\*\*[^*]+\*\*\s*','',new_caption,count=1)
        assert old_body in t
        t=t.replace(old_body,new_body,1)
        insert=(f'\n\n---\n\n**Figure {SUP} (supp). {SUP_TITLE}.**\n\n'
                f'[![{SUP}]({ROOT.as_posix()}/{SUP_PATH}.png)]({ROOT.as_posix()}/{SUP_PATH}.pdf)\n\n'
                f'{supp_body}\n\n---')
        pattern=rf'(---\n\n\*\*Figure {re.escape(MAIN)} \(main\).*?\n\n---)'
        t,n=re.subn(pattern,lambda m:m[0]+insert,t,count=1,flags=re.S)
        assert n==1
    else:
        insert=(f'\n\n[![Figure {SUP}. {SUP_TITLE}.]({SUP_PATH}.png)]({SUP_PATH}.pdf)\n\n'
                f'**Figure {SUP} (supp).** {SUP_TITLE}. [Full caption](Chapter3_figure_captions.md#{SUP_ANCHOR}).')
        pattern=rf'^\*\*Figure {re.escape(MAIN)} \(main\).*?$'
        t,n=re.subn(pattern,lambda m:m[0]+insert,t,count=1,flags=re.M)
        assert n==1
    texts[name]=t

texts['Chapter3_source_notes.md']+=(
    '\n\n## DepMap LOY frequency panel — 5 October 2026\n\n'
    f'At the author’s request, {MAIN}c now shows the percentage of Y-negative cell lines by OncoTree lineage. '
    'The denominator contains only definitive Y-negative and Y-positive calls among non-engineered, annotated male cancer models, '
    'before restriction by assay availability. The frozen `cohort.csv` contains 1,076 eligible models, of which 910 are classified '
    '(377 negative, 533 positive). Seventeen lineages meet the display threshold of ten classified models, accounting for 875 models. '
    'All 24 lineage counts, including those below the display threshold, are retained in the source table. '
    f'The former call-evidence panel is now {SUP}; its measurements and counts are unchanged. '
    'No dependency model was refitted. '
    '[Source tables and calculation](figures/chapter3/source_data/depmap_y_status_dependencies/README.md).\n')

for stem,selected in [(MAIN_STEM,MAIN_SELECTED),(SUP_STEM,SUP_SELECTED)]:
    folder=FIG/'selected'/selected
    folder.mkdir(exist_ok=True)
    for ext in ['pdf','svg','png']:
        shutil.copy2(FIG/f'{stem}.{ext}',folder/f'{selected}.{ext}')
for name in files:
    assert (ROOT/name).read_bytes()==original[name],f'Intervening edit: {name}'
    out=texts[name]
    if b'\r\n' in original[name]:out=out.replace('\n','\r\n')
    (ROOT/name).write_bytes(out.encode('utf-8'))

for selected,title,code in [(MAIN_SELECTED,'Combined DepMap main figure',MAIN),(SUP_SELECTED,'DepMap call-evidence supplement',SUP)]:
    (FIG/'selected'/selected/'selection_notes.md').write_text(
        f'# {title}\n\nActive working identifier: {code}. Updated 5 October 2026. '
        f'The main figure shows Y-status profiles (a–b), Y-negative percentages by lineage (c), and dependency results (d–e). '
        f'The former main panel c is retained as {SUP}. Main panel c uses classified eligible cancer models; '
        'the call-evidence supplement uses all annotated male models.\n\n'
        '- [Active captions](../../../../Chapter3_figure_captions.md)\n'
        f'- [Reproduction script](../../code/{MAIN_STEM}.py)\n'
        '- [Source tables and provenance](../../source_data/depmap_y_status_dependencies/README.md)\n\n'
        f'Previous exports and captions: `{archive.relative_to(ROOT).as_posix()}`. '
        'Dependency statistics are unchanged. No plot titles are used.\n',encoding='utf-8')
(DATA/'README.md').write_text(
    '# Source data — DepMap Y-status, lineage frequencies and dependencies\n\n'
    f'Main figure: [{MAIN}](../../{MAIN_STEM}.pdf). Call-evidence supplement: [{SUP}](../../{SUP_STEM}.pdf). '
    '[Active captions](../../../../Chapter3_figure_captions.md). '
    f'[Reproduction script](../../code/{MAIN_STEM}.py).\n\n'
    'The original five `fig*.csv` files are exact copies of the saved 5 October 2026 analysis exports. '
    'The call-tier table now supplies supplementary S01. Main panels a, b, d and e retain the original profiles and dependency results.\n\n'
    '`lineage_prevalence_models.csv` records the 1,076 eligible model identifiers, donor identifiers, lineages and integrated Y calls '
    'from the frozen analysis cohort. `lineage_prevalence_all.csv` contains all 24 lineages; '
    '`lineage_prevalence_displayed.csv` contains the 17 lineages with at least 10 definitive calls. '
    'LOY% = 100 × Y-negative/(Y-negative + Y-positive); set-aside models are excluded from that denominator. '
    'Counts are model counts, not unique-donor counts. No assay-availability restriction is applied. '
    'The display threshold is for readability and is not a statistical significance criterion.\n\n'
    '`provenance.json` records the exact source paths and hashes, sample counts, mapping and display rule. '
    'Scientific estimates in the dependency panels were not refitted. Outputs use the bundled Open Sans font, '
    'editable vector PDF/SVG and 300-dpi PNG.\n',encoding='utf-8')

def read_tsv(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t');return reader.fieldnames,list(reader)
def write_tsv(path,columns,rows):
    out=io.StringIO(newline='');w=csv.DictWriter(out,fieldnames=columns,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    path.write_bytes(out.getvalue().encode('utf-8'))
index=FIG/'selected/selected_figure_index.tsv'
columns,rows=read_tsv(index)
assert sum(r['working_id']==MAIN and 'depmap_y_status_dependencies' in r['selected_pdf'] for r in rows)==1
for r in rows:
    if r['working_id']==MAIN and 'depmap_y_status_dependencies' in r['selected_pdf']:
        r['sha256_pdf']=hashlib.sha256((ROOT/f'{MAIN_PATH}.pdf').read_bytes()).hexdigest()
rows.append(dict(working_id=SUP,chapter='3',section='3.4.1',role='supp',content=SUP_TITLE,status='selected',
                 selected_pdf=f'{SUP_SELECTED}/{SUP_SELECTED}.pdf',source_stem=SUP_STEM,
                 sha256_pdf=hashlib.sha256((ROOT/f'{SUP_PATH}.pdf').read_bytes()).hexdigest()))
write_tsv(index,columns,rows)
manifest=FIG/'figure_manifest.tsv'
columns,rows=read_tsv(manifest)
for stem,role,title in [(MAIN_STEM,'main','Y-status classification and associated gene dependencies in DepMap'),(SUP_STEM,'supp',SUP_TITLE)]:
    for ext in ['pdf','svg','png']:
        filename=f'{stem}.{ext}'
        matching=[r for r in rows if r['file']==filename]
        if matching:r=matching[0]
        else:r={k:'' for k in columns};r['file']=filename;rows.append(r)
        r.update(role=role,title=f'{role} - {title}',status='selected; updated 2026-10-05',
                 original_source=f'figures/chapter3/code/{MAIN_STEM}.py',
                 sha256=hashlib.sha256((FIG/filename).read_bytes()).hexdigest())
write_tsv(manifest,columns,rows)

review_manifest=ROOT/'exports/Results_edit_manifest.json'
m=json.loads(review_manifest.read_text(encoding='utf-8'))
position=next(i for i,f in enumerate(m['figures']) if f['code']==MAIN)
m['figures'].insert(position+1,dict(ref=f'{SUP_PATH}.pdf',code=SUP,role='supp'))
m['source_sha256']=hashlib.sha256((ROOT/'Chapter3_results.md').read_bytes()).hexdigest()
m['review_copy_sha256']=hashlib.sha256((ROOT/'Results_edit.md').read_bytes()).hexdigest()
review_manifest.write_text(json.dumps(m,indent=2),encoding='utf-8')

for name in ['Chapter3_results.md','Results_edit.md']:
    out=(ROOT/name).read_text(encoding='utf-8')
    assert out.count(f'**Figure {SUP} (supp)')==1
    assert out.split('### 3.4.2 ',1)[1] == original[name].decode('utf-8').replace('\r\n','\n').split('### 3.4.2 ',1)[1]
    for destination in re.findall(r'\]\(([^\n)]+)\)',out):
        d=destination.strip('<>').split('#',1)[0]
        if re.match(r'(?:https?:|mailto:)',d) or not d:continue
        p=Path(d)
        if not p.is_absolute():p=ROOT/p
        assert p.exists(),(name,d)
print(json.dumps({'main':MAIN,'supplement':SUP,'lineages_displayed':17,'classified_models_displayed':875,
                  'all_lineages_in_source_table':24,'later_results_preserved':True,'archive':str(archive)},indent=2))
